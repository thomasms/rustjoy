//! Port of `endfio.f` (PREPRO 2026-1): the ENDF record routines that all
//! PREPRO codes share, the number editing (`OUT9`, `OUT9G`, `OUT10`,
//! `OUT7`, `IN9`), the `INCORE` rounding, `ZAHOL`, `TERPIT`, `MAXIE` and
//! the timer.
//!
//! The COMMON blocks `/ENDFIO/`, `/HEADER/`, `/LEADER/`, `/LASTE/`,
//! `/FLAGS/`, `/COPC/`, `/COPI/`, `/NAMECODE/` are fields of [`Endf`].

#![allow(dead_code)]

use crate::ffmt;
use crate::fio::{self, Rd, Units, L};
use std::sync::OnceLock;

// ---------------------------------------------------------------------
// Tables built as the Fortran builds them on the first call.

struct Tabs {
    /// TENS(-99:99), index i+99.
    tens: [f64; 199],
    /// ROUNDER(-99:99) of OUT9/OUT9G/OUT10/OUT7.
    rounder: [f64; 199],
}

fn tabs() -> &'static Tabs {
    static T: OnceLock<Tabs> = OnceLock::new();
    T.get_or_init(|| {
        let mut tens = [0.0f64; 199];
        let mut rounder = [0.0f64; 199];
        tens[99] = 1.0;
        for i in 1..=99usize {
            tens[99 + i] = tens[99 + i - 1] * 10.0;
            tens[99 - i] = tens[99 - i + 1] / 10.0;
            rounder[99 + i] = 5.001e-05;
        }
        for i in -99..=-10i32 {
            rounder[(i + 99) as usize] = 5.001e-05;
        }
        for i in -9..=-4i32 {
            rounder[(i + 99) as usize] = 5.001e-06;
        }
        rounder[99 - 3] = 5.001e-07;
        rounder[99 - 2] = 5.001e-08;
        rounder[99 - 1] = 5.001e-09;
        rounder[99] = 5.001e-09;
        for i in 1..=8usize {
            rounder[99 + i] = 5.001e-09;
        }
        rounder[99 + 9] = 5.001e-06;
        Tabs { tens, rounder }
    })
}

/// TENS(i). Outside -99..99 the Fortran reads outside the array (only
/// INCORE with numbers below 1e-90 can do that); here it is 10**i.
#[inline]
pub fn tens(i: i32) -> f64 {
    if (-99..=99).contains(&i) {
        tabs().tens[(i + 99) as usize]
    } else {
        10f64.powi(i)
    }
}

#[inline]
fn rounder(i: i32) -> f64 {
    tabs().rounder[(i + 99) as usize]
}

const DIGITS: &[u8; 10] = b"0123456789";
const ZEROH: [u8; 11] = *b" 0.0       ";

/// Fortran INT() of a REAL*8 to INTEGER*4 (truncation; gfortran gives the
/// conversion result of cvttsd2si, which saturates to i32::MIN).
#[inline]
pub fn fint(x: f64) -> i32 {
    if x.is_nan() || x >= 2147483648.0 || x < -2147483648.0 {
        i32::MIN
    } else {
        x as i32
    }
}

/// Fortran conversion of a REAL*8 to INTEGER*8.
#[inline]
pub fn fint8(x: f64) -> i64 {
    if x.is_nan() || x >= 9.223372036854775808e18 || x < -9.223372036854775808e18 {
        i64::MIN
    } else {
        x as i64
    }
}

/// The mantissa digits of the E forms (labels 180..250 of OUT9).
fn put_e(field: &mut [u8], mut inmant: i64, mut iexp: i32, istart: usize, id: usize, kdot: usize) {
    let mut iexps = id + 1;
    let mut ii = id;
    for _ in istart..=id {
        if ii != kdot {
            let innext = inmant / 10;
            let i3 = inmant - 10 * innext;
            field[ii - 1] = DIGITS[i3 as usize];
            inmant = innext;
        }
        ii -= 1;
    }
    field[iexps - 1] = b'E';
    iexps += 1;
    if iexp >= 0 {
        field[iexps - 1] = b'+';
    } else {
        iexp = -iexp;
        field[iexps - 1] = b'-';
    }
    if iexp >= 10 {
        let kexp = iexp / 10;
        field[9] = DIGITS[kexp as usize];
        iexp %= 10;
    }
    field[10] = DIGITS[iexp as usize];
}

/// E form that leaves column 1 for the sign (labels 180-190 of OUT9).
fn e_col3(field: &mut [u8], zn: f64, mut iexp: i32) {
    field[2] = b'.';
    let (kdot, istart) = (3usize, 2usize);
    let id;
    let mut inmant;
    if iexp.abs() < 10 {
        inmant = fint8(1.0e5 * zn);
        if inmant < 1_000_000 {
            put_e(field, inmant, iexp, istart, 8, kdot);
            return;
        }
        inmant = 100_000;
        iexp += 1;
        if iexp.abs() < 10 {
            put_e(field, inmant, iexp, istart, 8, kdot);
            return;
        }
    }
    id = 7;
    inmant = fint8(1.0e4 * zn);
    if inmant >= 100_000 {
        inmant = 10_000;
        iexp += 1;
    }
    put_e(field, inmant, iexp, istart, id, kdot);
}

/// E form that uses column 1 (positive numbers in OUT10 and OUT7).
fn e_col2(field: &mut [u8], zn: f64, mut iexp: i32) {
    field[1] = b'.';
    let (kdot, istart) = (2usize, 1usize);
    let mut inmant;
    if iexp.abs() < 10 {
        inmant = fint8(1.0e6 * zn);
        if inmant < 10_000_000 {
            put_e(field, inmant, iexp, istart, 8, kdot);
            return;
        }
        inmant = 1_000_000;
        iexp += 1;
        if iexp.abs() < 10 {
            put_e(field, inmant, iexp, istart, 8, kdot);
            return;
        }
    }
    inmant = fint8(1.0e5 * zn);
    if inmant >= 1_000_000 {
        inmant = 100_000;
        iexp += 1;
    }
    put_e(field, inmant, iexp, istart, 7, kdot);
}

/// Shared start of OUT9/OUT9G/OUT7: sign, exponent, rounding. Returns
/// None for the zero field, else (z, zn, iexp).
fn out_start(zin: f64, field: &mut [u8], extra: bool) -> Option<(f64, f64, i32)> {
    let a = zin.abs();
    if a <= 0.0 || a < 1.0e-99 || a > 1.0e99 || zin.is_nan() {
        field[..11].copy_from_slice(&ZEROH);
        return None;
    }
    let mut z;
    if zin < 0.0 {
        field[0] = b'-';
        z = -zin;
    } else {
        field[0] = b' ';
        z = zin;
    }
    let mut iexp = fint(z.log10());
    if z < 1.0 {
        iexp -= 1;
    }
    if iexp.abs() > 97 {
        field[..11].copy_from_slice(&ZEROH);
        return None;
    }
    let mut zn = if zin < 0.0 {
        z * tens(-iexp) + rounder(iexp)
    } else if extra && (z <= 1.0e-3 || z >= 1.0e9) {
        z * tens(-iexp) + 0.1 * rounder(iexp)
    } else {
        z * tens(-iexp) + rounder(iexp)
    };
    if zn == 1.0 {
    } else if zn > 1.0 {
        if zn >= 10.0 {
            iexp += 1;
            zn /= 10.0;
        }
    } else {
        iexp -= 1;
    }
    z = zn * tens(iexp);
    Some((z, zn, iexp))
}

/// The 9 digit F form of OUT9/OUT9G (labels 130-150).
fn f9(field: &mut [u8], z: f64, mut iexp: i32) {
    let mut ipower = 8 - iexp;
    if iexp < 0 {
        ipower = 8;
    }
    let mut inmant = fint8(z * tens(ipower));
    if inmant >= 1_000_000_000 {
        inmant = 100_000_000;
        iexp += 1;
    }
    let mut idot = 3 + iexp;
    if idot <= 2 {
        idot = 2;
        inmant = fint8(z * 10.0 * tens(ipower));
    }
    field[(idot - 1) as usize] = b'.';
    let mut ii = 11i32;
    for _ in 2..=11 {
        if ii != idot {
            let innext = inmant / 10;
            let i3 = inmant - 10 * innext;
            field[(ii - 1) as usize] = DIGITS[i3 as usize];
            inmant = innext;
        }
        ii -= 1;
    }
}

/// `OUT9`: 11 column form with up to 9 digits.
pub fn out9(zin: f64, field: &mut [u8]) {
    let Some((z, zn, iexp)) = out_start(zin, field, true) else { return };
    if z <= 1.0e-3 || z >= 1.0e9 {
        // 2026/2/21 - always leave column 1 blank
        e_col3(field, zn, iexp);
        return;
    }
    f9(field, z, iexp);
}

/// `OUT9G`: as OUT9, column 1 is always '-' or blank.
pub fn out9g(zin: f64, field: &mut [u8]) {
    let Some((z, zn, iexp)) = out_start(zin, field, false) else { return };
    if z <= 1.0e-3 || z >= 1.0e9 {
        e_col3(field, zn, iexp);
        return;
    }
    f9(field, z, iexp);
}

/// `OUT7`: always the E form, column 1 used for positive numbers.
pub fn out7(zin: f64, field: &mut [u8]) {
    let Some((_z, zn, iexp)) = out_start(zin, field, true) else { return };
    if zin <= 0.0 {
        e_col3(field, zn, iexp);
    } else {
        e_col2(field, zn, iexp);
    }
}

/// `OUT10`: 11 column form with up to 10 digits (positive numbers).
pub fn out10(zin: f64, field: &mut [u8]) {
    let a = zin.abs();
    if a <= 0.0 || a < 1.0e-99 || a > 1.0e99 || zin.is_nan() {
        field[..11].copy_from_slice(&ZEROH);
        return;
    }
    if zin < 0.0 {
        out9(zin, field);
        return;
    }
    field[0] = b' ';
    let mut z = zin;
    let mut iexp = fint(z.log10());
    if iexp.abs() > 97 {
        field[..11].copy_from_slice(&ZEROH);
        return;
    }
    if z < 1.0 {
        iexp -= 1;
    }
    let mut zn = z * tens(-iexp) + 0.1 * rounder(iexp);
    if zn == 1.0 {
    } else if zn > 1.0 {
        if zn >= 10.0 {
            iexp += 1;
            zn /= 10.0;
        }
    } else {
        iexp -= 1;
    }
    z = zn * tens(iexp);
    if z <= 1.0e-3 || z >= 1.0e9 {
        e_col2(field, zn, iexp);
        return;
    }
    let mut ipower = 9 - iexp;
    if iexp < 0 {
        ipower = 9;
    }
    let mut inmant = fint8(z * tens(ipower));
    if inmant >= 10_000_000_000 {
        inmant = 1_000_000_000;
        iexp += 1;
    }
    let mut idot = 2 + iexp;
    if idot <= 1 {
        idot = 1;
        inmant = fint8(z * 10.0 * tens(ipower));
    }
    field[(idot - 1) as usize] = b'.';
    let mut ii = 11i32;
    for _ in 1..=11 {
        if ii != idot {
            let innext = inmant / 10;
            let i3 = inmant - 10 * innext;
            field[(ii - 1) as usize] = DIGITS[i3 as usize];
            inmant = innext;
        }
        ii -= 1;
    }
}

/// `OUT9` into a new 11 byte field (callers that start from blanks).
pub fn f9s(z: f64) -> [u8; 11] {
    let mut f = [b' '; 11];
    out9(z, &mut f);
    f
}

pub fn f9g(z: f64) -> [u8; 11] {
    let mut f = [b' '; 11];
    out9g(z, &mut f);
    f
}

pub fn f10(z: f64) -> [u8; 11] {
    let mut f = [b' '; 11];
    out10(z, &mut f);
    f
}

/// `IN9`: an 11 column field to a number. Only the first 11 bytes are
/// used (the Fortran can look at the byte after the field after a final
/// 'E', but the result does not depend on it).
pub fn in9(field: &[u8]) -> f64 {
    let fc = |j: usize| -> u8 {
        // j is 1-based
        if j >= 1 && j <= field.len() && j <= 11 {
            field[j - 1]
        } else {
            b' '
        }
    };
    let mut i = 1usize;
    while i <= 11 {
        if fc(i) != b' ' {
            break;
        }
        i += 1;
    }
    if i > 11 {
        return 0.0;
    }
    let mut sign = 1.0f64;
    if fc(i) == b'-' {
        sign = -1.0;
        i += 1;
    } else if fc(i) == b'+' {
        i += 1;
    }
    let mut x9in = 0.0f64;
    let mut ipt: i32 = -20;
    let mut imzero = 0;
    let mut j = i;
    // exponent state
    enum Next {
        Done,
        Exp(usize),     // label 160 with J
        ExpPlus(usize), // label 170
        ExpMinus(usize), // label 190
        Err(usize),
    }
    let mut next = Next::Done;
    'scan: {
        while j <= 11 {
            let c = fc(j);
            if c.is_ascii_digit() {
                let k = (c - b'0') as i32;
                if ipt >= 0 {
                    if k == 0 {
                        imzero += 1;
                        j += 1;
                        continue;
                    }
                    while imzero > 0 {
                        x9in *= 10.0;
                        ipt += 1;
                        imzero -= 1;
                    }
                }
                x9in = 10.0 * x9in + k as f64;
                ipt += 1;
                j += 1;
                continue;
            }
            if c == b'.' {
                ipt = 0;
                j += 1;
                continue;
            }
            if c == b' ' {
                // label 130: scan the rest for an exponent
                let i2 = j + 1;
                if i2 > 11 {
                    break 'scan;
                }
                let mut jj = i2;
                while jj <= 11 {
                    let c2 = fc(jj);
                    if c2 == b' ' {
                        jj += 1;
                        continue;
                    }
                    next = match c2 {
                        b'e' | b'E' | b'd' | b'D' => Next::Exp(jj),
                        b'-' => Next::ExpMinus(jj),
                        b'+' => Next::ExpPlus(jj),
                        _ => Next::Err(jj),
                    };
                    break 'scan;
                }
                break 'scan;
            }
            next = match c {
                b'e' | b'E' | b'd' | b'D' => Next::Exp(j),
                b'-' => Next::ExpMinus(j),
                b'+' => Next::ExpPlus(j),
                _ => Next::Err(j),
            };
            break 'scan;
        }
    }
    let (jstart, isign) = match next {
        Next::Done => {
            let mut e = x9in;
            if ipt > 0 {
                e /= tens(ipt);
            }
            return sign * e;
        }
        Next::Err(jj) => return in9_err(field, jj),
        Next::Exp(jj) => {
            let j2 = jj + 1;
            let c = fc(j2);
            if c == b'-' {
                (j2 + 1, -1)
            } else if c == b'+' {
                (j2 + 1, 1)
            } else {
                (j2, 1)
            }
        }
        Next::ExpPlus(jj) => (jj + 1, 1),
        Next::ExpMinus(jj) => (jj + 1, -1),
    };
    let mut iexp: i32 = 0;
    let mut ii = jstart;
    while ii <= 11 {
        let c = fc(ii);
        if c == b' ' {
            ii += 1;
            continue;
        }
        if c.is_ascii_digit() {
            iexp = iexp.wrapping_mul(10).wrapping_add((c - b'0') as i32);
            ii += 1;
            continue;
        }
        return in9_err(field, jstart);
    }
    let e = x9in;
    iexp *= isign;
    if ipt > 0 {
        iexp -= ipt;
    }
    if iexp.abs() <= 99 {
        sign * e * tens(iexp)
    } else {
        0.0
    }
}

fn in9_err(field: &[u8], j: usize) -> f64 {
    let mut mess = [b' '; 11];
    if (1..=11).contains(&j) {
        mess[j - 1] = b'*';
    }
    let mut f = [b' '; 11];
    for k in 0..11.min(field.len()) {
        f[k] = field[k];
    }
    fio::out(L::new().s(" ").b(&f));
    fio::out(L::new().s(" ").b(&mess));
    fio::out(" ERROR in Input Data...Translated as 0.");
    0.0
}

/// `INT11`: integer right adjusted in 11 columns.
pub fn int11(iout: i32) -> [u8; 11] {
    let mut f = [b' '; 11];
    if iout.unsigned_abs() == 0 {
        f[10] = b'0';
        return f;
    }
    let mut ii = iout.wrapping_abs() as i64;
    let mut id = 11usize;
    loop {
        let kk = (ii % 10) as usize;
        f[id - 1] = DIGITS[kk];
        ii /= 10;
        if ii <= 0 {
            break;
        }
        if id == 1 {
            id = 2;
            break;
        }
        id -= 1;
    }
    if iout < 0 {
        id -= 1;
        f[id - 1] = b'-';
    }
    f
}

// ---------------------------------------------------------------------
// INCORE rounding

/// `INCORE10`.
pub fn incore10(zin: &mut f64) {
    if *zin == 0.0 {
        return;
    }
    if !(*zin > 0.0) {
        incore9(zin);
        return;
    }
    let mut iexp = fint(zin.log10());
    if *zin < 1.0 {
        iexp -= 1;
    }
    if iexp.abs() > 99 {
        return;
    }
    let mut znorm = *zin * tens(-iexp) + 5.0e-10;
    if znorm >= 10.0 {
        iexp += 1;
        znorm = 1.0;
    }
    let z10 = znorm * tens(9);
    let i10 = fint8(z10);
    let z10 = i10 as f64;
    *zin = z10 * tens(iexp - 9);
}

/// `INCORE9`.
pub fn incore9(zin: &mut f64) {
    if *zin == 0.0 {
        return;
    }
    let pzin = zin.abs();
    let mut iexp = fint(pzin.log10());
    if pzin < 1.0 {
        iexp -= 1;
    }
    if iexp.abs() > 99 {
        return;
    }
    let mut znorm = pzin * tens(-iexp) + 5.0e-09;
    if znorm >= 10.0 {
        iexp += 1;
        znorm = 1.0;
    }
    let z9 = znorm * tens(8);
    let i9 = fint8(z9);
    let mut z9 = i9 as f64;
    if *zin < 0.0 {
        z9 = -z9;
    }
    *zin = z9 * tens(iexp - 8);
}

/// `INCORE7`.
pub fn incore7(zin: &mut f64) {
    if *zin == 0.0 {
        return;
    }
    let pzin = zin.abs();
    let mut iexp = fint(pzin.log10());
    if pzin < 1.0 {
        iexp -= 1;
    }
    if iexp.abs() > 99 {
        return;
    }
    let mut znorm = pzin * tens(-iexp) + 5.0e-07;
    if znorm >= 10.0 {
        iexp += 1;
        znorm = 1.0;
    }
    let z7 = znorm * tens(6);
    let i7 = fint8(z7);
    let mut z7 = i7 as f64;
    if *zin < 0.0 {
        z7 = -z7;
    }
    *zin = z7 * tens(iexp - 6);
}

// ---------------------------------------------------------------------
// ZAHOL, ZAHOLM, ZANAME

const ZATAB: &[u8; 236] = b"H HeLiBeB C N O F NeNaMgAlSiP S ClArK CaScTiV CrMnFeCoNiCuZnGaGeAsSeBrKrRbSrY ZrNbMoTcRuRhPdAgCdInSnSbTeI XeCsBaLaCePrNdPmSmEuGdTbDyHoErTmYbLuHfTaW ReOsIrPtAuHgTlPbBiPoAtRnFrRaAcThPaU NpPuAmCmBkCfEsFmMdNoLrRfDbSgBhHsMtDsRgCnNhFlMcLvTsOg";

/// `ZAHOL`: ZA as 10 characters ZZZ-SS-AAA.
pub fn zahol(za: i32) -> [u8; 10] {
    let mut f = [b' '; 10];
    if za == 1 {
        return *b" Neutron  ";
    }
    if za == 99120 || za == 99125 {
        return *b"Fiss.Prod.";
    }
    if za == 0 {
        return *b" Photon   ";
    }
    let mut z = za / 1000;
    let mut a = za - 1000 * z;
    f[3] = b'-';
    f[6] = b'-';
    let mut skip_z = false;
    if z > 0 && z <= 118 {
        f[4] = ZATAB[(2 * (z - 1)) as usize];
        f[5] = ZATAB[(2 * (z - 1) + 1) as usize];
    } else {
        f[4] = b'?';
        f[5] = b'?';
        if z < 0 || z > 999 {
            skip_z = true;
        }
    }
    if !skip_z {
        let mut ii = 3usize;
        for _ in 0..3 {
            let nextz = z / 10;
            let kz = z - 10 * nextz;
            f[ii - 1] = DIGITS[kz as usize];
            z = nextz;
            if z <= 0 {
                break;
            }
            ii -= 1;
        }
    }
    if a <= 0 {
        f[7] = b'N';
        f[8] = b'a';
        f[9] = b't';
        return f;
    }
    let mut idiv = 100;
    let mut imon = 0;
    let mut ii = 7usize;
    for _ in 0..3 {
        let ia = a / idiv;
        if !(ia == 0 && imon == 0) {
            imon = 1;
            ii += 1;
            f[ii - 1] = DIGITS[ia as usize];
        }
        a -= idiv * ia;
        idiv /= 10;
    }
    f
}

/// `ZAHOLM`: ZAHOL plus the isomer flag, 12 characters.
pub fn zaholm(za: i32, liso: i32) -> [u8; 12] {
    let mut f = [b' '; 12];
    f[..10].copy_from_slice(&zahol(za));
    if liso <= 0 {
        return f;
    }
    let mut i = 10usize;
    while i >= 1 {
        if f[i - 1] != b' ' {
            break;
        }
        i -= 1;
    }
    if i == 0 {
        return f;
    }
    i += 1;
    f[i - 1] = b'm';
    if liso == 2 {
        f[i - 1] = b'n';
    }
    if liso > 9 {
        f[i] = b'?';
    }
    if (3..=9).contains(&liso) {
        f[i] = DIGITS[liso as usize];
    }
    f
}

/// `ZANAME`: characters 1-8 of NAME1 become ZAzzzaaa.
pub fn zaname(iza: i32, name1: &mut [u8]) {
    name1[0] = b'Z';
    name1[1] = b'A';
    let mut kkza = iza;
    for kk in (3..=8).rev() {
        let idigit = kkza % 10;
        name1[kk - 1] = DIGITS[idigit.rem_euclid(10) as usize];
        kkza /= 10;
    }
}

// ---------------------------------------------------------------------
// TERPIT, HIGHEND

/// `TERPIT`: ENDF interpolation laws 1 to 6.
pub fn terpit(x: f64, x1: f64, x2: f64, y1: f64, y2: f64, interp: i32) -> f64 {
    if x == x1 {
        return y1;
    }
    if x == x2 {
        return y2;
    }
    if x1 == x2 || y1 == y2 {
        return y1;
    }
    let lin = || {
        let wt2 = (x - x1) / (x2 - x1);
        let wt1 = 1.0 - wt2;
        wt2 * y2 + wt1 * y1
    };
    match interp {
        1 => y1,
        2 => lin(),
        3 => {
            if x <= 0.0 || x1 <= 0.0 || x2 <= 0.0 {
                return lin();
            }
            let wt2 = (x / x1).ln() / (x2 / x1).ln();
            let wt1 = 1.0 - wt2;
            wt2 * y2 + wt1 * y1
        }
        4 => {
            if y1 <= 0.0 || y2 <= 0.0 {
                return lin();
            }
            let wt2 = (x - x1) / (x2 - x1);
            let wt1 = 1.0 - wt2;
            (wt2 * y2.ln() + wt1 * y1.ln()).exp()
        }
        5 => {
            if x <= 0.0 || x1 <= 0.0 || x2 <= 0.0 {
                return lin();
            }
            if y1 <= 0.0 || y2 <= 0.0 {
                return lin();
            }
            let wt2 = (x / x1).ln() / (x2 / x1).ln();
            let wt1 = 1.0 - wt2;
            (wt2 * y2.ln() + wt1 * y1.ln()).exp()
        }
        6 => {
            if x <= 0.0 || x1 <= 0.0 || x2 <= 0.0 {
                return lin();
            }
            if y1 <= 0.0 || y2 <= 0.0 {
                return lin();
            }
            let wt2 = (1.0 / x.sqrt() - 1.0 / x1.sqrt()) / (1.0 / x2.sqrt() - 1.0 / x1.sqrt());
            let wt1 = 1.0 - wt2;
            (wt2 * (x2 * y2).ln() + wt1 * (x1 * y1).ln()).exp() / x
        }
        // computed GO TO out of range falls through to label 10
        _ => y1,
    }
}

/// `HIGHEND`: make a spectrum end with Y = 0. Arrays are 1-based in the
/// Fortran; here x[0] is X(1).
pub fn highend(x: &mut [f64], y: &mut [f64], nbt: &mut [i32], int: &[i32], n1: usize, n2: &mut usize) {
    if y[*n2 - 1] == 0.0 {
        return;
    }
    if int[n1 - 1] == 1 {
        if x[*n2 - 1] != x[*n2 - 2] {
            y[*n2 - 1] = 0.0;
        } else {
            *n2 -= 1;
            y[*n2 - 1] = 0.0;
            nbt[n1 - 1] = *n2 as i32;
        }
        return;
    }
    *n2 += 1;
    x[*n2 - 1] = x[*n2 - 2];
    y[*n2 - 1] = 0.0;
    nbt[n1 - 1] = *n2 as i32;
}

// ---------------------------------------------------------------------
// The ENDF I/O state.

/// COMMON/HEADER/ or /LEADER/.
#[derive(Clone, Copy, Default, Debug)]
pub struct Head {
    pub c1: f64,
    pub c2: f64,
    pub l1: i32,
    pub l2: i32,
    pub n1: i32,
    pub n2: i32,
    pub mat: i32,
    pub mf: i32,
    pub mt: i32,
}

/// /ERRORCOM/ of SETERROR.
#[derive(Clone, Copy, Default, Debug)]
pub struct ErrCom {
    pub errxc3: f64,
    pub erminmax: f64,
    pub ermt2: f64,
    pub ermt102: f64,
    pub ersoft: f64,
    pub xclow: f64,
}

impl ErrCom {
    /// `SETERROR(ICODE)`.
    pub fn seterror(&mut self, icode: i32) {
        if self.errxc3 <= 0.0 {
            self.errxc3 = 1.0e-3;
        }
        if self.xclow <= 0.0 {
            self.xclow = 1.0e-30;
        }
        if icode == 1 {
            self.erminmax = self.errxc3;
            self.ermt2 = self.errxc3;
            self.ermt102 = self.errxc3;
            self.ersoft = self.errxc3;
            return;
        }
        if icode > 1 {
            self.erminmax = 0.1 * self.errxc3;
            self.ermt2 = 0.8 * self.errxc3;
            self.ermt102 = 0.5 * self.errxc3;
            self.ersoft = 1.0e1 * self.errxc3;
            if self.errxc3 <= 1.0e-4 {
                self.ermt2 = self.errxc3;
                self.ermt102 = self.errxc3;
                if self.erminmax < 1.0e-5 {
                    self.erminmax = 1.0e-5;
                }
                if self.ersoft > 1.0e-3 {
                    self.ersoft = 1.0e-3;
                }
                if self.ersoft < self.errxc3 {
                    self.ersoft = self.errxc3;
                }
            }
        }
    }
}

/// /MAXIECOM/ state.
pub struct Maxie {
    xmaxie: Vec<f64>,
    ymaxie: Vec<f64>,
    matm: Vec<i32>,
    mfm: Vec<i32>,
    mtm: Vec<i32>,
    maxie: usize,
    mypass: i32,
    mygoof: i32,
    pub xminimum: f64,
    xsmallest: f64,
    pub xmaximum: f64,
}

impl Default for Maxie {
    fn default() -> Self {
        Maxie {
            xmaxie: vec![0.0; 1000],
            ymaxie: vec![0.0; 1000],
            matm: vec![0; 1000],
            mfm: vec![0; 1000],
            mtm: vec![0; 1000],
            maxie: 0,
            mypass: 0,
            mygoof: 0,
            xminimum: 0.0,
            xsmallest: 0.0,
            xmaximum: 0.0,
        }
    }
}

pub struct Endf {
    pub inp: i32,
    pub outp: i32,
    pub itape: i32,
    pub otape: i32,
    /// COMMON/HEADER/.
    pub h: Head,
    /// COMMON/LEADER/.
    pub ld: Head,
    /// COMMON/LASTE/ELAST.
    pub elast: f64,
    /// COMMON/FLAGS/MINUS3,IMPLUS.
    pub minus3: i32,
    pub implus: i32,
    pointo_pass: bool,
    /// COMMON/COPC/CARD(17) as 66 bytes.
    pub card: [u8; 66],
    /// COMMON/COPI/MFIELD(3).
    pub mfield: [i32; 3],
    /// COMMON/NAMECODE/CODENAME (8 characters).
    pub codename: String,
    pub u: Units,
    tstart: f32,
    tpass: i32,
    pub maxie: Maxie,
    /// /MATMF2COM/MATMF2.
    pub matmf2: i32,
    pub err: ErrCom,
}

/// A record with its three trailing fields MAT, MF, MT.
pub fn ctl(mat: i32, mf: i32, mt: i32) -> String {
    let mut s = String::with_capacity(9);
    s.push_str(&ffmt::i(mat, 4));
    s.push_str(&ffmt::i(mf, 2));
    s.push_str(&ffmt::i(mt, 3));
    s
}

impl Endf {
    pub fn new(codename: &str) -> Endf {
        Endf {
            inp: 0,
            outp: 0,
            itape: 0,
            otape: 0,
            h: Head::default(),
            ld: Head::default(),
            elast: 0.0,
            minus3: 0,
            implus: 0,
            pointo_pass: false,
            card: [b' '; 66],
            mfield: [0; 3],
            codename: format!("{:<8}", codename),
            u: Units::new(),
            tstart: 0.0,
            tpass: 0,
            maxie: Maxie::default(),
            matmf2: 0,
            err: ErrCom::default(),
        }
    }

    /// Write to the listing unit OUTP.
    pub fn lst<S: AsRef<[u8]>>(&mut self, s: S) {
        let u = self.outp;
        self.u.write(u, s);
    }

    /// Write to the listing and to the screen.
    pub fn both<S: AsRef<[u8]>>(&mut self, s: S) {
        let u = self.outp;
        self.u.write(u, s.as_ref());
        fio::out(s);
    }

    fn rd_itape(&mut self) -> Vec<u8> {
        let t = self.itape;
        self.u.read_req(t)
    }

    // ---------------- control records ----------------

    /// `CONTI`: read a control line into /HEADER/ (EOF = TEND).
    pub fn conti(&mut self) {
        let t = self.itape;
        match self.u.read(t) {
            Rd::Eof => {
                self.h.mat = -1;
                self.h.mf = 0;
                self.h.mt = 0;
            }
            Rd::Line(l) => {
                self.h.c1 = in9(&fio::rd_a(&l, 0, 11));
                self.h.c2 = in9(&fio::rd_a(&l, 11, 11));
                self.h.l1 = fio::rd_i(&l, 22, 11, t);
                self.h.l2 = fio::rd_i(&l, 33, 11, t);
                self.h.n1 = fio::rd_i(&l, 44, 11, t);
                self.h.n2 = fio::rd_i(&l, 55, 11, t);
                self.h.mat = fio::rd_i(&l, 66, 4, t);
                self.h.mf = fio::rd_i(&l, 70, 2, t);
                self.h.mt = fio::rd_i(&l, 72, 3, t);
            }
        }
    }

    /// `CONTO`: write /HEADER/ as a control line.
    pub fn conto(&mut self) {
        if self.otape <= 0 {
            return;
        }
        if self.h.mt <= 0 {
            let (m, f) = (self.h.mat, self.h.mf);
            self.outs(m, f);
            return;
        }
        let h = self.h;
        self.cont_line(h.c1, h.c2, h.l1, h.l2, h.n1, h.n2);
    }

    fn cont_line(&mut self, c1: f64, c2: f64, l1: i32, l2: i32, n1: i32, n2: i32) {
        let mut f = [b' '; 22];
        out9(c1, &mut f[0..11]);
        out9(c2, &mut f[11..22]);
        let mut l = L::new().b(&f);
        for v in [l1, l2, n1, n2] {
            l.push(&ffmt::i(if v.unsigned_abs() == 0 { 0 } else { v }, 11));
        }
        l.push(&ctl(self.h.mat, self.h.mf, self.h.mt));
        let o = self.otape;
        self.u.write(o, l);
    }

    /// `CARDI`: read a line; MAT/MF/MT go to /LEADER/.
    pub fn cardi(&mut self) -> (f64, f64, i32, i32, i32, i32) {
        let t = self.itape;
        let l = self.rd_itape();
        let c1 = in9(&fio::rd_a(&l, 0, 11));
        let c2 = in9(&fio::rd_a(&l, 11, 11));
        let l1 = fio::rd_i(&l, 22, 11, t);
        let l2 = fio::rd_i(&l, 33, 11, t);
        let n1 = fio::rd_i(&l, 44, 11, t);
        let n2 = fio::rd_i(&l, 55, 11, t);
        self.ld.mat = fio::rd_i(&l, 66, 4, t);
        self.ld.mf = fio::rd_i(&l, 70, 2, t);
        self.ld.mt = fio::rd_i(&l, 72, 3, t);
        (c1, c2, l1, l2, n1, n2)
    }

    /// `CARDO`: write a line with /HEADER/ MAT/MF/MT.
    pub fn cardo(&mut self, c1: f64, c2: f64, l1: i32, l2: i32, n1: i32, n2: i32) {
        if self.otape <= 0 {
            return;
        }
        if self.h.mt <= 0 {
            let (m, f) = (self.h.mat, self.h.mf);
            self.outs(m, f);
            return;
        }
        self.cont_line(c1, c2, l1, l2, n1, n2);
    }

    /// `CARDIO`.
    pub fn cardio(&mut self) -> (f64, f64, i32, i32, i32, i32) {
        let r = self.cardi();
        self.cardo(r.0, r.1, r.2, r.3, r.4, r.5);
        r
    }

    // ---------------- interpolation tables ----------------

    /// `TERPI`: read N1 interpolation pairs.
    pub fn terpi(&mut self, nbt: &mut [i32], int: &mut [i32], n1: usize) {
        let t = self.itape;
        let mut i = 0;
        while i < n1 {
            let l = self.rd_itape();
            for k in 0..3 {
                if i >= n1 {
                    break;
                }
                nbt[i] = fio::rd_i(&l, 22 * k, 11, t);
                int[i] = fio::rd_i(&l, 22 * k + 11, 11, t);
                i += 1;
            }
        }
    }

    /// `TERPO`.
    pub fn terpo(&mut self, nbt: &[i32], int: &[i32], n1: usize) {
        if self.otape <= 0 {
            return;
        }
        let mut i1 = 0;
        while i1 < n1 {
            let i2 = (i1 + 3).min(n1);
            let mut l = String::with_capacity(80);
            for ii in i1..i2 {
                l.push_str(&ffmt::i(nbt[ii], 11));
                l.push_str(&ffmt::i(int[ii], 11));
            }
            for _ in (i2 - i1)..3 {
                l.push_str("                      ");
            }
            l.push_str(&ctl(self.h.mat, self.h.mf, self.h.mt));
            let o = self.otape;
            self.u.write(o, l);
            i1 += 3;
        }
    }

    /// `TERPTEST`.
    pub fn terptest(&mut self, nbt: &[i32], int: &[i32], n1: i32, n2: i32, mypath: i32) {
        let (mat, mf, mt) = (self.ld.mat, self.ld.mf, self.ld.mt);
        let err_hdr = format!(" {}ERROR==", "ERROR===".repeat(9));
        let pairs = |s: &mut Vec<String>| {
            s.push("      NBT  INT".to_string());
            for k in 0..n1.max(0) as usize {
                s.push(format!("{}{}", ffmt::i(nbt[k], 9), ffmt::i(int[k], 5)));
            }
        };
        let mut lines: Vec<String> = Vec::new();
        let head2 = |n2: i32| format!(" MAT/MF/MT/N2={}{}{}{}", ffmt::i(mat, 5), ffmt::i(mf, 3), ffmt::i(mt, 4), ffmt::i(n2, 10));
        let mut warn_only = false;
        let mut linmsg = false;
        if n1 <= 0 {
            lines.push(err_hdr.clone());
            lines.push(format!(" MAT/MF/MT/N1/N2={}{}{}{}{}", ffmt::i(mat, 5), ffmt::i(mf, 3), ffmt::i(mt, 4), ffmt::i(n1, 10), ffmt::i(n2, 10)));
            lines.push(" No Interpolation Law - N1 MUST be Positive".to_string());
        } else if n1 > 100 {
            lines.push(err_hdr.clone());
            lines.push(format!(" MAT/MF/MT/N1/N2={}{}{}{}{}", ffmt::i(mat, 5), ffmt::i(mf, 3), ffmt::i(mt, 4), ffmt::i(n1, 10), ffmt::i(n2, 10)));
            lines.push(format!(" Interpolation Law N1 Exceeds Dimension{}", ffmt::i(100, 5)));
        } else if n2 < 2 {
            lines.push(err_hdr.clone());
            lines.push(head2(n2));
            lines.push(" Point Count (N2) MUST be 2 or More".to_string());
            pairs(&mut lines);
        } else if n2 != nbt[n1 as usize - 1] {
            lines.push(err_hdr.clone());
            lines.push(head2(n2));
            lines.push(" Point Count (N2) and Interpolation Law do not Agree".to_string());
            pairs(&mut lines);
        } else {
            let mut nbtlast = 0;
            let mut bad = 0;
            for k in 0..n1 as usize {
                let ik = int[k];
                if mypath == 0 {
                    if (1..=6).contains(&ik) || (11..=15).contains(&ik) || (21..=25).contains(&ik) {
                        continue;
                    }
                    bad = 100;
                    break;
                }
                if mypath == 1 && !(1..=6).contains(&ik) {
                    bad = 120;
                    break;
                }
                if mypath == 2 && !(1..=2).contains(&ik) {
                    bad = 140;
                    break;
                }
                if mypath == 3 && ik != 2 {
                    bad = 160;
                    break;
                }
                if nbt[k] <= nbtlast {
                    bad = 190;
                    break;
                }
                nbtlast = nbt[k];
            }
            match bad {
                0 => return,
                100 => {
                    lines.push(err_hdr.clone());
                    lines.push(head2(n2));
                    lines.push(" Interpolation Law (INT) MUST be 1 to 6, 11 to 15, or 21 to 25".to_string());
                    pairs(&mut lines);
                }
                120 => {
                    lines.push(format!(" {}WARNING", "WARNING=".repeat(9)));
                    lines.push(head2(n2));
                    lines.push(" Interpolation Law (INT) Expected 1 to 6".to_string());
                    pairs(&mut lines);
                    warn_only = true;
                }
                140 => {
                    lines.push(err_hdr.clone());
                    lines.push(head2(n2));
                    lines.push(" Interpolation Law (INT) MUST be 1 or 2".to_string());
                    pairs(&mut lines);
                    linmsg = true;
                }
                160 => {
                    lines.push(err_hdr.clone());
                    lines.push(head2(n2));
                    lines.push(" Interpolation Law (INT) MUST be 2".to_string());
                    pairs(&mut lines);
                    linmsg = true;
                }
                _ => {
                    lines.push(err_hdr.clone());
                    lines.push(head2(n2));
                    lines.push(" Interpolation Boundary (NBT) MUST be in Ascending Order".to_string());
                    pairs(&mut lines);
                }
            }
        }
        let txt = lines.join("\n");
        self.both(&txt);
        if warn_only {
            return;
        }
        if linmsg {
            self.both(" You MUST Linearize this data before using this code.\n We recommend to linearize ENDF formatted data use,\n LINEAR  for MF=3 Cross Sections.\n LEGEND  for MF=4 Angular Distributions.\n SPECTRA for MF=5 Energy Spectra.");
        }
        self.both(" Execution Terminated");
        self.enderror();
    }

    // ---------------- tables of points ----------------

    /// `POINTI`: read IXY points, check ascending x (x and y 0-based).
    pub fn pointi(&mut self, x: &mut [f64], y: &mut [f64], ixy: usize) {
        let mut i = 0;
        while i < ixy {
            let l = fio::pad(self.rd_itape(), 66);
            let n = (ixy - i).min(3);
            for k in 0..n {
                x[i + k] = in9(&l[22 * k..22 * k + 11]);
                y[i + k] = in9(&l[22 * k + 11..22 * k + 22]);
            }
            i += 3;
        }
        for i in 0..ixy {
            if x[i] < self.elast {
                if (self.elast - x[i]).abs() > 1.0e-09 * self.elast {
                    let a = f9s(self.elast);
                    let b = f9s(x[i]);
                    let s = L::new()
                        .s(&format!("  {}\n", "-".repeat(78)))
                        .s(&format!("{}{}{}\n", ffmt::i(self.h.mat, 5), ffmt::i(self.h.mf, 3), ffmt::i(self.h.mt, 4)))
                        .s(" Energies Not in Ascending Energy Order\n")
                        .s("  Index      Energy\n")
                        .s(&ffmt::i(i as i32, 7))
                        .s(" ")
                        .b(&a)
                        .s("\n")
                        .s(&ffmt::i(i as i32 + 1, 7))
                        .s(" ")
                        .b(&b)
                        .s("\n")
                        .s(&format!("{} Execution Terminated.\n", " ".repeat(19)))
                        .s(&format!("  {}", "-".repeat(78)));
                    self.lst(s);
                }
                x[i] = self.elast;
            }
            self.elast = x[i];
        }
    }

    fn point_line(&mut self, x: &[f64], y: &[f64], i1: usize, i2: usize, mf3_10: bool) {
        let mut f = [b' '; 66];
        let mut k = 0;
        for ii in i1..i2 {
            if !(self.h.mf == 27 && self.h.mt == 506) && y[ii] < 0.0 {
                self.minus3 += 1;
            }
            if y[ii] > 0.0 {
                self.implus = 1;
            }
            if mf3_10 && self.h.mf == 3 {
                out10(x[ii], &mut f[11 * k..11 * k + 11]);
            } else {
                out9g(x[ii], &mut f[11 * k..11 * k + 11]);
            }
            k += 1;
            out9g(y[ii], &mut f[11 * k..11 * k + 11]);
            k += 1;
        }
        let l = L::new().b(&f).s(&ctl(self.h.mat, self.h.mf, self.h.mt));
        let o = self.otape;
        self.u.write(o, l);
    }

    /// `POINTO`.
    pub fn pointo(&mut self, x: &[f64], y: &[f64], ixy: usize) {
        if !self.pointo_pass {
            self.minus3 = 0;
            self.implus = self.minus3;
            self.minus3 = self.implus;
            self.pointo_pass = true;
        }
        if self.otape <= 0 || ixy == 0 {
            return;
        }
        let mut i1 = 0;
        while i1 < ixy {
            let i2 = (i1 + 3).min(ixy);
            self.point_line(x, y, i1, i2, true);
            i1 += 3;
        }
    }

    /// `POINTO9`.
    pub fn pointo9(&mut self, x: &[f64], y: &[f64], ixy: usize) {
        if self.otape <= 0 || ixy == 0 {
            return;
        }
        let mut i1 = 0;
        while i1 < ixy {
            let i2 = (i1 + 3).min(ixy);
            self.point_line(x, y, i1, i2, false);
            i1 += 3;
        }
    }

    /// `POINTIV`: read points from unit iu (no order check).
    pub fn pointiv(&mut self, iu: i32, x: &mut [f64], y: &mut [f64], ixy: usize) {
        let mut i = 0;
        while i < ixy {
            let l = fio::pad(self.u.read_req(iu), 66);
            let n = (ixy - i).min(3);
            for k in 0..n {
                x[i + k] = in9(&l[22 * k..22 * k + 11]);
                y[i + k] = in9(&l[22 * k + 11..22 * k + 22]);
            }
            i += 3;
        }
    }

    // ---------------- lists ----------------

    /// `LISTI`.
    pub fn listi(&mut self, x: &mut [f64], ix: usize) {
        let mut i1 = 0;
        while i1 < ix {
            let l = fio::pad(self.rd_itape(), 66);
            let n = (ix - i1).min(6);
            for k in 0..n {
                x[i1 + k] = in9(&l[11 * k..11 * k + 11]);
            }
            i1 += 6;
        }
    }

    /// `LISTIV`.
    pub fn listiv(&mut self, iu: i32, x: &mut [f64], ix: usize) {
        let mut i1 = 0;
        while i1 < ix {
            let l = fio::pad(self.u.read_req(iu), 66);
            let n = (ix - i1).min(6);
            for k in 0..n {
                x[i1 + k] = in9(&l[11 * k..11 * k + 11]);
            }
            i1 += 6;
        }
    }

    /// `LISTSKIP`.
    pub fn listskip(&mut self, ix: usize) {
        let mut i = 0;
        while i < ix {
            self.rd_itape();
            i += 6;
        }
    }

    /// `LISTO`.
    pub fn listo(&mut self, x: &[f64], ix: usize) {
        if self.otape <= 0 || ix == 0 {
            return;
        }
        let mut i1 = 0;
        while i1 < ix {
            let i2 = (i1 + 6).min(ix);
            let mut f = [b' '; 66];
            for (k, l) in (i1..i2).enumerate() {
                if self.h.mf == 2 && k == 0 {
                    out10(x[l], &mut f[11 * k..11 * k + 11]);
                } else {
                    out9(x[l], &mut f[11 * k..11 * k + 11]);
                }
            }
            let l = L::new().b(&f).s(&ctl(self.h.mat, self.h.mf, self.h.mt));
            let o = self.otape;
            self.u.write(o, l);
            i1 += 6;
        }
    }

    /// `LISTO9`.
    pub fn listo9(&mut self, x: &[f64], ix: usize) {
        if self.otape <= 0 || ix == 0 {
            return;
        }
        let mut i1 = 0;
        while i1 < ix {
            let i2 = (i1 + 6).min(ix);
            let mut f = [b' '; 66];
            for (k, l) in (i1..i2).enumerate() {
                out9g(x[l], &mut f[11 * k..11 * k + 11]);
            }
            let l = L::new().b(&f).s(&ctl(self.h.mat, self.h.mf, self.h.mt));
            let o = self.otape;
            self.u.write(o, l);
            i1 += 6;
        }
    }

    pub fn listio(&mut self, x: &mut [f64], ix: usize) {
        self.listi(x, ix);
        self.listo(x, ix);
    }

    pub fn listio9(&mut self, x: &mut [f64], ix: usize) {
        self.listi(x, ix);
        self.listo9(x, ix);
    }

    // ---------------- whole lines ----------------

    /// `LINEIN`: CARD and MAT/MF/MT into /HEADER/.
    pub fn linein(&mut self) {
        let t = self.itape;
        let l = fio::pad(self.rd_itape(), 66);
        self.card.copy_from_slice(&l[..66]);
        self.h.mat = fio::rd_i(&l, 66, 4, t);
        self.h.mf = fio::rd_i(&l, 70, 2, t);
        self.h.mt = fio::rd_i(&l, 72, 3, t);
    }

    /// `LINEOUT`.
    pub fn lineout(&mut self) {
        if self.otape <= 0 {
            return;
        }
        if self.h.mt <= 0 {
            let (m, f) = (self.h.mat, self.h.mf);
            self.outs(m, f);
            return;
        }
        let l = L::new().b(&self.card).s(&ctl(self.h.mat, self.h.mf, self.h.mt));
        let o = self.otape;
        self.u.write(o, l);
    }

    /// One line read by the COPY entries into CARD/MFIELD.
    fn copy_read(&mut self) {
        let t = self.itape;
        let l = fio::pad(self.rd_itape(), 66);
        self.card.copy_from_slice(&l[..66]);
        self.mfield[0] = fio::rd_i(&l, 66, 4, t);
        self.mfield[1] = fio::rd_i(&l, 70, 2, t);
        self.mfield[2] = fio::rd_i(&l, 72, 3, t);
    }

    fn copy_write(&mut self) {
        if self.otape <= 0 {
            return;
        }
        if self.mfield[2] > 0 {
            let l = L::new().b(&self.card).s(&ctl(self.mfield[0], self.mfield[1], self.mfield[2]));
            let o = self.otape;
            self.u.write(o, l);
        } else {
            let (a, b) = (self.mfield[0], self.mfield[1]);
            self.outs(a, b);
        }
    }

    /// `COPYT`: copy to the TEND line.
    pub fn copyt(&mut self) {
        loop {
            self.copy_read();
            self.copy_write();
            if self.mfield[0] < 0 {
                return;
            }
        }
    }

    /// `COPYM`: copy to the MEND line.
    pub fn copym(&mut self) {
        loop {
            self.copy_read();
            self.copy_write();
            if self.mfield[0] <= 0 {
                return;
            }
        }
    }

    /// `COPYF`: copy to the FEND line.
    pub fn copyf(&mut self) {
        loop {
            self.copy_read();
            self.copy_write();
            if self.mfield[1] <= 0 {
                return;
            }
        }
    }

    /// `COPYS`: copy to the SEND line.
    pub fn copys(&mut self) {
        loop {
            self.copy_read();
            self.copy_write();
            if self.mfield[2] <= 0 {
                return;
            }
        }
    }

    /// `COPYL`: copy the tape label.
    pub fn copyl(&mut self) {
        self.copy_read();
        if self.otape <= 0 {
            return;
        }
        if self.mfield[0] == 0 {
            self.mfield[0] = 7000;
        }
        self.mfield[1] = 0;
        self.mfield[2] = 0;
        let l = L::new().b(&self.card).s(&ctl(self.mfield[0], 0, 0));
        let o = self.otape;
        self.u.write(o, l);
    }

    /// `COPY1`: copy one line.
    pub fn copy1(&mut self) {
        self.copy_read();
        self.copy_write();
    }

    // ---------------- end lines ----------------

    pub fn outt(&mut self) {
        if self.otape <= 0 {
            return;
        }
        let o = self.otape;
        self.u.write(o, format!("{}  -1 0  0", " ".repeat(66)));
    }

    pub fn outm(&mut self) {
        if self.otape <= 0 {
            return;
        }
        let o = self.otape;
        self.u.write(o, format!("{}   0 0  0", " ".repeat(66)));
    }

    pub fn outf(&mut self, mat: i32) {
        if self.otape <= 0 {
            return;
        }
        let o = self.otape;
        self.u.write(o, format!("{}{} 0  0", " ".repeat(66), ffmt::i(mat, 4)));
    }

    pub fn outs(&mut self, mat: i32, mf: i32) {
        if self.otape <= 0 {
            return;
        }
        let o = self.otape;
        if mf <= 0 {
            self.u.write(o, format!("{}{} 0  0", " ".repeat(66), ffmt::i(mat, 4)));
        } else {
            self.u.write(o, format!("{}{}{}  0", " ".repeat(66), ffmt::i(mat, 4), ffmt::i(mf, 2)));
        }
    }

    // ---------------- skipping ----------------

    fn skip_read(&mut self) {
        let t = self.itape;
        let l = self.rd_itape();
        self.mfield[0] = fio::rd_i(&l, 66, 4, t);
        self.mfield[1] = fio::rd_i(&l, 70, 2, t);
        self.mfield[2] = fio::rd_i(&l, 72, 3, t);
    }

    pub fn skipt(&mut self) {
        loop {
            self.skip_read();
            if self.mfield[0] < 0 {
                return;
            }
        }
    }
    pub fn skipm(&mut self) {
        loop {
            self.skip_read();
            if self.mfield[0] <= 0 {
                return;
            }
        }
    }
    pub fn skipf(&mut self) {
        loop {
            self.skip_read();
            if self.mfield[1] <= 0 {
                return;
            }
        }
    }
    pub fn skips(&mut self) {
        loop {
            self.skip_read();
            if self.mfield[2] <= 0 {
                return;
            }
        }
    }
    pub fn skipl(&mut self) {
        self.skip_read();
    }
    pub fn skip1(&mut self) {
        self.skip_read();
    }

    // ---------------- hollerith ----------------

    /// `HOLLYI`: 66 characters.
    pub fn hollyi(&mut self) -> [u8; 66] {
        let l = fio::pad(self.rd_itape(), 66);
        let mut r = [b' '; 66];
        r.copy_from_slice(&l[..66]);
        r
    }

    /// `HOLLYO`.
    pub fn hollyo(&mut self, line66: &[u8]) {
        if self.otape <= 0 {
            return;
        }
        let l = L::new().b(&line66[..66]).s(&ctl(self.h.mat, self.h.mf, self.h.mt));
        let o = self.otape;
        self.u.write(o, l);
    }

    // ---------------- MAXIE ----------------

    pub fn maxie0(&mut self) {
        let m = &mut self.maxie;
        m.mypass = 1;
        m.maxie = 0;
        m.xmaximum = 0.0;
        m.xminimum = 0.0;
    }

    pub fn maxie1(&mut self, mat: i32, mf: i32, mt: i32) {
        let m = &mut self.maxie;
        if m.mypass <= 0 {
            m.mypass = 1;
            m.maxie = 0;
            m.xmaximum = 0.0;
            m.xminimum = 0.0;
        }
        if m.maxie < 999 {
            m.maxie += 1;
        }
        let k = m.maxie;
        m.matm[k] = mat;
        m.mfm[k] = mf;
        m.mtm[k] = mt;
        m.xmaxie[k] = 0.0;
        m.ymaxie[k] = 0.0;
    }

    pub fn maxie2(&mut self, xm: f64, ym: f64) {
        if xm <= 0.0 {
            return;
        }
        let m = &mut self.maxie;
        if m.mypass <= 0 {
            m.mypass = 1;
            m.maxie = 1;
            m.xmaximum = 0.0;
            m.xminimum = 0.0;
            m.matm[1] = 0;
            m.mfm[1] = 0;
            m.mtm[1] = 0;
        }
        let k = m.maxie;
        m.xmaxie[k] = xm;
        m.ymaxie[k] = ym;
        if m.mfm[k] == 27 {
            m.ymaxie[k] = 0.0;
        }
        if m.mfm[k] == 26 {
            m.ymaxie[k] = 0.0;
        }
    }

    /// `MAXIE3`: summary on unit 3 (and the screen if itty > 0).
    pub fn maxie3(&mut self, itty: i32) {
        if self.maxie.mypass <= 0 || self.maxie.maxie == 0 {
            return;
        }
        let n = self.maxie.maxie;
        {
            let m = &mut self.maxie;
            m.xmaximum = 0.0;
            m.xminimum = 0.0;
            for imt in 1..=n {
                if m.xmaxie[imt] > m.xmaximum {
                    m.xmaximum = m.xmaxie[imt];
                }
                if m.ymaxie[imt] > 0.0 {
                    if m.xminimum <= 0.0 {
                        m.xminimum = m.xmaxie[imt];
                    } else if m.xmaxie[imt] < m.xminimum {
                        m.xminimum = m.xmaxie[imt];
                    }
                }
            }
            if m.xminimum <= 0.0 {
                for imt in 1..=n {
                    if m.xmaxie[imt] > m.xminimum {
                        m.xminimum = m.xmaxie[imt];
                    }
                }
            }
            if m.xminimum <= 0.0 || m.xmaximum <= 0.0 {
                return;
            }
        }
        let dash = format!(" {}", "-".repeat(78));
        let fx = f9s(self.maxie.xmaximum);
        let s = L::new()
            .s(&dash)
            .s("\n All Cross Sections Read. Maximum Tabulated Energy = ")
            .b(&fx)
            .s(" eV\n")
            .s(&dash);
        self.u.write(3, &s);
        if itty > 0 {
            fio::out(&s);
        }
        if self.maxie.xminimum >= self.maxie.xmaximum {
            return;
        }
        self.maxie.mygoof = 1;
        let s = format!(" WARNING - The Following MTs Stop below the Maximum Tabulated Energy\n{}\n  MAT MF  MT  This MT     Maximum E\n{}", dash, dash);
        self.u.write(3, &s);
        if itty > 0 {
            fio::out(&s);
        }
        for imt in 1..=n {
            let m = &self.maxie;
            if m.ymaxie[imt] > 0.0 && m.xmaxie[imt] < m.xmaximum {
                let a = f9s(m.xmaxie[imt]);
                let s = L::new()
                    .s(&format!("{}{}{} ", ffmt::i(m.matm[imt], 5), ffmt::i(m.mfm[imt], 3), ffmt::i(m.mtm[imt], 4)))
                    .b(&a)
                    .s(" ")
                    .b(&fx);
                self.u.write(3, &s);
                if itty > 0 {
                    fio::out(&s);
                }
            }
        }
        let a = f9s(self.maxie.xminimum);
        let s = L::new()
            .s(&dash)
            .s("\n WARNING - Data between")
            .b(&a)
            .s(" and ")
            .b(&fx)
            .s(" are NOT Reliable\n")
            .s(&dash);
        self.u.write(3, &s);
        if itty > 0 {
            fio::out(&s);
        }
        let m = &mut self.maxie;
        if m.xsmallest <= 0.0 {
            m.xsmallest = m.xminimum;
        } else if m.xminimum < m.xsmallest {
            m.xsmallest = m.xminimum;
        }
    }

    /// `MAXIE4`.
    pub fn maxie4(&mut self, jtty: i32) {
        let m = &self.maxie;
        if m.mypass <= 0 || m.mygoof <= 0 || m.xsmallest <= 0.0 {
            return;
        }
        let a = f9s(m.xsmallest);
        let dash = format!(" {}", "-".repeat(78));
        let s = L::new()
            .s(&dash)
            .s("\n WARNING - The ENDF Output is NOT Relialble above")
            .b(&a)
            .s(" eV,\n")
            .s("           because it is NOT Consistently Defined for ALL MTs.\n")
            .s("           It is STRONGLY recommended that you NOT use this\n")
            .s("           data above this energy. If you need data up to\n")
            .s("           energies higher than this PLEASE correct the data\n")
            .s("           and only then try re-running this code.\n")
            .s(&dash);
        self.u.write(3, &s);
        if jtty > 0 {
            fio::out(&s);
        }
    }

    // ---------------- resonance limits ----------------

    fn list2_text(res1: f64, res2: f64, ures1: f64, ures2: f64) -> String {
        let eq = format!(" {}", "=".repeat(79));
        format!(
            "{eq}\n Resolved   Resonance Region Lowes Energy Limit{} eV\n Resolved   Resonance Region Upper Energy Limit{} eV\n Unresolved Resonance Region Lowes Energy Limit{} eV\n Unresolved Resonance Region Upper Energy Limit{} eV\n{eq}",
            dp(res1, 11, 4),
            dp(res2, 11, 4),
            dp(ures1, 11, 4),
            dp(ures2, 11, 4)
        )
    }

    /// `LIST2`.
    pub fn list2(&mut self, res1: f64, res2: f64, ures1: f64, ures2: f64) {
        let s = Self::list2_text(res1, res2, ures1, ures2);
        if self.outp > 0 {
            self.lst(&s);
        }
        if self.otape > 0 {
            fio::out(&s);
        }
    }

    /// `LIST2G`.
    pub fn list2g(&mut self, res1: f64, res2: f64, ures1: f64, ures2: f64) {
        let s = Self::list2_text(res1, res2, ures1, ures2);
        if self.outp > 0 {
            self.lst(&s);
        }
        fio::out(&s);
    }

    // ---------------- timer and stop ----------------

    fn time_lines(&mut self, label: &str, t: f32) {
        let eq = format!(" {}", "=".repeat(79));
        let s = format!("{eq}\n {} Total Execution Time{} Seconds\n{eq}", label, ffmt::f(t as f64, 20, 2));
        self.lst(&s);
        if self.otape > 0 {
            fio::out(&s);
        }
    }

    /// `TIMER`.
    pub fn timer(&mut self) {
        let tnow = cpu_seconds();
        if self.tpass == 0 {
            self.tstart = tnow;
        }
        self.tpass += 1;
        if self.tpass <= 1 {
            return;
        }
        let name = self.codename.clone();
        self.time_lines(&name, tnow - self.tstart);
    }

    /// `TIMER1`.
    pub fn timer1(&mut self) -> f32 {
        cpu_seconds() - self.tstart
    }

    /// `TIMEMAT`.
    pub fn timemat(&mut self) {
        let t = cpu_seconds() - self.tstart;
        self.time_lines("     MAT", t);
    }

    /// `TIMEERR`.
    pub fn timeerr(&mut self) {
        let t = cpu_seconds() - self.tstart;
        self.time_lines("   ERROR", t);
    }

    /// `ENDIT`: normal end.
    pub fn endit(&mut self) -> ! {
        self.timer();
        fio::stop(&mut self.u)
    }

    /// `ENDERROR`: stop on an error.
    pub fn enderror(&mut self) -> ! {
        self.timeerr();
        fio::stop(&mut self.u)
    }
}

/// `1PDw.d`.
pub fn dp(x: f64, w: usize, d: usize) -> String {
    ffmt::ep(x, w, d, 1).replace('E', "D")
}

/// `1PEw.d`.
pub fn ep1(x: f64, w: usize, d: usize) -> String {
    ffmt::ep(x, w, d, 1)
}

/// ETIME: user plus system CPU time of the process, in REAL*4.
pub fn cpu_seconds() -> f32 {
    #[cfg(unix)]
    {
        extern "C" {
            fn clock() -> std::os::raw::c_long;
        }
        // POSIX: CLOCKS_PER_SEC is 1000000.
        let c = unsafe { clock() };
        (c as f64 / 1.0e6) as f32
    }
    #[cfg(not(unix))]
    {
        static T0: OnceLock<std::time::Instant> = OnceLock::new();
        T0.get_or_init(std::time::Instant::now).elapsed().as_secs_f32()
    }
}
