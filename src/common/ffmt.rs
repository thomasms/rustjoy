//! Fortran formatted-I/O editing, matching libgfortran.
//!
//! Output: `Iw`, `Fw.d` (with `kP` scale factor), `Ew.d` (with `kP`), `Aw`.
//! Input:  `Iw`, `Fw.d`/`Ew.d` (blank-null mode, which is gfortran's default).
//!
//! All decimal conversions are exact and round-half-even on exact ties, which
//! is what both libgfortran (via glibc printf) and Rust's formatter do.

#![allow(dead_code)]

/// Small stack buffer used to build one edited field without heap
/// allocation. `push` reports overflow; callers then use the `_ref` code.
struct Buf {
    b: [u8; 96],
    n: usize,
    over: bool,
}

impl Buf {
    #[inline]
    fn new() -> Buf {
        Buf { b: [0; 96], n: 0, over: false }
    }
    #[inline]
    fn push(&mut self, c: u8) {
        if self.n < self.b.len() {
            self.b[self.n] = c;
            self.n += 1;
        } else {
            self.over = true;
        }
    }
    #[inline]
    fn bytes(&self) -> &[u8] {
        &self.b[..self.n]
    }
}

impl std::fmt::Write for Buf {
    fn write_str(&mut self, s: &str) -> std::fmt::Result {
        for &c in s.as_bytes() {
            self.push(c);
        }
        if self.over { Err(std::fmt::Error) } else { Ok(()) }
    }
}

// ------------------------------------------------------------------
// Exact decimal digits with integer arithmetic. These give the same
// digits as Rust's `{:.*e}` and `{:.*}` (exact value, round half to even)
// without its general-purpose machinery; when a value needs more than 128
// bits they return None and the callers use the formatter.

/// ax = m*2^e for a finite ax >= 0.
#[inline]
fn decompose(ax: f64) -> (u64, i32) {
    let b = ax.to_bits();
    let ex = ((b >> 52) & 0x7ff) as i32;
    let fr = b & ((1u64 << 52) - 1);
    if ex == 0 {
        (fr, -1074)
    } else {
        (fr | (1u64 << 52), ex - 1075)
    }
}

const POW10: [u128; 39] = {
    let mut t = [1u128; 39];
    let mut i = 1;
    while i < 39 {
        t[i] = t[i - 1] * 10;
        i += 1;
    }
    t
};

/// floor(m*2^e*10^s) and where the remainder r lies: 0 if r=0, 1 if
/// 0<r<1/2, 2 if r=1/2, 3 if r>1/2. None if it does not fit in u128.
#[inline]
fn scaled(m: u64, e: i32, s: i32) -> Option<(u128, u8)> {
    if m == 0 {
        return Some((0, 0));
    }
    let cmp = |r: u128, den: u128| -> u8 {
        if r == 0 {
            0
        } else {
            // den < 2^127, so 2r does not overflow
            match (2 * r).cmp(&den) {
                std::cmp::Ordering::Less => 1,
                std::cmp::Ordering::Equal => 2,
                std::cmp::Ordering::Greater => 3,
            }
        }
    };
    if s >= 0 {
        let p = *POW10.get(s as usize)?;
        let num = (m as u128).checked_mul(p)?;
        if e >= 0 {
            if e >= 127 || num.leading_zeros() <= e as u32 {
                return None;
            }
            return Some((num << e, 0));
        }
        let sh = (-e) as u32;
        if sh >= 127 {
            return None;
        }
        let den = 1u128 << sh;
        Some((num >> sh, cmp(num & (den - 1), den)))
    } else {
        let p = *POW10.get((-s) as usize)?;
        let (num, den) = if e >= 0 {
            if e >= 127 || (m as u128).leading_zeros() <= e as u32 {
                return None;
            }
            ((m as u128) << e, p)
        } else {
            let sh = (-e) as u32;
            if sh >= 127 || p.leading_zeros() <= sh + 1 {
                return None;
            }
            (m as u128, p << sh)
        };
        if den >= 1u128 << 126 {
            return None;
        }
        Some((num / den, cmp(num % den, den)))
    }
}

#[inline]
fn round_even(q: u128, f: u8) -> u128 {
    if f == 3 || (f == 2 && q & 1 == 1) {
        q + 1
    } else {
        q
    }
}

/// The n significant digits of a finite ax > 0 as Rust's `{:.(n-1)e}`
/// gives them: (q, k) with 10^(n-1) <= q < 10^n and exponent k.
fn sci_digits(ax: f64, n: usize) -> Option<(u64, i32)> {
    if !(1..=19).contains(&n) || !(ax > 0.0) || !ax.is_finite() {
        return None;
    }
    let (m, e) = decompose(ax);
    let lo = POW10[n - 1];
    let hi = POW10[n];
    let mut k = ax.log10().floor() as i32;
    for _ in 0..4 {
        let s = n as i32 - 1 - k;
        let (q, f) = scaled(m, e, s)?;
        if q < lo {
            k -= 1;
            continue;
        }
        if q >= hi {
            k += 1;
            continue;
        }
        let mut q = round_even(q, f);
        if q == hi {
            q = lo;
            k += 1;
        }
        return Some((q as u64, k));
    }
    None
}

/// Rust's `{:.dd}` of a finite ax >= 0 into b.
fn fixed_digits(ax: f64, dd: usize, b: &mut Buf) -> bool {
    if !ax.is_finite() || dd > 38 {
        return false;
    }
    let (m, e) = decompose(ax);
    let q = match scaled(m, e, dd as i32) {
        Some((q, f)) => round_even(q, f),
        None => return false,
    };
    let mut dg = [0u8; 40];
    let mut nd = 0;
    let mut v = q;
    while v > 0 {
        dg[nd] = b'0' + (v % 10) as u8;
        v /= 10;
        nd += 1;
    }
    // digits most significant first: dg[nd-1] .. dg[0]
    let digit = |j: usize| -> u8 { if j < nd { dg[j] } else { b'0' } };
    let nint = if nd > dd { nd - dd } else { 1 };
    for j in (dd..dd + nint).rev() {
        b.push(digit(j));
    }
    if dd > 0 {
        b.push(b'.');
        for j in (0..dd).rev() {
            b.push(digit(j));
        }
    }
    !b.over
}

/// Right-justify `s` in a field of width w (stars if it does not fit).
#[inline]
fn rjust_bytes(s: &[u8], w: usize) -> String {
    let mut out = String::with_capacity(w);
    if s.len() > w {
        for _ in 0..w {
            out.push('*');
        }
    } else {
        for _ in s.len()..w {
            out.push(' ');
        }
        // s is ASCII
        out.push_str(std::str::from_utf8(s).unwrap());
    }
    out
}

fn stars(w: usize) -> String {
    "*".repeat(w)
}

/// `Iw`.
pub fn i(v: i32, w: usize) -> String {
    let mut d = [0u8; 12];
    let mut n = 0;
    let mut u = (v as i64).unsigned_abs();
    loop {
        d[n] = b'0' + (u % 10) as u8;
        n += 1;
        u /= 10;
        if u == 0 {
            break;
        }
    }
    let mut s = [0u8; 12];
    let mut m = 0;
    if v < 0 {
        s[0] = b'-';
        m = 1;
    }
    for k in (0..n).rev() {
        s[m] = d[k];
        m += 1;
    }
    rjust_bytes(&s[..m], w)
}

/// `I0`.
pub fn i0(v: i32) -> String {
    v.to_string()
}

/// `Aw` output of a string (Fortran semantics: if the string is longer than w,
/// the leftmost w characters are written; if shorter, it is right-justified).
pub fn a(s: &str, w: usize) -> String {
    a_ref(s, w)
}

/// `kP,Fw.d`. Same rules as `fp_ref` (the digits come from Rust's exact
/// decimal conversion in both); only the string handling differs.
fn fp_core(x: f64, w: usize, d: usize, k: i32) -> Option<Buf> {
    use std::fmt::Write;
    let neg = x.is_sign_negative();
    let ax = x.abs();
    let dd = (d as i32 + k).max(0) as usize;
    let mut raw = Buf::new();
    if !fixed_digits(ax, dd, &mut raw) {
        raw = Buf::new();
        if write!(raw, "{:.*}", dd, ax).is_err() {
            return None;
        }
    }
    let raw = raw.bytes();
    // body: the decimal string with the point moved right by k
    let mut body = Buf::new();
    if k == 0 {
        for &c in raw {
            body.push(c);
        }
    } else {
        let p = raw.iter().position(|&c| c == b'.');
        let (ip, fpart) = match p {
            Some(p) => (&raw[..p], &raw[p + 1..]),
            None => (raw, &raw[raw.len()..]),
        };
        let nd = ip.len() + fpart.len();
        let digit = |j: usize| if j < ip.len() { ip[j] } else { fpart[j - ip.len()] };
        let mut point = ip.len() as i32 + k;
        // digits padded with zeros at the end (point > nd) or the front
        let lead = if point < 0 { (-point) as usize } else { 0 };
        if point < 0 {
            point = 0;
        }
        let total = (nd + lead).max(point as usize);
        let get = |j: usize| -> u8 {
            if j < lead {
                b'0'
            } else if j - lead < nd {
                digit(j - lead)
            } else {
                b'0'
            }
        };
        let point = point as usize;
        // integer part without leading zeros ("0" if empty)
        let mut j0 = 0;
        while j0 < point && get(j0) == b'0' {
            j0 += 1;
        }
        if j0 == point {
            body.push(b'0');
        } else {
            for j in j0..point {
                body.push(get(j));
            }
        }
        body.push(b'.');
        for j in point..total {
            body.push(get(j));
        }
    }
    if body.over {
        return None;
    }
    if d == 0 && !body.bytes().ends_with(b".") && !body.bytes().contains(&b'.') {
        body.push(b'.');
    }
    let bb = body.bytes();
    let sl = bb.len() + neg as usize;
    let drop0 = sl > w && bb.starts_with(b"0.");
    let mut s = Buf::new();
    if neg {
        s.push(b'-');
    }
    for &c in if drop0 { &bb[1..] } else { bb } {
        s.push(c);
    }
    if s.over {
        return None;
    }
    Some(s)
}

pub fn fp(x: f64, w: usize, d: usize, k: i32) -> String {
    match fp_core(x, w, d, k) {
        Some(s) => rjust_bytes(s.bytes(), w),
        None => fp_ref(x, w, d, k),
    }
}

/// Right-justify `s` into dst (stars if it does not fit).
#[inline]
fn rjust_into(s: &[u8], dst: &mut [u8]) {
    let w = dst.len();
    if s.len() > w {
        dst.fill(b'*');
    } else {
        dst[..w - s.len()].fill(b' ');
        dst[w - s.len()..].copy_from_slice(s);
    }
}

/// `fp` written into dst (w = dst.len()) without a heap string.
pub fn fp_into(x: f64, d: usize, k: i32, dst: &mut [u8]) {
    let w = dst.len();
    match fp_core(x, w, d, k) {
        Some(s) => rjust_into(s.bytes(), dst),
        None => {
            let r = fp_ref(x, w, d, k);
            rjust_into(r.as_bytes(), dst)
        }
    }
}

/// `i` written into dst (w = dst.len()) without a heap string.
pub fn i_into(v: i32, dst: &mut [u8]) {
    let mut d = [0u8; 12];
    let mut n = 0;
    let mut u = (v as i64).unsigned_abs();
    loop {
        d[n] = b'0' + (u % 10) as u8;
        n += 1;
        u /= 10;
        if u == 0 {
            break;
        }
    }
    let mut s = [0u8; 12];
    let mut m = 0;
    if v < 0 {
        s[0] = b'-';
        m = 1;
    }
    for k in (0..n).rev() {
        s[m] = d[k];
        m += 1;
    }
    rjust_into(&s[..m], dst)
}

/// `Fw.d` (scale factor 0).
pub fn f(x: f64, w: usize, d: usize) -> String {
    fp(x, w, d, 0)
}

/// `kP,Ew.d` with the default exponent field (gfortran). Same rules as
/// `ep_ref`, without heap strings.
pub fn ep(x: f64, w: usize, d: usize, k: i32) -> String {
    use std::fmt::Write;
    let neg = x.is_sign_negative();
    let ax = x.abs();
    if !ax.is_finite() {
        return ep_ref(x, w, d, k);
    }
    let nsig: usize = if k > 0 { d + 1 } else { d };
    let mut mant = Buf::new();
    let mut e10: i32;
    if ax == 0.0 {
        e10 = 0;
    } else if let Some((q, k10)) = sci_digits(ax, nsig.saturating_sub(1) + 1) {
        let n = nsig.saturating_sub(1) + 1;
        let mut dg = [0u8; 20];
        let mut v = q;
        for j in (0..n).rev() {
            dg[j] = b'0' + (v % 10) as u8;
            v /= 10;
        }
        for &c in &dg[..n] {
            mant.push(c);
        }
        e10 = k10;
    } else {
        let mut raw = Buf::new();
        if write!(raw, "{:.*e}", nsig.saturating_sub(1), ax).is_err() {
            return ep_ref(x, w, d, k);
        }
        let r = raw.bytes();
        let p = r.iter().position(|&c| c == b'e').unwrap();
        for &c in &r[..p] {
            if c != b'.' {
                mant.push(c);
            }
        }
        let (eneg, ed) = if r[p + 1] == b'-' { (true, &r[p + 2..]) } else { (false, &r[p + 1..]) };
        let mut ev = 0i32;
        for &c in ed {
            ev = ev * 10 + (c - b'0') as i32;
        }
        e10 = if eneg { -ev } else { ev };
    }
    let mut body = Buf::new();
    if ax == 0.0 {
        e10 = 0;
        if k > 0 {
            let kk = k as usize;
            if kk > d + 1 {
                return ep_ref(x, w, d, k);
            }
            for _ in 0..kk {
                body.push(b'0');
            }
            body.push(b'.');
            for _ in 0..(d + 1 - kk) {
                body.push(b'0');
            }
        } else {
            body.push(b'0');
            body.push(b'.');
            for _ in 0..d {
                body.push(b'0');
            }
        }
    } else if k > 0 {
        let kk = k as usize;
        let m = mant.bytes();
        if kk > m.len() {
            return ep_ref(x, w, d, k);
        }
        for &c in &m[..kk] {
            body.push(c);
        }
        body.push(b'.');
        for &c in &m[kk..] {
            body.push(c);
        }
        e10 = e10 - (k - 1);
    } else {
        body.push(b'0');
        body.push(b'.');
        for &c in mant.bytes() {
            body.push(c);
        }
        e10 += 1;
    }
    let mut ex = Buf::new();
    let ae = e10.unsigned_abs();
    if ae <= 99 {
        ex.push(b'E');
        ex.push(if e10 < 0 { b'-' } else { b'+' });
        ex.push(b'0' + (ae / 10) as u8);
        ex.push(b'0' + (ae % 10) as u8);
    } else if ae <= 999 {
        ex.push(if e10 < 0 { b'-' } else { b'+' });
        ex.push(b'0' + (ae / 100) as u8);
        ex.push(b'0' + (ae / 10 % 10) as u8);
        ex.push(b'0' + (ae % 10) as u8);
    } else {
        return stars(w);
    }
    if body.over || mant.over {
        return ep_ref(x, w, d, k);
    }
    let bb = body.bytes();
    let sl = neg as usize + bb.len() + ex.n;
    let drop0 = sl > w && bb.starts_with(b"0.");
    let mut s = Buf::new();
    if neg {
        s.push(b'-');
    }
    for &c in if drop0 { &bb[1..] } else { bb } {
        s.push(c);
    }
    for &c in ex.bytes() {
        s.push(c);
    }
    if s.over {
        return ep_ref(x, w, d, k);
    }
    rjust_bytes(s.bytes(), w)
}

/// `Ew.d` (scale factor 0).
pub fn e(x: f64, w: usize, d: usize) -> String {
    ep(x, w, d, 0)
}

/// `1P,Ew.d`, the common NJOY form.
pub fn e1p(x: f64, w: usize, d: usize) -> String {
    ep(x, w, d, 1)
}

// Reference versions (the code before the allocation-free rewrite, checked
// against gfortran). The fast versions above fall back to them for fields
// longer than the stack buffer, and `examples/ffmt_check.rs` compares the two.

fn stars_ref(w: usize) -> String {
    "*".repeat(w)
}

fn rjust_ref(s: &str, w: usize) -> String {
    if s.len() > w {
        stars_ref(w)
    } else {
        format!("{:>w$}", s, w = w)
    }
}

/// `Iw`.
pub fn i_ref(v: i32, w: usize) -> String {
    rjust_ref(&v.to_string(), w)
}

/// `I0`.
pub fn i0_ref(v: i32) -> String {
    v.to_string()
}

/// `Aw` output of a string (Fortran semantics: if the string is longer than w,
/// the leftmost w characters are written; if shorter, it is right-justified).
pub fn a_ref(s: &str, w: usize) -> String {
    if s.len() >= w {
        s[..w].to_string()
    } else {
        format!("{:>w$}", s, w = w)
    }
}

/// Shift the decimal point of a plain decimal string ("123.4567") right by k.
fn shift_point_ref(s: &str, k: i32) -> String {
    if k == 0 {
        return s.to_string();
    }
    let (ip, fp) = match s.find('.') {
        Some(p) => (&s[..p], &s[p + 1..]),
        None => (s, ""),
    };
    let mut digits: Vec<u8> = ip.bytes().chain(fp.bytes()).collect();
    let mut point = ip.len() as i32 + k;
    while point > digits.len() as i32 {
        digits.push(b'0');
    }
    while point < 0 {
        digits.insert(0, b'0');
        point += 1;
    }
    let point = point as usize;
    let mut ipart: String = String::from_utf8(digits[..point].to_vec()).unwrap();
    let fpart: String = String::from_utf8(digits[point..].to_vec()).unwrap();
    let trimmed = ipart.trim_start_matches('0').to_string();
    ipart = if trimmed.is_empty() { "0".to_string() } else { trimmed };
    if fpart.is_empty() {
        format!("{}.", ipart)
    } else {
        format!("{}.{}", ipart, fpart)
    }
}

/// `kP,Fw.d`.
pub fn fp_ref(x: f64, w: usize, d: usize, k: i32) -> String {
    let neg = x.is_sign_negative();
    let ax = x.abs();
    let dd = (d as i32 + k).max(0) as usize;
    let body = format!("{:.*}", dd, ax);
    let mut body = shift_point_ref(&body, k);
    if d == 0 && !body.ends_with('.') {
        // shift_point always leaves a point; with d=0 Rust gives no point.
        if !body.contains('.') {
            body.push('.');
        }
    }
    // gfortran prints a minus sign for negative values, including -0 and
    // values that round to zero.
    let sign = if neg { "-" } else { "" };
    let mut s = format!("{}{}", sign, body);
    if s.len() > w && body.starts_with("0.") {
        body = body[1..].to_string();
        s = format!("{}{}", sign, body);
    }
    rjust_ref(&s, w)
}

/// `Fw.d` (scale factor 0).
pub fn f_ref(x: f64, w: usize, d: usize) -> String {
    fp_ref(x, w, d, 0)
}

/// `kP,Ew.d` with the default exponent field (gfortran).
pub fn ep_ref(x: f64, w: usize, d: usize, k: i32) -> String {
    let neg = x.is_sign_negative();
    let ax = x.abs();
    // number of significant digits
    let nsig: usize = if k > 0 { d + 1 } else { d };
    let (mant, mut e10) = if ax == 0.0 {
        ("0".repeat(nsig.max(1)), 0i32)
    } else {
        let s = format!("{:.*e}", nsig.saturating_sub(1), ax);
        let p = s.find('e').unwrap();
        let m: String = s[..p].chars().filter(|c| *c != '.').collect();
        let e: i32 = s[p + 1..].parse().unwrap();
        (m, e)
    };
    let body;
    if ax == 0.0 {
        e10 = 0;
        if k > 0 {
            let kk = k as usize;
            body = format!("{}.{}", &"0".repeat(kk), "0".repeat(d + 1 - kk));
        } else {
            body = format!("0.{}", "0".repeat(d));
        }
    } else if k > 0 {
        let kk = k as usize;
        body = format!("{}.{}", &mant[..kk], &mant[kk..]);
        e10 = e10 - (k - 1);
    } else {
        // k == 0: 0.dddd
        body = format!("0.{}", mant);
        e10 += 1;
    }
    let exp = if e10.abs() <= 99 {
        format!("E{}{:02}", if e10 < 0 { '-' } else { '+' }, e10.abs())
    } else if e10.abs() <= 999 {
        format!("{}{:03}", if e10 < 0 { '-' } else { '+' }, e10.abs())
    } else {
        return stars_ref(w);
    };
    let sign = if neg { "-" } else { "" };
    let mut s = format!("{}{}{}", sign, body, exp);
    if s.len() > w && body.starts_with("0.") {
        s = format!("{}{}{}", sign, &body[1..], exp);
    }
    rjust_ref(&s, w)
}

/// `Ew.d` (scale factor 0).
pub fn e_ref(x: f64, w: usize, d: usize) -> String {
    ep_ref(x, w, d, 0)
}

/// `1P,Ew.d`, the common NJOY form.
pub fn e1p_ref(x: f64, w: usize, d: usize) -> String {
    ep_ref(x, w, d, 1)
}

// ----------------------------------------------------------------- input --

/// Return bytes [start, start+w) of a record, padded with blanks.
pub fn field(rec: &[u8], start: usize, w: usize) -> Vec<u8> {
    let mut v = Vec::with_capacity(w);
    for j in start..start + w {
        v.push(if j < rec.len() { rec[j] } else { b' ' });
    }
    v
}

#[derive(Debug)]
pub struct ReadError(pub String);

/// A field of a record, blank padded, held on the stack (for the numeric
/// fields of ENDF lines; `field` is the heap version).
pub struct Fixed {
    b: [u8; 128],
    n: usize,
}

impl std::ops::Deref for Fixed {
    type Target = [u8];
    #[inline]
    fn deref(&self) -> &[u8] {
        &self.b[..self.n]
    }
}

/// `field` without heap allocation (w <= 128).
#[inline]
pub fn fieldx(rec: &[u8], start: usize, w: usize) -> Fixed {
    assert!(w <= 128);
    let mut f = Fixed { b: [b' '; 128], n: w };
    if start < rec.len() {
        let e = (start + w).min(rec.len());
        f.b[..e - start].copy_from_slice(&rec[start..e]);
    }
    f
}

/// `Iw` input (blank-null). Fast path for plain integers; anything else is
/// handed to `read_i_ref`, so results and error messages are the same.
pub fn read_i(fld: &[u8]) -> Result<i32, ReadError> {
    let mut neg = false;
    let mut v: i64 = 0;
    let mut nd = 0;
    let mut seen = false;
    for &c in fld {
        if c == b' ' || c == b'\t' {
            continue;
        }
        if !seen && (c == b'+' || c == b'-') {
            neg = c == b'-';
            seen = true;
            continue;
        }
        seen = true;
        if c.is_ascii_digit() && nd < 18 {
            v = v * 10 + (c - b'0') as i64;
            nd += 1;
        } else {
            return read_i_ref(fld);
        }
    }
    if !seen {
        return Ok(0);
    }
    if nd == 0 {
        return read_i_ref(fld);
    }
    Ok((if neg { -v } else { v }) as i32)
}

/// `Fw.d` / `Ew.d` input (blank-null). Builds the same decimal text as
/// `read_f_ref` on the stack and parses it the same way; any field that is
/// not a plain number goes to `read_f_ref`.
pub fn read_f(fld: &[u8], d: usize) -> Result<f64, ReadError> {
    let mut s = [0u8; 128];
    let mut n = 0;
    for &c in fld {
        if c != b' ' && c != b'\t' {
            if n == s.len() {
                return read_f_ref(fld, d);
            }
            s[n] = c;
            n += 1;
        }
    }
    let s = &s[..n];
    if s.is_empty() {
        return Ok(0.0);
    }
    let mut i = 0;
    let mut neg = false;
    if s[i] == b'+' || s[i] == b'-' {
        neg = s[i] == b'-';
        i += 1;
    }
    let ip0 = i;
    while i < s.len() && s[i].is_ascii_digit() {
        i += 1;
    }
    let ipart = &s[ip0..i];
    let mut seen_point = false;
    let mut fpart: &[u8] = &s[i..i];
    if i < s.len() && s[i] == b'.' {
        seen_point = true;
        i += 1;
        let f0 = i;
        while i < s.len() && s[i].is_ascii_digit() {
            i += 1;
        }
        fpart = &s[f0..i];
    }
    if ipart.is_empty() && fpart.is_empty() && !seen_point {
        return read_f_ref(fld, d);
    }
    let mut expo: i64 = 0;
    if i < s.len() {
        let c = s[i];
        let mut j = i;
        let letter = matches!(c, b'e' | b'E' | b'd' | b'D' | b'q' | b'Q');
        if letter {
            j += 1;
        }
        let mut eneg = false;
        if j < s.len() && (s[j] == b'+' || s[j] == b'-') {
            eneg = s[j] == b'-';
            j += 1;
        } else if !letter {
            return read_f_ref(fld, d);
        }
        if j >= s.len() || s.len() - j > 9 {
            return read_f_ref(fld, d);
        }
        let mut ev: i64 = 0;
        while j < s.len() {
            if !s[j].is_ascii_digit() {
                return read_f_ref(fld, d);
            }
            ev = ev * 10 + (s[j] - b'0') as i64;
            j += 1;
        }
        expo = if eneg { -ev } else { ev };
    }
    if !seen_point && d > 0 {
        expo -= d as i64;
    }
    // text: [-]ipart[.fpart]e<expo>, as read_f_ref builds it
    let mut t = [0u8; 200];
    let mut m = 0;
    if neg {
        t[m] = b'-';
        m += 1;
    }
    if ipart.is_empty() {
        t[m] = b'0';
        m += 1;
    } else {
        t[m..m + ipart.len()].copy_from_slice(ipart);
        m += ipart.len();
    }
    if !fpart.is_empty() {
        t[m] = b'.';
        m += 1;
        t[m..m + fpart.len()].copy_from_slice(fpart);
        m += fpart.len();
    }
    t[m] = b'e';
    m += 1;
    let mut ed = [0u8; 24];
    let mut ne = 0;
    let neg_e = expo < 0;
    let mut u = expo.unsigned_abs();
    loop {
        ed[ne] = b'0' + (u % 10) as u8;
        ne += 1;
        u /= 10;
        if u == 0 {
            break;
        }
    }
    if neg_e {
        t[m] = b'-';
        m += 1;
    }
    for k in (0..ne).rev() {
        t[m] = ed[k];
        m += 1;
    }
    match std::str::from_utf8(&t[..m]).ok().and_then(|x| x.parse::<f64>().ok()) {
        Some(v) => Ok(v),
        None => read_f_ref(fld, d),
    }
}

/// `Iw` input (blank-null).
pub fn read_i_ref(fld: &[u8]) -> Result<i32, ReadError> {
    let s: String = fld
        .iter()
        .filter(|c| **c != b' ' && **c != b'\t')
        .map(|c| *c as char)
        .collect();
    if s.is_empty() {
        return Ok(0);
    }
    let t = s.strip_prefix('+').unwrap_or(&s);
    t.parse::<i64>()
        .map(|v| v as i32)
        .map_err(|_| ReadError(format!("Bad value during integer read: '{}'", s)))
}

/// `Fw.d` / `Ew.d` input (blank-null), for the given d.
pub fn read_f_ref(fld: &[u8], d: usize) -> Result<f64, ReadError> {
    let s: Vec<u8> = fld
        .iter()
        .copied()
        .filter(|c| *c != b' ' && *c != b'\t')
        .collect();
    if s.is_empty() {
        return Ok(0.0);
    }
    let bad = || ReadError(format!("Bad value during floating point read: '{}'", String::from_utf8_lossy(fld)));
    let mut i = 0;
    let mut neg = false;
    if s[i] == b'+' || s[i] == b'-' {
        neg = s[i] == b'-';
        i += 1;
    }
    let mut ipart = String::new();
    let mut fpart = String::new();
    let mut seen_point = false;
    let mut ndig = 0;
    while i < s.len() {
        let c = s[i];
        if c.is_ascii_digit() {
            if seen_point { fpart.push(c as char) } else { ipart.push(c as char) }
            ndig += 1;
            i += 1;
        } else if c == b'.' && !seen_point {
            seen_point = true;
            i += 1;
        } else {
            break;
        }
    }
    // gfortran accepts special values
    if ndig == 0 && !seen_point {
        let rest: String = s[i..].iter().map(|c| (*c as char).to_ascii_lowercase()).collect();
        if rest == "inf" || rest == "infinity" {
            return Ok(if neg { f64::NEG_INFINITY } else { f64::INFINITY });
        }
        if rest.starts_with("nan") {
            return Ok(f64::NAN);
        }
        return Err(bad());
    }
    let mut expo: i64 = 0;
    if i < s.len() {
        let c = s[i];
        let mut j = i;
        if matches!(c, b'e' | b'E' | b'd' | b'D' | b'q' | b'Q') {
            j += 1;
            if j >= s.len() {
                // gfortran treats a missing exponent after the letter as 0
                return Err(bad());
            }
        }
        let mut eneg = false;
        if j < s.len() && (s[j] == b'+' || s[j] == b'-') {
            eneg = s[j] == b'-';
            j += 1;
        } else if !matches!(c, b'e' | b'E' | b'd' | b'D' | b'q' | b'Q') {
            return Err(bad());
        }
        if j >= s.len() {
            return Err(bad());
        }
        let mut ev: i64 = 0;
        while j < s.len() {
            if !s[j].is_ascii_digit() {
                return Err(bad());
            }
            ev = ev.saturating_mul(10).saturating_add((s[j] - b'0') as i64);
            j += 1;
        }
        expo = if eneg { -ev } else { ev };
    } else {
        // no exponent
    }
    // implied decimal point
    if !seen_point && d > 0 {
        expo -= d as i64;
    }
    if ipart.is_empty() {
        ipart.push('0');
    }
    let txt = format!(
        "{}{}{}e{}",
        if neg { "-" } else { "" },
        ipart,
        if fpart.is_empty() { String::new() } else { format!(".{}", fpart) },
        expo
    );
    txt.parse::<f64>().map_err(|_| bad())
}
