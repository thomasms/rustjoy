//! Fortran formatted-I/O editing, matching libgfortran.
//!
//! Output: `Iw`, `Fw.d` (with `kP` scale factor), `Ew.d` (with `kP`), `Aw`.
//! Input:  `Iw`, `Fw.d`/`Ew.d` (blank-null mode, which is gfortran's default).
//!
//! All decimal conversions are exact and round-half-even on exact ties, which
//! is what both libgfortran (via glibc printf) and Rust's formatter do.

#![allow(dead_code)]

fn stars(w: usize) -> String {
    "*".repeat(w)
}

fn rjust(s: &str, w: usize) -> String {
    if s.len() > w {
        stars(w)
    } else {
        format!("{:>w$}", s, w = w)
    }
}

/// `Iw`.
pub fn i(v: i32, w: usize) -> String {
    rjust(&v.to_string(), w)
}

/// `I0`.
pub fn i0(v: i32) -> String {
    v.to_string()
}

/// `Aw` output of a string (Fortran semantics: if the string is longer than w,
/// the leftmost w characters are written; if shorter, it is right-justified).
pub fn a(s: &str, w: usize) -> String {
    if s.len() >= w {
        s[..w].to_string()
    } else {
        format!("{:>w$}", s, w = w)
    }
}

/// Shift the decimal point of a plain decimal string ("123.4567") right by k.
fn shift_point(s: &str, k: i32) -> String {
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
pub fn fp(x: f64, w: usize, d: usize, k: i32) -> String {
    let neg = x.is_sign_negative();
    let ax = x.abs();
    let dd = (d as i32 + k).max(0) as usize;
    let body = format!("{:.*}", dd, ax);
    let mut body = shift_point(&body, k);
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
    rjust(&s, w)
}

/// `Fw.d` (scale factor 0).
pub fn f(x: f64, w: usize, d: usize) -> String {
    fp(x, w, d, 0)
}

/// `kP,Ew.d` with the default exponent field (gfortran).
pub fn ep(x: f64, w: usize, d: usize, k: i32) -> String {
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
        return stars(w);
    };
    let sign = if neg { "-" } else { "" };
    let mut s = format!("{}{}{}", sign, body, exp);
    if s.len() > w && body.starts_with("0.") {
        s = format!("{}{}{}", sign, &body[1..], exp);
    }
    rjust(&s, w)
}

/// `Ew.d` (scale factor 0).
pub fn e(x: f64, w: usize, d: usize) -> String {
    ep(x, w, d, 0)
}

/// `1P,Ew.d`, the common NJOY form.
pub fn e1p(x: f64, w: usize, d: usize) -> String {
    ep(x, w, d, 1)
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

/// `Iw` input (blank-null).
pub fn read_i(fld: &[u8]) -> Result<i32, ReadError> {
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
pub fn read_f(fld: &[u8], d: usize) -> Result<f64, ReadError> {
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
