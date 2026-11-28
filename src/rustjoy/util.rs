//! Port of NJOY2016 `util.f90` and `phys.f90`.

#![allow(dead_code)]

use crate::endf::Ctx;
use crate::farr::F;
use crate::ffmt;
use crate::fio::{RecW, Status};
use crate::fmath::{int, log10, nint, powi, ipow};

pub mod physics {
    pub const PI: f64 = 3.141592653589793238e0;
    pub const EULER: f64 = 0.57721566490153286e0;
    pub const BK: f64 = 8.617333262e-5;
    pub const EV: f64 = 1.602176634e-12;
    pub const CLIGHT: f64 = 2.99792458e10;
    pub const AMU: f64 = 931.49410242e6 * EV / (CLIGHT * CLIGHT);
    pub const HBAR: f64 = 6.582119569e-16 * EV;
    pub const FINSTRI: f64 = 1.0e16 * HBAR / (EV * EV * CLIGHT);
    pub const AMASSN: f64 = 1.00866491595e0;
    pub const AMASSP: f64 = 1.007276466621e0;
    pub const AMASSD: f64 = 2.013553212745e0;
    pub const AMASST: f64 = 3.01550071621e0;
    pub const AMASSH: f64 = 3.014932247175e0;
    pub const AMASSA: f64 = 4.001506179127e0;
    pub const AMASSE: f64 = 5.48579909065e-4;
    pub const PNRATIO: f64 = AMASSP / AMASSN;
    pub const DNRATIO: f64 = AMASSD / AMASSN;
    pub const TNRATIO: f64 = AMASST / AMASSN;
    pub const HNRATIO: f64 = AMASSH / AMASSN;
    pub const ANRATIO: f64 = AMASSA / AMASSN;
    pub const EPAIR: f64 = AMASSE * AMU * CLIGHT * CLIGHT / EV;
}

/// `sigfig(x,ndig,idig)`: round x to ndig significant figures and add idig in
/// the last place.
pub fn sigfig(x: f64, ndig: i32, idig: i32) -> f64 {
    let bias = 1.0000000000001e0f64;
    let zero = 0.0f64;
    let ten = 10.0f64;
    let mut xx = 0.0f64;
    if x != zero {
        let ipwr0 = match floor_log10(x.abs()) {
            Some(k) => k,
            None => return sigfig_ref(x, ndig, idig),
        };
        let mut ipwr = ndig - 1 - ipwr0;
        let p = |n: i32| pow10_tab(n).unwrap_or_else(|| powi(ten, n));
        let mut ii = nint(x * p(ipwr) + p(ndig - 11));
        const IP10: [i32; 10] = [1, 10, 100, 1000, 10_000, 100_000, 1_000_000, 10_000_000, 100_000_000, 1_000_000_000];
        let top = if (0..10).contains(&ndig) { IP10[ndig as usize] } else { ipow(10, ndig) };
        if ii >= top {
            ii /= 10;
            ipwr -= 1;
        }
        ii += idig;
        xx = ii as f64 * p(-ipwr);
    }
    xx = xx * bias;
    xx
}

/// The `ipwr` of `sigfig` before `ndig` (`int(log10(ax))`, less 1 when
/// the log is negative) for finite, normal ax > 0, without the log. It is
/// floor(log10(ax)), found from a table of powers of ten. None when ax is
/// within a relative 1e-9 of a power of ten (where the rounded log10 could
/// be the integer itself) or is not normal: the caller then uses the log.
/// Away from the powers of ten, the log10 of glibc (error below 1 ulp, so
/// below 1e-13 here) cannot cross an integer, so both give the same value.
#[inline]
fn floor_log10(ax: f64) -> Option<i32> {
    if !(ax.is_normal()) {
        return None;
    }
    let t = tables();
    let e2 = ((ax.to_bits() >> 52) & 0x7ff) as i32 - 1023;
    // e2*log10(2) <= log10(ax) < (e2+1)*log10(2), which spans less than 1
    // floor(e2*log10(2)) for |e2| < 1700 (78913/2^18 is log10(2) to 1e-7)
    let mut k = (e2 * 78913) >> 18;
    let idx = |k: i32| (k + P10_OFF) as usize;
    if k < -P10_OFF || k + 1 > P10_OFF {
        return None;
    }
    if ax >= t.p10[idx(k + 1)] {
        k += 1;
        if k + 1 > P10_OFF {
            return None;
        }
    }
    let lo = t.p10[idx(k)];
    let hi = t.p10[idx(k + 1)];
    if ax < lo || ax >= hi {
        return None;
    }
    const REL: f64 = 1.0e-9;
    if ax - lo <= REL * lo || hi - ax <= REL * hi {
        return None;
    }
    Some(k)
}

const P10_OFF: i32 = 330;
const POW_OFF: i32 = 400;

struct Tabs {
    /// 10^k, correctly rounded, for k = -P10_OFF..=P10_OFF
    p10: Vec<f64>,
    /// powi(10, n) (the Fortran `ten**n`), for n = -POW_OFF..=POW_OFF
    pw: Vec<f64>,
}

fn tables() -> &'static Tabs {
    static T: std::sync::OnceLock<Tabs> = std::sync::OnceLock::new();
    T.get_or_init(|| Tabs {
        p10: (-P10_OFF..=P10_OFF).map(|k| format!("1e{}", k).parse::<f64>().unwrap()).collect(),
        pw: (-POW_OFF..=POW_OFF).map(|n| powi(10.0, n)).collect(),
    })
}

/// powi(10, n) from the table (the same bits), or None outside it.
#[inline]
fn pow10_tab(n: i32) -> Option<f64> {
    if (-POW_OFF..=POW_OFF).contains(&n) {
        Some(tables().pw[(n + POW_OFF) as usize])
    } else {
        None
    }
}

/// The original `sigfig`, as in the Fortran (log10 and `ten**n`). It is the
/// reference for `sigfig` (see tools/checks/ffmt_check.rs) and its path
/// for the inputs that the table does not cover.
pub fn sigfig_ref(x: f64, ndig: i32, idig: i32) -> f64 {
    let bias = 1.0000000000001e0f64;
    let zero = 0.0f64;
    let ten = 10.0f64;
    let mut xx = 0.0f64;
    if x != zero {
        let aa = log10(x.abs());
        let mut ipwr = int(aa);
        if aa < zero {
            ipwr -= 1;
        }
        ipwr = ndig - 1 - ipwr;
        let mut ii = nint(x * powi(ten, ipwr) + powi(ten, ndig - 11));
        if ii >= ipow(10, ndig) {
            ii /= 10;
            ipwr -= 1;
        }
        ii += idig;
        xx = ii as f64 * powi(ten, -ipwr);
    }
    xx = xx * bias;
    xx
}

/// `a10`: 10-character representation of a real.
pub fn a10(x: f64) -> String {
    let zero = 0.0f64;
    let one = 1.0f64;
    let ten = 10.0f64;
    let up = 0.000000001e0f64;
    let top6 = 9.999995e0f64;
    let top5 = 9.99995e0f64;
    let top4 = 9.9995e0f64;
    let f;
    let mut s;
    let mut n;
    if x == zero {
        f = 0.0;
        s = '+';
        n = 0;
    } else {
        let xx = x;
        n = int(log10(xx.abs()));
        let mut ff;
        if xx.abs() < one {
            n = 1 - n;
            ff = xx * powi(ten, n);
            s = '-';
            if (n.abs() < 10 && ff.abs() < top6)
                || (n.abs() < 100 && n.abs() >= 10 && ff.abs() < top5)
                || (n.abs() >= 100 && ff.abs() < top4)
            {
            } else {
                ff = ff / ten;
                ff = ff + up;
                n = n - 1;
                if n <= 0 {
                    s = '+';
                }
            }
        } else {
            ff = xx / powi(ten, n);
            s = '+';
            if (n.abs() < 10 && ff.abs() < top6)
                || (n.abs() < 100 && n.abs() >= 10 && ff.abs() < top5)
                || (n.abs() >= 100 && ff.abs() < top4)
            {
            } else {
                ff = ff / ten;
                ff = ff + up;
                n = n + 1;
            }
        }
        f = ff;
    }
    let r = if n.abs() < 10 {
        format!("{}{}{}", ffmt::f(f, 8, 5), s, ffmt::i(n, 1))
    } else if n.abs() < 100 {
        format!("{}{}{}", ffmt::f(f, 7, 4), s, ffmt::i(n, 2))
    } else {
        format!("{}{}{}", ffmt::f(f, 6, 3), s, ffmt::i(n, 3))
    };
    format!("{:<10}", r)
}

#[repr(C)]
struct Tm {
    tm_sec: i32,
    tm_min: i32,
    tm_hour: i32,
    tm_mday: i32,
    tm_mon: i32,
    tm_year: i32,
    tm_wday: i32,
    tm_yday: i32,
    tm_isdst: i32,
    tm_gmtoff: i64,
    tm_zone: *const u8,
}

extern "C" {
    fn time(t: *mut i64) -> i64;
    fn localtime_r(t: *const i64, tm: *mut Tm) -> *mut Tm;
}

fn now_tm() -> Tm {
    unsafe {
        let mut t: i64 = 0;
        time(&mut t);
        let mut tm: Tm = std::mem::zeroed();
        localtime_r(&t, &mut tm);
        tm
    }
}

fn zero_fill(s: String) -> String {
    s.replace(' ', "0")
}

/// `dater`: mm/dd/yy
pub fn dater() -> String {
    let tm = now_tm();
    zero_fill(format!(
        "{}/{}/{}",
        ffmt::i(tm.tm_mon + 1, 2),
        ffmt::i(tm.tm_mday, 2),
        ffmt::i((tm.tm_year + 1900) % 100, 2)
    ))
}

/// `wclock`: hh:mm:ss
pub fn wclock() -> String {
    let tm = now_tm();
    zero_fill(format!(
        "{}:{}:{}",
        ffmt::i(tm.tm_hour, 2),
        ffmt::i(tm.tm_min, 2),
        ffmt::i(tm.tm_sec, 2)
    ))
}

extern "C" {
    fn clock() -> i64;
}

impl Ctx {
    /// `error`: print message and `stop 77`.
    pub fn error(&mut self, from: &str, mess1: &str, mess2: &str) -> ! {
        if self.io.quiet {
            // worker thread: give up this value; the main program computes
            // it again and stops there with the real message
            self.io.touched = true;
            std::panic::resume_unwind(Box::new(crate::fio::WorkerAbort));
        }
        let m1 = mess1.trim_end();
        let m2 = mess2.trim_end();
        self.io.so(&format!("\n ***error in {}***{}", from, m1));
        if !m2.is_empty() {
            self.io.so(&format!("{}{}", " ".repeat(22), m2));
        }
        self.io.so(" ");
        self.io.se(&format!("\n ***error in {}***{}", from, m1));
        if !m2.is_empty() {
            self.io.se(&format!("{}{}", " ".repeat(22), m2));
        }
        self.io.se(" ");
        self.io.flush_all();
        eprintln!("STOP 77");
        std::process::exit(77);
    }

    /// `mess`: print a message and continue.
    pub fn mess(&mut self, from: &str, mess1: &str, mess2: &str) {
        let m1 = mess1.trim_end();
        let m2 = mess2.trim_end();
        self.io.so(&format!("\n ---message from {}---{}", from, m1));
        if !m2.is_empty() {
            self.io.so(&format!("{}{}", " ".repeat(26), m2));
        }
        self.io.se(&format!("\n ---message from {}---{}", from, m1));
        if !m2.is_empty() {
            self.io.se(&format!("{}{}", " ".repeat(26), m2));
        }
    }

    /// `timer`: processor time in seconds.
    pub fn timer(&self) -> f64 {
        let c = unsafe { clock() };
        c as f64 / 1.0e6
    }

    pub fn repoz(&mut self, ntape: i32) {
        if ntape == 0 {
            return;
        }
        self.io.rewind(ntape.abs());
    }

    pub fn skiprz(&mut self, nunit: i32, nrt: i32) {
        if nunit == 0 || nrt == 0 {
            return;
        }
        let nun = nunit.abs();
        let nr = nrt.abs();
        if nrt > 0 {
            for _ in 0..nr {
                if nunit > 0 {
                    let _ = self.io.read_line(nun);
                }
                if nunit < 0 {
                    let _ = self.io.read_rec(nun);
                }
            }
        } else {
            for _ in 0..nr {
                self.io.backspace(nun);
            }
        }
    }

    pub fn openz(&mut self, lun: i32, new: i32) {
        let nun = lun.abs();
        if (5..=7).contains(&nun) || nun == 0 {
            return;
        }
        if nun > 999 {
            self.error("openz", "illegal unit number.", " ");
        }
        let fname = format!("tape{}", nun);
        let formatted = lun >= 0;
        if (10..=19).contains(&nun) {
            self.io.open(nun, formatted, Status::Scratch, None);
        } else if new != 1 {
            self.io.open(nun, formatted, Status::Old, Some(&fname));
        } else {
            let _ = std::fs::remove_file(&fname);
            self.io.open(nun, formatted, Status::New, Some(&fname));
        }
    }

    pub fn closz(&mut self, lun: i32) {
        let nun = lun.abs();
        if nun < 10 {
            return;
        }
        if nun > 999 {
            self.error("closz", "illegal unit number.", " ");
        }
        self.io.close(nun, false);
    }

    /// `loada`: buffered write of na words to a scratch tape.
    pub fn loada(&mut self, i: i32, a: &[f64], na: i32, ntape: i32, buf: &mut [f64], nbuf: i32) {
        let nl = nbuf / na;
        let j = i.abs();
        let mut x = j as f64;
        x = x / nl as f64;
        let ix = j / nl;
        let xnow = (x - ix as f64) * nl as f64;
        let mut inow = nint(xnow);
        if inow == 0 {
            inow = nl;
        }
        let mut k = na * (inow - 1);
        for jj in 0..na as usize {
            k += 1;
            buf[(k - 1) as usize] = a[jj];
        }
        if i == 1 {
            self.repoz(-ntape);
        }
        if inow == nl || i < 0 {
            let mut r = RecW::new();
            r.rs(&buf[..nbuf as usize]);
            self.io.write_rec(ntape, r);
        }
    }

    /// `finda`: buffered read of na words from a scratch tape.
    pub fn finda(&mut self, i: i32, a: &mut [f64], na: i32, ntape: i32, buf: &mut [f64], nbuf: i32) {
        let nl = nbuf / na;
        let mut inow = i % nl;
        if inow == 0 {
            inow = nl;
        }
        if i == 1 {
            self.repoz(-ntape);
        }
        if inow == 1 {
            self.read_buf(ntape, buf, nbuf);
        }
        let mut k = na * (inow - 1);
        for jj in 0..na as usize {
            k += 1;
            a[jj] = buf[(k - 1) as usize];
        }
    }

    fn read_buf(&mut self, ntape: i32, buf: &mut [f64], nbuf: i32) {
        self.io.read_rec_with(ntape, |r| {
            for k in 0..nbuf as usize {
                buf[k] = r.r();
            }
        });
    }

    /// `scana`: find the position of e in a buffered scratch tape.
    pub fn scana(&mut self, e: f64, ip: &mut i32, np: i32, na: i32, ntape: i32, buf: &mut [f64], nbuf: i32) {
        if *ip != 0 {
            self.error("scana", "initial ip ne 0.", " ");
        }
        let nl = nbuf / na;
        let mut j = 1;
        let mut last = 0;
        let b = |buf: &[f64], k: i32| buf[(k - 1) as usize];
        let mut lbl = 110;
        loop {
            match lbl {
                110 => {
                    let ef = b(buf, 1);
                    let mut ibufl = na * (nl - 1) + 1;
                    if np < *ip + nl {
                        ibufl = (np - *ip - 1) * na + 1;
                    }
                    let el = b(buf, ibufl);
                    if e < ef && j == 1 {
                        // 130
                        *ip = 1;
                        return;
                    }
                    if e < ef && j > 1 {
                        // 140
                        self.io.backspace(ntape);
                        self.io.backspace(ntape);
                        self.read_buf(ntape, buf, nbuf);
                        *ip -= 1;
                        return;
                    }
                    if e >= ef && e < el {
                        lbl = 120;
                        continue;
                    }
                    if last > 0 {
                        let s = format!("did not find energy {}", ffmt::e1p(e, 12, 4));
                        self.error("scana", &s, " ");
                    }
                    self.read_buf(ntape, buf, nbuf);
                    j += 1;
                    *ip += nl;
                    if np <= *ip + nl {
                        last = 1;
                    }
                }
                120 => {
                    *ip += 1;
                    let jp = *ip - nl * (j - 1);
                    let indx = na * (jp - 1) + 1;
                    let en = b(buf, indx);
                    if en == e {
                        lbl = 145;
                        continue;
                    }
                    if en < e {
                        continue;
                    }
                    *ip -= 1;
                    lbl = 145;
                }
                _ => {
                    // 145
                    let ipn = *ip % nl;
                    if ipn != 1 {
                        return;
                    }
                    self.io.backspace(ntape);
                    if *ip < 1 {
                        *ip = 1;
                    }
                    return;
                }
            }
        }
    }
}

/// Helper: view of a 1-based slice as F.
pub fn fv(v: &mut [f64]) -> F<'_> {
    F(v)
}

/// Number of threads for the parallel prefetch (RECONR `resxs`, BROADR
/// `broadn`): the environment variable `NJOY_THREADS` if set (1 turns the
/// prefetch off), else the number of CPUs. The output does not depend on it.
pub fn threads() -> usize {
    static N: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
    *N.get_or_init(|| {
        if let Ok(v) = std::env::var("NJOY_THREADS") {
            if let Ok(n) = v.trim().parse::<usize>() {
                return n.max(1);
            }
        }
        std::thread::available_parallelism().map(|n| n.get()).unwrap_or(1)
    })
}

/// `NJOY_PF_STATS` set: print prefetch hit/miss counts to standard error.
pub fn pf_stats() -> bool {
    static S: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *S.get_or_init(|| std::env::var_os("NJOY_PF_STATS").is_some())
}
