//! Port of NJOY2016 `endf.f90` (ENDF/PENDF tape I/O and interpolation), plus
//! the run-time context (`Ctx`) that holds the I/O system and the `endf`
//! module variables.
//!
//! Unit-number conventions are those of NJOY: a positive unit is a formatted
//! (ASCII) tape, a negative unit is an unformatted (binary) tape, 0 is none.

#![allow(clippy::too_many_arguments)]

use crate::farr::F;
use crate::ffmt;
use crate::fio::{Io, RecW};
use crate::fmath::{exp, log, log10, nint, int, powi, pow, sqrt};

pub struct Ctx {
    pub io: Box<Io>,
    // endf module variables
    pub npage: i32,
    pub iverf: i32,
    pub c1h: f64,
    pub c2h: f64,
    pub l1h: i32,
    pub l2h: i32,
    pub n1h: i32,
    pub n2h: i32,
    pub math: i32,
    pub mfh: i32,
    pub mth: i32,
    pub nsh: i32,
    pub nsp: i32,
    pub nsc: i32,
    pub thr6: f64,
    nah: i32,
    ah: [f64; 6],
    ihol: i32,
    /// `rhol(17)` real array; `hol(17)` character(4) is its first 68 bytes.
    rhol: [u8; 136],
    nr: i32,
    nbt: [i32; 24],
    jnt: [i32; 24],
    /// contio's local `field(2)` (uninitialised in Fortran; kept here)
    cfield: [[u8; 11]; 2],
    g1: GetySave,
    g2: GetySave,
    pub t0: std::time::Instant,
}

#[derive(Default, Clone, Copy)]
struct GetySave {
    nwtot: i32,
    nr: i32,
    np: i32,
    lt: i32,
    ip1: i32,
    ip2: i32,
    ip: i32,
    ir: i32,
    nb: i32,
    nbt: i32,
    xlast: f64,
    ylast: f64,
}

#[inline]
fn hrec_fields(line: &[u8]) -> (i32, i32, i32, i32) {
    // '(i4,i2,i3,i5)' starting at column 67
    let math = rd_i(&ffmt::field(line, 66, 4));
    let mfh = rd_i(&ffmt::field(line, 70, 2));
    let mth = rd_i(&ffmt::field(line, 72, 3));
    let nsp = rd_i(&ffmt::field(line, 75, 5));
    (math, mfh, mth, nsp)
}

pub fn rd_i(f: &[u8]) -> i32 {
    match ffmt::read_i(f) {
        Ok(v) => v,
        Err(e) => crate::fio::fatal_runtime(&format!("Fortran runtime error: {}", e.0)),
    }
}

pub fn rd_e(f: &[u8]) -> f64 {
    match ffmt::read_f(f, 0) {
        Ok(v) => v,
        Err(e) => crate::fio::fatal_runtime(&format!("Fortran runtime error: {}", e.0)),
    }
}

fn tail(math: i32, mfh: i32, mth: i32, ns: i32) -> String {
    format!("{}{}{}{}", ffmt::i(math, 4), ffmt::i(mfh, 2), ffmt::i(mth, 3), ffmt::i(ns, 5))
}

fn bstr(b: &[u8]) -> String {
    String::from_utf8_lossy(b).to_string()
}

/// `a11`: 11-character ENDF representation of a real. `hx` is in/out because
/// one rare path leaves it unchanged in the Fortran.
pub fn a11(x: f64, hx: &mut [u8; 11]) {
    let zero = 0.0f64;
    let tenth = 0.1e0f64;
    let onem = 0.999999999e0f64;
    let one = 1.0f64;
    let ten = 10.0f64;
    let tmil = 0.999999999e7f64;
    let up = 0.000000001e0f64;
    let top7 = 9.9999995e0f64;
    let top6 = 9.999995e0f64;
    let top5 = 9.99995e0f64;
    let top9 = 9.999999995e0f64;
    let bot9 = 9.99999995e0f64;
    let set = |hx: &mut [u8; 11], s: String| {
        let b = s.as_bytes();
        // internal write: the record is the whole variable, blank padded
        for j in 0..11 {
            hx[j] = if j < b.len() { b[j] } else { b' ' };
        }
    };
    let f;
    let s;
    let mut n;
    if x == zero {
        f = 0.0;
        s = '+';
        n = 0;
    } else if x > tenth && x < tmil {
        // label 130
        let mut f: f64;
        let s: char;
        let mut n = int(log10(x.abs()));
        if x.abs() < onem {
            // label 140
            f = x;
            s = '-';
            n = -1;
        } else {
            f = x / powi(ten, n);
            s = '+';
            if n.abs() < 10 && f.abs() < top9 {
            } else if n.abs() >= 10 && f.abs() < bot9 {
            } else {
                f = f / ten;
                n = n + 1;
            }
        }
        // label 150
        if n <= 0 {
            set(hx, ffmt::f(f, 11, 8));
        }
        if n == 1 {
            set(hx, ffmt::fp(f, 11, 7, 1));
        }
        if n == 2 {
            set(hx, ffmt::fp(f, 11, 6, 2));
        }
        if n == 3 {
            set(hx, ffmt::fp(f, 11, 5, 3));
        }
        if n == 4 {
            set(hx, ffmt::fp(f, 11, 4, 4));
        }
        if n == 5 {
            set(hx, ffmt::fp(f, 11, 3, 5));
        }
        if n == 6 {
            set(hx, ffmt::fp(f, 11, 2, 6));
        }
        if n == -1 {
            n = 1;
        }
        if f > onem && &hx[9..11] == b"00" {
            set(hx, format!("{}{}{}", ffmt::f(f, 9, 6), s, ffmt::i(n, 1)));
        }
        if f > tenth && f < onem && hx[10] == b'0' {
            set(hx, format!("{}{}{}", ffmt::fp(f, 9, 6, 1), s, ffmt::i(n, 1)));
        }
        return;
    } else {
        let xx = x;
        n = int(log10(xx.abs()));
        let ff;
        if xx.abs() < one {
            // label 110
            n = 1 - n;
            let mut f2 = xx * powi(ten, n);
            let mut s2 = '-';
            if n.abs() < 10 && f2.abs() < top7 {
            } else if n.abs() < 100 && n.abs() >= 10 && f2.abs() < top6 {
            } else if n.abs() >= 100 && f2.abs() < top5 {
            } else {
                f2 = f2 / ten;
                f2 = f2 + up;
                n = n - 1;
                if n <= 0 {
                    s2 = '+';
                }
            }
            ff = f2;
            s = s2;
        } else {
            let mut f2 = xx / powi(ten, n);
            if n.abs() < 10 && f2.abs() < top7 {
            } else if n.abs() < 100 && n.abs() >= 10 && f2.abs() < top6 {
            } else if n.abs() >= 100 && f2.abs() < top5 {
            } else {
                f2 = f2 / ten;
                f2 = f2 + up;
                n = n + 1;
            }
            ff = f2;
            s = '+';
        }
        f = ff;
    }
    // label 160
    if n.abs() < 10 {
        set(hx, format!("{}{}{}", ffmt::f(f, 9, 6), s, ffmt::i(n, 1)));
    } else if n.abs() < 100 {
        set(hx, format!("{}{}{}", ffmt::f(f, 8, 5), s, ffmt::i(n, 2)));
    } else {
        set(hx, format!("{}{}{}", ffmt::f(f, 7, 4), s, ffmt::i(n, 3)));
    }
}

/// a11 into a fresh blank field (the common case).
pub fn a11s(x: f64) -> String {
    let mut h = [b' '; 11];
    a11(x, &mut h);
    bstr(&h)
}

impl Ctx {
    pub fn new(io: Box<Io>) -> Self {
        Ctx {
            io,
            npage: 306,
            iverf: -1,
            c1h: 0.0,
            c2h: 0.0,
            l1h: 0,
            l2h: 0,
            n1h: 0,
            n2h: 0,
            math: 0,
            mfh: 0,
            mth: 0,
            nsh: 0,
            nsp: 0,
            nsc: 0,
            thr6: 0.0,
            nah: 0,
            ah: [0.0; 6],
            ihol: 0,
            rhol: [0u8; 136],
            nr: 0,
            nbt: [0; 24],
            jnt: [0; 24],
            cfield: [[b' '; 11]; 2],
            g1: GetySave::default(),
            g2: GetySave::default(),
            t0: std::time::Instant::now(),
        }
    }

    fn rhol_get(&self, j: usize) -> f64 {
        f64::from_le_bytes(self.rhol[8 * (j - 1)..8 * j].try_into().unwrap())
    }
    fn rhol_set(&mut self, j: usize, v: f64) {
        self.rhol[8 * (j - 1)..8 * j].copy_from_slice(&v.to_le_bytes());
    }

    // ------------------------------------------------------------ helpers
    fn wbin(&mut self, unit: i32, nb: i32, nw: i32, a: &F, n: i32) {
        let mut r = RecW::new();
        r.i(self.math).i(self.mfh).i(self.mth).i(nb).i(nw);
        for i in 1..=n {
            r.r(a[i]);
        }
        self.io.write_rec(unit, r);
    }

    /// Binary read of `math,mfh,mth,nb,nw,(a(i),i=1,cnt(nw))`.
    fn rbin(&mut self, unit: i32, a: &mut F, cnt: impl Fn(i32) -> i32) -> (i32, i32) {
        let (math, mfh, mth, nb, nw, vals) = self.io.read_rec_with(unit, |r| {
            let math = r.i();
            let mfh = r.i();
            let mth = r.i();
            let nb = r.i();
            let nw = r.i();
            let n = cnt(nw);
            let mut v = Vec::with_capacity(n.max(0) as usize);
            for _ in 0..n {
                v.push(r.r());
            }
            (math, mfh, mth, nb, nw, v)
        });
        self.math = math;
        self.mfh = mfh;
        self.mth = mth;
        for (i, x) in vals.iter().enumerate() {
            a[i as i32 + 1] = *x;
        }
        (nb, nw)
    }

    // ------------------------------------------------------------ contio
    pub fn contio(&mut self, nin: i32, nout: i32, nscr: i32, a: &mut F, nb: &mut i32, nw: &mut i32) {
        if nin < 0 {
            let (b, w) = self.rbin(nin.abs(), a, |_| 6);
            *nb = b;
            *nw = w;
        } else if nin > 0 {
            let line = self.io.read_line(nin);
            for i in 0..6 {
                a[i as i32 + 1] = rd_e(&ffmt::field(&line, 11 * i, 11));
            }
            let (m1, m2, m3, ns) = hrec_fields(&line);
            self.math = m1;
            self.mfh = m2;
            self.mth = m3;
            self.nsp = ns;
        }
        *nb = 0;
        *nw = 6;
        self.c1h = a[1];
        self.c2h = a[2];
        self.l1h = nint(a[3]);
        self.l2h = nint(a[4]);
        self.n1h = nint(a[5]);
        self.n2h = nint(a[6]);
        a[3] = self.l1h as f64;
        a[4] = self.l2h as f64;
        a[5] = self.n1h as f64;
        a[6] = self.n2h as f64;
        if nout == 0 && nscr == 0 {
            return;
        }
        let mut inout = nout.abs();
        if nout < 0 {
            self.wbin(inout, *nb, *nw, a, 6);
            inout = 0;
        }
        let mut inscr = nscr.abs();
        if nscr < 0 {
            self.wbin(inscr, *nb, *nw, a, 6);
            inscr = 0;
        }
        if nout <= 0 && nscr <= 0 {
            return;
        }
        if self.math == -1 {
            self.atend(inout, inscr);
            return;
        }
        if self.math == 0 {
            self.amend(inout, inscr);
            return;
        }
        if self.mfh == 0 {
            self.afend(inout, inscr);
            return;
        }
        if self.mth == 0 {
            self.asend(inout, inscr);
            return;
        }
        let mut f1 = self.cfield[0];
        let mut f2 = self.cfield[1];
        a11(self.c1h, &mut f1);
        a11(self.c2h, &mut f2);
        self.cfield = [f1, f2];
        let body = format!(
            "{}{}{}{}{}{}",
            bstr(&f1),
            bstr(&f2),
            ffmt::i(self.l1h, 11),
            ffmt::i(self.l2h, 11),
            ffmt::i(self.n1h, 11),
            ffmt::i(self.n2h, 11)
        );
        if nscr > 0 {
            let s = format!("{}{}", body, tail(self.math, self.mfh, self.mth, self.nsc));
            self.io.write_line(nscr, &s);
            self.nsc += 1;
            if self.nsc > 99999 {
                self.nsc = 1;
            }
        }
        if nout > 0 {
            let s = format!("{}{}", body, tail(self.math, self.mfh, self.mth, self.nsh));
            self.io.write_line(nout, &s);
            self.nsh += 1;
            if self.nsh > 99999 {
                self.nsh = 1;
            }
        }
    }

    // ------------------------------------------------------------ listio
    pub fn listio(&mut self, nin: i32, nout: i32, nscr: i32, a: &mut F, nb: &mut i32, nw: &mut i32) {
        self.ihol = 0;
        if nin < 0 {
            let (b, w) = self.rbin(nin.abs(), a, |nw| nw);
            *nb = b;
            *nw = w;
            self.c1h = a[1];
            self.c2h = a[2];
            self.l1h = nint(a[3]);
            self.l2h = nint(a[4]);
            self.n1h = nint(a[5]);
            self.n2h = nint(a[6]);
            a[3] = self.l1h as f64;
            a[4] = self.l2h as f64;
            a[5] = self.n1h as f64;
            a[6] = self.n2h as f64;
        } else if nin > 0 {
            self.contio(nin, 0, 0, a, nb, nw);
            *nb = 6 + self.n1h;
            let mut iend = self.npage;
            if self.n1h < self.npage {
                iend = self.n1h;
            }
            *nw = iend + *nw;
            let mut k = 6;
            let mut i = 1;
            while i <= iend {
                self.lineio(nin, 0, 0);
                let mut l = self.n1h - k + 6;
                if l != 0 {
                    if l > 6 {
                        l = 6;
                    }
                    for j in 1..=l {
                        k += 1;
                        a[k] = self.ah[(j - 1) as usize];
                    }
                }
                i += 6;
            }
            a[1] = self.c1h;
            a[2] = self.c2h;
            a[3] = self.l1h as f64;
            a[4] = self.l2h as f64;
            a[5] = self.n1h as f64;
            a[6] = self.n2h as f64;
        }
        self.n1h = nint(a[5]);
        if nin == 0 && self.n1h <= self.npage {
            *nw = 6 + self.n1h;
        }
        if nin == 0 && self.n1h > self.npage {
            *nw = 6 + self.npage;
        }
        *nb = 6 + self.n1h - *nw;
        if nout == 0 && nscr == 0 {
            return;
        }
        self.c1h = a[1];
        self.c2h = a[2];
        self.l1h = nint(a[3]);
        self.l2h = nint(a[4]);
        self.n2h = nint(a[6]);
        let mut inout = nout.abs();
        if nout < 0 {
            self.wbin(inout, *nb, *nw, a, *nw);
            inout = 0;
        }
        let mut inscr = nscr.abs();
        if nscr < 0 {
            self.wbin(inscr, *nb, *nw, a, *nw);
            inscr = 0;
        }
        let (mut nbx, mut nwx) = (0, 0);
        self.contio(0, inout, inscr, a, &mut nbx, &mut nwx);
        let mut k = 6;
        let iend = *nw - 6;
        let mut i = 1;
        while i <= iend {
            let mut l = iend - k + 6;
            if l != 0 {
                if l > 6 {
                    l = 6;
                }
                for j in 1..=l {
                    k += 1;
                    self.ah[(j - 1) as usize] = a[k];
                }
                self.nah = l;
                self.lineio(0, inout, inscr);
            }
            i += 6;
        }
    }

    // ------------------------------------------------------------ tab1io
    pub fn tab1io(&mut self, nin: i32, nout: i32, nscr: i32, a: &mut F, nb: &mut i32, nw: &mut i32) {
        self.ihol = 0;
        if nin < 0 {
            let (b, w) = self.rbin(nin.abs(), a, |nw| nw);
            *nb = b;
            *nw = w;
            self.c1h = a[1];
            self.c2h = a[2];
            self.l1h = nint(a[3]);
            self.l2h = nint(a[4]);
            self.n1h = nint(a[5]);
            self.n2h = nint(a[6]);
            if self.n1h <= 0 {
                let s = format!(" illegal TAB1, nr<=0 for mf/mt = {}/{}", ffmt::i(self.mfh, 2), ffmt::i(self.mth, 3));
                self.error("endf", &s, "");
            }
            a[3] = self.l1h as f64;
            a[4] = self.l2h as f64;
            a[5] = self.n1h as f64;
            a[6] = self.n2h as f64;
            for i in 1..=self.n1h {
                self.nbt[i as usize] = nint(a[2 * i + 5]);
                self.jnt[i as usize] = nint(a[2 * i + 6]);
            }
        } else if nin > 0 {
            self.contio(nin, 0, 0, a, nb, nw);
            if self.n1h <= 0 {
                let s = format!(" illegal TAB1, nr<=0 for mf/mt = {}/{}", ffmt::i(self.mfh, 2), ffmt::i(self.mth, 3));
                self.error("endf", &s, "");
            }
            *nb = 6 + 2 * self.n1h + 2 * self.n2h;
            self.tablio(nin, 0, 0);
            let mut k = 6;
            for i in 1..=self.nr {
                k += 1;
                a[k] = self.nbt[i as usize] as f64;
                k += 1;
                a[k] = self.jnt[i as usize] as f64;
            }
            let mut lim = self.npage;
            if 2 * self.n2h < self.npage {
                lim = 2 * self.n2h;
            }
            let mut j = 0;
            let mut i = 1;
            while i <= lim {
                self.lineio(nin, 0, 0);
                let mut l = 6;
                if j + 6 > lim {
                    l = lim - j;
                }
                for m in 1..=l {
                    j += 1;
                    a[k + j] = self.ah[(m - 1) as usize];
                }
                i += 6;
            }
            a[1] = self.c1h;
            a[2] = self.c2h;
            a[3] = self.l1h as f64;
            a[4] = self.l2h as f64;
            a[5] = self.n1h as f64;
            a[6] = self.n2h as f64;
            *nw = k + j;
        }
        if nin != 0 && self.nbt[self.n1h as usize] != self.n2h {
            let s = format!(" illegal TAB1, nbt(nr)/=np for mf/mt = {}/{}", ffmt::i(self.mfh, 2), ffmt::i(self.mth, 3));
            self.error("endf", &s, "");
        }
        self.n1h = nint(a[5]);
        self.n2h = nint(a[6]);
        if nin == 0 && 2 * self.n2h <= self.npage {
            *nw = 6 + 2 * self.n1h + 2 * self.n2h;
        }
        if nin == 0 && 2 * self.n2h > self.npage {
            *nw = 6 + 2 * self.n1h + self.npage;
        }
        *nb = 6 + 2 * self.n1h + 2 * self.n2h - *nw;
        if nout == 0 && nscr == 0 {
            return;
        }
        self.c1h = a[1];
        self.c2h = a[2];
        self.l1h = nint(a[3]);
        self.l2h = nint(a[4]);
        for i in 1..=self.n1h {
            self.nbt[i as usize] = nint(a[2 * i + 5]);
            self.jnt[i as usize] = nint(a[2 * i + 6]);
        }
        let mut inout = nout.abs();
        if nout < 0 {
            self.wbin(inout, *nb, *nw, a, *nw);
            inout = 0;
        }
        let mut inscr = nscr.abs();
        if nscr < 0 {
            self.wbin(inscr, *nb, *nw, a, *nw);
            inscr = 0;
        }
        if nout <= 0 && nscr <= 0 {
            return;
        }
        let (mut nbx, mut nwx) = (0, 0);
        self.contio(0, inout, inscr, a, &mut nbx, &mut nwx);
        for i in 1..=self.n1h {
            self.nbt[i as usize] = nint(a[2 * i + 5]);
            self.jnt[i as usize] = nint(a[2 * i + 6]);
        }
        self.nr = self.n1h;
        self.tablio(0, inout, inscr);
        let mut k = 6 + 2 * self.nr;
        let kp1 = k + 1;
        let mut i = kp1;
        while i <= *nw {
            let mut l = 6;
            if k + 6 > *nw {
                l = *nw - k;
            }
            for j in 1..=l {
                k += 1;
                self.ah[(j - 1) as usize] = a[k];
            }
            self.nah = l;
            self.lineio(0, inout, inscr);
            i += 6;
        }
    }

    // ------------------------------------------------------------ tab2io
    pub fn tab2io(&mut self, nin: i32, nout: i32, nscr: i32, a: &mut F, nb: &mut i32, nw: &mut i32) {
        if nin < 0 {
            let (b, w) = self.rbin(nin.abs(), a, |nw| nw);
            *nb = b;
            *nw = w;
            self.c1h = a[1];
            self.c2h = a[2];
            self.l1h = nint(a[3]);
            self.l2h = nint(a[4]);
            self.n1h = nint(a[5]);
            self.n2h = nint(a[6]);
            a[3] = self.l1h as f64;
            a[4] = self.l2h as f64;
            a[5] = self.n1h as f64;
            a[6] = self.n2h as f64;
        } else if nin > 0 {
            self.contio(nin, 0, 0, a, nb, nw);
            self.tablio(nin, 0, 0);
            a[1] = self.c1h;
            a[2] = self.c2h;
            a[3] = self.l1h as f64;
            a[4] = self.l2h as f64;
            a[5] = self.n1h as f64;
            a[6] = self.n2h as f64;
            let mut k = 6;
            for i in 1..=self.n1h {
                k += 1;
                a[k] = self.nbt[i as usize] as f64;
                k += 1;
                a[k] = self.jnt[i as usize] as f64;
            }
            *nw = *nw + 2 * self.n1h;
        }
        *nb = 0;
        if nout == 0 && nscr == 0 {
            return;
        }
        self.c1h = a[1];
        self.c2h = a[2];
        self.l1h = nint(a[3]);
        self.l2h = nint(a[4]);
        self.n1h = nint(a[5]);
        self.n2h = nint(a[6]);
        for i in 1..=self.n1h {
            self.nbt[i as usize] = nint(a[2 * i + 5]);
            self.jnt[i as usize] = nint(a[2 * i + 6]);
        }
        let mut inout = nout.abs();
        if nout < 0 {
            self.wbin(inout, *nb, *nw, a, *nw);
            inout = 0;
        }
        let mut inscr = nscr.abs();
        if nscr < 0 {
            self.wbin(inscr, *nb, *nw, a, *nw);
            inscr = 0;
        }
        if nout <= 0 && nscr <= 0 {
            return;
        }
        let (mut nbx, mut nwx) = (0, 0);
        self.contio(0, inout, inscr, a, &mut nbx, &mut nwx);
        for i in 1..=self.n1h {
            self.nbt[i as usize] = nint(a[2 * i + 5]);
            self.jnt[i as usize] = nint(a[2 * i + 6]);
        }
        self.nr = self.n1h;
        self.tablio(0, inout, inscr);
    }

    // ------------------------------------------------------------ moreio
    pub fn moreio(&mut self, nin: i32, nout: i32, nscr: i32, a: &mut F, nb: &mut i32, nw: &mut i32) {
        if nin < 0 {
            let (b, w) = self.rbin(nin.abs(), a, |nw| nw);
            *nb = b;
            *nw = w;
        } else if nin > 0 {
            *nw = *nb;
            if *nw > self.npage {
                *nw = self.npage;
            }
            *nb = *nb - *nw;
            let mut k = 0;
            if self.ihol != 1 {
                let mut i = 1;
                while i <= *nw {
                    self.lineio(nin, 0, 0);
                    let mut l = 6;
                    if k + 6 > *nw {
                        l = *nw - k;
                    }
                    for j in 1..=l {
                        k += 1;
                        a[k] = self.ah[(j - 1) as usize];
                    }
                    i += 6;
                }
            } else {
                let mut i = 1;
                while i <= *nw {
                    self.hollio(nin, 0, 0);
                    for j in 1..=17 {
                        k += 1;
                        a[k] = self.rhol_get(j);
                    }
                    i += 17;
                }
            }
        }
        let mut inout = nout.abs();
        let mut inscr = nscr.abs();
        if nin == 0 {
            *nw = *nb;
        }
        if *nw > self.npage {
            *nw = self.npage;
        }
        if nin == 0 {
            *nb = *nb - *nw;
        }
        if nout < 0 {
            self.wbin(inout, *nb, *nw, a, *nw);
            inout = 0;
        }
        if nscr < 0 {
            self.wbin(inscr, *nb, *nw, a, *nw);
            inscr = 0;
        }
        if nout <= 0 && nscr <= 0 {
            return;
        }
        if self.ihol != 1 {
            let mut k = 0;
            let mut i = 1;
            while i <= *nw {
                let mut l = 6;
                if l + k > *nw {
                    l = *nw - k;
                }
                for j in 1..=l {
                    k += 1;
                    self.ah[(j - 1) as usize] = a[k];
                }
                self.nah = l;
                self.lineio(0, inout, inscr);
                i += 6;
            }
        } else {
            let mut k = 0;
            let mut i = 1;
            while i <= *nw {
                for j in 1..=17 {
                    k += 1;
                    self.rhol_set(j, a[k]);
                }
                self.hollio(0, inout, inscr);
                i += 17;
            }
        }
    }

    // ------------------------------------------------------------ tpidio
    pub fn tpidio(&mut self, nin: i32, nout: i32, nscr: i32, a: &mut F, nb: &mut i32, nw: &mut i32) {
        // rb(17) real equivalenced with hb(17) character(4): local, so we
        // start it from the a() contents each time it is used.
        if nin < 0 {
            let (b, w) = self.rbin(nin.abs(), a, |_| 17);
            *nb = b;
            *nw = w;
        } else if nin > 0 {
            let line = self.io.read_line(nin);
            let mut text = ffmt::field(&line, 0, 66);
            text.push(b' ');
            text.push(b' ');
            let (m1, m2, m3, _) = hrec_fields(&line);
            self.math = m1;
            self.mfh = m2;
            self.mth = m3;
            *nw = 17;
            // rb is a fresh local; its tail bytes are undefined in Fortran.
            let mut rb = [0.0f64; 17];
            crate::farr::pack_text(&text, &mut rb);
            for i in 1..=17 {
                a[i] = rb[(i - 1) as usize];
            }
        }
        *nb = 0;
        let inout = nout.abs();
        let inscr = nscr.abs();
        if nout < 0 {
            self.wbin(inout, *nb, *nw, a, 17);
        }
        if nscr < 0 {
            self.wbin(inscr, *nb, *nw, a, 17);
        }
        if nout <= 0 && nscr <= 0 {
            return;
        }
        let mut rb = [0.0f64; 17];
        for i in 1..=17 {
            rb[(i - 1) as usize] = a[i];
        }
        let t = crate::farr::unpack_text(&rb);
        let text = bstr(&t[..66]);
        if nout > 0 {
            let s = format!("{}{}", text, tail(self.math, self.mfh, self.mth, self.nsh));
            self.io.write_line(nout, &s);
            self.nsh += 1;
        }
        if nscr > 0 {
            let s = format!("{}{}", text, tail(self.math, self.mfh, self.mth, self.nsc));
            self.io.write_line(nscr, &s);
            self.nsc += 1;
        }
    }

    // ------------------------------------------------------------ hdatio
    pub fn hdatio(&mut self, nin: i32, nout: i32, nscr: i32, a: &mut F, nb: &mut i32, nw: &mut i32) {
        self.ihol = 1;
        if nin < 0 {
            let (b, w) = self.rbin(nin.abs(), a, |nw| nw);
            *nb = b;
            *nw = w;
            self.c1h = a[1];
            self.c2h = a[2];
            self.l1h = nint(a[3]);
            self.l2h = nint(a[4]);
            self.n1h = nint(a[5]) / 17;
            self.n2h = nint(a[6]);
            a[3] = self.l1h as f64;
            a[4] = self.l2h as f64;
            a[5] = (self.n1h * 17) as f64;
            a[6] = self.n2h as f64;
        } else if nin > 0 {
            self.contio(nin, 0, 0, a, nb, nw);
            a[5] = a[5] * 17.0;
            *nb = 6 + nint(a[5]);
            let mut iend = nint(a[5]);
            if iend > self.npage {
                iend = self.npage;
            }
            let mut k = 6;
            let mut i = 1;
            while i <= iend {
                self.hollio(nin, 0, 0);
                for j in 1..=17 {
                    k += 1;
                    a[k] = self.rhol_get(j);
                }
                i += 17;
            }
            *nw = k;
        }
        self.n1h = nint(a[5]) / 17;
        if nin == 0 && a[5] <= self.npage as f64 {
            *nw = 6 + nint(a[5]);
        }
        if nin == 0 && a[5] > self.npage as f64 {
            *nw = 6 + self.npage;
        }
        *nb = 6 + nint(a[5]) - *nw;
        let mut inout = nout.abs();
        let mut inscr = nscr.abs();
        a[5] = (self.n1h * 17) as f64;
        if nout < 0 {
            self.wbin(inout, *nb, *nw, a, *nw);
            inout = 0;
        }
        if nscr < 0 {
            self.wbin(inscr, *nb, *nw, a, *nw);
            inscr = 0;
        }
        if nout <= 0 && nscr <= 0 {
            return;
        }
        a[5] = self.n1h as f64;
        let ncds = (*nw - 6) / 17;
        let (mut nbx, mut nwx) = (0, 0);
        self.contio(0, inout, inscr, a, &mut nbx, &mut nwx);
        let mut j = 6;
        for _ in 1..=ncds {
            for k in 1..=17 {
                j += 1;
                self.rhol_set(k, a[j]);
            }
            self.hollio(0, inout, inscr);
        }
    }

    // ------------------------------------------------------------ dictio
    pub fn dictio(&mut self, nin: i32, nout: i32, nscr: i32, a: &mut F, nb: &mut i32, nw: &mut i32) {
        let nx = *nw;
        if nin < 0 {
            let inin = nin.abs();
            let mut k = 0;
            for _ in 1..=nx {
                let mut sub = a.at(k + 1);
                let (b, w) = self.rbin(inin, &mut sub, |_| 6);
                *nb = b;
                *nw = w;
                k += 6;
            }
            *nw = k;
        } else if nin > 0 {
            let mut k = 0;
            for _ in 1..=nx {
                let line = self.io.read_line(nin);
                for i in 0..6 {
                    a[k + i + 1] = rd_e(&ffmt::field(&line, 11 * i as usize, 11));
                }
                let (m1, m2, m3, ns) = hrec_fields(&line);
                self.math = m1;
                self.mfh = m2;
                self.mth = m3;
                self.nsp = ns;
                k += 6;
            }
            *nw = k;
        }
        *nb = 0;
        if nout == 0 && nscr == 0 {
            return;
        }
        let inout = nout.abs();
        let inscr = nscr.abs();
        let n = 6;
        let mut k = 0;
        for _ in 1..=nx {
            self.l1h = nint(a[k + 3]);
            self.l2h = nint(a[k + 4]);
            self.n1h = nint(a[k + 5]);
            self.n2h = nint(a[k + 6]);
            a[k + 3] = self.l1h as f64;
            a[k + 4] = self.l2h as f64;
            a[k + 5] = self.n1h as f64;
            a[k + 6] = self.n2h as f64;
            let body = format!(
                "{}{}{}{}{}",
                " ".repeat(22),
                ffmt::i(self.l1h, 11),
                ffmt::i(self.l2h, 11),
                ffmt::i(self.n1h, 11),
                ffmt::i(self.n2h, 11)
            );
            if nout < 0 {
                let sub = a.at(k + 1);
                self.wbin(inout, *nb, n, &sub, n);
            }
            if nout > 0 {
                let s = format!("{}{}", body, tail(self.math, self.mfh, self.mth, self.nsh));
                self.io.write_line(nout, &s);
                self.nsh += 1;
                if self.nsh > 99999 {
                    self.nsh = 1;
                }
            }
            if nscr < 0 {
                let sub = a.at(k + 1);
                self.wbin(inscr, *nb, n, &sub, n);
            }
            if nscr > 0 {
                let s = format!("{}{}", body, tail(self.math, self.mfh, self.mth, self.nsc));
                self.io.write_line(nscr, &s);
                self.nsc += 1;
                if self.nsc > 99999 {
                    self.nsc = 1;
                }
            }
            k += 6;
        }
        *nw = k;
    }

    // ------------------------------------------------------------ lineio
    fn lineio(&mut self, nin: i32, nout: i32, nscr: i32) {
        if nin > 0 {
            let line = self.io.read_line(nin);
            for i in 0..6 {
                self.ah[i] = rd_e(&ffmt::field(&line, 11 * i, 11));
            }
            let (m1, m2, m3, ns) = hrec_fields(&line);
            self.math = m1;
            self.mfh = m2;
            self.mth = m3;
            self.nsp = ns;
            self.nah = 6;
        }
        if nout != 0 || nscr != 0 {
            let mut body = String::with_capacity(66);
            for k in 1..=6 {
                if k <= self.nah {
                    body.push_str(&a11s(self.ah[(k - 1) as usize]));
                } else {
                    body.push_str("           ");
                }
            }
            if nout != 0 {
                let s = format!("{}{}", body, tail(self.math, self.mfh, self.mth, self.nsh));
                self.io.write_line(nout, &s);
                self.nsh += 1;
                if self.nsh > 99999 {
                    self.nsh = 1;
                }
            }
            if nscr != 0 {
                let s = format!("{}{}", body, tail(self.math, self.mfh, self.mth, self.nsc));
                self.io.write_line(nscr, &s);
                self.nsc += 1;
                if self.nsc > 99999 {
                    self.nsc = 1;
                }
            }
        }
    }

    // ------------------------------------------------------------ tablio
    fn tablio(&mut self, nin: i32, nout: i32, nscr: i32) {
        if nin != 0 {
            self.nr = self.n1h;
        }
        let mut i = 1;
        while i <= self.nr {
            let ip2 = i + 2;
            let mut ipend = ip2;
            if ip2 > self.nr {
                ipend = self.nr;
            }
            let nterm = ipend - i + 1;
            if nin != 0 {
                let line = self.io.read_line(nin);
                let mut c = 0usize;
                for j in i..=ip2 {
                    let v1 = rd_i(&ffmt::field(&line, 11 * c, 11));
                    let v2 = rd_i(&ffmt::field(&line, 11 * (c + 1), 11));
                    if (j as usize) < self.nbt.len() {
                        self.nbt[j as usize] = v1;
                        self.jnt[j as usize] = v2;
                    }
                    c += 2;
                }
                let (m1, m2, m3, ns) = hrec_fields(&line);
                self.math = m1;
                self.mfh = m2;
                self.mth = m3;
                self.nsp = ns;
            }
            if nout != 0 || nscr != 0 {
                let mut body = String::new();
                for j in i..=ipend {
                    body.push_str(&ffmt::i(self.nbt[j as usize], 11));
                    body.push_str(&ffmt::i(self.jnt[j as usize], 11));
                }
                body.push_str(&" ".repeat((22 * (3 - nterm)) as usize));
                if nout != 0 {
                    let s = format!("{}{}", body, tail(self.math, self.mfh, self.mth, self.nsh));
                    self.io.write_line(nout, &s);
                    self.nsh += 1;
                    if self.nsh > 99999 {
                        self.nsh = 1;
                    }
                }
                if nscr != 0 {
                    let s = format!("{}{}", body, tail(self.math, self.mfh, self.mth, self.nsc));
                    self.io.write_line(nscr, &s);
                    self.nsc += 1;
                    if self.nsc > 99999 {
                        self.nsc = 1;
                    }
                }
            }
            i += 3;
        }
    }

    // ------------------------------------------------------------ hollio
    fn hollio(&mut self, nin: i32, nout: i32, nscr: i32) {
        if nin != 0 {
            let line = self.io.read_line(nin);
            let t = ffmt::field(&line, 0, 66);
            self.rhol[..66].copy_from_slice(&t);
            self.rhol[66] = b' ';
            self.rhol[67] = b' ';
            let (m1, m2, m3, ns) = hrec_fields(&line);
            self.math = m1;
            self.mfh = m2;
            self.mth = m3;
            self.nsp = ns;
        }
        let text = bstr(&self.rhol[..66]);
        if nout != 0 {
            let s = format!("{}{}", text, tail(self.math, self.mfh, self.mth, self.nsh));
            self.io.write_line(nout, &s);
            if self.nsh > 99999 {
                self.nsh = 1;
            }
            self.nsh += 1;
        }
        if nscr != 0 {
            let s = format!("{}{}", text, tail(self.math, self.mfh, self.mth, self.nsc));
            self.io.write_line(nscr, &s);
            if self.nsc > 99999 {
                self.nsc = 1;
            }
            self.nsc += 1;
        }
    }

    // ------------------------------------------------- tosend/tofend/...
    fn copy_to(&mut self, nin: i32, nout: i32, nscr: i32, a: &mut F, which: u8, name: &str) {
        if nin > 0 && (nout < 0 || nscr < 0) {
            self.error(name, "mode conversion not allowed.", " ");
        }
        if nin < 0 && (nout > 0 || nscr > 0) {
            self.error(name, "mode conversion not allowed.", " ");
        }
        loop {
            let (mat, mf, mt);
            if nin > 0 {
                let line = self.io.read_line(nin);
                let text = bstr(&ffmt::field(&line, 0, 66));
                let (m1, m2, m3, _ns) = hrec_fields(&line);
                mat = m1;
                mf = m2;
                mt = m3;
                if nscr != 0 {
                    let s = format!("{}{}", text, tail(mat, mf, mt, self.nsc));
                    self.io.write_line(nscr, &s);
                    self.nsc += 1;
                    if self.nsc > 99999 {
                        self.nsc = 1;
                    }
                }
                if nout != 0 {
                    let s = format!("{}{}", text, tail(mat, mf, mt, self.nsh));
                    self.io.write_line(nout, &s);
                    self.nsh += 1;
                    if self.nsh > 99999 {
                        self.nsh = 1;
                    }
                }
            } else {
                let inin = nin.abs();
                let rec = self.io.read_rec(inin);
                let mut rr = crate::fio::RecR { b: &rec, off: 0 };
                mat = rr.i();
                mf = rr.i();
                mt = rr.i();
                let nb = rr.i();
                let nw = rr.i();
                let mut v = Vec::with_capacity(nw.max(0) as usize);
                for _ in 0..nw {
                    v.push(rr.r());
                }
                for (i, x) in v.iter().enumerate() {
                    a[i as i32 + 1] = *x;
                }
                let mk = || {
                    let mut r = RecW::new();
                    r.i(mat).i(mf).i(mt).i(nb).i(nw);
                    for x in &v {
                        r.r(*x);
                    }
                    r
                };
                if nscr != 0 {
                    self.io.write_rec(nscr.abs(), mk());
                }
                if nout != 0 {
                    self.io.write_rec(nout.abs(), mk());
                }
            }
            let cont = match which {
                b's' => mt > 0,
                b'f' => mf > 0,
                b'm' => mat > 0,
                _ => mat > -1,
            };
            if !cont {
                break;
            }
        }
    }

    pub fn tosend(&mut self, nin: i32, nout: i32, nscr: i32, a: &mut F) {
        self.copy_to(nin, nout, nscr, a, b's', "tosend");
    }
    pub fn tofend(&mut self, nin: i32, nout: i32, nscr: i32, a: &mut F) {
        self.copy_to(nin, nout, nscr, a, b'f', "tofend");
    }
    pub fn tomend(&mut self, nin: i32, nout: i32, nscr: i32, a: &mut F) {
        self.copy_to(nin, nout, nscr, a, b'm', "tomend");
    }
    pub fn totend(&mut self, nin: i32, nout: i32, nscr: i32, a: &mut F) {
        self.copy_to(nin, nout, nscr, a, b't', "totend");
    }

    fn zrec(&mut self, unit: i32) {
        let mut r = RecW::new();
        r.i(self.math).i(self.mfh).i(self.mth).i(0).i(6);
        for _ in 0..6 {
            r.r(0.0);
        }
        self.io.write_rec(unit, r);
    }

    pub fn asend(&mut self, nout: i32, nscr: i32) {
        self.mth = 0;
        if nout > 0 {
            self.nsh = 99999;
            let s = format!("{}{}", " ".repeat(66), tail(self.math, self.mfh, self.mth, self.nsh));
            self.io.write_line(nout, &s);
            self.nsh = 1;
        }
        if nout < 0 {
            self.zrec(nout.abs());
        }
        if nscr > 0 {
            self.nsc = 99999;
            let s = format!("{}{}", " ".repeat(66), tail(self.math, self.mfh, self.mth, self.nsc));
            self.io.write_line(nscr, &s);
            self.nsc = 1;
        }
        if nscr < 0 {
            self.zrec(nscr.abs());
        }
    }

    fn end_common(&mut self, nout: i32, nscr: i32) {
        if nout > 0 {
            self.nsh = 0;
            let s = format!("{}{}", " ".repeat(66), tail(self.math, self.mfh, self.mth, self.nsh));
            self.io.write_line(nout, &s);
            self.nsh += 1;
        }
        if nout < 0 {
            self.zrec(nout.abs());
        }
        if nscr > 0 {
            self.nsc = 0;
            let s = format!("{}{}", " ".repeat(66), tail(self.math, self.mfh, self.mth, self.nsc));
            self.io.write_line(nscr, &s);
            self.nsc += 1;
        }
        if nscr < 0 {
            self.zrec(nscr.abs());
        }
    }

    pub fn afend(&mut self, nout: i32, nscr: i32) {
        self.mfh = 0;
        self.mth = 0;
        self.end_common(nout, nscr);
    }
    pub fn amend(&mut self, nout: i32, nscr: i32) {
        self.math = 0;
        self.mfh = 0;
        self.mth = 0;
        self.end_common(nout, nscr);
    }
    pub fn atend(&mut self, nout: i32, nscr: i32) {
        self.math = -1;
        self.mfh = 0;
        self.mth = 0;
        self.end_common(nout, nscr);
    }

    // ------------------------------------------------------------ findf
    fn findf_read(&mut self, ntape: i32) -> (i32, i32, i32) {
        let itape = ntape.abs();
        if ntape > 0 {
            let line = self.io.read_line(itape);
            let (m1, m2, m3, _) = hrec_fields(&line);
            (m1, m2, m3)
        } else {
            self.io.read_rec_with(itape, |r| (r.i(), r.i(), r.i()))
        }
    }

    pub fn findf(&mut self, mat: i32, mf: i32, mt: i32, ntape: i32) {
        let itape = ntape.abs();
        let mut ff: (i32, i32, i32) = (0, 0, 0);
        let mut matl = 0;
        let mut mfl = 0;
        let mut igo: i32 = 0;
        // states correspond to the Fortran statement labels
        let mut lbl = 100;
        loop {
            match lbl {
                100 => {
                    ff = self.findf_read(ntape);
                    let (math, _, _) = ff;
                    if math == 0 {
                        lbl = 100;
                        continue;
                    }
                    if math == mat {
                        lbl = 120;
                        continue;
                    }
                    igo = 1;
                    if math == -1 {
                        igo = -1;
                    }
                    if math > mat {
                        igo = -1;
                    }
                                        lbl = 200;
                }
                120 => {
                    let (_math, mfh, _mth) = ff;
                    if mf == 0 {
                        lbl = 300;
                        continue;
                    }
                    if mfh == mf {
                        lbl = 130;
                        continue;
                    }
                    if mfh == 0 {
                        lbl = 100;
                        continue;
                    }
                    igo = 1;
                    if mfh > mf {
                        igo = -1;
                    }
                                        lbl = 200;
                }
                130 => {
                    let (_math, _mfh, mth) = ff;
                    if mt == 0 {
                        lbl = 300;
                        continue;
                    }
                    if mth == mt {
                        lbl = 300;
                        continue;
                    }
                    if mth == 0 {
                        self.skiprz(ntape, -2);
                        lbl = 100;
                        continue;
                    }
                    igo = 1;
                    if mth > mt {
                        igo = -1;
                    }
                                        lbl = 200;
                }
                200 => {
                    let (math0, mfh0, _) = ff;
                    if igo < 0 {
                        self.skiprz(ntape, -2);
                        matl = math0;
                        mfl = mfh0;
                    }
                    let r = self.findf_read(ntape);
                    ff = r;
                    let (math, _mfh, _mth) = r;
                    if math == mat {
                        lbl = 230;
                        continue;
                    }
                    if igo < 0 {
                        lbl = 220;
                        continue;
                    }
                    if math > mat || math == -1 {
                        lbl = 400;
                        continue;
                    }
                    lbl = 200;
                }
                220 => {
                    let (math, mfh, _) = ff;
                    if math < mat && math != 0 {
                        lbl = 400;
                        continue;
                    }
                    if math == matl && mfl == 0 && mfh == 0 {
                        igo = 1;
                    }
                    lbl = 200;
                }
                230 => {
                    let (_, mfh, _) = ff;
                    if mf == 0 {
                        lbl = 300;
                        continue;
                    }
                    if mfh == mf {
                        lbl = 240;
                        continue;
                    }
                    if mfh == 0 {
                        lbl = 200;
                        continue;
                    }
                    if igo < 0 && mfh < mf {
                        lbl = 400;
                        continue;
                    }
                    if igo > 0 && mfh > mf {
                        lbl = 400;
                        continue;
                    }
                    lbl = 200;
                }
                240 => {
                    let (_, _, mth) = ff;
                    if mt == 0 {
                        lbl = 300;
                        continue;
                    }
                    if mth == mt {
                        lbl = 300;
                        continue;
                    }
                    lbl = 200;
                }
                300 => {
                    self.skiprz(ntape, -2);
                    let r = self.findf_read(ntape);
                    ff = r;
                    let (math, mfh, mth) = r;
                    if mth > 0 {
                        lbl = 300;
                        continue;
                    }
                    if mt == 0 && mfh > 0 {
                        lbl = 300;
                        continue;
                    }
                    if mf == 0 && math == mat {
                        lbl = 300;
                        continue;
                    }
                    return;
                }
                _ => {
                    // 400
                    let s = format!(
                        " mat{} mf{} mt{} not on tape {}",
                        ffmt::i(mat, 4),
                        ffmt::i(mf, 2),
                        ffmt::i(mt, 3),
                        ffmt::i(itape, 3)
                    );
                    self.error("findf", &s, " ");
                }
            }
        }
    }

    // ------------------------------------------------------------ terp1
    pub fn terp1(&self, x1: f64, y1: f64, x2: f64, y2: f64, x: f64, y: &mut f64, i: i32) {
        terp1(x1, y1, x2, y2, x, y, i, self.thr6)
    }

    // ------------------------------------------------------------ terpa
    pub fn terpa(&self, y: &mut f64, x: f64, xnext: &mut f64, idis: &mut i32, a: &F, ip: &mut i32, ir: &mut i32) {
        let shade = 1.00001e0f64;
        let xbig = 1.0e12f64;
        let zero = 0.0f64;
        let nr = nint(a[5]);
        let np = nint(a[6]);
        if *ir > nr {
            *ir = nr;
        }
        if *ip > np {
            *ip = np;
        }
        let mut jr = 5 + 2 * *ir;
        let mut jp = 5 + 2 * nr + 2 * *ip;
        *idis = 0;
        let mut lbl = 110;
        loop {
            match lbl {
                110 => {
                    if x < a[jp] {
                        lbl = 120;
                        continue;
                    }
                    if *ip == np {
                        lbl = 150;
                        continue;
                    }
                    jp += 2;
                    *ip += 1;
                    let it = nint(a[jr]);
                    if *ip <= it {
                        lbl = 110;
                        continue;
                    }
                    jr += 2;
                    *ir += 1;
                    lbl = 110;
                }
                120 => {
                    if x > a[jp - 2] {
                        lbl = 130;
                        continue;
                    }
                    if x == a[jp - 2] {
                        lbl = 140;
                        continue;
                    }
                    if *ip == 2 {
                        lbl = 170;
                        continue;
                    }
                    jp -= 2;
                    *ip -= 1;
                    if *ir == 1 {
                        lbl = 110;
                        continue;
                    }
                    let it = nint(a[jr - 2]);
                    if *ip > it {
                        lbl = 110;
                        continue;
                    }
                    jr -= 2;
                    *ir -= 1;
                    lbl = 110;
                }
                130 => {
                    let int = nint(a[jr + 1]);
                    self.terp1(a[jp - 2], a[jp - 1], a[jp], a[jp + 1], x, y, int);
                    *xnext = a[jp];
                    if int == 1 {
                        *idis = 1;
                    }
                    if *ip == np {
                        return;
                    }
                    if a[jp + 2] == *xnext {
                        *idis = 1;
                    }
                    return;
                }
                140 => {
                    *y = a[jp - 1];
                    let int = nint(a[jr + 1]);
                    *xnext = a[jp];
                    if int == 1 {
                        *idis = 1;
                    }
                    if *ip == np {
                        return;
                    }
                    if a[jp + 2] == *xnext {
                        *idis = 1;
                    }
                    return;
                }
                150 => {
                    if x < shade * a[jp] {
                        // 160
                        *y = a[jp + 1];
                        *xnext = xbig;
                        if *y > zero {
                            *xnext = shade * shade * a[jp];
                        }
                        return;
                    }
                    *y = 0.0;
                    *xnext = xbig;
                    return;
                }
                _ => {
                    // 170
                    *y = 0.0;
                    *xnext = a[jp - 2];
                    *idis = 1;
                    return;
                }
            }
        }
    }

    // ------------------------------------------------------------ intega
    pub fn intega(&mut self, f: &mut f64, x1: f64, x2: f64, a: &F, ip: &mut i32, ir: &mut i32) {
        *f = 0.0;
        let nr = nint(a[5]);
        let np = nint(a[6]);
        if *ip > np {
            *ip = np;
        }
        if *ir > nr {
            *ir = nr;
        }
        let mut jr = 5 + 2 * *ir;
        let mut jp = 5 + 2 * nr + 2 * *ip;
        let mut xlo = 0.0;
        let mut lbl = 110;
        loop {
            match lbl {
                110 => {
                    if x1 < a[jp] {
                        lbl = 120;
                        continue;
                    }
                    if *ip == np {
                        return;
                    }
                    jp += 2;
                    *ip += 1;
                    let it = nint(a[jr]);
                    if *ip <= it {
                        continue;
                    }
                    jr += 2;
                    *ir += 1;
                }
                120 => {
                    if x1 >= a[jp - 2] {
                        xlo = x1;
                        lbl = 150;
                        continue;
                    }
                    if *ip > 2 {
                        // 125
                        jp -= 2;
                        *ip -= 1;
                        if *ir == 1 {
                            lbl = 110;
                            continue;
                        }
                        let it = nint(a[jr - 2]);
                        if *ip > it {
                            lbl = 110;
                            continue;
                        }
                        jr -= 2;
                        *ir -= 1;
                        lbl = 110;
                        continue;
                    }
                    if x2 <= a[jp - 2] {
                        return;
                    }
                    xlo = a[jp - 2];
                    lbl = 150;
                }
                _ => {
                    // 150
                    let mut xhi = a[jp];
                    if xhi > x2 {
                        xhi = x2;
                    }
                    let int = nint(a[jr + 1]);
                    *f = *f + self.gral(a[jp - 2], a[jp - 1], a[jp], a[jp + 1], xlo, xhi, int);
                    if xhi == x2 {
                        return;
                    }
                    if *ip == np {
                        return;
                    }
                    xlo = xhi;
                    jp += 2;
                    *ip += 1;
                    let it = nint(a[jr]);
                    if *ip <= it {
                        continue;
                    }
                    jr += 2;
                    *ir += 1;
                }
            }
        }
    }

    // ------------------------------------------------------------ gral
    pub fn gral(&mut self, xl: f64, yl: f64, xh: f64, yh: f64, x1: f64, x2: f64, int: i32) -> f64 {
        let brk = 0.1e0f64;
        let zero = 0.0f64;
        let one = 1.0f64;
        let mut gral = 0.0;
        if x2 == x1 {
            return gral;
        }
        if x2 < x1 {
            self.error("gral", "x2 lt x1.", " ");
        }
        let lin = |xl: f64, yl: f64, xh: f64, yh: f64| {
            let b = (yh - yl) / (xh - xl);
            let a = yl - b * xl;
            (x2 - x1) * (a + b * (x2 + x1) / 2.0)
        };
        let linlog = |xl: f64, yl: f64, xh: f64, yh: f64| {
            let b = (yh - yl) / log(xh / xl);
            let z = (x2 - x1) / x1;
            if z.abs() <= brk {
                (x2 - x1) * (yl + b * log(x1 / xl))
                    + b * x1 * z * z * (1.0 + z * (-one / 3.0 + z * (one / 6.0 - z / 10.0))) / 2.0
            } else {
                (x2 - x1) * (yl + b * log(x1 / xl)) + b * x1 * (1.0 + (x2 / x1) * (log(x2 / x1) - 1.0))
            }
        };
        let loglin = |xl: f64, yl: f64, xh: f64, yh: f64| {
            let b = log(yh / yl) / (xh - xl);
            let a = log(yl) - b * xl;
            let z = (x2 - x1) * b;
            if z.abs() <= brk {
                exp(a + b * x1) * (x2 - x1) * (1.0 + z * (one / 2.0 + z / 6.0))
            } else {
                exp(a + b * x1) * (exp(z) - 1.0) / b
            }
        };
        if int == 1 {
            gral = (x2 - x1) * yl;
        } else if int == 2 {
            gral = lin(xl, yl, xh, yh);
        } else if int == 3 {
            if xl <= zero || xh <= zero {
                gral = lin(xl, yl, xh, yh);
            } else {
                gral = linlog(xl, yl, xh, yh);
            }
        } else if int == 4 {
            if yl < zero || yh < zero {
                gral = lin(xl, yl, xh, yh);
            } else {
                gral = loglin(xl, yl, xh, yh);
            }
        } else if int == 5 {
            if xl <= zero || xh <= zero {
                if yl < zero || yh < zero {
                    gral = lin(xl, yl, xh, yh);
                } else {
                    gral = loglin(xl, yl, xh, yh);
                }
            } else if yl < zero || yh < zero {
                if xl <= zero || xh <= zero {
                    gral = lin(xl, yl, xh, yh);
                } else {
                    gral = linlog(xl, yl, xh, yh);
                }
            } else {
                let b = log(yh / yl) / log(xh / xl);
                let z = (b + 1.0) * log(x2 / x1);
                if z.abs() <= brk {
                    gral = yl * x1 * pow(x1 / xl, b) * log(x2 / x1) * (1.0 + z * (one / 2.0 + z / 6.0));
                } else {
                    gral = yl * x1 * pow(x1 / xl, b) * (pow(x2 / x1, b + 1.0) - 1.0) / (b + 1.0);
                }
            }
        }
        gral
    }

    // ------------------------------------------------------------ gety1/2
    pub fn gety1(&mut self, x: f64, xnext: &mut f64, idis: &mut i32, y1: &mut f64, itape: i32, a: &mut F) {
        let mut s = self.g1;
        self.gety_impl(true, &mut s, x, xnext, idis, y1, itape, a);
        self.g1 = s;
    }
    pub fn gety2(&mut self, x: f64, xnext: &mut f64, idis: &mut i32, y2: &mut f64, itape: i32, a: &mut F) {
        let mut s = self.g2;
        self.gety_impl(false, &mut s, x, xnext, idis, y2, itape, a);
        self.g2 = s;
    }

    fn gety_impl(
        &mut self,
        first: bool,
        s: &mut GetySave,
        x: f64,
        xnext: &mut f64,
        idis: &mut i32,
        y1: &mut f64,
        itape: i32,
        a: &mut F,
    ) {
        let down = 0.999999e0f64;
        let xbig = 1.0e12f64;
        let zero = 0.0f64;
        *idis = 0;
        let mut ln = 0;
        let mut lbl;
        if x > zero {
            lbl = 110;
        } else {
            let (mut nb, mut nw) = (0, 0);
            self.tab1io(itape, 0, 0, a, &mut nb, &mut nw);
            s.nb = nb;
            s.nwtot = nw;
            s.nr = nint(a[5]);
            s.np = nint(a[6]);
            s.lt = 6 + 2 * s.nr;
            s.ip1 = 1;
            s.ip2 = (nw - s.lt) / 2;
            if s.nb == 0 {
                s.ip2 += 2;
            }
            s.ir = 1;
            s.ip = 1;
            if !first {
                // gety2
                s.xlast = a[s.lt + 1];
                s.ylast = a[s.lt + 2];
                *xnext = s.xlast;
                *y1 = s.ylast;
                if s.ylast != zero {
                    s.xlast = down * s.xlast;
                }
                return;
            }
            lbl = 100;
        }
        loop {
            match lbl {
                100 => {
                    let loc = s.lt + 2 * (s.ip - s.ip1) + 1;
                    s.xlast = a[loc];
                    s.ylast = a[loc + 1];
                    *xnext = s.xlast;
                    *y1 = s.ylast;
                    if s.ylast != zero {
                        s.xlast = down * s.xlast;
                    }
                    if s.ylast != zero {
                        return;
                    }
                    if s.ip >= s.np - 1 {
                        return;
                    }
                    if a[loc + 3] != zero && s.ip > 1 {
                        *xnext = down * *xnext;
                    }
                    if a[loc + 3] != zero {
                        return;
                    }
                    s.ip += 1;
                    s.nbt = nint(a[5 + 2 * s.ir]);
                    if s.ip > s.nbt {
                        s.ir += 1;
                    }
                    if s.ip + 2 <= s.ip2 {
                        continue;
                    }
                    self.gety_more(s, a, itape);
                    continue;
                }
                110 => {
                    if x < s.xlast {
                        // 200
                        *y1 = 0.0;
                        *xnext = s.xlast;
                        return;
                    }
                    ln = 2 * (s.ip - s.ip1) + s.lt;
                    if x < a[ln + 1] {
                        lbl = 120;
                        continue;
                    }
                    if s.ip == s.np {
                        // 300
                        *y1 = a[ln + 2];
                        s.xlast = xbig;
                        *xnext = xbig;
                        return;
                    }
                    s.xlast = a[ln + 1];
                    s.ylast = a[ln + 2];
                    s.ip += 1;
                    s.nbt = nint(a[5 + 2 * s.ir]);
                    if s.ip > s.nbt {
                        s.ir += 1;
                    }
                    if s.ip + 2 <= s.ip2 {
                        continue;
                    }
                    if s.nb == 0 {
                        lbl = 120;
                        continue;
                    }
                    self.gety_more(s, a, itape);
                    continue;
                }
                _ => {
                    // 120
                    let int = nint(a[2 * s.ir + 6]);
                    if int == 1 {
                        *idis = 1;
                    }
                    self.terp1(s.xlast, s.ylast, a[ln + 1], a[ln + 2], x, y1, int);
                    *xnext = a[ln + 1];
                    if ln + 3 > s.nwtot && s.nb == 0 {
                        return;
                    }
                    let xn = a[ln + 3];
                    if xn == *xnext {
                        *idis = 1;
                    }
                    return;
                }
            }
        }
    }

    fn gety_more(&mut self, s: &mut GetySave, a: &mut F, itape: i32) {
        let lt = s.lt;
        let nwtot = s.nwtot;
        a[lt + 1] = a[nwtot - 3];
        a[lt + 2] = a[nwtot - 2];
        a[lt + 3] = a[nwtot - 1];
        a[lt + 4] = a[nwtot];
        let (mut nb, mut nw) = (s.nb, 0);
        {
            let mut sub = a.at(lt + 5);
            self.moreio(itape, 0, 0, &mut sub, &mut nb, &mut nw);
        }
        s.nb = nb;
        s.nwtot = nw + lt + 4;
        s.ip1 = s.ip;
        s.ip2 = s.ip1 + nw / 2 + 1;
        if s.nb == 0 {
            s.ip2 += 2;
        }
    }
}

/// terp1 as a free function (thr6 passed explicitly).
pub fn terp1(x1: f64, y1: f64, x2: f64, y2: f64, x: f64, y: &mut f64, i: i32, thr6: f64) {
    let zero = 0.0f64;
    if x2 == x1 {
        *y = y1;
        return;
    }
    if i == 1 || y2 == y1 || x == x1 {
        *y = y1;
    } else if i == 2 {
        *y = y1 + (x - x1) * (y2 - y1) / (x2 - x1);
    } else if i == 3 {
        *y = y1 + log(x / x1) * (y2 - y1) / log(x2 / x1);
    } else if i == 4 {
        *y = y1 * exp((x - x1) * log(y2 / y1) / (x2 - x1));
    } else if i == 5 {
        if y1 == zero {
            *y = y1;
        } else {
            *y = y1 * exp(log(x / x1) * log(y2 / y1) / log(x2 / x1));
        }
    } else if i == 6 {
        if y1 == zero {
            *y = y1;
        } else {
            let t = sqrt(x1 - thr6);
            let mut b = log((x2 * y2) / (x1 * y1));
            b = b / (1.0 / t - 1.0 / sqrt(x2 - thr6));
            let a = exp(b / t) * x1 * y1;
            *y = (a / x) * exp(-b / sqrt(x - thr6));
        }
    }
}
