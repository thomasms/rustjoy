//! MIXR: a PENDF tape of mixed cross sections (linear combinations of
//! the MF3 sections of several materials), ported from NJOY2016
//! `mixr.f90` (module mixm).
//!
//! The routine keeps the Fortran structure and indices. One flat array
//! `a` holds the loada buffer, the output scratch area and one TAB1 page
//! for each input material, at the same offsets as the Fortran, so that
//! a page that grows past its area acts the same way.

#![allow(clippy::too_many_lines)]
#![allow(clippy::collapsible_if)]
#![allow(clippy::collapsible_else_if)]

use crate::endf::Ctx;
use crate::farr::{pack_text, F, RV};
use crate::ffmt;
use crate::fio::{Kind, Val};
use crate::fmath::nint;
use crate::util::sigfig;

const NNINMX: usize = 10;
const NMTMX: usize = 20;
const NMATMX: usize = NNINMX;
const MXLNA: i32 = 300000;
const BIG: f64 = 1.0e10;
const ZERO: f64 = 0.0;
/// Extra words past the Fortran array size.
const SLACK: usize = 20000;

/// Module variables of `mixm` and the saved local of `gety`.
#[derive(Default)]
pub struct Mixr {
    ntape: i32,
    jtape: [i32; NNINMX + 2],
    nrt: [i32; NNINMX + 2],
    npt: [i32; NNINMX + 2],
    irt: [i32; NNINMX + 2],
    ipt: [i32; NNINMX + 2],
    ip1t: [i32; NNINMX + 2],
    ip2t: [i32; NNINMX + 2],
    nbt: [i32; NNINMX + 2],
    nwt: [i32; NNINMX + 2],
    /// `save lt` in gety: one value for all tapes.
    lt: i32,
}

fn geti(v: &Option<Val>, d: i32) -> i32 {
    match v {
        Some(Val::Int(i)) => *i,
        Some(Val::Real(x)) => *x as i32,
        _ => d,
    }
}

fn getr(v: &Option<Val>, d: f64) -> f64 {
    match v {
        Some(Val::Real(x)) => *x,
        Some(Val::Int(i)) => *i as f64,
        _ => d,
    }
}

impl Mixr {
    pub fn new() -> Self {
        Self::default()
    }

    /// Construct a PENDF tape of mixed cross sections.
    pub fn mixr(&mut self, c: &mut Ctx) {
        let (mut nb, mut nw) = (0i32, 0i32);
        let mut a = RV::new(MXLNA as usize + SLACK);
        let mut cc = [0.0f64; 2];
        let mut des = [0.0f64; 17];
        let lrp = -1i32;

        let time = c.timer();
        c.io.so(&format!(
            "\n mixr...calculation of mixed cross sections{}{}s",
            " ".repeat(26),
            ffmt::f(time, 8, 1)
        ));
        c.io.se(&format!("\n mixr...{}{}s", " ".repeat(61), ffmt::f(time, 8, 1)));
        let nscr = 10;
        let nbuf = 5000i32;
        let nc = 2;
        c.openz(-nscr, 0);

        //--read user input
        let mut nin = [0i32; NNINMX + 1];
        let mut k: Vec<Kind> = vec![Kind::Int];
        k.extend(std::iter::repeat(Kind::Int).take(NNINMX));
        let v = c.io.input.read(&k);
        // nout is an uninitialised local in Fortran; 0 here
        let nout = geti(&v[0], 0);
        for i in 1..=NNINMX {
            nin[i] = geti(&v[i], nin[i]);
        }
        let mut nnin = 0usize;
        for i in 1..=NNINMX {
            if nin[i] != 0 {
                nnin = i;
            }
        }
        let mut mtn = [0i32; NMTMX + 1];
        let k: Vec<Kind> = std::iter::repeat(Kind::Int).take(NMTMX).collect();
        let v = c.io.input.read(&k);
        for i in 1..=NMTMX {
            mtn[i] = geti(&v[i - 1], mtn[i]);
        }
        let mut nmt = 0usize;
        for i in 1..=NMTMX {
            if mtn[i] != 0 {
                nmt = i;
            }
        }
        let mut matn = [0i32; NMATMX + 1];
        let mut wtn = [0.0f64; NMATMX + 1];
        let mut k: Vec<Kind> = Vec::new();
        for _ in 0..NMATMX {
            k.push(Kind::Int);
            k.push(Kind::Real);
        }
        let v = c.io.input.read(&k);
        for i in 1..=NMATMX {
            matn[i] = geti(&v[2 * i - 2], matn[i]);
            wtn[i] = getr(&v[2 * i - 1], wtn[i]);
        }
        let mut nmat = 0usize;
        for i in 1..=NMATMX {
            if matn[i] != 0 {
                nmat = i;
            }
        }
        let v = c.io.input.read(&[Kind::Real]);
        let temp = getr(&v[0], 0.0);
        let v = c.io.input.read(&[Kind::Int, Kind::Real, Kind::Real]);
        let matd = geti(&v[0], 0);
        let za = getr(&v[1], 0.0);
        let awr = getr(&v[2], 0.0);
        let mut strng = vec![b' '; 66];
        let v = c.io.input.read(&[Kind::Char(66)]);
        if let Some(Val::Char(s)) = &v[0] {
            let b = s.as_bytes();
            for (j, x) in b.iter().take(66).enumerate() {
                strng[j] = *x;
            }
        }
        // read(strng,'(16a4,a2)') str: 64 characters, then 2 and 2 blanks
        let mut t68 = strng[..66].to_vec();
        t68.push(b' ');
        t68.push(b' ');
        pack_text(&t68, &mut des);
        let ndes = 17;
        c.io.so(&format!("\n unit for output tape ................. {}", ffmt::i(nout, 10)));
        for i in 1..=nnin {
            c.io.so(&format!(" unit for input tape .................. {}", ffmt::i(nin[i], 10)));
        }
        let mut s = String::from(" mt numbers for output tape ........... \n          ");
        for i in 1..=nmt {
            s.push_str(&ffmt::i(mtn[i], 6));
        }
        c.io.so(&s);
        let mut s = String::from(" mat,weight values for mixes .......... ");
        for i in 1..=nmat {
            s.push_str(&format!("\n          {}{}", ffmt::i(matn[i], 6), ffmt::f(wtn[i], 10, 5)));
        }
        c.io.so(&s);

        //--allocate storage
        let ibuf = 1i32;
        let iscr = ibuf + nbuf;
        let mut jscr = [0i32; NMATMX + 1];
        jscr[1] = iscr + c.npage + 50;
        for i in 2..=nmat {
            jscr[i] = jscr[i - 1] + c.npage + 50;
        }

        //--open the output and input units
        //--and locate the desired materials and temperature.
        c.openz(nout, 1);
        let mut awi = 1.0f64;
        let mut emax = 20.0e6f64;
        let mut nsub = 10i32;
        for i in 1..=nmat {
            c.openz(nin[i], 0);
            c.tpidio(nin[i], 0, 0, &mut a.at(iscr), &mut nb, &mut nw);
            let mut idone = 0;
            while idone == 0 {
                c.contio(nin[i], 0, 0, &mut a.at(iscr), &mut nb, &mut nw);
                if c.math < 0 {
                    c.error("mixr", "mat and temp not found", " ");
                }
                if c.math == matn[i] {
                    c.contio(nin[i], 0, 0, &mut a.at(iscr), &mut nb, &mut nw);
                    if c.n1h != 0 {
                        c.iverf = 4;
                    } else if c.n2h == 0 {
                        c.iverf = 5;
                    } else {
                        c.iverf = 6;
                    }
                    c.skiprz(nin[i], -1);
                    if c.iverf > 4 {
                        c.contio(nin[i], 0, 0, &mut a.at(iscr), &mut nb, &mut nw);
                    }
                    if c.iverf > 5 {
                        c.contio(nin[i], 0, 0, &mut a.at(iscr), &mut nb, &mut nw);
                        if a[iscr + 1] > emax {
                            emax = a[iscr + 1];
                        }
                        if i == 1 {
                            awi = a[iscr];
                            if awi < 0.1 {
                                nsub = 0;
                            }
                            if awi > 0.9980 && awi < 0.9990 {
                                nsub = 10010;
                            }
                            if awi > 1.9950 && awi < 1.9970 {
                                nsub = 10020;
                            }
                            if awi > 2.9895 && awi < 2.9897 {
                                nsub = 10030;
                            }
                            if awi > 2.9890 && awi < 2.9891 {
                                nsub = 20030;
                            }
                            if awi > 3.9670 && awi < 3.9680 {
                                nsub = 20040;
                            }
                        }
                    }
                    c.hdatio(nin[i], 0, 0, &mut a.at(iscr), &mut nb, &mut nw);
                    let t = a[iscr];
                    let test = 1.0f64;
                    if (t - temp).abs() < test {
                        idone = 1;
                    } else {
                        c.tomend(nin[i], 0, 0, &mut a.at(iscr));
                    }
                } else {
                    c.tomend(nin[i], 0, 0, &mut a.at(iscr));
                }
            }
            c.findf(matn[i], 3, 0, nin[i]);
        }

        //--write the tape id record and file 1 for the output tape
        c.math = 1;
        c.afend(nout, 0);
        c.math = matd;
        c.mfh = 1;
        c.mth = 451;
        let nx = nmt as i32;
        a[iscr] = za;
        a[iscr + 1] = awr;
        a[iscr + 2] = lrp as f64;
        a[iscr + 3] = 0.0;
        a[iscr + 4] = 0.0;
        a[iscr + 5] = 0.0;
        if c.iverf <= 4 {
            a[iscr + 5] = nx as f64;
        }
        c.contio(0, nout, 0, &mut a.at(iscr), &mut nb, &mut nw);
        if c.iverf > 4 {
            for j in 0..6 {
                a[iscr + j] = 0.0;
            }
            if c.iverf == 6 {
                a[iscr + 5] = 6.0;
            }
            c.contio(0, nout, 0, &mut a.at(iscr), &mut nb, &mut nw);
            if c.iverf >= 6 {
                a[iscr] = awi;
                a[iscr + 1] = emax;
                a[iscr + 2] = 0.0;
                a[iscr + 3] = 0.0;
                a[iscr + 4] = nsub as f64;
                a[iscr + 5] = 0.0;
                c.contio(0, nout, 0, &mut a.at(iscr), &mut nb, &mut nw);
            }
        }
        a[iscr] = temp;
        a[iscr + 1] = 0.0;
        a[iscr + 2] = 0.0;
        if c.iverf >= 6 {
            a[iscr + 2] = 1.0;
        }
        a[iscr + 3] = 0.0;
        a[iscr + 4] = ndes as f64;
        a[iscr + 5] = 0.0;
        if c.iverf >= 5 {
            a[iscr + 5] = nx as f64;
        }
        for i in 1..=ndes {
            a[iscr + 5 + i] = des[(i - 1) as usize];
        }
        c.hdatio(0, nout, 0, &mut a.at(iscr), &mut nb, &mut nw);
        let mut inow = iscr;
        for i in 1..=nmt {
            a[inow] = 0.0;
            a[inow + 1] = 0.0;
            a[inow + 2] = 3.0;
            a[inow + 3] = mtn[i] as f64;
            a[inow + 4] = 0.0;
            a[inow + 5] = 0.0;
            inow += 6;
        }
        nw = nmt as i32;
        c.dictio(0, nout, 0, &mut a.at(iscr), &mut nb, &mut nw);
        c.afend(nout, 0);

        //--loop over the desired reactions.
        for imt in 1..=nmt {
            let mtd = mtn[imt];

            //--find the section for mt on each input tape.
            let mut enext = BIG;
            let reset = -1.0f64;
            let (mut enxt, mut y) = (0.0f64, 0.0f64);
            let mut idis = 0i32;
            self.gety(c, reset, &mut enxt, &mut idis, &mut y, 0, &mut a.at(iscr));
            for i in 1..=nmat {
                let mut idone = 0;
                while idone == 0 {
                    c.contio(nin[i], 0, 0, &mut a.at(iscr), &mut nb, &mut nw);
                    if c.mfh == 0 || c.mth == mtd {
                        idone = 1;
                    } else {
                        c.tosend(nin[i], 0, 0, &mut a.at(iscr));
                    }
                }
                if c.mth == mtd {
                    let e = 0.0f64;
                    self.gety(c, e, &mut enxt, &mut idis, &mut y, nin[i], &mut a.at(jscr[i]));
                    if enxt < enext {
                        enext = enxt;
                    }
                    c.io.so(&format!(
                        " found section for mt={} and mat={}",
                        ffmt::i(mtd, 4),
                        ffmt::i(matn[i], 4)
                    ));
                } else {
                    let strng = format!("mt={} not present for mat={}", ffmt::i(mtd, 3), ffmt::i(matn[i], 4));
                    c.mess("mixr", &strng, " ");
                    c.findf(matn[i], 3, 0, nin[i]);
                }
            }

            //--loop over the energies to get a union grid
            let mut je = 0i32;
            while enext < BIG {
                let e = enext;

                //--mix the cross sections at this energy
                enext = BIG;
                let mut sig = 0.0f64;
                for i in 1..=nmat {
                    self.gety(c, e, &mut enxt, &mut idis, &mut y, nin[i], &mut a.at(jscr[i]));
                    if enxt < enext {
                        enext = enxt;
                    }
                    sig = sig + wtn[i] * y;
                }

                //--save the mixed value on the scratch tape
                je += 1;
                if enext >= BIG {
                    je = -je;
                }
                cc[0] = e;
                cc[1] = sig;
                let lo = (ibuf - 1) as usize;
                c.loada(je, &cc, nc, nscr, &mut a.0[lo..lo + nbuf as usize], nbuf);
            }
            let ne = -je;
            for i in 1..=nmat {
                c.tosend(nin[i], 0, 0, &mut a.at(iscr));
            }

            //--write head record for this mt
            c.math = matd;
            c.mfh = 3;
            c.mth = mtd;
            a[iscr] = za;
            a[iscr + 1] = awr;
            for j in 2..6 {
                a[iscr + j] = 0.0;
            }
            c.contio(0, nout, 0, &mut a.at(iscr), &mut nb, &mut nw);

            //--write tab1 record for this mt
            a[iscr] = 0.0;
            a[iscr + 1] = 0.0;
            a[iscr + 2] = 0.0;
            a[iscr + 3] = 0.0;
            a[iscr + 4] = 1.0;
            a[iscr + 5] = ne as f64;
            a[iscr + 6] = ne as f64;
            a[iscr + 7] = 2.0;
            let mut inow = iscr + 8;
            nb = 0;
            let mut jend = c.npage / 2;
            if jend > ne {
                jend = ne;
            }
            let mut je = 0i32;
            let mut idone = 0;
            while je < ne && idone == 0 {
                je += 1;
                let lo = (ibuf - 1) as usize;
                {
                    let (buf, _) = a.0[lo..].split_at_mut(nbuf as usize);
                    c.finda(je, &mut cc, nc, nscr, buf, nbuf);
                }
                a[inow] = cc[0];
                a[inow + 1] = sigfig(cc[1], 7, 0);
                inow += 2;
                if inow > MXLNA {
                    c.error("mixr", "mxlna array limit exceeded", " ");
                }
                if je >= jend {
                    if nb == 0 {
                        c.tab1io(0, nout, 0, &mut a.at(iscr), &mut nb, &mut nw);
                    } else {
                        c.moreio(0, nout, 0, &mut a.at(iscr), &mut nb, &mut nw);
                    }
                    if nb == 0 {
                        idone = 1;
                    } else {
                        inow = iscr;
                        jend += c.npage / 2;
                        if jend > ne {
                            jend = ne;
                        }
                    }
                }
            }
            c.asend(nout, 0);
        }

        //--finish off the output tape
        c.afend(nout, 0);
        c.amend(nout, 0);
        c.atend(nout, 0);

        c.closz(nout);
        for i in 1..=nmat {
            c.closz(nin[i]);
        }
        c.closz(-nscr);
        let time = c.timer();
        c.io.so(&format!("\n{}{}s\n{}", " ".repeat(69), ffmt::f(time, 8, 1), "*".repeat(77)));
    }

    /// `gety`: y(x) from a TAB1 record read page by page from `itape`,
    /// with the reading state kept for each tape. `x=0` starts a tape,
    /// `x<0` resets all tapes.
    #[allow(clippy::too_many_arguments)]
    fn gety(&mut self, c: &mut Ctx, x: f64, xnext: &mut f64, idis: &mut i32, y: &mut f64, itape: i32, a: &mut F) {
        let (nr, np, mut ir, mut ip, mut ip1, mut ip2, mut nb, mut nwtot);
        let mut nw = 0i32;
        let ktape: usize;

        *idis = 0;
        if x == ZERO {
            self.ntape += 1;
            let nt = self.ntape as usize;
            self.jtape[nt] = itape;
            nb = 0;
            c.tab1io(itape, 0, 0, a, &mut nb, &mut nw);
            nwtot = nw;
            nr = nint(a[5]);
            np = nint(a[6]);
            self.lt = 6 + 2 * nr;
            ip1 = 1;
            ip2 = (nw - self.lt) / 2;
            if nb == 0 {
                ip2 += 2;
            }
            ir = 1;
            ip = 2;
            *xnext = a[self.lt + 1];

            self.nrt[nt] = nr;
            self.npt[nt] = np;
            self.irt[nt] = ir;
            self.ipt[nt] = ip;
            self.ip1t[nt] = ip1;
            self.ip2t[nt] = ip2;
            self.nbt[nt] = nb;
            self.nwt[nt] = nwtot;
            return;
        }
        if !(x > ZERO) {
            self.ntape = 0;
            return;
        }

        if self.ntape == 0 {
            c.error("gety", "not properly initialized", " ");
        }
        let mut found = 0usize;
        for i in 1..=self.ntape as usize {
            if self.jtape[i] == itape {
                found = i;
                break;
            }
        }
        if found == 0 {
            *y = 0.0;
            *xnext = BIG;
            return;
        }
        ktape = found;
        nr = self.nrt[ktape];
        np = self.npt[ktape];
        ir = self.irt[ktape];
        ip = self.ipt[ktape];
        ip1 = self.ip1t[ktape];
        ip2 = self.ip2t[ktape];
        nb = self.nbt[ktape];
        nwtot = self.nwt[ktape];
        let lt = self.lt;

        let mut ln;
        // 0: label 130, 1: label 135, 2: label 140
        let go;
        loop {
            // label 125
            ln = 2 * (ip - ip1) + lt;
            if x < a[ln - 1] {
                go = 1;
                break;
            }
            if x < a[ln + 1] {
                go = 0;
                break;
            }
            if ip == np {
                go = 2;
                break;
            }

            //--no.  move up to next range.
            //--read in new page of data if needed.
            ip += 1;
            let nbx = nint(a[5 + 2 * ir]);
            if ip > nbx {
                ir += 1;
            }
            if ip < ip2 {
                continue;
            }
            if nb == 0 {
                go = 0;
                break;
            }
            a[lt + 1] = a[nwtot - 3];
            a[lt + 2] = a[nwtot - 2];
            a[lt + 3] = a[nwtot - 1];
            a[lt + 4] = a[nwtot];
            c.moreio(itape, 0, 0, &mut a.at(lt + 5), &mut nb, &mut nw);
            nwtot = nw + lt + 4;
            ip1 = ip - 1;
            ip2 = ip1 + nw / 2 + 1;
            if nb == 0 {
                ip2 += 2;
            }
        }
        match go {
            0 => {
                // label 130
                let int = nint(a[6 + 2 * ir]);
                if int == 1 {
                    *idis = 1;
                }
                c.terp1(a[ln - 1], a[ln], a[ln + 1], a[ln + 2], x, y, int);
                *xnext = a[ln + 1];
                if (ln + 3) > nwtot && nb == 0 {
                    return;
                }
                let xn = a[ln + 3];
                if xn == *xnext {
                    *idis = 1;
                }
            }
            1 => {
                // label 135
                *y = 0.0;
                *xnext = a[ln - 1];
            }
            _ => {
                // label 140
                *y = a[ln + 2];
                *xnext = BIG;
            }
        }

        // label 150
        self.nrt[ktape] = nr;
        self.npt[ktape] = np;
        self.irt[ktape] = ir;
        self.ipt[ktape] = ip;
        self.ip1t[ktape] = ip1;
        self.ip2t[ktape] = ip2;
        self.nbt[ktape] = nb;
        self.nwt[ktape] = nwtot;
    }
}
