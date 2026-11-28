//! POWR lib=3, part 2: subfile 4 (cross sections), subfile 5 (burnup
//! data), subfile 6 (P1 matrices) and the CLIB input on nout.

use super::cpm::EFLAG;
use super::{alloc_rv, dealloc, getr, pe, rows, rowsj, Powr, SLACK};
use crate::endf::Ctx;
use crate::farr::{IV, RV};
use crate::ffmt;
use crate::fio::{fatal_runtime, Kind, RecW};
use crate::fmath::nint;
use crate::util::physics::AMASSN;

const ZERO: f64 = 0.0;

/// The real(kr) word that `read('end ','(a4)') flag` gives: "end "
/// followed by four blanks.
pub(super) fn end_flag() -> f64 {
    f64::from_le_bytes(*b"end     ")
}

/// The first four characters of a real(kr) word holding text (`a4`).
fn a4_of(x: f64) -> String {
    let b = x.to_le_bytes();
    String::from_utf8_lossy(&b[..4]).to_string()
}

fn wr_reals(c: &mut Ctx, nun: i32, v: &[f64]) {
    let mut r = RecW::new();
    r.rs(v);
    c.io.write_rec(nun, r);
}

fn wr_ints(c: &mut Ctx, nun: i32, v: &[i32]) {
    let mut r = RecW::new();
    for &x in v {
        r.i(x);
    }
    c.io.write_rec(nun, r);
}

fn rd_reals(c: &mut Ctx, nun: i32, n: i32) -> Vec<f64> {
    c.io.read_rec_with(nun, |rr| (0..n.max(0)).map(|_| rr.r()).collect())
}

fn rd_ints(c: &mut Ctx, nun: i32, n: i32) -> Vec<i32> {
    c.io.read_rec_with(nun, |rr| (0..n.max(0)).map(|_| rr.i()).collect())
}

impl Powr {
    /// Process cross sections.
    pub(super) fn sfile4(&mut self, c: &mut Ctx, jtemp: &mut i32, nt: i32) {
        let dilinf = 1.0e10f64;
        let lz = 6;
        let ng = self.ngendf;
        let ngnd = self.ngnd;
        let (mscr, kscr) = (self.mscr, self.kscr);
        let (mut nb, mut nw) = (0i32, 0i32);

        //--allocate storage.
        if self.if4 != 0 {
            c.openz(-mscr, 1);
            c.repoz(-mscr);
            c.openz(-kscr, 1);
            c.repoz(-kscr);
            if self.nres == 0 {
                c.io.so("\n");
            }
            if self.nres > 0 {
                c.io.so("\n\n");
            }
            c.io.so("\n ***subfile 4***cross sections");
            alloc_rv(&mut self.al_fiss, &mut self.fiss, self.newmat, "fiss");
            if self.iu == 1 {
                alloc_rv(&mut self.al_uff, &mut self.uff, self.nfsg, "uff");
            }
            if self.ipu == 1 {
                alloc_rv(&mut self.al_puff, &mut self.puff, self.nfsg, "puff");
            }
            let g = ngnd as usize + SLACK;
            let mut snu = RV::new(g);
            let mut l1 = IV::new(g);
            let mut l2 = IV::new(g);
            let mut l1e = IV::new(g);
            let mut l2e = IV::new(g);
            let ngndsq = ngnd * ngnd;
            alloc_rv(&mut self.al_elas, &mut self.elas, ngndsq, "elas");
            let mut cspc = RV::new(self.ngn.max(0) as usize + 200 + SLACK);
            let mut abs1 = RV::new(g);
            let mut abs2 = RV::new(g);
            let mut sn2n = RV::new(g);
            let mut tot = RV::new(g);
            let mut xtr = RV::new(g);
            nw = nt + 1;
            if self.nres == 0 || cfg!(feature = "powr-cards") {
                alloc_rv(&mut self.al_tempt, &mut self.tempt, nw, "tempt");
            }
            let nwloc = ngnd * 6 + ngndsq;
            let mut a = RV::new(nwloc as usize + SLACK);
            let iloc = 1;
            let locab0 = iloc;
            let locsf0 = locab0 + ngnd;
            let locchi = locsf0 + ngnd;
            let locnus = locchi + ngnd;
            let locsg = locnus + ngnd;
            let locxtr = locsg + ngnd;
            let locxs = locxtr + ngnd;
            self.nfiss[1] = 0;
            self.nfiss[2] = 0;
            let il = 1;
            let lim = locxs + ngndsq - 1;
            let ngnd5 = ngnd * 5;

            //--search for desired materials.
            c.repoz(ng);
            c.tpidio(ng, 0, 0, &mut self.scr.at(1), &mut nb, &mut nw);
            let mut nmat = 0;
            let (mut matd, mut ident, mut ip1opt, mut mti, mut mtc, mut ntemp) = (0, 0, 0, 0, 0, 0);
            let mut sgref = 0.0f64;
            let (mut isg, mut mtabs1, mut mtabs2, mut i318, mut i618, mut mttot, mut mtel, mut mtex) =
                (0, 0, 0, 0, 0, 0, 0, 0);
            let (mut isof, mut nina, mut mf3mt2, mut jfiss, mut nth1, mut nth, mut if6) = (0, 0, 0, 0, 0, 0, 0);
            let (mut cnorm, mut dnorm) = (0.0f64, 0.0f64);
            let (mut nl, mut nz, mut ntw, mut l, mut iz, mut jic) = (0, 0, 0, 1, 0, 0);
            let (mut ng2, mut ig2lo, mut ig, mut igc, mut jg) = (0, 0, 0, 0, 0);
            let mut lbl = 120;
            loop {
                match lbl {
                    120 => {
                        nmat += 1;
                        matd = self.mat[nmat];
                        ident = self.idnt[nmat];
                        ip1opt = self.ip1tab[nmat];
                        mti = self.mtit[nmat];
                        mtc = self.mtct[nmat];
                        ntemp = self.ntmp[nmat];
                        sgref = self.s0rf[nmat];
                        isg = 0;
                        if (sgref - dilinf).abs() >= sgref / 100.0 {
                            isg = 1;
                        }
                        mtabs1 = 0;
                        mtabs2 = 0;
                        i318 = 0;
                        i618 = 0;
                        mttot = 0;
                        mtel = 0;
                        mtex = 0;
                        for i in 1..=ngnd {
                            abs2[i] = 0.0;
                            sn2n[i] = 0.0;
                            a[i - 1 + locxtr] = 0.0;
                            l1e[i] = ngnd;
                            l2e[i] = 1;
                            l1[i] = ngnd;
                            l2[i] = 1;
                        }
                        isof = 0;
                        if ident == self.iu235 {
                            isof = 1;
                        }
                        if ident == self.ipu239 {
                            isof = 2;
                        }
                        nina = self.ninat[nmat];
                        mf3mt2 = 0;
                        cnorm = 0.0;
                        dnorm = 0.0;
                        jfiss = 0;
                        *jtemp = -1;
                        nth1 = 0;
                        nth = 0;
                        if6 = 0;
                        for i in 1..=nwloc {
                            a[i - 1 + iloc] = 0.0;
                        }
                        if nina < 2 {
                            lbl = 140;
                            continue;
                        }
                        // read in absorption
                        let nzr = ngnd + 2;
                        let kinds = vec![Kind::Real; nzr as usize];
                        let v = c.io.input.read(&kinds);
                        for i in 1..=nzr {
                            self.scr[i] = getr(&v[(i - 1) as usize], self.scr[i]);
                        }
                        self.awrt[nmat] = self.scr[1];
                        self.tempt[1] = self.scr[2];
                        for i in 3..=nzr {
                            a[i - 3 + locab0] = self.scr[i];
                        }
                        if nina >= 3 {
                            let kinds = vec![Kind::Real; ngnd as usize];
                            for base in [locsf0, locnus, locxtr] {
                                let v = c.io.input.read(&kinds);
                                for i in 1..=ngnd {
                                    a[base + i - 1] = getr(&v[(i - 1) as usize], a[base + i - 1]);
                                }
                            }
                            for ia in 1..=ngnd {
                                let v = c.io.input.read(&[Kind::Real; 3]);
                                for i in 1..=3 {
                                    self.scr[i] = getr(&v[(i - 1) as usize], self.scr[i]);
                                }
                                let ib = nint(self.scr[1]);
                                if ib != 0 {
                                    l1[ia] = nint(self.scr[2]);
                                    l2[ia] = nint(self.scr[3]);
                                    let lone = nint(self.scr[2]);
                                    let ltwo = nint(self.scr[3]);
                                    let nzs = ltwo - lone + 1;
                                    let kinds = vec![Kind::Real; nzs.max(0) as usize];
                                    let v = c.io.input.read(&kinds);
                                    for i in 1..=nzs {
                                        self.scr[i] = getr(&v[(i - 1) as usize], self.scr[i]);
                                    }
                                    let loc = locxs + ib - 1 - ngnd;
                                    let mut k = 0;
                                    for ja in lone..=ltwo {
                                        k += 1;
                                        a[ngnd * ja + loc] = self.scr[k];
                                    }
                                } else {
                                    l1[ia] = ngnd;
                                    l2[ia] = 1;
                                }
                            }
                            *jtemp = 1;
                        } else {
                            *jtemp = -1;
                            matd = -1;
                        }
                        lbl = 500;
                    }
                    140 => {
                        c.contio(ng, 0, 0, &mut self.scr.at(1), &mut nb, &mut nw);
                        if c.math == matd {
                            lbl = 150;
                            continue;
                        }
                        if c.math == 0 {
                            lbl = 140;
                            continue;
                        }
                        if c.math > 0 {
                            c.tomend(ng, 0, 0, &mut self.scr.at(1));
                            lbl = 140;
                            continue;
                        }
                        c.repoz(ng);
                        c.tpidio(ng, 0, 0, &mut self.scr.at(1), &mut nb, &mut nw);
                        lbl = 140;
                    }
                    //--search for required reactions.
                    150 => {
                        nl = c.l1h;
                        nz = c.l2h;
                        ntw = c.n2h;
                        if c.mfh == 3 && c.mth >= 2 && mf3mt2 <= 0 {
                            // invert flux arrays
                            mf3mt2 = 1;
                            for i in 1..=ngnd {
                                if self.cflux[i] != ZERO {
                                    self.cflux[i] = 1.0 / self.cflux[i];
                                }
                                if self.cflux[ngnd + i] != ZERO {
                                    self.cflux[ngnd + i] = 1.0 / self.cflux[ngnd + i];
                                }
                                let mut loc = i;
                                if self.flux[loc] != ZERO {
                                    self.flux[loc] = 1.0 / self.flux[loc];
                                }
                                loc += ngnd;
                                if self.flux[loc] != ZERO {
                                    self.flux[loc] = 1.0 / self.flux[loc];
                                }
                            }
                        }
                        lbl = 160;
                    }
                    160 => {
                        l = 1;
                        c.listio(ng, 0, 0, &mut self.scr.at(l), &mut nb, &mut nw);
                        if c.mfh != 1 {
                            lbl = 170;
                            continue;
                        }
                        let temp = c.c1h;
                        self.ngn = c.l1h;
                        if self.ngn != ngnd && self.ngn != self.ngnf {
                            c.error("sfile4", "storage exceeded.", " ");
                        }
                        let mut idone = 0;
                        let mut i = 0;
                        while i < nz && idone == 0 {
                            i += 1;
                            iz = i;
                            if (sgref - self.scr[l + 5 + ntw + i]).abs() < sgref / 100.0 {
                                idone = 1;
                            }
                        }
                        *jtemp += 1;
                        jic = 0;
                        if *jtemp == 0 {
                            lbl = 195;
                            continue;
                        }
                        if ntemp > 0 && *jtemp > ntemp {
                            lbl = 320;
                            continue;
                        }
                        self.tempt[*jtemp] = temp;
                        lbl = 195;
                    }
                    170 => {
                        ng2 = c.l1h;
                        ig2lo = c.l2h;
                        ig = c.n2h;
                        igc = 1;
                        if ig != 0 {
                            igc = self.icgrp(ig);
                        }
                        // check for presence of sometimes temperature-dependent data
                        if *jtemp != 0 {
                            let mth = c.mth;
                            let mut skip = false;
                            if !(mth < 102 || mth > 150) {
                                if mth == 102 {
                                    if mtabs1 > 0 {
                                        skip = true;
                                    } else {
                                        mtabs1 = 1;
                                        for i in 1..=ngnd {
                                            abs1[i] = 0.0;
                                        }
                                    }
                                } else {
                                    // 172
                                    if mtabs2 > 0 {
                                        skip = true;
                                    } else {
                                        mtabs2 = 1;
                                        for i in 1..=ngnd {
                                            abs2[i] = 0.0;
                                        }
                                    }
                                }
                            }
                            if !skip {
                                // 174
                                if mth == 2 {
                                    if c.mfh != 6 {
                                        skip = true;
                                    } else {
                                        if mtel <= 0 {
                                            mtel = 1;
                                            for i in 1..=ngnd {
                                                l1e[i] = ngnd;
                                                l2e[i] = 1;
                                            }
                                            for i in 1..=ngndsq {
                                                self.elas[i] = 0.0;
                                            }
                                        }
                                        // 177
                                        if mtex > 0 {
                                            skip = true;
                                        } else {
                                            mtex = 1;
                                            for i in 1..=ngnd {
                                                xtr[i] = 0.0;
                                            }
                                        }
                                    }
                                }
                            }
                            if !skip {
                                // 179
                                if mth == 1 && mttot <= 0 {
                                    mttot = 1;
                                    for i in 1..=ngnd {
                                        tot[i] = 0.0;
                                    }
                                }
                            }
                        }
                        // 185
                        while nb != 0 {
                            l += nw;
                            c.moreio(ng, 0, 0, &mut self.scr.at(l), &mut nb, &mut nw);
                            if l + nw - 1 > self.nwscr {
                                c.error("sfile4", "storage exceeded.", " ");
                            }
                        }
                        if igc == 0 || igc > ngnd {
                            lbl = 218;
                            continue;
                        }
                        jg = ngnd - igc + 1;
                        l = 1;
                        if c.mth == 455 && *jtemp > 0 {
                            lbl = 240;
                            continue;
                        }
                        if c.mfh == 3 {
                            lbl = 205;
                            continue;
                        }
                        if c.mfh == 6 {
                            lbl = 255;
                            continue;
                        }
                        lbl = 195;
                    }
                    195 => {
                        c.tosend(ng, 0, 0, &mut self.scr.at(l));
                        lbl = 200;
                    }
                    200 => {
                        c.contio(ng, 0, 0, &mut self.scr.at(l), &mut nb, &mut nw);
                        if c.math == 0 {
                            lbl = 330;
                            continue;
                        }
                        if c.math == -1 {
                            lbl = 350;
                            continue;
                        }
                        if c.math != matd {
                            lbl = 350;
                            continue;
                        }
                        if c.mfh == 0 || c.mth == 0 {
                            lbl = 200;
                            continue;
                        }
                        lbl = 150;
                    }
                    //--file 3 cross sections.
                    205 => {
                        if *jtemp == 0 {
                            lbl = 220;
                            continue;
                        }
                        if c.mth < 102 || c.mth > 150 {
                            lbl = 214;
                            continue;
                        }
                        // absorption
                        let mut loca = l + lz + nl * nz;
                        if isg > 0 && nz >= iz {
                            loca = l + lz + nl * (iz - 1 + nz);
                        }
                        let mut locf = l + lz;
                        if isg > 0 && nz >= iz {
                            locf = l + lz + nl * (iz - 1);
                        }
                        let mut jgc = jg;
                        if isg > 0 && nz >= iz {
                            jgc = jg + ngnd;
                        }
                        if c.mth != 102 {
                            abs2[jg] = abs2[jg] + self.scr[loca] * self.scr[locf] * self.cflux[jgc];
                        } else {
                            abs1[jg] = abs1[jg] + self.scr[loca] * self.scr[locf] * self.cflux[jgc];
                        }
                        if ig < self.ngn {
                            lbl = 160;
                            continue;
                        }
                        if nina > 0 {
                            c.tomend(ng, 0, 0, &mut self.scr.at(1));
                        }
                        if nina > 0 {
                            lbl = 400;
                            continue;
                        }
                        lbl = 200;
                    }
                    214 => {
                        let mth = c.mth;
                        if mth == 1 {
                            lbl = 220;
                            continue;
                        }
                        if nina > 0 {
                            lbl = 195;
                            continue;
                        }
                        if mth == 16 {
                            lbl = 220;
                            continue;
                        }
                        if (mth < 18 || mth > 21) && mth != 38 {
                            lbl = 220;
                            continue;
                        }
                        if mth == 18 {
                            i318 = 1;
                        }
                        if !(mth == 18 || i318 == 0) {
                            lbl = 195;
                            continue;
                        }
                        // fission (216)
                        let loc = locsf0 + jg - 1;
                        let loca = l + lz + nl * nz;
                        let mut locf = l + lz;
                        a[loc] = a[loc] + self.scr[loca] * self.scr[locf] * self.cflux[jg];
                        let loc = locab0 + jg - 1;
                        if isg > 0 && nz >= iz {
                            locf = l + lz + nl * (iz - 1);
                        }
                        let mut jgc = jg;
                        if isg > 0 && nz >= iz {
                            jgc = jg + ngnd;
                        }
                        a[loc] = a[loc] + self.scr[loca] * self.scr[locf] * self.cflux[jgc];
                        lbl = 218;
                    }
                    218 => {
                        lbl = if ig < self.ngn { 160 } else { 200 };
                    }
                    220 => {
                        if c.mth != 1 && *jtemp > 0 {
                            lbl = 195;
                            continue;
                        }
                        if c.mth != 1 && *jtemp == 0 {
                            lbl = 236;
                            continue;
                        }
                        if ig == 1 {
                            for i in 1..=ngnd {
                                self.cflux[i] = 0.0;
                                self.cflux[i + ngnd] = 0.0;
                                self.flux[i] = 0.0;
                                self.flux[i + ngnd] = 0.0;
                            }
                            mf3mt2 = 0;
                        }
                        // p1 flux, for infinite dilution
                        let mut loc = jg;
                        let mut loca = l + lz;
                        if nl > 1 {
                            loca += 1;
                        }
                        self.flux[loc] = self.flux[loc] + self.scr[loca];
                        if isg != 0 && nz >= iz {
                            // for reference sigma zero
                            loc = jg + ngnd;
                            loca = l + lz + nl * (iz - 1);
                            if nl > 1 {
                                loca += 1;
                            }
                            self.flux[loc] = self.flux[loc] + self.scr[loca];
                        }
                        // p0 flux, for infinite dilution
                        let loca = l + lz;
                        self.cflux[jg] = self.cflux[jg] + self.scr[loca];
                        // for reference sigma zero
                        if !(isg == 0 || nz < iz) {
                            self.cflux[jg + ngnd] = self.cflux[jg + ngnd] + self.scr[l + lz + nl * (iz - 1)];
                        }
                        lbl = 250;
                    }
                    // n2n
                    236 => {
                        if *jtemp > 0 {
                            lbl = 195;
                            continue;
                        }
                        if c.mth == 16 {
                            let loca = l + lz + nl * nz;
                            let locf = l + lz;
                            sn2n[jg] = sn2n[jg] + self.scr[loca] * self.scr[locf] * self.cflux[jg];
                            lbl = 250;
                            continue;
                        }
                        // 238
                        if c.mth != mti {
                            lbl = 195;
                            continue;
                        }
                        if igc > nth1 && self.scr[l + lz + nz * nl] != 0.0 {
                            nth1 = igc;
                        }
                        nth = ngnd - nth1 + 1;
                        lbl = 250;
                    }
                    240 => {
                        if c.mth != 455 || isof == 0 {
                            lbl = 195;
                            continue;
                        }
                        jfiss = 1;
                        if c.mfh != 5 {
                            // delayed neutrons
                            let loc = locnus + jg - 1;
                            let loca = l + lz;
                            let locf = l + lz;
                            a[loc] = a[loc]
                                + self.scr[loca + 1] * self.scr[loca + 2] * self.scr[locf] * self.cflux[jg];
                            dnorm = dnorm + self.scr[loca] * self.scr[loca + 1] * self.scr[loca + 2];
                            lbl = if ig < self.ngn { 160 } else { 200 };
                            continue;
                        }
                        // 244
                        if isof == 0 {
                            lbl = 195;
                            continue;
                        }
                        for i in 2..=ng2 {
                            let jg2 = self.ngn - ig2lo - i + 3;
                            let jg2c = self.icgrp(jg2);
                            if jg2c != 0 && jg2c <= ngnd {
                                let loc = locchi + jg2c - 1;
                                let loca = l + lz - 1 + nl * (i - 1);
                                for j in 1..=nl {
                                    a[loc] = self.scr[j + loca] * dnorm + a[loc];
                                    cnorm = self.scr[j + loca] * dnorm + cnorm;
                                }
                            }
                        }
                        lbl = 250;
                    }
                    250 => {
                        lbl = if ig < self.ngn { 160 } else { 200 };
                    }
                    //--file 6 matrices.
                    255 => {
                        if *jtemp == 0 {
                            lbl = 300;
                            continue;
                        }
                        let mth = c.mth;
                        if (mth < 18 || mth > 21) && mth != 38 {
                            lbl = 270;
                            continue;
                        }
                        if mth == 18 {
                            i618 = 1;
                        }
                        if !(mth == 18 || i618 == 0) {
                            lbl = 195;
                            continue;
                        }
                        // fission (260)
                        jfiss = 1;
                        if ig != 0 {
                            let locf = l + lz;
                            for i in 2..=ng2 {
                                let loca = l + lz + (il - 1) + nl * nz * (i - 1);
                                let loc = locnus + jg - 1;
                                a[loc] = a[loc] + self.scr[loca] * self.scr[locf] * self.cflux[jg];
                                if isof != 0 {
                                    if ig2lo != 0 {
                                        // matrix part
                                        let ig2 = ig2lo + i - 2;
                                        let ig2c = self.icgrp(ig2);
                                        if ig2c != 0 && ig2c <= ngnd {
                                            let jg2c = ngnd - ig2c + 1;
                                            let loc = locchi + jg2c - 1;
                                            a[loc] = a[loc] + self.scr[loca] * self.scr[locf];
                                            cnorm = cnorm + self.scr[loca] * self.scr[locf];
                                        }
                                    } else {
                                        // spectrum part
                                        for j in 1..=self.ngn {
                                            let ig2 = ig2lo + j - 1;
                                            let ig2c = self.icgrp(ig2);
                                            if ig2c != 0 && ig2c <= ngnd {
                                                let jg2c = ngnd - ig2c + 1;
                                                let loc = locchi + jg2c + 1;
                                                a[loc] = a[loc] + self.scr[loca] * self.scr[locf] * cspc[j];
                                                cnorm = cnorm + self.scr[loca] * self.scr[locf] * cspc[j];
                                            }
                                        }
                                    }
                                }
                            }
                        } else {
                            // save constant spectrum
                            for i in 1..=ng2 {
                                cspc[ig2lo + i - 1] = self.scr[l + lz + i - 1];
                            }
                        }
                        lbl = 295;
                    }
                    // add temperature-dependent reactions to matrix.
                    270 => {
                        let mth = c.mth;
                        if mth != 2 && mth != mti && mth != mtc {
                            lbl = 295;
                            continue;
                        }
                        if mth != 2 {
                            jic = 1;
                        }
                        let mut max = 0;
                        let mut min = ngnd;
                        let mut locf = l + lz;
                        if isg > 0 && nz >= iz {
                            locf = l + lz + nl * (iz - 1);
                        }
                        let mut jgc = jg;
                        if isg > 0 && nz >= iz {
                            jgc = jg + ngnd;
                        }
                        if mth == 2 && jg >= nth {
                            lbl = 295;
                            continue;
                        }
                        for i in 2..=ng2 {
                            let ig2 = ig2lo + i - 2;
                            let ig2c = self.icgrp(ig2);
                            if ig2c != 0 && ig2c <= ngnd {
                                let jg2c = ngnd - ig2c + 1;
                                let mut loca = l + lz + (il - 1) + nl * nz * (i - 1);
                                if isg > 0 && nz >= iz {
                                    loca += nl * (iz - 1);
                                }
                                if self.scr[loca] != ZERO {
                                    if mth != 2 {
                                        let loc = locxs + jg - 1 + ngnd * (jg2c - 1);
                                        if (jg2c - jg) >= -1 || self.scr[loca] >= ZERO {
                                            a[loc] = a[loc] + self.scr[loca] * self.scr[l + lz] * self.cflux[jg];
                                            let loc = locxtr + jg - 1;
                                            a[loc] = a[loc] + self.scr[loca] * self.scr[l + lz] * self.cflux[jg];
                                            if jg2c > max {
                                                max = jg2c;
                                            }
                                            if jg2c < min {
                                                min = jg2c;
                                            }
                                        }
                                        if nl != 1 {
                                            if jg2c >= nth {
                                                if igc <= nth1 {
                                                    let index = jg2c;
                                                    let loc = locxtr + jg2c - 1;
                                                    a[loc] = a[loc]
                                                        - self.scr[loca + 1] * self.scr[l + lz + 1] * self.flux[index];
                                                    if ip1opt != 0 {
                                                        let loc = locxs + jg2c - 1 + ngnd * (jg2c - 1);
                                                        a[loc] = a[loc]
                                                            - self.scr[loca + 1]
                                                                * self.scr[l + lz + 1]
                                                                * self.flux[index];
                                                    }
                                                }
                                            }
                                        }
                                    } else {
                                        let loc = jg + ngnd * (jg2c - 1);
                                        self.elas[loc] =
                                            self.elas[loc] + self.scr[loca] * self.scr[locf] * self.cflux[jgc];
                                        if jg2c > max {
                                            max = jg2c;
                                        }
                                        if jg2c < min {
                                            min = jg2c;
                                        }
                                        xtr[jg] = xtr[jg] + self.scr[loca] * self.scr[locf] * self.cflux[jgc];
                                        if nl != 1 {
                                            let mut index = jg2c;
                                            if isg > 0 && nz >= iz {
                                                index += ngnd;
                                            }
                                            xtr[jg2c] = xtr[jg2c]
                                                - self.scr[loca + 1] * self.scr[locf + 1] * self.flux[index];
                                            if ip1opt != 0 {
                                                let loc = jg2c + ngnd * (jg2c - 1);
                                                self.elas[loc] = self.elas[loc]
                                                    - self.scr[loca + 1] * self.scr[locf + 1] * self.flux[index];
                                            }
                                        }
                                    }
                                }
                            }
                        }
                        let lt = l1[jg];
                        if min < lt && mth != 2 {
                            l1[jg] = min;
                        }
                        let lt = l1e[jg];
                        if min < lt && mth == 2 {
                            l1e[jg] = min;
                        }
                        let lt = l2[jg];
                        if max > lt && mth != 2 {
                            l2[jg] = max;
                        }
                        let lt = l2e[jg];
                        if max > lt && mth == 2 {
                            l2e[jg] = max;
                        }
                        lbl = 295;
                    }
                    295 => {
                        lbl = if ig < self.ngn { 160 } else { 200 };
                    }
                    //--store non-temperature dependent matrix .
                    300 => {
                        if *jtemp > 0 {
                            lbl = 205;
                            continue;
                        }
                        if c.mfh != 6 {
                            lbl = 195;
                            continue;
                        }
                        if if6 <= 0 {
                            if nth1 == 0 {
                                nth1 = ngnd - self.nfg - self.nrg;
                            }
                            nth = ngnd - nth1 + 1;
                            if6 = 1;
                        }
                        let mth = c.mth;
                        if (mth >= 18 && mth <= 21) || mth == 38 || mth == 2 || (mth >= 201 && mth <= 255) {
                            lbl = 195;
                            continue;
                        }
                        let mut mult = 1;
                        if mth == 16 {
                            mult = 2;
                        }
                        let locf = l + lz;
                        let mut max = 0;
                        let mut min = ngnd;
                        for i in 2..=ng2 {
                            let ig2 = ig2lo + i - 2;
                            let ig2c = self.icgrp(ig2);
                            if ig2c != 0 && ig2c <= ngnd {
                                let jg2c = ngnd - ig2c + 1;
                                let loca = l + lz + (il - 1) + nl * nz * (i - 1);
                                let loc = locxs + jg - 1 + ngnd * (jg2c - 1);
                                a[loc] = a[loc] + self.scr[loca] * self.scr[locf] * self.cflux[jg];
                                let loc = locxtr + jg - 1;
                                a[loc] = a[loc] + (self.scr[loca] * mult as f64) * self.scr[locf] * self.cflux[jg];
                                if jg2c > max {
                                    max = jg2c;
                                }
                                if jg2c < min {
                                    min = jg2c;
                                }
                                if nl != 1 {
                                    let loc = locxtr + jg2c - 1;
                                    let index = jg2c;
                                    a[loc] = a[loc] - self.scr[loca + 1] * self.scr[locf + 1] * self.flux[index];
                                    if ip1opt != 0 {
                                        let loc = locxs + jg2c - 1 + ngnd * (jg2c - 1);
                                        a[loc] = a[loc] - self.scr[loca + 1] * self.scr[locf + 1] * self.flux[index];
                                    }
                                }
                            }
                        }
                        let lt = l1[jg];
                        if min < lt {
                            l1[jg] = min;
                        }
                        let lt = l2[jg];
                        if max > lt {
                            l2[jg] = max;
                        }
                        lbl = if ig < self.ngn { 160 } else { 200 };
                    }
                    //--write non-temperature dependent matrix on scratch tape.
                    320 => {
                        *jtemp -= 1;
                        lbl = 330;
                    }
                    330 => {
                        if *jtemp > 0 {
                            lbl = 360;
                            continue;
                        }
                        c.repoz(-kscr);
                        let g = ngnd as usize;
                        wr_ints(c, kscr, &l1.0[..g]);
                        wr_ints(c, kscr, &l2.0[..g]);
                        wr_reals(c, kscr, &a.0[(locxtr - 1) as usize..(locxtr - 1) as usize + g]);
                        wr_reals(c, kscr, &a.0[(locxs - 1) as usize..lim as usize]);
                        c.skiprz(ng, -3);
                        c.findf(matd, 1, 451, ng);
                        if nth1 > 0 && nth > 0 {
                            lbl = 140;
                            continue;
                        }
                        if nth1 == 0 {
                            nth1 = ngnd - self.nfg - self.nrg;
                        }
                        nth = ngnd - nth1 + 1;
                        lbl = 140;
                    }
                    //--desired data is loaded.
                    350 => {
                        matd = -1;
                        if c.math != 0 {
                            c.skiprz(ng, -3);
                        }
                        lbl = 500;
                    }
                    360 => {
                        if jic > 0 {
                            lbl = 365;
                            continue;
                        }
                        *jtemp -= 1;
                        let strng = format!("use only {} temperatures for mat {}", ffmt::i(*jtemp, 2), ffmt::i(matd, 4));
                        c.mess("sfile4", &strng, "mti missing from higher temperatures");
                        lbl = 350;
                    }
                    365 => {
                        if *jtemp == 1 {
                            // calculate nu at lowest temperature
                            for i in 1..=ngnd {
                                if a[locsf0 + i - 1] == ZERO {
                                    snu[i] = 0.0;
                                }
                                if a[locsf0 + i - 1] != ZERO {
                                    snu[i] = a[locnus + i - 1] / a[locsf0 + i - 1];
                                }
                            }
                        } else {
                            // calculate nu*sigma f at higher temperatures
                            for i in 1..=ngnd {
                                a[i - 1 + locnus] = snu[i] * a[i - 1 + locsf0];
                            }
                        }
                        lbl = 400;
                    }
                    // calculate total absorption
                    400 => {
                        for i in 1..=ngnd {
                            a[i - 1 + locab0] = a[i - 1 + locab0] + abs1[i] + abs2[i] - sn2n[i];
                        }
                        // correct scattering matrix
                        let inn = nth - 1;
                        for i in 1..=inn {
                            if l1e[i] < l1[i] {
                                l1[i] = l1e[i];
                            }
                            if l2e[i] > l2[i] {
                                l2[i] = l2e[i];
                            }
                        }
                        for i in 1..=ngndsq {
                            a[i - 1 + locxs] = a[i - 1 + locxs] + self.elas[i];
                        }
                        // correct transport
                        for i in 1..=ngnd {
                            a[i - 1 + locxtr] = a[i - 1 + locxtr] + xtr[i];
                        }
                        lbl = 500;
                    }
                    500 => {
                        self.sf4out(c, matd, nmat, ntemp, jfiss, nina, locab0, locsf0, locnus, locxtr, locxs, &l1, &l2, &a);
                        if nina != 3 {
                            lbl = 505;
                            continue;
                        }
                        matd = -1;
                        nina = 0;
                        lbl = 500;
                    }
                    505 => {
                        if matd < 0 {
                            lbl = 540;
                            continue;
                        }
                        mtabs1 = 0;
                        mtabs2 = 0;
                        mttot = 0;
                        mtel = 0;
                        mtex = 0;
                        i318 = 0;
                        i618 = 0;
                        c.repoz(-kscr);
                        let g = ngnd;
                        let v = rd_ints(c, kscr, g);
                        for i in 1..=g {
                            l1[i] = v[(i - 1) as usize];
                        }
                        let v = rd_ints(c, kscr, g);
                        for i in 1..=g {
                            l2[i] = v[(i - 1) as usize];
                        }
                        let v = rd_reals(c, kscr, g);
                        for i in 1..=g {
                            a[locxtr + i - 1] = v[(i - 1) as usize];
                        }
                        let v = rd_reals(c, kscr, lim - locxs + 1);
                        for i in locxs..=lim {
                            a[i] = v[(i - locxs) as usize];
                        }
                        self.fiss[nmat] = jfiss as f64;
                        if isof != 0 && *jtemp <= 1 {
                            // normalize chi
                            if cnorm != ZERO {
                                cnorm = 1.0 / cnorm;
                            }
                            for i in 1..=self.nfsg {
                                let loc = locchi + i - 1;
                                if cnorm != ZERO {
                                    a[loc] = a[loc] * cnorm;
                                }
                                let is = isof as usize;
                                if self.nfiss[is] != 0 && a[loc] != ZERO {
                                    self.nfiss[is] += 1;
                                }
                                if self.nfiss[is] == 0 && a[loc] != ZERO {
                                    self.nfiss[is] = 1;
                                }
                                if isof == 1 {
                                    self.uff[i] = a[loc];
                                }
                                if isof == 2 {
                                    self.puff[i] = a[loc];
                                }
                            }
                        }
                        for i in 1..=ngnd5 {
                            a[i - 1 + iloc] = 0.0;
                        }
                        lbl = 200;
                    }
                    _ => {
                        // 540
                        if ntemp == 0 {
                            self.ntmp[nmat] = *jtemp;
                        }
                        if nmat == self.newmat {
                            break;
                        }
                        lbl = 120;
                    }
                }
            }

            //--sfile4 is finished.
            dealloc(&mut self.al_tempt, "tempt");
            dealloc(&mut self.al_flux, "flux");
            // re-invert coarse group p0 flux for infinite dilution
            for i in 1..=ngnd {
                if self.cflux[i] != ZERO {
                    self.cflux[i] = 1.0 / self.cflux[i];
                }
            }
            let _ = (mttot, tot, ident, i318, i618, mtel, mtex);
        }
        let mut r = RecW::new();
        r.r(EFLAG);
        c.io.write_rec(self.nscr, r);
    }

    /// Write subfile 4 on nscr.
    fn sf4out(
        &mut self,
        c: &mut Ctx,
        matd: i32,
        nmat: i32,
        ntemp: i32,
        ifiss: i32,
        nina: i32,
        locab0: i32,
        locsf0: i32,
        locnus: i32,
        locxtr: i32,
        locxs: i32,
        l1: &IV,
        l2: &IV,
        a: &RV,
    ) {
        let ncol = 5;
        let ncol1 = ncol - 1;
        let izero = 0;
        let ngnd = self.ngnd;
        let (mscr, nscr) = (self.mscr, self.nscr);
        let g = ngnd as usize;
        let sl = |a: &RV, i0: i32, n: i32| -> Vec<f64> { a.0[(i0 - 1) as usize..(i0 - 1 + n) as usize].to_vec() };
        let mut lbl;
        if nina == 2 {
            lbl = 120;
        } else if matd == self.s4_matl {
            lbl = 160;
        } else if self.s4_matl == 0 {
            lbl = 150;
        } else {
            lbl = 0;
        }
        let mut ident = 0;
        let mut jtemp = 0;
        let mut ktemp = 0;
        let (mut nc, mut ic) = (0, 0);
        loop {
            match lbl {
                0 => {
                    // material has been completed.  copy mscr to nscr.
                    let flag = end_flag();
                    wr_reals(c, mscr, &[flag]);
                    ident = self.idnt[nmat];
                    if self.iprint == 1 && nmat > 1 {
                        c.io.so("\n\n");
                    }
                    if self.iprint == 1 {
                        c.io.so(&ffmt::i(ident, 1));
                    }
                    c.repoz(-mscr);
                    jtemp = 1;
                    let tv: Vec<f64> = (1..=ntemp).map(|i| self.tempt[i]).collect();
                    wr_reals(c, nscr, &tv);
                    if self.iprint != 0 {
                        c.io.so("\n temperatures\n ------------");
                        c.io.so(&rowsj(&tv, 6, |x| pe(x, 12, 4)));
                    }
                    lbl = 100;
                }
                100 => {
                    ktemp = nint(self.tempt[jtemp]);
                    nc = 4;
                    ic = 1;
                    if self.iprint == 1 {
                        c.io.so(&format!(
                            "\n absorption for nuclide {} at {} deg k\n {}",
                            ffmt::i(ident, 6),
                            ffmt::i(ktemp, 4),
                            "-".repeat(43)
                        ));
                    }
                    lbl = 105;
                }
                105 => {
                    let b = rd_reals(c, mscr, ngnd);
                    wr_reals(c, nscr, &b);
                    let ifissl = self.s4_ifissl;
                    if self.iprint == 1 && ic == 1 {
                        c.io.so(&rowsj(&b, 6, |x| pe(x, 12, 4)));
                    }
                    if self.iprint == 1 && ifissl == 1 && ic > 1 && ic < nc {
                        c.io.so(&rowsj(&b, 6, |x| pe(x, 12, 4)));
                    }
                    if self.iprint == 1 && ic == nc {
                        c.io.so(&rowsj(&b, 6, |x| pe(x, 12, 4)));
                    }
                    ic += 1;
                    if self.s4_ninal == 1 {
                        lbl = 135;
                        continue;
                    }
                    if ic > nc {
                        lbl = 125;
                        continue;
                    }
                    if !(self.iprint == 0 || (ifissl == 0 && ic > 1 && ic < nc) || ic == 1) {
                        if ic == 2 {
                            c.io.so(&format!(
                                "\n fission for nuclide {} at {} deg k\n {}",
                                ffmt::i(ident, 6),
                                ffmt::i(ktemp, 4),
                                "-".repeat(40)
                            ));
                        } else if ic == 3 {
                            c.io.so(&format!(
                                "\n nu*fission for nuclide {} at {} deg k\n {}",
                                ffmt::i(ident, 6),
                                ffmt::i(ktemp, 4),
                                "-".repeat(43)
                            ));
                        } else if ic == 4 {
                            c.io.so(&format!(
                                "\n transport corrected total for nuclide {} at {} deg k\n {}",
                                ffmt::i(ident, 6),
                                ffmt::i(ktemp, 4),
                                "-".repeat(58)
                            ));
                        }
                    }
                    // 115
                    lbl = 105;
                }
                120 => {
                    ident = self.idnt[nmat];
                    ktemp = 0;
                    wr_reals(c, nscr, &[self.tempt[1]]);
                    wr_reals(c, nscr, &sl(a, locab0, ngnd));
                    if self.iprint == 1 {
                        if nmat > 1 {
                            c.io.so("\n\n");
                        }
                        c.io.so(&format!("\n nuclide {}\n -------", ffmt::i(ident, 6)));
                        c.io.so("\n temperatures\n ------------");
                        c.io.so(&pe(self.tempt[1], 12, 4));
                        c.io.so(&format!(
                            "\n absorption for nuclide {} at {} deg k\n {}",
                            ffmt::i(ident, 6),
                            ffmt::i(ktemp, 4),
                            "-".repeat(43)
                        ));
                        c.io.so(&rowsj(&sl(a, locab0, ngnd), 6, |x| pe(x, 12, 4)));
                    }
                    self.s4_ifissl = 0;
                    self.s4_ninal = 1;
                    lbl = 200;
                }
                125 => {
                    if self.iprint == 1 {
                        let mut s = format!(
                            "\n p0 scattering matrix for nuclide {} at {} deg k\n {}\n\n ia l1 l2    l1+0",
                            ffmt::i(ident, 6),
                            ffmt::i(ktemp, 4),
                            "-".repeat(53)
                        );
                        for i in 1..=ncol1 {
                            s.push_str(&format!("{}+{}", " ".repeat(10), ffmt::i(i, 1)));
                        }
                        s.push_str(&format!("\n {}----", "----------".repeat(7)));
                        c.io.so(&s);
                    }
                    let mut stop135 = false;
                    for _ig in 1..=ngnd {
                        let v = rd_ints(c, mscr, 3);
                        let (ia, lone, ltwo) = (v[0], v[1], v[2]);
                        if ia == 0 {
                            stop135 = true;
                            break;
                        }
                        wr_ints(c, nscr, &[ia, lone, ltwo]);
                        let nb = ltwo - lone + 1;
                        let b = rd_reals(c, mscr, nb);
                        wr_reals(c, nscr, &b);
                        let i1 = 1;
                        let mut lim = nb;
                        if lim - i1 + 1 > ncol {
                            lim = ncol;
                        }
                        let mut s = format!("{}{}{}", ffmt::i(ia, 3), ffmt::i(lone, 3), ffmt::i(ltwo, 3));
                        if lim >= i1 {
                            s.push_str("     ");
                        }
                        for i in i1..=lim {
                            s.push_str(&pe(b[(i - 1) as usize], 12, 4));
                        }
                        c.io.so(&s);
                        if lim != nb && cfg!(feature = "powr-cards") {
                            // tools/njoy-powr-cards.patch: '(14x,1p,5e12.4)'
                            let mut lim = lim;
                            while lim != nb {
                                let i1 = lim + 1;
                                lim = nb;
                                if lim - i1 + 1 > ncol {
                                    lim = i1 + ncol - 1;
                                }
                                let mut s = " ".repeat(14);
                                for i in i1..=lim {
                                    s.push_str(&pe(b[(i - 1) as usize], 12, 4));
                                }
                                c.io.so(&s);
                            }
                        } else if lim != nb {
                            // '(9x,i4,2x,1p,5e12.4)' with a real for the i4
                            // item: a run-time error of the Fortran library
                            c.io.flush_all();
                            fatal_runtime(
                                "Fortran runtime error: Expected INTEGER for item 1 in formatted transfer, got REAL",
                            );
                        }
                    }
                    if !stop135 {
                        // 130
                        let _ = rd_ints(c, mscr, 1);
                    }
                    lbl = 135;
                }
                135 => {
                    jtemp += 1;
                    lbl = if jtemp <= ntemp { 100 } else { 150 };
                }
                150 => {
                    if matd < 0 {
                        lbl = 200;
                        continue;
                    }
                    c.repoz(-mscr);
                    self.s4_matl = matd;
                    self.s4_ifissl = ifiss;
                    self.s4_ninal = nina;
                    lbl = 160;
                }
                // write data on mscr.
                160 => {
                    wr_reals(c, mscr, &sl(a, locab0, ngnd));
                    if nina <= 0 || nina == 3 {
                        wr_reals(c, mscr, &sl(a, locsf0, ngnd));
                        wr_reals(c, mscr, &sl(a, locnus, ngnd));
                        wr_reals(c, mscr, &sl(a, locxtr, ngnd));
                        for ia in 1..=ngnd {
                            let lone = l1[ia];
                            let ltwo = l2[ia];
                            if lone != ngnd || ltwo != 1 {
                                let mut b = Vec::new();
                                for i in lone..=ltwo {
                                    let loc = locxs + ia - 1 + ngnd * (i - 1);
                                    b.push(a[loc]);
                                }
                                wr_ints(c, mscr, &[ia, lone, ltwo]);
                                wr_reals(c, mscr, &b);
                            } else {
                                wr_ints(c, mscr, &[ia, ia, ia]);
                                wr_reals(c, mscr, &[ZERO]);
                            }
                        }
                        wr_ints(c, mscr, &[izero]);
                    }
                    lbl = 200;
                }
                _ => {
                    //--sf4out is finished.
                    if matd < 0 {
                        self.s4_matl = 0;
                    }
                    let _ = g;
                    return;
                }
            }
        }
    }

    /// Write subfile 5 burnup data.
    pub(super) fn sf5out(&mut self, c: &mut Ctx) {
        c.io.so("\n ***subfile 5***burnup data");
        let nfis1 = self.nfis + 1;
        let nscr = self.nscr;
        let mut r = RecW::new();
        r.i(self.ntis);
        for i in 1..=self.ntis {
            r.i(self.idat[i]);
        }
        c.io.write_rec(nscr, r);
        let mut r = RecW::new();
        r.i(self.nfis);
        for i in 1..=self.nfis {
            r.i(self.idbt[i]);
        }
        c.io.write_rec(nscr, r);
        for i in 1..=self.ntis {
            let index = 1 + (i - 1) * nfis1;
            let v: Vec<f64> = (1..=nfis1).map(|j| self.yld[index + j - 1]).collect();
            wr_reals(c, nscr, &v);
        }
        wr_reals(c, nscr, &[EFLAG]);
        dealloc(&mut self.al_idat, "idat");
        dealloc(&mut self.al_idbt, "idbt");
        dealloc(&mut self.al_yld, "yld");
    }

    /// Process P1 scattering matrix.
    pub(super) fn sfile6(&mut self, c: &mut Ctx, ntemp: &mut i32) {
        let lz = 6;
        let ng = self.ngendf;
        let ngnd = self.ngnd;
        let (mut nb, mut nw) = (0i32, 0i32);

        //--allocate storage.
        if self.if4 != 0 {
            c.io.so("\n\n");
            c.io.so("\n ***subfile 6***p1 matrices");
            if self.nnina != 0 {
                let nther = ngnd - self.nfg - self.nrg;
                let newmat = self.newmat;
                let mut ids = IV::new(newmat.max(0) as usize + SLACK);
                let ngndsq = ngnd * ngnd;
                let mut s = RV::new(ngndsq as usize + SLACK);
                nw = newmat * 6;
                let mut tempt = RV::new(nw.max(0) as usize + SLACK);
                let g = ngnd as usize + SLACK;
                let mut l1 = IV::new(g);
                let mut l2 = IV::new(g);
                let mut l1e = IV::new(g);
                let mut l2e = IV::new(g);
                let mut ntp = [0i32; 7 + 16];
                let mut tempr = [0.0f64; 8 + 16];

                //--search for desired materials.
                c.repoz(ng);
                c.tpidio(ng, 0, 0, &mut self.scr.at(1), &mut nb, &mut nw);
                let il = 2;
                let mut np1 = 0;
                let mut nmat = 0;
                let (mut matd, mut mti, mut mtc, mut ntemps, mut mtelas, mut itd) = (0, 0, 0, 0, 0, 0);
                let (mut nth1, mut nth, mut ip1, mut iza) = (0, 0, 0, 0);
                let (mut nl, mut nz, mut l, mut ng2, mut ig2lo, mut ig, mut igc, mut jg) = (0, 0, 1, 0, 0, 0, 0, 0);
                let mut temp = 0.0f64;
                let mut lbl = 100;
                loop {
                    match lbl {
                        100 => {
                            nmat += 1;
                            let ip1opt = self.ip1tab[nmat];
                            if ip1opt > 0 {
                                lbl = 320;
                                continue;
                            }
                            let nina = self.ninat[nmat];
                            if nina > 0 {
                                lbl = 320;
                                continue;
                            }
                            matd = self.mat[nmat];
                            mti = self.mtit[nmat];
                            mtc = self.mtct[nmat];
                            ntemps = self.ntmp[nmat];
                            mtelas = 0;
                            itd = 0;
                            for i in 1..=ngndsq {
                                s[i] = 0.0;
                            }
                            for i in 1..=ngnd {
                                l1[i] = ngnd;
                                l2[i] = 1;
                            }
                            *ntemp = -1;
                            nth1 = 0;
                            nth = nther;
                            ip1 = 0;
                            c.math = 1;
                            while c.math != matd {
                                c.contio(ng, 0, 0, &mut self.scr.at(1), &mut nb, &mut nw);
                                if c.math != matd {
                                    if c.math < 0 {
                                        c.repoz(ng);
                                        c.tpidio(ng, 0, 0, &mut self.scr.at(1), &mut nb, &mut nw);
                                    } else if c.math > 0 {
                                        c.tomend(ng, 0, 0, &mut self.scr.at(1));
                                    }
                                }
                            }
                            lbl = 130;
                        }
                        //--search for required reactions.
                        125 => {
                            c.contio(ng, 0, 0, &mut self.scr.at(1), &mut nb, &mut nw);
                            if c.math == 0 {
                                lbl = 200;
                                continue;
                            }
                            if c.math < 0 || c.math != matd {
                                lbl = 300;
                                continue;
                            }
                            if c.mfh == 0 || c.mth == 0 {
                                lbl = 125;
                                continue;
                            }
                            if c.mfh == 1 && *ntemp == ntemps {
                                lbl = 300;
                                continue;
                            }
                            lbl = if c.mfh == 1 { 130 } else { 135 };
                        }
                        130 => {
                            iza = nint(c.c1h);
                            lbl = 135;
                        }
                        135 => {
                            nl = c.l1h;
                            if nl < il && c.mfh != 1 {
                                lbl = 155;
                                continue;
                            }
                            nz = c.l2h;
                            lbl = 140;
                        }
                        140 => {
                            l = 1;
                            c.listio(ng, 0, 0, &mut self.scr.at(l), &mut nb, &mut nw);
                            if c.mfh == 1 {
                                temp = c.c1h;
                                *ntemp += 1;
                            }
                            if c.mfh == 3 {
                                lbl = 160;
                                continue;
                            }
                            if c.mfh != 6 {
                                lbl = 155;
                                continue;
                            }
                            if c.mth == 2 && mtelas <= 0 {
                                mtelas = 1;
                                for i in 1..=ngnd {
                                    l1e[i] = ngnd;
                                    l2e[i] = 1;
                                }
                                for i in 1..=ngndsq {
                                    self.elas[i] = 0.0;
                                }
                            }
                            // 148
                            ng2 = c.l1h;
                            ig2lo = c.l2h;
                            ig = c.n2h;
                            igc = self.icgrp(ig);
                            while nb != 0 {
                                l += nw;
                                c.moreio(ng, 0, 0, &mut self.scr.at(l), &mut nb, &mut nw);
                                if l + nw - 1 > self.nwscr {
                                    c.error("sfile6", "storage exceeded.", " ");
                                }
                            }
                            lbl = 155;
                        }
                        155 => {
                            c.tosend(ng, 0, 0, &mut self.scr.at(1));
                            lbl = 125;
                        }
                        160 => {
                            if igc == 0 || igc > ngnd {
                                lbl = 225;
                                continue;
                            }
                            jg = ngnd - igc + 1;
                            l = 1;

                            //--store non-temperature dependent scattering matrix
                            if *ntemp > 0 {
                                lbl = 210;
                                continue;
                            }
                            let mth = c.mth;
                            if c.mfh == 3 {
                                // 230: coarse group flux
                                if mth != 1 {
                                    lbl = 155;
                                    continue;
                                }
                                self.cflux[jg] = self.cflux[jg] + self.scr[l + lz];
                                lbl = if ig < self.ngn { 140 } else { 125 };
                                continue;
                            }
                            if (mth >= 18 && mth <= 21) || mth == 38 || mth == 2 {
                                lbl = 155;
                                continue;
                            }
                            if mth == mti || mth == mtc {
                                // 180
                                if igc > nth1 && igc < nth {
                                    nth1 = igc;
                                }
                                lbl = if igc < nth { 140 } else { 155 };
                                continue;
                            }
                            if mth >= 201 && mth <= 250 {
                                lbl = 155;
                                continue;
                            }
                            for i in 2..=ng2 {
                                let ig2 = ig2lo + i - 3;
                                let ig2c = self.icgrp(ig2);
                                if ig2c != 0 && ig2c <= ngnd {
                                    let jg2c = ngnd - ig2c + 1;
                                    let loca = l + lz + (il - 1) + nl * nz * (i - 1);
                                    let loc = 1 + jg - 1 + ngnd * (jg2c - 1);
                                    s[loc] = s[loc] + self.scr[loca];
                                }
                            }
                            ip1 = 1;
                            let lt = l1[jg];
                            let l1t = self.ngn - ig2lo - ng2 + 3;
                            let l1tc = self.icgrp(l1t);
                            if l1tc < lt && l1tc != 0 && l1tc <= ngnd {
                                l1[jg] = l1tc;
                            }
                            let lt = l2[jg];
                            let l2t = l1t + ng2 - 2;
                            let l2tc = self.icgrp(l2t);
                            if l2tc > lt && l2tc != l1t && l2tc > 0 && l2tc <= ngnd {
                                l2[jg] = l2tc;
                            }
                            lbl = if ig < self.ngn { 140 } else { 125 };
                        }
                        //--write non-temperature dependent matrix on scratch tape.
                        200 => {
                            if *ntemp > 0 {
                                lbl = 250;
                                continue;
                            }
                            let kscr = self.kscr;
                            c.repoz(-kscr);
                            let gg = ngnd as usize;
                            wr_ints(c, kscr, &l1.0[..gg]);
                            wr_ints(c, kscr, &l2.0[..gg]);
                            wr_reals(c, kscr, &s.0[..ngndsq as usize]);
                            c.repoz(-kscr);
                            c.skiprz(ng, -3);
                            c.findf(matd, 1, 451, ng);
                            nth = nth1;
                            lbl = 125;
                        }
                        //--add temperature-dependent reaction to matrix.
                        210 => {
                            let mth = c.mth;
                            if !(mth == 2 || mth == mti || mth == mtc) {
                                lbl = 155;
                                continue;
                            }
                            // 215
                            if mth == 2 && igc <= nth {
                                lbl = 140;
                                continue;
                            }
                            if mth == mti && igc > nth {
                                lbl = 155;
                                continue;
                            }
                            if mth == mtc && igc > nth && mtc != 0 {
                                lbl = 155;
                                continue;
                            }
                            itd = 1;
                            for i in 2..=ng2 {
                                let ig2 = ig2lo + i - 2;
                                let ig2c = self.icgrp(ig2);
                                if ig2c != 0 && ig2c <= ngnd {
                                    let jg2c = ngnd - ig2c + 1;
                                    let loca = l + lz + (il - 1) + nl * nz * (i - 1);
                                    let loc = jg + ngnd * (jg2c - 1);
                                    s[loc] = s[loc] + self.scr[loca];
                                    if mth == 2 {
                                        let loc = jg + ngnd * (jg2c - 1);
                                        self.elas[loc] = self.elas[loc] + self.scr[loca];
                                    }
                                }
                            }
                            ip1 = 1;
                            let lt = l1[jg];
                            let l1t = self.ngn - ig2lo - ng2 + 3;
                            let l1tc = self.icgrp(l1t);
                            if l1tc < lt && l1tc != 0 {
                                l1[jg] = l1tc;
                            }
                            if mth == 2 {
                                l1e[jg] = l1tc;
                            }
                            let lt = l2[jg];
                            let l2t = l1t + ng2 - 2;
                            let l2tc = self.icgrp(l2t);
                            if l2tc > lt && l2tc != l1t && l2tc > 0 && l2tc <= ngnd {
                                l2[jg] = l2tc;
                            }
                            if mth == 2 && l2tc != l1t {
                                l2e[jg] = l2tc;
                            }
                            lbl = 225;
                        }
                        225 => {
                            lbl = if ig < self.ngn { 140 } else { 125 };
                        }
                        //--desired data is loaded
                        250 => {
                            if ip1 == 0 {
                                let strng = format!("no required p1 matrices found for mat {}", ffmt::i(matd, 4));
                                c.error("sfile6", &strng, " ");
                            }
                            if mti > 0 && itd == 0 {
                                let strng = format!("no temperature-dependent reactions for mat {}", ffmt::i(matd, 4));
                                c.error("sfile6", &strng, "check mti");
                            }
                            if mtelas <= 0 {
                                for i in 1..=ngnd {
                                    if l1e[i] < l1[i] {
                                        l1[i] = l1e[i];
                                    }
                                    if l2e[1] > l2[i] {
                                        l2[i] = l2e[i];
                                    }
                                }
                                for i in 1..=ngndsq {
                                    s[i] = s[i] + self.elas[i];
                                }
                            }
                            let math = c.math;
                            self.sf6out(c, math, np1, &ntp, &ids, &tempt, &l1, &l2, &s);
                            if *ntemp >= 1 {
                                tempr[*ntemp as usize] = temp;
                            }
                            let kscr = self.kscr;
                            c.repoz(-kscr);
                            let gg = ngnd;
                            let v = rd_ints(c, kscr, gg);
                            for i in 1..=gg {
                                l1[i] = v[(i - 1) as usize];
                            }
                            let v = rd_ints(c, kscr, gg);
                            for i in 1..=gg {
                                l2[i] = v[(i - 1) as usize];
                            }
                            let v = rd_reals(c, kscr, ngndsq);
                            for i in 1..=ngndsq {
                                s[i] = v[(i - 1) as usize];
                            }
                            mtelas = 0;
                            lbl = 125;
                        }
                        //--this material is finished.
                        300 => {
                            c.skiprz(ng, -1);
                            np1 += 1;
                            ids[np1] = iza;
                            for i in 1..=*ntemp {
                                tempt[np1 + newmat * (i - 1)] = tempr[i as usize];
                            }
                            ntp[np1 as usize] = *ntemp;
                            lbl = 320;
                        }
                        _ => {
                            // 320
                            if nmat < newmat {
                                lbl = 100;
                                continue;
                            }
                            break;
                        }
                    }
                }

                //--sfile6 is finished.
                self.sf6out(c, -1, np1, &ntp, &ids, &tempt, &l1, &l2, &s);
                let _ = (jg, ip1);
            }
            // 400
            dealloc(&mut self.al_mtit, "mtit");
            dealloc(&mut self.al_mtct, "mtct");
            c.closz(self.mscr);
            c.closz(self.kscr);
        }
        // 410
        let nscr = self.nscr;
        wr_reals(c, nscr, &[EFLAG]);
        wr_reals(c, nscr, &[end_flag()]);
    }

    /// Write subfile 6.
    fn sf6out(&mut self, c: &mut Ctx, matd: i32, np1: i32, ntp: &[i32], ids: &IV, temp: &RV, l1: &IV, l2: &IV, s: &RV) {
        let ncol = 5;
        let ngnd = self.ngnd;
        let (mscr, nscr) = (self.mscr, self.nscr);

        //--store matrix on mscr.
        let ncol1 = ncol - 1;
        if self.s6_jtem == 0 {
            c.repoz(-mscr);
        }
        self.s6_jtem += 1;
        if matd > 0 {
            for ia in 1..=ngnd {
                let lone = l1[ia];
                let ltwo = l2[ia];
                if ltwo < lone {
                    let strng = format!(
                        "for group {} l1={} is lt l2={} at temp no. {}",
                        ffmt::i(ia, 3),
                        ffmt::i(lone, 2),
                        ffmt::i(ltwo, 2),
                        ffmt::i(self.s6_jtem, 2)
                    );
                    c.error("sf6out", &strng, " ");
                }
                let mut b = Vec::new();
                for i in lone..=ltwo {
                    let loc = ia + ngnd * (i - 1);
                    b.push(s[loc]);
                }
                wr_ints(c, mscr, &[ia, lone, ltwo]);
                wr_reals(c, mscr, &b);
            }

        //--write subfile 6 on nscr.
        } else if np1 != 0 {
            let k: Vec<i32> = (1..=np1).map(|i| ids[i]).collect();
            let mut v = vec![np1];
            v.extend_from_slice(&k);
            wr_ints(c, nscr, &v);
            if self.iprint != 0 {
                let mut s = format!("\n np1 {}\n ---\n\n nuclides\n --------", ffmt::i(np1, 2));
                for ch in k.chunks(11) {
                    s.push('\n');
                    for &x in ch {
                        s.push_str(&ffmt::i(x, 7));
                    }
                }
                c.io.so(&s);
                c.io.so("\n ntemp    temperatures\n -----    ------------");
            }
            for i in 1..=np1 {
                let ntemps = ntp[i as usize];
                let b: Vec<f64> = (1..=ntemps).map(|j| temp[i + self.newmat * (j - 1)]).collect();
                let mut r = RecW::new();
                r.i(ntemps).rs(&b);
                c.io.write_rec(nscr, r);
                if self.iprint == 1 {
                    // '(i5,5x,1p,5e12.4/(10x,5e12.4))'
                    let mut s = ffmt::i(ntemps, 5);
                    if !b.is_empty() {
                        // (a trailing 5x with no value after it writes nothing)
                        s.push_str("     ");
                    }
                    let first = b.len().min(5);
                    for &x in &b[..first] {
                        s.push_str(&pe(x, 12, 4));
                    }
                    for ch in b[first..].chunks(5) {
                        s.push('\n');
                        s.push_str(&" ".repeat(10));
                        for &x in ch {
                            s.push_str(&pe(x, 12, 4));
                        }
                    }
                    c.io.so(&s);
                }
            }
            c.repoz(-mscr);
            for ip1 in 1..=np1 {
                let ident = ids[ip1];
                let ntemps = ntp[ip1 as usize];
                for it in 1..=ntemps {
                    let ktemp = nint(temp[ip1 + self.newmat * (it - 1)]);
                    if self.iprint == 1 {
                        let mut s = format!(
                            "\n p1 scattering matrix for nuclide {} ({}) at {} deg k\n {}--------\n ia l1 l2    l1+0",
                            ffmt::i(ip1, 2),
                            ffmt::i(ident, 6),
                            ffmt::i(ktemp, 4),
                            "----------".repeat(5)
                        );
                        for i in 1..=ncol1 {
                            s.push_str(&format!("{}+{}", " ".repeat(10), ffmt::i(i, 1)));
                        }
                        s.push_str(&format!("\n {}----", "----------".repeat(7)));
                        c.io.so(&s);
                    }
                    for _ig in 1..=ngnd {
                        let v = rd_ints(c, mscr, 3);
                        let (ia, lone, ltwo) = (v[0], v[1], v[2]);
                        wr_ints(c, nscr, &[ia, lone, ltwo]);
                        let nb = ltwo - lone + 1;
                        let b = rd_reals(c, mscr, nb);
                        wr_reals(c, nscr, &b);
                        if self.iprint != 0 {
                            let mut i1 = 1;
                            let mut lim = nb;
                            if lim > ncol {
                                lim = ncol;
                            }
                            let mut s = format!("{}{}{}", ffmt::i(ia, 3), ffmt::i(lone, 3), ffmt::i(ltwo, 3));
                        if lim >= i1 {
                            s.push_str("     ");
                        }
                            for i in i1..=lim {
                                s.push_str(&pe(b[(i - 1) as usize], 12, 4));
                            }
                            c.io.so(&s);
                            while lim != nb {
                                i1 = lim + 1;
                                lim = nb;
                                if lim - i1 + 1 > ncol {
                                    lim = i1 + ncol - 1;
                                }
                                let mut s = " ".repeat(14);
                                for i in i1..=lim {
                                    s.push_str(&pe(b[(i - 1) as usize], 12, 4));
                                }
                                c.io.so(&s);
                            }
                        }
                    }
                }
            }
        }
    }

    /// Write the input for CLIB on nout.
    pub(super) fn cpmout(&mut self, c: &mut Ctx, nlib: i32, idat: i32) {
        let nold = 65;
        const IDNTO: [i32; 80] = [
            1001, 1002, 5000, 5010, 6000, 8000, 8001, 13000, 14000, 24000, 25000, 26000, 28000, 29063, 47000, 48000,
            49000, 64154, 64155, 64156, 64157, 64158, 66164, 71176, 36083, 45103, 45105, 47109, 54131, 54135, 55133,
            55134, 55135, 60143, 60145, 61147, 61148, 61248, 62147, 62149, 62150, 62151, 62152, 63153, 63154, 63155,
            401, 402, 92234, 92235, 92236, 92238, 93237, 94238, 94239, 94240, 94241, 94242, 95241, 95242, 95243,
            96242, 96244, 1, 302, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        ];
        let nreso = 4;
        const IDRESO: [i32; 20] = [92235, 92236, 92238, 94329, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
        let nf2 = 1;
        const NFR2: [i32; 1] = [92238];
        let np1o = 4;
        const IDP1O: [i32; 10] = [1001, 8000, 8001, 1002, 0, 0, 0, 0, 0, 0];
        let fmt = "(6e12.0)";
        let izero = 0;
        let ione = 1;
        let ithree = 3;
        let eps = 1.0e-5f64;
        let nout = self.nout;
        let nscr = self.nscr;
        let ngnd = self.ngnd;
        let newmat = self.newmat;
        let idnto = |ip: i32| -> i32 { if ip >= 1 && ip <= 80 { IDNTO[(ip - 1) as usize] } else { 0 } };
        let wl = |c: &mut Ctx, s: &str| c.io.write_line(nout, s);
        let wrows = |c: &mut Ctx, v: &[f64]| {
            for s in rows(v, 6, |x| pe(x, 12, 5)) {
                c.io.write_line(nout, &s);
            }
        };
        let flags = |v: &[i32]| -> String { v.iter().map(|&x| ffmt::i(x, 1)).collect() };

        //--initialize.
        let time = c.timer();
        c.io.so(&format!("\n ***output***{}{}s", " ".repeat(56), ffmt::f(time, 8, 1)));

        //--write input for subfile 1.
        wl(c, &format!("pri   {}{}", ffmt::i(izero, 6), ffmt::i(izero, 6)));
        for _ in 0..4 {
            wl(c, fmt);
        }
        if self.mode <= 1 {
            wl(c, &format!("{}{}", ffmt::i(nold, 6), ffmt::i(izero, 6)));
        }
        if self.mode == 2 {
            wl(c, &format!("{}{}", ffmt::i(izero, 6), ffmt::i(ione, 6)));
        }
        if self.mode <= 1 {
            wl(c, &ffmt::i(ione, 6));
        }
        let mut jfiss = self.nfiss[1];
        if self.nfiss[2] > jfiss {
            jfiss = self.nfiss[2];
        }
        if jfiss == 0 {
            jfiss = 22;
        }
        let new = nold + newmat;
        let l7 = |n7: i32| -> String {
            [nlib, idat, ngnd, self.nfg, self.nrg, jfiss, n7].iter().map(|&x| ffmt::i(x, 6)).collect()
        };
        if self.mode == 0 {
            wl(c, &l7(nold));
        }
        if self.mode == 1 {
            wl(c, &l7(new));
        }
        if self.mode == 2 {
            wl(c, &l7(newmat));
        }
        let mut sec = 0; // 100, 102 or 106
        if self.mode == 2 {
            sec = 100;
        } else {
            wl(c, &ffmt::i(izero, 6));
            wl(c, &ffmt::i(izero, 6));
            if self.mode == 0 {
                sec = 102;
            } else if self.mode == 1 {
                sec = 106;
            }
        }
        if sec == 100 {
            let ngnd1 = ngnd + 1;
            wrows(c, &self.egb.0[..ngnd1 as usize]);
            wrows(c, &self.cflux.0[..ngnd as usize]);
            let zs = vec![ZERO; jfiss.max(0) as usize];
            if self.iu == 0 || self.if4 == 0 {
                wrows(c, &zs);
            }
            if self.iu == 1 && self.if4 > 0 {
                wrows(c, &self.uff.0[..jfiss as usize]);
            }
            if self.ipu == 0 || self.if4 == 0 {
                wrows(c, &zs);
            }
            if self.ipu == 1 && self.if4 > 0 {
                wrows(c, &self.puff.0[..jfiss as usize]);
            }
        } else if sec != 0 {
            if sec == 102 {
                for i in 1..=newmat {
                    let ip = self.ipost[i];
                    if ip >= 0 && ip <= 80 {
                        self.co_iflgs[ip as usize] = 1;
                    }
                }
                if self.mode == 0 {
                    wl(c, &flags(&self.co_iflgs[1..=nold as usize]));
                }
            }
            // 106
            if self.mode > 0 {
                wl(c, &flags(&vec![izero; nold as usize]));
            }
        }
        // 108
        let mut np1 = 0;
        let mut inf2 = 1;
        let mut mp1o = 0;
        let mut ip = 0;
        for i in 1..=newmat {
            if self.mode == 0 {
                ip = self.ipost[i];
            }
            if self.mode > 0 {
                ip = 0;
            }
            let mut nina = self.ninat[i];
            let mut aw = self.awrt[i];
            if nina < 2 {
                aw = aw * AMASSN;
            }
            if nina == 2 {
                nina = 1;
            }
            if nina == 3 {
                nina = 0;
            }
            self.ninat[i] = nina;
            let mut ident = self.idnt[i];
            if self.mode == 0 {
                ident = idnto(ip);
            }
            for j in 1..=np1o {
                if ident == IDP1O[(j - 1) as usize] {
                    mp1o += 1;
                }
            }
            let mut indfis = self.indf[i];
            let jres = self.irest[i];
            if jres > 0 && indfis == 0 {
                indfis = 1;
            }
            let mut nf = 0;
            if jres > 0 {
                nf = 1;
            }
            if ident == NFR2[(inf2 - 1) as usize] && self.mode == 0 {
                nf = 2;
            }
            if inf2 < nf2 && ident == NFR2[(inf2 - 1) as usize] && self.mode == 0 {
                inf2 += 1;
            }
            let ntemps = self.ntmp[i];
            let mut ip1 = 0;
            let ip1opt = self.ip1tab[i];
            if ip1opt == 0 {
                ip1 = 1;
            }
            if ip1 == 1 {
                np1 += 1;
            }
            let mut s = ffmt::f(aw, 12, 4);
            for x in [ident, indfis, jres, nf, ntemps, nina, izero, ip1] {
                s.push_str(&ffmt::i(x, 6));
            }
            wl(c, &s);
        }
        let _ = ip;

        //--write input for subfile 2.
        let mut ifile;
        ifile = 2;
        if self.mode == 0 {
            wl(c, &flags(&self.co_iflgs[1..=nold as usize]));
        }
        if self.mode == 1 {
            wl(c, &flags(&vec![izero; nold as usize]));
        }
        c.repoz(-nscr);
        if self.nnina != 0 {
            for _imat in 1..=self.nnina {
                for _ig in 1..=self.nrg {
                    let b = rd_reals(c, nscr, 4);
                    wrows(c, &b);
                }
            }
        }
        let mut flag = rd_reals(c, nscr, 1)[0];
        let mut ok = (flag - self.co_eflag).abs() <= eps;

        //--write input for subfile 3.
        let mut input = 0;
        if ok {
            ifile = 3;
            let nresot = nreso - nf2;
            let z3 = |c: &mut Ctx| {
                c.io.write_line(nout, &ffmt::i(izero, 6));
                c.io.write_line(nout, &ffmt::i(izero, 6));
                c.io.write_line(nout, &ffmt::i(izero, 6));
            };
            let mut go145 = false;
            let mut go170 = false;
            if self.mode == 1 {
                // 140
            } else if self.mode == 2 && self.nnina == 0 {
                go170 = true;
            } else if self.mode == 2 || self.nres > 0 {
                go145 = true;
            }
            if !go170 && !go145 {
                // 140
                for _i in 1..=nreso {
                    z3(c);
                }
                if self.nres == 0 {
                    go170 = true;
                } else {
                    go145 = true;
                }
            }
            if go145 {
                let mut kres = 1;
                let mut jnf2 = 0;
                for nmat in 1..=newmat {
                    let jres = self.irest[nmat];
                    if jres <= 0 {
                        continue;
                    }
                    let ident = self.idnt[nmat];
                    let mut iwrite = 0;
                    if self.mode <= 0 {
                        // 152
                        loop {
                            if ident < IDRESO[(kres - 1) as usize] {
                                break;
                            }
                            if ident <= IDRESO[(kres - 1) as usize] {
                                iwrite = 1;
                                let mut inf2 = 0;
                                while inf2 < nf2 {
                                    inf2 += 1;
                                    if ident == NFR2[(inf2 - 1) as usize] {
                                        z3(c);
                                    }
                                    if ident == NFR2[(inf2 - 1) as usize] {
                                        jnf2 += 1;
                                    }
                                }
                                if kres < nreso {
                                    kres += 1;
                                }
                                break;
                            }
                            // 154
                            z3(c);
                            let mut inf2 = 0;
                            while inf2 < nf2 {
                                inf2 += 1;
                                if IDRESO[(kres - 1) as usize] == NFR2[(inf2 - 1) as usize] {
                                    z3(c);
                                }
                                if IDRESO[(kres - 1) as usize] == NFR2[(inf2 - 1) as usize] {
                                    jnf2 += 1;
                                }
                            }
                            if kres == nreso {
                                break;
                            }
                            kres += 1;
                        }
                    }
                    // 156
                    if iwrite > 0 {
                        wl(c, &ffmt::i(ione, 6));
                    }
                    let xid = rd_reals(c, nscr, 1)[0];
                    wrows(c, &[xid]);
                    if iwrite > 0 {
                        wl(c, &ffmt::i(ione, 6));
                    }
                    let v = rd_ints(c, nscr, 2);
                    let (ntemps, nx) = (v[0], v[1]);
                    wl(c, &format!("{}{}", ffmt::i(ntemps, 6), ffmt::i(nx, 6)));
                    let b = rd_reals(c, nscr, ntemps);
                    wrows(c, &b);
                    let b = rd_reals(c, nscr, nx);
                    wrows(c, &b);
                    if iwrite > 0 {
                        wl(c, &ffmt::i(ione, 6));
                    }
                    let nw = ntemps * nx;
                    let mut indfis = self.indf[nmat];
                    loop {
                        // 158
                        for _j in 1..=self.nrg {
                            let b = rd_reals(c, nscr, nw);
                            wrows(c, &b);
                        }
                        if indfis < 2 {
                            break;
                        }
                        indfis = 1;
                    }
                }
                // 150 (after the loop)
                if !(self.mode > 0 || kres >= nresot) {
                    let need = nreso - kres + nf2 - jnf2;
                    for _i in 1..=need {
                        z3(c);
                    }
                }
            }
            // 170
            flag = rd_reals(c, nscr, 1)[0];
            ok = (flag - self.co_eflag).abs() <= eps;
        }

        //--write input for subfile 4.
        if ok {
            ifile = 4;
            if self.if4 != 0 {
                if self.mode != 2 {
                    wl(c, &flags(&self.co_iflgs[1..=nold as usize]));
                    let mut low = 1;
                    let mut lim = nold / 8;
                    if lim * 8 < nold {
                        lim += 1;
                    }
                    for _i in 1..=lim {
                        let ihi = low + 7;
                        let lr = low / 8;
                        for j in 1..=80 {
                            self.co_iflgs[j] = 0;
                        }
                        if self.mode != 1 {
                            for j in 1..=newmat {
                                let ip = self.ipost[j];
                                if ip >= low && ip <= ihi {
                                    let ipf = (ip - lr * 8 - 1) * 10;
                                    for q in 1..=5 {
                                        let k = ipf + q;
                                        if k >= 0 && k <= 80 {
                                            self.co_iflgs[k as usize] = 1;
                                        }
                                    }
                                }
                            }
                        }
                        wl(c, &flags(&self.co_iflgs[1..=80]));
                        low += 8;
                    }
                }
                let mut nmat = 0;
                while nmat < newmat {
                    nmat += 1;
                    let ntemps = self.ntmp[nmat];
                    let tempr = rd_reals(c, nscr, ntemps);
                    let mut s = String::new();
                    for (k, &x) in tempr.iter().enumerate() {
                        if k > 0 && k % 8 == 0 {
                            wl(c, &s);
                            s.clear();
                        }
                        s.push_str(&ffmt::f(x, 10, 0));
                    }
                    wl(c, &s);
                    let mut itemp = 0;
                    let nina = self.ninat[nmat];
                    let mut lim = 4;
                    if nina > 0 {
                        lim = 1;
                    }
                    while itemp < ntemps {
                        wl(c, &format!("{}{}", ffmt::i(ione, 6), ffmt::i(izero, 6)));
                        for _i in 1..=lim {
                            let b = rd_reals(c, nscr, ngnd);
                            wrows(c, &b);
                        }
                        if nina <= 0 {
                            let ind = 3;
                            input = 0;
                            wl(c, &format!("{}{}", ffmt::i(ind, 6), ffmt::i(input, 6)));
                            for _j in 1..=ngnd {
                                let v = rd_ints(c, nscr, 3);
                                let (ia, l1, l2) = (v[0], v[1], v[2]);
                                wl(c, &format!("{}{}{}", ffmt::i(ia, 6), ffmt::i(l1, 6), ffmt::i(l2, 6)));
                                let nw = l2 - l1 + 1;
                                let b = rd_reals(c, nscr, nw);
                                wrows(c, &b);
                            }
                        }
                        itemp += 1;
                    }
                }
            }
            flag = rd_reals(c, nscr, 1)[0];
            ok = (flag - self.co_eflag).abs() <= eps;
        }

        //--write input for subfile 5.
        if ok {
            ifile = 5;
            if self.if5 == 0 {
                wl(c, &ffmt::i(izero, 6));
            }
            if self.if5 > 0 {
                wl(c, &ffmt::i(ione, 6));
            }
            if self.if5 != 0 {
                let v = c.io.read_rec_with(nscr, |rr| {
                    let n = rr.i();
                    let mut k = vec![n];
                    for _ in 0..n.max(0) {
                        k.push(rr.i());
                    }
                    k
                });
                let nti = v[0];
                let mut s = vec![ffmt::i(nti, 6)];
                for ch in v[1..].chunks(12) {
                    s.push(ch.iter().map(|&x| ffmt::i(x, 6)).collect());
                }
                if v.len() == 1 {
                    s.push(String::new());
                }
                for l in s {
                    wl(c, &l);
                }
                let v = c.io.read_rec_with(nscr, |rr| {
                    let n = rr.i();
                    let mut k = vec![n];
                    for _ in 0..n.max(0) {
                        k.push(rr.i());
                    }
                    k
                });
                let nfi = v[0];
                let mut s = vec![ffmt::i(nfi, 6)];
                for ch in v[1..].chunks(12) {
                    s.push(ch.iter().map(|&x| ffmt::i(x, 6)).collect());
                }
                if v.len() == 1 {
                    s.push(String::new());
                }
                for l in s {
                    wl(c, &l);
                }
                for _i in 1..=nti {
                    let v = rd_reals(c, nscr, 1 + nfi);
                    wrows(c, &v[..1]);
                    wrows(c, &v[1..]);
                }
                self.co_eflag = rd_reals(c, nscr, 1)[0];
                ok = (flag - self.co_eflag).abs() <= eps;
            }
        }

        //--write input for subfile 6.
        if ok {
            ifile = 6;
            if self.if4 != 0 {
                let mut ind = 0;
                if mp1o > 0 {
                    ind = 1;
                }
                wl(c, &format!("{}{}{}", ffmt::i(ind, 6), ffmt::i(ithree, 6), ffmt::i(input, 6)));
                if np1 != 0 {
                    if self.mode != 2 && ind != 1 {
                        wl(c, &ffmt::i(ione, 6));
                    }
                    let v = c.io.read_rec_with(nscr, |rr| {
                        let n = rr.i();
                        let mut k = vec![n];
                        for _ in 0..n.max(0) {
                            k.push(rr.i());
                        }
                        k
                    });
                    np1 = v[0];
                    // '(7i6)'
                    for ch in v.chunks(7) {
                        wl(c, &ch.iter().map(|&x| ffmt::i(x, 6)).collect::<String>());
                    }
                    let mut ip1 = 0;
                    while ip1 < np1 {
                        ip1 += 1;
                        let (ntemps, b) = c.io.read_rec_with(nscr, |rr| {
                            let n = rr.i();
                            let b: Vec<f64> = (0..n.max(0)).map(|_| rr.r()).collect();
                            (n, b)
                        });
                        self.ntmp[ip1] = ntemps;
                        // '(i6,1p,6e12.5)'
                        let mut s = ffmt::i(ntemps, 6);
                        let first = b.len().min(6);
                        for &x in &b[..first] {
                            s.push_str(&pe(x, 12, 5));
                        }
                        wl(c, &s);
                        for ch in b[first..].chunks(6) {
                            wl(c, &ch.iter().map(|&x| pe(x, 12, 5)).collect::<String>());
                        }
                    }
                    let mut ip1 = 0;
                    while ip1 < np1 {
                        ip1 += 1;
                        let mut itemp = 0;
                        let ntemps = self.ntmp[ip1];
                        while itemp < ntemps {
                            for _j in 1..=ngnd {
                                let v = rd_ints(c, nscr, 3);
                                let (ia, l1, l2) = (v[0], v[1], v[2]);
                                wl(c, &format!("{}{}{}", ffmt::i(ia, 6), ffmt::i(l1, 6), ffmt::i(l2, 6)));
                                let nw = l2 - l1 + 1;
                                let b = rd_reals(c, nscr, nw);
                                wrows(c, &b);
                            }
                            itemp += 1;
                        }
                    }
                } else {
                    wl(c, &ffmt::i(izero, 6));
                    wl(c, &ffmt::i(izero, 6));
                }
            }
            flag = rd_reals(c, nscr, 1)[0];
            ok = (flag - self.co_eflag).abs() <= eps;
        }

        if ok {
            //--write end of input flag.
            let flag = rd_reals(c, nscr, 1)[0];
            wl(c, &a4_of(flag));

            //--cpmout is finished.
            c.closz(nscr);
            return;
        }

        //--end flag not found on input tape
        let strng = format!("expected end flag not found on nscr for file {}", ffmt::i(ifile, 2));
        c.error("cpmout", &strng, " ");
    }
}
