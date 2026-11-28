//! Port of LEGEND (PREPRO 2026-1): linearized angular distributions
//! (MF=4) from Legendre coefficients or tables, with thinning and
//! corrections of negative distributions. Input `LEGEND.INP`, listing
//! `LEGEND.LST`, work file `LEGEND.TMP`.
//!
//! Arrays are 1-based as in the Fortran.

use crate::endfio::{f9s, fint, in9, out9, Endf};
use crate::ffmt::{f, i};
use crate::fio::{self, Rd, L};
use crate::linear::{blank, rdi};

const MAXPOINT: usize = 120000;
const MAXMAT: usize = 100;
const MAXLEGS: usize = 500;
const ISCR: i32 = 12;

fn dash() -> String {
    format!(" {}", "-".repeat(78))
}

struct Lg {
    io: Endf,
    // /TABLET/
    etab1: Vec<f64>,
    etab2: Vec<f64>,
    sigtab: Vec<f64>,
    dfltab: Vec<f64>,
    matab1: Vec<i32>,
    mftab1: Vec<i32>,
    mtab1: Vec<i32>,
    matab2: Vec<i32>,
    mftab2: Vec<i32>,
    mtab2: Vec<i32>,
    itab: usize,
    // /LEGCONF/
    flmult: Vec<f64>,
    pcon1: Vec<f64>,
    pcon2: Vec<f64>,
    values: Vec<f64>,
    dfl: Vec<f64>,
    sumn: Vec<f64>,
    pn: Vec<f64>,
    fl: Vec<f64>,
    flold: Vec<f64>,
    flkeep: Vec<[f64; 4]>,
    xmubase: Vec<f64>,
    ekeep: [f64; 4],
    elast4: f64,
    sum: f64,
    dflchn: f64,
    dflmax: f64,
    // /LEGCONI/
    lvalue: Vec<i32>,
    newfl: Vec<i32>,
    lkeep: [i32; 4],
    legin: i32,
    leguse: i32,
    kpass: i32,
    nbt: [i32; 101],
    int: [i32; 101],
    nbtc: [i32; 201],
    intc: [i32; 201],
    xmuin: Vec<f64>,
    ptin: Vec<f64>,
    xmuout: Vec<f64>,
    ptout: Vec<f64>,
    // other COMMON
    errmax: f64,
    errok: f64,
    erruse: f64,
    loope: i32,
    nptmax: usize,
    legmax: i32,
    e: f64,
    izalst: i32,
    matlst: i32,
    mtlst: i32,
    ne: i32,
    je: i32,
    typed: [i32; 7],
    imtabs: i32,
    imlegs: i32,
    imfix: i32,
    imvar: i32,
    n2in: usize,
    n2out: usize,
    sigmin: f64,
    field6: [[u8; 11]; 6],
    // SAVEd locals
    leglst: i32,
    th_dsign: f64,
    th_slpmax: f64,
    th_slpmin: f64,
    th_ithin: usize,
    contol_lstmat: i32,
    contol_lstmf: i32,
}

impl Lg {
    fn lst<S: AsRef<[u8]>>(&mut self, s: S) {
        self.io.lst(s);
    }
    fn both<S: AsRef<[u8]>>(&mut self, s: S) {
        self.io.both(s);
    }

    fn mywant(&self, mat: i32, mf: i32, mt: i32) -> i32 {
        let mut needed = 0;
        for k in 1..=self.itab {
            if mat > self.matab2[k] {
                continue;
            }
            needed = 1;
            if mat < self.matab1[k] {
                continue;
            }
            if mf < self.mftab1[k] || mf > self.mftab2[k] {
                continue;
            }
            if mt < self.mtab1[k] || mt > self.mtab2[k] {
                continue;
            }
            return 1;
        }
        if needed <= 0 {
            -1
        } else {
            0
        }
    }

    fn contil(&mut self) {
        loop {
            self.io.conti();
            if self.io.h.mt > 0 {
                let m = self.mywant(self.io.h.mat, self.io.h.mf, self.io.h.mt);
                if m < 0 {
                    break;
                }
                if m > 0 {
                    return;
                }
                self.io.skips();
                continue;
            }
            if self.io.h.mat < 0 {
                break;
            }
        }
        self.io.h.mat = -1;
        self.io.h.mf = 0;
        self.io.h.mt = 0;
    }

    fn contol(&mut self) {
        if self.contol_lstmat != -9999 && !(self.contol_lstmf == self.io.h.mf && self.contol_lstmat == self.io.h.mat) {
            let m = self.contol_lstmat;
            self.io.outf(m);
            if self.contol_lstmat != self.io.h.mat {
                self.io.outm();
            }
        }
        self.contol_lstmat = self.io.h.mat;
        self.contol_lstmf = self.io.h.mf;
        self.io.conto();
    }

    /// CARDI into LEADER C1, C2, L1, L2 only (N1, N2 are other variables).
    fn cardi4(&mut self) -> (i32, i32) {
        let (c1, c2, l1, l2, n1, n2) = self.io.cardi();
        self.io.ld.c1 = c1;
        self.io.ld.c2 = c2;
        self.io.ld.l1 = l1;
        self.io.ld.l2 = l2;
        (n1, n2)
    }

    fn cardo_ld(&mut self, n1: i32, n2: i32) {
        let ld = self.io.ld;
        self.io.cardo(ld.c1, ld.c2, ld.l1, ld.l2, n1, n2);
    }

    fn file1(&mut self) {
        let texts = [
            " ***************** Program LEGEND (VERSION 2026-1) ***********",
            " Linearized Angular Distributions from Tables or Coefficients ",
            " Allowable Error of12345678901 per-cent and123456 Max. Points ",
            " Tabulated Data...Copied                                      ",
            " Tabulated Data...Linearized                                  ",
            " Tabulated Data...Linearized/Thinned                          ",
            " Legendre Data....Copied                                      ",
            " Legendre Data....Reconstructed Table/Output Table            ",
            " Legendre Data....Reconstructed Table/Output Coefficients     ",
            " Correction.......None                                        ",
            " Correction.......Modify Coefficients                         ",
            " Correction.......Force Distributions to be Positive          ",
        ];
        let mut pd = vec![[b' '; 66]; 13];
        for (k, s) in texts.iter().enumerate() {
            pd[k + 1][..s.len()].copy_from_slice(s.as_bytes());
        }
        for k in 62..66 {
            pd[1][k] = b'*';
        }
        let (c1, c2, l1, l2, n1, n2) = self.io.cardi_ld();
        let (mut c1, mut c2, mut l1, mut l2, mut n1, mut n2) = (c1, c2, l1, l2, n1, n2);
        if !(n1 > 0 && n2 == 0) {
            let n2x = n2;
            self.io.cardo(c1, c2, l1, l2, n1, n2);
            (c1, c2, l1, l2, n1, n2) = self.io.cardi_ld();
            if n2x > 0 {
                self.io.cardo(c1, c2, l1, l2, n1, n2);
                (c1, c2, l1, l2, n1, n2) = self.io.cardi_ld();
            }
        }
        let n1in = n1;
        n1 += 6;
        self.io.ld.n1 = n1;
        self.io.cardo(c1, c2, l1, l2, n1, n2);
        for _ in 0..n1in.max(0) {
            self.io.copy1();
        }
        self.io.hollyo(&pd[1]);
        self.io.hollyo(&pd[2]);
        out9(self.errmax, &mut pd[3][19..30]);
        intout(self.nptmax as i32, &mut pd[3][43..49]);
        self.io.hollyo(&pd[3]);
        self.io.hollyo(&pd[4 + self.imtabs as usize]);
        self.io.hollyo(&pd[7 + self.imlegs as usize]);
        self.io.hollyo(&pd[10 + self.imfix as usize]);
    }

    fn xint(x: &[f64], y: &[f64], n: usize) -> f64 {
        let mut s = 0.0;
        if n >= 2 {
            for j in 2..=n {
                s += (x[j] - x[j - 1]) * (y[j] + y[j - 1]);
            }
        }
        s / 2.0
    }

    /// `ANGLE1`: the distribution at a cosine (SUM, PN, SUMN).
    fn angle1(&mut self, xmudum: f64) {
        let mut xmu = xmudum;
        if xmu.abs() < 1.0e-10 {
            xmu = 0.0;
        }
        self.sum = 0.5;
        if self.leguse < 1 {
            return;
        }
        self.pn[1] = xmu;
        self.sumn[1] = self.flmult[1] * self.fl[1] * self.pn[1];
        self.sum += self.sumn[1];
        if self.leguse == 1 {
            return;
        }
        self.pn[2] = self.pcon1[2] * xmu * self.pn[1] + self.pcon2[2];
        self.sumn[2] = self.flmult[2] * self.fl[2] * self.pn[2];
        self.sum += self.sumn[2];
        if self.leguse == 2 {
            return;
        }
        for l in 3..=self.leguse as usize {
            self.pn[l] = self.pcon1[l] * xmu * self.pn[l - 1] + self.pcon2[l] * self.pn[l - 2];
            self.sumn[l] = self.flmult[l] * self.fl[l] * self.pn[l];
            self.sum += self.sumn[l];
        }
    }

    fn anglen(&mut self) {
        if self.leguse <= 1 {
            self.n2out = 2;
            self.xmuout[1] = -1.0;
            self.xmuout[2] = 1.0;
            self.angle1(-1.0);
            self.ptout[1] = self.sum;
            self.angle1(1.0);
            self.ptout[2] = self.sum;
            return;
        }
        let lu = self.leguse as usize;
        let nmid = 2 * lu + 1;
        let nstart = nmid + 2 * lu;
        let dmu = 1.0 / ((nmid - 1) as f32 as f64);
        self.xmubase[1] = -1.0;
        for k in 2..=nstart {
            self.xmubase[k] = self.xmubase[k - 1] + dmu;
            if k == nmid {
                self.xmubase[k] = 0.0;
            }
        }
        self.xmubase[nstart] = 1.0;
        self.erruse = self.errok;
        self.xmuout[1] = -1.0;
        self.angle1(-1.0);
        self.ptout[1] = self.sum;
        'l30: loop {
            let mut nbase = 1usize;
            self.n2out = 1;
            'l40: loop {
                nbase += 1;
                let xb = self.xmubase[nbase];
                self.angle1(xb);
                let ptbase = self.sum;
                'l50: loop {
                    self.n2out += 1;
                    if self.n2out > MAXPOINT {
                        self.erruse *= 2.0;
                        continue 'l30;
                    }
                    while nbase < nstart && self.xmubase[nbase] < self.xmuout[self.n2out - 1] {
                        nbase += 1;
                    }
                    let n = self.n2out;
                    self.xmuout[n] = self.xmubase[nbase];
                    self.ptout[n] = ptbase;
                    loop {
                        let xmumid = 0.5 * (self.xmuout[n] + self.xmuout[n - 1]);
                        self.angle1(xmumid);
                        let ptmid = self.sum;
                        let ptav = 0.5 * (self.ptout[n] + self.ptout[n - 1]);
                        if (ptmid - ptav).abs() <= self.erruse * ptmid.abs() {
                            break;
                        }
                        if ptmid < 0.0 {
                            break;
                        }
                        self.xmuout[n] = xmumid;
                        self.ptout[n] = ptmid;
                    }
                    if self.xmuout[n] < self.xmubase[nbase] {
                        continue 'l50;
                    }
                    if nbase < nstart {
                        continue 'l40;
                    }
                    return;
                }
            }
        }
    }

    fn minnie(&mut self) {
        if self.n2out < 3 {
            return;
        }
        self.n2in = 2;
        self.xmuin[1] = self.xmuout[1];
        self.ptin[1] = self.ptout[1];
        self.xmuin[2] = self.xmuout[2];
        self.ptin[2] = self.ptout[2];
        let mut dsign = 1.0;
        if self.ptout[2] < self.ptout[1] {
            dsign = -1.0;
        }
        let n2out = self.n2out;
        let (mut xmumin, mut ptmin) = (0.0, 0.0);
        for i2 in 3..=n2out {
            let im1 = i2 - 1;
            let dp = self.ptout[i2] - self.ptout[im1];
            let mut go70 = true;
            if !(dp * dsign >= 0.0) {
                dsign = -dsign;
                if dp > 0.0 {
                    go70 = false;
                }
            }
            if !go70 {
                let mut xe = [0.0, self.xmuout[i2 - 2], self.xmuout[im1], self.xmuout[i2]];
                let mut pe = [0.0, self.ptout[i2 - 2], self.ptout[im1], self.ptout[i2]];
                xmumin = xe[2];
                ptmin = pe[2];
                for _ in 0..10 {
                    let mut hit = 0;
                    let mut xmumid = 0.0;
                    let mut ptmid = 0.0;
                    for j in 1..=2 {
                        xmumid = (xe[j] + xe[j + 1]) / 2.0;
                        self.angle1(xmumid);
                        ptmid = self.sum;
                        if ptmid < ptmin {
                            hit = j;
                            break;
                        }
                    }
                    if hit == 0 {
                        break;
                    }
                    let j = hit;
                    xe[1] = xe[j];
                    xe[3] = xe[j + 1];
                    pe[1] = pe[j];
                    pe[3] = pe[j + 1];
                    xe[2] = xmumid;
                    pe[2] = ptmid;
                    xmumin = xmumid;
                    ptmin = ptmid;
                }
                // label 40
                let mut go50 = self.n2out >= MAXPOINT;
                let mut skip = false;
                if !go50 {
                    let dmu = xmumin - self.xmuin[self.n2in];
                    if dmu.abs() < 0.001 {
                        go50 = true;
                    } else if dmu < 0.0 {
                        let n2inm1 = self.n2in;
                        self.n2in += 1;
                        let n = self.n2in;
                        self.xmuin[n] = self.xmuin[n2inm1];
                        self.ptin[n] = self.ptin[n2inm1];
                        self.xmuin[n2inm1] = xmumin;
                        self.ptin[n2inm1] = ptmin;
                        skip = true;
                    } else if dmu == 0.0 {
                        skip = true;
                    } else {
                        self.n2in += 1;
                        go50 = true;
                    }
                }
                if go50 && !skip {
                    let n = self.n2in;
                    self.xmuin[n] = xmumin;
                    self.ptin[n] = ptmin;
                }
            }
            // label 70
            self.n2in += 1;
            let n = self.n2in;
            self.xmuin[n] = self.xmuout[i2];
            self.ptin[n] = self.ptout[i2];
        }
        let _ = (xmumin, ptmin);
        self.n2out = self.n2in;
        for k in 1..=self.n2out {
            self.xmuout[k] = self.xmuin[k];
            self.ptout[k] = self.ptin[k];
        }
    }

    fn thin(&mut self) {
        if self.n2out < 3 || self.errok <= 0.0 {
            return;
        }
        let nstart = self.n2out;
        self.th_dsign = 1.0;
        if self.ptout[1] > self.ptout[2] {
            self.th_dsign = -1.0;
        }
        self.n2out = 1;
        self.th_ithin = 2;
        #[derive(Clone, Copy, PartialEq)]
        enum G {
            L10,
            L20,
            L30,
            L40,
            L50,
            L70,
            L80,
            L90,
            L100,
            L110,
        }
        for m in 2..=nstart {
            let mm1 = m - 1;
            let mut dy = self.ptout[m] - self.ptout[mm1];
            let mut g = G::L10;
            if !(dy * self.th_dsign >= 0.0) {
                self.th_dsign = -self.th_dsign;
                if self.xmuout[m] <= self.xmuout[mm1] {
                    g = G::L30;
                } else if m > self.th_ithin {
                    g = G::L70;
                }
            }
            loop {
                match g {
                    G::L10 => {
                        g = if self.ptout[m] * self.ptout[mm1] >= 0.0 {
                            G::L20
                        } else if m <= self.th_ithin {
                            G::L90
                        } else {
                            G::L80
                        };
                    }
                    G::L20 => g = if self.xmuout[m] > self.xmuout[mm1] { G::L50 } else { G::L30 },
                    G::L30 => {
                        dy = dy.abs();
                        g = if m != self.th_ithin {
                            G::L40
                        } else if dy <= 0.0 {
                            G::L100
                        } else {
                            G::L90
                        };
                    }
                    G::L40 => g = if dy <= 0.0 { G::L90 } else { G::L80 },
                    G::L50 => {
                        if self.xmuout[m] == 0.0 {
                            g = G::L90;
                            continue;
                        }
                        let n = self.n2out;
                        let dx = self.xmuout[m] - self.xmuout[n];
                        let slope = (self.ptout[m] - self.ptout[n]) / dx;
                        if m == self.th_ithin {
                            let dslope = self.errok * self.ptout[m] / dx;
                            self.th_slpmax = slope + dslope;
                            self.th_slpmin = slope - dslope;
                            g = G::L110;
                            continue;
                        }
                        if slope > self.th_slpmax || slope < self.th_slpmin {
                            g = G::L70;
                            continue;
                        }
                        let dslope = self.errok * self.ptout[m] / dx;
                        let slp1 = slope + dslope;
                        if slp1 < self.th_slpmax {
                            self.th_slpmax = slp1;
                        }
                        let slp2 = slope - dslope;
                        if slp2 > self.th_slpmin {
                            self.th_slpmin = slp2;
                        }
                        g = G::L110;
                    }
                    G::L70 => {
                        self.n2out += 1;
                        let n = self.n2out;
                        self.xmuout[n] = self.xmuout[mm1];
                        self.ptout[n] = self.ptout[mm1];
                        self.th_ithin = m;
                        g = G::L50;
                    }
                    G::L80 => {
                        self.n2out += 1;
                        let n = self.n2out;
                        self.xmuout[n] = self.xmuout[mm1];
                        self.ptout[n] = self.ptout[mm1];
                        g = G::L90;
                    }
                    G::L90 => {
                        self.n2out += 1;
                        let n = self.n2out;
                        self.xmuout[n] = self.xmuout[m];
                        self.ptout[n] = self.ptout[m];
                        g = G::L100;
                    }
                    G::L100 => {
                        self.th_ithin = m + 1;
                        g = G::L110;
                    }
                    G::L110 => break,
                }
            }
        }
        if self.th_ithin > nstart {
            return;
        }
        self.n2out += 1;
        let n = self.n2out;
        self.xmuout[n] = self.xmuout[nstart];
        self.ptout[n] = self.ptout[nstart];
    }

    fn postiv(&self) -> (f64, f64, i32) {
        let mut xmumin = self.xmuout[1];
        let mut ptmin = self.ptout[1];
        let mut imneg = 0;
        for k in 1..=self.n2out {
            if self.ptout[k] < 0.0 {
                imneg += 1;
            }
            if self.ptout[k] >= ptmin {
                continue;
            }
            xmumin = self.xmuout[k];
            ptmin = self.ptout[k];
        }
        (xmumin, ptmin, imneg)
    }

    fn sort(&mut self, l: usize) {
        if l <= 1 {
            return;
        }
        let mut ltop = l;
        for _ in 2..=l {
            let mut iswtch = 0;
            let mut jm1 = 1;
            for j in 2..=ltop {
                if !(self.values[j] >= self.values[jm1]) {
                    self.values.swap(j, jm1);
                    self.lvalue.swap(j, jm1);
                    iswtch = 1;
                }
                jm1 = j;
            }
            if iswtch <= 0 {
                return;
            }
            ltop -= 1;
        }
    }

    fn lowleg(&mut self, xmu: f64, ptmin: &mut f64) -> i32 {
        let mut newleg = 0;
        self.angle1(xmu);
        let lu = self.leguse.max(0) as usize;
        for l in 1..=lu {
            self.lvalue[l] = l as i32;
            self.values[l] = self.sumn[l];
        }
        self.sort(lu);
        if self.kpass == 1 {
            let s = format!(
                " Requested Minimum Value={}\n Negative Contribution of Each Legendre Order and per-cent\n Changed to Make Positive by Changing One Coefficient\n Index   Coefficient  P-n(Cosine)  Contribution   Change  Allowed",
                f(self.sigmin, 11, 6)
            );
            self.lst(&s);
        }
        loop {
            // label 20
            let mut ilow = 0usize;
            let mut ratlow = 100.0;
            for k in 1..=lu {
                let ii = self.lvalue[k] as usize;
                if self.newfl[ii] > 0 {
                    continue;
                }
                if self.values[k] >= 0.0 {
                    break;
                }
                self.dfl[k] = (self.sigmin - *ptmin) / (self.flmult[ii] * self.pn[ii]);
                let ratio = self.dfl[k] / self.fl[ii];
                let absrat = ratio.abs();
                if absrat > 10.0 {
                    break;
                }
                if absrat < ratlow {
                    ilow = k;
                    ratlow = absrat;
                }
                if self.kpass == 1 {
                    let s = format!(
                        "{}{}{}{}{}{}",
                        i(ii as i32, 6),
                        f(self.fl[ii], 14, 6),
                        f(self.pn[ii], 14, 6),
                        f(self.values[k], 14, 6),
                        f(100.0 * ratio, 9, 4),
                        f(100.0 * self.dflmax, 9, 4)
                    );
                    self.lst(&s);
                }
            }
            // label 50
            if ilow == 0 {
                return newleg;
            }
            self.kpass = 2;
            let ii = self.lvalue[ilow] as usize;
            let mut dfluse = self.dfl[ilow];
            let fltry = self.fl[ii] + dfluse;
            if (fltry - self.flold[ii]).abs() > self.dflmax * self.flold[ii].abs() {
                dfluse = -self.dflmax * self.flold[ii];
            }
            self.fl[ii] += dfluse;
            self.newfl[ii] = 1;
            *ptmin += self.flmult[ii] * dfluse * self.pn[ii];
            newleg += 1;
            if *ptmin < 0.99 * self.sigmin {
                continue;
            }
            return -1;
        }
    }

    fn lowtab(&mut self, ptmin: f64, ltt: i32) {
        let addit = (self.sigmin - ptmin) / (1.0 - 2.0 * self.sigmin);
        if self.n2out > 0 {
            for k in 1..=self.n2out {
                self.ptout[k] += addit;
            }
            if ltt == 2 {
                let s = format!("{} Added to Each Point to Make Distribution Positive", f(addit, 10, 5));
                self.lst(&s);
            }
        }
        if self.leguse <= 0 {
            return;
        }
        let sumx = 1.0 + 2.0 * addit;
        for l in 1..=self.leguse as usize {
            self.flold[l] = self.fl[l];
            self.fl[l] /= sumx;
            self.newfl[l] = 1;
        }
        let s = format!(" Each Legendre Coefficient Divided by{} to Make Distribution Positive", f(sumx, 10, 5));
        self.lst(&s);
    }

    fn tab4(&mut self) {
        let sumx = Self::xint(&self.xmuout, &self.ptout, self.n2out);
        for k in 1..=self.n2out {
            self.ptout[k] /= sumx;
        }
        self.io.ld.n1 = 1;
        self.io.ld.n2 = self.n2out as i32;
        self.cardo_ld(1, self.n2out as i32);
        self.io.terpo(&[self.n2out as i32], &[2], 1);
        self.io.minus3 = 0;
        self.io.implus = 0;
        let x = std::mem::take(&mut self.xmuout);
        let y = std::mem::take(&mut self.ptout);
        self.io.pointo(&x[1..], &y[1..], self.n2out);
        self.xmuout = x;
        self.ptout = y;
    }

    fn list4(&mut self) {
        let mut n1 = self.legin;
        if n1 <= 0 {
            n1 = 1;
            self.fl[1] = 0.0;
        } else if self.leguse < n1 {
            for k in (self.leguse + 1) as usize..=n1 as usize {
                self.fl[k] = 0.0;
            }
        }
        self.io.ld.n1 = n1;
        let n2 = self.io.ld.n2;
        self.cardo_ld(n1, n2);
        let fl = std::mem::take(&mut self.fl);
        self.io.listo9(&fl[1..], n1 as usize);
        self.fl = fl;
    }

    fn evar1(&mut self) {
        let d = dash();
        if self.je <= 1 {
            return;
        }
        out9(self.ekeep[2], &mut self.field6[0]);
        out9(self.ekeep[3], &mut self.field6[1]);
        if self.lkeep[2] > self.lkeep[3] {
            let s = L::new()
                .s(&d)
                .s("\n WARNING - ")
                .s(&i(self.lkeep[2], 3))
                .s(" Coefficients at E=")
                .b(&self.field6[0])
                .s(" eV\n       Only")
                .s(&i(self.lkeep[3], 3))
                .s(" Coefficients at E=")
                .b(&self.field6[1])
                .s(" eV\n")
                .s(&d);
            self.lst(&s);
        }
        let ltop = self.lkeep[2].max(self.lkeep[3]);
        if ltop > 0 {
            let mut ierr = 0;
            for l in 1..=ltop as usize {
                if self.lkeep[2] < l as i32 || self.lkeep[3] < l as i32 {
                    continue;
                }
                let (a, b) = (self.flkeep[l][2], self.flkeep[l][3]);
                if !(a.abs() > b.abs() || (a > 0.0 && b < 0.0) || (a < 0.0 && b > 0.0)) {
                    continue;
                }
                if ierr == 0 {
                    self.lst(format!("{d}\n WARNING - Coefficient Has Decreased or Changed Sign"));
                }
                ierr = 1;
                out9(a, &mut self.field6[2]);
                out9(b, &mut self.field6[3]);
                let s = L::new()
                    .s(" L=")
                    .s(&i(l as i32, 3))
                    .s(" ")
                    .b(&self.field6[2])
                    .s(" at E=")
                    .b(&self.field6[0])
                    .s(" eV\n       ")
                    .b(&self.field6[3])
                    .s(" at E=")
                    .b(&self.field6[1])
                    .s(" eV");
                self.lst(&s);
            }
            if ierr != 0 {
                self.lst(&d);
            }
        }
        let ltop = self.lkeep[3];
        if ltop <= 0 {
            return;
        }
        let mut flast = 1.0f64;
        let mut ierr = 0;
        for l in 1..=ltop as usize {
            if self.flkeep[l][3].abs() > flast.abs() {
                if ierr == 0 {
                    self.lst(format!("{d}\n WARNING - Coefficients do NOT Decrease with L"));
                }
                out9(flast, &mut self.field6[2]);
                out9(self.flkeep[l][3], &mut self.field6[3]);
                let s = L::new()
                    .s(" L=")
                    .s(&i(l as i32 - 1, 3))
                    .s(" ")
                    .b(&self.field6[2])
                    .s(" L=")
                    .s(&i(l as i32, 3))
                    .s(" ")
                    .b(&self.field6[3]);
                self.lst(&s);
                ierr = 1;
            }
            flast = self.flkeep[l][3];
        }
        if ierr != 0 {
            self.lst(&d);
        }
    }

    fn find_jtab(&self) -> usize {
        for j in 1..=self.itab {
            if self.matlst < self.matab1[j] || self.matlst > self.matab2[j] {
                continue;
            }
            if self.mtlst < self.mtab1[j] || self.mtlst > self.mtab2[j] {
                continue;
            }
            if self.e < self.etab1[j] || self.e > self.etab2[j] {
                continue;
            }
            return j;
        }
        self.itab
    }

    fn angleg(&mut self) {
        let dfltim = [0.0, 0.25, 0.5, 1.0];
        let d = dash();
        self.legin = self.io.ld.n1;
        self.leguse = self.legin;
        if self.legin >= 1 {
            let mut p = std::mem::take(&mut self.ptout);
            self.io.listi(&mut p[1..], self.legin as usize);
            self.ptout = p;
            if self.leguse > self.legmax {
                self.leguse = self.legmax;
            }
            let legtop = self.leguse;
            for _ in 0..legtop {
                if self.ptout[self.leguse as usize] != 0.0 {
                    break;
                }
                self.leguse -= 1;
            }
            if self.leguse <= 0 {
                self.leguse = 1;
            }
            for l in 1..=self.leguse as usize {
                self.fl[l] = self.ptout[l];
                self.flold[l] = self.fl[l];
                self.newfl[l] = 0;
            }
        }
        // label 40
        if self.je > 1 {
            for l in 2..=3 {
                let lm1 = l - 1;
                let lego = self.lkeep[l];
                self.lkeep[lm1] = self.lkeep[l];
                self.ekeep[lm1] = self.ekeep[l];
                if lego > 0 {
                    for j in 1..=lego as usize {
                        self.flkeep[j][lm1] = self.flkeep[j][l];
                    }
                }
            }
        }
        self.lkeep[3] = self.leguse;
        self.ekeep[3] = self.e;
        if self.leguse > 0 {
            for j in 1..=self.leguse as usize {
                self.flkeep[j][3] = self.flold[j];
            }
        } else {
            self.leguse = 0;
            self.fl[1] = 0.0;
        }
        self.anglen();
        self.minnie();
        let sumx = Self::xint(&self.xmuout, &self.ptout, self.n2out);
        out9(self.e, &mut self.field6[0]);
        let head = format!("{}{}{}", i(self.izalst, 7), i(self.matlst, 5), i(self.mtlst, 4));
        let s = if self.legin == self.leguse {
            L::new().s(&head).b(&self.field6[0]).s(&i(self.legin, 5)).s(&" ".repeat(17)).s(&i(self.n2out as i32, 11)).s(&f(sumx, 11, 6))
        } else {
            L::new()
                .s(&head)
                .b(&self.field6[0])
                .s(&i(self.legin, 5))
                .s(" (")
                .s(&i(self.leguse, 4))
                .s(")")
                .s(&" ".repeat(10))
                .s(&i(self.n2out as i32, 11))
                .s(&f(sumx, 11, 6))
        };
        self.both(&s);
        if self.imvar > 0 {
            if self.je > 1 && self.leguse > self.leglst {
                let mut lout = 0;
                for l in (self.leglst + 1) as usize..=self.leguse as usize {
                    if self.flkeep[l][3].abs() == 0.0 {
                        continue;
                    }
                    if lout == 0 {
                        self.lst(format!("{d}\n WARNING - Coefficients do NOT Start Equal to Zero.\n Interpolation to Lower Energies May Cause Errors.\n   L  Coefficient"));
                    }
                    lout = 1;
                    out9(self.flkeep[l][3], &mut self.field6[0]);
                    let s = L::new().s(&i(l as i32, 4)).s("  ").b(&self.field6[0]);
                    self.lst(&s);
                }
                if lout != 0 {
                    self.lst(&d);
                }
            }
            self.leglst = self.leguse;
            self.evar1();
        }
        if (sumx - 1.0).abs() > 0.05 {
            self.lst(format!("{d}\n WARNING - Angular Distribution NOT Normalized to within 5%\n{d}"));
        }
        if !(self.e > self.elast4) {
            out9(self.elast4, &mut self.field6[0]);
            out9(self.e, &mut self.field6[1]);
            let s = L::new()
                .s(&d)
                .s("\n WARNING - Energy is NOT in Ascending Order\n Last E =")
                .b(&self.field6[0])
                .s(" eV, Current E =")
                .b(&self.field6[1])
                .s(" eV");
            self.lst(&s);
            self.typed[6] += 1;
        }
        self.elast4 = self.e;
        if self.legin > self.legmax {
            let s = format!(
                "{d}\n WARNING - Too Many Legendre Coefficients\n Read{} Only Used First{} All Others Ignored\n{d}",
                i(self.legin, 4),
                i(self.legmax, 4)
            );
            self.lst(&s);
        }
        let (xmumin, mut ptmin, imneg) = self.postiv();
        'l270: {
            if imneg <= 0 {
                break 'l270;
            }
            let s = format!(
                "{d}\n WARNING - Negative Angular Distribution at{} Points.\n Minimum at Cosine={} Value={}\n{d}",
                i(imneg, 5),
                f(xmumin, 11, 6),
                f(ptmin, 11, 7)
            );
            self.lst(&s);
            self.typed[5] += 1;
            if self.imfix <= 0 {
                break 'l270;
            }
            let jtab = self.find_jtab();
            self.sigmin = self.sigtab[jtab];
            self.dflchn = self.dfltab[jtab];
            if self.imfix == 2 {
                self.lowtab(ptmin, 1);
                break 'l270;
            }
            self.anglen();
            self.minnie();
            let ptminx = ptmin;
            self.kpass = 0;
            let mut xmumin = xmumin;
            let mut ok250 = false;
            'l220: for looppc in 1..=3 {
                self.dflmax = self.dflchn * dfltim[looppc];
                ptmin = ptminx;
                for l in 1..=self.leguse as usize {
                    self.fl[l] = self.flold[l];
                    self.newfl[l] = 0;
                }
                for _loopl in 1..=self.leguse {
                    self.kpass += 1;
                    let newleg = self.lowleg(xmumin, &mut ptmin);
                    self.anglen();
                    if newleg >= 0 {
                        continue 'l220;
                    }
                    let (xm, pm, im) = self.postiv();
                    xmumin = xm;
                    ptmin = pm;
                    if im <= 0 && ptmin >= self.sigmin {
                        ok250 = true;
                        break 'l220;
                    }
                }
            }
            if !ok250 {
                if ptmin >= 0.0 {
                    let s = format!(
                        "{d}\n WARNING - Cannot Modify Coefficients to Make\n Distribution More Than Minimum={}\n However, Modified Distribution is Positive with Minimum={}\n Will Output Modified Coefficients or Angular Distribution\n{d}",
                        f(self.sigmin, 11, 6),
                        f(ptmin, 11, 6)
                    );
                    self.lst(&s);
                } else {
                    self.lst(format!("{d}\n WARNING - Cannot Modify Coefficients to Make\n Distribution Positive\n Will Output Unmodified Coefficients or Angular Distribution\n{d}"));
                    for k in 1..=self.leguse as usize {
                        self.fl[k] = self.flold[k];
                        self.newfl[k] = 0;
                    }
                    self.anglen();
                }
            }
            // label 250
            let mut ii = 0;
            for l in 1..=self.leguse.max(0) as usize {
                if self.newfl[l] <= 0 {
                    continue;
                }
                if ii == 0 {
                    self.lst(" Coefficients Modified to Eliminate Negative ValueS\n Index   Old Value   New Value    per-cent Change");
                }
                ii = 1;
                let pc = 100.0 * (self.fl[l] - self.flold[l]) / self.flold[l];
                let s = format!("{}{}{}{}", i(l as i32, 6), f(self.flold[l], 12, 6), f(self.fl[l], 12, 6), f(pc, 12, 6));
                self.lst(&s);
            }
            if ii != 0 {
                self.lst(&d);
            }
        }
        // label 270
        if self.imlegs > 1 {
            self.list4();
        } else {
            self.tab4();
        }
        if self.leguse > 0 {
            for k in 1..=self.leguse as usize {
                if self.fl[k].abs() > 1.0 {
                    let s = format!("{d}\n WARNING - Legendre Coefficient={}{} (Expect -1 TO 1)", i(k as i32, 4), f(self.fl[k], 11, 6));
                    self.lst(&s);
                }
            }
        }
    }

    fn angtab(&mut self) {
        let d = dash();
        let terp = |x: f64, x1: f64, x2: f64, y1: f64, y2: f64| ((x - x1) * y2 + (x2 - x) * y1) / (x2 - x1);
        let n1 = self.io.ld.n1;
        let n2 = self.io.ld.n2;
        self.n2in = n2.max(0) as usize;
        let n1u = n1.max(0) as usize;
        let mut nbt = self.nbt;
        let mut int = self.int;
        self.io.terpi(&mut nbt[1..], &mut int[1..], n1u.min(100));
        self.nbt = nbt;
        self.int = int;
        self.io.terptest(&self.nbt[1..], &self.int[1..], n1, n2, 1);
        if self.n2in > MAXPOINT {
            let s = format!(
                "{}{}{} Angular Distribution Contains Over {} Points\n{} N2 ={} (Execution Terminated)",
                i(self.izalst, 7),
                i(self.matlst, 5),
                i(self.mtlst, 4),
                i(MAXPOINT as i32, 6),
                " ".repeat(16),
                i(self.n2in as i32, 6)
            );
            self.both(&s);
            self.io.enderror();
        }
        self.io.elast = -1.0;
        let mut x = std::mem::take(&mut self.xmuin);
        let mut y = std::mem::take(&mut self.ptin);
        self.io.pointi(&mut x[1..], &mut y[1..], self.n2in);
        self.xmuin = x;
        self.ptin = y;
        self.erruse = self.errok;
        self.loope = 0;
        'l30: loop {
            self.n2out = 1;
            self.xmuout[1] = self.xmuin[1];
            self.ptout[1] = self.ptin[1];
            let mut nr2 = 1i64;
            for ir in 1..=n1u {
                let nr1 = nr2 + 1;
                nr2 = self.nbt[ir] as i64;
                let intype = self.int[ir];
                if intype == 2 {
                    let mut npt = nr1;
                    while npt <= nr2 {
                        self.n2out += 1;
                        if self.n2out > self.nptmax {
                            self.erruse *= 2.0;
                            self.loope += 1;
                            continue 'l30;
                        }
                        let n = self.n2out;
                        self.xmuout[n] = self.xmuin[npt as usize];
                        self.ptout[n] = self.ptin[npt as usize];
                        npt += 1;
                    }
                    continue;
                }
                if intype < 2 {
                    let mut npt = nr1;
                    while npt <= nr2 {
                        let km2 = self.n2out;
                        let km1 = self.n2out + 1;
                        self.n2out += 2;
                        if self.n2out > self.nptmax {
                            self.erruse *= 2.0;
                            self.loope += 1;
                            continue 'l30;
                        }
                        let n = self.n2out;
                        self.xmuout[n] = self.xmuin[npt as usize];
                        self.ptout[n] = self.ptin[npt as usize];
                        self.xmuout[km1] = self.xmuout[n];
                        self.ptout[km1] = self.ptout[km2];
                        npt += 1;
                    }
                    continue;
                }
                let mut npt = nr1;
                while npt <= nr2 {
                    let xend = self.xmuin[npt as usize];
                    let yend = self.ptin[npt as usize];
                    loop {
                        // label 80
                        let km2 = self.n2out;
                        let km1 = self.n2out + 1;
                        self.n2out += 2;
                        if self.n2out > self.nptmax {
                            self.erruse *= 2.0;
                            self.loope += 1;
                            continue 'l30;
                        }
                        let n = self.n2out;
                        self.xmuout[n] = xend;
                        self.ptout[n] = yend;
                        if self.xmuout[n] <= self.xmuout[km2] {
                            self.xmuout[km1] = self.xmuout[n];
                            self.ptout[km1] = self.ptout[n];
                            self.n2out = km1;
                        } else {
                            loop {
                                // label 90
                                let n = self.n2out;
                                let (xnkm1, ynkm1);
                                match intype {
                                    1 | 2 => break,
                                    3 => {
                                        xnkm1 = (self.xmuout[km2] * self.xmuout[n]).sqrt();
                                        ynkm1 = 0.5 * (self.ptout[km2] + self.ptout[n]);
                                    }
                                    4 => {
                                        xnkm1 = 0.5 * (self.xmuout[km2] + self.xmuout[n]);
                                        ynkm1 = (self.ptout[km2] * self.ptout[n]).sqrt();
                                    }
                                    5 => {
                                        xnkm1 = (self.xmuout[km2] * self.xmuout[n]).sqrt();
                                        ynkm1 = (self.ptout[km2] * self.ptout[n]).sqrt();
                                    }
                                    _ => {
                                        // computed GO TO out of range: label 100
                                        self.xmuout[km1] = self.xmuout[n];
                                        self.ptout[km1] = self.ptout[n];
                                        self.n2out = km1;
                                        break;
                                    }
                                }
                                let yaprox = terp(xnkm1, self.xmuout[km2], self.xmuout[n], self.ptout[km2], self.ptout[n]);
                                if (ynkm1 - yaprox).abs() <= (self.erruse * ynkm1).abs() {
                                    self.xmuout[km1] = xnkm1;
                                    self.ptout[km1] = ynkm1;
                                    break;
                                }
                                self.xmuout[n] = xnkm1;
                                self.ptout[n] = ynkm1;
                            }
                        }
                        // label 160
                        if self.xmuout[self.n2out] < xend {
                            continue;
                        }
                        break;
                    }
                    npt += 1;
                }
            }
            break;
        }
        if self.imtabs == 2 {
            self.thin();
        }
        let sumx = Self::xint(&self.xmuout, &self.ptout, self.n2out);
        out9(self.e, &mut self.field6[1]);
        let s = L::new()
            .s(&format!("{}{}{}", i(self.izalst, 7), i(self.matlst, 5), i(self.mtlst, 4)))
            .b(&self.field6[1])
            .s(&" ".repeat(11))
            .s(&i(self.n2in as i32, 11))
            .s(&i(self.n2out as i32, 11))
            .s(&f(sumx, 11, 6));
        self.both(&s);
        let (xmumin, ptmin, imneg) = self.postiv();
        if imneg > 0 {
            let s = format!(
                "{d}\n Negative Angular Distribution at{} Points.\n Mimimum at Cosine={} Value={}\n{d}",
                i(imneg, 5),
                f(xmumin, 11, 6),
                f(ptmin, 11, 7)
            );
            self.lst(&s);
            self.typed[5] += 1;
            if self.imfix == 2 {
                let jtab = self.find_jtab();
                self.sigmin = self.sigtab[jtab];
                self.lowtab(ptmin, 2);
            }
        }
        if (sumx - 1.0).abs() > 0.05 {
            self.lst(format!("{d}\n WARNING - Angular Distribution NOT Normalized to within 5%\n{d}"));
        }
        if !(self.e > self.elast4) {
            out9(self.elast4, &mut self.field6[0]);
            let s = L::new()
                .s(&d)
                .s("\n WARNING - Energy is NOT in Ascending Order\n Last E =")
                .b(&self.field6[0])
                .s(" eV, Current E =")
                .b(&self.field6[1])
                .s(" eV\n")
                .s(&d);
            self.lst(&s);
            self.typed[6] += 1;
        }
        self.elast4 = self.e;
        self.tab4();
    }

    fn file4(&mut self) {
        let d = dash();
        let ntape1 = self.io.itape;
        let ntape2 = self.io.otape;
        let iza = fint(self.io.h.c1);
        if self.io.h.mat == self.matlst {
            if self.io.h.mt != self.mtlst {
                self.lst(&d);
            }
        } else if self.matlst > 0 {
            self.lst(format!(" {}", "=".repeat(78)));
        }
        self.izalst = iza;
        self.matlst = self.io.h.mat;
        self.mtlst = self.io.h.mt;
        if self.io.h.l2 == 3 {
            self.io.h.l2 = 4;
        }
        let mut ltt = self.io.h.l2;
        let lvt = self.io.h.l1;
        let idl = |s: &Lg| format!("{}{}{}", i(s.izalst, 7), i(s.matlst, 5), i(s.mtlst, 4));
        #[derive(PartialEq)]
        enum G {
            L110,
            L120,
            L130,
        }
        let mut n1combo = 0;
        let mut n2combo = 0;
        let g: G;
        if !(0..=4).contains(&ltt) {
            let s = format!("{}{}{} LTT={} Expect 0 to 4...Section Copied.", i(iza, 7), i(self.io.h.mat, 5), i(self.io.h.mt, 4), i(ltt, 6));
            self.lst(&s);
            g = G::L110;
        } else {
            match ltt {
                0 => {
                    self.typed[1] += 1;
                    let s = format!("{} Isotropic...Copied", idl(self));
                    self.lst(&s);
                    g = G::L110;
                }
                1 => {
                    if self.imlegs > 0 {
                        g = G::L120;
                    } else {
                        let s = format!("{} Legendre Coefficients...Copied", idl(self));
                        self.lst(&s);
                        self.typed[2] += 1;
                        g = G::L110;
                    }
                }
                2 => {
                    if self.imtabs > 0 {
                        g = G::L130;
                    } else {
                        let s = format!("{} Legendre Coefficients...Copied", idl(self));
                        self.lst(&s);
                        self.typed[3] += 1;
                        g = G::L110;
                    }
                }
                _ => {
                    if self.imlegs == 1 && self.imtabs > 0 {
                        // label 80: copy the section to LEGEND.TMP
                        self.io.u.rewind(ISCR);
                        self.io.otape = ISCR;
                        let (n1, n2) = self.cardi4();
                        self.io.ld.n1 = n1;
                        self.io.ld.n2 = n2;
                        self.cardo_ld(n1, n2);
                        let (n1l, n2l) = self.cardi4();
                        self.cardo_ld(n1l, n2l);
                        let n1lu = n1l.max(0) as usize;
                        let mut nbtc = self.nbtc;
                        let mut intc = self.intc;
                        self.io.terpi(&mut nbtc[1..], &mut intc[1..], n1lu.min(200));
                        self.nbtc = nbtc;
                        self.intc = intc;
                        self.io.terptest(&self.nbtc[1..], &self.intc[1..], n1l, n2l, 1);
                        self.io.terpo(&self.nbtc[1..], &self.intc[1..], n1lu);
                        for _ in 0..n2l.max(0) {
                            let (n1, n2) = self.cardi4();
                            self.io.ld.n1 = n1;
                            self.io.ld.n2 = n2;
                            self.cardo_ld(n1, n2);
                            let mut fl = std::mem::take(&mut self.fl);
                            self.io.listi(&mut fl[1..], n1.max(0) as usize);
                            self.io.listo(&fl[1..], n1.max(0) as usize);
                            self.fl = fl;
                        }
                        let (n1t, n2t) = self.cardi4();
                        self.cardo_ld(n1t, n2t);
                        let n1tu = n1t.max(0) as usize;
                        let mut nbtc = self.nbtc;
                        let mut intc = self.intc;
                        self.io.terpi(&mut nbtc[n1lu + 1..], &mut intc[n1lu + 1..], n1tu.min(200 - n1lu));
                        self.nbtc = nbtc;
                        self.intc = intc;
                        self.io.terptest(&self.nbtc[n1lu + 1..], &self.intc[n1lu + 1..], n1t, n2t, 1);
                        self.io.terpo(&self.nbtc[n1lu + 1..], &self.intc[n1lu + 1..], n1tu);
                        for _ in 0..n2t.max(0) {
                            let (n1, n2) = self.cardi4();
                            self.io.ld.n1 = n1;
                            self.io.ld.n2 = n2;
                            self.cardo_ld(n1, n2);
                            let n1u = n1.max(0) as usize;
                            let mut nbt = self.nbt;
                            let mut int = self.int;
                            self.io.terpi(&mut nbt[1..], &mut int[1..], n1u.min(100));
                            self.nbt = nbt;
                            self.int = int;
                            self.io.terptest(&self.nbt[1..], &self.int[1..], n1, n2, 1);
                            self.io.terpo(&self.nbt[1..], &self.int[1..], n1u);
                            self.io.elast = -1.0;
                            let n2u = n2.max(0) as usize;
                            let mut x = std::mem::take(&mut self.xmuin);
                            let mut y = std::mem::take(&mut self.ptin);
                            self.io.pointi(&mut x[1..], &mut y[1..], n2u);
                            self.io.pointo(&x[1..], &y[1..], n2u);
                            self.xmuin = x;
                            self.ptin = y;
                        }
                        n1combo = n1l + n1t;
                        n2combo = n2l + n2t;
                        for k in (n1lu + 1)..=(n1combo.max(0) as usize) {
                            self.nbtc[k] += self.nbtc[n1lu];
                        }
                        self.io.u.endfile(ISCR);
                        self.io.u.rewind(ISCR);
                        self.io.itape = ISCR;
                        self.io.otape = ntape2;
                        g = G::L120;
                    } else {
                        let s = format!("{} Legendre Coefficients + Tabulated...Copied", idl(self));
                        self.lst(&s);
                        self.typed[4] += 1;
                        g = G::L110;
                    }
                }
            }
        }
        if g == G::L110 {
            self.contol();
            self.io.copys();
            return;
        }
        if g == G::L120 && self.imlegs == 1 {
            self.io.h.l2 = 2;
        }
        self.contol();
        let (n1, n2) = self.cardi4();
        self.io.ld.n1 = n1;
        self.io.ld.n2 = n2;
        self.cardo_ld(n1, n2);
        if lvt == 1 {
            let mut l = 1;
            while l <= n1 {
                self.io.copy1();
                l += 6;
            }
        }
        let loopltt = if ltt == 4 { 2 } else { 1 };
        for koop in 1..=loopltt {
            if loopltt == 2 && koop == 1 {
                ltt = 1;
            }
            if loopltt == 2 && koop == 2 {
                ltt = 2;
            }
            let (n1, n2) = self.cardi4();
            self.io.ld.n1 = n1;
            self.io.ld.n2 = n2;
            let n1u = n1.max(0) as usize;
            let mut nbt = self.nbt;
            let mut int = self.int;
            self.io.terpi(&mut nbt[1..], &mut int[1..], n1u.min(100));
            self.nbt = nbt;
            self.int = int;
            self.io.terptest(&self.nbt[1..], &self.int[1..], n1, n2, 1);
            if koop != 2 {
                if loopltt == 1 {
                    self.cardo_ld(n1, n2);
                    self.io.terpo(&self.nbt[1..], &self.int[1..], n1u);
                } else {
                    self.cardo_ld(n1combo, n2combo);
                    self.io.terpo(&self.nbtc[1..], &self.intc[1..], n1combo.max(0) as usize);
                }
                if ltt == 1 && (1..=n1u).any(|l| self.int[l] >= 4) {
                    // format 330 with no list: an empty line
                    self.lst("");
                }
            }
            self.ne = n2;
            if loopltt == 1 {
                self.typed[(ltt + 1) as usize] += self.ne;
            } else {
                self.typed[4] += self.ne;
            }
            self.elast4 = 0.0;
            if ltt == 1 {
                self.n2out = usize::MAX;
            }
            if ltt == 2 {
                self.leguse = -1;
            }
            let ne = self.ne;
            self.je = 1;
            while self.je <= ne {
                let (c1, c2, l1, l2, n1, n2) = self.io.cardi_ld();
                let _ = (c1, l1, l2, n1, n2);
                self.e = c2;
                if ltt == 1 {
                    self.angleg();
                } else {
                    self.angtab();
                    if self.loope > 0 {
                        out9(self.errok, &mut self.field6[0]);
                        let s = L::new().s(&d).s("\n WARNING...Allowable Error=").b(&self.field6[0]).s("\n").s(&d);
                        self.lst(&s);
                    }
                }
                let mut coslst = -1.0;
                let n2o = self.n2out;
                for j in 1..=n2o {
                    if self.xmuout[j] < coslst {
                        let s = format!(
                            "{d}\n WARNING...Cosines NOT in Ascending Order{}{}{}\n{d}",
                            i(j as i32, 6),
                            f(self.xmuout[j], 16, 8),
                            f(coslst, 16, 8)
                        );
                        self.lst(&s);
                    } else {
                        coslst = self.xmuout[j];
                    }
                }
                if (self.xmuout[1] + 1.0).abs() <= 0.001 {
                    self.xmuout[1] = -1.0;
                }
                if (self.xmuout[n2o] - 1.0).abs() <= 0.001 {
                    self.xmuout[n2o] = 1.0;
                }
                if (self.xmuout[1] + 1.0).abs() > 0.001 || (self.xmuout[n2o] - 1.0).abs() > 0.001 {
                    let s = format!(
                        "{d}\n WARNING...Cosines Do NOT Span -1 to +1Range={}{}\n{d}",
                        f(self.xmuout[1], 16, 8),
                        f(self.xmuout[n2o], 16, 8)
                    );
                    self.lst(&s);
                }
                self.je += 1;
            }
        }
        let (m, mf, mt) = (self.io.h.mat, self.io.h.mf, self.io.h.mt);
        self.io.maxie1(m, mf, mt);
        let e = self.e;
        self.io.maxie2(e, 1.0);
        self.io.itape = ntape1;
        self.io.otape = ntape2;
        self.io.copys();
    }
}

/// `INTOUT` of LEGEND (stops at the first zero digit value).
fn intout(v: i32, field: &mut [u8]) {
    for c in field.iter_mut() {
        *c = b' ';
    }
    let mut ii = v;
    for k in (0..field.len()).rev() {
        if ii <= 0 {
            break;
        }
        let kk = ii / 10;
        let ll = ii - 10 * kk;
        field[k] = b'0' + ll as u8;
        ii = kk;
    }
}

pub fn run() -> ! {
    let mut io = Endf::new("LEGEND");
    io.timer();
    io.inp = 2;
    io.outp = 3;
    io.itape = 10;
    io.otape = 11;
    io.u.open(3, "LEGEND.LST");
    io.u.open(ISCR, "LEGEND.TMP");
    let mut istat1 = if io.u.open_old(2, "LEGEND.INP") { 0 } else { 1 };
    let mut s = Lg {
        io,
        etab1: vec![0.0; MAXMAT + 1],
        etab2: vec![0.0; MAXMAT + 1],
        sigtab: vec![0.0; MAXMAT + 1],
        dfltab: vec![0.0; MAXMAT + 1],
        matab1: vec![0; MAXMAT + 1],
        mftab1: vec![0; MAXMAT + 1],
        mtab1: vec![0; MAXMAT + 1],
        matab2: vec![0; MAXMAT + 1],
        mftab2: vec![0; MAXMAT + 1],
        mtab2: vec![0; MAXMAT + 1],
        itab: 0,
        flmult: vec![0.0; MAXLEGS + 1],
        pcon1: vec![0.0; MAXLEGS + 1],
        pcon2: vec![0.0; MAXLEGS + 1],
        values: vec![0.0; MAXLEGS + 1],
        dfl: vec![0.0; MAXLEGS + 1],
        sumn: vec![0.0; MAXLEGS + 1],
        pn: vec![0.0; MAXLEGS + 1],
        fl: vec![0.0; MAXPOINT + 2],
        flold: vec![0.0; MAXLEGS + 1],
        flkeep: vec![[0.0; 4]; MAXLEGS + 1],
        xmubase: vec![0.0; 5 * MAXLEGS + 1],
        ekeep: [0.0; 4],
        elast4: 0.0,
        sum: 0.0,
        dflchn: 0.0,
        dflmax: 0.0,
        lvalue: vec![0; MAXLEGS + 1],
        newfl: vec![0; MAXLEGS + 1],
        lkeep: [0; 4],
        legin: 0,
        leguse: 0,
        kpass: 0,
        nbt: [0; 101],
        int: [0; 101],
        nbtc: [0; 201],
        intc: [0; 201],
        xmuin: vec![0.0; MAXPOINT + 3],
        ptin: vec![0.0; MAXPOINT + 3],
        xmuout: vec![0.0; MAXPOINT + 3],
        ptout: vec![0.0; MAXPOINT + 3],
        errmax: 0.0,
        errok: 0.0,
        erruse: 0.0,
        loope: 0,
        nptmax: MAXPOINT,
        legmax: MAXLEGS as i32,
        e: 0.0,
        izalst: -999999,
        matlst: -9999,
        mtlst: -9999,
        ne: 0,
        je: 0,
        typed: [0; 7],
        imtabs: 0,
        imlegs: 0,
        imfix: 0,
        imvar: 0,
        n2in: 0,
        n2out: 0,
        sigmin: 0.0,
        field6: [[0u8; 11]; 6],
        leglst: 0,
        th_dsign: 0.0,
        th_slpmax: 0.0,
        th_slpmin: 0.0,
        th_ithin: 0,
        contol_lstmat: -9999,
        contol_lstmf: -999,
    };
    s.both(format!(" Linearized ENDF/B Angular Distributions (LEGEND 2026-1)\n{}", dash()));
    s.flmult[1] = 1.5;
    let (mut con1, mut con2) = (3.0f64, 1.0f64);
    for l in 2..=MAXLEGS {
        s.flmult[l] = s.flmult[l - 1] + 1.0;
        s.pcon1[l] = con1 / (con2 + 1.0);
        s.pcon2[l] = -con2 / (con2 + 1.0);
        con1 += 2.0;
        con2 += 1.0;
    }
    let mut nptmu = MAXPOINT as i32;
    let mut have = false;
    if istat1 != 1 {
        if let Rd::Line(l) = s.io.u.read(2) {
            let fld = fio::rd_a(&l, 0, 11);
            let v = [rdi(&l, 11, 11), rdi(&l, 22, 11), rdi(&l, 33, 11), rdi(&l, 44, 11), rdi(&l, 55, 11)];
            if v.iter().all(|x| x.is_some()) {
                s.field6[0].copy_from_slice(&fld);
                s.errok = in9(&s.field6[0]);
                nptmu = v[0].unwrap();
                s.imtabs = v[1].unwrap();
                s.imlegs = v[2].unwrap();
                s.imfix = v[3].unwrap();
                s.imvar = v[4].unwrap();
                have = true;
            }
        }
    }
    if !have {
        istat1 = 1;
        s.errok = 1.0e-03;
        nptmu = MAXPOINT as i32;
        s.imtabs = 1;
        s.imlegs = 1;
        s.imfix = 0;
        s.imvar = 0;
    }
    if s.errok <= 0.0 {
        s.errok = 1.0e-03;
    }
    if nptmu <= 0 {
        nptmu = s.nptmax as i32;
    }
    if nptmu < 11 {
        nptmu = 11;
    }
    if nptmu > s.nptmax as i32 {
        nptmu = s.nptmax as i32;
    }
    s.nptmax = nptmu as usize;
    if !(0..=2).contains(&s.imtabs) {
        s.imtabs = 0;
    }
    if !(0..=2).contains(&s.imlegs) {
        s.imlegs = 0;
    }
    if !(0..=2).contains(&s.imfix) {
        s.imfix = 0;
    }
    if s.imvar != 1 {
        s.imvar = 0;
    }
    s.errmax = 100.0 * s.errok;
    out9(s.errok, &mut s.field6[0]);
    let l = L::new()
        .s(" Allowable Accuracy----------")
        .b(&s.field6[0])
        .s(" (")
        .s(&f(s.errmax, 9, 4))
        .s(" per-cent)\n Maximum Number of Cosines---")
        .s(&i(s.nptmax as i32, 11));
    s.lst(&l);
    let tabs = [
        " Tabulated Data--------------Copy",
        " Tabulated Data--------------Linearize      - Output Tables",
        " Tabulated Data--------------Linearize/Thin - Output Tables",
    ];
    let legs = [
        " Legendre Data---------------Copy",
        " Legendre Data---------------Tabulate       - Output Tables",
        " Legendre Data---------------Tabulate       - Output Coefficients",
    ];
    let fixs = [
        " Negative Data Correction----None",
        " Negative Table Correction---None\n Legendre Coef. Correction---Modify up to 1 per-cent",
        " Negative Data Correction----Force to be Positive",
    ];
    s.lst(tabs[s.imtabs as usize]);
    s.lst(legs[s.imlegs as usize]);
    s.lst(fixs[s.imfix as usize]);
    s.lst(if s.imvar == 0 { " Coefficient Tests-----------Off" } else { " Coefficient Tests-----------On" });
    let mut namein = None;
    let mut nameout = None;
    if istat1 != 1 {
        if let Rd::Line(l) = s.io.u.read(2) {
            let mut n = fio::rd_a(&l, 0, 60);
            if blank(&n) {
                n = fio::rd_a(b"ENDFB.IN", 0, 60);
            }
            namein = Some(n);
            if let Rd::Line(l) = s.io.u.read(2) {
                let mut n = fio::rd_a(&l, 0, 60);
                if blank(&n) {
                    n = fio::rd_a(b"ENDFB.OUT", 0, 60);
                }
                nameout = Some(n);
            }
        }
    }
    let namein = namein.unwrap_or_else(|| fio::rd_a(b"ENDFB.IN", 0, 60));
    let nameout = match nameout {
        Some(n) => n,
        None => {
            istat1 = 1;
            fio::rd_a(b"ENDFB.OUT", 0, 60)
        }
    };
    let l = L::new().s(&dash()).s("\n ENDF/B Input and Output Data Filenames\n ").b(&namein).s("\n ").b(&nameout);
    s.both(&l);
    let ni = String::from_utf8_lossy(&namein).into_owned();
    if !s.io.u.open_old(10, &ni) {
        let l = L::new().s("\n\n ERROR - Opening ENDF/B formatted file\n ").b(&namein).s("\n\n");
        s.both(&l);
        s.io.enderror();
    }
    let no = String::from_utf8_lossy(&nameout).into_owned();
    s.io.u.open(11, &no);
    s.lst(format!(
        "{d}\n Minimum Allowable Cross Section and Maximum Allowable Coefficient Change\n{d}\n  MAT1 MF1  MT1  MAT2 MF2  MT2         E1         E2  Sigma-Low   Delta-Fl\n{}(eV)       (eV)     (b/mu)\n{d}",
        " ".repeat(38),
        d = dash()
    ));
    let mut iend = 0;
    let mut itab = 1usize;
    let mut done = false;
    while itab <= MAXMAT {
        let mut go120 = istat1 == 1;
        if !go120 {
            let rec = match s.io.u.read(2) {
                Rd::Line(l) => {
                    let v = [rdi(&l, 0, 6), rdi(&l, 6, 2), rdi(&l, 8, 3), rdi(&l, 11, 6), rdi(&l, 17, 2), rdi(&l, 19, 3)];
                    if v.iter().all(|x| x.is_some()) {
                        Some((v.map(|x| x.unwrap()), fio::rd_a(&l, 22, 44)))
                    } else {
                        None
                    }
                }
                Rd::Eof => None,
            };
            match rec {
                None => go120 = true,
                Some((v, flds)) => {
                    s.matab1[itab] = v[0];
                    s.mftab1[itab] = v[1];
                    s.mtab1[itab] = v[2];
                    s.matab2[itab] = v[3];
                    s.mftab2[itab] = v[4];
                    s.mtab2[itab] = v[5];
                    for k in 0..4 {
                        s.field6[k].copy_from_slice(&flds[11 * k..11 * k + 11]);
                    }
                    s.etab1[itab] = in9(&s.field6[0]);
                    s.etab2[itab] = in9(&s.field6[1]);
                    s.sigtab[itab] = in9(&s.field6[2]);
                    s.dfltab[itab] = in9(&s.field6[3]);
                    let any = v.iter().any(|&x| x != 0)
                        || s.etab1[itab] != 0.0
                        || s.etab2[itab] != 0.0
                        || s.sigtab[itab] != 0.0
                        || s.dfltab[itab] != 0.0;
                    if any {
                        if s.matab2[itab] < s.matab1[itab] {
                            s.matab2[itab] = s.matab1[itab];
                        }
                        if s.mftab2[itab] < s.mftab1[itab] {
                            s.mftab2[itab] = s.mftab1[itab];
                        }
                        if s.mtab2[itab] < s.mtab1[itab] {
                            s.mtab2[itab] = s.mtab1[itab];
                        }
                        if s.etab1[itab] <= 0.0 {
                            s.etab1[itab] = 0.0;
                        }
                        if s.etab2[itab] <= 0.0 {
                            s.etab2[itab] = 1.0e9;
                        }
                        if s.sigtab[itab] <= 0.0 {
                            s.sigtab[itab] = 0.001;
                        }
                        if s.dfltab[itab] <= 0.01 {
                            s.dfltab[itab] = 0.01;
                        }
                    } else if itab > 1 {
                        itab -= 1;
                        done = true;
                        break;
                    } else {
                        go120 = true;
                    }
                }
            }
        }
        if go120 {
            istat1 = 1;
            iend = 1;
            if itab > 1 {
                itab -= 1;
                done = true;
                break;
            }
            s.matab1[itab] = 1;
            s.matab2[itab] = 9999;
            s.mftab1[itab] = 1;
            s.mftab2[itab] = 99;
            s.mtab1[itab] = 1;
            s.mtab2[itab] = 999;
            s.etab1[itab] = 0.0;
            s.etab2[itab] = 1.0e9;
            s.sigtab[itab] = 0.001;
            s.dfltab[itab] = 0.01;
        }
        // label 140
        out9(s.etab1[itab], &mut s.field6[0]);
        out9(s.etab2[itab], &mut s.field6[1]);
        out9(s.sigtab[itab], &mut s.field6[2]);
        let l = L::new()
            .s(&format!(
                "{}{}{}{}{}{}",
                i(s.matab1[itab], 6),
                i(s.mftab1[itab], 4),
                i(s.mtab1[itab], 5),
                i(s.matab2[itab], 6),
                i(s.mftab2[itab], 4),
                i(s.mtab2[itab], 5)
            ))
            .b(&s.field6[0])
            .b(&s.field6[1])
            .b(&s.field6[2])
            .s(&f(s.dfltab[itab], 11, 5));
        s.lst(&l);
        if iend > 0 {
            done = true;
            break;
        }
        itab += 1;
    }
    if !done {
        itab = MAXMAT;
    }
    s.itab = itab;
    s.typed = [0; 7];
    s.io.copyl();
    let l = L::new().s(&dash()).s("\n Tape Label\n").s(&dash()).s("\n ").b(&s.io.card).s(&i(s.io.h.mat, 4));
    s.both(&l);
    let mut lastmat = -10000;
    let hdr = format!(
        "{d}\n     ZA  MAT  MT  Energy-eV   Legendre     Points     Points\n{}Order           In        Out   Integral Messages\n{}Read (Used)\n{d}",
        " ".repeat(31),
        " ".repeat(28),
        d = dash()
    );
    loop {
        s.contil();
        if s.io.h.mat < 0 {
            if lastmat > 0 {
                s.io.maxie3(1);
                s.io.timemat();
            }
            break;
        }
        if s.io.h.mat != lastmat {
            if lastmat > 0 {
                s.io.maxie3(1);
                s.io.timemat();
            }
            s.io.maxie0();
            s.both(&hdr);
            lastmat = s.io.h.mat;
        }
        if s.io.h.mf == 1 && s.io.h.mt == 451 {
            s.contol();
            s.file1();
            s.io.copys();
            continue;
        }
        if s.io.h.mf == 4 {
            s.file4();
            continue;
        }
        s.contol();
        s.io.copys();
    }
    s.contol();
    let t = s.typed;
    let l = format!(
        "{d}\n Summary of Data as Read\n{d}\n Isotropic Sections----------------------{}\n Energies with Legendre Coefficients-----{}\n Energies with Tabulations---------------{}\n Energies with Combined------------------{}\n Energies with NEGATIVE distributions----{}\n Repeated Energies-----------------------{}",
        i(t[1], 8),
        i(t[2], 8),
        i(t[3], 8),
        i(t[4], 8),
        i(t[5], 8),
        i(t[6], 8),
        d = dash()
    );
    s.both(&l);
    s.io.maxie4(1);
    let _ = f9s;
    s.io.endit();
}
