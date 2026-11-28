//! GROUPR angular and energy distributions (getdis, getfle, getaed, aedi,
//! getgfl, getco, getgyl, getsed, anased, hnab).

use super::*;
use crate::fmath::{cos, ipow, powi, sin};
use crate::mathm::{erfc, legndr};

/// Module variables that getco reads.
pub(super) struct CoPar {
    pub awr: f64,
    pub awrp: f64,
    pub izap: i32,
    pub izat: i32,
}

impl Groupr {
    fn copar(&self) -> CoPar {
        CoPar { awr: self.awr, awrp: self.awrp, izap: self.izap, izat: self.izat }
    }

    /// Feed function for elastic or discrete inelastic scattering using
    /// data from File 4 or File 6.
    pub(super) fn getdis(
        &mut self,
        c: &mut Ctx,
        e: f64,
        enext: &mut f64,
        idisc: &mut i32,
        yld: &mut f64,
        ff: &mut V2,
        nl: &mut i32,
        ng: &mut i32,
        iglo: &mut i32,
        nq: &mut i32,
        nin: i32,
        _nlg: i32,
    ) {
        const QP4: [f64; 4] = [-0.8611363116e0, -0.3399810436e0, 0.3399810436e0, 0.8611363116e0];
        const QW4: [f64; 4] = [0.3478548451e0, 0.6521451549e0, 0.6521451549e0, 0.3478548451e0];
        const QP8: [f64; 8] = [
            -0.9602898565e0,
            -0.7966664774e0,
            -0.5255324099e0,
            -0.1834346425e0,
            0.1834346425e0,
            0.5255324099e0,
            0.7966664774e0,
            0.9602898565e0,
        ];
        const QW8: [f64; 8] = [
            0.1012285362e0,
            0.2223810345e0,
            0.3137066459e0,
            0.3626837834e0,
            0.3626837834e0,
            0.3137066459e0,
            0.2223810345e0,
            0.1012285363e0,
        ];
        const QP12: [f64; 12] = [
            -0.9815606342e0,
            -0.9041172564e0,
            -0.7699026742e0,
            -0.5873179543e0,
            -0.3678314990e0,
            -0.1252334085e0,
            0.1252334085e0,
            0.3678314990e0,
            0.5873179543e0,
            0.7699026742e0,
            0.9041172564e0,
            0.9815606342e0,
        ];
        const QW12: [f64; 12] = [
            0.0471753364e0,
            0.1069393260e0,
            0.1600783285e0,
            0.2031674267e0,
            0.2334925365e0,
            0.2491470458e0,
            0.2491470458e0,
            0.2334925365e0,
            0.2031674267e0,
            0.1600783285e0,
            0.1069393260e0,
            0.0471753364e0,
        ];
        const QP20: [f64; 20] = [
            -0.9931286e0,
            -0.9639719e0,
            -0.9122344e0,
            -0.839117e0,
            -0.7463319e0,
            -0.6360537e0,
            -0.510867e0,
            -0.3737061e0,
            -0.2277859e0,
            -0.0765265e0,
            0.0765265e0,
            0.2277859e0,
            0.3737061e0,
            0.5108670e0,
            0.6360537e0,
            0.7463319e0,
            0.8391170e0,
            0.9122344e0,
            0.9639719e0,
            0.9931286e0,
        ];
        const QW20: [f64; 20] = [
            0.0176140e0,
            0.0406014e0,
            0.0626720e0,
            0.0832767e0,
            0.1019301e0,
            0.1181945e0,
            0.1316886e0,
            0.1420961e0,
            0.1491730e0,
            0.1527534e0,
            0.1527534e0,
            0.1491730e0,
            0.1420961e0,
            0.1316886e0,
            0.1181945e0,
            0.1019301e0,
            0.0832767e0,
            0.0626720e0,
            0.0406014e0,
            0.0176140e0,
        ];
        let fm = 1.0e-12f64;
        let emax = 1.0e10f64;
        let small = 1.0e-10f64;
        let wcut = 0.94e0f64;
        let shade = 0.99999995e0f64;
        let zero = 0.0f64;
        let ten = 10.0f64;
        let mut fle = [0.0f64; 22];
        let mut p = [0.0f64; 22];
        let mut flt = [0.0f64; 21];
        let mfd = self.mfd;
        let mtd = self.mtd;
        let ngn = self.ngn;

        let mut nld = 21;
        let ngn1 = ngn + 1;
        *enext = emax;
        *idisc = 0;
        let mut lidp = 0;
        if self.izap == self.izat {
            lidp = 1;
        }
        *yld = 1.0;
        let lr = self.lrflag;
        if lr == 16 || lr == 21 || lr == 24 || lr == 26 || lr == 30 {
            *yld = 2.0;
        }
        if lr == 17 || lr == 25 || lr == 38 {
            *yld = 3.0;
        }
        if lr == 37 {
            *yld = 4.0;
        }
        let lcd = 2;
        let mut mft = 4;
        if mfd == 8 {
            mft = 8;
        }
        if mfd >= 21 && mfd <= 26 {
            mft = 8;
        }
        if mfd == 3 && mtd == 2 && self.izap > 1 {
            mft = 8;
        }
        let mut mtt = mtd;
        if mtd == 251 || mtd == 252 || mtd == 253 {
            mtt = 2;
        }
        let mut en = 0.0f64;
        let mut idis = 0i32;
        let mut lcdm = lcd;
        self.getfle(c, e, &mut en, &mut idis, &mut F(&mut fle[1..]), &mut nld, &mut lcdm, self.matd, mft, mtt, nin);
        if en < *enext * (1.0 - small) {
            *idisc = idis;
        }
        if en < *enext * (1.0 - small) {
            *enext = en;
        }
        let awr = self.awr;
        let aprime = self.aprime;
        let awr2 = awr * (awr + 1.0 - aprime) / aprime;
        if e == 0.0 {
            // 600
            self.sv.di_iecl = 0;
            self.sv.di_iech = 0;
            self.sv.di_ecl = 0.0;
            self.sv.di_ech = 0.0;
            self.sv.di_ecn = 0.0;
        } else {
            let ld = nld - 1;
            let mut npo = ld + *nl + (log(300.0 / awr)) as i32;
            if mfd == 3 {
                npo += self.lord;
            }
            if self.izap > 1 {
                npo += 6;
            }
            self.sv.di_nqp0 = 4;
            if npo > 8 {
                self.sv.di_nqp0 = 8;
            }
            if npo > 12 {
                self.sv.di_nqp0 = 12;
            }
            if npo > 16 {
                self.sv.di_nqp0 = 20;
            }
            let mut ast = 0.0f64;
            if e > self.thresh {
                ast = sqrt(awr2) * sqrt(1.0 - self.thresh / e);
            }
            *ng = 0;
            let mut ig = 0;
            let mut ngnow = 1;
            let mut ignow = 0;
            *iglo = self.sv.di_iech - 1;
            if *iglo < 1 {
                *iglo = 1;
            }
            let mut wmin = -1.0f64;
            let mut wmax = 1.0f64;
            if self.izap > 1 && mtd == 2 {
                wmax = wcut;
                if lidp == 1 {
                    wmin = -wcut;
                }
            }
            let mut whi = wmin;
            let mut lbl = 410;
            loop {
                match lbl {
                    410 => {
                        *ng += 1;
                        if *ng > ngn1 {
                            *ng = ngn1;
                            ig = ignow;
                            for il in 1..=*nl {
                                ff[(il, *ng)] = 0.0;
                            }
                            lbl = 465;
                            continue;
                        }
                        ngnow = *ng;
                        for il in 1..=*nl {
                            ff[(il, *ng)] = 0.0;
                        }
                        if ast == 0.0 {
                            lbl = 465;
                            continue;
                        }
                        lbl = 425;
                    }
                    425 => {
                        let wlo = whi;
                        ig += 1;
                        if *iglo + ig > ngn1 {
                            *ng = ngnow;
                            ig -= 1;
                            lbl = 465;
                            continue;
                        }
                        ignow = ig;
                        let ep = self.egn[*iglo + ig];
                        if (mtd == 252 || mtd == 253) && ep >= e {
                            lbl = 465;
                            continue;
                        }
                        whi = ((ep / e) * powc(1.0 + awr, 2) / aprime - (1.0 + ast * ast)) / (2.0 * ast);
                        if whi < shade * wmin {
                            whi = wmin;
                        }
                        if whi > shade * wmax {
                            whi = wmax;
                        }
                        if (mtd == 251 || mtd == 252 || mtd == 253) && wlo == whi {
                            lbl = 425;
                            continue;
                        }
                        if mfd == 3 && wlo == whi {
                            lbl = 425;
                            continue;
                        }
                        if wlo == whi {
                            lbl = 410;
                            continue;
                        }
                        if (mtd == 252 || mtd == 253) && whi >= wmax {
                            // 464
                            ff[(1, 1)] = 1.0;
                            ig += 1;
                            lbl = 465;
                            continue;
                        }

                        //--integrate between wlo and whi using gaussian quadrature
                        let aa = (whi + wlo) / 2.0;
                        let b = (whi - wlo) / 2.0;
                        let mut nqp = (npo as f64 * 2.0 * b) as i32;
                        let mut wqp: f64;
                        if nqp > npo {
                            nqp = npo;
                        }
                        if nqp <= 8 {
                            nqp = 4;
                        }
                        if nqp > 8 && nqp <= 12 {
                            nqp = 8;
                        }
                        if nqp > 12 && nqp <= 16 {
                            nqp = 12;
                        }
                        if nqp > 16 {
                            nqp = 20;
                        }
                        for il in 1..=*nl as usize {
                            flt[il] = 0.0;
                        }
                        for iqp in 1..=nqp {
                            let k = (iqp - 1) as usize;
                            let wqw;
                            if nqp == 4 {
                                wqp = aa + b * QP4[k];
                                wqw = b * QW4[k];
                            } else if nqp == 8 {
                                wqp = aa + b * QP8[k];
                                wqw = b * QW8[k];
                            } else if nqp == 12 {
                                wqp = aa + b * QP12[k];
                                wqw = b * QW12[k];
                            } else if nqp == 20 {
                                wqp = aa + b * QP20[k];
                                wqw = b * QW20[k];
                            } else {
                                c.error("getdis", "illegal nqp", " ");
                            }
                            if self.law == 4 {
                                wqp = -wqp;
                            }

                            //--compute scattering probability at this cm cosine
                            legndr(wqp, &mut p, ld);
                            let mut prob = 0.0f64;
                            for il in 1..=nld {
                                prob = prob + fle[il as usize] * p[(il - 1) as usize] * (2 * il - 1) as f64 / 2.0;
                            }

                            //--add in coulomb part for charged particles
                            if mtd == 2 && self.izap > 1 {
                                let spi = self.spi;
                                let i2s = nint(2.0 * spi);
                                let ai = self.awrp * AMASSN;
                                let at = awr * AMASSN;
                                let zt = (self.izat / 1000) as f64;
                                let zi = (self.izap / 1000) as f64;
                                let cc1 = 2.0 * AMU * EV * powc(fm, 2) / powc(HBAR, 2);
                                let ee = (EV / 10000000.0) * (CLIGHT / 10.0);
                                let cc2 = powc(ee, 4) * AMU / (2.0 * powc(HBAR, 2) * EV);
                                let wn = at * sqrt(cc1 * e * ai) / (ai + at);
                                let eta = zt * zi * sqrt(cc2 * ai / e);
                                let mut sigc = 0.0f64;
                                if lidp == 0 {
                                    sigc = (powc(eta, 2) / powc(wn, 2)) / powc(1.0 - wqp, 2);
                                }
                                if lidp == 1 {
                                    sigc = ((2.0 * powc(eta, 2) / powc(wn, 2)) / (1.0 - powc(wqp, 2)))
                                        * ((1.0 + powc(wqp, 2)) / (1.0 - powc(wqp, 2))
                                            + (ipow(-1, i2s) as f64) * cos(eta * log((1.0 + wqp) / (1.0 - wqp)))
                                                / (2.0 * spi + 1.0));
                                }
                                if lidp == 0 {
                                    prob = prob / (1.0 - wqp);
                                }
                                if lidp == 1 {
                                    prob = prob / (1.0 - wqp * wqp);
                                }
                                prob = prob + sigc;
                            }

                            //--multiply scattering probability by legendre polynomials at
                            //--corresponding lab cosine and accumulate integrals
                            let wlab = (1.0 + ast * wqp) / sqrt(1.0 + ast * ast + 2.0 * ast * wqp);
                            let l = *nl - 1;
                            legndr(wlab, &mut p, l);
                            prob = *yld * wqw * prob;
                            if mtd == 251 {
                                prob = wlab * prob;
                            }
                            for il in 1..=*nl as usize {
                                flt[il] = flt[il] + prob * p[il - 1];
                            }
                        }
                        for il in 1..=*nl {
                            ff[(il, *ng)] = ff[(il, *ng)] + flt[il as usize];
                        }
                        if whi >= wmax {
                            lbl = 465;
                            continue;
                        }
                        if *iglo + ig > ngn {
                            lbl = 465;
                            continue;
                        }
                        if mtd == 251 || mtd == 252 || mtd == 253 {
                            lbl = 425;
                            continue;
                        }
                        if mfd == 3 {
                            lbl = 425;
                            continue;
                        }
                        lbl = 410;
                    }
                    _ => break,
                }
            }
            // 465
            //--reduce significant figures for small values
            let ndig = 7;
            let fact = powi(ten, ndig);
            for il in 1..=*nl {
                for ii in 1..=*ng {
                    let iii = nint(fact * ff[(il, ii)] + powi(ten, ndig - 11));
                    ff[(il, ii)] = iii as f64 / fact;
                }
            }
            if mfd == 3 && mtd == 2 {
                *iglo = 1;
            }
            if mtd == 251 || mtd == 252 || mtd == 253 {
                for il in 1..=*nl {
                    if mtd == 251 || mtd == 252 {
                        ff[(il, *ng + 1)] = 1.0;
                    }
                    if *iglo + ig - 1 >= 1 {
                        if mtd == 252 {
                            ff[(il, *ng)] = ff[(il, *ng)] * log(self.egn[*iglo + ig] / self.egn[*iglo + ig - 1]);
                        }
                        if mtd == 253 {
                            ff[(il, *ng + 1)] = ff[(il, *ng)];
                            ff[(il, *ng)] = ff[(il, *ng)] * log(e / self.egn[*iglo + ig - 1]);
                        }
                    }
                }
                *ng += 1;
                *iglo = 1;
            }
        }

        //--check next critical point for discrete channel scattering
        // 610
        let s = &mut self.sv;
        let thresh = self.thresh;
        let q = self.q;
        let alpha = self.alpha;
        if !(e < shade * s.di_ecn) {
            // 620
            while !(e < shade * s.di_ecl) {
                s.di_ecl = emax;
                if s.di_iecl == ngn + 1 {
                    break;
                }
                s.di_iecl += 1;
                if !(thresh > zero) {
                    let test = 1.0;
                    if !(aprime > test) {
                        s.di_ecl = self.egn[s.di_iecl];
                        continue;
                    }
                    // 625
                    s.di_ecl = self.egn[s.di_iecl] * powc(awr + 1.0, 2) / (4.0 * awr);
                    continue;
                }
                // 630
                let ef = (awr + 1.0) * self.egn[s.di_iecl] / (awr + 1.0 - aprime);
                let aa = 1.0 + ef / (-q);
                let mut disc = (awr2 * aa - 1.0) * ef / (-q);
                if disc < zero {
                    continue;
                }
                disc = sqrt(disc);
                let af = (1.0 - disc) / aa;
                if powc(af, 2) == awr2 {
                    s.di_ecl = emax;
                } else {
                    s.di_ecl = thresh / (1.0 - powc(af, 2) / awr2);
                }
            }
            // 640
            while !(e < shade * s.di_ech) {
                s.di_ech = emax;
                if s.di_iech == ngn + 1 {
                    break;
                }
                s.di_iech += 1;
                if !(thresh > zero) {
                    let mut test = 1.0;
                    if aprime > test {
                        // 645
                        s.di_ech = emax;
                        continue;
                    }
                    test = test / 10.0;
                    if alpha < test {
                        continue;
                    }
                    s.di_ech = self.egn[s.di_iech] / alpha;
                    continue;
                }
                // 650
                let ef = (awr + 1.0) * self.egn[s.di_iech] / (awr + 1.0 - aprime);
                let aa = 1.0 + ef / (-q);
                let mut disc = (awr2 * aa - 1.0) * ef / (-q);
                if disc < zero {
                    continue;
                }
                disc = sqrt(disc);
                let af = (1.0 + disc) / aa;
                if powc(af, 2) == awr2 {
                    s.di_ech = emax;
                } else {
                    s.di_ech = thresh / (1.0 - powc(af, 2) / awr2);
                }
            }
            // 670
            s.di_ecn = s.di_ecl;
            if s.di_ech < s.di_ecn * (1.0 - small) {
                s.di_ecn = s.di_ech;
            }
        }
        // 680
        if s.di_ecn < *enext * (1.0 - small) {
            *idisc = 1;
        }
        if s.di_ecn < *enext * (1.0 - small) {
            *enext = s.di_ecn;
        }
        *nq = s.di_nqp0;
        if s.di_iecl == s.di_iech {
            *nq = 0;
        }
    }

    /// Retrieve or compute Legendre coefficients at e. Return isotropic
    /// distribution below the first point.
    pub(super) fn getfle(
        &mut self,
        c: &mut Ctx,
        e: f64,
        enext: &mut f64,
        idis: &mut i32,
        fle: &mut F,
        nle: &mut i32,
        lcd: &mut i32,
        matd: i32,
        mfd: i32,
        mtd: i32,
        nin: i32,
    ) {
        let emax = 1.0e10f64;
        let small = 1.0e-10f64;
        let over = 1.01e0f64;
        let zero = 0.0f64;
        let (mut nb, mut nwc) = (0i32, 0i32);
        let mut fls = std::mem::take(&mut self.sv.fl_fls);

        //--initialize
        *idis = 0;
        if !(e > zero) {
            self.sv.fl_iso = 0;
            if mfd != 8 {
                c.findf(matd, mfd, mtd, nin);
                c.contio(nin, 0, 0, &mut fls.at(1), &mut nb, &mut nwc);
                self.awr = fls[2];
                let lvt = nint(fls[3]);
                self.sv.fl_ltt = nint(fls[4]);
                self.sv.fl_ltt3 = self.sv.fl_ltt;
                if self.sv.fl_ltt == 3 {
                    self.sv.fl_ltt = 1;
                    self.sv.fl_lttn = 1;
                }
                if lvt == 0 {
                    c.contio(nin, 0, 0, &mut fls.at(1), &mut nb, &mut nwc);
                } else {
                    c.listio(nin, 0, 0, &mut fls.at(1), &mut nb, &mut nwc);
                    while nb != 0 {
                        c.moreio(nin, 0, 0, &mut fls.at(7), &mut nb, &mut nwc);
                    }
                }
                self.sv.fl_iso = nint(fls[3]);
                self.lct = nint(fls[4]);
            } else {
                if self.law.abs() == 3 {
                    self.sv.fl_iso = 1;
                }
                if self.law == -4 {
                    self.sv.fl_iso = 1;
                }
                self.sv.fl_ltt = 1;
            }
            // check for lab distribution for two-body scattering.
            if self.lct == 1 && self.awr >= 10.0 {
                if mtd >= 51 && mtd <= 90 {
                    let s = format!("lab distribution changed to cm for mt={}", ffmt::i(mtd, 3));
                    c.mess("getfle", &s, " ");
                    self.lct = 2;
                }
            }
            self.sv.fl_elo = emax;
            if self.sv.fl_iso != 1 {
                c.tab2io(nin, 0, 0, &mut fls.at(1), &mut nb, &mut nwc);
                self.sv.fl_lidp = -1;
                if self.law.abs() == 5 {
                    self.sv.fl_lidp = nint(fls[3]);
                }
                self.sv.fl_ne = nint(fls[6]);
                // read in raw data at first two incident energies.
                // retrieve or compute legendre coefficients.
                self.sv.fl_iraw = 1 + nwc;
                let iraw = self.sv.fl_iraw;
                let ltt = self.sv.fl_ltt;
                let lidp = self.sv.fl_lidp;
                let cp = self.copar();
                let lct = self.lct;
                for k in 0..2 {
                    let mut ll = iraw;
                    if ltt == 1 {
                        c.listio(nin, 0, 0, &mut fls.at(ll), &mut nb, &mut nwc);
                    }
                    if ltt == 2 {
                        c.tab1io(nin, 0, 0, &mut fls.at(ll), &mut nb, &mut nwc);
                    }
                    while nb != 0 {
                        ll += nwc;
                        c.moreio(nin, 0, 0, &mut fls.at(ll), &mut nb, &mut nwc);
                    }
                    if k == 0 {
                        self.sv.fl_elo = fls[iraw + 1];
                        self.sv.fl_nlo = *nle;
                    } else {
                        self.sv.fl_ehi = fls[iraw + 1];
                        self.sv.fl_nhi = *nle;
                    }
                    if lidp >= 0 {
                        fls[iraw + 3] = lidp as f64;
                    }
                    let (fl, nlx) = if k == 0 {
                        (&mut self.sv.fl_flo, &mut self.sv.fl_nlo)
                    } else {
                        (&mut self.sv.fl_fhi, &mut self.sv.fl_nhi)
                    };
                    getco(c, &cp, &mut fl[1..], nlx, *lcd, &mut fls.0[(iraw - 1) as usize..], lct, ltt, idis);
                }
                // check for total isotropy.
                if self.sv.fl_ne == 2 && self.sv.fl_nlo == 1 && self.sv.fl_nhi == 1 {
                    self.sv.fl_iso = 1;
                }
                self.sv.fl_int = nint(fls[8]);
            }
            if self.sv.fl_iso == 1 {
                *nle = 1;
            }
            *enext = self.sv.fl_elo;
            self.sv.fl_nne = 2;
            self.sv.fl_ir = 1;
            self.sv.fl_fls = fls;
            return;
        }

        //--normal entry
        //--is desired energy in current panel
        let mut lbl;
        if *nle == 1 || self.sv.fl_iso == 1 {
            lbl = 440;
        } else if e < self.sv.fl_ehi * (1.0 - small) {
            lbl = 300;
        } else if self.sv.fl_nne == self.sv.fl_ne && e < over * self.sv.fl_ehi {
            if self.sv.fl_ltt3 == 3 && self.sv.fl_lttn == 1 {
                lbl = 210;
            } else {
                lbl = 300;
            }
        } else {
            lbl = 210;
        }
        if lbl == 210 {
            //--no.  slide raw high energy data into low energy positions.
            loop {
                let s = &mut self.sv;
                for i in 1..=s.fl_nhi as usize {
                    s.fl_flo[i] = s.fl_fhi[i];
                }
                s.fl_nlo = s.fl_nhi;
                s.fl_elo = s.fl_ehi;

                //--read in new raw data at high energy.
                if s.fl_nne == s.fl_ne && s.fl_ltt3 == 3 && s.fl_lttn == 1 {
                    c.tab2io(nin, 0, 0, &mut fls.at(1), &mut nb, &mut nwc);
                    s.fl_ne = nint(fls[6]);
                    s.fl_nne = 0;
                    s.fl_ir = 1;
                    s.fl_ltt = 2;
                    s.fl_lttn = 2;
                } else if s.fl_nne == s.fl_ne {
                    c.error("getfle", "desired energy above highest given.", " ");
                }
                let iraw = s.fl_iraw;
                let ltt = s.fl_ltt;
                let mut ll = iraw;
                if ltt == 1 {
                    c.listio(nin, 0, 0, &mut fls.at(ll), &mut nb, &mut nwc);
                }
                if ltt == 2 {
                    c.tab1io(nin, 0, 0, &mut fls.at(ll), &mut nb, &mut nwc);
                }
                while nb != 0 {
                    ll += nwc;
                    c.moreio(nin, 0, 0, &mut fls.at(ll), &mut nb, &mut nwc);
                }
                s.fl_ehi = fls[iraw + 1];
                s.fl_nhi = *nle;
                if s.fl_lidp >= 0 {
                    fls[iraw + 3] = s.fl_lidp as f64;
                }
                let cp = CoPar { awr: self.awr, awrp: self.awrp, izap: self.izap, izat: self.izat };
                let lct = self.lct;
                let s = &mut self.sv;
                getco(c, &cp, &mut s.fl_fhi[1..], &mut s.fl_nhi, *lcd, &mut fls.0[(iraw - 1) as usize..], lct, ltt, idis);
                s.fl_nne += 1;
                if s.fl_nne as f64 > fls[5 + 2 * s.fl_ir] {
                    s.fl_ir += 1;
                }
                s.fl_int = nint(fls[6 + 2 * s.fl_ir]);
                if s.fl_ehi <= e * (1.0 - small) && s.fl_nne < s.fl_ne {
                    continue;
                }
                break;
            }
            lbl = 300;
        }
        self.sv.fl_fls = fls;
        if lbl == 300 {
            //--yes.  interpolate for coefficients at desired energy.
            let s = &self.sv;
            if e >= s.fl_elo * (1.0 - small) {
                // normal energies
                let mut nlmax = s.fl_nlo;
                if nlmax < s.fl_nhi {
                    nlmax = s.fl_nhi;
                }
                for i in 1..=*nle {
                    if i <= nlmax {
                        c.terp1(s.fl_elo, s.fl_flo[i as usize], s.fl_ehi, s.fl_fhi[i as usize], e, &mut fle[i], s.fl_int);
                    } else {
                        fle[i] = 0.0;
                    }
                }
                *nle = nlmax;
                *enext = s.fl_ehi;
                if s.fl_int == 1 {
                    *idis = 1;
                }
            } else {
                //--return isotropic distribution below first point.
                fle[1] = 1.0;
                for il in 2..=*nle {
                    fle[il] = 0.0;
                }
                *nle = 1;
                *enext = s.fl_ehi;
                *idis = 1;
            }
            return;
        }

        //--isotropic distribution or request.
        // 440
        fle[1] = 1.0;
        if *nle > 1 {
            for il in 2..=*nle {
                fle[il] = 0.0;
            }
            *nle = 1;
        }
        *enext = emax;
        *idis = 0;
    }

    /// Legendre coefficients at e for thermal angle-energy distributions
    /// in File 6.
    pub(super) fn getaed(
        &mut self,
        c: &mut Ctx,
        aed: &mut V2,
        e: f64,
        enext: &mut f64,
        idis: &mut i32,
        nl: &mut i32,
        nq: &mut i32,
        mat: i32,
        mf: i32,
        mt: i32,
        nin: i32,
        _nlg: i32,
    ) {
        let emin = 1.0e-5f64;
        let emax = 10.0e0f64;
        let etop = 1.0e10f64;
        let eps = 1.0e-5f64;
        let shade = 0.999999e0f64;
        let step = 1.1e0f64;
        let up = 1.00001e0f64;
        let small = 1.0e-10f64;
        let zero = 0.0f64;
        let maxaes = 200000i32;
        let (mut nb, mut nw) = (0i32, 0i32);
        let ngn = self.ngn;
        let nlv = *nl;
        let mut p = [0.0f64; 21];

        //--initialize.
        if e == zero {
            let mut pp = RV::new(20 + SLACK);
            c.findf(mat, mf, mt, nin);
            c.contio(nin, 0, 0, &mut pp.at(1), &mut nb, &mut nw);
            c.tab1io(nin, 0, 0, &mut pp.at(1), &mut nb, &mut nw);
            self.sv.ae_itt = -c.l1h;
            self.sv.ae_law = c.l2h;
            while nb != 0 {
                c.moreio(nin, 0, 0, &mut pp.at(1), &mut nb, &mut nw);
            }

            if self.sv.ae_law == 1 {
                //--incoherent elastic or inelastic with E-E'-mu ordering
                self.aes = RV::new(maxaes as usize + SLACK);
                let aes = &mut self.aes;
                let s = &mut self.sv;
                let l1 = 1;
                c.tab2io(nin, 0, 0, &mut aes.at(l1), &mut nb, &mut nw);
                s.ae_ir = 1;
                s.ae_ne = nint(aes[l1 + 5]);
                s.ae_nwt = nw;
                s.ae_l2 = l1 + s.ae_nwt;
                c.listio(nin, 0, 0, &mut aes.at(s.ae_l2), &mut nb, &mut nw);
                s.ae_elo = aes[s.ae_l2 + 1];
                s.ae_ncyc = nint(aes[s.ae_l2 + 5]);
                s.ae_nu = s.ae_ncyc - 2;
                s.ae_nlo = nint(aes[s.ae_l2 + 4]) / s.ae_ncyc;
                let mut loc = s.ae_l2 + nw;
                while nb != 0 {
                    c.moreio(nin, 0, 0, &mut aes.at(loc), &mut nb, &mut nw);
                    loc += nw;
                    if loc > maxaes {
                        c.error("getaed", "storage exceeded.", " ");
                    }
                }
                s.ae_llo = loc;
                if loc > maxaes {
                    c.error("getaed", "storage exceeded.", " ");
                }
                s.ae_lhi = s.ae_llo;
                s.ae_nw1 = loc - s.ae_l2;
                s.ae_ehi = 0.0;
                *enext = s.ae_elo;
                s.ae_ie = 1;
            } else if self.sv.ae_law == 7 {
                c.error("getaed", "thermal mf6/law7 not coded", " ");
            } else {
                //--coherent elastic
                self.sv.ae_nbrag = self.sv.ae_itt;
                let nbrag = self.sv.ae_nbrag;
                self.aes = alloc(2 * nbrag);
                *enext = emin;
                let mut clast = 0.0f64;
                self.sv.ae_jbrag = 0;

                //--store bragg edges
                let mut idone = 0;
                let mut xs = 0.0f64;
                while idone == 0 {
                    let ebrag = *enext;
                    c.gety1(ebrag, enext, idis, &mut xs, self.npend, &mut self.sigma.at(1));
                    let cnow = ebrag * xs;
                    if *enext > emax * (1.0 + small) {
                        idone = 1;
                    } else if (*enext - emax).abs() <= eps {
                        idone = 1;
                    } else if (cnow - clast) > (eps * clast + eps) {
                        self.sv.ae_jbrag += 1;
                        let jb = self.sv.ae_jbrag;
                        self.aes[jb] = ebrag;
                        self.aes[nbrag + jb] = cnow - clast;
                        clast = cnow;
                    }
                }

                //--reposition npend for getsig
                c.findf(mat, 3, mt, self.npend);
                c.contio(self.npend, 0, 0, &mut self.sigma.at(1), &mut nb, &mut nw);
                c.gety1(e, enext, idis, &mut xs, self.npend, &mut self.sigma.at(1));
                *enext = etop;
            }
            return;
        }

        let mut aes = std::mem::take(&mut self.aes);
        let done = self.getaed_normal(c, &mut aes, aed, e, enext, idis, nlv, nq, nin, &mut p, maxaes);
        self.aes = aes;
        let _ = (done, ngn, up, shade, step);
    }

    /// getaed, normal entries.
    fn getaed_normal(
        &mut self,
        c: &mut Ctx,
        aes: &mut RV,
        aed: &mut V2,
        e: f64,
        enext: &mut f64,
        idis: &mut i32,
        nl: i32,
        nq: &mut i32,
        nin: i32,
        p: &mut [f64; 21],
        maxaes: i32,
    ) -> bool {
        let etop = 1.0e10f64;
        let shade = 0.999999e0f64;
        let step = 1.1e0f64;
        let up = 1.00001e0f64;
        let small = 1.0e-10f64;
        let zero = 0.0f64;
        let (mut nb, mut nw) = (0i32, 0i32);
        let ngn = self.ngn;
        let s = &mut self.sv;

        if s.ae_itt == 1 {
            //--normal entry for incoherent inelastic
            let mut lbl;
            if e < s.ae_ehi * (1.0 - small) {
                lbl = 300;
            } else if s.ae_ie == s.ae_ne && e > s.ae_ehi {
                lbl = 410;
            } else if s.ae_ie == s.ae_ne {
                lbl = 300;
            } else {
                let l1 = 1;
                s.ae_l2 = l1 + s.ae_nwt;
                s.ae_l3 = s.ae_l2 + s.ae_nw1;
                lbl = if s.ae_ie == 1 { 220 } else { 205 };
                loop {
                    if lbl == 205 {
                        //--slide high data into low positions.
                        for i in 1..=s.ae_nw2 {
                            aes[i - 1 + s.ae_l2] = aes[i - 1 + s.ae_l3];
                        }
                        s.ae_llo = s.ae_l2 + s.ae_lhi - s.ae_l3;
                        s.ae_nw1 = s.ae_nw2;
                        s.ae_l3 = s.ae_l2 + s.ae_nw1;
                        s.ae_elo = s.ae_ehi;
                        s.ae_nlo = s.ae_nhi;
                    }
                    // 220: read new high data.
                    c.listio(nin, 0, 0, &mut aes.at(s.ae_l3), &mut nb, &mut nw);
                    s.ae_ehi = aes[s.ae_l3 + 1];
                    s.ae_nhi = nint(aes[s.ae_l3 + 4]) / s.ae_ncyc;
                    let mut loc = s.ae_l3 + nw;
                    while nb != 0 {
                        c.moreio(nin, 0, 0, &mut aes.at(loc), &mut nb, &mut nw);
                        loc += nw;
                        if loc > maxaes {
                            c.error("getaed", "storage exceeded.", " ");
                        }
                    }
                    if loc > maxaes {
                        c.error("getaed", "storage exceeded.", " ");
                    }
                    s.ae_nw2 = loc - s.ae_l3;
                    s.ae_ie += 1;
                    if e >= s.ae_ehi * (1.0 - small) && s.ae_ie < s.ae_ne {
                        lbl = 205;
                        continue;
                    }
                    if !(s.ae_ie as f64 <= aes[l1 + 4 + 2 * s.ae_ir]) {
                        s.ae_ir += 1;
                    }
                    break;
                }
                lbl = 300;
            }
            if lbl == 410 {
                for il in 1..=nl {
                    for ig in 1..=ngn {
                        aed[(il, ig)] = 0.0;
                    }
                }
                *enext = etop;
                *nq = 0;
                return true;
            }

            //--begin loop over energy groups.
            // 300
            for il in 1..=nl {
                for ig in 1..=ngn {
                    aed[(il, ig)] = 0.0;
                }
            }
            let (l2, l3, ncyc, nlo, nhi, elo, ehi) = (s.ae_l2, s.ae_l3, s.ae_ncyc, s.ae_nlo, s.ae_nhi, s.ae_elo, s.ae_ehi);
            let mut k1 = 1;
            let mut k2 = 1;
            let mut fl = [0.0f64; 21];
            let mut fl1 = [0.0f64; 21];
            let mut fl2 = [0.0f64; 21];
            let mut fi = [0.0f64; 21];
            let mut eg = 0.0f64;
            let eb = elo / 10.0;
            let mut i = 0;
            let mut egp: f64;
            let mut egp1 = 0.0f64;
            let mut egp2 = 0.0f64;
            let (mut eg1, mut eg2): (f64, f64);
            let (mut ek1, mut ei1, mut ei2): (f64, f64, f64);
            let mut ek2 = 0.0f64;
            // 330
            'groups: loop {
                i += 1;
                egp = self.egn[i + 1];
                if i == ngn {
                    egp = egp + 1.0;
                }
                eg1 = egp1;
                eg2 = egp2;
                egp1 = egp - e + elo;
                egp2 = egp + ehi - e;
                if egp1 < eb * (1.0 - small) {
                    egp1 = egp * eb / (e - elo + eb);
                    egp2 = egp * (ehi - elo + eb) / (e - elo + eb);
                }
                // 345
                loop {
                    // get next point projected from low side
                    loop {
                        if k1 > nlo {
                            // 360
                            ek1 = aes[l3 + 6 + ncyc * (nhi - 1)];
                            break;
                        }
                        ek1 = aes[l2 + 6 + ncyc * (k1 - 1)];
                        if ek1 > up * eg1 {
                            break;
                        }
                        k1 += 1;
                    }
                    // 355
                    ei1 = ek1 + e - elo;
                    if ek1 < elo / 10.0 {
                        ei1 = ek1 * (e - elo + eb) / eb;
                    }
                    // 375: get next point projected from high side
                    let mut to450 = false;
                    loop {
                        if k2 > nhi {
                            to450 = true;
                            break;
                        }
                        ek2 = aes[l3 + 6 + ncyc * (k2 - 1)];
                        if ek2 > up * eg2 {
                            break;
                        }
                        k2 += 1;
                    }
                    if to450 {
                        break;
                    }
                    // 385
                    ei2 = ek2 - ehi + e;
                    if ek2 < ehi - elo + eb {
                        ei2 = ek2 * (e - elo + eb) / (ehi - elo + eb);
                    }

                    //--do integrals to next point
                    let mut ei = egp;
                    if ei1 < ei {
                        ei = ei1;
                    }
                    if ei2 < ei {
                        ei = ei2;
                    }
                    if (ei - egp).abs() < egp * small {
                        aedi(egp1, &mut fl1, nl, &aes.0[(l2 - 1) as usize..]);
                        aedi(egp2, &mut fl2, nl, &aes.0[(l3 - 1) as usize..]);
                    } else if (ei - ei1).abs() < ei1 * small {
                        aedi(ek1, &mut fl1, nl, &aes.0[(l2 - 1) as usize..]);
                        let mut ee = ek1 + ehi - elo;
                        if ek1 < eb * (1.0 - small) {
                            ee = ek1 * (ehi - elo + eb) / eb;
                        }
                        aedi(ee, &mut fl2, nl, &aes.0[(l3 - 1) as usize..]);
                        eg1 = ek1;
                        eg2 = ee;
                    } else if (ei - ei2).abs() < ei2 * small {
                        let mut ee = ek2 - ehi + elo;
                        if ee < eb * (1.0 - small) {
                            ee = ek2 * eb / (ehi - elo + eb);
                        }
                        aedi(ee, &mut fl1, nl, &aes.0[(l2 - 1) as usize..]);
                        aedi(ek2, &mut fl2, nl, &aes.0[(l3 - 1) as usize..]);
                        eg2 = ek2;
                        eg1 = ee;
                    }
                    let f1 = (ehi - e) / (ehi - elo);
                    let f2 = (e - elo) / (ehi - elo);
                    for il in 1..=nl as usize {
                        fi[il] = f1 * fl1[il] + f2 * fl2[il];
                    }
                    for il in 1..=nl {
                        let iu = il as usize;
                        aed[(il, i)] = aed[(il, i)] + (fi[iu] + fl[iu]) * (ei - eg) / 2.0;
                        fl[iu] = fi[iu];
                    }
                    eg = ei;
                    if ei < egp * (1.0 - small) {
                        continue;
                    }
                    break;
                }
                // 450: close loop over energy groups
                if i < ngn && (k1 < nlo || k2 < nhi) {
                    continue 'groups;
                }
                break;
            }
            *enext = ehi;
            if *enext <= e * (1.0 + small) {
                *enext = etop;
            }
            *nq = 0;
            return true;
        }

        //--normal entry for elastic
        //--find group index of in-group
        // 490
        let mut ig = 1;
        while e >= self.egn[ig + 1] * (1.0 - small) && ig != ngn {
            ig += 1;
        }
        for il in 1..=nl {
            for i in 1..=ngn {
                aed[(il, i)] = 0.0;
            }
        }
        if s.ae_itt > 2 {
            // 600: coherent elastic
            if s.ae_jbrag != 0 {
                let mut ib = 1;
                while e >= shade * aes[ib] && ib <= s.ae_jbrag {
                    let u = 1.0 - 2.0 * aes[ib] / e;
                    legndr(u, &mut p[..], nl);
                    for il in 1..=nl {
                        aed[(il, ig)] = aed[(il, ig)] + aes[s.ae_nbrag + ib] * p[(il - 1) as usize] / e;
                    }
                    ib += 1;
                }
                if aed[(1, ig)] != zero {
                    if nl >= 2 {
                        for il in 2..=nl {
                            aed[(il, ig)] = aed[(il, ig)] / aed[(1, ig)];
                        }
                    }
                }
            }
            aed[(1, ig)] = 1.0;
            *enext = self.egn[ig + 1];
            *idis = 1;
            *nq = 0;
            return true;
        }

        //--incoherent elastic
        let mut lbl = 0;
        if e < s.ae_ehi * (1.0 - small) {
            lbl = 560;
        } else if s.ae_ie == s.ae_ne && e > step * s.ae_ehi {
            for il in 1..=nl {
                for ig in 1..=ngn {
                    aed[(il, ig)] = 0.0;
                }
            }
            *enext = etop;
            *nq = 0;
            return true;
        } else if s.ae_ie == s.ae_ne {
            lbl = 560;
        }
        if lbl == 0 {
            if s.ae_ie != 1 {
                //--slide high data into low positions
                for i in 1..=s.ae_nw2 {
                    aes[i - 1 + s.ae_l2] = aes[i - 1 + s.ae_l3];
                }
                s.ae_nw1 = s.ae_nw2;
                s.ae_l3 = s.ae_l2 + s.ae_nw1;
                s.ae_elo = s.ae_ehi;
                s.ae_llo = s.ae_lhi;
            }
            // 530: read new high data
            s.ae_ehi = 0.0;
            while e >= s.ae_ehi * (1.0 - small) && s.ae_ie < s.ae_ne {
                c.listio(nin, 0, 0, &mut aes.at(s.ae_l3), &mut nb, &mut nw);
                s.ae_ehi = aes[s.ae_l3 + 1];
                s.ae_nhi = 1;
                let mut loc = s.ae_l3 + nw;
                while nb != 0 {
                    c.moreio(nin, 0, 0, &mut aes.at(loc), &mut nb, &mut nw);
                    loc += nw;
                }
                s.ae_nw2 = loc - s.ae_l3;
                s.ae_ie += 1;
            }
        }

        //--compute legendre components
        // 560
        let (l2, l3, elo, ehi, nu) = (s.ae_l2, s.ae_l3, s.ae_elo, s.ae_ehi, s.ae_nu);
        let x1 = aes[l2 + 7];
        let x2 = aes[l3 + 7];
        let mut xs = x1 + (x2 - x1) * (e - elo) / (ehi - elo);
        xs = xs / nu as f64;
        for iu in 1..=nu {
            let u1 = aes[l2 + 7 + iu];
            let u2 = aes[l3 + 7 + iu];
            let u = u1 + (u2 - u1) * (e - elo) / (ehi - elo);
            legndr(u, &mut p[..], nl);
            for il in 1..=nl {
                aed[(il, ig)] = aed[(il, ig)] + xs * p[(il - 1) as usize];
            }
        }
        *enext = ehi;
        if *enext > self.egn[ig + 1] * (1.0 + small) {
            *enext = self.egn[ig + 1];
        }
        *idis = 1;
        if (*enext - ehi).abs() < ehi * small {
            *idis = 0;
        }
        if *enext <= e * (1.0 + small) {
            *enext = etop;
        }
        *nq = 0;
        true
    }

    /// Retrieve gamma angular distributions as Legendre coefficients for
    /// all gamma rays. Initialize if ed=0.
    pub(super) fn getgfl(
        &mut self,
        c: &mut Ctx,
        ed: f64,
        enext: &mut f64,
        idis: &mut i32,
        gfl: &mut V2,
        nl: &mut i32,
        _nlg: i32,
        ng: i32,
        mat: i32,
        mf: i32,
        mt: i32,
        nin: i32,
    ) {
        let emax = 1.0e10f64;
        let small = 1.0e-10f64;
        let up = 1.1e0f64;
        let zero = 0.0f64;
        let maxgfl = 500;
        let (mut nb, mut nw) = (0i32, 0i32);
        let mut b = RV::new(6 + SLACK);

        if !(ed > zero) {
            //--initialize
            self.lct = 1;
            c.findf(mat, mf, mt, nin);
            c.contio(nin, 0, 0, &mut b.at(1), &mut nb, &mut nw);
            self.sv.gf_li = nint(b[3]);
            self.sv.gf_ltt = nint(b[4]);
            let ni = nint(b[6]);
            let ntmp = 10000;
            let mut tmp = RV::new(ntmp as usize + SLACK);
            *enext = emax;

            //--reduce legendre order if possible
            if self.sv.gf_li == 1 {
                *nl = 1;
            }
            if self.sv.gf_li == 1 {
                self.flgn = alloc(1);
                return;
            }

            //--find starting location of each sub-section
            let mut ig = 0;
            let mut il = 1;
            let loca = &mut self.sv.gf_loca;
            if ni > 0 {
                while ig < ni {
                    c.contio(nin, 0, 0, &mut tmp.at(il), &mut nb, &mut nw);
                    ig += 1;
                    if ig > maxgfl {
                        c.error("getgfl", "too many gammas", " ");
                    }
                    loca[ig as usize] = il - 1;
                    il += nw;
                }
            }
            let ltt = self.sv.gf_ltt;
            while ig < ng {
                c.tab2io(nin, 0, 0, &mut tmp.at(il), &mut nb, &mut nw);
                ig += 1;
                if ig > maxgfl {
                    c.error("getgfl", "too many gammas", " ");
                }
                loca[ig as usize] = il - 1;
                let np = nint(tmp[il + 5]);
                let mut nnp = 0;
                il += nw;
                let mut loc = il;
                while nnp < np {
                    if loc + c.npage > ntmp {
                        c.error("getgfl", "storage exceeded", " ");
                    }
                    if ltt == 1 {
                        c.listio(nin, 0, 0, &mut tmp.at(loc), &mut nb, &mut nw);
                    }
                    if ltt == 2 {
                        c.tab1io(nin, 0, 0, &mut tmp.at(loc), &mut nb, &mut nw);
                    }
                    loc += nw;
                    while nb != 0 {
                        if loc + c.npage > ntmp {
                            c.error("getgfl", "storage exceeded", " ");
                        }
                        c.moreio(nin, 0, 0, &mut tmp.at(loc), &mut nb, &mut nw);
                        loc += nw;
                    }
                    nnp += 1;
                    if nnp == 1 && tmp[il + 1] < *enext {
                        *enext = tmp[il + 1];
                    }
                    il = loc;
                }
            }
            let na = il - 1;
            self.flgn = alloc(na);
            for i in 1..=na {
                self.flgn[i] = tmp[i];
            }

            //--sort locations into descending order
            if ng > 1 {
                let nh = ng - ni;
                for j in 1..=nh {
                    let mut nipj = ni + j;
                    if ni == 0 {
                        nipj = ng;
                    }
                    for i in 1..=nipj {
                        let il = 1 + loca[i as usize];
                        let im = 1 + loca[nipj as usize];
                        if self.flgn[il] <= self.flgn[im] {
                            let mut m = -1;
                            let ip1 = i + 1;
                            for _k in ip1..=nipj {
                                m += 1;
                                loca[(nipj - m) as usize] = loca[(nipj - m - 1) as usize];
                            }
                            loca[i as usize] = im - 1;
                        }
                    }
                }
            }
            return;
        }

        //--normal entry
        *enext = emax;
        *idis = 0;
        let li = self.sv.gf_li;
        let ltt = self.sv.gf_ltt;
        let mut flgn = std::mem::take(&mut self.flgn);
        let cp = self.copar();
        let lct = self.lct;
        let mut alo = [0.0f64; 66];
        let mut ahi = [0.0f64; 66];
        for j in 1..=ng {
            let mut iso = li == 1;
            let mut l = 0;
            let mut lint = 0;
            let mut lnow = 0;
            if !iso {
                l = 1 + self.sv.gf_loca[j as usize];
                lint = l;
                lnow = nint(flgn[l + 4]);
                if !(lnow > 0) {
                    iso = true;
                }
            }
            if !iso {
                // 220: for anisotropic photons
                let np = nint(flgn[l + 5]);
                let mut im = 1;
                let mut inn = 1;
                let mut lnext = lnow * 2;
                let mut elo: f64;
                let mut ehi = 0.0f64;
                loop {
                    // 240
                    l = l + lnext + 6;
                    elo = flgn[l + 1];
                    if elo > ed * (1.0 + small) && im == 1 {
                        iso = true;
                        break;
                    }
                    if *enext < elo * (1.0 - small) {
                        *enext = elo;
                    }
                    lnext = nint(flgn[l + 4]);
                    if ltt == 2 {
                        lnext = 2 * nint(flgn[l + 4]) + 2 * nint(flgn[l + 5]);
                    }
                    ehi = flgn[l + lnext + 7];
                    im += 1;
                    if im > np {
                        c.error("getgfl", "desired energy at highest given energy.", " ");
                    }
                    let nbt = nint(flgn[lint + 4 + 2 * inn]);
                    if im > nbt {
                        inn += 1;
                    }
                    if ehi <= ed * (1.0 + small) && im < np {
                        continue;
                    }
                    break;
                }
                if !iso {
                    if im == np && ed > up * ehi {
                        c.error("getgfl", "desired energy at highest given energy.", " ");
                    }
                    let int = nint(flgn[lint + 5 + 2 * inn]);
                    let mut nld = *nl;
                    getco(c, &cp, &mut alo[1..], &mut nld, 1, &mut flgn.0[(l - 1) as usize..], lct, ltt, idis);
                    l = l + lnext + 6;
                    nld = *nl;
                    getco(c, &cp, &mut ahi[1..], &mut nld, 1, &mut flgn.0[(l - 1) as usize..], lct, ltt, idis);
                    for i in 1..=*nl {
                        let mut y = 0.0;
                        c.terp1(elo, alo[i as usize], ehi, ahi[i as usize], ed, &mut y, int);
                        gfl[(i, j)] = y;
                    }
                    if ehi > *enext * (1.0 + small) {
                        continue;
                    }
                    *idis = 0;
                    if int == 1 {
                        *idis = 1;
                    }
                    *enext = ehi;
                    continue;
                }
            }
            // 205: for isotropic photons and energies less than threshold
            gfl[(1, j)] = 1.0;
            if *nl == 1 {
                continue;
            }
            for i in 2..=*nl {
                gfl[(i, j)] = 0.0;
            }
        }
        self.flgn = flgn;
    }

    /// Retrieve gamma yields from mf12 or compute them from mf13.
    pub(super) fn getgyl(
        &mut self,
        c: &mut Ctx,
        ed: f64,
        enext: &mut f64,
        idis: &mut i32,
        gyl: &mut F,
        eyl: &mut F,
        nyl: &mut i32,
        matd: i32,
        mfd: i32,
        mtd: i32,
        nin: i32,
    ) {
        let emax = 1.0e10f64;
        let zero = 0.0f64;
        let nylmax = 550;
        let (mut nb, mut nw) = (0i32, 0i32);

        //--initialize
        if ed == zero {
            let ntmp = 99000;
            let mut tmp = RV::new(ntmp as usize + SLACK);
            c.findf(matd, mfd, mtd, nin);
            c.contio(nin, 0, 0, &mut tmp.at(1), &mut nb, &mut nw);
            self.sv.gy_awr = tmp[2];
            let lo = nint(tmp[3]);
            if lo == 2 {
                c.error("getgyl", "lo=2 not coded.", " ");
            }
            *nyl = nint(tmp[5]);
            if *nyl > nylmax {
                c.error("getgyl", "too many gammas.", " ");
            }
            // find starting location for each sub section
            let mut loc = 1;
            c.tab1io(nin, 0, 0, &mut tmp.at(loc), &mut nb, &mut nw);
            if *nyl != 1 || mfd != 13 {
                loc += nw;
                while nb != 0 {
                    c.moreio(nin, 0, 0, &mut tmp.at(loc), &mut nb, &mut nw);
                    loc += nw;
                }
            }
            let mut il = loc;
            let mut ik = 0;
            eyl[1] = tmp[1];
            let mut nr = nint(tmp[5]);
            *enext = tmp[7 + 2 * nr];
            let loca = &mut self.sv.gy_loca;
            loca[1] = 0;
            let na: i32;
            if *nyl == 1 {
                let mut n = loc - 1;
                if mfd == 13 {
                    n = 1;
                }
                na = n;
            } else {
                while ik < *nyl {
                    if loc + c.npage > ntmp {
                        c.error("getgyl", "storage exceeded.", " ");
                    }
                    c.tab1io(nin, 0, 0, &mut tmp.at(loc), &mut nb, &mut nw);
                    loc += nw;
                    while nb != 0 {
                        if loc + c.npage > ntmp {
                            c.error("getgyl", "storage exceeded.", " ");
                        }
                        c.moreio(nin, 0, 0, &mut tmp.at(loc), &mut nb, &mut nw);
                        loc += nw;
                    }
                    ik += 1;
                    loca[ik as usize] = il - 1;
                    nr = nint(tmp[il + 4]);
                    let en = tmp[il + 6 + 2 * nr];
                    if en < *enext {
                        *enext = en;
                    }
                    il = loc;
                }
                il = 1 + loca[ik as usize];
                eyl[1] = tmp[il];
                self.sv.gy_ipd = 2;
                self.sv.gy_ird = 1;
                na = loc - 1;
            }
            // check for constant range at low energies
            if mfd == 12 {
                let mut idone = 0;
                let mut i = 0;
                while i < *nyl && idone == 0 {
                    i += 1;
                    let l = 1 + self.sv.gy_loca[i as usize];
                    let ey = tmp[l];
                    let lp = nint(tmp[l + 2]);
                    let nr = nint(tmp[l + 4]);
                    let np = nint(tmp[l + 5]);
                    let ll = l + 5 + 2 * nr;
                    let mut e = 0.0f64;
                    //--if histogram interpolation, start with e equal
                    //  to the second energy point
                    if nint(tmp[l + 7]) == 1 {
                        e = tmp[ll + 3];
                    }
                    for j in 2..=np {
                        if tmp[ll + 2 * j] == tmp[ll + 2] {
                            e = tmp[ll + 2 * j - 1];
                        }
                    }
                    if lp == 2 && ey / 100.0 < e {
                        e = ey / 100.0;
                    }
                    if e < self.econst {
                        self.econst = e;
                    }
                    if self.econst == zero {
                        idone = 1;
                    }
                }
            } else {
                self.econst = zero;
            }
            self.gyln = alloc(na);
            for i in 1..=na {
                self.gyln[i] = tmp[i];
            }
            return;
        }

        //--normal entry
        *enext = emax;
        *idis = 0;
        if *nyl != 1 || mfd != 13 {
            let awr = self.sv.gy_awr;
            let mut idisc = 0i32;
            let mut en = 0.0f64;
            for i in 1..=*nyl {
                let l = 1 + self.sv.gy_loca[i as usize];
                eyl[i] = self.gyln[l];
                let lp = nint(self.gyln[l + 2]);
                if lp == 2 {
                    eyl[i] = eyl[i] + ed * awr / (awr + 1.0);
                }
                let mut ip = 2;
                let mut ir = 1;
                c.terpa(&mut gyl[i], ed, &mut en, &mut idisc, &self.gyln.at(l), &mut ip, &mut ir);
                if en == *enext && idisc > *idis {
                    *idis = idisc;
                }
                if en < *enext {
                    *idis = idisc;
                }
                if en < *enext {
                    *enext = en;
                }
                if lp == 2 {
                    let ngi = self.ngg + 1;
                    if self.gyln[l] != zero {
                        for j in 1..=ngi {
                            let e = (self.egg[j] - self.gyln[l]) * (awr + 1.0) / awr;
                            if e >= ed && e < *enext {
                                *idis = 1;
                                *enext = e;
                            }
                        }
                    }
                }
            }
            if mfd != 12 {
                // divide by cross section to get yield
                let mut sig = 0.0f64;
                c.terpa(&mut sig, ed, &mut en, idis, &self.gyln.at(1), &mut self.sv.gy_ipd, &mut self.sv.gy_ird);
                if en == *enext && idisc > *idis {
                    *idis = idisc;
                }
                if en < *enext {
                    *idis = idisc;
                }
                if en < *enext {
                    *enext = en;
                }
                if sig != zero {
                    let rsig = 1.0 / sig;
                    for i in 1..=*nyl {
                        gyl[i] = gyl[i] * rsig;
                    }
                }
            }
        } else {
            gyl[1] = 1.0;
            eyl[1] = self.gyln[1];
        }
    }

    /// Secondary energy distribution for all sink groups. Initialize if
    /// ed=0 and return the number of subsections in nk. On a normal
    /// entry, nk=1 sums the subsections; nk>1 returns them separately.
    pub(super) fn getsed(
        &mut self,
        c: &mut Ctx,
        ed: f64,
        enext: &mut f64,
        idis: &mut i32,
        sed: &mut V2,
        which: i32,
        ng: i32,
        nk: &mut i32,
        matd: i32,
        mfd: i32,
        mtd: i32,
        nin: i32,
    ) {
        let emax = 1.0e10f64;
        let eps = 0.001e0f64;
        let small = 1.0e-10f64;
        let ebig = 1.0e8f64;
        let tenth = 0.1e0f64;
        let zero = 0.0f64;
        let nkmax = 20;
        let (mut nb, mut nw) = (0i32, 0i32);
        let mut ethi = 0.0f64;

        //--initialize
        if ed == zero {
            let ntmp: i32 = 300000;
            let mut tmp = RV::new(ntmp as usize + SLACK);
            c.findf(matd, mfd, mtd, nin);
            c.contio(nin, 0, 0, &mut tmp.at(1), &mut nb, &mut nw);
            self.sv.sd_nktot = nint(tmp[5]);
            if self.sv.sd_nktot > nkmax {
                c.error("getsed", "too many subsections.", " ");
            }
            *nk = self.sv.sd_nktot;
            let mut l = 1;
            let mut il: i32;
            let mut ik = 0;
            let mut m = 1;
            *enext = emax;
            while ik < *nk {
                ik += 1;
                il = l - 1;
                self.sv.sd_loc[ik as usize] = il;
                c.tab1io(nin, 0, 0, &mut tmp.at(l), &mut nb, &mut nw);
                l += nw;
                while nb != 0 {
                    if l + nw > ntmp {
                        c.error("getsed", "storage for tmp exceeded1", " ");
                    }
                    c.moreio(nin, 0, 0, &mut tmp.at(l), &mut nb, &mut nw);
                    l += nw;
                }
                let lf = nint(tmp[il + 4]);
                let mut nr = nint(tmp[il + 5]);
                let elo = tmp[il + 1 + 6 + 2 * nr];
                if elo < *enext {
                    *enext = elo;
                }

                if lf > 1 {
                    //--analytic subsection.
                    if lf != 3 {
                        c.tab1io(nin, 0, 0, &mut tmp.at(l), &mut nb, &mut nw);
                        nr = nint(tmp[l + 4]);
                        let mut np = nint(tmp[l + 5]);
                        let mut ln = l + 2 * nr + 5;
                        l += nw;
                        while nb != 0 {
                            c.moreio(nin, 0, 0, &mut tmp.at(l), &mut nb, &mut nw);
                            l += nw;
                        }
                        for ip in 2..=np {
                            let delta = (tmp[ln + 2 * ip] - tmp[ln + 2 * ip - 2]).abs();
                            if delta >= eps * tmp[ln + 2 * ip - 2] {
                                let ee = tmp[ln + 2 * ip - 3]
                                    + eps * tmp[ln + 2 * ip - 2] * (tmp[ln + 2 * ip - 1] - tmp[ln + 2 * ip - 3]) / delta;
                                if ee < self.econst {
                                    self.econst = ee;
                                }
                            }
                        }
                        if lf == 5 || lf == 11 {
                            let l1 = l;
                            c.tab1io(nin, 0, 0, &mut tmp.at(l), &mut nb, &mut nw);
                            nr = nint(tmp[l + 4]);
                            np = nint(tmp[l + 5]);
                            ln = l + 2 * nr + 5;
                            l += nw;
                            while nb != 0 {
                                c.moreio(nin, 0, 0, &mut tmp.at(l), &mut nb, &mut nw);
                                l += nw;
                            }
                            //extend lowest delayed bin using sqrt(e) shape
                            if self.ismooth > 0 && mtd == 455 && nint(tmp[l1 + 7]) == 1 {
                                l = self.sqrt_bin(c, &mut tmp, l1, ntmp, "storage for tmp exceeded3");
                            }
                            if lf != 5 {
                                for ip in 2..=np {
                                    let delta = (tmp[ln + 2 * ip] - tmp[ln + 2 * ip - 2]).abs();
                                    if delta >= eps * tmp[ln + 2 * ip - 2] {
                                        let ee = tmp[ln + 2 * ip - 3]
                                            + eps * tmp[ln + 2 * ip - 2] * (tmp[ln + 2 * ip - 1] - tmp[ln + 2 * ip - 3]) / delta;
                                        if ee < self.econst {
                                            self.econst = ee;
                                        }
                                    }
                                }
                            }
                        }
                    }
                    m = l;
                } else {
                    //--tabulated subsection
                    c.tab2io(nin, 0, 0, &mut tmp.at(l), &mut nb, &mut nw);
                    let mut ne = nint(tmp[l + 5]);
                    l += nw;
                    let iraw = l + ne * (ng + 1);
                    m = iraw;
                    let mut nne = 0;

                    //--read and average spectrum for each incident energy
                    while nne < ne {
                        let m1 = m;
                        c.tab1io(nin, 0, 0, &mut tmp.at(m), &mut nb, &mut nw);
                        m += nw;
                        while nb != 0 {
                            if m > ntmp {
                                c.error("getsed", "storage for tmp exceeded2", " ");
                            }
                            c.moreio(nin, 0, 0, &mut tmp.at(m), &mut nb, &mut nw);
                            m += nw;
                        }
                        // extend lowest delayed bin using sqrt(e) shape
                        if self.ismooth > 0 && mtd == 455 && nint(tmp[m1 + 7]) == 1 {
                            m = self.sqrt_bin(c, &mut tmp, m1, ntmp, "storage for tmp exceeded4");
                        }
                        nne += 1;
                        tmp[l] = tmp[m1 + 1];
                        l += 1;
                        let mut ip = 2;
                        let mut ir = 1;
                        for ig in 1..=ng {
                            let mut e1 = self.egv(which, ig);
                            if ig == 1 {
                                e1 = 0.0;
                            }
                            let mut e2 = self.egv(which, ig + 1);
                            if ig == ng {
                                e2 = ebig;
                            }
                            let mut f = 0.0f64;
                            c.intega(&mut f, e1, e2, &tmp.at(m1), &mut ip, &mut ir);
                            tmp[l] = f;
                            l += 1;
                        }
                    }

                    //--determine econst
                    l = self.sv.sd_loc[ik as usize] + 1;
                    let mut nr = nint(tmp[l + 4]);
                    let np = nint(tmp[l + 5]);
                    l = l + 6 + 2 * nr + 2 * np;
                    nr = nint(tmp[l + 4]);
                    ne = nint(tmp[l + 5]);
                    let mut ir = 1;
                    let mut nbt = nint(tmp[l + 4 + 2 * ir]);
                    let mut int = nint(tmp[l + 5 + 2 * ir]);
                    if int == 2 {
                        int = 22;
                    }
                    if int >= 11 && int <= 15 {
                        c.mess("getsed", "corresponding point interpolation not available", " ");
                    }
                    let lnow = l + 6 + 2 * nr;
                    for i in 2..=ne {
                        if i > nbt {
                            ir += 1;
                            nbt = nint(tmp[l + 4 + 2 * ir]);
                            int = nint(tmp[l + 5 + 2 * ir]);
                            if int == 2 {
                                int = 22;
                            }
                        }
                        let llo = lnow + (i - 2) * (ng + 1);
                        let elo = tmp[llo];
                        let lhi = lnow + (i - 1) * (ng + 1);
                        let ehi = tmp[lhi];
                        for ig in 1..=ng {
                            let delta = (tmp[lhi + ig] - tmp[llo + ig]).abs();
                            if delta > eps * tmp[llo + ig] {
                                if int == 1 && ehi < self.econst {
                                    self.econst = ehi;
                                }
                                if int != 1 {
                                    let ee = elo + eps * tmp[llo + ig] * (ehi - elo) / delta;
                                    if ee < self.econst {
                                        self.econst = ee;
                                    }
                                }
                            }
                        }
                    }
                    l = m;
                }
            }

            //--initialization complete
            let nc = m - 1;
            self.sedist = alloc(nc);
            for i in 1..=nc {
                self.sedist[i] = tmp[i];
            }
            self.sv.sd_nupm = 0;
            return;
        }

        //--normal entry.
        *enext = emax;
        *idis = 0;
        for ik in 1..=*nk {
            for ig in 1..=ng {
                sed[(ik, ig)] = 0.0;
            }
        }
        let mut sd = std::mem::take(&mut self.sedist);
        let mut iout = 0;
        let mut ik = 0;
        let mut int = 0i32;
        let mut ehi = 0.0f64;
        let mut nne = 0;
        let mut ne = 0;
        let mut idisc = 0i32;
        while ik < self.sv.sd_nktot && iout == 0 {
            ik += 1;
            let mut ikt = ik;
            if *nk == 1 {
                ikt = 1;
            }
            let lnow = self.sv.sd_loc[ik as usize] + 1;
            let lf = nint(sd[lnow + 3]);

            //--interpolate for fractional probability.
            let mut ip = 2;
            let mut ir = 1;
            let mut pe = 0.0f64;
            let mut eihi = 0.0f64;
            c.terpa(&mut pe, ed, &mut eihi, &mut idisc, &sd.at(lnow), &mut ip, &mut ir);
            if (eihi - *enext).abs() < *enext * small && idisc > *idis {
                *idis = idisc;
            }
            if eihi < *enext * (1.0 - small) {
                *idis = idisc;
            }
            if eihi < *enext * (1.0 - small) {
                *enext = eihi;
            }
            if pe > zero {
                let nr = nint(sd[lnow + 4]);
                let np = nint(sd[lnow + 5]);
                let mnow = lnow + 6 + 2 * nr + 2 * np;

                if lf > 1 {
                    //--analytic subsection.  compute sed at e.
                    ip = 2;
                    ir = 1;
                    for ig in 1..=ng {
                        let mut e1 = self.egv(which, ig);
                        if ig == 1 {
                            e1 = 0.0;
                        }
                        let mut e2 = self.egv(which, ig + 1);
                        if ig == ng {
                            e2 = ebig;
                        }
                        let mut s = 0.0f64;
                        self.anased(c, &mut s, ed, &mut ethi, &mut idisc, e1, e2, &mut sd.0[(lnow - 1) as usize..], &mut ip, &mut ir);
                        s = s * pe;
                        sed[(ikt, ig)] = sed[(ikt, ig)] + s;
                    }
                    if (ethi - *enext).abs() < *enext * small && idisc > *idis {
                        *idis = idisc;
                    }
                    if ethi < *enext * (1.0 - small) {
                        *idis = idisc;
                    }
                    if ethi < *enext * (1.0 - small) {
                        *enext = ethi;
                    }
                } else {
                    //--tabulated subsection.  interpolate for sed at e.
                    nne = 0;
                    let nr = nint(sd[mnow + 4]);
                    ne = nint(sd[mnow + 5]);
                    let nnow = mnow + 6 + 2 * nr;
                    let mut nnt = 6;
                    let mut nbt = nint(sd[mnow + nnt]);
                    int = nint(sd[mnow + nnt + 1]);
                    if int == 2 {
                        int = 22;
                    }
                    ehi = 0.0;
                    let (mut llo, mut lhi) = (0i32, 0i32);
                    let mut elo = 0.0f64;
                    while nne < ne && ed > ehi * (1.0 + small) && iout == 0 {
                        nne += 1;
                        if nne > ne {
                            iout = 1;
                        } else {
                            if nne > nbt {
                                nnt += 2;
                                nbt = nint(sd[mnow + nnt]);
                                int = nint(sd[mnow + nnt + 1]);
                                if int == 2 {
                                    int = 22;
                                }
                            }
                            llo = nnow + (ng + 1) * (nne - 1);
                            elo = sd[llo];
                            lhi = nnow + (ng + 1) * nne;
                            ehi = sd[lhi];
                        }
                    }
                    let (mut klo, mut khi) = (0i32, 0i32);
                    let (mut xlo, mut xhi) = (0.0f64, 0.0f64);
                    if int > 5 {
                        klo = nnow + (ng + 1) * ne;
                        elo = sd[klo + 1];
                        let mut nrlo = nint(sd[klo + 4]);
                        let mut nplo = nint(sd[klo + 5]);
                        xlo = sd[klo + 4 + 2 * nrlo + 2 * nplo];
                        khi = klo + 6 + 2 * nrlo + 2 * nplo;
                        ehi = sd[khi + 1];
                        let mut nrhi = nint(sd[khi + 4]);
                        let mut nphi = nint(sd[khi + 5]);
                        xhi = sd[khi + 4 + 2 * nrhi + 2 * nphi];
                        while nne < ne && ed > ehi * (1.0 + small) && iout == 0 {
                            klo = khi;
                            elo = ehi;
                            nrlo = nrhi;
                            nplo = nphi;
                            xlo = xhi;
                            khi = klo + 6 + 2 * nrlo + 2 * nplo;
                            ehi = sd[khi + 1];
                            nrhi = nint(sd[khi + 4]);
                            nphi = nint(sd[khi + 5]);
                            xhi = sd[khi + 4 + 2 * nrhi + 2 * nphi];
                        }
                    }
                    if iout == 0 {
                        if int >= 21 {
                            //--unit base.
                            let jnt = int - 20;
                            let mut xend = 0.0f64;
                            c.terp1(elo, xlo, ehi, xhi, ed, &mut xend, jnt);
                            for ig in 1..=ng {
                                let mut e1 = self.egv(which, ig);
                                if ig == 1 {
                                    e1 = 0.0;
                                }
                                let e1lo = e1 * xlo / xend;
                                let e1hi = e1 * xhi / xend;
                                let mut e2 = self.egv(which, ig + 1);
                                if ig == ng {
                                    e2 = ebig;
                                }
                                let e2lo = e2 * xlo / xend;
                                let e2hi = e2 * xhi / xend;
                                let (mut flo, mut fhi, mut fe) = (0.0f64, 0.0f64, 0.0f64);
                                ip = 2;
                                ir = 1;
                                c.intega(&mut flo, e1lo, e2lo, &sd.at(klo), &mut ip, &mut ir);
                                ip = 2;
                                ir = 1;
                                c.intega(&mut fhi, e1hi, e2hi, &sd.at(khi), &mut ip, &mut ir);
                                c.terp1(elo, flo, ehi, fhi, ed, &mut fe, jnt);
                                sed[(ikt, ig)] = sed[(ikt, ig)] + pe * fe;
                            }
                        } else if int >= 11 {
                            //--corresponding points.  not implemented.
                            for ig in 1..=ng {
                                sed[(ikt, ig)] = 0.0;
                            }
                        } else {
                            //--cartesion.
                            for ig in 1..=ng {
                                let mut s = 0.0f64;
                                c.terp1(elo, sd[llo + ig], ehi, sd[lhi + ig], ed, &mut s, int);
                                sed[(ikt, ig)] = sed[(ikt, ig)] + s * pe;
                            }
                        }
                        //--upscatter is not allowed in secondary energy
                        //--spectra, except for fission and gamma-ray
                        //--production. put the upscatters into the in-group.
                        let mut sigup = 0.0f64;
                        if mfd != 15 && mtd != 455 {
                            if mtd < 18 || (mtd > 21 && mtd != 38) {
                                let mut ig = ng;
                                while ed <= self.egv(which, ig) * (1.0 + eps) && ig > 1 {
                                    sed[(ikt, ig - 1)] = sed[(ikt, ig - 1)] + sed[(ikt, ig)];
                                    sigup = sed[(ikt, ig)];
                                    sed[(ikt, ig)] = 0.0;
                                    ig -= 1;
                                }
                            }
                        }
                        let test = tenth;
                        if sigup > test && self.sv.sd_nupm < 10 {
                            self.sv.sd_nupm += 1;
                            let s = format!("upscatter correction {}", ffmt::e1p(sigup, 12, 4));
                            c.mess("getsed", &s, " ");
                            if self.sv.sd_nupm == 10 {
                                c.mess("getsed", "additional messages supressed", " ");
                            }
                        }
                        if int == 1 && ehi < *enext {
                            idisc = 1;
                        }
                        if ehi < *enext * (1.0 - small) {
                            *enext = ehi;
                        }
                    }
                }
            }
        }
        self.sedist = sd;

        //--return the normal results
        if iout == 0 {
            return;
        }
        //--return zeroes outside range of table
        if int != 1 {
            for ig in 1..=ng {
                sed[(1, ig)] = 0.0;
            }
            *enext = ehi;
            if nne == ne {
                *enext = emax;
            }
            *idis = 1;
        } else {
            for ig in 1..=ng {
                sed[(1, ig)] = ehi;
            }
            *enext = emax;
        }
    }

    /// getsed: extend the lowest delayed bin using a sqrt(E) shape.
    /// Returns the new end of the record.
    fn sqrt_bin(&mut self, c: &mut Ctx, tmp: &mut RV, l1: i32, ntmp: i32, msg: &str) -> i32 {
        let ex = 40.0f64;
        let fx = super::mf6::f32lit(0.8409);
        c.io.so(" extending lowest delayed bin using sqrt(E)");
        let mut mm = nint(tmp[l1 + 5]);
        while tmp[l1 + 10] > ex {
            if l1 + 9 + 2 * mm > ntmp {
                c.error("getsed", msg, " ");
            }
            let mut ix = 2 * mm;
            while ix >= 1 {
                tmp[l1 + 9 + ix] = tmp[l1 + 7 + ix];
                ix -= 1;
            }
            tmp[l1 + 10] = fx * tmp[l1 + 12];
            let val = tmp[l1 + 9];
            tmp[l1 + 9] = sqrt(fx) * val;
            tmp[l1 + 9] = sigfig(tmp[l1 + 9], 7, 0);
            tmp[l1 + 11] = (1.0 - fx * sqrt(fx)) * val / (1.0 - fx);
            tmp[l1 + 11] = sigfig(tmp[l1 + 11], 7, 0);
            mm += 1;
        }
        tmp[l1 + 5] = mm as f64;
        tmp[l1 + 6] = mm as f64;
        l1 + 8 + 2 * mm
    }

    /// Analytic secondary energy distributions: the integral between
    /// sink energies ep1 and ep2 for source energy e. a is a 0-based
    /// slice of a(*).
    pub(super) fn anased(
        &mut self,
        c: &mut Ctx,
        g: &mut f64,
        e: f64,
        enext: &mut f64,
        idis: &mut i32,
        ep1: f64,
        ep2: f64,
        a: &mut [f64],
        ip: &mut i32,
        ir: &mut i32,
    ) {
        let xmin = 1.0e-3f64;
        let brk = 0.01e0f64;
        let small = 1.0e-10f64;
        let emax = 1.0e10f64;
        let etest = 1.1e-5f64;
        let zero = 0.0f64;
        let rp4 = sqrt(PI / 4.0);
        let epl = ep1;
        let eph = ep2;
        let av = |a: &[f64], i: i32| a[(i - 1) as usize];
        self.sv.an_new = 0;
        if ep1 < etest {
            self.sv.an_new = 1;
        }
        let new = self.sv.an_new;

        // check limit on integration.
        *enext = emax;
        *idis = 0;
        *g = 0.0;
        let lf = nint(av(a, 4));
        let mut de = 0.0f64;
        if lf != 12 {
            let u = av(a, 1);
            de = e - u;
            if epl > de * (1.0 + small) || de <= zero {
                return;
            }
        }

        //--retrieve theta by interpolation.
        if new != 0 {
            let nr = nint(av(a, 5));
            let np = nint(av(a, 6));
            let mut ip2 = 2;
            let mut ir2 = 1;
            self.sv.an_loct = 6 + 2 * nr + 2 * np;
            let loct = self.sv.an_loct;
            let mut temp1 = 0.0f64;
            let mut en = 0.0f64;
            let mut idisc = 0i32;
            c.terpa(&mut temp1, e, &mut en, &mut idisc, &F(&mut a[loct as usize..]), &mut ip2, &mut ir2);
            self.sv.an_theta = temp1;
            if lf != 12 {
                self.sv.an_xc = de / self.sv.an_theta;
            }
            if en < *enext * (1.0 - small) {
                *idis = idisc;
            }
            if en < *enext * (1.0 - small) {
                *enext = en;
            }
        }
        let theta = self.sv.an_theta;
        let xc = self.sv.an_xc;
        let loct = self.sv.an_loct;

        //--compute common quantities and integration limits.
        let (mut xlo, mut xhi) = (0.0f64, 0.0f64);
        if lf != 12 {
            xlo = epl / theta;
            if eph <= de * (1.0 + small) {
                xhi = eph / theta;
            }
            if eph > de * (1.0 + small) {
                xhi = xc;
            }
        }

        if lf == 5 {
            //--law 5.  general evaporation spectrum.
            let nr = nint(av(a, loct + 5));
            let np = nint(av(a, loct + 6));
            let locg = loct + 6 + 2 * nr + 2 * np;
            c.intega(g, xlo, xhi, &F(&mut a[locg as usize..]), ip, ir);
            *g = *g * theta;
        } else if lf == 7 {
            //--law 7.  simple fission spectrum.
            if new != 0 {
                self.sv.an_rc = sqrt(xc);
                let r4 = -xc;
                let rc = self.sv.an_rc;
                self.sv.an_bot = rp4 * (1.0 - erfc(rc)) - rc * exp(r4);
            }
            let rl = sqrt(xlo);
            let rh = sqrt(xhi);
            let r3 = -xlo;
            let r4 = -xhi;
            let top = rl * exp(r3) - rh * exp(r4) - rp4 * (erfc(rh) - erfc(rl));
            *g = top / self.sv.an_bot;
        } else if lf == 9 {
            //--law 9.  evaporation spectrum.
            self.sv.an_rc = xc;
            let rc = xc;
            if new != 0 {
                if xc < xmin {
                    self.sv.an_bot = rc * rc / 2.0;
                }
                if xc >= xmin {
                    self.sv.an_bot = 1.0 - exp(-rc) * (1.0 + rc);
                }
            }
            let r1 = xlo;
            let r2 = xhi;
            let mut top = 0.0f64;
            if xhi < xmin {
                top = (r2 * r2 - r1 * r1) / 2.0;
            }
            if xhi >= xmin {
                top = exp(-r1) * (r1 + 1.0) - exp(-r2) * (r2 + 1.0);
            }
            *g = top / self.sv.an_bot;
        } else if lf == 11 {
            //--law 11. energy dependent watt spectrum.
            if new != 0 {
                let nr = nint(av(a, loct + 5));
                let np = nint(av(a, loct + 6));
                *ip = 2;
                *ir = 1;
                let locb = loct + 6 + 2 * nr + 2 * np;
                let mut b = 0.0f64;
                let mut en = 0.0f64;
                let mut idisc = 0i32;
                c.terpa(&mut b, e, &mut en, &mut idisc, &F(&mut a[locb as usize..]), ip, ir);
                let mut rc = theta * b / 4.0;
                rc = sqrt(rc);
                self.sv.an_rc = rc;
                self.sv.an_ca = 4.0 * rc * exp(-rc * rc) / (3.0 * rp4);
            }
            let rc = self.sv.an_rc;
            let mut h = [0.0f64; 2];
            let mut rh = sqrt(xc);
            let r1 = -rc;
            let r2 = rh - rc;
            let r3 = rc;
            let r4 = rh + rc;
            self.hnab(&mut h, r1, r2);
            let mut bot = h[1] + rc * h[0];
            self.hnab(&mut h, r3, r4);
            bot = bot - h[1] + rc * h[0];
            self.sv.an_bot = bot;
            let rl = sqrt(xlo);
            rh = sqrt(xhi);
            let top;
            if rh >= brk {
                let r1 = rl - rc;
                let r2 = rh - rc;
                let r3 = rl + rc;
                let r4 = rh + rc;
                self.hnab(&mut h, r1, r2);
                let mut t = h[1] + rc * h[0];
                self.hnab(&mut h, r3, r4);
                t = t - h[1] + rc * h[0];
                top = t;
            } else {
                top = self.sv.an_ca * (powc(rh, 3) - powc(rl, 3));
            }
            *g = top / bot;
        } else if lf == 12 {
            //--law 12.  madland-nix fission spectrum
            let _ = (c, &mut *g);
            crate::fio::fatal_runtime("groupr: anased law 12 (Madland-Nix) is not ported");
        } else {
            c.error("anased", "illegal lf.", " ");
        }
    }

    /// Integral from a to b of (1/sqrt(pi))*(u**n)*exp(-u**2) for n=0
    /// and n=1. hh[0] is hh(1). Writes the module variable beta, as the
    /// Fortran does.
    pub(super) fn hnab(&mut self, hh: &mut [f64; 2], aa: f64, bb: f64) {
        let toler = 1.0e-5f64;
        const POW2: [f64; 5] = [1.4142135623731e0, 2.0e0, 2.8284271247462e0, 4.0e0, 5.6568542494924e0];
        let one = 1.0f64;
        let epsrel = 1.0e-9f64;
        let epsabs = 1.0e-30f64;
        let zero = 0.0f64;
        let resqpi = 1.0 / sqrt(PI);
        let mut cm = [0.0f64; 51];
        let mut cmstar = [0.0f64; 51];
        let mut fa = [0.0f64; 3];
        let mut fb = [0.0f64; 3];

        let expa = exp(-aa * aa);
        let expb = exp(-bb * bb);
        fa[1] = 1.0 - erfc(aa.abs());
        if aa < zero {
            fa[1] = -fa[1];
        }
        fb[1] = 1.0 - erfc(bb.abs());
        if bb < zero {
            fb[1] = -fb[1];
        }
        fa[2] = (1.0 - expa) * resqpi;
        fb[2] = (1.0 - expb) * resqpi;
        hh[0] = fb[1] - fa[1];
        let mut n: i32;
        let mut lbl = 0;
        if (hh[0].abs() > toler * fb[1].abs() && hh[0].abs() > toler * fa[1].abs()) || (aa - bb).abs() >= one {
            lbl = 110;
        }
        n = 0;
        loop {
            if lbl == 110 {
                hh[1] = fb[2] - fa[2];
                if (hh[1].abs() > toler * fb[2].abs() && hh[1].abs() > toler * fa[2].abs()) || (aa - bb).abs() >= one {
                    return;
                }
                n = 1;
            }
            //--compute as a taylor expansion
            // 120
            let mut sgn = 1.0f64;
            let (a, b);
            if bb < aa {
                a = bb.abs();
                b = aa.abs();
                sgn = -sgn;
            } else {
                a = aa.abs();
                b = bb.abs();
            }
            if bb < 0.0 && n % 2 != 0 {
                sgn = -sgn;
            }
            let h = (b - a) * POW2[0];
            let x = POW2[0] * a;
            let xx = x * x;
            let asq = a * a;
            let con = exp(-asq) * resqpi / POW2[n as usize];
            let mut mflag = 0;
            let mut k = n;
            let mut kd = 0;
            cm[1] = 1.0;
            let mut xk = 1.0f64;
            if k != 0 {
                xk = powi(x, k);
            }
            let mut s = h * xk;
            let mut fact = h;
            let mut idone = 0;
            let mut m = 1;
            while idone == 0 && m < 50 {
                m += 1;
                fact = fact * h / m as f64;
                let kstar = k;
                let kdstar = kd + 1;
                for j in 1..=kdstar as usize {
                    cmstar[j] = cm[j];
                }
                k = n - m + 1;
                if k < 0 {
                    let kk = k % 2;
                    k = 0;
                    if kk != 0 {
                        k = 1;
                    }
                }
                kd = (n + m - 1 - k) / 2;
                cm[(kd + 1) as usize] = -cmstar[kdstar as usize];
                let mut qmn = cm[(kd + 1) as usize];
                if kd != 0 {
                    for j in 1..=kd {
                        let jalpha = (2 * j + k - 1 - kstar) / 2;
                        self.beta = 0.0;
                        if jalpha != 0 {
                            self.beta = cmstar[jalpha as usize];
                        }
                        cm[j as usize] = (2 * j + k - 1) as f64 * cmstar[(jalpha + 1) as usize] - self.beta;
                    }
                    for j in 1..=kd {
                        qmn = qmn * xx + cm[(kd + 1 - j) as usize];
                    }
                }
                xk = 1.0;
                if k != 0 {
                    xk = powi(x, k);
                }
                let term = fact * xk * qmn;
                s = s + term;
                let xn1 = (n + 1) as f64;
                let xn2 = (h * x).abs();
                if m as f64 >= xn1.max(xn2) {
                    let test = epsabs + epsrel * s.abs();
                    if term.abs() <= test {
                        if mflag == 1 {
                            idone = 1;
                        }
                        mflag = 1;
                    } else {
                        mflag = 0;
                    }
                }
            }
            hh[n as usize] = 2.0 * con * s * sgn;
            if n == 0 {
                lbl = 110;
                continue;
            }
            return;
        }
    }
}

/// Retrieve or compute Legendre coefficients from raw ENDF data (see the
/// Fortran comments). fl is a 0-based slice of fl(*), cdat of c(*).
pub(super) fn getco(
    c: &mut Ctx,
    cp: &CoPar,
    fl: &mut [f64],
    nl: &mut i32,
    lcd: i32,
    cdat: &mut [f64],
    lct: i32,
    ltt: i32,
    idis: &mut i32,
) {
    let toler = 1.0e-8f64;
    let nqp = 64;
    let nlmax = 65;
    let zero = 0.0f64;
    let mut p = [0.0f64; 66];
    let l = *nl - 1;
    if *nl > nlmax {
        c.error("getco", "limited to 64 legendre coefficients.", " ");
    }
    let cv = |cd: &[f64], i: i32| cd[(i - 1) as usize];
    let mut lang = -1;
    if c.mfh == 6 {
        lang = c.l1h;
    }
    let direct = !(lcd != lct || lang == 1 || lang > 2 || ltt == 2);
    if direct {
        //--retrieve coefficients from raw data
        if lang != 2 {
            let np = nint(cv(cdat, 5)) + 1;
            fl[0] = 1.0;
            let mut nlz = 1;
            for il in 2..=*nl {
                if il <= np {
                    fl[(il - 1) as usize] = cv(cdat, il + 5);
                    if fl[(il - 1) as usize] != zero {
                        nlz = il;
                    }
                } else {
                    fl[(il - 1) as usize] = 0.0;
                }
            }
            *nl = nlz;
        } else {
            let np = nint(cv(cdat, 5));
            let mut nlz = 1;
            for il in 1..=*nl {
                if il <= np {
                    fl[(il - 1) as usize] = cv(cdat, il + 6);
                    if fl[(il - 1) as usize] != zero {
                        nlz = il;
                    }
                } else {
                    fl[(il - 1) as usize] = 0.0;
                }
            }
            *nl = nlz;
        }
        return;
    }

    //--integrate for fl in desired system
    // 300
    for il in 1..=*nl {
        fl[(il - 1) as usize] = 0.0;
    }
    let mut ipc = 2;
    let mut irc = 1;
    let mut ipm = 0;
    let mut intn = 0;
    if lang > 2 {
        let nmu = c.n2h;
        ipm = 6 + 2 * nmu;
        intn = c.l1h - 10;
    }
    let mut xnext = 0.0f64;
    for iq in 1..=nqp {
        let x = GQP[(iq - 1) as usize];
        let mut y = 0.0f64;
        if ltt == 2 || lang > 2 {
            // 350: retrieve scattering probability from raw tabular data
            if c.mfh != 6 {
                c.terpa(&mut y, x, &mut xnext, idis, &F(cdat), &mut ipc, &mut irc);
            } else {
                // 352
                loop {
                    if x < cv(cdat, 7) {
                        y = 0.0;
                        break;
                    }
                    if ipc > ipm {
                        y = 0.0;
                        break;
                    }
                    if x <= cv(cdat, ipc + 7) {
                        c.terp1(cv(cdat, ipc + 5), cv(cdat, ipc + 6), cv(cdat, ipc + 7), cv(cdat, ipc + 8), x, &mut y, intn);
                        break;
                    }
                    ipc += 2;
                }
            }
        } else if lang == 1 {
            // 360: coulomb nuclear+interference by amplitude expansion
            if !(cp.izap <= 1 || c.mth != 2) && lang == 1 {
                y = coulomb_y(cp, cdat, x, &mut p);
            }
        } else {
            //--compute scattering probability from raw coefficients
            let np = nint(cv(cdat, 5));
            legndr(x, &mut p, np);
            y = 1.0;
            y = y / 2.0;
            for ip in 1..=np {
                y = y + (2 * ip + 1) as f64 * p[ip as usize] * cv(cdat, ip + 6) / 2.0;
            }
        }

        //--multiply scattering probability by legendre polynomials in
        //--desired system and accumulate
        // 370
        let mut xn = x;
        if lcd == 2 && lct == 1 {
            c.error("getco", "lab to cm conversion not coded.", " ");
        }
        if lcd == 1 && lct == 2 {
            xn = (1.0 + cp.awr * x) / sqrt(1.0 + cp.awr * cp.awr + 2.0 * cp.awr * x);
        }
        legndr(xn, &mut p, l);
        y = y * GQW[(iq - 1) as usize];
        for il in 1..=*nl {
            fl[(il - 1) as usize] = fl[(il - 1) as usize] + y * p[(il - 1) as usize];
        }
    }

    //--reduce number of significant figures
    //--and count number of coefficients required
    let mut nlz = 1;
    for il in 2..=*nl {
        let j = nint(fl[(il - 1) as usize] / toler);
        fl[(il - 1) as usize] = j as f64 * toler;
        if j != 0 {
            nlz = il;
        }
    }
    *nl = nlz;
}

/// Complex number helpers for getco (gfortran -fcx-fortran-rules: plain
/// formulas for multiplication).
#[derive(Clone, Copy)]
struct Cx(f64, f64);
impl Cx {
    fn add(self, o: Cx) -> Cx {
        Cx(self.0 + o.0, self.1 + o.1)
    }
    fn mul(self, o: Cx) -> Cx {
        Cx(self.0 * o.0 - self.1 * o.1, self.0 * o.1 + self.1 * o.0)
    }
    fn mulr(self, r: f64) -> Cx {
        Cx(self.0 * r, self.1 * r)
    }
    fn divr(self, r: f64) -> Cx {
        Cx(self.0 / r, self.1 / r)
    }
    fn exp(self) -> Cx {
        let e = exp(self.0);
        Cx(e * cos(self.1), e * sin(self.1))
    }
}

/// getco label 360: charged-particle elastic nuclear+interference term.
/// `cmplx(a,b)` without KIND gives a default (single precision) complex.
fn coulomb_y(cp: &CoPar, cdat: &[f64], x: f64, p: &mut [f64; 66]) -> f64 {
    let cv = |i: i32| cdat[(i - 1) as usize];
    let c4 = |a: f64, b: f64| Cx(a as f32 as f64, b as f32 as f64);
    let zero = 0.0f64;
    let uno = 1.0f64;
    let e = cv(2);
    let lidp = nint(cv(4));
    let ai = cp.awrp * AMASSN;
    let zt = (cp.izat / 1000) as f64;
    let zi = (cp.izap / 1000) as f64;
    let ee = (EV / 10000000.0) * (CLIGHT / 10.0);
    let cc2 = powc(ee, 4) * AMU / (2.0 * powc(HBAR, 2) * EV);
    let eta = zt * zi * sqrt(cc2 * ai / e);
    let nt = nint(cv(6));
    let np = 2 * nt;
    legndr(x, &mut p[..], np);
    let iunit = c4(zero, uno);
    if lidp != 1 {
        let mut sigr = cv(7) / 2.0;
        for ip in 1..=np {
            sigr = sigr + (2 * ip + 1) as f64 * p[ip as usize] * cv(ip + 7) / 2.0;
        }
        let mut cs1 = c4(cv(8 + np), cv(9 + np)).divr(2.0);
        for it in 1..=nt {
            cs1 = cs1.add(c4(cv(8 + np + 2 * it), cv(9 + np + 2 * it)).mulr((2 * it + 1) as f64 * p[it as usize]).divr(2.0));
        }
        let carg1 = iunit.mul(Cx(eta, 0.0)).mul(Cx(log((1.0 - x) / 2.0), 0.0));
        let sigi = (-2.0 * eta / (1.0 - x)) * cs1.mul(carg1.exp()).0;
        (1.0 - x) * (sigr + sigi)
    } else {
        let mut sigr = cv(7) / 2.0;
        for it in 1..=nt {
            sigr = sigr + (4 * it + 1) as f64 * p[(2 * it) as usize] * cv(it + 7) / 2.0;
        }
        let mut cs1 = c4(cv(8 + nt), cv(9 + nt)).divr(2.0);
        let mut cs2 = cs1;
        let mut sgn = -1.0f64;
        for it in 1..=nt {
            let z = c4(cv(8 + nt + 2 * it), cv(9 + nt + 2 * it));
            cs1 = cs1.add(z.mulr((2 * it + 1) as f64 * p[it as usize]).divr(2.0));
            cs2 = cs2.add(z.mulr(sgn * (2 * it + 1) as f64 * p[it as usize]).divr(2.0));
            sgn = -sgn;
        }
        let carg1 = iunit.mul(Cx(eta, 0.0)).mul(Cx(log((1.0 - x) / 2.0), 0.0));
        let carg2 = iunit.mul(Cx(eta, 0.0)).mul(Cx(log((1.0 + x) / 2.0), 0.0));
        let t = cs1.mul(Cx(1.0 + x, 0.0)).mul(carg1.exp()).add(cs2.mul(Cx(1.0 - x, 0.0)).mul(carg2.exp()));
        let sigi = (-2.0 * eta / (1.0 - x * x)) * t.0;
        (1.0 - x * x) * (sigr + sigi)
    }
}

/// Interpolate for Legendre components of thermal scattering at ee. aa
/// is a 0-based slice of aa(*); fl[il] is fl(il).
pub(super) fn aedi(ee: f64, fl: &mut [f64; 21], nl: i32, aa: &[f64]) {
    let av = |i: i32| aa[(i - 1) as usize];
    let mut p = [0.0f64; 21];
    let ncyc = nint(av(6));
    let nw = nint(av(5));
    let np = nw / ncyc;
    let nu = ncyc - 2;
    for il in 1..=nl as usize {
        fl[il] = 0.0;
    }
    let mut ip = 0;
    let mut i = 1;
    while ip == 0 && i < np {
        if ee <= av(7 + ncyc * i) {
            ip = i;
        }
        i += 1;
    }
    if ip > 0 {
        let f1 = (av(7 + ncyc * ip) - ee) / (av(7 + ncyc * ip) - av(7 + ncyc * (ip - 1)));
        let f2 = (ee - av(7 + ncyc * (ip - 1))) / (av(7 + ncyc * ip) - av(7 + ncyc * (ip - 1)));
        for iu in 1..=nu {
            let u = av(8 + iu + ncyc * (ip - 1));
            legndr(u, &mut p, nl);
            for il in 1..=nl as usize {
                fl[il] = fl[il] + f1 * av(8 + ncyc * (ip - 1)) * p[il - 1] / nu as f64;
            }
            let u = av(8 + iu + ncyc * ip);
            legndr(u, &mut p, nl);
            for il in 1..=nl as usize {
                fl[il] = fl[il] + f2 * av(8 + ncyc * ip) * p[il - 1] / nu as f64;
            }
        }
    }
}

const GQP: [f64; 64] = [
    -9.99305042e-01,
    -9.96340117e-01,
    -9.91013371e-01,
    -9.83336254e-01,
    -9.73326828e-01,
    -9.61008800e-01,
    -9.46411375e-01,
    -9.29569172e-01,
    -9.10522137e-01,
    -8.89315446e-01,
    -8.65999398e-01,
    -8.40629296e-01,
    -8.13265315e-01,
    -7.83972359e-01,
    -7.52819907e-01,
    -7.19881850e-01,
    -6.85236313e-01,
    -6.48965471e-01,
    -6.11155355e-01,
    -5.71895646e-01,
    -5.31279464e-01,
    -4.89403146e-01,
    -4.46366017e-01,
    -4.02270158e-01,
    -3.57220158e-01,
    -3.11322872e-01,
    -2.64687162e-01,
    -2.17423644e-01,
    -1.69644420e-01,
    -1.21462819e-01,
    -7.29931218e-02,
    -2.43502927e-02,
    2.43502927e-02,
    7.29931218e-02,
    1.21462819e-01,
    1.69644420e-01,
    2.17423644e-01,
    2.64687162e-01,
    3.11322872e-01,
    3.57220158e-01,
    4.02270158e-01,
    4.46366017e-01,
    4.89403146e-01,
    5.31279464e-01,
    5.71895646e-01,
    6.11155355e-01,
    6.48965471e-01,
    6.85236313e-01,
    7.19881850e-01,
    7.52819907e-01,
    7.83972359e-01,
    8.13265315e-01,
    8.40629296e-01,
    8.65999398e-01,
    8.89315446e-01,
    9.10522137e-01,
    9.29569172e-01,
    9.46411375e-01,
    9.61008800e-01,
    9.73326828e-01,
    9.83336254e-01,
    9.91013371e-01,
    9.96340117e-01,
    9.99305042e-01,
];

const GQW: [f64; 64] = [
    1.78328072e-03,
    4.14703326e-03,
    6.50445797e-03,
    8.84675983e-03,
    1.11681395e-02,
    1.34630479e-02,
    1.57260305e-02,
    1.79517158e-02,
    2.01348232e-02,
    2.22701738e-02,
    2.43527026e-02,
    2.63774697e-02,
    2.83396726e-02,
    3.02346571e-02,
    3.20579284e-02,
    3.38051618e-02,
    3.54722133e-02,
    3.70551285e-02,
    3.85501532e-02,
    3.99537411e-02,
    4.12625632e-02,
    4.24735151e-02,
    4.35837245e-02,
    4.45905582e-02,
    4.54916279e-02,
    4.62847966e-02,
    4.69681828e-02,
    4.75401657e-02,
    4.79993886e-02,
    4.83447622e-02,
    4.85754674e-02,
    4.86909570e-02,
    4.86909570e-02,
    4.85754674e-02,
    4.83447622e-02,
    4.79993886e-02,
    4.75401657e-02,
    4.69681828e-02,
    4.62847966e-02,
    4.54916279e-02,
    4.45905582e-02,
    4.35837245e-02,
    4.24735151e-02,
    4.12625632e-02,
    3.99537411e-02,
    3.85501532e-02,
    3.70551285e-02,
    3.54722133e-02,
    3.38051618e-02,
    3.20579284e-02,
    3.02346571e-02,
    2.83396726e-02,
    2.63774697e-02,
    2.43527026e-02,
    2.22701738e-02,
    2.01348232e-02,
    1.79517158e-02,
    1.57260305e-02,
    1.34630479e-02,
    1.11681395e-02,
    8.84675983e-03,
    6.50445797e-03,
    4.14703326e-03,
    1.78328072e-03,
];

#[allow(dead_code)]
fn _unused() -> i32 {
    ipow(1, 1)
}
