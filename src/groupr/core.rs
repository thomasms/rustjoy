//! GROUPR reaction set-up, panel integration, printing, and retrieval of
//! flux, yields and cross sections (init, panel, displa, getflx, getyld,
//! getsig, stounr, getunr, terpu, getff, parts).

use super::*;

impl Groupr {
    /// Select parameters appropriate for this reaction type.
    pub(super) fn init(&mut self, c: &mut Ctx, ng: &mut i32, nl: &mut i32, nz: &mut i32) {
        let naed = 14;
        let ebig = 1.0e8f64;
        let (mfd, mtd) = (self.mfd, self.mtd);

        //--select control parameters
        *ng = 2;
        if mfd == 6 {
            *ng = self.ngn + 1;
        }
        if mfd == 8 {
            *ng = self.ngn + 1;
        }
        if mtd == 251 || mtd == 252 || mtd == 253 {
            *ng = 3;
        }
        if mtd == 452 || mtd == 455 || mtd == 456 {
            *ng = 3;
        }
        if mfd == 5 {
            *ng = self.ngn + 1;
        }
        if mfd == 16 || mfd == 17 {
            *ng = self.ngg + 1;
        }
        if mfd == 18 {
            *ng = self.ngg + 1;
        }
        if mfd >= 21 && mfd <= 26 {
            *ng = self.ngn + 1;
        }
        if mfd >= 31 && mfd <= 36 {
            *ng = self.ngn + 1;
        }
        *nl = 1;
        if mtd == 1 && self.nsigz > 1 {
            *nl = 2;
        }
        if mfd == 5 && mtd == 455 {
            *nl = self.ndelg;
        }
        if mfd == 6 || mfd == 8 || mfd == 16 || mfd == 17 || mfd == 18 {
            *nl = self.lord + 1;
        }
        if mfd >= 21 && mfd <= 26 {
            *nl = self.lord + 1;
        }
        if mfd >= 31 && mfd <= 36 {
            *nl = self.lord + 1;
        }
        if mtd == 18 {
            *nl = 1;
        }
        if mtd == 19 || mtd == 20 || mtd == 21 || mtd == 38 {
            *nl = 1;
        }
        *nz = 1;
        let zl3 = [1, 2, 18, 19, 51, 102, 301, 302, 318, 402, 443, 444];
        if mfd == 3 && zl3.contains(&mtd) {
            *nz = self.nsigz;
        }
        if mfd == 6 && (mtd == 2 || mtd == 18) {
            *nz = self.nsigz;
        }
        self.econst = 0.0;
        if mfd != 3 && mfd != 8 && mfd != 18 && mfd < 10000000 {
            if mtd == 18 || mtd == 19 {
                self.econst = ebig;
            }
            if mtd == 452 || mtd == 455 || mtd == 456 {
                self.econst = ebig;
            }
            if mtd == 102 {
                self.econst = ebig;
            }
        }

        //--copy file 6 to npend2 for use by retrieval routines.
        if mfd == 6 && mtd >= 221 && mtd <= 250 {
            let mut tmp = RV::new((c.npage + 50) as usize + SLACK);
            c.findf(self.matb, 6, mtd, self.npend);
            self.npend2 = naed;
            if self.npend < 0 {
                self.npend2 = -self.npend2;
            }
            c.repoz(self.npend2);
            c.nsc = 1;
            c.math = 1;
            c.afend(0, self.npend2);
            c.tosend(self.npend, 0, self.npend2, &mut tmp.at(1));
            c.amend(0, self.npend2);
            c.atend(0, self.npend2);
            c.repoz(self.npend2);
            c.skiprz(self.npend, -4);
        }
    }

    /// Generalized group constant integrals for one panel.
    pub(super) fn panel(
        &mut self,
        c: &mut Ctx,
        elo: &mut f64,
        ehi: &mut f64,
        ans: &mut R3,
        ff: &mut [f64],
        nl: &mut i32,
        nz: i32,
        ng: &mut i32,
        iglo: &mut i32,
        nlg: i32,
        ng2g: i32,
    ) {
        const QP2: [f64; 2] = [-1.0e0, 1.0e0];
        const QW2: [f64; 2] = [1.0e0, 1.0e0];
        const QP6: [f64; 6] = [-1.0e0, -0.76505532e0, -0.28523152e0, 0.28523152e0, 0.76505532e0, 1.0e0];
        const QW6: [f64; 6] = [0.06666667e0, 0.37847496e0, 0.55485838e0, 0.55485838e0, 0.37847496e0, 0.06666667e0];
        const QP10: [f64; 10] = [
            -1.0e0,
            -0.9195339082e0,
            -0.7387738651e0,
            -0.4779249498e0,
            -0.1652789577e0,
            0.1652789577e0,
            0.4779249498e0,
            0.7387738651e0,
            0.9195339082e0,
            1.0e0,
        ];
        const QW10: [f64; 10] = [
            0.0222222222e0,
            0.1333059908e0,
            0.2248893420e0,
            0.2920426836e0,
            0.3275397612e0,
            0.3275397612e0,
            0.2920426836e0,
            0.2248893420e0,
            0.1333059908e0,
            0.0222222222e0,
        ];
        let rndoff = 1.000002e0f64;
        let delta = 0.999995e0f64;
        let emax = 1.0e10f64;
        let small = 1.0e-10f64;
        let zero = 0.0f64;
        let smin = 1.0e-9f64;
        let ttest = 0.1e0f64;
        let trange = 1.0e-6f64;
        let mut yld = 0.0f64;
        let nlv = *nl;

        //--retrieve factors in integrands at lower boundary
        if self.matd != self.p_matl || self.mfd != self.p_mfl || self.mtd != self.p_mtl || self.ipan == 0 {
            self.ipan = 1;
            self.p_matl = self.matd;
            self.p_mfl = self.mfd;
            self.p_mtl = self.mtd;
            self.p_ig1w = 0;
            self.p_iglow = 0;
            self.p_eqw = 0.0;
            self.p_nq = 0;
            self.p_elast = 0.0;
            self.p_idisc = 0;
            //--reset these arrays whenever matd, mfd or mtd changes
            self.stmp = RM::new(nz as usize, nlv as usize);
            self.slst = RM::new(nz as usize, nlv as usize);
            self.ftmp = RM::new(nz as usize, nlv as usize);
            self.flst = RM::new(nz as usize, nlv as usize);
            self.stmp_alloc = true;
        }
        let elow = *elo;
        if *elo > *ehi * (1.0 + small) {
            c.error("panel", "elo.gt.ehi.", " ");
        }
        let mut go115 = false;
        if (*elo - self.p_elast).abs() >= self.p_elast * small {
            if *elo * rndoff < *ehi {
                *elo = *elo * rndoff;
            }
            self.p_elast = *elo;
            let mut fl = std::mem::take(&mut self.flst);
            let mut en = 0.0f64;
            let mut idiscf = 0i32;
            self.getflx(c, *elo, &mut en, &mut idiscf, &mut fl, nlv, nz);
            self.flst = fl;
            let mut sl = std::mem::take(&mut self.slst);
            let (mut enext, mut idisc) = (self.p_enext, self.p_idisc);
            self.getsig(c, *elo, &mut enext, &mut idisc, &mut sl, nlv, nz);
            self.slst = sl;
            self.p_enext = enext;
            self.p_idisc = idisc;
            self.p_en = en;
            self.p_idiscf = idiscf;
            if self.p_enext >= emax * (1.0 - small) {
                go115 = true;
            } else {
                self.merge_next(small);
                let (mut ng1, mut ig1, mut nq) = (self.p_ng1, self.p_ig1, self.p_nq);
                let mut en = self.p_en;
                let mut idiscf = self.p_idiscf;
                self.getff(
                    c,
                    *elo,
                    &mut en,
                    &mut idiscf,
                    &mut yld,
                    &mut V2::new(ff, nlg as usize),
                    nl,
                    &mut ng1,
                    &mut ig1,
                    &mut nq,
                    nlg,
                    ng2g,
                );
                self.p_ng1 = ng1;
                self.p_ig1 = ig1;
                self.p_nq = nq;
                self.p_en = en;
                self.p_idiscf = idiscf;
                self.merge_next(small);
                self.p_nq += 2;
                if self.p_nq > 10 {
                    self.p_nq = 10;
                }
            }
        }
        let nlv = *nl;

        if !go115 {
            //--retrieve cross section and flux at upper boundary
            if self.p_enext < delta * *ehi {
                *ehi = self.p_enext;
                self.p_ehigh = *ehi;
                if self.p_idisc > 0 && *ehi * delta > *elo {
                    self.p_ehigh = *ehi * delta;
                }
            } else {
                self.p_ehigh = delta * *ehi;
            }
            let mut ft = std::mem::take(&mut self.ftmp);
            let (mut en, mut idiscf) = (self.p_en, self.p_idiscf);
            self.getflx(c, self.p_ehigh, &mut en, &mut idiscf, &mut ft, nlv, nz);
            self.ftmp = ft;
            self.p_en = en;
            self.p_idiscf = idiscf;
            let mut st = std::mem::take(&mut self.stmp);
            let (mut enext, mut idisc) = (self.p_enext, self.p_idisc);
            self.getsig(c, self.p_ehigh, &mut enext, &mut idisc, &mut st, nlv, nz);
            self.stmp = st;
            self.p_enext = enext;
            self.p_idisc = idisc;
            if self.p_enext > emax * (1.0 - small) {
                // 110: check for zero cross section over entire panel
                if self.stmp[(1, 1)] == zero {
                    go115 = true;
                }
            } else {
                self.merge_next(small);
            }
        }

        if go115 {
            // 115
            for iz in 1..=nz {
                for il in 1..=nlv {
                    ans[(il, iz, 1)] = 1.0;
                    ans[(il, iz, 2)] = 0.0;
                }
            }
            *ng = 2;
            *iglo = 1;
            self.p_en = emax;
            self.p_nqp = 0;
        } else {
            // 125
            //--compute group fluxes assuming flux is linear over panel.
            let aq = (*ehi + elow) / 2.0;
            let bq = (*ehi - elow) / 2.0;
            for iz in 1..=nz {
                for il in 1..=nlv {
                    ans[(il, iz, 1)] = ans[(il, iz, 1)] + (self.ftmp[(iz, il)] + self.flst[(iz, il)]) * bq;
                }
            }

            //--integrate over panel using lobatto quadrature.
            let nq = self.p_nq;
            for iq in 1..=nq {
                let (mut eq, wq);
                let k = (iq - 1) as usize;
                if nq == 2 {
                    eq = aq + bq * QP2[k];
                    wq = bq * QW2[k];
                } else if nq == 6 {
                    eq = aq + bq * QP6[k];
                    wq = bq * QW6[k];
                } else if nq == 10 {
                    eq = aq + bq * QP10[k];
                    wq = bq * QW10[k];
                } else {
                    c.error("panel", "bad nq in panel", " ");
                }
                eq = sigfig(eq, 9, 0);
                let t1 = (eq - elow) / (*ehi - elow);

                //--retrieve the feed function at the quadrature point.
                //--first point was last point of previous panel.
                if iq > 1 {
                    if eq > self.p_ehigh {
                        eq = self.p_ehigh;
                    }
                    let (mut ng1, mut ig1, mut nqp) = (self.p_ng1, self.p_ig1, self.p_nqp);
                    let (mut en, mut idiscf) = (self.p_en, self.p_idiscf);
                    self.getff(
                        c,
                        eq,
                        &mut en,
                        &mut idiscf,
                        &mut yld,
                        &mut V2::new(ff, nlg as usize),
                        nl,
                        &mut ng1,
                        &mut ig1,
                        &mut nqp,
                        nlg,
                        ng2g,
                    );
                    self.p_ng1 = ng1;
                    self.p_ig1 = ig1;
                    self.p_nqp = nqp;
                    self.p_en = en;
                    self.p_idiscf = idiscf;
                }

                //--accumulate the ng*nl*nz integrals simultaneously.
                if *iglo == 0 {
                    *iglo = self.p_ig1;
                }
                let ffv = V2::new(ff, nlg as usize);
                let ig1 = self.p_ig1;
                for iz in 1..=nz {
                    for il in 1..=*nl {
                        let a = self.stmp[(iz, il)] * self.ftmp[(iz, il)];
                        let b = self.slst[(iz, il)] * self.flst[(iz, il)];
                        let mut rr = b + (a - b) * t1;
                        rr = rr * wq;
                        for ig in 1..=self.p_ng1 {
                            if ffv[(il, ig)] != zero {
                                let igt = ig1 + ig - *iglo + 1;
                                if igt <= 1 {
                                    if ig1 != self.p_ig1w || *iglo != self.p_iglow || (eq - self.p_eqw).abs() >= trange {
                                        let s1 = format!("thermal range problem at{}", ffmt::e1p(eq, 12, 4));
                                        let s2 = format!(
                                            "scattering to group {}.  expected {}",
                                            ffmt::i(ig1, 3),
                                            ffmt::i(*iglo, 3)
                                        );
                                        c.mess("panel", &s1, &s2);
                                        self.p_ig1w = ig1;
                                        self.p_iglow = *iglo;
                                        self.p_eqw = eq;
                                    }
                                }
                                if igt > 1 {
                                    ans[(il, iz, igt)] = ans[(il, iz, igt)] + rr * ffv[(il, ig)];
                                    if igt > *ng {
                                        *ng = igt;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            //--remove partial thermal group.
            if self.mtd >= 221 && self.mtd <= 250 {
                if *elo >= ttest {
                    if self.stmp[(1, 1)].abs() <= smin {
                        for iz in 1..=nz {
                            for il in 1..=*nl {
                                ans[(il, iz, 2)] = 0.0;
                            }
                        }
                        *ng = 2;
                    }
                }
            }
        }

        //--save cross section and flux for next panel.
        //--determine next point and quadrature order from last ff.
        // 220
        self.p_elast = self.p_ehigh;
        for iz in 1..=nz {
            for il in 1..=*nl {
                self.slst[(iz, il)] = self.stmp[(iz, il)];
                self.flst[(iz, il)] = self.ftmp[(iz, il)];
            }
        }
        self.merge_next(small);
        if self.p_enext <= *ehi {
            self.p_enext = rndoff * *ehi;
        }
        self.p_nq = self.p_nqp + 2;
        if self.p_nq > 10 {
            self.p_nq = 10;
        }
    }

    /// panel: take the smaller of en and enext (with its discontinuity flag).
    fn merge_next(&mut self, small: f64) {
        let (en, idiscf) = (self.p_en, self.p_idiscf);
        if (en - self.p_enext).abs() < self.p_enext * small && idiscf > self.p_idisc {
            self.p_idisc = idiscf;
        }
        if en < self.p_enext * (1.0 - small) {
            self.p_idisc = idiscf;
        }
        if en < self.p_enext * (1.0 - small) {
            self.p_enext = en;
        }
    }

    /// Display generalized group constants generated by panel. Remove
    /// small values and transform ratio records. Write a heading if
    /// called with ig<0.
    pub(super) fn displa(
        &mut self,
        c: &mut Ctx,
        ig: i32,
        ans: &mut R3,
        nl: i32,
        nz: i32,
        ng2: &mut i32,
        ig2lo: &mut i32,
        igzero: &mut i32,
    ) {
        let smin = 1.0e-9f64;
        let big = 1.0e10f64;
        let zero = 0.0f64;
        *igzero = 0;
        let nlnz10 = nl.max(nz).max(10);
        let mut result = RV::new(nlnz10 as usize + 1);
        let mut field: Vec<String> = vec![" ".repeat(10); nlnz10 as usize + 1];
        let mfd = self.mfd;
        let mtd = self.mtd;
        let iprint = self.iprint;
        let fl = |field: &Vec<String>, a: i32, b: i32| -> Vec<String> {
            let mut v = Vec::new();
            let mut i = a;
            while i <= b {
                v.push(field[i as usize].clone());
                i += 1;
            }
            v
        };
        // total format reversion of `(head, 1p, 9a11)`
        let revall = |head: &str, f: &[String], n: usize| -> String {
            if f.is_empty() {
                return head.to_string();
            }
            let mut out: Vec<String> = Vec::new();
            for ch in f.chunks(n) {
                let mut s = head.to_string();
                for x in ch {
                    s.push_str(&a11(x));
                }
                out.push(s);
            }
            out.join("\n")
        };

        //--write appropriate heading.
        if ig < 0 {
            c.io.so(" ");
            self.d_itr = 0;
            let mut lbl = 0;
            if mfd == 6
                || mfd == 8
                || mfd == 16
                || mfd == 17
                || mfd == 18
                || (mfd >= 21 && mfd <= 26)
                || (mfd >= 31 && mfd <= 36)
            {
                lbl = 130;
            } else if iprint != 1 {
                lbl = 150;
            } else if mfd == 5 {
                lbl = 140;
            } else if nz > 1 {
                lbl = 120;
            }
            if lbl == 0 {
                if nl == 1 {
                    c.io.so(" enrgy  group constants at\n group  infinite dilution");
                }
                if nl > 1 {
                    c.io.so(" enrgy  lgend  group constants at\n group  order  infinite dilution");
                }
            } else if lbl == 120 {
                for i in 2..=nz {
                    field[i as usize] = a10(self.sigz[i]);
                }
                let f = fl(&field, 2, nz);
                if nl == 1 {
                    c.io.so(&rev_lines(
                        " enrgy  group constants vs sigma zero\n group  infinity ",
                        &f,
                        9,
                        &" ".repeat(17),
                        9,
                    ));
                }
                if nl > 1 {
                    c.io.so(&rev_lines(
                        " enrgy  lgend  group constants vs sigma zero\n group  order  infinity ",
                        &f,
                        9,
                        &" ".repeat(24),
                        9,
                    ));
                }
            } else if lbl == 130 {
                self.d_itr = 1;
                if iprint == 1 {
                    if nl == 1 && nz == 1 {
                        // 135
                        c.io.so(" initl  final  isotropic matrix vs final group\n group  group  +0         +1         +2         +3");
                    } else {
                        let mut lt = nl - 1;
                        if lt > 4 {
                            lt = 4;
                        }
                        if nz == 1 {
                            let mut s = String::from(" initl  final  group constants vs legendre order\n group  group  0");
                            for l in 1..=lt {
                                s.push_str(&" ".repeat(10));
                                s.push_str(&ffmt::i(l, 1));
                            }
                            c.io.so(&s);
                        }
                        if nz > 1 {
                            for i in 2..=nz {
                                field[i as usize] = a10(self.sigz[i]);
                            }
                            let f = fl(&field, 2, nz);
                            c.io.so(&revall(
                                " initl  final  lgend  group constants vs sigma zero\n group  group  order  infinity ",
                                &f,
                                9,
                            ));
                        }
                    }
                }
            } else if lbl == 140 {
                let nd = nl;
                if nd == 1 {
                    c.io.so(" normalized fission spectrum");
                }
                if nd > 1 {
                    for i in 1..=nd {
                        field[i as usize] = a10(self.dntc[i as usize]);
                    }
                    let f = fl(&field, 1, nd);
                    c.io.so(&revall(" normalized delayed neutron spectra vs time constant\n group", &f, 9));
                }
            }
            // 150
            c.io.so(" ");
            self.d_irat = 0;
            if mtd == 251 || mtd == 252 || mtd == 253 {
                self.d_irat = 1;
            }
            if mtd == 452 || mtd == 455 || mtd == 456 {
                self.d_irat = 1;
            }
            return;
        }

        //--write out results in appropriate format.
        let irat = self.d_irat;
        let itr = self.d_itr;
        if mfd == 5 {
            //--fission spectrum (mfd=5/mtd=18)
            if mtd != 455 {
                // prompt
                for j in 1..=*ng2 {
                    ans[(1, 1, j)] = sigfig(ans[(1, 1, j)], 7, 0);
                }
                *igzero = 1;
                if iprint != 1 {
                    return;
                }
                let mut j = 1;
                while j <= *ng2 {
                    let mut lim = j + 9;
                    if lim > *ng2 {
                        lim = *ng2;
                    }
                    let mut s = format!(" {}  ", ffmt::i(j, 5));
                    for i in j..=lim {
                        s.push_str(&a11(&a10(ans[(1, 1, i)])));
                    }
                    c.io.so(&s);
                    j += 10;
                }
            } else {
                //--delayed
                *igzero = 1;
                for j in 1..=*ng2 {
                    for il in 1..=nl {
                        ans[(il, 1, j)] = sigfig(ans[(il, 1, j)], 7, 0);
                    }
                }
                if iprint != 1 {
                    return;
                }
                for j in 2..=*ng2 {
                    let i2 = j - 1;
                    let mut s = format!(" {}  ", ffmt::i(i2, 5));
                    for il in 1..=nl {
                        s.push_str(&a11(&a10(ans[(il, 1, j)])));
                    }
                    c.io.so(&s);
                }
            }
            return;
        }
        if ig == 0 {
            //--constant spectra
            for i in 1..=*ng2 {
                ans[(1, 1, i)] = sigfig(ans[(1, 1, i)], 7, 0);
            }
            *igzero = 1;
            if iprint != 1 {
                return;
            }
            let mut i = 1;
            while i <= *ng2 {
                let mut lim = i + 5;
                if lim > *ng2 {
                    lim = *ng2;
                }
                let i2 = *ig2lo - 1 + i;
                let mut s = format!(" spec   {}  ", ffmt::i(i2, 5));
                for j in i..=lim {
                    s.push_str(&a11(&a10(ans[(1, 1, j)])));
                }
                c.io.so(&s);
                i += 6;
            }
            c.io.so(" ");
            return;
        }
        if *ig2lo == 0 {
            //--constant productions
            *igzero = 1;
            for i in 1..=*ng2 {
                ans[(1, 1, i)] = sigfig(ans[(1, 1, i)], 7, 0);
            }
            if iprint != 1 {
                return;
            }
            let mut i = 1;
            while i <= *ng2 {
                let mut lim = i + 5;
                if lim > *ng2 {
                    lim = *ng2;
                }
                let mut s = format!(" {}   prod  ", ffmt::i(i, 5));
                for j in i..=lim {
                    s.push_str(&a11(&a10(ans[(1, 1, j)])));
                }
                c.io.so(&s);
                i += 6;
            }
            c.io.so(" ");
            return;
        }
        if nz <= 1 && itr != 0 && nl != 1 {
            //--infinite dilution transfer matrices (iz=1).
            let mut ilo = 0;
            let mut ihi = *ng2;
            for ig2 in 2..=*ng2 {
                let igt = *ig2lo + ig2 - 2;
                let mut izero = 0;
                for il in 1..=nl {
                    result[il] = 0.0;
                    if ans[(il, 1, ig2)] != zero {
                        if ans[(il, 1, 1)] == zero {
                            ans[(il, 1, 1)] = big;
                        }
                        ans[(il, 1, 1)] = sigfig(ans[(il, 1, 1)], 7, 0);
                        result[il] = ans[(il, 1, ig2)] / ans[(il, 1, 1)];
                        result[il] = sigfig(result[il], 7, 0);
                        if result[il].abs() >= smin {
                            izero = 1;
                        }
                    }
                }
                if ilo == 0 && (izero > 0 || ig2 == *ng2) {
                    ilo = ig2;
                }
                if ilo != 0 {
                    for il in 1..=nl {
                        ans[(il, 1, ig2 - ilo + 2)] = result[il];
                    }
                    if izero > 0 {
                        ihi = ig2;
                    }
                    if izero == 1 {
                        *igzero = 1;
                    }
                    if iprint == 1 && izero > 0 {
                        let mut s = format!(" {}  {} ", ffmt::i(ig, 5), ffmt::i(igt, 5));
                        for i in 1..=nl {
                            s.push_str(&a11(&a10(result[i])));
                        }
                        c.io.so(&s);
                    }
                }
            }
            *ig2lo = *ig2lo + ilo - 2;
            *ng2 = ihi - ilo + 2;
            return;
        }
        if nz <= 1 && itr != 0 && nl == 1 {
            //--isotropic infinite dilution transfer matrix (il=1,iz=1).
            let mut ilo = 0;
            let mut ihi = *ng2;
            let mut j = 2;
            for ig2 in 2..=*ng2 {
                if ilo > 0 {
                    j += 1;
                }
                if ans[(1, 1, 1)] == zero {
                    ans[(1, 1, 1)] = big;
                }
                ans[(1, 1, 1)] = sigfig(ans[(1, 1, 1)], 7, 0);
                ans[(1, 1, j)] = ans[(1, 1, ig2)] / ans[(1, 1, 1)];
                ans[(1, 1, j)] = sigfig(ans[(1, 1, j)], 7, 0);
                if ans[(1, 1, j)].abs() < smin / 1000.0 {
                    ans[(1, 1, j)] = 0.0;
                }
                if ans[(1, 1, j)] != zero {
                    if ilo == 0 {
                        ilo = ig2;
                    }
                    ihi = ig2;
                    *igzero = 1;
                }
            }
            if *igzero == 0 {
                ilo = *ng2;
            }
            *ig2lo = *ig2lo + ilo - 2;
            *ng2 = ihi - ilo + 2;
            if *igzero == 0 {
                return;
            }
            if iprint != 1 {
                return;
            }
            let mut ig2 = 2;
            while ig2 <= *ng2 {
                let igt = *ig2lo + ig2 - 2;
                let mut igj = ig2 + 5;
                if igj > *ng2 {
                    igj = *ng2;
                }
                let mut s = format!(" {}  {} ", ffmt::i(ig, 5), ffmt::i(igt, 5));
                for j in ig2..=igj {
                    s.push_str(&a11(&a10(ans[(1, 1, j)])));
                }
                c.io.so(&s);
                ig2 += 6;
            }
            return;
        }

        //--self-shielded cross-sections and transfer matrices.
        // 310
        let mut ilo = 0;
        let mut ihi = *ng2;
        for ig2 in 2..=*ng2 {
            let igt = *ig2lo + ig2 - 2;
            let mut izero = 0;
            for il in 1..=nl {
                for iz in 1..=nz {
                    result[iz] = 0.0;
                    if ig2 == 2 {
                        ans[(il, iz, 1)] = sigfig(ans[(il, iz, 1)], 7, 0);
                        if irat == 1 {
                            ans[(il, iz, 3)] = sigfig(ans[(il, iz, 3)], 7, 0);
                        }
                    }
                    if ans[(il, iz, ig2)] != zero {
                        let mut locd = 1;
                        if irat == 1 && ig2 == 2 {
                            locd = 3;
                        }
                        if ans[(il, iz, locd)] == zero {
                            ans[(il, iz, locd)] = big;
                        }
                        result[iz] = ans[(il, iz, ig2)] / ans[(il, iz, locd)];
                        result[iz] = sigfig(result[iz], 7, 0);
                        izero = 1;
                    }
                }
                if ilo == 0 && (izero > 0 || ig2 == *ng2) {
                    ilo = ig2;
                }
                if ilo == 0 && irat == 1 {
                    ilo = ig2;
                }
                if ilo != 0 {
                    for iz in 1..=nz {
                        ans[(il, iz, ig2 - ilo + 2)] = result[iz];
                    }
                    if izero > 0 {
                        ihi = ig2;
                    }
                    if izero == 1 {
                        *igzero = 1;
                    }
                    if izero != 0 && iprint == 1 {
                        if irat != 1 || ig2 != 3 {
                            let l = il - 1;
                            if itr == 1 {
                                let f: Vec<String> = (1..=nz).map(|i| a10(result[i])).collect();
                                let head = format!(
                                    " {}  {}      {}  ",
                                    ffmt::i(ig, 5),
                                    ffmt::i(igt, 5),
                                    ffmt::i(l, 1)
                                );
                                c.io.so(&rev_lines(&head, &f, 10, &" ".repeat(22), 10));
                            } else if nl > 1 {
                                let f: Vec<String> = (1..=nz).map(|i| a10(result[i])).collect();
                                let head = format!(" {}      {}  ", ffmt::i(ig, 5), ffmt::i(l, 1));
                                c.io.so(&rev_lines(&head, &f, 10, &" ".repeat(15), 10));
                                if mfd == 3 && mtd == 1 {
                                    let f: Vec<String> = (1..=nz).map(|i| a10(ans[(il, i, 1)])).collect();
                                    let head = format!("     flux {}  ", ffmt::i(l, 1));
                                    c.io.so(&rev_lines(&head, &f, 10, &" ".repeat(13), 10));
                                }
                            } else {
                                let f: Vec<String> = (1..=nz).map(|i| a10(ans[(il, i, 2)])).collect();
                                let head = format!(" {}  ", ffmt::i(ig, 5));
                                c.io.so(&rev_lines(&head, &f, 10, &" ".repeat(8), 10));
                                if mfd == 3 && mtd == 1 {
                                    let f: Vec<String> = (1..=nz).map(|i| a10(ans[(1, i, 1)])).collect();
                                    c.io.so(&rev_lines("  flx ", &f, 10, &" ".repeat(6), 10));
                                }
                            }
                        }
                    }
                }
            }
        }
        *ig2lo = *ig2lo + ilo - 2;
        *ng2 = ihi - ilo + 2;
        let _ = &mut field;
    }

    /// Retrieve the weighting flux.
    pub(super) fn getflx(&mut self, c: &mut Ctx, e: f64, enext: &mut f64, idis: &mut i32, flux: &mut RM, nl: i32, nz: i32) {
        let emax = 1.0e10f64;
        let zero = 0.0f64;
        let mut wtf = 0.0f64;

        //--initialize.
        if e == zero {
            if self.nsigz > 1 {
                self.falo = alloc(self.nfv);
                self.fahi = alloc(self.nfv);
                c.finda(1, &mut self.falo.0[..], self.nfv, self.nflx, &mut self.wtbuf.0[..], NBUF);
                c.finda(2, &mut self.fahi.0[..], self.nfv, self.nflx, &mut self.wtbuf.0[..], NBUF);
                *enext = self.falo[1];
                self.gf_ip = 0;
            } else {
                self.getwtf(e, enext, idis, 0, &mut wtf, c);
            }
            return;
        }

        //--normal entry
        if self.nsigz > 1 {
            //--branch for nz>1 (self shielded)
            if self.gf_ip == 0 {
                c.scana(e, &mut self.gf_ip, self.nfp, self.nfv, self.nflx, &mut self.wtbuf.0[..], NBUF);
                c.finda(self.gf_ip, &mut self.falo.0[..], self.nfv, self.nflx, &mut self.wtbuf.0[..], NBUF);
                c.finda(self.gf_ip + 1, &mut self.fahi.0[..], self.nfv, self.nflx, &mut self.wtbuf.0[..], NBUF);
            }

            //--check energy range
            //--read next block of data if necessary
            while e >= self.fahi[1] && self.gf_ip < self.nfp {
                for i in 1..=self.nfv {
                    self.falo[i] = self.fahi[i];
                }
                self.gf_ip += 1;
                c.finda(self.gf_ip, &mut self.fahi.0[..], self.nfv, self.nflx, &mut self.wtbuf.0[..], NBUF);
            }

            //--in the right energy range
            //--interpolate for desired flux components.
            let mut llord = self.lord;
            if self.lord == 0 && self.nsigz > 1 {
                llord = 1;
            }
            for iz in 1..=nz {
                for il in 1..=nl {
                    let l = 1 + il + (llord + 1) * (iz - 1);
                    let mut y = 0.0;
                    c.terp1(self.falo[1], self.falo[l], self.fahi[1], self.fahi[l], e, &mut y, 2);
                    flux[(iz, il)] = y;
                }
            }
            let nfv = self.nfv;
            let mut y = 0.0;
            c.terp1(self.falo[1], self.falo[nfv], self.fahi[1], self.fahi[nfv], e, &mut y, 2);
            self.xtot = y;
            *enext = self.fahi[1];
            *idis = 0;
            if self.gf_ip == self.nfp {
                *enext = emax;
            }
        } else {
            //--branch for nz=1 (infinite dilution only)
            self.getwtf(e, enext, idis, 0, &mut wtf, c);
            for il in 1..=nl {
                flux[(1, il)] = wtf;
            }
        }
    }

    /// Retrieve or compute the yield for mat, mf, mt on itape. If
    /// mt=455, retrieve the delayed neutron time constants. Initialize
    /// if e=0.
    pub(super) fn getyld(
        &mut self,
        c: &mut Ctx,
        e: f64,
        enext: &mut f64,
        idis: &mut i32,
        yld: &mut f64,
        mat: i32,
        mf: i32,
        mt: i32,
        lfs: i32,
        itape: i32,
    ) {
        let emax = 1.0e10f64;
        let zero = 0.0f64;
        let (mut nb, mut nw) = (0i32, 0i32);

        if e > zero {
            //--tabulated data
            if self.gy_lnu == 1 {
                //--polynomial data
                let loc = 1;
                *yld = self.yield_[loc + 6];
                let mut term = 1.0;
                let nc = nint(self.yield_[loc + 4]);
                for i in 2..=nc {
                    term = term * e;
                    *yld = *yld + self.yield_[loc + 5 + i] * term;
                }
                *enext = emax;
                *idis = 0;
                return;
            }
            c.terpa(yld, e, enext, idis, &self.yield_.at(1), &mut self.gy_ip, &mut self.gy_ir);
            return;
        }

        //--initialize
        let ntmp = 20000;
        let mut tmp = RV::new(ntmp + SLACK);
        let mut mft = mf;
        if mft >= 40000000 {
            mft = 10;
        }
        if mft >= 30000000 {
            mft = 9;
        }
        if mft >= 20000000 {
            mft = 6;
        }
        c.findf(mat, mft, mt, itape);
        c.contio(itape, 0, 0, &mut tmp.at(1), &mut nb, &mut nw);
        let mut nk = 0;
        if mf > 99 {
            nk = nint(tmp[5]);
        }
        let mut ik = 0;
        self.gy_lnu = 0;
        if mt >= 452 && mt <= 456 {
            self.gy_lnu = nint(tmp[4]);
        }
        let mut iza = 0;
        if mf >= 10000000 {
            iza = (mf / 10) % 1000000;
        }
        let mut loc: i32 = 1;
        let na: i32;
        let mut lbl;
        if mt == 455 {
            // 110
            let mut ll = loc;
            c.listio(itape, 0, 0, &mut tmp.at(ll), &mut nb, &mut nw);
            while nb != 0 {
                ll += nw;
                c.moreio(itape, 0, 0, &mut tmp.at(ll), &mut nb, &mut nw);
            }
            let lnd = nint(tmp[loc + 4]);
            if lnd > 8 {
                c.error("getyld", "illegal lnd.", " ");
            }
            for i in 1..=lnd {
                self.dntc[i as usize] = tmp[i + 5 + loc];
            }
            lbl = 130;
        } else if self.gy_lnu != 1 {
            lbl = 130;
        } else {
            let mut ll = loc;
            c.listio(itape, 0, 0, &mut tmp.at(ll), &mut nb, &mut nw);
            let mut n = nw;
            while nb != 0 {
                ll += nw;
                c.moreio(itape, 0, 0, &mut tmp.at(ll), &mut nb, &mut nw);
                n += nw;
            }
            *enext = emax;
            na = n;
            // 190
            *idis = 0;
            self.yield_ = alloc(na);
            for i in 1..=na {
                self.yield_[i] = tmp[i];
            }
            return;
        }
        loop {
            match lbl {
                130 => {
                    ik += 1;
                    loc = 1;
                    c.tab1io(itape, 0, 0, &mut tmp.at(loc), &mut nb, &mut nw);
                    loc += nw;
                    while nb != 0 {
                        c.moreio(itape, 0, 0, &mut tmp.at(loc), &mut nb, &mut nw);
                        loc += nw;
                    }
                    if nk == 0 {
                        lbl = 180;
                        continue;
                    }
                    let mut izn = 0;
                    let mut lfn = 0;
                    if mft == 9 || mft == 10 {
                        izn = nint(tmp[3]);
                    }
                    if mft == 6 {
                        izn = nint(tmp[1]);
                    }
                    if mft == 9 || mft == 10 {
                        lfn = nint(tmp[4]);
                    }
                    if mft == 6 {
                        lfn = nint(tmp[3]);
                    }
                    if iza == izn && lfs == lfn {
                        lbl = 180;
                        continue;
                    }
                    if mft > 6 && izn == 0 && lfs == lfn {
                        lbl = 180;
                        continue;
                    }
                    if !(mft == 9 || mft == 10) {
                        self.law = nint(tmp[4]);
                        c.skip6(itape, 0, 0, &mut tmp.at(loc), self.law);
                    }
                    // 170
                    if ik < nk {
                        lbl = 130;
                        continue;
                    }
                    let s = format!("unable to find nuclide for iza={} lfs={}", ffmt::i(iza, 7), ffmt::i(lfs, 3));
                    c.error("getyld", &s, " ");
                }
                _ => {
                    // 180
                    let na = loc - 1;
                    self.gy_ip = 2;
                    self.gy_ir = 1;
                    let nr = nint(tmp[5]);
                    *enext = tmp[7 + 2 * nr];
                    // 190
                    *idis = 0;
                    self.yield_ = alloc(na);
                    for i in 1..=na {
                        self.yield_[i] = tmp[i];
                    }
                    return;
                }
            }
        }
    }

    /// Retrieve the reaction cross section defined by mfd and mtd.
    /// Use unresolved cross sections if needed. Initialize if e=0.
    pub(super) fn getsig(&mut self, c: &mut Ctx, e: f64, enext: &mut f64, idis: &mut i32, sig: &mut RM, nl: i32, nz: i32) {
        let emin = 1.0e-5f64;
        let step = 1.01e0f64;
        let vc = 1.919e8f64;
        let ezero = 1.0e7f64;
        let step2 = 1.10e0f64;
        let zero = 0.0f64;
        let mut nb = 0i32;
        let mut nw: i32;
        let mut s = 0.0f64;
        let mfd = self.mfd;
        let mtd = self.mtd;

        //--initialize
        if e == zero {
            self.gs_nsig = self.npend;
            if mfd == 5 {
                // spectra
                *enext = emin;
                self.thresh = 0.0;
                self.lrflag = 0;
                self.sigma = alloc(1);
                return;
            }
            let mut mf = 3;
            if mfd == 13 || mfd == 17 {
                mf = 13;
            }
            if mfd >= 10000000 {
                mf = 3;
            }
            if mfd >= 20000000 {
                mf = 3;
            }
            if mfd >= 30000000 {
                mf = 3;
            }
            if mfd >= 40000000 {
                mf = 10;
            }
            let mut mt = 0;
            if mtd == 259 || mtd == 258 || mtd == 257 {
                // non cross section quantities
                *enext = emin;
                self.thresh = 0.0;
                self.lrflag = 0;
                self.sigma = alloc(1);
                self.gs_mt = mt;
                return;
            }
            let iverf = c.iverf;
            if mtd <= 150 {
                mt = mtd;
            }
            if mtd >= 152 && mtd <= 200 {
                mt = mtd;
            }
            if mtd >= 201 && mtd <= 250 {
                mt = mtd;
            }
            if mtd > 299 && mtd < 451 {
                mt = mtd;
            }
            if iverf >= 6 && mtd == 500 {
                mt = mtd;
            }
            if iverf < 6 && mtd >= 700 && mtd <= 799 {
                mt = mtd;
            }
            if iverf >= 6 && mtd >= 600 && mtd <= 849 {
                mt = mtd;
            }
            if iverf >= 6 && mtd >= 850 && mtd <= 874 {
                mt = mtd;
            }
            if iverf >= 6 && mtd >= 875 && mtd <= 891 {
                mt = mtd;
            }
            if iverf >= 6 && mtd == 599 {
                mt = mtd;
            }
            if mtd == 251 || mtd == 252 || mtd == 253 {
                mt = 2;
            }
            if mtd == 261 {
                mt = 261;
            }
            if mtd == 452 || mtd == 455 || mtd == 456 {
                mt = 18;
            }
            self.gs_mt = mt;
            if mt == 0 {
                let s = format!("{} invalid in endf", ffmt::i(mtd, 4));
                c.error("getsig", "illegal mt.", &s);
            }
            let nsig = self.gs_nsig;
            c.findf(self.matd, mf, mt, nsig);
            nw = c.npage + 50;
            self.sigma = RV::new(nw as usize + SLACK);
            c.contio(nsig, 0, 0, &mut self.sigma.at(1), &mut nb, &mut nw);
            self.awr = self.sigma[2];
            if self.awrp != zero {
                self.awr = self.awr / self.awrp;
            }
            if mf == 10 {
                let nfs = c.n1h;
                let mut jfs = -1;
                for i in 1..=nfs {
                    c.tab1io(nsig, 0, 0, &mut self.sigma.at(1), &mut nb, &mut nw);
                    if c.l1h == self.izar && c.l2h == self.lfs {
                        jfs = i;
                    }
                    while nb != 0 {
                        c.moreio(nsig, 0, 0, &mut self.sigma.at(1), &mut nb, &mut nw);
                    }
                }
                if jfs < 0 {
                    let s = format!(
                        "can't find mf,mt,izar,lfs = {}{}{}{}",
                        ffmt::i(mf, 9),
                        ffmt::i(mt, 9),
                        ffmt::i(self.izar, 9),
                        ffmt::i(self.lfs, 5)
                    );
                    c.error("getsig", &s, " ");
                }
                let nskip = jfs - 1;
                c.skiprz(nsig, -1);
                c.findf(self.matd, mf, mt, nsig);
                c.contio(nsig, 0, 0, &mut self.sigma.at(1), &mut nb, &mut nw);
                if nskip > 0 {
                    for _i in 1..=nskip {
                        c.tab1io(nsig, 0, 0, &mut self.sigma.at(1), &mut nb, &mut nw);
                        while nb != 0 {
                            c.moreio(nsig, 0, 0, &mut self.sigma.at(1), &mut nb, &mut nw);
                        }
                    }
                }
            }
            c.gety1(e, enext, idis, &mut s, nsig, &mut self.sigma.at(1));
            self.q = c.c2h;
            self.lrflag = c.l2h;
            if mf != 3 {
                self.lrflag = 0;
            }
            self.thresh = (self.awr + 1.0) * (-self.q) / self.awr;
            self.alpha = powc(self.awr - 1.0, 2) / powc(self.awr + 1.0, 2);
            return;
        }

        //--retrieve point cross sections.
        //--check for unresolved cross sections.
        if mtd == 259 {
            *idis = 0;
            *enext = step * e;
            sig[(1, 1)] = 1.0 / sqrt(vc * e);
        } else if mtd == 258 {
            *idis = 0;
            *enext = step * e;
            sig[(1, 1)] = log(ezero / e);
        } else if mtd == 257 {
            *idis = 0;
            *enext = step * e;
            sig[(1, 1)] = e;
        } else if mtd == 2 && self.izap > 1 {
            // cp elastic
            *idis = 0;
            *enext = step2 * e;
            for il in 1..=nl {
                for iz in 1..=nz {
                    sig[(iz, il)] = 1.0;
                }
            }
        } else {
            // other reactions
            let nsig = self.gs_nsig;
            c.gety1(e, enext, idis, &mut s, nsig, &mut self.sigma.at(1));
            let mut en = 0.0f64;
            for il in 1..=nl {
                for iz in 1..=nz {
                    sig[(iz, il)] = s;
                }
                let mut mtt = self.gs_mt;
                if il > 1 && self.gs_mt == 1 {
                    mtt = 261;
                }
                en = *enext;
                if nz > 1 {
                    let k = (il - 1) as usize * sig.n1;
                    self.getunr(c, mtt, e, &mut en, &mut F(&mut sig.v[k..]));
                }
            }
            if en < zero {
                return;
            }
            if en < *enext {
                *idis = 0;
            }
            if en < *enext {
                *enext = en;
            }
        }
    }

    /// Read self-shielded unresolved resonance cross sections from ntape
    /// and store them for retrieval with getunr.
    pub(super) fn stounr(&mut self, c: &mut Ctx, mat: i32, temp: f64, ntape: i32) {
        let (mut nb, mut nw) = (0i32, 0i32);
        //--initialize
        self.nunr = 0;
        if ntape == 0 {
            return;
        }
        //--locate mat and temp
        if self.lrp == -1 {
            self.unr = alloc(1);
            return;
        }
        let ntmp = 10000i32;
        let nsig0: i32;
        let mut l: i32;
        let mut tmp = RV::new(ntmp as usize + SLACK);
        loop {
            // 100
            let mfd = 2;
            let mtd = 152;
            c.findf(mat, 2, 0, ntape);
            c.contio(ntape, 0, 0, &mut tmp.at(1), &mut nb, &mut nw);
            let mut found = false;
            if !(c.mfh > mfd) {
                if c.mth == mtd {
                    found = true;
                } else {
                    c.tosend(ntape, 0, 0, &mut tmp.at(1));
                    c.contio(ntape, 0, 0, &mut tmp.at(1), &mut nb, &mut nw);
                    if !(c.mfh > mfd) && c.mth == mtd {
                        found = true;
                    }
                }
            }
            if !found {
                // 105
                self.unr = alloc(1);
                return;
            }
            //--read in desired data
            // 110
            l = 1;
            self.intunr = c.n2h;
            self.lssf = c.l1h;
            c.listio(ntape, 0, 0, &mut tmp.at(l), &mut nb, &mut nw);
            let temz = tmp[l];
            if (temz - temp).abs() <= 1.0 {
                break;
            }
            if temz > temp {
                let s = format!("cannot find temp={}{}", ffmt::e1p(temp, 12, 3), ffmt::e1p(temz, 12, 3));
                c.error("stounr", &s, " ");
            }
            c.tomend(ntape, 0, 0, &mut tmp.at(l));
        }
        // 120
        nsig0 = nint(tmp[4]);
        self.nunr = nint(tmp[6]);
        l += nw;
        while nb != 0 {
            if l > ntmp {
                c.error("stounr", "storage exceeded.", " ");
            }
            c.moreio(ntape, 0, 0, &mut tmp.at(l), &mut nb, &mut nw);
            l += nw;
        }
        if nsig0 > ntmp {
            c.error("stounr", "storage exceeded.", " ");
        }
        self.unr = alloc(l);
        for i in 1..=l {
            self.unr[i] = tmp[i];
        }

        //--check for sigmazero conflicts
        if self.nsigz == 1 {
            return;
        }
        if nsig0 <= 1 {
            c.mess("stounr", "no unresolved sigma zero data", "will use infinite dilution only");
            return;
        }
        for i in 1..=self.nsigz {
            let mut ok = false;
            for j in 1..=nsig0 {
                if (self.sigz[i] - self.unr[6 + j]).abs() < self.sigz[i] / 100.0 {
                    ok = true;
                    break;
                }
            }
            if !ok {
                c.mess("stounr", "sigma zero grids do not match", "will interpolate for desired values");
                return;
            }
        }
    }

    /// Retrieve unresolved resonance cross sections stored by stounr and
    /// correct sig for self-shielding. enext=-1 when e is outside the
    /// unresolved range or there are no unresolved data.
    pub(super) fn getunr(&mut self, c: &mut Ctx, mt: i32, e: f64, enext: &mut f64, sig: &mut F) {
        let eps = 1.0e-4f64;
        let del = 1.0e-10f64;
        let zero = 0.0f64;

        //--check nunr and energy range
        *enext = -1.0;
        let mut iovl = 0;
        if self.nunr == 0 {
            return;
        }
        let unr = &self.unr;
        let nunr = self.nunr;
        let nx = nint(unr[3]);
        let nsig0 = nint(unr[4]);
        if nsig0 <= 1 {
            return;
        }
        let ncyc = 1 + nx * nsig0;
        let ibase = 7 + nsig0;
        let e1 = unr[ibase].abs();
        if e < e1 - del * e1 {
            return;
        }
        let enunr = unr[ibase + (nunr - 1) * ncyc].abs();
        if e > enunr + del * enunr {
            return;
        }

        //--determine reaction type
        let mut ix = 0;
        if mt == 1 {
            ix = 1;
        }
        if mt == 2 {
            ix = 2;
        }
        if mt == 18 {
            ix = 3;
        }
        if mt == 19 {
            ix = 3;
        }
        if mt == 102 {
            ix = 4;
        }
        if mt == 261 {
            ix = 5;
        }
        if ix == 0 {
            return;
        }
        let locl = ibase + (nunr - 1) * ncyc + nsig0 * (ix - 1);
        let locs = 7;

        //--find interpolation range
        let lim = nunr - 1;
        let mut iflag = 0;
        let mut i = 0;
        let (mut loc, mut loc1) = (0, 0);
        while iflag == 0 && i < lim {
            i += 1;
            loc = ibase + (i - 1) * ncyc;
            loc1 = loc + ncyc;
            let mut elo = unr[loc].abs();
            elo = elo - del * elo;
            let ehi = unr[loc1].abs();
            if e > elo && e < ehi {
                iflag = 1;
            }
        }
        let sl_ = |k: i32| &unr.0[(k - 1) as usize..];
        let mut s = 0.0f64;

        if iflag == 0 {
            //--special exit for last point
            let sinf = unr[locl + 1];
            if unr[locl] < zero {
                iovl = 1;
            }
            if ix == 1 && iovl == 1 {
                self.xtot = sig[1] - sinf;
            }
            for is in 1..=self.nsigz {
                let mut sigb = self.sigz[is];
                if iovl == 1 {
                    sigb = sigb + self.xtot;
                }
                terpu(&mut s, sigb, sl_(locl + 1), sl_(locs), nsig0);
                if self.lssf == 0 {
                    sig[is] = sig[is] + s - sinf;
                }
                if self.lssf > 0 && sinf != zero {
                    sig[is] = sig[is] * s / sinf;
                }
            }
            *enext = enunr * (1.0 + eps);
        } else {
            //--interpolate for cross-sections requested
            let el = unr[loc].abs();
            let en = unr[loc1].abs();
            if unr[loc] < zero && unr[loc1] < zero {
                iovl = 1;
            }
            let ib = nsig0 * (ix - 1) + 1;
            let mut sinf = 0.0f64;
            c.terp1(el, unr[loc + ib], en, unr[loc1 + ib], e, &mut sinf, self.intunr);
            if ix == 1 && iovl == 1 {
                self.xtot = sig[1] - sinf;
            }
            for is in 1..=self.nsigz {
                let mut sigb = self.sigz[is];
                if iovl == 1 {
                    sigb = sigb + self.xtot;
                }
                let (mut sl, mut sn) = (0.0f64, 0.0f64);
                terpu(&mut sl, sigb, sl_(loc + ib), sl_(locs), nsig0);
                terpu(&mut sn, sigb, sl_(loc1 + ib), sl_(locs), nsig0);
                if sl < zero || sn < zero {
                    c.mess("getunr", " Warning, negative URR cross sections found, check unresr", " ");
                }
                c.terp1(el, sl, en, sn, e, &mut s, self.intunr);
                if self.lssf == 0 {
                    sig[is] = sig[is] + s - sinf;
                }
                if self.lssf > 0 && sinf != zero {
                    sig[is] = sig[is] * s / sinf;
                }
            }
            *enext = en;
        }
    }

    /// Compute the feed function or yield for the desired reaction type.
    pub(super) fn getff(
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
        nlg: i32,
        _ng2g: i32,
    ) {
        let emax = 1.0e10f64;
        let small = 1.0e-10f64;
        let smin = 1.0e-9f64;
        let zero = 0.0f64;
        let mfd = self.mfd;
        let mtd = self.mtd;
        let iverf = c.iverf;
        let mut en: f64;
        let mut idis = 0i32;

        //--select option according to reaction type requested
        let opt: i32 = if mtd == 251 || mtd == 252 || mtd == 253 {
            400
        } else if self.izap > 1 && mfd == 3 && mtd == 2 {
            800
        } else if mfd == 3 {
            100
        } else if mfd == 5 {
            200
        } else if mfd == 12 || mfd == 13 {
            100
        } else if mfd >= 10000000 {
            100
        } else if mfd == 6 && mtd == 2 {
            400
        } else if mfd == 6 && mtd >= 221 && mtd <= 250 {
            700
        } else if mfd == 6 && (mtd >= 51 && mtd <= 90) {
            400
        } else if mfd == 6 && (mtd >= 6 && mtd <= 9) {
            400
        } else if mfd == 6 {
            200
        } else if mfd == 16 || mfd == 17 {
            300
        } else if mfd == 8 {
            800
        } else if mfd == 18 {
            810
        } else if mfd >= 21 && mfd <= 26 {
            800
        } else if mfd >= 31 && mfd <= 36 {
            let r = |a: i32, b: i32| mtd >= a && mtd < b;
            if mtd == 2 || r(51, 91) {
                400
            } else if iverf < 6 && (r(700, 719) || r(720, 739) || r(740, 759) || r(760, 779) || r(780, 799)) {
                400
            } else if iverf >= 6
                && (r(600, 649) || r(650, 699) || r(700, 749) || r(750, 799) || r(800, 849) || r(875, 891))
            {
                400
            } else {
                0
            }
        } else {
            0
        };
        if opt == 0 {
            let s = format!("do not know how to handle mf,mt: {},{}", ffmt::i(mfd, 2), ffmt::i(mtd, 3));
            c.error("getff", &s, " ");
        }

        match opt {
            100 => {
                //--self-shielded group cross sections and ratios.
                *enext = emax;
                *idisc = 0;
                en = *enext;
                *nq = 0;
                *ng = 1;
                if mtd == 452 || mtd == 455 || mtd == 456 {
                    *ng = 2;
                }
                *iglo = 1;
                *yld = 1.0;
                if mfd == 12 || (mfd >= 30000000 && mfd < 40000000) {
                    self.getyld(c, e, &mut en, &mut idis, yld, self.matd, mfd, mtd, self.lfs, self.nend3);
                }
                if mfd >= 20000000 && mfd < 30000000 {
                    self.getyld(c, e, &mut en, &mut idis, yld, self.matd, mfd, mtd, self.isom, self.nend3);
                }
                if mtd == 452 || mtd == 455 || mtd == 456 {
                    self.getyld(c, e, &mut en, &mut idis, yld, self.matd, 1, mtd, 0, self.nend3);
                }
                if en < *enext * (1.0 - small) {
                    *idisc = idis;
                }
                if en < *enext * (1.0 - small) {
                    *enext = en;
                }
                if e == zero {
                    return;
                }
                for il in 1..=*nl {
                    ff[(il, 1)] = *yld;
                    if *ng == 2 {
                        ff[(il, 2)] = 1.0;
                    }
                }
            }
            200 => self.getff_cont(c, e, enext, idisc, yld, ff, nl, ng, iglo, nq),
            300 => self.getff_gam(c, e, enext, idisc, yld, ff, nl, ng, iglo, nq, nlg),
            400 => {
                //--discrete channel scattering by the center-of-mass method
                self.parts(c, mfd, mtd);
                self.getdis(c, e, enext, idisc, yld, ff, nl, ng, iglo, nq, self.nend2, nlg);
            }
            700 => {
                //--angle-energy distribution matrices (mf6), thermal data
                *idisc = 0;
                self.getaed(c, ff, e, enext, idisc, nl, nq, self.matd, mfd, mtd, self.npend2, nlg);
                *yld = 1.0;
                *ng = self.ngn;
                if e == zero {
                    return;
                }
                *iglo = self.ngn;
                let mut ighi = 1;
                for ig in 1..=self.ngn {
                    for il in 1..=*nl {
                        if ff[(il, ig)].abs() > small {
                            if ig < *iglo {
                                *iglo = ig;
                            }
                            if ig > ighi {
                                ighi = ig;
                            }
                        }
                    }
                }
                if *iglo > ighi {
                    *iglo = ighi;
                }
                *iglo = 1;
                if *iglo < 1 {
                    *iglo = 1;
                }
                *ng = ighi - *iglo + 1;
            }
            800 | 810 => {
                //--energy-angle distributions from endf-6 file 6
                let (which, ngx) = if opt == 800 { (1, self.ngn) } else { (2, self.ngg) };
                self.getmf6(c, ff, e, enext, idisc, yld, which, ngx, nl, iglo, ng, nq, self.matd, mfd, mtd, self.nend3, nlg);
                if opt == 800 && ((mtd >= 18 && mtd <= 21) || mtd == 38) {
                    en = 0.0;
                    self.getyld(c, e, &mut en, &mut idis, yld, self.matd, 1, 456, 0, self.nend3);
                    for ig in 1..=*ng {
                        for il in 1..=*nl {
                            ff[(il, ig)] = ff[(il, ig)] * *yld;
                        }
                    }
                    if en < *enext * (1.0 - small) {
                        *enext = en;
                        *idisc = idis;
                    }
                }
            }
            _ => unreachable!(),
        }
        let _ = smin;
    }

    /// getff, neutron continuum transfer matrices and spectra (label 200).
    fn getff_cont(
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
    ) {
        let emax = 1.0e10f64;
        let small = 1.0e-10f64;
        let smin = 1.0e-9f64;
        let zero = 0.0f64;
        let mfd = self.mfd;
        let mtd = self.mtd;
        let mut idis = 0i32;
        let mut fle = [0.0f64; 22];
        *enext = emax;
        *idisc = 0;
        let mut en = *enext;
        *nq = 0;
        *yld = 1.0;
        if e == zero {
            self.ff_ifirst = 0;
        }
        if e > zero {
            self.ff_ifirst += 1;
        }
        let ifirst = self.ff_ifirst;
        let mut to202 = mfd == 5;
        if !to202 && e < self.econst * (1.0 - small) && ifirst == 1 {
            to202 = true;
        }
        if !to202 {
            if mtd == 16 || mtd == 24 || mtd == 26 || mtd == 30 {
                *yld = 2.0;
            }
            if mtd == 11 || mtd == 41 {
                *yld = 2.0;
            }
            if mtd == 17 || mtd == 25 {
                *yld = 3.0;
            }
            if mtd == 42 {
                *yld = 3.0;
            }
            if mtd == 37 {
                *yld = 4.0;
            }
            if (mtd >= 18 && mtd <= 21) || mtd == 38 {
                self.getyld(c, e, &mut en, &mut idis, yld, self.matd, 1, 456, 0, self.nend3);
            }
            if en < *enext * (1.0 - small) {
                *idisc = idis;
            }
            if en < *enext * (1.0 - small) {
                *enext = en;
            }
            let mut nle = *nl;
            let mut lcd = 1;
            self.getfle(c, e, &mut en, &mut idis, &mut F(&mut fle[1..]), &mut nle, &mut lcd, self.matd, 4, mtd, self.nend2);
            if e == 0.0 && nle == 1 {
                *nl = 1;
            }
            Self::merge3(en, idis, enext, idisc, small);
        }
        // 202
        let mut lbl;
        let mut nk = 1;
        if e < self.econst * (1.0 - small) && ifirst > 1 {
            lbl = 230;
        } else {
            if mtd == 455 {
                nk = self.ndelg;
            }
            if e == zero {
                self.sede = alloc(nk * self.ngn);
                self.sede_ld = nk as usize;
            }
            let mut sed = std::mem::take(&mut self.sede);
            let ld = self.sede_ld;
            self.getsed(c, e, &mut en, &mut idis, &mut V2::new(&mut sed.0, ld), 1, self.ngn, &mut nk, self.matd, 5, mtd, self.nendf);
            self.sede = sed;
            Self::merge3(en, idis, enext, idisc, small);
            if !(e > zero) {
                if self.econst == zero {
                    return;
                }
                // find group where e-dependence begins
                let ngn1 = self.ngn + 1;
                for ig in 1..=ngn1 {
                    if self.econst > self.egn[ig] * (1.0 + small) {
                        self.jconst = ig - 1;
                    }
                }
                self.econst = self.egn[self.jconst + 1];
                if self.jconst <= 2 {
                    self.econst = 0.0;
                }
                return;
            }
            // 205: for e.gt.econst, return matrix
            if mfd == 5 {
                lbl = 260;
            } else if e < self.econst * (1.0 - small) {
                lbl = 230;
            } else {
                *iglo = 1;
                let mut ighi = 1;
                let ld = self.sede_ld;
                for ig in 1..=self.ngn {
                    let temp = *yld * self.sede[1 + (ig - 1) * ld as i32];
                    if temp > 0.0 {
                        ighi = ig;
                    }
                    for il in 1..=*nl {
                        ff[(il, ig)] = temp * fle[il as usize];
                    }
                }
                *ng = ighi - *iglo + 1;
                return;
            }
        }
        if lbl == 230 {
            // for e.lt.econst, return yield only
            if ifirst != 1 {
                *iglo = 1;
                *ng = 1;
                ff[(1, 1)] = *yld;
                *enext = self.econst;
                *idisc = 1;
                return;
            }
            lbl = 260;
        }
        let _ = lbl;
        // 260: spectrum calculation (e.lt.econst or mfd=5)
        *iglo = 0;
        if mfd == 5 {
            *iglo = 1;
        }
        let mut ighi = 0i32;
        let ld = self.sede_ld as i32;
        let sd = |s: &RV, ik: i32, ig: i32| s[ik + (ig - 1) * ld];
        for ig in 1..=self.ngn {
            let mut go = *iglo > 0;
            if !go {
                for ik in 1..=nk {
                    if sd(&self.sede, ik, ig) > smin / 1000.0 || mtd == 18 {
                        go = true;
                        break;
                    }
                }
            }
            if !go {
                continue;
            }
            // 280
            if *iglo == 0 {
                *iglo = ig;
            }
            if mfd != 5 {
                ighi = ig;
            }
            for ik in 1..=nk {
                let jg = ig - *iglo + 1;
                ff[(ik, jg)] = sd(&self.sede, ik, ig);
                if mfd == 5 {
                    if ff[(ik, jg)] != zero {
                        ighi = ig;
                    }
                }
            }
        }
        if *iglo == 0 {
            ighi = 1;
        }
        if *iglo == 0 {
            *iglo = 1;
        }
        *ng = ighi - *iglo + 1;
    }

    /// Take en as the next point if it is lower (getff and others).
    fn merge3(en: f64, idis: i32, enext: &mut f64, idisc: &mut i32, small: f64) {
        if (en - *enext).abs() < *enext * small && idis > *idisc {
            *idisc = idis;
        }
        if en < *enext * (1.0 - small) {
            *idisc = idis;
        }
        if en < *enext * (1.0 - small) {
            *enext = en;
        }
    }

    /// getff, gamma production matrices (label 300).
    fn getff_gam(
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
        nlg: i32,
    ) {
        let emax = 1.0e10f64;
        let small = 1.0e-10f64;
        let zero = 0.0f64;
        let mfd = self.mfd;
        let mtd = self.mtd;
        let mut idis = 0i32;
        let mut en = 0.0f64;
        *enext = emax;
        *idisc = 0;
        *nq = 0;
        let mut mft = 0;
        if mfd == 16 {
            mft = 12;
        }
        if mfd == 17 {
            mft = 13;
        }
        // initialize retrieval routines and assign storage.
        if e == zero {
            let mut y = [0.0f64; 8];
            let mut eccn = [0.0f64; 8];
            let mut nyl = self.ff_nyl;
            self.getgyl(c, e, &mut en, &mut idis, &mut F(&mut y[..]), &mut F(&mut eccn[..]), &mut nyl, self.matd, mft, mtd, self.nend2);
            self.ff_nyl = nyl;
            let econ = eccn[0];
            if en < *enext * (1.0 - small) {
                *enext = en;
            }
            self.gyle = alloc(nyl);
            self.eyle = alloc(nyl);
            self.gfle = alloc(nlg * nyl);
            self.gfle_ld = nlg as usize;
            let mut gf = std::mem::take(&mut self.gfle);
            self.getgfl(c, e, &mut en, &mut idis, &mut V2::new(&mut gf.0, nlg as usize), nl, nlg, nyl, self.matd, 14, mtd, self.nend2);
            self.gfle = gf;
            if econ == zero {
                let mut nk = 1;
                let mut s = [0.0f64; 8];
                self.getsed(c, e, &mut en, &mut idis, &mut V2::new(&mut s[..], 1), 2, self.ngg, &mut nk, self.matd, 15, mtd, self.nendf);
                if en < *enext * (1.0 - small) {
                    *enext = en;
                }
                self.sede = alloc(nk * self.ngg);
                self.sede_ld = nk as usize;
            }
            // find group containing lowest energy photon
            self.ff_igmin = 1;
            while econ >= self.egg[self.ff_igmin + 1] {
                self.ff_igmin += 1;
            }
            if self.econst != zero {
                // find group where e dependence starts
                let ngn1 = self.ngn + 1;
                for ig in 1..=ngn1 {
                    if self.econst > self.egn[ig] {
                        self.jconst = ig - 1;
                    }
                }
                self.econst = self.egn[self.jconst + 1];
                if self.jconst <= 2 {
                    self.econst = 0.0;
                }
            }
            self.ff_ifirst = 0;
            return;
        }
        // retrieve factors of feed function.
        self.ff_ifirst += 1;
        let ifirst = self.ff_ifirst;
        let nyl0 = self.ff_nyl;
        let mut nyl = nyl0;
        let mut gy = std::mem::take(&mut self.gyle);
        let mut ey = std::mem::take(&mut self.eyle);
        self.getgyl(c, e, &mut en, &mut idis, &mut gy.at(1), &mut ey.at(1), &mut nyl, self.matd, mft, mtd, self.nend2);
        self.gyle = gy;
        self.eyle = ey;
        self.ff_nyl = nyl;
        Self::merge3(en, idis, enext, idisc, small);
        if e < self.econst * (1.0 - small) && ifirst > 1 {
            // 395
            *iglo = -1;
            *ng = 1;
            *yld = 0.0;
            for ik in 1..=nyl {
                *yld = *yld + self.gyle[ik];
            }
            ff[(1, 1)] = *yld;
            *enext = self.econst;
            *idisc = 1;
            return;
        }
        let mut gf = std::mem::take(&mut self.gfle);
        self.getgfl(c, e, &mut en, &mut idis, &mut V2::new(&mut gf.0, nlg as usize), nl, nlg, nyl, self.matd, 14, mtd, self.nend2);
        self.gfle = gf;
        Self::merge3(en, idis, enext, idisc, small);
        if self.eyle[nyl] == zero {
            let mut nk = 1;
            let mut sed = std::mem::take(&mut self.sede);
            self.getsed(c, e, &mut en, &mut idis, &mut V2::new(&mut sed.0, 1), 2, self.ngg, &mut nk, self.matd, 15, mtd, self.nendf);
            self.sede = sed;
            Self::merge3(en, idis, enext, idisc, small);
            if !(e > zero) {
                return;
            }
        }
        // 330
        *iglo = self.ff_igmin;
        let mut ighi = 1;
        let gld = self.gfle_ld;
        for ig in 1..=self.ngg {
            let epl = self.egg[ig];
            let eph = self.egg[ig + 1];
            for il in 1..=*nl {
                ff[(il, ig)] = 0.0;
            }
            for ik in 1..=nyl {
                if self.gyle[ik] == zero {
                    continue;
                }
                if self.eyle[ik] != zero {
                    if self.eyle[ik] < epl * (1.0 - small) {
                        continue;
                    }
                    if self.eyle[ik] >= eph * (1.0 - small) {
                        continue;
                    }
                }
                // 340
                let mut tempo = self.gyle[ik];
                if self.eyle[ik] == zero {
                    tempo = tempo * self.sede[ig];
                }
                for il in 1..=*nl {
                    ff[(il, ig)] = ff[(il, ig)] + tempo * self.gfle[il + (ik - 1) * gld as i32];
                    if ig > ighi {
                        ighi = ig;
                    }
                }
            }
        }
        *yld = 0.0;
        for ik in 1..=nyl {
            *yld = self.gyle[ik] + *yld;
        }
        *ng = ighi - *iglo + 1;
        if e >= self.econst || ifirst != 1 {
            if *iglo == 1 {
                return;
            }
            let iglo1 = *iglo - 1;
            for ig in 1..=*ng {
                for il in 1..=*nl {
                    ff[(il, ig)] = ff[(il, ig + iglo1)];
                }
            }
            return;
        }
        *iglo = 0;
        ighi = 0;
        for ig in 1..=self.ngg {
            ff[(1, ig)] = ff[(1, ig)] / *yld;
            if ff[(1, ig)] >= small {
                if *iglo == 0 {
                    *iglo = ig;
                }
                ighi = ig;
            }
        }
        *ng = ighi - *iglo + 1;
        for ig in 1..=*ng {
            ff[(1, ig)] = ff[(1, *iglo - 1 + ig)];
        }
    }

    /// Set up particles for reactions that use file 4.
    pub(super) fn parts(&mut self, c: &Ctx, mfd: i32, mtd: i32) {
        self.law = 0;
        self.zap = 1.0;
        self.aprime = 1.0;
        let awr = self.awr;
        let set = |s: &mut Self, zap: f64, ratio: f64, mfp: i32| {
            s.zap = zap;
            s.aprime = ratio;
            if mfd != mfp {
                s.aprime = awr + 1.0 - s.aprime;
                s.law = 4;
            }
        };
        if c.iverf < 6 {
            if mfd == 3 {
                return;
            }
            if mtd == 2 || (mtd >= 51 && mtd <= 91) {
                self.zap = 1.0;
                self.aprime = 1.0;
                if mfd != 6 {
                    self.aprime = awr;
                    self.law = 4;
                }
            } else if mtd >= 700 && mtd < 720 {
                set(self, 1001.0, PNRATIO, 31);
            } else if mtd >= 720 && mtd < 740 {
                set(self, 1002.0, DNRATIO, 32);
            } else if mtd >= 740 && mtd < 760 {
                set(self, 1003.0, TNRATIO, 33);
            } else if mtd >= 760 && mtd < 780 {
                set(self, 2003.0, HNRATIO, 34);
            } else if mtd >= 780 && mtd < 800 {
                set(self, 2004.0, ANRATIO, 35);
            }
        } else {
            if mfd == 3 && mtd != 2 {
                return;
            }
            if mtd == 2 {
                self.zap = self.izap as f64;
                self.aprime = 1.0;
                if mfd != 6 {
                    self.aprime = awr;
                    self.law = 4;
                }
            } else if mtd >= 50 && mtd <= 91 {
                self.zap = 1.0;
                self.aprime = 1.0;
                if mfd != 6 {
                    self.aprime = awr;
                    self.law = 4;
                }
            } else if mtd >= 600 && mtd < 650 {
                set(self, 1001.0, PNRATIO, 31);
            } else if mtd >= 650 && mtd < 700 {
                set(self, 1002.0, DNRATIO, 32);
            } else if mtd >= 700 && mtd < 750 {
                set(self, 1003.0, TNRATIO, 33);
            } else if mtd >= 750 && mtd < 800 {
                set(self, 2003.0, HNRATIO, 34);
            } else if mtd >= 800 && mtd < 850 {
                set(self, 2004.0, ANRATIO, 35);
            }
        }
    }
}

/// Interpolate for the cross section at a desired sigma zero value for
/// unresolved data. xs and sigs are 0-based slices of the Fortran
/// arrays xs(ns) and sigs(ns).
pub(super) fn terpu(sig: &mut f64, sigz: f64, xs: &[f64], sigs: &[f64], ns: i32) {
    let rlo = 0.99e0f64;
    let rhi = 10.1e0f64;
    let sg = |i: i32| sigs[(i - 1) as usize];
    *sig = 0.0;
    for is in 1..=ns {
        let mut terp = 1.0;
        let sm = sg(is);
        if sigz < sg(ns) && sm > rhi * sg(ns) {
            continue;
        }
        if sigz < sg(2) && sm > rhi * sg(2) {
            continue;
        }
        if sigz >= sg(2) && sm < rlo * sg(2) {
            continue;
        }
        for js in 1..=ns {
            let st = sg(js);
            if (sm - st).abs() < sm / 1000.0 {
                continue;
            }
            if sigz >= sg(ns) {
                // 130
                if sigz >= sg(2) {
                    // 140
                    if st < rlo * sg(2) {
                        continue;
                    }
                    terp = terp * ((st / sigz) - 1.0) / ((st / sm) - 1.0);
                    continue;
                }
                if st > rhi * sg(2) {
                    continue;
                }
                terp = terp * log(sigz / st) / log(sm / st);
                continue;
            }
            if st > rhi * sg(ns) {
                continue;
            }
            terp = terp * (powc(sigz, 2) - powc(st, 2)) / (powc(sm, 2) - powc(st, 2));
        }
        *sig = *sig + terp * xs[(is - 1) as usize];
    }
}
