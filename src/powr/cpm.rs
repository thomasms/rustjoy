//! POWR lib=3: the EPRI-CPM library (CLIB input), part 1: input,
//! pinit, subfile 2 (resonance data) and subfile 3 (resonance
//! integrals).

use super::{alloc_iv, alloc_rv, dealloc, geti, getr, pe, rowsj, Powr, SLACK};
use crate::endf::Ctx;
use crate::farr::{IV, RV};
use crate::ffmt;
use crate::fio::{fatal_runtime, Kind, RecW};
use crate::fmath::{log, nint, powc};

const ZERO: f64 = 0.0;
pub(super) const EFLAG: f64 = 1.0e35;

/// A 1-based real array whose reads below index 1 give 0 (the Fortran
/// reads memory before the array there).
fn rd0(v: &RV, i: i32) -> f64 {
    if i < 1 { 0.0 } else { v[i] }
}

impl Powr {
    /// Format multigroup cross sections from groupr for the EPRI-CPM
    /// library code CLIB.
    pub(super) fn cpm(&mut self, c: &mut Ctx, iprnt: i32, iclaps: i32) {
        let nrg69 = 13;
        let nrg185 = 15;
        self.ngnd = 69;
        self.nfg = 14;
        self.nrg = nrg185;
        if iclaps > 0 {
            self.nrg = nrg69;
        }
        self.nfsg = 22;
        self.ngref = 60;
        if self.nrg == nrg69 {
            self.ngref = 55;
        }
        let nt = 6;
        self.nzmax = 10;
        self.iprint = iprnt;
        // arrays im and rleth are available for temporary storage
        let immax = self.ngnd + 1 + self.nrg + 1;
        let iscr = 1;
        let mut im = IV::new(70 + SLACK);
        let mut z = RV::new(70 + SLACK);
        let mut km = IV::new(70 + SLACK);

        //--read and write user input.
        c.io.so("\n ***epri-cpm library***");
        c.io.se("\n ***epri-cpm library***");
        self.mode = 0;
        self.if5 = 0;
        self.if4 = 1;
        let (mut nlib, mut idatt, mut iopt) = (0i32, 0i32, 0i32);
        let v = c.io.input.read(&[Kind::Int; 7]);
        nlib = geti(&v[0], nlib);
        idatt = geti(&v[1], idatt);
        self.newmat = geti(&v[2], self.newmat);
        iopt = geti(&v[3], iopt);
        self.mode = geti(&v[4], self.mode);
        self.if5 = geti(&v[5], self.if5);
        self.if4 = geti(&v[6], self.if4);
        let newmat = self.newmat;
        for i in 1..=newmat {
            im[i] = 0;
        }
        if iopt != 1 {
            if newmat > immax {
                c.error("cpm", "exceeded storage in input array z.", " ");
            }
            let kinds = vec![Kind::Int; newmat.max(0) as usize];
            let v = c.io.input.read(&kinds);
            for i in 1..=newmat {
                im[i] = geti(&v[(i - 1) as usize], im[i]);
            }
        }
        c.io.so(&format!(
            " number of library .................... {}\n \
             date ................................. {}\n \
             no. of materials to be added ......... {}\n \
             add option (0 only mats read, 1 all) . {}\n \
             mode (0 replace, 1 add, 2 create ) ... {}\n \
             burn option (0 no, 1 yes, 2 only) .... {}\n \
             xsec data option (0 no, 1 yes) ....... {}",
            ffmt::i(nlib, 10),
            ffmt::i(idatt, 10),
            ffmt::i(newmat, 10),
            ffmt::i(iopt, 10),
            ffmt::i(self.mode, 10),
            ffmt::i(self.if5, 10),
            ffmt::i(self.if4, 10)
        ));
        if iopt == 0 {
            let mut s = String::from(" desired mats ......................... ");
            for ch in im.0[..(newmat.max(0) as usize)].chunks(12) {
                s.push_str("\n   ");
                for &x in ch {
                    s.push_str(&ffmt::i(x, 6));
                }
            }
            if newmat <= 0 {
                s.push_str("\n");
            }
            c.io.so(&s);
        }
        if iopt == 1 {
            c.io.so(&format!(" first {} mat(s) on ngendf will be processed.", ffmt::i(newmat, 2)));
        }
        c.io.so(
            "\n nina ntemp nsigz ref. sig0 ires pot. xsec mti mtc ip1 norf pos posr\n \
             ---- ----- ----- --------- ---- --------- --- --- --- ---- --- ----",
        );
        self.nres = 0;
        self.nnina = 0;
        let mut k = 0;
        let mut ntemp = 0i32;
        for _i in 1..=newmat {
            let mut nina = 0;
            ntemp = 0;
            let mut nsigz = 0;
            let mut sgref = 0.0;
            let mut ires = 0;
            let mut sigp = 0.0;
            let mut mti = 0;
            let mut mtc = 0;
            let mut ip1opt = 0;
            let mut norf = 0;
            let mut ipos = 0;
            let mut iposr = 0;
            let v = c.io.input.read(&[
                Kind::Int,
                Kind::Int,
                Kind::Int,
                Kind::Real,
                Kind::Int,
                Kind::Real,
                Kind::Int,
                Kind::Int,
                Kind::Int,
                Kind::Int,
                Kind::Int,
                Kind::Int,
            ]);
            nina = geti(&v[0], nina);
            ntemp = geti(&v[1], ntemp);
            nsigz = geti(&v[2], nsigz);
            sgref = getr(&v[3], sgref);
            ires = geti(&v[4], ires);
            sigp = getr(&v[5], sigp);
            mti = geti(&v[6], mti);
            mtc = geti(&v[7], mtc);
            ip1opt = geti(&v[8], ip1opt);
            norf = geti(&v[9], norf);
            ipos = geti(&v[10], ipos);
            iposr = geti(&v[11], iposr);
            if nina == 0 || nina == 3 {
                self.nnina += 1;
            }
            if ires == 1 {
                self.nres += 1;
            }
            if nina != 0 {
                ip1opt = 1;
            }
            c.io.so(&format!(
                "{}{}{}{}{}{}{}{}{}{}{}{}",
                ffmt::i(nina, 4),
                ffmt::i(ntemp, 5),
                ffmt::i(nsigz, 6),
                pe(sgref, 12, 2),
                ffmt::i(ires, 4),
                pe(sigp, 11, 3),
                ffmt::i(mti, 4),
                ffmt::i(mtc, 4),
                ffmt::i(ip1opt, 4),
                ffmt::i(norf, 5),
                ffmt::i(ipos, 5),
                ffmt::i(iposr, 4)
            ));
            if cfg!(feature = "powr-cards") {
                // tools/njoy-powr-cards.patch: keep the values where they
                // are read back below
                let vals = [
                    nina as f64,
                    ntemp as f64,
                    nsigz as f64,
                    sgref,
                    ires as f64,
                    sigp,
                    mti as f64,
                    mtc as f64,
                    ip1opt as f64,
                    norf as f64,
                    ipos as f64,
                    iposr as f64,
                ];
                for (j, x) in vals.iter().enumerate() {
                    self.scr[iscr + k + j as i32] = *x;
                }
            }
            k += 12;
        }
        if self.if5 != 0 {
            // file 5 input
            let iz = newmat + 1;
            let v = c.io.input.read(&[Kind::Int, Kind::Int]);
            self.ntis = geti(&v[0], self.ntis);
            self.nfis = geti(&v[1], self.nfis);
            let ntis = self.ntis;
            let nfis = self.nfis;
            let lspace = immax - newmat;
            let nfis1 = nfis + 1;
            let nfis2 = nfis + 2;
            let mut scr2 = RV::new(0);
            if lspace < nfis2 {
                scr2 = RV::new(nfis2.max(0) as usize + SLACK);
            }
            alloc_iv(&mut self.al_idat, &mut self.idat, ntis, "idat");
            alloc_iv(&mut self.al_idbt, &mut self.idbt, ntis, "idbt");
            let nw = ntis * nfis1;
            alloc_rv(&mut self.al_yld, &mut self.yld, nw, "yld");
            for i in 1..=nw {
                self.yld[i] = 0.0;
            }
            let kinds = vec![Kind::Int; nfis.max(0) as usize];
            let v = c.io.input.read(&kinds);
            for i in 1..=nfis {
                self.idbt[i] = geti(&v[(i - 1) as usize], self.idbt[i]);
            }
            for i in 1..=nfis {
                if lspace >= nfis2 {
                    im[iz - 1 + i] = self.idbt[i];
                }
                if lspace < nfis2 {
                    km[i] = self.idbt[i];
                }
            }
            let ids: Vec<i32> = if lspace >= nfis2 {
                (1..=nfis).map(|i| im[iz - 1 + i]).collect()
            } else {
                (1..=nfis).map(|i| km[i]).collect()
            };
            let mut s = format!(
                " no. time-dependent nuclides .......... {}\n \
                 no. fissionable burnup nuclides ...... {}\n \
                 fissionable burnup nuclide ids ....... \n",
                ffmt::i(ntis, 10),
                ffmt::i(nfis, 10)
            );
            if ids.is_empty() {
                s.push_str("\n");
            }
            for ch in ids.chunks(9) {
                s.push('\n');
                s.push_str(&" ".repeat(16));
                for &x in ch {
                    s.push_str(&ffmt::i(x, 12));
                }
            }
            c.io.so(&s);
            c.io.so("\n  ident    decay     yields\n ------  ----------  ------");
            for i in 1..=ntis {
                for _j in 1..=nfis2 {
                    if lspace >= nfis2 {
                        z[iz - 1 + _j] = 0.0;
                    }
                    if lspace < nfis2 {
                        scr2[i] = 0.0;
                    }
                }
                let nz = nfis2;
                let kinds = vec![Kind::Real; nz.max(0) as usize];
                let v = c.io.input.read(&kinds);
                for j in 1..=nz {
                    if lspace >= nfis2 {
                        z[iz + j - 1] = getr(&v[(j - 1) as usize], z[iz + j - 1]);
                    } else {
                        scr2[j] = getr(&v[(j - 1) as usize], scr2[j]);
                    }
                }
                let identa = if lspace >= nfis2 { nint(z[iz]) } else { nint(scr2[1]) };
                self.idat[i] = identa;
                let index = 1 + (i - 1) * nfis1;
                let mut izero = 0;
                for j in 2..=nz {
                    if lspace >= nfis2 {
                        self.yld[index + j - 2] = z[iz + j - 1];
                    }
                    if lspace < nfis2 {
                        self.yld[index + j - 2] = scr2[j];
                    }
                    if j > 2 && self.yld[index + j - 2] != ZERO {
                        izero = 1;
                    }
                }
                let vals: Vec<f64> = if izero > 0 {
                    (1..=nfis1).map(|j| self.yld[index + j - 1]).collect()
                } else {
                    vec![self.yld[index]]
                };
                // '(i7,1p,10e12.4/(19x,9e12.4))'
                let mut s = ffmt::i(identa, 7);
                let first = vals.len().min(10);
                for &x in &vals[..first] {
                    s.push_str(&pe(x, 12, 4));
                }
                for ch in vals[first..].chunks(9) {
                    s.push('\n');
                    s.push_str(&" ".repeat(19));
                    for &x in ch {
                        s.push_str(&pe(x, 12, 4));
                    }
                }
                c.io.so(&s);
            }
            if self.if5 != 1 {
                let nw = 3 * ntis;
                let mut f1p = RV::new(nw.max(0) as usize + SLACK);
                c.io.so("   ident     aw     indfis  ntemp\n   -----  --------  ------  -----");
                for i in 1..=ntis {
                    let index = 1 + (i - 1) * 3;
                    let v = c.io.input.read(&[Kind::Real; 3]);
                    for j in 1..=3 {
                        f1p[index + j - 1] = getr(&v[(j - 1) as usize], f1p[index + j - 1]);
                    }
                    let identa = self.idat[i];
                    // '(3x,i6,f10.4,i6,i7)' with a real for the i6 item
                    let _ = (identa, f1p[index]);
                    c.io.flush_all();
                    fatal_runtime(
                        "Fortran runtime error: Expected INTEGER for item 3 in formatted transfer, got REAL",
                    );
                }
            }
        }
        c.io.so("\n");

        //--allocate storage.
        alloc_iv(&mut self.al_mat, &mut self.mat, newmat, "mat");
        alloc_iv(&mut self.al_idnt, &mut self.idnt, newmat, "idnt");
        alloc_iv(&mut self.al_ninat, &mut self.ninat, newmat, "ninat");
        alloc_iv(&mut self.al_irest, &mut self.irest, newmat, "irest");
        alloc_rv(&mut self.al_s0rf, &mut self.s0rf, newmat, "s0rf");
        alloc_iv(&mut self.al_ip1tab, &mut self.ip1tab, newmat, "ip1tab");
        alloc_iv(&mut self.al_norft, &mut self.norft, newmat, "norft");
        alloc_iv(&mut self.al_mtit, &mut self.mtit, newmat, "mtit");
        alloc_iv(&mut self.al_mtct, &mut self.mtct, newmat, "mtct");
        alloc_rv(&mut self.al_awrt, &mut self.awrt, newmat, "awrt");
        alloc_iv(&mut self.al_indf, &mut self.indf, newmat, "indf");
        alloc_iv(&mut self.al_ntmp, &mut self.ntmp, newmat, "ntmp");
        alloc_iv(&mut self.al_nszt, &mut self.nszt, newmat, "nszt");
        let mut nwflx = self.nrg * nt * self.nzmax;
        if nwflx < self.ngnd * 2 {
            nwflx = self.ngnd * 2;
        }
        alloc_rv(&mut self.al_flux, &mut self.flux, nwflx, "flux");
        alloc_rv(&mut self.al_sigpt, &mut self.sigpt, newmat, "sigpt");
        let mut nn = self.nnina;
        if nn == 0 {
            nn = self.nres;
        }
        let nwglam = nn * self.nrg;
        if nwglam > 0 {
            alloc_rv(&mut self.al_glamt, &mut self.glamt, nwglam, "glamt");
        }
        if self.mode == 0 {
            alloc_iv(&mut self.al_ipost, &mut self.ipost, newmat, "ipost");
            alloc_iv(&mut self.al_iposrt, &mut self.iposrt, newmat, "iposrt");
        }
        if self.nnina != 0 || self.nres != 0 {
            let jscr = iscr + k + 1;
            for i in 1..=nn {
                for j in 1..=self.nrg {
                    self.glamt[i + nn * (j - 1)] = 0.0;
                }
                let kinds = vec![Kind::Real; self.nrg as usize];
                let v = c.io.input.read(&kinds);
                for kk in 1..=self.nrg {
                    self.scr[jscr + kk - 1] = getr(&v[(kk - 1) as usize], self.scr[jscr + kk - 1]);
                }
                if cfg!(feature = "powr-cards") {
                    for j in 1..=self.nrg {
                        self.glamt[i + nn * (j - 1)] = self.scr[jscr + j - 1];
                    }
                }
            }
        }
        let mut k = 0;
        for i in 1..=newmat {
            self.mat[i] = im[i];
            self.ninat[i] = nint(self.scr[iscr + k]);
            self.ntmp[i] = nint(self.scr[iscr + k + 1]);
            self.nszt[i] = nint(self.scr[iscr + k + 2]);
            self.s0rf[i] = self.scr[iscr + k + 3];
            self.irest[i] = nint(self.scr[iscr + k + 4]);
            self.sigpt[i] = self.scr[iscr + k + 5];
            self.mtit[i] = nint(self.scr[iscr + k + 6]);
            self.mtct[i] = nint(self.scr[iscr + k + 7]);
            self.ip1tab[i] = nint(self.scr[iscr + k + 8]);
            self.norft[i] = nint(self.scr[iscr + k + 9]);
            k += 12;
        }
        if self.mode <= 0 {
            let mut k = 0;
            for i in 1..=newmat {
                self.ipost[i] = nint(self.scr[iscr + k + 10]);
                self.iposrt[i] = nint(self.scr[iscr + k + 11]);
                k += 12;
            }
        }
        self.nscr = 10;
        self.kscr = 11;
        self.mscr = 12;
        c.repoz(self.nout);
        self.pinit(c, iopt);

        //--process data for resonance calculation.
        self.sfile2(c);

        //--process effective resonance integrals.
        self.sfile3(c, &mut ntemp, nt);

        //--process cross sections.
        self.sfile4(c, &mut ntemp, nt);

        //--process burnup data.
        if self.if5 > 0 {
            self.sf5out(c);
        }

        //--process p1 scattering matrix.
        self.sfile6(c, &mut ntemp);

        //--write output tape.
        if self.nout != 0 {
            self.cpmout(c, nlib, idatt);
        }

        //--cpm is finished.
        c.repoz(self.nout);
    }

    /// Store materials from ngendf.
    fn pinit(&mut self, c: &mut Ctx, iopt: i32) {
        let ng = self.ngendf;
        let ngnd = self.ngnd;
        let (mut nb, mut nw) = (0i32, 0i32);
        let ngnd1 = ngnd + 1;
        alloc_rv(&mut self.al_egb, &mut self.egb, ngnd1, "egb");
        self.iu235 = 92235;
        self.ipu239 = 94239;
        self.iu = 0;
        self.nfiss[1] = 0;
        self.ipu = 0;
        self.nfiss[2] = 0;
        if iopt == 1 {
            c.io.so("\n materials to be processed ............ ");
        }
        c.repoz(ng);
        c.tpidio(ng, 0, 0, &mut self.scr.at(1), &mut nb, &mut nw);
        let mut matl = 0;
        let mut nmat = 0;
        let mut matd = 0i32;
        let (mut icont, mut mat1, mut irepoz) = (0, 0, 0);
        let (mut i318, mut mess3, mut i618, mut mess6) = (0, 0, 0, 0);
        let mut lbl = 100;
        if self.nnina <= 0 {
            let mut n2 = 0;
            for i in 1..=self.newmat {
                if self.ninat[i] == 2 {
                    n2 += 1;
                }
            }
            if n2 == self.newmat {
                lbl = 110;
            }
        }
        loop {
            match lbl {
                100 => {
                    nmat += 1;
                    if nmat > self.newmat {
                        lbl = 180;
                        continue;
                    }
                    icont = 0;
                    mat1 = 0;
                    irepoz = 0;
                    self.idnt[nmat] = 0;
                    if self.mat[nmat] < 0 {
                        self.idnt[nmat] = 100;
                    }
                    if self.idnt[nmat] > 0 {
                        self.mat[nmat] = -self.mat[nmat];
                    }
                    if iopt == 0 {
                        matd = self.mat[nmat];
                    }
                    let nina = self.ninat[nmat];
                    if nina < 2 {
                        lbl = 110;
                        continue;
                    }
                    self.idnt[nmat] = matd + self.idnt[nmat];
                    lbl = 100;
                }
                110 => {
                    c.contio(ng, 0, 0, &mut self.scr.at(1), &mut nb, &mut nw);
                    if c.math == 0 {
                        lbl = 110;
                        continue;
                    }
                    if c.math < 0 {
                        lbl = 170;
                        continue;
                    }
                    icont += 1;
                    if iopt == 0 && icont == 1 {
                        mat1 = c.math;
                    }
                    if c.math == matl {
                        lbl = 160;
                        continue;
                    }
                    matl = c.math;
                    if iopt == 0 && irepoz > 0 && c.math == mat1 {
                        lbl = 170;
                        continue;
                    }
                    if iopt == 0 && c.math != matd {
                        lbl = 160;
                        continue;
                    }
                    if iopt == 1 {
                        self.mat[nmat] = c.math;
                    }
                    if iopt == 1 {
                        c.io.so(&format!("        {}", ffmt::i(c.math, 5)));
                    }
                    let iza = nint(c.c1h);
                    self.idnt[nmat] = iza + self.idnt[nmat];
                    self.awrt[nmat] = c.c2h;
                    if iza == self.iu235 {
                        self.iu = 1;
                    }
                    if iza == self.ipu239 {
                        self.ipu = 1;
                    }
                    if nmat <= 1 {
                        let nz = c.l2h;
                        let ntw = c.n2h;
                        c.listio(ng, 0, 0, &mut self.scr.at(1), &mut nb, &mut nw);
                        let i1 = 7 + ntw + nz;
                        self.ngn = c.l1h;
                        if self.ngn == ngnd {
                            for i in 1..=ngnd1 {
                                self.egb[i] = self.scr[ngnd1 - i + i1];
                            }
                        } else {
                            // calculate group bounds for collapsed structure
                            if self.ngn != self.ngnf {
                                c.error("init", "wrong group structure.", " ");
                            }
                            let ngn1 = self.ngn + 1;
                            let lastg = 0;
                            for i in 1..=ngn1 {
                                let ig = self.icgrp(i);
                                if ig != lastg && ig != 0 {
                                    let jg = ngnd - ig + 1;
                                    self.egb[1 + jg] = self.scr[i1 + i - 1];
                                }
                            }
                        }
                    }
                    // check for presence of mf3, mt18 and mf6, mt18.
                    i318 = 0;
                    mess3 = 0;
                    i618 = 0;
                    mess6 = 0;
                    lbl = 125;
                }
                125 => {
                    c.tosend(ng, 0, 0, &mut self.scr.at(1));
                    lbl = 130;
                }
                130 => {
                    c.contio(ng, 0, 0, &mut self.scr.at(1), &mut nb, &mut nw);
                    if c.math == 0 {
                        lbl = 150;
                        continue;
                    }
                    if c.math == -1 {
                        lbl = 170;
                        continue;
                    }
                    if c.math != matd {
                        lbl = 145;
                        continue;
                    }
                    if c.mfh == 0 || c.mth == 0 {
                        lbl = 130;
                        continue;
                    }
                    if c.mfh == 3 && c.mth >= 18 && c.mth <= 38 {
                        // 135
                        if c.mth > 21 && c.mth < 38 {
                            lbl = 125;
                            continue;
                        }
                        if c.mth == 18 {
                            i318 = 1;
                        }
                        if c.mth == 18 || i318 == 0 || mess3 > 0 {
                            lbl = 125;
                            continue;
                        }
                        let strng = format!(
                            "mat={} mf={} has both mt18 and mt{}",
                            ffmt::i(c.math, 4),
                            ffmt::i(c.mfh, 2),
                            ffmt::i(c.mth, 2)
                        );
                        c.mess("pinit", &strng, "mt18 used in all subroutines");
                        mess3 = 1;
                        lbl = if mess3 > 0 && mess6 > 0 { 150 } else { 125 };
                        continue;
                    }
                    if c.mfh == 6 && c.mth >= 18 && c.mth <= 38 {
                        // 140
                        if c.mth > 21 && c.mth < 38 {
                            lbl = 125;
                            continue;
                        }
                        if c.mth == 18 {
                            i618 = 1;
                        }
                        if c.mth == 18 || i618 == 0 || mess6 > 0 {
                            lbl = 125;
                            continue;
                        }
                        let strng = format!(
                            "mat={} mf={} has both mt18 and mt{}",
                            ffmt::i(c.math, 4),
                            ffmt::i(c.mfh, 2),
                            ffmt::i(c.mth, 2)
                        );
                        c.mess("pinit", &strng, "mt18 used in all subroutines");
                        mess6 = 1;
                        lbl = if mess3 > 0 && mess6 > 0 { 150 } else { 125 };
                        continue;
                    }
                    lbl = 125;
                }
                145 => {
                    c.skiprz(ng, -3);
                    c.math = matd;
                    lbl = 150;
                }
                150 => {
                    if nmat == self.newmat {
                        lbl = 180;
                        continue;
                    }
                    if c.math == 0 {
                        lbl = 100;
                        continue;
                    }
                    c.tomend(ng, 0, 0, &mut self.scr.at(1));
                    lbl = 100;
                }
                160 => {
                    c.tomend(ng, 0, 0, &mut self.scr.at(1));
                    lbl = 110;
                }
                170 => {
                    if irepoz > 0 || iopt > 0 {
                        // 175
                        if iopt == 0 {
                            let strng = format!("mat={} not found on gendf", ffmt::i(matd, 4));
                            c.error("pinit", &strng, " ");
                        }
                        nmat -= 1;
                        let strng = format!("only {} mats found on gendf", ffmt::i(nmat, 2));
                        c.mess("pinit", &strng, " ");
                        self.newmat = nmat;
                        lbl = 180;
                        continue;
                    }
                    c.repoz(ng);
                    c.tpidio(ng, 0, 0, &mut self.scr.at(1), &mut nb, &mut nw);
                    irepoz = 1;
                    matl = 0;
                    lbl = 110;
                }
                _ => {
                    // 180
                    let eg = &self.egb.0[..ngnd1 as usize];
                    if self.ngn == ngnd {
                        c.io.so(&format!(
                            "\n group structure\n ---------------\n{}",
                            rowsj(eg, 6, |x| pe(x, 12, 4))
                        ));
                    }
                    if self.ngn != ngnd {
                        c.io.so(&format!(
                            "\n group structure (collapsed from {} group structure)\n \
                             ----------------------------------------------------\n{}",
                            ffmt::i(self.ngn, 3),
                            rowsj(eg, 6, |x| pe(x, 12, 4))
                        ));
                    }
                    // find lethargy factors for resonance groups
                    let il = self.nfg + 1;
                    let ih = self.nfg + self.nrg;
                    for i in il..=ih {
                        let aleth = log(self.egb[i] / self.egb[i + 1]);
                        self.rleth[i - self.nfg] = 1.0 / aleth;
                    }
                    // find lethargy factor for ngref
                    let jgref = ngnd1 - self.ngref;
                    let aleth = log(self.egb[jgref] / self.egb[1 + jgref]);
                    self.reflth = aleth;
                    if self.mode < 2 {
                        dealloc(&mut self.al_egb, "egb");
                    }
                    let _ = (i318, i618);
                    return;
                }
            }
        }
    }

    /// Process data for resonance calculation.
    fn sfile2(&mut self, c: &mut Ctx) {
        let eps = 1.0e-5f64;
        let dilinf = 1.0e10f64;
        let small = 1.0e-9f64;
        let lz = 6;
        let ng = self.ngendf;
        let ngnd = self.ngnd;
        let nscr = self.nscr;
        let (mut nb, mut nw) = (0i32, 0i32);

        //--allocate storage.
        c.openz(-nscr, 1);
        if self.nnina != 0 {
            c.io.so("\n\n");
            c.io.so("\n ***subfile 2***resonance data");
            let nrg = self.nrg;
            let nwloc = nrg * 3;
            let mut a = RV::new(nwloc as usize + SLACK);
            let iloc = 1;
            let loctot = iloc;
            let locnus = loctot + nrg;
            let locsf0 = locnus + nrg;
            let nghi = ngnd - self.nfg;
            let nglo = nghi - nrg + 1;

            //--search for desired materials.
            c.repoz(ng);
            c.tpidio(ng, 0, 0, &mut self.scr.at(1), &mut nb, &mut nw);
            let mut nmat = 0;
            let mut jnina = 0;
            let (mut matd, mut isg, mut mf3mt2, mut i318, mut i618) = (0, 0, 0, 0, 0);
            let mut sgref = 0.0f64;
            let (mut nl, mut nz, mut ntw, mut l, mut ig, mut jg, mut isgz) = (0, 0, 0, 1, 0, 0, 0);
            let mut ng2 = 0;
            let mut lbl = 100;
            loop {
                match lbl {
                    100 => {
                        nmat += 1;
                        matd = self.mat[nmat];
                        let nina = self.ninat[nmat];
                        if nina == 3 {
                            lbl = 315;
                            continue;
                        }
                        if nina > 0 {
                            lbl = 320;
                            continue;
                        }
                        jnina += 1;
                        sgref = self.s0rf[nmat];
                        isg = 0;
                        if (sgref - dilinf).abs() >= sgref / 100.0 {
                            isg = 1;
                        }
                        mf3mt2 = 0;
                        i318 = 0;
                        i618 = 0;
                        for i in 1..=nwloc {
                            a[i - 1 + iloc] = 0.0;
                        }
                        for i in 1..=ngnd {
                            self.cflux[i] = 0.0;
                            self.cflux[i + ngnd] = 0.0;
                        }
                        lbl = 110;
                    }
                    110 => {
                        c.contio(ng, 0, 0, &mut self.scr.at(1), &mut nb, &mut nw);
                        if c.math == matd {
                            lbl = 130;
                            continue;
                        }
                        if c.math == 0 {
                            lbl = 110;
                            continue;
                        }
                        if c.math > 0 {
                            c.tomend(ng, 0, 0, &mut self.scr.at(1));
                            lbl = 110;
                            continue;
                        }
                        c.repoz(ng);
                        c.tpidio(ng, 0, 0, &mut self.scr.at(1), &mut nb, &mut nw);
                        lbl = 110;
                    }
                    //--search for required reactions.
                    120 => {
                        c.contio(ng, 0, 0, &mut self.scr.at(1), &mut nb, &mut nw);
                        lbl = if c.math == 0 { 300 } else { 130 };
                    }
                    130 => {
                        nl = c.l1h;
                        nz = c.l2h;
                        ntw = c.n2h;
                        if c.mfh == 0 || c.mth == 0 {
                            lbl = 120;
                            continue;
                        }
                        if c.math != matd {
                            lbl = 300;
                            continue;
                        }
                        if !(c.mfh != 3 || c.mth < 2 || mf3mt2 > 0) {
                            mf3mt2 = 1;
                            for i in 1..=ngnd {
                                if self.cflux[i] != ZERO {
                                    self.cflux[i] = 1.0 / self.cflux[i];
                                }
                                if self.cflux[ngnd + i] != ZERO {
                                    self.cflux[ngnd + i] = 1.0 / self.cflux[ngnd + i];
                                }
                            }
                        }
                        lbl = 140;
                    }
                    140 => {
                        l = 1;
                        c.listio(ng, 0, 0, &mut self.scr.at(l), &mut nb, &mut nw);
                        if c.mfh == 1 {
                            // 165
                            isgz = 0;
                            for i in 1..=nz {
                                let sigz = self.scr[6 + ntw + i];
                                if (sigz - sgref).abs() <= eps * sigz {
                                    isgz = i;
                                }
                            }
                            if isgz == 0 {
                                c.error("sfile2", "reference sigma zero not found.", " ");
                            }
                            lbl = 170;
                            continue;
                        }
                        ng2 = c.l1h;
                        ig = c.n2h;
                        while nb != 0 {
                            l += nw;
                            c.moreio(ng, 0, 0, &mut self.scr.at(l), &mut nb, &mut nw);
                            if l + nw - 1 > self.nwscr {
                                c.error("sfile2", "storage exceeded.", " ");
                            }
                        }
                        if ig == 0 {
                            lbl = 140;
                            continue;
                        }
                        let igc = self.icgrp(ig);
                        if igc == 0 || igc > ngnd {
                            lbl = 250;
                            continue;
                        }
                        if igc < nglo || igc > nghi {
                            lbl = 250;
                            continue;
                        }
                        let kg = ngnd - igc + 1;
                        jg = kg - self.nfg;
                        l = 1;
                        if c.mth == 455 {
                            // 240
                            if c.mfh == 5 {
                                lbl = 170;
                                continue;
                            }
                            // delayed neutrons
                            let loc = locnus + jg - 1;
                            let loca = l + lz;
                            let locf = l + lz;
                            a[loc] = a[loc]
                                + self.scr[loca + 1] * self.scr[loca + 2] * self.scr[locf] * self.cflux[jg];
                            lbl = 250;
                            continue;
                        }
                        if c.mfh == 3 {
                            lbl = 200;
                            continue;
                        }
                        if c.mfh == 6 {
                            lbl = 220;
                            continue;
                        }
                        lbl = 170;
                    }
                    170 => {
                        c.tosend(ng, 0, 0, &mut self.scr.at(1));
                        lbl = 120;
                    }
                    //--file 3 cross sections.
                    200 => {
                        let mth = c.mth;
                        if mth == 18 {
                            i318 = 1;
                        }
                        if mth > 18 && i318 == 1 {
                            lbl = 170;
                            continue;
                        }
                        if mth != 1 && mth != 2 && (mth < 18 || mth > 21) && mth != 38 {
                            lbl = 170;
                            continue;
                        }
                        if (mth >= 18 && mth <= 21) || mth == 38 {
                            // fission
                            let loc = locsf0 + jg - 1;
                            let loca = l + lz + nl * nz;
                            let locf = l + lz;
                            a[loc] = a[loc] + self.scr[loca] * self.scr[locf] * self.cflux[jg];
                        } else if mth == 2 {
                            // elastic
                            let loc = loctot + jg - 1;
                            let mut loca = l + lz + nl * nz;
                            if isg > 0 && nz >= isgz {
                                loca = l + lz + nl * (isgz - 1 + nz);
                            }
                            let mut locf = l + lz;
                            if isg > 0 && nz >= isgz {
                                locf = l + lz + nl * (isgz - 1);
                            }
                            let mut jgc = jg;
                            if isg > 0 && nz >= isgz {
                                jgc = jg + ngnd;
                            }
                            a[loc] = a[loc] + self.scr[loca] * self.scr[locf] * self.cflux[jgc];
                        } else if mth == 1 {
                            // total: coarse group flux for infinite dilution
                            self.cflux[jg] = self.cflux[jg] + self.scr[l + lz];
                            if isg != 0 && nz >= isgz {
                                // for reference sigma zero
                                self.cflux[jg + ngnd] = self.cflux[jg + ngnd] + self.scr[l + lz + nl * (isgz - 1)];
                            }
                        }
                        lbl = 250;
                    }
                    //--file 6 matrices.
                    220 => {
                        let mth = c.mth;
                        if (mth < 18 || mth > 21) && mth != 38 {
                            lbl = 170;
                            continue;
                        }
                        if mth == 18 {
                            i618 = 1;
                        }
                        if mth > 18 && i618 == 1 {
                            lbl = 170;
                            continue;
                        }
                        // fission
                        let loc = locnus + jg - 1;
                        let locf = l + lz;
                        let loca = l + lz - nl * nz;
                        for i in 2..=ng2 {
                            a[loc] = self.scr[nl * nz * i + loca] * self.scr[locf] * self.cflux[jg] + a[loc];
                        }
                        lbl = 250;
                    }
                    250 => {
                        lbl = if ig < self.ngn { 140 } else { 120 };
                    }
                    //--desired data is loaded.
                    //--find nu from nu times sigma f.
                    300 => {
                        for ig in 1..=nrg {
                            let loc = locsf0 + ig - 1;
                            let loca = locnus + ig - 1;
                            if a[loc] != ZERO {
                                a[loca] = a[loca] / a[loc];
                            }
                            if a[loca] < small {
                                a[loca] = 0.0;
                            }
                        }
                        lbl = 319;
                    }
                    //--read in nus and tot for nuclides having nina=3
                    315 => {
                        let kinds = vec![Kind::Real; nrg as usize];
                        let v = c.io.input.read(&kinds);
                        for j in 1..=nrg {
                            a[locnus + j - 1] = getr(&v[(j - 1) as usize], a[locnus + j - 1]);
                        }
                        let v = c.io.input.read(&kinds);
                        for j in 1..=nrg {
                            a[loctot + j - 1] = getr(&v[(j - 1) as usize], a[loctot + j - 1]);
                        }
                        lbl = 319;
                    }
                    319 => {
                        self.sf2out(c, nmat, jnina, loctot, locnus, &a);
                        lbl = 320;
                    }
                    _ => {
                        // 320
                        if nmat < self.newmat {
                            lbl = 100;
                            continue;
                        }
                        break;
                    }
                }
            }

            //--sfile2 is finished.
            if self.iprint == 1 {
                c.io.so("\n");
            }
        }
        let mut r = RecW::new();
        r.r(EFLAG);
        c.io.write_rec(nscr, r);
    }

    /// Write subfile 2.
    fn sf2out(&mut self, c: &mut Ctx, nmat: i32, jnina: i32, loctot: i32, locnus: i32, a: &RV) {
        let ident = self.idnt[nmat];
        if self.iprint == 1 {
            c.io.so(&format!("\n nuclide {}\n -------", ffmt::i(ident, 6)));
        }
        if self.iprint == 1 {
            c.io.so(&format!(
                "\n  group{}lambda{}sigp{}nu{}sigsctot\n {}-",
                " ".repeat(7),
                " ".repeat(8),
                " ".repeat(10),
                " ".repeat(9),
                "----------".repeat(6)
            ));
        }
        let sigp = self.sigpt[nmat];
        for i in 1..=self.nrg {
            let ig = i + self.nfg;
            let glam = rd0(&self.glamt, jnina + self.nnina * (i - 1));
            let mut r = RecW::new();
            r.r(glam).r(sigp).r(a[locnus - 1 + i]).r(a[loctot - 1 + i]);
            c.io.write_rec(self.nscr, r);
            if self.iprint == 1 {
                c.io.so(&format!(
                    "{}     {}{}{}{}",
                    ffmt::i(ig, 5),
                    pe(glam, 13, 4),
                    pe(sigp, 13, 4),
                    pe(a[locnus - 1 + i], 13, 4),
                    pe(a[loctot - 1 + i], 13, 4)
                ));
            }
        }
    }

    /// Process effective resonance integrals.
    fn sfile3(&mut self, c: &mut Ctx, jtemp: &mut i32, nt: i32) {
        let lz = 6;
        let ng = self.ngendf;
        let ngnd = self.ngnd;
        let (mut nb, mut nw) = (0i32, 0i32);

        //--allocate storage.
        if self.nres != 0 {
            c.io.so("\n\n");
            c.io.so("\n ***subfile 3***resonance integrals");
            let nzmax = self.nzmax;
            let nrg = self.nrg;
            nw = nt + 1;
            let mut tempt = RV::new(nw as usize + SLACK);
            let nwloc = nzmax * nt * nrg * 3;
            let mut a = RV::new(nwloc as usize + SLACK);
            let iloc = 1;
            let nwflxr = nzmax * nt;
            let mut flxr = RV::new(nwflxr as usize + SLACK);
            let mut sigz = RV::new(nzmax as usize + SLACK);
            let mut snu = RV::new(nrg as usize + SLACK);
            let mut absv = RV::new(nt as usize + SLACK);
            let mut nwfa = 0;
            if nzmax * nt > nrg + 1 {
                nwfa = nzmax * nt;
            }
            let mut fa = RV::new((nwfa.max(70)) as usize + SLACK);
            let locabs = iloc;
            let locnus = locabs + nzmax * nt * nrg;
            let locsf = locnus + nzmax * nt * nrg;
            let locs = [0, locabs, locnus, locsf];
            let mut nn = self.nnina;
            if nn == 0 {
                nn = self.nres;
            }
            let _ = nn;
            let nghi = ngnd - self.nfg;
            let nglo = nghi - nrg + 1;

            //--search for desired materials.
            c.repoz(ng);
            c.tpidio(ng, 0, 0, &mut self.scr.at(1), &mut nb, &mut nw);
            let il = 1;
            let mut nmat = 0;
            let mut irmat = 0;
            let mut jnina = 0;
            let (mut matd, mut nina, mut inorf, mut jres, mut ntnmat, mut msigz) = (0, 0, 0, 0, 0, 0);
            let (mut jfiss, mut nl, mut nz, mut ntw, mut nsigz, mut i318, mut i618) = (0, 0, 0, 0, 0, 0, 0);
            let (mut l, mut ig, mut igc, mut jg, mut lim) = (1, 0, 0, 0, 0);
            let mut xid = 0.0f64;
            let mut lbl = 110;
            loop {
                match lbl {
                    110 => {
                        nmat += 1;
                        matd = self.mat[nmat];
                        nina = self.ninat[nmat];
                        if nina == 0 {
                            jnina += 1;
                        }
                        self.indf[nmat] = 0;
                        inorf = self.norft[nmat];
                        jres = self.irest[nmat];
                        if jres <= 0 {
                            lbl = 470;
                            continue;
                        }
                        if self.nnina == 0 {
                            jnina += 1;
                        }
                        ntnmat = self.ntmp[nmat];
                        msigz = self.nszt[nmat];
                        irmat += 1;
                        for i in 1..=nwflxr {
                            flxr[i] = 0.0;
                        }
                        let nwflx = nzmax * nt * nrg;
                        for i in 1..=nwflx {
                            self.flux[i] = 0.0;
                        }
                        for i in 1..=nwloc {
                            a[i - 1 + iloc] = 0.0;
                        }
                        for i in 1..=nrg {
                            snu[i] = 0.0;
                        }
                        *jtemp = 0;
                        jfiss = 1;
                        lbl = 120;
                    }
                    120 => {
                        c.contio(ng, 0, 0, &mut self.scr.at(1), &mut nb, &mut nw);
                        if c.math == matd {
                            lbl = 130;
                            continue;
                        }
                        if c.math == 0 {
                            lbl = 120;
                            continue;
                        }
                        if c.math > 0 {
                            c.tomend(ng, 0, 0, &mut self.scr.at(1));
                            lbl = 120;
                            continue;
                        }
                        c.repoz(ng);
                        c.tpidio(ng, 0, 0, &mut self.scr.at(1), &mut nb, &mut nw);
                        lbl = 120;
                    }
                    //--search for required reactions.
                    130 => {
                        xid = self.idnt[nmat] as f64;
                        lbl = 135;
                    }
                    135 => {
                        nl = c.l1h;
                        nz = c.l2h;
                        ntw = c.n2h;
                        lbl = 140;
                    }
                    140 => {
                        l = 1;
                        c.listio(ng, 0, 0, &mut self.scr.at(l), &mut nb, &mut nw);
                        if c.mfh != 1 {
                            lbl = 160;
                            continue;
                        }
                        let temp = c.c1h;
                        *jtemp += 1;
                        tempt[*jtemp] = temp;
                        nsigz = nz;
                        if nz > msigz && msigz != 0 {
                            nsigz = msigz;
                        }
                        for i in 1..=nsigz {
                            let index = 1 + nsigz - i;
                            sigz[index] = self.scr[6 + ntw + i];
                        }
                        absv[*jtemp] = 0.0;
                        i318 = 0;
                        i618 = 0;
                        lbl = 300;
                    }
                    160 => {
                        let (mfh, mth) = (c.mfh, c.mth);
                        if mth >= 18 && mth <= 21 && inorf > 0 {
                            lbl = 300;
                            continue;
                        }
                        if mth == 38 && inorf > 0 {
                            lbl = 300;
                            continue;
                        }
                        if !(mfh == 6 && *jtemp == 1) {
                            if mfh != 3 {
                                lbl = 300;
                                continue;
                            }
                            if !(mth == 1 || (mth == 455 && *jtemp == 1)) {
                                if mth < 18 || mth > 150 {
                                    lbl = 300;
                                    continue;
                                }
                                if mth > 21 && mth < 102 && mth != 38 {
                                    lbl = 300;
                                    continue;
                                }
                                if mth == 18 {
                                    i318 = 1;
                                }
                                if mth > 18 && mth < 39 && i318 > 0 {
                                    lbl = 300;
                                    continue;
                                }
                            }
                        }
                        // 162
                        ig = c.n2h;
                        while nb != 0 {
                            l += nw;
                            c.moreio(ng, 0, 0, &mut self.scr.at(l), &mut nb, &mut nw);
                            if l + nw - 1 > self.nwscr {
                                c.error("sfile3", "storage exceeded.", " ");
                            }
                        }
                        if ig == 0 {
                            lbl = 140;
                            continue;
                        }
                        igc = self.icgrp(ig);
                        if !(c.mth == 1 && igc == self.ngref) {
                            if igc < nglo || igc > nghi {
                                lbl = 185;
                                continue;
                            }
                        }
                        // 172
                        let kg = ngnd - igc + 1;
                        jg = kg - self.nfg;
                        lim = nsigz;
                        if nsigz > nz {
                            lim = nz;
                        }
                        if c.mfh == 6 {
                            lbl = 230;
                            continue;
                        }
                        if c.mth == 1 {
                            lbl = 190;
                            continue;
                        }
                        if c.mth == 455 {
                            lbl = 225;
                            continue;
                        }
                        if c.mth == 18 {
                            lbl = 175;
                            continue;
                        }
                        if c.mth < 102 && i318 > 0 {
                            lbl = 215;
                            continue;
                        }
                        lbl = 175;
                    }
                    // absorption
                    175 => {
                        let iadd = nsigz + nzmax * ((*jtemp - 1) + nt * (jg - 1));
                        let loca = l + lz + nl * (nz - 1);
                        let locb = l + lz - nl;
                        for jz in 1..=lim {
                            a[iadd - jz + locabs] = a[iadd - jz + locabs]
                                + self.scr[nl * jz + loca] * self.scr[nl * jz + locb] / self.flux[iadd - jz + 1];
                        }
                        if c.mth < 102 {
                            lbl = 215;
                            continue;
                        }
                        absv[*jtemp] = 1.0;
                        lbl = 185;
                    }
                    185 => {
                        lbl = if ig < self.ngn { 140 } else { 310 };
                    }
                    // flux
                    190 => {
                        if igc == self.ngref && (self.ngref > nghi || self.ngref < nglo) {
                            lbl = 205;
                            continue;
                        }
                        // 195
                        let loc = 1 + nsigz + nzmax * ((*jtemp - 1) + nt * (jg - 1));
                        let loca = l + lz - nl;
                        for jz in 1..=lim {
                            self.flux[loc - jz] = self.scr[nl * jz + loca] + self.flux[loc - jz];
                        }
                        lbl = if igc == self.ngref { 205 } else { 185 };
                    }
                    // reference group flux
                    205 => {
                        for jz in 1..=lim {
                            let loca = l + lz + nl * (jz - 1);
                            let loc = 1 + nsigz - jz + nzmax * (*jtemp - 1);
                            flxr[loc] = flxr[loc] + self.scr[loca];
                        }
                        lbl = 185;
                    }
                    // nu sigma f
                    215 => {
                        let iadd = nsigz + nzmax * ((*jtemp - 1) + nt * (jg - 1));
                        let loca = l + lz + nl * (nz - 1);
                        let locb = l + lz - nl;
                        for jz in 1..=lim {
                            a[iadd - jz + locsf] = a[iadd - jz + locsf]
                                + self.scr[nl * jz + loca] * self.scr[nl * jz + locb] / self.flux[iadd - jz + 1];
                        }
                        jfiss = 2;
                        lbl = 185;
                    }
                    // delayed neutrons
                    225 => {
                        let loca = l + lz;
                        let locf = 1 + nsigz - 1 + nzmax * nt * (jg - 1);
                        snu[jg] = snu[jg]
                            + self.scr[loca + 1] * self.scr[loca + 2] * self.scr[loca] / self.flux[locf];
                        lbl = 185;
                    }
                    // save nu from file 6
                    230 => {
                        let mth = c.mth;
                        if (mth < 18 || mth > 22) && mth != 38 {
                            lbl = 300;
                            continue;
                        }
                        if mth == 18 {
                            i618 = 1;
                        }
                        if !(mth == 18 || i618 == 0) {
                            lbl = 300;
                            continue;
                        }
                        // 235
                        let ng2 = c.l1h;
                        let loca = l + lz + (il - 1) - nl * nz;
                        let locf = 1 + nsigz - 1 + nzmax * nt * (jg - 1);
                        for i in 2..=ng2 {
                            snu[jg] = snu[jg] + self.scr[nl * nz * i + loca] * self.scr[l + lz] / self.flux[locf];
                        }
                        lbl = 185;
                    }
                    300 => {
                        c.tosend(ng, 0, 0, &mut self.scr.at(1));
                        lbl = 310;
                    }
                    310 => {
                        c.contio(ng, 0, 0, &mut self.scr.at(1), &mut nb, &mut nw);
                        if c.math < 0 {
                            lbl = 350;
                            continue;
                        }
                        if c.math == 0 || c.mfh == 0 || c.mth == 0 {
                            lbl = 310;
                            continue;
                        }
                        if c.math != matd {
                            lbl = 350;
                            continue;
                        }
                        if c.mfh == 1 {
                            lbl = 360;
                            continue;
                        }
                        if c.mfh != 3 && c.mfh != 6 {
                            lbl = 300;
                            continue;
                        }
                        lbl = 135;
                    }
                    350 => {
                        lbl = if *jtemp > 1 || jfiss == 1 { 400 } else { 370 };
                    }
                    360 => {
                        lbl = if *jtemp > 1 || jfiss == 1 { 130 } else { 370 };
                    }
                    //--correct nu
                    370 => {
                        for jg in 1..=nrg {
                            let iterm = nsigz - 1 + nzmax * nt * (jg - 1);
                            let loca = locsf + iterm;
                            if a[loca] == ZERO {
                                snu[jg] = 0.0;
                            }
                            if a[loca] != ZERO {
                                snu[jg] = snu[jg] / a[loca];
                            }
                        }
                        if c.math != matd {
                            lbl = 400;
                            continue;
                        }
                        if ntnmat == 0 {
                            lbl = 130;
                            continue;
                        }
                        if *jtemp == ntnmat {
                            lbl = 400;
                            continue;
                        }
                        lbl = 130;
                    }
                    //--desired data is loaded.
                    400 => {
                        if jfiss != 1 {
                            // calculate nu times sigma f
                            for jg in 1..=nrg {
                                for jz in 1..=nsigz {
                                    let iterm = jz - 1 + nzmax * nt * (jg - 1) - nzmax;
                                    for jtem in 1..=*jtemp {
                                        a[nzmax * jtem + iterm + locnus] = a[nzmax * jtem + iterm + locsf] * snu[jg];
                                    }
                                }
                            }
                        }
                        // check flux arrays for temperature dependence.
                        for jz in 1..=nsigz {
                            for it in 1..=*jtemp {
                                if flxr[jz + nzmax * (it - 1)] == 0.0 {
                                    flxr[jz + nzmax * (it - 1)] = rd0(&flxr, jz + nzmax * (it - 2));
                                }
                            }
                        }
                        for jg in 1..=nrg {
                            for jz in 1..=nsigz {
                                for it in 1..=*jtemp {
                                    let index = jz + nzmax * ((it - 1) + nt * (jg - 1));
                                    if self.flux[index] == ZERO {
                                        let indexl = jz + nzmax * ((it - 2) + nt * (jg - 1));
                                        self.flux[index] = rd0(&self.flux, indexl);
                                    }
                                }
                            }
                        }
                        // correct non-temperature dependent absorption if necessary
                        if *jtemp != 1 {
                            for it in 2..=*jtemp {
                                if absv[it] <= ZERO {
                                    for jg in 1..=nrg {
                                        let loc = locabs + nzmax * ((it - 2) + nt * (jg - 1));
                                        let loca = loc + nzmax;
                                        for jz in 1..=nsigz {
                                            a[jz + loca] = a[jz + loca] + a[jz + loc];
                                        }
                                    }
                                }
                            }
                        }
                        if (self.nnina > 0 && nina > 0) || (self.nnina == 0 && jres == 0) {
                            // convert njoy cross sections to cpm resonance integrals
                            for jg in 1..=nrg {
                                let amda = rd0(&self.glamt, jnina + self.nnina * (jg - 1));
                                let sigp = self.sigpt[nmat];
                                for ix in 1..=2 {
                                    for it in 1..=*jtemp {
                                        let loc = locs[ix] + nzmax * (it - 1 + nt * (jg - 1));
                                        cpmize(nsigz, &sigz, &mut a.0[(loc - 1) as usize..], amda, sigp);
                                    }
                                }
                            }
                            // write out results
                            let _ = nwfa;
                            self.sf3out(c, xid, nmat, *jtemp, nsigz, nt, jfiss, &tempt, &sigz, &flxr, &locs, &mut fa, &a);
                            self.indf[nmat] = jfiss;
                        }
                        lbl = 470;
                    }
                    _ => {
                        // 470
                        if nmat == self.newmat {
                            break;
                        }
                        if irmat == 0 {
                            lbl = 110;
                            continue;
                        }
                        if c.math != matd {
                            c.skiprz(ng, -1);
                        }
                        lbl = 110;
                    }
                }
            }
        }

        //--sfile3 is finished.
        dealloc(&mut self.al_nszt, "nszt");
        dealloc(&mut self.al_norft, "norft");
        dealloc(&mut self.al_sigpt, "sigpt");
        dealloc(&mut self.al_glamt, "glamt");
        let mut r = RecW::new();
        r.r(EFLAG);
        c.io.write_rec(self.nscr, r);
    }

    /// Write subfile 3.
    fn sf3out(
        &mut self,
        c: &mut Ctx,
        xid: f64,
        nmat: i32,
        ntemp: i32,
        nsigz: i32,
        nt: i32,
        ifiss: i32,
        temp: &RV,
        sigz: &RV,
        flxr: &RV,
        locs: &[i32; 4],
        b: &mut RV,
        a: &RV,
    ) {
        let nmax = 5;
        let ncol = 8;
        let undl = "-------------";
        let nzmax = self.nzmax;
        let nscr = self.nscr;
        let ncol1 = ncol + 1;
        let mut lim = nsigz;
        if nsigz > ncol {
            lim = ncol;
        }
        let ident = self.idnt[nmat];
        let mut n = nsigz;
        if n > nmax {
            n = nmax;
        }
        let mut r = RecW::new();
        r.r(xid);
        c.io.write_rec(nscr, r);
        let mut r = RecW::new();
        r.i(ntemp).i(nsigz);
        c.io.write_rec(nscr, r);
        let mut r = RecW::new();
        for i in 1..=ntemp {
            r.r(temp[i]);
        }
        c.io.write_rec(nscr, r);
        let mut r = RecW::new();
        for i in 1..=nsigz {
            r.r(sigz[i]);
        }
        c.io.write_rec(nscr, r);
        c.io.so(&format!("\n xid {}\n ---", ffmt::f(xid, 8, 2)));
        c.io.so(&format!(
            "\n ntemp {}  nsigp {}\n -----     -----\n\n temperatures (down)\n -------------------",
            ffmt::i(ntemp, 2),
            ffmt::i(nsigz, 2)
        ));
        let tv: Vec<f64> = (1..=ntemp).map(|i| temp[i]).collect();
        c.io.so(&rows10(&tv, 6, " ", |x| pe(x, 11, 3)));
        c.io.so("\n potential cross sections (across)\n ---------------------------------");
        let sv: Vec<f64> = (1..=nsigz).map(|i| sigz[i]).collect();
        c.io.so(&rows10(&sv, 9, &" ".repeat(10), |x| pe(x, 13, 5)));
        let jgref = self.ngnd - self.ngref + 1;
        let mut m = n;
        if m < 4 {
            m = 4;
        }
        c.io.so(&format!(
            "\n group      flux per unit lethargy normalized at group {} for nuclide {}\n ----------{}",
            ffmt::i(jgref, 2),
            ffmt::i(ident, 6),
            undl.repeat(m as usize)
        ));
        let mut cc = [0.0f64; 11];
        for ig in 1..=self.nrg {
            for it in 1..=ntemp {
                let locf = nzmax * ((it - 1) + nt * (ig - 1));
                let locr = nzmax * (it - 1);
                for jz in 1..=nsigz {
                    cc[jz as usize] = self.flux[jz + locf] * self.rleth[ig] * self.reflth / flxr[jz + locr];
                }
                let jg = ig + self.nfg;
                let cv = |i0: i32, i1: i32| -> Vec<f64> { (i0..=i1).map(|i| cc[i as usize]).collect() };
                if it == 1 {
                    let mut s = ffmt::i(jg, 5);
                    if lim >= 1 {
                        s.push_str("     ");
                    }
                    for x in cv(1, lim) {
                        s.push_str(&pe(x, 13, 5));
                    }
                    c.io.so(&s);
                }
                if it == 1 && nsigz > ncol {
                    c.io.so(&rows10(&cv(ncol1, nsigz), 9, &" ".repeat(10), |x| pe(x, 13, 5)));
                }
                if it > 1 {
                    c.io.so(&rows10(&cv(1, nsigz), 9, &" ".repeat(10), |x| pe(x, 13, 5)));
                }
            }
        }
        for il in 1..=2usize {
            let loc = locs[il];
            if self.iprint == 1 && il == 1 {
                c.io.so(&format!(
                    "\n group      absorption for nuclide {}\n ----------{}",
                    ffmt::i(ident, 6),
                    undl.repeat(n as usize)
                ));
            }
            if il != 2 || ifiss >= 2 {
                if self.iprint == 1 && il == 2 {
                    c.io.so(&format!(
                        "\n group      fission for nuclide {}\n ----------{}",
                        ffmt::i(ident, 6),
                        undl.repeat(n as usize)
                    ));
                }
                for ig in 1..=self.nrg {
                    let mut nb = 0;
                    let jg = ig + self.nfg;
                    for it in 1..=ntemp {
                        let loca = loc - 1 + nzmax * ((it - 1) + nt * (ig - 1));
                        for jz in 1..=nsigz {
                            nb += 1;
                            b[nb] = a[jz + loca];
                            cc[jz as usize] = b[nb];
                        }
                        let cv = |i0: i32, i1: i32| -> Vec<f64> { (i0..=i1).map(|i| cc[i as usize]).collect() };
                        if self.iprint != 0 {
                            if il != 2 || ifiss != 1 {
                                if it == 1 {
                                    let mut s = format!("\n{}", ffmt::i(jg, 5));
                                    if lim >= 1 {
                                        s.push_str("     ");
                                    }
                                    for x in cv(1, lim) {
                                        s.push_str(&pe(x, 13, 5));
                                    }
                                    c.io.so(&s);
                                }
                                if it == 1 && nsigz > ncol {
                                    c.io.so(&rows10(&cv(ncol1, nsigz), 9, &" ".repeat(10), |x| pe(x, 13, 5)));
                                }
                                if it > 1 {
                                    c.io.so(&rows10(&cv(1, lim), 9, &" ".repeat(10), |x| pe(x, 13, 5)));
                                }
                                if it > 1 && nsigz > ncol {
                                    c.io.so(&rows10(&cv(ncol1, nsigz), 9, &" ".repeat(10), |x| pe(x, 13, 5)));
                                }
                            }
                        }
                    }
                    let mut r = RecW::new();
                    for ib in 1..=nb {
                        r.r(b[ib]);
                    }
                    c.io.write_rec(nscr, r);
                }
            }
        }
    }
}

/// Records of `'(<lead>,1p,<per>e..)'` (the lead is repeated on each
/// record by reversion; one empty-but-lead record when there are no
/// values).
pub(super) fn rows10(v: &[f64], per: usize, lead: &str, f: impl Fn(f64) -> String) -> String {
    if v.is_empty() {
        return String::new();
    }
    let mut out: Vec<String> = Vec::new();
    for ch in v.chunks(per) {
        let mut s = lead.to_string();
        for &x in ch {
            s.push_str(&f(x));
        }
        out.push(s);
    }
    out.join("\n")
}

/// Convert NJOY shielded cross sections to CPM resonance integrals using
/// CPM values for intermediate-resonance lambda.
fn cpmize(nsigz: i32, sigz: &RV, sign: &mut [f64], amda: f64, sigp: f64) {
    let mut sig = [0.0f64; 11 + 16];
    let sg = |i: i32| sigz[i];
    for iz in 1..=nsigz {
        sig[iz as usize] = sign[(iz - 1) as usize];
    }
    for iz in 1..=nsigz {
        let mut sm = sg(iz);
        let test = 1001.0f64;
        if sm <= test {
            sign[(iz - 1) as usize] = 0.0;
            let mut sn = sm - amda * sigp;
            if sn < ZERO {
                sn = 0.0;
            }
            let mut jl = iz - 1;
            if jl < 1 {
                jl = 1;
            }
            let mut jh = iz + 1;
            if jh > nsigz {
                jh = nsigz;
            }
            for jz in jl..=jh {
                let mut temp = 1.0f64;
                sm = sg(jz);
                for kz in jl..=jh {
                    let st = sg(kz);
                    if kz != jz {
                        if sn >= sg(1) {
                            temp = temp * log(sn / st) / log(sm / st);
                        } else {
                            temp = temp * (powc(sn, 2) - powc(st, 2)) / (powc(sm, 2) - powc(st, 2));
                        }
                    }
                }
                sign[(iz - 1) as usize] = sign[(iz - 1) as usize] + temp * sig[jz as usize];
            }
        }
    }
    for iz in 1..=nsigz {
        let sm = sg(iz);
        let k = (iz - 1) as usize;
        sign[k] = sm * sign[k] / (sm + sign[k]);
    }
}
