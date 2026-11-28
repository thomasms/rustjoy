//! GROUPR File 6 feed functions (getmf6, cm2lab, f6cm, f6ddx, f6dis,
//! bach, ll2lab, f6lab, gam102, f6psp).

use super::*;
use crate::fmath::{cosh, ipow, powi, sinh};
use crate::mathm::legndr;

const MXLG: usize = 65;

/// Saved locals of the File 4/5/6 routines.
pub(super) struct Saved {
    // getmf6
    pub m6_nne: i32,
    pub m6_ne: i32,
    pub m6_int: i32,
    pub m6_jlo: i32,
    pub m6_elo: f64,
    pub m6_jhi: i32,
    pub m6_ehi: f64,
    pub m6_terml: [f64; MXLG + 1],
    pub m6_pspmax: f64,
    pub m6_langn: i32,
    pub m6_lepn: i32,
    pub m6_disc102: f64,
    pub m6_idis: i32,
    pub m6_iyss: [i32; 10],
    pub m6_izss: [i32; 10],
    pub m6_jjss: [i32; 10],
    pub m6_jloss: Vec<i32>,
    pub m6_nss: i32,
    pub m6_jzap: i32,
    pub m6_lct3: i32,
    // f6cm
    pub cm_na: i32,
    pub cm_eps: f64,
    pub cm_xc: f64,
    pub cm_ndnow: i32,
    pub cm_npnow: i32,
    pub cm_ncnow: i32,
    pub cm_elmax: f64,
    pub cm_e: f64,
    pub cm_epmax: f64,
    // f6ddx
    pub dx_enow: f64,
    pub dx_efirst: f64,
    pub dx_nl: i32,
    pub dx_inow: i32,
    pub dx_lnow: i32,
    pub dx_mnow: i32,
    pub dx_ncnow: i32,
    pub dx_na: i32,
    pub dx_illdef: i32,
    pub dx_x1: f64,
    pub dx_x2: f64,
    // f6dis
    pub ds_nl: i32,
    pub ds_inow: i32,
    pub ds_lnow: i32,
    pub ds_mnow: i32,
    pub ds_ncnow: i32,
    pub ds_na: i32,
    pub ds_enow: f64,
    // f6lab
    pub lb_xend: f64,
    pub lb_ilo: i32,
    pub lb_llo: i32,
    pub lb_mlo: i32,
    pub lb_nclo: i32,
    pub lb_xlo: f64,
    pub lb_ihi: i32,
    pub lb_lhi: i32,
    pub lb_mhi: i32,
    pub lb_nchi: i32,
    pub lb_xhi: f64,
    pub lb_elo: f64,
    pub lb_ehi: f64,
    pub lb_intt: i32,
    // getdis
    pub di_iecl: i32,
    pub di_iech: i32,
    pub di_ecl: f64,
    pub di_ech: f64,
    pub di_ecn: f64,
    pub di_nqp0: i32,
    // getfle
    pub fl_iso: i32,
    pub fl_lidp: i32,
    pub fl_ltt: i32,
    pub fl_iraw: i32,
    pub fl_ir: i32,
    pub fl_nne: i32,
    pub fl_ne: i32,
    pub fl_int: i32,
    pub fl_elo: f64,
    pub fl_ehi: f64,
    pub fl_nlo: i32,
    pub fl_nhi: i32,
    pub fl_flo: Vec<f64>,
    pub fl_fhi: Vec<f64>,
    pub fl_ltt3: i32,
    pub fl_lttn: i32,
    pub fl_fls: RV,
    // getaed
    pub ae_nwt: i32,
    pub ae_ir: i32,
    pub ae_ncyc: i32,
    pub ae_nu: i32,
    pub ae_ne: i32,
    pub ae_ie: i32,
    pub ae_nw1: i32,
    pub ae_nw2: i32,
    pub ae_l2: i32,
    pub ae_l3: i32,
    pub ae_llo: i32,
    pub ae_nlo: i32,
    pub ae_elo: f64,
    pub ae_lhi: i32,
    pub ae_nhi: i32,
    pub ae_ehi: f64,
    pub ae_itt: i32,
    pub ae_law: i32,
    pub ae_jbrag: i32,
    pub ae_nbrag: i32,
    // getgfl
    pub gf_li: i32,
    pub gf_loca: Vec<i32>,
    pub gf_ltt: i32,
    // getgyl
    pub gy_awr: f64,
    pub gy_loca: Vec<i32>,
    pub gy_ird: i32,
    pub gy_ipd: i32,
    // getsed
    pub sd_nktot: i32,
    pub sd_nupm: i32,
    pub sd_loc: [i32; 21],
    // anased
    pub an_new: i32,
    pub an_theta: f64,
    pub an_xc: f64,
    pub an_rc: f64,
    pub an_bot: f64,
    pub an_ca: f64,
    pub an_loct: i32,
    // f6psp
    pub ps_cn: f64,
    pub ps_eimax: f64,
    pub ps_ex: f64,
}

impl Default for Saved {
    fn default() -> Self {
        Saved {
            m6_nne: 0,
            m6_ne: 0,
            m6_int: 0,
            m6_jlo: 0,
            m6_elo: 0.0,
            m6_jhi: 0,
            m6_ehi: 0.0,
            m6_terml: [0.0; MXLG + 1],
            m6_pspmax: 0.0,
            m6_langn: 0,
            m6_lepn: 0,
            m6_disc102: 0.0,
            m6_idis: 0,
            m6_iyss: [0; 10],
            m6_izss: [0; 10],
            m6_jjss: [0; 10],
            m6_jloss: vec![0; 502],
            m6_nss: 0,
            m6_jzap: 0,
            m6_lct3: 0,
            cm_na: 0,
            cm_eps: 0.0,
            cm_xc: 0.0,
            cm_ndnow: 0,
            cm_npnow: 0,
            cm_ncnow: 0,
            cm_elmax: 0.0,
            cm_e: 0.0,
            cm_epmax: 0.0,
            dx_enow: 0.0,
            dx_efirst: 0.0,
            dx_nl: 0,
            dx_inow: 0,
            dx_lnow: 0,
            dx_mnow: 0,
            dx_ncnow: 0,
            dx_na: 0,
            dx_illdef: 0,
            dx_x1: 0.0,
            dx_x2: 0.0,
            ds_nl: 0,
            ds_inow: 0,
            ds_lnow: 0,
            ds_mnow: 0,
            ds_ncnow: 0,
            ds_na: 0,
            ds_enow: 0.0,
            lb_xend: 0.0,
            lb_ilo: 0,
            lb_llo: 0,
            lb_mlo: 0,
            lb_nclo: 0,
            lb_xlo: 0.0,
            lb_ihi: 0,
            lb_lhi: 0,
            lb_mhi: 0,
            lb_nchi: 0,
            lb_xhi: 0.0,
            lb_elo: 0.0,
            lb_ehi: 0.0,
            lb_intt: 0,
            di_iecl: 0,
            di_iech: 0,
            di_ecl: 0.0,
            di_ech: 0.0,
            di_ecn: 0.0,
            di_nqp0: 0,
            fl_iso: 0,
            fl_lidp: 0,
            fl_ltt: 0,
            fl_iraw: 0,
            fl_ir: 0,
            fl_nne: 0,
            fl_ne: 0,
            fl_int: 0,
            fl_elo: 0.0,
            fl_ehi: 0.0,
            fl_nlo: 0,
            fl_nhi: 0,
            fl_flo: vec![0.0; MXLG + 1],
            fl_fhi: vec![0.0; MXLG + 1],
            fl_ltt3: 0,
            fl_lttn: 0,
            fl_fls: RV::new(1000 + SLACK),
            ae_nwt: 0,
            ae_ir: 0,
            ae_ncyc: 0,
            ae_nu: 0,
            ae_ne: 0,
            ae_ie: 0,
            ae_nw1: 0,
            ae_nw2: 0,
            ae_l2: 0,
            ae_l3: 0,
            ae_llo: 0,
            ae_nlo: 0,
            ae_elo: 0.0,
            ae_lhi: 0,
            ae_nhi: 0,
            ae_ehi: 0.0,
            ae_itt: 0,
            ae_law: 0,
            ae_jbrag: 0,
            ae_nbrag: 0,
            gf_li: 0,
            gf_loca: vec![0; 502],
            gf_ltt: 0,
            gy_awr: 0.0,
            gy_loca: vec![0; 552],
            gy_ird: 0,
            gy_ipd: 0,
            sd_nktot: 0,
            sd_nupm: 0,
            sd_loc: [0; 21],
            an_new: 0,
            an_theta: 0.0,
            an_xc: 0.0,
            an_rc: 0.0,
            an_bot: 0.0,
            an_ca: 0.0,
            an_loct: 0,
            ps_cn: 0.0,
            ps_eimax: 0.0,
            ps_ex: 0.0,
        }
    }
}

/// `fx=.8409`: a default-real literal assigned to a real(kr) variable.
pub(super) fn f32lit(x: f32) -> f64 {
    x as f64
}

impl Groupr {
    /// Group bound j of the structure `which` (1 neutron, 2 gamma).
    /// Past the last bound the Fortran reads whatever follows the array
    /// (the jgmax search of getmf6 for incident energies above the top
    /// group); here that is taken as +infinity, which ends the search.
    pub(super) fn egv(&self, which: i32, j: i32) -> f64 {
        let (n, e) = if which == 1 { (self.ngn, &self.egn) } else { (self.ngg, &self.egg) };
        if j > n + 1 {
            return f64::INFINITY;
        }
        e[j]
    }

    /// Compute the feed function from ENDF-6 File 6 data.
    pub(super) fn getmf6(
        &mut self,
        c: &mut Ctx,
        ans: &mut V2,
        ed: f64,
        enext: &mut f64,
        idisc: &mut i32,
        yld: &mut f64,
        which: i32,
        ng: i32,
        nl: &mut i32,
        iglo: &mut i32,
        ng2: &mut i32,
        nq: &mut i32,
        matd: i32,
        mfd: i32,
        mtd: i32,
        nin: i32,
        nlg: i32,
    ) {
        let ncmax = 1000;
        let emax = 1.0e10f64;
        let small = 1.0e-10f64;
        let shade = 1.1999e0f64;
        let step = 1.2e0f64;
        let up = 1.00001e0f64;
        let dn = 0.99999e0f64;
        let eps = 0.02e0f64;
        let zero = 0.0f64;
        let alight = 5.0f64;
        let ntmp: i32 = 2000000;
        let maxss = 500;
        let nssm = 9;
        let (mut nb, mut nw) = (0i32, 0i32);
        let mut yldd = 0.0f64;

        if !(ed > zero) {
            //--initialize
            self.sv.m6_idis = 0;
            let mut tmp = RV::new(ntmp as usize + SLACK);
            let max = ntmp;
            let mfn = 6;
            c.findf(matd, mfn, mtd, nin);
            c.contio(nin, 0, 0, &mut tmp.at(1), &mut nb, &mut nw);
            let nk = nint(tmp[5]);
            self.lct = nint(tmp[4]);
            self.sv.m6_lct3 = self.lct;
            self.izat = nint(tmp[1]);
            let mut zad = 1.0f64;
            if mfd == 3 {
                zad = self.izap as f64;
            }
            if mfd == 18 {
                zad = 0.0;
            }
            if mfd == 18 {
                self.lct = 1;
            }
            if mfd == 21 {
                zad = 1001.0;
            }
            if mfd == 22 {
                zad = 1002.0;
            }
            if mfd == 23 {
                zad = 1003.0;
            }
            if mfd == 24 {
                zad = 2003.0;
            }
            if mfd == 25 {
                zad = 2004.0;
            }
            if mfd == 26 {
                zad = 9999.0;
            }
            let jzad = nint(zad);
            let iy = 7;
            self.sv.m6_nss = 1;
            self.sv.m6_jjss[1] = 1;
            self.sv.m6_iyss[1] = iy;
            self.sv.m6_elo = emax;
            self.sv.m6_ehi = 0.0;
            let mut ik = 0;
            let mut l: i32;
            // 100
            loop {
                ik += 1;
                if ik > nk && mfd == 18 {
                    // 199: special exit. no mf6 gammas.
                    *enext = -1.0;
                    *idisc = 0;
                    return;
                }
                if ik > nk {
                    c.error("getmf6", "desired particle not found.", " ");
                }
                l = iy;
                c.tab1io(nin, 0, 0, &mut tmp.at(l), &mut nb, &mut nw);
                l += nw;
                while nb != 0 {
                    c.moreio(nin, 0, 0, &mut tmp.at(l), &mut nb, &mut nw);
                    l += nw;
                }
                self.sv.m6_disc102 = 0.0;
                self.sv.m6_jzap = nint(tmp[iy]);
                let jzap = self.sv.m6_jzap;
                self.zap = jzap as f64;
                self.awp = tmp[iy + 1];
                self.lip = nint(tmp[iy + 2]);
                self.law = nint(tmp[iy + 3]);
                if jzap == 0 && self.law == 2 {
                    self.sv.m6_disc102 = self.awp;
                    self.awp = 0.0;
                }
                if self.sv.m6_lct3 == 3 {
                    self.lct = 1;
                    if jzap >= 1 && jzap <= 2004 {
                        self.lct = 2;
                    }
                }
                self.aprime = self.awp;
                if self.awrp != 0.0 {
                    self.aprime = self.awp / self.awrp;
                }
                if jzap == jzad {
                    break;
                }
                if mfd == 26 && jzap > 2004 {
                    break;
                }
                c.skip6(nin, 0, 0, &mut tmp.at(l), self.law);
            }
            // 130
            //--for discrete recoil, back up to particle distribution,
            //--which is assumed to be the first subsection.
            if self.law == 4 {
                c.skiprz(nin, -2);
                c.findf(matd, mfn, mtd, nin);
                l = 1;
                c.contio(nin, 0, 0, &mut tmp.at(l), &mut nb, &mut nw);
                l += nw;
                c.tab1io(nin, 0, 0, &mut tmp.at(l), &mut nb, &mut nw);
                tmp[l] = self.zap;
                tmp[l + 1] = self.awp;
                let lf = nint(tmp[l + 3]);
                tmp[l + 3] = self.law as f64;
                if lf == 3 {
                    self.law = -4;
                }
                l += nw;
                while nb != 0 {
                    c.moreio(nin, 0, 0, &mut tmp.at(l), &mut nb, &mut nw);
                    l += nw;
                }
            }

            //--check law
            // 140
            let lbl;
            loop {
                if self.sv.m6_disc102 > zero {
                    lbl = 195;
                    break;
                }
                if self.law >= 2 && self.law <= 5 {
                    lbl = 194;
                    break;
                }
                if self.law == -4 {
                    lbl = 194;
                    break;
                }
                if self.law > 7 {
                    c.error("getmf6", "illegal law.", " ");
                }
                let nss = self.sv.m6_nss as usize;
                let jzap = self.sv.m6_jzap;
                let awr = self.awr;
                if self.law == 1 {
                    //--read in data for law 1
                    self.sv.m6_izss[nss] = l;
                    let iz = l;
                    c.tab2io(nin, 0, 0, &mut tmp.at(l), &mut nb, &mut nw);
                    let lang = nint(tmp[l + 2]);
                    let lep = nint(tmp[l + 3]);
                    self.sv.m6_ne = nint(tmp[l + 5]);
                    self.sv.m6_int = nint(tmp[l + 7]);
                    // force unit base interpolation for smoother scattering source
                    if self.sv.m6_int == 2 {
                        self.sv.m6_int = 22;
                    }
                    let mut nn: i32;
                    l += nw;
                    let mut nnn = 0;
                    for ie in 1..=self.sv.m6_ne {
                        let ilo = l;
                        let mut jlo = ilo;
                        let ii = self.sv.m6_jjss[nss] + ie - 1;
                        if ii > maxss {
                            c.error("getmf6", "too many subsection energy points", " ");
                        }
                        self.sv.m6_jloss[ii as usize] = jlo;
                        c.listio(nin, 0, 0, &mut tmp.at(l), &mut nb, &mut nw);
                        if tmp[l + 1] < self.sv.m6_elo {
                            self.sv.m6_elo = tmp[l + 1];
                        }
                        if tmp[l + 1] > self.sv.m6_ehi {
                            self.sv.m6_ehi = tmp[l + 1];
                        }
                        nn = nint(tmp[l + 5]);
                        let mut ncyc = nint(tmp[l + 4] / nn as f64);
                        let int = self.sv.m6_int;
                        if ie > 1 && nn != nnn && int >= 11 && int <= 20 {
                            self.sv.m6_int = int + 10;
                            tmp[iz + 7] = self.sv.m6_int as f64;
                            c.mess(
                                "getmf6",
                                "bad grids for corresponding-point interpolation-",
                                "changing to unit-base interpolation",
                            );
                        }
                        nnn = nn;
                        l += nw;
                        while nb != 0 {
                            if l > ntmp {
                                c.error("getmf6", "storage exceeded", " ");
                            }
                            c.moreio(nin, 0, 0, &mut tmp.at(l), &mut nb, &mut nw);
                            l += nw;
                        }
                        if nn == 2 && lep == 2 {
                            let mut j = 0;
                            if tmp[ilo + 6 + ncyc] < tmp[ilo + 1] / powc(awr + 1.0, 2) {
                                j = 1;
                                tmp[ilo + 7] = 0.0;
                                tmp[ilo + 8] = tmp[ilo + 1] / powc(awr + 1.0, 2);
                                tmp[ilo + 9] = 2.0 / tmp[ilo + 8];
                            } else if tmp[ilo + 7] != 0.0 && tmp[ilo + 9] == 0.0 {
                                j = 1;
                                tmp[ilo + 9] = tmp[ilo + 7];
                                tmp[ilo + 7] = 0.0;
                            }
                            if j != 0 {
                                c.io.so(&format!(" patching low-energy distribution at{}", ffmt::e1p(tmp[ilo + 1], 10, 3)));
                            }
                        }
                        if self.ismooth > 0 && jzap == 1 && lep == 1 {
                            let fx = f32lit(0.8409);
                            let ex = 40.0f64;
                            ncyc = nint(tmp[ilo + 3]) + 2;
                            let mut cx = tmp[ilo + 6 + ncyc] * tmp[ilo + 7];
                            let mut nx = nint(tmp[ilo + 4]);
                            let mut n = nint(tmp[ilo + 5]);
                            while n > 2 {
                                let cxx = cx + tmp[ilo + 7 + ncyc] * (tmp[ilo + 6 + 2 * ncyc] - tmp[ilo + 6 + ncyc]);
                                if (cxx / pow(tmp[ilo + 6 + 2 * ncyc], 1.5) - cx / pow(tmp[ilo + 6 + ncyc], 1.5)).abs()
                                    > cx / pow(tmp[ilo + 6 + ncyc], 1.5) / 50.0
                                {
                                    break;
                                }
                                tmp[ilo + 7] = (tmp[ilo + 7] * tmp[ilo + 6 + ncyc]
                                    + tmp[ilo + 7 + ncyc] * (tmp[ilo + 6 + 2 * ncyc] - tmp[ilo + 6 + ncyc]))
                                    / tmp[ilo + 6 + 2 * ncyc];
                                for ix in 1..=(nx - 2 * ncyc) {
                                    tmp[ilo + 5 + ix + ncyc] = tmp[ilo + 5 + ix + 2 * ncyc];
                                }
                                cx = cxx;
                                nx -= ncyc;
                                n -= 1;
                            }
                            c.io.so(&format!(
                                " extending histogram as sqrt(E) below{} eV for E={} eV",
                                ffmt::e1p(tmp[ilo + 6 + ncyc], 10, 2),
                                ffmt::e1p(tmp[ilo + 1], 10, 2)
                            ));
                            while tmp[ilo + 6 + ncyc] > ex {
                                let mut ix = nx;
                                while ix >= 1 {
                                    tmp[ilo + 5 + ncyc + ix] = tmp[ilo + 5 + ix];
                                    tmp[ilo + 5 + ncyc + ix] = sigfig(tmp[ilo + 5 + ncyc + ix], 7, 0);
                                    ix -= 1;
                                }
                                tmp[ilo + 6 + ncyc] = fx * tmp[ilo + 6 + 2 * ncyc];
                                tmp[ilo + 6 + ncyc] = sigfig(tmp[ilo + 6 + ncyc], 7, 0);
                                let val = tmp[ilo + 7];
                                tmp[ilo + 7] = sqrt(fx) * val;
                                tmp[ilo + 7] = sigfig(tmp[ilo + 7], 7, 0);
                                tmp[ilo + 7 + ncyc] = (1.0 - fx * sqrt(fx)) * val / (1.0 - fx);
                                tmp[ilo + 7 + ncyc] = sigfig(tmp[ilo + 7 + ncyc], 7, 0);
                                nx += ncyc;
                                n += 1;
                                tmp[ilo + 4] = nx as f64;
                                tmp[ilo + 5] = n as f64;
                            }
                            l = ilo + 6 + nx;
                        } else if self.ismooth > 0 && jzap == 1 && lep == 2 {
                            ncyc = nint(tmp[ilo + 3]) + 2;
                            let mut nx = nint(tmp[ilo + 4]);
                            let mut n = nint(tmp[ilo + 5]);
                            c.io.so(&format!(
                                " extending lin-lin as sqrt(E) below{} eV for E={} eV",
                                ffmt::e1p(tmp[ilo + 6 + ncyc], 10, 2),
                                ffmt::e1p(tmp[ilo + 1], 10, 2)
                            ));
                            let ex = 40.0f64;
                            let fx = 0.50f64;
                            nn = 0;
                            let cx = (tmp[ilo + 6 + ncyc] - tmp[ilo + 6]) * (tmp[ilo + 7 + ncyc] + tmp[ilo + 7]) / 2.0;
                            let mut cxx = 0.0f64;
                            for i in 1..=(ncyc - 1) {
                                tmp[ilo + 6 + i] = 0.0;
                            }
                            while tmp[ilo + 6 + ncyc] > ex {
                                nn += 1;
                                let mut ix = nx;
                                while ix >= ncyc {
                                    tmp[ilo + 5 + ncyc + ix] = tmp[ilo + 5 + ix];
                                    ix -= 1;
                                }
                                tmp[ilo + 6 + ncyc] = fx * tmp[ilo + 6 + 2 * ncyc];
                                tmp[ilo + 6 + ncyc] = sigfig(tmp[ilo + 6 + ncyc], 6, 0);
                                for i in 1..=(ncyc - 1) {
                                    tmp[ilo + 6 + ncyc + i] =
                                        tmp[ilo + 6 + 2 * ncyc + i] * sqrt(tmp[ilo + 6 + ncyc] / tmp[ilo + 6 + 2 * ncyc]);
                                }
                                if nn > 1 {
                                    cxx = cxx
                                        + (tmp[ilo + 6 + 2 * ncyc] - tmp[ilo + 6 + ncyc])
                                            * (tmp[ilo + 7 + 2 * ncyc] + tmp[ilo + 7 + ncyc])
                                            / 2.0;
                                }
                                nx += ncyc;
                                n += 1;
                                tmp[ilo + 4] = nx as f64;
                                tmp[ilo + 5] = n as f64;
                            }
                            cxx = cxx + tmp[ilo + 6 + ncyc] * tmp[ilo + 7 + ncyc] / 2.0;
                            let dx = tmp[ilo + 6 + (nn + 1) * ncyc] - tmp[ilo + 6 + nn * ncyc];
                            let mut rn = 1.0f64;
                            if cxx + dx * tmp[ilo + 7 + nn * ncyc] / 2.0 != zero {
                                rn = (cx - dx * tmp[ilo + 7 + (nn + 1) * ncyc] / 2.0)
                                    / (cxx + dx * tmp[ilo + 7 + nn * ncyc] / 2.0);
                            }
                            for i in 1..=nn {
                                for j in 1..=(ncyc - 1) {
                                    tmp[ilo + 6 + i * ncyc + j] = rn * tmp[ilo + 6 + i * ncyc + j];
                                }
                            }
                            l = ilo + 6 + nx;
                        }
                        if self.lct == 2 || (self.lct == 3 && self.awp <= alight) {
                            self.cm2lab(c, ilo, &mut jlo, &mut l, &mut tmp, *nl, lang, lep, max);
                            for i in jlo..=l {
                                tmp[ilo + i - jlo] = tmp[i];
                            }
                            l = ilo + l - jlo;
                        }
                        self.sv.m6_jlo = jlo;
                    }
                } else if self.law == 6 {
                    //--read in data for law 6.  convert them to law 1.
                    c.contio(nin, 0, 0, &mut tmp.at(l), &mut nb, &mut nw);
                    let apsx = tmp[l];
                    let npsx = nint(tmp[l + 5]);
                    tmp[l] = 0.0;
                    tmp[l + 1] = 0.0;
                    tmp[l + 2] = 0.0;
                    tmp[l + 3] = 0.0;
                    let mut lang = 0;
                    let lep = 2;
                    self.lct = 2;
                    let irr = nint(tmp[iy + 4]);
                    let npp = nint(tmp[iy + 5]);
                    let mut enow = tmp[iy + 6 + 2 * irr] / shade;
                    self.sv.m6_pspmax = tmp[iy + 4 + 2 * irr + 2 * npp];
                    self.sv.m6_izss[nss] = l;
                    l += 8;
                    self.sv.m6_ne = 0;
                    while enow < self.sv.m6_pspmax {
                        self.sv.m6_ne += 1;
                        let ii = self.sv.m6_jjss[nss] + self.sv.m6_ne - 1;
                        if ii > maxss {
                            c.error("getmf6", "too many subsection energy points", " ");
                        }
                        let ilo = l;
                        self.sv.m6_jloss[ii as usize] = ilo;
                        tmp[ilo] = apsx;
                        enow = enow * step;
                        if enow > self.sv.m6_pspmax {
                            enow = self.sv.m6_pspmax;
                        }
                        tmp[ilo + 1] = enow;
                        if enow < self.sv.m6_elo {
                            self.sv.m6_elo = enow;
                        }
                        tmp[ilo + 2] = 0.0;
                        tmp[ilo + 3] = 0.0;
                        tmp[ilo + 4] = 0.0;
                        tmp[ilo + 5] = npsx as f64;
                        l += 6;
                        let mut jlo = 0;
                        self.cm2lab(c, ilo, &mut jlo, &mut l, &mut tmp, *nl, lang, lep, max);
                        for i in jlo..=l {
                            tmp[ilo + i - jlo] = tmp[i];
                        }
                        l = ilo + l - jlo;
                        self.sv.m6_jlo = jlo;
                    }
                    let i = self.sv.m6_izss[nss];
                    lang = 1;
                    tmp[i + 2] = lang as f64;
                    tmp[i + 3] = lep as f64;
                    tmp[i + 4] = 1.0;
                    tmp[i + 5] = self.sv.m6_ne as f64;
                    tmp[i + 6] = self.sv.m6_ne as f64;
                    tmp[i + 7] = 22.0;
                    self.law = 1;
                    tmp[iy + 3] = self.law as f64;
                } else if self.law == 7 {
                    //--read in data for law 7
                    self.sv.m6_izss[nss] = l;
                    c.tab2io(nin, 0, 0, &mut tmp.at(l), &mut nb, &mut nw);
                    self.sv.m6_ne = nint(tmp[l + 5]);
                    l += nw;
                    for ie in 1..=self.sv.m6_ne {
                        let ilo = l;
                        let mut jlo = ilo;
                        let ii = self.sv.m6_jjss[nss] + ie - 1;
                        if ii > maxss {
                            c.error("getmf6", "too many subsection energy points", " ");
                        }
                        self.sv.m6_jloss[ii as usize] = ilo;
                        c.tab2io(nin, 0, 0, &mut tmp.at(l), &mut nb, &mut nw);
                        if tmp[l + 1] < self.sv.m6_elo {
                            self.sv.m6_elo = tmp[l + 1];
                        }
                        if tmp[l + 1] > self.sv.m6_ehi {
                            self.sv.m6_ehi = tmp[l + 1];
                        }
                        let nmu = nint(tmp[l + 5]);
                        l += nw;
                        for _i in 1..=nmu {
                            c.tab1io(nin, 0, 0, &mut tmp.at(l), &mut nb, &mut nw);
                            l += nw;
                            while nb != 0 {
                                if l > ntmp {
                                    c.error("getmf6", "storage exceeded", " ");
                                }
                                c.moreio(nin, 0, 0, &mut tmp.at(l), &mut nb, &mut nw);
                                l += nw;
                            }
                        }
                        self.ll2lab(c, ilo, &mut jlo, &mut l, &mut tmp, *nl, max);
                        for i in jlo..=l {
                            tmp[ilo + i - jlo] = tmp[i];
                        }
                        l = ilo + l - jlo;
                        self.sv.m6_jlo = jlo;
                    }
                }

                //--check for multiple particle subsections
                if ik == nk {
                    lbl = 195;
                    break;
                }
                ik += 1;
                c.tab1io(nin, 0, 0, &mut tmp.at(l), &mut nb, &mut nw);
                let l1 = l;
                l += nw;
                while nb != 0 {
                    if l > ntmp {
                        c.error("getmf6", "storage exceeded", " ");
                    }
                    c.moreio(nin, 0, 0, &mut tmp.at(l), &mut nb, &mut nw);
                    l += nw;
                }
                if nint(tmp[l1]) != self.sv.m6_jzap {
                    lbl = 195;
                    break;
                }
                self.law = nint(tmp[l1 + 3]);
                self.sv.m6_nss += 1;
                if self.sv.m6_nss > nssm {
                    c.error("getmf6", "too many subsections for one particle", " ");
                }
                let nss = self.sv.m6_nss as usize;
                self.sv.m6_iyss[nss] = l1;
                self.sv.m6_jjss[nss] = self.sv.m6_jjss[nss - 1] + self.sv.m6_ne;
            }
            if lbl == 194 {
                //--discrete two-body scattering
                let nss = self.sv.m6_nss as usize;
                self.sv.m6_izss[nss] = l;
                let mut elo = self.sv.m6_elo;
                self.getdis(c, ed, &mut elo, idisc, &mut yldd, ans, nl, ng2, iglo, nq, nin, nlg);
                self.sv.m6_elo = elo;
                l += ncmax;
            }
            //--initialization complete
            // 195
            if self.sv.m6_nss > 1 {
                c.mess("getmf6", "there are multiple subsections in mf6", "for this emitted particle");
            }
            let nc = l - 1;
            self.ddmf6 = alloc(nc);
            for i in 1..=nc {
                self.ddmf6[i] = tmp[i];
            }
            *enext = self.sv.m6_elo;
            *idisc = 0;
            return;
        }

        //--normal entry
        // 200
        for j in 1..=ng {
            for l in 1..=*nl {
                ans[(l, j)] = 0.0;
            }
        }
        *enext = emax;
        *yld = 0.0;
        let mut iss = 1;
        let mut dd = std::mem::take(&mut self.ddmf6);
        let mut jgmax = 0i32;
        let mut term = [0.0f64; MXLG + 1];
        loop {
            // 207
            let iy = self.sv.m6_iyss[iss as usize];
            let iz = self.sv.m6_izss[iss as usize];
            self.lip = nint(dd[iy + 2]);
            self.law = nint(dd[iy + 3]);
            if self.sv.m6_disc102 > zero {
                // 600: discrete relativistic capture gamma
                self.ddmf6 = dd;
                let d102 = self.sv.m6_disc102;
                self.gam102(ans, ed, enext, d102, self.law, *nl, iglo, ng2, nq);
                return;
            }
            if self.law >= 2 && self.law <= 5 {
                // 500: laws 2 thru 5 -- discrete two-body scattering
                self.ddmf6 = dd;
                self.getdis(c, ed, enext, idisc, &mut yldd, ans, nl, ng2, iglo, nq, nin, nlg);
                return;
            }

            //--interpolate for fractional probability
            let mut ip = 2;
            let mut ir = 1;
            let mut pe = 0.0f64;
            let mut eihi = 0.0f64;
            c.terpa(&mut pe, ed, &mut eihi, idisc, &dd.at(iy), &mut ip, &mut ir);
            if (eihi - *enext).abs() < *enext * small && *idisc > self.sv.m6_idis {
                self.sv.m6_idis = *idisc;
            }
            if eihi < *enext * (1.0 - small) {
                self.sv.m6_idis = *idisc;
            }
            if eihi < *enext * (1.0 - small) {
                *enext = eihi;
            }
            let mut go450 = pe == zero;
            if !go450 {
                *yld = *yld + pe;
                let lang = nint(dd[iz + 2]);
                let lep = nint(dd[iz + 3]);
                let law = self.law;
                let lct = self.lct;
                self.sv.m6_lepn = lep;
                if law == 1 && lct == 2 {
                    self.sv.m6_lepn = 2;
                }
                if law == 7 {
                    self.sv.m6_lepn = 2;
                }
                self.sv.m6_langn = lang;
                if law == 1 && lct == 2 {
                    self.sv.m6_langn = 1;
                }
                if law == 1 && lct == 3 && self.awp < alight {
                    self.sv.m6_langn = 1;
                }
                if law == 7 {
                    self.sv.m6_langn = 1;
                }

                //--find jgmax
                jgmax = ng + 1;
                if self.q <= zero && self.sv.m6_jzap != 0 {
                    jgmax = 1;
                    while self.egv(which, jgmax + 1) < ed * (1.0 - small) {
                        jgmax += 1;
                    }
                }

                //--find desired energy panel
                self.sv.m6_nne = nint(dd[iz + 5]);
                self.sv.m6_int = nint(dd[iz + 7]);
                if self.sv.m6_int == 2 {
                    self.sv.m6_int = 22;
                }
                let jj = self.sv.m6_jjss[iss as usize];
                let mut found = false;
                for ie in 1..=self.sv.m6_nne {
                    self.sv.m6_jlo = self.sv.m6_jloss[(jj + ie - 1) as usize];
                    self.sv.m6_jhi = self.sv.m6_jloss[(jj + ie) as usize];
                    let (jlo, jhi) = (self.sv.m6_jlo, self.sv.m6_jhi);
                    if ed >= dd[jlo + 1] * (1.0 - small) && ed <= dd[jhi + 1] * (1.0 + small) {
                        self.sv.m6_elo = dd[jlo + 1];
                        self.sv.m6_ehi = dd[jhi + 1];
                        found = true;
                        break;
                    }
                }
                if !found {
                    go450 = true;
                }
            }
            if !go450 {
                //--set up loop over secondary energy
                // 300
                let (jlo, jhi) = (self.sv.m6_jlo, self.sv.m6_jhi);
                let (law, int, langn, lepn) = (self.law, self.sv.m6_int, self.sv.m6_langn, self.sv.m6_lepn);
                let mut ep = 0.0f64;
                let mut epnext = 0.0f64;
                let mut terml = self.sv.m6_terml;
                {
                    let clo = &dd.0[(jlo - 1) as usize..];
                    let chi = &dd.0[(jhi - 1) as usize..];
                    self.f6lab(c, ep, &mut epnext, &mut terml[1..], *nl, law, int, langn, lepn, clo, chi, ed);
                }
                let mut jg = 1;
                // 320
                loop {
                    let mut en = epnext;
                    if en >= emax * (1.0 - small) {
                        break;
                    }
                    if jg < ng && ep >= self.egv(which, jg + 1) * (1.0 - small) {
                        jg += 1;
                    }
                    if jg < ng && self.egv(which, jg + 1) < en * (1.0 - small) {
                        en = self.egv(which, jg + 1);
                    }
                    {
                        let clo = &dd.0[(jlo - 1) as usize..];
                        let chi = &dd.0[(jhi - 1) as usize..];
                        self.f6lab(c, en, &mut epnext, &mut term[1..], *nl, law, int, langn, lepn, clo, chi, ed);
                    }
                    //--trapazoidal integration over panel.
                    for l in 1..=*nl {
                        ans[(l, jg)] = ans[(l, jg)] + pe * (en - ep) * (term[l as usize] + terml[l as usize]) / 2.0;
                        terml[l as usize] = term[l as usize];
                    }
                    ep = en;
                }
                self.sv.m6_terml = terml;

                //--add in discrete energies for lab, if any.
                // 420
                if self.lct != 2 {
                    if self.lct != 3 || self.awp >= alight {
                        let ndlo = nint(dd[jlo + 2]);
                        let nplo = nint(dd[jlo + 5]);
                        let nclo = nint(dd[jlo + 4]) / nplo;
                        let nphi = nint(dd[jhi + 5]);
                        let nchi = nint(dd[jhi + 4]) / nphi;
                        if ndlo > 0 {
                            let awr = self.awr;
                            let (elo, ehi) = (self.sv.m6_elo, self.sv.m6_ehi);
                            for i in 1..=ndlo {
                                let llo = jlo + 6 + nclo * (i - 1);
                                let lhi = jhi + 6 + nchi * (i - 1);
                                let mut el = dd[llo];
                                let mut eh = dd[lhi];
                                if el < zero {
                                    el = ed * awr / (awr + 1.0) - el;
                                }
                                if eh < zero {
                                    eh = ed * awr / (awr + 1.0) - eh;
                                }
                                let (mut e0, mut g0) = (0.0f64, 0.0f64);
                                c.terp1(elo, el, ehi, eh, ed, &mut e0, 2);
                                c.terp1(elo, dd[llo + 1], ehi, dd[lhi + 1], ed, &mut g0, 2);
                                for j in 1..=ng {
                                    let mut e1 = self.egv(which, j);
                                    if j == 1 {
                                        e1 = 0.0;
                                    }
                                    let mut e2 = self.egv(which, j + 1);
                                    if j == ng {
                                        e2 = emax;
                                    }
                                    if e0 >= e1 * (1.0 - small) && e0 < e2 * (1.0 - small) {
                                        ans[(1, j)] = ans[(1, j)] + pe * g0;
                                    }
                                }
                            }
                        }
                    }
                }
            }

            //--finished with energy loop.
            // 450
            *iglo = 1;
            *ng2 = ng;
            *nq = 0;
            let ehi = self.sv.m6_ehi;
            if ed < ehi * (1.0 - small) && ehi < *enext * (1.0 - small) {
                *enext = ehi;
            }

            //--loop back to process additional particles, if any.
            iss += 1;
            if iss > self.sv.m6_nss {
                break;
            }
        }
        self.ddmf6 = dd;

        //--check normalization
        //--move illegal upscatters to ingroup
        // 700
        if ed <= zero {
            return;
        }
        if *yld == zero {
            return;
        }
        let mut test = 0.0f64;
        let jj = jgmax - *iglo + 1;
        for i in 1..=*ng2 {
            test = test + ans[(1, i)];
            if i > jj && nint(self.zap) == self.izap && self.q <= zero {
                for l in 1..=*nl {
                    ans[(l, jj)] = ans[(l, jj)] + ans[(l, i)];
                    ans[(l, i)] = 0.0;
                }
            }
        }
        if test == zero {
            return;
        }
        let (elo, ehi) = (self.sv.m6_elo, self.sv.m6_ehi);
        if (test - *yld).abs() > eps && (ed < up * elo || ed > dn * ehi) {
            c.io.so(&format!("     normalization {}{}", ffmt::e1p(ed, 12, 4), ffmt::e1p(test, 15, 6)));
        }
        test = *yld / test;
        for l in 1..=*nl {
            for i in 1..=*ng2 {
                ans[(l, i)] = test * ans[(l, i)];
            }
        }
    }

    /// Convert the cm distribution starting at inow into the laboratory
    /// frame. The lab distribution starts at jnow.
    pub(super) fn cm2lab(
        &mut self,
        c: &mut Ctx,
        inow: i32,
        jnow: &mut i32,
        lnow: &mut i32,
        cc: &mut RV,
        nl: i32,
        lang: i32,
        lep: i32,
        max: i32,
    ) {
        const IMAX: usize = 15;
        let tol = 0.01e0f64;
        let eps = 1.0e-8f64;
        let small = 1.0e-10f64;
        let emax = 1.0e10f64;
        let zero = 0.0f64;
        let mut term = [0.0f64; MXLG + 1];
        let mut x = [0.0f64; IMAX + 2];
        let mut y = vec![[0.0f64; MXLG + 1]; IMAX + 2];
        let nlu = nl as usize;

        //--jnow points to the new lab distribution
        *jnow = *lnow;
        cc[*jnow] = 0.0;
        cc[*jnow + 1] = cc[inow + 1];
        *lnow = *jnow + 6;

        //--prime the stack with first point
        let mut sum = 0.0f64;
        x[2] = 0.0;
        let mut epnext = 0.0f64;
        let iu = (inow - 1) as usize;
        self.f6cm(c, x[2], &mut epnext, &mut term[1..], nl, lang, lep, &cc.0[iu..]);
        for l in 1..=nlu {
            y[2][l] = term[l];
        }

        //--prime the stack with the top of the next panel
        let mut xm = 0.0f64;
        while epnext < emax * (1.0 - small) {
            x[1] = epnext;
            self.f6cm(c, x[1], &mut epnext, &mut term[1..], nl, lang, lep, &cc.0[iu..]);
            for l in 1..=nlu {
                y[1][l] = term[l];
            }

            //--adaptive integration over panel for cm.
            let mut i = 2usize;
            while i > 1 {
                //--test for convergence.
                let mut iflag = 0;
                if i > 1 && i < IMAX {
                    let da = (x[i - 1] - x[i]) * (y[i - 1][1] + y[i][1]) / 2.0;
                    let mut rat = 2.0f64;
                    if x[i] != zero {
                        rat = x[i - 1] / x[i];
                    }
                    if rat > 1.0 + tol {
                        if da >= eps {
                            xm = (x[i - 1] + x[i]) / 2.0;
                            if xm != x[i - 1] && xm != x[i] {
                                let mut epn = x[i];
                                self.f6cm(c, xm, &mut epn, &mut term[1..], nl, lang, lep, &cc.0[iu..]);
                                for l in 1..=nlu {
                                    let ym = (y[i - 1][l] + y[i][l]) / 2.0;
                                    let test = l as f64 * tol * (ym.abs());
                                    let dely = (term[l] - ym).abs();
                                    if dely > test {
                                        iflag = 1;
                                    }
                                }
                            }
                        }
                    }
                }

                if iflag == 1 {
                    //--fails tests. add midpoint to stack and continue.
                    i += 1;
                    x[i] = x[i - 1];
                    for l in 1..=nlu {
                        y[i][l] = y[i - 1][l];
                    }
                    x[i - 1] = xm;
                    for l in 1..=nlu {
                        y[i - 1][l] = term[l];
                    }
                } else {
                    //--passes tests. take top point off of the stack.
                    cc[*lnow] = x[i];
                    for l in 1..=nlu {
                        cc[*lnow + l as i32] = y[i][l];
                    }
                    *lnow = *lnow + 1 + nl;
                    if *lnow + nl + 1 > max {
                        c.error("cm2lab", "storage exceeded.", " ");
                    }
                    sum = sum + (x[i - 1] - x[i]) * (y[i - 1][1] + y[i][1]) / 2.0;
                    i -= 1;
                }
            }

            //--continue looping over panels.
            x[2] = x[1];
            for l in 1..=nlu {
                y[2][l] = y[1][l];
            }
        }
        if (sum - 1.0).abs() > tol {
            let s = format!(
                "mt={} e={} lab sum={}",
                ffmt::i(c.mth, 3),
                ffmt::e1p(cc[inow + 1], 10, 3),
                ffmt::f(sum, 6, 3)
            );
            c.mess("cm2lab", "lab normalization problem ", &s);
        }

        //--linearization and conversion complete.
        let nw = *lnow - *jnow - 6;
        let np = nw / (nl + 1);
        cc[*jnow + 4] = nw as f64;
        cc[*jnow + 5] = np as f64;
        cc[*jnow + 2] = 0.0;
        cc[*jnow + 3] = (nl - 1) as f64;
    }

    /// Legendre coefficients of the double-differential cross section at
    /// lab secondary energy ep from File 6 cm data in cnow (cnow[0] is
    /// cnow(1)). Call with ep=0 to initialize.
    pub(super) fn f6cm(
        &mut self,
        c: &mut Ctx,
        ep: f64,
        epnext: &mut f64,
        term: &mut [f64],
        nl: i32,
        lang: i32,
        lep: i32,
        cnow: &[f64],
    ) {
        const IMAX: usize = 10;
        let tol = 0.005e0f64;
        let emax = 1.0e10f64;
        let check = 0.99999e0f64;
        let rdn = 0.999995e0f64;
        let rup = 1.000005e0f64;
        let amil = 1.0e-6f64;
        let athou = 1.0e-3f64;
        let dn = 0.999999e0f64;
        let up = 1.000001e0f64;
        let step = 0.05e0f64;
        let one = 1.0e0f64;
        let tiny = 1.0e-4f64;
        let small = 1.0e-10f64;
        let zero = 0.0f64;
        let cn = |i: i32| cnow[(i - 1) as usize];
        let nlu = nl as usize;
        // term is 0-based here: term[l-1] is term(l)
        let mut p = [0.0f64; MXLG + 1];
        let mut x = [0.0f64; IMAX + 2];
        let mut y = vec![[0.0f64; MXLG + 1]; IMAX + 2];
        let mut yt = [0.0f64; MXLG + 1];
        let (mut epn, mut epm, mut w) = (0.0f64, 0.0f64, 0.0f64);

        //--initialize.
        if nl > MXLG as i32 {
            c.error("f6cm", "nl>mxlg", " ");
        }
        if ep == zero {
            let s = &mut self.sv;
            s.cm_e = cn(2);
            s.cm_ndnow = nint(cn(3));
            s.cm_npnow = nint(cn(6));
            s.cm_ncnow = nint(cn(5)) / s.cm_npnow;
            s.cm_xc = self.aprime / powc(self.awr + 1.0, 2);
            *epnext = emax;
            self.sv.cm_elmax = 0.0;
            let mut epnn = emax;
            self.sv.cm_epmax = 0.0;
            let e = self.sv.cm_e;
            let xc = self.sv.cm_xc;
            if self.sv.cm_npnow != self.sv.cm_ndnow {
                if lang > 0 {
                    self.f6ddx(c, ep, &mut epn, &mut epm, w, cnow, lang, lep);
                }
                if lang == 0 {
                    self.f6psp(c, ep, &mut epn, &mut epm, w, e, cnow);
                }
                if epn < epnn * (1.0 - small) {
                    epnn = epn;
                }
                if epm > self.sv.cm_epmax * (1.0 + small) {
                    self.sv.cm_epmax = epm;
                }
                let epmax = self.sv.cm_epmax;
                let mut elmax = e * powc(sqrt(epmax / e) + sqrt(xc), 2);
                elmax = up * elmax;
                self.sv.cm_elmax = elmax;
                epn = step * elmax;
                if epmax * up < xc * e {
                    epn = e * powc(sqrt(epmax / e) - sqrt(xc), 2);
                }
                if epn < *epnext {
                    *epnext = epn;
                }
                self.sv.cm_eps = tiny / elmax;
                *epnext = dn * *epnext;
            }
            if self.sv.cm_ndnow != 0 {
                let _ = self.f6dis(c, 0, &mut epn, &mut epm, w, cnow, lang);
                for i in 1..=self.sv.cm_ndnow {
                    let _ = self.f6dis(c, i, &mut epn, &mut epm, w, cnow, lang);
                    if epn < epnn * (1.0 - small) {
                        epnn = epn;
                    }
                    if epn > self.sv.cm_epmax * (1.0 + small) {
                        self.sv.cm_epmax = epn;
                    }
                    let mut epx = dn * e * powc(sqrt(epn / e) - sqrt(xc), 2);
                    if epx < *epnext * (1.0 - small) {
                        *epnext = epx;
                    }
                    epx = up * e * powc(sqrt(epn / e) + sqrt(xc), 2);
                    if epx > self.sv.cm_elmax * (1.0 + small) {
                        self.sv.cm_elmax = epx;
                    }
                }
            }
            self.sv.cm_na = nl - 1;
            for l in 0..nlu {
                term[l] = 0.0;
            }
            return;
        }

        //--normal entry.
        for l in 0..nlu {
            term[l] = 0.0;
        }
        let (e, xc, na, elmax, epmax, eps) =
            (self.sv.cm_e, self.sv.cm_xc, self.sv.cm_na, self.sv.cm_elmax, self.sv.cm_epmax, self.sv.cm_eps);
        let (ndnow, npnow, ncnow) = (self.sv.cm_ndnow, self.sv.cm_npnow, self.sv.cm_ncnow);
        let epnxt: f64;
        if ndnow == npnow {
            // 330: check for contributions from delta functions
            let mut ep_nxt: f64;
            if ndnow != 0 {
                ep_nxt = emax;
                for i in 1..=ndnow {
                    let ll = 7 + ncnow * (i - 1);
                    let epp = cn(ll);
                    let w = (ep - epp - xc * e) / (2.0 * sqrt(xc * e) * sqrt(epp));
                    let u = (xc * e + ep - epp) / (2.0 * sqrt(xc * e) * sqrt(ep));
                    let wmin = -1.0;
                    let wmax = 1.0;
                    if w < wmin {
                        epn = e * powc(sqrt(epp / e) - sqrt(xc), 2);
                        if ep < check * epn {
                            epn = rdn * epn;
                        } else {
                            epn = rup * epn;
                        }
                        if epn < ep_nxt * (1.0 - small) {
                            ep_nxt = epn;
                        }
                    } else if w <= wmax {
                        legndr(u, &mut p, na);
                        let mut f = self.f6dis(c, i, &mut epn, &mut epm, w, cnow, lang);
                        f = f / (2.0 * sqrt(xc * e * epp));
                        epn = e * powc(sqrt(epp / e) + sqrt(xc), 2);
                        if ep < check * epn {
                            epn = rdn * epn;
                        } else {
                            epn = rup * epn;
                        }
                        if epn < ep_nxt * (1.0 - small) {
                            ep_nxt = epn;
                        }
                        for l in 0..nlu {
                            term[l] = term[l] + f * p[l];
                        }
                    }
                }
            } else {
                ep_nxt = self.cm_epnxt_garbage();
            }
            epnxt = ep_nxt;
        } else {
            let xx = ep / e;
            let mut cc = xc / xx;
            let mut cv = sqrt(cc);
            let mut umin = (1.0 + cc - epmax / ep) / (2.0 * cv);
            if !(umin < check) {
                if cv <= one {
                    // 420
                    *epnext = emax;
                    return;
                }
                // 320
                epnxt = (*epnext + elmax) / 2.0;
                self.f6cm_tail(term, nl, tol, epnext, epnxt);
                return;
            }
            // 215
            if umin < -one {
                umin = -one;
            }
            let dm = (1.0 - umin) * athou + amil;
            let mut epnn = emax;

            //--for the continuum part of the distribution,
            //--integrate over the trajectory ep=constant
            //--using adaptive integration on lab cosine.
            let umax = 1.0f64;
            let mut un = umax;
            x[2] = un;
            let mut s = 0.0f64;
            let mut u = 0.0f64;
            while x[2] > umin {
                //--load the boundaries for this cosine panel.
                let mut j = 2;
                while j == 2 {
                    u = un;
                    let mut yy = 1.0 + cc - 2.0 * cv * u;
                    if yy < amil {
                        yy = amil;
                        cv = u - athou;
                        cc = cv * cv;
                    }
                    let epp = yy * ep;
                    w = (u - cv) / sqrt(yy);
                    if lang > 0 {
                        s = self.f6ddx(c, epp, &mut epn, &mut epm, w, cnow, lang, lep);
                    }
                    if lang == 0 {
                        s = self.f6psp(c, epp, &mut epn, &mut epm, w, e, cnow);
                    }
                    if u == umax && epn < epnn {
                        epnn = epn;
                    }
                    legndr(u, &mut p, na);
                    j = 1;
                    if u == umax {
                        j = 2;
                    }
                    x[j] = u;
                    for l in 1..=nlu {
                        y[j][l] = p[l - 1] * s / sqrt(yy);
                    }
                    un = u - (epn - epp) / (2.0 * cv * ep);
                    if un > u - tiny {
                        un = u - tiny;
                    }
                    if un < umin + tiny / 10.0 {
                        un = umin;
                    }
                }

                //--do an adaptive reconstruction for this panel.
                let mut i = 2usize;
                while i > 1 {
                    //--test for convergence.
                    let mut iflag = 0;
                    if i > 1 && i < IMAX {
                        let dx = x[i] - x[i - 1];
                        if dx >= dm {
                            let da = dx * (y[i][1] + y[i - 1][1]) / 2.0;
                            if da >= eps {
                                u = (x[i - 1] + x[i]) / 2.0;
                                let yy = 1.0 + cc - 2.0 * cv * u;
                                let epp = yy * ep;
                                w = (u - cv) / sqrt(yy);
                                if lang > 0 {
                                    s = self.f6ddx(c, epp, &mut epn, &mut epm, w, cnow, lang, lep);
                                }
                                if lang == 0 {
                                    s = self.f6psp(c, epp, &mut epn, &mut epm, w, e, cnow);
                                }
                                legndr(u, &mut p, na);
                                for l in 1..=nlu {
                                    yt[l] = p[l - 1] * s / sqrt(yy);
                                }
                                for l in 1..=nlu {
                                    let ym = (y[i - 1][l] + y[i][l]) / 2.0;
                                    let dy = (yt[l] - ym).abs();
                                    if dy > l as f64 * tol * ym.abs() + eps {
                                        iflag = 1;
                                    }
                                }
                            }
                        }
                    }

                    if iflag == 1 {
                        //--fails test. add midpoint to stack.
                        i += 1;
                        x[i] = x[i - 1];
                        for l in 1..=nlu {
                            y[i][l] = y[i - 1][l];
                        }
                        x[i - 1] = u;
                        for l in 1..=nlu {
                            y[i - 1][l] = yt[l];
                        }
                    } else {
                        //--passes test. take top point off of stack.
                        for l in 1..=nlu {
                            term[l - 1] = term[l - 1] + (x[i] - x[i - 1]) * (y[i][l] + y[i - 1][l]) / 2.0;
                        }
                        i -= 1;
                    }
                }

                //--continue loop over cosine panels
                x[2] = x[1];
                for l in 1..=nlu {
                    y[2][l] = y[1][l];
                }
            }

            //--select next energy point
            let mut ep_nxt = e * powc(sqrt(epnn / e) - sqrt(xc), 2);
            if ep_nxt <= ep * (1.0 + small) {
                ep_nxt = e * powc(sqrt(epnn / e) + sqrt(xc), 2);
            }
            if ep_nxt <= ep * (1.0 + small) {
                ep_nxt = 3.0 * ep / 2.0;
            }
            if ep_nxt > elmax * (1.0 + small) && elmax > ep * (1.0 + small) {
                ep_nxt = elmax;
            }
            epnxt = ep_nxt;
        }
        self.f6cm_tail(term, nl, tol, epnext, epnxt);
    }

    /// f6cm label 330 with no discrete lines: epnxt is not set there.
    fn cm_epnxt_garbage(&self) -> f64 {
        0.0
    }

    /// f6cm 380: remove small, inaccurate moments and select epnext.
    fn f6cm_tail(&self, term: &mut [f64], nl: i32, tol: f64, epnext: &mut f64, epnxt: f64) {
        if nl > 1 {
            let test = tol * term[0].abs();
            for l in 2..=nl as usize {
                if term[l - 1].abs() < test {
                    term[l - 1] = 0.0;
                }
            }
        }
        *epnext = epnxt;
    }

    /// Double differential cross section for secondary energy ep and cm
    /// cosine w from a File 6 subsection in cnow (continuum part only).
    /// Call with ep=0 to initialize for each new e.
    pub(super) fn f6ddx(
        &mut self,
        c: &mut Ctx,
        ep: f64,
        epnext: &mut f64,
        epmax: &mut f64,
        w: f64,
        cnow: &[f64],
        lang: i32,
        lep: i32,
    ) -> f64 {
        let up = 1.00001e0f64;
        let dn = 0.99999e0f64;
        let off = 0.999995e0f64;
        let emax = 1.0e10f64;
        let step = 0.05e0f64;
        let small = 1.0e-10f64;
        let zero = 0.0f64;
        let cn = |i: i32| cnow[(i - 1) as usize];
        let mut p = [0.0f64; MXLG + 1];

        //--initialize
        if ep == zero {
            let s = &mut self.sv;
            s.dx_enow = cn(2);
            let ndnow = nint(cn(3));
            let npnow = nint(cn(6));
            s.dx_ncnow = nint(cn(5)) / npnow;
            s.dx_mnow = 7 + (npnow - 1) * s.dx_ncnow;
            s.dx_inow = 7 + s.dx_ncnow * ndnow;
            s.dx_lnow = s.dx_inow;
            if cn(s.dx_lnow) <= zero {
                s.dx_lnow += s.dx_ncnow;
            }
            *epnext = cn(s.dx_lnow);
            *epmax = cn(s.dx_mnow);
            if *epnext > step * *epmax {
                *epnext = step * *epmax;
            }
            s.dx_nl = s.dx_ncnow - 1;
            if s.dx_nl > MXLG as i32 {
                c.error("f6ddx", "nl>mxlg", " ");
            }
            let s = &mut self.sv;
            s.dx_na = s.dx_nl - 1;
            s.dx_efirst = 0.0;
            s.dx_illdef = 0;
            return 0.0;
        }

        let mut t: f64;
        if ep < self.sv.dx_efirst * (1.0 - small) {
            //--use analytic form below efirst.
            t = 0.0;
        } else {
            //--check for desired interpolation range
            let (inow, mnow, ncnow) = (self.sv.dx_inow, self.sv.dx_mnow, self.sv.dx_ncnow);
            let mut lnow = self.sv.dx_lnow;
            let mut idone = 0;
            while idone == 0 {
                *epnext = cn(lnow);
                if ep < off * *epnext {
                    let eplast = cn(lnow - ncnow);
                    if ep >= off * eplast {
                        idone = 1;
                    } else if lnow <= inow + ncnow {
                        idone = 1;
                    } else {
                        lnow -= ncnow;
                    }
                } else if lnow >= mnow {
                    *epnext = emax;
                    idone = 1;
                } else {
                    lnow += ncnow;
                }
            }
            self.sv.dx_lnow = lnow;
            let (nl, na) = (self.sv.dx_nl, self.sv.dx_na);

            if lang == 1 {
                //--legendre coefficients
                legndr(w, &mut p, na);
                t = 0.0;
                for l in 1..=nl {
                    if l <= ncnow - 1 && ep >= cn(inow) * (1.0 - small) && ep <= cn(mnow) * (1.0 + small) {
                        let lll = lnow + l;
                        let x1 = cn(lnow - ncnow);
                        let x2 = cn(lnow);
                        let y1 = cn(lll - ncnow);
                        let y2 = cn(lll);
                        let mut tt = 0.0f64;
                        if x1 == x2 {
                            tt = y1;
                        } else {
                            c.terp1(x1, y1, x2, y2, ep, &mut tt, lep);
                        }
                        t = t + (2 * l - 1) as f64 * tt * p[(l - 1) as usize] / 2.0;
                    }
                }
                if t < zero {
                    t = 0.0;
                }
            } else if lang == 2 {
                //--kalbach systematics
                let mut s = 0.0f64;
                let mut r = 0.0f64;
                let (mut x1, mut x2) = (self.sv.dx_x1, self.sv.dx_x2);
                if ep >= cn(inow) * (1.0 - small) && ep <= cn(mnow) * (1.0 + small) {
                    x1 = cn(lnow - ncnow);
                    x2 = cn(lnow);
                    if x1 == x2 && lep > 1 {
                        x2 = sigfig(x2, 6, 1);
                        if self.sv.dx_illdef == 0 {
                            c.mess("f6ddx", "vertical segment(s) in distribution", "y(x) is ill defined");
                            self.sv.dx_illdef = 1;
                        }
                    }
                    let y1 = cn(lnow - ncnow + 1);
                    let y2 = cn(lnow + 1);
                    c.terp1(x1, y1, x2, y2, ep, &mut s, lep);
                    let y1 = cn(lnow - ncnow + 2);
                    let y2 = cn(lnow + 2);
                    c.terp1(x1, y1, x2, y2, ep, &mut r, lep);
                }
                self.sv.dx_x1 = x1;
                self.sv.dx_x2 = x2;
                let iza2 = nint(self.zap);
                let mut aa = 0.0f64;
                if na == 2 {
                    let y1 = cn(lnow - ncnow + 3);
                    let y2 = cn(lnow + 3);
                    c.terp1(x1, y1, x2, y2, ep, &mut aa, lep);
                } else {
                    aa = bach(c, self.izap, iza2, self.izat, self.sv.dx_enow, ep);
                }
                t = aa * (cosh(aa * w) + r * sinh(aa * w)) / (2.0 * sinh(aa));
                t = t * s;
                if t < zero {
                    t = 0.0;
                }
            } else {
                //--tabulated angular distribution
                if lang < 11 || lang > 15 {
                    c.error("f6ddx", "illegal lang.", " ");
                }
                let int = lang - 10;
                let ii = lnow - ncnow;
                let jj = lnow;
                let x1 = cn(ii);
                let mut x2 = cn(jj);
                if x1 == x2 && lep > 1 {
                    x2 = sigfig(x2, 6, 1);
                    if self.sv.dx_illdef == 0 {
                        c.mess("f6ddx", "vertical segment(s) in distribution", "y(x) is ill defined");
                        self.sv.dx_illdef = 1;
                    }
                }
                let y1 = cn(ii + 1);
                let y2 = cn(jj + 1);
                let mut s = 0.0f64;
                c.terp1(x1, y1, x2, y2, ep, &mut s, lep);
                let mut tii = 0.0f64;
                let mut ia = 1;
                while ia <= na {
                    if w >= cn(ii + 1 + ia) && w <= cn(ii + 3 + ia) {
                        c.terp1(cn(ii + 1 + ia), cn(ii + 2 + ia), cn(ii + 3 + ia), cn(ii + 4 + ia), w, &mut tii, int);
                    }
                    ia += 2;
                }
                let mut tjj = 0.0f64;
                ia = 1;
                while ia <= na {
                    if w >= cn(jj + 1 + ia) && w <= cn(jj + 3 + ia) {
                        c.terp1(cn(jj + 1 + ia), cn(jj + 2 + ia), cn(jj + 3 + ia), cn(jj + 4 + ia), w, &mut tjj, int);
                    }
                    ia += 2;
                }
                t = 0.0;
                c.terp1(x1, tii, x2, tjj, ep, &mut t, lep);
                t = t * s;
            }
        }

        //--return final value
        if lep > 1 {
        } else if (*epnext - emax).abs() < emax * small {
        } else if ep >= dn * dn * *epnext {
            *epnext = up * *epnext;
        } else {
            *epnext = dn * *epnext;
        }
        t
    }

    /// Double differential cross section for discrete energy i and cm
    /// cosine w from a File 6 subsection in cnow. Call with i=0 to
    /// initialize for each new e.
    pub(super) fn f6dis(
        &mut self,
        c: &mut Ctx,
        i: i32,
        epnext: &mut f64,
        epmax: &mut f64,
        w: f64,
        cnow: &[f64],
        lang: i32,
    ) -> f64 {
        let zero = 0.0f64;
        let cn = |i: i32| cnow[(i - 1) as usize];
        let mut p = [0.0f64; MXLG + 1];

        //--initialize
        if i == 0 {
            let s = &mut self.sv;
            s.ds_enow = cn(2);
            let ndnow = nint(cn(3));
            let npnow = nint(cn(6));
            s.ds_ncnow = nint(cn(5)) / npnow;
            s.ds_inow = 7;
            s.ds_lnow = s.ds_inow;
            s.ds_mnow = s.ds_lnow + (ndnow - 1) * s.ds_ncnow;
            *epnext = cn(s.ds_lnow);
            *epmax = cn(s.ds_mnow);
            s.ds_nl = s.ds_ncnow - 1;
            s.ds_na = s.ds_nl - 1;
            return 0.0;
        }

        //--select desired discrete line
        let ncnow = self.sv.ds_ncnow;
        self.sv.ds_inow = 7 + (i - 1) * ncnow;
        let inow = self.sv.ds_inow;
        *epnext = cn(inow);
        let (nl, na) = (self.sv.ds_nl, self.sv.ds_na);
        let mut t: f64;

        if lang == 1 {
            legndr(w, &mut p, na);
            t = 0.0;
            for l in 1..=nl {
                if l <= ncnow - 1 {
                    t = t + (2 * l - 1) as f64 * p[(l - 1) as usize] * cn(inow + 1 + l) / 2.0;
                }
            }
            t = t * cn(inow + 1);
            if t < zero {
                t = 0.0;
            }
        } else if lang == 2 {
            let s = cn(inow + 1);
            let r = cn(inow + 2);
            let iza2 = nint(self.zap);
            let aa = if na == 2 {
                cn(inow + 3)
            } else {
                bach(c, self.izap, iza2, self.izat, self.sv.ds_enow, *epnext)
            };
            t = aa * (cosh(aa * w) + r * sinh(aa * w)) / (2.0 * sinh(aa));
            t = t * s;
            if t < zero {
                t = 0.0;
            }
        } else if lang >= 11 && lang <= 15 {
            let int = lang - 10;
            t = 0.0;
            let mut ia = 1;
            while ia <= na - 2 {
                if w >= cn(inow + 1 + ia) && w <= cn(inow + 3 + ia) {
                    c.terp1(cn(inow + 1 + ia), cn(inow + 2 + ia), cn(inow + 3 + ia), cn(inow + 4 + ia), w, &mut t, int);
                }
                ia += 2;
            }
            t = t * cn(inow + 1);
        } else {
            c.error("f6dis", "illegal lang.", " ");
        }
        t
    }

    /// Change a law 7 lab angle-energy distribution into a lab Legendre
    /// coefficient representation in law 1 format.
    pub(super) fn ll2lab(&mut self, c: &mut Ctx, inow: i32, jnow: &mut i32, lnow: &mut i32, cc: &mut RV, nl: i32, max: i32) {
        const QP: [f64; 8] = [
            -0.9602898565e0,
            -0.7966664774e0,
            -0.5255324099e0,
            -0.1834346425e0,
            0.1834346425e0,
            0.5255324099e0,
            0.7966664774e0,
            0.9602898565e0,
        ];
        const QW: [f64; 8] = [
            0.1012285362e0,
            0.2223810345e0,
            0.3137066459e0,
            0.3626837834e0,
            0.3626837834e0,
            0.3137066459e0,
            0.2223810345e0,
            0.1012285363e0,
        ];
        let tol = 0.005e0f64;
        let shade = 0.99999e0f64;
        let emax = 1.0e10f64;
        let small = 1.0e-10f64;
        let zero = 0.0f64;
        let mut term = [0.0f64; MXLG + 1];
        let mut p = [0.0f64; MXLG + 1];
        let mut amu = [0.0f64; 52];
        let mut fmu = [0.0f64; 52];

        //--initialize output record
        *jnow = *lnow;
        cc[*jnow] = 0.0;
        cc[*jnow + 1] = cc[inow + 1];
        cc[*jnow + 2] = 0.0;
        cc[*jnow + 3] = (nl - 1) as f64;
        let nll = nl - 1;
        let mut ie = 0;
        *lnow = *jnow + 6;

        //--loop over union grid of secondary energies
        let mut epnext = 0.0f64;
        let mut ip = 2;
        let mut ir = 1;
        let mut idis = 0;
        while epnext < emax * (1.0 - small) {
            let ep = epnext;
            epnext = emax;
            if *lnow + nl + 1 > max {
                c.error("ll2lab", "storage exceeded.", " ");
            }

            //--retrieve the angular distribution
            let nmu = nint(cc[inow + 5]);
            let mut nr = nint(cc[inow + 4]);
            let mut next = inow + 6 + 2 * nr;
            for i in 1..=nmu {
                let ll = next;
                amu[i as usize] = cc[ll + 1];
                let np = nint(cc[ll + 5]);
                nr = nint(cc[ll + 4]);
                next = ll + 6 + 2 * nr + 2 * np;
                let mut epn = 0.0f64;
                let mut fm = 0.0f64;
                c.terpa(&mut fm, ep, &mut epn, &mut idis, &cc.at(ll), &mut ip, &mut ir);
                fmu[i as usize] = fm;
                let test = shade * epn;
                if idis > 0 && ep < test {
                    epn = test;
                }
                if epn < epnext * (1.0 - small) {
                    epnext = epn;
                }
            }

            //--compute legendre coefficients for this distribution.
            for l in 1..=nl as usize {
                term[l] = 0.0;
            }
            for i in 0..8 {
                legndr(QP[i], &mut p, nll);
                let mut j = 1usize;
                while QP[i] >= amu[j] {
                    j += 1;
                }
                for l in 1..=nl as usize {
                    let f = (QP[i] - amu[j - 1]) / (amu[j] - amu[j - 1]);
                    let ff = (1.0 - f) * fmu[j - 1] + f * fmu[j];
                    term[l] = term[l] + ff * p[l - 1] * QW[i];
                }
            }

            //--check normalization
            let nmm = nmu - 1;
            let mut sum = 0.0f64;
            for i in 1..=nmm as usize {
                sum = sum + (amu[i + 1] - amu[i]) * (fmu[i + 1] + fmu[i]) / 2.0;
            }
            let mut fact = 1.0f64;
            if term[1] != zero {
                fact = sum / term[1];
            }
            cc[*lnow] = ep;
            for i in 1..=nl {
                cc[*lnow + i] = term[i as usize] * fact;
            }

            //--is the previous point still needed?
            if ie < 3 {
                ie += 1;
                *lnow = *lnow + nl + 1;
            } else {
                let ln = *lnow;
                let f1 = (cc[ln - nl - 1] - cc[ln - 2 * nl - 2]) / (cc[ln] - cc[ln - 2 * nl - 2]);
                let f2 = 1.0 - f1;
                let mut iflag = 0;
                for i in 1..=nl {
                    let tmid = f1 * cc[ln - 2 * nl - 2 + i] + f2 * cc[ln + i];
                    let test = tol * tmid.abs() + small / 100.0;
                    if (tmid - cc[ln - nl - 1 + i]).abs() > test {
                        iflag = 1;
                    }
                }
                if iflag == 0 {
                    cc[ln - nl - 1] = cc[ln];
                    for i in 1..=nl {
                        cc[ln - nl - 1 + i] = cc[ln + i];
                    }
                } else {
                    ie += 1;
                    *lnow = *lnow + nl + 1;
                }
            }
        }

        //--finished.
        let nw = *lnow - *jnow - 6;
        cc[*jnow + 4] = nw as f64;
        cc[*jnow + 5] = ie as f64;
    }

    /// Legendre coefficients of the double-differential cross section at
    /// lab ep from File 6 law 1 data (continuum part). clo and chi are
    /// 0-based slices of clo(*) and chi(*). Call with ep=0 to initialize.
    pub(super) fn f6lab(
        &mut self,
        c: &mut Ctx,
        ep: f64,
        epnext: &mut f64,
        term: &mut [f64],
        nl: i32,
        _law: i32,
        int: i32,
        lang: i32,
        lep: i32,
        clo: &[f64],
        chi: &[f64],
        e: f64,
    ) {
        const QP: [f64; 8] = [
            -0.9602898565e0,
            -0.7966664774e0,
            -0.5255324099e0,
            -0.1834346425e0,
            0.1834346425e0,
            0.5255324099e0,
            0.7966664774e0,
            0.9602898565e0,
        ];
        const QW: [f64; 8] = [
            0.1012285362e0,
            0.2223810345e0,
            0.3137066459e0,
            0.3626837834e0,
            0.3626837834e0,
            0.3137066459e0,
            0.2223810345e0,
            0.1012285363e0,
        ];
        let up = 1.00001e0f64;
        let dn = 0.99999e0f64;
        let emax = 1.0e10f64;
        let small = 1.0e-10f64;
        let emin = 1.0e-5f64;
        let zero = 0.0f64;
        let cl = |i: i32| clo[(i - 1) as usize];
        let ch = |i: i32| chi[(i - 1) as usize];
        let nlu = nl as usize;
        let mut term1 = [0.0f64; MXLG + 1];
        let mut term2 = [0.0f64; MXLG + 1];
        let mut p = [0.0f64; MXLG + 1];

        //--initialization
        if ep == zero {
            let s = &mut self.sv;
            s.lb_xlo = 1.0;
            s.lb_xhi = 1.0;
            s.lb_xend = 1.0;
            s.lb_elo = cl(2);
            s.lb_ehi = ch(2);
            let ndlo = nint(cl(3));
            let ndhi = nint(ch(3));
            s.lb_intt = int;
            if s.lb_intt > 20 {
                s.lb_intt = int - 20;
            }
            if s.lb_intt > 10 {
                s.lb_intt = int - 10;
            }
            let nplo = nint(cl(6));
            s.lb_nclo = nint(cl(5)) / nplo;
            s.lb_mlo = 7 + (nplo - 1) * s.lb_nclo;
            s.lb_ilo = 7 + s.lb_nclo * ndlo;
            s.lb_llo = s.lb_ilo;
            if cl(s.lb_llo) <= zero {
                s.lb_llo += s.lb_nclo;
            }
            if int > 20 {
                s.lb_xlo = cl(s.lb_mlo);
            }
            let nphi = nint(ch(6));
            s.lb_nchi = nint(ch(5)) / nphi;
            s.lb_mhi = 7 + (nphi - 1) * s.lb_nchi;
            s.lb_ihi = 7 + s.lb_nchi * ndhi;
            s.lb_lhi = s.lb_ihi;
            if ch(s.lb_lhi) <= zero {
                s.lb_lhi += s.lb_nchi;
            }
            if int > 20 {
                s.lb_xhi = ch(s.lb_mhi);
            }
            if int < 11 || int > 15 {
                s.lb_xend = s.lb_xlo + (e - s.lb_elo) * (s.lb_xhi - s.lb_xlo) / (s.lb_ehi - s.lb_elo);
                if cl(s.lb_ilo + 1) > 0.0 || ch(s.lb_ihi + 1) > zero {
                    *epnext = emin;
                } else {
                    *epnext = cl(s.lb_llo) / s.lb_xlo;
                    let epn = ch(s.lb_lhi) / s.lb_xhi;
                    if epn < *epnext * (1.0 - small) {
                        *epnext = epn;
                    }
                    *epnext = *epnext * s.lb_xend;
                }
            } else if cl(s.lb_ilo + 1) > 0.0 || ch(s.lb_ihi + 1) > zero {
                *epnext = emin;
            } else {
                let (elo, ehi, llo, lhi) = (s.lb_elo, s.lb_ehi, s.lb_llo, s.lb_lhi);
                c.terp1(elo, cl(llo), ehi, ch(lhi), e, epnext, 2);
            }
            for l in 0..nlu {
                term[l] = 0.0;
            }
            return;
        }

        //--normal entry
        //--law 1.  find energy range.
        let mut f1 = 0.0f64;
        let mut f2 = 0.0f64;
        let (ilo, mlo, nclo, xlo) = (self.sv.lb_ilo, self.sv.lb_mlo, self.sv.lb_nclo, self.sv.lb_xlo);
        let (ihi, mhi, nchi, xhi) = (self.sv.lb_ihi, self.sv.lb_mhi, self.sv.lb_nchi, self.sv.lb_xhi);
        let (elo, ehi, intt, xend) = (self.sv.lb_elo, self.sv.lb_ehi, self.sv.lb_intt, self.sv.lb_xend);
        let mut llo = self.sv.lb_llo;
        let mut lhi = self.sv.lb_lhi;
        let mut epn = 0.0f64;
        if int < 11 || int > 15 {
            let mut idone = 0;
            while idone == 0 {
                if llo > mlo {
                    *epnext = emax;
                    idone = 1;
                } else {
                    *epnext = cl(llo) * xend / xlo;
                    if ep < *epnext * (1.0 - small) {
                        idone = 1;
                    } else {
                        llo += nclo;
                    }
                }
            }
            idone = 0;
            while idone == 0 {
                if lhi > mhi {
                    epn = emax;
                    idone = 1;
                } else {
                    epn = ch(lhi) * xend / xhi;
                    if ep < epn * (1.0 - small) {
                        idone = 1;
                    } else {
                        lhi += nchi;
                    }
                }
            }
            if epn < *epnext * (1.0 - small) {
                *epnext = epn;
            }
        } else {
            let mut idone = 0;
            while idone == 0 {
                c.terp1(elo, cl(llo), ehi, ch(lhi), e, epnext, 2);
                if ep < *epnext * (1.0 - small) {
                    idone = 1;
                } else if llo >= mlo {
                    *epnext = emax;
                    idone = 1;
                } else {
                    llo += nclo;
                    lhi += nchi;
                }
            }
            if llo != ilo && llo <= mlo {
                let mut eplast = 0.0f64;
                c.terp1(elo, cl(llo - nclo), ehi, ch(lhi - nchi), e, &mut eplast, 2);
                f1 = (*epnext - ep) / (*epnext - eplast);
                f2 = 1.0 - f1;
            }
        }
        self.sv.lb_llo = llo;
        self.sv.lb_lhi = lhi;

        //--interpolate for coefficients.
        if lang == 1 {
            for l in 1..=nl {
                let lu = l as usize;
                if int < 11 || int > 15 {
                    term1[lu] = 0.0;
                    if llo != ilo && llo <= mlo {
                        if l <= nclo - 1 {
                            c.terp1(cl(llo - nclo), cl(llo - nclo + l), cl(llo), cl(llo + l), ep * xlo / xend, &mut term1[lu], lep);
                            term1[lu] = term1[lu] * xlo;
                        }
                    }
                    term2[lu] = 0.0;
                    if lhi != ihi && lhi <= mhi {
                        if l <= nchi - 1 {
                            c.terp1(ch(lhi - nchi), ch(lhi - nchi + l), ch(lhi), ch(lhi + l), ep * xhi / xend, &mut term2[lu], lep);
                            term2[lu] = term2[lu] * xhi;
                        }
                    }
                    c.terp1(elo, term1[lu], ehi, term2[lu], e, &mut term[lu - 1], intt);
                    term[lu - 1] = term[lu - 1] / xend;
                } else {
                    term1[lu] = 0.0;
                    if llo != ihi && llo <= mlo {
                        if l <= nclo - 1 {
                            let epp = f1 * cl(llo) + f2 * cl(llo - nclo);
                            c.terp1(cl(llo - nclo), cl(llo - nclo + l), cl(llo), cl(llo + l), epp, &mut term1[lu], lep);
                        }
                    }
                    term2[lu] = 0.0;
                    if lhi != ihi && lhi <= mhi {
                        if l <= nchi - 1 {
                            let epp = f1 * ch(lhi) + f2 * ch(lhi - nchi);
                            c.terp1(ch(lhi - nchi), ch(lhi - nchi + l), ch(lhi), ch(lhi + l), epp, &mut term2[lu], lep);
                        }
                    }
                    c.terp1(elo, term1[lu], ehi, term2[lu], e, &mut term[lu - 1], intt);
                }
            }
        } else {
            //--compute coefficients from tabulated distribution.
            if lang < 11 || lang > 15 {
                c.error("f6lab", "illegal lang.", " ");
            }
            let int2 = lang - 10;
            let nll = nl - 1;
            let na = (nclo - 2) / 2;
            for l in 1..=nlu {
                term1[l] = 0.0;
            }
            for i in 0..8 {
                legndr(QP[i], &mut p, nll);
                let mut t0 = 0.0f64;
                if llo != ilo && llo <= mlo {
                    let (mut t1, mut t2) = (0.0f64, 0.0f64);
                    let mut j = 1;
                    let mut jj = llo - nclo + 2;
                    while cl(jj) <= QP[i] && j < na {
                        j += 1;
                        jj += 2;
                    }
                    c.terp1(cl(jj - 2), cl(jj - 1), cl(jj), cl(jj + 1), QP[i], &mut t1, int2);
                    j = 1;
                    jj = llo + 2;
                    while cl(jj) <= QP[i] && j < na {
                        j += 1;
                        jj += 2;
                    }
                    c.terp1(cl(jj - 2), cl(jj - 1), cl(jj), cl(jj + 1), QP[i], &mut t2, int2);
                    c.terp1(cl(llo - nclo), t1, cl(llo), t2, ep * xlo / xhi, &mut t0, lep);
                }
                t0 = t0 * QW[i];
                for l in 1..=nlu {
                    term1[l] = term1[l] + t0 * p[l - 1];
                }
            }
            for l in 1..=nlu {
                term2[l] = 0.0;
            }
            for i in 0..8 {
                legndr(QP[i], &mut p, nll);
                let mut t0 = 0.0f64;
                if lhi != ihi && lhi <= mhi {
                    let (mut t1, mut t2) = (0.0f64, 0.0f64);
                    let mut j = 1;
                    let mut jj = lhi - nchi + 2;
                    while ch(jj) <= QP[i] && j < na {
                        j += 1;
                        jj += 2;
                    }
                    c.terp1(ch(jj - 2), ch(jj - 1), ch(jj), ch(jj + 1), QP[i], &mut t1, int2);
                    j = 1;
                    jj = lhi + 2;
                    while ch(jj) <= QP[i] && j < na {
                        j += 1;
                        jj += 2;
                    }
                    c.terp1(ch(jj - 2), ch(jj - 1), ch(jj), ch(jj + 1), QP[i], &mut t2, int2);
                    c.terp1(ch(lhi - nchi), t1, ch(lhi), t2, ep * xhi / xhi, &mut t0, lep);
                }
                t0 = t0 * QW[i];
                for l in 1..=nlu {
                    term2[l] = term2[l] + t0 * p[l - 1];
                }
            }
            for l in 1..=nlu {
                c.terp1(elo, term1[l], ehi, term2[l], e, &mut term[l - 1], intt);
            }
        }

        //--finished with law 1.
        if lep > 1 {
            return;
        }
        if (*epnext - emax).abs() < small * emax {
            return;
        }
        if ep >= dn * dn * *epnext {
            *epnext = up * *epnext;
        } else {
            *epnext = dn * *epnext;
        }
    }

    /// The relativistic discrete gamma or its recoil of mf6/mt102.
    /// ans has leading dimension nl in the Fortran; only column 1 is
    /// used, which is the same memory for any leading dimension.
    pub(super) fn gam102(
        &mut self,
        ans: &mut V2,
        ed: f64,
        enext: &mut f64,
        disc102: f64,
        law: i32,
        nl: i32,
        iglo: &mut i32,
        ng2: &mut i32,
        nq: &mut i32,
    ) {
        let big = 1.0e10f64;
        let awr = self.awr;
        if law == 2 {
            let edis = disc102 + ed * awr / (awr + 1.0);
            for i in 1..=self.ngg {
                if edis >= self.egg[i] && edis < self.egg[i + 1] {
                    *iglo = i;
                }
            }
            *ng2 = 1;
            for i in 1..=nl {
                ans.v[(i - 1) as usize] = 0.0;
            }
            ans.v[0] = 1.0;
            *nq = 0;
            if ed < self.egg[*iglo + 1] {
                *enext = self.egg[*iglo + 1];
            } else if ed < self.egg[*iglo] {
                *enext = self.egg[*iglo];
            } else {
                *enext = big;
            }
        } else if law == 4 {
            let edis = ed / (awr + 1.0);
            for i in 1..=self.ngg {
                if edis >= self.egn[i] && edis < self.egn[i + 1] {
                    *iglo = i;
                }
            }
            *ng2 = 1;
            for i in 1..=nl {
                ans.v[(i - 1) as usize] = 1.0;
            }
            *nq = 0;
            if ed < self.egn[*iglo + 1] {
                *enext = self.egn[*iglo + 1];
            } else if ed < self.egn[*iglo] {
                *enext = self.egn[*iglo];
            } else {
                *enext = big;
            }
        }
    }

    /// n-body phase-space double differential cross section (lab e, cm
    /// secondary energy ep and cosine w). Call with ep=0 to initialize.
    pub(super) fn f6psp(
        &mut self,
        c: &mut Ctx,
        ep: f64,
        epnext: &mut f64,
        epmax: &mut f64,
        _w: f64,
        e: f64,
        cc: &[f64],
    ) -> f64 {
        let c3 = 1.2732e0f64;
        let c4 = 3.2813e0f64;
        let c5 = 5.8205e0f64;
        let thrhaf = 1.5e0f64;
        let sevhaf = 3.5e0f64;
        let step = 0.05e0f64;
        let emax = 1.0e10f64;
        let up = 1.000001e0f64;
        let small = 1.0e-10f64;
        let zero = 0.0f64;

        if ep == zero {
            let apsx = cc[0];
            let npsx = nint(cc[5]);
            let f1 = (apsx - self.awp) / apsx;
            let f2 = self.awr / (self.awr + 1.0);
            self.sv.ps_ex = thrhaf * npsx as f64 - 4.0;
            let mut eimax = f1 * (f2 * e + self.q);
            if eimax <= zero {
                eimax = 1.0;
            }
            self.sv.ps_eimax = eimax;
            if npsx == 3 {
                self.sv.ps_cn = c3 / powc(eimax, 2);
            } else if npsx == 4 {
                self.sv.ps_cn = c4 / pow(eimax, sevhaf);
            } else if npsx == 5 {
                self.sv.ps_cn = c5 / powc(eimax, 5);
            } else {
                c.error("f6psp", "3, 4, or 5 particles only.", " ");
            }
            *epnext = step * eimax;
            *epmax = eimax;
            0.0
        } else {
            let (cn, eimax, ex) = (self.sv.ps_cn, self.sv.ps_eimax, self.sv.ps_ex);
            let mut s = 0.0f64;
            if ep < eimax * (1.0 - small) {
                s = cn * sqrt(ep) * pow(eimax - ep, ex);
            }
            *epnext = ep + step * eimax;
            if *epnext > eimax * (1.0 + small) {
                *epnext = up * eimax;
            }
            if ep >= eimax * (1.0 - small) {
                *epnext = emax;
            }
            s
        }
    }
}

/// Compute the Kalbach a parameter.
pub(super) fn bach(c: &mut Ctx, iza1i: i32, iza2: i32, izat: i32, e: f64, ep: f64) -> f64 {
    let third = 0.333333333e0f64;
    let twoth = 0.666666667e0f64;
    let fourth = 1.33333333e0f64;
    let c1 = 15.68e0f64;
    let c2 = -28.07e0f64;
    let c3 = -18.56e0f64;
    let c4 = 33.22e0f64;
    let c5 = -0.717e0f64;
    let c6 = 1.211e0f64;
    let s2 = 2.22e0f64;
    let s3 = 8.48e0f64;
    let s4 = 7.72e0f64;
    let s5 = 28.3e0f64;
    let brk1 = 130.0e0f64;
    let brk2 = 41.0e0f64;
    let half = 0.5e0f64;
    let b1 = 0.04e0f64;
    let b2 = 1.8e-6f64;
    let b3 = 6.7e-7f64;
    let d1 = 9.3e0f64;
    let tomev = 1.0e-6f64;
    let zero = 0.0f64;
    let emc2 = tomev * AMASSN * AMU * CLIGHT * CLIGHT / EV;

    let mut iza1 = iza1i;
    if iza1i == 0 {
        iza1 = 1;
    }
    let mut iza = izat;
    let map: [(i32, i32); 23] = [
        (6000, 6012),
        (12000, 12024),
        (14000, 14028),
        (16000, 16032),
        (17000, 17035),
        (19000, 19039),
        (20000, 20040),
        (22000, 22048),
        (23000, 23051),
        (24000, 24052),
        (26000, 26056),
        (28000, 28058),
        (29000, 29063),
        (31000, 31069),
        (40000, 40090),
        (42000, 42096),
        (48000, 48112),
        (49000, 49115),
        (50000, 50120),
        (63000, 63151),
        (72000, 72178),
        (74000, 74184),
        (82000, 82208),
    ];
    for (a, b) in map.iter() {
        if iza == *a {
            iza = *b;
        }
    }
    let aa = (iza % 1000) as f64;
    if aa == zero {
        let s = format!("dominant isotope not known for {}", ffmt::i(iza, 8));
        c.error("bach", &s, " ");
    }
    let za = (iza / 1000) as f64;
    let ac = aa + (iza1 % 1000) as f64;
    let zc = za + (iza1 / 1000) as f64;
    let ab = ac - (iza2 % 1000) as f64;
    let zb = zc - (iza2 / 1000) as f64;
    let na = nint(aa - za);
    let nb = nint(ab - zb) as f64;
    let nc = nint(ac - zc) as f64;
    let naf = na as f64;
    let mut sa = c1 * (ac - aa)
        + c2 * (powc(nc - zc, 2) / ac - powc(naf - za, 2) / aa)
        + c3 * (pow(ac, twoth) - pow(aa, twoth))
        + c4 * (powc(nc - zc, 2) / pow(ac, fourth) - powc(naf - za, 2) / pow(aa, fourth))
        + c5 * (powc(zc, 2) / pow(ac, third) - powc(za, 2) / pow(aa, third))
        + c6 * (powc(zc, 2) / ac - powc(za, 2) / aa);
    if iza1 == 1002 {
        sa = sa - s2;
    }
    if iza1 == 1003 {
        sa = sa - s3;
    }
    if iza1 == 2003 {
        sa = sa - s4;
    }
    if iza1 == 2004 {
        sa = sa - s5;
    }
    let mut sb = c1 * (ac - ab)
        + c2 * (powc(nc - zc, 2) / ac - powc(nb - zb, 2) / ab)
        + c3 * (pow(ac, twoth) - pow(ab, twoth))
        + c4 * (powc(nc - zc, 2) / pow(ac, fourth) - powc(nb - zb, 2) / pow(ab, fourth))
        + c5 * (powc(zc, 2) / pow(ac, third) - powc(zb, 2) / pow(ab, third))
        + c6 * (powc(zc, 2) / ac - powc(zb, 2) / ab);
    if iza2 == 1002 {
        sb = sb - s2;
    }
    if iza2 == 1003 {
        sb = sb - s3;
    }
    if iza2 == 2003 {
        sb = sb - s4;
    }
    if iza2 == 2004 {
        sb = sb - s5;
    }
    let ecm = aa * e / ac;
    let ea = ecm * tomev + sa;
    let eb = tomev * ep * ac / ab + sb;
    let mut x1 = eb;
    if ea > brk1 {
        x1 = brk1 * eb / ea;
    }
    let mut x3 = eb;
    if ea > brk2 {
        x3 = brk2 * eb / ea;
    }
    let mut fa = 1.0f64;
    if iza1 == 2004 {
        fa = 0.0;
    }
    let mut fb = 1.0f64;
    if iza2 == 1 {
        fb = half;
    }
    if iza2 == 2004 {
        fb = 2.0;
    }
    let mut bb = b1 * x1 + b2 * powc(x1, 3) + b3 * fa * fb * powc(x3, 4);
    if iza1i == 0 {
        let mut fact = d1 / sqrt(ep * tomev);
        let mut test = 1.0f64;
        if fact < test {
            fact = test;
        }
        test = 4.0;
        if fact > test {
            fact = test;
        }
        bb = bb * sqrt((tomev * e) / (2.0 * emc2)) * fact;
    }
    let _ = (ipow(1, 1), powi(1.0, 1));
    bb
}
