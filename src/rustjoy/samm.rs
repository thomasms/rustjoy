//! Port of NJOY2016 `samm.f90`: the SAMMY R-matrix routines that RECONR
//! uses for LRF=7 (R-matrix limited) resonance parameters.
//!
//! RECONR sets `Want_Partial_Derivs`, `Want_Angular_Dist`, `Want_SAMRML_BW`
//! and `Want_SAMRML_RM` to false, and nothing changes them. All code that is
//! only reachable with one of those flags set (derivatives, Legendre
//! coefficients) is therefore dead in RECONR and is not ported; the places
//! are marked "flag-guarded".
//!
//! Module arrays are flat column-major memory with the declared shape (plus
//! slack), so out-of-shape index arithmetic behaves as in the Fortran.

#![allow(dead_code, unused_variables, unused_assignments, clippy::too_many_arguments, clippy::needless_range_loop)]

use std::ops::{Index, IndexMut};

use crate::endf::Ctx;
use crate::farr::{F, RM, RV};
use crate::ffmt;
use crate::fmath::{atan, atan2, cos, exp, log, log10, nint, pow, powi, sin, sqrt};
use crate::util::physics::{AMASSA, AMASSD, AMASSH, AMASSN, AMASSP, AMASST, AMU, CLIGHT, EULER, EV, FINSTRI, HBAR, PI};
use crate::util::sigfig;

const SLACK: usize = 64;

/// Flat column-major real array with up to 4 dimensions (1-based).
#[derive(Clone, Default, Debug)]
pub struct FA {
    pub v: Vec<f64>,
    pub d: [i32; 4],
}

impl FA {
    pub fn new(d: &[i32]) -> FA {
        let mut dd = [1i32; 4];
        let mut n: usize = 1;
        for (i, x) in d.iter().enumerate() {
            dd[i] = (*x).max(0);
            n *= dd[i] as usize;
        }
        FA { v: vec![0.0; n + SLACK], d: dd }
    }
    #[inline]
    pub fn off(&self, i: i32, j: i32, k: i32, l: i32) -> usize {
        let d = &self.d;
        ((i - 1) + d[0] * ((j - 1) + d[1] * ((k - 1) + d[2] * (l - 1)))) as usize
    }
    /// View starting at an element, as Fortran passes `a(i,j,k,l)` to `a(*)`.
    pub fn view(&mut self, i: i32, j: i32, k: i32, l: i32) -> F<'_> {
        let o = self.off(i, j, k, l);
        F(&mut self.v[o..])
    }
}

macro_rules! idx_impl {
    ($t:ty, $e:ty) => {
        impl Index<(i32, i32)> for $t {
            type Output = $e;
            #[inline]
            fn index(&self, (i, j): (i32, i32)) -> &$e {
                &self.v[self.off(i, j, 1, 1)]
            }
        }
        impl IndexMut<(i32, i32)> for $t {
            #[inline]
            fn index_mut(&mut self, (i, j): (i32, i32)) -> &mut $e {
                let o = self.off(i, j, 1, 1);
                &mut self.v[o]
            }
        }
        impl Index<(i32, i32, i32)> for $t {
            type Output = $e;
            #[inline]
            fn index(&self, (i, j, k): (i32, i32, i32)) -> &$e {
                &self.v[self.off(i, j, k, 1)]
            }
        }
        impl IndexMut<(i32, i32, i32)> for $t {
            #[inline]
            fn index_mut(&mut self, (i, j, k): (i32, i32, i32)) -> &mut $e {
                let o = self.off(i, j, k, 1);
                &mut self.v[o]
            }
        }
        impl Index<(i32, i32, i32, i32)> for $t {
            type Output = $e;
            #[inline]
            fn index(&self, (i, j, k, l): (i32, i32, i32, i32)) -> &$e {
                &self.v[self.off(i, j, k, l)]
            }
        }
        impl IndexMut<(i32, i32, i32, i32)> for $t {
            #[inline]
            fn index_mut(&mut self, (i, j, k, l): (i32, i32, i32, i32)) -> &mut $e {
                let o = self.off(i, j, k, l);
                &mut self.v[o]
            }
        }
    };
}

/// Flat column-major integer array.
#[derive(Clone, Default, Debug)]
pub struct FI {
    pub v: Vec<i32>,
    pub d: [i32; 4],
}

impl FI {
    pub fn new(d: &[i32]) -> FI {
        let mut dd = [1i32; 4];
        let mut n: usize = 1;
        for (i, x) in d.iter().enumerate() {
            dd[i] = (*x).max(0);
            n *= dd[i] as usize;
        }
        FI { v: vec![0; n + SLACK], d: dd }
    }
    #[inline]
    pub fn off(&self, i: i32, j: i32, k: i32, l: i32) -> usize {
        let d = &self.d;
        ((i - 1) + d[0] * ((j - 1) + d[1] * ((k - 1) + d[2] * (l - 1)))) as usize
    }
}

idx_impl!(FA, f64);
idx_impl!(FI, i32);

#[inline]
fn sq(x: f64) -> f64 {
    x * x
}

#[derive(Default, Clone)]
pub struct Samm {
    pub want_partial_derivs: bool,
    pub want_angular_dist: bool,
    pub want_samrml_bw: bool,
    pub want_samrml_rm: bool,

    ngroup: i32,
    mchan: i32,
    nres: i32,
    npp: i32,
    ntriag: i32,
    npar: i32,
    lllmax: i32,
    nback: i32,
    nrext: i32,
    nparr: i32,
    nerm: i32,
    ngroupm: [i32; 13],
    nppm: [i32; 13],
    nresm: [i32; 13],
    awr: f64,
    su: f64,

    nchan: FI,
    ema: FA,
    emb: FA,
    spina: FA,
    spinb: FA,
    qqq: FA,
    kza: FI,
    kzb: FI,
    ishift: FI,
    lpent: FI,
    mt: FI,
    pa: FA,
    pb: FA,
    sspin: FA,
    parity: FA,
    nresg: FI,
    ipp: FI,
    lspin: FI,
    backgr: FI,
    backgrdata: FA,
    chspin: FA,
    bound: FA,
    rdeff: FA,
    rdtru: FA,
    eres: FA,
    gamgam: FA,
    gamma: FA,
    nent: FI,
    next: FI,
    goj: FA,
    zke: FA,
    zkfe: FA,
    zkte: FA,
    zeta: FA,
    echan: FA,
    betapr: FA,
    gbetpr: FA,
    beta: FA,
    uuuu: FA,
    duuu: FA,
    iduu: FI,
    crss: FA,
    sigmas: FA,
    alphar: FA,
    alphai: FA,
    difen: FA,
    xden: FA,
    sinsqr: FA,
    sin2ph: FA,
    dphi: FA,
    dpdr: FA,
    dsdr: FA,
    sinphi: FA,
    cosphi: FA,
    rootp: FA,
    elinvr: FA,
    elinvi: FA,
    psmall: FA,
    xqr: FA,
    xqi: FA,
    xxxxr: FA,
    xxxxi: FA,
    qr: FA,
    qi: FA,
    rmat: FA,
    ymat: FA,
    yinv: FA,
    cscs: FA,
    ddddd: FA,
    par: FA,

    // SAVEd locals
    findsp_ifneg: i32,
}

impl Samm {
    // ============================================================ cssammy
    pub fn cssammy(
        &mut self,
        c: &mut Ctx,
        e: f64,
        sigp: &mut RV,
        siga: &mut RM,
        sigd: &mut RM,
        mmtres: &[i32; 11],
        nmtres: i32,
        ncoef: i32,
        nresp: i32,
        ier: i32,
        _res: &RV,
    ) {
        self.su = e;
        self.abpart(ier);
        self.crosss(c, ier);
        sigp[1] = self.sigmas[(1, ier)] + self.sigmas[(2, ier)];
        sigp[2] = self.sigmas[(1, ier)];
        sigp[3] = 0.0;
        sigp[4] = self.sigmas[(2, ier)];
        if nmtres > 2 {
            for i in 3..=nmtres {
                sigp[4] = sigp[4] - self.sigmas[(i, ier)];
                if mmtres[i as usize] == 18 {
                    sigp[3] = self.sigmas[(i, ier)];
                } else {
                    sigp[2 + i] = self.sigmas[(i, ier)];
                }
            }
        }
        // Want_Angular_Dist / Want_Partial_Derivs: flag-guarded (false)
    }

    // ============================================================ s2sammy
    pub fn s2sammy(&mut self, c: &mut Ctx, nin: i32, res: &mut RV, maxres: i32, mmtres: &mut [i32; 11], nmtres: &mut i32) {
        let zero = 0.0f64;
        let half = 0.5e0f64;
        let (mut nb, mut nw) = (0i32, 0i32);
        self.nchan = FI::new(&[100, 12]);
        *nmtres = 0;
        self.nrext = 0;

        let mut jnow = 1;
        c.contio(nin, 0, 0, &mut res.at(jnow), &mut nb, &mut nw);
        if c.mfh != 2 {
            c.io.so("\n mat has no file 2.");
            c.tosend(nin, 0, 0, &mut res.all());
            return;
        }
        self.awr = c.c2h;
        let nis = c.n1h;
        self.nres = 0;
        self.ngroup = 0;
        self.npp = 0;
        self.mchan = 0;
        self.npar = 0;
        self.nback = 0;
        let mut kres;
        let mut kpp;
        let mut kchan;
        let mut kback;
        let mut jj: i32 = 0;

        if nis > 1 {
            c.mess("s2sammy", "multiple isotopes do not work with sammy ", " method--reverting to normal RM processing");
            *nmtres = 0;
        } else {
            c.io.so(" resonance range information\n ----------------------------------\n ier     energy-range     lru lrf  method");
            for _iso in 1..=nis {
                c.contio(nin, 0, 0, &mut res.at(jnow), &mut nb, &mut nw);
                jnow += nw;
                if jnow > maxres {
                    c.error("s2sammy", "res storage exceeded", " ");
                }
                let lfw = c.l2h;
                let ner = c.n1h;
                self.nerm = ner;
                for ier in 1..=ner {
                    c.contio(nin, 0, 0, &mut res.at(jnow), &mut nb, &mut nw);
                    let e1 = c.c1h;
                    let e2 = c.c2h;
                    let lru = c.l1h;
                    let lrf = c.l2h;
                    if lrf == 4 {
                        c.mess("s2sammy", "adler-adler format does not work with sammy", " method, reverting to normal a-a processing");
                        *nmtres = 0;
                        return;
                    }
                    let nro = c.n1h;
                    if lru == 1 {
                        res[jnow + 4] = nro as f64;
                    } else {
                        res[jnow + 4] = lfw as f64;
                    }
                    let naps = c.n2h;
                    res[jnow + 5] = naps as f64;
                    let mut mode = lrf + 10 * (lru - 1);
                    if lru == 0 {
                        mode = 0;
                    }
                    if nro != 0 {
                        if lru == 1 && lrf != 2 {
                            c.error("s2sammy", "energy-dep scattering length", "only works for mlbw");
                        }
                        c.tab1io(nin, 0, 0, &mut res.at(jnow), &mut nb, &mut nw);
                        jnow += nw;
                        if jnow > maxres {
                            c.error("s2sammy", "res storage exceeded", " ");
                        }
                    }
                    if mode <= 3 {
                        // scanning SLBW/MLBW (mode<=2) and RM (mode=3)
                        c.contio(nin, 0, 0, &mut res.at(jnow), &mut nb, &mut nw);
                        jnow += nw;
                        if jnow > maxres {
                            c.error("s2sammy", "res storage exceeded", " ");
                        }
                        let spin = c.c1h;
                        let nls = c.n1h;
                        let mut parity = 1.0;
                        if spin < zero {
                            parity = -1.0;
                        }
                        let mut s = [0.0f64; 3];
                        s[1] = spin.abs() - half;
                        s[2] = spin.abs() + half;
                        let mut is = 2;
                        if spin == zero {
                            s[1] = half;
                            is = 1;
                        }
                        let mut ng = 0;
                        let mut parl = -parity;
                        for ll in 1..=nls {
                            let el = (ll - 1) as f64;
                            parl = -parl;
                            for iis in 1..=is {
                                let mut capj = (s[iis].abs() - el).abs();
                                let capjmx = s[iis].abs() + el;
                                loop {
                                    ng += 1;
                                    capj = capj + 1.0;
                                    if capj > capjmx {
                                        break;
                                    }
                                }
                            }
                        }
                        self.ngroupm[ier as usize] = ng;
                        if ng > self.ngroup {
                            self.ngroup = ng;
                        }
                        kchan = 1;
                        kres = 0;
                        kpp = 2;
                        let mut kpar = 0;
                        let mut kki = 0;
                        let mut kkf = 0;
                        for _i in 1..=nls {
                            c.listio(nin, 0, 0, &mut res.at(jnow), &mut nb, &mut nw);
                            jj = jnow + nw;
                            while nb != 0 {
                                c.moreio(nin, 0, 0, &mut res.at(jj), &mut nb, &mut nw);
                                jj += nw;
                                if jj > maxres {
                                    c.error("s2sammy", "res storage exceeded", " ");
                                }
                            }
                            let lrx = nint(res[jnow + 3]);
                            let nrs = nint(res[jnow + 5]);
                            kres += nrs;
                            jj = jnow;
                            if mode <= 2 {
                                let mut kf = 0;
                                let mut ki = 0;
                                for _j in 1..=nrs {
                                    jj += 6;
                                    let gamf = res[jj + 5];
                                    let mut gamx = 0.0;
                                    if lrx > 0 {
                                        gamx = res[jj + 2] - res[jj + 3] - res[jj + 4] - res[jj + 5];
                                    }
                                    if gamx < 1.0e7 {
                                        gamx = 0.0;
                                    }
                                    if gamf != zero {
                                        kf = 1;
                                    }
                                    if gamx != zero {
                                        ki = 1;
                                    }
                                }
                                if kf == 1 {
                                    kkf = 1;
                                }
                                if ki == 1 {
                                    kki = 1;
                                }
                                kchan = 1 + kf + ki;
                                if kchan > 1 {
                                    kpp = 3;
                                }
                            } else {
                                let mut kk = 1;
                                for _j in 1..=nrs {
                                    jj += 6;
                                    let gamf = res[jj + 4];
                                    let gamf2 = res[jj + 5];
                                    if gamf != zero && kk < 2 {
                                        kk = 2;
                                    }
                                    if gamf2 != zero && kk < 3 {
                                        kk = 3;
                                    }
                                }
                                if kk > 1 {
                                    kchan = kk;
                                    kpp = 3;
                                }
                            }
                            kpar += self.nres * (kpp + 2);
                        }
                        if kpar > self.npar {
                            self.npar = kpar;
                        }
                        if kchan > self.mchan {
                            self.mchan = kchan;
                        }
                        self.nresm[ier as usize] = kres;
                        if kres > self.nres {
                            self.nres = kres;
                        }
                        self.nppm[ier as usize] = kpp;
                        if kpp > self.npp {
                            self.npp = kpp;
                        }
                        mmtres[1] = 2;
                        mmtres[2] = 102;
                        *nmtres = 2;
                        if mode <= 2 {
                            if kkf > 0 && kki == 0 {
                                mmtres[3] = 18;
                                *nmtres = 3;
                            }
                            if kkf == 0 && kki > 0 {
                                mmtres[3] = 51;
                                *nmtres = 3;
                            }
                            if kkf > 0 && kki > 0 {
                                mmtres[3] = 18;
                                mmtres[4] = 51;
                                *nmtres = 4;
                            }
                            if !self.want_samrml_bw {
                                *nmtres = 0;
                            }
                        } else {
                            if self.npp > 2 {
                                mmtres[3] = 18;
                                *nmtres = 3;
                            }
                            if !self.want_samrml_rm {
                                *nmtres = 0;
                            }
                        }
                    } else if mode == 7 {
                        if nis > 1 {
                            c.error("s2sammy", "multiple isotopes do not work with sammy method", " ");
                        }
                        if nro != 0 {
                            if mode != 2 {
                                c.error("s2sammy", "energy-dep scattering length", "only works for mlbw");
                            }
                        }
                        c.contio(nin, 0, 0, &mut res.at(jnow), &mut nb, &mut nw);
                        jnow += nw;
                        if jnow > maxres {
                            c.error("s2sammy", "res storage exceeded", " ");
                        }
                        let nls = c.n1h;
                        self.ngroupm[ier as usize] = nls;
                        if nls > self.ngroup {
                            self.ngroup = nls;
                        }
                        c.listio(nin, 0, 0, &mut res.at(jnow), &mut nb, &mut nw);
                        kpp = c.l1h;
                        for ipp in 1..=kpp {
                            mmtres[ipp as usize] = nint(res[jnow + 15 + 12 * (ipp - 1)]);
                        }
                        self.nppm[ier as usize] = kpp;
                        if kpp > self.npp {
                            self.npp = kpp;
                        }
                        mmtres[1] = 2;
                        mmtres[2] = 102;
                        *nmtres = self.npp;
                        jnow += nw;
                        kchan = 0;
                        kback = 0;
                        kres = 0;
                        let mut kpar = 0;
                        for igroup in 1..=self.ngroup {
                            c.listio(nin, 0, 0, &mut res.at(jnow), &mut nb, &mut nw);
                            jnow += nw;
                            while nb != 0 {
                                // note: the Fortran uses jj here before it is set
                                c.moreio(nin, 0, 0, &mut res.at(jj.max(1)), &mut nb, &mut nw);
                                jj += nw;
                                if jj > maxres {
                                    c.error("s2sammy", "res storage exceeded", " ");
                                }
                            }
                            let kbk = c.l1h;
                            let nch = c.n2h;
                            let ich = nch - 1;
                            self.nchan[(igroup, ier)] = ich;
                            if ich > kchan {
                                kchan = ich;
                            }
                            c.listio(nin, 0, 0, &mut res.at(jnow), &mut nb, &mut nw);
                            jj = jnow + nw;
                            while nb != 0 {
                                c.moreio(nin, 0, 0, &mut res.at(jj), &mut nb, &mut nw);
                                jj += nw;
                                if jj > maxres {
                                    c.error("s2sammy", "res storage exceeded", " ");
                                }
                            }
                            let nrs = nint(res[jnow + 3]);
                            kpar += nrs * (ich + 2);
                            kres += nrs;
                            if kbk > 0 {
                                for _i in 1..=kbk {
                                    c.contio(nin, 0, 0, &mut res.at(jnow), &mut nb, &mut nw);
                                    let lbk = c.l2h;
                                    if lbk == 1 {
                                        c.tab1io(nin, 0, 0, &mut res.at(jj), &mut nb, &mut nw);
                                        jj = jnow + nw;
                                        while nb != 0 {
                                            c.moreio(nin, 0, 0, &mut res.at(jj), &mut nb, &mut nw);
                                            jj += nw;
                                            if jj > maxres {
                                                c.error("s2sammy", "res storage exceeded", " ");
                                            }
                                        }
                                        c.tab1io(nin, 0, 0, &mut res.at(jj), &mut nb, &mut nw);
                                        jj = jnow + nw;
                                        while nb != 0 {
                                            c.moreio(nin, 0, 0, &mut res.at(jj), &mut nb, &mut nw);
                                            jj += nw;
                                            if jj > maxres {
                                                c.error("s2sammy", "res storage exceeded", " ");
                                            }
                                        }
                                        if jj - jnow + 1 > kback {
                                            kback = jj - jnow + 1;
                                        }
                                    } else if lbk == 2 || lbk == 3 {
                                        c.listio(nin, 0, 0, &mut res.at(jj), &mut nb, &mut nw);
                                        jj = jnow + nw;
                                        while nb != 0 {
                                            c.moreio(nin, 0, 0, &mut res.at(jj), &mut nb, &mut nw);
                                            jj += nw;
                                            if jj > maxres {
                                                c.error("s2sammy", "res storage exceeded", " ");
                                            }
                                        }
                                        if lbk == 2 {
                                            if 7 > kback {
                                                kback = 7;
                                            }
                                        } else if 5 > kback {
                                            kback = 5;
                                        }
                                    }
                                }
                            }
                        }
                        if kpar > self.npar {
                            self.npar = kpar;
                        }
                        if kchan > self.mchan {
                            self.mchan = kchan;
                        }
                        if kback > self.nback {
                            self.nback = kback;
                        }
                        self.nresm[ier as usize] = kres;
                        if kres > self.nres {
                            self.nres = kres;
                        }
                    }
                    let method = if *nmtres == 0 { "   normal" } else { "   sammy" };
                    c.io.so(&format!(
                        "{} {}{}{}{}{}",
                        ffmt::i(ier, 3),
                        ffmt::e1p(e1, 10, 3),
                        ffmt::e1p(e2, 10, 3),
                        ffmt::i(lru, 4),
                        ffmt::i(lrf, 4),
                        method
                    ));
                }
            }
        }
        // label 10
        if *nmtres > 0 {
            self.lllmax = 10;
            self.allo();
        }
    }

    // ============================================================ rdsammy
    pub fn rdsammy(
        &mut self,
        c: &mut Ctx,
        nin: i32,
        ier: i32,
        jnow: &mut i32,
        nro: i32,
        naps: i32,
        mode: i32,
        el: f64,
        eh: f64,
        enode: &mut RV,
        nodes: &mut i32,
        nodmax: i32,
        spin: &mut f64,
        ascat: &mut f64,
        res: &mut RV,
        maxres: i32,
    ) {
        let zero = 0.0f64;
        let one = 1.0f64;
        let three = 3.0f64;
        let half = 0.5e0f64;
        let zz8 = 0.08e0f64;
        let z123 = 0.123e0f64;
        let (mut nb, mut nw) = (0i32, 0i32);
        if nro != 0 {
            c.error("rdsammy", "energy dependent scattering radius", "not currently allowed for sammy rml or rm methods");
        }
        let mut idum = vec![0i32; self.nres.max(0) as usize + 2];
        let mut igamma = 0;
        c.contio(nin, 0, 0, &mut res.at(*jnow), &mut nb, &mut nw);
        *jnow += nw;
        if *jnow > maxres {
            c.error("rdsammy", "res storage exceeded", " ");
        }
        let nls = c.n1h;
        *spin = c.c1h;
        *ascat = c.c2h;

        let push_nodes = |c: &mut Ctx, enode: &mut RV, nodes: &mut i32, e0: f64, hw: f64| {
            let mut ndig = 5;
            if e0 > zero {
                ndig = 2 + nint(log10(e0 / (hw / 10.0)));
            }
            if ndig < 5 {
                ndig = 5;
            }
            if ndig > 9 {
                ndig = 9;
            }
            if e0 > el && e0 < eh {
                *nodes += 1;
                enode[*nodes] = sigfig(e0, ndig, 0);
            }
            if (e0 + hw) > el && (e0 + hw) < eh {
                if *nodes >= nodmax {
                    c.error("rdsammy", "storage in enode exceeded.", " ");
                }
                let ehalf = e0 + hw;
                *nodes += 1;
                enode[*nodes] = sigfig(ehalf, ndig, 0);
            }
            if (e0 - hw) > el && (e0 - hw) < eh {
                if *nodes >= nodmax {
                    c.error("rdsammy", "storage in enode exceeded.", " ");
                }
                let ehalf = e0 - hw;
                *nodes += 1;
                enode[*nodes] = sigfig(ehalf, ndig, 0);
            }
        };

        if mode <= 3 {
            for i in 1..=self.nppm[ier as usize] {
                self.ishift[(i, ier)] = 0;
                self.lpent[(i, ier)] = 0;
                self.kza[(i, ier)] = 0;
                self.kzb[(i, ier)] = 0;
                self.mt[(i, ier)] = 0;
                self.spina[(i, ier)] = 0.0;
                self.spinb[(i, ier)] = 0.0;
                self.ema[(i, ier)] = 0.0;
                self.emb[(i, ier)] = 0.0;
                self.qqq[(i, ier)] = 0.0;
            }
            self.spinb[(2, ier)] = *spin;
            self.lpent[(2, ier)] = 1;
            self.mt[(1, ier)] = 102;
            self.mt[(2, ier)] = 2;
            if self.nppm[ier as usize] == 3 {
                self.mt[(3, ier)] = 18;
            }
            self.ppdefs(ier);
            for i in 1..=self.mchan {
                for j in 1..=self.ngroupm[ier as usize] {
                    self.lspin[(i, j, ier)] = 0;
                    self.chspin[(i, j, ier)] = 0.0;
                    self.backgr[(i, j, ier)] = 0;
                    // backgrdata(i,igroup,ier,:)=0 with an undefined igroup;
                    // with nback=0 (no LRF=7 background) this does nothing
                }
            }
            let mut pari = 1.0;
            if *spin < zero {
                pari = -1.0;
            }
            let mut s = [0.0f64; 3];
            s[1] = spin.abs() - half;
            s[2] = spin.abs() + half;
            let mut is = 2;
            if *spin == zero {
                s[1] = half;
                is = 1;
            }
            let s1 = s[1];
            let s2 = s[2];
            let mut ng = 0;
            let mut parl = -pari;
            for ll in 1..=nls {
                let ell = (ll - 1) as f64;
                let l = ll - 1;
                parl = -parl;
                for iis in 1..=is {
                    let mut capj = (s[iis].abs() - ell).abs();
                    let capjmx = s[iis].abs() + ell;
                    loop {
                        let cc = capj * parl;
                        ng += 1;
                        self.sspin[(ng, ier)] = cc;
                        self.parity[(ng, ier)] = parl;
                        self.ipp[(1, ng, ier)] = 2;
                        self.lspin[(1, ng, ier)] = l;
                        self.chspin[(1, ng, ier)] = s[iis];
                        capj = capj + 1.0;
                        if capj > capjmx {
                            break;
                        }
                    }
                }
            }
            if mode <= 2 {
                self.ngroupm[ier as usize] = ng;
            } else {
                self.ngroup = ng;
            }
            for i in 1..=self.ngroupm[ier as usize] {
                self.nchan[(i, ier)] = 1;
            }
            let mut ires = 0;
            for i in 1..=nls {
                c.listio(nin, 0, 0, &mut res.at(*jnow), &mut nb, &mut nw);
                let mut jj = *jnow + nw;
                while nb != 0 {
                    c.moreio(nin, 0, 0, &mut res.at(jj), &mut nb, &mut nw);
                    jj += nw;
                    if jj > maxres {
                        c.error("rdsammy", "res storage exceeded", " ");
                    }
                }
                let mut lrx = 0;
                let awri;
                let aptru;
                let apeff;
                if mode <= 2 {
                    let _qx = res[*jnow + 1];
                    lrx = nint(res[*jnow + 3]);
                    awri = res[*jnow];
                    self.emb[(2, ier)] = awri;
                    let mut ap = z123 * pow(awri * AMASSN, one / three) + zz8;
                    apeff = *ascat;
                    if naps == 1 {
                        ap = *ascat;
                    }
                    aptru = ap;
                } else {
                    awri = res[*jnow];
                    self.emb[(2, ier)] = awri;
                    let apl = res[*jnow + 1];
                    if apl == zero {
                        aptru = z123 * pow(awri * AMASSN, one / three) + zz8;
                        apeff = *ascat;
                    } else {
                        apeff = apl;
                        if naps == 1 {
                            aptru = apl;
                        } else {
                            aptru = z123 * pow(awri * AMASSN, one / three) + zz8;
                        }
                    }
                }
                let nrs = nint(res[*jnow + 5]);
                let llll = nint(res[*jnow + 2]);
                let mut igx = [0i32; 21];
                let mut igxm = 0;
                let ngl = if mode <= 2 { self.ngroupm[ier as usize] } else { self.ngroup };
                for j in 1..=ngl {
                    if llll == self.lspin[(1, j, ier)] {
                        igxm += 1;
                        igx[igxm as usize] = j;
                    }
                }
                for ig in 1..=igxm {
                    if igx[ig as usize] != 0 {
                        self.rdeff[(1, igx[ig as usize], ier)] = apeff;
                        self.rdtru[(1, igx[ig as usize], ier)] = aptru;
                    }
                }
                jj = *jnow + 6;
                for _j in 1..=nrs {
                    ires += 1;
                    let mut hw = res[jj + 2] / 2.0;
                    hw = hw + (res[jj + 3] + res[jj + 4].abs() + res[jj + 5].abs()) / 2.0;
                    push_nodes(c, enode, nodes, res[jj], hw);
                    let spinjj = res[jj + 1];
                    let ig = self.findsp(c, llll, spinjj, s1, s2, ier);
                    idum[ires as usize] = ig;
                    self.eres[(ires, ier)] = res[jj];
                    if mode <= 2 {
                        self.gamma[(1, ires, ier)] = res[jj + 3];
                        self.gamgam[(ires, ier)] = res[jj + 4];
                        let mut gamf = res[jj + 5];
                        gamf = gamf * crate::fmath::ipow(-1, ires) as f64;
                        let mut gamx = 0.0;
                        if lrx > 0 {
                            gamx = res[jj + 2] - res[jj + 3] - res[jj + 4] - res[jj + 5];
                        }
                        if gamx < 1.0e-7 {
                            gamx = 0.0;
                        }
                        if gamf != zero && gamx == zero {
                            self.gamma[(2, ires, ier)] = gamf;
                            if self.nchan[(ig, ier)] < 2 {
                                self.nchan[(ig, ier)] = 2;
                                self.ipp[(2, ig, ier)] = 3;
                            }
                        }
                        if gamf == zero && gamx != 0.0 {
                            self.gamma[(2, ires, ier)] = gamx;
                            if self.nchan[(ig, ier)] < 2 {
                                self.nchan[(ig, ier)] = 2;
                                self.ipp[(3, ig, ier)] = 3;
                            }
                        }
                        if gamf != zero && gamx != 0.0 {
                            self.gamma[(2, ires, ier)] = gamf;
                            self.gamma[(3, ires, ier)] = gamx;
                            if self.nchan[(i, ier)] < 3 {
                                self.nchan[(ig, ier)] = 3;
                                self.ipp[(3, ig, ier)] = 3;
                            }
                        }
                    } else {
                        // the mode=3 (RM) branch is not reachable from
                        // RECONR (Want_SAMRML_RM is false)
                        c.error("rdsammy", "sammy RM branch is not ported (unreachable in RECONR)", " ");
                    }
                    jj += 6;
                }
            }
            self.rearrange(c, &mut idum, ier);
        } else if mode == 7 {
            self.ngroup = nls;
            let ifg = c.l1h;
            let krm = c.l2h;
            if ifg != 0 {
                c.error("rdsammy", "reduced resonance widths are not supported (IFG=1)", " ");
            }
            if krm != 3 {
                c.error("rdsammy", "LRF=7 currently only supports Reich-Moore (KRM=3)", " ");
            }
            c.listio(nin, 0, 0, &mut res.at(*jnow), &mut nb, &mut nw);
            let mut jj = *jnow + nw;
            while nb != 0 {
                c.moreio(nin, 0, 0, &mut res.at(jj), &mut nb, &mut nw);
                jj += nw;
                if jj > maxres {
                    c.error("rdsammy", "res storage exceeded", " ");
                }
            }
            self.npp = nint(res[*jnow + 2]);
            jj = *jnow + 6;
            for i in 1..=self.npp {
                self.ema[(i, ier)] = res[jj];
                self.emb[(i, ier)] = res[jj + 1];
                self.kza[(i, ier)] = nint(res[jj + 2]);
                self.kzb[(i, ier)] = nint(res[jj + 3]);
                self.spina[(i, ier)] = res[jj + 4];
                self.spinb[(i, ier)] = res[jj + 5];
                self.qqq[(i, ier)] = res[jj + 6];
                self.lpent[(i, ier)] = nint(res[jj + 7]);
                self.ishift[(i, ier)] = nint(res[jj + 8]);
                self.mt[(i, ier)] = nint(res[jj + 9]);
                if self.mt[(i, ier)] == 2 {
                    *spin = self.spinb[(i, ier)];
                }
                self.pa[(i, ier)] = res[jj + 10];
                self.pb[(i, ier)] = res[jj + 11];
                jj += 12;
            }
            self.ppdefs(ier);
            let mut kres = 0;
            for igroup in 1..=self.ngroup {
                c.listio(nin, 0, 0, &mut res.at(*jnow), &mut nb, &mut nw);
                jj = *jnow + nw;
                while nb != 0 {
                    c.moreio(nin, 0, 0, &mut res.at(jj), &mut nb, &mut nw);
                    jj += nw;
                    if jj > maxres {
                        c.error("rdsammy", "res storage exceeded", " ");
                    }
                }
                self.sspin[(igroup, ier)] = c.c1h;
                self.parity[(igroup, ier)] = c.c2h;
                let kbk = c.l1h;
                let ichp1 = c.n2h;
                let ichan = ichp1 - 1;
                let mut ix = 0;
                let mut i = 0;
                jj = *jnow + 6;
                for ich in 1..=ichp1 {
                    let ippx = nint(res[jj]);
                    if ippx == 1 {
                        igamma = ich;
                    } else {
                        i += 1;
                        self.ipp[(i, igroup, ier)] = ippx;
                        if self.ipp[(i, igroup, ier)] == 2 {
                            ix += 1;
                        }
                        self.lspin[(i, igroup, ier)] = nint(res[jj + 1]);
                        self.chspin[(i, igroup, ier)] = res[jj + 2];
                        self.bound[(i, igroup, ier)] = res[jj + 3];
                        self.rdeff[(i, igroup, ier)] = res[jj + 4];
                        let ip = self.ipp[(i, igroup, ier)];
                        if self.mt[(ip, ier)] == 2 {
                            *ascat = res[jj + 4];
                        }
                        self.rdtru[(i, igroup, ier)] = res[jj + 5];
                        self.backgr[(i, igroup, ier)] = 0;
                        for k in 1..=self.nback {
                            self.backgrdata[(i, igroup, ier, k)] = 0.0;
                        }
                    }
                    jj += 6;
                }
                c.listio(nin, 0, 0, &mut res.at(*jnow), &mut nb, &mut nw);
                jj = *jnow + nw;
                while nb != 0 {
                    c.moreio(nin, 0, 0, &mut res.at(jj), &mut nb, &mut nw);
                    jj += nw;
                    if jj > maxres {
                        c.error("rdsammy", "res storage exceeded", " ");
                    }
                }
                self.nresg[(igroup, ier)] = c.l2h;
                if self.nresg[(igroup, ier)] > 0 {
                    jj = *jnow + 6;
                    let mut nx = 6;
                    if 2 + ichan > 6 {
                        nx = 12;
                    }
                    for _ires in 1..=self.nresg[(igroup, ier)] {
                        kres += 1;
                        let mut hw = res[jj + 1] / 2.0;
                        for j in 1..=ichan {
                            hw = hw + res[jj + 1 + j].abs() / 2.0;
                        }
                        push_nodes(c, enode, nodes, res[jj], hw);
                        self.eres[(kres, ier)] = res[jj];
                        self.gamgam[(kres, ier)] = res[jj + 1];
                        for j in 1..=ichan {
                            self.gamma[(j, kres, ier)] = res[jj + 1 + j];
                        }
                        if igamma != 1 {
                            let x = self.gamma[(igamma + 1, kres, ier)];
                            if igamma > 1 {
                                let mut i2 = igamma;
                                while i2 >= 2 {
                                    self.gamma[(i2, kres, ier)] = self.gamma[(i2 - 1, kres, ier)];
                                    i2 -= 1;
                                }
                            }
                            self.gamma[(1, kres, ier)] = self.gamgam[(kres, ier)];
                            self.gamgam[(kres, ier)] = x;
                        }
                        jj += nx;
                    }
                }
                if kbk > 0 {
                    self.nrext = 1;
                    if krm != 3 {
                        c.error(
                            "rdsammy",
                            "Background rmatrix elements in LRF=7 are only allowed for Reich-Moore (KRM=3)",
                            " ",
                        );
                    }
                    for _l in 1..=kbk {
                        c.contio(nin, 0, 0, &mut res.at(*jnow), &mut nb, &mut nw);
                        let lch = c.l1h;
                        let lbk = c.l2h;
                        if lbk != 0 && lch == 1 {
                            let strng = format!(
                                "channel {} in spin group {} is the eliminated capture channel",
                                ffmt::i(lch, 3),
                                ffmt::i(igroup, 3)
                            );
                            c.error("rdsammy", &strng, "Background rmatrix data is not allowed for this channel");
                        }
                        if lbk != 0 && lch != 1 {
                            self.backgr[(lch - 1, igroup, ier)] = lbk;
                        }
                        if lbk == 1 {
                            let mut jjb = 3;
                            self.backgrdata[(lch - 1, igroup, ier, 1)] = jjb as f64;
                            {
                                let mut v = self.backgrdata.view(lch - 1, igroup, ier, jjb);
                                c.tab1io(nin, 0, 0, &mut v, &mut nb, &mut nw);
                            }
                            jjb += nw;
                            while nb != 0 {
                                let mut v = self.backgrdata.view(lch - 1, igroup, ier, jjb);
                                c.moreio(nin, 0, 0, &mut v, &mut nb, &mut nw);
                                jjb += nw;
                            }
                            self.backgrdata[(lch - 1, igroup, ier, 2)] = jjb as f64;
                            {
                                let mut v = self.backgrdata.view(lch - 1, igroup, ier, jjb);
                                c.tab1io(nin, 0, 0, &mut v, &mut nb, &mut nw);
                            }
                            jjb = *jnow + nw;
                            while nb != 0 {
                                let mut v = self.backgrdata.view(lch - 1, igroup, ier, jjb);
                                c.moreio(nin, 0, 0, &mut v, &mut nb, &mut nw);
                                jjb += nw;
                            }
                        } else if lbk == 2 || lbk == 3 {
                            c.listio(nin, 0, 0, &mut res.at(*jnow), &mut nb, &mut nw);
                            jj = *jnow + nw;
                            while nb != 0 {
                                c.moreio(nin, 0, 0, &mut res.at(jj), &mut nb, &mut nw);
                                jj += nw;
                                if jj > maxres {
                                    c.error("rdsammy", "res storage exceeded", " ");
                                }
                            }
                            self.backgrdata[(lch - 1, igroup, ier, 1)] = c.c1h;
                            self.backgrdata[(lch - 1, igroup, ier, 2)] = c.c2h;
                            self.backgrdata[(lch - 1, igroup, ier, 3)] = res[*jnow + 6];
                            self.backgrdata[(lch - 1, igroup, ier, 4)] = res[*jnow + 7];
                            self.backgrdata[(lch - 1, igroup, ier, 5)] = res[*jnow + 8];
                            if lbk == 2 {
                                self.backgrdata[(lch - 1, igroup, ier, 6)] = res[*jnow + 9];
                                self.backgrdata[(lch - 1, igroup, ier, 7)] = res[*jnow + 10];
                            }
                        }
                    }
                }
            }
            if self.ngroupm[ier as usize] > 1 {
                for igroup in 2..=self.ngroupm[ier as usize] {
                    self.nresg[(igroup, ier)] = self.nresg[(igroup - 1, ier)] + self.nresg[(igroup, ier)];
                }
            }
        } else {
            c.error("rdsammy", "only RML and RM currently coded", " ");
        }
        crate::reconr::order(enode, nodes);
        // imf2 = 0: the debug listing of the parameters is never written
    }

    fn findsp(&mut self, c: &mut Ctx, llll: i32, spinjj: f64, s1: f64, s2: f64, ier: i32) -> i32 {
        let zero = 0.0f64;
        let ig = 0;
        for i in 1..=self.ngroup {
            // the Fortran passes sspin as a 1-D array: spin(i) is sspin(i,1)
            let a = self.sspin.v[(i - 1) as usize].abs();
            if llll == self.lspin[(1, i, ier)] {
                if spinjj == a {
                    if self.chspin[(1, i, ier)] == s2 || spinjj == zero || self.findsp_ifneg == 0 {
                        return i;
                    }
                } else if spinjj == -a && self.chspin[(1, i, ier)] == s1 {
                    self.findsp_ifneg = 1;
                    return i;
                }
            }
        }
        if self.findsp_ifneg == 1 && ig == 0 {
            c.error("findsp", "quantum numbers in file 2 do not make sense", " ");
        }
        ig
    }

    fn ppdefs(&mut self, ier: i32) {
        let pnratio = AMASSP / AMASSN;
        let dnratio = AMASSD / AMASSN;
        let tnratio = AMASST / AMASSN;
        let hnratio = AMASSH / AMASSN;
        let anratio = AMASSA / AMASSN;
        for i in 1..=self.nppm[ier as usize] {
            let m = self.mt[(i, ier)];
            if m == 102 {
            } else if m == 2 {
                self.ema[(i, ier)] = 1.0;
                self.spina[(i, ier)] = 0.5;
                self.kza[(i, ier)] = 0;
                self.lpent[(i, ier)] = 1;
            } else if m == 18 {
                self.lpent[(i, ier)] = 0;
            } else if m > 50 && m < 99 {
                self.ema[(i, ier)] = 1.0;
                self.spina[(i, ier)] = 0.5;
                self.kza[(i, ier)] = 0;
                self.lpent[(i, ier)] = 1;
            } else if m == 103 {
                self.ema[(i, ier)] = pnratio;
                self.spina[(i, ier)] = 0.5;
                self.kza[(i, ier)] = 1;
                self.lpent[(i, ier)] = 1;
            } else if m == 104 {
                self.ema[(i, ier)] = dnratio;
                self.spina[(i, ier)] = 1.0;
                self.kza[(i, ier)] = 1;
                self.lpent[(i, ier)] = 1;
            } else if m == 105 {
                self.ema[(i, ier)] = tnratio;
                self.spina[(i, ier)] = 0.5;
                self.kza[(i, ier)] = 1;
                self.lpent[(i, ier)] = 1;
            } else if m == 106 {
                self.ema[(i, ier)] = hnratio;
                self.spina[(i, ier)] = 0.5;
                self.kza[(i, ier)] = 2;
                self.lpent[(i, ier)] = 1;
            } else if m == 107 {
                self.ema[(i, ier)] = anratio;
                self.spina[(i, ier)] = 0.0;
                self.pa[(i, ier)] = 0.0;
                self.kza[(i, ier)] = 2;
                self.lpent[(i, ier)] = 1;
            }
        }
    }

    fn rearrange(&mut self, c: &mut Ctx, idum: &mut [i32], ier: i32) {
        let nr = self.nresm[ier as usize];
        for _i in 1..=nr {
            let mut isort = 0;
            for j in 2..=nr {
                let (j0, j1) = (j as usize, (j - 1) as usize);
                if idum[j0] < idum[j1] {
                    isort = 1;
                    idum.swap(j0, j1);
                    let a = self.eres[(j, ier)];
                    self.eres[(j, ier)] = self.eres[(j - 1, ier)];
                    self.eres[(j - 1, ier)] = a;
                    let a = self.gamgam[(j, ier)];
                    self.gamgam[(j, ier)] = self.gamgam[(j - 1, ier)];
                    self.gamgam[(j - 1, ier)] = a;
                    for k in 1..=self.mchan {
                        let a = self.gamma[(k, j, ier)];
                        self.gamma[(k, j, ier)] = self.gamma[(k, j - 1, ier)];
                        self.gamma[(k, j - 1, ier)] = a;
                    }
                }
            }
            if isort == 0 {
                break;
            }
        }
        for i in 1..=self.ngroupm[ier as usize] {
            self.nresg[(i, ier)] = 0;
        }
        for j in 1..=nr {
            let g = idum[j as usize];
            self.nresg[(g, ier)] += 1;
        }
        if self.ngroupm[ier as usize] > 1 {
            for i in 2..=self.ngroupm[ier as usize] {
                self.nresg[(i, ier)] = self.nresg[(i - 1, ier)] + self.nresg[(i, ier)];
            }
        }
        if self.nresm[ier as usize] != self.nresg[(self.ngroup, ier)] {
            c.error("rearrange", "nres fault", " ");
        }
    }

    // ============================================================ ppsammy
    pub fn ppsammy(&mut self, c: &mut Ctx, ier: i32, ncoef: &mut i32, nresp: &mut i32) {
        for i in 1..=self.mchan {
            for j in 1..=self.ngroup {
                self.echan[(i, j, ier)] = 0.0;
                self.zeta[(i, j, ier)] = 0.0;
            }
        }
        self.checkqn(c, ier);
        self.fxradi(ier);
        self.betset(c, ier);
        if self.npar > 0 {
            *nresp = self.npar;
        }
        *ncoef = 1;
        // Want_Angular_Dist / Want_Partial_Derivs: flag-guarded (false)
    }

    fn allo(&mut self) {
        let (npp, nerm, ng, mc, nres) = (self.npp, self.nerm, self.ngroup, self.mchan, self.nres);
        self.ntriag = (mc * (mc + 1)) / 2;
        let nt = self.ntriag;
        self.ema = FA::new(&[npp, nerm]);
        self.emb = FA::new(&[npp, nerm]);
        self.kza = FI::new(&[npp, nerm]);
        self.kzb = FI::new(&[npp, nerm]);
        self.spina = FA::new(&[npp, nerm]);
        self.spinb = FA::new(&[npp, nerm]);
        self.qqq = FA::new(&[npp, nerm]);
        self.ishift = FI::new(&[npp, nerm]);
        self.lpent = FI::new(&[npp, nerm]);
        self.mt = FI::new(&[npp, nerm]);
        self.pa = FA::new(&[npp, nerm]);
        self.pb = FA::new(&[npp, nerm]);
        self.sspin = FA::new(&[ng, nerm]);
        self.parity = FA::new(&[ng, nerm]);
        self.nresg = FI::new(&[ng, nerm]);
        self.ipp = FI::new(&[mc, ng, nerm]);
        self.lspin = FI::new(&[mc, ng, nerm]);
        self.backgr = FI::new(&[mc, ng, nerm]);
        self.backgrdata = FA::new(&[mc, ng, nerm, self.nback]);
        self.chspin = FA::new(&[mc, ng, nerm]);
        self.bound = FA::new(&[mc, ng, nerm]);
        self.rdeff = FA::new(&[mc, ng, nerm]);
        self.rdtru = FA::new(&[mc, ng, nerm]);
        self.eres = FA::new(&[nres, nerm]);
        self.gamgam = FA::new(&[nres, nerm]);
        self.gamma = FA::new(&[mc, nres, nerm]);
        self.nent = FI::new(&[nres, nerm]);
        self.next = FI::new(&[nres, nerm]);
        self.goj = FA::new(&[nres, nerm]);
        self.zke = FA::new(&[mc, ng, nerm]);
        self.zkfe = FA::new(&[mc, ng, nerm]);
        self.zkte = FA::new(&[mc, ng, nerm]);
        self.zeta = FA::new(&[mc, ng, nerm]);
        self.echan = FA::new(&[mc, ng, nerm]);
        self.betapr = FA::new(&[mc, nres, nerm]);
        self.gbetpr = FA::new(&[3, nres, nerm]);
        self.beta = FA::new(&[nt, nres, nerm]);
        self.uuuu = FA::new(&[(mc + 2) * nres, nerm]);
        self.duuu = FA::new(&[(mc + 2) * nres, nerm]);
        self.iduu = FI::new(&[(mc + 2) * nres, nerm]);
        self.crss = FA::new(&[npp, nerm]);
        self.sigmas = FA::new(&[npp, nerm]);
        self.alphar = FA::new(&[nres, nerm]);
        self.alphai = FA::new(&[nres, nerm]);
        self.difen = FA::new(&[nres, nerm]);
        self.xden = FA::new(&[nres, nerm]);
        self.sinsqr = FA::new(&[mc, nerm]);
        self.sin2ph = FA::new(&[mc, nerm]);
        self.dphi = FA::new(&[mc, nerm]);
        self.dpdr = FA::new(&[mc, nerm]);
        self.dsdr = FA::new(&[mc, nerm]);
        self.cosphi = FA::new(&[mc, nerm]);
        self.sinphi = FA::new(&[mc, nerm]);
        self.rootp = FA::new(&[mc, nerm]);
        self.elinvr = FA::new(&[mc, nerm]);
        self.elinvi = FA::new(&[mc, nerm]);
        self.psmall = FA::new(&[mc, nerm]);
        self.xqr = FA::new(&[mc, mc, nerm]);
        self.xqi = FA::new(&[mc, mc, nerm]);
        self.xxxxr = FA::new(&[nt, nerm]);
        self.xxxxi = FA::new(&[nt, nerm]);
        self.qr = FA::new(&[nt, nt, nerm]);
        self.qi = FA::new(&[nt, nt, nerm]);
        self.rmat = FA::new(&[2, nt, nerm]);
        self.ymat = FA::new(&[2, nt, nerm]);
        self.yinv = FA::new(&[2, nt, nerm]);
        self.cscs = FA::new(&[2, nt, nerm]);
        self.ddddd = FA::new(&[npp, nerm]);
        self.par = FA::new(&[self.npar]);
    }

    fn checkqn(&mut self, c: &mut Ctx, ier: i32) {
        let zero = 0.0f64;
        let half = 0.5e0f64;
        for j in 1..=self.ngroupm[ier as usize] {
            if crate::fmath::fmod(self.sspin[(j, ier)], half) != zero {
                let strng = format!("spin group number {} with spin ={}", ffmt::i(j, 3), ffmt::f(self.sspin[(j, ier)], 10, 5));
                c.error("checkqn", "error in quantum numbers", &strng);
            }
            self.goj[(j, ier)] = 0.0;
            self.nent[(j, ier)] = 0;
            self.next[(j, ier)] = 0;
            let nchanj = self.nchan[(j, ier)];
            for n in 1..=nchanj {
                let ip = self.ipp[(n, j, ier)];
                if ip == 2 {
                    self.nent[(j, ier)] += 1;
                    if self.goj[(j, ier)] == zero {
                        self.goj[(j, ier)] = (2.0 * self.sspin[(j, ier)].abs() + 1.0)
                            / ((2.0 * self.spina[(ip, ier)].abs() + 1.0) * (2.0 * self.spinb[(ip, ier)].abs() + 1.0));
                    }
                } else if ip > 2 {
                    self.next[(j, ier)] += 1;
                }
                if self.lspin[(n, j, ier)] < 0 {
                    wrongi(c, "lspin", self.lspin[(n, j, ier)]);
                }
                if self.chspin[(n, j, ier)] == zero || self.sspin[(j, ier)] == zero {
                } else {
                    let mut x = 1.0;
                    let i = self.lspin[(n, j, ier)] % 2;
                    if i != 0 {
                        x = -x;
                    }
                    if self.chspin[(n, j, ier)] < zero {
                        x = -x;
                    }
                    if self.sspin[(j, ier)] < zero {
                        x = -x;
                    }
                    if x < zero {
                        c.io.so(&format!(
                            " *** Parity problem ***\n     Group and channel #{}{}\n     Spin, L, Chspin ={}{}{}",
                            ffmt::i(j, 5),
                            ffmt::i(n, 5),
                            ffmt::f(self.sspin[(j, ier)], 7, 1),
                            ffmt::i(self.lspin[(n, j, ier)], 4),
                            ffmt::f(self.chspin[(n, n, ier)], 7, 1)
                        ));
                    }
                }
                let sa = self.spina[(ip, ier)].abs();
                let sb = self.spinb[(ip, ier)].abs();
                let mut smin = (sa - sb).abs();
                let mut smax = sa + sb;
                let ch = self.chspin[(n, j, ier)];
                let hdr = format!(
                    "\n ****** Error in quantum numbers for Group Number{} and Channel Number{}",
                    ffmt::i(j, 3),
                    ffmt::i(n, 2)
                );
                if ch.abs() < smin {
                    c.io.so(&hdr);
                    c.io.so(&format!(
                        " ******* |Chspin|<|Spina-Spinb|=Smin   {}{}{}{}",
                        ffmt::f(ch, 5, 1),
                        ffmt::f(self.spina[(ip, ier)], 5, 1),
                        ffmt::f(self.spinb[(ip, ier)], 5, 1),
                        ffmt::f(smin, 5, 1)
                    ));
                } else if ch.abs() > smax {
                    c.io.so(&hdr);
                    c.io.so(&format!(
                        " ******* |Chspin|>|Spina+Spinb|=Smax   {}{}{}{}",
                        ffmt::f(ch, 5, 1),
                        ffmt::f(self.spina[(ip, ier)], 5, 1),
                        ffmt::f(self.spinb[(ip, ier)], 5, 1),
                        ffmt::f(smax, 5, 1)
                    ));
                }
                if self.lpent[(ip, ier)] > 0 {
                    let ls = self.lspin[(n, j, ier)];
                    smin = ((ls as f32) as f64 - ch.abs()).abs().abs();
                    smax = ls as f64 + ch.abs();
                    let sj = self.sspin[(j, ier)];
                    if sj.abs() < smin {
                        c.io.so(&hdr);
                        c.io.so(&format!(
                            " ****** |Spin|<|Lspin-Chspin|=Smin   {}{}{}{}",
                            ffmt::f(sj, 5, 1),
                            ffmt::i(ls, 5),
                            ffmt::f(ch, 5, 1),
                            ffmt::f(smin, 5, 1)
                        ));
                    } else if sj.abs() > smax {
                        c.io.so(&hdr);
                        c.io.so(&format!(
                            " ****** |Spin|>|Lspin-chspin|=Smin   {}{}{}{}",
                            ffmt::f(sj, 5, 1),
                            ffmt::i(ls, 5),
                            ffmt::f(ch, 5, 1),
                            ffmt::f(smin, 5, 1)
                        ));
                    }
                }
            }
        }
    }

    fn fxradi(&mut self, ier: i32) {
        let hbarrr = HBAR / EV;
        let amuevv = AMU * CLIGHT * CLIGHT / EV;
        let cspeed = CLIGHT / 100.0;
        let fininv = 1.0e0 / FINSTRI;
        let zero = 0.0f64;
        let ff = 1.0e15f64;
        let twomhb = sqrt(2.0 * AMASSN * amuevv) / (hbarrr * ff * cspeed);
        let etac = fininv * amuevv / (hbarrr * ff * cspeed) * AMASSN;
        for ippx in 1..=self.nppm[ier as usize] {
            if self.ema[(ippx, ier)] == zero {
                self.ema[(ippx, ier)] = 1.0;
            }
            if self.emb[(ippx, ier)] == zero {
                self.emb[(ippx, ier)] = self.awr;
            }
        }
        for kgroup in 1..=self.ngroupm[ier as usize] {
            let mut factor = self.emb[(2, ier)] + self.ema[(2, ier)];
            let alabcm = self.emb[(2, ier)] / factor;
            factor = alabcm / self.ema[(2, ier)];
            for ichan in 1..=self.nchan[(kgroup, ier)] {
                let ippx = self.ipp[(ichan, kgroup, ier)];
                let mut aa = self.emb[(ippx, ier)] + self.ema[(ippx, ier)];
                aa = self.emb[(ippx, ier)] / aa;
                if self.qqq[(ippx, ier)] != zero && self.echan[(ichan, kgroup, ier)] == zero {
                    self.echan[(ichan, kgroup, ier)] = -self.qqq[(ippx, ier)] / alabcm;
                }
                let redmas = aa * self.ema[(ippx, ier)];
                let z = twomhb * sqrt(redmas * factor);
                self.zke[(ichan, kgroup, ier)] = z;
                self.zkfe[(ichan, kgroup, ier)] = z * self.rdeff[(ichan, kgroup, ier)] * 10.0;
                self.zkte[(ichan, kgroup, ier)] = z * self.rdtru[(ichan, kgroup, ier)] * 10.0;
                if self.kzb[(ippx, ier)] != 0 && self.kza[(ippx, ier)] != 0 {
                    self.zeta[(ichan, kgroup, ier)] = etac * self.kzb[(ippx, ier)] as f64 * self.kza[(ippx, ier)] as f64
                        * redmas
                        / self.zke[(ichan, kgroup, ier)];
                }
            }
        }
    }

    fn betset(&mut self, c: &mut Ctx, ier: i32) {
        let zero = 0.0f64;
        let half = 0.5e0f64;
        let knpar = 0;
        if self.nres > 0 {
            let mut minres = 1;
            for ig in 1..=self.ngroupm[ier as usize] {
                let mmaxc = self.nchan[(ig, ier)];
                for j in 1..=mmaxc {
                    let lsp = self.lspin[(j, ig, ier)];
                    let ipx = self.ipp[(j, ig, ier)];
                    for ires in minres..=self.nresg[(ig, ier)] {
                        self.betapr[(j, ires, ier)] = 0.0;
                        let mut p = 1.0;
                        let g = self.gamma[(j, ires, ier)];
                        if g != zero {
                            if self.lpent[(ipx, ier)] <= 0 {
                                self.betapr[(j, ires, ier)] = sqrt(half * g.abs());
                                if g < zero {
                                    self.betapr[(j, ires, ier)] = -self.betapr[(j, ires, ier)];
                                }
                            } else {
                                let mut ex = (self.eres[(ires, ier)] - self.echan[(j, ig, ier)]).abs();
                                if ex != zero {
                                    ex = sqrt(ex);
                                    let rho = self.zkte[(j, ig, ier)] * ex;
                                    if self.zeta[(j, ig, ier)] == zero {
                                        p = pf(rho, lsp);
                                        let (mut hr, mut hi, mut dp, mut ds, mut iffy) = (0.0, 0.0, 0.0, 0.0, 0);
                                        pgh(rho, lsp, self.bound[(j, ig, ier)], &mut hr, &mut hi, &mut p, &mut dp, &mut ds, self.ishift[(ipx, ier)], &mut iffy);
                                    } else {
                                        let eta = self.zeta[(j, ig, ier)] / ex;
                                        let jdopha = 0;
                                        let (mut hr, mut hi, mut dp, mut ds, mut iffy) = (0.0, 0.0, 0.0, 0.0, 0);
                                        let (mut sp, mut cp, mut dph) = (0.0, 0.0, 0.0);
                                        pghcou(
                                            c, rho, lsp, self.bound[(j, ig, ier)], &mut hr, &mut hi, &mut p, &mut dp, &mut ds,
                                            self.ishift[(ipx, ier)], &mut iffy, eta, &mut sp, &mut cp, &mut dph, jdopha,
                                        );
                                        if p < zero {
                                            p = 1.0;
                                            let s = format!("p<0 set to 1 at eres={}", ffmt::e1p(self.eres[(ires, ier)], 12, 4));
                                            c.mess("betset", &s, " ");
                                        } else if p == zero {
                                            p = 1.0;
                                            let s = format!("p=0 set to 1 at eres={}", ffmt::e1p(self.eres[(ires, ier)], 12, 4));
                                            c.mess("betset", &s, " ");
                                        }
                                    }
                                }
                            }
                            self.betapr[(j, ires, ier)] = sqrt(half * g.abs() / p);
                            if g < zero {
                                self.betapr[(j, ires, ier)] = -self.betapr[(j, ires, ier)];
                            }
                            // dum(j,ires) is only used for derivatives
                        }
                    }
                }
                for ires in minres..=self.nresg[(ig, ier)] {
                    let mut kl = 0;
                    for k in 1..=mmaxc {
                        for l in 1..=k {
                            kl += 1;
                            self.beta[(kl, ires, ier)] = self.betapr[(l, ires, ier)] * self.betapr[(k, ires, ier)];
                        }
                    }
                }
                for ires in minres..=self.nresg[(ig, ier)] {
                    self.gbetpr[(2, ires, ier)] = self.gamgam[(ires, ier)].abs() / 2.0;
                    self.gbetpr[(1, ires, ier)] = sqrt(self.gbetpr[(2, ires, ier)]);
                    if self.gamgam[(ires, ier)] < zero {
                        self.gbetpr[(1, ires, ier)] = -self.gbetpr[(1, ires, ier)];
                    }
                    self.gbetpr[(3, ires, ier)] = sq(self.gbetpr[(2, ires, ier)]);
                }
                // Want_Partial_Derivs: flag-guarded (false)
                minres = self.nresg[(ig, ier)] + 1;
            }
        }
        self.nparr = knpar;
        self.npar = knpar;
    }

    fn abpart(&mut self, ier: i32) {
        for n in 1..=self.nresm[ier as usize] {
            self.xden[(n, ier)] = 0.0;
            self.alphar[(n, ier)] = 0.0;
            self.alphai[(n, ier)] = 0.0;
            self.difen[(n, ier)] = self.eres[(n, ier)] - self.su;
            let aa = sq(self.difen[(n, ier)]) + self.gbetpr[(3, n, ier)];
            self.xden[(n, ier)] = 1.0 / aa;
            self.alphar[(n, ier)] = self.difen[(n, ier)] * self.xden[(n, ier)];
            self.alphai[(n, ier)] = self.gbetpr[(2, n, ier)] * self.xden[(n, ier)];
        }
        // Want_Partial_Derivs: flag-guarded (false)
    }

    fn crosss(&mut self, c: &mut Ctx, ier: i32) {
        let fourpi = 4.0 * PI / 100.0;
        for i in 1..=self.nppm[ier as usize] {
            self.sigmas[(i, ier)] = 0.0;
        }
        let mut maxr = 0;
        for n in 1..=self.ngroupm[ier as usize] {
            for i in 1..=self.nppm[ier as usize] {
                self.crss[(i, ier)] = 0.0;
            }
            let minr = maxr + 1;
            maxr = self.nresg[(n, ier)];
            let mut nchann = self.nchan[(n, ier)];
            let nentnn = self.nent[(n, ier)];
            let nextnn = self.next[(n, ier)];
            let mut lrmat = 0;
            self.setr(c, n, &mut lrmat, minr, maxr, nentnn, &mut nchann, ier);
            if lrmat == 1 {
                self.zeror(nchann, ier);
            } else {
                self.yinvrs(c, nchann, ier);
                self.setxqx(nchann, ier);
            }
            self.sectio(n, nentnn, nextnn, nchann, ier);
            for ip in 1..=self.nppm[ier as usize] {
                self.sigmas[(ip, ier)] = self.crss[(ip, ier)] + self.sigmas[(ip, ier)];
            }
        }
        for ip in 1..=self.nppm[ier as usize] {
            self.sigmas[(ip, ier)] = self.sigmas[(ip, ier)] * fourpi / self.su;
        }
    }

    fn setr(&mut self, c: &mut Ctx, n: i32, lrmat: &mut i32, min: i32, max: i32, _nent: i32, nchan: &mut i32, ier: i32) {
        let zero = 0.0f64;
        let one = 1.0f64;
        let tiny = 1.0e-8f64;
        let su = self.su;
        let mut nnntot = *nchan + 1;
        for _i in 1..=*nchan {
            nnntot -= 1;
            if su >= zero && su <= self.echan[(nnntot, n, ier)] {
                *nchan = nnntot - 1;
            } else {
                break;
            }
        }
        let nchan = *nchan;
        let mut kl = 0;
        for k in 1..=nchan {
            for l in 1..=k {
                kl += 1;
                self.rmat[(1, kl, ier)] = 0.0;
                self.rmat[(2, kl, ier)] = 0.0;
                if l == k && self.backgr[(k, n, ier)] > 0 {
                    let bg = self.backgr[(k, n, ier)];
                    let bd = |s: &Samm, m: i32| s.backgrdata[(k, n, ier, m)];
                    if bg == 1 {
                        let mut ip = 2;
                        let mut ir = 1;
                        let (mut v, mut enext, mut idis) = (0.0, 0.0, 0);
                        {
                            let o = self.backgrdata.off(k, n, ier, 3);
                            let a = F(&mut self.backgrdata.v[o..]);
                            c.terpa(&mut v, su, &mut enext, &mut idis, &a, &mut ip, &mut ir);
                        }
                        self.rmat[(1, kl, ier)] = v;
                        ip = 2;
                        ir = 1;
                        let m2 = crate::fmath::int(bd(self, 2));
                        {
                            let o = self.backgrdata.off(k, n, ier, m2);
                            let a = F(&mut self.backgrdata.v[o..]);
                            c.terpa(&mut v, su, &mut enext, &mut idis, &a, &mut ip, &mut ir);
                        }
                        self.rmat[(2, kl, ier)] = v;
                    } else if bg == 2 {
                        self.rmat[(1, kl, ier)] = bd(self, 3) - bd(self, 7) * (bd(self, 2) - bd(self, 1))
                            + (bd(self, 4) + bd(self, 5) * su) * su
                            - (bd(self, 6) + bd(self, 7) * su) * log((bd(self, 2) - su) / (su - bd(self, 1)));
                    } else if bg == 3 {
                        let esum = bd(self, 1) + bd(self, 2);
                        let ediff = bd(self, 2) - bd(self, 1);
                        self.rmat[(1, kl, ier)] = bd(self, 3) + 2.0 * bd(self, 4) * atanh((2.0 * su - esum) / ediff);
                        self.rmat[(2, kl, ier)] = bd(self, 5) / ediff / (1.0 - sq((2.0 * su - esum) / ediff));
                    }
                }
            }
        }
        if max >= min && min > 0 {
            for ires in min..=max {
                kl = 0;
                for k in 1..=nchan {
                    for l in 1..=k {
                        kl += 1;
                        if su > self.echan[(k, n, ier)] && su > self.echan[(l, n, ier)] && self.beta[(kl, ires, ier)] != zero {
                            self.rmat[(1, kl, ier)] = self.rmat[(1, kl, ier)] + self.alphar[(ires, ier)] * self.beta[(kl, ires, ier)];
                            self.rmat[(2, kl, ier)] = self.rmat[(2, kl, ier)] + self.alphai[(ires, ier)] * self.beta[(kl, ires, ier)];
                        }
                    }
                }
            }
        }
        let mut nonzero = false;
        kl = 0;
        'outer: for k in 1..=nchan {
            for _l in 1..=k {
                kl += 1;
                if self.rmat[(1, kl, ier)] != zero || self.rmat[(2, kl, ier)] != zero {
                    nonzero = true;
                    break 'outer;
                }
            }
        }
        if !nonzero {
            *lrmat = 1;
        }
        kl = 0;
        for k in 1..=nchan {
            for _l in 1..=k {
                kl += 1;
                self.ymat[(1, kl, ier)] = -self.rmat[(1, kl, ier)];
                self.ymat[(2, kl, ier)] = -self.rmat[(2, kl, ier)];
            }
        }
        for i in 1..=nchan {
            self.psmall[(i, ier)] = 0.0;
        }
        let mut ii = 0;
        for i in 1..=nchan {
            let ipx = self.ipp[(i, n, ier)];
            ii += i;
            self.rootp[(i, ier)] = 1.0;
            self.elinvr[(i, ier)] = 0.0;
            self.elinvi[(i, ier)] = -1.0;
            self.dpdr[(i, ier)] = 0.0;
            self.dsdr[(i, ier)] = 0.0;
            let mut iffy = 0;
            if su > self.echan[(i, n, ier)] {
                if self.lpent[(ipx, ier)] <= 0 {
                    self.ymat[(2, ii, ier)] = self.ymat[(2, ii, ier)] - 1.0;
                } else {
                    let lsp = self.lspin[(i, n, ier)];
                    let ex = sqrt(su - self.echan[(i, n, ier)]);
                    let rho = self.zkte[(i, n, ier)] * ex;
                    let rhof = self.zkfe[(i, n, ier)] * ex;
                    let (mut hr, mut hi, mut p, mut dp, mut ds) = (0.0, 0.0, 0.0, 0.0, 0.0);
                    if self.zeta[(i, n, ier)] == zero {
                        let (mut a, mut b, mut cc, mut d, mut e) = (0.0, 0.0, 0.0, 0.0, 0.0);
                        sinsix(rhof, lsp, &mut a, &mut b, &mut cc, &mut d, &mut e);
                        self.sinsqr[(i, ier)] = a;
                        self.sin2ph[(i, ier)] = b;
                        self.dphi[(i, ier)] = cc;
                        self.sinphi[(i, ier)] = d;
                        self.cosphi[(i, ier)] = e;
                        pgh(rho, lsp, self.bound[(i, n, ier)], &mut hr, &mut hi, &mut p, &mut dp, &mut ds, self.ishift[(ipx, ier)], &mut iffy);
                    } else {
                        let eta = self.zeta[(i, n, ier)] / ex;
                        let (mut sp, mut cp, mut dph) = (self.sinphi[(i, ier)], self.cosphi[(i, ier)], self.dphi[(i, ier)]);
                        if rho == rhof {
                            let jdopha = 1;
                            pghcou(
                                c, rho, lsp, self.bound[(i, n, ier)], &mut hr, &mut hi, &mut p, &mut dp, &mut ds,
                                self.ishift[(ipx, ier)], &mut iffy, eta, &mut sp, &mut cp, &mut dph, jdopha,
                            );
                        } else {
                            let jdopha = 0;
                            let (mut sinphx, mut cosphx, mut dphix) = (0.0, 0.0, 0.0);
                            pghcou(
                                c, rho, lsp, self.bound[(i, n, ier)], &mut hr, &mut hi, &mut p, &mut dp, &mut ds,
                                self.ishift[(ipx, ier)], &mut iffy, eta, &mut sinphx, &mut cosphx, &mut dphix, jdopha,
                            );
                            let jdopha = 1;
                            let (mut hrx, mut hix, mut px, mut dpx, mut dsx) = (0.0, 0.0, 0.0, 0.0, 0.0);
                            pghcou(
                                c, rhof, lsp, self.bound[(i, n, ier)], &mut hrx, &mut hix, &mut px, &mut dpx, &mut dsx,
                                self.ishift[(ipx, ier)], &mut iffy, eta, &mut sp, &mut cp, &mut dph, jdopha,
                            );
                        }
                        self.sinphi[(i, ier)] = sp;
                        self.cosphi[(i, ier)] = cp;
                        self.dphi[(i, ier)] = dph;
                        self.sinsqr[(i, ier)] = sq(sp);
                        self.sin2ph[(i, ier)] = 2.0 * sp * cp;
                    }
                    self.rootp[(i, ier)] = sqrt(p);
                    self.dpdr[(i, ier)] = dp;
                    self.dsdr[(i, ier)] = ds;
                    let ish = self.ishift[(ipx, ier)];
                    if iffy == 0 && !(ish <= 0 && (1.0 - p * self.rmat[(2, ii, ier)] == one || p < tiny)) {
                        self.elinvr[(i, ier)] = hr;
                        self.elinvi[(i, ier)] = hi;
                        self.ymat[(1, ii, ier)] = hr + self.ymat[(1, ii, ier)];
                        self.ymat[(2, ii, ier)] = hi + self.ymat[(2, ii, ier)];
                    } else {
                        self.psmall[(i, ier)] = self.rootp[(i, ier)];
                        self.ymat[(1, ii, ier)] = p * self.ymat[(1, ii, ier)];
                        self.ymat[(2, ii, ier)] = p * self.ymat[(2, ii, ier)] - one;
                        self.rmat[(1, ii, ier)] = p * self.rmat[(1, ii, ier)];
                        self.rmat[(2, ii, ier)] = p * self.rmat[(2, ii, ier)];
                        if nchan > 1 {
                            let rp = self.rootp[(i, ier)];
                            if i > 1 {
                                for j in 1..=i - 1 {
                                    let ji = (i * (i - 1)) / 2 + j;
                                    self.ymat[(1, ji, ier)] = rp * self.ymat[(1, ji, ier)];
                                    self.ymat[(2, ji, ier)] = rp * self.ymat[(2, ji, ier)];
                                    self.rmat[(1, ji, ier)] = rp * self.rmat[(1, ji, ier)];
                                    self.rmat[(2, ji, ier)] = rp * self.rmat[(2, ji, ier)];
                                }
                            }
                            if i < nchan {
                                for j in i + 1..=nchan {
                                    let ji = (j * (j - 1)) / 2 + i;
                                    self.ymat[(1, ji, ier)] = rp * self.ymat[(1, ji, ier)];
                                    self.ymat[(2, ji, ier)] = rp * self.ymat[(2, ji, ier)];
                                    self.rmat[(1, ji, ier)] = rp * self.rmat[(1, ji, ier)];
                                    self.rmat[(2, ji, ier)] = rp * self.rmat[(2, ji, ier)];
                                }
                            }
                        }
                        self.rootp[(i, ier)] = 1.0;
                        self.elinvr[(i, ier)] = 0.0;
                        self.elinvi[(i, ier)] = -1.0;
                        self.dpdr[(i, ier)] = 0.0;
                        self.dsdr[(i, ier)] = 0.0;
                    }
                }
            } else {
                self.rootp[(i, ier)] = 0.0;
                self.elinvr[(i, ier)] = 1.0;
                self.elinvi[(i, ier)] = 0.0;
                self.dpdr[(i, ier)] = 0.0;
                self.dsdr[(i, ier)] = 0.0;
            }
        }
        if *lrmat != 1 {
            kl = 0;
            for k in 1..=nchan {
                kl += k;
                if self.ymat[(1, kl, ier)] == zero && self.ymat[(2, kl, ier)] == zero {
                    self.ymat[(1, kl, ier)] = 1.0;
                }
            }
        }
    }

    fn zeror(&mut self, nchan: i32, ier: i32) {
        let mut jk = 0;
        for j in 1..=nchan {
            for _k in 1..=j {
                jk += 1;
                self.xxxxr[(jk, ier)] = 0.0;
                self.xxxxi[(jk, ier)] = 0.0;
            }
        }
        for j in 1..=nchan {
            for k in 1..=nchan {
                self.xqr[(k, j, ier)] = 0.0;
                self.xqi[(k, j, ier)] = 0.0;
            }
        }
    }

    fn yinvrs(&mut self, c: &mut Ctx, nchan: i32, ier: i32) {
        if nchan == 1 {
            onech(&self.ymat, &mut self.yinv, ier);
        } else if nchan == 2 {
            twoch(&mut self.ymat, &mut self.yinv, ier);
        } else if nchan == 3 {
            threech(&mut self.ymat, &mut self.yinv, ier);
        } else {
            yfour(c, &mut self.ymat, &mut self.yinv, nchan);
        }
    }

    fn setxqx(&mut self, nchan: i32, ier: i32) {
        for i in 1..=nchan {
            for j in 1..=nchan {
                self.xqr[(i, j, ier)] = 0.0;
                self.xqi[(i, j, ier)] = 0.0;
            }
        }
        for i in 1..=nchan {
            for j in 1..=nchan {
                let ij = ijkl(j, i);
                for k in 1..=nchan {
                    let jk = ijkl(k, j);
                    self.xqr[(k, i, ier)] = self.xqr[(k, i, ier)] + self.yinv[(1, ij, ier)] * self.rmat[(1, jk, ier)]
                        - self.yinv[(2, ij, ier)] * self.rmat[(2, jk, ier)];
                    self.xqi[(k, i, ier)] = self.xqi[(k, i, ier)] + self.yinv[(1, ij, ier)] * self.rmat[(2, jk, ier)]
                        + self.yinv[(2, ij, ier)] * self.rmat[(1, jk, ier)];
                }
            }
        }
        let mut ij = 0;
        for i in 1..=nchan {
            let plr = self.rootp[(i, ier)] * self.elinvr[(i, ier)];
            let pli = self.rootp[(i, ier)] * self.elinvi[(i, ier)];
            for j in 1..=i {
                ij += 1;
                self.xxxxr[(ij, ier)] = self.rootp[(j, ier)] * (self.xqr[(j, i, ier)] * plr - self.xqi[(j, i, ier)] * pli);
                self.xxxxi[(ij, ier)] = self.rootp[(j, ier)] * (self.xqi[(j, i, ier)] * plr + self.xqr[(j, i, ier)] * pli);
            }
        }
    }

    fn sectio(&mut self, n: i32, nent: i32, next: i32, nchan: i32, ier: i32) {
        for jj in 1..=self.npp {
            self.crss[(jj, ier)] = 0.0;
        }
        let mut ii = 0;
        let mut ij = 0;
        for i in 1..=nent {
            let zz = sq(self.zke[(i, n, ier)]);
            ii += i;
            let mut termn = self.sinsqr[(i, ier)] * (1.0 - 2.0 * self.xxxxi[(ii, ier)]) - self.sin2ph[(i, ier)] * self.xxxxr[(ii, ier)];
            termn = termn / zz;
            for j in 1..=i {
                ij += 1;
                let mut ar = (sq(self.xxxxr[(ij, ier)]) + sq(self.xxxxi[(ij, ier)])) / zz;
                if i != j {
                    ar = ar + ar;
                }
                termn = termn + ar;
            }
            self.crss[(1, ier)] = termn + self.crss[(1, ier)];
        }
        self.crss[(1, ier)] = self.crss[(1, ier)] * self.goj[(n, ier)];
        ii = 0;
        ij = 0;
        for i in 1..=nent {
            ii += i;
            let zz = sq(self.zke[(i, n, ier)]);
            let mut terma = self.xxxxi[(ii, ier)] / zz;
            for j in 1..=i {
                ij += 1;
                let mut ar = (-sq(self.xxxxr[(ij, ier)]) - sq(self.xxxxi[(ij, ier)])) / zz;
                if i != j {
                    ar = ar + ar;
                }
                terma = terma + ar;
            }
            self.crss[(2, ier)] = terma + self.crss[(2, ier)];
        }
        self.crss[(2, ier)] = self.crss[(2, ier)] * self.goj[(n, ier)];
        for jj in 1..=next {
            let j = jj + nent;
            if j <= nchan {
                let ip = self.ipp[(j, n, ier)];
                for i in 1..=nent {
                    let zz = sq(self.zke[(i, n, ier)]);
                    let ij = (j * (j - 1)) / 2 + i;
                    self.crss[(ip, ier)] = self.crss[(ip, ier)] + (sq(self.xxxxr[(ij, ier)]) + sq(self.xxxxi[(ij, ier)])) / zz;
                }
            }
        }
        if self.npp > 2 {
            for ip in 3..=self.npp {
                self.crss[(ip, ier)] = self.crss[(ip, ier)] * self.goj[(n, ier)];
            }
        }
        // Want_Angular_Dist: flag-guarded (false)
    }

    // ============================================================ desammy
    pub fn desammy(&mut self) {
        let keep = (self.want_partial_derivs, self.want_angular_dist, self.want_samrml_bw, self.want_samrml_rm);
        let ifneg = self.findsp_ifneg;
        let (ngroup, mchan, nres, npp, ntriag, npar, lllmax, nback, nrext, nparr, nerm) = (
            self.ngroup, self.mchan, self.nres, self.npp, self.ntriag, self.npar, self.lllmax, self.nback, self.nrext,
            self.nparr, self.nerm,
        );
        let (ngroupm, nppm, nresm, awr, su) = (self.ngroupm, self.nppm, self.nresm, self.awr, self.su);
        *self = Samm::default();
        self.want_partial_derivs = keep.0;
        self.want_angular_dist = keep.1;
        self.want_samrml_bw = keep.2;
        self.want_samrml_rm = keep.3;
        self.findsp_ifneg = ifneg;
        self.ngroup = ngroup;
        self.mchan = mchan;
        self.nres = nres;
        self.npp = npp;
        self.ntriag = ntriag;
        self.npar = npar;
        self.lllmax = lllmax;
        self.nback = nback;
        self.nrext = nrext;
        self.nparr = nparr;
        self.nerm = nerm;
        self.ngroupm = ngroupm;
        self.nppm = nppm;
        self.nresm = nresm;
        self.awr = awr;
        self.su = su;
    }
}

fn wrongi(c: &mut Ctx, name: &str, i: i32) {
    let s = format!("\n $$$ Variable \"{}\" probably has wrong value{} $$$", ffmt::a(name, 10), ffmt::i(i, 15));
    c.io.se(&s);
    c.io.so(&s);
}

extern "C" {
    #[link_name = "atanh"]
    fn c_atanh(x: f64) -> f64;
}
fn atanh(x: f64) -> f64 {
    unsafe { c_atanh(x) }
}

pub fn ijkl(m: i32, n: i32) -> i32 {
    if m <= n {
        (n * (n - 1)) / 2 + m
    } else {
        (m * (m - 1)) / 2 + n
    }
}

// ---------------------------------------------------------------- phases

pub fn sinsix(rho: f64, lspin: i32, sinsqr: &mut f64, sin2ph: &mut f64, dphi: &mut f64, sinphi: &mut f64, cosphi: &mut f64) {
    let a = rho;
    let l = lspin;
    let c = cos(a);
    let s = sin(a);
    let (x, y, f);
    if l == 0 {
        *sinsqr = s * s;
        *sin2ph = 2.0 * c * s;
        *cosphi = c;
        *sinphi = s;
        *dphi = 1.0;
        return;
    } else if l == 1 {
        let x = a;
        let g = 1.0 + sq(x);
        let d = (s - c * x) / g;
        *sinsqr = (s - c * x) * d;
        *sin2ph = 2.0 * (c + s * x) * d;
        *dphi = x * x / g;
        let dd = sqrt(g);
        *cosphi = (c + s * x) / dd;
        *sinphi = (s - c * x) / dd;
        return;
    } else if l == 2 {
        let a2 = a * a;
        let x0 = 3.0;
        y = 3.0 - a2;
        f = x0 * y + 2.0 * a2 * x0;
        x = x0 * a;
    } else if l == 3 {
        let a2 = a * a;
        let x0 = 15.0 - a2;
        y = 15.0 - a2 * 6.0;
        f = x0 * y - 2.0 * a2 * y + 12.0 * a2 * x0;
        x = x0 * a;
    } else if l == 4 {
        let a2 = a * a;
        let a4 = a2 * a2;
        let x0 = 105.0 - 10.0 * a2;
        y = 105.0 - 45.0 * a2 + a4;
        f = x0 * y - 20.0 * a2 * y + a2 * x0 * (90.0 - 4.0 * a2);
        x = x0 * a;
    } else {
        let (mut p, mut dp, mut ss, mut ds, mut xx, mut ff, mut yy) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        genpsf(a, l, &mut p, &mut dp, &mut ss, &mut ds, &mut xx, &mut ff, &mut yy);
        x = xx;
        f = ff;
        y = yy;
    }
    let g = sq(y) + sq(x);
    let d = (s * y - c * x) / g;
    *sinsqr = (s * y - c * x) * d;
    *sin2ph = 2.0 * (c * y + s * x) * d;
    *dphi = 1.0 - f / g;
    let dd = sqrt(g);
    *cosphi = (c * y + s * x) / dd;
    *sinphi = (s * y - c * x) / dd;
}

pub fn pgh(rho: f64, l: i32, bound: f64, g: &mut f64, h: &mut f64, p: &mut f64, dp: &mut f64, ds: &mut f64, ishift: i32, iffy: &mut i32) {
    let zero = 0.0f64;
    *iffy = 0;
    let r = rho;
    let mut gg = 0.0;
    *dp = 0.0;
    *ds = 0.0;
    let mut hh;
    if l == 0 {
        hh = r;
        *dp = 1.0;
        if ishift > 0 {
            gg = -bound;
        }
    } else if l == 1 {
        let r2 = r * r;
        let d = 1.0 + r2;
        hh = r * r2 / d;
        *dp = (r2 / d) * (3.0 + r2) / d;
        if ishift > 0 {
            gg = -1.0 / d - bound;
            *ds = 2.0 * (r / d) / d;
        }
    } else if l == 2 {
        let r2 = r * r;
        let r4 = r2 * r2;
        let d = 9.0 + r2 * (3.0 + r2);
        hh = r * r4 / d;
        *dp = r2 * (r2 / d) * (45.0 + r2 * (9.0 + r2)) / d;
        if ishift > 0 {
            *ds = (18.0 + 3.0 * r2) / d;
            gg = -*ds - bound;
            *ds = -6.0 * r / d + *ds * 4.0 * r * (1.5e0 + r2) / d;
        }
    } else if l == 3 {
        let r2 = r * r;
        let r4 = r2 * r2;
        let r6 = r4 * r2;
        let d = 225.0 + r2 * (45.0 + r2 * (6.0 + r2));
        hh = r * r6 / d;
        *dp = (r6 / d) * (1575.0 + (225.0 + r2 * (18.0 + r2)) * r2) / d;
        if ishift > 0 {
            *ds = (675.0 + r2 * (90.0 + r2 * 6.0)) / d;
            gg = -*ds - bound;
            *ds = -24.0 * (7.5e0 + r2) * r / d + *ds * 6.0 * r * (15.0 + r2 * (4.0 + r2)) / d;
        }
    } else if l == 4 {
        let r2 = r * r;
        let r4 = r2 * r2;
        let r8 = r4 * r4;
        let d = 11025.0 + r2 * (1575.0 + r2 * (135.0 + r2 * (10.0 + r2)));
        hh = r * r8 / d;
        *dp = (r8 / d) * (9.0 - (1575.0 + r2 * (270.0 + r2 * (30.0 + r2 * 4.0))) * 2.0 * (r2 / d));
        if ishift > 0 {
            *ds = (44100.0 + r2 * (4725.0 + r2 * (270.0 + 10.0 * r2))) / d;
            gg = -*ds - bound;
            *ds = -60.0 * r * (157.5e0 + r2 * (18.0 + r2)) / d + *ds * r * (3150.0 + r2 * (540.0 + r2 * (60.0 + r2 * 8.0))) / d;
        }
    } else {
        let (mut s, mut x, mut db, mut y) = (0.0, 0.0, 0.0, 0.0);
        hh = 0.0;
        genpsf(rho, l, &mut hh, dp, &mut s, ds, &mut x, &mut db, &mut y);
        if ishift > 0 {
            gg = s - bound;
        }
    }
    if hh <= 1.0e-35 {
        hh = zero;
    }
    if gg == zero && hh == zero {
        *g = 0.0;
        *h = 0.0;
        *p = 0.0;
        *ds = 0.0;
        *iffy = 1;
    } else if gg == zero {
        *g = 0.0;
        *h = -1.0 / hh;
        *p = hh;
    } else if hh == zero {
        *g = 1.0 / gg;
        *h = 0.0;
        *p = 0.0;
    } else if hh + gg != hh {
        if hh + gg != gg {
            let d = sq(hh) + sq(gg);
            *g = gg / d;
            *h = -hh / d;
            *p = hh;
        } else {
            *g = 1.0 / gg;
            *h = -(hh / gg) / gg;
            *p = hh;
        }
    } else {
        *g = (gg / hh) / hh;
        *h = -1.0 / hh;
        *p = hh;
    }
}

pub fn genpsf(rho: f64, ll: i32, pp: &mut f64, dp: &mut f64, ss: &mut f64, ds: &mut f64, qq: &mut f64, dq: &mut f64, yy: &mut f64) {
    let lwant = ll;
    let a = rho;
    let a2 = a * a;
    let a4 = a2 * a2;
    let a8 = a4 * a4;
    let mut l = 4;
    let mut yyy = 105.0 - 45.0 * a2 + a4;
    let mut qqq = 105.0 - 10.0 * a2;
    let mut dqq = (qqq - 20.0 * a2) * yyy + a2 * qqq * (90.0 - 4.0 * a2);
    qqq = qqq * a;
    let mut d = 11025.0 + a2 * (1575.0 + a2 * (135.0 + a2 * (10.0 + a2)));
    let mut ppp = a * a8 / d;
    let mut dpp = (a8 / d) * (9.0 - (1575.0 + a2 * (270.0 + a2 * (30.0 + a2 * 4.0))) * 2.0 * (a2 / d));
    // note: the Fortran uses dss here before it is set
    let mut dss = 0.0f64;
    let mut sss = -dss;
    dss = -60.0 * a * (157.5e0 + a2 * (18.0 + a2)) / d + dss * a * (3150.0 + a2 * (540.0 + a2 * (60.0 + a2 * 8.0))) / d;
    loop {
        l += 1;
        let el = l as f64;
        let els = el - sss;
        let x = ppp / els;
        let dx = (dpp + x * dss) / els;
        let yyx = yyy - qqq * x;
        let qqx = qqq + x * yyy;
        d = qqx * (qqq * dx * yyy + dqq * x);
        let dqx = (dqq + dx * sq(yyy)) * yyx + d;
        d = sq(els) + sq(ppp);
        let ppx = a2 * ppp / d;
        let c = 1.0 + a * (els * dss - ppp * dpp) / d;
        let mut dpx = 2.0 * ppp * c + a * dpp;
        dpx = a * dpx / d;
        let ssx = a2 * els / d - el;
        let mut dsx = 2.0 * els * c - a * dss;
        dsx = a * dsx / d;
        if l != lwant {
            ppp = ppx;
            dpp = dpx;
            sss = ssx;
            dss = dsx;
            qqq = qqx;
            dqq = dqx;
            yyy = yyx;
        } else {
            *pp = ppx;
            *dp = dpx;
            *ss = ssx;
            *ds = dsx;
            *qq = qqx;
            *dq = dqx;
            *yy = yyx;
            return;
        }
    }
}

pub fn pf(rho: f64, l: i32) -> f64 {
    let a = rho;
    if l == 0 {
        a
    } else if l == 1 {
        let a2 = a * a;
        a * a2 / (1.0 + a2)
    } else if l == 2 {
        let a2 = a * a;
        let a4 = a2 * a2;
        a * a4 / (9.0 + a2 * (3.0 + a2))
    } else if l == 3 {
        let a2 = a * a;
        let a4 = a2 * a2;
        let a6 = a2 * a4;
        let d = 225.0 + a2 * (45.0 + a2 * (6.0 + a2));
        a * a6 / d
    } else if l == 4 {
        let a2 = a * a;
        let a4 = a2 * a2;
        let a8 = a4 * a4;
        let d = 11025.0 + a2 * (1575.0 + a2 * (135.0 + a2 * (10.0 + a2)));
        a * a8 / d
    } else {
        let (mut p, mut dp, mut s, mut ds, mut qq, mut db, mut yy) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
        genpsf(rho, l, &mut p, &mut dp, &mut s, &mut ds, &mut qq, &mut db, &mut yy);
        p
    }
}

pub fn pghcou(
    c: &mut Ctx,
    rho: f64,
    l: i32,
    bound: f64,
    g: &mut f64,
    h: &mut f64,
    p: &mut f64,
    dp: &mut f64,
    ds: &mut f64,
    ishift: i32,
    iffy: &mut i32,
    eta: f64,
    sinphi: &mut f64,
    cosphi: &mut f64,
    dphi: &mut f64,
    jdopha: i32,
) {
    let zero = 0.0f64;
    *iffy = 0;
    let mut gg = 0.0;
    *dp = 0.0;
    *ds = 0.0;
    let jdoder = 1;
    let (mut ifail, mut pent, mut shift, mut dpent, mut dshift) = (0, 0.0, 0.0, 0.0, 0.0);
    pspcou(c, rho, l, eta, ishift, jdopha, jdoder, &mut ifail, &mut pent, &mut shift, &mut dpent, &mut dshift, sinphi, cosphi, dphi);
    let mut hh = pent;
    *dp = dpent;
    let s = shift;
    *ds = dshift;
    if ishift > 0 {
        gg = s - bound;
    } else {
        *ds = 0.0;
    }
    if hh <= 1.0e-35 {
        hh = 0.0;
    }
    *p = 0.0;
    *h = 0.0;
    *g = 0.0;
    if gg == zero && hh == zero {
        *dp = 0.0;
        *iffy = 1;
        return;
    }
    if hh == zero {
        *g = 1.0 / gg;
        return;
    }
    *p = hh;
    if gg == zero {
        *h = -1.0 / hh;
        return;
    }
    if hh + gg == hh {
        *g = (gg / hh) / hh;
        *h = -1.0 / hh;
        return;
    }
    if hh + gg == gg {
        *g = 1.0 / gg;
        *h = -(hh / gg) / gg;
        return;
    }
    let d = sq(hh) + sq(gg);
    *g = gg / d;
    *h = -hh / d;
}

/// Coulomb penetrability, shift and phase.
pub fn pspcou(
    c: &mut Ctx,
    rho: f64,
    lll: i32,
    eta: f64,
    ishift: i32,
    jdopha: i32,
    jdoder: i32,
    ifail: &mut i32,
    pencoul: &mut f64,
    shiftcoul: &mut f64,
    dpencoul: &mut f64,
    dshiftcoul: &mut f64,
    sinphi: &mut f64,
    cosphi: &mut f64,
    dphi: &mut f64,
) {
    *ifail = 0;
    *pencoul = 0.0;
    *shiftcoul = 0.0;
    *dpencoul = 0.0;
    *dshiftcoul = 0.0;
    *sinphi = 0.0;
    *cosphi = 1.0;
    *dphi = 0.0;
    if lll < 0 {
        *ifail = 10;
    }
    if *ifail > 0 {
        return;
    }
    let l = lll + 1;
    let llmin = 0;
    let mut llmax = lll + 2;
    if llmax > 100 {
        c.error("pspcou", "llmax larger than 100", " ");
    }
    let mut f = [0.0f64; 101];
    let mut fpr = [0.0f64; 101];
    let mut g = [0.0f64; 101];
    let mut gpr = [0.0f64; 101];
    let mut dummy = [0.0f64; 101];
    let (mut p, mut s, mut dp, mut sp, mut cp, mut dphix) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    if rho < 1.02e0 {
        coulx(c, eta, rho, lll, &mut llmax, &mut f, &mut fpr, &mut g, &mut gpr, &mut dummy, pencoul, shiftcoul, dpencoul, sinphi, cosphi, dphi, jdopha, jdoder, ishift);
        if ishift > 0 && jdoder > 0 {
            let rho2 = 1.01e0 * rho;
            let eta2 = eta * rho / rho2;
            coulx(c, eta2, rho2, lll, &mut llmax, &mut f, &mut fpr, &mut g, &mut gpr, &mut dummy, &mut p, &mut s, &mut dp, &mut sp, &mut cp, &mut dphix, jdopha, jdoder, ishift);
            *dshiftcoul = (s - *shiftcoul) / (rho2 - rho);
        }
    } else {
        let mut iexp = 0;
        coulfg(c, rho, eta, lll, llmax, &mut f, &mut g, &mut fpr, &mut gpr, pencoul, shiftcoul, dpencoul, sinphi, cosphi, dphi, jdopha, jdoder, ishift, ifail, &mut iexp);
        if *ifail != 0 {
            if ishift > 0 && jdoder > 0 {
                let rho2 = 1.01e0 * rho;
                let eta2 = eta * rho / rho2;
                coulfg(c, rho2, eta2, lll, llmax, &mut f, &mut g, &mut fpr, &mut gpr, &mut p, &mut s, &mut dp, &mut sp, &mut cp, &mut dphix, jdopha, jdoder, ishift, ifail, &mut iexp);
                if *ifail > 0 {
                    let msg = format!(
                        " Ifail,Rho,Eta,L,F,G={}{}{}{}{}{}",
                        ffmt::i(*ifail, 3),
                        ffmt::f(rho, 8, 5),
                        ffmt::f(eta, 8, 5),
                        ffmt::i(llmin, 5),
                        ffmt::e1p(f[l as usize], 12, 4),
                        ffmt::e1p(g[l as usize], 12, 4)
                    );
                    c.io.se(&msg);
                    c.io.so(&msg);
                    return;
                }
                *dshiftcoul = (s - *shiftcoul) / (rho2 - rho);
            }
        }
        // the Fortran's "else if (ifail.gt.0)" branch cannot be reached
    }
}

#[allow(clippy::too_many_arguments)]
pub fn coulfg(
    c: &mut Ctx,
    xx: f64,
    eta1: f64,
    lll: i32,
    llmax: i32,
    fc: &mut [f64; 101],
    gc: &mut [f64; 101],
    fcp: &mut [f64; 101],
    gcp: &mut [f64; 101],
    pcoul: &mut f64,
    scoul: &mut f64,
    dpcoul: &mut f64,
    sinphi: &mut f64,
    cosphi: &mut f64,
    dphi: &mut f64,
    jdopha: i32,
    jdoder: i32,
    ishift: i32,
    ifail: &mut i32,
    iexp: &mut i32,
) {
    let zero = 0.0f64;
    let one = 1.0f64;
    let two = 2.0f64;
    let ten = 10.0f64;
    let ten2 = 100.0f64;
    let half = 0.5e0f64;
    let accur = 1.0e-16f64;
    let tm30 = 1.0e-30f64;
    let abort = 2.0e4f64;
    // arrays are 1-based in the Fortran: fc(1) is fc[1]
    *ifail = 0;
    *iexp = 1;
    let eta = eta1;
    let mut gjwkb = 0.0;
    let mut fjwkb = 0.0;
    let mut paccq = 1.0;
    let acc = accur;
    let acc4 = acc * ten2 * ten2 / ten;
    let acch = sqrt(acc);
    if xx <= acch {
        *ifail = -1;
        c.io.se(&format!(
            " FOR Xx = {} TRY SMALL-X SOLUTIONS OR X NEGATIVE\n SQUARE ROOT ACCURACY PARAMETER = {}\n",
            ffmt::e1p(xx, 12, 3),
            ffmt::e1p(acch, 12, 3)
        ));
        return;
    }
    let x = xx;
    let xlm = 0.0f64;
    let e2mm1 = eta * eta + xlm * xlm + xlm;
    let mut xlturn = x * (x - 2.0 * eta) < xlm * xlm + xlm;
    let xll = llmax as f64;
    let l1 = (llmax + 1) as usize;
    let xi = 1.0 / x;
    let mut fcl = 1.0;
    let mut pk = xll + 1.0;
    let px = pk + abort;
    let mut ek;
    let mut f;
    let mut pk1;
    loop {
        // 10
        ek = eta / pk;
        f = (ek + pk * xi) * fcl + (fcl - 1.0) * xi;
        pk1 = pk + 1.0;
        if (eta * x + pk * pk1).abs() <= acc {
            fcl = (1.0 + ek * ek) / (1.0 + sq(eta / pk1));
            pk = 2.0 + pk;
            continue;
        }
        break;
    }
    let mut d = 1.0 / ((pk + pk1) * (xi + ek / pk1));
    let mut df = -fcl * (1.0 + ek * ek) * d;
    if fcl != one {
        fcl = -1.0;
    }
    if d < zero {
        fcl = -fcl;
    }
    f = f + df;
    let mut p = 1.0;
    let fail110 = |c: &mut Ctx, f: f64, df: f64, pk: f64, px: f64, ifail: &mut i32| {
        *ifail = 1;
        c.io.se(&format!(
            " CF1 HAS FAILED TO CONVERGE AFTER {} ITERATIONS\n F,Df,Pk,Px,Accur =  {}{}{}{}{}\n\n",
            ffmt::f(abort, 10, 0),
            ffmt::e1p(f, 12, 3),
            ffmt::e1p(df, 12, 3),
            ffmt::e1p(pk, 12, 3),
            ffmt::e1p(px, 12, 3),
            ffmt::e1p(acc, 12, 3)
        ));
    };
    loop {
        // 20
        pk = pk1;
        pk1 = pk1 + 1.0;
        ek = eta / pk;
        let tk = (pk + pk1) * (xi + ek / pk1);
        d = tk - d * (1.0 + ek * ek);
        if d.abs() <= acch {
            let msg = format!(
                "\n CF1 ACCURACY LOSS: D,Df,Acch,K,Eta/K,Eta,X = {}{}{}{}{}{}{}\n",
                ffmt::e1p(d, 9, 2),
                ffmt::e1p(df, 9, 2),
                ffmt::e1p(acch, 9, 2),
                ffmt::e1p(pk, 9, 2),
                ffmt::e1p(ek, 9, 2),
                ffmt::e1p(eta, 9, 2),
                ffmt::e1p(x, 9, 2)
            );
            c.io.se(&msg);
            c.io.so(&msg);
            p = p + 1.0;
            if p > two {
                fail110(c, f, df, pk, px, ifail);
                return;
            }
        }
        d = 1.0 / d;
        if d < zero {
            fcl = -fcl;
        }
        df = df * (d * tk - 1.0);
        f = f + df;
        if pk > px {
            fail110(c, f, df, pk, px, ifail);
            return;
        }
        if df.abs() >= f.abs() * acc {
            continue;
        }
        break;
    }
    let mut xl = 0.0;
    if llmax > 0 {
        fcl = fcl * tm30;
        let mut fpl = fcl * f;
        fcp[l1] = fpl;
        fc[l1] = fcl;
        xl = xll;
        for lp in 1..=llmax {
            let el = eta / xl;
            let rl = sqrt(1.0 + el * el);
            let sl = el + xl * xi;
            let l = l1 - lp as usize;
            let fcl1 = (fcl * sl + fpl) / rl;
            fpl = fcl1 * sl - fcl * rl;
            fcl = fcl1;
            fc[l] = fcl;
            fcp[l] = fpl;
            gc[l + 1] = rl;
            xl = xl - 1.0;
        }
        if fcl == zero {
            fcl = acc;
        }
        f = fpl / fcl;
    }
    if xlturn {
        jwkb(x, eta, xlm.max(zero), &mut fjwkb, &mut gjwkb, iexp);
    }
    let w;
    let gam;
    let q;
    if *iexp > 1 || gjwkb > 1.0 / (acch * ten2) {
        w = fjwkb;
        gam = gjwkb * w;
        p = f;
        q = 1.0;
    } else {
        xlturn = false;
        let ta = 2.0 * abort;
        pk = 0.0;
        let wi = eta + eta;
        p = 0.0;
        let mut qv = 1.0 - eta * xi;
        let mut ar = -e2mm1;
        let mut ai = eta;
        let br = 2.0 * (x - eta);
        let mut bi = 2.0;
        let mut dr = br / (br * br + bi * bi);
        let mut di = -bi / (br * br + bi * bi);
        let mut dp = -xi * (ar * di + ai * dr);
        let mut dq = xi * (ar * dr - ai * di);
        loop {
            p = p + dp;
            qv = qv + dq;
            pk = pk + 2.0;
            ar = ar + pk;
            ai = ai + wi;
            bi = bi + 2.0;
            let dd = ar * dr - ai * di + br;
            di = ai * dr + ar * di + bi;
            let cc = 1.0 / (dd * dd + di * di);
            dr = cc * dd;
            di = -cc * di;
            let a = br * dr - bi * di - 1.0;
            let b = bi * dr + br * di;
            let cc2 = dp * a - dq * b;
            dq = dp * b + dq * a;
            dp = cc2;
            if pk > ta {
                *ifail = 2;
                c.io.se(&format!(
                    " CF2 HAS FAILED TO CONVERGE AFTER {}ITERATIONS\n P,Q,Dp,Dq,Accur =  {}{}{}{}{}\n\n",
                    ffmt::f(abort, 7, 0),
                    ffmt::e1p(p, 17, 7),
                    ffmt::e1p(qv, 17, 7),
                    ffmt::e1p(dp, 17, 7),
                    ffmt::e1p(dq, 17, 7),
                    ffmt::e1p(acc, 12, 3)
                ));
                return;
            }
            if dp.abs() + dq.abs() >= (p.abs() + qv.abs()) * acc {
                continue;
            }
            break;
        }
        paccq = half * acc;
        if qv.abs() < one {
            paccq = half * acc / qv.abs();
        }
        if p.abs() > qv.abs() {
            paccq = paccq * p.abs();
        }
        if qv <= acc4 * p.abs() {
            *ifail = 3;
            let dpacc = p.abs() * acc * 10000.0;
            c.io.se(&format!(
                " Final Q.le.abs(P)*Acc*10**4\n       Q, ... , P, Acc = {}{}{}{}\n       Llmax = {}\n",
                ffmt::e1p(qv, 12, 3),
                ffmt::e1p(dpacc, 12, 3),
                ffmt::e1p(p, 12, 3),
                ffmt::e1p(acc, 12, 3),
                ffmt::i(llmax, 5)
            ));
            return;
        }
        gam = (f - p) / qv;
        w = 1.0 / sqrt((f - p) * gam + qv);
        q = qv;
    }
    let fcm = crate::fmath::sign(w, fcl);
    fc[1] = fcm;
    let mut gcl = if !xlturn { fcm * gam } else { gjwkb };
    gc[1] = gcl;
    let mut gpl = gcl * (p - q / gam);
    gcp[1] = gpl;
    fcp[1] = fcm * f;
    let w = w / fcl.abs();
    for l in 1..=llmax as usize {
        xl = xl + 1.0;
        let el = eta / xl;
        let rl = gc[l + 1];
        let sl = el + xl * xi;
        let gcl1 = (sl * gcl - gpl) / rl;
        gpl = rl * gcl - sl * gcl1;
        gcl = gcl1;
        gc[l + 1] = gcl1;
        gcp[l + 1] = gpl;
        fcp[l + 1] = w * fcp[l + 1];
        fc[l + 1] = w * fc[l + 1];
    }
    let l = (lll + 1) as usize;
    if *iexp > 1 {
        if *iexp < 150 {
            let asq = sq(gc[l]);
            let aaa = powi(ten, -*iexp * 2);
            *pcoul = xx / asq * aaa;
            if jdopha > 0 {
                *sinphi = fc[l] / gc[l] * aaa;
                *cosphi = 1.0 - sq(*sinphi);
                if jdoder > 0 {
                    *dphi = gcp[l] / gc[l] * *sinphi - fcp[l] / gc[l] * aaa;
                }
            }
        }
        if ishift > 0 {
            let sss = xx * gcp[l] / gc[l];
            *scoul = sss;
        }
    } else {
        let asq = sq(fc[l]) + sq(gc[l]);
        *pcoul = xx / asq;
        let sss = xx * (fc[l] * fcp[l] + gc[l] * gcp[l]) / asq;
        if ishift > 0 {
            *scoul = sss;
        }
        if jdoder > 0 {
            *dpcoul = (1.0 - 2.0 * sss) / asq;
        }
        if jdopha > 0 {
            let a = sqrt(asq);
            *sinphi = fc[l] / a;
            *cosphi = gc[l] / a;
            if jdoder > 0 {
                *dphi = (gcp[l] * fc[l] - gc[l] * fcp[l]) / asq;
            }
        }
    }
    let _ = paccq;
}

pub fn jwkb(xx: f64, eta1: f64, xl: f64, fjwkb: &mut f64, gjwkb: &mut f64, iexp: &mut i32) {
    let zero = 0.0f64;
    let half = 0.5e0f64;
    let ten = 10.0f64;
    let aloge = 0.434294481903251816667932e0f64;
    let six35 = 0.1714285714285714285714285714285714285714285e0f64;
    let x = xx;
    let eta = eta1;
    let gh2 = x * (eta + eta - x);
    let xll1 = (xl * xl + xl).max(zero);
    if gh2 + xll1 <= zero {
        return;
    }
    let hll = xll1 + six35;
    let hl = sqrt(hll);
    let sl = eta / hl + hl / x;
    let rl2 = 1.0 + eta * eta / hll;
    let gh = sqrt(gh2 + hll) / x;
    let mut phi = x * gh - half * (hl * log(sq(gh + sl) / rl2) - log(gh));
    if eta != zero {
        phi = phi - eta * atan2(x * gh, x - eta);
    }
    let phi10 = -phi * aloge;
    *iexp = crate::fmath::int(phi10);
    if *iexp > 70 {
        *gjwkb = pow(ten, phi10 - *iexp as f64);
    } else {
        *gjwkb = exp(-phi);
        *iexp = 0;
    }
    *fjwkb = half / (gh * *gjwkb);
}

#[allow(clippy::too_many_arguments)]
pub fn coulx(
    c: &mut Ctx,
    eeta: f64,
    rrho: f64,
    lll: i32,
    llmax: &mut i32,
    f: &mut [f64; 101],
    fpr: &mut [f64; 101],
    g: &mut [f64; 101],
    gpr: &mut [f64; 101],
    sigma: &mut [f64; 101],
    p: &mut f64,
    s: &mut f64,
    dp: &mut f64,
    sinphi: &mut f64,
    cosphi: &mut f64,
    dphi: &mut f64,
    jdopha: i32,
    jdoder: i32,
    ishift: i32,
) {
    let five = 5.0f64;
    let eta = eeta;
    let rho = rrho;
    let lmax = *llmax;
    if eta > 10.0 * rho && eta > five {
        bigeta(c, eta, rho, lll, llmax, f, fpr, g, gpr, p, s, dp, sinphi, cosphi, dphi, jdopha, jdoder, ishift);
    } else {
        let mut rhoi = 0.0;
        let mut u = 0.0;
        let mut upr = 0.0;
        if eta >= five {
            asymp1(eta, &mut rhoi, &mut u, &mut upr);
            taylor(eta, rho, &mut u, &mut upr, &mut rhoi);
        } else {
            let mut sigma0 = 0.0;
            xsigll(eta, sigma, &mut sigma0, lmax);
            asymp2(eta, rho, &mut u, &mut upr, &mut rhoi, sigma0);
            if u.abs() <= 1.0e25 {
                taylor(eta, rho, &mut u, &mut upr, &mut rhoi);
            }
        }
        let g0 = u;
        let g0pr = upr;
        if g0.abs() > 1.0e25 {
            end1(llmax, f, fpr, g, gpr, g0, g0pr);
        } else {
            getfg(eta, rho, llmax, lll, f, fpr, g, gpr, g0, g0pr);
        }
        getps(rho, lll, f, fpr, g, gpr, p, s, dp, sinphi, cosphi, dphi, jdopha, jdoder, ishift);
    }
}

pub fn asymp1(eeta: f64, rhoi: &mut f64, u: &mut f64, upr: &mut f64) {
    let eta = eeta;
    let ceta = pow(eta, 0.3333333333333333e0);
    let seta = sqrt(ceta);
    let temp = 1.0 / sq(ceta);
    *u = 1.223404016e0
        * seta
        * (1.0
            + sq(temp)
                * (0.04959570165e0
                    + temp * (-0.008888888889e0 + sq(temp) * (0.002455199181e0 + temp * (-0.0009108958061e0 + sq(temp) * 0.0002534684115e0)))));
    *upr = -0.7078817734e0
        * (1.0
            + temp
                * (-0.1728260369e0
                    + sq(temp) * (0.0003174603174e0 + temp * (-0.003581214850e0 + sq(temp) * (0.0003117824680e0 - temp * 0.0009073966427e0)))))
        / seta;
    *rhoi = 2.0 * eta;
}

pub fn xsigll(eeta: f64, sigma: &mut [f64; 101], sigma0: &mut f64, lmax: i32) {
    let small = 0.000001e0f64;
    let ber = [
        0.0,
        0.1666666666666666666666666666666666667e0,
        -0.0333333333333333333333333333333333333e0,
        0.0238095238095238095238095238095238095e0,
        -0.0333333333333333333333333333333333333e0,
        0.0757575757575757575757575757575757576e0,
    ];
    let zero = 0.0f64;
    let mmmxxx = 100000;
    let eta = eeta;
    let peta = eta.abs();
    if peta >= 3.0e0 {
        let mut sum = 0.0;
        for i in 1..=5 {
            let xi = i as f64;
            let m = 2 * i - 1;
            let xm = m as f64;
            sum = sum + ber[i as usize] / (2.0 * xi * xm * powi(peta, m));
        }
        *sigma0 = PI / 4.0 + peta * (log(peta) - 1.0) - sum;
    } else {
        let mut sumas = 0.0;
        for is in 1..=mmmxxx {
            let s = is as f64;
            let temp1 = peta / s;
            let mut as_;
            if s <= 2.0 * peta {
                as_ = temp1 - atan(temp1);
            } else {
                as_ = 0.0;
                let mut k = 0;
                for j in 1..=mmmxxx {
                    let m = j + j + 1;
                    let xm = m as f64;
                    let add = powi(temp1, m) / xm;
                    if k == 0 {
                        as_ = as_ + add;
                        k = 1;
                    } else {
                        as_ = as_ - add;
                        k = 0;
                    }
                    if (add / as_).abs() <= small {
                        break;
                    }
                }
            }
            sumas = sumas + as_;
            if (as_ / sumas).abs() <= small {
                break;
            }
        }
        *sigma0 = -EULER * peta + sumas;
    }
    if eta < zero {
        *sigma0 = -*sigma0;
    }
    sigma[1] = *sigma0;
    if lmax > 0 {
        for ll in 1..=lmax as usize {
            let xl = ll as f64;
            sigma[ll + 1] = sigma[ll] + atan(eta / xl);
        }
    }
}

pub fn asymp2(eeta: f64, rrho: f64, u: &mut f64, upr: &mut f64, rhoi: &mut f64, sigma0: f64) {
    let del = 100.0f64;
    let epslon = 0.000001e0f64;
    let ten = 10.0f64;
    let eta = eeta;
    let rho = rrho;
    *rhoi = (rho * 2.0).max(ten).max(ten * eta);
    let mut z = [0.0f64; 5];
    'l10: loop {
        let mut xn = 0.0f64;
        let mut zold = [0.0f64; 5];
        let mut znew = [0.0f64; 5];
        let mut bigz = [0.0f64; 5];
        zold[1] = 1.0;
        zold[2] = 0.0;
        zold[3] = 0.0;
        zold[4] = 1.0 - eta / *rhoi;
        for i in 1..=4 {
            z[i] = zold[i];
            bigz[i] = z[i].abs();
        }
        let mut jcheck = 0;
        for _n in 1..=100 {
            let temp = 2.0 * (xn + 1.0) * *rhoi;
            let an = (2.0 * xn + 1.0) * eta / temp;
            let bn = (eta * eta - xn * (xn + 1.0)) / temp;
            xn = xn + 1.0;
            znew[1] = an * zold[1] - bn * zold[2];
            znew[2] = an * zold[2] + bn * zold[1];
            znew[3] = an * zold[3] - bn * zold[4] - znew[1] / *rhoi;
            znew[4] = an * zold[4] + bn * zold[3] - znew[2] / *rhoi;
            let mut icheck = 0;
            for i in 1..=4 {
                z[i] = z[i] + znew[i];
                zold[i] = znew[i];
                let temp2 = z[i].abs();
                bigz[i] = bigz[i].max(z[i].abs());
                if bigz[i] / temp2 > del {
                    *rhoi = *rhoi * 2.0;
                    continue 'l10;
                }
                if (znew[i] / z[i]).abs() <= epslon {
                    icheck += 1;
                }
            }
            let w = z[1] * z[4] - z[2] * z[3];
            if w.abs() > ten {
                *rhoi = *rhoi * 2.0;
                continue 'l10;
            }
            if icheck == 4 {
                jcheck += 1;
                if jcheck >= 4 {
                    break 'l10;
                }
            } else {
                jcheck = 0;
            }
        }
        // 50
        *rhoi = *rhoi * 2.0;
    }
    // 60
    let phi = *rhoi - eta * log(2.0 * *rhoi) + sigma0;
    let cosphi = cos(phi);
    let sinphi = sin(phi);
    let g0 = z[1] * cosphi - z[2] * sinphi;
    let g0pr = z[3] * cosphi - z[4] * sinphi;
    *u = g0;
    *upr = g0pr;
}

pub fn taylor(eeta: f64, rrho: f64, u: &mut f64, upr: &mut f64, rhoi: &mut f64) {
    let epslon = 1.0e-6f64;
    let bigger = 1.0e10f64;
    let biggst = 1.0e30f64;
    let del = 100.0f64;
    let zero = 0.0f64;
    let eta = eeta;
    let rho = rrho;
    let mut delta = rho - *rhoi;
    if delta == zero {
        return;
    }
    let mut a = [0.0f64; 101];
    'l10: loop {
        a[1] = *u;
        a[2] = delta * *upr;
        a[3] = -delta * delta / 2.0 * (1.0 - 2.0 * eta / *rhoi) * a[1];
        let mut nstart: i32 = 4;
        let (sum, sumpr);
        'l20: loop {
            let mut jcheck = 0;
            let mut s = 0.0f64;
            let mut spr = 0.0f64;
            let mut big = 0.0f64;
            let mut bigpr = 0.0f64;
            let mut n = 1;
            let mut go40 = false;
            let mut go60 = false;
            while n <= 100 {
                let xn = (n - 1) as f64;
                let nu = n as usize;
                if n >= nstart {
                    a[nu] = -(delta * (xn - 1.0) * (xn - 2.0) * a[nu - 1]
                        + (*rhoi - 2.0 * eta) * sq(delta) * a[nu - 2]
                        + (delta * delta * delta) * a[nu - 3])
                        / (*rhoi * (xn - 1.0) * xn);
                    if a[nu] > bigger {
                        go40 = true;
                        break;
                    }
                }
                s = s + a[nu];
                spr = spr + xn * a[nu];
                if s >= biggst || spr >= biggst {
                    go40 = true;
                    break;
                }
                big = big.max(s.abs());
                bigpr = bigpr.max(spr.abs());
                if s == zero || spr == zero {
                    jcheck = 0;
                } else {
                    if (big / s).abs() >= del || (bigpr / spr).abs() >= del {
                        go40 = true;
                        break;
                    }
                    if (a[nu] / s).abs() >= epslon || (xn * a[nu] / spr).abs() >= epslon {
                        jcheck = 0;
                    } else {
                        jcheck += 1;
                        if jcheck >= 4 {
                            go60 = true;
                            break;
                        }
                    }
                }
                n += 1;
            }
            if go60 {
                sum = s;
                sumpr = spr;
                break 'l20;
            }
            if !go40 {
                n = 100;
            }
            // 40
            nstart = nstart.max(n + 1);
            let m = nstart - 1;
            delta = delta / 2.0;
            let mut temp = 2.0f64;
            for nn in 1..=m as usize {
                temp = temp / 2.0;
                a[nn] = a[nn] * temp;
            }
        }
        // 60
        *u = sum;
        *upr = sumpr / delta;
        *rhoi = *rhoi + delta;
        delta = rho - *rhoi;
        if delta.abs() >= epslon {
            continue 'l10;
        }
        return;
    }
}

pub fn end1(llmax: &mut i32, f: &mut [f64; 101], fpr: &mut [f64; 101], g: &mut [f64; 101], gpr: &mut [f64; 101], g0: f64, g0pr: f64) {
    *llmax = 0;
    g[1] = g0;
    gpr[1] = g0pr;
    f[1] = 0.0;
    fpr[1] = 0.0;
}

pub fn getfg(eta: f64, rho: f64, llmax: &mut i32, lll: i32, f: &mut [f64; 101], fpr: &mut [f64; 101], g: &mut [f64; 101], gpr: &mut [f64; 101], g0: f64, g0pr: f64) {
    let big = 1.0e12f64;
    let lmax = *llmax;
    let mut limit = 3.max(lmax + 1);
    g[1] = g0;
    gpr[1] = g0pr;
    g[2] = ((eta + 1.0 / rho) * g[1] - gpr[1]) / sqrt(sq(eta) + 1.0);
    for l in 3..=limit {
        let lu = l as usize;
        let xl = (l - 1) as f64;
        let temp1 = sqrt(xl * xl + eta * eta);
        g[lu] = (2.0 * xl - 1.0) / temp1 * (eta / (xl - 1.0) + xl / rho) * g[lu - 1]
            - xl / temp1 * sqrt(1.0 + sq(eta / (xl - 1.0))) * g[lu - 2];
        if g[lu].abs() > big && l > lll {
            limit = l;
            break;
        }
    }
    // the recursion below runs far past the arrays' 100 elements in the
    // Fortran only through scalars, so plain scalars are used here too
    let mut gm2 = g[(limit - 1) as usize];
    let mut gm1 = g[limit as usize];
    let mut il = -1;
    let mut j = limit;
    let mut gm = 0.0;
    while j <= 10000 {
        let xl = j as f64;
        let temp1 = sqrt(xl * xl + eta * eta);
        gm = (2.0 * xl - 1.0) / temp1 * (eta / (xl - 1.0) + xl / rho) * gm1 - xl / temp1 * sqrt(1.0 + sq(eta / (xl - 1.0))) * gm2;
        if (g[limit as usize] / gm).abs() > 1.0e-4 {
            il = -2;
        }
        if il > 0 {
            break;
        }
        il += 1;
        gm2 = gm1;
        gm1 = gm;
        j += 1;
    }
    let mut xl = j as f64;
    let mut fp1 = xl / gm / sqrt(sq(xl) + sq(eta));
    let mut fp2 = 0.0;
    let mut l = j - 1;
    for _ll in 1..=(j - 3 - limit) {
        l -= 1;
        xl = l as f64;
        let temp2 = sqrt(sq(xl + 1.0) + sq(eta));
        let fp = ((2.0 * xl + 3.0) * (eta / (xl + 2.0) + (xl + 1.0) / rho) * fp1 - (xl + 1.0) * sqrt(1.0 + sq(eta / (xl + 2.0))) * fp2) / temp2;
        fp2 = fp1;
        fp1 = fp;
    }
    for _ll in 1..=limit + 2 {
        l -= 1;
        xl = l as f64;
        let temp2 = sqrt(sq(xl + 1.0) + sq(eta));
        let fp = ((2.0 * xl + 3.0) * (eta / (xl + 2.0) + (xl + 1.0) / rho) * fp1 - (xl + 1.0) * sqrt(1.0 + sq(eta / (xl + 2.0))) * fp2) / temp2;
        f[(l + 1) as usize] = fp;
        fp2 = fp1;
        fp1 = fp;
    }
    fpr[1] = (1.0 / rho + eta) * f[1] - sqrt(1.0 + sq(eta)) * f[2];
    for l in 2..=limit {
        let lu = l as usize;
        let xl = l as f64;
        fpr[lu] = (xl / rho + eta / xl) * f[lu] - sqrt(1.0 + sq(eta / xl)) * f[lu + 1];
        let temp1 = eta / (xl - 1.0);
        gpr[lu] = sqrt(1.0 + sq(temp1)) * g[lu - 1] - ((xl - 1.0) / rho + temp1) * g[lu];
    }
    if limit <= lmax {
        *llmax = limit - 1;
    }
}

#[allow(clippy::too_many_arguments)]
pub fn bigeta(
    c: &mut Ctx,
    eeta: f64,
    rrho: f64,
    lll: i32,
    llmax: &mut i32,
    f: &mut [f64; 101],
    fpr: &mut [f64; 101],
    g: &mut [f64; 101],
    gpr: &mut [f64; 101],
    p: &mut f64,
    s: &mut f64,
    dp: &mut f64,
    sinphi: &mut f64,
    cosphi: &mut f64,
    dphi: &mut f64,
    jdopha: i32,
    jdoder: i32,
    ishift: i32,
) {
    let half = 0.5e0f64;
    let eta = eeta;
    let rho = rrho;
    let q = 2.0 * rho * eta;
    let zhalf = sqrt(q);
    let z = zhalf * 2.0;
    // 1./(k+1) is a default (single precision) real division
    let inv1 = |k: i32| (1.0f32 / ((k + 1) as f32)) as f64;
    let ksq = |k: i32| ((k + 1) * (k + 1)) as f64;

    let mut sum = 1.0;
    let mut a = q;
    let mut ok = false;
    for k in 1..=100 {
        if sum + a == sum {
            ok = true;
            break;
        }
        sum = sum + a;
        a = a * q / ksq(k);
    }
    if !ok {
        c.error("bigeta", "I0 sum failed", " ");
    }
    let ai0 = sum;

    sum = -(log(zhalf) + EULER) * ai0;
    a = q;
    let mut b = 1.0;
    ok = false;
    for k in 1..=100 {
        if sum + a * b == sum {
            ok = true;
            break;
        }
        sum = sum + a * b;
        b = b + inv1(k);
        a = a * q / ksq(k);
    }
    if !ok {
        c.error("bigeta", "K0 sum failed", " ");
    }
    let ak0 = sum;

    sum = 1.0;
    a = q / 2.0;
    ok = false;
    for k in 1..=100 {
        if sum + a == sum {
            ok = true;
            break;
        }
        sum = sum + a;
        a = a * q / ((k + 1) * (k + 2)) as f64;
    }
    if !ok {
        c.error("bigeta", "L1 sum failed", " ");
    }
    let ai1 = sum * zhalf;

    sum = 1.0 / z + (log(zhalf) + EULER) * ai1 - zhalf * half;
    a = zhalf * q / 2.0;
    b = 1.0;
    let mut cc = half * half;
    ok = false;
    for k in 1..=100 {
        if sum - a * (b + cc) == sum {
            ok = true;
            break;
        }
        sum = sum - a * (b + cc);
        cc = half / (k + 2) as f64;
        b = b + inv1(k);
        a = a * q / ((k + 1) * (k + 2)) as f64;
    }
    if !ok {
        c.error("bigeta", "K1 sum failed", " ");
    }
    let ak1 = sum;

    let mut c1 = sqrt(PI * rho);
    let mut d = sqrt(2.0 * PI * eta);
    f[1] = c1 * ai1;
    fpr[1] = d * ai0;
    c1 = 2.0 * c1 / PI;
    d = 2.0 * d / PI;
    g[1] = c1 * ak1;
    gpr[1] = -d * ak0;
    let g0 = g[1];
    let g0pr = gpr[1];
    if lll > 0 {
        getfg(eta, rho, llmax, lll, f, fpr, g, gpr, g0, g0pr);
    }
    let a = PI * eta;
    let b = exp(-a);
    let n = (lll + 1) as usize;
    if f[n] * b * b + g[n] == g[n] {
        *p = (rho * b / g[n]) / g[n] * b;
        if ishift > 0 {
            *s = rho * gpr[n] / g[n];
        }
        if jdopha != 0 {
            *sinphi = f[n] / g[n] * b * b;
        }
    } else {
        f[n] = f[n] * b;
        g[n] = g[n] / b;
        fpr[n] = fpr[n] * b;
        gpr[n] = gpr[n] / b;
        getps(rho, lll, f, fpr, g, gpr, p, s, dp, sinphi, cosphi, dphi, jdopha, jdoder, ishift);
    }
}

#[allow(clippy::too_many_arguments)]
pub fn getps(
    rho: f64,
    lll: i32,
    f: &[f64; 101],
    fpr: &[f64; 101],
    g: &[f64; 101],
    gpr: &[f64; 101],
    p: &mut f64,
    s: &mut f64,
    dp: &mut f64,
    sinphi: &mut f64,
    cosphi: &mut f64,
    dphi: &mut f64,
    jdopha: i32,
    jdoder: i32,
    ishift: i32,
) {
    let n = (lll + 1) as usize;
    let asq = sq(f[n]) + sq(g[n]);
    let a = sqrt(asq);
    *p = rho / asq;
    let mut ss = 0.0;
    if jdoder != 0 || ishift > 0 {
        ss = rho * (f[n] * fpr[n] + g[n] * gpr[n]) / asq;
    }
    if jdoder > 0 {
        *dp = (1.0 - 2.0 * ss) / asq;
    }
    if ishift > 0 {
        *s = ss;
    }
    if jdopha > 0 {
        *sinphi = f[n] / a;
        *cosphi = g[n] / a;
        if jdoder > 0 {
            *dphi = (gpr[n] * f[n] - g[n] * fpr[n]) / asq;
        }
    }
}

// ---------------------------------------------------- matrix inversion
// rmat/rinv here are the (2,ntriag,nerm) arrays; `r(i,k)` is rmat(i,k,ier).

fn onech(rmat: &FA, rinv: &mut FA, ier: i32) {
    let zero = 0.0f64;
    let r1 = rmat[(1, 1, ier)];
    let r2 = rmat[(2, 1, ier)];
    if r1 + r2 == r1 {
        rinv[(1, 1, ier)] = 1.0 / r1;
        rinv[(2, 1, ier)] = -(r2 / r1) / r1;
    } else if r1 + r2 == r2 {
        rinv[(1, 1, ier)] = (r1 / r2) / r2;
        rinv[(2, 1, ier)] = -1.0 / r2;
    } else if r1 == zero {
        rinv[(1, 1, ier)] = 0.0;
        rinv[(2, 1, ier)] = -1.0 / r2;
    } else {
        let aa = sq(r1) + sq(r2);
        rinv[(1, 1, ier)] = r1 / aa;
        rinv[(2, 1, ier)] = -r2 / aa;
    }
}

fn twoch(rmat: &mut FA, rinv: &mut FA, ier: i32) {
    let zero = 0.0f64;
    let (k11, k21, k22) = (1, 2, 3);
    if rmat[(1, k11, ier)] != zero || rmat[(1, k21, ier)] != zero || rmat[(1, k22, ier)] != zero {
        let mut a = 0.0;
        let mut k = 0;
        for i in 1..=3 {
            if rmat[(1, i, ier)].abs() > a {
                k = i;
            }
            if rmat[(1, i, ier)].abs() > a {
                a = rmat[(1, i, ier)].abs();
            }
            if rmat[(2, i, ier)].abs() > a {
                k = i;
            }
            if rmat[(2, i, ier)].abs() > a {
                a = rmat[(2, i, ier)].abs();
            }
        }
        if a > zero {
            if k == 2 {
                for i in 1..=3 {
                    rmat[(1, i, ier)] = rmat[(1, i, ier)] / a;
                    rmat[(2, i, ier)] = rmat[(2, i, ier)] / a;
                }
            } else {
                rmat[(1, k, ier)] = rmat[(1, k, ier)] / a;
                rmat[(2, k, ier)] = rmat[(2, k, ier)] / a;
                rmat[(1, 2, ier)] = rmat[(1, 2, ier)] / sqrt(a);
                rmat[(2, 2, ier)] = rmat[(2, 2, ier)] / sqrt(a);
            }
        }
        let r = |i: i32, j: i32| rmat[(i, j, ier)];
        let bbr = r(1, k11) * r(1, k22) - r(2, k11) * r(2, k22) - sq(r(1, k21)) + sq(r(2, k21));
        let bbi = r(1, k11) * r(2, k22) + r(2, k11) * r(1, k22) - 2.0 * r(1, k21) * r(2, k21);
        let (aar, aai);
        if bbr + bbi != bbi {
            if bbr + bbi != bbr {
                let aa = 1.0 / (bbr * bbr + bbi * bbi);
                aar = bbr * aa;
                aai = -bbi * aa;
            } else {
                aar = 1.0 / bbr;
                aai = -(bbi / bbr) / bbr;
            }
        } else {
            aar = (bbr / bbi) / bbi;
            aai = -1.0 / bbi;
        }
        rinv[(1, k11, ier)] = aar * r(1, k22) - aai * r(2, k22);
        rinv[(2, k11, ier)] = aar * r(2, k22) + aai * r(1, k22);
        rinv[(1, k21, ier)] = -aar * r(1, k21) + aai * r(2, k21);
        rinv[(2, k21, ier)] = -aar * r(2, k21) - aai * r(1, k21);
        rinv[(1, k22, ier)] = aar * r(1, k11) - aai * r(2, k11);
        rinv[(2, k22, ier)] = aar * r(2, k11) + aai * r(1, k11);
        if a != zero {
            if k == 2 {
                for i in 1..=3 {
                    rinv[(1, i, ier)] = rinv[(1, i, ier)] / a;
                    rinv[(2, i, ier)] = rinv[(2, i, ier)] / a;
                }
            } else {
                rinv[(1, k, ier)] = rinv[(1, k, ier)] / a;
                rinv[(2, k, ier)] = rinv[(2, k, ier)] / a;
                rinv[(1, 2, ier)] = rinv[(1, 2, ier)] / sqrt(a);
                rinv[(2, 2, ier)] = rinv[(2, 2, ier)] / sqrt(a);
            }
        }
    } else if rmat[(2, k21, ier)] != zero {
        rinv[(1, k11, ier)] = 0.0;
        rinv[(1, k21, ier)] = 0.0;
        rinv[(1, k22, ier)] = 0.0;
        let a = rmat[(2, k11, ier)] * rmat[(2, k22, ier)] - sq(rmat[(2, k21, ier)]);
        rinv[(2, k11, ier)] = -rmat[(2, k22, ier)] / a;
        rinv[(2, k21, ier)] = rmat[(2, k21, ier)] / a;
        rinv[(2, k22, ier)] = -rmat[(2, k11, ier)] / a;
    } else {
        rinv[(1, k11, ier)] = 0.0;
        rinv[(1, k21, ier)] = 0.0;
        rinv[(1, k22, ier)] = 0.0;
        rinv[(2, k11, ier)] = -1.0 / rmat[(2, k11, ier)];
        rinv[(2, k21, ier)] = 0.0;
        rinv[(2, k22, ier)] = -1.0 / rmat[(2, k22, ier)];
    }
}

fn threech(rmat: &mut FA, rinv: &mut FA, ier: i32) {
    let zero = 0.0f64;
    let (k11, k21, k31, k22, k32, k33) = (1, 2, 4, 3, 5, 6);
    let (mut a1, mut a2, mut a3) = (0.0, 0.0, 0.0);
    // scale3 receives rmat(2,*): the first nerm slice, whatever ier is
    scale3(&mut a1, &mut a2, &mut a3, &mut rmat.v[..]);
    let mut izz = 0;
    for i in 1..=6 {
        if rmat[(1, i, ier)] != zero {
            izz = 1;
        }
    }
    if rmat[(2, 2, ier)] != zero || rmat[(2, 4, ier)] != zero || rmat[(2, 5, ier)] != zero {
        izz = 1;
    }
    if izz == 0 {
        for i in 1..=6 {
            rinv[(1, i, ier)] = 0.0;
            rinv[(2, i, ier)] = 0.0;
        }
        rinv[(2, 1, ier)] = -1.0 / rmat[(2, 1, ier)];
        rinv[(2, 3, ier)] = -1.0 / rmat[(2, 3, ier)];
        rinv[(2, 6, ier)] = -1.0 / rmat[(2, 6, ier)];
        return;
    }
    let r = |i: i32, j: i32| rmat[(i, j, ier)];
    let fc1r = r(1, k22) * r(1, k33) - r(2, k22) * r(2, k33) - sq(r(1, k32)) + sq(r(2, k32));
    let fc1i = r(2, k22) * r(1, k33) + r(1, k22) * r(2, k33) - 2.0 * r(1, k32) * r(2, k32);
    let fc2r = r(1, k21) * r(1, k33) - r(2, k21) * r(2, k33) - r(1, k31) * r(1, k32) + r(2, k31) * r(2, k32);
    let fc2i = r(2, k21) * r(1, k33) + r(1, k21) * r(2, k33) - r(1, k31) * r(2, k32) - r(2, k31) * r(1, k32);
    let fc3r = r(1, k21) * r(1, k32) - r(2, k21) * r(2, k32) - r(1, k31) * r(1, k22) + r(2, k31) * r(2, k22);
    let fc3i = r(2, k21) * r(1, k32) + r(1, k21) * r(2, k32) - r(2, k31) * r(1, k22) - r(1, k31) * r(2, k22);
    let aar = r(1, k11) * fc1r - r(1, k21) * fc2r + r(1, k31) * fc3r - (r(2, k11) * fc1i - r(2, k21) * fc2i + r(2, k31) * fc3i);
    let aai = r(1, k11) * fc1i - r(1, k21) * fc2i + r(1, k31) * fc3i + (r(2, k11) * fc1r - r(2, k21) * fc2r + r(2, k31) * fc3r);
    let (dcr, dci);
    if aar + aai != aai {
        if aar + aai != aar {
            let aa = sq(aar) + sq(aai);
            dcr = aar / aa;
            dci = -aai / aa;
        } else {
            dcr = 1.0 / aar;
            dci = -(aai / aar) / aar;
        }
    } else {
        dcr = (aar / aai) / aai;
        dci = -1.0 / aai;
    }
    let v11r = fc1r * dcr - fc1i * dci;
    let v11i = fc1r * dci + fc1i * dcr;
    let v21r = -fc2r * dcr + fc2i * dci;
    let v21i = -fc2r * dci - fc2i * dcr;
    let t22a = r(1, k11) * r(1, k33) - sq(r(1, k31)) - r(2, k11) * r(2, k33) + sq(r(2, k31));
    let t22b = r(2, k11) * r(1, k33) - 2.0 * r(1, k31) * r(2, k31) + r(1, k11) * r(2, k33);
    let v22r = t22a * dcr - t22b * dci;
    let v22i = t22a * dci + t22b * dcr;
    let v31r = fc3r * dcr - fc3i * dci;
    let v31i = fc3r * dci + fc3i * dcr;
    let t32a = r(1, k11) * r(1, k32) - r(1, k21) * r(1, k31) - r(2, k11) * r(2, k32) + r(2, k21) * r(2, k31);
    let t32b = r(1, k11) * r(2, k32) - r(1, k21) * r(2, k31) + r(2, k11) * r(1, k32) - r(2, k21) * r(1, k31);
    let v32r = -t32a * dcr + t32b * dci;
    let v32i = -t32a * dci - t32b * dcr;
    let t33a = r(1, k11) * r(1, k22) - sq(r(1, k21)) - r(2, k11) * r(2, k22) + sq(r(2, k21));
    let t33b = r(2, k11) * r(1, k22) - 2.0 * r(1, k21) * r(2, k21) + r(1, k11) * r(2, k22);
    let v33r = t33a * dcr - t33b * dci;
    let v33i = t33a * dci + t33b * dcr;
    rinv[(1, k11, ier)] = v11r;
    rinv[(2, k11, ier)] = v11i;
    rinv[(1, k21, ier)] = v21r;
    rinv[(2, k21, ier)] = v21i;
    rinv[(1, k22, ier)] = v22r;
    rinv[(2, k22, ier)] = v22i;
    rinv[(1, k31, ier)] = v31r;
    rinv[(2, k31, ier)] = v31i;
    rinv[(1, k32, ier)] = v32r;
    rinv[(2, k32, ier)] = v32i;
    rinv[(1, k33, ier)] = v33r;
    rinv[(2, k33, ier)] = v33i;
    // unscale3 receives rmat(2,*) and rinv(2,*): the first slices
    let (rm, ri) = (&mut rmat.v[..], &mut rinv.v[..]);
    unscale3(a1, a2, a3, rm, ri);
}

/// rmat(2,*) flat: element (i,k) is r[(i-1)+2*(k-1)].
#[inline]
fn p2(i: i32, k: i32) -> usize {
    ((i - 1) + 2 * (k - 1)) as usize
}

fn scale3(a1: &mut f64, a2: &mut f64, a3: &mut f64, r: &mut [f64]) {
    let aa = 1.0e10f64;
    *a1 = 0.0;
    *a2 = 0.0;
    *a3 = 0.0;
    if r[p2(1, 1)].abs() >= aa || r[p2(2, 1)].abs() >= aa {
        let mut bb = r[p2(1, 1)].abs();
        let cc = r[p2(2, 1)].abs();
        if cc > bb {
            bb = cc;
        }
        *a1 = sqrt(bb);
        r[p2(1, 1)] /= bb;
        r[p2(2, 1)] /= bb;
        r[p2(1, 2)] /= *a1;
        r[p2(2, 2)] /= *a1;
        r[p2(1, 4)] /= *a1;
        r[p2(2, 4)] /= *a1;
    }
    if r[p2(1, 3)].abs() >= aa || r[p2(2, 3)].abs() >= aa {
        let mut bb = r[p2(1, 3)].abs();
        let cc = r[p2(2, 3)].abs();
        if cc > bb {
            bb = cc;
        }
        *a2 = sqrt(bb);
        r[p2(1, 2)] /= *a2;
        r[p2(2, 2)] /= *a2;
        r[p2(1, 3)] /= bb;
        r[p2(2, 3)] /= bb;
        r[p2(1, 5)] /= *a2;
        r[p2(2, 5)] /= *a2;
    }
    if r[p2(1, 6)].abs() >= aa || r[p2(2, 6)].abs() >= aa {
        let mut bb = r[p2(1, 6)].abs();
        let cc = r[p2(2, 6)].abs();
        if cc > bb {
            bb = cc;
        }
        *a3 = sqrt(bb);
        r[p2(1, 4)] /= *a3;
        r[p2(2, 4)] /= *a3;
        r[p2(1, 5)] /= *a3;
        r[p2(2, 5)] /= *a3;
        r[p2(1, 6)] /= bb;
        r[p2(2, 6)] /= bb;
    }
}

fn unscale3(a1: f64, a2: f64, a3: f64, r: &mut [f64], ri: &mut [f64]) {
    let zero = 0.0f64;
    if a1 > zero {
        let bb = sq(a1);
        r[p2(1, 1)] *= bb;
        r[p2(2, 1)] *= bb;
        r[p2(1, 2)] *= a1;
        r[p2(2, 2)] *= a1;
        r[p2(1, 4)] *= a1;
        r[p2(2, 4)] *= a1;
        ri[p2(1, 1)] /= bb;
        ri[p2(2, 1)] /= bb;
        ri[p2(1, 2)] /= a1;
        ri[p2(2, 2)] /= a1;
        ri[p2(1, 4)] /= a1;
        ri[p2(2, 4)] /= a1;
    }
    if a2 > zero {
        let bb = sq(a2);
        r[p2(1, 2)] *= a2;
        r[p2(2, 2)] *= a2;
        r[p2(1, 3)] *= bb;
        r[p2(2, 3)] *= bb;
        r[p2(1, 5)] *= a2;
        r[p2(2, 5)] *= a2;
        ri[p2(1, 2)] /= a2;
        ri[p2(2, 2)] /= a2;
        ri[p2(1, 3)] /= bb;
        ri[p2(2, 3)] /= bb;
        ri[p2(1, 5)] /= a2;
        ri[p2(2, 5)] /= a2;
    }
    if a3 > zero {
        let bb = sq(a3);
        r[p2(1, 4)] *= a3;
        r[p2(2, 4)] *= a3;
        r[p2(1, 5)] *= a3;
        r[p2(2, 5)] *= a3;
        r[p2(1, 6)] *= bb;
        r[p2(2, 6)] *= bb;
        ri[p2(1, 4)] /= a3;
        ri[p2(2, 4)] /= a3;
        ri[p2(1, 5)] /= a3;
        ri[p2(2, 5)] /= a3;
        ri[p2(1, 6)] /= bb;
        ri[p2(2, 6)] /= bb;
    }
}

/// yfour receives rmat(2,*) and rinv(2,*) (first slices).
fn yfour(c: &mut Ctx, rmat: &mut FA, rinv: &mut FA, nchan: i32) {
    let mut iii = vec![0i32; nchan as usize + 1];
    let mut info = 0;
    xspfa(&mut rmat.v[..], nchan, &mut iii, &mut info);
    if info != 0 {
        let s = format!(" Problem in xspfa with info={}", ffmt::i(info, 5));
        c.io.se(&s);
        c.io.so(&s);
    }
    let mut kj = 0;
    for k in 1..=nchan {
        for _j in 1..=k {
            kj += 1;
        }
    }
    let mut dummy = vec![0.0f64; 2 * kj as usize + 2 * nchan as usize + SLACK];
    kj = 0;
    for k in 1..=nchan {
        for j in 1..=nchan {
            dummy[p2(1, j)] = 0.0;
            dummy[p2(2, j)] = 0.0;
        }
        dummy[p2(1, k)] = 1.0;
        xspsl(&rmat.v[..], nchan, &iii, &mut dummy);
        for j in 1..=k {
            kj += 1;
            rinv.v[p2(1, kj)] = dummy[p2(1, j)];
            rinv.v[p2(2, kj)] = dummy[p2(2, j)];
        }
    }
}

// LINPACK-style complex symmetric packed factorisation. `ap` is ap(2,*)
// flat; a sub-array argument `ap(1,m)` is the slice starting at p2(1,m).

fn ixamax(n: i32, sx: &[f64], incx: i32) -> i32 {
    if n < 1 {
        return 0;
    }
    if n == 1 {
        return 1;
    }
    let mut ix_max = 1;
    if incx != 1 {
        let mut ix = 1;
        let mut smax = sq(sx[p2(1, 1)]) + sq(sx[p2(2, 1)]);
        ix += incx;
        for i in 2..=n {
            let aa = sq(sx[p2(1, ix)]) + sq(sx[p2(2, ix)]);
            if aa > smax {
                ix_max = i;
                smax = aa;
            }
            ix += incx;
        }
    } else {
        let mut smax = sq(sx[p2(1, 1)]) + sq(sx[p2(2, 1)]);
        for i in 2..=n {
            let aa = sq(sx[p2(1, i)]) + sq(sx[p2(2, i)]);
            if aa > smax {
                ix_max = i;
                smax = aa;
            }
        }
    }
    ix_max
}

/// y := y + (sa + i sai) x, with x at offset xo and y at offset yo of the
/// same buffer `a` (both are views into ap or b in the Fortran).
fn xaxpy2(n: i32, sa: f64, sai: f64, x: &[f64], xo: usize, y: &mut [f64], yo: usize) {
    let zero = 0.0f64;
    if n <= 0 {
        return;
    }
    if sa == zero && sai == zero {
        return;
    }
    let step = |i: i32, y: &mut [f64]| {
        let (x1, x2) = (x[xo + p2(1, i)], x[xo + p2(2, i)]);
        let aa = y[yo + p2(1, i)] + sa * x1 - sai * x2;
        y[yo + p2(2, i)] = y[yo + p2(2, i)] + sa * x2 + sai * x1;
        y[yo + p2(1, i)] = aa;
    };
    let m = n % 4;
    if m != 0 {
        for i in 1..=m {
            step(i, y);
        }
        if n < 4 {
            return;
        }
    }
    let mut i = m + 1;
    while i <= n {
        step(i, y);
        step(i + 1, y);
        step(i + 2, y);
        step(i + 3, y);
        i += 4;
    }
}

/// xaxpy where x and y are both views into the same array `a`.
fn xaxpy_same(n: i32, sa: f64, sai: f64, a: &mut [f64], xo: usize, yo: usize) {
    let zero = 0.0f64;
    if n <= 0 || (sa == zero && sai == zero) {
        return;
    }
    let step = |i: i32, a: &mut [f64]| {
        let (x1, x2) = (a[xo + p2(1, i)], a[xo + p2(2, i)]);
        let aa = a[yo + p2(1, i)] + sa * x1 - sai * x2;
        a[yo + p2(2, i)] = a[yo + p2(2, i)] + sa * x2 + sai * a[xo + p2(1, i)];
        a[yo + p2(1, i)] = aa;
    };
    let m = n % 4;
    if m != 0 {
        for i in 1..=m {
            step(i, a);
        }
        if n < 4 {
            return;
        }
    }
    let mut i = m + 1;
    while i <= n {
        step(i, a);
        step(i + 1, a);
        step(i + 2, a);
        step(i + 3, a);
        i += 4;
    }
}

fn xdot(xdoti: &mut f64, n: i32, sx: &[f64], xo: usize, sy: &[f64], yo: usize) -> f64 {
    let mut stemp = 0.0f64;
    let mut stempi = 0.0f64;
    *xdoti = 0.0;
    if n <= 0 {
        return 0.0;
    }
    let x1 = |i: i32| sx[xo + p2(1, i)];
    let x2 = |i: i32| sx[xo + p2(2, i)];
    let y1 = |i: i32| sy[yo + p2(1, i)];
    let y2 = |i: i32| sy[yo + p2(2, i)];
    let m = n % 5;
    if m != 0 {
        for i in 1..=m {
            stemp = stemp + x1(i) * y1(i) - x2(i) * y2(i);
            stempi = stempi + x2(i) * y1(i) + x1(i) * y2(i);
        }
        if n < 5 {
            *xdoti = stempi;
            return stemp;
        }
    }
    let mut i = m + 1;
    while i <= n {
        stemp = stemp + x1(i) * y1(i) - x2(i) * y2(i) + x1(i + 1) * y1(i + 1) - x2(i + 1) * y2(i + 1)
            + x1(i + 2) * y1(i + 2)
            - x2(i + 2) * y2(i + 2)
            + x1(i + 3) * y1(i + 3)
            - x2(i + 3) * y2(i + 3)
            + x1(i + 4) * y1(i + 4)
            - x2(i + 4) * y2(i + 4);
        stempi = stempi + x2(i) * y1(i) + x1(i) * y2(i) + x2(i + 1) * y1(i + 1) + x1(i + 1) * y2(i + 1)
            + x2(i + 2) * y1(i + 2)
            + x1(i + 2) * y2(i + 2)
            + x2(i + 3) * y1(i + 3)
            + x1(i + 3) * y2(i + 3)
            + x2(i + 4) * y1(i + 4)
            + x1(i + 4) * y2(i + 4);
        i += 5;
    }
    *xdoti = stempi;
    stemp
}

fn xswap_same(n: i32, a: &mut [f64], xo: usize, yo: usize) {
    if n <= 0 {
        return;
    }
    // the unrolled loops in the Fortran do the same element swaps in order
    for i in 1..=n {
        for k in 1..=2 {
            let (u, v) = (xo + p2(k, i), yo + p2(k, i));
            a.swap(u, v);
        }
    }
}

fn xspfa(ap: &mut [f64], n: i32, kpvt: &mut [i32], info: &mut i32) {
    let zero = 0.0f64;
    let alpha = (1.0 + sqrt(17.0e0)) / 8.0;
    let a = |ap: &[f64], i: i32, m: i32| ap[p2(i, m)];
    *info = 0;
    let mut im = 0;
    let mut k = n;
    let mut ik = (n * (n - 1)) / 2;
    loop {
        if k == 0 {
            break;
        }
        if k <= 1 {
            kpvt[1] = 1;
            if a(ap, 1, 1) == zero && a(ap, 2, 1) == zero {
                *info = 1;
            }
            break;
        }
        let km1 = k - 1;
        let kk = ik + k;
        let absakk = sq(a(ap, 1, kk)) + sq(a(ap, 2, kk));
        let imax = ixamax(k - 1, &ap[p2(1, ik + 1)..], 1);
        let imk = ik + imax;
        let colmax = sq(a(ap, 1, imk)) + sq(a(ap, 2, imk));
        let kstep;
        let swap;
        if absakk >= alpha * colmax {
            kstep = 1;
            swap = false;
        } else {
            let mut rowmax = 0.0f64;
            let imaxp1 = imax + 1;
            im = (imax * (imax - 1)) / 2;
            let mut imj = im + 2 * imax;
            for j in imaxp1..=k {
                let aa = sq(a(ap, 1, imj)) + sq(a(ap, 2, imj));
                rowmax = rowmax.max(aa);
                imj += j;
            }
            if imax != 1 {
                let jmax = ixamax(imax - 1, &ap[p2(1, im + 1)..], 1);
                let jmim = jmax + im;
                let aa = sq(a(ap, 1, jmim)) + sq(a(ap, 2, jmim));
                rowmax = rowmax.max(aa);
            }
            let imim = imax + im;
            let aa = sq(a(ap, 1, imim)) + sq(a(ap, 2, imim));
            if aa >= alpha * rowmax {
                kstep = 1;
                swap = true;
            } else if absakk >= alpha * colmax * (colmax / rowmax) {
                kstep = 1;
                swap = false;
            } else {
                kstep = 2;
                swap = imax != km1;
            }
        }
        if absakk.max(colmax) == zero {
            kpvt[k as usize] = k;
            *info = k;
        } else if kstep != 2 {
            if swap {
                xswap_same(imax, ap, p2(1, im + 1), p2(1, ik + 1));
                let mut imj = ik + imax;
                for jj in imax..=k {
                    let j = k + imax - jj;
                    let jk = ik + j;
                    ap.swap(p2(1, jk), p2(1, imj));
                    ap.swap(p2(2, jk), p2(2, imj));
                    imj -= j - 1;
                }
            }
            let mut ij = ik - (k - 1);
            for jj in 1..=km1 {
                let j = k - jj;
                let jk = ik + j;
                let aa = sq(a(ap, 1, kk)) + sq(a(ap, 2, kk));
                let dmulk = -(a(ap, 1, jk) * a(ap, 1, kk) + a(ap, 2, jk) * a(ap, 2, kk)) / aa;
                let dmulki = (a(ap, 1, jk) * a(ap, 2, kk) - a(ap, 2, jk) * a(ap, 1, kk)) / aa;
                xaxpy_same(j, dmulk, dmulki, ap, p2(1, ik + 1), p2(1, ij + 1));
                ap[p2(1, jk)] = dmulk;
                ap[p2(2, jk)] = dmulki;
                ij -= j - 1;
            }
            kpvt[k as usize] = k;
            if swap {
                kpvt[k as usize] = imax;
            }
        } else {
            let km1k = ik + k - 1;
            let ikm1 = ik - (k - 1);
            if swap {
                xswap_same(imax, ap, p2(1, im + 1), p2(1, ikm1 + 1));
                let mut imj = ikm1 + imax;
                for jj in imax..=km1 {
                    let j = km1 + imax - jj;
                    let jkm1 = ikm1 + j;
                    ap.swap(p2(1, jkm1), p2(1, imj));
                    ap.swap(p2(2, jkm1), p2(2, imj));
                    imj -= j - 1;
                }
                ap.swap(p2(1, km1k), p2(1, imk));
                ap.swap(p2(2, km1k), p2(2, imk));
            }
            let km2 = k - 2;
            if km2 != 0 {
                let aa = sq(a(ap, 1, km1k)) + sq(a(ap, 2, km1k));
                let ak = (a(ap, 1, kk) * a(ap, 1, km1k) + a(ap, 2, kk) * a(ap, 2, km1k)) / aa;
                let aki = (a(ap, 2, kk) * a(ap, 1, km1k) - a(ap, 1, kk) * a(ap, 2, km1k)) / aa;
                let km1km1 = ikm1 + k - 1;
                let akm1 = (a(ap, 1, km1km1) * a(ap, 1, km1k) + a(ap, 2, km1km1) * a(ap, 2, km1k)) / aa;
                let akm1i = (a(ap, 2, km1km1) * a(ap, 1, km1k) - a(ap, 1, km1km1) * a(ap, 2, km1k)) / aa;
                let denom = 1.0 - (ak * akm1 - aki * akm1i);
                let denomi = -(ak * akm1i + aki * akm1);
                let dd = sq(denom) + sq(denomi);
                let mut ij = ik - (k - 1) - (k - 2);
                for jj in 1..=km2 {
                    let j = km1 - jj;
                    let jk = ik + j;
                    let bk = (a(ap, 1, jk) * a(ap, 1, km1k) + a(ap, 2, jk) * a(ap, 2, km1k)) / aa;
                    let bki = (a(ap, 2, jk) * a(ap, 1, km1k) - a(ap, 1, jk) * a(ap, 2, km1k)) / aa;
                    let jkm1 = ikm1 + j;
                    let bkm1 = (a(ap, 1, jkm1) * a(ap, 1, km1k) + a(ap, 2, jkm1) * a(ap, 2, km1k)) / aa;
                    let bkm1i = (a(ap, 2, jkm1) * a(ap, 1, km1k) - a(ap, 1, jkm1) * a(ap, 2, km1k)) / aa;
                    let mut xx = akm1 * bk - akm1i * bki - bkm1;
                    let mut xxi = akm1 * bki + akm1i * bk - bkm1i;
                    let dmulk = (xx * denom + xxi * denomi) / dd;
                    let dmulki = (xxi * denom - xx * denomi) / dd;
                    xx = ak * bkm1 - aki * bkm1i - bk;
                    xxi = ak * bkm1i + aki * bkm1 - bki;
                    let dmlkm1 = (xx * denom + xxi * denomi) / dd;
                    let dmlkmi = (xxi * denom - xx * denomi) / dd;
                    xaxpy_same(j, dmulk, dmulki, ap, p2(1, ik + 1), p2(1, ij + 1));
                    xaxpy_same(j, dmlkm1, dmlkmi, ap, p2(1, ikm1 + 1), p2(1, ij + 1));
                    ap[p2(1, jk)] = dmulk;
                    ap[p2(2, jk)] = dmulki;
                    ap[p2(1, jkm1)] = dmlkm1;
                    ap[p2(2, jkm1)] = dmlkmi;
                    ij -= j - 1;
                }
            }
            kpvt[k as usize] = 1 - k;
            if swap {
                kpvt[k as usize] = -imax;
            }
            kpvt[(k - 1) as usize] = kpvt[k as usize];
        }
        ik -= k - 1;
        if kstep == 2 {
            ik -= k - 2;
        }
        k -= kstep;
    }
}

fn xspsl(ap: &[f64], n: i32, kpvt: &[i32], b: &mut [f64]) {
    let a = |i: i32, m: i32| ap[p2(i, m)];
    let mut k = n;
    let mut ik = (n * (n - 1)) / 2;
    while k != 0 {
        let kk = ik + k;
        if kpvt[k as usize] >= 0 {
            if k != 1 {
                let kp = kpvt[k as usize];
                if kp != k {
                    b.swap(p2(1, k), p2(1, kp));
                    b.swap(p2(2, k), p2(2, kp));
                }
                let (s1, s2) = (b[p2(1, k)], b[p2(2, k)]);
                xaxpy2(k - 1, s1, s2, ap, p2(1, ik + 1), b, 0);
            }
            let aa = sq(a(1, kk)) + sq(a(2, kk));
            let ab = (b[p2(1, k)] * a(1, kk) + b[p2(2, k)] * a(2, kk)) / aa;
            b[p2(2, k)] = (b[p2(2, k)] * a(1, kk) - b[p2(1, k)] * a(2, kk)) / aa;
            b[p2(1, k)] = ab;
            k -= 1;
            ik -= k;
        } else {
            let ikm1 = ik - (k - 1);
            if k != 2 {
                let kp = kpvt[k as usize].abs();
                if kp != k - 1 {
                    b.swap(p2(1, k - 1), p2(1, kp));
                    b.swap(p2(2, k - 1), p2(2, kp));
                }
                let (s1, s2) = (b[p2(1, k)], b[p2(2, k)]);
                xaxpy2(k - 2, s1, s2, ap, p2(1, ik + 1), b, 0);
                let (s1, s2) = (b[p2(1, k - 1)], b[p2(2, k - 1)]);
                xaxpy2(k - 2, s1, s2, ap, p2(1, ikm1 + 1), b, 0);
            }
            let km1k = ik + k - 1;
            let kk = ik + k;
            let aa = sq(a(1, km1k)) + sq(a(2, km1k));
            let ak = (a(1, kk) * a(1, km1k) + a(2, kk) * a(2, km1k)) / aa;
            let aki = (a(2, kk) * a(1, km1k) - a(1, kk) * a(2, km1k)) / aa;
            let km1km1 = ikm1 + k - 1;
            let akm1 = (a(1, km1km1) * a(1, km1k) + a(2, km1km1) * a(2, km1k)) / aa;
            let akm1i = (a(2, km1km1) * a(1, km1k) - a(1, km1km1) * a(2, km1k)) / aa;
            let bk = (b[p2(1, k)] * a(1, km1k) + b[p2(2, k)] * a(2, km1k)) / aa;
            let bki = (b[p2(2, k)] * a(1, km1k) - b[p2(1, k)] * a(2, km1k)) / aa;
            let bkm1 = (b[p2(1, k - 1)] * a(1, km1k) + b[p2(2, k - 1)] * a(2, km1k)) / aa;
            let bkm1i = (b[p2(2, k - 1)] * a(1, km1k) - b[p2(1, k - 1)] * a(2, km1k)) / aa;
            let denom = ak * akm1 - aki * akm1i - 1.0;
            let denomi = ak * akm1i + aki * akm1;
            let dd = sq(denom) + sq(denomi);
            let mut ab = akm1 * bk - akm1i * bki - bkm1;
            let mut abi = akm1i * bk + akm1 * bki - bkm1i;
            b[p2(1, k)] = (ab * denom + abi * denomi) / dd;
            b[p2(2, k)] = (abi * denom - ab * denomi) / dd;
            ab = ak * bkm1 - aki * bkm1i - bk;
            abi = aki * bkm1 + ak * bkm1i - bki;
            b[p2(1, k - 1)] = (ab * denom + abi * denomi) / dd;
            b[p2(2, k - 1)] = (abi * denom - ab * denomi) / dd;
            k -= 2;
            ik = ik - (k + 1) - k;
        }
    }
    k = 1;
    ik = 0;
    while k <= n {
        let mut xdoti = 0.0;
        if kpvt[k as usize] >= 0 {
            if k != 1 {
                let d = xdot(&mut xdoti, k - 1, ap, p2(1, ik + 1), b, 0);
                b[p2(1, k)] += d;
                b[p2(2, k)] += xdoti;
                let kp = kpvt[k as usize];
                if kp != k {
                    b.swap(p2(1, k), p2(1, kp));
                    b.swap(p2(2, k), p2(2, kp));
                }
            }
            ik += k;
            k += 1;
        } else {
            if k != 1 {
                let d = xdot(&mut xdoti, k - 1, ap, p2(1, ik + 1), b, 0);
                b[p2(1, k)] += d;
                b[p2(2, k)] += xdoti;
                let ikp1 = ik + k;
                let d = xdot(&mut xdoti, k - 1, ap, p2(1, ikp1 + 1), b, 0);
                b[p2(1, k + 1)] += d;
                b[p2(2, k + 1)] += xdoti;
                let kp = kpvt[k as usize].abs();
                if kp != k {
                    b.swap(p2(1, k), p2(1, kp));
                    b.swap(p2(2, k), p2(2, kp));
                }
            }
            ik = ik + k + k + 1;
            k += 2;
        }
    }
}
