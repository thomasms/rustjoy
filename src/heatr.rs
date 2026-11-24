//! HEATR: prompt kerma (heating) and radiation damage energy production,
//! ported from NJOY2016 `heatr.f90` (module heatm).
//!
//! The routines keep the Fortran structure: module variables are fields of
//! `Heatr`, `save`d locals of each routine are kept in a small struct of
//! their own, and labels become `match` states where the control flow
//! needs them. Arrays keep the Fortran (1-based) indices.

#![allow(clippy::too_many_arguments)]
#![allow(clippy::collapsible_if)]
#![allow(clippy::collapsible_else_if)]

use crate::endf::Ctx;
use crate::farr::{F, RV};
use crate::ffmt;
use crate::fio::{Kind, Val};
use crate::fmath::{cosh, exp, int, log, log10, nint, pow, powi, sinh, sqrt};
use crate::mathm::legndr;
use crate::util::physics::{AMASSN, AMU, ANRATIO, CLIGHT, DNRATIO, EV, HNRATIO, PI, PNRATIO, TNRATIO};
use crate::util::sigfig;

const NBUF: i32 = 1000;
const NQAMAX: usize = 30;
const ILMAX: i32 = 200;
const MAXMF6: usize = 320;
/// Extra words past the Fortran array sizes, so that writes the Fortran
/// does past an array end (without checking) do not stop the port first.
const SLACK: usize = 4096;

/// 64-point Gauss-Legendre abscissae and weights (disbar, getsix, hgetco).
const QP64: [f64; 64] = [
    -9.99305042E-01, -9.96340117E-01, -9.91013371E-01, -9.83336254E-01, -9.73326828E-01, -9.61008800E-01,
    -9.46411375E-01, -9.29569172E-01, -9.10522137E-01, -8.89315446E-01, -8.65999398E-01, -8.40629296E-01,
    -8.13265315E-01, -7.83972359E-01, -7.52819907E-01, -7.19881850E-01, -6.85236313E-01, -6.48965471E-01,
    -6.11155355E-01, -5.71895646E-01, -5.31279464E-01, -4.89403146E-01, -4.46366017E-01, -4.02270158E-01,
    -3.57220158E-01, -3.11322872E-01, -2.64687162E-01, -2.17423644E-01, -1.69644420E-01, -1.21462819E-01,
    -7.29931218E-02, -2.43502927E-02, 2.43502927E-02, 7.29931218E-02, 1.21462819E-01, 1.69644420E-01,
    2.17423644E-01, 2.64687162E-01, 3.11322872E-01, 3.57220158E-01, 4.02270158E-01, 4.46366017E-01,
    4.89403146E-01, 5.31279464E-01, 5.71895646E-01, 6.11155355E-01, 6.48965471E-01, 6.85236313E-01,
    7.19881850E-01, 7.52819907E-01, 7.83972359E-01, 8.13265315E-01, 8.40629296E-01, 8.65999398E-01,
    8.89315446E-01, 9.10522137E-01, 9.29569172E-01, 9.46411375E-01, 9.61008800E-01, 9.73326828E-01,
    9.83336254E-01, 9.91013371E-01, 9.96340117E-01, 9.99305042E-01,
];
const QW64: [f64; 64] = [
    1.78328072E-03, 4.14703326E-03, 6.50445797E-03, 8.84675983E-03, 1.11681395E-02, 1.34630479E-02,
    1.57260305E-02, 1.79517158E-02, 2.01348232E-02, 2.22701738E-02, 2.43527026E-02, 2.63774697E-02,
    2.83396726E-02, 3.02346571E-02, 3.20579284E-02, 3.38051618E-02, 3.54722133E-02, 3.70551285E-02,
    3.85501532E-02, 3.99537411E-02, 4.12625632E-02, 4.24735151E-02, 4.35837245E-02, 4.45905582E-02,
    4.54916279E-02, 4.62847966E-02, 4.69681828E-02, 4.75401657E-02, 4.79993886E-02, 4.83447622E-02,
    4.85754674E-02, 4.86909570E-02, 4.86909570E-02, 4.85754674E-02, 4.83447622E-02, 4.79993886E-02,
    4.75401657E-02, 4.69681828E-02, 4.62847966E-02, 4.54916279E-02, 4.45905582E-02, 4.35837245E-02,
    4.24735151E-02, 4.12625632E-02, 3.99537411E-02, 3.85501532E-02, 3.70551285E-02, 3.54722133E-02,
    3.38051618E-02, 3.20579284E-02, 3.02346571E-02, 2.83396726E-02, 2.63774697E-02, 2.43527026E-02,
    2.22701738E-02, 2.01348232E-02, 1.79517158E-02, 1.57260305E-02, 1.34630479E-02, 1.11681395E-02,
    8.84675983E-03, 6.50445797E-03, 4.14703326E-03, 1.78328072E-03,
];
/// 4-point Gauss-Legendre (capdam, anadam, tabdam, tabsq6, tabsqr).
const QP4: [f64; 4] = [-0.86114e0, -0.33998e0, 0.33998e0, 0.86114e0];
const QW4: [f64; 4] = [0.34785e0, 0.65215e0, 0.65215e0, 0.34785e0];

// ---------------------------------------------------------------------
// saved locals of the routines
// ---------------------------------------------------------------------

#[derive(Default, Clone)]
struct CapdamS {
    en: f64,
    damn: f64,
    el: f64,
    daml: f64,
    zx: f64,
    ax: f64,
    denom: f64,
    z: f64,
    aw1fac: f64,
}

#[derive(Default, Clone)]
struct DisbarS {
    imiss: i32,
    thresh: f64,
    awp: f64,
    afact: f64,
    arat: f64,
    en: f64,
    cn: f64,
    el: f64,
    cl: f64,
    z: f64,
    daml: f64,
    damn: f64,
    enx: f64,
    enext: f64,
}

#[derive(Default, Clone)]
struct ConbarS {
    nktot: i32,
    awr: f64,
    awfac: f64,
    aw1fac: f64,
    z: f64,
    za: f64,
    mtt: i32,
    loc: [i32; 13],
    d1: [f64; 13],
    d2: [f64; 13],
    e1: [f64; 13],
    e2: [f64; 13],
    mf1: i32,
    mt1: i32,
    // t1, t2 are not saved in the Fortran, but they are set only when the
    // routine is initialized (e=0) and read on later calls: the values of
    // the initializing call are kept here.
    t1: f64,
    t2: f64,
}

#[derive(Default, Clone)]
struct HgtyldS {
    lnu: i32,
    ip: i32,
    ir: i32,
}

#[derive(Default, Clone)]
struct SixbarS {
    iraw: i32,
    nne: i32,
    ne: i32,
    law: i32,
    lang: i32,
    lep: i32,
    intl: i32,
    elo: f64,
    ehi: f64,
    flo: f64,
    fhi: f64,
    dlo: f64,
    dhi: f64,
    disc102: f64,
    zp: f64,
    ap: f64,
    zt: f64,
    at: f64,
}

#[derive(Default, Clone)]
struct H6cmS {
    na: i32,
    eps: f64,
    xc: f64,
    ndnow: i32,
    npnow: i32,
    ncnow: i32,
    elmax: f64,
    e: f64,
    epmax: f64,
}

#[derive(Default, Clone)]
struct H6ddxS {
    efirst: f64,
    enow: f64,
    nl: i32,
    inow: i32,
    mnow: i32,
    lnow: i32,
    ncnow: i32,
    na: i32,
    illdef: i32,
}

#[derive(Default, Clone)]
struct H6disS {
    nl: i32,
    inow: i32,
    lnow: i32,
    mnow: i32,
    ncnow: i32,
    na: i32,
    enow: f64,
}

#[derive(Default, Clone)]
struct H6pspS {
    cn: f64,
    ex: f64,
    eimax: f64,
}

#[derive(Clone)]
struct HgtfleS {
    iso: i32,
    iraw: i32,
    ir: i32,
    nne: i32,
    ne: i32,
    inn: i32,
    elo: f64,
    ehi: f64,
    nlo: i32,
    nhi: i32,
    flo: [f64; 66],
    fhi: [f64; 66],
    ltt: i32,
    ltt3: i32,
    lttn: i32,
    lct: i32,
    awr: f64,
}

impl Default for HgtfleS {
    fn default() -> Self {
        HgtfleS {
            iso: 0,
            iraw: 0,
            ir: 0,
            nne: 0,
            ne: 0,
            inn: 0,
            elo: 0.0,
            ehi: 0.0,
            nlo: 0,
            nhi: 0,
            flo: [0.0; 66],
            fhi: [0.0; 66],
            ltt: 0,
            ltt3: 0,
            lttn: 0,
            lct: 0,
            awr: 0.0,
        }
    }
}

#[derive(Default, Clone)]
struct GambarS {
    lf: i32,
    z: f64,
    awr: f64,
    jstart: i32,
    jend: i32,
    istart: i32,
    ilo: i32,
    ihi: i32,
    nnt: i32,
    ne: i32,
    nne: i32,
    inn: i32,
    nbt: i32,
    elo: f64,
    ehi: f64,
    flo: f64,
    fhi: f64,
    glo: f64,
    ghi: f64,
    hlo: f64,
    hhi: f64,
}

/// HEATR module variables.
pub struct Heatr {
    nendf: i32,
    nin: i32,
    nout: i32,
    nplot: i32,
    matd: i32,
    mgam: i32,
    npk: i32,
    mtp: [i32; 29],
    nqa: i32,
    mta: [i32; NQAMAX + 1],
    mt303: i32,
    mt19: i32,
    ne: i32,
    kchk: i32,
    iprint: i32,
    lqs: [i32; NQAMAX + 1],
    qbar: RV,
    qa: [f64; NQAMAX + 1],
    efirst: f64,
    elast: f64,
    za: f64,
    awr: f64,
    elist: [f64; ILMAX as usize + 2],
    brk: f64,
    qdel: f64,
    etabmax: f64,
    mt103: i32,
    mt104: i32,
    mt105: i32,
    mt106: i32,
    mt107: i32,
    mt16: i32,
    miss4: [i32; 251],
    nmiss4: i32,
    i6: i32,
    mt6: [i32; MAXMF6 + 1],
    i6g: i32,
    mt6no: [i32; MAXMF6 + 1],
    i6p: i32,
    mt6yp: [i32; MAXMF6 + 1],
    jp: i32,
    jpn: i32,
    jpp: i32,
    q: f64,
    zat: f64,
    awrt: f64,
    zap: f64,
    awp: f64,
    lct: i32,
    idame: i32,
    ebot: f64,
    etop: f64,
    mt458: i32,
    nply: i32,
    lfc: i32,
    ifc1: i32,
    ifc2: i32,
    ifc3: i32,
    ifc4: i32,
    ifc5: i32,
    ifc6: i32,
    c458: RV,
    /// cpoly(0:nply), hpoly(0:nply): index i is element i
    cpoly: Vec<f64>,
    hpoly: Vec<f64>,
    afr: RV,
    anp: RV,
    agp: RV,
    emc2: f64,
    tm: f64,
    rtm: f64,
    izat: i32,
    izap: i32,
    // saved locals
    cd: CapdamS,
    db: DisbarS,
    df_rel: f64,
    df_fl: f64,
    cb: ConbarS,
    hy: HgtyldS,
    sb: SixbarS,
    h6: H6cmS,
    hd: H6ddxS,
    hs: H6disS,
    hp: H6pspS,
    hf: HgtfleS,
    gb: GambarS,
}

impl Default for Heatr {
    fn default() -> Self {
        Self::new()
    }
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

/// Real array with the Fortran length n plus slack.
fn rv(n: i32) -> RV {
    RV::new(n.max(0) as usize + SLACK)
}

/// `write(nsyso,'(1x,1p,ke14.4)') vals` for up to k values per record
/// (format reversion starts a new record with `1x`).
fn e14_lines(vals: &[f64], k: usize) -> String {
    let mut s = String::new();
    for (i, v) in vals.iter().enumerate() {
        if i > 0 && i % k == 0 {
            s.push('\n');
        }
        if i % k == 0 {
            s.push(' ');
        }
        s.push_str(&ffmt::e1p(*v, 14, 4));
    }
    s
}

impl Heatr {
    pub fn new() -> Self {
        Heatr {
            nendf: 0,
            nin: 0,
            nout: 0,
            nplot: 0,
            matd: 0,
            mgam: 0,
            npk: 0,
            mtp: [0; 29],
            nqa: 0,
            mta: [0; NQAMAX + 1],
            mt303: 0,
            mt19: 0,
            ne: 0,
            kchk: 0,
            iprint: 0,
            lqs: [0; NQAMAX + 1],
            qbar: RV::new(1),
            qa: [0.0; NQAMAX + 1],
            efirst: 0.0,
            elast: 0.0,
            za: 0.0,
            awr: 0.0,
            elist: [0.0; ILMAX as usize + 2],
            brk: 0.0,
            qdel: 0.0,
            etabmax: 0.0,
            mt103: 0,
            mt104: 0,
            mt105: 0,
            mt106: 0,
            mt107: 0,
            mt16: 0,
            miss4: [0; 251],
            nmiss4: 0,
            i6: 0,
            mt6: [0; MAXMF6 + 1],
            i6g: 0,
            mt6no: [0; MAXMF6 + 1],
            i6p: 0,
            mt6yp: [0; MAXMF6 + 1],
            jp: 0,
            jpn: 0,
            jpp: 0,
            q: 0.0,
            zat: 0.0,
            awrt: 0.0,
            zap: 0.0,
            awp: 0.0,
            lct: 0,
            idame: 0,
            ebot: 0.0,
            etop: 0.0,
            mt458: 0,
            nply: 0,
            lfc: 0,
            ifc1: 0,
            ifc2: 0,
            ifc3: 0,
            ifc4: 0,
            ifc5: 0,
            ifc6: 0,
            c458: RV::new(1),
            cpoly: Vec::new(),
            hpoly: Vec::new(),
            afr: RV::new(1),
            anp: RV::new(1),
            agp: RV::new(1),
            emc2: 0.0,
            tm: 0.0,
            rtm: 0.0,
            izat: 0,
            izap: 0,
            cd: Default::default(),
            db: Default::default(),
            df_rel: 0.0,
            df_fl: 0.0,
            cb: Default::default(),
            hy: Default::default(),
            sb: Default::default(),
            h6: Default::default(),
            hd: Default::default(),
            hs: Default::default(),
            hp: Default::default(),
            hf: Default::default(),
            gb: Default::default(),
        }
    }

    // ================================================================ heatr
    pub fn heatr(&mut self, c: &mut Ctx) {
        const NPKMAX: i32 = 28;
        let qtest = 99.0e6f64;
        let qflag = -1.0e-9f64;
        let zero = 0.0f64;
        let maxqbar: i32 = 10000;
        let (mut nb, mut nw) = (0i32, 0i32);
        let mut mtk = [0i32; 2 + NPKMAX as usize + 1];
        let mut z = RV::new(17 + SLACK);

        //--start
        let time = c.timer();
        c.io.so(&format!("\n heatr...prompt kerma{}{}s", " ".repeat(48), ffmt::f(time, 8, 1)));
        c.io.se(&format!("\n heatr...{}{}s", " ".repeat(60), ffmt::f(time, 8, 1)));

        //--read user input.
        self.nplot = 0;
        let v = c.io.input.read(&[Kind::Int, Kind::Int, Kind::Int, Kind::Int]);
        self.nendf = geti(&v[0], self.nendf);
        self.nin = geti(&v[1], self.nin);
        self.nout = geti(&v[2], self.nout);
        self.nplot = geti(&v[3], self.nplot);
        c.openz(self.nendf, 0);
        c.openz(self.nin, 0);
        c.openz(self.nout, 1);
        if self.nplot != 0 {
            c.openz(self.nplot, 1);
        }
        self.npk = 0;
        self.nqa = 0;
        let mut ntemp = 0;
        let mut local = 0;
        self.iprint = 0;
        self.brk = 0.0;
        let v = c.io.input.read(&[Kind::Int, Kind::Int, Kind::Int, Kind::Int, Kind::Int, Kind::Int, Kind::Real]);
        self.matd = geti(&v[0], self.matd);
        self.npk = geti(&v[1], self.npk);
        self.nqa = geti(&v[2], self.nqa);
        ntemp = geti(&v[3], ntemp);
        local = geti(&v[4], local);
        self.iprint = geti(&v[5], self.iprint);
        self.brk = getr(&v[6], self.brk);
        self.kchk = 0;
        if self.iprint == 2 {
            self.kchk = 1;
        }
        if self.iprint == 2 {
            self.iprint = 1;
        }
        let mut npkk = self.npk + 3;
        if self.kchk == 1 {
            npkk = 3 * self.npk + 7;
        }
        if npkk > NPKMAX {
            c.error("heatr", "requested too many kerma mt-s (6+mt301 allowed).", " ");
        }
        if self.nqa > NQAMAX as i32 {
            c.error("heatr", "requested too many q values.", " ");
        }
        if self.npk > 0 {
            let kinds: Vec<Kind> = (0..self.npk).map(|_| Kind::Int).collect();
            let v = c.io.input.read(&kinds);
            for i in 1..=self.npk as usize {
                mtk[i] = geti(&v[i - 1], mtk[i]);
            }
        }
        let mut tmp = RV::new(maxqbar as usize + 1000 + SLACK);
        let mut loc = 1i32;
        if self.nqa != 0 {
            let kinds: Vec<Kind> = (0..self.nqa).map(|_| Kind::Int).collect();
            let v = c.io.input.read(&kinds);
            for i in 1..=self.nqa as usize {
                self.mta[i] = geti(&v[i - 1], self.mta[i]);
            }
            let kinds: Vec<Kind> = (0..self.nqa).map(|_| Kind::Real).collect();
            let v = c.io.input.read(&kinds);
            for i in 1..=self.nqa as usize {
                self.qa[i] = getr(&v[i - 1], self.qa[i]);
            }
            for i in 1..=self.nqa as usize {
                self.lqs[i] = 0;
                if self.qa[i] >= qtest {
                    self.lqs[i] = loc;
                    let nz = 1000;
                    let flag = qflag;
                    for j in 1..=nz {
                        tmp[loc + j - 1] = flag;
                    }
                    let kinds: Vec<Kind> = (0..nz).map(|_| Kind::Real).collect();
                    let v = c.io.input.read(&kinds);
                    for j in 1..=nz {
                        tmp[loc + j - 1] = getr(&v[(j - 1) as usize], tmp[loc + j - 1]);
                    }
                    let mut nz0 = 0;
                    for j in 1..=nz {
                        if tmp[loc + j - 1] > flag {
                            nz0 = j;
                        }
                    }
                    loc += nz0;
                    if loc > maxqbar {
                        c.error("heatr", "too much energy-dependent q data", " ");
                    }
                }
            }
        }
        let nzq = loc - 1;
        if nzq > 0 {
            self.qbar = rv(nzq);
            for i in 1..=nzq {
                self.qbar[i] = tmp[i];
            }
        }
        drop(tmp);
        if (self.nin < 0 && self.nout > 0) || (self.nin > 0 && self.nout < 0) {
            c.error("heatr", "mode conversion not allowed between nin and nout.", " ");
        }

        //--check the partial kermas
        if self.npk > 0 {
            let mut isave = 0;
            for i in 1..=self.npk as usize {
                if mtk[i] == 443 && self.kchk == 0 {
                    self.kchk = 2;
                }
                if mtk[i] == 301 {
                    isave = i as i32;
                }
            }
            if isave > 0 {
                c.mess("heatr", "mt301 always calculated", "--you do not need to ask for it.");
                self.npk -= 1;
                if isave <= self.npk {
                    for i in isave..=self.npk {
                        mtk[i as usize] = mtk[(i + 1) as usize];
                    }
                }
            }
            for i in 1..=self.npk {
                self.mtp[(self.npk + 3 - i) as usize] = mtk[i as usize];
            }
        }
        self.mtp[1] = 0;
        self.mtp[2] = 301;
        self.npk += 2;
        if self.npk > 3 {
            let n = self.npk as usize;
            horder(&mut self.mtp[1..=n]);
        }

        //--list all input parameters
        c.io.so(&format!(
            "\n input endf unit ...................... {}\n input pendf unit ..................... {}\n output pendf unit .................... {}",
            ffmt::i(self.nendf, 10),
            ffmt::i(self.nin, 10),
            ffmt::i(self.nout, 10)
        ));
        c.io.so(&format!(" mat to be processed .................. {}", ffmt::i(self.matd, 10)));
        c.io.so(&format!(
            " no. temperatures (0=all) ............. {}\n gamma heat (0 nonlocal, 1 local) ..... {}\n print option (0 min, 1 more, 2 chk) .. {}",
            ffmt::i(ntemp, 10),
            ffmt::i(local, 10),
            ffmt::i(self.iprint, 10)
        ));
        if self.brk == zero {
            c.io.so(" damage displacement energy ...........    default");
        } else {
            c.io.so(&format!(" damage displacement energy ........... {} ev", ffmt::f(self.brk, 10, 1)));
        }
        if ntemp == 0 {
            ntemp = 100;
        }
        if self.npk >= 3 {
            c.io.so(&format!(" partial kerma mt-s desired ........... {}", ffmt::i(self.mtp[3], 10)));
        }
        if self.npk >= 4 {
            let mut s = String::new();
            for i in 4..=self.npk as usize {
                if i > 4 {
                    s.push('\n');
                }
                s.push_str(&" ".repeat(40));
                s.push_str(&ffmt::i(self.mtp[i], 10));
            }
            c.io.so(&s);
        }
        if self.nqa != 0 {
            c.io.so(&format!(" auxiliary reactions .................. {}", ffmt::i(self.mta[1], 10)));
            if self.nqa > 1 {
                let mut s = String::new();
                for i in 2..=self.nqa as usize {
                    if i > 2 {
                        s.push('\n');
                    }
                    s.push_str(&" ".repeat(40));
                    s.push_str(&ffmt::i(self.mta[i], 10));
                }
                c.io.so(&s);
            }
            c.io.so(&format!(" auxiliary q values ................... {}", ffmt::e1p(self.qa[1], 12, 4)));
            if self.nqa > 1 {
                let mut s = String::new();
                for i in 2..=self.nqa as usize {
                    if i > 2 {
                        s.push('\n');
                    }
                    s.push_str(&" ".repeat(40));
                    s.push_str(&ffmt::e1p(self.qa[i], 12, 4));
                }
                c.io.so(&s);
            }
            if nzq != 0 {
                for i in 1..=self.nqa as usize {
                    if self.qa[i] >= qtest {
                        c.nsh = 0;
                        c.math = 1;
                        c.mfh = 1;
                        c.mth = self.mta[i];
                        c.io.so(" input q ----");
                        // nsyso is unit 0 (ERROR_UNIT, connected to 'output'),
                        // and tab1io writes nothing for nout=0
                        c.tab1io(0, 0, 0, &mut self.qbar.all(), &mut nb, &mut nw);
                    }
                }
            }
        }

        //--assign scratch units
        let mut nscr = 10;
        if self.nendf < 0 {
            nscr = -nscr;
        }
        let mut nend4 = 13;
        if self.nendf < 0 {
            nend4 = -nend4;
        }
        let mut nend6 = 14;
        if self.nendf < 0 {
            nend6 = -nend6;
        }
        let mut iold = 11;
        let mut inew = 12;
        c.openz(-iold, 1);
        c.openz(-inew, 1);
        c.openz(nscr, 1);
        c.openz(nend4, 1);
        c.openz(nend6, 1);

        //--process input pendf tape
        c.repoz(self.nin);
        c.repoz(self.nout);
        c.nsh = 0;
        c.tpidio(self.nin, self.nout, 0, &mut z.all(), &mut nb, &mut nw);
        c.findf(self.matd, 1, 451, self.nin);
        let mut itemp = 0;

        //--loop over desired temperatures for this mat
        let mut idone = 0;
        while idone == 0 {
            c.contio(self.nin, 0, 0, &mut z.all(), &mut nb, &mut nw);
            if c.math != self.matd {
                idone = 1;
            } else {
                itemp += 1;
                if itemp > ntemp {
                    idone = 1;
                } else {
                    self.za = c.c1h;
                    self.awr = c.c2h;
                    self.emc2 = AMASSN * AMU * CLIGHT * CLIGHT / EV;
                    self.tm = self.emc2 * (self.awr + 1.0);
                    self.rtm = 1.0 / self.tm;

                    //--check on default damage displacement energy
                    if self.brk == zero {
                        let iz = nint(self.za / 1000.0);
                        self.brk = if iz == 4 {
                            31.0
                        } else if iz == 6 {
                            31.0
                        } else if iz == 12 {
                            25.0
                        } else if iz == 13 {
                            27.0
                        } else if iz == 14 {
                            25.0
                        } else if iz == 20 {
                            40.0
                        } else if (22..=29).contains(&iz) {
                            40.0
                        } else if iz == 40 {
                            40.0
                        } else if iz == 41 {
                            40.0
                        } else if iz == 42 {
                            60.0
                        } else if iz == 47 {
                            60.0
                        } else if iz == 73 {
                            90.0
                        } else if iz == 74 {
                            90.0
                        } else if iz == 79 {
                            30.0
                        } else if iz == 82 {
                            25.0
                        } else {
                            25.0
                        };
                        c.io.so(&format!("\n default damage energy ={} ev", ffmt::f(self.brk, 5, 1)));
                    }

                    //--analyze dictionary and read energy grids.
                    self.hinit(c, iold, nend4, nend6, nscr);

                    //--calculate heating from neutron files.
                    let time = c.timer();
                    c.io.so(&format!("{}temp{}  {}s", " ".repeat(61), ffmt::i(itemp, 2), ffmt::f(time, 8, 1)));
                    self.nheat(c, &mut iold, &mut inew, nscr, nend4, nend6, local);

                    //--correct for gamma production.
                    if self.mgam > 0 && self.mgam != 10 && local == 0 {
                        let time = c.timer();
                        c.io.so(&format!("{}temp{}  {}s", " ".repeat(61), ffmt::i(itemp, 2), ffmt::f(time, 8, 1)));
                        self.gheat(c, &mut iold, &mut inew, nscr);
                    }

                    //--write kermas on pendf tape.
                    self.hout(c, iold);
                }
            }
        }

        //--heatr is finished.
        self.c458 = RV::new(1);
        self.cpoly = Vec::new();
        self.hpoly = Vec::new();
        self.afr = RV::new(1);
        self.anp = RV::new(1);
        self.agp = RV::new(1);
        c.atend(self.nout, 0);
        c.repoz(self.nout);
        c.repoz(self.nin);
        c.repoz(self.nendf);
        c.closz(self.nin);
        c.closz(self.nout);
        c.closz(self.nendf);
        c.closz(-iold);
        c.closz(-inew);
        c.closz(nscr);
        c.closz(nend4);
        c.closz(nend6);
        let time = c.timer();
        c.io.so(&format!("\n{}{}s\n {}*******", " ".repeat(69), ffmt::f(time, 8, 1), "**********".repeat(7)));
    }
}

/// `horder`: arrange the elements of iarray in ascending order (iarray is
/// the Fortran array, element i at index i-1).
fn horder(iarray: &mut [i32]) {
    let n = iarray.len();
    for i in 0..n {
        for j in i..n {
            if iarray[i] > iarray[j] {
                iarray.swap(i, j);
            }
        }
    }
}

impl Heatr {
    // ================================================================ hinit
    /// Analyze the dictionary and set flags; read the energy grid of the
    /// total cross section onto scratch unit iold.
    fn hinit(&mut self, c: &mut Ctx, iold: i32, nend4: i32, nend6: i32, nscr: i32) {
        let fact = 0.99999e0f64;
        let rup = 1.000000001e0f64;
        let big = 1.0e10f64;
        let zero = 0.0f64;
        let mevev = 1.0e-6f64;
        let tol = 1.0e-5f64;
        let (mut nb, mut nw) = (0i32, 0i32);
        let mut cc = RV::new(30 + SLACK);
        let nendf = self.nendf;
        let matd = self.matd;

        //--analyze dictionary
        let mut lrel = 0;
        let mut nmod = 0;
        self.i6 = 0;
        self.mt103 = 0;
        self.mt104 = 0;
        self.mt105 = 0;
        self.mt106 = 0;
        self.mt107 = 0;
        self.mt16 = 0;
        self.nmiss4 = 0;
        self.mgam = 0;
        let mut ifiss = 0;
        self.mt19 = 0;
        let mut i519 = 0;
        self.mt458 = 0;
        self.qdel = 0.0;
        self.jp = 0;
        self.jpn = 0;
        self.jpp = 0;
        self.etop = 20000000.0;
        self.ebot = 1.0;
        self.ebot /= 100000.0;
        self.etabmax = -1.0;
        c.repoz(nendf);
        c.tpidio(nendf, 0, 0, &mut cc.all(), &mut nb, &mut nw);
        c.findf(matd, 1, 451, nendf);
        c.contio(nendf, 0, 0, &mut cc.all(), &mut nb, &mut nw);
        let n6 = c.n2h;
        c.contio(nendf, 0, 0, &mut cc.all(), &mut nb, &mut nw);
        let mut nx = 0;
        if c.n1h != 0 {
            c.iverf = 4;
            nx = n6;
        } else if c.n2h == 0 {
            c.iverf = 5;
        } else {
            c.iverf = 6;
        }
        c.skiprz(nendf, -1);
        if c.iverf >= 5 {
            c.contio(nendf, 0, 0, &mut cc.all(), &mut nb, &mut nw);
        }
        if c.iverf == 6 {
            c.contio(nendf, 0, 0, &mut cc.all(), &mut nb, &mut nw);
        }
        if c.iverf == 6 {
            if c.c2h > self.etop {
                self.etop = c.c2h;
            }
            lrel = c.l1h;
            nmod = c.n2h;
        }
        let mut scr = rv(c.npage + 50);
        c.hdatio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
        while nb != 0 {
            c.moreio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
        }
        if c.iverf >= 5 {
            nx = c.n2h;
        }
        nw = 6 * nx;
        if nw < c.npage + 50 {
            nw = c.npage + 50;
        }
        scr = rv(nw);
        nw = nx;
        if nw != 0 {
            c.dictio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
            let mut i = 1;
            let mut ifini = 0;
            while i <= nw && ifini == 0 {
                let mfd = nint(scr[i + 2]);
                let mtd = nint(scr[i + 3]);
                if mfd != 2 {
                    if mfd == 12 && mtd != 460 {
                        self.mgam = 1;
                    }
                    if mfd == 13 {
                        self.mgam = 1;
                    }
                    if mfd >= 13 {
                        ifini = 1;
                    } else {
                        if mfd <= 3 {
                            if mfd == 1 && mtd == 458 {
                                self.mt458 = 1;
                            }
                            if mtd == 19 {
                                self.mt19 = 1;
                            }
                            if mtd == 18 || mtd == 19 {
                                ifiss = 1;
                            }
                            let (lo, hi) = if c.iverf >= 6 { (600, 850) } else { (700, 800) };
                            if c.iverf >= 6 {
                                if (600..650).contains(&mtd) {
                                    self.mt103 = 1;
                                }
                                if (650..700).contains(&mtd) {
                                    self.mt104 = 1;
                                }
                                if (700..750).contains(&mtd) {
                                    self.mt105 = 1;
                                }
                                if (750..800).contains(&mtd) {
                                    self.mt106 = 1;
                                }
                                if (800..850).contains(&mtd) {
                                    self.mt107 = 1;
                                }
                                if (875..891).contains(&mtd) {
                                    self.mt16 = 1;
                                }
                            } else {
                                if (700..720).contains(&mtd) {
                                    self.mt103 = 1;
                                }
                                if (720..740).contains(&mtd) {
                                    self.mt104 = 1;
                                }
                                if (740..760).contains(&mtd) {
                                    self.mt105 = 1;
                                }
                                if (760..780).contains(&mtd) {
                                    self.mt106 = 1;
                                }
                                if (780..800).contains(&mtd) {
                                    self.mt107 = 1;
                                }
                            }
                            if mtd >= lo && mtd <= hi {
                                let mut j = i;
                                let mut idone = 0;
                                while j <= nw && idone == 0 {
                                    let mf4 = nint(scr[j + 2]);
                                    let mt4 = nint(scr[j + 3]);
                                    if mf4 == 4 && mt4 == mtd {
                                        idone = 1;
                                    } else if mf4 == 6 && mt4 == mtd {
                                        idone = 1;
                                    } else {
                                        j += 6;
                                    }
                                }
                                if idone == 0 {
                                    let strng = format!("mf4 and 6 missing, isotropy assumed for mt {}", ffmt::i(mtd, 3));
                                    c.mess("hinit", &strng, " ");
                                    self.nmiss4 += 1;
                                    self.miss4[self.nmiss4 as usize] = mtd;
                                }
                            }
                        } else {
                            if mfd == 5 {
                                if mtd == 19 {
                                    i519 = 1;
                                }
                            } else if mfd == 6 {
                                self.i6 += 1;
                                if self.i6 > MAXMF6 as i32 {
                                    c.error("hinit", "too many mf6 reactions", " ");
                                }
                                self.mt6[self.i6 as usize] = mtd;
                            }
                        }
                    }
                }
                i += 6;
            }
        }
        if self.mt19 == 1 && i519 == 0 {
            self.mt19 = 2;
        }
        if self.mt19 == 1 {
            c.mess("hinit", "mt18 is redundant", "19, 20, 21, ... will be used.");
        }
        if self.mt19 == 2 {
            c.mess("hinit", "mt19 has no spectrum", "mt18 spectrum will be used.");
        }

        //--calculate q change for fissionables in mt 458
        if ifiss != 0 {
            nw = 10000;
            scr = rv(nw);
            self.c458 = RV::new(1);
            self.cpoly = Vec::new();
            self.hpoly = Vec::new();
            self.afr = RV::new(1);
            self.anp = RV::new(1);
            self.agp = RV::new(1);
            if self.mt458 == 1 {
                c.findf(matd, 1, 458, nendf);
                c.contio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                self.lfc = nint(scr[4]);
                let nfc = nint(scr[6]);
                self.ifc1 = 0;
                self.ifc2 = 0;
                self.ifc3 = 0;
                self.ifc4 = 0;
                self.ifc5 = 0;
                self.ifc6 = 0;
                let mut l = 1;
                c.listio(nendf, 0, 0, &mut scr.at(l), &mut nb, &mut nw);
                l += nw;
                while nb != 0 {
                    c.moreio(nendf, 0, 0, &mut scr.at(l), &mut nb, &mut nw);
                    l += nw;
                }
                self.nply = nint(scr[4]);
                if self.lfc == 0 {
                    self.cpoly = vec![0.0; self.nply.max(0) as usize + 1];
                    self.hpoly = vec![0.0; self.nply.max(0) as usize + 1];
                }
                self.c458 = rv(nint(scr[5]));
                // c0 and c1 terms are fine; for ENDF/B-VII.1 the c2 and
                // higher terms are given as if energy were in MeV
                c.io.so("\n fission energy components\n");
                for k in 1..=18 {
                    self.c458[k] = scr[k + 6];
                }
                if self.lfc == 0 {
                    self.qdel = self.c458[5] + self.c458[9] + self.c458[11];
                    if self.nply >= 1 {
                        for k in 19..=36 {
                            self.c458[k] = scr[k + 6];
                        }
                        if self.nply > 1 {
                            for n in 2..=self.nply {
                                let efix = if nmod != 7 || lrel != 1 { 1.0 } else { powi(mevev, n - 1) };
                                for k in 0..18 {
                                    self.c458[1 + 18 * n + k] = scr[7 + 18 * n + k] * efix;
                                }
                            }
                        }
                        c.io.so(&format!(
                            "   fission products: polynomial of order {}\n   prompt neutrons : polynomial of order {}\n   prompt gammas   : polynomial of order {}",
                            ffmt::i(self.nply, 2),
                            ffmt::i(self.nply, 2),
                            ffmt::i(self.nply, 2)
                        ));
                    } else {
                        c.io.so("   fission products: thermal point\n   prompt neutrons : thermal point\n   prompt gammas   : thermal point");
                    }
                    for i in 0..=self.nply {
                        self.cpoly[i as usize] = self.c458[1 + i * 18] + self.c458[3 + i * 18] + self.c458[7 + i * 18];
                        self.hpoly[i as usize] = self.c458[1 + i * 18] + self.c458[7 + i * 18];
                    }
                } else if self.lfc == 1 {
                    self.qdel = 0.0;
                    for _i in 1..=nfc {
                        let mut l = 1;
                        c.tab1io(nendf, 0, 0, &mut scr.at(l), &mut nb, &mut nw);
                        while nb != 0 {
                            l += nw;
                            c.moreio(nendf, 0, 0, &mut scr.at(l), &mut nb, &mut nw);
                        }
                        let ifc = c.l2h;
                        let nl = 6 + 2 * c.n1h + 2 * c.n2h;
                        if self.etabmax < 0.0 {
                            self.etabmax = scr[nl - 1];
                        } else if ((self.etabmax - scr[nl - 1]).abs() / self.etabmax) > tol {
                            let strng = format!("upper energy mismatch for ifc={} in mt=458.", ffmt::i(ifc, 2));
                            c.error("hinit", &strng, " ");
                        }
                        let copy = |scr: &RV| {
                            let mut a = rv(nl);
                            for k in 1..=nl {
                                a[k] = scr[k];
                            }
                            a
                        };
                        if ifc == 1 {
                            self.ifc1 = 1;
                            self.afr = copy(&scr);
                        } else if ifc == 2 {
                            self.ifc2 = 1;
                            self.anp = copy(&scr);
                        } else if ifc == 3 {
                            self.ifc3 = 1;
                            self.qdel += scr[6 + 2 * c.n1h + 2];
                        } else if ifc == 4 {
                            self.ifc4 = 1;
                            self.agp = copy(&scr);
                        } else if ifc == 5 {
                            self.ifc5 = 1;
                            self.qdel += scr[6 + 2 * c.n1h + 2];
                        } else if ifc == 6 {
                            self.ifc6 = 1;
                            self.qdel += scr[6 + 2 * c.n1h + 2];
                        }
                    }
                    if self.ifc3 == 0 {
                        self.qdel += self.c458[5];
                    }
                    if self.ifc5 == 0 {
                        self.qdel += self.c458[9];
                    }
                    if self.ifc6 == 0 {
                        self.qdel += self.c458[11];
                    }
                    c.io.so(if self.ifc1 == 1 {
                        "   fission products: tabulated values"
                    } else {
                        "   fission products: thermal point"
                    });
                    c.io.so(if self.ifc2 == 1 {
                        "   prompt neutrons : tabulated values"
                    } else {
                        "   prompt neutrons : thermal point"
                    });
                    c.io.so(if self.ifc4 == 1 {
                        "   prompt gammas   : tabulated values"
                    } else {
                        "   prompt gammas   : thermal point"
                    });
                    if self.etabmax < 0.0 {
                        c.error("nheat", "no tabulated fission q components found in mt=458.", " ");
                    }
                } else {
                    c.error("hinit", "bad LFC in mt=458.", " ");
                }
            } else {
                self.lfc = 0;
                self.nply = 0;
                self.cpoly = vec![0.0; 1];
                c.mess("hinit", "mt458 is missing for this mat", " ");
            }
        }

        //--prepare a scratch tape for file 6
        //--check for missing recoil subsections
        //--and check for photon data
        self.i6p = 0;
        self.mt6yp = [0; MAXMF6 + 1];
        // set before use whenever a section is processed
        let mut ii6 = 0usize;
        if self.i6 > 0 {
            c.repoz(nend6);
            c.nsp = 1;
            c.math = 1;
            c.afend(nend6, 0);
            c.findf(matd, 6, 0, nendf);
            ii6 = 0;
        }
        c.mfh = 6;
        while c.mfh == 6 {
            c.contio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
            if c.mth == 18 {
                self.jp = c.l1h;
                self.jpn = self.jp % 10;
                self.jpp = self.jp - 10 * self.jpn;
            }
            if c.mth != 18 || (c.mth == 18 && self.jp == 0) {
                c.contio(0, nend6, 0, &mut scr.all(), &mut nb, &mut nw);
            } else {
                //--delete mt=18 from mt6 array and skip this section (for now)
                let i6t = self.i6;
                for i in 1..=i6t as usize {
                    if self.mt6[i] < 18 {
                        continue;
                    }
                    if self.mt6[i] == 18 {
                        self.i6 -= 1;
                    }
                    if self.mt6[i] > 18 && self.i6 != i6t {
                        self.mt6[i - 1] = self.mt6[i];
                    }
                }
                c.tosend(nendf, 0, 0, &mut scr.all());
                continue;
            }
            if c.mfh != 0 && c.mth != 0 {
                let nk = c.n1h;
                ii6 += 1;
                self.mt6no[ii6] = 0;
                let mut zar = c.c1h + 1.0;
                let mut awrr = c.c2h + 1.0;
                let ielem = nint(c.c1h) % 1000;
                for ik in 1..=nk {
                    let mut l = 1;
                    c.tab1io(nendf, nend6, 0, &mut scr.at(l), &mut nb, &mut nw);
                    l += nw;
                    while nb != 0 {
                        c.moreio(nendf, nend6, 0, &mut scr.at(l), &mut nb, &mut nw);
                        l += nw;
                    }
                    let nr = c.n1h;
                    self.zap = c.c1h;
                    let law = c.l2h;
                    if nint(self.zap) == 0 {
                        self.i6p += 1;
                        self.mt6yp[self.i6p as usize] = c.mth;
                    }
                    let lf = c.l1h;
                    if c.mth == 5 {
                        zar = 0.0;
                    }
                    if c.mth != 5 {
                        let yld = scr[6 + 2 * nr + 2];
                        zar -= self.zap * yld;
                        awrr -= c.c2h * yld;
                        if self.zap == zero && (ik > 1 || nk == 1) {
                            if nk != 1 || self.zap != zero {
                                if zar != zero {
                                    let test = 3000.0;
                                    if ielem != 0 || zar.abs() >= test {
                                        let strng = format!(
                                            "mf6, mt{} does not give recoil za={}",
                                            ffmt::i(c.mth, 3),
                                            ffmt::i(nint(zar), 6)
                                        );
                                        if c.mth == 102 {
                                            c.mess("hinit", &strng, "photon momentum recoil used.");
                                        } else if c.mth != 10 {
                                            c.mess("hinit", &strng, "one-particle recoil approx. used.");
                                        }
                                        self.mt6no[ii6] = ik - 1;
                                        zar = 0.0;
                                    }
                                }
                            }
                        }
                    }
                    if c.mth == 102 && self.zap > 0.0 && lf == 0 {
                        let strng = format!("mf6, mt{} has recoil with no spectrum", ffmt::i(c.mth, 3));
                        c.mess("hinit", &strng, "photon momentum recoil used.");
                        self.mt6no[ii6] = nk;
                    }
                    if self.zap == zero {
                        self.mgam = 10 + self.mgam % 10;
                    }
                    if law == 6 {
                        c.contio(nendf, nend6, 0, &mut scr.all(), &mut nb, &mut nw);
                    } else if law == 1 || law == 2 || law == 5 {
                        c.tab2io(nendf, nend6, 0, &mut scr.all(), &mut nb, &mut nw);
                        self.ne = c.n2h;
                        for _ie in 1..=self.ne {
                            c.listio(nendf, nend6, 0, &mut scr.all(), &mut nb, &mut nw);
                            while nb != 0 {
                                c.moreio(nendf, nend6, 0, &mut scr.all(), &mut nb, &mut nw);
                            }
                        }
                    } else if law == 7 {
                        c.tab2io(nendf, nend6, 0, &mut scr.all(), &mut nb, &mut nw);
                        self.ne = c.n2h;
                        for _ie in 1..=self.ne {
                            c.tab2io(nendf, nend6, 0, &mut scr.all(), &mut nb, &mut nw);
                            let nmu = c.n2h;
                            for _imu in 1..=nmu {
                                c.tab1io(nendf, nend6, 0, &mut scr.all(), &mut nb, &mut nw);
                                while nb != 0 {
                                    c.moreio(nendf, nend6, 0, &mut scr.all(), &mut nb, &mut nw);
                                }
                            }
                        }
                    } else if law < 0 {
                        c.skip6(nendf, nend6, 0, &mut scr.all(), law);
                    }
                }
                let _ = awrr;
                if zar != zero {
                    let test = 3000.0;
                    if ielem != 0 || zar.abs() >= test {
                        let strng = format!("mf6, mt{} does not give recoil za={}", ffmt::i(c.mth, 3), ffmt::i(nint(zar), 6));
                        if c.mth == 102 {
                            c.mess("hinit", &strng, "photon momentum recoil used.");
                        } else if c.mth != 10 {
                            c.mess("hinit", &strng, "one-particle recoil approx. used.");
                        }
                        self.mt6no[ii6] = nk;
                    }
                }
            }
        }
        c.amend(nend6, 0);
        c.atend(nend6, 0);
        c.repoz(nend6);
        if self.mgam == 0 {
            c.io.so("\n no photon production files...all photon energy will be deposited locally.");
        }

        //--copy endf tape to nscr for use in conbar and gheat
        //--and to nend4 for use in disbar
        //--convert photon transition probability arrays if present
        c.repoz(nendf);
        c.repoz(nend4);
        c.repoz(nscr);
        c.nsp = 1;
        c.nsc = 1;
        c.math = 1;
        c.afend(nend4, nscr);
        c.findf(matd, 1, 451, nendf);
        self.hconvr(c, nendf, nend4, nscr);
        c.amend(nend4, nscr);
        c.atend(nend4, nscr);
        c.repoz(nend4);
        c.repoz(nscr);
        c.repoz(nendf);

        //--store energy grid of the total cross section.
        c.findf(matd, 3, 1, self.nin);
        let mut npkk = self.npk;
        if self.kchk == 1 {
            npkk = 3 * self.npk - 2;
        }
        if self.mgam > 0 {
            npkk += 3;
        }
        if self.mgam == 0 && self.i6 > 0 {
            npkk += 1;
        }
        for i in 1..=npkk {
            cc[i] = 0.0;
        }
        c.contio(self.nin, 0, 0, &mut scr.all(), &mut nb, &mut nw);
        let mut e = 0.0f64;
        let mut enext = 0.0f64;
        let mut idis = 0i32;
        let mut y = 0.0f64;
        c.gety1(e, &mut enext, &mut idis, &mut y, self.nin, &mut scr.all());
        enext *= rup;
        self.efirst = enext;
        let mut idnx = 0;
        let mut bufo = vec![0.0f64; NBUF as usize];
        self.ne = 0;
        let mut nen = 1;
        while nen > 0 {
            self.ne += 1;
            let test = enext * fact;
            if idnx > 0 && test > e {
                enext = test;
            }
            e = enext;
            c.gety1(e, &mut enext, &mut idis, &mut y, self.nin, &mut scr.all());
            nen = self.ne;
            if enext > big - big / 100.0 {
                nen = -nen;
            }
            cc[1] = e;
            c.loada(nen, &cc.0, npkk, iold, &mut bufo, NBUF);
            if idis > 0 && idnx > 0 {
                idnx = 0;
            }
            if idis > 0 && idnx == 0 {
                idnx = 1;
            }
            if idis == 0 {
                idnx = 0;
            }
        }
        self.elast = e;
        c.tosend(self.nin, 0, 0, &mut scr.all());
        self.mt303 = 0;
        if self.npk < 3 {
            return;
        }
        for i in 3..=self.npk as usize {
            if self.mtp[i] == 303 {
                self.mt303 = i as i32;
            }
        }
    }
}

/// Fission-type reaction (mt 18-21 or 38).
fn is_fis(mt: i32) -> bool {
    (18..=21).contains(&mt) || mt == 38
}

impl Heatr {
    // ================================================================ nheat
    /// Contributions to heating from the neutron files. iold and inew are
    /// exchanged as in the Fortran (they are arguments by reference).
    fn nheat(&mut self, c: &mut Ctx, iold: &mut i32, inew: &mut i32, nscr: i32, nend4: i32, nend6: i32, local: i32) {
        let small = 1.0e-10f64;
        let big = 1.0e10f64;
        let efis = 15.0e6f64;
        let fq1 = 8.07e6f64;
        let fq2 = 1.307e0f64;
        let one = 1.0f64;
        let ten = 10.0f64;
        let qtest = 99.0e6f64;
        let smin = 1.0e-9f64;
        let break1 = 1.0e6f64;
        let break2 = 2.0e6f64;
        let break3 = 5.0e6f64;
        let break4 = 20.0e6f64;
        let step1 = 0.2e6f64;
        let step2 = 0.5e6f64;
        let step3 = 1.0e6f64;
        let step4 = 5.0e6f64;
        let up = 1.1e0f64;
        let zero = 0.0f64;
        let tol = 1.0e-5f64;
        let (mut nb, mut nw) = (0i32, 0i32);
        let matd = self.matd;
        let nin = self.nin;
        let nendf = self.nendf;

        //--allocate storage
        let mut b = rv(c.npage + 50);
        let nwmax: i32 = 250000;
        let mut d = rv(nwmax);
        let na: i32 = 500000;
        let mut a = rv(na);
        let mut scr = rv(c.npage + 50);
        let mut bufo = vec![0.0f64; NBUF as usize];
        let mut bufn = vec![0.0f64; NBUF as usize];
        let mut cv = RV::new(30 + SLACK);
        let mut imt = [0i32; 31];
        let awfac = 1.0 / self.awr;
        let aw1fac = 1.0 / (self.awr + 1.0);
        let mut npkk = self.npk;
        if self.kchk == 1 {
            npkk = 3 * self.npk - 2;
        }
        if self.mgam > 0 {
            npkk += 3;
        }
        if self.mgam == 0 && self.i6 > 0 {
            npkk += 1;
        }
        let mut npkkk = 0;
        for ipk in 3..=self.npk {
            if self.mtp[ipk as usize] == 442 {
                npkkk = ipk;
            }
        }
        // check for requested dame energy calculation
        self.idame = 0;
        let mut ipk = 2;
        while self.idame == 0 && ipk < self.npk {
            ipk += 1;
            if self.mtp[ipk as usize] >= 444 {
                self.idame = 1;
            }
        }
        let mut n6 = 0i32;
        let mut j6 = 0i32;
        self.i6g = 0;
        let mut irec = 0i32;
        let mut jrec = 0i32;
        let mut last6 = 0;
        let mut new6 = 0;
        let mut qsave = 0.0f64;
        let mut pnue = 0.0f64;

        // locals of the Fortran routine (kept across the labels)
        let mut e = 0.0f64;
        let mut thresh = 0.0f64;
        let mut idis = 0i32;
        let mut y = 0.0f64;
        let mut t = 0.0f64;
        let mut qendf = 0.0f64;
        let mut lr = 0i32;
        let mut qs = 0.0f64;
        let mut q0 = 0.0f64;
        let mut h0 = 0.0f64;
        let mut icon = 0i32;
        let mut ebar = 0.0f64;
        let mut dame = 0.0f64;
        let mut mtd = 0i32;
        let mut iimt = 0usize;
        let mut yld = 0.0f64;
        let mut nwa = na;
        let mut nwm = nwmax;
        let mut iflag = 0i32;
        let mut nmt = 0i32;
        let mut ipx = 2;
        let mut irx = 1;
        let mut ipfr = 2;
        let mut irfr = 1;
        let mut ipnp = 2;
        let mut irnp = 1;
        let mut ipgp = 2;
        let mut irgp = 1;
        let mut elst = 0.0f64;
        let mut ilist = 1i32;
        let mut yld0 = 0.0f64;
        let mut nn = 0i32;
        let mut enext = 0.0f64;
        let mut enx = 0.0f64;
        let mut idx = 0i32;
        let mut h = 0.0f64;
        let mut ebal6 = 0.0f64;
        let mut hmin = 0.0f64;
        let mut hmax = 0.0f64;
        let mut qfr = 0.0f64;
        let mut qnp = 0.0f64;
        let mut qgp = 0.0f64;

        let mut lbl = 105;
        loop {
            match lbl {
                105 => {
                    //--loop over non-redundant mt-s in file 3.
                    c.contio(nin, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                    if c.mfh == 0 {
                        lbl = 300;
                        continue;
                    }
                    if c.mth == 0 {
                        continue;
                    }
                    let mth = c.mth;
                    let v5 = c.iverf <= 5;
                    let skip = mth == 3
                        || mth == 4
                        || mth == 10
                        || mth == 26
                        || mth == 27
                        || (mth == 18 && self.mt19 == 1)
                        || ((19..=21).contains(&mth) && self.mt19 == 2)
                        || (mth == 38 && self.mt19 == 2)
                        || mth == 101
                        || (mth == 103 && self.mt103 > 0)
                        || (mth == 104 && self.mt104 > 0)
                        || (mth == 105 && self.mt105 > 0)
                        || (mth == 106 && self.mt106 > 0)
                        || (mth == 107 && self.mt107 > 0)
                        || (mth == 16 && self.mt16 > 0)
                        || (v5 && mth == 719)
                        || (v5 && mth == 739)
                        || (v5 && mth == 759)
                        || (v5 && mth == 779)
                        || (v5 && mth == 799)
                        || (mth > 120 && mth < 152)
                        || (mth > 200 && mth < 600);
                    if skip {
                        c.tosend(nin, 0, 0, &mut scr.all());
                        continue;
                    }
                    self.zat = c.c1h;
                    self.awrt = c.c2h;
                    self.izat = nint(self.zat);
                    lbl = 120;
                }
                120 => {
                    //--compute heating from this section.
                    e = 0.0;
                    c.gety1(e, &mut thresh, &mut idis, &mut y, nin, &mut b.all());
                    thresh = sigfig(thresh, 7, 0);
                    t = c.c1h;
                    self.q = c.c2h;
                    qendf = self.q;
                    lr = c.l2h;
                    qs = 0.0;
                    q0 = 0.0;
                    h0 = 0.0;
                    icon = 0;
                    ebar = 0.0;
                    dame = 0.0;
                    mtd = c.mth;

                    //--check for an endf-6 file 6
                    let mut go123 = j6 > 0;
                    if !go123 {
                        let mut found = false;
                        for i in 1..=self.i6 as usize {
                            if mtd == self.mt6[i] {
                                found = true;
                                break;
                            }
                        }
                        if found {
                            // 121
                            c.findf(matd, 6, mtd, nend6);
                            c.contio(nend6, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                            n6 = nint(scr[5]);
                            self.lct = nint(scr[4]);
                            self.jp = nint(scr[3]);
                            self.jpn = self.jp % 10;
                            self.jpp = self.jp - 10 * self.jpn;
                            if c.mth == 18 && self.jpn != 0 {
                                c.tosend(nend6, 0, 0, &mut scr.all());
                                lbl = 105;
                                continue;
                            }
                            self.i6g += 1;
                            go123 = true;
                        }
                    }
                    if go123 {
                        // 123
                        j6 += 1;
                        if irec > 0 && self.mt6no[self.i6g as usize] > 0 {
                            j6 -= 1;
                        }
                        last6 = 0;
                        if j6 == n6 && irec > 0 {
                            last6 = 1;
                        }
                        if j6 == n6 && self.mt6no[self.i6g as usize] == 0 {
                            last6 = 1;
                        }
                        new6 = 0;
                        if j6 == 1 && irec == 0 {
                            new6 = 1;
                        }
                    }
                    // 124
                    //--adjust fission q to e(in)=0 prompt value
                    if is_fis(c.mth) {
                        if self.mt458 == 1 {
                            qendf = self.c458[15];
                            self.q = qendf - self.qdel;
                            let strng1 = format!("changed q from {} to {}", ffmt::e1p(qendf, 14, 6), ffmt::e1p(self.q, 14, 6));
                            let strng2 = format!("for mt {} by taking out delayed components", ffmt::i(c.mth, 3));
                            c.mess("nheat", &strng1, &strng2);
                        } else {
                            self.q = qendf;
                            self.cpoly[0] = self.q;
                        }
                    }

                    //--choose appropriate ground state q value.
                    iimt = 0;
                    if self.nqa != 0 {
                        for i in 1..=self.nqa as usize {
                            if self.mta[i] != c.mth {
                                continue;
                            }
                            qs = self.qa[i];
                            iimt = i;
                            break;
                        }
                    }
                    // 140
                    lbl = 140;
                }
                140 => {
                    if j6 > 0 {
                        lbl = 175;
                        continue;
                    }
                    if c.mth > 15 && c.mth < 51 {
                        // 150
                        q0 = self.q;
                        if iimt != 0 {
                            q0 = qs;
                        }
                        lbl = 160;
                        continue;
                    }
                    if c.mth > 100 {
                        lbl = 170;
                        continue;
                    }
                    icon = 1;
                    q0 = 0.0;
                    if lr != 0 && lr != 31 {
                        q0 = t;
                    }
                    if iimt != 0 {
                        q0 = qs;
                    }
                    if c.mth == 91 {
                        lbl = 160;
                        continue;
                    }
                    yld = 1.0;
                    // accomodate neutron producing pseudo-levels
                    if lr == 16 || lr == 21 || lr == 24 || lr == 26 || lr == 30 {
                        yld = 2.0;
                    }
                    if lr == 17 || lr == 25 || lr == 38 {
                        yld = 3.0;
                    }
                    if lr == 37 {
                        yld = 4.0;
                    }
                    nwa = na;
                    let (q, za, awr) = (self.q, self.za, self.awr);
                    self.disbar(c, e, &mut ebar, &mut dame, q, za, awr, nend4, matd, mtd, &mut a, nwa);
                    lbl = 180;
                }
                160 => {
                    icon = 2;
                    nwm = nwmax;
                    nwa = na;
                    let q = self.q;
                    self.conbar(c, e, &mut ebar, &mut dame, nendf, nscr, matd, mtd, &mut a, nwa, &mut d, nwm, &mut yld, q);
                    lbl = 180;
                }
                170 => {
                    let iverf = c.iverf;
                    if (iverf <= 5 && (700..=799).contains(&c.mth)) || (iverf >= 6 && (600..=849).contains(&c.mth)) {
                        // 171
                        let m = c.mth;
                        let go174 = (iverf <= 5 && (m == 718 || m == 738 || m == 758 || m == 778 || m == 798))
                            || (iverf >= 6 && (m == 649 || m == 699 || m == 749 || m == 799 || m == 849));
                        if go174 {
                            // 174: continuum charged-particle emission
                            yld = 1.0;
                            if iverf <= 5 {
                                q0 = qsave;
                            }
                            if iverf >= 6 {
                                q0 = t;
                            }
                            if self.idame > 0 {
                                let (q, za, awr) = (self.q, self.za, self.awr);
                                self.capdam(c, e, &mut dame, q, za, awr, mtd);
                            }
                            lbl = 180;
                            continue;
                        }
                        // discrete charged-particle levels
                        icon = 1;
                        if iverf <= 5 && (m == 700 || m == 720 || m == 740 || m == 760 || m == 780) {
                            qsave = self.q;
                        }
                        if iverf <= 5 {
                            q0 = qsave;
                        }
                        if iverf >= 6 {
                            q0 = t;
                        }
                        yld = 1.0;
                        nwa = na;
                        let (q, za, awr) = (self.q, self.za, self.awr);
                        self.disbar(c, e, &mut ebar, &mut dame, q, za, awr, nend4, matd, mtd, &mut a, nwa);
                        lbl = 180;
                        continue;
                    }
                    q0 = self.q;
                    if self.idame > 0 {
                        let (q, za, awr) = (self.q, self.za, self.awr);
                        self.capdam(c, e, &mut dame, q, za, awr, mtd);
                    }
                    lbl = 180;
                }
                175 => {
                    // file 6 reactions
                    nwm = nwmax;
                    nwa = na;
                    self.sixbar(
                        c, e, &mut ebar, &mut yld, &mut dame, nend6, &mut a, nwa, nscr, &mut d, nwm, n6, j6, &mut irec,
                        &mut jrec, &mut iflag,
                    );
                    if ebar < zero {
                        lbl = 291;
                        continue;
                    }
                    icon = -1;
                    q0 = t;
                    if iimt != 0 {
                        q0 = qs;
                    }
                    if is_fis(mtd) {
                        q0 = self.q;
                    }
                    lbl = 180;
                }
                180 => {
                    nmt = 0;
                    if self.npk > 2 {
                        indx(self.npk, &self.mtp, self.mt19, mtd, &mut imt, &mut nmt, c.iverf);
                    }
                    let dmg = |s: &str| -> String { s.to_string() };
                    let _ = dmg;
                    if icon < 0 {
                        // 179
                        self.izap = nint(self.zap);
                        if irec > 0 {
                            if self.izat % 1000 == 0 {
                                self.izap = 1000 * (self.izat / 1000 - self.izap / 1000);
                            } else {
                                self.izap = nint(self.zat + 1.0 - self.zap);
                            }
                        }
                        if self.iprint == 1 {
                            let tail = if self.idame > 0 { format!("{}damage", " ".repeat(8)) } else { String::new() };
                            if qs >= qtest {
                                // 178
                                c.io.so(&format!(
                                    "\n file six heating for mt{}, particle ={}{}q0 = variable\n{}e{}qbar{}ebar{}yield{}xsec{}heating{}",
                                    ffmt::i(mtd, 3),
                                    ffmt::i(self.izap, 6),
                                    " ".repeat(5),
                                    " ".repeat(14),
                                    " ".repeat(10),
                                    " ".repeat(10),
                                    " ".repeat(9),
                                    " ".repeat(10),
                                    " ".repeat(7),
                                    tail
                                ));
                            } else {
                                c.io.so(&format!(
                                    "\n file six heating for mt{}, particle ={}{}q = {}\n{}e{}ebar{}yield{}xsec{}heating{}",
                                    ffmt::i(mtd, 3),
                                    ffmt::i(self.izap, 6),
                                    " ".repeat(5),
                                    ffmt::e1p(q0, 12, 4),
                                    " ".repeat(14),
                                    " ".repeat(10),
                                    " ".repeat(9),
                                    " ".repeat(10),
                                    " ".repeat(7),
                                    tail
                                ));
                            }
                        }
                    } else if self.iprint == 1 {
                        let tail = if self.idame > 0 { format!("{}damage", " ".repeat(8)) } else { String::new() };
                        if qs >= qtest {
                            // 181
                            c.io.so(&format!(
                                "\n neutron heating for mt{}   q0= variable{}q ={}\n{}e{}qbar{}ebar{}yield{}xsec{}heating{}",
                                ffmt::i(mtd, 3),
                                " ".repeat(5),
                                ffmt::e(qendf, 12, 4),
                                " ".repeat(14),
                                " ".repeat(10),
                                " ".repeat(10),
                                " ".repeat(9),
                                " ".repeat(10),
                                " ".repeat(7),
                                tail
                            ));
                        } else if is_fis(mtd) {
                            // 183
                            c.io.so(&format!(
                                "\n neutron heating for mt{}   q0= variable fission q\n{}e{}q0{}ebar{}yield{}xsec{}heating{}",
                                ffmt::i(mtd, 3),
                                " ".repeat(14),
                                " ".repeat(12),
                                " ".repeat(10),
                                " ".repeat(9),
                                " ".repeat(10),
                                " ".repeat(7),
                                tail
                            ));
                            if irec > 0 && icon < 0 {
                                self.izap = 100;
                            }
                        } else {
                            c.io.so(&format!(
                                "\n neutron heating for mt{}   q0 ={}{}q ={}\n{}e{}ebar{}yield{}xsec{}heating{}",
                                ffmt::i(mtd, 3),
                                ffmt::e1p(q0, 12, 4),
                                " ".repeat(5),
                                ffmt::e1p(qendf, 12, 4),
                                " ".repeat(14),
                                " ".repeat(10),
                                " ".repeat(9),
                                " ".repeat(10),
                                " ".repeat(7),
                                tail
                            ));
                        }
                    }
                    // 182
                    ipx = 2;
                    irx = 1;
                    ipfr = 2;
                    irfr = 1;
                    ipnp = 2;
                    irnp = 1;
                    ipgp = 2;
                    irgp = 1;
                    elst = 0.0;
                    ilist = 1;
                    if icon == 0 {
                        yld = 0.0;
                    }
                    yld0 = 0.0;

                    //--compute heating for energies in union grid.
                    nn = 0;
                    lbl = 190;
                }
                190 => {
                    nn += 1;
                    c.finda(nn, &mut cv.0, npkk, *iold, &mut bufo, NBUF);
                    e = cv[1];
                    e = sigfig(e, 9, 0);
                    if new6 == 1 {
                        cv[npkk] = 0.0;
                    }
                    if mtd == 2 && !(e < elst && nn < self.ne) {
                        if ilist > ILMAX {
                            c.error("nheat", "storage exceeded.", " ");
                        }
                        self.elist[ilist as usize] = e;
                        self.elist[(ilist + 1) as usize] = big;
                        if e < break1 {
                            let mut iii = int(log10(e * (1.0 + small)) + small);
                            if log10(e * (1.0 + small)) < zero {
                                iii -= 1;
                            }
                            let xxx = powi(ten, iii);
                            if 2.0 * xxx > e * (1.0 + small) && e >= one * (1.0 - small) {
                                elst = 2.0 * xxx;
                            } else if 5.0 * xxx > e * (1.0 + small) && e >= one * (1.0 - small) {
                                elst = 5.0 * xxx;
                            } else {
                                elst = 10.0 * xxx;
                            }
                        } else {
                            let mut test = break4;
                            test = sigfig(test, 9, 1);
                            if e >= test {
                                elst = step4 * int(e / step4 + up) as f64;
                            }
                            if e < test {
                                elst = step3 * int(e / step3 + up) as f64;
                            }
                            test = break3;
                            test = sigfig(test, 9, 1);
                            if e < test {
                                elst = step2 * int(e / step2 + up) as f64;
                            }
                            test = break2;
                            test = sigfig(test, 9, 1);
                            if e < test {
                                elst = step1 * int(e / step1 + up) as f64;
                            }
                        }
                        elst -= elst / 10000000.0;
                    }
                    // 193
                    if e < thresh {
                        lbl = 290;
                        continue;
                    }
                    c.gety1(e, &mut enext, &mut idis, &mut y, nin, &mut b.all());
                    // check for energy-dependent q
                    if iimt > 0 && self.qa[iimt] >= qtest {
                        let lq0 = self.lqs[iimt];
                        c.terpa(&mut q0, e, &mut enx, &mut idx, &self.qbar.at(lq0), &mut ipx, &mut irx);
                        self.q = q0;
                    }
                    let (q, za, awr) = (self.q, self.za, self.awr);
                    if icon == 1 {
                        self.disbar(c, e, &mut ebar, &mut dame, q, za, awr, nend4, matd, mtd, &mut a, nwa);
                    }
                    if icon == 2 {
                        self.conbar(c, e, &mut ebar, &mut dame, nendf, nscr, matd, mtd, &mut a, nwa, &mut d, nwm, &mut yld, q);
                    }
                    if icon == 0 && self.idame > 0 {
                        self.capdam(c, e, &mut dame, q, za, awr, mtd);
                    }
                    if icon < 0 {
                        self.sixbar(
                            c, e, &mut ebar, &mut yld, &mut dame, nend6, &mut a, nwa, nscr, &mut d, nwm, n6, j6, &mut irec,
                            &mut jrec, &mut iflag,
                        );
                    }
                    if yld0 == zero {
                        yld0 = yld;
                    }
                    // correct fission q for current incident energy
                    if is_fis(mtd) {
                        if self.lfc == 0 {
                            if self.nply == 0 {
                                q0 = self.cpoly[0] - (yld - yld0) * fq1 + fq2 * e;
                                if self.mt458 == 1 {
                                    h0 = self.hpoly[0];
                                } else {
                                    h0 = q0 - ebar * yld;
                                }
                            } else {
                                let np = self.nply as usize;
                                q0 = self.cpoly[np];
                                h0 = self.hpoly[np];
                                let mut i = np as i32 - 1;
                                while i >= 0 {
                                    q0 = q0 * e + self.cpoly[i as usize];
                                    h0 = h0 * e + self.hpoly[i as usize];
                                    i -= 1;
                                }
                            }
                            pnue = q0 - h0;
                        } else {
                            if self.etabmax < e && ((self.etabmax - e).abs() / self.etabmax) > tol {
                                c.error("nheat", "upper energy tabulated fission q components is too low.", " ");
                            }
                            if self.ifc1 == 1 {
                                c.terpa(&mut qfr, e, &mut enx, &mut idx, &self.afr.all(), &mut ipfr, &mut irfr);
                            } else {
                                qfr = self.c458[1];
                            }
                            if self.ifc2 == 1 {
                                c.terpa(&mut qnp, e, &mut enx, &mut idx, &self.anp.all(), &mut ipnp, &mut irnp);
                            } else {
                                qnp = self.c458[3] - (yld - yld0) * fq1 + fq2 * e;
                            }
                            if self.ifc4 == 1 {
                                c.terpa(&mut qgp, e, &mut enx, &mut idx, &self.agp.all(), &mut ipgp, &mut irgp);
                            } else {
                                qgp = self.c458[7];
                            }
                            pnue = qnp;
                            h0 = qfr + qgp;
                            q0 = h0 + pnue;
                        }
                        q0 -= e;
                    }
                    if icon >= 0 {
                        ebal6 = 0.0;
                        dame *= y;
                        if !is_fis(mtd) {
                            h = (e + q0 - ebar * yld) * y;
                        } else {
                            h = h0 * y;
                        }
                    } else {
                        if !is_fis(mtd) {
                            h = ebar * yld * y;
                        } else {
                            h = h0 * y;
                        }
                        dame *= y;
                        cv[npkk] += h;
                        if self.izap == 0 {
                            cv[npkk - 1] += h;
                        }
                        if self.izap == 0 && npkkk > 0 {
                            cv[npkkk] += h;
                        }
                        if self.izap == 1 {
                            h = 0.0;
                        }
                        if self.izap == 0 && local == 0 {
                            h = 0.0;
                        }
                        if self.izap <= 1 {
                            dame = 0.0;
                        }
                        ebal6 = 0.0;
                        if last6 == 1 && c.mth != 5 && c.mth != 102 {
                            ebal6 = (e + q0) * y - cv[npkk];
                        }
                        if is_fis(mtd) {
                            ebal6 = 0.0;
                        }
                    }
                    if (46..=49).contains(&mtd) {
                        // treat second neutron from sequential (n,2n), mt=46-49.
                        if qs >= zero {
                            c.error("nheat", "binding energy for sequential n2n needed.", " ");
                        }
                        h = (qs - ebar) * y;
                    }

                    //--accumulate total and partial heating factors.
                    cv[2] = cv[2] + h + ebal6;
                    if is_fis(mtd) {
                        ebar = pnue / yld;
                    }
                    let ebarp = sigfig(ebar, 9, 0);
                    let yldp = sigfig(yld, 9, 0);
                    let yp = sigfig(y, 9, 0);
                    let hp = sigfig(h, 9, 0);
                    let damep = sigfig(dame, 9, 0);
                    let ebal6p = sigfig(ebal6, 9, 0);
                    let el = self.elist[ilist as usize];
                    if !(self.iprint != 1 || (e - el).abs() > small * e || y <= smin) {
                        if qs >= qtest || is_fis(mtd) {
                            // 198
                            if self.idame == 0 {
                                c.io.so(&e14_lines(&[e, q0, ebarp, yldp, yp, hp], 7));
                            } else {
                                c.io.so(&e14_lines(&[e, q0, ebarp, yldp, yp, hp, damep], 7));
                            }
                        } else if self.idame == 0 {
                            c.io.so(&e14_lines(&[e, ebarp, yldp, yp, hp], 7));
                        } else {
                            c.io.so(&e14_lines(&[e, ebarp, yldp, yp, hp, damep], 7));
                        }
                    }
                    // 199
                    if ebal6.abs() >= 100.0 * y {
                        if self.iprint == 1 && (e - el).abs() <= small * e && y > smin {
                            if qs < qtest {
                                c.io.so(&format!("{}ebal{}", " ".repeat(53), ffmt::e1p(ebal6p, 14, 4)));
                            } else {
                                c.io.so(&format!("{}ebal{}", " ".repeat(67), ffmt::e1p(ebal6p, 14, 4)));
                            }
                        }
                    }
                    if nmt != 0 {
                        // partial kermas
                        for ii in 1..=nmt as usize {
                            let indxx = imt[ii];
                            if self.mtp[indxx as usize] >= 444 {
                                cv[indxx] += dame;
                            }
                            if self.mtp[indxx as usize] < 442 {
                                cv[indxx] = cv[indxx] + h + ebal6;
                            }
                        }
                    }

                    //--compute kinematic limits on heating, if desired.
                    if self.kchk == 0 {
                        lbl = 290;
                        continue;
                    }
                    let iverf = c.iverf;
                    if icon < 0 {
                        // 275
                        hmin = h;
                        hmax = h;
                    } else if mtd == 2 {
                        hmin = h;
                        hmax = hmin;
                    } else if mtd == 102 {
                        // 215
                        hmin = y * e * aw1fac;
                        let tt = self.q + self.awr * e * aw1fac;
                        hmax = hmin + y * tt * tt * self.rtm / 2.0;
                    } else if (51..91).contains(&mtd)
                        || (iverf <= 5
                            && ((700..718).contains(&mtd)
                                || (720..738).contains(&mtd)
                                || (740..758).contains(&mtd)
                                || (760..778).contains(&mtd)
                                || (780..798).contains(&mtd)))
                        || (iverf >= 6
                            && ((600..649).contains(&mtd)
                                || (650..699).contains(&mtd)
                                || (700..749).contains(&mtd)
                                || (750..799).contains(&mtd)
                                || (800..849).contains(&mtd)))
                    {
                        // 221
                        hmin = (e + self.q - ebar * yld) * y;
                        hmin = (e + self.q - ebar * yld) * y;
                        hmax = hmin;
                        if lr != 0 && lr != 31 {
                            hmax = (e + q0 - ebar * yld) * y;
                        }
                    } else if mtd == 91 {
                        // 230
                        hmin = y * (e + ebar) * awfac;
                        hmax = hmin;
                        if lr != 0 && lr != 31 {
                            hmax = (e + q0 - ebar) * y;
                        }
                    } else if mtd >= 103 {
                        // 240
                        hmin = y * e * awfac;
                        hmax = h;
                    } else if mtd == 18 || mtd == 20 || mtd == 21 || mtd == 38 {
                        // 250
                        hmin = y * (e + q0 - ebar * yld / 2.0 - efis);
                        hmax = h;
                    } else if mtd == 16 && !(self.awr < 10.0) {
                        // 260
                        hmin = 0.0;
                        hmax = y * (e + ebar) / (self.awr - 1.0);
                        if hmax > h {
                            hmax = h;
                        }
                    } else if mtd == 17 && !(self.awr < 10.0) {
                        // 265
                        hmin = 0.0;
                        hmax = y * (e + 2.0 * ebar) / (self.awr - 2.0);
                        if hmax > h {
                            hmax = h;
                        }
                    } else {
                        // 270
                        hmin = 0.0;
                        hmax = h;
                    }
                    // 285
                    let npkt = self.npk - 1;
                    if self.kchk == 1 {
                        cv[2 + npkt] += hmin;
                        cv[2 + 2 * npkt] += hmax;
                    }
                    // 283
                    if nmt != 0 {
                        for ii in 1..=nmt as usize {
                            let indxx = imt[ii];
                            let m = self.mtp[indxx as usize];
                            if m == 442 {
                                continue;
                            }
                            if m == 443 {
                                cv[indxx] += hmax;
                                continue;
                            }
                            if m >= 444 {
                                continue;
                            }
                            if self.kchk != 1 {
                                continue;
                            }
                            cv[indxx + npkt] += hmin;
                            cv[indxx + 2 * npkt] += hmax;
                        }
                    }
                    lbl = 290;
                }
                290 => {
                    //--store results on scratch and go to next energy.
                    let mut intt = nn;
                    if intt == self.ne {
                        intt = -intt;
                    }
                    c.loada(intt, &cv.0, npkk, *inew, &mut bufn, NBUF);
                    if e == self.elist[ilist as usize] {
                        ilist += 1;
                    }
                    if nn < self.ne {
                        lbl = 190;
                        continue;
                    }

                    //--for file 6, loop over subsections
                    std::mem::swap(iold, inew);
                    lbl = 291;
                }
                291 => {
                    if j6 >= n6 {
                        // make sure we're at the end of section (will not be
                        // so if there was no particle distribution for the
                        // last subsection in this section).
                        if ebar < zero {
                            c.tosend(nin, 0, 0, &mut b.all());
                            c.skiprz(nin, -1);
                        }
                        // 296
                        if j6 == 0 {
                            lbl = 105;
                            continue;
                        }
                        let i6g = self.i6g as usize;
                        if !(irec > 0 && self.mt6no[i6g] > 0) {
                            if self.iprint == 1 && self.i6g > 0 && self.izap != 0 && iflag == 0 {
                                c.io.so(&format!("\n   no explicit file 6 photon production for mt{}", ffmt::i(mtd, 3)));
                            }
                            if self.mt6no[i6g] != 0 {
                                if c.mth != 102 {
                                    c.io.so("\n   generating recoil with one-particle approx.");
                                }
                                if c.mth == 102 {
                                    c.io.so("\n   generating recoil from photon momentum.");
                                }
                                irec = self.mt6no[i6g];
                                c.skiprz(nin, -2);
                                c.findf(matd, 3, mtd, nin);
                                lbl = 105;
                                continue;
                            }
                        }
                        // 293
                        irec = 0;
                        jrec = 0;
                        j6 = 0;
                        n6 = 0;
                        lbl = 105;
                        continue;
                    }
                    c.skiprz(nin, -2);
                    c.findf(matd, 3, mtd, nin);
                    lbl = 105;
                }
                _ => {
                    // 300: nheat is finished.
                    c.skiprz(nin, -4);
                    let _ = (hmin, hmax, jrec);
                    return;
                }
            }
        }
    }
}

/// `indx`: indices in mtp of the requested kermas that reaction mt
/// contributes to.
fn indx(npk: i32, mtp: &[i32; 29], mt19: i32, mt: i32, imt: &mut [i32; 31], nmt: &mut i32, iverf: i32) {
    *nmt = 0;
    if npk < 3 {
        return;
    }
    for i in 3..=npk {
        let mtpi = mtp[i as usize];
        let mut iflag = 0;
        if mtpi == 442 {
            iflag = 1;
        }
        if mtpi == 443 {
            iflag = 1;
        }
        if mtpi == 444 {
            iflag = 1;
        }
        if mtpi == 445 && mt == 2 {
            iflag = 1;
        }
        if mtpi == 446 && (51..=91).contains(&mt) {
            iflag = 1;
        }
        if mtpi == 447 && (102..=120).contains(&mt) {
            iflag = 1;
        }
        if iverf < 6 {
            if mtpi == 447 && (700..=799).contains(&mt) {
                iflag = 1;
            }
        } else if mtpi == 447 && (600..=849).contains(&mt) {
            iflag = 1;
        }
        if mtpi == 303 && mt != 2 && mt != 3 {
            iflag = 1;
        }
        if mtpi == 304 && (51..=91).contains(&mt) {
            iflag = 1;
        }
        if mtpi == 318 {
            if mt19 == 0 && mt == 18 {
                iflag = 1;
            }
            if mt19 == 1 && (mt == 19 || mt == 20 || mt == 21 || mt == 38) {
                iflag = 1;
            }
            if mt19 == 2 && mt == 18 {
                iflag = 1;
            }
        }
        if mtpi == 401 {
            if (102..=120).contains(&mt) {
                iflag = 1;
            }
            if iverf < 6 && (700..800).contains(&mt) {
                iflag = 1;
            }
            if iverf >= 6 && (600..850).contains(&mt) {
                iflag = 1;
            }
        }
        if iverf < 6 {
            if mtpi == 403 && (700..=719).contains(&mt) {
                iflag = 1;
            }
            if mtpi == 404 && (720..=739).contains(&mt) {
                iflag = 1;
            }
            if mtpi == 405 && (740..=759).contains(&mt) {
                iflag = 1;
            }
            if mtpi == 406 && (760..=779).contains(&mt) {
                iflag = 1;
            }
            if mtpi == 407 && (780..=799).contains(&mt) {
                iflag = 1;
            }
        } else {
            if mtpi == 403 && (600..=649).contains(&mt) {
                iflag = 1;
            }
            if mtpi == 404 && (650..=699).contains(&mt) {
                iflag = 1;
            }
            if mtpi == 405 && (700..=749).contains(&mt) {
                iflag = 1;
            }
            if mtpi == 406 && (750..=799).contains(&mt) {
                iflag = 1;
            }
            if mtpi == 407 && (800..=849).contains(&mt) {
                iflag = 1;
            }
        }
        if mtpi == mt + 300 {
            iflag = 1;
        }
        if iflag == 1 {
            *nmt += 1;
            imt[*nmt as usize] = i;
        }
    }
}

impl Heatr {
    // =============================================================== capdam
    /// Damage energy for capture reactions (call with ee=0 to initialize).
    fn capdam(&mut self, c: &Ctx, ee: f64, dame: &mut f64, q: f64, za: f64, awr: f64, mtd: i32) {
        let third = 0.333333333e0f64;
        let econ = 1.029e6f64;
        let eps = 1.0e-10f64;
        let step = 1.1e0f64;
        let zero = 0.0f64;
        let iverf = c.iverf;

        //--initialize
        let mut e = ee;
        if e == zero {
            let s = &mut self.cd;
            s.zx = 0.0;
            if mtd > 102 {
                s.zx = 1.0;
            }
            if mtd == 106 || mtd == 107 {
                s.zx = 2.0;
            }
            if iverf <= 5 && (mtd == 778 || mtd == 798) {
                s.zx = 2.0;
            }
            s.ax = 0.0;
            if mtd > 102 {
                s.ax = 1.0;
            }
            if mtd == 104 {
                s.ax = 2.0;
            }
            if iverf <= 5 && mtd == 738 {
                s.ax = 2.0;
            }
            if mtd == 105 || mtd == 106 {
                s.ax = 3.0;
            }
            if iverf <= 5 && (mtd == 758 || mtd == 778) {
                s.ax = 3.0;
            }
            if mtd == 107 {
                s.ax = 4.0;
            }
            if iverf <= 5 && mtd == 798 {
                s.ax = 4.0;
            }
            s.denom = 1.0 / pow(awr, third);
            if s.ax != zero {
                s.denom = 1.0 / (pow(s.ax, third) + pow(awr, third));
            }
            let iz = nint(za / 1000.0);
            s.z = iz as f64;
            s.en = 0.0;
            let (en, z, zx, ax) = (s.en, s.z, s.zx, s.ax);
            self.cd.damn = self.df(en, z - zx, awr + 1.0 - ax, z, awr);
            *dame = 0.0;
            self.cd.aw1fac = 1.0 / (awr + 1.0);
            return;
        }

        //--normal entry
        //--if energy not in current range,
        //--construct new interpolation range.
        if ee >= self.cd.en * (1.0 - eps) {
            self.cd.el = self.cd.en;
            self.cd.daml = self.cd.damn;
            e = step * self.cd.el;
            if e < ee * (1.0 - eps) {
                e = (1.0 - eps) * ee;
            }
            *dame = 0.0;
            let (z, zx, ax, denom, aw1fac) = (self.cd.z, self.cd.zx, self.cd.ax, self.cd.denom, self.cd.aw1fac);
            if mtd == 102 {
                // capture followed by gamma emission
                let mut er = e * aw1fac;
                let ea = awr * er + q;
                er += ea * ea * self.rtm / 2.0;
                *dame = self.df(er, z, awr + 1.0, z, awr);
            } else if ax != zero {
                // capture followed by particle emission
                let ec = econ * zx * z * denom;
                let mut ea = q + awr * e * aw1fac;
                if ea >= zero {
                    let et = (awr + 1.0 - ax) * e * aw1fac;
                    if ea > ec * (1.0 + eps) {
                        ea = ec;
                    }
                    for iq in 0..4 {
                        let er = (et - 2.0 * sqrt(et * ax * ea) * QP4[iq] + ax * ea) * aw1fac;
                        *dame += QW4[iq] * self.df(er, z - zx, awr + 1.0 - ax, z, awr) / 2.0;
                    }
                }
            }
            self.cd.en = e;
            self.cd.damn = *dame;
        }

        //--interpolate for damage energy
        let s = &self.cd;
        c.terp1(s.el, s.daml, s.en, s.damn, ee, dame, 2);
    }

    // =============================================================== disbar
    /// Ebar and damage energy for discrete scattering.
    fn disbar(
        &mut self,
        c: &mut Ctx,
        ee: f64,
        ebar: &mut f64,
        dame: &mut f64,
        q: f64,
        za: f64,
        awr: f64,
        ntape: i32,
        matd: i32,
        mtd: i32,
        a: &mut RV,
        nwa: i32,
    ) {
        let edis = 25.0e0f64;
        let step = 1.1e0f64;
        let small = 1.0e-10f64;
        let zero = 0.0f64;
        let mut fl = [0.0f64; 66];
        let mut p = [0.0f64; 66];
        let mut idis = 0i32;
        let iverf = c.iverf;

        //--initialize if e=0.
        let mut e = ee;
        if e == zero {
            self.db.imiss = 0;
            for i in 1..=self.nmiss4 as usize {
                if self.miss4[i] == mtd {
                    self.db.imiss = 1;
                }
            }
            if self.db.imiss != 1 {
                let mfd = 4;
                let mut nld = 60;
                let lcd = 2;
                let mut enext = self.db.enext;
                self.hgtfle(c, e, &mut enext, &mut idis, &mut fl, &mut nld, lcd, matd, mfd, mtd, ntape, a, nwa);
                self.db.enext = enext;
            } else {
                self.db.enext = self.etop;
            }
            let s = &mut self.db;
            s.thresh = ((awr + 1.0) / awr) * (-q);
            s.awp = 1.0;
            if iverf <= 5 && (700..718).contains(&mtd) {
                s.awp = PNRATIO;
            }
            if iverf <= 5 && (720..738).contains(&mtd) {
                s.awp = DNRATIO;
            }
            if iverf <= 5 && (740..758).contains(&mtd) {
                s.awp = TNRATIO;
            }
            if iverf <= 5 && (760..778).contains(&mtd) {
                s.awp = HNRATIO;
            }
            if iverf <= 5 && (780..798).contains(&mtd) {
                s.awp = ANRATIO;
            }
            if iverf >= 6 && (600..649).contains(&mtd) {
                s.awp = PNRATIO;
            }
            if iverf >= 6 && (650..699).contains(&mtd) {
                s.awp = DNRATIO;
            }
            if iverf >= 6 && (700..749).contains(&mtd) {
                s.awp = TNRATIO;
            }
            if iverf >= 6 && (750..799).contains(&mtd) {
                s.awp = HNRATIO;
            }
            if iverf >= 6 && (800..849).contains(&mtd) {
                s.awp = ANRATIO;
            }
            s.afact = s.awp / powi(awr + 1.0, 2);
            s.arat = s.awp / (awr + 1.0 - s.awp);
            s.en = 0.0;
            s.cn = 0.0;
            s.el = 0.0;
            s.cl = 0.0;
            s.enx = 0.0;
            if self.idame != 0 {
                let iz = nint(za / 1000.0);
                self.db.z = iz as f64;
                let z = self.db.z;
                self.db.damn = self.df(e, z, awr, z, awr);
                self.db.daml = self.db.damn;
                if mtd == 2 {
                    self.db.daml = 0.0;
                    self.db.enx = edis * self.db.arat / (4.0 * self.db.afact);
                    self.db.damn = 0.0;
                }
            }
            return;
        }

        //--normal entry.
        let mut go120 = false;
        let mut go130 = ee <= self.db.en * (1.0 + small);
        if !go130 {
            self.db.el = self.db.en;
            self.db.cl = self.db.cn;
            if self.idame != 0 {
                self.db.daml = self.db.damn;
            }
            e = step * self.db.el;
            if self.db.enext < e * (1.0 - small) {
                e = self.db.enext;
            }
            if ee > e * (1.0 + small) {
                e = ee;
            }
            let mut r = 1.0f64;
            if mtd != 2 {
                if self.db.thresh >= e * (1.0 - small) {
                    r = 0.0;
                } else {
                    r = sqrt(1.0 - self.db.thresh / e);
                }
            }
            let nld;
            if self.db.imiss != 1 {
                let mfd = 4;
                let mut n = 60;
                let lcd = 2;
                let mut enext = self.db.enext;
                self.hgtfle(c, e, &mut enext, &mut idis, &mut fl, &mut n, lcd, matd, mfd, mtd, ntape, a, nwa);
                self.db.enext = enext;
                nld = n;
            } else {
                nld = 1;
                fl[1] = 1.0;
                fl[2] = 0.0;
            }
            self.db.en = e;
            let b = r * sqrt(awr / self.db.arat);
            let g = b * self.db.arat;
            self.db.cn = 0.0;
            let mut test = 1.0f64;
            test /= 1000000.0;
            if (self.db.awp - 1.0).abs() > test {
                go120 = true;
            } else {
                //--compute ebar
                let mut wbar = fl[2];
                if wbar > QP64[63] {
                    wbar = QP64[63];
                }
                self.db.cn = (1.0 + 2.0 * b * wbar + b * b) * self.db.afact;
                if self.idame == 0 {
                    go130 = true;
                } else if e < self.db.enx * (1.0 - small) {
                    go130 = true;
                } else {
                    go120 = true;
                }
            }
            if go120 {
                //--angle integration by gauss-legendre method.
                *dame = 0.0;
                let z = self.db.z;
                let (afact, arat) = (self.db.afact, self.db.arat);
                for iq in 0..64 {
                    let u = QP64[iq];
                    legndr(u, &mut p[1..], nld);
                    let mut f = 0.0f64;
                    for il in 1..=nld as usize {
                        f += (2 * il - 1) as f64 * fl[il] * p[il] / 2.0;
                    }
                    let e2 = e * (1.0 - 2.0 * g * u + g * g) * afact / arat;
                    *dame += QW64[iq] * f * self.df(e2, z, awr, z, awr);
                }
                if *dame < zero {
                    *dame = zero;
                }
                self.db.damn = *dame;
            }
        }

        //--interpolate for results.
        // 130
        let mut ce = 0.0f64;
        let s = &self.db;
        c.terp1(s.el, s.cl, s.en, s.cn, ee, &mut ce, 2);
        *ebar = ee * ce;
        if self.idame > 0 {
            c.terp1(s.el, s.daml, s.en, s.damn, ee, dame, 2);
        }
    }

    // =================================================================== df
    /// Lindhard damage function (call with e=0 for each reaction to set the
    /// constants).
    fn df(&mut self, e: f64, zr: f64, ar: f64, zl: f64, al: f64) -> f64 {
        let twothd = 0.666666667e0f64;
        let threeq = 0.75e0f64;
        let sixth = 0.166666667e0f64;
        let onep5 = 1.5e0f64;
        let c1 = 30.724e0f64;
        let c2 = 0.0793e0f64;
        let c3 = 3.4008e0f64;
        let c4 = 0.40244e0f64;
        let zero = 0.0f64;
        if zr == zero {
            0.0
        } else if e <= zero {
            let el = c1 * zr * zl * sqrt(pow(zr, twothd) + pow(zl, twothd)) * (ar + al) / al;
            self.df_rel = 1.0 / el;
            let denom = pow(pow(zr, twothd) + pow(zl, twothd), threeq) * pow(ar, onep5) * sqrt(al);
            self.df_fl = c2 * pow(zr, twothd) * sqrt(zl) * pow(ar + al, onep5) / denom;
            0.0
        } else if e < self.brk {
            0.0
        } else {
            let ep = e * self.df_rel;
            e / (1.0 + self.df_fl * (c3 * pow(ep, sixth) + c4 * pow(ep, threeq) + ep))
        }
    }

    // =============================================================== hgtfle
    /// Retrieve or compute Legendre coefficients at e (isotropic below the
    /// first point). fle is 1-based (fle[i] is fle(i)).
    fn hgtfle(
        &mut self,
        c: &mut Ctx,
        e: f64,
        enext: &mut f64,
        idis: &mut i32,
        fle: &mut [f64; 66],
        nle: &mut i32,
        lcd: i32,
        matd: i32,
        mfd: i32,
        mtd: i32,
        nin: i32,
        cc: &mut RV,
        _nc: i32,
    ) {
        let small = 1.0e-10f64;
        let zero = 0.0f64;
        let mut nb = 0i32;
        let mut nwc = 0i32;

        //--initialize
        *idis = 0;
        if e <= zero {
            let s = &mut self.hf;
            s.iso = 0;
            c.findf(matd, mfd, mtd, nin);
            c.contio(nin, 0, 0, &mut cc.all(), &mut nb, &mut nwc);
            s.awr = cc[2];
            let lvt = nint(cc[3]);
            s.ltt = nint(cc[4]);
            s.ltt3 = s.ltt;
            if s.ltt == 3 {
                s.ltt = 1;
                s.lttn = 1;
            }
            if lvt == 0 {
                c.contio(nin, 0, 0, &mut cc.all(), &mut nb, &mut nwc);
            } else {
                c.listio(nin, 0, 0, &mut cc.all(), &mut nb, &mut nwc);
                while nb != 0 {
                    c.moreio(nin, 0, 0, &mut cc.at(7), &mut nb, &mut nwc);
                }
            }
            s.iso = nint(cc[3]);
            s.lct = nint(cc[4]);
            // check for lab distribution for two-body scattering.
            let test = 10.0;
            if s.lct == 1 && s.awr >= test {
                if (51..=90).contains(&mtd) {
                    let strng = format!("lab distribution changed to cm for mt {}", ffmt::i(mtd, 3));
                    c.mess("hgtfle", &strng, " ");
                    s.lct = 2;
                }
            }
            if s.iso != 1 {
                c.tab2io(nin, 0, 0, &mut cc.at(1), &mut nb, &mut nwc);
                s.ne = nint(cc[6]);
                // read in raw data at first incident energy point
                // retrieve or compute legendre coefficients.
                s.iraw = 1 + nwc;
                let iraw = s.iraw;
                if s.ltt == 1 {
                    c.listio(nin, 0, 0, &mut cc.at(iraw), &mut nb, &mut nwc);
                }
                if s.ltt == 2 {
                    c.tab1io(nin, 0, 0, &mut cc.at(iraw), &mut nb, &mut nwc);
                }
                if nb != 0 {
                    let mut l = iraw;
                    while nb != 0 {
                        l += nwc;
                        c.moreio(nin, 0, 0, &mut cc.at(l), &mut nb, &mut nwc);
                    }
                }
                s.elo = cc[iraw + 1];
                s.nlo = *nle;
                let (lct, ltt, awr) = (s.lct, s.ltt, s.awr);
                hgetco(c, &mut s.flo, &mut s.nlo, lcd, &cc.at(iraw), lct, ltt, awr, idis);
                if s.ltt == 1 {
                    c.listio(nin, 0, 0, &mut cc.at(iraw), &mut nb, &mut nwc);
                }
                if s.ltt == 2 {
                    c.tab1io(nin, 0, 0, &mut cc.at(iraw), &mut nb, &mut nwc);
                }
                if nb != 0 {
                    let mut l = iraw;
                    while nb != 0 {
                        l += nwc;
                        c.moreio(nin, 0, 0, &mut cc.at(l), &mut nb, &mut nwc);
                    }
                }
                s.ehi = cc[iraw + 1];
                s.nhi = *nle;
                hgetco(c, &mut s.fhi, &mut s.nhi, lcd, &cc.at(iraw), lct, ltt, awr, idis);
                // check for total isotropy.
                if s.ne == 2 && s.nlo == 1 && s.nhi == 1 {
                    s.iso = 1;
                }
            }
            if s.iso == 1 {
                *nle = 1;
                *enext = self.etop;
            } else {
                *enext = s.elo;
            }
            let s = &mut self.hf;
            s.nne = 2;
            s.ir = 1;
            s.inn = nint(cc[6 + 2 * s.ir]);
            return;
        }

        //--normal entry
        //--is desired energy in current panel
        let s = &mut self.hf;
        let mut lbl = if *nle == 1 || s.iso == 1 {
            150
        } else if e < s.ehi * (1.0 - small) {
            130
        } else if s.nne == s.ne && e < s.ehi + s.ehi / 100.0 {
            130
        } else {
            120
        };
        loop {
            match lbl {
                120 => {
                    //--no.  slide raw high energy data into low energy positions.
                    for i in 1..=s.nhi as usize {
                        s.flo[i] = s.fhi[i];
                    }
                    s.nlo = s.nhi;
                    s.elo = s.ehi;

                    //--read in new raw data at high energy.
                    if s.nne == s.ne && s.ltt3 == 3 && s.lttn == 1 {
                        c.tab2io(nin, 0, 0, &mut cc.at(1), &mut nb, &mut nwc);
                        s.ne = nint(cc[6]);
                        s.nne = 0;
                        s.ir = 1;
                        s.ltt = 2;
                        s.lttn = 2;
                    } else if s.nne == s.ne {
                        c.error("hgtfle", "desired energy above highest given.", " ");
                    }
                    let iraw = s.iraw;
                    if s.ltt == 1 {
                        c.listio(nin, 0, 0, &mut cc.at(iraw), &mut nb, &mut nwc);
                    }
                    if s.ltt == 2 {
                        c.tab1io(nin, 0, 0, &mut cc.at(iraw), &mut nb, &mut nwc);
                    }
                    if nb != 0 {
                        let mut l = iraw;
                        while nb != 0 {
                            l += nwc;
                            c.moreio(nin, 0, 0, &mut cc.at(l), &mut nb, &mut nwc);
                        }
                    }
                    s.ehi = cc[iraw + 1];
                    s.nhi = *nle;
                    let (lct, ltt, awr) = (s.lct, s.ltt, s.awr);
                    hgetco(c, &mut s.fhi, &mut s.nhi, lcd, &cc.at(iraw), lct, ltt, awr, idis);
                    s.nne += 1;
                    if s.nne > nint(cc[5 + 2 * s.ir]) {
                        s.ir += 1;
                    }
                    s.inn = nint(cc[6 + 2 * s.ir]);
                    if s.ehi <= e * (1.0 + small) && s.nne < s.ne {
                        continue;
                    }
                    lbl = 130;
                }
                130 => {
                    //--yes.  interpolate for coefficients at desired energy.
                    if e < s.elo * (1.0 - small) {
                        lbl = 140;
                        continue;
                    }
                    let mut nlmax = s.nlo;
                    if nlmax < s.nhi {
                        nlmax = s.nhi;
                    }
                    for i in 1..=*nle as usize {
                        if i as i32 <= nlmax {
                            let mut innt = s.inn;
                            if innt == 4 || innt == 5 {
                                if s.flo[i] * s.fhi[i] <= zero {
                                    innt -= 2;
                                }
                            }
                            c.terp1(s.elo, s.flo[i], s.ehi, s.fhi[i], e, &mut fle[i], innt);
                        } else {
                            fle[i] = 0.0;
                        }
                    }
                    *nle = nlmax;
                    *enext = s.ehi;
                    if s.inn == 1 {
                        *idis = 1;
                    }
                    return;
                }
                140 => {
                    //--return isotropic distribution below first point.
                    fle[1] = 1.0;
                    for il in 2..=*nle as usize {
                        fle[il] = 0.0;
                    }
                    *nle = 1;
                    *enext = s.ehi;
                    *idis = 1;
                    return;
                }
                _ => {
                    //--isotropic distribution or request.
                    fle[1] = 1.0;
                    if *nle > 1 {
                        for il in 2..=*nle as usize {
                            fle[il] = 0.0;
                        }
                        *nle = 1;
                    }
                    *enext = self.etop;
                    *idis = 0;
                    return;
                }
            }
        }
    }
}

/// `hgetco`: Legendre coefficients from raw ENDF data (fl is 1-based).
fn hgetco(c: &mut Ctx, fl: &mut [f64; 66], nl: &mut i32, lcd: i32, cv: &F, lct: i32, ltt: i32, awr: f64, idis: &mut i32) {
    let toler = 1.0e-6f64;
    let rtoler = 1.0e6f64;
    let half = 0.5e0f64;
    let nlmax = 65;
    let zero = 0.0f64;
    let mut p = [0.0f64; 66];

    //--determine representation
    let l = *nl - 1;
    if *nl > nlmax {
        c.error("hgetco", "limited to 64 legendre coefficients.", " ");
    }
    let integrate = ltt == 2 || (lcd == 1 && lct != 1) || (lcd == 2 && lct < 2);
    if !integrate {
        //--retrieve coefficients from raw data
        let np = nint(cv[5]) + 1;
        fl[1] = 1.0;
        for il in 2..=*nl {
            if il <= np {
                fl[il as usize] = cv[il + 5];
            } else {
                fl[il as usize] = 0.0;
            }
        }
        if np < *nl {
            *nl = np;
        }
        return;
    }

    //--integrate for fl in desired system
    fl[1] = 1.0;
    for il in 2..=*nl as usize {
        fl[il] = 0.0;
    }
    let mut ipc = 2;
    let mut irc = 1;
    let mut xnext = 0.0f64;
    for iq in 0..64 {
        let x = QP64[iq];

        //--compute scattering probability
        let mut y;
        if ltt != 2 {
            // from coefficients
            let np = nint(cv[5]);
            legndr(x, &mut p[1..], np);
            y = half;
            for ip in 1..=np {
                y = (2 * ip + 1) as f64 * p[(ip + 1) as usize] * cv[ip + 6] / 2.0 + y;
            }
        } else {
            // from tabular data
            y = 0.0;
            c.terpa(&mut y, x, &mut xnext, idis, cv, &mut ipc, &mut irc);
        }

        //--multiply scattering probability by legendre polynomials
        //--in desired system and accumulate
        let mut xn = x;
        if lcd == 2 && lct == 1 {
            c.error("hgetco", "lab to cm conversion not coded.", " ");
        }
        if lcd == 1 && lct >= 2 {
            xn = (1.0 + awr * x) / sqrt(1.0 + awr * awr + 2.0 * awr * x);
        }
        legndr(xn, &mut p[1..], l);
        y *= QW64[iq];
        for il in 2..=*nl as usize {
            fl[il] += y * p[il];
        }
    }

    //--reduce number of significant figures
    //--and count number of coefficients required
    let mut nlz = 1;
    for il in 2..=*nl {
        let j = nint(fl[il as usize] * rtoler);
        fl[il as usize] = j as f64 * toler;
        if fl[il as usize].abs() > zero {
            nlz = il;
        }
    }
    *nl = nlz;
}

impl Heatr {
    // =============================================================== conbar
    /// Neutron ebar and damage energy for continuum scattering (File 5).
    /// nbmax is changed by hgtyld, as the Fortran argument is.
    fn conbar(
        &mut self,
        c: &mut Ctx,
        e: f64,
        ebar: &mut f64,
        dame: &mut f64,
        nin: i32,
        nscr: i32,
        matd: i32,
        mtd: i32,
        cc: &mut RV,
        ncmax: i32,
        b: &mut RV,
        mut nbmax: i32,
        yld: &mut f64,
        q: f64,
    ) {
        const NKMAX: i32 = 12;
        let emax = 1.0e10f64;
        let step = 1.5e0f64;
        let small = 1.0e-10f64;
        let zero = 0.0f64;
        let (mut nb, mut nw) = (0i32, 0i32);
        let mut enext = 0.0f64;
        let mut idis = 0i32;

        //--initialize if e=0.
        if e == zero {
            let mut l = 1;
            if !is_fis(mtd) {
                *yld = 1.0;
                if mtd == 16 || mtd == 30 {
                    *yld = 2.0;
                }
                if mtd == 24 || mtd == 26 {
                    *yld = 2.0;
                }
                if mtd == 11 || mtd == 41 {
                    *yld = 2.0;
                }
                if mtd == 17 || mtd == 25 {
                    *yld = 3.0;
                }
                if mtd == 38 || mtd == 42 {
                    *yld = 3.0;
                }
                if mtd == 37 {
                    *yld = 4.0;
                }
            } else {
                self.cb.mf1 = 1;
                self.cb.mt1 = 452;
                if self.nply > 0 {
                    self.cb.mt1 = 456;
                }
                // delayed neutrons are treated the same as fast neutrons.
                // but when nply.ne.0, use prompt neutron data only
                let (mf1, mt1) = (self.cb.mf1, self.cb.mt1);
                self.hgtyld(c, e, &mut enext, &mut idis, yld, matd, mf1, mt1, nscr, b, &mut nbmax);
            }
            let mfd = 5;
            c.findf(matd, mfd, mtd, nin);
            c.contio(nin, 0, 0, &mut cc.at(l), &mut nb, &mut nw);
            self.cb.nktot = nint(cc[5]);
            if self.cb.nktot > NKMAX {
                c.error("conbar", "nktot gt nkmax.", " ");
            }
            self.cb.awr = cc[2];
            self.cb.awfac = 1.0 / self.cb.awr;
            self.cb.aw1fac = 1.0 / (self.cb.awr + 1.0);
            let iz = nint(cc[1] / 1000.0);
            self.cb.z = iz as f64;
            // treat damage for (n,np) as (n,p),etc.
            self.cb.za = cc[1];
            self.cb.mtt = 0;
            if mtd == 22 {
                self.cb.mtt = 107;
            }
            if mtd == 28 {
                self.cb.mtt = 103;
            }
            if mtd == 32 {
                self.cb.mtt = 104;
            }
            if mtd == 33 {
                self.cb.mtt = 105;
            }
            if mtd == 34 {
                self.cb.mtt = 106;
            }
            if self.cb.mtt > 0 {
                let (za, awr, mtt) = (self.cb.za, self.cb.awr, self.cb.mtt);
                self.capdam(c, e, dame, q, za, awr, mtt);
            }
            let mut il = 1;
            let mut ik = 0;
            while ik < self.cb.nktot {
                ik += 1;
                let iku = ik as usize;
                self.cb.d2[iku] = 0.0;
                self.cb.e2[iku] = 0.0;
                self.cb.loc[iku] = il;
                l = il;
                c.tab1io(nin, 0, 0, &mut cc.at(l), &mut nb, &mut nw);
                l += nw;
                while nb != 0 {
                    c.moreio(nin, 0, 0, &mut cc.at(l), &mut nb, &mut nw);
                    l += nw;
                }
                let lf = nint(cc[il + 3]);
                il = l;

                //--analytic subsection.
                if lf > 1 {
                    if lf == 3 {
                        let (awr, aw1fac) = (self.cb.awr, self.cb.aw1fac);
                        self.cb.t1 = (awr * awr + 1.0) * aw1fac * aw1fac;
                        self.cb.t2 = awr * aw1fac;
                    } else {
                        c.tab1io(nin, 0, 0, &mut cc.at(l), &mut nb, &mut nw);
                        l += nw;
                        while nb != 0 {
                            if l > ncmax {
                                c.error("conbar", "insufficient storage for raw endf data.", " ");
                            }
                            c.moreio(nin, 0, 0, &mut cc.at(l), &mut nb, &mut nw);
                            l += nw;
                        }
                        il = l;
                        if lf == 5 || lf == 11 {
                            c.tab1io(nin, 0, 0, &mut cc.at(l), &mut nb, &mut nw);
                            l += nw;
                            while nb != 0 {
                                if l > ncmax {
                                    c.error("conbar", "insufficient storage for raw endf data.", " ");
                                }
                                c.moreio(nin, 0, 0, &mut cc.at(l), &mut nb, &mut nw);
                                l += nw;
                            }
                            il = l;
                        }
                    }

                //--tabulated subsection
                } else {
                    c.tab2io(nin, 0, 0, &mut cc.at(l), &mut nb, &mut nw);
                    let ne = nint(cc[l + 5]);
                    l += nw;
                    let mut nne = 0;
                    let lnow = l;
                    let iraw = lnow + 3 * ne;

                    //--read spectrum for each incident energy
                    //--compute heating and damage contributions
                    while nne < ne {
                        l = iraw;
                        c.tab1io(nin, 0, 0, &mut cc.at(l), &mut nb, &mut nw);
                        l += nw;
                        while nb != 0 {
                            if l > ncmax {
                                c.error("conbar", "insufficient storage for raw endf data.", " ");
                            }
                            c.moreio(nin, 0, 0, &mut cc.at(l), &mut nb, &mut nw);
                            l += nw;
                        }
                        nne += 1;
                        let ehi = cc[iraw + 1];
                        let mut fhi = 0.0f64;
                        tabbar(c, &mut fhi, &mut cc.at(iraw), lf);
                        let mut dhi = 0.0f64;
                        if self.idame != 0 {
                            if self.cb.mtt != 0 {
                                let (za, awr, mtt) = (self.cb.za, self.cb.awr, self.cb.mtt);
                                self.capdam(c, ehi, &mut dhi, q, za, awr, mtt);
                            } else {
                                let (z, awr, awfac) = (self.cb.z, self.cb.awr, self.cb.awfac);
                                self.tabdam(ehi, &mut dhi, z, awr, awfac, &cc.at(iraw));
                            }
                        }
                        cc[lnow + 3 * (nne - 1)] = ehi;
                        cc[lnow + 3 * (nne - 1) + 1] = fhi;
                        cc[lnow + 3 * (nne - 1) + 2] = dhi;
                    }
                    il = lnow + 3 * ne;
                }
            }
            return;
        }

        //--normal entry.
        enext = emax;
        idis = 0;
        *ebar = 0.0;
        *dame = 0.0;

        //--get energy dependent yld.
        if mtd == 18 || mtd == 19 || mtd == 20 || mtd == 21 || mtd == 38 {
            self.cb.mf1 = 1;
            let (mf1, mt1) = (self.cb.mf1, self.cb.mt1);
            self.hgtyld(c, e, &mut enext, &mut idis, yld, matd, mf1, mt1, nscr, b, &mut nbmax);
        }

        //--loop over subsections.
        for ik in 1..=self.cb.nktot {
            let iku = ik as usize;
            let mut lnow = self.cb.loc[iku];
            let lf = nint(cc[lnow + 3]);

            //--interpolate for fractional probability.
            let mut ip = 2;
            let mut ir = 1;
            let mut pe = 0.0f64;
            let mut eihi = 0.0f64;
            let mut idisc = 0i32;
            c.terpa(&mut pe, e, &mut eihi, &mut idisc, &cc.at(lnow), &mut ip, &mut ir);
            if (eihi - enext).abs() < small * enext && idisc > idis {
                idis = idisc;
            }
            if eihi < enext * (1.0 - small) {
                idis = idisc;
            }
            if eihi < enext * (1.0 - small) {
                enext = eihi;
            }
            let nr = nint(cc[lnow + 4]);
            let np = nint(cc[lnow + 5]);
            let ll = lnow + 6 + 2 * nr + 2 * (np - 1);
            let etp = cc[ll];
            if pe > zero {
                //--tabulated subsection.
                //--interpolate for results
                if lf == 1 {
                    lnow = lnow + 6 + 2 * nr + 2 * np;
                    let ne = nint(cc[lnow + 5]);
                    let nnr = nint(cc[lnow + 4]);
                    let inn = 2;
                    lnow = lnow + 6 + 2 * nnr;
                    let mut ie = 0;
                    let mut idone = 0;
                    let mut ii = 0;
                    while ie < ne && idone == 0 {
                        ie += 1;
                        if cc[lnow + 3 * ie] >= e * (1.0 - small) {
                            ii = lnow + 3 * (ie - 1);
                            idone = 1;
                        }
                    }
                    if idone == 0 {
                        ii = lnow + 3 * (ne - 1);
                    }
                    let mut f = 0.0f64;
                    c.terp1(cc[ii], cc[ii + 1], cc[ii + 3], cc[ii + 4], e, &mut f, inn);
                    *ebar += f * pe;
                    if self.idame != 0 {
                        let mut d = 0.0f64;
                        c.terp1(cc[ii], cc[ii + 2], cc[ii + 3], cc[ii + 5], e, &mut d, inn);
                        *dame += d * pe;
                    }

                //--analytic subsection.
                } else {
                    let mut s = 0.0f64;
                    let (t1, t2) = (self.cb.t1, self.cb.t2);
                    anabar(c, &mut s, e, &mut cc.at(lnow), t1, t2);
                    *ebar += s * pe;
                    if self.idame != 0 {
                        if self.cb.mtt != 0 {
                            let (za, awr, mtt) = (self.cb.za, self.cb.awr, self.cb.mtt);
                            self.capdam(c, e, dame, q, za, awr, mtt);
                        } else {
                            if e >= self.cb.e2[iku] * (1.0 - small) {
                                let mut idone = 0;
                                while idone == 0 {
                                    self.cb.d1[iku] = self.cb.d2[iku];
                                    self.cb.e1[iku] = self.cb.e2[iku];
                                    self.cb.e2[iku] = step * self.cb.e1[iku];
                                    if self.cb.e2[iku] > etp * (1.0 + small) {
                                        self.cb.e2[iku] = etp;
                                    }
                                    if self.cb.e1[iku] == zero {
                                        self.cb.e2[iku] = e;
                                    }
                                    let (e2, z, awr) = (self.cb.e2[iku], self.cb.z, self.cb.awr);
                                    let mut d2 = 0.0f64;
                                    self.anadam(c, e2, &mut d2, z, awr, &mut cc.at(lnow));
                                    self.cb.d2[iku] = d2;
                                    if (self.cb.e2[iku] - etp).abs() < small * etp {
                                        idone = 1;
                                    } else if self.cb.e1[iku] != zero && e <= self.cb.e2[iku] * (1.0 + small) {
                                        idone = 1;
                                    }
                                }
                            }
                            let mut d = 0.0f64;
                            let cb = &self.cb;
                            c.terp1(cb.e1[iku], cb.d1[iku], cb.e2[iku], cb.d2[iku], e, &mut d, 2);
                            *dame += d * pe;
                        }
                    }
                }
            }
        }
    }

    // =============================================================== hgtyld
    /// Retrieve or compute the yield for mat, mf, mt on itape (e=0 to
    /// initialize; na is changed as the Fortran argument is).
    fn hgtyld(
        &mut self,
        c: &mut Ctx,
        e: f64,
        enext: &mut f64,
        idis: &mut i32,
        yld: &mut f64,
        mat: i32,
        mf: i32,
        mt: i32,
        itape: i32,
        a: &mut RV,
        na: &mut i32,
    ) {
        let emax = 1.0e10f64;
        let zero = 0.0f64;
        let (mut nb, mut nw) = (0i32, 0i32);

        //--initialize.
        if e == zero {
            c.findf(mat, mf, mt, itape);
            c.contio(itape, 0, 0, &mut a.all(), &mut nb, &mut nw);
            self.hy.lnu = nint(a[4]);
            if c.mth == 455 {
                let mut loc = 1;
                c.listio(itape, 0, 0, &mut a.at(loc), &mut nb, &mut nw);
                loc += nw;
                while nb != 0 {
                    c.moreio(itape, 0, 0, &mut a.at(loc), &mut nb, &mut nw);
                    loc += nw;
                }
                let lnd = nint(a[5]);
                if lnd != 6 && lnd != 8 {
                    c.error("hgtyld", "illegal lnd, must be 6 or 8", " ");
                }
                loc = 1;
                c.tab1io(itape, 0, 0, &mut a.at(loc), &mut nb, &mut nw);
                loc += nw;
                if loc > *na {
                    c.error("hgtyld", "storage exceeded.", " ");
                }
                while nb != 0 {
                    c.moreio(itape, 0, 0, &mut a.at(loc), &mut nb, &mut nw);
                    loc += nw;
                    if loc > *na {
                        c.error("hgtyld", "storage exceeded.", " ");
                    }
                }
                *na = loc - 1;
                self.hy.ip = 2;
                self.hy.ir = 1;
                let nr = nint(a[5]);
                *enext = a[7 + 2 * nr];
            } else {
                let mut loc = 1;
                if self.hy.lnu != 1 {
                    c.tab1io(itape, 0, 0, &mut a.at(loc), &mut nb, &mut nw);
                    loc += nw;
                    if loc > *na {
                        c.error("hgtyld", "storage exceeded.", " ");
                    }
                    while nb != 0 {
                        c.moreio(itape, 0, 0, &mut a.at(loc), &mut nb, &mut nw);
                        loc += nw;
                        if loc > *na {
                            c.error("hgtyld", "storage exceeded.", " ");
                        }
                    }
                    *na = loc - 1;
                    self.hy.ip = 2;
                    self.hy.ir = 1;
                    let nr = nint(a[5]);
                    *enext = a[7 + 2 * nr];
                } else {
                    loc = 1;
                    c.listio(itape, 0, 0, &mut a.at(loc), &mut nb, &mut nw);
                    loc += nw;
                    *na = nw;
                    while nb != 0 {
                        c.moreio(itape, 0, 0, &mut a.at(loc), &mut nb, &mut nw);
                        loc += nw;
                        *na += nw;
                    }
                    *enext = emax;
                }
            }
            *idis = 0;
            return;
        }

        //--tabulated data.
        if self.hy.lnu != 1 {
            let (mut ip, mut ir) = (self.hy.ip, self.hy.ir);
            c.terpa(yld, e, enext, idis, &a.all(), &mut ip, &mut ir);
            self.hy.ip = ip;
            self.hy.ir = ir;

        //--polynomial data.
        } else {
            *yld = a[7];
            let mut term = 1.0f64;
            let nc = nint(a[5]);
            for i in 2..=nc {
                term *= e;
                *yld += a[i + 6] * term;
            }
            *enext = emax;
            *idis = 0;
        }
    }

    // =============================================================== anadam
    /// Damage energy for an analytic section of File 5 (evaporation only).
    fn anadam(&mut self, c: &Ctx, e: f64, d: &mut f64, z: f64, awr: f64, a: &mut F) {
        const IMAX: usize = 10;
        let eps = 0.05e0f64;
        let delta = 1.0e-7f64;
        let mut x = [0.0f64; IMAX + 1];
        let mut y = [0.0f64; IMAX + 1];

        //--retrieve spectrum temperature
        *d = 0.0;
        let u = a[1];
        let nr = nint(a[5]);
        let np = nint(a[6]);
        let lf = nint(a[4]);
        let awfac = 1.0 / awr;
        if lf != 9 {
            return;
        }
        if e <= u {
            return;
        }
        if (e - u) < delta * e {
            return;
        }
        let test = 10.0;
        if (e - u) < test {
            return;
        }
        let lnext = 7 + 2 * nr + 2 * np;
        let mut ip = 2;
        let mut ir = 1;
        let mut theta = 0.0f64;
        let mut enext = 0.0f64;
        let mut idis = 0i32;
        c.terpa(&mut theta, e, &mut enext, &mut idis, &a.at(lnext), &mut ip, &mut ir);

        //--adaptive energy integration
        //--gauss-legendre angle integration
        let mut i = 4usize;
        let mut ep = 0.0f64;
        while i > 1 {
            i -= 1;
            if i == 3 {
                ep = 1.0;
            }
            if i == 2 {
                ep = (e - u) / 2.0;
            }
            if i == 2 && theta < e - u {
                ep = theta;
            }
            if i == 1 {
                ep = e - u;
            }
            x[i] = ep;
            y[i] = 0.0;
            for iq in 0..4 {
                let er = (e - 2.0 * sqrt(e * ep) * QP4[iq] + ep) * awfac;
                y[i] += QW4[iq] * self.df(er, z, awr, z, awr) / 2.0;
            }
            y[i] *= sed(e, ep, theta, u);
        }
        i = 3;
        let mut j = 0;
        let mut xlast = 0.0f64;
        let mut ylast = 0.0f64;
        while i > 1 {
            // test for convergence
            let mut iflag = 0;
            let mut xm = 0.0f64;
            let mut yt = 0.0f64;
            if i > 1 && i < IMAX {
                xm = (x[i - 1] + x[i]) / 2.0;
                let ym = (y[i - 1] + y[i]) / 2.0;
                yt = 0.0;
                for iq in 0..4 {
                    let er = (e - 2.0 * sqrt(e * xm) * QP4[iq] + xm) * awfac;
                    yt += QW4[iq] * self.df(er, z, awr, z, awr) / 2.0;
                }
                yt *= sed(e, xm, theta, u);
                if (yt - ym).abs() > eps * yt {
                    iflag = 1;
                }
            }
            if iflag == 1 {
                // fails test: add midpoint to stack and continue.
                i += 1;
                x[i] = x[i - 1];
                x[i - 1] = xm;
                y[i] = y[i - 1];
                y[i - 1] = yt;
            } else {
                // passes test: use top point off the stack.
                j += 1;
                if j > 1 {
                    *d += (x[i] - xlast) * (y[i] + ylast) / 2.0;
                }
                xlast = x[i];
                ylast = y[i];
                i -= 1;
            }
        }
    }

    // =============================================================== tabdam
    /// Damage energy for a tabulated section of File 5.
    fn tabdam(&mut self, e: f64, d: &mut f64, z: f64, awr: f64, awfac: f64, a: &F) {
        let zero = 0.0f64;

        //--initialize
        if e == zero {
            *d = self.df(e, z, awr, z, awr);
            return;
        }

        //--trapazoidal energy integration
        //--gauss-legendre angle integration
        *d = 0.0;
        let nr = nint(a[5]);
        let np = nint(a[6]);
        let ibase = 6 + 2 * nr;
        let mut xh = a[ibase + 1];
        let mut _yh = a[ibase + 2];
        let mut fl = 0.0f64;
        for i in 1..=np {
            let xl = xh;
            xh = a[ibase + 2 * i - 1];
            let yh = a[ibase + 2 * i];
            _yh = yh;
            let mut f = 0.0f64;
            for iq in 0..4 {
                let er = (e - 2.0 * sqrt(e * xh) * QP4[iq] + xh) * awfac;
                f += QW4[iq] * self.df(er, z, awr, z, awr) / 2.0;
            }
            f *= yh;
            if i > 1 {
                *d += (xh - xl) * (f + fl) / 2.0;
            }
            fl = f;
        }
    }
}

/// `anabar`: mean energy for an analytic section of File 5.
fn anabar(c: &Ctx, s: &mut f64, e: f64, a: &mut F, t1: f64, t2: f64) {
    let thrhaf = 1.5e0f64;
    let zero = 0.0f64;
    let rpi2 = sqrt(PI) / 2.0;
    let r = |z: f64| -> f64 {
        let a0 = 0.3275911e0f64;
        1.0 / (1.0 + a0 * z.abs())
    };
    let rerfc = |z: f64| -> f64 {
        let a1 = 0.254829592e0f64;
        let a2 = -0.284496736e0f64;
        let a3 = 1.421413741e0f64;
        let a4 = -1.453152027e0f64;
        let a5 = 1.061405429e0f64;
        r(z) * (a1 + r(z) * (a2 + r(z) * (a3 + r(z) * (a4 + a5 * r(z)))))
    };

    //--retrieve common factors.
    *s = 0.0;
    let u = a[1];
    if e <= u {
        return;
    }
    let lf = nint(a[4]);
    let nr = nint(a[5]);
    let np = nint(a[6]);
    let lnext = 7 + 2 * nr + 2 * np;
    let mut ip = 2;
    let mut ir = 1;
    let mut theta = 0.0f64;
    let mut enext = 0.0f64;
    let mut idis = 0i32;
    c.terpa(&mut theta, e, &mut enext, &mut idis, &a.at(lnext), &mut ip, &mut ir);

    if lf == 3 {
        //--law 3--discrete energy spectrum.
        *s = e * t1 - t2 * theta;
    } else if lf == 5 {
        //--law 5--general evaporation spectrum.
        let nrn = nint(a[lnext + 4]);
        let npn = nint(a[lnext + 5]);
        let lnxt = lnext + 2 * nrn + 2 * npn + 6;
        tabbar_ro(s, &mut a.at(lnxt), lf);
        *s *= theta;
    } else if lf == 7 {
        //--law 7--simple fission spectrum.
        if theta != zero {
            let b1 = (e - u) / theta;
            let b2 = sqrt(b1);
            let b3 = exp(-b1);
            let expa = exp(-b2 * b2);
            *s = theta * (thrhaf - b1 * b2 * b3 / (rpi2 * (1.0 - expa * rerfc(b2)) - b2 * b3));
        }
    } else if lf == 9 {
        //--law 9--evaporation spectrum.
        if theta != zero {
            let b1 = (e - u) / theta;
            let _b2 = sqrt(b1);
            let b3 = exp(-b1);
            let mut test = 1.0f64;
            test /= 1000.0;
            if b1 >= test {
                *s = theta * (2.0 - b1 * b1 * b3 / (1.0 - (b1 + 1.0) * b3));
            } else {
                *s = 4.0 * (e - u) / 3.0;
            }
        }
    } else if lf == 11 {
        //--law 11--energy-dependent watt spectrum.
        let nrn = nint(a[lnext + 4]);
        let npn = nint(a[lnext + 5]);
        let lnxt = lnext + 6 + 2 * nrn + 2 * npn;
        ip = 2;
        ir = 1;
        let mut bb = 0.0f64;
        c.terpa(&mut bb, e, &mut enext, &mut idis, &a.at(lnxt), &mut ip, &mut ir);
        *s = theta * (thrhaf + theta * bb / 4.0);
    } else if lf == 12 {
        //--law 12--energy dependent fission neutron spectrum (madland-nix)
        *s = ((a[lnext] + a[lnext + 1]) / 2.0) + (4.0 * theta / 3.0);
    }
}

/// `sed`: secondary energy distribution for anadam.
fn sed(e: f64, ep: f64, theta: f64, u: f64) -> f64 {
    let xeu = (e - u) / theta;
    let aeu = xeu.abs();
    let mut test = 1.0f64;
    test /= 10000.0;
    if aeu >= test {
        ep * exp(-ep / theta) / (theta * theta * (1.0 - exp(-xeu) * (1.0 + xeu)))
    } else {
        ep * exp(-ep / theta) / (theta * theta * xeu * xeu * (1.0 - xeu) / 2.0)
    }
}

/// `tabbar` (stops with the Fortran message for an illegal law).
fn tabbar(c: &mut Ctx, f: &mut f64, a: &mut F, law: i32) {
    if law > 0 && law != 1 && law != 5 {
        c.error("tabbar", "coded for lf=1 and lf=5 only.", " ");
    }
    tabbar_ro(f, a, law);
}

/// `tabbar` without the law check: mean energy of a tabulated section of
/// File 5 (law<0: a File 6 subsection with interpolation |law|).
fn tabbar_ro(f: &mut f64, a: &mut F, law: i32) {
    let nr;
    let np;
    let ibase;
    let mut ir = 0;
    let mut nbt = 0;
    let mut inn;
    let ncyc;
    if law >= 0 {
        *f = 0.0;
        nr = nint(a[5]);
        np = nint(a[6]);
        ibase = 6 + 2 * nr;
        ir = 1;
        nbt = nint(a[7]);
        inn = nint(a[8]);
        ncyc = 2;
    } else {
        *f = 0.0;
        np = nint(a[6]);
        inn = law.abs();
        ncyc = nint(a[4]) + 2;
        ibase = 6;
    }
    let mut xh = a[ibase + 1];
    let mut xhsq = xh * xh;
    let mut yh = a[ibase + 2];

    //--accumulate contributions to integral.
    for i in 2..=np {
        let xl = xh;
        xh = a[ibase + ncyc * (i - 1) + 1];
        let xlsq = xhsq;
        xhsq = xh * xh;
        let delxsq = xhsq - xlsq;
        let yl = yh;
        yh = a[ibase + ncyc * (i - 1) + 2];
        if law > 0 && i > nbt {
            ir += 1;
            nbt = nint(a[6 + 2 * ir - 1]);
            inn = nint(a[6 + 2 * ir]);
        }
        if xl != xh {
            if inn == 1 {
                //--y is a constant.
                *f += yl * delxsq / 2.0;
            } else if inn == 2 {
                //--y is linear in x.
                let b = (yh - yl) / (xh - xl);
                *f += (yl - xl * b) * delxsq / 2.0 + b * (xh * xhsq - xl * xlsq) / 3.0;
            } else if inn == 3 {
                //--y is linear in ln(x)
                let alxl = log(xl);
                let alxh = log(xh);
                let b = (yh - yl) / (alxh - alxl);
                *f = *f + (yl - b * alxl + b / 2.0) * delxsq / 2.0 + b * (xhsq * alxh - xlsq * alxl) / 2.0;
            } else if inn == 4 {
                //--ln(y) is linear in x
                let alyl = log(yl);
                let alyh = log(yh);
                let b = (alyh - alyl) / (xh - xl);
                *f += (yh * (b * xh - 1.0) - yl * (b * xl - 1.0)) / powi(b, 2);
            } else if inn == 5 {
                //--ln(y) is linear in ln(x)
                let alxl = log(xl);
                let alxh = log(xh);
                let alyl = log(yl);
                let alyh = log(yh);
                let b = (alyh - alyl) / (alxh - alxl);
                *f += yl * xlsq * (pow(xh / xl, b + 2.0) - 1.0) / (b + 2.0);
            }
        }
    }
}

impl Heatr {
    // =============================================================== sixbar
    /// Charged-particle ebar and damage energy for reactions in ENDF-6
    /// File 6 (e=0 to initialize each subsection). irec is the index of the
    /// driving particle for recoils (0 for normal particles).
    fn sixbar(
        &mut self,
        c: &mut Ctx,
        e: f64,
        ebar: &mut f64,
        yld: &mut f64,
        dame: &mut f64,
        nin: i32,
        cc: &mut RV,
        _ncmax: i32,
        nscr: i32,
        b: &mut RV,
        mut nbmax: i32,
        _n2: i32,
        j6: i32,
        irec: &mut i32,
        jrec: &mut i32,
        iflag: &mut i32,
    ) {
        let small = 1.0e-10f64;
        let emax = 1.0e10f64;
        let zero = 0.0f64;
        let (mut nb, mut nw) = (0i32, 0i32);
        let mut l = 1i32;
        let mut enext = 0.0f64;
        let mut idis = 0i32;

        let mut lbl;
        if e > zero {
            lbl = 300;
        } else if *irec == 0 {
            lbl = 110;
        } else {
            if *jrec > 0 {
                *irec = *jrec;
            }
            lbl = 100;
        }
        loop {
            match lbl {
                100 => {
                    //--for recoils, back up to the driving particle,
                    //--or restore position after doing recoil.
                    let matd = c.math;
                    let mfd = 6;
                    let mtd = c.mth;
                    c.skiprz(nin, -2);
                    c.findf(matd, mfd, mtd, nin);
                    l = 1;
                    c.contio(nin, 0, 0, &mut cc.at(l), &mut nb, &mut nw);
                    let mut ik = 0;
                    while ik < *irec - 1 {
                        ik += 1;
                        c.tab1io(nin, 0, 0, &mut cc.at(l), &mut nb, &mut nw);
                        l += nw;
                        while nb != 0 {
                            c.moreio(nin, 0, 0, &mut cc.at(l), &mut nb, &mut nw);
                            l += nw;
                        }
                        self.sb.law = c.l2h;
                        let law = self.sb.law;
                        c.skip6(nin, 0, 0, &mut cc.at(l), law);
                    }
                    if *jrec == *irec {
                        *irec = 0;
                        *jrec = 0;
                    }
                    lbl = 110;
                }
                110 => {
                    //--start reading in the desired subsection
                    l = 1;
                    c.tab1io(nin, 0, 0, &mut cc.at(l), &mut nb, &mut nw);
                    l += nw;
                    while nb != 0 {
                        c.moreio(nin, 0, 0, &mut cc.at(l), &mut nb, &mut nw);
                        l += nw;
                    }
                    self.zap = c.c1h;
                    self.awp = c.c2h;
                    self.sb.law = c.l2h;

                    //--error for awp=0. for non photons
                    if self.zap != 0.0 && self.awp == 0.0 {
                        let strng = format!("awp is zero for for mt{} particle {}", ffmt::i(c.mth, 3), ffmt::i(nint(self.zap), 5));
                        c.error("sixbar", &strng, " ");
                    }
                    *iflag = 0;
                    self.sb.disc102 = 0.0;
                    if self.zap == zero {
                        *iflag = 1;
                        self.sb.disc102 = self.awp;
                        self.awp = 0.0;
                    }
                    *ebar = 0.0;
                    *yld = 0.0;
                    if self.sb.law > 0 {
                        lbl = 150;
                        continue;
                    }

                    //--warning for law=0.
                    let strng = format!("no distribution for mt{} particle {}", ffmt::i(c.mth, 3), ffmt::i(nint(self.zap), 5));
                    c.mess("sixbar", &strng, " ");
                    *ebar = -2.0;
                    return;
                }
                150 => {
                    //--watch for 2-body recoils
                    if self.sb.law != 4 {
                        lbl = 210;
                        continue;
                    }
                    *irec = j6 - 1;
                    *jrec = j6 + 1;
                    lbl = 100;
                }
                210 => {
                    //--continue processing this subsection
                    let s = &mut self.sb;
                    s.zp = int(self.zap / 1000.0) as f64;
                    s.zt = int(self.zat / 1000.0) as f64;
                    if *irec > 0 {
                        s.zp = s.zt - s.zp;
                    }
                    s.ap = self.awp;
                    if *irec > 0 {
                        s.ap = self.awr + 1.0 - self.awp;
                    }
                    s.at = self.awr;
                    if self.zap == zero {
                        s.ap = self.awr + 1.0;
                        s.zp = s.zt;
                    }
                    let (zp, ap, zt, at) = (s.zp, s.ap, s.zt, s.at);
                    *dame = self.df(e, zp, ap, zt, at);
                    if self.sb.disc102 > zero {
                        // 295
                        *yld = 1.0;
                        let law = self.sb.law;
                        c.skip6(nin, 0, 0, &mut cc.at(l), law);
                        return;
                    }
                    let law = self.sb.law;
                    if law != 3 && law != 6 {
                        c.tab2io(nin, 0, 0, &mut cc.at(l), &mut nb, &mut nw);
                        self.sb.lang = nint(cc[l + 2]);
                        self.sb.lep = nint(cc[l + 3]);
                        self.sb.ne = c.n2h;
                        let nnt = l + 6;
                        self.sb.intl = nint(cc[nnt + 1]);
                        self.sb.intl = 2;
                        l += nw;
                    }
                    self.sb.nne = 0;
                    self.sb.iraw = l;
                    if law == 3 {
                        //--law 3 -- isotropic distribution
                        self.sb.ehi = emax;
                        self.sb.lang = 0;
                        return;
                    }
                    if law == 6 {
                        //--law 6 -- phase-space distribution
                        c.contio(nin, 0, 0, &mut cc.at(l), &mut nb, &mut nw);
                        self.sb.lang = 0;
                        self.sb.lep = 2;
                        l += nw;
                        self.sb.ehi = emax;
                        self.sb.dhi = 0.0;
                        self.sb.fhi = 0.0;
                        return;
                    }
                    if law == 7 {
                        //--law 7 -- laboratory angle-energy law
                        c.tab2io(nin, 0, 0, &mut cc.at(l), &mut nb, &mut nw);
                        let nmu = c.n2h;
                        l += nw;
                        for _imu in 1..=nmu {
                            c.tab1io(nin, 0, 0, &mut cc.at(l), &mut nb, &mut nw);
                            l += nw;
                            while nb != 0 {
                                c.moreio(nin, 0, 0, &mut cc.at(l), &mut nb, &mut nw);
                                l += nw;
                            }
                        }
                    } else {
                        //--laws 1, 2, and 5.
                        c.listio(nin, 0, 0, &mut cc.at(l), &mut nb, &mut nw);
                        l += nw;
                        while nb != 0 {
                            c.moreio(nin, 0, 0, &mut cc.at(l), &mut nb, &mut nw);
                            l += nw;
                        }
                    }
                    let iraw = self.sb.iraw;
                    self.sb.nne += 1;
                    self.sb.elo = cc[iraw + 1];
                    self.sb.ehi = 0.0;
                    self.sb.dhi = 0.0;
                    self.sb.fhi = 0.0;

                    //--compute low values of ebar and dame
                    // 290
                    let elo = self.sb.elo;
                    let (mut flo, mut dlo) = (0.0f64, 0.0f64);
                    if c.mth != 102 || (self.zap == 0.0 && *irec == 0) {
                        let (law, mut lang, lep) = (self.sb.law, self.sb.lang, self.sb.lep);
                        self.getsix(c, elo, &mut flo, &mut dlo, &mut cc.at(iraw), law, &mut lang, lep, *irec);
                        self.sb.lang = lang;
                    } else {
                        let ztt = int(self.zat / 1000.0) as f64;
                        let (law, awrt) = (self.sb.law, self.awrt);
                        self.tabsq6(c, &mut flo, &mut dlo, cc, iraw, law, ztt, awrt, elo);
                    }
                    self.sb.flo = flo;
                    self.sb.dlo = dlo;
                    if is_fis(c.mth) {
                        let matd = c.math;
                        let mf1 = 1;
                        let mut mt1 = 452;
                        if self.nply > 0 {
                            mt1 = 456;
                        }
                        // delayed neutrons are treated the same as fast
                        // neutrons; when nply.ne.0, use prompt data only.
                        self.hgtyld(c, e, &mut enext, &mut idis, yld, matd, mf1, mt1, nscr, b, &mut nbmax);
                    }
                    return;
                }
                300 => {
                    //--normal entry
                    if self.sb.disc102 > zero {
                        // 430: discrete relativistic capture gamma
                        let s = &self.sb;
                        let (disc102, zp, ap, zt, at) = (s.disc102, s.zp, s.ap, s.zt, s.at);
                        self.hgam102(e, ebar, dame, disc102, *irec, zp, ap, zt, at);
                        *yld = 1.0;
                        return;
                    }

                    //--interpolate for particle yield
                    let mut ip = 2;
                    let mut ir = 1;
                    let mut pe = 0.0f64;
                    let mut eihi = 0.0f64;
                    let mut idisc = 0i32;
                    c.terpa(&mut pe, e, &mut eihi, &mut idisc, &cc.at(1), &mut ip, &mut ir);
                    if *irec > 0 {
                        pe = 1.0;
                    }
                    loop {
                        // 305: is desired energy in current panel
                        let s = &self.sb;
                        let go400 = if e < s.ehi * (1.0 - small) {
                            true
                        } else if s.nne == 1 {
                            false
                        } else if s.nne == s.ne && e <= s.ehi * (1.0 + small) {
                            true
                        } else if s.nne == s.ne {
                            // 450: return zeros outside range of table
                            *ebar = 0.0;
                            *dame = 0.0;
                            return;
                        } else {
                            //--no. slide high-energy data into low-energy
                            //--positions and read in new data for high energy
                            let s = &mut self.sb;
                            s.elo = s.ehi;
                            s.flo = s.fhi;
                            s.dlo = s.dhi;
                            false
                        };
                        if go400 {
                            break;
                        }
                        // 310
                        let law = self.sb.law;
                        let iraw = self.sb.iraw;
                        if law == 1 || law == 2 || law == 5 {
                            l = iraw;
                            c.listio(nin, 0, 0, &mut cc.at(l), &mut nb, &mut nw);
                            l += nw;
                            while nb != 0 {
                                c.moreio(nin, 0, 0, &mut cc.at(l), &mut nb, &mut nw);
                                l += nw;
                            }
                            self.sb.nne += 1;
                            self.sb.ehi = cc[iraw + 1];
                        } else if law == 6 {
                            l = iraw;
                            c.contio(nin, 0, 0, &mut cc.at(l), &mut nb, &mut nw);
                        } else if law == 7 {
                            l = iraw;
                            c.tab2io(nin, 0, 0, &mut cc.at(l), &mut nb, &mut nw);
                            let nmu = c.n2h;
                            l += nw;
                            for _imu in 1..=nmu {
                                c.tab1io(nin, 0, 0, &mut cc.at(l), &mut nb, &mut nw);
                                l += nw;
                                while nb != 0 {
                                    c.moreio(nin, 0, 0, &mut cc.at(l), &mut nb, &mut nw);
                                    l += nw;
                                }
                            }
                            self.sb.ehi = cc[iraw + 1];
                            self.sb.nne += 1;
                        }

                        //--compute high values of ebar and dame
                        let ehi = self.sb.ehi;
                        let (mut fhi, mut dhi) = (0.0f64, 0.0f64);
                        if c.mth != 102 || (self.zap == 0.0 && *irec == 0) {
                            let (law, mut lang, lep) = (self.sb.law, self.sb.lang, self.sb.lep);
                            self.getsix(c, ehi, &mut fhi, &mut dhi, &mut cc.at(iraw), law, &mut lang, lep, *irec);
                            self.sb.lang = lang;
                        } else {
                            let ztt = int(self.zat / 1000.0) as f64;
                            let awrt = self.awrt;
                            self.tabsq6(c, &mut fhi, &mut dhi, cc, iraw, law, ztt, awrt, ehi);
                        }
                        self.sb.fhi = fhi;
                        self.sb.dhi = dhi;
                    }

                    //--yes.  compute result.
                    // 400
                    let law = self.sb.law;
                    let iraw = self.sb.iraw;
                    if law == 3 || law == 6 {
                        cc[iraw + 1] = e;
                        let (mut f, mut d) = (0.0f64, 0.0f64);
                        let (mut lang, lep) = (self.sb.lang, self.sb.lep);
                        self.getsix(c, e, &mut f, &mut d, &mut cc.at(iraw), law, &mut lang, lep, *irec);
                        self.sb.lang = lang;
                        *ebar = f;
                        *yld = pe;
                        *dame = pe * d;
                        return;
                    }

                    //--yes. interpolate for results
                    // 420
                    if self.sb.nne == 1 {
                        *ebar = 0.0;
                        *dame = 0.0;
                        return;
                    }
                    let mut s = 0.0f64;
                    let sb = &self.sb;
                    c.terp1(sb.elo, sb.flo, sb.ehi, sb.fhi, e, &mut s, sb.intl);
                    *ebar = s;
                    *yld = pe;
                    if is_fis(c.mth) {
                        let matd = c.math;
                        let mf1 = 1;
                        let mt1 = 452;
                        self.hgtyld(c, e, &mut enext, &mut idis, yld, matd, mf1, mt1, nscr, b, &mut nbmax);
                    }
                    let mut d = 0.0f64;
                    let sb = &self.sb;
                    c.terp1(sb.elo, sb.dlo, sb.ehi, sb.dhi, e, &mut d, sb.intl);
                    *dame = pe * d;
                    return;
                }
                _ => unreachable!(),
            }
        }
    }

    // =============================================================== getsix
    /// Ebar and damage for one incident energy from File 6 data in cv (the
    /// emitted particle for irec=0, the recoil for irec>0). lang is changed
    /// as the Fortran argument is.
    fn getsix(
        &mut self,
        c: &mut Ctx,
        e: f64,
        ebar: &mut f64,
        dame: &mut f64,
        cv: &mut F,
        law: i32,
        lang: &mut i32,
        lep: i32,
        irec: i32,
    ) {
        const IMAX: usize = 10;
        let tol = 0.02e0f64;
        let small = 1.0e-8f64;
        let emax = 1.0e10f64;
        let zero = 0.0f64;
        let one = 1.0f64;
        let mut p = [0.0f64; 66];

        //--branch to appropriate law
        let mut zp = int(self.zap / 1000.0) as f64;
        let zt = int(self.zat / 1000.0) as f64;
        if irec > 0 {
            zp = zt - zp;
        }
        let mut ap = self.awp;
        if irec > 0 {
            ap = self.awrt + 1.0 - self.awp;
        }
        let at = self.awrt;
        let awp = self.awp;
        let awrt = self.awrt;
        let zap = self.zap;

        if law == 7 {
            //--law 7: trapazoidal integration over energy for each angle,
            //--then trapazoidal integration over angle
            *ebar = 0.0;
            *dame = 0.0;
            let mut dl = 0.0f64;
            let mut hl = 0.0f64;
            let mut ul = 0.0f64;
            let nmu = nint(cv[6]);
            let mut l = 7 + 2 * nint(cv[5]);
            let mut xl = 0.0f64;
            let mut el = 0.0f64;
            let mut fl = 0.0f64;
            for imu in 1..=nmu {
                let u = cv[l + 1];
                let nr = nint(cv[l + 4]);
                let np = nint(cv[l + 5]);
                let iint = nint(cv[l + 7]);
                let next = l + 6 + 2 * nr + 2 * np;
                let ibase = l + 5 + 2 * nr;
                let mut h = 0.0f64;
                let mut d = 0.0f64;
                for i in 1..=np {
                    let mut xx = cv[ibase + 2 * i - 1];
                    let yy = cv[ibase + 2 * i];
                    let en = yy * xx;
                    if irec > 0 {
                        xx = (e - 2.0 * sqrt(e * awp * xx) * u + awp * xx) / (awrt + 1.0 - awp);
                    }
                    if i > 1 {
                        if iint == 1 {
                            h += (xx - xl) * el;
                        } else {
                            h += (xx - xl) * (en + el) / 2.0;
                        }
                    }
                    let f = yy * self.df(xx, zp, ap, zt, at);
                    if i > 1 {
                        if iint == 1 {
                            d += (xx - xl) * fl;
                        } else {
                            d += (xx - xl) * (f + fl) / 2.0;
                        }
                    }
                    xl = xx;
                    el = en;
                    fl = f;
                }
                l = next;
                if imu > 1 {
                    *ebar += (u - ul) * (h + hl) / 2.0;
                }
                if imu > 1 {
                    *dame += (u - ul) * (d + dl) / 2.0;
                }
                ul = u;
                hl = h;
                dl = d;
            }
            return;
        }
        let lab = (law == 1 && self.lct == 1) || (law == 1 && self.lct == 3 && zap > 2004.0) || (law == 1 && zap == 0.0);
        if lab {
            //--law 1 for lab system.
            //--trapazoidal energy integration, and
            //--gauss-legendre angle integration, if needed.
            *ebar = 0.0;
            *dame = 0.0;
            let nd = nint(cv[3]);
            let na = nint(cv[4]);
            let nep = nint(cv[6]);
            let ncyc = na + 2;
            let mut h = 0.0f64;
            let mut d = 0.0f64;
            let mut xl = 0.0f64;
            let mut yl = 0.0f64;
            let mut el = 0.0f64;
            let mut fl = 0.0f64;
            // read before set if lep is neither 1 nor >1 (Fortran locals)
            let mut t1 = 0.0f64;
            let mut t2 = 0.0f64;
            for i in 1..=nep {
                let l = 7 + ncyc * (i - 1);
                let xx = cv[l].abs();
                let yy = cv[l + 1];
                let mut x2 = xx;
                let en;
                let f;
                let f2;
                if irec == 0 {
                    en = yy * xx;
                    let zpp = int(zap / 1000.0) as f64;
                    let ztt = int(self.zat / 1000.0) as f64;
                    f2 = self.df(xx, zpp, awp, ztt, awrt);
                    f = yy * f2;
                } else {
                    let mut enn = 0.0f64;
                    let mut ff = 0.0f64;
                    for iq in 0..64 {
                        let mut bb = 0.0f64;
                        legndr(QP64[iq], &mut p[1..], na);
                        for ia in 1..=na + 1 {
                            bb += (2 * ia - 1) as f64 * cv[l + ia] * p[ia as usize] / 2.0;
                        }
                        if yy != zero {
                            bb /= yy;
                        }
                        let er = (e - 2.0 * sqrt(e * awp * xx) * QP64[iq] + awp * xx) / (awrt + 1.0 - awp);
                        enn += QW64[iq] * bb * er;
                        ff += QW64[iq] * bb * self.df(er, zp, ap, zt, at);
                    }
                    x2 = enn;
                    en = enn * yy;
                    f2 = ff;
                    f = ff * yy;
                }
                if i <= nd {
                    h += en;
                    d += f;
                } else if i == nd + 1 {
                    xl = xx;
                    yl = yy;
                    el = en;
                    fl = f;
                } else {
                    if lep == 1 {
                        t1 = x2 * yl + el;
                    }
                    if lep > 1 {
                        t1 = en + el;
                    }
                    h += (xx - xl) * t1 / 2.0;
                    if lep == 1 {
                        t2 = f2 * yl + fl;
                    }
                    if lep > 1 {
                        t2 = f + fl;
                    }
                    d += (xx - xl) * t2 / 2.0;
                    xl = xx;
                    yl = yy;
                    el = en;
                    fl = f;
                }
            }
            *ebar = h;
            *dame = d;
            return;
        }
        if law == 2 || law == 3 {
            //--laws 2 and 3: discrete-level scattering
            *ebar = 0.0;
            *dame = 0.0;
            *lang = 0;
            if law == 2 {
                *lang = nint(cv[3]);
            }
            let mut nld = 0;
            if law == 2 {
                nld = nint(cv[6]);
            }
            let thresh = ((awrt + 1.0) / awrt) * (-self.q);
            if e < thresh {
                return;
            }
            let mut beta = sqrt(awrt * (awrt + 1.0 - awp) * (1.0 - thresh / e) / awp);
            let mut afact = awp / powi(1.0 + awrt, 2);
            let arec = awp / (awrt + 1.0 - awp);
            if *lang == 0 {
                // legendre coefficients; gauss-legendre integration for damage
                let mut f1 = 0.0f64;
                if law != 3 {
                    f1 = cv[7];
                }
                if irec > 0 {
                    f1 = -f1;
                    beta *= arec;
                    afact /= arec;
                }
                *ebar = afact * e * (1.0 + 2.0 * beta * f1 + beta * beta);
                if self.idame == 0 {
                    return;
                }
                for iq in 0..64 {
                    let mut u = QP64[iq];
                    if irec > 0 {
                        u = -u;
                    }
                    let mut f = one / 2.0;
                    if law != 3 {
                        legndr(u, &mut p[1..], nld);
                        for il in 1..=nld {
                            f += (2 * il + 1) as f64 * cv[6 + il] * p[(il + 1) as usize] / 2.0;
                        }
                    }
                    let e2 = afact * e * (1.0 + 2.0 * beta * u + beta * beta);
                    *dame += QW64[iq] * f * self.df(e2, zp, ap, zt, at);
                }
            } else {
                // tabulated distribution; trapazoidal integration
                let nmu = nint(cv[6]);
                if irec > 0 {
                    beta *= arec;
                    afact /= arec;
                }
                let mut el = 0.0f64;
                let mut fl = 0.0f64;
                let mut ul = 0.0f64;
                for imu in 1..=nmu {
                    let mut l = 6 + 2 * imu - 1;
                    if irec > 0 {
                        l = 6 + 2 * nmu - 2 * imu + 1;
                    }
                    let mut u = cv[l];
                    if irec > 0 {
                        u = -u;
                    }
                    let yy = cv[l + 1];
                    let e2 = afact * e * (1.0 + 2.0 * beta * u + beta * beta);
                    let en = yy * e2;
                    let fn_ = yy * self.df(e2, zp, ap, zt, at);
                    if imu > 1 {
                        *ebar += (u - ul) * (en + el) / 2.0;
                        *dame += (u - ul) * (fn_ + fl) / 2.0;
                    }
                    ul = u;
                    el = en;
                    fl = fn_;
                }
            }
            return;
        }
        if law == 6 {
            *lang = 0;
        }

        //--law 1 in the cm system.
        //--adaptive integration over lab energy.
        let mut nl = 1;
        if irec > 0 {
            nl = 2;
        }
        let nlu = nl as usize;
        let mut summ = 0.0f64;
        *ebar = 0.0;
        *dame = 0.0;
        let mut x = [0.0f64; IMAX + 1];
        let mut y = [[0.0f64; 3]; IMAX + 1];
        let mut term = [0.0f64; 3];
        let mut epnext = 0.0f64;
        x[2] = 0.0;
        let lng = *lang;
        self.h6cm(c, x[2], &mut epnext, &mut term, nl, lng, lep, cv);
        for l in 1..=nlu {
            y[2][l] = term[l];
        }
        let mut el = x[2] * y[2][1];
        let mut fl = 0.0f64;
        let mut i;
        'l120: loop {
            x[1] = epnext;
            let x1 = x[1];
            self.h6cm(c, x1, &mut epnext, &mut term, nl, lng, lep, cv);
            for l in 1..=nlu {
                y[1][l] = term[l];
            }
            i = 2usize;
            let dy = (y[1][1] + y[2][1]) / 10000.0;
            // 140 / 160 / 190
            let mut lbl = 140;
            loop {
                if lbl == 140 {
                    if i == IMAX {
                        lbl = 190;
                        continue;
                    }
                    let da = (x[i - 1] - x[i]) * (y[i - 1][1] + y[i][1]) / 2.0;
                    if da < small {
                        lbl = 190;
                        continue;
                    }
                    let xm = (x[i - 1] + x[i]) / 2.0;
                    if xm >= x[i - 1] || xm <= x[i] {
                        lbl = 190;
                        continue;
                    }
                    let mut epn = x[i];
                    self.h6cm(c, xm, &mut epn, &mut term, nl, lng, lep, cv);
                    let ym = (y[i - 1][1] + y[i][1]) / 2.0;
                    let test = 2.0 * tol * ym.abs() + 2.0 * dy;
                    if (term[1] - ym).abs() > test {
                        // 160: fails test
                        i += 1;
                        x[i] = x[i - 1];
                        for l in 1..=nlu {
                            y[i][l] = y[i - 1][l];
                        }
                        x[i - 1] = xm;
                        for l in 1..=nlu {
                            y[i - 1][l] = term[l];
                        }
                        continue;
                    }
                    lbl = 190;
                    continue;
                }
                // 190: passes test
                if i == 1 {
                    return;
                }
                summ += (x[i - 1] - x[i]) * (y[i - 1][1] + y[i][1]) / 2.0;
                let en;
                let f;
                if irec <= 0 {
                    en = y[i - 1][1] * x[i - 1];
                    *ebar += (x[i - 1] - x[i]) * (en + el) / 2.0;
                    f = y[i - 1][1] * self.df(x[i - 1], zp, ap, zt, at);
                    *dame += (x[i - 1] - x[i]) * (f + fl) / 2.0;
                } else {
                    let mut f1 = 0.0f64;
                    if y[i - 1][1] != zero {
                        f1 = y[i - 1][2] / y[i - 1][1];
                    }
                    let mut enn = (e - 2.0 * sqrt(e * awp * x[i - 1]) * f1 + awp * x[i - 1]) / (awrt + 1.0 - awp);
                    f = y[i - 1][1] * self.df(enn, zp, ap, zt, at);
                    enn *= y[i - 1][1];
                    en = enn;
                    *ebar += (x[i - 1] - x[i]) * (en + el) / 2.0;
                    *dame += (x[i - 1] - x[i]) * (f + fl) / 2.0;
                }
                el = en;
                fl = f;
                i -= 1;
                if i > 1 {
                    lbl = 140;
                    continue;
                }
                let test = emax - emax / 100.0;
                if epnext >= test {
                    // back to 190 with i=1: finished
                    let _ = summ;
                    return;
                }
                x[2] = x[1];
                for l in 1..=nlu {
                    y[2][l] = y[1][l];
                }
                continue 'l120;
            }
        }
    }

    // ================================================================= h6cm
    /// Legendre coefficients of the double-differential cross section at
    /// lab secondary energy ep, from File 6 cm data in cnow (ep=0 to
    /// initialize). term is 1-based.
    fn h6cm(&mut self, c: &mut Ctx, ep: f64, epnext: &mut f64, term: &mut [f64; 3], nl: i32, lang: i32, lep: i32, cnow: &mut F) {
        const IMAX: usize = 10;
        let tol = 0.02e0f64;
        let emax = 1.0e10f64;
        let step = 0.05e0f64;
        let check = 0.99999e0f64;
        let rdn = 0.999995e0f64;
        let rup = 1.000005e0f64;
        let tiny = 1.0e-9f64;
        let dn = 0.999999e0f64;
        let up = 1.000001e0f64;
        let thou = 1.0e-3f64;
        let amil = 1.0e-6f64;
        let one = 1.0e0f64;
        let onep5 = 1.5e0f64;
        let du = 0.0001e0f64;
        let small = 1.0e-10f64;
        let zero = 0.0f64;
        let nlu = nl as usize;
        let mut p = [0.0f64; 66];
        let mut epn = 0.0f64;
        let mut epm = 0.0f64;
        let mut w;

        //--initialize.
        if ep <= zero {
            self.h6.e = cnow[2];
            self.h6.ndnow = nint(cnow[3]);
            self.h6.npnow = nint(cnow[6]);
            self.h6.ncnow = nint(cnow[5]) / self.h6.npnow;
            let aprime = self.awp;
            self.h6.xc = aprime / powi(self.awrt + 1.0, 2);
            *epnext = emax;
            self.h6.elmax = 0.0;
            let mut epnn = emax;
            self.h6.epmax = 0.0;
            epm = 0.0;
            epn = 0.0;
            w = 0.0;
            let e = self.h6.e;
            if self.h6.npnow != self.h6.ndnow {
                let _t = if lang > 0 {
                    self.h6ddx(c, ep, &mut epn, &mut epm, w, cnow, lang, lep)
                } else {
                    self.h6psp(c, ep, &mut epn, &mut epm, w, e, cnow)
                };
                if epn < epnn * (1.0 - small) {
                    epnn = epn;
                }
                if epm > self.h6.epmax * (1.0 + small) {
                    self.h6.epmax = epm;
                }
                let s = &mut self.h6;
                s.elmax = e * powi(sqrt(s.epmax / e) + sqrt(s.xc), 2);
                s.elmax *= up;
                epn = step * s.elmax;
                if s.epmax * up < s.xc * e {
                    epn = e * powi(sqrt(s.epmax / e) - sqrt(s.xc), 2);
                }
                if epn < *epnext {
                    *epnext = epn;
                }
                s.eps = tiny / s.elmax;
                *epnext *= dn;
            }
            if self.h6.ndnow != 0 {
                let mut i = 0;
                let _t = self.h6dis(c, i, &mut epn, &mut epm, w, cnow, lang);
                for ii in 1..=self.h6.ndnow {
                    i = ii;
                    let _t = self.h6dis(c, i, &mut epn, &mut epm, w, cnow, lang);
                    if epn < epnn * (1.0 - small) {
                        epnn = epn;
                    }
                    let s = &mut self.h6;
                    if epn > s.epmax * (1.0 + small) {
                        s.epmax = epn;
                    }
                    let mut epx = dn * e * powi(sqrt(epn / e) - sqrt(s.xc), 2);
                    if epx < *epnext * (1.0 - small) {
                        *epnext = epx;
                    }
                    epx = up * e * powi(sqrt(epn / e) + sqrt(s.xc), 2);
                    if epx > s.elmax * (1.0 + small) {
                        s.elmax = epx;
                    }
                }
            }
            self.h6.na = nl - 1;
            for l in 1..=nlu {
                term[l] = 0.0;
            }
            return;
        }

        //--normal entry.
        for l in 1..=nlu {
            term[l] = 0.0;
        }
        let (e, xc, epmax, elmax, eps, na) = (self.h6.e, self.h6.xc, self.h6.epmax, self.h6.elmax, self.h6.eps, self.h6.na);
        // epnxt is read unset when ndnow = npnow = 0 (no data at all)
        let mut epnxt = 0.0f64;
        if self.h6.ndnow == self.h6.npnow {
            // 330: check for contributions from delta functions
            if self.h6.ndnow != 0 {
                epnxt = emax;
                let ncnow = self.h6.ncnow;
                for i in 1..=self.h6.ndnow {
                    let ll = 7 + ncnow * (i - 1);
                    let epp = cnow[ll];
                    w = (ep - epp - xc * e) / (2.0 * sqrt(xc * e) * sqrt(epp));
                    let u = (xc * e + ep - epp) / (2.0 * sqrt(xc * e) * sqrt(ep));
                    let wlo = -1.0;
                    let whi = 1.0;
                    if w >= wlo && w <= whi {
                        legndr(u, &mut p[1..], na);
                        let mut f = self.h6dis(c, i, &mut epn, &mut epm, w, cnow, lang);
                        f /= 2.0 * sqrt(xc * e * epp);
                        epn = e * powi(sqrt(epp / e) + sqrt(xc), 2);
                        if ep < check * epn {
                            epn *= rdn;
                        } else {
                            epn *= rup;
                        }
                        if epn < epnxt {
                            epnxt = epn;
                        }
                        for l in 1..=nlu {
                            term[l] += f * p[l];
                        }
                    } else if w < wlo {
                        epn = e * powi(sqrt(epp / e) - sqrt(xc), 2);
                        if ep < check * epn {
                            epn *= rdn;
                        } else {
                            epn *= rup;
                        }
                        if epn < epnxt {
                            epnxt = epn;
                        }
                    }
                }
            }
            // 400
            *epnext = epnxt;
            return;
        }
        let xx = ep / e;
        let mut cc = xc / xx;
        let mut cq = sqrt(cc);
        let mut umin = (1.0 + cc - epmax / ep) / (2.0 * cq);
        if !(umin < check) {
            if cq <= one {
                // 420
                *epnext = emax;
                return;
            }
            if cq > one {
                // 320
                epnxt = (*epnext + elmax) / 2.0;
                *epnext = epnxt;
                return;
            }
            // (c is NaN: the Fortran falls through to 215)
        }
        // 215
        if umin < -one {
            umin = -one;
        }
        let mut dm = (1.0 - umin) / 1000.0 + amil;
        if dm < thou / 5.0 {
            dm = thou / 5.0;
        }
        let mut epnn = emax;

        //--for the continuum part of the distribution,
        //--integrate over the trajectory ep=constant
        //--using adaptive integration on lab cosine.
        let umax = 1.0f64;
        let mut un = umax;
        let mut x = [0.0f64; IMAX + 1];
        let mut y = [[0.0f64; 66]; IMAX + 1];
        let mut yt = [0.0f64; 66];
        x[2] = un;
        // 220
        while !(x[2] <= umin) {
            let mut j = 2;
            while j == 2 {
                let u = un;
                let mut yy = 1.0 + cc - 2.0 * cq * u;
                let test = amil;
                if yy < test {
                    yy = amil;
                    cq = u - thou;
                    cc = cq * cq;
                }
                let epp = yy * ep;
                w = (u - cq) / sqrt(yy);
                let s = if lang > 0 {
                    self.h6ddx(c, epp, &mut epn, &mut epm, w, cnow, lang, lep)
                } else {
                    self.h6psp(c, epp, &mut epn, &mut epm, w, e, cnow)
                };
                if u == umax && epn < epnn {
                    epnn = epn;
                }
                legndr(u, &mut p[1..], na);
                j = 1;
                if u == umax {
                    j = 2;
                }
                x[j] = u;
                for l in 1..=nlu {
                    y[j][l] = p[l] * s / sqrt(yy);
                }
                un = u - (epn - epp) / (2.0 * cq * ep);
                if un > u - du {
                    un = u - du;
                }
                if un < umin + du / 10.0 {
                    un = umin;
                }
            }
            let mut i = 2usize;
            loop {
                // 240
                let mut pass = true;
                if i != IMAX {
                    let dx = x[i] - x[i - 1];
                    if !(dx < dm) {
                        let da = dx * (y[i][1] + y[i - 1][1]) / 2.0;
                        if !(da < eps) {
                            let u = (x[i - 1] + x[i]) / 2.0;
                            let yy = 1.0 + cc - 2.0 * cq * u;
                            let epp = yy * ep;
                            w = (u - cq) / sqrt(yy);
                            let mut s = 0.0f64;
                            if lang > 0 {
                                s = self.h6ddx(c, epp, &mut epn, &mut epm, w, cnow, lang, lep);
                            }
                            if lang == 0 {
                                s = self.h6psp(c, epp, &mut epn, &mut epm, w, e, cnow);
                            }
                            legndr(u, &mut p[1..], na);
                            for l in 1..=nlu {
                                yt[l] = p[l] * s / sqrt(yy);
                            }
                            for l in 1..=nlu {
                                let ym = (y[i - 1][l] + y[i][l]) / 2.0;
                                if (yt[l] - ym).abs() > l as f64 * tol * ym.abs() + eps {
                                    pass = false;
                                    break;
                                }
                            }
                            if !pass {
                                // 260: fails test
                                i += 1;
                                x[i] = x[i - 1];
                                for l in 1..=nlu {
                                    y[i][l] = y[i - 1][l];
                                }
                                x[i - 1] = u;
                                for l in 1..=nlu {
                                    y[i - 1][l] = yt[l];
                                }
                                continue;
                            }
                        }
                    }
                }
                // 290: passes test
                for l in 1..=nlu {
                    term[l] += (x[i] - x[i - 1]) * (y[i][l] + y[i - 1][l]) / 2.0;
                }
                i -= 1;
                if i > 1 {
                    continue;
                }
                break;
            }
            x[2] = x[1];
            for l in 1..=nlu {
                y[2][l] = y[1][l];
            }
        }
        // 310
        epnxt = e * powi(sqrt(epnn / e) - sqrt(xc), 2);
        if epnxt <= ep * (1.0 + small) {
            epnxt = e * powi(sqrt(epnn / e) + sqrt(xc), 2);
        }
        if epnxt <= ep * (1.0 + small) {
            epnxt = onep5 * ep;
        }
        if epnxt > elmax * (1.0 + small) && elmax > ep * (1.0 + small) {
            epnxt = elmax;
        }
        // 400
        *epnext = epnxt;
    }

    // ================================================================ h6ddx
    /// Continuum double-differential cross section at secondary energy ep
    /// and cm cosine w from a File 6 subsection (ep=0 to initialize).
    fn h6ddx(&mut self, c: &mut Ctx, ep: f64, epnext: &mut f64, epmax: &mut f64, w: f64, cnow: &mut F, lang: i32, lep: i32) -> f64 {
        const NLMAX: i32 = 65;
        let up = 1.00001e0f64;
        let dn = 0.99999e0f64;
        let off = 0.999995e0f64;
        let emax = 1.0e10f64;
        let small = 1.0e-10f64;
        let zero = 0.0f64;
        let mut p = [0.0f64; 66];

        //--initialize
        if ep == zero {
            let s = &mut self.hd;
            s.enow = cnow[2];
            let ndnow = nint(cnow[3]);
            let npnow = nint(cnow[6]);
            s.ncnow = nint(cnow[5]) / npnow;
            s.mnow = 7 + (npnow - 1) * s.ncnow;
            s.inow = 7 + s.ncnow * ndnow;
            s.lnow = s.inow;
            if cnow[s.lnow] <= zero {
                s.lnow += s.ncnow;
            }
            *epnext = cnow[s.lnow];
            *epmax = cnow[s.mnow];
            s.nl = s.ncnow - 1;
            if s.nl > NLMAX {
                c.error("h6ddx", "too many legendre terms", " ");
            }
            s.na = s.nl - 1;
            s.efirst = 0.0;
            s.illdef = 0;
            return 0.0;
        }

        //--check for desired interpolation range
        let s = &mut self.hd;
        let mut t;
        let lbl;
        if ep < s.efirst * (1.0 - small) {
            lbl = 200;
        } else if ep > cnow[s.mnow] {
            lbl = 190;
        } else {
            let mut l = 130;
            loop {
                if l == 130 {
                    *epnext = cnow[s.lnow];
                    if ep < off * *epnext {
                        l = 140;
                        continue;
                    }
                    if s.lnow >= s.mnow {
                        break;
                    }
                    s.lnow += s.ncnow;
                    continue;
                }
                // 140
                if s.lnow <= s.ncnow {
                    l = 200;
                    break;
                }
                let eplast = cnow[s.lnow - s.ncnow];
                if ep >= off * eplast {
                    break;
                }
                if s.lnow <= s.inow + s.ncnow {
                    break;
                }
                s.lnow -= s.ncnow;
                l = 130;
            }
            lbl = if l == 200 { 200 } else { 210 };
        }
        if lbl == 190 {
            //--return zero above table
            t = 0.0;
            *epnext = emax;
        } else if lbl == 200 {
            //--use analytic form below efirst.
            t = 0.0;
        } else if lang == 1 {
            //--legendre coefficients
            let (nl, na, ncnow, lnow, inow, mnow) = (s.nl, s.na, s.ncnow, s.lnow, s.inow, s.mnow);
            legndr(w, &mut p[1..], na);
            t = 0.0;
            for l in 1..=nl {
                let mut tt = 0.0f64;
                if l <= ncnow - 1 && ep >= cnow[inow] * (1.0 - small) && ep <= cnow[mnow] * (1.0 + small) {
                    let lll = lnow + l;
                    let x1 = cnow[lnow - ncnow];
                    let mut x2 = cnow[lnow];
                    if x1 == x2 && lep > 1 {
                        x2 = sigfig(x2, 6, 1);
                        if self.hd.illdef == 0 {
                            self.hd.illdef = 1;
                        }
                    }
                    let y1 = cnow[lll - ncnow];
                    let y2 = cnow[lll];
                    c.terp1(x1, y1, x2, y2, ep, &mut tt, lep);
                }
                t += p[l as usize] * (2 * l - 1) as f64 * tt / 2.0;
            }
            if t < zero {
                t = 0.0;
            }
        } else if lang == 2 {
            //--kalbach systematics
            let (ncnow, lnow, inow, mnow, enow) = (s.ncnow, s.lnow, s.inow, s.mnow, s.enow);
            let mut ss = 0.0f64;
            let mut r = 0.0f64;
            if ep >= cnow[inow] * (1.0 - small) && ep <= cnow[mnow] * (1.0 + small) {
                let x1 = cnow[lnow - ncnow];
                let mut x2 = cnow[lnow];
                if x1 == x2 && lep > 1 {
                    x2 = sigfig(x2, 6, 1);
                    if self.hd.illdef == 0 {
                        self.hd.illdef = 1;
                    }
                }
                let y1 = cnow[lnow - ncnow + 1];
                let y2 = cnow[lnow + 1];
                c.terp1(x1, y1, x2, y2, ep, &mut ss, lep);
                let y1 = cnow[lnow - ncnow + 2];
                let y2 = cnow[lnow + 2];
                c.terp1(x1, y1, x2, y2, ep, &mut r, lep);
            }
            let iza2 = nint(self.zap);
            let aa = bacha(c, self.izap, iza2, self.izat, enow, ep);
            t = aa * (cosh(aa * w) + r * sinh(aa * w)) / (2.0 * sinh(aa));
            t *= ss;
            if t < zero {
                t = 0.0;
            }
        } else if (11..=15).contains(&lang) {
            //--tabulated angular distribution
            let (ncnow, lnow, na) = (s.ncnow, s.lnow, s.na);
            let inn = lang - 10;
            let ii = lnow - ncnow;
            let jj = lnow;
            let x1 = cnow[ii];
            let mut x2 = cnow[jj];
            if x1 == x2 && lep > 1 {
                x2 = sigfig(x2, 6, 1);
                if self.hd.illdef == 0 {
                    self.hd.illdef = 1;
                }
            }
            let y1 = cnow[ii + 1];
            let y2 = cnow[jj + 1];
            let mut ss = 0.0f64;
            c.terp1(x1, y1, x2, y2, ep, &mut ss, lep);
            let mut tii = 0.0f64;
            let mut ia = 1;
            while ia <= na {
                if w >= cnow[ii + 1 + ia] && w <= cnow[ii + 3 + ia] {
                    c.terp1(cnow[ii + 1 + ia], cnow[ii + 2 + ia], cnow[ii + 3 + ia], cnow[ii + 4 + ia], w, &mut tii, inn);
                }
                ia += 2;
            }
            let mut tjj = 0.0f64;
            ia = 1;
            while ia <= na {
                if w >= cnow[jj + 1 + ia] && w <= cnow[jj + 3 + ia] {
                    c.terp1(cnow[jj + 1 + ia], cnow[jj + 2 + ia], cnow[jj + 3 + ia], cnow[jj + 4 + ia], w, &mut tjj, inn);
                }
                ia += 2;
            }
            t = 0.0;
            c.terp1(x1, tii, x2, tjj, ep, &mut t, lep);
            t *= ss;
        } else {
            //--bad value for lang
            c.error("h6ddx", "illegal lang.", " ");
        }

        //--return final value
        // 350
        if lep > 1 {
            return t;
        }
        if *epnext == emax {
            return t;
        }
        if ep < dn * dn * *epnext {
            *epnext *= dn;
        } else {
            *epnext *= up;
        }
        t
    }

    // ================================================================ h6dis
    /// Double-differential cross section for discrete energy i and cm
    /// cosine w from a File 6 subsection (i=0 to initialize).
    fn h6dis(&mut self, c: &mut Ctx, i: i32, epnext: &mut f64, epmax: &mut f64, w: f64, cnow: &mut F, lang: i32) -> f64 {
        let zero = 0.0f64;
        let mut p = [0.0f64; 66];

        //--initialize
        if i == 0 {
            let s = &mut self.hs;
            s.enow = cnow[2];
            let ndnow = nint(cnow[3]);
            let npnow = nint(cnow[6]);
            s.ncnow = nint(cnow[5]) / npnow;
            s.inow = 7;
            s.lnow = s.inow;
            s.mnow = s.lnow + (ndnow - 1) * s.ncnow;
            *epnext = cnow[s.lnow];
            *epmax = cnow[s.mnow];
            s.nl = s.ncnow - 1;
            s.na = s.nl - 1;
            return 0.0;
        }

        //--select desired discrete line
        let s = &mut self.hs;
        s.inow = 7 + (i - 1) * s.ncnow;
        *epnext = cnow[s.inow];
        let (inow, ncnow, nl, na, enow) = (s.inow, s.ncnow, s.nl, s.na, s.enow);
        let mut t;

        if lang == 1 {
            //--legendre coefficients
            legndr(w, &mut p[1..], na);
            t = 0.0;
            for l in 1..=nl {
                if l <= ncnow - 1 {
                    t += (2 * l - 1) as f64 * p[l as usize] * cnow[inow + 1 + l] / 2.0;
                }
            }
            t *= cnow[inow + 1];
            if t < zero {
                t = 0.0;
            }
        } else if lang == 2 {
            //--kalbach systematics
            let ss = cnow[inow + 1];
            let r = cnow[inow + 2];
            let iza2 = nint(self.zap);
            let aa = bacha(c, self.izap, iza2, self.izat, enow, *epnext);
            t = aa * (cosh(aa * w) + r * sinh(aa * w)) / (2.0 * sinh(aa));
            t *= ss;
            if t < zero {
                t = 0.0;
            }
        } else if (11..=15).contains(&lang) {
            //--tabulated angular distribution
            let inn = lang - 10;
            t = 0.0;
            let mut ia = 1;
            while ia <= na - 2 {
                if w >= cnow[inow + 1 + ia] && w <= cnow[inow + 3 + ia] {
                    c.terp1(cnow[inow + 1 + ia], cnow[inow + 2 + ia], cnow[inow + 3 + ia], cnow[inow + 4 + ia], w, &mut t, inn);
                }
                ia += 2;
            }
            t *= cnow[inow + 1];
        } else {
            //--illegal lang
            c.error("h6dis", "illegal lang.", " ");
        }
        t
    }

    // ================================================================ h6psp
    /// n-body phase-space double-differential cross section (ep=0 to
    /// initialize for each new e).
    fn h6psp(&mut self, c: &mut Ctx, ep: f64, epnext: &mut f64, epmax: &mut f64, _w: f64, e: f64, cv: &mut F) -> f64 {
        let c3 = 1.2732e0f64;
        let c4 = 3.2813e0f64;
        let c5 = 5.8205e0f64;
        let thrhaf = 1.5e0f64;
        let sevhaf = 3.5e0f64;
        let step = 0.05e0f64;
        let rndup = 1.000001e0f64;
        let emax = 1.0e10f64;
        let small = 1.0e-10f64;
        let zero = 0.0f64;

        //--initialize
        if ep == zero {
            let apsx = cv[1];
            let npsx = nint(cv[6]);
            let f1 = (apsx - self.awp) / apsx;
            let f2 = self.awrt / (self.awrt + 1.0);
            let s = &mut self.hp;
            s.ex = thrhaf * npsx as f64 - 4.0;
            s.eimax = f1 * (f2 * e + self.q);
            if s.eimax <= zero {
                s.eimax = 1.0;
            }
            if npsx == 3 {
                s.cn = c3 / powi(s.eimax, 2);
            } else if npsx == 4 {
                s.cn = c4 / pow(s.eimax, sevhaf);
            } else if npsx == 5 {
                s.cn = c5 / powi(s.eimax, 5);
            } else {
                c.error("h6psp", "3, 4, or 5 particles only.", " ");
            }
            let s = &self.hp;
            *epnext = step * s.eimax;
            let test = 1.0;
            if s.eimax == test {
                *epnext = emax;
            }
            *epmax = s.eimax;
            return 0.0;
        }

        //--compute cross section
        let s = &self.hp;
        let mut ss = 0.0f64;
        if ep < s.eimax * (1.0 - small) {
            ss = s.cn * sqrt(ep) * pow(s.eimax - ep, s.ex);
        }
        *epnext = ep + step * s.eimax;
        if *epnext > s.eimax * (1.0 + small) {
            *epnext = rndup * s.eimax;
        }
        if ep >= s.eimax * (1.0 - small) {
            *epnext = emax;
        }
        ss
    }

    // =============================================================== tabsq6
    /// Average photon recoil and damage energy for File 6 capture photons.
    /// The data are cc(iraw...), the photon yield is the TAB1 at cc(1).
    fn tabsq6(&mut self, c: &mut Ctx, g: &mut f64, h: &mut f64, cc: &mut RV, iraw: i32, _law: i32, z: f64, awr: f64, e: f64) {
        let a = |k: i32| iraw + k - 1;
        let ein = 2.0 * self.tm;
        let rein = 1.0 / ein;
        *g = 0.0;
        *h = 0.0;
        let mut s = 0.0f64;
        let nd = nint(cc[a(3)]);
        let np = nint(cc[a(6)]);
        let ncyc = nint(cc[a(5)]) / np;
        let mut ibase = 6;
        let inn = 2;

        //--interpolate the photon yield
        let mut ip = 2;
        let mut ir = 1;
        let mut yield_ = 0.0f64;
        let mut enext = 0.0f64;
        let mut idis = 0i32;
        c.terpa(&mut yield_, e, &mut enext, &mut idis, &cc.at(1), &mut ip, &mut ir);

        //--accumulate contributions from discrete levels
        if nd != 0 {
            for i in 1..=nd {
                let x = cc[a(ibase + ncyc * (i - 1) + 1)];
                let y = cc[a(ibase + ncyc * (i - 1) + 2)];
                let xr = x * x * rein;
                *g += xr * y;
                let awc = awr + 1.0;
                *h += self.df(xr, z, awc, z, awr) * y;
                s += y;
            }
        }

        //--accumulate contributions from continuum
        if np > nd {
            let nc = np - nd;
            ibase += ncyc * nd;
            let mut xh = cc[a(ibase + 1)];
            let mut yh = cc[a(ibase + 2)];
            for i in 2..=nc {
                let xl = xh;
                xh = cc[a(ibase + ncyc * (i - 1) + 1)];
                let yl = yh;
                yh = cc[a(ibase + ncyc * (i - 1) + 2)];
                if xl != xh {
                    let dx = xh - xl;
                    for j in 0..4 {
                        let x = xl + (1.0 + QP4[j]) * dx / 2.0;
                        let mut y = 0.0f64;
                        c.terp1(xl, yl, xh, yh, x, &mut y, inn);
                        let xr = x * x * rein;
                        *g += QW4[j] * xr * dx * y;
                        *h += QW4[j] * self.df(xr, z, awr + 1.0, z, awr) * dx * y;
                        s += QW4[j] * dx * y;
                    }
                }
            }
        }

        //--finished
        *g = yield_ * *g / s;
        *h = yield_ * *h / s;
    }

    // ============================================================== hgam102
    /// The relativistic discrete gamma (or its recoil) of mf6/mt102.
    fn hgam102(&mut self, e: f64, ebar: &mut f64, dame: &mut f64, disc102: f64, irec: i32, zp: f64, ap: f64, zt: f64, at: f64) {
        if irec == 0 {
            *ebar = disc102 + e * self.awr / (1.0 + self.awr);
            *dame = 0.0;
        } else {
            //--include recoil energy plus photon "kick" energy
            let er = e / (self.awr + 1.0);
            let eg2 = disc102 * disc102 * self.rtm / 2.0;
            *ebar = er + eg2;
            *dame = self.df(*ebar, zp, ap, zt, at);
        }
    }
}

/// `bacha`: the Kalbach a parameter.
fn bacha(c: &mut Ctx, iza1i: i32, iza2: i32, izat: i32, e: f64, ep: f64) -> f64 {
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
    let tomev = 1.0e-6f64;
    let zero = 0.0f64;

    let mut iza1 = iza1i;
    if iza1i == 0 {
        iza1 = 1;
    }
    let iza = match izat {
        6000 => 6012,
        12000 => 12024,
        14000 => 14028,
        16000 => 16032,
        17000 => 17035,
        19000 => 19039,
        20000 => 20040,
        22000 => 22048,
        23000 => 23051,
        24000 => 24052,
        26000 => 26056,
        28000 => 28058,
        29000 => 29063,
        31000 => 31069,
        40000 => 40090,
        42000 => 42096,
        48000 => 48112,
        49000 => 49115,
        50000 => 50120,
        63000 => 63151,
        72000 => 72178,
        74000 => 74184,
        82000 => 82208,
        x => x,
    };
    let aa = (iza % 1000) as f64;
    if aa == zero {
        let strng = format!("dominant isotope not known for {}", ffmt::i(iza, 8));
        c.error("bacha", &strng, " ");
    }
    let za = (iza / 1000) as f64;
    let ac = aa + (iza1 % 1000) as f64;
    let zc = za + (iza1 / 1000) as f64;
    let ab = ac - (iza2 % 1000) as f64;
    let zb = zc - (iza2 / 1000) as f64;
    let na = nint(aa - za) as f64;
    let nb = nint(ab - zb) as f64;
    let nc = nint(ac - zc) as f64;
    let mut sa = c1 * (ac - aa)
        + c2 * (powi(nc - zc, 2) / ac - powi(na - za, 2) / aa)
        + c3 * (pow(ac, twoth) - pow(aa, twoth))
        + c4 * (powi(nc - zc, 2) / pow(ac, fourth) - powi(na - za, 2) / pow(aa, fourth))
        + c5 * (powi(zc, 2) / pow(ac, third) - powi(za, 2) / pow(aa, third))
        + c6 * (powi(zc, 2) / ac - powi(za, 2) / aa);
    if iza1 == 1002 {
        sa -= s2;
    }
    if iza1 == 1003 {
        sa -= s3;
    }
    if iza1 == 2003 {
        sa -= s4;
    }
    if iza1 == 2004 {
        sa -= s5;
    }
    let mut sb = c1 * (ac - ab)
        + c2 * (powi(nc - zc, 2) / ac - powi(nb - zb, 2) / ab)
        + c3 * (pow(ac, twoth) - pow(ab, twoth))
        + c4 * (powi(nc - zc, 2) / pow(ac, fourth) - powi(nb - zb, 2) / pow(ab, fourth))
        + c5 * (powi(zc, 2) / pow(ac, third) - powi(zb, 2) / pow(ab, third))
        + c6 * (powi(zc, 2) / ac - powi(zb, 2) / ab);
    if iza2 == 1002 {
        sb -= s2;
    }
    if iza2 == 1003 {
        sb -= s3;
    }
    if iza2 == 2003 {
        sb -= s4;
    }
    if iza2 == 2004 {
        sb -= s5;
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
    b1 * x1 + b2 * powi(x1, 3) + b3 * fa * fb * powi(x3, 4)
}

impl Heatr {
    // =============================================================== hconvr
    /// Convert photon transition probability arrays (LO=2) to yields
    /// (LO=1), add MT456 if necessary, and copy all other sections.
    fn hconvr(&mut self, c: &mut Ctx, nin: i32, nout: i32, nscr: i32) {
        const NQMX: usize = 450;
        let zero = 0.0f64;
        let (mut nb, mut nw) = (0i32, 0i32);
        let imax: i32 = 50;
        let lmax: i32 = 500;
        let mut ngam = [0i32; 301];
        let mut mtq = [0i32; NQMX + 1];
        let mut eeq = [0.0f64; NQMX + 1];
        let mut scr = rv(c.npage + 50);
        // allocated at the first LO=2 section
        let mut istor = 0;
        let mut ev = RV::new(1);
        let mut eg = RV::new(1);
        let mut es = RV::new(1);
        let mut y = RV::new(1);
        let mut aa = RV::new(1);
        let mut r = RV::new(1);
        let mut nu = RV::new(1);
        let mut nnu = 0i32;
        let mut lnu = 0i32;

        //--initialize
        let mut mt0 = 0;
        let mut mt0old = 0;
        let mut l2flg = 0;
        c.nsh = 0;
        c.nsc = 0;

        //--check for consistancy between mf3 mt's and, if present, mf12
        //  mt's (see the Fortran for the rules).
        let mut mtnow = 0;
        let mut mtmess = 0;
        let mut i10 = 0;
        let mut nmf3 = 0usize;
        let mut maths = 0;
        let in_first = |mth: i32, iverf: i32| -> bool {
            mth == 51
                || (iverf >= 6 && (mth == 601 || mth == 651 || mth == 701 || mth == 751 || mth == 801 || mth == 876))
                || (iverf <= 5 && (mth == 701 || mth == 721 || mth == 741 || mth == 761 || mth == 781))
        };
        let in_next = |mth: i32, iverf: i32| -> bool {
            (mth > 51 && mth < 91)
                || (iverf >= 6 && mth > 601 && mth < 648)
                || (iverf >= 6 && mth > 651 && mth < 698)
                || (iverf >= 6 && mth > 701 && mth < 748)
                || (iverf >= 6 && mth > 751 && mth < 798)
                || (iverf >= 6 && mth > 801 && mth < 848)
                || (iverf >= 6 && mth > 876 && mth < 890)
                || (iverf <= 5 && mth > 701 && mth < 718)
                || (iverf <= 5 && mth > 721 && mth < 738)
                || (iverf <= 5 && mth > 741 && mth < 758)
                || (iverf <= 5 && mth > 761 && mth < 778)
                || (iverf <= 5 && mth > 781 && mth < 798)
        };
        loop {
            // 10
            c.contio(nin, 0, 0, &mut scr.all(), &mut nb, &mut nw);
            if i10 == 0 {
                i10 = 1;
                maths = c.math;
            }
            if c.math == 0 {
                break;
            }
            if c.mfh < 3 || (c.mfh > 3 && c.mfh < 12) {
                c.tofend(nin, 0, 0, &mut scr.all());
                continue;
            }
            if c.mfh > 12 {
                break;
            }
            if c.mfh == 3 {
                if in_first(c.mth, c.iverf) {
                    mtnow = c.mth;
                } else if in_next(c.mth, c.iverf) {
                    if mtnow + 1 != c.mth {
                        mtmess = 1;
                        let strng = format!("mf3, mt{} is missing", ffmt::i(c.mth - 1, 2));
                        c.mess("hconvr", &strng, "");
                    }
                    mtnow = c.mth;
                }
                c.contio(nin, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                nmf3 += 1;
                mtq[nmf3] = c.mth;
                eeq[nmf3] = scr[2];
                c.tosend(nin, 0, 0, &mut scr.all());
                continue;
            } else if c.mfh == 12 {
                if in_first(c.mth, c.iverf) {
                    mtnow = c.mth;
                } else if in_next(c.mth, c.iverf) {
                    if mtnow + 1 != c.mth {
                        let strng = format!("mf12, mt{} may be missing", ffmt::i(c.mth - 1, 2));
                        c.mess("hconvr", &strng, "discrete photon data may be incomplete");
                    }
                    mtnow = c.mth;
                }
                c.tosend(nin, 0, 0, &mut scr.all());
                continue;
            }
            if mtmess != 0 {
                c.error("hconvr", "missing mf3 mt's, probable endf error", "");
            }
        }

        //--mt checks are done, rewind the endf tape and position
        //  at the first record for this material
        // 20
        c.repoz(nin);
        c.findf(maths, 1, 451, nin);

        //--loop over all sections of this material
        let mut lbl = 110;
        loop {
            match lbl {
                110 => {
                    c.contio(nin, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                    if c.mfh == 0 && c.math != 0 {
                        // 130
                        c.contio(0, nout, nscr, &mut scr.all(), &mut nb, &mut nw);
                        continue;
                    }
                    if c.mfh == 1 {
                        lbl = 550;
                        continue;
                    }
                    if c.mfh == 12 && c.mth != 460 {
                        // 120
                        if c.l1h == 2 {
                            lbl = 150;
                            continue;
                        }
                        c.contio(0, nout, nscr, &mut scr.all(), &mut nb, &mut nw);
                        c.tosend(nin, nout, nscr, &mut scr.all());
                        continue;
                    }
                    if (c.math == 1149 && c.mfh == 13 && c.mth >= 51 && c.iverf < 5)
                        || (c.math == 1150 && c.mfh == 13 && c.mth <= 54 && c.iverf < 5)
                    {
                        // 140
                        c.tosend(nin, 0, 0, &mut scr.all());
                        let strng = format!("gamma prod patch made for mt {}", ffmt::i(c.mth, 3));
                        c.mess("hconvr", &strng, " ");
                        continue;
                    }
                    if c.mfh == 14 && l2flg == 1 {
                        lbl = 500;
                        continue;
                    }
                    c.contio(0, nout, nscr, &mut scr.all(), &mut nb, &mut nw);
                    if c.math == 0 {
                        lbl = 600;
                        continue;
                    }
                    c.tosend(nin, nout, nscr, &mut scr.all());
                }
                150 => {
                    //--allocate storage
                    if istor == 0 {
                        istor = 1;
                        let imax2 = imax * imax;
                        ev = RV::new((imax as usize) + SLACK);
                        eg = RV::new(lmax as usize + SLACK);
                        es = RV::new(lmax as usize + SLACK);
                        y = RV::new(lmax as usize + SLACK);
                        aa = RV::new(imax2 as usize + SLACK);
                        r = RV::new(imax2 as usize + SLACK);
                        ev[1] = 0.0;
                        for i in 1..=imax2 {
                            r[i] = 0.0;
                            aa[i] = 0.0;
                        }
                        for i in 1..=imax {
                            r[imax * (i - 1) + i] = 1.0;
                        }
                    }

                    //--convert transition probability array to yields
                    self.za = c.c1h;
                    self.awr = c.c2h;
                    let lg = c.l2h;
                    l2flg = 1;
                    let mut l = 1;
                    c.listio(nin, 0, 0, &mut scr.at(l), &mut nb, &mut nw);
                    l += nw;
                    while nb != 0 {
                        c.moreio(nin, 0, 0, &mut scr.at(l), &mut nb, &mut nw);
                        l += nw;
                    }

                    //--set base value for mt0
                    let mth = c.mth;
                    if (51..=91).contains(&mth) && mt0 != 49 {
                        mt0 = 49;
                    }
                    if c.iverf >= 6 {
                        if (601..=649).contains(&mth) && mt0 != 599 {
                            mt0 = 599;
                        }
                        if (651..=699).contains(&mth) && mt0 != 649 {
                            mt0 = 649;
                        }
                        if (701..=749).contains(&mth) && mt0 != 699 {
                            mt0 = 699;
                        }
                        if (751..=799).contains(&mth) && mt0 != 749 {
                            mt0 = 749;
                        }
                        if (801..=849).contains(&mth) && mt0 != 799 {
                            mt0 = 799;
                        }
                        if (876..=891).contains(&mth) && mt0 != 874 {
                            mt0 = 874;
                        }
                    } else if c.iverf <= 5 {
                        if (701..=719).contains(&mth) && mt0 != 699 {
                            mt0 = 699;
                        }
                        if (721..=739).contains(&mth) && mt0 != 719 {
                            mt0 = 719;
                        }
                        if (741..=759).contains(&mth) && mt0 != 739 {
                            mt0 = 739;
                        }
                        if (761..=779).contains(&mth) && mt0 != 759 {
                            mt0 = 759;
                        }
                        if (781..=799).contains(&mth) && mt0 != 779 {
                            mt0 = 779;
                        }
                    }

                    //--load the e(_) array with mf3 -q (overwritten with
                    //  mf12 data when possible).
                    if mt0 != mt0old {
                        for k in 1..=imax {
                            ev[k] = 0.0;
                        }
                        mt0old = mt0;
                        let m1 = mt0 + 2;
                        // m2 is set on every path the data can take
                        let mut m2 = 0;
                        if mt0 == 49 {
                            m2 = 91;
                        } else if c.iverf >= 6 && mt0 != 874 {
                            m2 = m1 + 48;
                        } else if c.iverf >= 6 && mt0 == 874 {
                            m2 = m1 + 15;
                        } else if c.iverf <= 5 {
                            m2 = m1 + 17;
                        }
                        let mut nn = 1usize;
                        let mut mttst = -1;
                        while mttst < m2 && nn <= nmf3 {
                            mttst = mtq[nn];
                            if mttst >= m1 && mttst <= m2 {
                                ev[mttst - mt0] = -eeq[nn];
                            }
                            nn += 1;
                        }
                    }
                    let j = mth - mt0;
                    ev[j] = scr[1];
                    let n = c.n2h;
                    let jm1 = j - 1;
                    for kk in 1..=jm1 {
                        let k = jm1 - kk + 1;
                        let mut ii = 0;
                        let mut idone = 0;
                        let mut i = 0;
                        while ii < n && idone == 0 {
                            ii += 1;
                            i = ii;
                            let ei = scr[6 + (lg + 1) * i - lg];
                            if ei == zero && ev[k] == zero {
                                idone = 1;
                            } else if ei != zero {
                                if (ei - ev[k]).abs() / ei <= 0.0001 {
                                    idone = 1;
                                }
                            }
                        }
                        if idone == 0 {
                            aa[(k - 1) * imax + j] = 0.0;
                            r[(k - 1) * imax + j] = 0.0;
                        } else {
                            let p = scr[7 + (lg + 1) * i - lg];
                            let mut g = 1.0f64;
                            if lg == 2 {
                                g = scr[6 + (lg + 1) * i];
                            }
                            aa[(k - 1) * imax + j] = p * g;
                            r[(k - 1) * imax + j] = p;
                        }
                        if k != jm1 {
                            let kp1 = k + 1;
                            for ii in kp1..=jm1 {
                                r[(k - 1) * imax + j] = r[(k - 1) * imax + j] + r[(ii - 1) * imax + j] * aa[(k - 1) * imax + ii];
                            }
                        }
                    }
                    let mut l = 0;
                    let mut ysum = 0.0f64;
                    for i in 2..=j {
                        let ja = j + 2 - i;
                        let eja = ev[ja];
                        for ii in 1..=jm1 {
                            let jb = j - ii;
                            let ejb = ev[jb];
                            let yy = aa[(jb - 1) * imax + ja] * r[(ja - 1) * imax + j];
                            if yy > zero {
                                l += 1;
                                if l > lmax {
                                    c.error("hconvr", "too many lo=2 gammas", " ");
                                }
                                eg[l] = eja - ejb;
                                es[l] = eja;
                                y[l] = yy;
                                ysum += y[l];
                            }
                        }
                    }
                    ev[1] = 0.0;

                    //--arrange gamma energies in descending order
                    if l > 1 {
                        let lm1 = l - 1;
                        for i in 1..=lm1 {
                            let ip1 = i + 1;
                            for ii in ip1..=l {
                                if eg[i] < eg[ii] {
                                    let mut tsave = eg[i];
                                    eg[i] = eg[ii];
                                    eg[ii] = tsave;
                                    tsave = y[i];
                                    y[i] = y[ii];
                                    y[ii] = tsave;
                                    tsave = es[i];
                                    es[i] = es[ii];
                                    es[ii] = tsave;
                                }
                            }
                        }
                    }

                    //--output computed yields in endf lo=1 format
                    scr[1] = self.za;
                    scr[2] = self.awr;
                    scr[3] = 1.0;
                    scr[4] = 0.0;
                    scr[5] = l as f64;
                    scr[6] = 0.0;
                    let mut mtl = 49;
                    if c.iverf <= 5 && mth >= 700 {
                        mtl = 649;
                    }
                    if c.iverf >= 6 && mth >= 600 {
                        mtl = 549;
                    }
                    ngam[(mth - mtl) as usize] = l;
                    c.contio(0, nout, nscr, &mut scr.all(), &mut nb, &mut nw);
                    // output tab1 sum record
                    if l > 1 {
                        scr[1] = 0.0;
                        scr[2] = 0.0;
                        scr[3] = 0.0;
                        scr[5] = 1.0;
                        scr[6] = 2.0;
                        scr[7] = 2.0;
                        scr[8] = 2.0;
                        scr[9] = self.ebot;
                        scr[10] = ysum;
                        scr[11] = self.etop;
                        scr[12] = ysum;
                        nw = 12;
                        c.tab1io(0, nout, nscr, &mut scr.all(), &mut nb, &mut nw);
                    }
                    // output rest of tab1 records
                    for i in 1..=l {
                        scr[1] = eg[i];
                        scr[2] = es[i];
                        scr[3] = 0.0;
                        scr[4] = 2.0;
                        scr[5] = 1.0;
                        scr[6] = 2.0;
                        scr[7] = 2.0;
                        scr[8] = 2.0;
                        scr[9] = self.ebot;
                        scr[10] = y[i];
                        scr[11] = self.etop;
                        scr[12] = y[i];
                        nw = 12;
                        c.tab1io(0, nout, nscr, &mut scr.all(), &mut nb, &mut nw);
                    }

                    //--copy send record and loop over remaining sections.
                    c.contio(nin, nout, nscr, &mut scr.all(), &mut nb, &mut nw);
                    lbl = 110;
                }
                500 => {
                    //--angular distributions
                    if (51..=90).contains(&c.mth) {
                        scr[3] = 1.0;
                        scr[4] = 0.0;
                        scr[6] = 0.0;
                        scr[5] = ngam[(c.mth - 50) as usize] as f64;
                        c.contio(0, nout, nscr, &mut scr.all(), &mut nb, &mut nw);
                        c.tosend(nin, 0, 0, &mut scr.all());
                        c.asend(nout, nscr);
                    } else {
                        c.contio(0, nout, nscr, &mut scr.all(), &mut nb, &mut nw);
                        c.tosend(nin, nout, nscr, &mut scr.all());
                    }
                    lbl = 110;
                }
                550 => {
                    //--file1.  add 456=452 if necessary.
                    c.contio(0, nout, nscr, &mut scr.all(), &mut nb, &mut nw);
                    c.tosend(nin, nout, nscr, &mut scr.all());
                    let mut no455 = 0;
                    c.contio(nin, nout, nscr, &mut scr.all(), &mut nb, &mut nw);
                    if c.mfh == 0 {
                        lbl = 110;
                        continue;
                    }
                    if c.mth != 452 {
                        // 595
                        c.tofend(nin, nout, nscr, &mut scr.all());
                        lbl = 110;
                        continue;
                    }
                    let mut tmpnu = rv(c.npage + 46);
                    let mut l = 1;
                    lnu = c.l2h;
                    if lnu == 1 {
                        c.listio(nin, nout, nscr, &mut tmpnu.at(l), &mut nb, &mut nw);
                    } else if lnu == 2 {
                        c.tab1io(nin, nout, nscr, &mut tmpnu.at(l), &mut nb, &mut nw);
                    }
                    nnu = nb + nw + 1;
                    nu = rv(nnu);
                    for k in 1..=nw {
                        nu[k] = tmpnu[k];
                    }
                    while nb != 0 {
                        l += nw;
                        if l > nnu {
                            c.error("hconvr", "exceeded storage for nubar", " ");
                        }
                        c.moreio(nin, nout, nscr, &mut nu.at(l), &mut nb, &mut nw);
                    }
                    c.tosend(nin, nout, nscr, &mut scr.all());
                    // go570: the loop ended at label 570 (else at 590)
                    let mut go570 = false;
                    loop {
                        // 565
                        c.contio(nin, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                        if c.mth == 456 {
                            break;
                        }
                        if c.mth == 0 && no455 == 0 {
                            go570 = true;
                            break;
                        }
                        if c.mfh == 0 {
                            break;
                        }
                        if c.mth == 455 {
                            no455 = 1;
                        }
                        if c.mth > 456 && no455 == 0 {
                            go570 = true;
                            break;
                        }
                        c.contio(0, nout, nscr, &mut scr.all(), &mut nb, &mut nw);
                        c.tosend(nin, nout, nscr, &mut scr.all());
                    }
                    if go570 && no455 != 1 {
                        // 570: add mt456 as a copy of mt452
                        let mf = c.mfh;
                        let mt = c.mth;
                        let l1 = c.l1h;
                        let l2 = c.l2h;
                        let n1 = c.n1h;
                        let n2 = c.n2h;
                        c.mfh = 1;
                        c.mth = 456;
                        scr[3] = 0.0;
                        scr[4] = lnu as f64;
                        scr[5] = 0.0;
                        scr[6] = 0.0;
                        let mut l = 1;
                        c.contio(0, nout, nscr, &mut scr.all(), &mut nb, &mut nw);
                        if lnu == 1 {
                            c.listio(0, nout, nscr, &mut nu.at(l), &mut nb, &mut nw);
                        }
                        if lnu == 2 {
                            c.tab1io(0, nout, nscr, &mut nu.at(l), &mut nb, &mut nw);
                        }
                        while nb != 0 {
                            l += nw;
                            if l > nnu {
                                c.error("hconvr", "exceeded storage for nubar", " ");
                            }
                            c.moreio(0, nout, nscr, &mut nu.at(l), &mut nb, &mut nw);
                        }
                        c.asend(nout, nscr);
                        c.mfh = mf;
                        c.mth = mt;
                        scr[3] = l1 as f64;
                        scr[4] = l2 as f64;
                        scr[5] = n1 as f64;
                        scr[6] = n2 as f64;
                    }
                    // 590
                    c.contio(0, nout, nscr, &mut scr.all(), &mut nb, &mut nw);
                    if c.mfh == 0 {
                        lbl = 110;
                        continue;
                    }
                    // 595
                    c.tofend(nin, nout, nscr, &mut scr.all());
                    lbl = 110;
                }
                _ => {
                    //--conver is finished.
                    let _ = (istor, nnu);
                    return;
                }
            }
        }
    }
}

impl Heatr {
    // ================================================================ gheat
    /// Subtract the energy carried away by photons. iold and inew are
    /// exchanged as in the Fortran.
    fn gheat(&mut self, c: &mut Ctx, iold: &mut i32, inew: &mut i32, nscr: i32) {
        let small = 1.0e-10f64;
        let big = 1.0e10f64;
        let smin = 1.0e-9f64;
        let qtest = 99.0e6f64;
        let zero = 0.0f64;
        let (mut nb, mut nw) = (0i32, 0i32);
        let matd = self.matd;
        let nin = self.nin;
        let nendf = self.nendf;
        let mut cv = RV::new(30 + SLACK);
        let mut imt = [0i32; 31];

        //--allocate buffers for loada/finda
        let mut bufo = vec![0.0f64; NBUF as usize];
        let mut bufn = vec![0.0f64; NBUF as usize];

        //--find first section of mf12 or mf13
        let mut scr = rv(c.npage + 50);
        let mut mfd = 3;
        let afact = self.awr / (self.awr + 1.0);
        let aw1fac = 1.0 / (self.awr + 1.0);
        let z = nint(self.za / 1000.0) as f64;
        let mut e = 0.0f64;
        let awr = self.awr;
        let mut dame = self.df(e, z, awr, z, awr);
        // gheat's own mgam (a local that hides the module variable)
        let mut mgam = 2;
        let mut hk = 0.0f64;
        loop {
            // 100
            c.contio(nscr, 0, 0, &mut scr.all(), &mut nb, &mut nw);
            if c.mfh == 12 {
                mgam = 1;
            }
            if c.mfh >= 12 {
                break;
            }
            c.tofend(nscr, 0, 0, &mut scr.all());
        }
        // 102
        let mfc = c.mfh;
        let mtc = c.mth;
        c.findf(matd, mfc, mtc, nscr);
        let mut npkk = self.npk;
        if self.kchk == 1 {
            npkk = 3 * self.npk - 2;
        }
        if mgam > 0 {
            npkk += 3;
        }
        let npkt = self.npk - 1;

        //--sum all reactions found in mf12 and mf13.
        //--reconstruct redundant cross sections if needed.
        let mut b = rv(c.npage + 50);
        let nd: i32 = 10000;
        let mut d = rv(nd);

        // locals kept across the labels
        let mut mfx = 0;
        let mut mtd = 0;
        let mut mtx = 0;
        let mut nk = 0;
        let mut ik = 0;
        let mut cerr = 0.0f64;
        let mut elow = 0.0f64;
        let mut ehigh = 0.0f64;
        let mut thresh = 0.0f64;
        let mut jdis = 0i32;
        let mut y = 0.0f64;
        let mut ilist = 1i32;
        let mut egkr = 0.0f64;
        let mut nwd = nd;
        let mut ebar = 0.0f64;
        // egam, edam: set by disgam or gambar before they are used
        let mut egam = 0.0f64;
        let mut edam = 0.0f64;
        let mut lp = 0;
        let mut lqx = 0i32;
        let mut enxt = 0.0f64;
        let mut idis = 0i32;
        let mut x = 0.0f64;
        let mut damn = 0.0f64;
        let mut ipx = 2;
        let mut irx = 1;
        let mut nmt = 0i32;
        let mut h = 0.0f64;
        let mut enext = 0.0f64;
        let mut enx = 0.0f64;
        let mut idx = 0i32;

        let mut lbl = if mgam == 2 { 210 } else { 105 };
        loop {
            match lbl {
                105 => {
                    c.contio(nscr, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                    if c.mfh > 13 {
                        lbl = 250;
                        continue;
                    }
                    if c.mfh == 0 {
                        lbl = 220;
                        continue;
                    }
                    if c.mth == 0 {
                        continue;
                    }
                    if c.mfh == 12 && c.mth == 460 {
                        c.tosend(nscr, 0, 0, &mut scr.all());
                        continue;
                    }
                    // skip over this mf/mt if photon data were already found in mf6
                    if self.i6p > 0 {
                        let mut mf6flg = 0;
                        for i in 1..=self.i6p as usize {
                            if self.mt6yp[i] != c.mth {
                                continue;
                            }
                            mf6flg = 1;
                            break;
                        }
                        if mf6flg != 0 {
                            let strng1 = format!("skipping mf{}/mt = {}", ffmt::i(c.mfh, 2), ffmt::i(c.mth, 3));
                            let strng2 = "photons were already processed in mf6";
                            c.mess("gheat ", &strng1, strng2);
                            c.tosend(nscr, 0, 0, &mut scr.all());
                            continue;
                        }
                    }
                    let lo = c.l1h;
                    if lo == 2 {
                        c.error("gheat", "lo=2 not coded.", " ");
                    }
                    mfx = c.mfh;
                    mtd = c.mth;
                    mtx = c.mth;
                    if c.mth == 3 {
                        mtd = 1;
                    }
                    if self.iprint == 1 && c.mfh == 12 && self.idame > 0 && mtd == 102 {
                        c.io.so(&format!(
                            "\n photon energy (from yields) mf{}, mt{}\n{}e{}ebar/err{}egam{}edam{}xsec{}yield{}heating{}damage",
                            ffmt::i(c.mfh, 2),
                            ffmt::i(c.mth, 3),
                            " ".repeat(14),
                            " ".repeat(6),
                            " ".repeat(10),
                            " ".repeat(10),
                            " ".repeat(10),
                            " ".repeat(9),
                            " ".repeat(7),
                            " ".repeat(8)
                        ));
                    }
                    if self.iprint == 1 && c.mfh == 12 && (self.idame == 0 || mtd != 102) {
                        c.io.so(&format!(
                            "\n photon energy (from yields) mf{}, mt{}\n{}e{}ebar{}xsec{}yield{}energy{}heating",
                            ffmt::i(c.mfh, 2),
                            ffmt::i(c.mth, 3),
                            " ".repeat(14),
                            " ".repeat(10),
                            " ".repeat(10),
                            " ".repeat(9),
                            " ".repeat(8),
                            " ".repeat(7)
                        ));
                    }
                    if self.iprint == 1 && c.mfh == 13 {
                        c.io.so(&format!(
                            "\n photon energy (from xsecs) mf{}, mt{}\n{}e{}ebar{}xsec{}energy{}heating",
                            ffmt::i(c.mfh, 2),
                            ffmt::i(c.mth, 3),
                            " ".repeat(14),
                            " ".repeat(10),
                            " ".repeat(10),
                            " ".repeat(8),
                            " ".repeat(7)
                        ));
                    }
                    lbl = 110;
                }
                110 => {
                    nk = c.n1h;
                    ik = 0;
                    cerr = 0.0;
                    // check for leading or trailing zeroes to locate ranges
                    // where this partial is not defined
                    e = 0.0;
                    c.gety1(e, &mut enxt, &mut jdis, &mut y, nscr, &mut scr.all());
                    e = enxt;
                    let mut el = sigfig(e, 7, 0);
                    elow = big;
                    ehigh = 0.0;
                    let test = big - big / 100.0;
                    while e < test {
                        c.gety1(e, &mut enxt, &mut jdis, &mut y, nscr, &mut scr.all());
                        if y > zero && e < elow {
                            elow = el;
                        }
                        if y > zero {
                            ehigh = 0.0;
                        }
                        if ehigh == zero && y == zero {
                            ehigh = sigfig(e, 7, 0);
                        }
                        el = sigfig(e, 7, 0);
                        e = enxt;
                    }
                    if ehigh == zero {
                        ehigh = sigfig(e, 7, 0);
                    }
                    if nk == 1 {
                        c.findf(matd, mfx, mtx, nscr);
                        c.contio(nscr, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                    }
                    lbl = 130;
                }
                130 => {
                    e = 0.0;
                    c.gety1(e, &mut thresh, &mut jdis, &mut y, nscr, &mut scr.all());
                    ilist = 1;
                    egkr = c.c1h;
                    nwd = nd;
                    if egkr == zero {
                        self.gambar(c, e, &mut ebar, &mut egam, &mut edam, nendf, matd, mtx, &mut d, nwd);
                    }
                    lp = c.l1h;
                    ik += 1;
                    if self.iprint == 1 && egkr == zero {
                        c.io.so(&format!("{}  continuum gammas", ffmt::i(ik, 6)));
                    }
                    if self.iprint == 1 && egkr > zero {
                        c.io.so(&format!("{}{} ev gamma", ffmt::i(ik, 6), ffmt::e1p(egkr, 13, 4)));
                    }
                    lqx = 0;
                    if mfd != 13 {
                        c.findf(matd, mfd, mtd, nin);
                        c.contio(nin, 0, 0, &mut b.all(), &mut nb, &mut nw);
                        c.gety2(e, &mut enxt, &mut idis, &mut x, nin, &mut b.all());
                        if c.mth == 102 {
                            self.q = c.c2h;
                            if self.nqa > 0 {
                                for i in 1..=self.nqa as usize {
                                    if self.mta[i] == c.mth {
                                        self.q = self.qa[i];
                                        if self.q >= qtest {
                                            lqx = self.lqs[i];
                                        }
                                    }
                                }
                            }
                            let (q, za, awr, mth) = (self.q, self.za, self.awr, c.mth);
                            self.capdam(c, e, &mut damn, q, za, awr, mth);
                            disgam(self, e, &mut egam, &mut edam, z, awr);
                            ipx = 2;
                            irx = 1;
                        }
                        if enxt > thresh {
                            thresh = enxt;
                        }
                    }
                    // 150
                    if ik <= 1 {
                        nmt = 0;
                        if self.npk >= 3 {
                            let mut mty = mtx;
                            if mty == 3 {
                                mty = 1;
                            }
                            if mty == 4 {
                                mty = 51;
                            }
                            if mty == 18 && self.mt19 == 1 {
                                mty = 19;
                            }
                            indx(self.npk, &self.mtp, self.mt19, mty, &mut imt, &mut nmt, c.iverf);
                        }
                    }
                    // 160
                    let mut i = 0;
                    loop {
                        // 190
                        i += 1;
                        c.finda(i, &mut cv.0, npkk, *iold, &mut bufo, NBUF);
                        e = sigfig(cv[1], 9, 0);
                        let skip = e < elow * (1.0 - small) || e > ehigh * (1.0 + small);
                        if !skip {
                            c.gety1(e, &mut enext, &mut jdis, &mut y, nscr, &mut scr.all());
                            x = zero;
                            let mut zero_x = false;
                            if mfd == 3 {
                                c.gety2(e, &mut enext, &mut idis, &mut x, nin, &mut b.all());
                                if x == zero {
                                    zero_x = true;
                                }
                            }
                            if !zero_x {
                                if lqx != 0 {
                                    c.terpa(&mut self.q, e, &mut enx, &mut idx, &self.qbar.at(lqx), &mut ipx, &mut irx);
                                }
                                let mut egk = egkr;
                                if lp == 2 {
                                    egk = egkr + e * afact;
                                }
                                ebar = egk;
                                dame = 0.0;
                                if ik == 1 {
                                    cv[npkk] = 0.0;
                                }
                                if egkr == zero {
                                    self.gambar(c, e, &mut ebar, &mut egam, &mut edam, nendf, matd, mtx, &mut d, nwd);
                                }
                                let el = self.elist[ilist as usize];
                                if mfd == 13 {
                                    // 170: mf13 reactions
                                    h = -y * ebar;
                                    cv[npkk - 1] -= h;
                                    cv[npkk] += h;
                                    let hx = h;
                                    cv[2] += h;
                                    if nmt > 0 {
                                        for ii in 1..=nmt as usize {
                                            let indxx = imt[ii];
                                            if self.mtp[indxx as usize] < 442 {
                                                cv[indxx] += h;
                                            }
                                            if self.mtp[indxx as usize] == 442 {
                                                cv[indxx] -= h;
                                            }
                                        }
                                    }
                                    let mut prt = !(self.iprint == 0 || (e - el).abs() > small * e);
                                    if prt && ik < nk && y < smin {
                                        prt = false;
                                    }
                                    if prt && nk == 1 && y < smin {
                                        prt = false;
                                    }
                                    let mut test = -1.0f64;
                                    test /= 10000.0;
                                    if prt && nk > 1 && ik == nk && cv[npkk] > test {
                                        prt = false;
                                    }
                                    if prt {
                                        let ebarp = sigfig(ebar, 9, 0);
                                        let yp = sigfig(y, 9, 0);
                                        let hxp = sigfig(hx, 9, 0);
                                        let hp = sigfig(h, 9, 0);
                                        c.io.so(&e14_lines(&[e, ebarp, yp, -hxp, hp], 8));
                                        if ik == nk && nk > 1 {
                                            let mut subtot = -cv[npkk];
                                            subtot = sigfig(subtot, 9, 0);
                                            c.io.so(&format!("{}subtot   {}", " ".repeat(48), ffmt::e1p(subtot, 14, 4)));
                                        }
                                    }
                                } else {
                                    let mut go180 = false;
                                    if c.mth == 102 {
                                        // photon recoil correction
                                        if egkr != zero {
                                            disgam(self, egkr, &mut egam, &mut edam, z, awr);
                                        }
                                        h = egam * x * y;
                                        hk = h;
                                        cv[npkk - 1] += x * y * ebar;
                                        cv[npkk] += y * ebar;
                                        dame = edam * x * y;
                                        if ik == nk {
                                            if self.idame > 0 {
                                                let (q, za, mth) = (self.q, self.za, c.mth);
                                                self.capdam(c, e, &mut damn, q, za, awr, mth);
                                                dame -= damn * x;
                                            }
                                            hk -= x * self.rtm * powi(self.q + e * awr * aw1fac, 2) / 2.0;
                                            let cfix = (e + self.q - e * aw1fac) * x;
                                            h -= cfix;
                                            cv[npkk - 2] = cfix - x * cv[npkk];
                                            let eava = afact * e + self.q;
                                            cerr = 100.0 * (cv[npkk] - eava) / eava;
                                        }
                                        let ebarp = sigfig(ebar, 9, 0);
                                        let xp = sigfig(x, 9, 0);
                                        let yp = sigfig(y, 9, 0);
                                        let hp = sigfig(h, 9, 0);
                                        let egamp = sigfig(egam, 9, 0);
                                        let edamp = sigfig(edam, 9, 0);
                                        let damep = sigfig(dame, 9, 0);
                                        cerr = sigfig(cerr, 4, 0);
                                        if !(self.iprint != 1 || (e - el).abs() > small * e || y <= smin) {
                                            if self.idame == 0 {
                                                c.io.so(&e14_lines(&[e, ebarp, xp, yp, hp], 8));
                                            } else {
                                                c.io.so(&e14_lines(&[e, ebarp, egamp, edamp, xp, yp, hp, damep], 8));
                                            }
                                            if ik == nk {
                                                c.io.so(&format!(" {}{} pc", ffmt::e1p(e, 14, 4), ffmt::f(cerr, 10, 1)));
                                            }
                                        }
                                    } else {
                                        // 164: other mf12 reactions
                                        if mtd != 2 && mtd != 460 {
                                            h = -y * x * ebar;
                                        }
                                        if mtd == 2 {
                                            h = y * x * ebar;
                                        }
                                        cv[npkk - 1] -= h;
                                        cv[npkk] += h;
                                    }
                                    // 166
                                    let hx = h;
                                    cv[2] += h;
                                    if nmt > 0 {
                                        for ii in 1..=nmt as usize {
                                            let indxx = imt[ii];
                                            let m = self.mtp[indxx as usize];
                                            if m >= 444 {
                                                cv[indxx] += dame;
                                            }
                                            if m < 442 {
                                                cv[indxx] += h;
                                            }
                                            if m == 442 {
                                                cv[indxx] -= h;
                                            }
                                            if m == 443 && c.mth == 102 {
                                                cv[indxx] += hk;
                                            }
                                        }
                                    }
                                    if self.iprint == 0 || e != el {
                                        go180 = true;
                                    }
                                    if !go180 && c.mth == 102 {
                                        go180 = true;
                                    }
                                    if !go180 && ik < nk && y < smin {
                                        go180 = true;
                                    }
                                    if !go180 && nk == 1 && y < smin {
                                        go180 = true;
                                    }
                                    let mut test = -1.0f64;
                                    test /= 10000.0;
                                    if !go180 && nk > 1 && ik == nk && cv[npkk] > test {
                                        go180 = true;
                                    }
                                    if !go180 {
                                        let ebarp = sigfig(ebar, 9, 0);
                                        let xp = sigfig(x, 9, 0);
                                        let yp = sigfig(y, 9, 0);
                                        let hxp = sigfig(hx, 9, 0);
                                        let hp = sigfig(h, 9, 0);
                                        c.io.so(&e14_lines(&[e, ebarp, xp, yp, -hxp, hp], 8));
                                        if ik == nk && nk > 1 {
                                            let mut subtot = -cv[npkk];
                                            subtot = sigfig(subtot, 9, 0);
                                            c.io.so(&format!("{}subtot   {}", " ".repeat(48), ffmt::e1p(subtot, 14, 4)));
                                        }
                                    }
                                }
                            }
                        }
                        // 180, 195
                        let mut it = i;
                        if i == self.ne {
                            it = -it;
                        }
                        c.loada(it, &cv.0, npkk, *inew, &mut bufn, NBUF);
                        if (e - self.elist[ilist as usize]).abs() < small * e {
                            ilist += 1;
                        }
                        if i < self.ne {
                            continue;
                        }
                        break;
                    }

                    //--loop over all subsections and reactions.
                    std::mem::swap(iold, inew);
                    if ik < nk {
                        lbl = 130;
                        continue;
                    }
                    if mfd == 13 {
                        lbl = 105;
                        continue;
                    }
                    if mtx == 3 && mtd == 1 {
                        // 200
                        mtd = 2;
                        c.findf(matd, mfx, mtx, nscr);
                        c.contio(nscr, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                        lbl = 110;
                        continue;
                    }
                    lbl = 105;
                }
                210 => {
                    c.mess("gheat", "no file 12 for this material.", " ");
                    c.findf(matd, 13, 0, nscr);
                    lbl = 220;
                }
                220 => {
                    //--finished with present file.
                    if mfd == 13 {
                        lbl = 250;
                        continue;
                    }
                    mfd = 13;
                    mtd = 0;
                    lbl = 105;
                }
                _ => {
                    // 250
                    if self.kchk == 1 && self.mt303 != 0 {
                        //--print kinematic test for total photon production
                        ilist = 1;
                        c.io.so(&format!(
                            "\n photon energy production check\n{}e{}ev-barns{}min{}max",
                            " ".repeat(14),
                            " ".repeat(6),
                            " ".repeat(11),
                            " ".repeat(11)
                        ));
                        let mt303 = self.mt303;
                        for i in 1..=self.ne {
                            c.finda(i, &mut cv.0, npkk, *iold, &mut bufo, NBUF);
                            e = sigfig(cv[1], 9, 0);
                            if (e - self.elist[ilist as usize]).abs() <= small * e {
                                ilist += 1;
                                let elo = cv[npkk - 1] + cv[mt303] - cv[mt303 + 2 * npkt] + cv[npkk - 2];
                                let ehi = cv[npkk - 1] + cv[mt303] - cv[mt303 + npkt] + cv[npkk - 2];
                                let mut jerr = 0;
                                if cv[npkk - 1] > ehi + ehi / 10.0 {
                                    jerr = 1;
                                }
                                if cv[npkk - 1] < elo - elo / 10.0 {
                                    jerr = -1;
                                }
                                let line = e14_lines(&[e, cv[npkk - 1], elo, ehi], 8);
                                if jerr == 0 {
                                    c.io.so(&line);
                                }
                                if jerr > 0 {
                                    c.io.so(&format!("{}   ++++", line));
                                }
                                if jerr < 0 {
                                    c.io.so(&format!("{}   ----", line));
                                }
                            }
                        }
                    }
                    // 300
                    let _ = (dame, mgam);
                    return;
                }
            }
        }
    }

    // =============================================================== gambar
    /// Mean energy for continuous photon spectra (File 15) and the photon
    /// recoil correction for capture (e=0 to initialize).
    fn gambar(&mut self, c: &mut Ctx, e: f64, ebar: &mut f64, esqb: &mut f64, esqd: &mut f64, nin: i32, matd: i32, mtd: i32, a: &mut RV, nwamax: i32) {
        let small = 1.0e-10f64;
        let up = 1.000001e0f64;
        let zero = 0.0f64;
        let (mut nb, mut nw) = (0i32, 0i32);

        //--initialize if e=0.
        if e <= zero {
            let mfd = 15;
            c.findf(matd, mfd, mtd, nin);
            let mut l = 1;
            c.contio(nin, 0, 0, &mut a.at(l), &mut nb, &mut nw);
            let iz = nint(c.c1h / 1000.0);
            self.gb.z = iz as f64;
            self.gb.awr = c.c2h;
            c.tab1io(nin, 0, 0, &mut a.at(l), &mut nb, &mut nw);
            l += nw;
            while nb != 0 {
                if l > nwamax {
                    c.error("gambar", "storage exceeded in a.", " ");
                }
                c.moreio(nin, 0, 0, &mut a.at(l), &mut nb, &mut nw);
                l += nw;
            }
            self.gb.lf = 1;

            //--read in first incident energy point.
            l = 1;
            c.tab2io(nin, 0, 0, &mut a.at(l), &mut nb, &mut nw);
            let s = &mut self.gb;
            s.ne = nint(a[l + 5]);
            s.nnt = l + 6;
            s.nbt = nint(a[s.nnt]);
            s.inn = nint(a[s.nnt + 1]);
            s.inn = 2;
            s.istart = l + nw;
            l = s.istart;
            c.tab1io(nin, 0, 0, &mut a.at(l), &mut nb, &mut nw);
            l += nw;
            while nb != 0 {
                if l > nwamax {
                    c.error("gambar", "storage exceeded in a.", " ");
                }
                c.moreio(nin, 0, 0, &mut a.at(l), &mut nb, &mut nw);
                l += nw;
            }
            let s = &mut self.gb;
            let iend = l - 1;
            s.jstart = iend + 1;
            s.jend = 0;
            s.ehi = a[s.istart + 1];
            s.elo = s.ehi;
            s.nne = 1;
            s.ilo = 0;
            s.ihi = 0;
            s.fhi = 0.0;
            s.ghi = 0.0;
            s.hhi = 0.0;
            return;
        }

        //--normal entry.  is desired energy in current panel.
        let mut lbl = if e < self.gb.ehi * (1.0 - small) {
            140
        } else if self.gb.nne == 1 {
            130
        } else if self.gb.nne == self.gb.ne {
            160
        } else {
            120
        };
        loop {
            match lbl {
                120 => {
                    //--no.  slide high energy data into low energy positions.
                    let s = &mut self.gb;
                    let mut j = -1;
                    for il in s.jstart..=s.jend {
                        j += 1;
                        a[j + s.istart] = a[il];
                    }
                    let iend = s.istart + j;
                    s.jstart = iend + 1;
                    s.elo = s.ehi;
                    s.flo = s.fhi;
                    if mtd == 102 {
                        s.glo = s.ghi;
                    }
                    if mtd == 102 {
                        s.hlo = s.hhi;
                    }
                    s.ilo = s.ihi;
                    lbl = 130;
                }
                130 => {
                    //--read new high energy data record.
                    if self.gb.nne == self.gb.ne {
                        c.error("gambar", "requested energy gt highest given.", " ");
                    }
                    let mut l = self.gb.jstart;
                    c.tab1io(nin, 0, 0, &mut a.at(l), &mut nb, &mut nw);
                    l += nw;
                    while nb != 0 {
                        if l > nwamax {
                            c.error("gambar", "storage exceeded in a.", " ");
                        }
                        c.moreio(nin, 0, 0, &mut a.at(l), &mut nb, &mut nw);
                        l += nw;
                    }
                    let s = &mut self.gb;
                    s.jend = l - 1;
                    s.nne += 1;
                    s.ihi = 0;
                    s.ehi = a[s.jstart + 1];
                    if s.ehi < e * (1.0 - small) && s.nne < s.ne {
                        lbl = 120;
                        continue;
                    }

                    //--integrate at high and low energies.
                    let (istart, jstart, lf) = (s.istart, s.jstart, s.lf);
                    let (z, awr) = (s.z, s.awr);
                    if self.gb.ilo <= 0 {
                        let mut flo = 0.0f64;
                        tabbar(c, &mut flo, &mut a.at(istart), lf);
                        self.gb.flo = flo;
                        if mtd == 102 {
                            let (mut glo, mut hlo) = (0.0f64, 0.0f64);
                            self.tabsqr(c, &mut glo, &mut hlo, &mut a.at(istart), lf, z, awr);
                            self.gb.glo = glo;
                            self.gb.hlo = hlo;
                        }
                    }
                    let mut fhi = 0.0f64;
                    tabbar(c, &mut fhi, &mut a.at(jstart), lf);
                    self.gb.fhi = fhi;
                    if mtd == 102 {
                        let (mut ghi, mut hhi) = (0.0f64, 0.0f64);
                        self.tabsqr(c, &mut ghi, &mut hhi, &mut a.at(jstart), lf, z, awr);
                        self.gb.ghi = ghi;
                        self.gb.hhi = hhi;
                    }
                    self.gb.ihi = 1;
                    lbl = 140;
                }
                140 => {
                    //--yes.  interpolate for mean energy.
                    let s = &mut self.gb;
                    if s.nne == 1 {
                        //--return zero below threshold.
                        *ebar = 0.0;
                        *esqb = 0.0;
                        *esqd = 0.0;
                        return;
                    }
                    if s.nne > s.nbt {
                        s.nnt += 2;
                        s.nbt = nint(a[s.nnt]);
                        s.inn = nint(a[s.nnt + 1]);
                        s.inn = 2;
                    }
                    if s.elo == s.ehi {
                        s.elo = sigfig(s.elo, 7, -1);
                        s.ehi = sigfig(s.ehi, 7, 1);
                    }
                    let mut sen = 0.0f64;
                    c.terp1(s.elo, s.flo, s.ehi, s.fhi, e, &mut sen, s.inn);
                    if mtd == 102 {
                        c.terp1(s.elo, s.glo, s.ehi, s.ghi, e, esqb, s.inn);
                    }
                    if mtd == 102 {
                        c.terp1(s.elo, s.hlo, s.ehi, s.hhi, e, esqd, s.inn);
                    }
                    *ebar = sen;
                    return;
                }
                _ => {
                    //--return zero or fhi above table.
                    let s = &self.gb;
                    if s.inn != 1 && e >= up * s.ehi {
                        *ebar = 0.0;
                        *esqb = 0.0;
                        *esqd = 0.0;
                    } else {
                        *ebar = s.fhi;
                        *esqb = s.ghi;
                        *esqd = s.hhi;
                    }
                    return;
                }
            }
        }
    }

    // =============================================================== tabsqr
    /// Average photon recoil and damage energy for a tabulated section of
    /// File 15 (capture).
    fn tabsqr(&mut self, c: &Ctx, g: &mut f64, h: &mut f64, a: &mut F, _law: i32, z: f64, awr: f64) {
        let ein = 2.0 * self.tm;
        let rein = 1.0 / ein;
        *g = 0.0;
        *h = 0.0;
        let mut s = 0.0f64;
        let nr = nint(a[5]);
        let np = nint(a[6]);
        let ibase = 6 + 2 * nr;
        let mut ir = 1;
        let mut nbt = nint(a[7]);
        let mut inn = nint(a[8]);
        let mut xh = a[ibase + 1];
        let mut yh = a[ibase + 2];

        //--accumulate contributions to integrals
        for i in 2..=np {
            let xl = xh;
            xh = a[ibase + 2 * i - 1];
            let yl = yh;
            yh = a[ibase + 2 * i];
            if i > nbt {
                ir += 1;
                nbt = nint(a[6 + 2 * ir - 1]);
                inn = nint(a[6 + 2 * ir]);
            }
            if xl != xh {
                let dx = xh - xl;
                for j in 0..4 {
                    let x = xl + (1.0 + QP4[j]) * dx / 2.0;
                    let mut y = 0.0f64;
                    c.terp1(xl, yl, xh, yh, x, &mut y, inn);
                    let xr = x * x * rein;
                    *g += QW4[j] * xr * dx * y;
                    *h += QW4[j] * self.df(xr, z, awr + 1.0, z, awr) * dx * y;
                    s += QW4[j] * dx * y;
                }
            }
        }
        *g /= s;
        *h /= s;
    }
}

/// `disgam`: recoil and damage energy of discrete capture photons (e=0
/// does nothing).
fn disgam(ht: &mut Heatr, e: f64, egam: &mut f64, edam: &mut f64, z: f64, awr: f64) {
    if e == 0.0 {
        return;
    }
    *egam = e * e * ht.rtm / 2.0;
    *edam = ht.df(*egam, z, awr + 1.0, z, awr);
}

impl Heatr {
    // ================================================================= hout
    /// Write the output tape with the corrected dictionary and the kerma
    /// MT-s added; optionally write the plot file.
    fn hout(&mut self, c: &mut Ctx, iold: i32) {
        let small = 1.0e-10f64;
        let big = 1.0e10f64;
        let emax = 20.0e6f64;
        let ecut = 0.01e6f64;
        let zero = 0.0f64;
        let ten = 10.0f64;
        let nwscr: i32 = 10000;
        let (mut nb, mut nw) = (0i32, 0i32);
        let matd = self.matd;
        let nin = self.nin;
        let nout = self.nout;
        let npk = self.npk;
        let ne = self.ne;
        let mut cv = RV::new(30 + SLACK);
        let mut ncds = [0i32; 31];
        let mut klo = ["    "; 10];
        let mut khi = ["    "; 10];

        // allocate scratch storage for endf records and buffer
        let mut scr = rv(nwscr);
        let mut buf = vec![0.0f64; NBUF as usize];

        //--construct tab1 records for the kerma factors
        let mut n = 0;
        let mut nscr = 15;
        if nout < 0 {
            nscr = -nscr;
        }
        c.openz(nscr, 1);
        c.repoz(nscr);
        if self.iprint == 1 {
            // '(/'' final kerma factors''/14x,''e'',8(11x,i3))': format
            // reversion starts a new record with the group (11x,i3)
            let mut s = format!("\n final kerma factors\n{}e", " ".repeat(14));
            for (k, i) in (2..=npk as usize).enumerate() {
                if k > 0 && k % 8 == 0 {
                    s.push('\n');
                }
                s.push_str(&" ".repeat(11));
                s.push_str(&ffmt::i(self.mtp[i], 3));
            }
            c.io.so(&s);
        }
        let mut npkk = npk;
        if self.kchk == 1 {
            npkk = 3 * npk - 2;
        }
        if self.mgam > 0 {
            npkk += 3;
        }
        if self.mgam == 0 && self.i6 > 0 {
            npkk += 1;
        }
        let npkt = npk - 1;
        let mut npkd = npk;
        let mut npktd = npkt;
        for i in 2..=npk as usize {
            if self.mtp[i] >= 442 {
                npkd -= 1;
                npktd -= 1;
            }
        }
        let mut ilist = 1i32;
        c.nsc = 0;
        c.math = matd;
        c.mfh = 3;
        let mut inpk = 1i32;
        loop {
            // 105
            inpk += 1;
            c.mth = self.mtp[inpk as usize];
            scr[1] = self.za;
            scr[2] = self.awr;
            scr[3] = 0.0;
            scr[4] = 0.0;
            scr[5] = 0.0;
            scr[6] = 0.0;
            c.contio(0, 0, nscr, &mut scr.all(), &mut nb, &mut nw);
            let mut nn = 0;
            scr[1] = 0.0;
            scr[2] = 0.0;
            scr[5] = 1.0;
            scr[6] = 0.0;
            scr[7] = 0.0;
            scr[8] = 2.0;
            let mut ibase = 8;
            let mut i = 0;
            nw = 0;
            let mut allzero = false;
            loop {
                // 110
                nn += 1;
                if nn > ne {
                    allzero = true;
                    break;
                }
                c.finda(nn, &mut cv.0, npkk, iold, &mut buf, NBUF);
                for j in 1..=npkk {
                    if cv[j] > big {
                        cv[j] = sigfig(cv[j], 6, 0);
                    } else {
                        cv[j] = sigfig(cv[j], 7, 0);
                    }
                }
                let e = sigfig(cv[1], 9, 0);
                if self.iprint == 1 && (e - self.elist[ilist as usize]).abs() <= small * e {
                    let mut ihi = 0;
                    let mut ilo = 0;
                    if self.kchk == 1 {
                        ihi = 0;
                        ilo = 0;
                        for j in 2..=npkd {
                            klo[j as usize] = "    ";
                            let test = cv[j + npkt] - cv[j + npkt] / 10.0;
                            if cv[j] < test {
                                klo[j as usize] = "low ";
                                ilo = 1;
                            }
                            khi[j as usize] = "    ";
                            let test = cv[j + 2 * npkt] + cv[j + 2 * npkt] / 10.0;
                            if cv[j] > test {
                                khi[j as usize] = "high";
                                ihi = 1;
                            }
                        }
                    }
                    if ilo == 1 {
                        c.io.so(" ");
                    }
                    if self.kchk == 1 {
                        c.io.so(&a4_lines(&klo[2..=npkd.max(1) as usize]));
                    }
                    if self.kchk == 1 {
                        let v: Vec<f64> = (1..=npktd).map(|j| cv[npk + j]).collect();
                        c.io.so(&minmax_lines("min", &v));
                    }
                    let v: Vec<f64> = (1..=npk).map(|j| cv[j]).collect();
                    c.io.so(&e14_lines(&v, 9));
                    if self.kchk == 1 {
                        let v: Vec<f64> = (1..=npktd).map(|j| cv[npk + npkt + j]).collect();
                        c.io.so(&minmax_lines("max", &v));
                    }
                    if ihi == 1 {
                        c.io.so(&a4_lines(&khi[2..=npkd.max(1) as usize]));
                    }
                }
                if e >= self.elist[ilist as usize] {
                    ilist += 1;
                }
                if i <= 0 {
                    if cv[inpk] == zero && i == 0 && ibase == 8 {
                        continue;
                    }
                    if nw <= 0 {
                        n = ne - nn + 1;
                        scr[6] = n as f64;
                        scr[7] = n as f64;
                    }
                }
                // 118
                i += 2;
                scr[ibase + i - 1] = cv[1];
                scr[ibase + i] = cv[inpk];
                if i < c.npage && nn != ne {
                    continue;
                }
                if ibase == 0 {
                    // 120
                    c.moreio(0, 0, nscr, &mut scr.all(), &mut nb, &mut nw);
                    if nn == ne {
                        break;
                    }
                    i = 0;
                    continue;
                }
                c.tab1io(0, 0, nscr, &mut scr.all(), &mut nb, &mut nw);
                if nn == ne {
                    break;
                }
                ibase = 0;
                i = 0;
            }
            if allzero {
                // 125: all zeroes for this mt
                scr[6] = 2.0;
                scr[7] = 2.0;
                scr[9] = self.efirst;
                scr[10] = 0.0;
                scr[11] = self.elast;
                scr[12] = 0.0;
                nw = 12;
                c.tab1io(0, 0, nscr, &mut scr.all(), &mut nb, &mut nw);
                n = 2;
            }
            // 130
            c.asend(0, nscr);
            ncds[inpk as usize] = 3 + (n + 2) / 3;
            if inpk < npk {
                continue;
            }
            break;
        }
        c.afend(0, nscr);
        c.repoz(nscr);

        //--update dictionary in file 1
        c.findf(matd, 1, 451, nin);
        c.contio(nin, 0, 0, &mut scr.all(), &mut nb, &mut nw);
        let mut nx = c.n2h;
        let mut inow = 7;
        if c.iverf >= 5 {
            c.contio(nin, 0, 0, &mut scr.at(inow), &mut nb, &mut nw);
            inow += 6;
            if c.iverf == 6 {
                c.contio(nin, 0, 0, &mut scr.at(inow), &mut nb, &mut nw);
                inow += 6;
            }
        }
        c.hdatio(nin, 0, 0, &mut scr.at(inow), &mut nb, &mut nw);
        let nwd = nint(scr[inow + 4]);
        if c.iverf >= 5 {
            nx = nint(scr[inow + 5]);
        }
        inow += nw;
        while nb != 0 {
            c.moreio(nin, 0, 0, &mut scr.at(inow), &mut nb, &mut nw);
            inow += nw;
        }
        nw = nx;
        let idict = inow;
        let idict1 = idict - 1;
        c.dictio(nin, 0, 0, &mut scr.at(idict), &mut nb, &mut nw);
        nw = 6 * (nx + 20);
        if nw < c.npage + 50 {
            nw = c.npage + 50;
        }
        let mut b = rv(nw);
        let mut inpk = 2i32;
        let mut mtn = self.mtp[inpk as usize];
        let mut j = 0i32;

        //--loop over existing dictionary records
        for i in 1..=nx {
            let mut iowr = 0;
            let ia = 6 * (i - 1);
            let mfi = nint(scr[idict1 + ia + 3]);
            let mti = nint(scr[idict1 + ia + 4]);
            //--insert heating mt dictionary data (or overwrite
            //  original dictionary data with revised data).
            while (mfi == 3 && mtn <= mti) || (mfi > 3 && inpk <= npk) {
                if mfi == 3 && mtn == mti {
                    iowr = 1;
                }
                b[j + 1] = 0.0;
                b[j + 2] = 0.0;
                b[j + 3] = 3.0;
                b[j + 4] = self.mtp[inpk as usize] as f64;
                b[j + 5] = ncds[inpk as usize] as f64;
                b[j + 6] = 0.0;
                j += 6;
                inpk += 1;
                if inpk <= npk {
                    mtn = self.mtp[inpk as usize];
                } else {
                    mtn = 10000;
                }
            }
            //--insert original dictionary data (unless the heating
            //  mt entry above was an overwrite of existing data).
            if mfi != 3 || (mfi == 3 && iowr == 0) {
                b[j + 1] = 0.0;
                b[j + 2] = 0.0;
                b[j + 3] = mfi as f64;
                b[j + 4] = mti as f64;
                b[j + 5] = scr[idict1 + ia + 5];
                b[j + 6] = scr[idict1 + ia + 6];
                j += 6;
            }
        }
        //--insert rest of heating mt dictionary data, if any.
        for i in inpk..=npk {
            b[j + 1] = 0.0;
            b[j + 2] = 0.0;
            b[j + 3] = 3.0;
            b[j + 4] = self.mtp[i as usize] as f64;
            b[j + 5] = ncds[i as usize] as f64;
            b[j + 6] = 0.0;
            j += 6;
        }
        c.mth = 451;
        let nxn = j / 6;
        let mut inow = 1;
        if c.iverf <= 4 {
            scr[inow + 5] = nxn as f64;
        }
        c.contio(0, nout, 0, &mut scr.at(inow), &mut nb, &mut nw);
        inow += 6;
        if c.iverf == 6 {
            c.contio(0, nout, 0, &mut scr.at(inow), &mut nb, &mut nw);
        }
        if c.iverf == 6 {
            inow += 6;
        }
        if c.iverf > 4 {
            c.contio(0, nout, 0, &mut scr.at(inow), &mut nb, &mut nw);
            inow += 6;
            scr[inow + 5] = nxn as f64;
        }
        nw = 6 + nwd;
        c.hdatio(0, nout, 0, &mut scr.at(inow), &mut nb, &mut nw);
        inow += nw;
        while nb != 0 {
            c.moreio(0, nout, 0, &mut scr.at(inow), &mut nb, &mut nw);
            inow += nw;
        }
        nw = nxn;
        c.dictio(0, nout, 0, &mut b.all(), &mut nb, &mut nw);
        c.tofend(nin, nout, 0, &mut scr.all());
        c.tofend(nin, nout, 0, &mut scr.all());

        //--add new reactions to file 3.
        c.contio(nin, 0, 0, &mut scr.all(), &mut nb, &mut nw);
        if c.mfh != 3 || c.mth == 0 {
            let strng = format!("nin out of order.  read mfh,mth = {}{}", ffmt::i(c.mfh, 5), ffmt::i(c.mth, 5));
            c.error("hout", &strng, " ");
        }
        let mut mfnin = c.mfh;
        let mut mtnin = c.mth;
        let inow = 7;
        c.contio(nscr, 0, 0, &mut scr.at(inow), &mut nb, &mut nw);
        if c.mfh != 3 || c.mth == 0 {
            let strng = format!("nscr out of order.  read mfh,mth = {}{}", ffmt::i(c.mfh, 5), ffmt::i(c.mth, 5));
            c.error("hout", &strng, " ");
        }
        let mut mfnscr = c.mfh;
        let mut mtnscr = c.mth;
        let inow6 = inow + 6;
        while mfnin == 3 && mfnscr == 3 {
            if mtnin < mtnscr {
                c.mfh = mfnin;
                c.mth = mtnin;
                c.contio(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                c.tosend(nin, nout, 0, &mut scr.at(inow6));
                c.contio(nin, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                mfnin = c.mfh;
                mtnin = c.mth;
            } else if mtnin == mtnscr {
                c.mfh = mfnscr;
                c.mth = mtnscr;
                c.contio(0, nout, 0, &mut scr.at(inow), &mut nb, &mut nw);
                c.tosend(nscr, nout, 0, &mut scr.at(inow6));
                c.tosend(nin, 0, 0, &mut scr.at(inow6));
                c.contio(nin, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                mfnin = c.mfh;
                mtnin = c.mth;
                c.contio(nscr, 0, 0, &mut scr.at(inow), &mut nb, &mut nw);
                mfnscr = c.mfh;
                mtnscr = c.mth;
            } else {
                while mtnin >= mtnscr && mfnscr == 3 {
                    c.mfh = mfnscr;
                    c.mth = mtnscr;
                    c.contio(0, nout, 0, &mut scr.at(inow), &mut nb, &mut nw);
                    c.tosend(nscr, nout, 0, &mut scr.at(inow6));
                    c.contio(nscr, 0, 0, &mut scr.at(inow), &mut nb, &mut nw);
                    mfnscr = c.mfh;
                    mtnscr = c.mth;
                }
            }
        }
        if mfnin == 3 && mfnscr == 0 {
            c.mfh = mfnin;
            c.mth = mtnin;
            c.contio(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
        } else if mfnin == 0 && mfnscr == 3 {
            c.mfh = mfnscr;
            c.mth = mtnscr;
            c.contio(0, nout, 0, &mut scr.at(inow), &mut nb, &mut nw);
            c.tofend(nscr, nout, 0, &mut scr.at(inow6));
        } else if mfnin == 0 && mfnscr == 0 {
            c.mfh = mfnin;
            c.mth = mtnin;
            c.contio(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
        }
        c.tomend(nin, nout, 0, &mut scr.all());
        drop(b);
        c.closz(nscr);

        //--plot heating, photon energy production, and kinematic limits
        if self.nplot == 0 {
            return;
        }
        let nplot = self.nplot;
        let mt303 = self.mt303;
        let wl = |c: &mut Ctx, s: &str| c.io.write_line(nplot, s);
        let xy = |x: f64, y: f64| format!("{}{}/", ffmt::e1p(x, 14, 6), ffmt::e1p(y, 14, 6));
        wl(c, "/");

        // y values of the three curves (k=1 upper, 2 value, 3 lower)
        // for the heating plots (kind 0) and photon plots (kind 1)
        let yk = |cv: &RV, kind: i32, k: i32| -> f64 {
            if kind == 0 {
                match k {
                    1 => cv[2 + 2 * npkt],
                    2 => cv[2],
                    _ => cv[2 + npkt],
                }
            } else {
                match k {
                    1 => cv[npkk - 1] + cv[mt303] - cv[mt303 + npkt] + cv[npkk - 2],
                    2 => cv[npkk - 1],
                    _ => cv[npkk - 1] + cv[mt303] - cv[mt303 + 2 * npkt] + cv[npkk - 2],
                }
            }
        };
        let nkind = if mt303 == 0 { 1 } else { 2 };
        for kind in 0..nkind {
            for linlin in 0..2 {
                // log-log first, then lin-lin
                let elo0 = {
                    let mut elo = 1.0f64;
                    elo /= 1000.0;
                    elo
                };
                let (thin, xstart) = if linlin == 0 {
                    let rat = emax / elo0;
                    (powi_r(ten, log10(rat) / 2500.0), elo0)
                } else {
                    (emax / 2500.0, ecut)
                };
                let accept = |x: f64, xlast: f64| if linlin == 0 { x >= thin * xlast } else { x >= xlast + thin };
                let mut ylo = big;
                let mut yhi = -big;
                let mut xlast = xstart;
                for i in 1..=ne {
                    c.finda(i, &mut cv.0, npkk, iold, &mut buf, NBUF);
                    let x = cv[1];
                    if accept(x, xlast) {
                        if kind == 0 {
                            if cv[2] > yhi {
                                yhi = cv[2];
                            }
                            if cv[2] < ylo {
                                ylo = cv[2];
                            }
                            if cv[2 + npkt] > yhi {
                                yhi = cv[2 + npkt];
                            }
                            if cv[2 + npkt] < ylo {
                                ylo = cv[2 + npkt];
                            }
                            if cv[2 + 2 * npkt] > yhi {
                                yhi = cv[2 + 2 * npkt];
                            }
                            if cv[2 + 2 * npkt] < ylo {
                                ylo = cv[2 + 2 * npkt];
                            }
                        } else {
                            for k in [1, 2, 3] {
                                let y = yk(&cv, 1, k);
                                if y > yhi {
                                    yhi = y;
                                }
                                if y < ylo {
                                    ylo = y;
                                }
                            }
                        }
                        xlast = x;
                    }
                }
                if linlin == 0 {
                    if ylo <= zero {
                        ylo = small;
                    }
                    if yhi <= zero {
                        yhi = small;
                    }
                }
                let ylabel = if kind == 0 { "<h>eating (e<v>-barns)" } else { "<p>hoton energy prod (e<v>-barns)" };
                let vlabel = if kind == 0 { "<e>nergy-balance heating" } else { "<p>hoton energy production" };
                let upper = if kind == 0 && linlin == 1 { "<u>pper limit>" } else { "<u>pper limit" };
                let mut x = 0.0f64;
                for k in 1..=3 {
                    if k == 1 {
                        wl(c, "1/");
                        wl(c, "'Energy-Balance Check'/");
                        wl(c, "/");
                        wl(c, if linlin == 0 { "4 0 2 1/" } else { "1 0 2 1/" });
                        wl(c, "/");
                        wl(c, "'<e>nergy (e<v>)'/");
                        wl(c, "/");
                        wl(c, &format!("'{}'/", ylabel));
                        wl(c, "/");
                        wl(c, "0 0 4/");
                        wl(c, &format!("'{}'/", upper));
                    } else if k == 2 {
                        wl(c, "2/");
                        wl(c, "/");
                        wl(c, "/");
                        wl(c, &format!("'{}'/", vlabel));
                    } else {
                        wl(c, "3/");
                        wl(c, "/");
                        wl(c, "0 0 4/");
                        wl(c, "'<l>ower limit'/");
                    }
                    wl(c, "0/");
                    let mut xlast = xstart;
                    let mut jj = 0;
                    for i in 1..=ne {
                        c.finda(i, &mut cv.0, npkk, iold, &mut buf, NBUF);
                        x = cv[1];
                        let mut y = yk(&cv, kind, k);
                        if accept(x, xlast) && jj < 2500 {
                            jj += 1;
                            if linlin == 0 && y <= zero {
                                y = small;
                            }
                            wl(c, &xy(x, y));
                            xlast = x;
                        }
                    }
                    if k == 1 {
                        wl(c, &xy(x, ylo));
                    }
                    if k == 1 {
                        wl(c, &xy(x, yhi));
                    }
                    wl(c, "/ end");
                }
            }
        }

        //--plotting finished
        wl(c, "99/");
        c.closz(nplot);
    }
}

/// `ten**(r)` with a real exponent (libm pow, as gfortran).
fn powi_r(x: f64, r: f64) -> f64 {
    pow(x, r)
}

/// `'(15x,8(6x,a4,4x))'`: format reversion restarts at the group. An X
/// skip at the end of a record writes nothing (libgfortran), so the 4x
/// after the last item of each record does not appear.
fn a4_lines(items: &[&str]) -> String {
    let mut s = " ".repeat(15);
    for (k, it) in items.iter().enumerate() {
        if k > 0 && k % 8 == 0 {
            s.push('\n');
        } else if k > 0 {
            s.push_str(&" ".repeat(4));
        }
        s.push_str(&" ".repeat(6));
        s.push_str(it);
    }
    s
}

/// `'(12x,''min'',1p,8e14.4)'`: format reversion restarts the whole format.
fn minmax_lines(label: &str, vals: &[f64]) -> String {
    let mut s = String::new();
    for (k, v) in vals.iter().enumerate() {
        if k % 8 == 0 {
            if k > 0 {
                s.push('\n');
            }
            s.push_str(&" ".repeat(12));
            s.push_str(label);
        }
        s.push_str(&ffmt::e1p(*v, 14, 4));
    }
    if vals.is_empty() {
        s.push_str(&" ".repeat(12));
        s.push_str(label);
    }
    s
}
