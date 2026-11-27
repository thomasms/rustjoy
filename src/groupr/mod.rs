//! GROUPR: self-shielded multigroup cross sections, scattering matrices
//! and photon production matrices, ported from NJOY2016 `groupr.f90`
//! (module groupm).
//!
//! The routines keep the Fortran structure and 1-based indices. The
//! module variables and the saved locals (which keep their values from
//! one call, and one GROUPR run, to the next, as in the Fortran) are the
//! fields of `Groupr`. Statement labels that are jumped to become states
//! of small `loop { match lbl { .. } }` machines.

#![allow(clippy::too_many_lines)]
#![allow(clippy::too_many_arguments)]
#![allow(clippy::collapsible_if)]
#![allow(clippy::collapsible_else_if)]
#![allow(clippy::assign_op_pattern)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::manual_range_contains)]
#![allow(clippy::if_same_then_else)]
#![allow(clippy::neg_cmp_op_on_partial_ord)]

use std::ops::{Index, IndexMut};

mod conver;
mod core;
mod dist;
mod mf6;
mod weight;

pub use weight::{gengpg, gengpn};

use crate::endf::Ctx;
use crate::farr::{pack_text, unpack_text, F, R3, RM, RV};
use crate::ffmt;
use crate::fio::{Kind, RecW, Val};
use crate::fmath::{exp, log, nint, pow, powc, sqrt};
use crate::groupr_tables::*;
use crate::util::physics::*;
use crate::util::{a10, sigfig};

const NBUF: i32 = 1000;
const MAXR1: usize = 1000;
const MAXR2: usize = 10000;
/// Extra words past the Fortran array sizes.
const SLACK: usize = 4096;

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

/// The first n characters of a list-directed character value, blank padded.
fn getc(v: &Option<Val>, n: usize, d: &str) -> String {
    match v {
        Some(Val::Char(s)) => pad(s, n),
        _ => d.to_string(),
    }
}

/// `s` blank-padded to n characters and cut there.
fn pad(s: &str, n: usize) -> String {
    let mut t: Vec<u8> = s.bytes().take(n).collect();
    while t.len() < n {
        t.push(b' ');
    }
    String::from_utf8_lossy(&t).to_string()
}

/// `a11` editing of a 10-character field: one leading blank.
fn a11(s: &str) -> String {
    format!(" {}", s)
}

/// Lines written by a format `(head, n1 a11, (/, cont, n2 a11))` for the
/// given fields, as libgfortran writes them: a full first or inner group
/// is followed by the slash of the group (an empty record when no field
/// is left), and format reversion starts one more record.
fn rev_lines(head: &str, f: &[String], n1: usize, cont: &str, n2: usize) -> String {
    let mut s = String::from(head);
    let k1 = f.len().min(n1);
    for x in &f[..k1] {
        s.push_str(&a11(x));
    }
    if f.len() < n1 {
        return s;
    }
    s.push('\n');
    let mut k = k1;
    let mut first = true;
    loop {
        if k >= f.len() && !first {
            break;
        }
        if !first {
            // reversion: a new record, then the slash
            s.push('\n');
            s.push('\n');
        }
        first = false;
        let m = (f.len() - k).min(n2);
        if m == 0 {
            // the slash was done; nX before no item gives nothing
            break;
        }
        s.push_str(cont);
        for x in &f[k..k + m] {
            s.push_str(&a11(x));
        }
        k += m;
        if m < n2 {
            break;
        }
        if k >= f.len() {
            break;
        }
    }
    s
}

/// Column-major 2-D view with leading dimension `ld` on a flat slice
/// (`a(ld,*)` dummy arguments).
pub struct V2<'a> {
    pub v: &'a mut [f64],
    pub ld: usize,
}

impl<'a> V2<'a> {
    pub fn new(v: &'a mut [f64], ld: usize) -> Self {
        V2 { v, ld }
    }
    /// Flat 1-based view starting at element (i,j).
    pub fn at(&mut self, i: i32, j: i32) -> F<'_> {
        let k = (i - 1) as usize + (j - 1) as usize * self.ld;
        F(&mut self.v[k..])
    }
}

impl<'a> Index<(i32, i32)> for V2<'a> {
    type Output = f64;
    #[inline]
    fn index(&self, (i, j): (i32, i32)) -> &f64 {
        &self.v[(i - 1) as usize + (j - 1) as usize * self.ld]
    }
}
impl<'a> IndexMut<(i32, i32)> for V2<'a> {
    #[inline]
    fn index_mut(&mut self, (i, j): (i32, i32)) -> &mut f64 {
        &mut self.v[(i - 1) as usize + (j - 1) as usize * self.ld]
    }
}

/// 1-based integer list with a 2-D form `a(6,n)`.
#[derive(Clone)]
struct I2 {
    v: Vec<i32>,
}
impl I2 {
    fn new(n: usize) -> Self {
        I2 { v: vec![0; 6 * n] }
    }
    fn g(&self, i: i32, j: i32) -> i32 {
        self.v[(i - 1) as usize + 6 * (j - 1) as usize]
    }
    fn s(&mut self, i: i32, j: i32, x: i32) {
        self.v[(i - 1) as usize + 6 * (j - 1) as usize] = x;
    }
}

/// Module variables of `groupm` and the saved locals.
pub struct Groupr {
    // particle and gamma group structures
    ign: i32,
    ngn: i32,
    egn: RV,
    igg: i32,
    ngg: i32,
    egg: RV,
    // weighting function
    iwt: i32,
    jsigz: i32,
    jtemp: i32,
    wght: RV,
    wtbuf: RV,
    nfp: i32,
    nfv: i32,
    nfscr: i32,
    nflmax: i32,
    felo: f64,
    fehi: f64,
    sigpot: f64,
    alpha2: f64,
    alpha3: f64,
    beta: f64,
    sam: f64,
    gamma: f64,
    // temperatures and background cross sections
    temp: RV,
    sigz: RV,
    ntemp: i32,
    nsigz: i32,
    // input mat, Legendre order, print, smoothing, title
    matb: i32,
    lord: i32,
    iprint: i32,
    ismooth: i32,
    rtitle: [f64; 17],
    // units
    nendf: i32,
    npend: i32,
    nend2: i32,
    nend3: i32,
    nflx: i32,
    npend2: i32,
    ninwt: i32,
    ngout1: i32,
    ngout2: i32,
    ntw: i32,
    // auto reaction lists
    mf4: Vec<i32>,
    mf6: Vec<i32>,
    mf12: Vec<i32>,
    mf13: Vec<i32>,
    mf18: Vec<i32>,
    mf4r: I2,
    mf6p: I2,
    mf10f: Vec<i32>,
    mf10s: Vec<i32>,
    mf10i: Vec<i32>,
    lfs8: Vec<i32>,
    mlfs8: Vec<i32>,
    lastza: i32,
    izatest: i32,
    // reaction and kinematics parameters
    matd: i32,
    mfd: i32,
    mtd: i32,
    econst: f64,
    jconst: i32,
    ndelg: i32,
    dntc: [f64; 9],
    awr: f64,
    q: f64,
    thresh: f64,
    alpha: f64,
    lrflag: i32,
    awrp: f64,
    izap: i32,
    xtot: f64,
    zap: f64,
    awp: f64,
    aprime: f64,
    lip: i32,
    law: i32,
    lct: i32,
    izat: i32,
    spi: f64,
    emaxx: f64,
    ebeg: f64,
    lfs: i32,
    isom: i32,
    izar: i32,
    // unresolved resonance parameters
    nunr: i32,
    intunr: i32,
    lrp: i32,
    lssf: i32,
    unr: RV,
    // temporary arrays for reaction data
    yield_: RV,
    sigma: RV,
    fls: RV,
    flgn: RV,
    gyln: RV,
    gyle: RV,
    eyle: RV,
    gfle: RV,
    gfle_ld: usize,
    sedist: RV,
    sede: RV,
    sede_ld: usize,
    ddmf6: RV,
    aes: RV,
    stmp: RM,
    slst: RM,
    ftmp: RM,
    flst: RM,
    stmp_alloc: bool,
    falo: RV,
    fahi: RV,
    ipan: i32,
    // ---- saved locals
    nx_ir: i32,
    gw_ip: i32,
    gw_ir: i32,
    gw_ipl: i32,
    gw_step: f64,
    fw_el: f64,
    fw_en: f64,
    fw_xl: f64,
    fw_xn: f64,
    fw_wl: f64,
    fw_wn: f64,
    fw_templ: f64,
    p_elast: f64,
    p_idisc: i32,
    p_matl: i32,
    p_mfl: i32,
    p_mtl: i32,
    p_ig1w: i32,
    p_iglow: i32,
    p_eqw: f64,
    p_nq: i32,
    p_nqp: i32,
    p_ig1: i32,
    p_ng1: i32,
    p_enext: f64,
    /// panel local `ehigh` (not saved; a jump to 115 reaches 220 without
    /// setting it, so the value is the one from the last call)
    p_ehigh: f64,
    /// panel locals `en`, `idiscf` (same reason)
    p_en: f64,
    p_idiscf: i32,
    d_irat: i32,
    d_itr: i32,
    gf_ip: i32,
    gy_lnu: i32,
    gy_ip: i32,
    gy_ir: i32,
    gs_mt: i32,
    gs_nsig: i32,
    ff_ifirst: i32,
    ff_nyl: i32,
    ff_igmin: i32,
    /// main-program local `ig1` (kept from reaction to reaction)
    m_ig1: i32,
    sv: mf6::Saved,
}


impl Default for Groupr {
    fn default() -> Self {
        Groupr {
            ign: 0,
            ngn: 0,
            egn: RV::new(0),
            igg: 0,
            ngg: 0,
            egg: RV::new(0),
            iwt: 0,
            jsigz: 0,
            jtemp: 0,
            wght: RV::new(0),
            wtbuf: RV::new(0),
            nfp: 0,
            nfv: 0,
            nfscr: 0,
            nflmax: 0,
            felo: 0.0,
            fehi: 0.0,
            sigpot: 0.0,
            alpha2: 0.0,
            alpha3: 0.0,
            beta: 0.0,
            sam: 0.0,
            gamma: 0.0,
            temp: RV::new(0),
            sigz: RV::new(0),
            ntemp: 0,
            nsigz: 0,
            matb: 0,
            lord: 0,
            iprint: 0,
            ismooth: 0,
            rtitle: [0.0; 17],
            nendf: 0,
            npend: 0,
            nend2: 0,
            nend3: 0,
            nflx: 0,
            npend2: 0,
            ninwt: 0,
            ngout1: 0,
            ngout2: 0,
            ntw: 0,
            mf4: vec![0; MAXR1 + 1],
            mf6: vec![0; MAXR1 + 1],
            mf12: vec![0; MAXR1 + 1],
            mf13: vec![0; MAXR1 + 1],
            mf18: vec![0; MAXR1 + 1],
            mf4r: I2::new(MAXR1 + 1),
            mf6p: I2::new(MAXR1 + 1),
            mf10f: vec![0; MAXR2 + 1],
            mf10s: vec![0; MAXR2 + 1],
            mf10i: vec![0; MAXR2 + 1],
            lfs8: vec![0; MAXR2 + 1],
            mlfs8: vec![0; MAXR2 + 1],
            lastza: 0,
            izatest: 0,
            matd: 0,
            mfd: 0,
            mtd: 0,
            econst: 0.0,
            jconst: 0,
            ndelg: 0,
            dntc: [0.0; 9],
            awr: 0.0,
            q: 0.0,
            thresh: 0.0,
            alpha: 0.0,
            lrflag: 0,
            awrp: 0.0,
            izap: 0,
            xtot: 0.0,
            zap: 0.0,
            awp: 0.0,
            aprime: 0.0,
            lip: 0,
            law: 0,
            lct: 0,
            izat: 0,
            spi: 0.0,
            emaxx: 0.0,
            ebeg: 0.0,
            lfs: 0,
            isom: 0,
            izar: 0,
            nunr: 0,
            intunr: 0,
            lrp: 0,
            lssf: 0,
            unr: RV::new(0),
            yield_: RV::new(0),
            sigma: RV::new(0),
            fls: RV::new(0),
            flgn: RV::new(0),
            gyln: RV::new(0),
            gyle: RV::new(0),
            eyle: RV::new(0),
            gfle: RV::new(0),
            gfle_ld: 1,
            sedist: RV::new(0),
            sede: RV::new(0),
            sede_ld: 1,
            ddmf6: RV::new(0),
            aes: RV::new(0),
            stmp: RM::new(0, 0),
            slst: RM::new(0, 0),
            ftmp: RM::new(0, 0),
            flst: RM::new(0, 0),
            stmp_alloc: false,
            falo: RV::new(0),
            fahi: RV::new(0),
            ipan: 0,
            nx_ir: 0,
            gw_ip: 0,
            gw_ir: 0,
            gw_ipl: 0,
            gw_step: 0.0,
            fw_el: 0.0,
            fw_en: 0.0,
            fw_xl: 0.0,
            fw_xn: 0.0,
            fw_wl: 0.0,
            fw_wn: 0.0,
            fw_templ: 1.0e10,
            p_elast: 0.0,
            p_idisc: 0,
            p_matl: 0,
            p_mfl: 0,
            p_mtl: 0,
            p_ig1w: 0,
            p_iglow: 0,
            p_eqw: 0.0,
            p_nq: 0,
            p_nqp: 0,
            p_ig1: 0,
            p_ng1: 0,
            p_enext: 0.0,
            p_ehigh: 0.0,
            p_en: 0.0,
            p_idiscf: 0,
            d_irat: 0,
            d_itr: 0,
            gf_ip: 0,
            gy_lnu: 0,
            gy_ip: 0,
            gy_ir: 0,
            gs_mt: 0,
            gs_nsig: 0,
            ff_ifirst: 0,
            ff_nyl: 0,
            ff_igmin: 0,
            m_ig1: 0,
            sv: mf6::Saved::default(),
        }
    }
}

/// Fortran `allocate(a(n))` for a module array (the contents are not
/// set; a few words past the end are kept for safety).
fn alloc(n: i32) -> RV {
    RV::new(n.max(1) as usize + 8)
}

impl Groupr {
    pub fn new() -> Self {
        Self::default()
    }

    /// Compute self-shielded group-averaged cross sections.
    pub fn groupr(&mut self, c: &mut Ctx) {
        let small = 1.0e-6f64;
        let zero = 0.0f64;
        let (mut nb, mut nw) = (0i32, 0i32);
        let mut nwds: i32;
        let mut itend: i32;

        //--initialize
        let mut nwscr: i32 = 10000;
        self.nfscr = 10;
        self.nflx = 11;
        self.nend2 = 12;
        self.nend3 = 13;
        self.npend2 = 0;
        self.ninwt = 0;
        self.ipan = 0;
        c.openz(-self.nflx, 1);
        self.wtbuf = RV::new(NBUF as usize);
        self.ebeg = 10.0 * small;

        //--write heading and read user input
        let mut time = c.timer();
        c.io.so(&format!(
            "\n groupr...compute self-shielded group-averaged cross-sections{}{}s",
            " ".repeat(8),
            ffmt::f(time, 8, 1)
        ));
        c.io.se(&format!("\n groupr...{}{}s", " ".repeat(59), ffmt::f(time, 8, 1)));
        let mut iaddmt = 0i32;
        self.ruinb(c, &mut iaddmt);
        nwscr = nwscr.max(self.nsigz + self.ngn + self.ngg + 9);
        let mut scr = RV::new(nwscr.max(c.npage + 50) as usize + SLACK);
        if self.nendf < 0 {
            self.nend2 = -self.nend2;
        }
        if self.nendf < 0 {
            self.nend3 = -self.nend3;
        }
        c.openz(self.nend2, 1);
        c.openz(self.nend3, 1);
        if self.npend < 0 {
            self.nfscr = -self.nfscr;
        }
        c.openz(self.nfscr, 1);
        let mut naed = 14;
        if self.npend < 0 {
            naed = -naed;
        }
        c.openz(naed, 1);

        //--determine endf format being used
        c.tpidio(self.nendf, 0, 0, &mut scr.at(1), &mut nb, &mut nw);
        c.contio(self.nendf, 0, 0, &mut scr.at(1), &mut nb, &mut nw);
        c.contio(self.nendf, 0, 0, &mut scr.at(1), &mut nb, &mut nw);
        if c.n1h != 0 {
            c.iverf = 4;
        } else if c.n2h == 0 {
            c.iverf = 5;
        } else {
            c.iverf = 6;
        }
        c.io.so(&format!("\n using endf-{} format", ffmt::i(c.iverf, 1)));
        if c.iverf <= 4 {
            self.emaxx = 1.5e7;
        } else if c.iverf == 5 {
            self.emaxx = 2.0e7;
        } else {
            c.contio(self.nendf, 0, 0, &mut scr.at(1), &mut nb, &mut nw);
            self.emaxx = 2.0e7;
            if c.c2h > zero {
                self.emaxx = c.c2h;
            }
        }

        //--locate position for new material on old gout tape.
        c.repoz(self.ngout1);
        c.repoz(self.ngout2);
        c.math = 0;
        c.mfh = 0;
        c.mth = 0;
        c.nsh = 0;
        itend = 0;
        nwds = 17;
        if self.ngout1 == 0 {
            let mut rt = self.rtitle;
            c.tpidio(0, self.ngout2, 0, &mut F(&mut rt[..]), &mut nb, &mut nwds);
            self.rtitle = rt;
        }
        self.ntw = 1;
        self.rtitle[0] = 0.0;
        c.repoz(self.nendf);
        c.repoz(self.npend);
        c.tpidio(self.nendf, 0, 0, &mut scr.at(1), &mut nb, &mut nw);
        c.tpidio(self.npend, 0, 0, &mut scr.at(1), &mut nb, &mut nw);
        let mut lbl: i32 = 170;
        if !(self.ngout1 == 0 || self.ngout2 == 0) {
            c.tpidio(self.ngout1, self.ngout2, 0, &mut scr.at(1), &mut nb, &mut nwds);
            if iaddmt <= 0 {
                lbl = 120;
            }
        }

        let mut za = 0.0f64;
        let mut tempin: f64;
        let mut mtname = " ".repeat(60);
        let mut ig1 = self.m_ig1;
        let mut nz: i32 = 1;
        loop {
            match lbl {
                120 => {
                    c.contio(self.ngout1, 0, 0, &mut scr.at(1), &mut nb, &mut nw);
                    if c.math == -1 {
                        lbl = 160;
                    } else if c.math == 1 {
                        lbl = 120;
                    } else if c.math == self.matb {
                        // 140
                        c.tomend(self.ngout1, 0, 0, &mut scr.at(1));
                        lbl = 170;
                    } else if c.math > self.matb {
                        lbl = 160;
                    } else {
                        c.contio(0, self.ngout2, 0, &mut scr.at(1), &mut nb, &mut nw);
                        c.tomend(self.ngout1, self.ngout2, 0, &mut scr.at(1));
                        lbl = 120;
                    }
                }
                160 => {
                    c.skiprz(self.ngout1, -1);
                    lbl = 170;
                }
                170 => {
                    //--locate desired material on endf tape
                    //--make two copies for use by retrieval routines.
                    c.repoz(self.nendf);
                    if iaddmt < 0 {
                        c.contio(self.nendf, 0, 0, &mut scr.at(1), &mut nb, &mut nw);
                        if c.math < 0 {
                            lbl = 650;
                            continue;
                        }
                        self.matb = c.math;
                        c.skiprz(self.nendf, -1);
                    } else {
                        c.findf(self.matb, 1, 451, self.nendf);
                    }
                    c.repoz(self.nend2);
                    c.repoz(self.nend3);
                    c.nsh = 1;
                    c.nsc = 1;
                    c.math = 1;
                    c.afend(self.nend2, self.nend3);
                    self.conver(c, self.nendf, self.nend2, self.nend3);
                    c.atend(self.nend2, self.nend3);
                    c.repoz(self.nend2);
                    c.repoz(self.nend3);

                    //--search for desired mat and temperatures on pendf tape
                    //--store unresolved resonance cross sections if present.
                    for itemp in 1..=self.ntemp {
                        c.findf(self.matb, 1, 451, self.npend);
                        c.contio(self.npend, 0, 0, &mut scr.at(1), &mut nb, &mut nw);
                        self.lrp = c.l1h;
                        self.jtemp = itemp;
                        tempin = self.temp[itemp];
                        let eps = tempin / 10000.0;
                        let mut diff = 100000000.0f64;
                        while diff > 1.0 + eps {
                            za = c.c1h;
                            self.awr = c.c2h;
                            if c.iverf >= 5 {
                                c.contio(self.npend, 0, 0, &mut scr.at(1), &mut nb, &mut nw);
                            }
                            if c.iverf >= 6 {
                                c.contio(self.npend, 0, 0, &mut scr.at(1), &mut nb, &mut nw);
                            }
                            c.hdatio(self.npend, 0, 0, &mut scr.at(1), &mut nb, &mut nw);
                            diff = (c.c1h - tempin).abs();
                            if diff > 1.0 + eps {
                                if c.c1h > tempin {
                                    let s = format!(
                                        "unable to find mat={} t={}",
                                        ffmt::i(self.matb, 4),
                                        ffmt::e1p(tempin, 12, 4)
                                    );
                                    c.error("groupr", &s, " ");
                                }
                                c.tomend(self.npend, 0, 0, &mut scr.at(1));
                                c.contio(self.npend, 0, 0, &mut scr.at(1), &mut nb, &mut nw);
                                if c.math != self.matb {
                                    let s = format!(
                                        "unable to find mat={} t={}",
                                        ffmt::i(self.matb, 4),
                                        ffmt::e1p(tempin, 12, 4)
                                    );
                                    c.error("groupr", &s, " ");
                                }
                            } else {
                                if diff > small {
                                    let s = format!(
                                        "   difference between temperatures desired  and found is {}",
                                        ffmt::e1p(diff, 10, 2)
                                    );
                                    c.io.so(&s);
                                    c.io.se(&s);
                                }
                                let mut z = [0.0f64; 17];
                                for i in 1..=17 {
                                    z[i - 1] = scr[6 + i as i32];
                                }
                                let tz = unpack_text(&z);
                                c.io.so(&format!(
                                    "\n processing mat {}\n {}\n  {}",
                                    ffmt::i(self.matb, 6),
                                    "-".repeat(38),
                                    String::from_utf8_lossy(&tz)
                                ));
                                c.io.so(" ");
                                self.stounr(c, self.matb, tempin, self.npend);
                            }
                        }

                        //--generate weighting flux on nflx.
                        if self.nsigz > 1 {
                            let mut nwf = self.nflmax * (3 + self.nsigz);
                            if nwf == 0 {
                                nwf = 1;
                            }
                            let mut flxc = RV::new(nwf as usize + 8);
                            self.genflx(c, &mut flxc);
                        }

                        //--write head record for this material on gout tape.
                        if self.ngout2 != 0 {
                            if iaddmt == 0 {
                                nz = self.nsigz;
                                c.nsh = 1;
                                c.math = self.matb;
                                c.mfh = 1;
                                c.mth = 451;
                                scr[1] = za;
                                scr[2] = self.awr;
                                scr[3] = 0.0;
                                scr[4] = nz as f64;
                                scr[5] = -1.0;
                                scr[6] = self.ntw as f64;
                                c.contio(0, self.ngout2, 0, &mut scr.at(1), &mut nb, &mut nw);
                                scr[1] = tempin;
                                scr[2] = 0.0;
                                scr[3] = self.ngn as f64;
                                scr[4] = self.ngg as f64;
                                scr[5] = 0.0;
                                scr[6] = 0.0;
                                nw = 6;
                                for i in 1..=self.ntw {
                                    nw += 1;
                                    scr[nw] = self.rtitle[(i - 1) as usize];
                                }
                                for i in 1..=nz {
                                    nw += 1;
                                    scr[nw] = self.sigz[i];
                                }
                                let ngnp1 = self.ngn + 1;
                                for i in 1..=ngnp1 {
                                    nw += 1;
                                    scr[nw] = self.egn[i];
                                }
                                let nggp1 = self.ngg + 1;
                                for i in 1..=nggp1 {
                                    nw += 1;
                                    scr[nw] = self.egg[i];
                                }
                                if self.ngg == 0 {
                                    scr[nw] = 0.0;
                                }
                                let np = nw - 6;
                                scr[5] = np as f64;
                                let mut loc = 1;
                                c.listio(0, self.ngout2, 0, &mut scr.at(loc), &mut nb, &mut nw);
                                while nb != 0 {
                                    loc += nw;
                                    c.moreio(0, self.ngout2, 0, &mut scr.at(loc), &mut nb, &mut nw);
                                }
                                c.afend(self.ngout2, 0);
                            } else {
                                c.contio(self.ngout1, self.ngout2, 0, &mut scr.at(1), &mut nb, &mut nw);
                                c.tofend(self.ngout1, self.ngout2, 0, &mut scr.at(1));
                            }
                        }

                        //--main entry point for loop over reactions
                        let mut ngi = self.ngn;
                        if self.izap == 0 {
                            ngi = self.ngg;
                        }
                        self.matd = self.matb;
                        self.mtd = 1;
                        let mut mtdp: i32 = 1;
                        let mut iauto: i32 = 0;
                        let mut need_name = false;
                        let mut izam: i32 = 0;
                        let mut fzam: f64 = 0.0;
                        let mut rl: i32 = 360;
                        loop {
                            match rl {
                                360 => {
                                    self.lfs = 0;
                                    self.isom = 0;
                                    izam = 0;
                                    fzam = 0.0;
                                    if self.mtd + mtdp < 0 {
                                        rl = 400;
                                        continue;
                                    }
                                    if iauto == 0 {
                                        rl = 365;
                                        continue;
                                    }
                                    self.nextr(c, &mut iauto, self.matd, &mut mtdp, &mut scr);
                                    rl = if mtdp != 0 { 390 } else { 365 };
                                }
                                365 => {
                                    mtdp = -1000;
                                    let v = c.io.input.read(&[Kind::Int, Kind::Int, Kind::Char(60)]);
                                    self.mfd = geti(&v[0], self.mfd);
                                    mtdp = geti(&v[1], mtdp);
                                    let strng = getc(&v[2], 60, &" ".repeat(60));
                                    let mfd = self.mfd;
                                    if mfd < -1 || mfd == 1 || mfd == 2 || mfd == 4 {
                                        rl = 381;
                                        continue;
                                    }
                                    if mfd == 7 || mfd == 9 || mfd == 11 {
                                        rl = 381;
                                        continue;
                                    }
                                    if mfd == 14 || mfd == 15 {
                                        rl = 381;
                                        continue;
                                    }
                                    if mfd > 18 && mfd < 21 {
                                        rl = 381;
                                        continue;
                                    }
                                    if mfd > 26 && mfd < 31 {
                                        rl = 381;
                                        continue;
                                    }
                                    if mfd > 36 && mfd < 10000000 {
                                        rl = 381;
                                        continue;
                                    }
                                    if mfd >= 12 && mfd <= 18 && self.igg == 0 {
                                        c.error("groupr", "photons not allowed with igg=0.", " ");
                                    }
                                    if mfd == 0 {
                                        rl = 590;
                                        continue;
                                    }
                                    if mfd == -1 {
                                        rl = 383;
                                        continue;
                                    }
                                    if mtdp == -1000 {
                                        rl = 382;
                                        continue;
                                    }
                                    mtname = strng;
                                    rl = 390;
                                }
                                381 => {
                                    let s = format!("illegal mfd ={}", ffmt::i(self.mfd, 3));
                                    c.error("groupr", &s, " ");
                                }
                                382 => {
                                    iauto = 1;
                                    need_name = true;
                                    mtdp = 0;
                                    self.nextr(c, &mut iauto, self.matd, &mut mtdp, &mut scr);
                                    if mtdp != 0 {
                                        rl = 390;
                                        continue;
                                    }
                                    let s = format!("auto finds no reactions for mf={}", ffmt::i(self.mfd, 3));
                                    c.mess("groupr", &s, " ");
                                    rl = 365;
                                }
                                383 => {
                                    // read card9a format and down-convert
                                    let v = c.io.input.read(&[Kind::Int, Kind::Int, Kind::Int, Kind::Int]);
                                    self.mfd = geti(&v[0], self.mfd);
                                    mtdp = geti(&v[1], mtdp);
                                    self.izar = geti(&v[2], self.izar);
                                    self.lfs = geti(&v[3], self.lfs);
                                    if mtdp >= 0 {
                                        self.mtd = mtdp;
                                    }
                                    self.isom = self.lfs;
                                    match self.mfd {
                                        3 => self.mfd = 10000000,
                                        6 => self.mfd = 20000000,
                                        9 => self.mfd = 30000000,
                                        10 => self.mfd = 40000000,
                                        _ => {
                                            let s = format!("illegal file ={}", ffmt::i(self.mfd, 3));
                                            c.error("groupr", &s, " ");
                                        }
                                    }
                                    self.mfd = self.mfd + self.izar * 10;
                                    izam = self.izar * 1000 + self.lfs;
                                    fzam = self.izar as f64 + self.lfs as f64 / 1000.0;
                                    need_name = true;
                                    rl = 400;
                                }
                                390 => {
                                    if mtdp >= 0 {
                                        self.mtd = mtdp;
                                    }
                                    let mfd = self.mfd;
                                    if iauto == 0 && mfd >= 10000000 {
                                        if mfd == 40000000 {
                                            self.lfs = 0;
                                            self.isom = 0;
                                            izam = -1;
                                            self.izar = -1;
                                        } else {
                                            let mut itmp = mfd / 10000000;
                                            itmp = (mfd - 10000000 * itmp) / 10;
                                            self.izar = itmp;
                                            self.lfs = mfd - (10000000 * (mfd / 10000000) + 10 * itmp);
                                            self.isom = self.lfs;
                                            if self.lfs < 10 {
                                                izam = mfd % 10000000;
                                            } else {
                                                izam = 10 * (mfd % 10000000);
                                            }
                                        }
                                    } else if iauto == 1 && mfd >= 10000000 {
                                        if mfd == 40000000 {
                                            izam = -1;
                                            self.izar = -1;
                                        } else {
                                            self.izar = (mfd - ((mfd / 10000000) * 10000000)) / 10;
                                            if self.lfs < 10 {
                                                izam = mfd % 10000000 + self.lfs;
                                            } else {
                                                izam = 10 * (mfd % 10000000) + self.lfs;
                                            }
                                        }
                                    }
                                    rl = 400;
                                }
                                400 => {
                                    if mtdp < 0 {
                                        self.mtd += 1;
                                    }
                                    if iauto == 0 && (self.mfd == 13 || self.mfd == 16) {
                                        let (mut mf, mt) = (self.mfd, self.mtd);
                                        self.mfchk(&mut mf, mt);
                                        self.mfd = mf;
                                    }
                                    if iauto == 0 && self.mfd >= 21 && self.mfd <= 26 {
                                        let (mut mf, mt) = (self.mfd, self.mtd);
                                        self.mfchk2(&mut mf, mt);
                                        self.mfd = mf;
                                    }
                                    if need_name {
                                        mtname = namer(c, self.izap, self.izar, self.lfs, self.mfd, self.mtd);
                                    }
                                    rl = if iaddmt <= 0 { 490 } else { 410 };
                                }
                                410 => {
                                    c.contio(self.ngout1, 0, 0, &mut scr.at(1), &mut nb, &mut nw);
                                    if c.math == 0 && c.mfh == 0 && c.mth == 0 {
                                        c.skiprz(self.ngout1, -1);
                                        rl = 490;
                                    } else if c.mfh == 0 || c.mth == 0 {
                                        rl = 410;
                                    } else if c.mfh != self.mfd || c.mth != self.mtd {
                                        c.contio(0, self.ngout2, 0, &mut scr.at(1), &mut nb, &mut nw);
                                        c.tosend(self.ngout1, self.ngout2, 0, &mut scr.at(1));
                                        rl = 410;
                                    } else {
                                        c.tosend(self.ngout1, 0, 0, &mut scr.at(1));
                                        rl = 490;
                                    }
                                }
                                490 => {
                                    time = c.timer();
                                    c.skiprz(self.npend, -1);
                                    let (mut ng2g, mut nlg) = (0i32, 0i32);
                                    self.init(c, &mut ng2g, &mut nlg, &mut nz);
                                    let mut ans = R3::new(nlg as usize, nz as usize, ng2g as usize);
                                    let mut ng2 = ng2g;
                                    let mut nl = nlg;
                                    let mut ff = vec![0.0f64; (nlg * ng2g) as usize + SLACK];
                                    let mut prod = R3::new(nlg as usize, nz as usize, ngi.max(1) as usize);
                                    let mut flux = RM::new(nz as usize, nl as usize);
                                    let mut sig = RM::new(nz as usize, nl as usize);
                                    let ee = 0.0f64;
                                    let mut first = 0.0f64;
                                    let mut en = 0.0f64;
                                    let mut idis = 0i32;
                                    let mut yld = 0.0f64;
                                    let mut iglo = 0i32;
                                    let mut nq = 0i32;
                                    let mut ig2lo = 0i32;
                                    let mut igzero = 0i32;
                                    let mut enext = 0.0f64;
                                    let (mut elo, mut ehi): (f64, f64);
                                    if self.mfd != 5 {
                                        self.getsig(c, ee, &mut first, &mut idis, &mut sig, nl, nz);
                                        self.getflx(c, ee, &mut en, &mut idis, &mut flux, nl, nz);
                                    }
                                    self.getff(
                                        c,
                                        ee,
                                        &mut en,
                                        &mut idis,
                                        &mut yld,
                                        &mut V2::new(&mut ff, nlg as usize),
                                        &mut nl,
                                        &mut ng2,
                                        &mut iglo,
                                        &mut nq,
                                        nlg,
                                        ng2g,
                                    );
                                    let test = -1.0f64;
                                    let mut sl: i32;
                                    let mut ig: i32 = 0;
                                    if en == test {
                                        sl = 582;
                                    } else {
                                        ig1 = 0;
                                        ig = 0;
                                        sl = 0;
                                    }
                                    //--spectrum calculation (mfd=5 or constant spectra)
                                    let mut nll = 1i32;
                                    if sl == 0 {
                                        sl = if self.mfd == 5 {
                                            502
                                        } else if self.econst == zero {
                                            510
                                        } else {
                                            505
                                        };
                                    }
                                    loop {
                                        match sl {
                                            502 | 505 => {
                                                if sl == 502 && self.mtd == 455 {
                                                    for i in 1..=self.ndelg {
                                                        for iz in 1..=nz {
                                                            ans[(i, iz, 1)] = self.dntc[i as usize];
                                                        }
                                                    }
                                                    nll = self.ndelg;
                                                }
                                                // 505
                                                let en0 = en;
                                                self.getff(
                                                    c,
                                                    en0,
                                                    &mut enext,
                                                    &mut idis,
                                                    &mut yld,
                                                    &mut V2::new(&mut ff, nlg as usize),
                                                    &mut nll,
                                                    &mut ng2,
                                                    &mut ig2lo,
                                                    &mut nq,
                                                    nlg,
                                                    ng2g,
                                                );
                                                let f2 = V2::new(&mut ff, nlg as usize);
                                                for il in 1..=nll {
                                                    for i in 1..=ng2 {
                                                        let mut j = i;
                                                        if self.mtd == 455 {
                                                            j = i + 1;
                                                        }
                                                        for iz in 1..=nz {
                                                            ans[(il, iz, j)] = f2[(il, i)];
                                                        }
                                                    }
                                                }
                                                if self.mfd == 5 {
                                                    ig = ngi;
                                                }
                                                if self.mtd == 455 {
                                                    ng2 += 1;
                                                }
                                                sl = 540;
                                            }
                                            510 => {
                                                //--loop over initial energy groups
                                                ig += 1;
                                                if ngi == self.ngn {
                                                    elo = self.egn[ig];
                                                    ehi = self.egn[ig + 1];
                                                } else {
                                                    elo = self.egg[ig];
                                                    ehi = self.egg[ig + 1];
                                                }
                                                ig2lo = 0;
                                                ng2 = 2;
                                                if ehi <= first {
                                                    sl = 580;
                                                    continue;
                                                }
                                                if elo >= self.emaxx && ig != ngi {
                                                    sl = 580;
                                                    continue;
                                                }
                                                enext = ehi;
                                                for x in ans.v.iter_mut() {
                                                    *x = 0.0;
                                                }
                                                loop {
                                                    // 530
                                                    self.panel(
                                                        c,
                                                        &mut elo,
                                                        &mut enext,
                                                        &mut ans,
                                                        &mut ff,
                                                        &mut nl,
                                                        nz,
                                                        &mut ng2,
                                                        &mut ig2lo,
                                                        nlg,
                                                        ng2g,
                                                    );
                                                    if enext == ehi {
                                                        break;
                                                    }
                                                    elo = enext;
                                                    enext = ehi;
                                                }
                                                sl = 540;
                                            }
                                            540 => {
                                                //--accumulate production below econst, and print it
                                                if self.econst > 0.0 && self.mfd != 5 && ig != 0 && ig <= self.jconst {
                                                    for iz in 1..=nz {
                                                        for il in 1..=nl {
                                                            ans[(il, iz, 2)] = ans[(il, iz, 2)] / ans[(il, iz, 1)];
                                                            prod[(il, iz, ig)] = ans[(il, iz, 2)];
                                                        }
                                                    }
                                                    ig2lo = 0;
                                                    if ig >= self.jconst {
                                                        let mut jc = self.jconst;
                                                        self.displa(c, ig, &mut prod, nl, nz, &mut jc, &mut ig2lo, &mut igzero);
                                                        self.jconst = jc;
                                                    }
                                                } else {
                                                    //--print results for this initial energy group
                                                    if ig1 <= 0 {
                                                        ig1 = 1;
                                                        if tempin != 0.0 {
                                                            c.io.so(&format!(
                                                                " group constants at t={} deg k{}{}s",
                                                                ffmt::e1p(tempin, 9, 3),
                                                                " ".repeat(32),
                                                                ffmt::f(time, 8, 1)
                                                            ));
                                                        }
                                                        if tempin == 0.0 {
                                                            c.io.so(&format!(
                                                                "\n group constants at t=zero deg k{}{}s",
                                                                " ".repeat(37),
                                                                ffmt::f(time, 8, 1)
                                                            ));
                                                        }
                                                        if self.mfd < 1000000 {
                                                            c.io.so(&format!(
                                                                " for mf{} and mt{} {}",
                                                                ffmt::i(self.mfd, 3),
                                                                ffmt::i(self.mtd, 3),
                                                                mtname
                                                            ));
                                                        } else {
                                                            let mfdn = self.mfd / 10000000;
                                                            let jzam = izam;
                                                            let head = match mfdn {
                                                                1 => Some(" for mf3 mt"),
                                                                2 => Some(" for mf3*mf6 mt"),
                                                                3 => Some(" for mf3*mf9 mt"),
                                                                4 => Some(" for mf10 mt"),
                                                                _ => None,
                                                            };
                                                            if let Some(h) = head {
                                                                c.io.so(&format!(
                                                                    "{}{} zam{} {}",
                                                                    h,
                                                                    ffmt::i(self.mtd, 3),
                                                                    ffmt::i(jzam, 10),
                                                                    mtname
                                                                ));
                                                            }
                                                        }
                                                        if self.lrflag > 0 {
                                                            c.io.so(&format!(
                                                                "{}lr{} particle emission",
                                                                " ".repeat(14),
                                                                ffmt::i(self.lrflag, 3)
                                                            ));
                                                        }
                                                        if self.econst > zero {
                                                            c.io.so(&format!(
                                                                "\n spectrum constant below {} ev  ({} groups)",
                                                                ffmt::e1p(self.econst, 10, 4),
                                                                ffmt::i(self.jconst, 4)
                                                            ));
                                                        }
                                                        c.math = self.matb;
                                                        c.mfh = self.mfh_out();
                                                        c.mth = self.mtd;
                                                        scr[1] = za;
                                                        scr[2] = izam as f64;
                                                        if fzam != 0.0 {
                                                            scr[2] = fzam;
                                                        }
                                                        scr[3] = nl as f64;
                                                        scr[4] = nz as f64;
                                                        scr[5] = self.lrflag as f64;
                                                        scr[6] = ngi as f64;
                                                        c.contio(0, self.ngout2, 0, &mut scr.at(1), &mut nb, &mut nw);
                                                        self.displa(c, -1, &mut ans, nl, nz, &mut ng2, &mut ig2lo, &mut igzero);
                                                    }
                                                    self.displa(c, ig, &mut ans, nl, nz, &mut ng2, &mut ig2lo, &mut igzero);
                                                }

                                                //--write group constants on output tape
                                                let lim = nl * nz * ng2;
                                                if self.ngout2 != 0 {
                                                    if igzero != 0 || ig == ngi {
                                                        c.math = self.matb;
                                                        c.mfh = self.mfh_out();
                                                        c.mth = self.mtd;
                                                        scr[1] = tempin;
                                                        scr[2] = 0.0;
                                                        scr[3] = ng2 as f64;
                                                        scr[4] = ig2lo as f64;
                                                        scr[5] = lim as f64;
                                                        scr[6] = ig as f64;
                                                        let mut i = 0;
                                                        let mut j = 6;
                                                        let mut ibase = 6;
                                                        for it in 1..=ng2 {
                                                            for iz in 1..=nz {
                                                                for il in 1..=nl {
                                                                    i += 1;
                                                                    j += 1;
                                                                    scr[j] = ans[(il, iz, it)];
                                                                    if j >= (c.npage + ibase) || i == lim {
                                                                        if ibase != 0 {
                                                                            c.listio(0, self.ngout2, 0, &mut scr.at(1), &mut nb, &mut j);
                                                                            ibase = 0;
                                                                            j = 0;
                                                                        } else {
                                                                            c.moreio(0, self.ngout2, 0, &mut scr.at(1), &mut nb, &mut nw);
                                                                            j = 0;
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                                sl = 580;
                                            }
                                            580 => {
                                                if ig < ngi && self.mfd != 5 {
                                                    sl = 510;
                                                    continue;
                                                }
                                                if ig1 == 0 {
                                                    c.io.so(&format!(
                                                        "\n threshold is above highest energy bound for mt{}",
                                                        ffmt::i(self.mtd, 3)
                                                    ));
                                                }
                                                c.io.so(" ");
                                                break;
                                            }
                                            _ => {
                                                // 582
                                                c.io.so(&format!("\n no photons in mf6, mt{}", ffmt::i(c.mth, 3)));
                                                break;
                                            }
                                        }
                                    }
                                    //--loop over other desired reactions
                                    // 586
                                    self.sigma = RV::new(0);
                                    self.yield_ = RV::new(0);
                                    self.fls = RV::new(0);
                                    self.flgn = RV::new(0);
                                    self.gyln = RV::new(0);
                                    self.gyle = RV::new(0);
                                    self.eyle = RV::new(0);
                                    self.gfle = RV::new(0);
                                    self.sedist = RV::new(0);
                                    self.sede = RV::new(0);
                                    self.ddmf6 = RV::new(0);
                                    self.aes = RV::new(0);
                                    if c.mth > 0 {
                                        c.tosend(self.npend, 0, 0, &mut scr.at(1));
                                    }
                                    if ig1 != 0 {
                                        c.asend(self.ngout2, 0);
                                    }
                                    rl = 360;
                                }
                                590 => {
                                    //--loop over desired temperatures
                                    self.unr = RV::new(0);
                                    if iaddmt > 0 {
                                        // copy rest of material if necessary
                                        // 600
                                        c.contio(self.ngout1, 0, 0, &mut scr.at(1), &mut nb, &mut nw);
                                        if c.math > 0 && c.math != self.matd {
                                            // 620
                                            c.skiprz(self.ngout1, -2);
                                        } else {
                                            // 610
                                            c.contio(0, self.ngout2, 0, &mut scr.at(1), &mut nb, &mut nw);
                                            if c.math < 0 {
                                                itend = 1;
                                            }
                                            if c.math > 0 {
                                                c.tomend(self.ngout1, self.ngout2, 0, &mut scr.at(1));
                                            }
                                        }
                                        break;
                                    }
                                    c.amend(self.ngout2, 0);
                                    if itemp == self.ntemp {
                                        break;
                                    }
                                    c.tomend(self.npend, 0, 0, &mut scr.at(1));
                                    c.contio(self.npend, 0, 0, &mut scr.at(1), &mut nb, &mut nw);
                                    if c.math != self.matb {
                                        c.error("groupr", "unable to find next temp.", " ");
                                    }
                                    break;
                                }
                                _ => unreachable!(),
                            }
                        }
                        // 630
                    }
                    c.tomend(self.nendf, 0, 0, &mut scr.at(1));

                    //--loop over desired materials.
                    if iaddmt < 0 {
                        lbl = 160;
                        continue;
                    }
                    let v = c.io.input.read(&[Kind::Int]);
                    self.matb = geti(&v[0], self.matb);
                    iaddmt = 0;
                    if self.matb < 0 && self.ngout1 != 0 {
                        iaddmt = 1;
                    }
                    if self.matb < 0 && self.ngout1 == 0 {
                        iaddmt = -1;
                    }
                    self.matb = self.matb.abs();
                    if self.matb == 0 {
                        lbl = 650;
                    } else if self.ngout1 == 0 {
                        lbl = 160;
                    } else {
                        lbl = 120;
                    }
                }
                _ => {
                    // 650: groupr is complete
                    break;
                }
            }
        }
        self.m_ig1 = ig1;
        if self.ngout1 == 0 {
            c.atend(self.ngout2, 0);
        }
        if self.ngout1 != 0 && itend == 0 {
            c.totend(self.ngout1, self.ngout2, 0, &mut scr.at(1));
        }
        self.egn = RV::new(0);
        self.egg = RV::new(0);
        self.wtbuf = RV::new(0);
        self.wght = RV::new(0);
        if self.stmp_alloc {
            self.stmp_alloc = false;
            self.stmp = RM::new(0, 0);
            self.slst = RM::new(0, 0);
            self.ftmp = RM::new(0, 0);
            self.flst = RM::new(0, 0);
        }
        self.falo = RV::new(0);
        self.fahi = RV::new(0);
        c.repoz(self.nendf);
        c.repoz(self.npend);
        c.repoz(self.ngout1);
        c.repoz(self.ngout2);
        c.repoz(self.ninwt);
        c.closz(self.nendf);
        c.closz(self.npend);
        c.closz(self.ngout1);
        c.closz(self.ngout2);
        c.closz(self.nend2);
        c.closz(self.nend3);
        c.closz(self.ninwt);
        c.closz(naed);
        c.closz(self.nfscr);
        c.closz(-self.nflx);
        time = c.timer();
        c.io.so(&format!("{}{}s\n {}*******", " ".repeat(69), ffmt::f(time, 8, 1), "**********".repeat(7)));
    }

    /// The output MF for the current mfd.
    fn mfh_out(&self) -> i32 {
        let mut mfh = self.mfd;
        if self.mfd == 8 {
            mfh = 6;
        }
        if self.mfd == 17 {
            mfh = 16;
        }
        if self.mfd == 18 {
            mfh = 16;
        }
        if self.mfd >= 31 && self.mfd <= 36 {
            mfh = self.mfd - 10;
        }
        if self.mfd >= 10000000 {
            mfh = 3;
        }
        mfh
    }

    /// Read user input for group averaging.
    fn ruinb(&mut self, c: &mut Ctx, iaddmt: &mut i32) {
        let sigzmx = 1.0e10f64;
        let zero = 0.0f64;
        //--read and display user input
        self.ngout1 = 0;
        self.ngout2 = 0;
        let v = c.io.input.read(&[Kind::Int, Kind::Int, Kind::Int, Kind::Int]);
        self.nendf = geti(&v[0], self.nendf);
        self.npend = geti(&v[1], self.npend);
        self.ngout1 = geti(&v[2], self.ngout1);
        self.ngout2 = geti(&v[3], self.ngout2);
        c.openz(self.nendf, 0);
        c.openz(self.npend, 0);
        c.openz(self.ngout1, 0);
        c.openz(self.ngout2, 1);
        self.iprint = 1;
        self.ismooth = 1;
        self.ntemp = 1;
        self.nsigz = 1;
        let v = c.io.input.read(&[Kind::Int; 9]);
        self.matb = geti(&v[0], self.matb);
        self.ign = geti(&v[1], self.ign);
        self.igg = geti(&v[2], self.igg);
        self.iwt = geti(&v[3], self.iwt);
        self.lord = geti(&v[4], self.lord);
        self.ntemp = geti(&v[5], self.ntemp);
        self.nsigz = geti(&v[6], self.nsigz);
        self.iprint = geti(&v[7], self.iprint);
        self.ismooth = geti(&v[8], self.ismooth);
        *iaddmt = 0;
        if self.matb < 0 && self.ngout1 != 0 {
            *iaddmt = 1;
        }
        if self.matb < 0 && self.ngout1 == 0 {
            *iaddmt = -1;
        }
        self.ntw = 17;
        let v = c.io.input.read(&[Kind::Char(66)]);
        let text = getc(&v[0], 66, &" ".repeat(66));
        // read(text,'(16a4,a2)') title; title is equivalenced with rtitle
        let mut t68: Vec<u8> = text.as_bytes()[..66].to_vec();
        t68.push(b' ');
        t68.push(b' ');
        pack_text(&t68, &mut self.rtitle);
        c.io.so(&format!(
            "\n unit for endf tape ................... {}\n unit for pendf tape .................. {}\n unit for input gout tape ............. {}\n unit for output gout tape ............ {}",
            ffmt::i(self.nendf, 10),
            ffmt::i(self.npend, 10),
            ffmt::i(self.ngout1, 10),
            ffmt::i(self.ngout2, 10)
        ));
        c.io.so(&format!(
            " mat to be processed .................. {}\n neutron group structure option ....... {}\n gamma group option ................... {}\n weight function option ............... {}\n legendre order ....................... {}\n print option (0 min, 1 max) .......... {}\n smoothing option (0 off, 1 on) ....... {}",
            ffmt::i(self.matb, 10),
            ffmt::i(self.ign, 10),
            ffmt::i(self.igg, 10),
            ffmt::i(self.iwt, 10),
            ffmt::i(self.lord, 10),
            ffmt::i(self.iprint, 10),
            ffmt::i(self.ismooth, 10)
        ));
        if self.ismooth != 0 && self.ismooth != 1 {
            c.error("ruinb", "illegal ismooth.", " ");
        }
        let tt = unpack_text(&self.rtitle);
        c.io.so(&format!("\n run title\n {}\n      {}", "-".repeat(38), String::from_utf8_lossy(&tt[..66])));
        c.io.so(" ");
        self.temp = alloc(self.ntemp);
        self.sigz = alloc(self.nsigz);
        let v = c.io.input.read(&vec![Kind::Real; self.ntemp.max(0) as usize]);
        for i in 1..=self.ntemp {
            self.temp[i] = getr(&v[(i - 1) as usize], self.temp[i]);
        }
        if self.temp[1] == zero {
            c.io.so(" temperatures (kelvin) ................       zero");
        }
        if self.temp[1] != zero {
            c.io.so(&format!(" temperatures (kelvin) ................ {}", ffmt::e1p(self.temp[1], 10, 2)));
        }
        for i in 2..=self.ntemp {
            c.io.so(&format!("{}{}", " ".repeat(40), ffmt::e1p(self.temp[i], 10, 2)));
        }
        let v = c.io.input.read(&vec![Kind::Real; self.nsigz.max(0) as usize]);
        for i in 1..=self.nsigz {
            self.sigz[i] = getr(&v[(i - 1) as usize], self.sigz[i]);
        }
        self.sigz[1] = sigzmx;
        if self.nsigz > 0 {
            c.io.so(" sigma zeroes .........................   infinity");
        }
        for i in 2..=self.nsigz {
            c.io.so(&format!("{}{}", " ".repeat(40), ffmt::e1p(self.sigz[i], 10, 2)));
        }
        if *iaddmt > 0 {
            c.io.so("\n replacing and/or adding mts");
        }
        //--exclude ign=-1 (only used in errorr)
        if self.ign == -1 {
            c.error("gengpn", "illegal group structure.", " ");
        }
        let (mut ngn, mut egn) = (self.ngn, std::mem::take(&mut self.egn));
        gengpn(c, self.ign, &mut ngn, &mut egn);
        self.ngn = ngn;
        self.egn = egn;
        let (mut ngg, mut egg) = (self.ngg, std::mem::take(&mut self.egg));
        gengpg(c, self.igg, &mut ngg, &mut egg);
        self.ngg = ngg;
        self.egg = egg;
        self.genwtf(c);
    }

    /// Find the next reaction in this file of an ENDF or PENDF tape. If
    /// mtd=0, the first reaction is found. Returns iauto=0 at file end.
    /// The `mfd` dummy is the module mfd (the only actual argument).
    fn nextr(&mut self, c: &mut Ctx, iauto: &mut i32, matd: i32, mtd: &mut i32, scr: &mut RV) {
        let (mut nb, mut nw) = (0i32, 0i32);
        //--for mf3, process the reactions on pendf for this temp.
        //--exclude thermal reactions
        if self.mfd <= 3 {
            let nin = self.npend;
            if *mtd <= 0 {
                let mft = self.mfd;
                c.findf(matd, mft, *mtd, nin);
            }
            let mut idone = 0;
            while idone == 0 {
                c.contio(nin, 0, 0, &mut scr.at(1), &mut nb, &mut nw);
                if c.mfh <= 0 || c.mth != 0 {
                    *mtd = c.mth;
                    if *mtd <= 200 {
                        idone = 1;
                    }
                    if *mtd >= 203 && *mtd <= 207 {
                        idone = 1;
                    }
                    if *mtd > 300 {
                        idone = 1;
                    }
                    if idone == 0 {
                        c.tosend(nin, 0, 0, &mut scr.at(1));
                    }
                }
            }
            if *mtd <= 0 {
                *iauto = 0;
            }
            c.skiprz(nin, -1);
            return;
        }

        //--for other reactions, use lists from conver.
        if *mtd == 0 {
            self.nx_ir = 0;
        }
        self.nx_ir += 1;
        let mut lbl = 0;
        // neutrons
        if self.mfd == 6 {
            *mtd = self.mf4[self.nx_ir as usize];
            if *mtd != 0 && (*mtd).abs() <= 91 {
                return;
            }
            self.nx_ir = 1;
            self.mfd = 8;
        }
        // 220
        if self.mfd == 8 {
            *mtd = self.mf6[self.nx_ir as usize];
            if *mtd != 0 {
                return;
            }
            lbl = 280;
        }
        if lbl == 0 {
            // 225: nuclide production (first time)
            if self.mfd == 10 {
                self.nx_ir = 1;
                if self.mf10f[1] == 0 {
                    lbl = 280;
                } else {
                    self.next_prod();
                    *mtd = self.mf10s[self.nx_ir as usize];
                    return;
                }
            }
        }
        if lbl == 0 {
            // 230: photons
            if self.mfd == 16 {
                *mtd = self.mf12[self.nx_ir as usize];
                if *mtd != 0 {
                    return;
                }
                self.nx_ir = 1;
                self.mfd = 17;
            }
            if self.mfd == 17 {
                *mtd = self.mf13[self.nx_ir as usize];
                if *mtd != 0 {
                    return;
                }
                self.nx_ir = 1;
                self.mfd = 18;
            }
            if self.mfd == 18 {
                *mtd = self.mf18[self.nx_ir as usize];
                if *mtd != 0 {
                    return;
                }
                lbl = 280;
            }
        }
        if lbl == 0 {
            // charged particles: mf6 lists then mf4 recoil lists
            for k in 1..=6 {
                if self.mfd == 20 + k {
                    *mtd = self.mf6p.g(k, self.nx_ir);
                    if *mtd != 0 {
                        return;
                    }
                    self.nx_ir = 1;
                    self.mfd = 30 + k;
                }
                if self.mfd == 30 + k {
                    *mtd = self.mf4r.g(k, self.nx_ir);
                    if *mtd != 0 {
                        return;
                    }
                    lbl = 280;
                    break;
                }
            }
        }
        if lbl == 0 {
            // 430: production (after first 10/ entry)
            if self.mfd >= 10000000 && self.mf10f[self.nx_ir as usize] != 0 {
                self.next_prod();
                *mtd = self.mf10s[self.nx_ir as usize];
                return;
            }
        }
        //--return next mt, or zero when finished.
        // 280
        *mtd = 0;
        *iauto = 0;
        if self.mfd == 8 {
            self.mfd = 6;
        }
        if self.mfd == 17 || self.mfd == 18 {
            self.mfd = 16;
        }
        if self.mfd >= 31 && self.mfd <= 36 {
            self.mfd -= 10;
        }
    }

    /// nextr: set mfd, lfs and isom for production entry ir.
    fn next_prod(&mut self) {
        let ir = self.nx_ir as usize;
        self.mfd = self.mf10i[ir];
        self.lfs = self.lfs8[ir];
        self.isom = self.mlfs8[ir];
        if self.mf10f[ir] == 3 {
            self.mfd += 10000000;
        }
        if self.mf10f[ir] == 6 {
            self.mfd += 20000000;
        }
        if self.mf10f[ir] == 9 {
            self.mfd += 30000000;
        }
        if self.mf10f[ir] == 10 {
            self.mfd += 40000000;
        }
    }

    /// Check for 12 required when 13 is given, or for 17 required when
    /// 16 is given.
    fn mfchk(&self, mf: &mut i32, mt: i32) {
        let max = MAXR1 as i32;
        let list = if *mf == 13 {
            &self.mf12
        } else if *mf == 16 {
            &self.mf13
        } else {
            return;
        };
        let newmf = if *mf == 13 { 12 } else { 17 };
        let mut i = 1;
        while i <= max {
            let li = list[i as usize];
            if li == 0 {
                return;
            }
            if li == mt {
                *mf = newmf;
                return;
            }
            if i > 1 {
                if li <= 0 && list[(i - 1) as usize] <= mt && -li >= mt {
                    *mf = newmf;
                    return;
                }
            }
            i += 1;
        }
    }

    /// Check for charged particle production from mf4 recoils when
    /// mf=21-26 is requested.
    fn mfchk2(&self, mf: &mut i32, mt: i32) {
        let max = MAXR1 as i32;
        let mut i = 1;
        while i <= max {
            let mtd = self.mf4r.g(*mf - 20, i);
            if mtd == 0 {
                return;
            }
            if mt == mtd {
                *mf += 10;
                return;
            }
            i += 1;
        }
    }
}

/// Standard names for ENDF reactions (`namer`), as 60 characters. The
/// letter 'h' is used for He-3.
fn namer(c: &Ctx, izad: i32, izar: i32, lfs: i32, mfd: i32, mtd: i32) -> String {
    const IP: [i32; 7] = [0, 1, 1001, 1002, 1003, 2003, 2004];
    const NP: [&str; 7] = ["g", "n", "p", "d", "t", "h", "a"];
    const I2: [i32; 17] = [16, 17, 18, 6, 8, 21, 22, 23, 24, 25, 26, 31, 32, 33, 34, 35, 36];
    const N2: [&str; 17] = [
        "gamma   ", "gamma   ", "gamma   ", "neutron ", "neutron ", "proton  ", "deuteron", "triton  ", "he3     ",
        "alpha   ", "recoil  ", "proton  ", "deuteron", "triton  ", "he3     ", "alpha   ", "recoil  ",
    ];
    const IR: [i32; 113] = [
        11, 16, 17, 22, 23, 24, 25, 28, 29, 30, 32, 33, 34, 35, 36, 37, 41, 42, 44, 45, 108, 109, 111, 112, 113, 114,
        115, 116, 117, 152, 153, 154, 155, 156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170,
        171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181, 182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192,
        193, 194, 195, 196, 197, 198, 199, 200, 18, 19, 20, 21, 38, 102, 103, 104, 105, 106, 107, 201, 202, 203, 204,
        205, 206, 207, 501, 502, 504, 515, 516, 517, 1, 2, 3, 4, 5, 10, 101, 301, 443, 444, 599,
    ];
    const NR: [&str; 113] = [
        "2nd    ", "2n     ", "3n     ", "na     ", "n3a    ", "2na    ", "3na    ", "np     ", "n2a    ", "2n2a   ",
        "nd     ", "nt     ", "nh     ", "nd2a   ", "nt2a   ", "4n     ", "2np    ", "3np    ", "n2p    ", "npa    ",
        "2a     ", "3a     ", "2p     ", "pa     ", "t2a    ", "d2a    ", "pd     ", "pt     ", "da     ", "5n     ",
        "6n     ", "2nt    ", "ta     ", "4np    ", "3nd    ", "nda    ", "2npa   ", "7n     ", "8n     ", "5np    ",
        "6np    ", "7np    ", "4na    ", "5na    ", "6na    ", "7na    ", "4nd    ", "5nd    ", "6nd    ", "3nt    ",
        "4nt    ", "5nt    ", "6nt    ", "2nh    ", "3nh    ", "4nh    ", "3n2p   ", "3n2a   ", "3npa   ", "dt     ",
        "npd    ", "npt    ", "ndt    ", "nph    ", "ndh    ", "nth    ", "nta    ", "2n2p   ", "ph     ", "dh     ",
        "ha     ", "4n2p   ", "4n2a   ", "4npa   ", "3p     ", "n3p    ", "3n2pa  ", "5n2p   ", "fission", "f      ",
        "nf     ", "2nf    ", "3nf    ", "g      ", "p      ", "d      ", "t      ", "h      ", "a      ", "xn     ",
        "xg     ", "xp     ", "xd     ", "xt     ", "xh     ", "xa     ", "total  ", "coh    ", "incoh  ", "paire  ",
        "pair   ", "pairn  ", "total  ", "elastic", "nonel  ", "inel   ", "x      ", "cont   ", "dis    ", "heat   ",
        "kerma  ", "damage ", "stop   ",
    ];

    //--get projectile name
    let mut proj = "?";
    for i in 0..7 {
        if izad == IP[i] {
            proj = NP[i];
        }
    }

    //--get reaction name
    let mut reac: String = "?".to_string();
    let f1 = |p: &str, n: i32| format!("{}{}", p, ffmt::i(n, 1));
    let f2 = |p: &str, n: i32| format!("{}{}", p, ffmt::i(n, 2));
    let m = mtd;
    if c.iverf < 6 {
        if m >= 51 && m < 60 {
            reac = f1("n0", m - 50);
        } else if m >= 60 && m <= 90 {
            reac = f2("n", m - 50);
        } else if m == 91 {
            reac = "nc".into();
        } else if m >= 700 && m <= 709 {
            reac = f1("p0", m - 700);
        } else if m >= 710 && m <= 718 {
            reac = f2("p", m - 700);
        } else if m == 719 {
            reac = "pc".into();
        } else if m >= 720 && m <= 729 {
            reac = f1("d0", m - 720);
        } else if m >= 730 && m <= 738 {
            reac = f2("d", m - 720);
        } else if m == 739 {
            reac = "dc".into();
        } else if m >= 740 && m <= 749 {
            reac = f1("t0", m - 740);
        } else if m >= 750 && m <= 758 {
            reac = f2("t", m - 740);
        } else if m == 759 {
            reac = "tc".into();
        } else if m >= 760 && m <= 769 {
            reac = f1("h0", m - 760);
        } else if m >= 770 && m <= 778 {
            reac = f2("h", m - 760);
        } else if m == 779 {
            reac = "hc".into();
        } else if m >= 780 && m <= 789 {
            reac = f1("a0", m - 780);
        } else if m >= 790 && m <= 798 {
            reac = f2("a", m - 780);
        } else if m == 799 {
            reac = "ac".into();
        } else if m == 602 {
            reac = "photo".into();
        }
    } else {
        if m >= 50 && m <= 59 {
            reac = f1("n0", m - 50);
        } else if m >= 60 && m <= 90 {
            reac = f2("n", m - 50);
        } else if m == 91 {
            reac = "nc".into();
        } else if m >= 600 && m <= 609 {
            reac = f1("p0", m - 600);
        } else if m >= 610 && m <= 648 {
            reac = f2("p", m - 600);
        } else if m == 649 {
            reac = "pc".into();
        } else if m >= 650 && m <= 659 {
            reac = f1("d0", m - 650);
        } else if m >= 660 && m <= 698 {
            reac = f2("d", m - 650);
        } else if m == 699 {
            reac = "dc".into();
        } else if m >= 700 && m <= 709 {
            reac = f1("t0", m - 700);
        } else if m >= 710 && m <= 748 {
            reac = f2("t", m - 700);
        } else if m == 749 {
            reac = "tc".into();
        } else if m >= 750 && m <= 759 {
            reac = f1("h0", m - 750);
        } else if m >= 760 && m <= 798 {
            reac = f2("h", m - 750);
        } else if m == 799 {
            reac = "hc".into();
        } else if m >= 800 && m <= 809 {
            reac = f1("a0", m - 800);
        } else if m >= 810 && m <= 848 {
            reac = f2("a", m - 800);
        } else if m == 849 {
            reac = "ac".into();
        } else if m >= 875 && m <= 884 {
            reac = f1("2n0", m - 875);
        } else if m >= 885 && m <= 890 {
            reac = f2("2n", m - 875);
        } else if m == 891 {
            reac = "2nc".into();
        } else if m == 522 {
            reac = "photo".into();
        } else if m == 500 {
            reac = "stop".into();
        }
    }
    let mut reac = pad(&reac, 7);
    if reac.trim_end() == "?" {
        for i in 0..113 {
            if mtd == IR[i] {
                reac = NR[i].to_string();
            }
        }
    }

    //--get outgoing particle name
    let dummy: String = if mfd == 3 {
        format!("({},{})-cross-section.", proj, reac)
    } else if mfd >= 10000000 {
        let azam = if lfs == 0 {
            format!("{}", izar)
        } else if mfd < 30000000 {
            format!("{}m{}", izar, lfs)
        } else {
            format!("{}e{}", izar, lfs)
        };
        let azam = pad(&azam, 10);
        format!("({},{})-{}-production.", proj, reac, azam.trim_end())
    } else {
        let mut part = "?       ";
        for i in 0..17 {
            if mfd == I2[i] {
                part = N2[i];
            }
        }
        format!("({},{})-{}-matrix.", proj, reac, part)
    };
    let mut d: Vec<u8> = pad(&dummy, 60).into_bytes();

    //--strip unneeded spaces
    let mut i = 0usize;
    let mut res: Vec<u8> = vec![b' '; 60];
    let mut j = 0usize;
    while d[i] != b'.' {
        while d[i] == b' ' {
            i += 1;
        }
        j += 1;
        if d[i] == b'-' {
            d[i] = b' ';
        }
        res[j - 1] = d[i];
        i += 1;
    }
    String::from_utf8_lossy(&res).to_string()
}
