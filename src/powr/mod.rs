//! POWR: input for the EPRI-CELL codes GAMTAP (fast, lib=1) and LIBRAR
//! (thermal, lib=2), and for the EPRI-CPM code CLIB (lib=3), ported from
//! NJOY2016 `powr.f90` (module powm).
//!
//! The routines keep the Fortran structure, 1-based indices and the order
//! of every floating-point operation. The module variables of `powm` and
//! the saved locals (therm `valu`, sf4out `matl/ifissl/ninal`, sf6out
//! `jtem`, cpmout `iflgs/eflag`) are the fields of `Powr`, so that they
//! keep their values from one POWR run to the next, as in the Fortran.
//! The scratch array `scr` is one of them: cpm takes its card-5 values
//! from `scr`, not from the variables it reads them into (an NJOY bug,
//! kept), so they are what an earlier POWR run left there (zero in the
//! first run).
//!
//! Character data that the Fortran reads with `a4` into real variables
//! (nuclide descriptions, the `end ` and `pri ` flags) are kept as text
//! where they are only written out again with `a` edits; the `end ` flag
//! that goes through the binary scratch files is the real number with
//! the same bytes.

#![allow(clippy::too_many_lines)]
#![allow(clippy::too_many_arguments)]
#![allow(clippy::collapsible_if)]
#![allow(clippy::collapsible_else_if)]
#![allow(clippy::assign_op_pattern)]
#![allow(clippy::needless_range_loop)]
// the label machines keep the Fortran assignments, some of them dead
#![allow(unused_assignments)]

mod cpm;
mod cpm4;
mod therm;

use crate::endf::Ctx;
use crate::farr::{IV, RV};
use crate::ffmt;
use crate::fio::{fatal_runtime, Kind, Val};
use crate::fmath::{exp, log, nint, sqrt};

/// Extra words past the Fortran array sizes.
pub(crate) const SLACK: usize = 4096;
const ZERO: f64 = 0.0;

pub(crate) fn geti(v: &Option<Val>, d: i32) -> i32 {
    match v {
        Some(Val::Int(i)) => *i,
        Some(Val::Real(x)) => *x as i32,
        _ => d,
    }
}

pub(crate) fn getr(v: &Option<Val>, d: f64) -> f64 {
    match v {
        Some(Val::Real(x)) => *x,
        Some(Val::Int(i)) => *i as f64,
        _ => d,
    }
}

/// A list-directed character value of length n (unchanged if null).
pub(crate) fn getc(v: &Option<Val>, n: usize, d: &str) -> String {
    match v {
        Some(Val::Char(s)) => {
            let mut t: String = s.chars().take(n).collect();
            while t.len() < n {
                t.push(' ');
            }
            t
        }
        _ => d.to_string(),
    }
}

/// `read(s,'(na4)')` into n real(kr) words: four characters each, blank
/// padded past the end of the internal record. Only the four characters
/// are kept (an `a4` edit writes them back).
pub(crate) fn a4s(s: &str, n: usize) -> Vec<String> {
    let b = s.as_bytes();
    (0..n)
        .map(|i| {
            (0..4)
                .map(|j| {
                    let k = 4 * i + j;
                    if k < b.len() { b[k] as char } else { ' ' }
                })
                .collect()
        })
        .collect()
}

/// `1p,ew.d`.
#[inline]
pub(crate) fn pe(x: f64, w: usize, d: usize) -> String {
    ffmt::ep(x, w, d, 1)
}

/// The records of `write(u,fmt) v` for a format of the form `(n<edit>)`
/// with no inner parentheses (with `1p` in front, which is re-applied on
/// reversion): `per` values a record, and one empty record when there
/// are no values.
pub(crate) fn rows(v: &[f64], per: usize, f: impl Fn(f64) -> String) -> Vec<String> {
    let mut out = Vec::new();
    if v.is_empty() {
        out.push(String::new());
        return out;
    }
    for ch in v.chunks(per) {
        let mut s = String::new();
        for &x in ch {
            s.push_str(&f(x));
        }
        out.push(s);
    }
    out
}

/// Same as `rows`, joined for the listing.
pub(crate) fn rowsj(v: &[f64], per: usize, f: impl Fn(f64) -> String) -> String {
    rows(v, per, f).join("\n")
}

pub(crate) fn rowsi(v: &[i32], per: usize, w: usize) -> Vec<String> {
    let mut out = Vec::new();
    if v.is_empty() {
        out.push(String::new());
        return out;
    }
    for ch in v.chunks(per) {
        let mut s = String::new();
        for &x in ch {
            s.push_str(&ffmt::i(x, w));
        }
        out.push(s);
    }
    out
}

/// gfortran's messages for allocate/deallocate errors.
pub(crate) fn alloc_iv(flag: &mut bool, v: &mut IV, n: i32, name: &str) {
    if *flag {
        fatal_runtime(&format!(
            "Fortran runtime error: Attempting to allocate already allocated variable '{}'",
            name
        ));
    }
    *flag = true;
    *v = IV::new(n.max(0) as usize + SLACK);
}

pub(crate) fn alloc_rv(flag: &mut bool, v: &mut RV, n: i32, name: &str) {
    if *flag {
        fatal_runtime(&format!(
            "Fortran runtime error: Attempting to allocate already allocated variable '{}'",
            name
        ));
    }
    *flag = true;
    *v = RV::new(n.max(0) as usize + SLACK);
}

pub(crate) fn dealloc(flag: &mut bool, name: &str) {
    if !*flag {
        fatal_runtime(&format!("Fortran runtime error: Attempt to DEALLOCATE unallocated '{}'", name));
    }
    *flag = false;
}

/// The coarse-group maps of `icgrp` (fine 185-group index, 1..186).
const IGF: [i32; 186] = [
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 1, 1, 2, 2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5,
    6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12,
    12, 13, 13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18, 19, 19, 20, 20, 21, 21, 22,
    22, 23, 23, 24, 24, 25, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30, 31, 31, 32,
    32, 33, 33, 34, 34, 35, 35, 36, 36, 37, 37, 38, 38, 39, 39, 40, 40, 41, 41, 42,
    42, 43, 43, 44, 44, 45, 45, 46, 46, 47, 47, 48, 48, 49, 49, 50, 50, 51, 51, 52,
    52, 53, 53, 54, 54, 55, 55, 56, 56, 57, 57, 58, 58, 59, 59, 60, 60, 61, 61, 62,
    62, 63, 63, 64, 64, 65, 65, 66, 66, 67, 67, 68, 68, 69, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
];
const IGT: [i32; 186] = [
    0, 1, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
    19, 19, 20, 21, 21, 22, 22, 23, 23, 24, 24, 25, 26, 27, 28, 29, 30, 30, 31, 32,
    32, 33, 34, 35, 36, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
];
const IGC: [i32; 186] = [
    1, 2, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
    20, 20, 21, 22, 22, 23, 23, 24, 25, 26, 27, 28, 29, 29, 29, 30, 31, 32, 33, 34,
    34, 34, 35, 35, 36, 36, 36, 37, 37, 37, 37, 38, 39, 39, 40, 41, 41, 42, 42, 42,
    42, 43, 43, 44, 44, 44, 45, 45, 45, 45, 46, 46, 46, 46, 46, 47, 47, 47, 47, 48,
    48, 48, 48, 48, 49, 49, 49, 49, 49, 49, 49, 50, 50, 50, 50, 50, 50, 50, 51, 51,
    51, 51, 51, 52, 52, 52, 53, 53, 53, 54, 54, 54, 54, 55, 55, 55, 55, 56, 56, 56,
    56, 57, 57, 57, 57, 58, 58, 58, 58, 59, 59, 59, 59, 60, 60, 60, 60, 61, 61, 61,
    61, 62, 62, 62, 62, 63, 63, 63, 63, 64, 64, 64, 64, 65, 65, 65, 65, 66, 66, 66,
    66, 67, 67, 67, 67, 68, 68, 68, 68, 69, 69, 69, 69, 70, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0,
];

/// Module variables of `powm` and the saved locals.
pub struct Powr {
    // i/o units
    pub(crate) ngendf: i32,
    pub(crate) nout: i32,
    pub(crate) nscr: i32,
    pub(crate) kscr: i32,
    pub(crate) mscr: i32,
    // global
    pub(crate) lib: i32,
    pub(crate) ngnf: i32,
    pub(crate) ngn: i32,
    pub(crate) ngnd: i32,
    pub(crate) cflux: RV,
    pub(crate) nwscr: i32,
    pub(crate) scr: RV,
    // epri-cell
    pub(crate) tmpr: RV,
    pub(crate) sigz: RV,
    pub(crate) ntp: i32,
    pub(crate) nsigz: i32,
    pub(crate) jwf: i32,
    pub(crate) iwr: i32,
    pub(crate) izref: i32,
    pub(crate) rtemp: f64,
    pub(crate) matd: i32,
    pub(crate) nid: i32,
    pub(crate) nff: i32,
    pub(crate) iwa: i32,
    pub(crate) iwf: i32,
    pub(crate) nfs: i32,
    pub(crate) xla: [f64; 4],
    pub(crate) xld: [f64; 4],
    pub(crate) xlol: [f64; 4],
    // epri-cpm
    pub(crate) nfg: i32,
    pub(crate) nrg: i32,
    pub(crate) nfsg: i32,
    pub(crate) ngref: i32,
    pub(crate) newmat: i32,
    pub(crate) iprint: i32,
    pub(crate) nzmax: i32,
    pub(crate) mode: i32,
    pub(crate) if5: i32,
    pub(crate) ntis: i32,
    pub(crate) nfis: i32,
    pub(crate) if4: i32,
    pub(crate) nnina: i32,
    pub(crate) nres: i32,
    pub(crate) iu235: i32,
    pub(crate) iu: i32,
    pub(crate) ipu239: i32,
    pub(crate) ipu: i32,
    pub(crate) nfiss: [i32; 3],
    pub(crate) rleth: RV,
    pub(crate) reflth: f64,
    // allocatable module arrays (and their allocation status)
    pub(crate) mat: IV,
    pub(crate) al_mat: bool,
    pub(crate) flux: RV,
    pub(crate) al_flux: bool,
    pub(crate) idnt: IV,
    pub(crate) al_idnt: bool,
    pub(crate) ninat: IV,
    pub(crate) al_ninat: bool,
    pub(crate) irest: IV,
    pub(crate) al_irest: bool,
    pub(crate) s0rf: RV,
    pub(crate) al_s0rf: bool,
    pub(crate) ip1tab: IV,
    pub(crate) al_ip1tab: bool,
    pub(crate) norft: IV,
    pub(crate) al_norft: bool,
    pub(crate) mtit: IV,
    pub(crate) al_mtit: bool,
    pub(crate) mtct: IV,
    pub(crate) al_mtct: bool,
    pub(crate) awrt: RV,
    pub(crate) al_awrt: bool,
    pub(crate) indf: IV,
    pub(crate) al_indf: bool,
    pub(crate) ntmp: IV,
    pub(crate) al_ntmp: bool,
    pub(crate) nszt: IV,
    pub(crate) al_nszt: bool,
    pub(crate) glamt: RV,
    pub(crate) al_glamt: bool,
    pub(crate) ipost: IV,
    pub(crate) al_ipost: bool,
    pub(crate) iposrt: IV,
    pub(crate) al_iposrt: bool,
    pub(crate) egb: RV,
    pub(crate) al_egb: bool,
    pub(crate) sigpt: RV,
    pub(crate) al_sigpt: bool,
    pub(crate) uff: RV,
    pub(crate) al_uff: bool,
    pub(crate) puff: RV,
    pub(crate) al_puff: bool,
    pub(crate) fiss: RV,
    pub(crate) al_fiss: bool,
    pub(crate) yld: RV,
    pub(crate) al_yld: bool,
    pub(crate) idat: IV,
    pub(crate) al_idat: bool,
    pub(crate) idbt: IV,
    pub(crate) al_idbt: bool,
    pub(crate) elas: RV,
    pub(crate) al_elas: bool,
    pub(crate) tempt: RV,
    pub(crate) al_tempt: bool,
    // saved locals
    /// therm `valu(20)` (data statement: saved)
    pub(crate) valu: [f64; 21],
    /// sf4out `matl`, `ifissl`, `ninal`
    pub(crate) s4_matl: i32,
    pub(crate) s4_ifissl: i32,
    pub(crate) s4_ninal: i32,
    /// sf6out `jtem`
    pub(crate) s6_jtem: i32,
    /// cpmout `iflgs(80)` (index 0 is the word before it, which a
    /// position of 0 reaches) and `eflag`
    pub(crate) co_iflgs: [i32; 81],
    pub(crate) co_eflag: f64,
}

impl Default for Powr {
    fn default() -> Self {
        Powr {
            ngendf: 0,
            nout: 0,
            nscr: 0,
            kscr: 0,
            mscr: 0,
            lib: 0,
            ngnf: 0,
            ngn: 0,
            ngnd: 0,
            cflux: RV::new(140 + SLACK),
            nwscr: 5000,
            scr: RV::new(5000 + SLACK),
            tmpr: RV::new(10 + 16),
            sigz: RV::new(10 + 16),
            ntp: 0,
            nsigz: 0,
            jwf: 0,
            iwr: 0,
            izref: 0,
            rtemp: 0.0,
            matd: 0,
            nid: 0,
            nff: 0,
            iwa: 0,
            iwf: 0,
            nfs: 0,
            xla: [0.0; 4],
            xld: [0.0; 4],
            xlol: [0.0; 4],
            nfg: 0,
            nrg: 0,
            nfsg: 0,
            ngref: 0,
            newmat: 0,
            iprint: 0,
            nzmax: 0,
            mode: 0,
            if5: 0,
            ntis: 0,
            nfis: 0,
            if4: 0,
            nnina: 0,
            nres: 0,
            iu235: 0,
            iu: 0,
            ipu239: 0,
            ipu: 0,
            nfiss: [0; 3],
            rleth: RV::new(19 + 16),
            reflth: 0.0,
            mat: IV::default(),
            al_mat: false,
            flux: RV::default(),
            al_flux: false,
            idnt: IV::default(),
            al_idnt: false,
            ninat: IV::default(),
            al_ninat: false,
            irest: IV::default(),
            al_irest: false,
            s0rf: RV::default(),
            al_s0rf: false,
            ip1tab: IV::default(),
            al_ip1tab: false,
            norft: IV::default(),
            al_norft: false,
            mtit: IV::default(),
            al_mtit: false,
            mtct: IV::default(),
            al_mtct: false,
            awrt: RV::default(),
            al_awrt: false,
            indf: IV::default(),
            al_indf: false,
            ntmp: IV::default(),
            al_ntmp: false,
            nszt: IV::default(),
            al_nszt: false,
            glamt: RV::default(),
            al_glamt: false,
            ipost: IV::default(),
            al_ipost: false,
            iposrt: IV::default(),
            al_iposrt: false,
            egb: RV::default(),
            al_egb: false,
            sigpt: RV::default(),
            al_sigpt: false,
            uff: RV::default(),
            al_uff: false,
            puff: RV::default(),
            al_puff: false,
            fiss: RV::default(),
            al_fiss: false,
            yld: RV::default(),
            al_yld: false,
            idat: IV::default(),
            al_idat: false,
            idbt: IV::default(),
            al_idbt: false,
            elas: RV::default(),
            al_elas: false,
            tempt: RV::default(),
            al_tempt: false,
            valu: [0.0; 21],
            s4_matl: 0,
            s4_ifissl: 0,
            s4_ninal: 0,
            s6_jtem: 0,
            co_iflgs: [0; 81],
            co_eflag: 1.0e35,
        }
    }
}

impl Powr {
    pub fn new() -> Self {
        Self::default()
    }

    /// Produce input for the EPRI-CELL codes GAMTAP (fast) and LIBRAR
    /// (thermal), and the EPRI-CPM code CLIB.
    pub fn powr(&mut self, c: &mut Ctx) {
        self.ngnf = 185;
        let mut igprnt = 0;
        let sec = c.timer();
        c.io.so(&format!(
            "\n powr...produce input for epri-cell libraries{}{}s",
            " ".repeat(23),
            ffmt::f(sec, 8, 1)
        ));
        c.io.se(&format!("\n powr...{}{}s", " ".repeat(61), ffmt::f(sec, 8, 1)));

        //--read and write common input.
        let v = c.io.input.read(&[Kind::Int, Kind::Int]);
        self.ngendf = geti(&v[0], self.ngendf);
        self.nout = geti(&v[1], self.nout);
        let mut iprnt = 0;
        let mut iclaps = 0;
        let v = c.io.input.read(&[Kind::Int, Kind::Int, Kind::Int]);
        self.lib = geti(&v[0], self.lib);
        iprnt = geti(&v[1], iprnt);
        iclaps = geti(&v[2], iclaps);
        if self.nout < 0 {
            c.mess("powr", "nout must be positive.  changed.", " ");
        }
        if self.nout < 0 {
            self.nout = -self.nout;
        }
        c.openz(self.ngendf, 0);
        c.openz(self.nout, 1);
        c.io.so(&format!(
            "\n input gendf unit ..................... {}\n \
             output unit .......................... {}\n \
             library (1 fast, 2 thermal, 3 cpm) ... {}\n \
             print option (0 min, 1 max) .......... {}\n \
             grp collapsing option (0 yes, 1 no) .. {}",
            ffmt::i(self.ngendf, 10),
            ffmt::i(self.nout, 10),
            ffmt::i(self.lib, 10),
            ffmt::i(iprnt, 10),
            ffmt::i(iclaps, 10)
        ));
        c.repoz(self.nout);

        //--branch for desired library
        if self.lib == 1 {
            self.fast(c, iprnt, &mut igprnt);
        } else if self.lib == 2 {
            self.therm(c, iprnt, &mut igprnt);
        } else if self.lib == 3 {
            self.cpm(c, iprnt, iclaps);
        } else {
            c.error("powr", "illegal library requested.", " ");
        }

        c.repoz(self.ngendf);
        c.repoz(self.nout);
        c.closz(self.ngendf);
        c.closz(self.nout);
        let sec = c.timer();
        c.io.so(&format!("{}{}s\n {}*******", " ".repeat(69), ffmt::f(sec, 8, 1), "*".repeat(70)));
    }

    /// For the group ig in the 185-group fine structure, return the coarse
    /// group for the fast, thermal, or CPM group structure.
    pub(crate) fn icgrp(&self, ig: i32) -> i32 {
        if self.ngn == self.ngnd {
            return ig;
        }
        if ig == 0 {
            return 0;
        }
        let k = (ig - 1) as usize;
        if self.lib == 1 {
            IGF[k]
        } else if self.lib == 2 {
            IGT[k]
        } else if self.lib == 3 {
            IGC[k]
        } else {
            // icgrp is undefined here in the Fortran; not reached
            0
        }
    }

    /// Format multigroup cross sections from groupr for the EPRI-CELL
    /// fast library maintenance code GAMTAP.
    fn fast(&mut self, c: &mut Ctx, iprint: i32, igprnt: &mut i32) {
        let misc = 1;
        let izero = 0;
        let na = 10000;
        let mut a = RV::new(na as usize + SLACK);
        let mut z = RV::new(70 + SLACK);
        let mut word: Vec<String> = vec!["    ".to_string(); 4];
        let mut fsn: Vec<String> = vec!["    ".to_string(); 10];
        // loc(1..10): locab0 locsf0 locchi locnus loce0 loce1 locin loc2n
        // locabs locsf (an equivalence in the Fortran)
        let mut loc = [0i32; 11];
        let (mut xld1, mut xld2, mut xld3) = (0.0f64, 0.0f64, 0.0f64);

        // set up i/o units
        self.nscr = 10;
        self.kscr = 11;
        c.openz(self.nscr, 1);
        c.openz(self.kscr, 1);
        c.repoz(self.nscr);
        c.repoz(self.kscr);
        let nscr = self.nscr;
        let kscr = self.kscr;

        // extract the fast data
        let mut ifiss = 0;
        self.ngnd = 68;
        let ngnd = self.ngnd;
        let ngmin = 19;
        let ngmax = 62;
        c.io.so("\n ***epri-cell fast library***");
        self.nff = 1;
        self.matd = 1;
        while self.matd > 0 {
            self.rtemp = 300.0;
            let mut iff = 1;
            let mut nsgz = 0;
            self.izref = 1;
            let v = c.io.input.read(&[Kind::Int, Kind::Real, Kind::Int, Kind::Int, Kind::Int]);
            self.matd = geti(&v[0], self.matd);
            self.rtemp = getr(&v[1], self.rtemp);
            iff = geti(&v[2], iff);
            nsgz = geti(&v[3], nsgz);
            self.izref = geti(&v[4], self.izref);
            if self.matd > 0 {
                let mut iread = 0;
                if self.matd < 0 {
                    iread = 1;
                }
                self.matd = self.matd.abs();
                let v = c.io.input.read(&[Kind::Char(12)]);
                let hid = getc(&v[0], 12, &" ".repeat(12));
                word = a4s(&hid, 4);
                if iread != 1 {
                    let v = c.io.input.read(&[Kind::Char(40)]);
                    let fsp = getc(&v[0], 40, &" ".repeat(40));
                    fsn = a4s(&fsp, 10);
                } else {
                    self.nid = self.matd;
                    self.iwa = 1;
                    self.iwr = 0;
                    self.iwf = 0;
                    for i in 1..=3 {
                        self.xlol[i] = 0.0;
                        self.xla[i] = 0.0;
                    }
                    xld1 = 0.0;
                    xld2 = 0.0;
                    xld3 = 0.0;
                    for i in 1..=ngnd {
                        z[i] = 0.0;
                    }
                    let kinds = vec![Kind::Real; ngnd as usize];
                    let v = c.io.input.read(&kinds);
                    for i in 1..=ngnd {
                        z[i] = getr(&v[(i - 1) as usize], z[i]);
                    }
                    c.io.so(&format!(" mat to be processed .................. {}", ffmt::i(self.matd, 10)));
                    c.io.so(&format!(
                        "\n description of nuclide\n {}--\n      {}",
                        "----".repeat(9),
                        word.concat()
                    ));
                    c.io.so("\n absorption only for this material.  nid=matd.");
                }

                //--write out input parameters.
                if iread != 1 {
                    c.io.so(&format!(" mat to be processed .................. {}", ffmt::i(self.matd, 10)));
                    c.io.so(&format!(
                        " reference temperature ................ {}\n \
                         f-factor option (0 no, 1 yes) ........ {}",
                        pe(self.rtemp, 10, 3),
                        ffmt::i(iff, 10)
                    ));
                    if nsgz > 0 {
                        c.io.so(&format!(" max. no. of sigma zeroes ............. {}", ffmt::i(nsgz, 10)));
                    }
                    if nsgz == 0 {
                        c.io.so(&format!(" max. no. of sigma zeroes ............. {}all", " ".repeat(7)));
                    }
                    c.io.so(&format!(" index of reference sigma zero ........ {}", ffmt::i(self.izref, 10)));
                    c.io.so(&format!(
                        "\n description of nuclide\n {}--\n      {}",
                        "----".repeat(9),
                        word.concat()
                    ));

                    //--set up pointers for gam.
                    self.gamll(c, igprnt, &mut iff, nsgz);

                    //--format cross sections for gam.
                    self.gamxs(c, &mut loc, &mut a, na);

                    //--format shielding factors for gam.
                    if self.nff > 0 && self.iwr > 0 {
                        self.gamff(c, ngmin, ngmax, &mut loc, &mut a, na);
                    }

                    //--write nscr and print output, if desired.
                    xld1 = self.xld[1] + 1.0;
                    xld2 = self.xld[2] + 1.0;
                    xld3 = self.xld[3] + 1.0;
                    if self.nfs + self.iwf != 0 {
                        ifiss = self.iwf + self.nfs;
                        c.io.write_line(
                            kscr,
                            &format!("{}{}{}", ffmt::i(self.nid, 6), ffmt::i(izero, 2), fsn.concat()),
                        );
                        if ifiss == 1 {
                            c.io.so("\n fission spectrum\n ----------------");
                        }
                        if ifiss > 1 {
                            c.io.so("\n fission spectra\n ---------------");
                        }
                        c.io.so(&format!("      {}total\n", fsn.concat()));
                        let mut l = loc[3];
                        let mut lim = l + ngnd - 1;
                        for s in rows(&a.0[(l - 1) as usize..lim.max(l - 1) as usize], 6, |x| pe(x, 12, 5)) {
                            c.io.write_line(kscr, &s);
                        }
                        c.io.so(&rowsj(&a.0[(l - 1) as usize..lim.max(l - 1) as usize], 6, |x| pe(x, 12, 4)));
                        if self.nfs != 0 {
                            for ifs in 1..=self.nfs {
                                c.io.write_line(
                                    kscr,
                                    &format!("{}{}{}", ffmt::i(self.nid, 6), ffmt::i(ifs, 2), fsn.concat()),
                                );
                                c.io.so(&format!("\n      {}delayed{}\n", fsn.concat(), ffmt::i(ifs, 1)));
                                l += ngnd;
                                lim += ngnd;
                                let sl = &a.0[(l - 1) as usize..lim.max(l - 1) as usize];
                                c.io.so(&rowsj(sl, 6, |x| pe(x, 12, 4)));
                                for s in rows(sl, 6, |x| pe(x, 12, 5)) {
                                    c.io.write_line(kscr, &s);
                                }
                            }
                            self.nfs += self.iwf;
                        }
                    }
                }

                //--write more nout and print output
                c.io.write_line(
                    nscr,
                    &format!(
                        "{}{}{}{}",
                        ffmt::i(self.nid, 10),
                        ffmt::i(self.iwa, 10),
                        ffmt::i(self.iwf, 10),
                        ffmt::i(self.iwr, 10)
                    ),
                );
                c.io.write_line(nscr, &word.concat());
                c.io.write_line(
                    nscr,
                    &format!(
                        "{}{}{}{}",
                        " ".repeat(36),
                        pe(self.xlol[3], 12, 5),
                        pe(self.xla[3], 12, 5),
                        pe(xld3, 12, 5)
                    ),
                );
                let six = [self.xlol[1], self.xla[1], xld1, self.xlol[2], self.xla[2], xld2];
                c.io.write_line(nscr, &rows(&six, 6, |x| pe(x, 12, 5))[0]);
                let sl = |a: &RV, i0: i32, i1: i32| -> Vec<f64> {
                    if i1 < i0 { Vec::new() } else { a.0[(i0 - 1) as usize..i1 as usize].to_vec() }
                };
                if self.iwa != 0 {
                    let lim = loc[1] + ngnd - 1;
                    if iread == 0 {
                        for s in rows(&sl(&a, loc[1], lim), 6, |x| pe(x, 12, 5)) {
                            c.io.write_line(nscr, &s);
                        }
                    }
                    if iread > 0 {
                        for s in rows(&sl(&z, 1, ngnd), 6, |x| pe(x, 12, 5)) {
                            c.io.write_line(nscr, &s);
                        }
                    }
                    if iprint != 0 {
                        c.io.so("\n absorption cross section\n ------------------------");
                        if iread == 0 {
                            c.io.so(&rowsj(&sl(&a, loc[1], lim), 6, |x| pe(x, 12, 4)));
                        }
                        if iread > 0 {
                            c.io.so(&rowsj(&sl(&z, 1, ngnd), 6, |x| pe(x, 12, 4)));
                        }
                    }
                }
                if self.iwf != 0 {
                    let limsf = loc[2] + ngnd - 1;
                    for s in rows(&sl(&a, loc[2], limsf), 6, |x| pe(x, 12, 5)) {
                        c.io.write_line(nscr, &s);
                    }
                    let lim = loc[4] + ngnd - 1;
                    for s in rows(&sl(&a, loc[4], lim), 6, |x| pe(x, 12, 5)) {
                        c.io.write_line(nscr, &s);
                    }
                    if iprint != 0 {
                        c.io.so("\n fission cross section\n ---------------------");
                        c.io.so(&rowsj(&sl(&a, loc[2], limsf), 6, |x| pe(x, 12, 4)));
                        c.io.so("\n fission nu\n ----------");
                        c.io.so(&rowsj(&sl(&a, loc[4], lim), 6, |x| pe(x, 12, 4)));
                        if self.nfs != 0 {
                            let limc = loc[3] - 1 + ngnd;
                            c.io.so("\n fission chi\n -----------");
                            c.io.so(&rowsj(&sl(&a, loc[3], limc), 6, |x| pe(x, 12, 4)));
                        }
                    }
                }
                if self.xlol[3] != ZERO {
                    let lim = loc[5] + nint(self.xlol[3]) - 1;
                    for s in rows(&sl(&a, loc[5], lim), 6, |x| pe(x, 12, 5)) {
                        c.io.write_line(nscr, &s);
                    }
                    let lim = loc[6] + nint(self.xlol[3]) - 1;
                    for s in rows(&sl(&a, loc[6], lim), 6, |x| pe(x, 12, 5)) {
                        c.io.write_line(nscr, &s);
                    }
                    if iprint != 0 {
                        let la = nint(self.xla[3]);
                        let ld = nint(self.xld[3]);
                        c.io.so("\n\n p0 elastic scattering matrix\n ----------------------------");
                        pgam(c, &a.0[(loc[5] - 1) as usize..], ngnd, la, ld);
                        c.io.so("\n\n p1 elastic scattering matrix\n ----------------------------");
                        pgam(c, &a.0[(loc[6] - 1) as usize..], ngnd, la, ld);
                    }
                }
                if self.xlol[1] != ZERO {
                    let lim = loc[7] + nint(self.xlol[1]) - 1;
                    for s in rows(&sl(&a, loc[7], lim), 6, |x| pe(x, 12, 5)) {
                        c.io.write_line(nscr, &s);
                    }
                    if iprint != 0 {
                        let la = nint(self.xla[1]);
                        let ld = nint(self.xld[1]);
                        c.io.so("\n\n inelastic scattering matrix\n ---------------------------");
                        pgam(c, &a.0[(loc[7] - 1) as usize..], ngnd, la, ld);
                    }
                }
                if self.xlol[2] != ZERO {
                    let lim = loc[8] + nint(self.xlol[2]) - 1;
                    for s in rows(&sl(&a, loc[8], lim), 6, |x| pe(x, 12, 5)) {
                        c.io.write_line(nscr, &s);
                    }
                    if iprint != 0 {
                        let la = nint(self.xla[2]);
                        let ld = nint(self.xld[2]);
                        c.io.so("\n\n n2n scattering matrix\n ---------------------");
                        pgam(c, &a.0[(loc[8] - 1) as usize..], ngnd, la, ld);
                    }
                }
                if self.nff != 0 && self.iwr != 0 {
                    let nsigz = self.nsigz;
                    let ntp = self.ntp;
                    c.io.write_line(
                        nscr,
                        &format!(
                            "{}{}{}{}",
                            ffmt::i(nsigz, 6),
                            ffmt::i(ntp, 6),
                            ffmt::i(misc, 6),
                            ffmt::i(self.jwf, 6)
                        ),
                    );
                    for s in rows(&sl(&self.sigz, 1, nsigz), 6, |x| pe(x, 12, 5)) {
                        c.io.write_line(nscr, &s);
                    }
                    for s in rows(&sl(&self.tmpr, 1, ntp), 6, |x| pe(x, 12, 5)) {
                        c.io.write_line(nscr, &s);
                    }
                    let istart = loc[9] + (ngmax - self.iwr) * nsigz * ntp;
                    let lim = istart + self.iwr * nsigz * ntp - 1;
                    for s in rows(&sl(&a, istart, lim), 6, |x| pe(x, 12, 5)) {
                        c.io.write_line(nscr, &s);
                    }
                    if iprint != 0 {
                        c.io.so(
                            "\n\n absorption self-shielding factors\n ---------------------------------",
                        );
                        pgamff(c, &a.0[(loc[9] - 1) as usize..], ntp, nsigz, self.iwr, &self.sigz, &self.tmpr, ngmax);
                    }
                    if self.jwf != 0 {
                        let istart = loc[10] + (ngmax - self.iwr) * nsigz * ntp;
                        let lim = istart + self.iwr * nsigz * ntp - 1;
                        for s in rows(&sl(&a, istart, lim), 6, |x| pe(x, 12, 5)) {
                            c.io.write_line(nscr, &s);
                        }
                        if iprint != 0 {
                            c.io.so("\n\n fission self-shielding factors\n ------------------------------");
                            pgamff(
                                c,
                                &a.0[(loc[10] - 1) as usize..],
                                ntp,
                                nsigz,
                                self.iwr,
                                &self.sigz,
                                &self.tmpr,
                                ngmax,
                            );
                        }
                    }
                    // write a dummy amisc array
                    let nwmisc = self.iwr * 10;
                    let zs = vec![ZERO; nwmisc.max(0) as usize];
                    for s in rows(&zs, 6, |x| pe(x, 12, 5)) {
                        c.io.write_line(nscr, &s);
                    }
                }
            }
        }

        //--fast is finished.
        c.repoz(nscr);
        mergt(c, ifiss, kscr, nscr, self.nout);
        c.closz(nscr);
        c.closz(nscr);
    }

    /// Read through gendf tape and assign storage and pointers for GAM
    /// format.
    fn gamll(&mut self, c: &mut Ctx, igprnt: &mut i32, iff: &mut i32, nsgz: i32) {
        let (mut nb, mut nw) = (0i32, 0i32);
        let eps = 1.0e-4f64;
        let ntpmax = 5;
        let nszmax = 7;
        let ng = self.ngendf;
        let ngnd = self.ngnd;
        let (mut ldin, mut ldp, mut ld2n) = (0i32, 0i32, 0i32);
        let (mut lain, mut la2n, mut lap) = (0i32, 0i32, 0i32);
        let lz = 6;
        let mut ntw = 0;
        let (mut ig, mut igc, mut ig2loc) = (0i32, 0i32, 0i32);

        //--read data for desired material.
        c.repoz(ng);
        c.tpidio(ng, 0, 0, &mut self.scr.at(1), &mut nb, &mut nw);
        let mut lbl = 110;
        loop {
            match lbl {
                110 => {
                    c.contio(ng, 0, 0, &mut self.scr.at(1), &mut nb, &mut nw);
                    if c.math == -1 {
                        let strng = format!("mat={} not on tape {}", ffmt::i(self.matd, 4), ffmt::i(ng, 3));
                        c.error("gamll", &strng, " ");
                    }
                    if c.math == 0 || c.mfh == 0 || c.mth == 0 {
                        lbl = 110;
                        continue;
                    }
                    if c.math == self.matd {
                        lbl = 120;
                        continue;
                    }
                    c.tomend(ng, 0, 0, &mut self.scr.at(1));
                    lbl = 110;
                }
                120 => {
                    self.iwa = 1;
                    self.iwf = 0;
                    self.nfs = 0;
                    self.iwr = 0;
                    c.nsp = 0;
                    self.ntp = -1;
                    self.jwf = 0;
                    ldin = 0;
                    ldp = 0;
                    ld2n = 0;
                    let iza = nint(c.c1h);
                    self.nid = iza * 10;
                    lbl = 140;
                }
                130 => {
                    c.contio(ng, 0, 0, &mut self.scr.at(1), &mut nb, &mut nw);
                    lbl = if c.math != self.matd { 400 } else { 140 };
                }
                140 => {
                    self.nsigz = c.l2h;
                    if self.nsigz == 1 {
                        *iff = 0;
                    }
                    if self.nsigz > 1 && *iff > 0 {
                        self.iwr = 1;
                    }
                    if self.nsigz > nsgz && nsgz > 0 {
                        self.nsigz = nsgz;
                    }
                    if self.nsigz > nszmax {
                        let strng = format!("{} allowed", ffmt::i(nszmax, 2));
                        c.error("gamll", "reqested too many sigma zeroes", &strng);
                    }
                    ntw = c.n2h;
                    c.listio(ng, 0, 0, &mut self.scr.at(1), &mut nb, &mut nw);
                    let temp = c.c1h;
                    if (temp - self.rtemp).abs() <= eps {
                        lbl = 160;
                        continue;
                    }
                    if temp > self.rtemp && self.ntp < 0 {
                        c.error("gamll", "reference temp not on gendf tape.", " ");
                    }
                    if self.ntp < 0 {
                        lbl = 150;
                        continue;
                    }
                    self.ntp += 1;
                    if self.ntp > ntpmax {
                        let strng = format!("{} allowed", ffmt::i(ntpmax, 2));
                        c.error("gamll", "requested too many temperatures", &strng);
                    }
                    self.tmpr[self.ntp] = temp;
                    lbl = 220;
                }
                150 => {
                    c.tomend(ng, 0, 0, &mut self.scr.at(1));
                    lbl = 130;
                }
                160 => {
                    let temp = c.c1h;
                    self.ngn = c.l1h;
                    if self.ngn != ngnd && self.ngn != self.ngnf {
                        c.error("gamll", "incorrect group structure or group type.", " ");
                    }
                    if *igprnt <= 0 {
                        let mut i1 = 7 + ntw + self.nsigz + self.ngn;
                        if self.ngn == ngnd {
                            let mut v = vec![self.scr[i1]];
                            for i in 1..=self.ngn {
                                v.push(self.scr[i1 - i]);
                            }
                            // '(/'' neutron group structure''/'' ---...'',(1p,6e12.4))'
                            let mut s = String::from("\n neutron group structure\n -----------------------");
                            let first = v.len().min(6);
                            for &x in &v[..first] {
                                s.push_str(&pe(x, 12, 4));
                            }
                            for ch in v[first..].chunks(6) {
                                s.push('\n');
                                for &x in ch {
                                    s.push_str(&pe(x, 12, 4));
                                }
                            }
                            c.io.so(&s);
                        }
                        if self.ngn != ngnd {
                            let ngnd1 = ngnd + 1;
                            let ngn1 = self.ngn + 1;
                            let mut eg = RV::new(ngnd1 as usize + SLACK);
                            let mut lastg = 0;
                            i1 -= ngn1;
                            for i in 1..=ngn1 {
                                let ig = self.icgrp(i);
                                if ig != lastg {
                                    if ig != 0 {
                                        let jg = ngnd1 - ig;
                                        eg[1 + jg] = self.scr[i1 + i];
                                    }
                                    lastg = ig;
                                }
                            }
                            c.io.so(&format!(
                                "\n neutron group structure\n -----------------------\n{}",
                                rowsj(&eg.0[..ngnd1 as usize], 6, |x| pe(x, 12, 4))
                            ));
                        }
                        *igprnt = 1;
                    }
                    self.ntp = 0;
                    if self.rtemp > 0.0 {
                        self.ntp = 1;
                    }
                    if self.rtemp > 0.0 {
                        self.tmpr[1] = temp;
                    }
                    for i in 1..=self.nsigz {
                        self.sigz[i] = self.scr[i + ntw + lz];
                    }
                    if self.izref == 1 {
                        c.io.so("\n ref. sigma zero for elastic matrix = infinity");
                    }
                    if self.izref != 1 {
                        c.io.so(&format!(
                            "\n ref. sigma zero for elastic matrix ={}",
                            pe(self.sigz[self.izref], 10, 2)
                        ));
                    }
                    lain = ngnd;
                    la2n = ngnd;
                    lap = ngnd;
                    lbl = 220;
                }
                220 => {
                    c.tosend(ng, 0, 0, &mut self.scr.at(1));
                    lbl = 250;
                }
                //--process relevant sections for this material.
                250 => {
                    c.contio(ng, 0, 0, &mut self.scr.at(1), &mut nb, &mut nw);
                    if c.math == 0 && self.nsigz > 1 {
                        lbl = 130;
                        continue;
                    }
                    if c.math != self.matd {
                        lbl = 400;
                        continue;
                    }
                    if c.mfh == 0 || c.mth == 0 {
                        lbl = 250;
                        continue;
                    }
                    if c.mfh == 3 {
                        //--check for presence of fission.
                        if c.mth == 18 || c.mth == 19 {
                            self.iwf = 1;
                        }
                        c.tosend(ng, 0, 0, &mut self.scr.at(1));
                        lbl = 250;
                        continue;
                    }
                    if c.mfh == 5 {
                        //--check for presence of delayed neutron spectra (mf=5)
                        if c.mth == 455 {
                            let nl = c.l1h;
                            self.nfs += nl;
                        }
                        c.tosend(ng, 0, 0, &mut self.scr.at(1));
                        lbl = 250;
                        continue;
                    }
                    if c.mfh == 6 {
                        lbl = 300;
                        continue;
                    }
                    c.tosend(ng, 0, 0, &mut self.scr.at(1));
                    lbl = 250;
                }
                //--process transfer matrix data (mf=6).
                300 => {
                    c.listio(ng, 0, 0, &mut self.scr.at(1), &mut nb, &mut nw);
                    let mut l = 1;
                    while nb != 0 {
                        l += nw;
                        c.moreio(ng, 0, 0, &mut self.scr.at(l), &mut nb, &mut nw);
                    }
                    ig = c.n2h;
                    igc = self.icgrp(ig);
                    if igc == 0 || igc > ngnd {
                        lbl = 350;
                        continue;
                    }
                    let mut ig2lo = c.l2h;
                    let mut idone = 0;
                    while idone == 0 {
                        ig2loc = self.icgrp(ig2lo);
                        if ig2loc != 0 && ig2loc <= ngnd {
                            idone = 1;
                        } else {
                            if ig2lo == ig {
                                idone = 1;
                            } else {
                                ig2lo += 1;
                            }
                        }
                    }
                    if c.mth < 51 || c.mth > 91 {
                        lbl = 360;
                        continue;
                    }
                    // inelastic
                    if lain > igc {
                        lain = igc;
                    }
                    if ldin < igc - ig2loc {
                        ldin = igc - ig2loc;
                    }
                    lbl = 350;
                }
                350 => {
                    lbl = if ig < self.ngn { 300 } else { 250 };
                }
                360 => {
                    if c.mth != 2 {
                        lbl = 370;
                        continue;
                    }
                    // elastic
                    if lap > igc {
                        lap = igc;
                    }
                    if ldp < igc - ig2loc {
                        ldp = igc - ig2loc;
                    }
                    lbl = if ig < self.ngn { 300 } else { 250 };
                }
                370 => {
                    if c.mth == 16 || (c.mth > 5 && c.mth < 10) || (c.mth > 45 && c.mth < 50) {
                        // n2n
                        if la2n > igc {
                            la2n = igc;
                        }
                        if ld2n < igc - ig2loc {
                            ld2n = igc - ig2loc;
                        }
                        lbl = 250;
                        continue;
                    }
                    c.tosend(ng, 0, 0, &mut self.scr.at(1));
                    lbl = 250;
                }
                //--set up storage pointers.
                _ => {
                    self.xla[1] = 0.0;
                    self.xla[2] = 0.0;
                    self.xla[3] = 0.0;
                    self.xlol[1] = 0.0;
                    self.xlol[2] = 0.0;
                    self.xlol[3] = 0.0;
                    self.xld[1] = ldin as f64;
                    self.xld[2] = ld2n as f64;
                    self.xld[3] = ldp as f64;
                    let ngf = ngnd as f64;
                    if ldin != 0 {
                        self.xla[1] = (ngnd - lain + 1) as f64;
                        let mut lx = nint(self.xla[1]) + nint(self.xld[1]) - ngnd - 1;
                        if lx < 0 {
                            lx = 0;
                        }
                        self.xlol[1] = self.xla[1] * (self.xld[1] + 1.0) - ((lx * (lx - 1)) / 2) as f64;
                    }
                    if ld2n != 0 {
                        self.xla[2] = (ngnd - la2n + 1) as f64;
                        let mut lx = nint(self.xla[2]) + nint(self.xld[2]) - ngnd - 1;
                        if lx < 0 {
                            lx = 0;
                        }
                        self.xlol[2] = self.xla[2] * (self.xld[2] + 1.0) - ((lx * (lx - 1)) / 2) as f64;
                    }
                    if ldp != 0 {
                        self.xla[3] = (ngnd - lap + 1) as f64;
                        let mut lx = nint(self.xla[3]) + nint(self.xld[3]) - ngnd - 1;
                        if lx < 0 {
                            lx = 0;
                        }
                        self.xlol[3] = self.xla[3] * (self.xld[3] + 1.0) - ((lx * (lx - 1)) / 2) as f64;
                    }
                    let _ = ngf;
                    return;
                }
            }
        }
    }

    /// Read through ngendf tape again and accumulate desired cross
    /// sections and matrices in the array a.
    fn gamxs(&mut self, c: &mut Ctx, loc_: &mut [i32; 11], a: &mut RV, nloc: i32) {
        let eps = 1.0e-4f64;
        let ng = self.ngendf;
        let ngnd = self.ngnd;
        let (mut nb, mut nw) = (0i32, 0i32);

        //--initialize data pointers.
        let ngnd1 = ngnd + 1;
        let ldin = nint(self.xld[1]);
        let kdin = ngnd1 - ldin;
        let ld2n = nint(self.xld[2]);
        let kd2n = ngnd1 - ld2n;
        let ldp = nint(self.xld[3]);
        let kdp = ngnd1 - ldp;
        let lz = 6;
        let mut mf3mt2 = 0;
        let mut dnu = 0.0f64;
        let mut jgdnu = 0;
        let mut dnorm = 0.0f64;
        let mut cnorm = 0.0f64;
        let iloc = 1;
        let locab0 = iloc;
        let locsf0 = locab0 + ngnd;
        let locchi = locsf0 + ngnd;
        let locdla = locchi + ngnd;
        let locnus = locdla + ngnd * self.nfs;
        let loce0 = locnus + ngnd;
        let loce1 = loce0 + nint(self.xlol[3]);
        let locin = loce1 + nint(self.xlol[3]);
        let loc2n = locin + nint(self.xlol[1]);
        let locabs = loc2n + nint(self.xlol[2]);
        loc_[1] = locab0;
        loc_[2] = locsf0;
        loc_[3] = locchi;
        loc_[4] = locnus;
        loc_[5] = loce0;
        loc_[6] = loce1;
        loc_[7] = locin;
        loc_[8] = loc2n;
        loc_[9] = locabs;
        let na = locabs - 1;
        if na - iloc > nloc {
            c.error("gamxs", "storage exceeded.", " ");
        }
        for i in iloc..=na {
            a[i] = 0.0;
        }
        for i in 1..=ngnd {
            self.cflux[i] = 0.0;
            self.cflux[i + ngnd] = 0.0;
        }
        let mut cspc = RV::new(200 + SLACK);
        for i in 1..=self.ngn {
            cspc[i] = 0.0;
        }
        let mut sumd = [0.0f64; 7];
        let mut nldla = 0;
        let mut i318 = 0;
        let mut i618 = 0;
        let mut tfound = 0.0f64;
        let (mut nl, mut nz, mut ng2, mut ig2lo, mut ig, mut igc, mut jg) = (0, 0, 0, 0, 0, 0, 0);
        let mut l = 1;
        c.repoz(ng);
        c.tpidio(ng, 0, 0, &mut self.scr.at(1), &mut nb, &mut nw);
        let mut lbl = 120;
        loop {
            match lbl {
                120 => {
                    c.contio(ng, 0, 0, &mut self.scr.at(1), &mut nb, &mut nw);
                    if c.math == self.matd {
                        lbl = 130;
                        continue;
                    }
                    c.tomend(ng, 0, 0, &mut self.scr.at(1));
                    lbl = 120;
                }
                //--process relevant sections for reference temperature.
                130 => {
                    tfound = 0.0;
                    if c.mfh != 1 {
                        lbl = 150;
                        continue;
                    }
                    c.tosend(ng, 0, 0, &mut self.scr.at(1));
                    lbl = 140;
                }
                140 => {
                    l = 1;
                    c.contio(ng, 0, 0, &mut self.scr.at(l), &mut nb, &mut nw);
                    if c.math == 0 {
                        lbl = 500;
                        continue;
                    }
                    if c.mfh == 0 || c.mth == 0 {
                        lbl = 140;
                        continue;
                    }
                    if c.mfh != 1 {
                        lbl = 150;
                        continue;
                    }
                    c.tosend(ng, 0, 0, &mut self.scr.at(l));
                    lbl = 140;
                }
                150 => {
                    nl = c.l1h;
                    nz = c.l2h;
                    if c.mfh != 3 || c.mth < 2 || mf3mt2 > 0 {
                        lbl = 170;
                        continue;
                    }
                    mf3mt2 = 1;
                    for i in 1..=ngnd {
                        if self.cflux[i] != ZERO {
                            self.cflux[i] = 1.0 / self.cflux[i];
                        }
                        if self.cflux[ngnd + i] != ZERO {
                            self.cflux[ngnd + i] = 1.0 / self.cflux[ngnd + i];
                        }
                    }
                    lbl = 170;
                }
                170 => {
                    c.listio(ng, 0, 0, &mut self.scr.at(l), &mut nb, &mut nw);
                    let temp = c.c1h;
                    if (temp - self.rtemp).abs() > eps && tfound == 0.0 {
                        lbl = 210;
                        continue;
                    }
                    if (temp - self.rtemp).abs() <= eps {
                        tfound = 1.0;
                    }
                    ng2 = c.l1h;
                    ig2lo = c.l2h;
                    ig = c.n2h;
                    igc = 1;
                    if ig > 0 {
                        igc = self.icgrp(ig);
                    }
                    while nb != 0 {
                        l += nw;
                        c.moreio(ng, 0, 0, &mut self.scr.at(l), &mut nb, &mut nw);
                        if l + nw - 1 > self.nwscr {
                            c.error("gamxs", "storage exceeded.", " ");
                        }
                    }
                    if igc == 0 || igc > ngnd {
                        lbl = 260;
                        continue;
                    }
                    jg = ngnd - igc + 1;
                    l = 1;
                    if c.mth == 455 {
                        lbl = 330;
                        continue;
                    }
                    if c.mfh == 3 {
                        lbl = 250;
                        continue;
                    }
                    if c.mfh == 6 {
                        lbl = 300;
                        continue;
                    }
                    lbl = 200;
                }
                200 => {
                    c.tosend(ng, 0, 0, &mut self.scr.at(l));
                    lbl = 140;
                }
                210 => {
                    c.tomend(ng, 0, 0, &mut self.scr.at(l));
                    lbl = 140;
                }
                //--file 3 cross sections
                250 => {
                    if c.mth != 1 {
                        lbl = 270;
                        continue;
                    }
                    // coarse group p0 flux
                    self.cflux[jg] = self.cflux[jg] + self.scr[l + lz];
                    if nl == 1 {
                        lbl = 260;
                        continue;
                    }
                    // coarse group p1 flux
                    self.cflux[ngnd + jg] = self.cflux[ngnd + jg] + self.scr[l + lz + 1];
                    lbl = 260;
                }
                260 => {
                    lbl = if ig < self.ngn { 170 } else { 140 };
                }
                270 => {
                    if c.mth < 102 || c.mth > 150 {
                        lbl = 280;
                        continue;
                    }
                    // capture
                    let loc = locab0 + jg - 1;
                    let loca = l + lz + nl * nz;
                    a[loc] = a[loc] + self.scr[loca] * self.scr[l + lz] * self.cflux[jg];
                    lbl = if ig < self.ngn { 170 } else { 140 };
                }
                280 => {
                    if (c.mth < 18 || c.mth > 21) && c.mth != 38 {
                        lbl = 200;
                        continue;
                    }
                    // fission
                    if c.mth == 18 {
                        i318 = 1;
                    }
                    if c.mth > 18 && i318 > 0 {
                        lbl = 200;
                        continue;
                    }
                    let loc = locsf0 + jg - 1;
                    let loca = l + lz + nl * nz;
                    a[loc] = a[loc] + self.scr[loca] * self.scr[l + lz] * self.cflux[jg];
                    let loc = locab0 + jg - 1;
                    a[loc] = a[loc] + self.scr[loca] * self.scr[l + lz] * self.cflux[jg];
                    lbl = if ig < self.ngn { 170 } else { 140 };
                }
                //--file 6 matrices
                300 => {
                    if igc == 0 || igc > ngnd {
                        lbl = 320;
                        continue;
                    }
                    if (c.mth < 18 || c.mth > 21) && c.mth != 38 {
                        lbl = 370;
                        continue;
                    }
                    //--fission
                    if c.mth == 18 {
                        i618 = 1;
                    }
                    if c.mth > 18 && i618 > 0 {
                        lbl = 200;
                        continue;
                    }
                    if ig != 0 {
                        for k in 2..=ng2 {
                            let loca = l + lz + nl * nz * (k - 1);
                            let loc = locnus + jg - 1;
                            if ig2lo != 0 {
                                // matrix part
                                a[loc] = a[loc] + self.scr[loca] * self.scr[l + lz] * self.cflux[jg];
                                let locf = l + lz;
                                let ig2 = ig2lo + k - 2;
                                let ig2c = self.icgrp(ig2);
                                if ig2c != 0 && ig2c <= ngnd {
                                    let jg2c = ngnd - ig2c + 1;
                                    let loc = locchi + jg2c - 1;
                                    a[loc] = a[loc] + self.scr[locf] * self.scr[loca];
                                    cnorm = cnorm + self.scr[locf] * self.scr[loca];
                                }
                            } else {
                                // spectrum part
                                for i in 1..=self.ngn {
                                    a[loc] = a[loc]
                                        + cspc[i] * self.scr[loca] * self.scr[loca - 1] * self.cflux[jg];
                                    let ig2c = self.icgrp(i);
                                    if ig2c != 0 && ig2c <= ngnd {
                                        let jg2c = ngnd - ig2c - 1;
                                        let locc = locchi + jg2c - 1;
                                        a[locc] = a[locc] + cspc[i] * self.scr[loca] * self.scr[loca - 1];
                                        cnorm = cnorm + self.scr[loca] * self.scr[loca - 1];
                                    }
                                }
                            }
                        }
                    } else {
                        // save constant spectrum
                        for k in 1..=ng2 {
                            cspc[ig2lo + k - 1] = self.scr[l + lz + nl * nz * (k - 1)];
                        }
                    }
                    lbl = 320;
                }
                320 => {
                    lbl = if ig < self.ngn { 170 } else { 140 };
                }
                // delayed neutrons
                330 => {
                    if c.mfh != 5 {
                        let loc = locnus + jg - 1;
                        let loca = l + lz;
                        a[loc] = a[loc]
                            + self.scr[loca + 1] * self.scr[loca + 2] * self.scr[l + lz] * self.cflux[jg];
                        dnorm = dnorm + self.scr[loca] * self.scr[loca + 1] * self.scr[loca + 2];
                        let mut test = 1.0f64;
                        test = test / 10.0;
                        if self.scr[loca + 2] >= test && dnu == ZERO {
                            jgdnu = jg;
                        }
                        if self.scr[loca + 2] >= test && dnu == ZERO {
                            dnu = self.scr[loca + 1];
                        }
                    } else {
                        nldla = nl;
                        if nldla != 6 {
                            c.error("gamxs", "no. deleay groups not 6.", " ");
                        }
                        for il in 1..=nl {
                            sumd[il as usize] = 0.0;
                        }
                        let locb = l + lz - 1;
                        for k in 2..=ng2 {
                            let ig2 = ig2lo + k - 2;
                            let ig2c = self.icgrp(ig2);
                            if ig2c != 0 && ig2c <= ngnd {
                                let jgc = ngnd - ig2c + 1;
                                let loc = locchi + jgc - 1;
                                let loca = l + lz - 1 + nl * (k - 1);
                                let locd = locdla + jgc - 1;
                                for jgt in 1..=nl {
                                    a[loc] = self.scr[jgt + loca] * dnorm + a[loc];
                                    a[ngnd * (jgt - 1) + locd] =
                                        a[ngnd * (jgt - 1) + locd] + self.scr[jgt + loca] * dnu;
                                    if jgc == ngnd {
                                        a[ngnd * (jgt - 1) + locd] = a[jgt + locb];
                                    }
                                    if jgc < ngnd - 1 {
                                        sumd[jgt as usize] = sumd[jgt as usize] + self.scr[jgt + loca] * dnu;
                                    }
                                    cnorm = self.scr[jgt + loca] * dnorm + cnorm;
                                }
                            }
                        }
                    }
                    lbl = if ig < self.ngn { 170 } else { 140 };
                }
                //--inelastic
                370 => {
                    if c.mth < 51 || c.mth > 91 {
                        lbl = 410;
                        continue;
                    }
                    let mut jgp = jg;
                    if jg > kdin {
                        jgp = kdin;
                    }
                    let mut lc = (ldin + 1) * (jgp - 1);
                    if jg > kdin {
                        let lim = jg - kdin;
                        for jgt in 1..=lim {
                            lc = lc + ldin + 2 - jgt;
                        }
                    }
                    for k in 2..=ng2 {
                        let ig2 = ig2lo + k - 2;
                        let ig2c = self.icgrp(ig2);
                        if ig2c != 0 && ig2c <= ngnd {
                            let mut jg2c = ngnd - ig2c + 1;
                            if jg2c < jg {
                                jg2c = jg;
                            }
                            let loca = l + lz + nl * nz * (k - 1);
                            let loc = locin + lc + jg2c - jg;
                            a[loc] = a[loc] + self.scr[loca] * self.scr[l + lz] * self.cflux[jg];
                        }
                    }
                    lbl = if ig < self.ngn { 170 } else { 140 };
                }
                //--elastic
                410 => {
                    if c.mth != 2 {
                        lbl = 450;
                        continue;
                    }
                    let mut jgp = jg;
                    if jg > kdp {
                        jgp = kdp;
                    }
                    let mut lc = (ldp + 1) * (jgp - 1);
                    if jg > kdp {
                        let lim = jg - kdp;
                        for jgt in 1..=lim {
                            lc = lc + ldp + 2 - jgt;
                        }
                    }
                    for k in 2..=ng2 {
                        let ig2 = ig2lo + k - 2;
                        let ig2c = self.icgrp(ig2);
                        if ig2c != 0 && ig2c <= ngnd {
                            let mut jg2c = ngnd - ig2c + 1;
                            if jg2c < jg {
                                jg2c = jg;
                            }
                            let loc = loce0 + lc + jg2c - jg;
                            let mut loca = l + lz + nl * (self.izref - 1 + nz * (k - 1));
                            a[loc] = a[loc] + self.scr[loca] * self.scr[l + lz] * self.cflux[jg];
                            if nl != 1 {
                                loca += 1;
                                let loc = loce1 + lc + jg2c - jg;
                                let mut flx = self.cflux[ngnd + jg];
                                if self.cflux[ngnd + jg] == 0.0 {
                                    flx = self.cflux[jg];
                                }
                                a[loc] = a[loc] + 3.0 * self.scr[loca] * self.scr[l + lz + 1] * flx;
                            }
                        }
                    }
                    lbl = if ig < self.ngn { 170 } else { 140 };
                }
                450 => {
                    if c.mth == 16 {
                        lbl = 460;
                        continue;
                    }
                    if c.mth < 6 || c.mth > 49 {
                        lbl = 200;
                        continue;
                    }
                    if c.mth > 9 && c.mth < 46 {
                        lbl = 200;
                        continue;
                    }
                    lbl = 460;
                }
                //--n2n
                460 => {
                    let mut jgp = jg;
                    if jg > kd2n {
                        jgp = kd2n;
                    }
                    let mut lc = (ld2n + 1) * (jgp - 1);
                    if jg > kd2n {
                        let lim = jg - kd2n;
                        for jgt in 1..=lim {
                            lc = lc + ld2n + 2 - jgt;
                        }
                    }
                    for k in 2..=ng2 {
                        let ig2 = ig2lo + k - 2;
                        let ig2c = self.icgrp(ig2);
                        if ig2c != 0 && ig2c <= ngnd {
                            let mut jg2c = ngnd - ig2c + 1;
                            if jg2c < jg {
                                jg2c = jg;
                            }
                            let loca = l + lz + nl * nz * (k - 1);
                            let loc = loc2n + lc + jg2c - jg;
                            a[loc] = a[loc] + self.scr[loca] * self.scr[l + lz] * self.cflux[jg] / 2.0;
                        }
                    }
                    lbl = if ig < self.ngn { 170 } else { 140 };
                }
                //--finished loading a.
                //--find nu from nu times sigma f.
                _ => {
                    for ig in 1..=ngnd {
                        let loc1 = locsf0 + ig - 1;
                        if a[loc1] != ZERO {
                            let loc2 = locnus + ig - 1;
                            a[loc2] = a[loc2] / a[loc1];
                        }
                    }

                    //--normalize chi and delayed chi
                    if cnorm != ZERO {
                        cnorm = 1.0 / cnorm;
                        for ig in 1..=ngnd {
                            a[ig - 1 + locchi] = a[ig - 1 + locchi] * cnorm;
                        }
                        if dnu != ZERO {
                            let rnorm = 1.0 / a[jgdnu - 1 + locnus];
                            let lim = ngnd - 1;
                            for ig in 1..=lim {
                                let locd = locdla - 1 + ig;
                                for il in 1..=nldla {
                                    if ig < lim {
                                        a[ngnd * (il - 1) + locd] = a[ngnd * (il - 1) + locd] * rnorm;
                                    } else {
                                        a[ngnd * (il - 1) + locd] = sumd[il as usize] * rnorm;
                                    }
                                }
                            }
                            let fracd = dnu / a[jgdnu - 1 + locnus];
                            let mut v1 = String::new();
                            for il in 1..=nldla {
                                v1.push_str(&pe(a[ngnd * (il - 1) + locdla + ngnd - 1], 13, 5));
                            }
                            let mut v2 = String::new();
                            for il in 1..=nldla {
                                v2.push_str(&pe(a[ngnd * (il - 1) + locdla + ngnd - 2], 13, 5));
                            }
                            // nldla is 6 or (no delayed spectra read) 0; with
                            // no items left the format stops at the first
                            // e13.5, after '1x' (which writes nothing)
                            if nldla == 0 {
                                c.io.so(&format!(
                                    "\n nudlay ={}\n delayed fraction ={}\n time constants\n",
                                    pe(dnu, 13, 5),
                                    pe(fracd, 13, 5)
                                ));
                            } else {
                                c.io.so(&format!(
                                    "\n nudlay ={}\n delayed fraction ={}\n time constants\n {}\n time group fractions\n {}",
                                    pe(dnu, 13, 5),
                                    pe(fracd, 13, 5),
                                    v1,
                                    v2
                                ));
                            }
                        }
                    }

                    //--for reference temperatures greater than zero...
                    if self.rtemp > 0.0 {
                        c.skiprz(ng, -2);
                    }
                    if self.rtemp > 0.0 {
                        c.findf(self.matd, 1, 451, ng);
                    }
                    let _ = (nb, nw);
                    return;
                }
            }
        }
    }

    /// Calculate f-factors for absorption and fission.
    fn gamff(&mut self, c: &mut Ctx, ngmin: i32, ngmax: i32, loc_: &mut [i32; 11], a: &mut RV, nloc: i32) {
        let toler = 0.01f64;
        let ng = self.ngendf;
        let ngnd = self.ngnd;
        let (mut nb, mut nw) = (0i32, 0i32);
        let locab0 = loc_[1];
        let locsf0 = loc_[2];
        let locabs = loc_[9];

        //--initialize.
        let iloc = 1;
        let lz = 6;
        let ntp = self.ntp;
        let nsigz = self.nsigz;
        let locsf = locabs + ngnd * ntp * nsigz;
        loc_[10] = locsf;
        let na = locsf + ngnd * ntp * nsigz - 1;
        if na - iloc > nloc {
            c.error("gamff", "storage exceeded.", " ");
        }
        for i in locabs..=na {
            a[i] = 0.0;
        }
        let mut itp = 1;
        let mut mf3mt1 = 1;
        let mut i318 = 0;
        let (mut nl, mut nz, mut lim, mut ig, mut igc, mut jg) = (0, 0, 0, 0, 0, 0);
        let mut l = 1;

        //--read through rest of temperatures for this mat,
        //--and store absorption and fission cross-sections.
        let mut lbl = 110;
        loop {
            match lbl {
                110 => {
                    l = 1;
                    c.contio(ng, 0, 0, &mut self.scr.at(l), &mut nb, &mut nw);
                    if c.math > 0 {
                        lbl = 120;
                        continue;
                    }
                    if itp == ntp {
                        lbl = 300;
                        continue;
                    }
                    itp += 1;
                    i318 = 0;
                    lbl = 120;
                }
                120 => {
                    if c.mfh == 0 || c.mth == 0 {
                        lbl = 110;
                        continue;
                    }
                    if c.mfh == 3 {
                        lbl = 130;
                        continue;
                    }
                    c.tosend(ng, 0, 0, &mut self.scr.at(l));
                    lbl = 110;
                }
                130 => {
                    if mf3mt1 <= 0 && c.mfh == 3 && c.mth >= 2 {
                        mf3mt1 = 1;
                        for i in 1..=ngnd {
                            if self.cflux[i] != ZERO {
                                self.cflux[i] = 1.0 / self.cflux[i];
                            }
                        }
                    }
                    if c.mth > 101 && c.mth < 108 {
                        lbl = 160;
                        continue;
                    }
                    if (c.mth > 17 && c.mth < 22) || c.mth == 38 {
                        lbl = 160;
                        continue;
                    }
                    if c.mth == 1 {
                        lbl = 160;
                        continue;
                    }
                    lbl = 150;
                }
                150 => {
                    c.tosend(ng, 0, 0, &mut self.scr.at(l));
                    lbl = 110;
                }
                160 => {
                    nl = c.l1h;
                    nz = c.l2h;
                    lim = nz;
                    if nsigz < nz {
                        lim = nsigz;
                    }
                    lbl = 170;
                }
                170 => {
                    l = 1;
                    c.listio(ng, 0, 0, &mut self.scr.at(l), &mut nb, &mut nw);
                    ig = c.n2h;
                    igc = self.icgrp(ig);
                    while nb != 0 {
                        l += nw;
                        if l + nw - 1 > self.nwscr {
                            c.error("gamff", "storage exceeded.", " ");
                        }
                        c.moreio(ng, 0, 0, &mut self.scr.at(l), &mut nb, &mut nw);
                    }
                    if igc == 0 || igc > ngnd {
                        lbl = 210;
                        continue;
                    }
                    jg = ngnd - igc + 1;
                    if c.mth < 102 {
                        lbl = 220;
                        continue;
                    }
                    // capture
                    let loc = locabs - 1 + nsigz * ((itp - 1) + ntp * (jg - 1));
                    let loca = l + lz + nl * (nz - 1);
                    for jz in 1..=lim {
                        a[jz + loc] = self.scr[nl * jz + loca] * self.scr[l + lz] * self.cflux[jg] + a[jz + loc];
                    }
                    lbl = 210;
                }
                210 => {
                    lbl = if ig < self.ngn { 170 } else { 110 };
                }
                220 => {
                    if c.mth == 1 {
                        // flux
                        if igc <= 1 {
                            for i in 1..=ngnd {
                                self.cflux[i] = 0.0;
                            }
                            mf3mt1 = 0;
                        }
                        self.cflux[jg] = self.cflux[jg] + self.scr[l + lz];
                        lbl = 210;
                        continue;
                    }
                    // fission
                    if c.mth == 18 {
                        i318 = 1;
                    }
                    if c.mth > 18 && i318 > 0 {
                        lbl = 150;
                        continue;
                    }
                    self.jwf = 1;
                    let iadd = nsigz * ((itp - 1) + ntp * (jg - 1)) - 1;
                    let loca = l + lz + nl * (nz - 1);
                    for jz in 1..=lim {
                        a[jz + iadd + locsf] =
                            a[jz + iadd + locsf] + self.scr[nl * jz + loca] * self.scr[l + lz] * self.cflux[jg];
                        a[jz + iadd + locabs] =
                            a[jz + iadd + locabs] + self.scr[nl * jz + loca] * self.scr[l + lz] * self.cflux[jg];
                    }
                    lbl = 210;
                }
                //--finished loading a.
                //--calculate f-factors.
                _ => {
                    let mut iglo = ngmin;
                    for ig in ngmin..=ngmax {
                        let mut ab0 = a[locab0 + ig - 1];
                        if ab0 != ZERO {
                            ab0 = 1.0 / ab0;
                        }
                        let mut sf0 = a[locsf0 + ig - 1];
                        if sf0 != ZERO {
                            sf0 = 1.0 / sf0;
                        }
                        for jz in 1..=nsigz {
                            for itp in 1..=ntp {
                                let iadd = jz - 1 + nsigz * ((itp - 1) + ntp * (ig - 1));
                                let mut loc = locabs + iadd;
                                if ab0 == ZERO {
                                    a[loc] = 1.0;
                                }
                                if ab0 != ZERO {
                                    a[loc] = a[loc] * ab0;
                                }
                                if a[loc] <= ZERO {
                                    a[loc] = 1.0;
                                }
                                if (1.0 - a[loc]).abs() > toler && iglo == ngmin {
                                    iglo = ig;
                                }
                                a[loc] = log(self.sigz[jz]) - 2.0 * log(a[loc]);
                                if self.jwf != 0 {
                                    loc = locsf + iadd;
                                    if sf0 == ZERO {
                                        a[loc] = 1.0;
                                    }
                                    if sf0 != ZERO {
                                        a[loc] = a[loc] * sf0;
                                    }
                                    if a[loc] <= ZERO {
                                        a[loc] = 1.0;
                                    }
                                    a[loc] = log(self.sigz[jz]) - 2.0 * log(a[loc]);
                                }
                            }
                        }
                    }
                    self.iwr = ngmax - iglo + 1;
                    return;
                }
            }
        }
    }
}

/// Print out scattering matrices in packed GAM format.
fn pgam(c: &mut Ctx, a: &[f64], ng: i32, la: i32, ld: i32) {
    let ncol = 5;
    let ncol1 = ncol - 1;
    let a = |k: i32| a[(k - 1) as usize];
    let mut h = String::from("\n initl final    cross section vs final group\n group group    +0");
    for i in 1..=ncol1 {
        h.push_str(&format!("{}+{}", " ".repeat(10), ffmt::i(i, 1)));
    }
    h.push_str(&format!("\n {}--------", "----------".repeat(7)));
    c.io.so(&h);
    let mut k = 0;
    for ig in 1..=la {
        let mut n = ng - ig + 1;
        if n > ld + 1 {
            n = ld + 1;
        }
        if n <= ld {
            n += 1;
        }
        // do not print leading and trailing zeroes
        let mut i1 = 0;
        let mut i2 = 0;
        for j in 1..=n {
            if i1 == 0 && a(k + j) != 0.0 {
                i1 = j;
            }
            if a(k + j) != 0.0 {
                i2 = j;
            }
        }
        if i1 != 0 {
            let mut nw = i2 - i1 + 1;
            let mut m = nw;
            k = k + i1 - 1;
            if m > ncol {
                m = ncol;
            }
            let mut ig2 = ig + i1 - 1;
            let mut s = format!("{}{}  ", ffmt::i(ig, 5), ffmt::i(ig2, 6));
            for i in 1..=m {
                s.push_str(&pe(a(k + i), 12, 4));
            }
            c.io.so(&s);
            while nw > 0 {
                k += m;
                nw -= m;
                if nw > 0 {
                    ig2 += ncol;
                    m = nw;
                    if m > ncol {
                        m = ncol;
                    }
                    let mut s = format!("     {}  ", ffmt::i(ig2, 6));
                    for i in 1..=m {
                        s.push_str(&pe(a(k + i), 12, 4));
                    }
                    c.io.so(&s);
                }
            }
        }
        k = k + n - i2;
    }
}

/// Print out form factors in packed GAM format.
fn pgamff(c: &mut Ctx, a: &[f64], ntp: i32, nsigz: i32, iwr: i32, sigz: &RV, tmpr: &RV, ngmax: i32) {
    let jopt = 1;
    let a3 = |is: i32, it: i32, ig: i32| a[((is - 1) + nsigz * ((it - 1) + ntp * (ig - 1))) as usize];
    let mut s = String::from("     temperatures (down)");
    for (k, i) in (1..=ntp).enumerate() {
        if k > 0 && k % 10 == 0 {
            s.push_str("\n");
        }
        s.push_str(&ffmt::f(tmpr[i], 7, 0));
    }
    c.io.so(&s);
    let mut s = String::from("     sigzeroes (across)\n         ");
    for (k, i) in (1..=nsigz).enumerate() {
        if k > 0 && k % 10 == 0 {
            s.push_str("\n         ");
        }
        s.push_str(&pe(sigz[i], 12, 0));
    }
    c.io.so(&s);
    let n = ngmax + 1 - iwr;
    let mut f = [0.0f64; 11];
    for ig in n..=ngmax {
        let mut s = format!("\n group{}   ", ffmt::i(ig, 2));
        if jopt <= 0 {
            for is in 1..=nsigz {
                s.push_str(&pe(a3(is, 1, ig), 12, 4));
            }
        } else {
            for is in 1..=nsigz {
                f[is as usize] = exp(-a3(is, 1, ig) / 2.0) * sqrt(sigz[is]);
            }
            for is in 1..=nsigz {
                s.push_str(&pe(f[is as usize], 12, 4));
            }
        }
        c.io.so(&s);
        if ntp != 1 {
            for it in 2..=ntp {
                let mut s = " ".repeat(11);
                if jopt <= 0 {
                    for is in 1..=nsigz {
                        s.push_str(&pe(a3(is, it, ig), 12, 4));
                    }
                } else {
                    for is in 1..=nsigz {
                        f[is as usize] = exp(-a3(is, it, ig) / 2.0) * sqrt(sigz[is]);
                    }
                    for is in 1..=nsigz {
                        s.push_str(&pe(f[is as usize], 12, 4));
                    }
                }
                c.io.so(&s);
            }
        }
    }
}

/// Copy tape n1 if present and n2 to n3 (`20a4` records: the first 80
/// characters of each line, blank padded).
fn mergt(c: &mut Ctx, ifiss: i32, n1: i32, n2: i32, n3: i32) {
    c.repoz(n2);
    c.repoz(n3);
    if ifiss > 0 {
        c.repoz(n1);
    }
    let copy = |c: &mut Ctx, nin: i32| {
        while !c.io.line_eof(nin) {
            let r = c.io.read_line(nin);
            let mut s: Vec<u8> = r.iter().take(80).copied().collect();
            while s.len() < 80 {
                s.push(b' ');
            }
            c.io.write_line(n3, &String::from_utf8_lossy(&s));
        }
    };
    if ifiss != 0 {
        copy(c, n1);
    }
    copy(c, n2);
    c.repoz(n3);
}
