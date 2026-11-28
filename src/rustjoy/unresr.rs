//! UNRESR: self-shielded cross sections in the unresolved resonance range
//! (Bondarenko method, ETOX averages), ported from NJOY2016 `unresr.f90`
//! (module unresm).
//!
//! The routines keep the Fortran structure: module variables are fields of
//! `Unresr`, `save`d locals are fields too, and labels become loops or
//! `match` states where the control flow needs them. Arrays keep the
//! Fortran (1-based) indices.

#![allow(clippy::too_many_arguments)]
#![allow(clippy::collapsible_if)]
#![allow(clippy::collapsible_else_if)]
#![allow(clippy::needless_range_loop)]

use crate::endf::Ctx;
use crate::farr::{RM, RV};
use crate::ffmt;
use crate::fio::{Kind, Val};
use crate::fmath::{atan, cos, exp, int, log, nint, pow, powc, sin, sqrt};
use crate::util::physics::{AMASSN, AMU, BK, EV, HBAR, PI};
use crate::util::sigfig;

const NMTX: usize = 5;
const MAXSCR: usize = 1000;
const MAXEUNR: usize = 150;
const MAXARRY: i32 = 10000;
/// Extra words past the Fortran array sizes, so that writes the Fortran
/// does past an array end (without checking) do not stop the port first.
const SLACK: usize = 20000;
/// Room for more sections than the Fortran arrays (20) hold.
const NSECT: usize = 200;

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

/// Format `(..., 1p, e10.3 / (40x, e10.3))` for a list of values. A list
/// of one value leaves an empty record (the slash is done before the list
/// is found to be empty).
fn list40(head: &str, v: &[f64]) -> String {
    if v.is_empty() {
        return head.to_string();
    }
    let mut s = format!("{}{}", head, ffmt::e1p(v[0], 10, 3));
    for x in &v[1..] {
        s.push('\n');
        s.push_str(&" ".repeat(40));
        s.push_str(&ffmt::e1p(*x, 10, 3));
    }
    if v.len() == 1 {
        s.push('\n');
    }
    s
}

/// UNRESR module variables (and the `save`d local of unresl).
pub struct Unresr {
    nendf: i32,
    nin: i32,
    nout: i32,
    nscr: i32,
    nunr: i32,
    matd: i32,
    nsigz: i32,
    mtx: [i32; NMTX + 1],
    npnts: [i32; NMTX + 1],
    intunr: i32,
    lssf: i32,
    nro: i32,
    naps: i32,
    isot: [i32; NSECT + 1],
    modet: [i32; NSECT + 1],
    ibaset: [i32; NSECT + 1],
    nsect: i32,
    abnt: [f64; NSECT + 1],
    elt: [f64; NSECT + 1],
    eht: [f64; NSECT + 1],
    /// tr(62,62), ti(62,62), column-major.
    tr: Vec<f64>,
    ti: Vec<f64>,
    /// `noissue` in unresl (initialized, so saved).
    noissue: bool,
    // locals of unresr kept between materials (the Fortran keeps them
    // because the material loop stays inside one call)
    iprint: i32,
    ntemp: i32,
}

impl Default for Unresr {
    fn default() -> Self {
        Self::new()
    }
}

impl Unresr {
    pub fn new() -> Self {
        Unresr {
            nendf: 0,
            nin: 0,
            nout: 0,
            nscr: 0,
            nunr: 0,
            matd: 0,
            nsigz: 0,
            mtx: [0; NMTX + 1],
            npnts: [0; NMTX + 1],
            intunr: 0,
            lssf: 0,
            nro: 0,
            naps: 0,
            isot: [0; NSECT + 1],
            modet: [0; NSECT + 1],
            ibaset: [0; NSECT + 1],
            nsect: 0,
            abnt: [0.0; NSECT + 1],
            elt: [0.0; NSECT + 1],
            eht: [0.0; NSECT + 1],
            tr: vec![0.0; 62 * 62],
            ti: vec![0.0; 62 * 62],
            noissue: true,
            iprint: 0,
            ntemp: 0,
        }
    }

    #[inline]
    fn tr(&self, i: i32, j: i32) -> f64 {
        self.tr[((i - 1) + (j - 1) * 62) as usize]
    }
    #[inline]
    fn ti(&self, i: i32, j: i32) -> f64 {
        self.ti[((i - 1) + (j - 1) * 62) as usize]
    }

    // ============================================================== unresr
    /// Compute unresolved resonance cross sections.
    pub fn unresr(&mut self, c: &mut Ctx) {
        let (mut nb, mut nw) = (0i32, 0i32);
        let mut za = 0.0f64;
        let mut awr = 0.0f64;
        let mut bkgz = [0.0f64; 5];
        let mut scr = RV::new(MAXSCR + SLACK);
        let mut eunr = RV::new(0);
        let mut arry = RV::new(0);
        let mut sb = RM::new(0, 4);
        let mut temp: Vec<f64> = Vec::new();
        let mut sigz: Vec<f64> = Vec::new();
        let mut sigu: RM;

        //--initialize
        self.nscr = 10;
        c.nsh = 0;
        self.mtx[1] = 1;
        self.mtx[2] = 2;
        self.mtx[3] = 18;
        self.mtx[4] = 19;
        self.mtx[5] = 102;
        let time = c.timer();
        let v = c.io.input.read(&[Kind::Int, Kind::Int, Kind::Int]);
        self.nendf = geti(&v[0], self.nendf);
        self.nin = geti(&v[1], self.nin);
        self.nout = geti(&v[2], self.nout);
        if self.nin < 0 && self.nout > 0 {
            c.error("unresr", "mode conversion between nin and nout not allowed.", " ");
        }
        if self.nin > 0 && self.nout < 0 {
            c.error("unresr", "mode conversion between nin and nout not allowed.", " ");
        }
        if self.nin < 0 {
            self.nscr = -self.nscr;
        }
        let (nendf, nin, nout, nscr) = (self.nendf, self.nin, self.nout, self.nscr);
        c.openz(nendf, 0);
        c.openz(nin, 0);
        c.openz(nout, 1);
        c.openz(nscr, 1);
        c.io.so(&format!(
            "\n unresr...calculation of unresolved resonance cross sections{}{}s",
            " ".repeat(9),
            ffmt::f(time, 8, 1)
        ));
        c.io.se(&format!("\n unresr...{}{}s", " ".repeat(59), ffmt::f(time, 8, 1)));
        c.io.so(&format!(
            "\n unit for input endf tape ............. {}\n unit for input pendf tape ............ {}\n unit for output pendf tape ........... {}",
            ffmt::i(nendf, 10),
            ffmt::i(nin, 10),
            ffmt::i(nout, 10)
        ));
        self.uwtab();
        c.repoz(nin);
        c.repoz(nendf);
        c.repoz(nout);
        c.tpidio(nin, nout, 0, &mut scr.all(), &mut nb, &mut nw);
        c.tpidio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);

        //--loop over requested materials
        'mat: loop {
            self.iprint = 0;
            self.ntemp = 1;
            self.nsigz = 1;
            let v = c.io.input.read(&[Kind::Int; 4]);
            self.matd = geti(&v[0], self.matd);
            self.ntemp = geti(&v[1], self.ntemp);
            self.nsigz = geti(&v[2], self.nsigz);
            self.iprint = geti(&v[3], self.iprint);
            if self.matd == 0 {
                break 'mat;
            }
            let (matd, ntemp, nsigz, iprint) = (self.matd, self.ntemp, self.nsigz, self.iprint);
            let nt = ntemp.max(0) as usize;
            let ns = nsigz.max(0) as usize;
            // values not given keep the old ones (zero on a new array)
            let mut t2 = vec![0.0f64; nt];
            let kinds: Vec<Kind> = (0..nt).map(|_| Kind::Real).collect();
            let v = c.io.input.read(&kinds);
            for i in 0..nt {
                t2[i] = getr(&v[i], if i < temp.len() { temp[i] } else { 0.0 });
            }
            temp = t2;
            let mut s2 = vec![0.0f64; ns];
            let kinds: Vec<Kind> = (0..ns).map(|_| Kind::Real).collect();
            let v = c.io.input.read(&kinds);
            for i in 0..ns {
                s2[i] = getr(&v[i], if i < sigz.len() { sigz[i] } else { 0.0 });
            }
            sigz = s2;
            sigu = RM::new(5, ns.max(1));
            c.io.so(&list40("\n temperatures ......................... ", &temp));
            c.io.so(&list40(" sigma zero values .................... ", &sigz));
            c.io.so(&format!(" print option (0 min., 1 max.) ........ {}", ffmt::i(iprint, 10)));

            // 0: normal, 1: no resonance parameters, 2: no unresolved
            let mut skip = 0;

            //--loop over requested temperatures
            for it in 1..=ntemp {
                let temz = temp[(it - 1) as usize];
                let time = c.timer();
                c.io.so(&format!(
                    "\n mat = {}    temp = {}{}{}s",
                    ffmt::i(matd, 4),
                    ffmt::e1p(temz, 10, 3),
                    " ".repeat(37),
                    ffmt::f(time, 8, 1)
                ));

                //--read resonance parameters from file 2 on endf tape.
                //--read background cross sections from file 3
                for i in 1..=NMTX {
                    self.npnts[i] = 0;
                }
                if it <= 1 {
                    self.nunr = 0;
                    c.findf(matd, 1, 451, nendf);
                    c.contio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                    let lrp = c.l1h;
                    if lrp == 0 {
                        skip = 1;
                        break;
                    }
                    c.contio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                    if c.n1h != 0 {
                        c.iverf = 4;
                    } else if c.n2h == 0 {
                        c.iverf = 5;
                    } else {
                        c.iverf = 6;
                    }
                    c.findf(matd, 2, 151, nendf);
                    c.contio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                    za = c.c1h;
                    awr = c.c2h;
                    eunr = RV::new(MAXEUNR + SLACK);
                    arry = RV::new(MAXARRY as usize + SLACK);
                    self.rdunf2(c, &mut eunr, &mut scr, &mut arry, MAXARRY);
                    if self.nunr == 0 {
                        skip = 2;
                        break;
                    }
                    if self.lssf <= 0 {
                        c.findf(matd, 3, 1, nendf);
                        sb = RM::new(self.nunr as usize, 4);
                        self.rdunf3(c, &eunr, &mut sb, &mut scr);
                    }
                }
                // label 130
                let nunr = self.nunr;
                let nbb = 12 + nsigz + nunr * (1 + 5 * nsigz);
                let mut b = RV::new(nbb.max(1) as usize + 64);
                let mut l = 1i32;
                b[l] = za;
                b[l + 1] = awr;
                b[l + 2] = self.lssf as f64;
                b[l + 3] = 0.0;
                b[l + 4] = 0.0;
                b[l + 5] = self.intunr as f64;
                b[l + 6] = temz;
                b[l + 7] = 0.0;
                b[l + 8] = 5.0;
                b[l + 9] = nsigz as f64;
                b[l + 10] = (nsigz + nunr * (1 + 5 * nsigz)) as f64;
                b[l + 11] = nunr as f64;
                l += 11;
                // store sigz into next nsigz location of a.
                for i in 0..ns {
                    l += 1;
                    b[l] = sigz[i];
                }

                //--compute unresolved resonance cross-sections
                //--for all grid points and values of sigma zero.
                for ie in 1..=nunr {
                    let ez = eunr[ie];
                    l += 1;
                    b[l] = ez;
                    if self.lssf <= 0 {
                        bkgz[1] = sb[(ie, 1)];
                        bkgz[2] = sb[(ie, 2)];
                        bkgz[3] = sb[(ie, 3)];
                        bkgz[4] = sb[(ie, 4)];
                    } else {
                        bkgz[1] = 0.0;
                        bkgz[2] = 0.0;
                        bkgz[3] = 0.0;
                        bkgz[4] = 0.0;
                    }
                    if iprint == 1 {
                        c.io.so(&format!(" energy ={}", ffmt::e1p(ez, 12, 4)));
                    }
                    self.unresl(c, ez, temz, &sigz, nsigz, &bkgz, &mut sigu, &mut arry);
                    for ix in 1..=5 {
                        for is in 1..=nsigz {
                            l += 1;
                            if l as usize <= b.0.len() {
                                b[l] = sigfig(sigu[(ix, is)], 7, 0);
                            }
                        }
                        if l > nbb {
                            c.error("unresr", "storage exceeded.", " ");
                        }
                        if iprint == 1 {
                            let mut s = String::new();
                            for j in 1..=nsigz {
                                if j > 1 && (j - 1) % 8 == 0 {
                                    s.push('\n');
                                }
                                if (j - 1) % 8 == 0 {
                                    s.push(' ');
                                }
                                s.push_str(&ffmt::e1p(sigu[(ix, j)], 11, 3));
                            }
                            c.io.so(&s);
                        }
                    }
                }
                if nin != 0 && nout != 0 {
                    //--check for previous mt152
                    let mut new = 1;
                    c.findf(matd, 2, 151, nin);
                    c.tosend(nin, 0, 0, &mut scr.all());
                    c.contio(nin, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                    if c.mth == 152 {
                        new = 0;
                    }
                    c.findf(matd, 1, 451, nin);

                    //--write new pendf tape
                    let mut nx = 0i32;
                    c.contio(nin, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                    if c.iverf < 5 {
                        nx = c.n2h;
                        if nx > 0 && new > 0 {
                            scr[6] += 1.0;
                        }
                    }
                    c.contio(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                    if c.iverf >= 5 {
                        c.contio(nin, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                    }
                    if c.iverf == 6 {
                        c.contio(nin, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                    }
                    c.hdatio(nin, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                    scr[1] = temz;
                    if c.iverf >= 5 {
                        nx = c.n2h;
                        if nx > 0 && new > 0 {
                            scr[6] += 1.0;
                        }
                    }
                    c.hdatio(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                    while nb != 0 {
                        c.moreio(nin, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                    }
                    nw = nx;
                    if nw != 0 {
                        let ncc = 6 * nx + 6;
                        let mut cc = RV::new(ncc.max(1) as usize + 64);
                        let nd = 6 * nx;
                        let mut d = RV::new(nd.max(1) as usize + 64);
                        c.dictio(nin, 0, 0, &mut d.all(), &mut nb, &mut nw);
                        let mut j = 0i32;
                        let mut newmat = 0;
                        let ncds = 2 + (nint(b[11]) - 1) / 6;
                        let mut i = 1;
                        while i <= nw {
                            let mfd = nint(d[i + 2]);
                            let mtd = nint(d[i + 3]);
                            if mfd > 2 || mtd == 152 {
                                if newmat <= 0 {
                                    newmat = 1;
                                    cc[j + 1] = 0.0;
                                    cc[j + 2] = 0.0;
                                    cc[j + 3] = 2.0;
                                    cc[j + 4] = 152.0;
                                    cc[j + 5] = ncds as f64;
                                    cc[j + 6] = 0.0;
                                    if new > 0 {
                                        j += 6;
                                    }
                                }
                            }
                            for k in 1..=6 {
                                cc[k + j] = d[k + i - 1];
                            }
                            j += 6;
                            i += 6;
                        }
                        let mut nwd = j / 6;
                        c.dictio(0, nout, 0, &mut cc.all(), &mut nb, &mut nwd);
                    }
                    // copy through file 2, mt 151
                    c.tofend(nin, nout, 0, &mut scr.all());
                    c.tosend(nin, nout, 0, &mut scr.all());
                    if new == 0 {
                        c.tosend(nin, 0, 0, &mut scr.all());
                    }
                    // write new mt on output tape
                    c.math = matd;
                    c.mfh = 2;
                    c.mth = 152;
                    c.contio(0, nout, 0, &mut b.all(), &mut nb, &mut nw);
                    let mut lstart = 7;
                    let mut nwds = l - 6;
                    nw = nwds;
                    if nw > c.npage + 6 {
                        nw = c.npage + 6;
                    }
                    c.listio(0, nout, 0, &mut b.at(lstart), &mut nb, &mut nw);
                    loop {
                        lstart += nw;
                        nwds -= nw;
                        nw = nwds;
                        if nw > c.npage {
                            nw = c.npage;
                        }
                        if nb == 0 {
                            break;
                        }
                        c.moreio(0, nout, 0, &mut b.at(lstart), &mut nb, &mut nw);
                    }
                    c.asend(nout, 0);
                    // copy file end card and rest of material
                    c.tomend(nin, nout, 0, &mut scr.all());
                }

                //--write report of calculation
                let time = c.timer();
                c.io.so(&format!(
                    " generated cross sections at {} points{}{}s",
                    ffmt::i(nunr, 3),
                    " ".repeat(30),
                    ffmt::f(time, 8, 1)
                ));
            }
            if skip == 0 {
                continue 'mat;
            }

            if skip == 1 {
                //--mat has no resonance parameters. copy as is to nout.
                let strng1 = format!("mat{} has no resonance parameters", ffmt::i(matd, 5));
                c.mess("unresr", &strng1, "copy as is to nout");
            } else {
                //--mat has no unresolved.  copy as is to nout
                let strng1 = format!("mat{} has no unresolved parameters", ffmt::i(matd, 5));
                c.mess("unresr", &strng1, "copy as is to nout");
            }
            // label 355
            c.findf(matd, 1, 451, nin);
            let mut it = 0;
            loop {
                c.contio(nin, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                if c.math != matd {
                    continue 'mat;
                }
                c.contio(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                c.tomend(nin, nout, 0, &mut scr.all());
                it += 1;
                if it >= ntemp {
                    break;
                }
            }
        }

        //--finished with unresr.
        c.atend(nout, 0);
        c.repoz(nout);
        c.repoz(nin);
        c.repoz(nendf);
        c.closz(nendf);
        c.closz(nin);
        c.closz(nout);
        c.closz(nscr);
        let time = c.timer();
        c.io.so(&format!("\n{}{}s\n {}*******", " ".repeat(69), ffmt::f(time, 8, 1), "**********".repeat(7)));
    }

    // ============================================================== rdunf2
    /// Read and/or copy unresolved resonance parameters from file 2.
    fn rdunf2(&mut self, c: &mut Ctx, eunr: &mut RV, scr: &mut RV, arry: &mut RV, jx: i32) {
        const EGRIDU: [f64; 78] = [
            1.0e1, 1.25e1, 1.5e1, 1.7e1, 2.0e1, 2.5e1, 3.0e1, 3.5e1, 4.0e1, 5.0e1, 6.0e1, 7.2e1, 8.5e1, 1.0e2,
            1.25e2, 1.5e2, 1.7e2, 2.0e2, 2.5e2, 3.0e2, 3.5e2, 4.0e2, 5.0e2, 6.0e2, 7.2e2, 8.5e2, 1.0e3, 1.25e3,
            1.5e3, 1.7e3, 2.0e3, 2.5e3, 3.0e3, 3.5e3, 4.0e3, 5.0e3, 6.0e3, 7.2e3, 8.5e3, 1.0e4, 1.25e4, 1.5e4,
            1.7e4, 2.0e4, 2.5e4, 3.0e4, 3.5e4, 4.0e4, 5.0e4, 6.0e4, 7.2e4, 8.5e4, 1.0e5, 1.25e5, 1.5e5, 1.7e5,
            2.0e5, 2.5e5, 3.0e5, 3.5e5, 4.0e5, 5.0e5, 6.0e5, 7.2e5, 8.5e5, 1.0e6, 1.25e6, 1.5e6, 1.7e6, 2.0e6,
            2.5e6, 3.0e6, 3.5e6, 4.0e6, 5.0e6, 6.0e6, 7.2e6, 8.5e6,
        ];
        const NGRIDU: usize = 78;
        const WIDE: f64 = 1.26e0;
        const ONEMEV: f64 = 1.0e6;
        const SMALL: f64 = 1.0e-12;
        const STEP: f64 = 1.01e0;
        const ZERO: f64 = 0.0;
        let (mut nb, mut nw) = (0i32, 0i32);
        let nendf = self.nendf;
        let mut enow = 0.0f64;

        //--initialize
        let mut indep = 0;
        let mut elr = 0.0f64;
        let mut ehr = ONEMEV;
        self.nunr = 1;
        eunr[1] = ONEMEV;

        //--head card read in calling program
        let nis = c.n1h;
        self.nsect = 0;
        let mut inow = 1i32;

        //--loop over isotopes
        for is in 1..=nis {
            c.contio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
            let abn = c.c2h;
            let lfw = c.l2h;
            let ner = c.n1h;
            let mut ier = 0;

            //--loop over all energy ranges
            loop {
                // label 110
                c.contio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                ier += 1;
                let el = sigfig(c.c1h, 7, 0);
                let eh = sigfig(c.c2h, 7, 0);
                let lru = c.l1h;
                let lrf = c.l2h;
                self.nro = c.n1h;
                self.naps = c.n2h;
                if lru != 2 {
                    if eh > elr {
                        elr = eh;
                    }

                    //--read through this sub-section
                    c.contio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                    let nls = c.n1h;
                    if lrf == 4 || lrf == 7 {
                        c.listio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                        while nb != 0 {
                            c.moreio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                        }
                    }
                    for _l in 1..=nls {
                        if lrf != 4 {
                            c.listio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                            while nb != 0 {
                                c.moreio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                            }
                            if lrf == 7 {
                                c.listio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                                while nb != 0 {
                                    c.moreio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                                }
                            }
                        } else {
                            c.contio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                            let njs = c.n1h;
                            for _j in 1..=njs {
                                c.listio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                                while nb != 0 {
                                    c.moreio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                                }
                            }
                        }
                    }
                    if ier < ner {
                        continue;
                    }
                    break;
                }

                //--read and copy unresolved resonance parameters
                // label 140
                self.nsect += 1;
                let ns = self.nsect as usize;
                self.isot[ns] = is;
                self.abnt[ns] = abn;
                self.elt[ns] = el;
                self.eht[ns] = eh;
                let mode = lrf + 4 * (lru - 1);
                self.modet[ns] = mode;
                self.ibaset[ns] = inow;
                if eh < ehr {
                    ehr = eh;
                }
                // add el and eh to list of nodes
                // shade nodes around discontinuities
                let et = sigfig(el, 7, -1);
                ilist(et, eunr, &mut self.nunr);
                let et = sigfig(el, 7, 1);
                ilist(et, eunr, &mut self.nunr);
                let et = sigfig(eh, 7, -1);
                ilist(et, eunr, &mut self.nunr);
                let et = sigfig(eh, 7, 1);
                ilist(et, eunr, &mut self.nunr);

                // if present, read and store the energy-dependent
                // scattering radius tab1 data.
                if self.nro == 1 {
                    c.tab1io(nendf, 0, 0, &mut arry.at(inow), &mut nb, &mut nw);
                    let mut jj = inow + nw;
                    while nb != 0 {
                        c.moreio(nendf, 0, 0, &mut arry.at(jj), &mut nb, &mut nw);
                        jj += nw;
                        if jj > jx {
                            c.error("rdunf2", "storage in a exceeded", " ");
                        }
                    }
                    inow = jj;
                }

                //--branch to specified representation
                if lrf != 2 {
                    self.intunr = 2;
                    if lfw != 1 {
                        //--all parameters independent of energy
                        c.contio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                        let spi = c.c1h;
                        let ay = c.c2h;
                        self.lssf = c.l1h;
                        let nls = c.n1h;
                        // loop over l states
                        for l in 1..=nls {
                            c.listio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                            let awri = c.c1h;
                            let ll = c.l1h;
                            let njs = c.n2h;
                            if l == 1 {
                                arry[inow] = awri;
                                arry[inow + 1] = ay;
                                arry[inow + 2] = spi;
                                arry[inow + 3] = lfw as f64;
                                arry[inow + 4] = nls as f64;
                                arry[inow + 5] = 0.0;
                                inow += 6;
                            }
                            arry[inow] = njs as f64;
                            arry[inow + 1] = ll as f64;
                            inow += 2;
                            let mut jnow = 6;
                            for _n in 1..=njs {
                                for i in 1..=5 {
                                    arry[inow] = scr[jnow + i];
                                    inow += 1;
                                }
                                jnow += 6;
                            }
                        }
                        if inow > jx {
                            c.error("rdunf2", "storage exceeded.", " ");
                        }
                        indep = 1;
                    } else {
                        //--fission widths energy dependent
                        // label 190
                        c.listio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                        let spi = c.c1h;
                        let ay = c.c2h;
                        self.lssf = c.l1h;
                        let mut ne = c.n1h;
                        let nls = c.n2h;
                        arry[inow + 1] = ay;
                        arry[inow + 2] = spi;
                        arry[inow + 3] = lfw as f64;
                        arry[inow + 4] = nls as f64;
                        arry[inow + 5] = ne as f64;
                        // store fission width energy points
                        let inow1 = inow + 5;
                        for i in 1..=ne {
                            enow = sigfig(scr[i + 6], 7, 0);
                            arry[inow1 + i] = enow;
                            // add to list of energy nodes
                            if i > 1 && i <= ne {
                                ilist(enow, eunr, &mut self.nunr);
                            }
                        }
                        if enow < eh {
                            let strng = format!(
                                "between{} and{} eV",
                                ffmt::e1p(enow, 12, 4),
                                ffmt::e1p(eh, 12, 4)
                            );
                            c.error("rdunf2", "energy dependent data undefined", &strng);
                        }
                        // loop over l states
                        for l in 1..=nls {
                            c.contio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                            if l == 1 {
                                arry[inow] = c.c1h;
                                inow = inow + 6 + ne;
                            }
                            let njs = c.n1h;
                            arry[inow] = njs as f64;
                            arry[inow + 1] = c.l1h as f64;
                            inow += 2;
                            // loop over j states
                            for _j in 1..=njs {
                                c.listio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                                ne = c.n1h - 6;
                                arry[inow] = c.l2h as f64;
                                // store parameters
                                for ir in 1..=5 {
                                    arry[inow + ir] = scr[6 + ir];
                                }
                                let inow1 = inow + 5;
                                // store fission widths
                                for k in 1..=ne {
                                    arry[k + inow1] = scr[k + 12];
                                }
                                inow = inow + 6 + ne;
                            }
                        }
                        if inow > jx {
                            c.error("rdunf2", "storage exceeded.", " ");
                        }
                    }
                } else {
                    //--all parameters energy dependent
                    // label 250
                    let ist = inow;
                    // first cont record
                    c.contio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                    self.lssf = c.l1h;
                    let nls = c.n1h;
                    arry[inow + 1] = c.c2h;
                    arry[inow + 2] = c.c1h;
                    arry[inow + 3] = c.n1h as f64;
                    inow += 4;
                    // loop over l states
                    for l in 1..=nls {
                        c.contio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                        if l == 1 {
                            arry[ist] = c.c1h;
                        }
                        let njs = c.n1h;
                        arry[inow] = njs as f64;
                        arry[inow + 1] = c.l1h as f64;
                        inow += 2;
                        // loop over j states
                        for j in 1..=njs {
                            c.listio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                            let mut loc = 1 + nw;
                            while nb != 0 {
                                c.moreio(nendf, 0, 0, &mut scr.at(loc), &mut nb, &mut nw);
                                loc += nw;
                            }
                            let ne = c.n2h;
                            arry[inow] = c.c1h;
                            arry[inow + 1] = c.l1h as f64;
                            self.intunr = c.l1h;
                            arry[inow + 2] = ne as f64;
                            // store numbers of degrees of freedom
                            for k in 3..=6 {
                                arry[k + inow] = scr[k + 6];
                            }
                            inow += 7;
                            let mut jnow = 12;
                            // store parameters
                            for n in 1..=ne {
                                for k in 1..=6 {
                                    if scr[jnow + k] <= ZERO {
                                        scr[jnow + k] = SMALL;
                                    }
                                    arry[inow + k - 1] = scr[jnow + k];
                                }
                                inow += 6;
                                // add to list of energy nodes
                                enow = sigfig(scr[jnow + 1], 7, 0);
                                if n != 1 && n != ne && l == 1 && j == 1 {
                                    ilist(enow, eunr, &mut self.nunr);
                                }
                                jnow += 6;
                            }
                            if enow < eh {
                                let strng = format!(
                                    "between{} and{} eV",
                                    ffmt::e1p(enow, 12, 4),
                                    ffmt::e1p(eh, 12, 4)
                                );
                                c.error("rdunf2", "energy dependent data undefined", &strng);
                            }
                        }
                    }
                    if inow > jx {
                        c.error("rdunf2", "storage exceeded.", " ");
                    }
                }

                //--continue loop over isotopes and energy ranges
                // label 290
                if ier < ner {
                    continue;
                }
                break;
            }
        }
        // label 300
        self.nunr -= 1;
        if self.nunr == 0 {
            return;
        }

        //--add extra nodes if needed
        //--remove first and last energy nodes
        //--flag resolved-unresolved overlap energies
        self.nunr += 1;
        let mut i = 1i32;
        let mut elast = eunr[2];
        loop {
            // label 350
            i += 1;
            let enext = eunr[i + 1];
            if enext >= ONEMEV {
                break;
            }
            if !(enext < WIDE * elast && indep == 0) {
                let mut et = elast;
                loop {
                    // label 360
                    let mut enut = enext;
                    for ii in 0..NGRIDU {
                        if EGRIDU[ii] > STEP * et {
                            enut = EGRIDU[ii];
                            break;
                        }
                    }
                    et = enut;
                    if et >= enext {
                        break;
                    }
                    ilist(et, eunr, &mut self.nunr);
                    i += 1;
                }
            }
            // label 370
            elast = eunr[i + 1];
        }
        // label 380
        let mut iovlp = 0;
        let lim = self.nunr - 1;
        self.nunr = 0;
        let mut en = 0.0f64;
        for i in 2..=lim {
            let mut et = eunr[i];
            if et >= en {
                if et < elr {
                    et = -et;
                }
                if et < ZERO {
                    iovlp = 1;
                }
                self.nunr += 1;
                eunr[self.nunr] = et;
                en = sigfig(et.abs(), 7, 2);
            }
        }
        self.nunr -= 1;
        if iovlp == 1 {
            c.mess("rdunf2", "resolved-unresolved overlap energies", "marked with minus signs.");
        }
        if eunr[self.nunr] > ehr {
            let strng = format!("unresolved-smooth overlap above e={}", ffmt::e1p(ehr, 12, 4));
            c.mess("rdunf2", &strng, " ");
        }
    }

    // ============================================================== rdunf3
    /// Read unresolved region background cross sections from file 3.
    /// Compute the background cross sections on the unresolved energy grid
    /// and store for use in unresl.
    fn rdunf3(&mut self, c: &mut Ctx, eunr: &RV, sb: &mut RM, scr: &mut RV) {
        const UP: f64 = 1.00001e0;
        const DN: f64 = 0.99999e0;
        let (mut nb, mut nw) = (0i32, 0i32);
        let nendf = self.nendf;
        let nunr = self.nunr;
        let mut enext = 0.0f64;
        let mut idis = 0i32;
        let mut bkg = 0.0f64;
        let mut sig = 0.0f64;

        let mut ix = 1usize;
        // 100: read a cont record; 105: test it; 160: next reaction
        let mut lab = 100;
        loop {
            match lab {
                100 => {
                    c.contio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                    lab = 105;
                }
                105 => {
                    let mut kx = ix as i32;
                    if self.mtx[ix] > 18 {
                        kx = ix as i32 - 1;
                    }
                    if c.mth < self.mtx[ix] {
                        c.tosend(nendf, 0, 0, &mut scr.all());
                        lab = 100;
                        continue;
                    }
                    if c.mth == self.mtx[ix] {
                        // label 130
                        if self.mtx[ix] == 19 {
                            self.npnts[ix] = self.npnts[ix - 1];
                        } else {
                            let mut e = 0.0f64;
                            c.gety1(e, &mut enext, &mut idis, &mut bkg, nendf, &mut scr.all());
                            let np = nint(scr[6]);
                            self.npnts[ix] = np;
                            for ie in 1..=nunr {
                                e = eunr[ie].abs();
                                if ie == 1 {
                                    e *= UP;
                                }
                                if ie == nunr {
                                    e *= DN;
                                }
                                c.gety1(e, &mut enext, &mut idis, &mut sig, nendf, &mut scr.all());
                                sb[(ie, kx)] = sig;
                            }
                            c.tosend(nendf, 0, 0, &mut scr.all());
                        }
                    } else {
                        for ie in 1..=nunr {
                            sb[(ie, kx)] = 0.0;
                        }
                    }
                    lab = 160;
                }
                _ => {
                    // label 160
                    ix += 1;
                    if ix > NMTX {
                        return;
                    }
                    if self.mtx[ix] == c.mth {
                        lab = 105;
                    } else if self.mtx[ix] > c.mth {
                        lab = 100;
                    } else {
                        lab = 160;
                    }
                }
            }
        }
    }

    // ============================================================== unresl
    /// Computes effective cross-sections in the unresolved region using
    /// the etox method.
    fn unresl(
        &mut self,
        c: &mut Ctx,
        ee: f64,
        tt: f64,
        sig0: &[f64],
        nsig0: i32,
        sigbkg: &[f64; 5],
        sigu: &mut RM,
        arry: &mut RV,
    ) {
        const MXKS: i32 = 100;
        const NQP: i32 = 10;
        #[rustfmt::skip]
        const QW: [f64; 40] = [
            1.1120413e-1, 2.3546798e-1, 2.8440987e-1,
            2.2419127e-1, 0.10967668e0, 0.030493789e0,
            0.0042930874e0, 2.5827047e-4, 4.9031965e-6,
            1.4079206e-8, 0.033773418e0, 0.079932171e0,
            0.12835937e0, 0.17652616e0, 0.21347043e0,
            0.21154965e0, 0.13365186e0, 0.022630659e0,
            1.6313638e-5, 2.745383e-31, 3.3376214e-4,
            0.018506108e0, 0.12309946e0, 0.29918923e0,
            0.33431475e0, 0.17766657e0, 0.042695894e0,
            4.0760575e-3, 1.1766115e-4, 5.0989546e-7,
            1.7623788e-3, 0.021517749e0, 0.080979849e0,
            0.18797998e0, 0.30156335e0, 0.29616091e0,
            0.10775649e0, 2.5171914e-3, 8.9630388e-10,
            0.0,
        ];
        #[rustfmt::skip]
        const QP: [f64; 40] = [
            3.0013465e-3, 7.8592886e-2, 4.3282415e-1,
            1.3345267e0, 3.0481846e0, 5.8263198e0,
            9.9452656e0, 1.5782128e1, 23.996824e0,
            36.216208e0, 1.3219203e-2, 7.2349624e-2,
            0.19089473e0, 0.39528842e0, 0.74083443e0,
            1.3498293e0, 2.5297983e0, 5.2384894e0,
            13.821772e0, 75.647525e0, 1.0004488e-3,
            0.026197629e0, 0.14427472e0, 0.44484223e0,
            1.0160615e0, 1.9421066e0, 3.3150885e0,
            5.2607092e0, 7.9989414e0, 12.072069e0,
            0.013219203e0, 0.072349624e0, 0.19089473e0,
            0.39528842e0, 0.74083443e0, 1.3498293e0,
            2.5297983e0, 5.2384894e0, 13.821772e0,
            75.647525e0,
        ];
        let qw = |i: i32, j: i32| QW[((i - 1) + (j - 1) * 10) as usize];
        let qp = |i: i32, j: i32| QP[((i - 1) + (j - 1) * 10) as usize];
        const RC1: f64 = 0.123e0;
        const RC2: f64 = 0.08e0;
        const THIRD: f64 = 0.333333333e0;
        const ZERO: f64 = 0.0;
        const ONE: f64 = 1.0e0;
        const SMALL: f64 = 1.0e-8;

        let cwaven = sqrt(2.0 * AMASSN * AMU * EV) * 1.0e-12 / HBAR;

        //--initialize
        let c2 = BK;
        let ac = 1.0f64;
        let e = ee.abs();
        let ntempu = 1i32;
        let mut tempu = [0.0f64; 2];
        tempu[1] = tt;
        if tempu[1] < ONE {
            tempu[1] = ONE;
        }

        //--make a pass to compute the potential scattering cross section
        //--then a second pass to compute unresolved cross sections
        let ns0 = nsig0.max(0) as usize;
        let mut ispot = 0;
        let mut spot = 0.0f64;
        let mut sint = 0.0f64;
        // sigf(ntempu,nsig0,5), t(4,nsig0,ntempu): ntempu is 1
        let mut sigf = RM::new(ns0.max(1), 5);
        let mut sigm = RV::new(ns0.max(1));
        let mut t = RM::new(4, ns0.max(1));
        let mut abns = RV::new(1);
        let mut tj = RM::new(1, 1);
        let mut tk = RM::new(1, 1);
        let mut tl = RM::new(1, 1);
        let mut sigbt = 0.0f64;
        let mut ks;
        let mut del = [0.0f64; 4];
        let mut gg = [0.0f64; 6];
        let mut yy = [0.0f64; 4];
        let (mut ay, mut enx, mut idx) = (0.0f64, 0.0f64, 0i32);
        let mut lfw = 0i32;
        let mut ne = 0i32;
        let mut iest = 0i32;
        let (mut xj, mut xk) = (0.0f64, 0.0f64);
        let (mut vl, mut ps) = (0.0f64, 0.0f64);
        let mut gfx = 0.0f64;
        let (mut dx, mut gxx, mut gnox, mut ggx) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        let (mut amux, mut amuf, mut aj, mut amun);

        loop {
            //--find sections of resonance parameters which contribute
            // label 200
            ks = 0i32;
            for i in 1..=self.nsect {
                let iu = i as usize;
                let isect = i;
                let abn = self.abnt[iu];
                if e < self.elt[iu] || e > self.eht[iu] {
                    continue;
                }
                let mode = self.modet[iu];
                if !(5..=6).contains(&mode) {
                    continue;
                }
                let mut inow = self.ibaset[isect as usize];

                //--retrieve nuclide information
                if self.nro == 1 {
                    let iro = inow;
                    let mut ip = 2;
                    let mut ir = 1;
                    c.terpa(&mut ay, e, &mut enx, &mut idx, &arry.at(iro), &mut ip, &mut ir);
                    inow = inow + 6 + 2 * nint(arry[iro + 4]) + 2 * nint(arry[iro + 5]);
                } else {
                    ay = arry[inow + 1];
                }
                let awri = arry[inow];
                let mut aaa = 0.0f64;
                if self.naps == 2 && self.nro == 1 {
                    aaa = arry[inow + 1];
                }
                let spi = arry[inow + 2];
                let nls;
                if mode == 6 {
                    nls = nint(arry[inow + 3]);
                    inow += 4;
                } else {
                    lfw = nint(arry[inow + 3]);
                    nls = nint(arry[inow + 4]);
                    ne = nint(arry[inow + 5]);
                    inow += 6;
                    if lfw != 0 {
                        // save starting location for fission width energies
                        iest = inow;
                        inow += ne;
                    }
                }

                //--calculate nuclide dependent parameters
                let rat = awri / (awri + 1.0);
                let aw = awri * AMASSN;
                let aa;
                if self.naps == 0 {
                    aa = RC1 * pow(aw, THIRD) + RC2;
                } else if self.naps == 1 {
                    aa = ay;
                } else if self.naps == 2 && self.nro == 1 {
                    aa = aaa;
                } else {
                    c.error("unresl", "illegal naps", " ");
                }
                let e2 = sqrt(e);
                let k = cwaven * rat * e2;
                let ab = 4.0 * PI / powc(k, 2);
                let rho = k * aa;
                let rhoc = k * ay;
                for itp in 1..=ntempu as usize {
                    del[itp] = 2.0 * e2 * sqrt(tempu[itp] * c2 / awri);
                }

                //--loop over all sequences(l,j) in this section
                for _l in 1..=nls {
                    let njs = nint(arry[inow]);
                    let ll = nint(arry[inow + 1]);

                    //--calculate penetrability and phase shifts
                    uunfac(ll, rho, rhoc, ONE, &mut vl, &mut ps);
                    inow += 2;
                    for j in 1..=njs {
                        //--retrieve parameters for this sequence
                        ks += 1;
                        if ks > MXKS {
                            c.error("unresl", "storage exceeded, increase mxks", " ");
                        }
                        if ispot != 0 {
                            abns[ks] = abn;
                        }
                        if mode != 6 {
                            amux = 0.0;
                            gxx = 0.0;
                            if lfw == 0 {
                                amuf = 0.0;
                                dx = arry[inow];
                                aj = arry[inow + 1];
                                amun = arry[inow + 2];
                                gnox = arry[inow + 3];
                                ggx = arry[inow + 4];
                                gfx = 0.0;
                                inow += 5;
                            } else {
                                amuf = arry[inow];
                                dx = arry[inow + 1];
                                aj = arry[inow + 2];
                                amun = arry[inow + 3];
                                gnox = arry[inow + 4];
                                ggx = arry[inow + 5];
                                inow += 6;
                                let ifst = inow;
                                let intl = 2;
                                intrf(c, e, &mut gfx, iest, ifst, ne, intl, arry);
                                inow += ne;
                            }
                        } else {
                            aj = arry[inow];
                            let _int = nint(arry[inow + 1]);
                            let int = 2;
                            ne = nint(arry[inow + 2]);
                            amux = arry[inow + 3];
                            amun = arry[inow + 4];
                            amuf = arry[inow + 6];
                            inow += 7;
                            intr(c, e, &mut dx, &mut gxx, &mut gnox, &mut ggx, &mut gfx, inow, ne, int, arry);
                            inow += 6 * ne;
                        }
                        let gj = (2.0 * aj + 1.0) / (4.0 * spi + 2.0);
                        let gnx = gnox * vl * e2 * amun;

                        //--compute potential scattering and interference correction
                        if ispot != 1 {
                            if j == 1 {
                                spot += abn * ab * (2 * ll + 1) as f64 * powc(sin(ps), 2);
                            }
                            sint -= abn * ab * PI * gj * gnx * powc(sin(ps), 2) / dx;

                        //--compute etox statistical averages for this sequence
                        } else {
                            for is0 in 1..=nsig0 {
                                tk[(ks, is0)] = 0.0;
                                for kx in 1..=4 {
                                    t[(kx, is0)] = 0.0;
                                }
                            }
                            let d = dx;
                            let mut mu = nint(amuf);
                            if gfx <= SMALL {
                                mu = 0;
                            }
                            let nu = nint(amun);
                            let mut lu = nint(amux);
                            if gxx < SMALL {
                                lu = 0;
                            }
                            gg[2] = ggx;
                            let nqf = if mu > 0 { NQP } else { 1 };
                            let nqn = if nu > 0 { NQP } else { 1 };
                            let nqx = if lu > 0 { NQP } else { 1 };
                            for kf in 1..=nqf {
                                if mu == 0 {
                                    gg[1] = gfx;
                                }
                                if mu > 0 {
                                    gg[1] = qp(kf, mu) * gfx;
                                }
                                for kn in 1..=nqn {
                                    if nu == 0 {
                                        gg[3] = gnx;
                                    }
                                    if nu > 0 {
                                        gg[3] = qp(kn, nu) * gnx;
                                    }
                                    for kl in 1..=nqx {
                                        if lu == 0 {
                                            gg[4] = gxx;
                                        }
                                        if lu > 0 {
                                            gg[4] = qp(kl, lu) * gxx;
                                        }
                                        gg[5] = gg[1] + gg[2] + gg[3] + gg[4];
                                        let s0u = ab * gj * gg[3] / gg[5];
                                        for itp in 1..=ntempu as usize {
                                            let sti = gg[5] / del[itp];
                                            for is0 in 1..=nsig0 {
                                                let beta = sigm[is0] / s0u;
                                                if beta < ZERO && beta.abs() < ONE && self.noissue {
                                                    c.mess(
                                                        "unresl",
                                                        "Square root of negative number detected",
                                                        "Probably caused by negative background xs in urr",
                                                    );
                                                    self.noissue = false;
                                                }
                                                self.ajku(beta, sti, &mut xj, &mut xk);
                                                if mu > 0 {
                                                    xj *= qw(kf, mu);
                                                }
                                                if mu > 0 {
                                                    xk *= qw(kf, mu);
                                                }
                                                if nu > 0 {
                                                    xj *= qw(kn, nu);
                                                }
                                                if nu > 0 {
                                                    xk *= qw(kn, nu);
                                                }
                                                if lu > 0 {
                                                    xj *= qw(kl, lu);
                                                }
                                                if lu > 0 {
                                                    xk *= qw(kl, lu);
                                                }
                                                for kx in 1..=4 {
                                                    t[(kx, is0)] += xj * gg[kx as usize];
                                                }
                                                tk[(ks, is0)] += xk * gg[5];
                                            }
                                        }
                                    }
                                }
                            }
                            for is0 in 1..=nsig0 {
                                tk[(ks, is0)] = ac * tk[(ks, is0)] / d;
                                let mut ttj = 0.0f64;
                                for kx in 1..=4 {
                                    ttj += t[(kx, is0)];
                                    let kkx = kx + (ks - 1) * 3;
                                    tl[(kkx, is0)] = ac * t[(kx, is0)] / d;
                                }
                                tj[(ks, is0)] = ac * ttj / d;
                            }
                        }
                    }
                }
            }
            if ispot != 0 {
                break;
            }
            // label 110
            //--use current value of ks to allocate remaining arrays
            let ksz = ks.max(0) as usize;
            abns = RV::new(ksz.max(1));
            tj = RM::new(ksz.max(1), ns0.max(1));
            tk = RM::new(ksz.max(1), ns0.max(1));
            tl = RM::new(3 * ksz + 1, ns0.max(1));
            sigbt = sigbkg[1] + spot + sint;
            for is in 1..=nsig0 {
                sigm[is] = sigbt + sig0[(is - 1) as usize];
                if sigm[is] < ZERO {
                    c.mess("unresl", "Negative background xs in urr may cause issues", "Check the evaluation");
                }
            }
            ispot = 1;
        }

        //--compute average cross-sections by summing over sequences
        let ns = ks;
        for is0 in 1..=nsig0 {
            yy[1] = 0.0;
            yy[2] = 0.0;
            yy[3] = 0.0;
            let mut yj = 0.0f64;
            let mut yk = 0.0f64;
            for ks in 1..=ns {
                xk = 0.0;
                xj = 0.0;
                for ksp in 1..=ns {
                    if ksp != ks {
                        xj += tj[(ksp, is0)] * abns[ksp];
                        xk += tk[(ksp, is0)] * abns[ksp];
                    }
                }
                for kx in 1..=3 {
                    let knm1 = kx + (ks - 1) * 3;
                    yy[kx as usize] += tl[(knm1, is0)] * (ONE - xj) * abns[ks];
                }
                let ttt = tj[(ks, is0)] * (ONE - xj) * abns[ks];
                yj += ttt;
                yk = yk + (tk[(ks, is0)] - tj[(ks, is0)]) * (ONE - xk) * abns[ks] + ttt;
            }
            sigf[(is0, 4)] = 0.0;
            for i in 1..=3 {
                sigf[(is0, i)] = sigm[is0] * yy[i as usize] / (ONE - yj);
                sigf[(is0, 4)] += sigf[(is0, i)];
            }
            let mut fact = sig0[(is0 - 1) as usize];
            let term = (ONE - yj) / (ONE - yk);
            let term1 = fact * (yk - yj) / (1.0 - yk);
            fact = sigbt;
            let sigtr = fact * term + term1;
            sigf[(is0, 5)] = sigtr;
        }

        //--output cross section in desired order.
        for is in 1..=nsig0 {
            sigu[(1, is)] = sigf[(is, 4)] + sigbt;
            sigu[(2, is)] = sigf[(is, 3)] + sigbkg[2] + spot + sint;
            sigu[(3, is)] = sigf[(is, 1)] + sigbkg[3];
            sigu[(4, is)] = sigf[(is, 2)] + sigbkg[4];
            sigu[(5, is)] = sigf[(is, 5)];
        }
    }

    // ================================================================ quikw
    /// Used to calculate chi and xi line shape functions.
    fn quikw(&self, ax: f64, y: f64, rew: &mut f64, aimw: &mut f64, ki: i32) {
        const BREAK1: f64 = 36.0e0;
        const BREAK2: f64 = 144.0e0;
        const BREAK3: f64 = 10000.0e0;
        const C1: f64 = 0.2752551e0;
        const C2: f64 = 2.724745e0;
        const C3: f64 = 0.5124242e0;
        const C4: f64 = 0.05176536e0;
        const C5: f64 = 1.1283792e0;
        const ZERO: f64 = 0.0;

        let rpi = sqrt(PI);
        let mut aki = 1.0f64;
        if ax < ZERO {
            aki = -1.0;
        }
        let x = ax.abs();
        let test = x * x + y * y;
        if test < BREAK1 {
            let ii = int(x * 10.0);
            let jj = int(y * 10.0);
            let i = ii + 2;
            let j = jj + 2;
            let n = j - 1;
            let p = 10.0 * x - ii as f64;
            let q = 10.0 * y - jj as f64;
            let p2 = p * p;
            let q2 = q * q;
            let pq = p * q;
            let hp = p / 2.0;
            let hq = q / 2.0;
            let hq2 = q2 / 2.0;
            let hp2 = p2 / 2.0;
            let a1 = hq2 - hq;
            let a2 = hp2 - hp;
            let a3 = 1.0 + pq - p2 - q2;
            let a4 = hp2 - pq + hp;
            let a5 = hq2 - pq + hq;
            *rew = a1 * self.tr(i, n)
                + a2 * self.tr(i - 1, j)
                + a3 * self.tr(i, j)
                + a4 * self.tr(i + 1, j)
                + a5 * self.tr(i, j + 1)
                + pq * self.tr(i + 1, j + 1);
            if ki > 0 {
                *aimw = a1 * self.ti(i, n)
                    + a2 * self.ti(i - 1, j)
                    + a3 * self.ti(i, j)
                    + a4 * self.ti(i + 1, j)
                    + a5 * self.ti(i, j + 1)
                    + pq * self.ti(i + 1, j + 1);
                *aimw *= aki;
            }
        } else if test < BREAK2 {
            let a1 = powc(x, 2) - powc(y, 2);
            let a2 = 2.0 * x * y;
            let a3 = powc(a2, 2);
            let a4 = a1 - C1;
            let a5 = a1 - C2;
            let d1 = C3 / (powc(a4, 2) + a3);
            let d2 = C4 / (powc(a5, 2) + a3);
            *rew = d1 * (a2 * x - a4 * y) + d2 * (a2 * x - a5 * y);
            if ki > 0 {
                *aimw = d1 * (a4 * x + a2 * y) + d2 * (a5 * x + a2 * y);
                *aimw *= aki;
            }
        } else if test < BREAK3 {
            let a1 = (powc(x, 2) - powc(y, 2)) * 2.0;
            let a2 = 4.0 * x * y;
            let a4 = a1 - 1.0;
            let d1 = C5 / (powc(a4, 2) + powc(a2, 2));
            *rew = d1 * (a2 * x - a4 * y);
            if ki > 0 {
                *aimw = d1 * (a4 * x + a2 * y);
                *aimw *= aki;
            }
        } else {
            let a1 = 1.0 / (rpi * test);
            *rew = y * a1;
            if ki > 0 {
                *aimw = x * a1 * aki;
            }
        }
    }

    // ================================================================= ajku
    /// Calculates j and k integrals.
    fn ajku(&self, b: f64, st: f64, aj: &mut f64, ak: &mut f64) {
        const XG: [f64; 8] = [
            0.095012510e0,
            0.28160355e0,
            0.45801678e0,
            0.61787624e0,
            0.75540441e0,
            0.8656312e0,
            0.94457502e0,
            0.98940093e0,
        ];
        const WG: [f64; 8] = [
            0.18945061e0,
            0.18260342e0,
            0.16915652e0,
            0.14959600e0,
            0.12462897e0,
            0.09515851e0,
            0.06225352e0,
            0.02715246e0,
        ];
        const EPS: f64 = 0.0002e0;
        const ZERO: f64 = 0.0;

        let pi2 = PI / 2.0;
        let sqpi = sqrt(PI);
        let term = if (-st / 2.0) < log(f64::MIN_POSITIVE) { 0.0 } else { exp(-st / 2.0) };
        let ep = (1.0 - term) / EPS - b;
        if ep <= ZERO {
            *aj = pi2 / b;
            *ak = 2.0 * *aj;
        } else {
            let mut y = 0.0f64;
            let mut z = 0.0f64;
            let mut yk = 0.0f64;
            let mut zk = 0.0f64;
            let y1 = st / 2.0;
            let ki = 0;
            let b1 = b / (sqpi * y1);
            let a = sqrt((1.0 + b) / b);
            let a2 = a * a;
            let c = 200.0 / (a * st);
            let remj = (pi2 - atan(c)) / (b * a);
            let remk = ((1.0 + a2) * remj - (1.0 - a2) * c / ((c * c + 1.0) * b * a)) / (2.0 * a2);
            let (mut rew, mut aimw) = (0.0f64, 0.0f64);
            for n in 0..8 {
                let ax = 5.0 * XG[n] + 5.0;
                self.quikw(ax, y1, &mut rew, &mut aimw, ki);
                let x1 = 1.0 / (1.0 + b1 / rew);
                let z1 = x1 * x1 / rew;
                let ax = -5.0 * XG[n] + 5.0;
                self.quikw(ax, y1, &mut rew, &mut aimw, ki);
                let x2 = 1.0 / (1.0 + b1 / rew);
                let z2 = x2 * x2 / rew;
                let ax = 45.0 * XG[n] + 55.0;
                self.quikw(ax, y1, &mut rew, &mut aimw, ki);
                let x3 = 1.0 / (1.0 + b1 / rew);
                let z3 = x3 * x3 / rew;
                let ax = -45.0 * XG[n] + 55.0;
                self.quikw(ax, y1, &mut rew, &mut aimw, ki);
                let x4 = 1.0 / (1.0 + b1 / rew);
                let z4 = x4 * x4 / rew;
                y += WG[n] * (x3 + x4);
                yk += WG[n] * (z3 + z4);
                z += WG[n] * (x1 + x2);
                zk += WG[n] * (z1 + z2);
            }
            *aj = 5.0 * (z + 9.0 * y) / y1 + remj;
            *ak = *aj + 5.0 * b1 * (zk + 9.0 * yk) / y1 + remk;
        }
    }

    // ================================================================ uwtab
    /// Generate table of complex probability integral w.
    fn uwtab(&mut self) {
        const TENTH: f64 = 0.1e0;
        let mut x = [0.0f64; 63];
        let mut y = [0.0f64; 63];
        let x0 = -TENTH;
        let y0 = -TENTH;
        let dx = TENTH;
        let dy = TENTH;
        let nx = 62usize;
        let ny = 62usize;
        y[1] = y0;
        let nx1 = nx - 1;
        let ny1 = ny - 1;
        for i in 1..=ny1 {
            y[i + 1] = i as f64 * dy + y0;
        }
        x[1] = x0;
        for i in 1..=nx1 {
            x[i + 1] = i as f64 * dx + x0;
        }
        for i in 1..=nx {
            for j in 1..=ny {
                let (rwt, aimt) = uw(x[i], y[j]);
                self.tr[(i - 1) + (j - 1) * 62] = rwt;
                self.ti[(i - 1) + (j - 1) * 62] = aimt;
            }
        }
    }
}

/// Insert a new energy into a list of energies in ascending order. Omit
/// duplicate values. Initial list must be primed with one energy larger
/// than any others to be added.
fn ilist(e: f64, elist: &mut RV, nlist: &mut i32) {
    let mut jl = 0;
    let n = *nlist;
    let mut placed = false;
    for i in 1..=n {
        if e > elist[i] {
            continue;
        }
        if e == elist[i] {
            placed = true;
            break;
        }
        jl = n - i + 1;
        for j in 1..=jl {
            elist[n - j + 2] = elist[n - j + 1];
        }
        elist[i] = e;
        placed = true;
        break;
    }
    if !placed {
        jl = 1;
        elist[n + 1] = e;
    }
    if jl > 0 {
        *nlist += 1;
    }
}

/// Calculates the penetrability factor (vl) and phase shift (ps).
fn uunfac(l: i32, rho: f64, rhoc: f64, amun: f64, vl: &mut f64, ps: &mut f64) {
    let r2 = rho * rho;
    if l == 0 {
        *vl = amun;
        *ps = rhoc;
    } else if l == 1 {
        *vl = amun * r2 / (1.0 + r2);
        *ps = rhoc - atan(rhoc);
    } else {
        let r4 = r2 * r2;
        *vl = amun * r4 / (9.0 + 3.0 * r2 + r4);
        *ps = rhoc - atan(3.0 * rhoc / (3.0 - rhoc * rhoc));
    }
}

/// Interpolates fission widths for unresolved representation 1.
fn intrf(c: &Ctx, e: f64, gfx: &mut f64, iest: i32, ifst: i32, ne: i32, int: i32, a: &RV) {
    for i in 2..=ne {
        let i1 = iest + i - 1;
        let i2 = ifst + i - 1;
        if e >= a[i1 - 1] && e <= a[i1] {
            c.terp1(a[i1 - 1], a[i2 - 1], a[i1], a[i2], e, gfx, int);
        }
    }
}

/// Interpolates energy dependent unresolved resonance parameters.
fn intr(
    c: &Ctx,
    e: f64,
    dx: &mut f64,
    gxx: &mut f64,
    gnox: &mut f64,
    ggx: &mut f64,
    gfx: &mut f64,
    inow: i32,
    ne: i32,
    int: i32,
    a: &RV,
) {
    const SMALL: f64 = 1.0e-8;
    let mut i1 = inow;
    for i in 2..=ne {
        let i2 = i1 + 6;
        if e >= a[i1] && e <= a[i2] {
            c.terp1(a[i1], a[i1 + 1], a[i2], a[i2 + 1], e, dx, int);
            c.terp1(a[i1], a[i1 + 2], a[i2], a[i2 + 2], e, gxx, int);
            c.terp1(a[i1], a[i1 + 3], a[i2], a[i2 + 3], e, gnox, int);
            c.terp1(a[i1], a[i1 + 4], a[i2], a[i2 + 4], e, ggx, int);
            c.terp1(a[i1], a[i1 + 5], a[i2], a[i2 + 5], e, gfx, int);
        }
        if i < ne {
            i1 = i2;
        }
    }
    if *gxx < SMALL {
        *gxx = 0.0;
    }
    if *gfx < SMALL {
        *gfx = 0.0;
    }
}

/// Compute complex probability integral w.
fn uw(rez: f64, aim1: f64) -> (f64, f64) {
    const C1: f64 = 1.1283792e0;
    const C2: f64 = 1.5e0;
    const UP: f64 = 1.0e15;
    const DN: f64 = 1.0e-15;
    const BRK1: f64 = 1.25e0;
    const BRK2: f64 = 5.0e0;
    const BRK3: f64 = 1.863636e0;
    const BRK4: f64 = 4.1e0;
    const BRK5: f64 = 1.71e0;
    const BRK6: f64 = 2.89e0;
    const BRK7: f64 = 1.18e0;
    const BRK8: f64 = 5.76e0;
    const BRK9: f64 = 1.5e0;
    const EPS: f64 = 1.0e-7;
    const ZERO: f64 = 0.0;

    let rpi = sqrt(PI);
    let mut rew = 0.0f64;
    let mut aimw = 0.0f64;
    let mut aimz = aim1.abs();
    let abrez = rez.abs();
    if abrez + aimz == ZERO {
        return (1.0, 0.0);
    }

    let r2 = rez * rez;
    let ai2 = aimz * aimz;
    let go350 = if abrez + BRK1 * aimz - BRK2 > ZERO {
        true
    } else if abrez + BRK3 * aimz - BRK4 > ZERO {
        true
    } else if r2 + BRK5 * ai2 - BRK6 < ZERO {
        false
    } else if r2 + BRK7 * ai2 - BRK8 >= ZERO {
        false
    } else {
        aimz - BRK9 >= ZERO
    };
    let kw;
    let mut lab;
    if !go350 {
        kw = 2;
        aimz = aim1;
        lab = 420;
    } else if aim1 >= ZERO {
        kw = 1;
        lab = 370;
    } else {
        kw = 2;
        aimz = aim1;
        lab = 420;
    }

    let (mut ak, mut el, mut h, mut b, mut a) = (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
    let (mut tempm, mut temel, mut g, mut c, mut d) = (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
    let (mut am, mut aak, mut rv) = (0.0f64, 0.0f64, 0.0f64);
    let (mut temp4, mut ajp, mut tempc, mut tempd) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
    let (mut temp1, mut aj, mut ajsig, mut sig2p) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
    let (mut expc, mut exps) = (0.0f64, 0.0f64);
    let mut _k = 0i32;
    let mut _jsig = 0i32;
    loop {
        match lab {
            370 => {
                // wa is obtained from asymtotic series
                rv = 2.0 * (r2 - ai2);
                ak = 4.0 * rez * aimz;
                el = ak;
                h = 0.0;
                b = 0.0;
                a = 0.0;
                tempm = 0.0;
                temel = 0.0;
                g = 1.0;
                c = -(C1 * aimz);
                d = C1 * rez;
                am = rv - 1.0;
                aak = 1.0;
                _k = 0;
                lab = 380;
            }
            380 => {
                let ajtemp = 2.0 * aak;
                temp4 = (1.0 - ajtemp) * ajtemp;
                ajp = rv - (4.0 * aak + 1.0);
                lab = 480;
            }
            390 => {
                aak += 1.0;
                _k += 1;
                let pr = rew;
                let pim = aimw;
                let amagn = powc(tempm, 2) + powc(temel, 2);
                rew = (tempc * tempm + tempd * temel) / amagn;
                aimw = (tempm * tempd - temel * tempc) / amagn;
                if (rew - pr).abs() >= EPS || (aimw - pim).abs() >= EPS {
                    lab = 380;
                } else {
                    return (rew, aimw);
                }
            }
            420 => {
                // wt is obtained from taylor series
                temp1 = r2 + ai2;
                let temp2 = 2.0 * temp1 * temp1;
                aj = -(r2 - ai2) / temp2;
                ak = 2.0 * rez * aimz / temp2;
                c = 0.0;
                b = 0.0;
                ajsig = 0.0;
                d = 0.0;
                _jsig = 0;
                g = 0.0;
                h = 0.0;
                el = 0.0;
                a = 1.0;
                am = 1.0;
                let sigp = C2;
                let expon = exp(temp2 * aj);
                expc = expon * cos(temp2 * ak);
                exps = -(expon * sin(temp2 * ak));
                sig2p = 2.0 * sigp;
                lab = 430;
            }
            430 => {
                let aj4sig = 4.0 * ajsig;
                let aj4sm1 = aj4sig - 1.0;
                let temp3 = 1.0 / (aj4sm1 * (aj4sig + 3.0));
                let tt4 = sig2p * (2.0 * ajsig - 1.0);
                temp4 = tt4 / (aj4sm1 * (aj4sig + 1.0) * (aj4sig - 3.0) * aj4sm1);
                ajp = aj + temp3;
                lab = 480;
            }
            440 => {
                ajsig += 1.0;
                _jsig += 1;
                let temp7 = rpi * (am * am + el * el);
                let reff = (aimz * (c * am + d * el) - rez * (am * d - c * el)) / temp7 / temp1;
                let aimf = (aimz * (am * d - c * el) + rez * (c * am + d * el)) / temp7 / temp1;
                let pr = rew;
                let pim = aimw;
                rew = expc - reff;
                aimw = exps - aimf;
                if (rew - pr).abs() >= EPS || (aimw - pim).abs() >= EPS {
                    // label 470
                    sig2p = 2.0 * ajsig;
                    lab = 430;
                } else {
                    return (rew, aimw);
                }
            }
            _ => {
                // label 480
                tempc = ajp * c + temp4 * a - ak * d;
                tempd = ajp * d + temp4 * b + ak * c;
                temel = ajp * el + temp4 * h + ak * am;
                tempm = ajp * am + temp4 * g - ak * el;
                a = c;
                b = d;
                g = am;
                h = el;
                c = tempc;
                d = tempd;
                am = tempm;
                el = temel;
                if (tempm.abs() + temel.abs() - UP) < ZERO {
                    if (tempm.abs() + temel.abs() - DN) <= ZERO {
                        c *= UP;
                        d *= UP;
                        am *= UP;
                        el *= UP;
                        tempc *= UP;
                        tempd *= UP;
                        tempm *= UP;
                        temel *= UP;
                    }
                } else {
                    c *= DN;
                    d *= DN;
                    am *= DN;
                    el *= DN;
                    tempc *= DN;
                    tempd *= DN;
                    tempm *= DN;
                    temel *= DN;
                }
                // label 520
                lab = if kw == 1 { 390 } else { 440 };
            }
        }
    }
}
