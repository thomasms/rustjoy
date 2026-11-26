//! PURR: probability tables and Bondarenko moments in the unresolved
//! resonance range, by Monte Carlo sampling of resonance ladders. Ported
//! from NJOY2016 `purr.f90` (module purm).
//!
//! The routines keep the Fortran structure: module variables are fields of
//! `Purr`, `save`d locals are fields too, and labels become loops or
//! `match` states where the control flow needs them. Arrays keep the
//! Fortran (1-based) indices.

#![allow(clippy::too_many_arguments)]
#![allow(clippy::collapsible_if)]
#![allow(clippy::collapsible_else_if)]
#![allow(clippy::needless_range_loop)]


use crate::endf::Ctx;
use crate::farr::{R3, RM, RV};
use crate::ffmt;
use crate::fio::{fatal_runtime, Kind, RecW, Val};
use crate::fmath::{atan, cos, exp, int, log, nint, pow, powc, sin, sqrt};
use crate::util::physics::{AMASSN, AMU, EV, HBAR, PI};
use crate::util::sigfig;

const MXNS0: usize = 100;
const MEUNR: i32 = 150;
const JX: i32 = 10000;
const MAXSCR: i32 = 20000;
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

/// Format `(head,1p,e11.3,10e11.3/(16x,10e11.3))` for a value and a list.
fn list10(head: &str, t: f64, v: &[f64]) -> String {
    let mut s = format!("{}{}", head, ffmt::e1p(t, 11, 3));
    let n1 = v.len().min(10);
    for x in &v[..n1] {
        s.push_str(&ffmt::e1p(*x, 11, 3));
    }
    if v.len() == 10 {
        // the slash is done before the list is found to be empty
        s.push('\n');
        return s;
    }
    let mut k = n1;
    while k < v.len() {
        s.push('\n');
        s.push_str(&" ".repeat(16));
        let m = (v.len() - k).min(10);
        for x in &v[k..k + m] {
            s.push_str(&ffmt::e1p(*x, 11, 3));
        }
        k += m;
    }
    s
}

/// PURR module variables (and the `save`d locals of unrest and rann).
pub struct Purr {
    nendf: i32,
    nin: i32,
    nout: i32,
    init: i32,
    nunr: i32,
    lssf: i32,
    iinel: i32,
    iabso: i32,
    iprint: i32,
    nermax: i32,
    nladr: i32,
    nmode: i32,
    nbin: i32,
    /// tr(41,27), ti, trs, tis, column-major.
    tr: Vec<f64>,
    ti: Vec<f64>,
    trs: Vec<f64>,
    tis: Vec<f64>,
    nsect: i32,
    isot: [i32; NSECT + 1],
    modet: [i32; NSECT + 1],
    ibaset: [i32; NSECT + 1],
    abnt: [f64; NSECT + 1],
    elt: [f64; NSECT + 1],
    eht: [f64; NSECT + 1],
    nseq0: i32,
    nro: i32,
    naps: i32,
    intunr: i32,
    e: f64,
    cth: [f64; MXNS0 + 1],
    csz: [f64; MXNS0 + 1],
    cc2p: [f64; MXNS0 + 1],
    cs2p: [f64; MXNS0 + 1],
    cgn: [f64; MXNS0 + 1],
    cgg: [f64; MXNS0 + 1],
    cgf: [f64; MXNS0 + 1],
    cgx: [f64; MXNS0 + 1],
    cgt: [f64; MXNS0 + 1],
    dbar: [f64; MXNS0 + 1],
    spot: f64,
    dbarin: f64,
    sigi: [f64; 5],
    ndfn: [i32; MXNS0 + 1],
    ndff: [i32; MXNS0 + 1],
    ndfx: [i32; MXNS0 + 1],
    eunr: RV,
    sb: RV,
    bval: R3,
    tmin: Vec<f64>,
    tmax: Vec<f64>,
    tsum: Vec<f64>,
    arry: RV,
    arry_alloc: bool,
    ipl: i32,
    epl: [f64; 201],
    sigpl: R3,
    kk: i32,
    rr: f64,
    mflg1: i32,
    mflg2: i32,
    // saved locals of unrest
    navoid: i32,
    elow: f64,
    // saved locals of rann
    iff: i32,
    iy: i32,
    ir: [i32; 98],
}

impl Default for Purr {
    fn default() -> Self {
        Self::new()
    }
}

#[inline]
fn t41(v: &[f64], i: i32, j: i32) -> f64 {
    v[((i - 1) + (j - 1) * 41) as usize]
}

impl Purr {
    pub fn new() -> Self {
        Purr {
            nendf: 0,
            nin: 0,
            nout: 0,
            init: 0,
            nunr: 0,
            lssf: 0,
            iinel: 0,
            iabso: 0,
            iprint: 0,
            nermax: 0,
            nladr: 0,
            nmode: 0,
            nbin: 0,
            tr: vec![0.0; 41 * 27],
            ti: vec![0.0; 41 * 27],
            trs: vec![0.0; 41 * 27],
            tis: vec![0.0; 41 * 27],
            nsect: 0,
            isot: [0; NSECT + 1],
            modet: [0; NSECT + 1],
            ibaset: [0; NSECT + 1],
            abnt: [0.0; NSECT + 1],
            elt: [0.0; NSECT + 1],
            eht: [0.0; NSECT + 1],
            nseq0: 0,
            nro: 0,
            naps: 0,
            intunr: 0,
            e: 0.0,
            cth: [0.0; MXNS0 + 1],
            csz: [0.0; MXNS0 + 1],
            cc2p: [0.0; MXNS0 + 1],
            cs2p: [0.0; MXNS0 + 1],
            cgn: [0.0; MXNS0 + 1],
            cgg: [0.0; MXNS0 + 1],
            cgf: [0.0; MXNS0 + 1],
            cgx: [0.0; MXNS0 + 1],
            cgt: [0.0; MXNS0 + 1],
            dbar: [0.0; MXNS0 + 1],
            spot: 0.0,
            dbarin: 0.0,
            sigi: [0.0; 5],
            ndfn: [0; MXNS0 + 1],
            ndff: [0; MXNS0 + 1],
            ndfx: [0; MXNS0 + 1],
            eunr: RV::new(MEUNR as usize + SLACK),
            sb: RV::new(0),
            bval: R3::empty(),
            tmin: Vec::new(),
            tmax: Vec::new(),
            tsum: Vec::new(),
            arry: RV::new(0),
            arry_alloc: false,
            ipl: 0,
            epl: [0.0; 201],
            sigpl: R3::empty(),
            kk: 0,
            rr: 0.0,
            mflg1: 0,
            mflg2: 0,
            navoid: 0,
            elow: 0.0,
            iff: 0,
            iy: 0,
            ir: [0; 98],
        }
    }

    // ================================================================ purr
    /// Probabilistic unresolved calculation of Bondarenko moments and
    /// probability tables.
    pub fn purr(&mut self, c: &mut Ctx) {
        const TREF: f64 = 300.0;
        const SIGMIN: f64 = 1.0e-5;
        const ZERO: f64 = 0.0;
        let (mut nb, mut nw) = (0i32, 0i32);
        let mut bkgz = [0.0f64; 5];
        let mut temp: Vec<f64> = vec![0.0];
        let mut sigz: Vec<f64> = vec![0.0];
        let mut sigu = R3::empty();
        let mut tabl = R3::empty();
        let mut tval = RM::new(0, 0);
        let mut heat = R3::empty();
        let mut alloc_temp = false;

        //--initialize
        self.ipl = 0;
        self.nmode = 0;
        self.nermax = 1000;
        let nsamp = 10000i32;
        c.nsh = 0;
        let time = c.timer();
        c.io.so(&format!(
            "\n purr...probabalistic unresolved calculation{}{}s",
            " ".repeat(25),
            ffmt::f(time, 8, 1)
        ));
        c.io.se(&format!("\n purr...{}{}s", " ".repeat(61), ffmt::f(time, 8, 1)));
        let v = c.io.input.read(&[Kind::Int, Kind::Int, Kind::Int]);
        self.nendf = geti(&v[0], self.nendf);
        self.nin = geti(&v[1], self.nin);
        self.nout = geti(&v[2], self.nout);
        let (nendf, nin, nout) = (self.nendf, self.nin, self.nout);
        if (nin < 0 && nout > 0) || (nin > 0 && nout < 0) {
            c.error("purr", "mode conversion between nin and nout not allowed", " ");
        }
        c.openz(nendf, 0);
        c.openz(nin, 0);
        c.openz(nout, 1);
        let nscr = 10i32;
        c.openz(-nscr, 1);
        let mut a = RV::new(MAXSCR as usize + SLACK);
        c.io.so(&format!(
            " unit for input endf tape ............. {}\n unit for input pendf tape ............ {}\n unit for output pendf tape ........... {}",
            ffmt::i(nendf, 10),
            ffmt::i(nin, 10),
            ffmt::i(nout, 10)
        ));
        self.uwtab2();
        c.repoz(nin);
        c.repoz(nendf);
        c.repoz(nout);
        c.tpidio(nin, nout, 0, &mut a.all(), &mut nb, &mut nw);
        c.tpidio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);

        // set random number sequence
        self.kk = -101;
        let _ = self.rann();

        //--loop over requested materials
        'mat: loop {
            self.iprint = 1;
            let mut nunx = 0i32;
            let mut ntemp = 1i32;
            let mut nsigz = 1i32;
            self.init = 0;
            let v = c.io.input.read(&[Kind::Int; 7]);
            let matd = geti(&v[0], 0);
            ntemp = geti(&v[1], ntemp);
            nsigz = geti(&v[2], nsigz);
            self.nbin = geti(&v[3], self.nbin);
            self.nladr = geti(&v[4], self.nladr);
            self.iprint = geti(&v[5], self.iprint);
            nunx = geti(&v[6], nunx);
            if matd == 0 {
                break 'mat;
            }
            let nbin = self.nbin;
            c.repoz(-nscr);
            let nt = ntemp.max(0) as usize;
            let ns = nsigz.max(0) as usize;
            temp = vec![0.0; nt + 1];
            sigz = vec![0.0; ns + 1];
            sigu = R3::new(5, ns, nt);
            self.sigpl = R3::new(200, 5, ns);
            alloc_temp = true;
            let kinds: Vec<Kind> = (0..nt).map(|_| Kind::Real).collect();
            let v = c.io.input.read(&kinds);
            for i in 0..nt {
                temp[i + 1] = getr(&v[i], 0.0);
            }
            let kinds: Vec<Kind> = (0..ns).map(|_| Kind::Real).collect();
            let v = c.io.input.read(&kinds);
            for i in 0..ns {
                sigz[i + 1] = getr(&v[i], 0.0);
            }
            c.io.so(&format!(" temperatures ......................... {}", ffmt::e1p(temp[1], 10, 3)));
            for i in 2..=nt {
                c.io.so(&format!("{}{}", " ".repeat(40), ffmt::e1p(temp[i], 10, 3)));
            }
            c.io.so(" sigma zero values ....................   infinity");
            for i in 2..=ns {
                c.io.so(&format!("{}{}", " ".repeat(40), ffmt::e1p(sigz[i], 10, 3)));
            }
            c.io.so(&format!(
                " number of probability bins ........... {}\n number of resonance ladders .......... {}\n print option (0=min, 1=max) .......... {}\n no. of energy points (0=all) ......... {}",
                ffmt::i(nbin, 10),
                ffmt::i(self.nladr, 10),
                ffmt::i(self.iprint, 10),
                ffmt::i(nunx, 10)
            ));
            if nbin < 15 {
                c.error("purr", "nbin should be 15 or more", " ");
            }

            //--allocate storage for ladders and tables
            if self.arry_alloc {
                fatal_runtime("Fortran runtime error: Attempting to allocate already allocated variable 'arry'");
            }
            self.arry = RV::new(JX as usize + SLACK);
            self.arry_alloc = true;
            let nb_ = nbin.max(0) as usize;
            tabl = R3::new(nb_, 5, nt);
            tval = RM::new(nb_, nt);

            //--process this material
            let time = c.timer();
            let s = format!("\n mat = {}{}{}s", ffmt::i(matd, 4), " ".repeat(58), ffmt::f(time, 8, 1));
            c.io.so(&s);
            c.io.se(&s);

            //--read resonance parameters from file 2 on endf tape.
            //--read background cross sections from file 3
            self.nunr = 0;
            c.findf(matd, 1, 451, nendf);
            c.contio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
            let za = c.c1h;
            let awr = c.c2h;
            let lrp = c.l1h;
            // 0: normal, 1: no resonance parameters, 2: no unresolved
            let mut skip = 0;
            if lrp <= 0 {
                skip = 1;
            } else {
                c.contio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
                if c.n1h != 0 {
                    c.iverf = 4;
                } else if c.n2h == 0 {
                    c.iverf = 5;
                } else {
                    c.iverf = 6;
                }
                c.findf(matd, 2, 151, nendf);
                c.contio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
                self.rdf2un(c, &mut a);
                if self.nunr == 0 {
                    skip = 2;
                }
            }
            if skip != 0 {
                let strng1 = if skip == 1 {
                    format!("mat{} has no resonance parameters", ffmt::i(matd, 5))
                } else {
                    format!("mat{} has no unresolved parameters", ffmt::i(matd, 5))
                };
                c.mess("purr", &strng1, "copy as is to nout");
                // label 385
                c.findf(matd, 1, 451, nin);
                let mut it = 0;
                loop {
                    c.contio(nin, 0, 0, &mut a.all(), &mut nb, &mut nw);
                    if c.math != matd {
                        continue 'mat;
                    }
                    c.contio(0, nout, 0, &mut a.all(), &mut nb, &mut nw);
                    c.tomend(nin, nout, 0, &mut a.all());
                    it += 1;
                    if it >= ntemp {
                        continue 'mat;
                    }
                }
            }
            let nunr = self.nunr;
            let maxtst = (12 + nsigz + nunr * (1 + 5 * nsigz)).max(12 + (1 + 6 * nbin) * nunr);
            if maxtst > MAXSCR {
                let strng1 = format!("maxscr is too small, increase to at least {} words", ffmt::i(maxtst, 7));
                c.error("purr", &strng1, " ");
            }
            c.findf(matd, 3, 1, nendf);
            self.rdf3un(c, &mut a);

            //--read in the total and partial heating cross sections
            heat = R3::new(4, nunr as usize, nt);
            let mut ihave = 0;
            self.rdheat(c, &mut a, &mut heat, &temp, ntemp, nunr, &mut ihave, matd);
            if ihave == 0 {
                c.mess("purr", "no heating found on pendf", "ur heating set to zero");
            }
            if ihave == 1 {
                c.mess("purr", "no partial heating xsecs found on pendf", "ur heating will not selfshield");
            }

            //--compute unresolved resonance cross-sections
            //--for all grid points, temperatures, and sigma zero values.
            if nunx == 0 || nunx > nunr {
                nunx = nunr;
            }
            if nunx == 1 {
                nunx = 2;
            }
            let mut nstep = nunr / (nunx - 1);
            if nstep < 1 {
                nstep = 1;
            }
            let mut ie = 1 - nstep;
            self.mflg1 = 0;
            self.mflg2 = 0;
            let mut w = Work::new(self.nermax, nsamp, ntemp);
            while ie < nunr {
                ie += nstep;
                if ie > nunr {
                    ie = nunr;
                }
                let ez = self.eunr[ie];
                bkgz[1] = self.sb[ie];
                bkgz[2] = self.sb[nunr + ie];
                bkgz[3] = self.sb[nunr * 2 + ie];
                bkgz[4] = self.sb[nunr * 3 + ie];
                let mut sigx = bkgz[1] - bkgz[2] - bkgz[3] - bkgz[4];
                if sigx < SIGMIN {
                    sigx = 0.0;
                }
                let _ = sigx;
                self.unresx(c, ez, TREF);
                self.unrest(c, &bkgz, &sigz, &mut sigu, &mut tabl, &mut tval, &mut w, &temp, nsigz, ntemp, nbin, nsamp);

                //--purr can run for a long time ... keep the user informed
                let time = c.timer();
                c.io.se(&format!(
                    "  {} of {} loops done for all temps & sig0s.{}{}s",
                    ffmt::i(ie, 5),
                    ffmt::i(nunx, 5),
                    " ".repeat(19),
                    ffmt::f(time, 8, 1)
                ));

                //--write bondarenko data and probability table to scratch file.
                let mut r = RecW::new();
                r.r(ez);
                c.io.write_rec(nscr, r);
                let mut r = RecW::new();
                r.rs(&tabl.v);
                c.io.write_rec(nscr, r);
                let mut r = RecW::new();
                r.rs(&sigu.v);
                c.io.write_rec(nscr, r);
            }

            //--check for previous mt152 and mt153
            let mut new152 = 1;
            let mut new153 = 1;
            c.findf(matd, 2, 151, nin);
            c.tosend(nin, 0, 0, &mut a.all());
            c.contio(nin, 0, 0, &mut a.all(), &mut nb, &mut nw);
            if c.mth == 152 {
                new152 = 0;
                c.tosend(nin, 0, 0, &mut a.all());
                c.contio(nin, 0, 0, &mut a.all(), &mut nb, &mut nw);
            }
            if c.mth == 153 {
                new153 = 0;
            }
            c.findf(matd, 1, 451, nin);

            //--write data to new pendf tape for all temperatures
            for it in 1..=ntemp {
                let temz = temp[it as usize];
                let mut nx = 0;
                c.contio(nin, 0, 0, &mut a.all(), &mut nb, &mut nw);
                if c.iverf < 5 {
                    nx = c.n2h;
                    if nx > 0 && new152 > 0 {
                        a[6] += 1.0;
                    }
                    if nx > 0 && new153 > 0 {
                        a[6] += 1.0;
                    }
                }
                c.contio(0, nout, 0, &mut a.all(), &mut nb, &mut nw);
                if c.iverf >= 5 {
                    c.contio(nin, nout, 0, &mut a.all(), &mut nb, &mut nw);
                }
                if c.iverf >= 6 {
                    c.contio(nin, nout, 0, &mut a.all(), &mut nb, &mut nw);
                }
                c.hdatio(nin, 0, 0, &mut a.all(), &mut nb, &mut nw);
                a[1] = temz;
                if c.iverf >= 5 {
                    nx = c.n2h;
                    if nx > 0 && new152 > 0 {
                        a[6] += 1.0;
                    }
                    if nx > 0 && new153 > 0 {
                        a[6] += 1.0;
                    }
                }
                c.hdatio(0, nout, 0, &mut a.all(), &mut nb, &mut nw);
                while nb != 0 {
                    c.moreio(nin, nout, 0, &mut a.all(), &mut nb, &mut nw);
                }
                nw = nx;
                if nw != 0 {
                    let nc = 6 * (nx + new152 + new153);
                    let mut cc = RV::new(nc.max(1) as usize + 64);
                    let nd = 6 * nx;
                    let mut d = RV::new(nd.max(1) as usize + 64);
                    c.dictio(nin, 0, 0, &mut d.all(), &mut nb, &mut nw);
                    let mut j = 0i32;
                    let mut newmat = 0;
                    let mut ncds = nsigz + nunx * (1 + 5 * nsigz);
                    ncds = 2 + (ncds - 1) / 6;
                    let mut nc153 = (1 + 6 * nbin) * nunx;
                    nc153 = 2 + (nc153 - 1) / 6;
                    let mut i = 1;
                    while i <= nw {
                        let mut copy = true;
                        if newmat <= 0 {
                            let mfd = nint(d[i + 2]);
                            let mtd = nint(d[i + 3]);
                            if mfd < 2 || (mfd == 2 && mtd == 151) {
                                // label 220
                            } else if mfd == 2 && (mtd == 152 || mtd == 153) {
                                copy = false;
                            } else {
                                newmat = 1;
                                cc[1 + j] = 0.0;
                                cc[2 + j] = 0.0;
                                cc[3 + j] = 2.0;
                                cc[4 + j] = 152.0;
                                cc[5 + j] = ncds as f64;
                                cc[6 + j] = 0.0;
                                cc[7 + j] = 0.0;
                                cc[8 + j] = 0.0;
                                cc[9 + j] = 2.0;
                                cc[10 + j] = 153.0;
                                cc[11 + j] = nc153 as f64;
                                cc[12 + j] = 0.0;
                                j += 12;
                            }
                        }
                        if copy {
                            for k in 1..=6 {
                                cc[j + k] = d[k + i - 1];
                            }
                            j += 6;
                        }
                        i += 6;
                    }
                    let mut nwd = j / 6;
                    c.dictio(0, nout, 0, &mut cc.all(), &mut nb, &mut nwd);
                }
                // copy through file 2, mt 151
                // label 250
                c.tofend(nin, nout, 0, &mut a.all());
                c.tosend(nin, nout, 0, &mut a.all());
                c.contio(nin, 0, 0, &mut a.all(), &mut nb, &mut nw);
                if c.mth == 152 {
                    c.tosend(nin, 0, 0, &mut a.all());
                    c.contio(nin, 0, 0, &mut a.all(), &mut nb, &mut nw);
                }
                if c.mth == 153 {
                    c.tosend(nin, 0, 0, &mut a.all());
                    c.contio(nin, 0, 0, &mut a.all(), &mut nb, &mut nw);
                }

                //--write new mt152 on output tape
                //--for bondarenko self-shielded cross sections
                c.math = matd;
                c.mfh = 2;
                c.mth = 152;
                let mut l = 1i32;
                a[l] = za;
                a[l + 1] = awr;
                a[l + 2] = self.lssf as f64;
                a[l + 3] = 0.0;
                a[l + 4] = 0.0;
                a[l + 5] = self.intunr as f64;
                a[l + 6] = temz;
                a[l + 7] = 0.0;
                a[l + 8] = 5.0;
                a[l + 9] = nsigz as f64;
                a[l + 10] = (nsigz + nunx * (1 + 5 * nsigz)) as f64;
                a[l + 11] = nunx as f64;
                l += 11;
                for i in 1..=nsigz {
                    l += 1;
                    a[l] = sigz[i as usize];
                }
                c.repoz(-nscr);
                let mut ie = 1 - nstep;
                while ie < nunr {
                    ie += nstep;
                    if ie > nunr {
                        ie = nunr;
                    }
                    self.read_scr(c, nscr, &mut tabl, &mut sigu);
                    l += 1;
                    a[l] = self.eunr[ie];
                    for ix in 1..=5 {
                        for is in 1..=nsigz {
                            l += 1;
                            a[l] = sigfig(sigu[(ix, is, it)], 7, 0);
                        }
                    }
                }
                c.contio(0, nout, 0, &mut a.all(), &mut nb, &mut nw);
                let mut lstart = 7;
                let mut nwds = l - 1 - 5;
                nw = nwds;
                if nw > c.npage + 6 {
                    nw = c.npage + 6;
                }
                c.listio(0, nout, 0, &mut a.at(lstart), &mut nb, &mut nw);
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
                    c.moreio(0, nout, 0, &mut a.at(lstart), &mut nb, &mut nw);
                }
                c.asend(nout, 0);

                //--write new mt153 on output tape
                //--for probability tables
                c.math = matd;
                c.mfh = 2;
                c.mth = 153;
                let mut n = 1i32;
                a[n] = za;
                a[n + 1] = awr;
                a[n + 2] = self.iinel as f64;
                a[n + 3] = self.iabso as f64;
                a[n + 4] = self.intunr as f64;
                a[n + 5] = nbin as f64;
                a[n + 6] = temz;
                a[n + 7] = 0.0;
                a[n + 8] = self.lssf as f64;
                a[n + 9] = 0.0;
                a[n + 10] = ((1 + 6 * nbin) * nunx) as f64;
                a[n + 11] = nunx as f64;
                n += 11;
                c.repoz(-nscr);
                let lssf = self.lssf;
                let mut ie = 1 - nstep;
                while ie < nunr {
                    ie += nstep;
                    if ie > nunr {
                        ie = nunr;
                    }
                    self.read_scr(c, nscr, &mut tabl, &mut sigu);
                    n += 1;
                    a[n] = sigfig(self.eunr[ie], 7, 0);
                    let n1 = n;
                    for i in 1..=5 {
                        for j in 1..=nbin {
                            n += 1;
                            a[n] = sigfig(tabl[(j, i, it)], 7, 0);
                        }
                    }
                    for _j in 1..=nbin {
                        n += 1;
                        a[n] = 0.0;
                    }
                    for i in 2..=5 {
                        for j in 1..=nbin {
                            let l = n1 + j + nbin * (i - 1);
                            if lssf == 1 {
                                if sigu[(i - 1, 1, 1)] != 0.0 {
                                    a[l] /= sigu[(i - 1, 1, 1)];
                                } else {
                                    a[l] = 1.0;
                                }
                            }
                            a[l] = sigfig(a[l], 7, 0);
                        }
                    }
                    for j in 1..=nbin {
                        let l = n1 + j + 5 * nbin;
                        if ihave == 1 {
                            if lssf == 1 {
                                a[l] = 1.0;
                            } else {
                                a[l] = heat[(1, ie, it)] / sigu[(1, 1, 1)];
                            }
                        } else if ihave == 2 {
                            a[l] += heat[(1, ie, it)] - heat[(2, ie, it)] - heat[(3, ie, it)] - heat[(4, ie, it)];
                            let mut h = heat[(2, ie, it)];
                            if lssf == 1 {
                                h *= a[n1 + j + 2 * nbin];
                            } else if sigu[(2, 1, 1)] != ZERO {
                                h = h * a[n1 + j + 2 * nbin] / sigu[(2, 1, 1)];
                            }
                            a[l] += h;
                            let mut h = heat[(3, ie, it)];
                            if lssf == 1 {
                                h *= a[n1 + j + 3 * nbin];
                            } else if sigu[(3, 1, 1)] != ZERO {
                                h = h * a[n1 + j + 3 * nbin] / sigu[(3, 1, 1)];
                            }
                            a[l] += h;
                            let mut h = heat[(4, ie, it)];
                            if lssf == 1 {
                                h *= a[n1 + j + 4 * nbin];
                            } else if sigu[(4, 1, 1)] != ZERO {
                                h = h * a[n1 + j + 4 * nbin] / sigu[(4, 1, 1)];
                            }
                            a[l] += h;
                            if a[n1 + j + nbin] != ZERO {
                                a[l] /= a[n1 + j + nbin];
                            }
                            if lssf == 1 {
                                a[l] /= heat[(1, ie, it)];
                            }
                        }
                    }
                }
                c.contio(0, nout, 0, &mut a.all(), &mut nb, &mut nw);
                let mut lstart = 7;
                let mut nwds = 6 * nbin * nunx;
                nw = nwds;
                if nw > c.npage + 6 {
                    nw = c.npage + 6;
                }
                c.listio(0, nout, 0, &mut a.at(lstart), &mut nb, &mut nw);
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
                    c.moreio(0, nout, 0, &mut a.at(lstart), &mut nb, &mut nw);
                }
                c.asend(nout, 0);
                // write file end card for file 2
                c.afend(nout, 0);
                // copy rest of material to output tape
                c.contio(nin, nout, 0, &mut a.all(), &mut nb, &mut nw);
                c.tomend(nin, nout, 0, &mut a.all());
            }

            //--write report of calculation
            let time = c.timer();
            c.io.so(&format!(
                "\n generated cross sections at {} points{}{}s",
                ffmt::i(nunx, 3),
                " ".repeat(30),
                ffmt::f(time, 8, 1)
            ));
            self.sb = RV::new(0);
        }

        //--finished with purr
        // label 400
        if !self.arry_alloc {
            fatal_runtime("Fortran runtime error: Attempt to DEALLOCATE unallocated 'arry'");
        }
        self.arry_alloc = false;
        if !alloc_temp {
            fatal_runtime("Fortran runtime error: Attempt to DEALLOCATE unallocated 'temp'");
        }
        let _ = (&temp, &sigz, &sigu, &tabl, &tval, &heat);
        c.atend(nout, 0);
        c.repoz(nout);
        c.closz(nendf);
        c.closz(nin);
        c.closz(nout);
        c.closz(-nscr);
        let time = c.timer();
        c.io.so(&format!("\n{}{}s\n{}*******", " ".repeat(69), ffmt::f(time, 8, 1), "*".repeat(70)));
    }

    /// Read one energy point (ez, tabl, sigu) back from the scratch file.
    fn read_scr(&mut self, c: &mut Ctx, nscr: i32, tabl: &mut R3, sigu: &mut R3) {
        c.io.read_rec_with(nscr, |r| {
            let _ = r.r();
        });
        c.io.read_rec_with(nscr, |r| {
            for x in tabl.v.iter_mut() {
                *x = r.r();
            }
        });
        c.io.read_rec_with(nscr, |r| {
            for x in sigu.v.iter_mut() {
                *x = r.r();
            }
        });
    }

    // ============================================================== rdf2un
    /// Read and/or copy unresolved resonance parameters from File 2.
    fn rdf2un(&mut self, c: &mut Ctx, a: &mut RV) {
        const EGRIDU: [f64; 78] = [
            1.0e1, 1.25e1, 1.5e1, 1.7e1, 2.0e1, 2.5e1, 3.0e1, 3.5e1, 4.0e1, 5.0e1, 6.0e1, 7.2e1, 8.5e1, 1.0e2,
            1.25e2, 1.5e2, 1.7e2, 2.0e2, 2.5e2, 3.0e2, 3.5e2, 4.0e2, 5.0e2, 6.0e2, 7.2e2, 8.5e2, 1.0e3, 1.25e3,
            1.5e3, 1.7e3, 2.0e3, 2.5e3, 3.0e3, 3.5e3, 4.0e3, 5.0e3, 6.0e3, 7.2e3, 8.5e3, 1.0e4, 1.25e4, 1.5e4,
            1.7e4, 2.0e4, 2.5e4, 3.0e4, 3.5e4, 4.0e4, 5.0e4, 6.0e4, 7.2e4, 8.5e4, 1.0e5, 1.25e5, 1.5e5, 1.7e5,
            2.0e5, 2.5e5, 3.0e5, 3.5e5, 4.0e5, 5.0e5, 6.0e5, 7.2e5, 8.5e5, 1.0e6, 1.25e6, 1.5e6, 1.7e6, 2.0e6,
            2.5e6, 3.0e6, 3.5e6, 4.0e6, 5.0e6, 6.0e6, 7.2e6, 8.5e6,
        ];
        const NGRIDU: usize = 78;
        const ETOP: f64 = 5.0e6;
        const SMALL: f64 = 1.0e-12;
        const WIDE: f64 = 1.26e0;
        let (mut nb, mut nw) = (0i32, 0i32);
        let nendf = self.nendf;

        let mut indep = 0;
        let mut elr = 0.0f64;
        let mut ehr = ETOP;
        self.nunr = 1;
        self.eunr[1] = ETOP;

        //--head card read in calling program
        let nis = c.n1h;
        self.nsect = 0;
        let mut inow = 1i32;

        //--loop over isotopes
        'iso: for is in 1..=nis {
            c.contio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
            let abn = c.c2h;
            let lfw = c.l2h;
            let ner = c.n1h;
            let mut ier = 0;

            //--loop over all energy ranges
            loop {
                // label 115
                c.contio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
                if c.mth == 0 {
                    break 'iso;
                }
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
                    c.contio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
                    let nls = c.n1h;
                    if lrf == 4 || lrf == 7 {
                        c.listio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
                        while nb != 0 {
                            c.moreio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
                        }
                    }
                    for _l in 1..=nls {
                        if lrf != 4 {
                            c.listio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
                            while nb != 0 {
                                c.moreio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
                            }
                            if lrf == 7 {
                                c.listio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
                                while nb != 0 {
                                    c.moreio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
                                }
                            }
                        } else {
                            c.contio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
                            let njs = c.n1h;
                            for _j in 1..=njs {
                                c.listio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
                                while nb != 0 {
                                    c.moreio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
                                }
                            }
                        }
                    }
                    if ier == ner {
                        continue 'iso;
                    }
                    continue;
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
                ilist2(et, &mut self.eunr, &mut self.nunr);
                let et = sigfig(el, 7, 1);
                ilist2(et, &mut self.eunr, &mut self.nunr);
                let et = sigfig(eh, 7, -1);
                ilist2(et, &mut self.eunr, &mut self.nunr);
                let et = sigfig(eh, 7, 1);
                ilist2(et, &mut self.eunr, &mut self.nunr);

                // if present, read and store the energy-dependent
                // scattering radius tab1 data.
                if self.nro == 1 {
                    c.tab1io(nendf, 0, 0, &mut self.arry.at(inow), &mut nb, &mut nw);
                    let mut jj = inow + nw;
                    while nb != 0 {
                        c.moreio(nendf, 0, 0, &mut self.arry.at(jj), &mut nb, &mut nw);
                        jj += nw;
                        if jj > JX {
                            c.error("rdf2un", "storage in a exceeded", " ");
                        }
                    }
                    inow = jj;
                }

                //--branch to specified representation
                let arry = &mut self.arry;
                if lrf != 2 {
                    if lfw != 1 {
                        //--all parameters independent of energy
                        c.contio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
                        let spi = c.c1h;
                        let ay = c.c2h;
                        self.lssf = c.l1h;
                        let nls = c.n1h;
                        // loop over l states
                        for l in 1..=nls {
                            c.listio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
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
                            let mut jnow = 1 + 5;
                            for _n in 1..=njs {
                                for i in 1..=5 {
                                    arry[inow] = a[jnow + i];
                                    inow += 1;
                                }
                                jnow += 6;
                            }
                        }
                        if inow > JX {
                            c.error("rdf2un", "storage exceeded", " ");
                        }
                        indep = 1;
                        self.intunr = 5;
                    } else {
                        //--fission widths energy dependent
                        c.listio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
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
                            arry[inow1 + i] = a[1 + i + 5];
                            let enow = a[1 + i + 5];
                            arry[inow1 + i] = enow;
                            // use as energy nodes
                            if i > 1 && i < ne {
                                if self.nunr >= MEUNR {
                                    c.error("rdunf2", "too many ur energy points", " ");
                                }
                                if enow >= el && enow <= eh {
                                    ilist2(enow, &mut self.eunr, &mut self.nunr);
                                }
                            }
                        }
                        // loop over l states
                        for l in 1..=nls {
                            c.contio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
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
                                c.listio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
                                ne = c.n1h - 6;
                                arry[inow] = c.l2h as f64;
                                // store parameters
                                for ir in 1..=5 {
                                    arry[inow + ir] = a[1 + 5 + ir];
                                }
                                let inow1 = inow + 5;
                                // store fission widths
                                for k in 1..=ne {
                                    arry[inow1 + k] = a[1 + 11 + k];
                                }
                                inow = inow + 6 + ne;
                            }
                        }
                        if inow > JX {
                            c.error("rdf2un", "storage exceeded", " ");
                        }
                        self.intunr = 5;
                    }
                } else {
                    //--all parameters energy dependent
                    let ist = inow;
                    // first cont record
                    c.contio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
                    self.lssf = c.l1h;
                    let nls = c.n1h;
                    arry[inow + 1] = c.c2h;
                    arry[inow + 2] = c.c1h;
                    arry[inow + 3] = c.n1h as f64;
                    inow += 4;
                    // loop over l states
                    for l in 1..=nls {
                        c.contio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
                        if l == 1 {
                            arry[ist] = c.c1h;
                        }
                        let njs = c.n1h;
                        arry[inow] = njs as f64;
                        arry[inow + 1] = c.l1h as f64;
                        inow += 2;
                        // loop over j states
                        for j in 1..=njs {
                            c.listio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
                            let mut loc = 1 + nw;
                            while nb != 0 {
                                c.moreio(nendf, 0, 0, &mut a.at(loc), &mut nb, &mut nw);
                                loc += nw;
                            }
                            let ne = c.n2h;
                            self.intunr = c.l1h;
                            arry[inow] = c.c1h;
                            arry[inow + 1] = c.l1h as f64;
                            arry[inow + 2] = ne as f64;
                            // store numbers of degrees of freedom
                            for k in 3..=6 {
                                arry[inow + k] = a[1 + 5 + k];
                            }
                            inow += 7;
                            let mut jnow = 1 + 11;
                            // store parameters
                            for n in 1..=ne {
                                for k in 1..=6 {
                                    if a[jnow + k] <= 0.0 {
                                        a[jnow + k] = SMALL;
                                    }
                                    arry[inow + k - 1] = a[jnow + k];
                                }
                                inow += 6;
                                if n != 1 && n != ne && l == 1 && j == 1 {
                                    // store energy as a node
                                    let enow = a[jnow + 1];
                                    if self.nunr >= MEUNR {
                                        c.error("rdunf2", "too many ur energy points", " ");
                                    }
                                    if enow >= el && enow <= eh {
                                        ilist2(enow, &mut self.eunr, &mut self.nunr);
                                    }
                                }
                                jnow += 6;
                            }
                        }
                    }
                    if inow > JX {
                        c.error("rdf2un", "storage exceeded", " ");
                    }
                }
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

        //--remove first and last energy nodes
        //--flag resolved-unresolved overlap energies
        self.nunr += 1;
        let mut i = 1i32;
        let mut elast = self.eunr[2];
        loop {
            i += 1;
            let enext = self.eunr[1 + i];
            if enext >= ETOP {
                break;
            }
            if enext >= WIDE * elast || indep != 0 {
                let mut et = elast;
                loop {
                    let mut found = false;
                    let mut enut = 0.0f64;
                    for ii in 0..NGRIDU {
                        enut = EGRIDU[ii];
                        if enut > et + et / 100.0 {
                            found = true;
                            break;
                        }
                    }
                    if !found {
                        enut = enext;
                    }
                    et = enut;
                    if et >= enext {
                        break;
                    }
                    if self.nunr >= MEUNR {
                        c.error("rdunf2", "too many ur energy points", " ");
                    }
                    ilist2(et, &mut self.eunr, &mut self.nunr);
                    i += 1;
                }
            }
            elast = self.eunr[1 + i];
        }
        let mut iovlp = 0;
        self.nunr -= 2;
        for i in 2..=self.nunr {
            let mut et = self.eunr[i];
            if et < elr {
                et = -et;
            }
            if et > ehr {
                et = -et;
            }
            if et < 0.0 {
                iovlp = 1;
            }
            self.eunr[i - 1] = et;
        }
        self.nunr -= 1;
        if iovlp == 1 {
            c.mess("rdf2un", "resolved-unresolved overlap energies", "marked with minus signs");
        }
    }

    // ============================================================== rdf3un
    /// Read unresolved region background cross sections from File 3.
    /// Compute the background cross sections on the unresolved energy grid
    /// and store for use in unresx. Also, analyze the competitive
    /// reactions, if any.
    fn rdf3un(&mut self, c: &mut Ctx, a: &mut RV) {
        const MTX: [i32; 5] = [0, 1, 2, 18, 102];
        const UP: f64 = 1.00001e0;
        const DN: f64 = 0.99999e0;
        const SMALL: f64 = 1.0e-5;
        const TOL: f64 = 1.0e-6;
        let (mut nb, mut nw) = (0i32, 0i32);
        let (mut enext, mut idis) = (0.0f64, 0i32);
        let nendf = self.nendf;
        let nunr = self.nunr;

        //--allocate storage
        let mut thr: Vec<f64> = vec![0.0; 2 * 1000 + 1];
        self.sb = RV::new((nunr * 4).max(0) as usize);

        //--loop over reactions
        let mut jthr = 0i32;
        let mut ix = 1usize;
        let mtmax = 891;
        c.contio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
        while c.mfh == 3 && c.mth <= mtmax {
            let mth = c.mth;
            if mth < 6
                || (mth > 10 && mth < 12)
                || (mth > 15 && mth < 43)
                || (mth > 43 && mth < 46)
                || (mth > 49 && mth < 92)
                || (mth > 100 && mth < 110)
                || (mth > 110 && mth < 118)
                || (mth > 150 && mth <= 200)
                || mth >= 600
            {
                //--get reaction threshold and store data
                let mut bkg = 0.0f64;
                c.gety1(0.0, &mut enext, &mut idis, &mut bkg, nendf, &mut a.all());
                if mth >= 4 {
                    jthr += 1;
                    let k = (1 + 2 * (jthr - 1)) as usize;
                    if k + 1 >= thr.len() {
                        thr.resize(k + 2, 0.0);
                    }
                    thr[k] = mth as f64;
                    thr[k + 1] = enext;
                }
                //--skip unresolved resonance reaction if required
                while MTX[ix] != 102 && MTX[ix] < mth {
                    ix += 1;
                }
                //--store data is this is one we need
                if mth == MTX[ix] {
                    let ibase = 1 + nunr * (ix as i32 - 1) - 1;
                    for ie in 1..=nunr {
                        let mut e = self.eunr[ie].abs();
                        if ie == 1 {
                            e *= UP;
                        }
                        if ie == nunr {
                            e *= DN;
                        }
                        let mut y = 0.0f64;
                        c.gety1(e, &mut enext, &mut idis, &mut y, nendf, &mut a.all());
                        self.sb[ibase + ie] = y;
                    }
                }
            }
            c.tosend(nendf, 0, 0, &mut a.all());
            c.contio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
        }

        //--output lssf flag
        c.io.so(&format!(" evaluation lssf equal to{}", ffmt::i(self.lssf, 2)));

        //--check for competitive reactions
        let sb = &mut self.sb;
        let mut icx = 0;
        for ie in 1..=nunr {
            let sigx = sb[ie] - sb[nunr + ie] - sb[2 * nunr + ie] - sb[3 * nunr + ie];
            if icx == 0 && sigx > SMALL {
                icx = ie;
            }
        }
        self.iinel = -1;
        self.iabso = -1;
        let mut ecomp = 0.0f64;
        if icx > 0 {
            ecomp = self.eunr[icx];
            c.io.so(&format!("\n competition starts at{}", ffmt::e1p(ecomp, 12, 4)));
            for kthr in 1..=jthr {
                let mtc = nint(thr[(1 + 2 * (kthr - 1)) as usize]);
                if mtc > 4 && mtc != 18 && mtc != 19 && mtc != 102 {
                    if UP * thr[(2 + 2 * (kthr - 1)) as usize] < self.eunr[nunr] {
                        c.io.so(&format!("   ur competes with mt{}", ffmt::i(mtc, 3)));
                        if mtc >= 51 && mtc <= 91 {
                            if self.iinel < 0 {
                                self.iinel = mtc;
                            } else {
                                self.iinel = 4;
                            }
                        } else {
                            if self.iabso < 0 {
                                self.iabso = mtc;
                            } else {
                                self.iabso = 0;
                            }
                        }
                    }
                }
            }
            c.io.so(&format!(
                " inelastic competition flag set to {}\n absorption competition flag set to {}",
                ffmt::i(self.iinel, 3),
                ffmt::i(self.iabso, 3)
            ));
        }

        //--sanity check for lssf=1, derive corresponding sb(i) array
        if self.lssf > 0 {
            for ie in 1..=nunr {
                let total = sb[ie];
                sb[ie] = sb[ie] - sb[ie + nunr] - sb[ie + 2 * nunr] - sb[ie + 3 * nunr];
                if sb[ie] > TOL * total {
                    //--this can only happen if there is competition and above ecomp
                    if self.iinel < 0 && self.iabso < 0 {
                        let strng1 =
                            format!("total xs greater than its components at e={}", ffmt::e1p(self.eunr[ie], 12, 4));
                        c.mess("purr", &strng1, "check evaluation file");
                        sb[ie] = 0.0;
                    } else if self.eunr[ie] < ecomp {
                        let strng1 =
                            format!("total xs greater than its components at e={}", ffmt::e1p(self.eunr[ie], 12, 4));
                        c.mess("purr", &strng1, "check evaluation file below competition");
                        sb[ie] = 0.0;
                    }
                } else {
                    //--either roundoff or total is smaller than its components
                    if sb[ie] < -TOL * total {
                        let strng1 =
                            format!("total xs less than its components at e={}", ffmt::e1p(self.eunr[ie], 12, 4));
                        c.mess("purr", &strng1, "check evaluation file");
                    }
                    sb[ie] = 0.0;
                }
                sb[ie + nunr] = 0.0;
                sb[ie + 2 * nunr] = 0.0;
                sb[ie + 3 * nunr] = 0.0;
            }
        }
    }

    // ============================================================== rdheat
    /// Read the total heating (MT=301) and partial heating cross sections
    /// for elastic (302), fission (318), and capture (402) from the pendf
    /// tape on the unresolved energy grid.
    fn rdheat(
        &mut self,
        c: &mut Ctx,
        a: &mut RV,
        heat: &mut R3,
        _temp: &[f64],
        ntemp: i32,
        nunr: i32,
        ihave: &mut i32,
        matd: i32,
    ) {
        let (mut nb, mut nw) = (0i32, 0i32);
        let (mut enext, mut idis) = (0.0f64, 0i32);
        let nin = self.nin;

        //--loop over the temperatures
        c.repoz(nin);
        c.tpidio(nin, 0, 0, &mut a.all(), &mut nb, &mut nw);
        let mut itemp = 0;
        while itemp < ntemp {
            itemp += 1;
            c.findf(matd, 3, 1, nin);
            *ihave = 0;
            c.mth = 1;
            while c.mth != 0 {
                c.contio(nin, 0, 0, &mut a.all(), &mut nb, &mut nw);
                if c.mth != 0 {
                    let mut ix = 0;
                    if c.mth == 301 {
                        ix = 1;
                    }
                    if c.mth == 302 {
                        ix = 2;
                    }
                    if c.mth == 318 {
                        ix = 3;
                    }
                    if c.mth == 402 {
                        ix = 4;
                    }
                    if ix > 0 {
                        self.e = 0.0;
                        let mut h = 0.0f64;
                        c.gety1(self.e, &mut enext, &mut idis, &mut h, nin, &mut a.all());
                        if ix == 1 {
                            *ihave = 1;
                        }
                        if ix > 1 {
                            *ihave = 2;
                        }
                        for ie in 1..=nunr {
                            self.e = self.eunr[ie].abs();
                            c.gety1(self.e, &mut enext, &mut idis, &mut h, nin, &mut a.all());
                            heat[(ix, ie, itemp)] = h;
                        }
                    }
                    c.tosend(nin, 0, 0, &mut a.all());
                }
            }
            c.tomend(nin, 0, 0, &mut a.all());
        }
        c.repoz(nin);
    }

    // ============================================================== unresx
    /// Construct ladder parameters, infinite-dilution cross sections, and
    /// potential scattering for this energy.
    fn unresx(&mut self, c: &mut Ctx, ee: f64, tt: f64) {
        const RC1: f64 = 0.123e0;
        const RC2: f64 = 0.08e0;
        const THIRD: f64 = 0.333333333e0;
        const CON1: f64 = 2901.34e0;
        let cwaven = sqrt(2.0 * AMASSN * AMU * EV) * 1.0e-12 / HBAR;

        //--initialize
        self.e = ee.abs();
        let e = self.e;
        let t = tt;
        self.spot = 0.0;
        self.dbarin = 0.0;
        self.nseq0 = 0;
        self.sigi = [0.0; 5];
        let mut lfw = 0i32;
        let mut vl = 0.0f64;
        let mut ps = 0.0f64;
        let mut ne = 0i32;
        let mut iest = 0i32;
        let mut aa: f64;
        let (mut dx, mut gxx, mut gnox, mut ggx, mut gfx) = (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
        let (mut amux, mut amuf, mut aj, mut amun): (f64, f64, f64, f64);

        //--find sections of resonance parameters which contribute
        for i in 1..=self.nsect as usize {
            let isect = i;
            let abn = self.abnt[i];
            let mode = self.modet[i];
            if (e >= self.elt[i] && e <= self.eht[i]) && (mode == 5 || mode == 6) {
                //--retrieve starting location in storage for this section
                let mut inow = self.ibaset[isect];

                //--retrieve nuclide information
                let mut ay = 0.0f64;
                if self.nro == 1 {
                    let iro = inow;
                    let mut ip = 2;
                    let mut ir = 1;
                    let (mut enx, mut idx) = (0.0f64, 0i32);
                    c.terpa(&mut ay, e, &mut enx, &mut idx, &self.arry.at(iro), &mut ip, &mut ir);
                    inow = inow + 6 + 2 * nint(self.arry[iro + 4]) + 2 * nint(self.arry[iro + 5]);
                } else {
                    ay = self.arry[inow + 1];
                }
                let awri = self.arry[inow];
                let mut aaa = 0.0f64;
                if self.naps == 2 && self.nro == 1 {
                    aaa = self.arry[inow + 1];
                }
                let spi = self.arry[inow + 2];
                let nls;
                if mode != 6 {
                    lfw = nint(self.arry[inow + 3]);
                    nls = nint(self.arry[inow + 4]);
                    ne = nint(self.arry[inow + 5]);
                    inow += 6;
                    // save starting location for fission width energies
                    if lfw != 0 {
                        iest = inow;
                        inow += ne;
                    }
                } else {
                    nls = nint(self.arry[inow + 3]);
                    inow += 4;
                }

                //--calculate nuclide dependent parameters
                let rat = awri / (awri + 1.0);
                let aw = awri * AMASSN;
                if self.naps == 0 {
                    aa = RC1 * pow(aw, THIRD) + RC2;
                } else if self.naps == 1 {
                    aa = ay;
                } else if self.naps == 2 && self.nro == 1 {
                    aa = aaa;
                } else {
                    c.error("unresx", "illegal naps", " ");
                }
                let e2 = sqrt(e);
                let k = cwaven * rat * e2;
                let ab = 4.0 * PI / powc(k, 2);
                let rho = k * aa;
                let rhoc = k * ay;

                //--loop over all sequences (l,j) in this section
                for _l in 1..=nls {
                    let njs = nint(self.arry[inow]);
                    let ll = nint(self.arry[inow + 1]);
                    inow += 2;
                    for j in 1..=njs {
                        //--retrieve parameters for this sequence
                        if mode != 6 {
                            amux = 0.0;
                            gxx = 0.0;
                            if lfw == 0 {
                                amuf = 0.0;
                                dx = self.arry[inow];
                                aj = self.arry[inow + 1];
                                amun = self.arry[inow + 2];
                                gnox = self.arry[inow + 3];
                                ggx = self.arry[inow + 4];
                                gfx = 0.0;
                                inow += 5;
                            } else {
                                amuf = self.arry[inow];
                                dx = self.arry[inow + 1];
                                aj = self.arry[inow + 2];
                                amun = self.arry[inow + 3];
                                gnox = self.arry[inow + 4];
                                ggx = self.arry[inow + 5];
                                inow += 6;
                                // interpolate for fission width at this energy
                                let ifst = inow;
                                intrf2(c, e, &mut gfx, iest, ifst, ne, 2, &self.arry);
                                inow += ne;
                            }
                        } else {
                            aj = self.arry[inow];
                            let int_ = nint(self.arry[inow + 1]);
                            ne = nint(self.arry[inow + 2]);
                            amux = self.arry[inow + 3];
                            amun = self.arry[inow + 4];
                            amuf = self.arry[inow + 6];
                            inow += 7;
                            // interpolate for parameters at this energy
                            intr2(c, e, &mut dx, &mut gxx, &mut gnox, &mut ggx, &mut gfx, inow, ne, int_, &self.arry);
                            inow += 6 * ne;
                        }

                        //--calculate penetrability and phase shifts
                        unfac2(ll, rho, rhoc, amun, &mut vl, &mut ps);
                        let gj = (2.0 * aj + 1.0) / (4.0 * spi + 2.0);
                        let mut gnx = gnox * vl * e2 * amun;

                        //--compute potential scattering
                        if j == 1 {
                            self.spot += abn * ab * (2 * ll + 1) as f64 * powc(sin(ps), 2);
                        }

                        //--store ladder parameters
                        self.dbarin += 1.0 / dx;
                        let lu = nint(amux);
                        let mu = nint(amuf);
                        let nu = nint(amun);
                        gnx /= nu as f64;
                        if mu <= 0 {
                            gfx = 0.0;
                        }
                        if lu <= 0 {
                            gxx = 0.0;
                        }
                        self.nseq0 += 1;
                        if self.nseq0 > MXNS0 as i32 {
                            c.error("unresx", "too many sequences, increase mxns0", " ");
                        }
                        let n = self.nseq0 as usize;
                        self.csz[n] = abn * ab * gj;
                        let mut tp = t;
                        if tp == 0.0 {
                            tp = 1.0;
                        }
                        self.cth[n] = sqrt(CON1 * awri / (e * tp));
                        self.cc2p[n] = cos(2.0 * ps);
                        self.cs2p[n] = sin(2.0 * ps);
                        self.cgn[n] = gnx;
                        self.cgg[n] = ggx;
                        self.cgf[n] = gfx;
                        self.cgx[n] = gxx;
                        self.cgt[n] = gnx + ggx + gfx + gxx;
                        self.ndfn[n] = nu;
                        self.ndff[n] = mu;
                        self.ndfx[n] = lu;
                        self.dbar[n] = dx;

                        //--compute infinite dilution cross sections
                        let gs = gnrx(gnx, gfx, ggx, nu, mu, lu, gxx, 1);
                        let gc = gnrx(gnx, gfx, ggx, nu, mu, lu, gxx, 2);
                        let gf = gnrx(gnx, gfx, ggx, nu, mu, lu, gxx, 3);
                        let temp = abn * PI * ab * gj * gnx / (2.0 * dx);
                        self.sigi[2] += temp * (gs * gnx - 2.0 * powc(sin(ps), 2));
                        self.sigi[3] += temp * gf * gfx;
                        self.sigi[4] += temp * gc * ggx;
                    }
                }
            }
        }
    }

    // =============================================================== ladr2
    /// Generate ladders of statistically chosen resonances from average
    /// resonance data from ENDF data.
    fn ladr2(&mut self, c: &mut Ctx, nseqz: usize, elow: f64, ehigh: f64, w: &mut Work) -> i32 {
        const CHISQ: [f64; 80] = [
            1.31003e-3, 9.19501e-3, 0.0250905e0, 0.049254e0, 0.0820892e0, 0.124169e0, 0.176268e0, 0.239417e0,
            0.314977e0, 0.404749e0, 0.511145e0, 0.637461e0, 0.788315e0, 0.970419e0, 1.194e0, 1.47573e0,
            1.84547e0, 2.36522e0, 3.20371e0, 5.58201e0, 0.0508548e0, 0.156167e0, 0.267335e0, 0.38505e0,
            0.510131e0, 0.643564e0, 0.786543e0, 0.940541e0, 1.1074e0, 1.28947e0, 1.48981e0, 1.71249e0,
            1.96314e0, 2.24984e0, 2.58473e0, 2.98744e0, 3.49278e0, 4.17238e0, 5.21888e0, 7.99146e0,
            0.206832e0, 0.470719e0, 0.691933e0, 0.901674e0, 1.10868e0, 1.31765e0, 1.53193e0, 1.75444e0,
            1.98812e0, 2.23621e0, 2.50257e0, 2.79213e0, 3.11143e0, 3.46967e0, 3.88053e0, 4.36586e0,
            4.96417e0, 5.75423e0, 6.94646e0, 10.0048e0, 0.459462e0, 0.893735e0, 1.21753e0, 1.50872e0,
            1.78605e0, 2.05854e0, 2.33194e0, 2.61069e0, 2.89878e0, 3.20032e0, 3.51995e0, 3.86331e0,
            4.23776e0, 4.65345e0, 5.12533e0, 5.67712e0, 6.35044e0, 7.22996e0, 8.541e0, 11.8359e0,
        ];
        const START: f64 = 19.9999e0;
        const ZERO: f64 = 0.0;
        let chisq = |n: i32, ndf: i32| CHISQ[((n - 1) + (ndf - 1) * 20) as usize];

        //--select resonance parameters
        //--until resonance energy is above ehi
        let dcon = self.dbar[nseqz] * sqrt(4.0 / PI);
        let nr0 = 0;
        let mut ir = nr0;
        loop {
            ir += 1;
            if ir > self.nermax {
                c.error("ladr2", "too many resonances in ladder", " ");
            }
            let iru = ir as usize;
            if ir != nr0 + 1 {
                // select next resonance location from wigner distribution
                self.rr = self.rann();
                w.er[iru] = w.er[iru - 1] + dcon * sqrt(-log(self.rr));
            } else {
                // choose starting point from uniform distribution
                self.rr = self.rann();
                w.er[iru] = elow + dcon * self.rr;
            }

            //--select level widths
            let i = nseqz;
            let gg = self.cgg[i];
            let mut gf = self.cgf[i];
            let mut gx = self.cgx[i];
            let mut gn = self.cgn[i];
            // choose neutron width for ndn degrees of freedom
            let mut ndf = self.ndfn[i];
            if ndf > 0 {
                gn /= ndf as f64;
            }
            let n = int(1.0 + START * self.rann());
            let xr = chisq(n, ndf);
            gn *= xr;
            // choose fission width for ndf degrees of freedom
            if self.cgf[i] != ZERO && self.ndff[i] != 0 {
                ndf = self.ndff[i];
                if ndf > 0 {
                    gf /= ndf as f64;
                }
                let n = int(1.0 + START * self.rann());
                let xr = chisq(n, ndf);
                gf *= xr;
            }
            // choose competitive width for ndx degrees of freedom
            if self.cgx[i] != ZERO && self.ndfx[i] != 0 {
                ndf = self.ndfx[i];
                if ndf > 0 {
                    gx /= ndf as f64;
                }
                let mut dn = 20.0f64;
                dn -= dn / 10000.0;
                let n = int(1.0 + dn * self.rann());
                let xr = chisq(n, ndf);
                gx *= xr;
            }
            w.gt[iru] = gn + gf + gg + gx;
            w.gnr[iru] = gn / w.gt[iru];
            w.gfr[iru] = gf / w.gt[iru];
            w.ggr[iru] = gg / w.gt[iru];
            w.gxr[iru] = gx / w.gt[iru];
            if w.er[iru] > ehigh {
                break;
            }
        }
        ir - 1
    }

    // ================================================================ rann
    /// Random number generator. Set kk negative to reset the seed.
    fn rann(&mut self) -> f64 {
        const M: i32 = 714025;
        const IA: i32 = 1366;
        const IC: i32 = 150889;
        let rm: f64 = 1.0 / M as f64;
        loop {
            if self.kk < 0 || self.iff == 0 {
                self.iff = 1;
                self.kk = (IC - self.kk) % M;
                for j in 1..=97 {
                    self.kk = (IA * self.kk + IC) % M;
                    self.ir[j] = self.kk;
                }
                self.kk = (IA * self.kk + IC) % M;
                self.iy = self.kk;
            }
            let j = 1 + (97 * self.iy) / M;
            if !(1..=97).contains(&j) {
                crate::fio::fatal_runtime(" ***error in rann***failed");
            }
            self.iy = self.ir[j as usize];
            let r = self.iy as f64 * rm;
            self.kk = (IA * self.kk + IC) % M;
            self.ir[j as usize] = self.kk;
            if r != 0.0 {
                return r;
            }
        }
    }

    // ============================================================== unrest
    /// Probability tables and Bondarenko moments.
    fn unrest(
        &mut self,
        c: &mut Ctx,
        bkg: &[f64; 5],
        sig0: &[f64],
        sigf: &mut R3,
        tabl: &mut R3,
        tval: &mut RM,
        w: &mut Work,
        temp: &[f64],
        nsig0: i32,
        ntemp: i32,
        nbin: i32,
        nsamp: i32,
    ) {
        const NMR: [&str; 5] = ["", "tot", "els", "fis", "cap"];
        const TREF: f64 = 300.0e0;
        const BIG: f64 = 1.0e6;
        const ZERO: f64 = 0.0;
        let rpi = sqrt(PI);
        let mut izeroprob = 0;
        let mut ibadxs = 0;
        // data for optional bondarenko plots
        if self.ipl < 200 {
            self.ipl += 1;
        }
        self.epl[self.ipl as usize] = self.e;

        //--initialize
        if self.init == 0 {
            self.init = 1;
            self.navoid = 300;
            self.elow = 10.0;
            let nt = ntemp.max(0) as usize;
            self.bval = R3::new(8, nsig0.max(0) as usize, nt);
            self.tmin = vec![0.0; nt + 1];
            self.tmax = vec![0.0; nt + 1];
            self.tsum = vec![0.0; nt + 1];
        }
        let mut dmin = 100000.0f64;
        for i in 1..=self.nseq0 as usize {
            if self.dbar[i] < dmin {
                dmin = self.dbar[i];
            }
        }
        let erange = (9 * self.nermax) as f64 * dmin / 10.0;
        let nres = int(erange * self.dbarin);
        let nrest = nres - 2 * self.navoid;
        let elow = self.elow;
        let ehigh = elow + erange;
        let emin = elow + self.navoid as f64 / self.dbarin;
        let emax = ehigh - self.navoid as f64 / self.dbarin;
        let espan = emax - emin;
        if nres <= 0 || emax < emin {
            c.error("unrest", "bad value for nres or emin>emax, increase dmin", "");
        }
        let dbart = 1.0 / self.dbarin;
        let sigx = bkg[1] - bkg[2] - bkg[3] - bkg[4];
        for itemp in 1..=ntemp {
            for i in (1..=nsig0).rev() {
                for j in 1..=8 {
                    self.bval[(j, i, itemp)] = 0.0;
                }
            }
        }

        if self.iprint > 0 {
            c.io.so(&format!(
                "\ne={}  spot={}  dbar={}  sigx={}\n{}total{}elastic{}fission{}capture",
                ffmt::e1p(self.e, 11, 4),
                ffmt::e1p(self.spot, 11, 4),
                ffmt::e1p(dbart, 11, 4),
                ffmt::e1p(sigx, 11, 4),
                " ".repeat(13),
                " ".repeat(5),
                " ".repeat(5),
                " ".repeat(5)
            ));
        }

        let (mut tav, mut tvar, mut eav, mut evar) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        let (mut fav, mut fvar, mut cav, mut cvar) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
        let spot = self.spot;

        //--loop over ladders
        //--use first pass to set prob table limits
        for iladr in 1..=self.nladr {
            let ne = nsamp;

            //--choose random energy grid
            for ie in 1..=ne {
                self.rr = self.rann();
                w.es[ie as usize] = emin + espan * self.rr;
                for itemp in 1..=ntemp {
                    w.cap[(itemp, ie)] = 0.0;
                    w.fis[(itemp, ie)] = 0.0;
                    w.els[(itemp, ie)] = spot;
                }
            }
            fsort(&mut w.es, ne);

            //--loop over sequences
            for k in 1..=self.nseq0 as usize {
                let nr = self.ladr2(c, k, elow, ehigh, w);

                //--loop over temperatures
                for itemp in 1..=ntemp {
                    let ctx = self.cth[k] * sqrt(TREF / temp[itemp as usize]);

                    //--accumulate contributions to cross sections
                    for jr in 1..=nr as usize {
                        self.accum(w, k, jr, itemp, ctx, ne, rpi);
                    }
                }
            }

            //--eliminate negative elastic cross sections; reset to 1 microbarn
            for itemp in 1..=ntemp {
                for ie in 1..=ne {
                    if w.els[(itemp, ie)] < -bkg[2] {
                        w.els[(itemp, ie)] = -bkg[2] + 1.0 / BIG;
                    }
                }
            }

            //--compute infinitely-dilute cross sections
            if iladr == 1 {
                tav = 0.0;
                tvar = 0.0;
                eav = 0.0;
                evar = 0.0;
                fav = 0.0;
                fvar = 0.0;
                cav = 0.0;
                cvar = 0.0;
            }
            let mut totf = 0.0f64;
            let mut elsf = 0.0f64;
            let mut capf = 0.0f64;
            let mut fisf = 0.0f64;
            for ie in 1..=ne {
                totf = totf + w.els[(1, ie)] + w.fis[(1, ie)] + w.cap[(1, ie)] + bkg[1];
                elsf = elsf + w.els[(1, ie)] + bkg[2];
                capf = capf + w.cap[(1, ie)] + bkg[4];
                fisf = fisf + w.fis[(1, ie)] + bkg[3];
            }
            totf /= ne as f64;
            elsf /= ne as f64;
            capf /= ne as f64;
            fisf /= ne as f64;
            tav += totf;
            tvar += totf * totf;
            eav += elsf;
            evar += elsf * elsf;
            fav += fisf;
            fvar += fisf * fisf;
            cav += capf;
            cvar += capf * capf;
            if self.iprint > 0 {
                c.io.so(&format!(
                    "{}{}{}{}{}",
                    ffmt::i(iladr, 6),
                    ffmt::e1p(totf, 12, 4),
                    ffmt::e1p(elsf, 12, 4),
                    ffmt::e1p(fisf, 12, 4),
                    ffmt::e1p(capf, 12, 4)
                ));
            }

            //--renormalize reaction cross sections
            //--compute total cross section
            if self.nmode == 1 {
                let efact = (self.sigi[2] + bkg[2]) / (elsf - spot);
                let cfact = (self.sigi[4] + bkg[4]) / capf;
                let mut ffact = 1.0;
                if fisf != ZERO {
                    ffact = (self.sigi[3] + bkg[3]) / fisf;
                }
                for ie in 1..=ne {
                    for itemp in 1..=ntemp {
                        w.els[(itemp, ie)] = spot + (w.els[(itemp, ie)] - spot) * efact;
                        w.cap[(itemp, ie)] *= cfact;
                        w.fis[(itemp, ie)] *= ffact;
                    }
                }
            }

            //--loop over temperatures
            //--using a different total cross section bin structure for each
            for itemp in 1..=ntemp {
                let it = itemp as usize;
                //--choose probability table bounds
                //--and zero accumulators
                if iladr == 1 {
                    for ie in 1..=ne {
                        w.es[ie as usize] = w.els[(itemp, ie)] + w.fis[(itemp, ie)] + w.cap[(itemp, ie)] + bkg[1];
                    }
                    fsort(&mut w.es, ne);
                    self.tmin[it] = w.es[1];
                    self.tmax[it] = w.es[ne as usize];
                    let nebin = int(nsamp as f64 / ((nbin - 10) as f64 + 1.76));
                    let mut ibin = nebin / 200;
                    if ibin <= 0 {
                        if self.mflg1 == 0 {
                            self.mflg1 = 1;
                            c.mess("purr", "reset ibin=1, consider larger nsamp or smaller nbin", "");
                        }
                        ibin = 1;
                    }
                    for i in 1..=nbin - 1 {
                        if ibin > nsamp {
                            if self.mflg2 == 0 {
                                self.mflg2 = 1;
                                c.mess("purr", "reset ibin=nsamp, consider smaller nbin", "");
                            }
                            ibin = nsamp;
                        }
                        tval[(i, itemp)] = w.es[ibin as usize];
                        if i > 1 {
                            if tval[(i, itemp)] <= tval[(i - 1, itemp)] {
                                tval[(i, itemp)] = tval[(i - 1, itemp)] + tval[(i - 1, itemp)] / 20.0;
                            }
                        }
                        if i == 1 {
                            ibin += nebin / 40;
                        }
                        if i == 2 {
                            ibin += nebin / 10;
                        }
                        if i == 3 {
                            ibin += nebin / 4;
                        }
                        if i == 4 {
                            ibin += nebin / 2;
                        }
                        if i > 4 && i < nbin - 5 {
                            ibin += nebin;
                        }
                        if i == nbin - 5 {
                            ibin += nebin / 2;
                        }
                        if i == nbin - 4 {
                            ibin += nebin / 4;
                        }
                        if i == nbin - 3 {
                            ibin += nebin / 10;
                        }
                        if i == nbin - 2 {
                            ibin += nebin / 40;
                        }
                        if i == nbin - 1 {
                            ibin += nebin / 200;
                        }
                    }
                    tval[(nbin, itemp)] = BIG;
                    for i in 1..=nbin {
                        for j in 1..=5 {
                            tabl[(i, j, itemp)] = 0.0;
                        }
                    }
                    self.tsum[it] = 0.0;
                }

                //--accumulate cross sections into tables
                let nb = nbin as usize;
                for ie in 1..=ne {
                    let tot = w.els[(itemp, ie)] + w.fis[(itemp, ie)] + w.cap[(itemp, ie)] + bkg[1];
                    if tot < self.tmin[it] {
                        self.tmin[it] = tot;
                    }
                    if tot > self.tmax[it] {
                        self.tmax[it] = tot;
                    }
                    let (mut ii, mfl) = fsrch(tot, &tval.v[(it - 1) * nb..], nbin);
                    if mfl != 2 && ii < nbin {
                        ii += 1;
                    }
                    if ii < 1 {
                        ii = 1;
                    }
                    self.tsum[it] += 1.0;
                    tabl[(ii, 1, itemp)] += 1.0;
                    tabl[(ii, 2, itemp)] += tot;
                    tabl[(ii, 3, itemp)] = tabl[(ii, 3, itemp)] + w.els[(itemp, ie)] + bkg[2];
                    tabl[(ii, 4, itemp)] = tabl[(ii, 4, itemp)] + w.fis[(itemp, ie)] + bkg[3];
                    tabl[(ii, 5, itemp)] = tabl[(ii, 5, itemp)] + w.cap[(itemp, ie)] + bkg[4];
                    for i in 1..=nsig0 {
                        let s0 = sig0[i as usize];
                        let tem = s0 / (s0 + tot);
                        let bv = &mut self.bval;
                        bv[(1, i, itemp)] += tot * tem;
                        bv[(2, i, itemp)] += (w.els[(itemp, ie)] + bkg[2]) * tem;
                        bv[(3, i, itemp)] += (w.fis[(itemp, ie)] + bkg[3]) * tem;
                        bv[(4, i, itemp)] += (w.cap[(itemp, ie)] + bkg[4]) * tem;
                        bv[(5, i, itemp)] += tot * tem * tem;
                        bv[(6, i, itemp)] += tem;
                        bv[(7, i, itemp)] += tem * tem;
                    }
                }
            }
        }

        //--write overall average cross sections
        let nladr = self.nladr as f64;
        tav /= nladr;
        let mut argt = tvar / nladr - tav * tav;
        if argt < ZERO {
            argt = 0.0;
        }
        eav /= nladr;
        let mut arge = evar / nladr - eav * eav;
        if arge < ZERO {
            arge = 0.0;
        }
        fav /= nladr;
        let mut argf = fvar / nladr - fav * fav;
        if argf < ZERO {
            argf = 0.0;
        }
        cav /= nladr;
        let mut argc = cvar / nladr - cav * cav;
        if argc < ZERO {
            argc = 0.0;
        }
        tvar = 100.0 * sqrt(argt) / tav;
        evar = 100.0 * sqrt(arge) / eav;
        if fav != ZERO {
            fvar = 100.0 * sqrt(argf) / fav;
        }
        cvar = 100.0 * sqrt(argc) / cav;
        let sigi = &mut self.sigi;
        sigi[1] = sigi[2] + sigi[3] + sigi[4] + spot + bkg[1];
        sigi[2] = sigi[2] + spot + bkg[2];
        sigi[3] += bkg[3];
        sigi[4] += bkg[4];
        let sigi = self.sigi;
        if self.iprint > 0 {
            let e4 = |v: &[f64]| v.iter().map(|x| ffmt::e1p(*x, 12, 4)).collect::<String>();
            c.io.so(&format!(
                "  bkgd{}\n  infd{}\n  aver{}\n  pcsd{}{}{}{}\n  nres{}",
                e4(&bkg[1..5]),
                e4(&sigi[1..5]),
                e4(&[tav, eav, fav, cav]),
                ffmt::f(tvar, 12, 2),
                ffmt::f(evar, 12, 2),
                ffmt::f(fvar, 12, 2),
                ffmt::f(cvar, 12, 2),
                ffmt::i(nrest, 12)
            ));
        }

        //--compute and write bondarenko table
        let bhead = format!(
            "{}temp{}sig0{}p0 total{}elastic{}fission{}capture{}p1 total",
            " ".repeat(9),
            " ".repeat(6),
            " ".repeat(4),
            " ".repeat(5),
            " ".repeat(5),
            " ".repeat(5),
            " ".repeat(4)
        );
        if self.iprint > 0 {
            c.io.so(&format!("\n bondarenko cross sections by direct sampling\n{}", bhead));
        }
        let bv = &self.bval;
        for i in 1..=nsig0 {
            for itemp in 1..=ntemp {
                sigf[(1, i, itemp)] = bv[(1, i, itemp)] / bv[(6, i, itemp)];
                sigf[(2, i, itemp)] = bv[(2, i, itemp)] / bv[(6, i, itemp)];
                sigf[(3, i, itemp)] = bv[(3, i, itemp)] / bv[(6, i, itemp)];
                sigf[(4, i, itemp)] = bv[(4, i, itemp)] / bv[(6, i, itemp)];
                sigf[(5, i, itemp)] = bv[(5, i, itemp)] / bv[(7, i, itemp)];
            }
        }
        let bline = |c: &mut Ctx, sigf: &R3| {
            for itemp in 1..=ntemp {
                for i in 1..=nsig0 {
                    let mut s = format!(
                        "   {}{}",
                        ffmt::e1p(temp[itemp as usize], 10, 3),
                        ffmt::e1p(sig0[i as usize], 10, 3)
                    );
                    for j in 1..=5 {
                        s.push_str(&ffmt::e1p(sigf[(j, i, itemp)], 12, 4));
                    }
                    c.io.so(&s);
                }
            }
        };
        if self.iprint > 0 {
            bline(c, sigf);
        }

        //--normalize and write probability table
        if self.iprint > 0 {
            c.io.so("\n probability table");
        }
        let tnorm = 0.0f64;
        let enorm = 0.0f64;
        let fnorm = 0.0f64;
        let cnorm = 0.0f64;
        for i in 1..=nbin {
            for itemp in 1..=ntemp {
                let it = itemp as usize;
                self.tmin[it] += tnorm;
                tval[(nbin, itemp)] = self.tmax[it];
                let mut denom = tabl[(i, 1, itemp)];
                if denom == ZERO {
                    denom = 1.0;
                }
                tabl[(i, 1, itemp)] /= self.tsum[it];
                tabl[(i, 2, itemp)] /= denom;
                tabl[(i, 3, itemp)] /= denom;
                tabl[(i, 4, itemp)] /= denom;
                tabl[(i, 5, itemp)] /= denom;
                tabl[(i, 2, itemp)] -= tnorm;
                tabl[(i, 3, itemp)] -= enorm;
                tabl[(i, 4, itemp)] -= fnorm;
                tabl[(i, 5, itemp)] -= cnorm;
                if tabl[(i, 1, itemp)] <= ZERO {
                    izeroprob += 1;
                }
                for j in 2..=5 {
                    if tabl[(i, j, itemp)] < ZERO {
                        ibadxs += 1;
                    }
                }
            }
        }
        for itemp in 1..=ntemp {
            let it = itemp as usize;
            let t = temp[it];
            let col = |j: i32| -> Vec<f64> { (1..=nbin).map(|i| tabl[(i, j, itemp)]).collect() };
            for ixx in 1..=4 {
                if ixx == 1 {
                    c.io.so(&format!(" tmin{}{}", ffmt::e1p(t, 11, 3), ffmt::e1p(self.tmin[it], 11, 3)));
                    let tv: Vec<f64> = (1..=nbin).map(|i| tval[(i, itemp)]).collect();
                    c.io.so(&list10(" tmax", t, &tv));
                    c.io.so(&list10(" prob", t, &col(1)));
                }
                c.io.so(&list10(&format!(" {} ", NMR[ixx as usize]), t, &col(ixx + 1)));
            }
        }
        if izeroprob > 0 {
            let strng = format!("ptable has {} zero probability bins", ffmt::i(izeroprob, 3));
            c.mess("purr", &strng, " ");
        }
        if ibadxs > 0 {
            let strng = format!("ptable has {} negative xs values", ffmt::i(ibadxs, 3));
            c.mess("purr", &strng, " ");
        }

        //--compute bondarenko cross sections from prob. table
        if self.iprint > 0 {
            c.io.so(&format!("\n bondarenko cross sections from probability table\n{}", bhead));
        }
        let bv = &mut self.bval;
        for i in 1..=nsig0 {
            for j in 1..=7 {
                for itemp in 1..=ntemp {
                    bv[(j, i, itemp)] = 0.0;
                }
            }
            let s0 = sig0[i as usize];
            for itemp in 1..=ntemp {
                for j in 1..=nbin {
                    if tabl[(j, 1, itemp)] != ZERO {
                        let den = s0 / (s0 + tabl[(j, 2, itemp)]);
                        let ttt = tabl[(j, 1, itemp)];
                        bv[(1, i, itemp)] += ttt * tabl[(j, 2, itemp)] * den;
                        bv[(2, i, itemp)] += ttt * tabl[(j, 3, itemp)] * den;
                        bv[(3, i, itemp)] += ttt * tabl[(j, 4, itemp)] * den;
                        bv[(4, i, itemp)] += ttt * tabl[(j, 5, itemp)] * den;
                        bv[(5, i, itemp)] += ttt * tabl[(j, 2, itemp)] * den * den;
                        bv[(6, i, itemp)] += ttt * den;
                        bv[(7, i, itemp)] += ttt * den * den;
                    }
                }
                sigf[(1, i, itemp)] = bv[(1, i, itemp)] / bv[(6, i, itemp)];
                sigf[(2, i, itemp)] = bv[(2, i, itemp)] / bv[(6, i, itemp)];
                sigf[(3, i, itemp)] = bv[(3, i, itemp)] / bv[(6, i, itemp)];
                sigf[(4, i, itemp)] = bv[(4, i, itemp)] / bv[(6, i, itemp)];
                sigf[(5, i, itemp)] = bv[(5, i, itemp)] / bv[(7, i, itemp)];
            }
        }
        if self.iprint > 0 {
            bline(c, sigf);
        }

        //--renormalize probability table and bondarenko
        //--cross sections to the computed infinitely-dilute values
        for itemp in 1..=ntemp {
            for i in 2..=5 {
                for j in 1..=nbin {
                    if sigf[(i - 1, 1, itemp)] != 0.0 {
                        tabl[(j, i, itemp)] = tabl[(j, i, itemp)] * sigi[(i - 1) as usize] / sigf[(i - 1, 1, itemp)];
                    }
                }
            }
        }
        let ipl = self.ipl;
        for itemp in 1..=ntemp {
            for i in (1..=nsig0).rev() {
                for j in 1..=5 {
                    let mut k = j;
                    if j == 5 {
                        k = 1;
                    }
                    if sigf[(j, 1, itemp)] != 0.0 {
                        sigf[(j, i, itemp)] = sigf[(j, i, itemp)] * sigi[k as usize] / sigf[(j, 1, itemp)];
                    }
                    self.sigpl[(ipl, j, i)] = sigf[(j, i, 1)];
                }
            }
        }
    }

    /// The body of loop 210 of unrest: add the contribution of resonance
    /// jr of sequence k to the cross sections at temperature itemp.
    fn accum(&self, w: &mut Work, k: usize, jr: usize, itemp: i32, ctx: f64, ne: i32, rpi: f64) {
        const CON1: f64 = 31.83e0;
        const CON2: f64 = 20.0e0;
        const BREAK1: f64 = 0.5e0;
        const BREAK2: f64 = 3.9e0;
        let chek1 = CON1 * w.gt[jr];
        let chek2 = CON2 / ctx;
        let mut chekn = chek1;
        if chek2 > chek1 {
            chekn = chek2;
        }
        let delr = chek1 + chekn;
        let elo = w.er[jr] - delr;
        let ehi = w.er[jr] + delr;
        let (i7, _) = fsrch(ehi, &w.es[1..], ne);
        let (i0, _) = fsrch(elo, &w.es[1..], i7);
        if i0 == i7 {
            return;
        }
        let it = i7 - i0 + 1;
        let y = ctx * w.gt[jr] / 2.0;
        let yy = y * y;
        let szy = self.csz[k] * w.gnr[jr] * rpi * y;
        let cc2 = szy * (self.cc2p[k] - 1.0 + w.gnr[jr]);
        let cs2 = szy * self.cs2p[k];
        let ccg = szy * w.ggr[jr];
        let ccf = szy * w.gfr[jr];
        for ie in i0..=i7 {
            w.xs[ie as usize] = ctx * (w.es[ie as usize] - w.er[jr]);
        }
        let cap = &mut w.cap;
        let fis = &mut w.fis;
        let els = &mut w.els;
        let xs = &w.xs;
        let mut add = |ie: i32, rew: f64, aimw: f64| {
            cap[(itemp, ie)] += ccg * rew;
            fis[(itemp, ie)] += ccf * rew;
            els[(itemp, ie)] = els[(itemp, ie)] + cc2 * rew + cs2 * aimw;
        };
        // search xs(i0:i0+it-1) and shift the result as the Fortran does
        let srch = |test: f64| -> i32 {
            let (mut i, mfl) = fsrch(test, &xs[i0 as usize..], it);
            if i == 1 && mfl == 2 {
                i = 0;
            }
            i + i0 - 1
        };

        //--accumulate cross sections using an
        //--appropriate approximation in each range
        let mut i6 = i0 - 1;
        let (mut i1, mut i2, mut i3, mut i4, mut i5) = (0i32, 0i32, 0i32, 0i32, 0i32);
        let mut lab = if y > 100.0 { 340 } else { 0 };
        loop {
            match lab {
                0 => {
                    // asymptotic term for x.gt.100 or y.gt.100
                    i1 = srch(-100.0);
                    if i1 < i0 {
                        i1 = i0 - 1;
                    }
                    lab = 240;
                    if i1 <= i0 {
                        continue;
                    }
                    for ie in i0..=i1 {
                        let (rew, aimw) = wasym(xs[ie as usize], y, yy);
                        add(ie, rew, aimw);
                    }
                    if i1 >= i7 {
                        return;
                    }
                }
                240 => {
                    // rational approximation for x.gt.6 or y.gt.6 but x and y .le.100
                    i5 = i1;
                    if y > 6.0 {
                        lab = 330;
                        continue;
                    }
                    i2 = srch(-6.0);
                    if i2 < i1 {
                        i2 = i1;
                    }
                    lab = 250;
                    if i2 <= i1 {
                        continue;
                    }
                    for ie in i1 + 1..=i2 {
                        let (rew, aimw) = wrat2(xs[ie as usize], y, yy);
                        add(ie, rew, aimw);
                    }
                    if i2 >= i7 {
                        return;
                    }
                }
                250 => {
                    // rational approximation for x.gt.3.9 or y.gt.3.0 but x and y .le.6
                    i4 = i2;
                    if y > 3.0 {
                        lab = 320;
                        continue;
                    }
                    i3 = srch(-BREAK2);
                    if i3 < i2 {
                        i3 = i2;
                    }
                    lab = 260;
                    if i3 <= i2 {
                        continue;
                    }
                    for ie in i2 + 1..=i3 {
                        let (rew, aimw) = wrat3(xs[ie as usize], y, yy);
                        add(ie, rew, aimw);
                    }
                    if i3 >= i7 {
                        return;
                    }
                }
                260 => {
                    // six point table interpolation for 0.le.abs(x).le.3.9
                    i4 = srch(BREAK2);
                    if i4 < i3 {
                        i4 = i3;
                    }
                    lab = 320;
                    if i4 <= i3 {
                        continue;
                    }
                    let is = i3 + 1;
                    let (tr, ti, tempor, joff) = if y >= BREAK1 {
                        // y.ge.0.5
                        (&self.tr, &self.ti, 10.0 * y, -3)
                    } else {
                        // y.lt.0.5
                        (&self.trs, &self.tis, 50.0 * y, 2)
                    };
                    let jj = int(tempor);
                    let j = jj + joff;
                    let n = j - 1;
                    let q = tempor - jj as f64;
                    let q2 = q * q;
                    let hq = 0.5 * q;
                    let hq2 = 0.5 * q2;
                    let a1 = hq2 - hq;
                    for ie in is..=i4 {
                        let x = xs[ie as usize];
                        let ax = x.abs();
                        let mut aki = 1.0;
                        if x < 0.0 {
                            aki = -1.0;
                        }
                        let tempor = 10.0 * ax;
                        let ii = int(tempor);
                        let i = ii + 2;
                        let p = tempor - ii as f64;
                        let p2 = p * p;
                        let hp = 0.5 * p;
                        let hp2 = 0.5 * p2;
                        let a2 = hp2 - hp;
                        let pq = p * q;
                        let a3 = 1.0 + pq - p2 - q2;
                        let a4 = hp2 - pq + hp;
                        let a5 = hq2 - pq + hq;
                        let rew = a1 * t41(tr, i, n)
                            + a2 * t41(tr, i - 1, j)
                            + a3 * t41(tr, i, j)
                            + a4 * t41(tr, i + 1, j)
                            + a5 * t41(tr, i, j + 1)
                            + pq * t41(tr, i + 1, j + 1);
                        let mut aimw = a1 * t41(ti, i, n)
                            + a2 * t41(ti, i - 1, j)
                            + a3 * t41(ti, i, j)
                            + a4 * t41(ti, i + 1, j)
                            + a5 * t41(ti, i, j + 1)
                            + pq * t41(ti, i + 1, j + 1);
                        aimw *= aki;
                        add(ie, rew, aimw);
                    }
                    if i4 >= i7 {
                        return;
                    }
                }
                320 => {
                    // rational approximation for x.gt.3.9 or y.gt.3.0 but x and y .le.6
                    i5 = srch(6.0);
                    if i5 < i4 {
                        i5 = i4;
                    }
                    lab = 330;
                    if i5 <= i4 {
                        continue;
                    }
                    for ie in i4 + 1..=i5 {
                        let (rew, aimw) = wrat3(xs[ie as usize], y, yy);
                        add(ie, rew, aimw);
                    }
                    if i5 >= i7 {
                        return;
                    }
                }
                330 => {
                    // rational approximation for x.gt.6 or y.gt.6 but x and y .le.100
                    i6 = srch(100.0);
                    if i6 < i5 {
                        i6 = i5;
                    }
                    lab = 340;
                    if i6 <= i5 {
                        continue;
                    }
                    for ie in i5 + 1..=i6 {
                        let (rew, aimw) = wrat2(xs[ie as usize], y, yy);
                        add(ie, rew, aimw);
                    }
                    if i6 >= i7 {
                        return;
                    }
                }
                _ => {
                    // asymptotic term for x.gt.100 or y.gt.100
                    for ie in i6 + 1..=i7 {
                        let (rew, aimw) = wasym(xs[ie as usize], y, yy);
                        add(ie, rew, aimw);
                    }
                    return;
                }
            }
        }
    }

    // ============================================================== uwtab2
    /// Compute the coarse and fine, real and imaginary parts of the w table.
    fn uwtab2(&mut self) {
        const TENTH: f64 = 0.1e0;
        const OH2: f64 = 0.02e0;
        let mut ax = [0.0f64; 42];
        let mut ay = [0.0f64; 28];
        let ix = |i: i32, j: i32| ((i - 1) + (j - 1) * 41) as usize;
        ax[1] = -TENTH;
        ax[2] = 0.0;
        let delx = TENTH;
        for i in 3..=41 {
            ax[i] = ax[i - 1] + delx;
        }
        ay[1] = 4.0 * TENTH;
        ay[2] = 5.0 * TENTH;
        let dely = TENTH;
        for j in 3..=27 {
            ay[j] = ay[j - 1] + dely;
        }
        for i in 2..=41 {
            let x = ax[i as usize];
            for j in 1..=27 {
                let y = ay[j as usize];
                let (rew, aimw) = uw2(x, y);
                self.tr[ix(i, j)] = rew;
                self.ti[ix(i, j)] = aimw;
            }
        }
        for j in 1..=27 {
            self.tr[ix(1, j)] = self.tr[ix(3, j)];
            self.ti[ix(1, j)] = -self.ti[ix(3, j)];
            self.ti[ix(2, j)] = 0.0;
        }
        ay[1] = -OH2;
        ay[2] = 0.0;
        let dely = OH2;
        for j in 3..=27 {
            ay[j] = ay[j - 1] + dely;
        }
        for i in 2..=41 {
            let x = ax[i as usize];
            for j in 1..=27 {
                let y = ay[j as usize];
                let (rew, aimw) = uw2(x, y);
                self.trs[ix(i, j)] = rew;
                self.tis[ix(i, j)] = aimw;
            }
        }
        for j in 1..=27 {
            self.trs[ix(1, j)] = self.trs[ix(3, j)];
            self.tis[ix(1, j)] = -self.tis[ix(3, j)];
            self.tis[ix(2, j)] = 0.0;
        }
    }
}

/// Work arrays of unrest and ladr2 (allocated per material in purr).
struct Work {
    er: Vec<f64>,
    gnr: Vec<f64>,
    gfr: Vec<f64>,
    ggr: Vec<f64>,
    gxr: Vec<f64>,
    gt: Vec<f64>,
    es: Vec<f64>,
    xs: Vec<f64>,
    cap: RM,
    fis: RM,
    els: RM,
}

impl Work {
    fn new(nermax: i32, nsamp: i32, ntemp: i32) -> Self {
        let nr = nermax as usize + 2;
        let ns = nsamp as usize + 1;
        let nt = ntemp.max(1) as usize;
        Work {
            er: vec![0.0; nr],
            gnr: vec![0.0; nr],
            gfr: vec![0.0; nr],
            ggr: vec![0.0; nr],
            gxr: vec![0.0; nr],
            gt: vec![0.0; nr],
            es: vec![0.0; ns],
            xs: vec![0.0; ns],
            cap: RM::new(nt, nsamp as usize),
            fis: RM::new(nt, nsamp as usize),
            els: RM::new(nt, nsamp as usize),
        }
    }
}

/// Asymptotic form of w for x.gt.100 or y.gt.100.
#[inline]
fn wasym(x: f64, y: f64, yy: f64) -> (f64, f64) {
    const C1: f64 = 0.5641895835e0;
    let test = x * x + yy;
    let a1 = C1 / test;
    (y * a1, x * a1)
}

/// Rational approximation of w for x.gt.6 or y.gt.6.
#[inline]
fn wrat2(x: f64, y: f64, yy: f64) -> (f64, f64) {
    const C2: f64 = 0.2752551e0;
    const C3: f64 = 2.724745e0;
    const C4: f64 = 0.5124242e0;
    const C5: f64 = 0.05176536e0;
    let a1 = x * x - yy;
    let a2 = 2.0 * x * y;
    let a3 = a2 * a2;
    let temp1 = a2 * x;
    let temp2 = a2 * y;
    let a4 = a1 - C2;
    let a5 = a1 - C3;
    let f1 = C4 / (a4 * a4 + a3);
    let f2 = C5 / (a5 * a5 + a3);
    let rew = f1 * (temp1 - a4 * y) + f2 * (temp1 - a5 * y);
    let aimw = f1 * (a4 * x + temp2) + f2 * (a5 * x + temp2);
    (rew, aimw)
}

/// Rational approximation of w for x.gt.3.9 or y.gt.3.
#[inline]
fn wrat3(x: f64, y: f64, yy: f64) -> (f64, f64) {
    const D1: f64 = 0.4613135e0;
    const D2: f64 = 0.1901635e0;
    const D3: f64 = 0.09999216e0;
    const D4: f64 = 1.7844927e0;
    const D5: f64 = 0.002883894e0;
    const D6: f64 = 5.5253437e0;
    let a1 = x * x - yy;
    let a2 = 2.0 * x * y;
    let a3 = a2 * a2;
    let temp1 = a2 * x;
    let temp2 = a2 * y;
    let a4 = a1 - D2;
    let a5 = a1 - D4;
    let a6 = a1 - D6;
    let e1 = D1 / (a4 * a4 + a3);
    let e2 = D3 / (a5 * a5 + a3);
    let e3 = D5 / (a6 * a6 + a3);
    let rew = e1 * (temp1 - a4 * y) + e2 * (temp1 - a5 * y) + e3 * (temp1 - a6 * y);
    let aimw = e1 * (a4 * x + temp2) + e2 * (a5 * x + temp2) + e3 * (a6 * x + temp2);
    (rew, aimw)
}

/// Insert a new energy into a list of energies in ascending order. Omit
/// duplicate values.
fn ilist2(e: f64, elist: &mut RV, nlist: &mut i32) {
    let mut jl = 0;
    let n = *nlist;
    let mut idone = 0;
    let mut i = 0;
    while i < n && idone == 0 {
        i += 1;
        if e <= elist[i] {
            if e != elist[i] {
                jl = n - i + 1;
                for j in 1..=jl {
                    elist[n - j + 2] = elist[n - j + 1];
                }
                elist[i] = e;
            }
            idone = 1;
        }
    }
    if idone == 0 {
        jl = 1;
        elist[n + 1] = e;
    }
    if jl > 0 {
        *nlist += 1;
    }
}

/// Calculates the penetrability factor (vl) and phase shift (ps).
fn unfac2(l: i32, rho: f64, rhoc: f64, amun: f64, vl: &mut f64, ps: &mut f64) {
    if l == 0 {
        *vl = amun;
        *ps = rhoc;
    } else if l == 1 {
        let r2 = rho * rho;
        *vl = amun * r2 / (1.0 + r2);
        *ps = rhoc - atan(rhoc);
    } else if l == 2 {
        let r2 = rho * rho;
        let r4 = r2 * r2;
        *vl = amun * r4 / (9.0 + 3.0 * r2 + r4);
        *ps = rhoc - atan(3.0 * rhoc / (3.0 - rhoc * rhoc));
    }
}

/// Interpolates fission widths for unresolved representation 1.
fn intrf2(c: &Ctx, e: f64, gfx: &mut f64, iest: i32, ifst: i32, ne: i32, int: i32, a: &RV) {
    let mut i = 1;
    let mut idone = 0;
    let mut i1 = 0;
    let mut i2 = 0;
    while i < ne && idone == 0 {
        i += 1;
        i1 = iest + i - 1;
        i2 = ifst + i - 1;
        if e < a[i1] {
            idone = 1;
        }
    }
    c.terp1(a[i1 - 1], a[i2 - 1], a[i1], a[i2], e, gfx, int);
}

/// Interpolates energy dependent unresolved resonance parameters.
fn intr2(
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
    let mut i2 = 0;
    let mut i = 1;
    let mut idone = 0;
    while i < ne && idone == 0 {
        i += 1;
        i2 = i1 + 6;
        if e < a[i2] {
            idone = 1;
        } else if i < ne {
            i1 = i2;
        }
    }
    c.terp1(a[i1], a[i1 + 1], a[i2], a[i2 + 1], e, dx, int);
    c.terp1(a[i1], a[i1 + 2], a[i2], a[i2 + 2], e, gxx, int);
    c.terp1(a[i1], a[i1 + 3], a[i2], a[i2 + 3], e, gnox, int);
    c.terp1(a[i1], a[i1 + 4], a[i2], a[i2 + 4], e, ggx, int);
    c.terp1(a[i1], a[i1 + 5], a[i2], a[i2 + 5], e, gfx, int);
    if *gxx < SMALL {
        *gxx = 0.0;
    }
    if *gfx < SMALL {
        *gfx = 0.0;
    }
}

/// Calculates fluctuation integrals for unresolved resonances using the
/// MC2-2 quadrature scheme.
fn gnrx(galpha: f64, gbeta: f64, gamma: f64, mu: i32, nu: i32, lamda: i32, df: f64, id: i32) -> f64 {
    const QW: [f64; 40] = [
        1.1120413e-1, 2.3546798e-1, 2.8440987e-1, 2.2419127e-1, 0.10967668e0, 0.030493789e0, 0.0042930874e0,
        2.5827047e-4, 4.9031965e-6, 1.4079206e-8, 0.033773418e0, 0.079932171e0, 0.12835937e0, 0.17652616e0,
        0.21347043e0, 0.21154965e0, 0.13365186e0, 0.022630659e0, 1.6313638e-5, 2.745383e-31, 3.3376214e-4,
        0.018506108e0, 0.12309946e0, 0.29918923e0, 0.33431475e0, 0.17766657e0, 0.042695894e0, 4.0760575e-3,
        1.1766115e-4, 5.0989546e-7, 1.7623788e-3, 0.021517749e0, 0.080979849e0, 0.18797998e0, 0.30156335e0,
        0.29616091e0, 0.10775649e0, 2.5171914e-3, 8.9630388e-10, 0.0,
    ];
    const QP: [f64; 40] = [
        3.0013465e-3, 7.8592886e-2, 4.3282415e-1, 1.3345267e0, 3.0481846e0, 5.8263198e0, 9.9452656e0,
        1.5782128e+1, 23.996824e0, 36.216208e0, 1.3219203e-2, 7.2349624e-2, 0.19089473e0, 0.39528842e0,
        0.74083443e0, 1.3498293e0, 2.5297983e0, 5.2384894e0, 13.821772e0, 75.647525e0, 1.0004488e-3,
        0.026197629e0, 0.14427472e0, 0.44484223e0, 1.0160615e0, 1.9421066e0, 3.3150885e0, 5.2607092e0,
        7.9989414e0, 12.072069e0, 0.013219203e0, 0.072349624e0, 0.19089473e0, 0.39528842e0, 0.74083443e0,
        1.3498293e0, 2.5297983e0, 5.2384894e0, 13.821772e0, 75.647525e0,
    ];
    const ZERO: f64 = 0.0;
    let qw = |j: usize, m: i32| QW[(j - 1) + (m as usize - 1) * 10];
    let qp = |j: usize, m: i32| QP[(j - 1) + (m as usize - 1) * 10];
    let mut s = 0.0f64;
    if galpha <= ZERO {
        return s;
    }
    if gamma <= ZERO {
        return s;
    }
    if gbeta < ZERO {
        return s;
    }
    if gbeta <= ZERO && df < ZERO {
        return s;
    }
    if gbeta <= ZERO {
        if df <= ZERO {
            for j in 1..=10 {
                let xj = qp(j, mu);
                let wj = qw(j, mu);
                if id == 1 {
                    s += wj * xj * xj / (galpha * xj + gamma);
                } else if id == 2 {
                    s += wj * xj / (galpha * xj + gamma);
                }
            }
        } else {
            for j in 1..=10 {
                let xj = qp(j, mu);
                let wj = qw(j, mu);
                for k in 1..=10 {
                    let xk = qp(k, lamda);
                    let wk = qw(k, lamda);
                    if id == 1 {
                        s += wj * wk * xj * xj / (galpha * xj + gamma + df * xk);
                    } else if id == 2 {
                        s += wj * wk * xj / (galpha * xj + gamma + df * xk);
                    }
                }
            }
        }
    } else if df <= ZERO {
        if df >= ZERO {
            for j in 1..=10 {
                let xj = qp(j, mu);
                let wj = qw(j, mu);
                for k in 1..=10 {
                    let xk = qp(k, nu);
                    let wk = qw(k, nu);
                    if id == 1 {
                        s += wj * wk * xj * xj / (galpha * xj + gbeta * xk + gamma);
                    } else if id == 2 {
                        s += wj * wk * xj / (galpha * xj + gbeta * xk + gamma);
                    } else if id == 3 {
                        s += wj * wk * xj * xk / (galpha * xj + gbeta * xk + gamma);
                    }
                }
            }
        }
    } else {
        for j in 1..=10 {
            let xj = qp(j, mu);
            let wj = qw(j, mu);
            for k in 1..=10 {
                let xk = qp(k, nu);
                let wk = qw(k, nu);
                for l in 1..=10 {
                    let xl = qp(l, lamda);
                    let wl = qw(l, lamda);
                    if id == 1 {
                        s += wj * wk * wl * xj * xj / (galpha * xj + gbeta * xk + gamma + df * xl);
                    } else if id == 2 {
                        s += wj * wk * wl * xj / (galpha * xj + gbeta * xk + gamma + df * xl);
                    } else if id == 3 {
                        s += wj * wk * wl * xj * xk / (galpha * xj + gbeta * xk + gamma + df * xl);
                    }
                }
            }
        }
    }
    s
}

/// Floating-point sort of x(1..n) into increasing order.
fn fsort(x: &mut [f64], n: i32) {
    if n <= 1 {
        return;
    }
    quicksort(x, 1, n);
}

fn quicksort(a: &mut [f64], l: i32, r: i32) {
    let mut i = l;
    let mut j = r;
    let pivot = a[((l + r) / 2) as usize];
    loop {
        while a[i as usize] < pivot {
            i += 1;
        }
        while a[j as usize] > pivot {
            j -= 1;
        }
        if i <= j {
            a.swap(i as usize, j as usize);
            i += 1;
            j -= 1;
        }
        if i > j {
            break;
        }
    }
    if l < j {
        quicksort(a, l, j);
    }
    if i < r {
        quicksort(a, i, r);
    }
}

/// Floating-point search: find i with x.ge.xarray(i) and x.lt.xarray(i+1).
/// `xarray[0]` is xarray(1). Returns (i, k): k=2 if x is below the lower
/// limit (i=1), k=3 if above the upper limit (i=n), else k=1.
fn fsrch(x: f64, xarray: &[f64], n: i32) -> (i32, i32) {
    let xa = |i: i32| xarray[(i - 1) as usize];
    if x < xa(1) {
        return (1, 2);
    }
    if x > xa(n) {
        return (n, 3);
    }
    let mut i1 = 1;
    let mut i2 = n;
    loop {
        if i1 + 1 == i2 {
            break;
        }
        if x >= xa((i1 + i2) / 2) {
            i1 = (i1 + i2) / 2;
        } else {
            i2 = (i1 + i2) / 2;
        }
    }
    (i1, 1)
}

/// Compute complex probability integral w.
fn uw2(rez: f64, aim1: f64) -> (f64, f64) {
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
                if (rew - pr).abs() >= EPS {
                    lab = 380;
                } else if rez == ZERO {
                    return (rew, 0.0);
                } else if (aimw - pim).abs() >= EPS {
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
                if (rew - pr).abs() >= EPS || (rez != ZERO && (aimw - pim).abs() >= EPS) {
                    // label 470
                    sig2p = 2.0 * ajsig;
                    lab = 430;
                } else if rez == ZERO {
                    return (rew, 0.0);
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
