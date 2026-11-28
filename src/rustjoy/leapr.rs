//! LEAPR: thermal scattering law S(alpha,beta) in ENDF-6 File 7 format,
//! ported from NJOY2016 `leapr.f90` (module leapm).
//!
//! The port keeps the Fortran structure and the order of every floating
//! point operation, so that the results are the same bit for bit. Arrays
//! keep the Fortran (1-based) indices. The places where the Fortran reads
//! or writes outside an array, or uses a value that it did not set, are
//! marked in comments.

#![allow(clippy::too_many_arguments)]
#![allow(clippy::too_many_lines)]
#![allow(clippy::collapsible_if)]
#![allow(clippy::collapsible_else_if)]
#![allow(clippy::assign_op_pattern)]
#![allow(clippy::needless_range_loop)]

use crate::endf::Ctx;
use crate::farr::{IV, R3, RV};
use crate::ffmt::{a, e, e1p, f, i};
use crate::fio::{Kind, RecW, Val};
use crate::fmath::{cos, exp, int, ipow, log, nint, powc, powi, sin, sqrt};
use crate::util::physics::{AMASSN, AMU, BK, EV, HBAR, PI};
use crate::util::sigfig;

/// Min phonon expansion order for the time warning message.
const MAXNPHON: i32 = 250;
const THERM: f64 = 0.0253;
const ZERO: f64 = 0.0;
/// Extra words past the Fortran array sizes, for writes the Fortran does
/// past an array end without a check.
const SLACK: usize = 2000;

fn seti(v: &Option<Val>, x: &mut i32) {
    match v {
        Some(Val::Int(k)) => *x = *k,
        Some(Val::Real(r)) => *x = *r as i32,
        _ => {}
    }
}

fn setr(v: &Option<Val>, x: &mut f64) {
    match v {
        Some(Val::Real(r)) => *x = *r,
        Some(Val::Int(k)) => *x = *k as f64,
        _ => {}
    }
}

/// `read(nsysi,*) (x(i),i=1,n)`: values not given are left unchanged.
fn read_reals(c: &mut Ctx, x: &mut RV, n: i32) {
    if n <= 0 {
        let _ = c.io.input.read(&[]);
        return;
    }
    let v = c.io.input.read(&vec![Kind::Real; n as usize]);
    for k in 1..=n {
        setr(&v[(k - 1) as usize], &mut x[k]);
    }
}

/// The normalization and sum-rule checks of an S(alpha,beta) column
/// (the same loop is in trans, discre and skold).
fn moments(be: impl Fn(i32) -> f64, s: impl Fn(i32) -> f64, nbeta: i32) -> (f64, f64) {
    let mut sum0 = 0.0f64;
    let mut sum1 = 0.0f64;
    let mut ff1l = 0.0f64;
    let mut ff2l = 0.0f64;
    let mut bel = 0.0f64;
    for ibeta in 1..=nbeta {
        let b = be(ibeta);
        let ff2 = s(ibeta);
        let ff1 = s(ibeta) * exp(-b);
        if ibeta > 1 {
            sum0 = sum0 + (b - bel) * (ff1l + ff2l + ff1 + ff2) / 2.0;
            sum1 = sum1 + (b - bel) * (ff2l * bel + ff2 * b - ff1l * bel - ff1 * b) / 2.0;
            ff1l = ff1;
            ff2l = ff2;
            bel = b;
        } else {
            bel = b;
            ff1l = ff1;
            ff2l = ff2;
            sum0 = 0.0;
            sum1 = 0.0;
        }
    }
    (sum0, sum1)
}

#[derive(Default)]
pub struct Leapr {
    // user's input
    nout: i32,
    iprint: i32,
    nphon: i32,
    mat: i32,
    smin: f64,
    za: f64,
    awr: f64,
    spr: f64,
    npr: i32,
    iel: i32,
    nss: i32,
    ncold: i32,
    nsk: i32,
    b7: f64,
    aws: f64,
    sps: f64,
    mss: i32,
    nalpha: i32,
    nbeta: i32,
    lat: i32,
    alpha: RV,
    beta: RV,
    twt: f64,
    c: f64,
    tbeta: f64,
    tempr: RV,
    np1: i32,
    delta1: f64,
    p1: RV,
    nd: i32,
    bdel: RV,
    adel: RV,
    nka: i32,
    dka: f64,
    ska: RV,
    cfrac: f64,
    // other global variables for the module
    naint: i32,
    nbint: i32,
    f0: f64,
    tbar: f64,
    arat: f64,
    tev: f64,
    deltab: f64,
    ssm: R3,
    ssp: R3,
    dwpix: RV,
    dwp1: RV,
    tempf: RV,
    tempf1: RV,
    // saved locals
    ni: i32,
    temp: f64,
    hz: [f64; 17],
}

impl Leapr {
    pub fn new() -> Self {
        Leapr::default()
    }

    /// Calculate S(alpha,beta) and write it in ENDF-6 File 7 format.
    pub fn leapr(&mut self, c: &mut Ctx) {
        //--initialize
        let time = c.timer();
        c.io.so(&format!(
            "\n leapr...compute thermal scattering law{}{}s",
            " ".repeat(30),
            f(time, 8, 1)
        ));
        c.io.se(&format!("\n leapr...{}{}s", " ".repeat(60), f(time, 8, 1)));

        //--read user's global input
        let v = c.io.input.read(&[Kind::Int]);
        seti(&v[0], &mut self.nout);
        let mut text = " ".repeat(80);
        let v = c.io.input.read(&[Kind::Char(80)]);
        if let Some(Val::Char(s)) = &v[0] {
            text = s.clone();
        }
        c.io.so(&format!("\n     {}", text));
        let mut ntempr = 1i32;
        self.iprint = 1;
        self.nphon = 100;
        let v = c.io.input.read(&[Kind::Int, Kind::Int, Kind::Int]);
        seti(&v[0], &mut ntempr);
        seti(&v[1], &mut self.iprint);
        seti(&v[2], &mut self.nphon);
        let mut isabt = 0i32;
        let mut ilog = 0i32;
        self.smin = 1.0e-75;
        let v = c.io.input.read(&[Kind::Int, Kind::Real, Kind::Int, Kind::Int, Kind::Real]);
        seti(&v[0], &mut self.mat);
        setr(&v[1], &mut self.za);
        seti(&v[2], &mut isabt);
        seti(&v[3], &mut ilog);
        setr(&v[4], &mut self.smin);
        c.io.so(&format!(
            "\n no. of temperatures .................. {}\n \
             print flag ........................... {}\n \
             phonon-expansion order ............... {}\n \
             endf mat number ...................... {}\n \
             za ................................... {}\n \
             isabt ................................ {}\n \
             ilog ................................. {}\n \
             smin.................. ............... {}",
            i(ntempr, 10),
            i(self.iprint, 10),
            i(self.nphon, 10),
            i(self.mat, 10),
            i(nint(self.za), 10),
            i(isabt, 10),
            i(ilog, 10),
            e1p(self.smin, 10, 3)
        ));
        if isabt != 0 {
            c.io.so(
                "\n*** Warning.  isabt=1 pendf tapes CANNOT be processed by the NJOY THERMR module ***",
            );
            c.mess("leapr", "isabt=1 pendf tapes CANNOT be", "processed by thermr.");
        }
        self.iel = 0;
        self.ncold = 0;
        self.nsk = 0;
        let v = c.io.input.read(&[Kind::Real, Kind::Real, Kind::Int, Kind::Int, Kind::Int, Kind::Int]);
        setr(&v[0], &mut self.awr);
        setr(&v[1], &mut self.spr);
        seti(&v[2], &mut self.npr);
        seti(&v[3], &mut self.iel);
        seti(&v[4], &mut self.ncold);
        seti(&v[5], &mut self.nsk);
        c.io.so(&format!(
            "\n awr for principal scatterer .......... {}\n \
             free xsec for principal scatterer .... {}\n \
             number of principal atoms ............ {}\n \
             elastic option ....................... {}\n \
             cold moderator option ................ {}\n \
             s(kappa) option ...................... {}",
            f(self.awr, 10, 3),
            f(self.spr, 10, 3),
            i(self.npr, 10),
            i(self.iel, 10),
            i(self.ncold, 10),
            i(self.nsk, 10)
        ));
        self.b7 = 0.0;
        self.aws = 0.0;
        self.sps = 0.0;
        self.mss = 0;
        let v = c.io.input.read(&[Kind::Int, Kind::Real, Kind::Real, Kind::Real, Kind::Int]);
        seti(&v[0], &mut self.nss);
        setr(&v[1], &mut self.b7);
        setr(&v[2], &mut self.aws);
        setr(&v[3], &mut self.sps);
        seti(&v[4], &mut self.mss);
        if self.nss > 0 {
            c.io.so(&format!(
                "\n secondary scatterer type ............. {}\n \
                 awr for secondary scatterer .......... {}\n \
                 free xsec for secondary scatterer .... {}\n \
                 number of secondary atoms ............ {}",
                f(self.b7, 10, 3),
                f(self.aws, 10, 3),
                f(self.sps, 10, 3),
                i(self.mss, 10)
            ));
        }
        self.lat = 0;
        let v = c.io.input.read(&[Kind::Int, Kind::Int, Kind::Int]);
        seti(&v[0], &mut self.nalpha);
        seti(&v[1], &mut self.nbeta);
        seti(&v[2], &mut self.lat);
        self.naint = 1;
        self.nbint = 1;
        let (nalpha, nbeta) = (self.nalpha, self.nbeta);
        self.alpha = RV::new(nalpha.max(0) as usize);
        self.beta = RV::new(nbeta.max(0) as usize);
        read_reals(c, &mut self.alpha, nalpha);
        read_reals(c, &mut self.beta, nbeta);

        // warn for excessive computation time
        if self.nphon > MAXNPHON {
            let strng = format!(" phonon expansion order is larger than {}", i(MAXNPHON, 3));
            c.mess("leapr", &strng, "calculation time may be excessive");
        }

        //--open the output unit
        c.openz(self.nout, 1);

        //--allocate storage for ssm (and ssp if needed)
        let (nb_, na_, nt_) = (nbeta as usize, nalpha as usize, ntempr as usize);
        self.ssm = R3::new(nb_, na_, nt_);
        self.ssp = if self.ncold != 0 { R3::new(nb_, na_, nt_) } else { R3::empty() };

        //--allocate storage for tempr, dwpix, dwp1, tempf and tempf1
        self.tempr = RV::new(nt_);
        self.dwpix = RV::new(nt_);
        self.dwp1 = RV::new(nt_);
        self.tempf = RV::new(nt_);
        self.tempf1 = RV::new(nt_);

        //--loop over desired scatterers and temperatures
        let mut isecs = 0;
        self.arat = 1.0;
        let mut idone = 0;
        while idone == 0 {
            if isecs == 0 {
                c.io.so("\n principal scatterer...");
            }
            if isecs > 0 {
                self.arat = self.aws / self.awr;
                c.io.so(&format!(
                    "\n secondary scatterer...\n input alpha values divided by{}",
                    f(self.arat, 7, 3)
                ));
            }
            for itemp in 1..=ntempr {
                let v = c.io.input.read(&[Kind::Real]);
                setr(&v[0], &mut self.temp);
                let temp = self.temp;
                self.tempr[itemp] = temp.abs();
                c.io.so(&format!("\n doing temp ={}", f(temp, 10, 2)));
                self.tev = BK * temp.abs();
                if itemp == 1 || temp >= ZERO {
                    //--read in parameters for a continuous distribution
                    let v = c.io.input.read(&[Kind::Real, Kind::Int]);
                    setr(&v[0], &mut self.delta1);
                    seti(&v[1], &mut self.ni);
                    let ni = self.ni;
                    self.p1 = RV::new(ni.max(0) as usize);
                    read_reals(c, &mut self.p1, ni);
                    self.np1 = ni;
                    let v = c.io.input.read(&[Kind::Real, Kind::Real, Kind::Real]);
                    setr(&v[0], &mut self.twt);
                    setr(&v[1], &mut self.c);
                    setr(&v[2], &mut self.tbeta);

                    //--read in oscillator data
                    let v = c.io.input.read(&[Kind::Int]);
                    seti(&v[0], &mut self.nd);
                    if self.nd > 0 {
                        let nd = self.nd;
                        self.bdel = RV::new(nd as usize);
                        self.adel = RV::new(nd as usize);
                        read_reals(c, &mut self.bdel, nd);
                        read_reals(c, &mut self.adel, nd);
                    }

                    //--read in pair correlation function, s_kappa
                    if self.nsk > 0 || self.ncold > 0 {
                        let v = c.io.input.read(&[Kind::Int, Kind::Real]);
                        seti(&v[0], &mut self.nka);
                        setr(&v[1], &mut self.dka);
                        let nka = self.nka;
                        self.ska = RV::new(nka.max(0) as usize);
                        read_reals(c, &mut self.ska, nka);
                        if self.nsk == 1 {
                            c.io.so("\n s(kappa) for vineyard method");
                        }
                        if self.nsk == 2 {
                            c.io.so("\n s(kappa) for skold method");
                        }
                        for k in 1..=nka {
                            c.io.so(&format!(
                                "{}{}",
                                e1p(self.dka * k as f64, 12, 4),
                                e1p(self.ska[k], 12, 4)
                            ));
                        }
                    }

                    //--read in coherent fraction for skold method
                    if self.nsk > 0 {
                        let v = c.io.input.read(&[Kind::Real]);
                        setr(&v[0], &mut self.cfrac);
                    }
                }

                //--continuous part of distribution
                let (np1, nphon) = (self.np1, self.nphon);
                self.contin(c, temp, itemp, np1, nphon);

                //--translational part, if any
                if self.twt > ZERO {
                    self.trans(c, itemp);
                }

                //--discrete oscillators, if any
                if self.nd > 0 {
                    self.discre(c, itemp);
                }

                //--check for special hydrogen and deuterium options
                if self.ncold > 0 {
                    self.coldh(c, itemp, temp);
                }

                //--check for skold option for correlations
                if self.nsk == 2 && self.ncold == 0 {
                    self.skold(c, itemp, temp);
                }
            }
            self.adel = RV::default();
            self.bdel = RV::default();

            //--save ssm for principal scatterer on scratch file
            if self.nss == 0 || self.b7 > ZERO || isecs > 0 {
                idone = 1;
            } else {
                isecs += 1;
                self.copys(c, ntempr);
                for itemp in 1..=ntempr {
                    self.tempf1[itemp] = self.tempf[itemp];
                    self.dwp1[itemp] = self.dwpix[itemp];
                }
            }
        }

        //--construct bragg edges if coherent was requested.
        let mut nedge = 0i32;
        let mut bragg;
        if self.iel > 0 {
            let mut maxb = 60000i32;
            bragg = RV::new(maxb as usize + SLACK);
            let emax = 5.0f64;
            coher(c, self.iel, self.npr, &mut bragg, &mut nedge, &mut maxb, emax);
        } else {
            bragg = RV::new(1 + SLACK);
        }

        //--write output in endf format
        let mut isym = 0;
        if self.ncold != 0 {
            isym = 1;
        }
        if isabt == 1 {
            isym = isym + 2;
        }
        let mscr = (8 + 2 * nalpha).max(8 + 2 * nedge).max(4000);
        let mut scr = RV::new(mscr as usize + nbeta.max(0) as usize + SLACK);
        self.endout(c, ntempr, &bragg, nedge, &mut scr, mscr, isym, ilog);

        self.p1 = RV::default();
        self.ssm = R3::empty();
        self.ssp = R3::empty();

        //--finished
        c.closz(self.nout);
        let time = c.timer();
        c.io.so(&format!("{}{}s\n {}", " ".repeat(69), f(time, 8, 1), "*".repeat(77)));
    }

    /// Column (beta values) of ssm for alpha index j and temperature k.
    fn col(&self, j: i32, k: i32) -> std::ops::Range<usize> {
        let nb = self.ssm.n1;
        let off = nb * ((j - 1) as usize + self.ssm.n2 * (k - 1) as usize);
        off..off + nb
    }

    // ================================================================ contin
    /// S(alpha,beta) at temp for the continuous phonon frequency
    /// distribution (phonon expansion).
    fn contin(&mut self, c: &mut Ctx, _temp: f64, itemp: i32, np: i32, maxn: i32) {
        let tiny = 1.0e-30f64;
        let explim = -250.0f64;
        let (nalpha, nbeta, iprint) = (self.nalpha, self.nbeta, self.iprint);

        //--write heading
        c.io.so("\n solid-type contributions to scattering law");

        //--allocate temporary arrays
        let np1 = self.np1.max(0) as usize;
        let mut p = RV::new(np1);
        let ntot = (self.nphon.max(1) as usize) * np1;
        let mut tlast = RV::new(ntot);
        let mut tnow = RV::new(ntot);
        let mut xa = RV::new(nalpha as usize);
        let mut maxt = IV::new(nbeta as usize);

        //--calculate various parameters for this temperature
        self.start(c, itemp, &mut p, np);
        let mut sc = 1.0f64;
        if self.lat == 1 {
            sc = THERM / self.tev;
        }
        let (f0, deltab, arat) = (self.f0, self.deltab, self.arat);

        //--start the phonon expansion sum with t1
        for k in 1..=np {
            tlast[k] = p[k];
        }
        // terpt(t, be) depends on beta only: once per beta, not per alpha
        let mut stk = vec![0.0f64; nbeta.max(0) as usize + 1];
        for k in 1..=nbeta {
            stk[k as usize] = terpt(&p, np, deltab, self.beta[k] * sc);
        }
        for j in 1..=nalpha {
            let al = self.alpha[j] * sc / arat;
            xa[j] = log(al * f0);
            let ex = -f0 * al + xa[j];
            let mut exx = 0.0;
            if ex > explim {
                exx = exp(ex);
            }
            for k in 1..=nbeta {
                let st = stk[k as usize];
                let mut add = st * exx;
                if add < tiny {
                    add = 0.0;
                }
                self.ssm[(k, j, itemp)] = add;
            }
        }
        let mut npl = np;

        //--do the phonon expansion sum
        for j in 1..=nbeta {
            maxt[j] = nalpha + 1;
        }
        if iprint == 2 {
            c.io.so("\n normalization check for phonon expansion");
        }
        if maxn > MAXNPHON {
            let time = c.timer();
            c.io.se(&format!("\n performing phonon expansion sum{}{}s", " ".repeat(37), f(time, 8, 1)));
        }
        // exp(-m*delta) for the convolution (the same values the Fortran
        // computes again inside the loop)
        let ntab = ntot.max(np1) + 1;
        let mut etab = vec![0.0f64; ntab];
        for m in 0..ntab {
            etab[m] = exp(-(m as f64 * deltab));
        }
        for n in 2..=maxn {
            let npn = np + npl - 1;
            let ckk = convol(&p, &tlast, &mut tnow, np, npl, npn, deltab, &etab);
            if iprint == 2 {
                c.io.so(&format!("     {}{}", i(n, 5), f(ckk, 12, 5)));
            }
            for k in 1..=nbeta {
                stk[k as usize] = terpt(&tnow, npn, deltab, self.beta[k] * sc);
            }
            for j in 1..=nalpha {
                let al = self.alpha[j] * sc / arat;
                xa[j] = xa[j] + log(al * f0 / n as f64);
                let ex = -f0 * al + xa[j];
                let mut exx = 0.0;
                if ex > explim {
                    exx = exp(ex);
                }
                for k in 1..=nbeta {
                    let st = stk[k as usize];
                    let mut add = st * exx;
                    if add < tiny {
                        add = 0.0;
                    }
                    let s = self.ssm[(k, j, itemp)] + add;
                    self.ssm[(k, j, itemp)] = s;
                    if s != ZERO && n >= maxn {
                        if add > s / 1000.0 && j < maxt[k] {
                            maxt[k] = j;
                        }
                    }
                }
            }
            for k in 1..=npn {
                tlast[k] = tnow[k];
            }
            npl = npn;
            if n % MAXNPHON == 0 {
                let time = c.timer();
                c.io.se(&format!(
                    "  {} of {} loops done for phonon expansion sum{}{}s",
                    i(n, 5),
                    i(maxn, 5),
                    " ".repeat(17),
                    f(time, 8, 1)
                ));
            }
        }
        if maxn > MAXNPHON {
            let time = c.timer();
            c.io.se(&format!("\n done with phonon expansion sum{}{}s", " ".repeat(38), f(time, 8, 1)));
        }

        //--print out start of sct range for each beta
        if iprint != 0 {
            c.io.so("\n         beta   alpha sct");
            for k in 1..=nbeta {
                if k > 1 {
                    if maxt[k] > maxt[k - 1] {
                        maxt[k] = maxt[k - 1];
                    }
                }
                if maxt[k] > nalpha {
                    c.io.so(&format!(" {}      none", f(self.beta[k] * sc, 12, 4)));
                } else {
                    c.io.so(&format!(
                        " {}{}",
                        f(self.beta[k] * sc, 12, 4),
                        f(self.alpha[maxt[k]] * sc / arat, 12, 4)
                    ));
                }
            }
        }

        //---check the moments of s(alpha,beta)
        if iprint == 1 {
            c.io.so("\n sab checks\n      alpha      norm   sum rule");
        }
        let tbeta = self.tbeta;
        let tbar = self.tbar;
        for j in 1..=nalpha {
            let mut iprt = (j - 1) % self.naint + 1;
            if j == nalpha {
                iprt = 1;
            }
            if iprt == 1 {
                let al = self.alpha[j] * sc / arat;
                if iprint == 2 {
                    c.io.so(&format!("\n  alpha={}", f(al, 10, 4)));
                    c.io.so("     beta       s(alpha,beta)      ss(alpha,beta)     ss(alpha,-beta)");
                }
                let mut bel = 0.0f64;
                let mut ff1l = 0.0f64;
                let mut ff2l = 0.0f64;
                let mut sum0 = 0.0f64;
                let mut sum1 = 0.0f64;
                for k in 1..=nbeta {
                    let mut jprt = (k - 1) % self.nbint + 1;
                    if k == nbeta {
                        jprt = 1;
                    }
                    let be = self.beta[k] * sc;
                    let alw = al * tbeta;
                    let alp = alw * tbar;
                    let ex = -powc(alw - be, 2) / (4.0 * alp);
                    let mut ssct = 0.0;
                    if ex > explim {
                        ssct = exp(ex) / sqrt(4.0 * PI * alp);
                    }
                    let mut tag = "   ";
                    if j >= maxt[k] {
                        tag = "sct";
                        self.ssm[(k, j, itemp)] = ssct;
                    }
                    let s = self.ssm[(k, j, itemp)];
                    let ff2 = s;
                    let ff1 = s * exp(-be);
                    let ff0 = s * exp(-be / 2.0);
                    if jprt == 1 && iprint == 2 && ff2 > ZERO {
                        c.io.so(&format!(
                            "{}{}{}{}    {}",
                            f(be, 10, 4),
                            e1p(ff0, 18, 5),
                            e1p(ff1, 20, 5),
                            e1p(ff2, 20, 5),
                            a(tag, 3)
                        ));
                    }
                    if k > 1 {
                        sum0 = sum0 + (be - bel) * (ff1l + ff2l + ff1 + ff2) / 2.0;
                        sum1 = sum1 + (be - bel) * (ff2l * bel + ff2 * be - ff1l * bel - ff1 * be) / 2.0;
                        ff1l = ff1;
                        ff2l = ff2;
                        bel = be;
                    } else {
                        bel = be;
                        ff1l = ff1;
                        ff2l = ff2;
                        sum0 = 0.0;
                        sum1 = 0.0;
                    }
                }
                sum0 = sum0 / (1.0 - exp(-al * f0));
                sum1 = sum1 / al / tbeta;
                if iprint == 2 {
                    c.io.so(&format!("   normalization check ={}", f(sum0, 8, 4)));
                    c.io.so(&format!("        sum rule check ={}", f(sum1, 8, 4)));
                } else if iprint == 1 {
                    c.io.so(&format!(" {}{}{}", f(al, 10, 4), f(sum0, 10, 4), f(sum1, 10, 4)));
                }
            }
        }
    }

    // ================================================================= start
    /// Integral functions of the phonon frequency distribution.
    fn start(&mut self, c: &mut Ctx, itemp: i32, p: &mut RV, _np: i32) {
        let tev = self.tev;
        let tbeta = self.tbeta;

        //--copy input spectrum into p array
        let deltab = self.delta1 / tev;
        self.deltab = deltab;
        let np1 = self.np1;
        for k in 1..=np1 {
            p[k] = self.p1[k];
        }
        let npt = np1;
        let mut u = deltab;
        let v = exp(deltab / 2.0);
        p[1] = p[2] / powc(deltab, 2);
        let mut vv = v;
        for j in 2..=np1 {
            p[j] = p[j] / (u * (vv - 1.0 / vv));
            vv = v * vv;
            u = u + deltab;
        }

        //--calculate normalizing constant, an
        let mut tau = 1.0f64;
        tau = tau / 2.0;
        let mut an = fsum(1, p, npt, tau, deltab);
        an = an / tbeta;
        for k in 1..=npt {
            p[k] = p[k] / an;
        }

        //--calculate debye-waller lambda and effective temperature
        let f0 = fsum(0, p, npt, tau, deltab);
        let tbar = fsum(2, p, npt, tau, deltab) / (2.0 * tbeta);
        self.f0 = f0;
        self.tbar = tbar;

        //--convert p(beta) into t1(beta)
        for k in 1..=npt {
            let be = deltab * (k - 1) as f64;
            p[k] = p[k] * exp(be / 2.0) / f0;
        }

        //--print out the results
        c.io.so("\n frequency distribution");
        c.io.so("        e        rho(e)      beta     rho(beta)          bigp      t1(beta)");
        for k in 1..=npt {
            let be = deltab * (k - 1) as f64;
            let ee = be * tev;
            let z = exp(be / 2.0);
            let bigp = p[k] * f0 / exp(be / 2.0);
            let rho = bigp * be * (z - 1.0 / z);
            let rhoe = rho * an;
            c.io.so(&format!(
                "{}{}{}{}{}{}",
                f(ee, 9, 5),
                e1p(rhoe, 14, 4),
                f(be, 10, 4),
                e1p(rho, 14, 4),
                e1p(bigp, 14, 4),
                e1p(p[k], 14, 4)
            ));
        }
        self.dwpix[itemp] = f0;
        self.tempf[itemp] = tbar * self.tempr[itemp];
        c.io.so(&format!(
            "\n p(beta) is scaled to{}\n\n for p-bound only\n             effective temp = {}\n        debye-waller lambda={}",
            f(tbeta, 9, 5),
            f(self.tempf[itemp], 10, 3),
            f(f0, 10, 6)
        ));
        c.io.so(&format!("          factor_ph = {}", f(tbar, 10, 6)));
    }

    // ================================================================= trans
    /// Add a translational (diffusion or free-gas) part to S(alpha,beta).
    fn trans(&mut self, c: &mut Ctx, itemp: i32) {
        let c0 = 0.4f64;
        let c1 = 1.0f64;
        let c2 = 1.42f64;
        let c3 = 0.2f64;
        let c4 = 10.0f64;
        let (nalpha, nbeta, iprint) = (self.nalpha, self.nbeta, self.iprint);
        let (twt, cc, tbeta, f0) = (self.twt, self.c, self.tbeta, self.f0);

        //--write heading for translational calculation
        let time = c.timer();
        c.io.so(&format!("\n translational part of scattering law{}{}s", " ".repeat(32), f(time, 8, 1)));
        let mut sc = 1.0f64;
        if self.lat == 1 {
            sc = THERM / self.tev;
        }

        //---allocate scratch storage
        let ndmax = nbeta.max(1000000);
        let mut betan = RV::new(nbeta as usize);
        let mut ap = RV::new(ndmax as usize + SLACK);
        let mut sd = RV::new(ndmax as usize + SLACK);
        let mut sb = RV::new(ndmax as usize + SLACK);
        let mut etab: Vec<f64> = Vec::new();
        let mut etab_delta = f64::NAN;

        let (mut pf_hit, mut pf_miss) = (0u64, 0u64);
        //--alpha loop
        if iprint == 1 {
            c.io.so("\n sab checks\n      alpha      norm   sum rule");
        }
        for ialpha in 1..=nalpha {
            let mut iprt = (ialpha - 1) % self.naint + 1;
            if ialpha == nalpha {
                iprt = 1;
            }
            let al = self.alpha[ialpha] * sc / self.arat;
            if iprt == 1 && iprint == 2 {
                c.io.so(&format!("\n alpha={}", f(al, 10, 4)));
            }

            //--choose beta interval for convolution
            let mut ded = c0 * (twt * cc * al) / sqrt(c1 + c2 * (twt * cc * al) * cc);
            if ded == ZERO {
                ded = c3 * sqrt(twt * al);
            }
            let deb = c4 * al * self.deltab;
            let mut delta = deb;
            if ded < delta {
                delta = ded;
            }
            let nu = 1;
            if iprt == 1 && iprint == 2 {
                c.io.so(&format!(
                    "\n delta d={}     delta b={}          delta={}",
                    e(ded, 18, 5),
                    e(deb, 18, 5),
                    e(delta, 18, 5)
                ));
            }

            //--make table of s-diffusion or s-free on this interval
            let nsd = self.stable(c, &mut ap, &mut sd, al, delta, iprt, nu, ndmax);
            if nsd > 1 {
                //--copy original ss(-beta) to a temporary array
                for k in 1..=nbeta {
                    betan[k] = self.beta[k] * sc;
                    ap[k] = self.ssm[(k, ialpha, itemp)];
                }
                // log(ap(j)) (or slim) as sbfill computes it for each point
                let lap: Vec<f64> =
                    (0..=nbeta).map(|k| if k == 0 || ap[k] <= ZERO { -225.0 } else { log(ap[k]) }).collect();

                //--loop over beta values
                if iprt == 1 && iprint == 2 {
                    c.io.so(
                        "\n results after convolution \n     beta       s(alpha,beta)      ss(alpha,beta)     ss(alpha,-beta)",
                    );
                }
                // With more than one thread the betas are computed ahead in
                // parallel (see `trans_pre`); sb is then set at the end of
                // the alpha to what the serial loop leaves in it.
                let nbt = nsd;
                let (pre, info) = trans_pre(&sd, nbt, delta, &lap, &betan, nbeta, ndmax, al, f0, &sb);
                let mut sbw: Option<RV> = None;
                for ibeta in 1..=nbeta {
                    let mut jprt = (ibeta - 1) % self.nbint + 1;
                    if ibeta == nbeta {
                        jprt = 1;
                    }
                    let be = betan[ibeta];
                    let s = if pre.is_empty() {
                        trans_beta(
                            c, &mut sb, &mut etab, &mut etab_delta, &sd, nbt, &mut delta, be, &lap, &betan, nbeta,
                            ibeta, ndmax, al, f0, None,
                        )
                    } else {
                        match pre[ibeta as usize] {
                            Some((sv, d0, d1)) if delta.to_bits() == d0.to_bits() => {
                                delta = d1;
                                pf_hit += 1;
                                sv
                            }
                            _ => {
                                // computed here the way a worker does it
                                pf_miss += 1;
                                let w = sbw.get_or_insert_with(|| RV::new(sb.0.len()));
                                let pt = Some((info.as_slice(), &sb));
                                trans_beta(
                                    c, w, &mut etab, &mut etab_delta, &sd, nbt, &mut delta, be, &lap, &betan, nbeta,
                                    ibeta, ndmax, al, f0, pt,
                                )
                            }
                        }
                    };

                    //--store results
                    self.ssm[(ibeta, ialpha, itemp)] = s;
                    if s != ZERO {
                        let s1 = s * exp(-be / 2.0);
                        let s2 = s * exp(-be);
                        if iprt == 1 && jprt == 1 && iprint == 2 {
                            c.io.so(&format!("{}{}{}{}", f(be, 10, 4), e1p(s1, 18, 5), e1p(s2, 20, 5), e1p(s, 20, 5)));
                        }
                    }
                }

                if !pre.is_empty() {
                    sb_final(&mut sb, &info, nbt, &lap, &betan, nbeta);
                }

                //--check moments of calculated s(alpha,beta).
                if iprt == 1 {
                    let col = self.col(ialpha, itemp);
                    let v = &self.ssm.v[col];
                    let (sum0, mut sum1) = moments(|k| betan[k], |k| v[(k - 1) as usize], nbeta);
                    sum1 = sum1 / al / (tbeta + twt);
                    if iprint == 2 {
                        c.io.so(&format!("     normalization check ={}", f(sum0, 8, 4)));
                        c.io.so(&format!("          sum rule check ={}", f(sum1, 8, 4)));
                    } else if iprint == 1 {
                        c.io.so(&format!(" {}{}{}", f(al, 10, 4), f(sum0, 10, 4), f(sum1, 10, 4)));
                    }
                }
            }
        }

        if crate::util::pf_stats() {
            eprintln!("leapr trans: {} betas from workers, {} computed here", pf_hit, pf_miss);
        }

        //--update effective temperature
        self.tempf[itemp] = (tbeta * self.tempf[itemp] + twt * self.tempr[itemp]) / (tbeta + twt);
        c.io.so(&format!("\n     new effective temp = {}", f(self.tempf[itemp], 10, 3)));
    }

    // ================================================================ stable
    /// Table of S-diffusion or S-free (the -beta side) at intervals delta.
    /// Returns nsd.
    fn stable(
        &self,
        c: &mut Ctx,
        ap: &mut RV,
        sd: &mut RV,
        al: f64,
        delta: f64,
        iprt: i32,
        nu: i32,
        ndmax: i32,
    ) -> i32 {
        let quart = 0.25f64;
        let eps = 1.0e-7f64;
        let one = 1.0f64;
        let iprint = self.iprint;
        let twt = self.twt;
        let cc = self.c;
        let mut j;

        //--diffusion branch
        if cc != ZERO {
            if iprt == 1 && iprint == 2 {
                c.io.so(&format!("\n{}", "    beta    ssdiff(-beta)".repeat(4)));
            }
            let d = twt * cc;
            let c2 = sqrt(cc * cc + quart);
            let mut c3 = 2.0 * d * al;
            let c4 = c3 * c3;
            let c8 = c2 * c3 / PI;
            c3 = 2.0 * d * cc * al;
            let mut be = 0.0f64;
            j = 1i32;
            let mut idone = 0;
            let mut c5 = 0.0f64;
            while idone == 0 {
                let c6 = sqrt(be * be + c4);
                let c7 = c6 * c2;
                if c7 <= one {
                    c5 = c8 * exp(c3 + be / 2.0);
                }
                if c7 > one {
                    let ex = c3 - c7 + be / 2.0;
                    c5 = c8 * exp(ex);
                }
                sd[j] = c5 * besk1(c7) / c6;
                ap[j] = be;
                be = be + delta;
                j = j + 1;
                if j % 2 == 0 {
                    if j >= ndmax {
                        idone = 1;
                    }
                    if eps * sd[1] >= sd[j - 1] {
                        idone = 1;
                    }
                }
            }
            j = j - 1;

        //--free-gas branch
        } else {
            if iprt == 1 && iprint == 2 {
                c.io.so(&format!("\n{}", "    beta    ssfree(-beta)".repeat(4)));
            }
            let mut be = 0.0f64;
            j = 1;
            let wal = twt * al;
            let mut idone = 0;
            while idone == 0 {
                let ex = -powc(wal - be, 2) / (4.0 * wal);
                let sfree = exp(ex) / sqrt(4.0 * PI * wal);
                sd[j] = sfree;
                ap[j] = be;
                be = be + delta;
                j = j + 1;
                if j % 2 == 0 {
                    if j >= ndmax {
                        idone = 1;
                    }
                    if eps * sd[1] >= sd[j - 1] {
                        idone = 1;
                    }
                }
            }
            j = j - 1;
        }
        if iprt == 1 && iprint == 2 {
            let mut line = String::new();
            let mut cnt = 0;
            let mut k = 1;
            while k <= j {
                line.push_str(&f(ap[k], 10, 5));
                line.push_str(&e1p(sd[k], 15, 7));
                cnt += 1;
                if cnt == 4 {
                    c.io.so(&line);
                    line.clear();
                    cnt = 0;
                }
                k += nu;
            }
            if cnt > 0 {
                c.io.so(&line);
            }
        }
        // the moment check is off in the Fortran (icheck=0)
        j
    }

    // ================================================================ discre
    /// Convolve discrete oscillators with the continuous S(alpha,beta).
    fn discre(&mut self, c: &mut Ctx, itemp: i32) {
        let small = 1.0e-8f64;
        let vsmall = 1.0e-10f64;
        let tiny = 1.0e-20f64;
        let (nalpha, nbeta, iprint, nd) = (self.nalpha, self.nbeta, self.iprint, self.nd);
        let (twt, tbeta) = (self.twt, self.tbeta);
        let tev = self.tev;

        //--write heading
        let time = c.timer();
        c.io.so(&format!(
            "\n discrete-oscillator part of scattering law{}{}s",
            " ".repeat(26),
            f(time, 8, 1)
        ));
        let mut sc = 1.0f64;
        if self.lat == 1 {
            sc = THERM / tev;
        }

        //--allocate scratch storage
        let mut betan = RV::new(nbeta as usize);
        let mut exb = RV::new(nbeta as usize);
        let mut sexpb = RV::new(nbeta as usize + 1);
        let maxbb = 2 * nbeta + 1;
        let mut bex = RV::new(maxbb as usize);
        let mut rdbex = RV::new(maxbb as usize);
        let mut sex = RV::new(maxbb as usize);
        let maxdd = 500i32;
        let mut bes = RV::new(maxdd as usize);
        let mut wts = RV::new(maxdd as usize);
        let mut ben = RV::new(maxdd as usize);
        let mut wtn = RV::new(maxdd as usize);
        let mut bdeln = [0.0f64; 51];
        let mut eb = [0.0f64; 51];
        let mut dbw = [0.0f64; 51];
        let mut ar = [0.0f64; 51];
        let mut dist = [0.0f64; 51];
        let mut bplus = [0.0f64; 51];
        let mut bminus = [0.0f64; 51];
        let mut bzero = 0.0f64;

        //--set up oscillator parameters
        let mut dwt = 0.0f64;
        for k in 1..=nd as usize {
            bdeln[k] = self.bdel[k as i32] / tev;
            dwt = dwt + self.adel[k as i32];
        }
        let mut tsave = 0.0f64;
        let dw0 = self.dwpix[itemp];
        for k in 1..=nd as usize {
            let ad = self.adel[k as i32];
            eb[k] = exp(bdeln[k] / 2.0);
            let sn = (eb[k] - 1.0 / eb[k]) / 2.0;
            let cn = (eb[k] + 1.0 / eb[k]) / 2.0;
            ar[k] = ad / (sn * bdeln[k]);
            dist[k] = ad * self.bdel[k as i32] * cn / (2.0 * sn);
            tsave = tsave + dist[k] / BK;
            dbw[k] = ar[k] * cn;
            if self.dwpix[itemp] > ZERO {
                self.dwpix[itemp] = self.dwpix[itemp] + dbw[k];
            }
        }
        c.io.so("\n add delta functions\n");
        for k in 1..=nd as usize {
            c.io.so(&format!("     {}{}{}", i(k as i32, 3), e1p(bdeln[k], 14, 4), e1p(self.adel[k as i32], 14, 4)));
        }

        //--prepare functions of beta
        for k in 1..=nbeta {
            let be = self.beta[k] * sc;
            exb[k] = exp(-be / 2.0);
            betan[k] = be;
        }
        let nbx = bfill(&mut bex, &mut rdbex, &betan, nbeta);
        // wt and tbart are not reset for each alpha (as in the Fortran)
        let mut _wt = tbeta;
        let mut tbart = self.tempf[itemp] / self.tempr[itemp];

        //--main alpha loop
        if iprint == 1 {
            c.io.so("\n sab checks\n      alpha      norm   sum rule");
        }
        for nal in 1..=nalpha {
            let mut iprt = (nal - 1) % self.naint + 1;
            if nal == nalpha {
                iprt = 1;
            }
            let al = self.alpha[nal] * sc / self.arat;
            if iprt == 1 && iprint == 2 {
                c.io.so(&format!("\n   alpha={}", f(al, 10, 5)));
            }
            let dwf = exp(-al * dw0);
            if iprt == 1 && iprint == 2 {
                c.io.so(&format!("\n      debye-waller factor={}", e1p(dwf, 12, 4)));
            }
            let col = self.col(nal, itemp);
            exts(&self.ssm.v[col.clone()], &mut sex, &exb, &betan, nbeta);
            for j in 1..=nbeta {
                sexpb[j] = 0.0;
            }

            //---initialize for delta function calculation
            ben[1] = 0.0;
            wtn[1] = 1.0;
            let mut nn = 1i32;
            let mut n = 0i32;

            //--loop over all oscillators
            for k in 1..=nd as usize {
                let dwc = al * dbw[k];
                let x = al * ar[k];
                bfact(x, &mut bzero, &mut bplus, &mut bminus, dwc, bdeln[k]);

                //--do convolution for the delta functions
                //--n=0 term
                for m in 1..=nn {
                    let besn = ben[m];
                    let wtsn = wtn[m] * bzero;
                    if besn <= ZERO || wtsn >= small {
                        if n < maxdd {
                            n = n + 1;
                            bes[n] = besn;
                            wts[n] = wtsn;
                        }
                    }
                }

                //--negative n terms
                let mut kk = 0usize;
                let mut idone = 0;
                while kk < 50 && idone == 0 {
                    kk = kk + 1;
                    if bminus[kk] <= ZERO {
                        idone = 1;
                    } else {
                        for m in 1..=nn {
                            let besn = ben[m] - kk as f64 * bdeln[k];
                            let wtsn = wtn[m] * bminus[kk];
                            if wtsn >= small && n < maxdd {
                                n = n + 1;
                                bes[n] = besn;
                                wts[n] = wtsn;
                            }
                        }
                    }
                }

                //--positive n terms
                kk = 0;
                idone = 0;
                while kk < 50 && idone == 0 {
                    kk = kk + 1;
                    if bplus[kk] <= ZERO {
                        idone = 1;
                    } else {
                        for m in 1..=nn {
                            let besn = ben[m] + kk as f64 * bdeln[k];
                            let wtsn = wtn[m] * bplus[kk];
                            if wtsn >= small && n < maxdd {
                                n = n + 1;
                                bes[n] = besn;
                                wts[n] = wtsn;
                            }
                        }
                    }
                }

                //--continue oscillator loop
                nn = n;
                for m in 1..=nn {
                    ben[m] = bes[m];
                    wtn[m] = wts[m];
                }
                n = 0;
                _wt = _wt + self.adel[k as i32];
                tbart = tbart + dist[k] / BK / self.tempr[itemp];
            }
            n = nn;

            //--sort the discrete lines
            //--and throw out the smallest ones
            nn = n - 1;
            for i1 in 2..=n - 1 {
                for j in i1 + 1..=n {
                    if wts[j] >= wts[i1] {
                        let save = wts[j];
                        wts[j] = wts[i1];
                        wts[i1] = save;
                        let save = bes[j];
                        bes[j] = bes[i1];
                        bes[i1] = save;
                    }
                }
            }
            let mut i1 = 0;
            let mut idone = 0;
            while i1 < nn && idone == 0 {
                i1 = i1 + 1;
                n = i1;
                if wts[i1] < 100.0 * small && i1 > 5 {
                    idone = 1;
                }
            }

            //--report the discrete lines
            if iprint >= 2 {
                c.io.so("\n      discrete lines\n              beta      weight");
                let mut sumn = 0.0f64;
                let mut sumr = 0.0f64;
                for m in 1..=n {
                    c.io.so(&format!("      {}{}", f(bes[m], 12, 4), e1p(wts[m], 12, 4)));
                    sumn = sumn + wts[m];
                    sumr = sumr - bes[m] * wts[m];
                }
                sumr = sumr / al / dwt;
                c.io.so(&format!("        norm check{}\n        rule check{}", f(sumn, 10, 4), f(sumr, 10, 4)));
            }

            //--add the continuum part to the scattering law
            for m in 1..=n {
                for j in 1..=nbeta {
                    let be = -betan[j] - bes[m];
                    let st = sint(be, &bex, &rdbex, &sex, nbx, al, tbeta + twt, tbart, &betan, nbeta);
                    let add = wts[m] * st;
                    if add >= tiny {
                        sexpb[j] = sexpb[j] + add;
                    }
                }
            }

            //--add the delta functions to the scattering law
            //--delta(0.) is saved for the incoherent elastic
            if twt <= ZERO {
                let mut m = 0;
                let mut idone = 0;
                while m < n && idone == 0 {
                    m = m + 1;
                    if dwf < vsmall {
                        idone = 1;
                    } else {
                        if bes[m] < ZERO {
                            let be = -bes[m];
                            if be <= betan[nbeta - 1] {
                                let mut db = 1000.0f64;
                                idone = 0;
                                let mut j = 0;
                                let mut jj = 0;
                                while j < nbeta && idone == 0 {
                                    j = j + 1;
                                    jj = j;
                                    if (be - betan[j]).abs() > db {
                                        idone = 1;
                                    } else {
                                        db = (be - betan[j]).abs();
                                    }
                                }
                                let mut add;
                                if jj <= 2 {
                                    add = wts[m] / betan[jj];
                                } else {
                                    add = 2.0 * wts[m] / (betan[jj] - betan[jj - 2]);
                                }
                                add = add * dwf;
                                // (jj=1 would write sexpb(0) in the Fortran)
                                if add >= tiny && jj >= 2 {
                                    sexpb[jj - 1] = sexpb[jj - 1] + add;
                                }
                            }
                        }
                    }
                }
            }

            //--record the results
            for j in 1..=nbeta {
                self.ssm[(j, nal, itemp)] = sexpb[j];
            }
            if iprt == 1 && iprint == 2 {
                c.io.so("\n     beta       s(alpha,beta)       ss(alpha,beta)     ss(alpha,-beta)");
            }
            for k in 1..=nbeta {
                let be = self.beta[k] * sc;
                let ss = self.ssm[(k, nal, itemp)];
                let s1 = ss * exp(-be / 2.0);
                let s2 = ss * exp(-be);
                let mut jprt = (k - 1) % self.nbint + 1;
                if k == nbeta {
                    jprt = 1;
                }
                if iprt == 1 && jprt == 1 && iprint == 2 {
                    c.io.so(&format!("{}{}{}{}", f(betan[k], 10, 4), e1p(s1, 18, 5), e1p(s2, 20, 5), e1p(ss, 20, 5)));
                }
            }

            //--check moments of calculated s(alpha,beta).
            if iprt == 1 {
                let v = &self.ssm.v[col];
                let (mut sum0, mut sum1) = moments(|k| betan[k], |k| v[(k - 1) as usize], nbeta);
                if twt == ZERO {
                    sum0 = sum0 / (1.0 - exp(-al * self.dwpix[itemp]));
                }
                sum1 = sum1 / al;
                if iprint == 2 {
                    c.io.so(&format!("     normalization check ={}", f(sum0, 8, 4)));
                    c.io.so(&format!("          sum rule check ={}", f(sum1, 8, 4)));
                } else if iprint == 1 {
                    c.io.so(&format!(" {}{}{}", f(al, 10, 4), f(sum0, 10, 4), f(sum1, 10, 4)));
                }
            }
        }

        //--finished
        self.tempf[itemp] = (tbeta + twt) * self.tempf[itemp] + tsave;
        c.io.so(&format!(
            "\n       new effective temp = {}\n  new debye-waller lambda={}",
            f(self.tempf[itemp], 10, 3),
            f(self.dwpix[itemp], 10, 6)
        ));
        c.io.so(&format!(" discr.-oscill. part of eff. temp = {}", f(tsave, 10, 3)));
    }

    // ================================================================= coldh
    /// Convolve S(alpha,beta) with the rotational modes of ortho or para
    /// hydrogen or deuterium (Young and Koppel).
    fn coldh(&mut self, c: &mut Ctx, itemp: i32, temp: f64) {
        let pmass = 1.6726231e-24f64;
        let dmass = 3.343586e-24f64;
        let deh = 0.0147f64;
        let ded = 0.0074f64;
        let sampch = 0.356f64;
        let sampcd = 0.668f64;
        let sampih = 2.526f64;
        let sampid = 0.403f64;
        let small = 1.0e-6f64;
        let angst = 1.0e-8f64;
        let ifree = 0;
        let nokap = 0;
        let jterm = 3;
        let (nalpha, nbeta, iprint) = (self.nalpha, self.nbeta, self.iprint);

        //--write header
        let time = c.timer();
        c.io.so(&format!(
            "\n cold hydrogen or deuterium scattering{}{}s",
            " ".repeat(31),
            f(time, 8, 1)
        ));

        //--allocate scratch storage
        let mut betan = RV::new(nbeta as usize);
        let mut exb = RV::new(nbeta as usize);
        let maxbb = 2 * nbeta + 1;
        let mut bex = RV::new(maxbb as usize);
        let mut rdbex = RV::new(maxbb as usize);
        let mut sex = RV::new(maxbb as usize);
        let mut nbx = 0;

        //--set up constants
        let tev = BK * temp.abs();
        let mut sc = 1.0f64;
        if self.lat == 1 {
            sc = THERM / tev;
        }
        let law = self.ncold + 1;
        let mut de = deh;
        if law > 3 {
            de = ded;
        }
        let x = de / tev;
        let (amassm, sampc, bp, sampi);
        if law > 3 {
            amassm = 6.69e-24f64;
            sampc = sampcd;
            bp = HBAR / 2.0 * sqrt(2.0 / ded / EV / dmass) / angst;
            sampi = sampid;
        } else {
            amassm = 3.3464e-24f64;
            sampc = sampch;
            bp = HBAR / 2.0 * sqrt(2.0 / deh / EV / pmass) / angst;
            sampi = sampih;
        }
        let wt = self.twt + self.tbeta;
        let tbart = self.tempf[itemp] / self.tempr[itemp];
        let (mut swe, mut swo) = (0.0f64, 0.0f64);

        //--main alpha loop
        for nal in 1..=nalpha {
            let mut iprt = (nal - 1) % self.naint + 1;
            if nal == nalpha {
                iprt = 1;
            }
            let al = self.alpha[nal] * sc / self.arat;
            let alp = wt * al;
            let waven = angst * sqrt(amassm * tev * EV * al) / HBAR;
            let y = bp * waven;
            if iprt == 1 {
                c.io.so(&format!("\n\n{}   alpha={}", i(nal, 4), f(al, 10, 5)));
            }
            let mut sk = terpk(&self.ska, self.nka, self.dka, waven);
            if nokap == 1 {
                sk = 1.0;
            }
            if iprt == 1 && iprint == 2 {
                c.io.so(&format!("\n wave number             ={}", f(waven, 10, 4)));
                c.io.so(&format!(" static structure factor ={}", f(sk, 10, 4)));
                c.io.so(" oscillators:");
            }

            //--spin-correlation factors
            if law == 2 {
                swe = powc(sampi, 2) / 3.0;
            }
            if law == 2 {
                swo = sk * powc(sampc, 2) + 2.0 * powc(sampi, 2) / 3.0;
            }
            if law == 3 {
                swe = sk * powc(sampc, 2);
            }
            if law == 3 {
                swo = powc(sampi, 2);
            }
            if law == 4 {
                swe = sk * powc(sampc, 2) + 5.0 * powc(sampi, 2) / 8.0;
            }
            if law == 4 {
                swo = 3.0 * powc(sampi, 2) / 8.0;
            }
            if law == 5 {
                swe = 3.0 * powc(sampi, 2) / 4.0;
            }
            if law == 5 {
                swo = sk * powc(sampc, 2) + powc(sampi, 2) / 4.0;
            }
            let snorm = powc(sampi, 2) + powc(sampc, 2);
            swe = swe / snorm;
            swo = swo / snorm;

            //--prepare arrays for sint
            if nal == 1 {
                for k in 1..=nbeta {
                    let mut be = self.beta[k];
                    if self.lat == 1 {
                        be = be * THERM / tev;
                    }
                    exb[k] = exp(-be / 2.0);
                    betan[k] = be;
                }
                nbx = bfill(&mut bex, &mut rdbex, &betan, nbeta);
            }
            let col = self.col(nal, itemp);
            exts(&self.ssm.v[col], &mut sex, &exb, &betan, nbeta);

            //--loop over all beta values
            //    results for positive beta go into ssp
            //    results for negative beta go into ssm
            let jjmax = 2 * nbeta - 1;
            let mut k = 0;
            for jj in 1..=jjmax {
                if jj < nbeta {
                    k = nbeta - jj + 1;
                }
                if jj >= nbeta {
                    k = jj - nbeta + 1;
                }
                let mut be = betan[k];
                if jj < nbeta {
                    be = -be;
                }
                let mut sn = 0.0f64;
                let mut total = 0.0f64;

                //--loop over all oscillators
                // para-h2: j=0,2,....; ortho-h2: j=1,3,....
                // ortho-d2: j=0,2,....; para-d2: j=1,3,....
                let mut ipo = 1;
                if law == 2 || law == 5 {
                    ipo = 2;
                }
                let mut jt1 = 2 * jterm;
                if ipo == 2 {
                    jt1 = jt1 + 1;
                }
                let mut l = ipo;
                while l <= jt1 {
                    let j = l - 1;
                    let pj = bt(j, x);

                    //--sum over even values of j-prime
                    let mut snlg = 0.0f64;
                    let mut lp = 1;
                    while lp <= 10 {
                        let jp = lp - 1;
                        let betap = (-j * (j + 1) + jp * (jp + 1)) as f64 * x / 2.0;
                        let tmp = (2 * jp + 1) as f64 * pj * swe * 4.0 * sumh(c, j, jp, y);
                        if jj == 1 && tmp >= small {
                            c.io.so(&format!("     {}{}{}{}", f(betap, 10, 4), i(j, 4), i(jp, 4), f(tmp, 10, 6)));
                            total = total + tmp;
                        }
                        let bn = be + betap;
                        let add;
                        if ifree == 1 {
                            let mut ex = -powc(alp - bn.abs(), 2) / (4.0 * alp);
                            if bn > ZERO {
                                ex = ex - bn;
                            }
                            add = exp(ex) / sqrt(4.0 * PI * alp);
                        } else {
                            add = sint(bn, &bex, &rdbex, &sex, nbx, al, wt, tbart, &betan, nbeta);
                        }
                        snlg = snlg + tmp * add;
                        lp += 2;
                    }

                    //--sum over the odd values of j-prime
                    let mut snlk = 0.0f64;
                    let mut lp = 2;
                    while lp <= 10 {
                        let jp = lp - 1;
                        let betap = (-j * (j + 1) + jp * (jp + 1)) as f64 * x / 2.0;
                        let tmp = (2 * jp + 1) as f64 * pj * swo * 4.0 * sumh(c, j, jp, y);
                        if jj == 1 && tmp >= small {
                            c.io.so(&format!("     {}{}{}{}", f(betap, 10, 4), i(j, 4), i(jp, 4), f(tmp, 10, 6)));
                            total = total + tmp;
                        }
                        let bn = be + betap;
                        let add;
                        if ifree == 1 {
                            let mut ex = -powc(alp - bn.abs(), 2) / (4.0 * alp);
                            if bn > ZERO {
                                ex = ex - bn;
                            }
                            add = exp(ex) / sqrt(4.0 * PI * alp);
                        } else {
                            add = sint(bn, &bex, &rdbex, &sex, nbx, al, wt, tbart, &betan, nbeta);
                        }
                        snlk = snlk + tmp * add;
                        lp += 2;
                    }

                    //--continue the j loop
                    sn = sn + snlg + snlk;
                    l += 2;
                }
                if jj == 1 && iprt == 1 {
                    c.io.so(&format!("     total     {}", f(total, 10, 6)));
                }

                //--continue the beta loop
                if jj <= nbeta {
                    self.ssm[(k, nal, itemp)] = sn;
                }
                if jj >= nbeta {
                    self.ssp[(k, nal, itemp)] = sn;
                }
            }

            //--record the results
            if iprt == 1 {
                c.io.so(
                    "\n     beta       s(alpha,beta)       s(alpha,-beta)       ss(alpha,beta)     ss(alpha,-beta)",
                );
            }
            for k in 1..=nbeta {
                let mut jprt = (k - 1) % self.nbint + 1;
                if k == nbeta {
                    jprt = 1;
                }
                let down = self.ssm[(k, nal, itemp)] * exb[k];
                let mut up = 0.0;
                if exb[k] != ZERO {
                    up = self.ssp[(k, nal, itemp)] / exb[k];
                }
                if iprt == 1 && jprt == 1 {
                    c.io.so(&format!(
                        "{}{}{}{}{}",
                        f(betan[k], 10, 4),
                        e1p(up, 18, 5),
                        e1p(down, 20, 5),
                        e1p(self.ssp[(k, nal, itemp)], 20, 5),
                        e1p(self.ssm[(k, nal, itemp)], 20, 5)
                    ));
                }
            }

            //--check moments of calculated s(alpha,beta).
            let mut sum0 = 0.0f64;
            let mut bel = 0.0f64;
            let mut ff1l = 0.0f64;
            let mut ff2l = 0.0f64;
            for nbe in 1..=nbeta {
                let be = betan[nbe];
                let ff2 = self.ssm[(nbe, nal, itemp)];
                let ff1 = self.ssp[(nbe, nal, itemp)];
                if nbe != 1 {
                    sum0 = sum0 + (be - bel) * (ff1l + ff2l + ff1 + ff2) / 2.0;
                    ff1l = ff1;
                    ff2l = ff2;
                    bel = be;
                } else {
                    bel = be;
                    ff1l = ff1;
                    ff2l = ff2;
                    sum0 = 0.0;
                }
            }
            c.io.so(&format!("     normalization check ={}", f(sum0, 8, 4)));
        }
    }

    // ================================================================= copys
    /// Copy sab for the principal scatterer to the scratch tape on unit 10.
    fn copys(&mut self, c: &mut Ctx, ntempr: i32) {
        let nscr = -10;
        c.openz(nscr, 1);
        for k in 1..=ntempr {
            for j in 1..=self.nalpha {
                let col = self.col(j, k);
                let mut r = RecW::new();
                r.rs(&self.ssm.v[col]);
                c.io.write_rec(-nscr, r);
            }
        }
    }

    // ================================================================= skold
    /// Skold approximation for intermolecular coherence.
    fn skold(&mut self, c: &mut Ctx, itemp: i32, temp: f64) {
        let angst = 1.0e-8f64;
        let (nalpha, nbeta, iprint) = (self.nalpha, self.nbeta, self.iprint);
        let mut scoh = vec![0.0f64; 1001.max(nalpha as usize + 1)];

        //--apply the skold approximation
        let tev = BK * temp.abs();
        let mut sc = 1.0f64;
        if self.lat == 1 {
            sc = THERM / tev;
        }
        let amass = self.awr * AMASSN * AMU;
        for k in 1..=nbeta {
            for j in 1..=nalpha {
                let al = self.alpha[j] * sc / self.arat;
                let waven = angst * sqrt(2.0 * amass * tev * EV * al) / HBAR;
                let sk = terpk(&self.ska, self.nka, self.dka, waven);
                let ap = self.alpha[j] / sk;
                let mut kk = 0;
                for m in 1..=nalpha {
                    kk = m;
                    if ap < self.alpha[m] {
                        break;
                    }
                }
                if kk == 1 {
                    kk = 2;
                }
                let s1 = self.ssm[(k, kk - 1, itemp)];
                let s2 = self.ssm[(k, kk, itemp)];
                if s1 == ZERO || s2 == ZERO {
                    scoh[j as usize] = ZERO;
                } else {
                    let mut y = 0.0;
                    c.terp1(self.alpha[kk - 1], s1, self.alpha[kk], s2, ap, &mut y, 5);
                    scoh[j as usize] = y;
                }
                scoh[j as usize] = scoh[j as usize] * sk;
            }
            for j in 1..=nalpha {
                self.ssm[(k, j, itemp)] = (1.0 - self.cfrac) * self.ssm[(k, j, itemp)] + self.cfrac * scoh[j as usize];
            }
        }

        //--report the results
        if iprint == 2 {
            c.io.so("\n results after applying skold approximation");
        }
        for nal in 1..=nalpha {
            let mut iprt = (nal - 1) % self.naint + 1;
            if nal == nalpha {
                iprt = 1;
            }
            let al = self.alpha[nal] * sc / self.arat;
            if iprt == 1 && iprint == 2 {
                c.io.so(&format!("\n   alpha={}", f(al, 10, 5)));
            }
            if iprt == 1 && iprint == 2 {
                c.io.so("\n     beta       s(alpha,beta)       ss(alpha,beta)     ss(alpha,-beta)");
            }
            for k in 1..=nbeta {
                let be = self.beta[k] * sc;
                let ss = self.ssm[(k, nal, itemp)];
                let s1 = ss * exp(-be / 2.0);
                let s2 = ss * exp(-be);
                let mut jprt = (k - 1) % self.nbint + 1;
                if k == nbeta {
                    jprt = 1;
                }
                if iprt == 1 && jprt == 1 && iprint == 2 {
                    c.io.so(&format!("{}{}{}{}", f(self.beta[k], 10, 4), e1p(s1, 18, 5), e1p(s2, 20, 5), e1p(ss, 20, 5)));
                }
            }
            if iprt == 1 {
                let col = self.col(nal, itemp);
                let v = &self.ssm.v[col];
                let beta = &self.beta;
                let (sum0, mut sum1) = moments(|k| beta[k], |k| v[(k - 1) as usize], nbeta);
                sum1 = sum1 / al;
                if iprint == 2 {
                    c.io.so(&format!("     normalization check ={}", f(sum0, 8, 4)));
                    c.io.so(&format!("          sum rule check ={}", f(sum1, 8, 4)));
                } else if iprint == 1 {
                    c.io.so(&format!(" {}{}{}", f(al, 10, 4), f(sum0, 10, 4), f(sum1, 10, 4)));
                }
            }
        }
    }

    // ================================================================ endout
    /// ENDF output.
    fn endout(
        &mut self,
        c: &mut Ctx,
        ntempr: i32,
        bragg: &RV,
        nedge: i32,
        scr: &mut RV,
        mscr: i32,
        isym: i32,
        ilog: i32,
    ) {
        let small = 1.0e-9f64;
        let tiny = -999.0f64;
        let tol = 0.9e-7f64;
        let up = 1.01f64;
        let (nalpha, nbeta, nout) = (self.nalpha, self.nbeta, self.nout);
        let (awr, aws, nss, b7) = (self.awr, self.aws, self.nss, self.b7);
        let (mut nb, mut nw) = (0i32, 0i32);

        //--write header
        let time = c.timer();
        c.io.so(&format!("\n\n endf output{}{}s", " ".repeat(57), f(time, 8, 1)));
        if isym == 0 {
            c.io.so("\n endf output is s");
        } else if isym == 1 {
            c.io.so("\n endf output is s for +&- beta");
        } else if isym == 2 {
            c.io.so("\n endf output is ss for - beta");
        } else if isym == 3 {
            c.io.so("\n endf output is ss for +&- beta");
        }
        if ilog == 1 {
            c.io.so(" log10 s or ss is written");
        }

        //---compute bound scattering cross sections
        let sb = self.spr * powc((1.0 + awr) / awr, 2);
        let mut sbs = 0.0f64;
        if aws != ZERO {
            sbs = self.sps * powc((1.0 + aws) / aws, 2);
        }

        //--for mixed moderators, merge ssm results
        if nss != 0 && b7 <= ZERO {
            let srat = sbs / sb;
            let nscr = -10;
            c.repoz(nscr);
            for k in 1..=ntempr {
                for j in 1..=nalpha {
                    let vals = c.io.read_rec_with(-nscr, |r| {
                        let mut v = Vec::with_capacity(nbeta as usize);
                        for _ in 0..nbeta {
                            v.push(r.r());
                        }
                        v
                    });
                    for k2 in 1..=nbeta {
                        scr[k2] = vals[(k2 - 1) as usize];
                    }
                    for k2 in 1..=nbeta {
                        self.ssm[(k2, j, k)] = srat * self.ssm[(k2, j, k)] + scr[k2];
                    }
                }
            }
        }

        //---display endf t-effective and debye-waller integral
        for k in 1..=ntempr {
            if nss == 0 || b7 > ZERO {
                self.dwpix[k] = self.dwpix[k] / (awr * self.tempr[k] * BK);
            } else {
                self.dwpix[k] = self.dwpix[k] / (aws * self.tempr[k] * BK);
                self.dwp1[k] = self.dwp1[k] / (awr * self.tempr[k] * BK);
            }
        }
        if nss == 0 || b7 > ZERO {
            c.io.so("\n    i     tempr       tempf       awr*w");
            for k in 1..=ntempr {
                c.io.so(&format!(
                    "{}{}{}{}",
                    i(k, 5),
                    e1p(self.tempr[k], 12, 3),
                    e1p(self.tempf[k], 12, 3),
                    e1p(awr * self.dwpix[k], 12, 3)
                ));
            }
        } else {
            c.io.so(&format!(
                "\n{}principal{}secondary\n    i     tempr       tempf       awr*w       tempf       aws*w",
                " ".repeat(26),
                " ".repeat(15)
            ));
            for k in 1..=ntempr {
                c.io.so(&format!(
                    "{}{}{}{}{}{}",
                    i(k, 5),
                    e1p(self.tempr[k], 12, 3),
                    e1p(self.tempf1[k], 12, 3),
                    e1p(awr * self.dwp1[k], 12, 3),
                    e1p(self.tempf[k], 12, 3),
                    e1p(aws * self.dwpix[k], 12, 3)
                ));
            }
        }

        //--write endf-6 file 1.
        if self.iel == 0 && self.twt == ZERO {
            self.iel = -1;
        }
        let iel = self.iel;
        c.io.so("\n\n");
        let nprnt = 0;
        c.nsp = 0;
        c.nsc = 0;
        c.math = 1;
        c.mfh = 0;
        c.mth = 0;
        c.nsh = 0;
        let blank = [b' '; 68];
        crate::farr::pack_text(&blank, &mut self.hz);
        let mut z = RV(self.hz.to_vec());
        c.tpidio(0, nout, nprnt, &mut z.all(), &mut nb, &mut nw);
        c.math = self.mat;
        c.mfh = 1;
        c.mth = 451;
        scr[1] = self.za;
        scr[2] = awr;
        scr[3] = -1.0;
        scr[4] = 0.0;
        scr[5] = 0.0;
        scr[6] = 0.0;
        c.contio(0, nout, nprnt, &mut scr.at(1), &mut nb, &mut nw);
        scr[1] = 0.0;
        scr[2] = 0.0;
        scr[3] = 0.0;
        scr[4] = 0.0;
        scr[5] = 0.0;
        scr[6] = 6.0;
        c.contio(0, nout, nprnt, &mut scr.at(1), &mut nb, &mut nw);
        scr[1] = 1.0;
        scr[2] = sigfig(THERM * self.beta[nbeta], 7, 0);
        scr[3] = 0.0;
        scr[4] = 0.0;
        scr[5] = 12.0;
        scr[6] = 6.0;
        c.contio(0, nout, nprnt, &mut scr.at(1), &mut nb, &mut nw);
        scr[1] = 0.0;
        scr[2] = 0.0;
        scr[3] = 0.0;
        scr[4] = 0.0;
        scr[5] = 0.0;
        scr[6] = 2.0;
        if iel != 0 {
            scr[6] = 3.0;
        }
        let mut n = 6i32;
        let mut nc = 0i32;
        let mut idone = 0;
        while idone == 0 {
            let mut text = "$".to_string();
            let v = c.io.input.read(&[Kind::Char(66)]);
            if let Some(Val::Char(s)) = &v[0] {
                text = s.clone();
            }
            if text.trim_end() == "$" {
                idone = 1;
            } else {
                nc = nc + 1;
                let mut t68: Vec<u8> = text.as_bytes().to_vec();
                t68.resize(66, b' ');
                t68.extend_from_slice(b"  ");
                crate::farr::pack_text(&t68, &mut self.hz);
                for j in 1..=17 {
                    scr[n + j] = self.hz[(j - 1) as usize];
                }
                n = n + 17;
                if n > mscr {
                    c.error("endout", "scratch storage exceeded for hollerith data", " ");
                }
            }
        }
        scr[5] = (17 * nc) as f64;
        let mut l = 1;
        c.hdatio(0, nout, nprnt, &mut scr.at(l), &mut nb, &mut nw);
        while nb != 0 {
            l = l + nw;
            c.moreio(0, nout, nprnt, &mut scr.at(l), &mut nb, &mut nw);
        }
        scr[1] = 0.0;
        scr[2] = 0.0;
        scr[3] = 1.0;
        scr[4] = 451.0;
        scr[5] = (5 + nc + 1) as f64;
        if iel != 0 {
            scr[5] = scr[5] + 1.0;
        }
        scr[6] = 0.0;
        let mut ii = 6;
        let mut ncards = 0i32;
        if iel != 0 {
            scr[ii + 1] = 0.0;
            scr[ii + 2] = 0.0;
            scr[ii + 3] = 7.0;
            scr[ii + 4] = 2.0;
            if iel < 0 {
                ncards = 3 + (2 * ntempr + 4) / 6;
            }
            if iel > 0 {
                ncards = 3 + (2 * nedge + 4) / 6;
            }
            if iel > 0 && ntempr > 1 {
                ncards = ncards + (ntempr - 1) * (1 + (nedge + 5) / 6);
            }
            scr[ii + 5] = ncards as f64;
            scr[ii + 6] = 0.0;
            ii = ii + 6;
        }
        scr[ii + 1] = 0.0;
        scr[ii + 2] = 0.0;
        scr[ii + 3] = 7.0;
        scr[ii + 4] = 4.0;
        ncards = 2 + (2 * nalpha + 4) / 6;
        if ntempr > 1 {
            ncards = ncards + (ntempr - 1) * (1 + (nalpha + 5) / 6);
        }
        ncards = 5 + nbeta * ncards;
        scr[ii + 5] = ncards as f64;
        scr[ii + 6] = 0.0;
        nw = 2;
        if iel != 0 {
            nw = 3;
        }
        c.dictio(0, nout, nprnt, &mut scr.at(1), &mut nb, &mut nw);
        c.asend(nout, nprnt);
        c.afend(nout, nprnt);

        //--write incoherent elastic part
        if iel < 0 {
            c.mfh = 7;
            c.mth = 2;
            scr[1] = self.za;
            scr[2] = awr;
            scr[3] = 2.0;
            scr[4] = 0.0;
            scr[5] = 0.0;
            scr[6] = 0.0;
            c.contio(0, nout, nprnt, &mut scr.at(1), &mut nb, &mut nw);
            scr[1] = sb * self.npr as f64;
            scr[2] = 0.0;
            scr[3] = 0.0;
            scr[4] = 0.0;
            scr[5] = 1.0;
            let mut ndw = ntempr;
            if ndw == 1 {
                ndw = 2;
            }
            scr[6] = ndw as f64;
            scr[7] = ndw as f64;
            scr[8] = 2.0;
            for k in 1..=ndw {
                if k <= ntempr {
                    scr[2 * k + 7] = self.tempr[k];
                    scr[2 * k + 8] = sigfig(self.dwpix[k], 7, 0);
                } else {
                    scr[2 * k + 7] = scr[2 * k + 5];
                    scr[2 * k + 8] = scr[2 * k + 6];
                }
            }
            c.tab1io(0, nout, nprnt, &mut scr.at(1), &mut nb, &mut nw);
            c.asend(nout, nprnt);
        }

        //--write coherent elastic part
        if iel >= 1 {
            c.mfh = 7;
            c.mth = 2;
            scr[1] = self.za;
            scr[2] = awr;
            scr[3] = 1.0;
            scr[4] = 0.0;
            scr[5] = 0.0;
            scr[6] = 0.0;
            c.contio(0, nout, nprnt, &mut scr.at(1), &mut nb, &mut nw);

            //--thin out the 1/e part at high energies.
            let mut sum = 0.0f64;
            let mut suml = 0.0f64;
            let mut w = self.dwpix[1];
            if nss > 0 && b7 == ZERO {
                w = (self.dwpix[1] + self.dwp1[1]) / 2.0;
            }
            // jmax is not set if no edge adds anything (Fortran: undefined)
            let mut jmax = 0i32;
            for j in 1..=nedge {
                let e = bragg[1 + 2 * j - 2];
                sum = sum + exp(-4.0 * w * e) * bragg[1 + 2 * j - 1];
                if sum - suml > tol * sum {
                    jmax = j;
                    suml = sum;
                }
            }

            //--now output the records for each temperature
            let npage = c.npage;
            let mut lo = vec![0.0f64; 4 * npage as usize + 8];
            fn put(scr: &mut RV, lo: &mut [f64], k: i32, v: f64) {
                if k >= 1 {
                    scr[k] = v;
                } else {
                    lo[(-k) as usize] = v;
                }
            }
            fn get(scr: &RV, lo: &[f64], k: i32) -> f64 {
                if k >= 1 {
                    scr[k]
                } else {
                    lo[(-k) as usize]
                }
            }
            for k in 1..=ntempr {
                if k == 1 {
                    scr[1] = self.tempr[k];
                    scr[2] = 0.0;
                    scr[3] = (ntempr - 1) as f64;
                    scr[4] = 0.0;
                    scr[5] = 1.0;
                    scr[6] = jmax as f64;
                    scr[7] = jmax as f64;
                    scr[8] = 1.0;
                    let mut ii = 8;
                    let mut jj = 0;
                    let mut j1 = 0;
                    sum = 0.0;
                    let mut w = self.dwpix[k];
                    if nss > 0 && b7 == ZERO {
                        w = (self.dwpix[k] + self.dwp1[k]) / 2.0;
                    }
                    let mut j = 0;
                    while j < nedge {
                        j = j + 1;
                        let e = bragg[1 + 2 * j - 2];
                        if j <= jmax {
                            jj = jj + 2;
                        }
                        // when 2*jmax is a multiple of npage, the Fortran
                        // writes scr(-1) and scr(0) here (before the array);
                        // those two words are kept in `lo`
                        put(scr, &mut lo, ii + jj - 1 - j1, sigfig(e, 7, 0));
                        put(scr, &mut lo, ii + jj - j1, sum + exp(-4.0 * w * e) * bragg[1 + 2 * j - 1]);
                        sum = get(scr, &lo, ii + jj - j1);
                        put(scr, &mut lo, ii + jj - j1, sigfig(sum, 7, 0));
                        if j >= nedge || jj - j1 >= npage {
                            if ii >= 8 {
                                c.tab1io(0, nout, nprnt, &mut scr.at(1), &mut nb, &mut nw);
                                ii = 0;
                                j1 = npage;
                            } else {
                                c.moreio(0, nout, nprnt, &mut scr.at(1), &mut nb, &mut nw);
                                j1 = j1 + npage;
                            }
                        }
                    }
                } else {
                    scr[1] = self.tempr[k];
                    scr[2] = 0.0;
                    scr[3] = 2.0;
                    scr[4] = 0.0;
                    scr[5] = jmax as f64;
                    scr[6] = 0.0;
                    let mut ii = 6;
                    let mut jj = 0;
                    let mut j1 = 0;
                    sum = 0.0;
                    let mut w = self.dwpix[k];
                    if nss > 0 && b7 == ZERO {
                        w = (self.dwpix[k] + self.dwp1[k]) / 2.0;
                    }
                    for j in 1..=nedge {
                        if j <= jmax {
                            jj = jj + 1;
                        }
                        let e = sigfig(bragg[1 + 2 * jj - 2], 7, 0);
                        put(scr, &mut lo, ii + jj - j1, sum + exp(-4.0 * w * e) * bragg[1 + 2 * jj - 1]);
                        sum = get(scr, &lo, ii + jj - j1);
                        put(scr, &mut lo, ii + jj - j1, sigfig(sum, 7, 0));
                        if j >= nedge || jj - j1 >= npage {
                            if ii >= 6 {
                                c.listio(0, nout, nprnt, &mut scr.at(1), &mut nb, &mut nw);
                                ii = 0;
                                j1 = npage;
                            } else {
                                c.moreio(0, nout, nprnt, &mut scr.at(1), &mut nb, &mut nw);
                                j1 = j1 + npage;
                            }
                        }
                    }
                }
            }
            c.asend(nout, nprnt);
        }

        //--write inelastic part
        c.mfh = 7;
        c.mth = 4;
        scr[1] = self.za;
        scr[2] = awr;
        scr[3] = 0.0;
        scr[4] = self.lat as f64;
        scr[5] = isym as f64;
        scr[6] = 0.0;
        c.contio(0, nout, nprnt, &mut scr.at(1), &mut nb, &mut nw);
        scr[1] = 0.0;
        scr[2] = 0.0;
        scr[3] = 0.0;
        if ilog != 0 {
            scr[3] = 1.0;
        }
        scr[4] = 0.0;
        scr[5] = 6.0;
        if nss > 0 {
            scr[5] = (6 * (nss + 1)) as f64;
        }
        scr[6] = nss as f64;
        scr[7] = self.npr as f64 * self.spr;
        scr[8] = self.beta[nbeta];
        scr[9] = awr;
        scr[10] = sigfig(THERM * self.beta[nbeta], 7, 0);
        scr[11] = 0.0;
        scr[12] = self.npr as f64;
        if nss != 0 {
            scr[13] = b7;
            scr[14] = self.mss as f64 * self.sps;
            scr[15] = aws;
            scr[16] = 0.0;
            scr[17] = 0.0;
            scr[18] = self.mss as f64;
        }
        c.listio(0, nout, nprnt, &mut scr.at(1), &mut nb, &mut nw);
        scr[1] = 0.0;
        scr[2] = 0.0;
        scr[3] = 0.0;
        scr[4] = 0.0;
        scr[5] = 1.0;
        let mut nbt = nbeta;
        if isym == 1 || isym == 3 {
            nbt = 2 * nbeta - 1;
        }
        scr[6] = nbt as f64;
        scr[7] = nbt as f64;
        scr[8] = 4.0;
        c.tab2io(0, nout, nprnt, &mut scr.at(1), &mut nb, &mut nw);
        let smin = self.smin;
        for ib in 1..=nbt {
            for nt in 1..=ntempr {
                let mut sc = 1.0f64;
                if self.lat == 1 {
                    sc = THERM / (BK * self.tempr[nt]);
                }
                // value for alpha index j: off is the word offset in scr
                // (7+2*j for the TAB1, 6+j for the LIST)
                let tab = nt == 1;
                scr[1] = self.tempr[nt];
                if isym % 2 == 0 {
                    scr[2] = self.beta[ib];
                }
                if isym % 2 == 1 && ib < nbeta {
                    scr[2] = -self.beta[nbeta - ib + 1];
                }
                if isym % 2 == 1 && ib >= nbeta {
                    scr[2] = self.beta[ib - nbeta + 1];
                }
                let be = scr[2] * sc;
                if tab {
                    scr[3] = (ntempr - 1) as f64;
                    scr[4] = 0.0;
                    scr[5] = 1.0;
                    scr[6] = nalpha as f64;
                    scr[7] = nalpha as f64;
                    scr[8] = 4.0;
                } else {
                    scr[3] = 4.0;
                    scr[4] = 0.0;
                    scr[5] = nalpha as f64;
                    scr[6] = 0.0;
                }
                for j in 1..=nalpha {
                    let idx = if tab { 8 + 2 * j } else { 6 + j };
                    if tab {
                        scr[7 + 2 * j] = self.alpha[j];
                    }
                    // the LIST records use 0 (not tiny) for log of 0 when
                    // isym=0
                    let tlog = if tab || isym != 0 { tiny } else { 0.0 };
                    let mut val = scr[idx];
                    if isym == 0 {
                        let s = self.ssm[(ib, j, nt)];
                        if ilog == 0 {
                            val = s * exp(-be / 2.0);
                            val = sig76(val, small);
                        } else {
                            val = tlog;
                            if s > ZERO {
                                val = log(s) - be / 2.0;
                                val = sigfig(val, 7, 0);
                            }
                        }
                    } else if isym == 1 {
                        let s = if ib < nbeta { self.ssm[(nbeta - ib + 1, j, nt)] } else { self.ssp[(ib - nbeta + 1, j, nt)] };
                        if ilog == 0 {
                            val = s * exp(be / 2.0);
                            val = sig76(val, small);
                        } else {
                            val = tlog;
                            if s > ZERO {
                                val = log(s) + be / 2.0;
                                val = sigfig(val, 7, 0);
                            }
                        }
                    } else if isym == 2 {
                        let s = self.ssm[(ib, j, nt)];
                        if ilog == 0 {
                            val = s;
                            val = sig76(val, small);
                        } else {
                            val = tlog;
                            if s > ZERO {
                                val = log(s);
                                val = sigfig(val, 7, 0);
                            }
                        }
                    } else if isym == 3 {
                        if ib < nbeta {
                            let s = self.ssm[(nbeta - ib + 1, j, nt)];
                            if ilog == 0 {
                                val = s;
                                val = sig76(val, small);
                            } else {
                                val = tlog;
                                if s > ZERO {
                                    val = log(s);
                                    val = sigfig(val, 7, 0);
                                }
                            }
                        } else {
                            let s = self.ssp[(ib - nbeta + 1, j, nt)];
                            if ilog == 0 {
                                val = s;
                                val = sig76(val, small);
                            } else {
                                val = tlog;
                                // the Fortran tests ssp(-nbeta+1,j,nt) here
                                if self.ssp_bad(j, nt) > ZERO {
                                    val = log(s);
                                    val = sigfig(val, 7, 0);
                                }
                            }
                        }
                    }
                    scr[idx] = val;
                    if ilog == 0 && scr[idx] < smin {
                        scr[idx] = 0.0;
                    }
                }
                if tab {
                    c.tab1io(0, nout, nprnt, &mut scr.at(1), &mut nb, &mut nw);
                    let mut ll = 1 + nw;
                    while nb != 0 {
                        c.moreio(0, nout, nprnt, &mut scr.at(ll), &mut nb, &mut nw);
                        ll = ll + nw;
                    }
                } else {
                    c.listio(0, nout, nprnt, &mut scr.at(1), &mut nb, &mut nw);
                    let mut ll = 1;
                    while nb != 0 {
                        ll = ll + nw;
                        c.moreio(0, nout, nprnt, &mut scr.at(ll), &mut nb, &mut nw);
                    }
                }
            }
        }
        if nss != 0 && b7 <= ZERO {
            scr[1] = 0.0;
            scr[2] = 0.0;
            scr[3] = 0.0;
            scr[4] = 0.0;
            scr[5] = 1.0;
            let ntf = ntempr;
            scr[6] = ntf as f64;
            scr[7] = ntf as f64;
            scr[8] = 2.0;
            for k in 1..=ntf {
                if k <= ntempr {
                    scr[2 * k + 7] = sigfig(self.tempr[k], 7, 0);
                    scr[2 * k + 8] = sigfig(self.tempf1[k], 7, 0);
                } else {
                    scr[2 * k + 7] = sigfig(up * scr[2 * k + 5], 7, 0);
                    scr[2 * k + 8] = scr[2 * k + 6];
                }
            }
            c.tab1io(0, nout, nprnt, &mut scr.at(1), &mut nb, &mut nw);
        }
        scr[1] = 0.0;
        scr[2] = 0.0;
        scr[3] = 0.0;
        scr[4] = 0.0;
        scr[5] = 1.0;
        let ntf = ntempr;
        scr[6] = ntf as f64;
        scr[7] = ntf as f64;
        scr[8] = 2.0;
        for k in 1..=ntf {
            if k <= ntempr {
                scr[2 * k + 7] = sigfig(self.tempr[k], 7, 0);
                scr[2 * k + 8] = sigfig(self.tempf[k], 7, 0);
            } else {
                scr[2 * k + 7] = sigfig(up * scr[2 * k + 5], 7, 0);
                scr[2 * k + 8] = scr[2 * k + 6];
            }
        }
        c.tab1io(0, nout, nprnt, &mut scr.at(1), &mut nb, &mut nw);
        c.asend(nout, nprnt);
        c.afend(nout, nprnt);
        c.amend(nout, nprnt);
        c.atend(nout, nprnt);
    }

    /// `ssp(-nbeta+1,j,nt)` as the Fortran reads it: the word nbeta places
    /// before ssp(1,j,nt) in column-major order. For j=1, nt=1 this is
    /// before the start of the array; the value there is taken as positive.
    fn ssp_bad(&self, j: i32, nt: i32) -> f64 {
        let nb = self.ssp.n1 as i64;
        let flat = -nb + nb * (j as i64 - 1) + nb * self.ssp.n2 as i64 * (nt as i64 - 1);
        if flat < 0 {
            1.0
        } else {
            self.ssp.v[flat as usize]
        }
    }
}

/// `sigfig(x,7,0)` if x >= small, else `sigfig(x,6,0)`.
fn sig76(x: f64, small: f64) -> f64 {
    if x >= small {
        sigfig(x, 7, 0)
    } else {
        sigfig(x, 6, 0)
    }
}

/// Integrals over the phonon frequency of 2*p*beta**n*hyperbolic.
fn fsum(n: i32, p: &RV, np: i32, tau: f64, deltab: f64) -> f64 {
    let arg = deltab * tau / 2.0;
    let edsq = exp(arg);
    let mut v = 1.0f64;
    let an = (1 - 2 * (n % 2)) as f64;
    let mut be = 0.0f64;
    let mut fs = 0.0f64;
    let mut w = 1.0f64;
    for ij in 1..=np {
        if n > 0 {
            w = powi(be, n);
        }
        let mut ff = ((p[ij] * v) * v + (p[ij] * an / v) / v) * w;
        if ij == 1 || ij == np {
            ff = ff / 2.0;
        }
        fs = fs + ff;
        be = be + deltab;
        v = v * edsq;
    }
    fs * deltab
}

/// Interpolate in a table of t_n(beta).
fn terpt(tn: &RV, ntn: i32, delta: f64, be: f64) -> f64 {
    if be > ntn as f64 * delta {
        return 0.0;
    }
    let mut k = int(be / delta);
    if k < ntn - 1 {
        let bt = k as f64 * delta;
        let btp = bt + delta;
        k = k + 1;
        tn[k] + (be - bt) * (tn[k + 1] - tn[k]) / (btp - bt)
    } else {
        0.0
    }
}

/// Next term of the phonon expansion: convolve t1 with tlast into tnext.
/// `etab[m]` is exp(-m*delta). Returns the integral check ckk.
fn convol(t1: &RV, tlast: &RV, tnext: &mut RV, n1: i32, nl: i32, nn: i32, delta: f64, etab: &[f64]) -> f64 {
    let tiny = 1.0e-30f64;
    let mut ckk = 0.0f64;
    for k in 1..=nn {
        let mut tk = 0.0f64;
        for j in 1..=n1 {
            let t1j = t1[j];
            if t1j > ZERO {
                let i1 = k + j - 2;
                let i2 = k - j;
                let mut f1 = 0.0;
                if i1 + 1 <= nl {
                    f1 = tlast[i1 + 1] * etab[(j - 1) as usize];
                }
                let mut f2 = 0.0;
                if i2 >= 0 && i2 + 1 <= nl {
                    f2 = tlast[i2 + 1];
                } else if i2 < 0 && 1 - i2 <= nl {
                    f2 = tlast[1 - i2] * etab[(-i2) as usize];
                }
                let mut cc = t1j * (f1 + f2);
                if j == 1 || j == n1 {
                    cc = cc / 2.0;
                }
                tk = tk + cc;
            }
        }
        tk = tk * delta;
        if tk < tiny {
            tk = 0.0;
        }
        tnext[k] = tk;
        let mut cc = tk;
        let be = (k - 1) as f64 * delta;
        // etab[m] is exp(-(m*delta)): the same operations
        let eb = match etab.get((k - 1) as usize) {
            Some(&v) => v,
            None => exp(-be),
        };
        cc = cc + tk * eb;
        if k == 1 || k == nn {
            cc = cc / 2.0;
        }
        ckk = ckk + cc;
    }
    ckk * delta
}

/// Interpolate (log-linear) in a table of S(alpha,beta).
fn terps(sd: &RV, nsd: i32, delta: f64, be: f64) -> f64 {
    let slim = -225.0f64;
    if be > delta * nsd as f64 {
        return 0.0;
    }
    let mut k = int(be / delta);
    if k < nsd - 1 {
        let bt = k as f64 * delta;
        let btp = bt + delta;
        k = k + 1;
        let st = if sd[k] <= ZERO { slim } else { log(sd[k]) };
        let stp = if sd[k + 1] <= ZERO { slim } else { log(sd[k + 1]) };
        let stt = st + (be - bt) * (stp - st) / (btp - bt);
        let mut t = 0.0;
        if stt > slim {
            t = exp(stt);
        }
        return t;
    }
    0.0
}

/// The loop of sbfill for one beta: the points k = 1, 2, ... of the
/// table, with delta increased where the Fortran increases it. The value of
/// point k is computed, and given to set, only for lo < k <= hi; the search
/// for the beta range runs for every point (as in the Fortran) when hi > 0,
/// except for the points up to lo when mono is set: with betan strictly
/// increasing the search ends at the same range from any start.
/// Returns the number of points.
#[allow(clippy::too_many_arguments)]
#[inline]
fn sbfill_core(
    nbt: i32, delta: &mut f64, be: f64, ls: &[f64], betan: &RV, nbeta: i32, lo: i32, hi: i32, mono: bool,
    mut set: impl FnMut(i32, f64),
) -> i32 {
    let shade = 1.00001f64;
    let slim = -225.0f64;
    let bmin = -be - (nbt - 1) as f64 * *delta;
    let bmax = -be + (nbt - 1) as f64 * *delta + *delta / 100.0;
    let mut j = nbeta;
    let mut k = 0;
    let mut bet = bmin;
    while bet <= bmax {
        k = k + 1;
        if hi > 0 && !(mono && k <= lo) {
            let b = bet.abs();
            // search for correct beta range
            let mut idone = 0;
            while idone == 0 {
                if b > betan[j] {
                    if j == nbeta && b < shade * betan[j] {
                        idone = 1;
                    } else {
                        if j == nbeta {
                            idone = 2;
                        } else {
                            // move up
                            j = j + 1;
                        }
                    }
                } else {
                    if b > betan[j - 1] {
                        idone = 1;
                    } else {
                        if j == 2 {
                            idone = 1;
                        } else {
                            // move down
                            j = j - 1;
                        }
                    }
                }
            }
            if k > lo && k <= hi {
                // interpolate in this range
                if idone == 1 {
                    // ls(j) is log(s(j)), or slim if s(j) <= 0
                    let st = ls[j as usize];
                    let stm = ls[(j - 1) as usize];
                    let mut v = st + (b - betan[j]) * (stm - st) / (betan[j - 1] - betan[j]);
                    if bet > ZERO {
                        v = v - bet;
                    }
                    let arg = v;
                    v = 0.0;
                    if arg > slim {
                        v = exp(arg);
                    }
                    set(k, v);
                } else {
                    set(k, 0.0);
                }
            }
        }
        // if delta is too small for the current value of beta, increase it
        while bet == bet + *delta {
            *delta = *delta * 10.0;
        }
        bet = bet + *delta;
    }
    k
}

/// s(beta) on a new energy grid for convolution with a diffusion or
/// free-gas shape. delta can be increased (as in the Fortran).
/// Returns the number of points it set (sb(1..k)).
#[allow(clippy::too_many_arguments)]
fn sbfill(
    c: &mut Ctx,
    sb: &mut RV,
    nbt: i32,
    delta: &mut f64,
    be: f64,
    ls: &[f64],
    betan: &RV,
    nbeta: i32,
    _nbe: i32,
    ndmax: i32,
) -> i32 {
    let bmin = -be - (nbt - 1) as f64 * *delta;
    let bmax = -be + (nbt - 1) as f64 * *delta + *delta / 100.0;
    if 1 + int((bmax - bmin) / *delta) > ndmax {
        let strng = format!("ndmax needs to be at least {}", i(1 + int((bmax - bmin) / *delta), 6));
        c.error("sbfill", &strng, " ");
    }
    sbfill_core(nbt, delta, be, ls, betan, nbeta, 0, i32::MAX, false, |k, v| sb[k] = v)
}

/// The delta that sbfill leaves for (nbt, delta, be) and its number of
/// points, from the same loop without the table (no exp). None where
/// sbfill stops (ndmax too small).
fn sbfill_delta(nbt: i32, delta: f64, be: f64, ndmax: i32) -> Option<(f64, i32)> {
    let mut delta = delta;
    let bmin = -be - (nbt - 1) as f64 * delta;
    let bmax = -be + (nbt - 1) as f64 * delta + delta / 100.0;
    if 1 + int((bmax - bmin) / delta) > ndmax {
        return None;
    }
    let dummy = RV::new(0);
    let k = sbfill_core(nbt, &mut delta, be, &[], &dummy, 0, 0, 0, false, |_, _| {});
    Some((delta, k))
}

/// betan(1..nbeta) is strictly increasing (no NaN).
fn ascending_strict(betan: &RV, nbeta: i32) -> bool {
    let n = (nbeta.max(0) as usize).min(betan.0.len());
    betan.0[..n].windows(2).all(|w| w[0] < w[1])
}

/// sb(lo+1..hi) as the serial loop has it just before beta i of an alpha:
/// each point m has the value of the latest earlier beta of the alpha that
/// set it (info(ip) = (delta at the start of beta ip, number of points it
/// set)), else the value in base (sb at the start of the alpha).
#[allow(clippy::too_many_arguments)]
fn sb_patch(
    sb: &mut RV, lo: i32, hi: i32, i: i32, info: &[(f64, i32)], base: &RV, nbt: i32, ls: &[f64], betan: &RV,
    nbeta: i32, mono: bool,
) {
    let mut covered = lo;
    let mut ip = i;
    while covered < hi && ip > 1 {
        ip -= 1;
        let (d0, ks) = info[ip as usize];
        if ks > covered {
            let top = ks.min(hi);
            let mut d = d0;
            sbfill_core(nbt, &mut d, betan[ip], ls, betan, nbeta, covered, top, mono, |k, v| sb[k] = v);
            covered = top;
        }
    }
    for m in covered + 1..=hi {
        sb[m] = base[m];
    }
}

/// Set sb to what the serial beta loop of an alpha leaves in it: each point
/// from the latest beta that set it (the others are not changed).
fn sb_final(sb: &mut RV, info: &[(f64, i32)], nbt: i32, ls: &[f64], betan: &RV, nbeta: i32) {
    let mono = ascending_strict(betan, nbeta);
    let mut covered = 0;
    for ip in (1..info.len()).rev() {
        let (d0, ks) = info[ip];
        if ks > covered {
            let mut d = d0;
            sbfill_core(nbt, &mut d, betan[ip as i32], ls, betan, nbeta, covered, ks, mono, |k, v| sb[k] = v);
            covered = ks;
        }
    }
}

/// One beta of the convolution loop of `trans`: sbfill, the convolution
/// sum and the delta-function part; delta is left as the Fortran leaves
/// it. With pt = (info, base) the points that the convolution reads but
/// this sbfill does not set are set first with `sb_patch`, so that the
/// result does not depend on what sb held before.
#[allow(clippy::too_many_arguments)]
fn trans_beta(
    c: &mut Ctx, sb: &mut RV, etab: &mut Vec<f64>, etab_delta: &mut f64, sd: &RV, nbt: i32, delta: &mut f64, be: f64,
    lap: &[f64], betan: &RV, nbeta: i32, ibeta: i32, ndmax: i32, al: f64, f0: f64, pt: Option<(&[(f64, i32)], &RV)>,
) -> f64 {
    let tiny = 1.0e-30f64;
    let mut s = 0.0f64;

    //--prepare table of continuous ss on new interval
    let kset = sbfill(c, sb, nbt, delta, be, lap, betan, nbeta, ibeta, ndmax);
    if let Some((info, base)) = pt {
        let need = 2 * nbt - 1;
        if kset < need {
            let mono = ascending_strict(betan, nbeta);
            sb_patch(sb, kset, need, ibeta, info, base, nbt, lap, betan, nbeta, mono);
        }
    }
    let delta = *delta;

    //--convolve s-transport with s-continuous
    if etab_delta.to_bits() != delta.to_bits() || etab.len() < nbt as usize {
        *etab_delta = delta;
        *etab = (0..nbt.max(1) as usize).map(|m| exp(-(m as f64 * delta))).collect();
    }
    for k in 1..=nbt {
        let mut ff = (2 * ((k - 1) % 2 + 1)) as f64;
        if k == 1 || k == nbt {
            ff = 1.0;
        }
        s = s + ff * sd[k] * sb[nbt + k - 1];
        s = s + ff * sd[k] * sb[nbt - k + 1] * etab[(k - 1) as usize];
    }
    s = s * delta / 3.0;
    if s < tiny {
        s = 0.0;
    }
    let st = terps(sd, nbt, delta, be);
    if st > ZERO {
        s = s + exp(-al * f0) * st;
    }
    s
}

/// `trans_beta` for every beta of an alpha, split between threads, with
/// pt = (info, base): result(ibeta) = (s, delta at the start, delta at the
/// end), info(ibeta) = (delta at the start, number of sbfill points). The
/// delta at the start of each beta and the number of points come first
/// from `sbfill_delta`, so the betas do not depend on each other. base is
/// sb at the start of the alpha. A result is None if its work stopped or
/// tried to write (the caller computes it again). Both are empty with one
/// thread or little work: the caller then runs the serial loop.
#[allow(clippy::too_many_arguments, clippy::type_complexity)]
fn trans_pre(
    sd: &RV, nbt: i32, delta: f64, lap: &[f64], betan: &RV, nbeta: i32, ndmax: i32, al: f64, f0: f64, base: &RV,
) -> (Vec<Option<(f64, f64, f64)>>, Vec<(f64, i32)>) {
    let nth = crate::util::threads();
    if nth <= 1 || nbeta < 2 || (nbeta as i64) * (nbt as i64) < 200_000 {
        return (Vec::new(), Vec::new());
    }
    let nb = nbeta as usize;
    // the delta at the start of each beta and the points that sbfill sets,
    // as in the serial loop; (NaN, 0) after a beta where sbfill stops
    // (in blocks: each block is found in parallel with the delta at its
    // start; where delta changes, the next block starts after that beta)
    let mut info = vec![(f64::NAN, 0i32); nb + 1];
    let mut d = delta;
    let mut first = 1usize;
    let block = 64 * nth;
    'chain: while first <= nb {
        let last = (first + block - 1).min(nb);
        let mut r: Vec<Option<(f64, i32)>> = vec![None; last - first + 1];
        let per = r.len().div_ceil(nth);
        std::thread::scope(|sc| {
            for (ic, part) in r.chunks_mut(per).enumerate() {
                sc.spawn(move || {
                    for (j, slot) in part.iter_mut().enumerate() {
                        let ibeta = first + ic * per + j;
                        *slot = sbfill_delta(nbt, d, betan[ibeta as i32], ndmax);
                    }
                });
            }
        });
        for ibeta in first..=last {
            match r[ibeta - first] {
                Some((dn, k)) => {
                    info[ibeta] = (d, k);
                    if dn.to_bits() != d.to_bits() {
                        d = dn;
                        first = ibeta + 1;
                        continue 'chain;
                    }
                }
                None => {
                    info[ibeta] = (d, 0);
                    break 'chain;
                }
            }
        }
        first = last + 1;
    }
    let info_r = &info;
    let mut out: Vec<Option<(f64, f64, f64)>> = vec![None; nb + 1];
    let chunk = nb.div_ceil(nth);
    std::thread::scope(|sc| {
        for (ic, part) in out[1..].chunks_mut(chunk).enumerate() {
            sc.spawn(move || {
                let mut wc = Ctx::new(crate::fio::Io::detached());
                let mut sb = RV::new(base.0.len());
                let mut etab = Vec::new();
                let mut etab_delta = f64::NAN;
                for (j, slot) in part.iter_mut().enumerate() {
                    let ibeta = (ic * chunk + j + 1) as i32;
                    let be = betan[ibeta];
                    let d0 = info_r[ibeta as usize].0;
                    if d0.is_nan() {
                        continue;
                    }
                    let mut d = d0;
                    wc.io.touched = false;
                    let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                        trans_beta(
                            &mut wc, &mut sb, &mut etab, &mut etab_delta, sd, nbt, &mut d, be, lap, betan, nbeta,
                            ibeta, ndmax, al, f0, Some((info_r.as_slice(), base)),
                        )
                    }));
                    *slot = match r {
                        Ok(s) if !wc.io.touched => Some((s, d0, d)),
                        _ => None,
                    };
                }
            });
        }
    });
    (out, info)
}

/// Modified Bessel function K1 (the exponential part for x>1 is omitted).
fn besk1(x: f64) -> f64 {
    let c0 = 0.125f64;
    let c1 = 0.442850424f64;
    let c2 = 0.584115288f64;
    let c3 = 6.070134559f64;
    let c4 = 17.864913364f64;
    let c5 = 48.858995315f64;
    let c6 = 90.924600045f64;
    let c7 = 113.795967431f64;
    let c8 = 85.331474517f64;
    let c9 = 32.00008698f64;
    let c10 = 3.999998802f64;
    let c11 = 1.304923514f64;
    let c12 = 1.47785657f64;
    let c13 = 16.402802501f64;
    let c14 = 44.732901977f64;
    let c15 = 115.837493464f64;
    let c16 = 198.437197312f64;
    let c17 = 222.869709703f64;
    let c18 = 142.216613971f64;
    let c19 = 40.000262262f64;
    let c20 = 1.999996391f64;
    let c21 = 1.0f64;
    let c22 = 0.5f64;
    let c23 = 0.5772156649f64;
    let c24 = 1.0f64;
    let c25 = 0.0108241775f64;
    let c26 = 0.0788000118f64;
    let c27 = 0.2581303765f64;
    let c28 = 0.5050238576f64;
    let c29 = 0.663229543f64;
    let c30 = 0.6283380681f64;
    let c31 = 0.4594342117f64;
    let c32 = 0.2847618149f64;
    let c33 = 0.1736431637f64;
    let c34 = 0.1280426636f64;
    let c35 = 0.1468582957f64;
    let c36 = 0.4699927013f64;
    let c37 = 1.2533141373f64;
    let test = 1.0f64;
    if x <= test {
        let v = c0 * x;
        let u = v * v;
        let bi1 = (((((((((c1 * u + c2) * u + c3) * u + c4) * u + c5) * u + c6) * u + c7) * u + c8) * u + c9) * u
            + c10)
            * v;
        let bi3 = ((((((((c11 * u + c12) * u + c13) * u + c14) * u + c15) * u + c16) * u + c17) * u + c18) * u
            + c19)
            * u
            + c20;
        c21 / x + bi1 * (log(c22 * x) + c23) - v * bi3
    } else {
        let u = c24 / x;
        let mut bi3 = -c25 * u + c26;
        bi3 = bi3 * u - c27;
        bi3 = bi3 * u + c28;
        bi3 = bi3 * u - c29;
        bi3 = bi3 * u + c30;
        bi3 = bi3 * u - c31;
        bi3 = bi3 * u + c32;
        bi3 = bi3 * u - c33;
        bi3 = bi3 * u + c34;
        bi3 = bi3 * u - c35;
        bi3 = bi3 * u + c36;
        bi3 = bi3 * u + c37;
        sqrt(u) * bi3
    }
}

/// Bessel function terms for discrete oscillators.
fn bfact(x: f64, bzero: &mut f64, bplus: &mut [f64; 51], bminus: &mut [f64; 51], dwc: f64, betai: f64) {
    let c0 = 3.75f64;
    let c1 = 1.0f64;
    let c2 = 3.5156229f64;
    let c3 = 3.0899424f64;
    let c4 = 1.2067492f64;
    let c5 = 0.2659732f64;
    let c6 = 0.0360768f64;
    let c7 = 0.0045813f64;
    let c8 = 0.39894228f64;
    let c9 = 0.01328592f64;
    let c10 = 0.00225319f64;
    let c11 = 0.00157565f64;
    let c12 = 0.00916281f64;
    let c13 = 0.02057706f64;
    let c14 = 0.02635537f64;
    let c15 = 0.01647633f64;
    let c16 = 0.00392377f64;
    let c17 = 0.5f64;
    let c18 = 0.87890594f64;
    let c19 = 0.51498869f64;
    let c20 = 0.15084934f64;
    let c21 = 0.02658733f64;
    let c22 = 0.00301532f64;
    let c23 = 0.00032411f64;
    let c24 = 0.02282967f64;
    let c25 = 0.02895312f64;
    let c26 = 0.01787654f64;
    let c27 = 0.00420059f64;
    let c28 = 0.39894228f64;
    let c29 = 0.03988024f64;
    let c30 = 0.00362018f64;
    let c31 = 0.00163801f64;
    let c32 = 0.01031555f64;
    let big = 1.0e10f64;
    let tiny = 1.0e-30f64;
    let one = 1.0f64;
    let mut bn = [0.0f64; 51];

    //--compute bessi0
    let y = x / c0;
    let bessi0;
    if y <= one {
        let u = y * y;
        bessi0 = c1 + u * (c2 + u * (c3 + u * (c4 + u * (c5 + u * (c6 + u * c7)))));
    } else {
        let v = 1.0 / y;
        bessi0 = (c8 + v * (c9 + v * (c10 + v * (-c11 + v * (c12 + v * (-c13 + v * (c14 + v * (-c15 + v * c16))))))))
            / sqrt(x);
    }

    //--compute bessi1
    let mut bessi1;
    if y <= one {
        let u = y * y;
        bessi1 = (c17 + u * (c18 + u * (c19 + u * (c20 + u * (c21 + u * (c22 + u * c23)))))) * x;
    } else {
        let v = 1.0 / y;
        bessi1 = c24 + v * (-c25 + v * (c26 - v * c27));
        bessi1 = c28 + v * (-c29 + v * (-c30 + v * (c31 + v * (-c32 + v * bessi1))));
        bessi1 = bessi1 / sqrt(x);
    }

    //--generate higher orders by reverse recursion
    let imax = 50usize;
    bn[imax] = 0.0;
    bn[imax - 1] = 1.0;
    let mut k = imax - 1;
    while k > 1 {
        bn[k - 1] = bn[k + 1] + k as f64 * (2.0 / x) * bn[k];
        k = k - 1;
        if bn[k] >= big {
            for j in k..=imax {
                bn[j] = bn[j] / big;
            }
        }
    }
    let rat = bessi1 / bn[1];
    for k in 1..=imax {
        bn[k] = bn[k] * rat;
        if bn[k] < tiny {
            bn[k] = 0.0;
        }
    }

    //--apply exponential terms to bessel functions
    if y <= one {
        *bzero = bessi0 * exp(-dwc);
        for k in 1..=imax {
            bminus[k] = 0.0;
            bplus[k] = 0.0;
            if bn[k] != ZERO {
                let arg = -dwc - k as f64 * betai / 2.0;
                bplus[k] = exp(arg) * bn[k];
                if bplus[k] < tiny {
                    bplus[k] = 0.0;
                }
                let arg = -dwc + k as f64 * betai / 2.0;
                bminus[k] = exp(arg) * bn[k];
                if bminus[k] < tiny {
                    bminus[k] = 0.0;
                }
            }
        }
    } else {
        *bzero = bessi0 * exp(-dwc + x);
        for k in 1..=imax {
            bminus[k] = 0.0;
            bplus[k] = 0.0;
            if bn[k] != ZERO {
                let arg = -dwc - k as f64 * betai / 2.0 + x;
                bplus[k] = exp(arg) * bn[k];
                if bplus[k] < tiny {
                    bplus[k] = 0.0;
                }
                let arg = -dwc + k as f64 * betai / 2.0 + x;
                bminus[k] = exp(arg) * bn[k];
                if bminus[k] < tiny {
                    bminus[k] = 0.0;
                }
            }
        }
    }
}

/// Set up the arrays bex and rdbex used by sint. Returns nbx.
fn bfill(bex: &mut RV, rdbex: &mut RV, betan: &RV, nbeta: i32) -> i32 {
    let small = 1.0e-9f64;
    let mut k = nbeta;
    for m in 1..=nbeta {
        bex[m] = -betan[k];
        k = k - 1;
    }
    if betan[1] <= small {
        bex[nbeta] = 0.0;
        k = nbeta + 1;
    } else {
        k = nbeta + 2;
        bex[nbeta + 1] = betan[1];
    }
    for m in 2..=nbeta {
        bex[k] = betan[m];
        k = k + 1;
    }
    let nbxm = k - 2;
    let nbx = k - 1;
    for m in 1..=nbxm {
        rdbex[m] = 1.0 / (bex[m + 1] - bex[m]);
    }
    nbx
}

/// Set up the array sex for sint from one column of ssm (0-based slice).
fn exts(sexpb: &[f64], sex: &mut RV, exb: &RV, betan: &RV, nbeta: i32) {
    let small = 1.0e-9f64;
    let sp = |m: i32| sexpb[(m - 1) as usize];
    let mut k = nbeta;
    for m in 1..=nbeta {
        sex[m] = sp(k);
        k = k - 1;
    }
    if betan[1] <= small {
        sex[nbeta] = sp(1);
        k = nbeta + 1;
    } else {
        k = nbeta + 2;
        sex[nbeta + 1] = sp(1);
    }
    for m in 2..=nbeta {
        sex[k] = sp(m) * exb[m] * exb[m];
        k = k + 1;
    }
}

/// Interpolate in the scattering function, or use the SCT approximation
/// outside the range in memory.
fn sint(
    x: f64,
    bex: &RV,
    rdbex: &RV,
    sex: &RV,
    nbx: i32,
    alph: f64,
    wt: f64,
    tbart: f64,
    betan: &RV,
    nbeta: i32,
) -> f64 {
    let slim = -225.0f64;

    //--sct approx
    if x.abs() > betan[nbeta] {
        let sv;
        if alph <= ZERO {
            sv = 0.0;
        } else {
            let mut ex = -powc(wt * alph - x.abs(), 2) / (4.0 * wt * alph * tbart);
            if x > ZERO {
                ex = ex - x;
            }
            sv = exp(ex) / (4.0 * PI * wt * alph * tbart);
        }
        return sv;
    }

    //--interpolation
    let mut k1 = 1;
    let mut k2 = nbeta;
    let mut k3 = nbx;
    //  bisect for x
    let mut idone = 0;
    while idone == 0 {
        if x == bex[k2] {
            return sex[k2];
        } else if x > bex[k2] {
            k1 = k2;
            k2 = (k3 - k2) / 2 + k2;
            if k3 - k1 <= 1 {
                idone = 1;
            }
        } else {
            k3 = k2;
            k2 = (k2 - k1) / 2 + k1;
            if k3 - k1 <= 1 {
                idone = 1;
            }
        }
    }
    let ss1 = if sex[k1] <= ZERO { slim } else { log(sex[k1]) };
    let ss3 = if sex[k3] <= ZERO { slim } else { log(sex[k3]) };
    let ex = ((bex[k3] - x) * ss1 + (x - bex[k1]) * ss3) * rdbex[k1];
    let mut sv = 0.0;
    if ex > slim {
        sv = exp(ex);
    }
    sv
}

/// Statistical weight factor for cold hydrogen or deuterium.
fn bt(j: i32, x: f64) -> f64 {
    let half = 0.5f64;
    let mut yy = half * j as f64 * (j + 1) as f64;
    let a = (2 * j + 1) as f64 * exp(-yy * x);
    let mut b = 0.0f64;
    for m in 1..=10 {
        let mut k = 2 * m - 2;
        if j % 2 == 1 {
            k = k + 1;
        }
        yy = half * k as f64 * (k + 1) as f64;
        b = b + (2 * k + 1) as f64 * exp(-yy * x);
    }
    a / (2.0 * b)
}

/// Sum over Bessel functions and Clebsch-Gordon coefficients.
fn sumh(c: &mut Ctx, j: i32, jp: i32, y: f64) -> f64 {
    if j == 0 {
        powc(sjbes(c, jp, y) * cn(j, jp, jp), 2)
    } else if jp == 0 {
        powc(sjbes(c, j, y) * cn(j, 0, j), 2)
    } else {
        let mut sum1 = 0.0f64;
        let imk = (j - jp).abs() + 1;
        let ipk1 = j + jp + 1;
        let mpk = ipk1 - imk;
        let ipk = if mpk <= 9 { ipk1 } else { imk + 9 };
        for n in imk..=ipk {
            let n1 = n - 1;
            sum1 = sum1 + powc(sjbes(c, n1, y) * cn(j, jp, n1), 2);
        }
        sum1
    }
}

/// log(k!) by the sum of logs, then exp, as in cn.
fn lfact(k: i32) -> f64 {
    let mut s = 0.0f64;
    for m in 1..=k {
        let zi = m as f64;
        s = s + log(zi);
    }
    s
}

/// Clebsch-Gordon coefficients for cold hydrogen or deuterium.
fn cn(jj: i32, ll: i32, nn: i32) -> f64 {
    let kdet = (jj + ll + nn) / 2;
    let kdel = jj + ll + nn - 2 * kdet;
    if kdel == 0 {
        let ka1 = jj + ll + nn;
        let ka2 = jj + ll - nn;
        let ka3 = jj - ll + nn;
        let ka4 = ll - jj + nn;
        let kb1 = ka1 / 2;
        let kb2 = ka2 / 2;
        let kb3 = ka3 / 2;
        let kb4 = ka4 / 2;
        let fa = |k: i32| {
            let s = lfact(k);
            let mut fact = 1.0f64;
            if s > ZERO {
                fact = exp(s);
            }
            sqrt(fact)
        };
        let fb = |k: i32| {
            let s = lfact(k);
            let mut b = 1.0f64;
            if s > ZERO {
                b = exp(s);
            }
            b
        };
        let a1 = fa(ka1);
        let a2 = fa(ka2);
        let a3 = fa(ka3);
        let a4 = fa(ka4);
        let b1 = fb(kb1);
        let b2 = fb(kb2);
        let b3 = fb(kb3);
        let b4 = fb(kb4);
        let mut rat = (2 * nn + 1) as f64;
        rat = rat / (jj + ll + nn + 1) as f64;
        let iwign = (jj + ll - nn) / 2;
        let mut wign = ipow(-1, iwign) as f64;
        wign = wign * sqrt(rat) * b1 / a1 * a2 / b2 * a3 / b3 * a4 / b4;
        wign
    } else {
        0.0
    }
}

/// Spherical Bessel functions for cold hydrogen or deuterium.
fn sjbes(c: &mut Ctx, n: i32, x: f64) -> f64 {
    let huge = 1.0e25f64;
    let small = 2.0e-38f64;
    let break1 = 3.0e4f64;
    let break2 = 7.0e-4f64;
    let break3 = 0.2f64;
    let one = 1.0f64;
    let ten = 10.0f64;
    let hund = 100.0f64;

    //--check for large arguments
    if n >= 30000 || x > break1 {
        let strng = format!("value is not accurate  n = {}{}x = {}", i(n, 7), " ".repeat(10), e(x, 14, 7));
        c.mess("sjbes", &strng, " ");
        return 0.0;
    }

    //--check for bad arguments
    if x < ZERO || n < 0 {
        let strng = format!("argument is invalid  n = {}{}x = {}", i(n, 7), " ".repeat(10), e(x, 14, 7));
        c.error("sjbes", &strng, " ");
    }

    //--compute normal values
    let bessel;
    if x <= break2 {
        let w = 1.0f64;
        if n == 0 {
            bessel = w;
        } else if n > 10 {
            bessel = 0.0;
        } else {
            let mut t1 = 3.0f64;
            let mut t2 = 1.0f64;
            let mut t3 = 0.0f64;
            for _ in 1..=n {
                t3 = t2 * x / t1;
                t1 = t1 + 2.0;
                t2 = t3;
            }
            bessel = t3;
        }
    } else {
        let w;
        if x < break3 {
            let y = powc(x, 2);
            w = 1.0 - y * (1.0 - y / 20.0) / 6.0;
        } else {
            w = sin(x) / x;
        }
        if n == 0 {
            bessel = w;
        } else {
            let l;
            if x >= hund {
                l = int(x / 50.0 + 18.0);
            } else if x >= ten {
                l = int(x / 10.0 + 10.0);
            } else if x > one {
                l = int(x / 2.0 + 5.0);
            } else {
                l = 5;
            }
            let iii = int(x);
            let mut kmax = n;
            if iii > n {
                kmax = iii;
            }
            let nm = kmax + l;
            let z = 1.0 / x;
            let mut t3 = 0.0f64;
            let mut t2 = small;
            let mut t1 = 0.0f64;
            let mut sj = 0.0f64;
            for m in 1..=nm {
                let k = nm - m;
                t1 = (2 * k + 3) as f64 * z * t2 - t3;
                if n == k {
                    sj = t1;
                }
                if t1.abs() >= huge {
                    t1 = t1 / huge;
                    t2 = t2 / huge;
                    sj = sj / huge;
                }
                t3 = t2;
                t2 = t1;
            }
            bessel = w * sj / t1;
        }
    }
    bessel
}

/// Interpolate in a table of ska(kappa).
fn terpk(ska: &RV, nka: i32, delta: f64, be: f64) -> f64 {
    if be > nka as f64 * delta {
        return 1.0;
    }
    let mut k = int(be / delta);
    if k < nka - 1 {
        let bt = k as f64 * delta;
        let btp = bt + delta;
        k = k + 1;
        return ska[k] + (be - bt) * (ska[k + 1] - ska[k]) / (btp - bt);
    }
    1.0
}

/// Bragg energies and structure factors for coherent elastic scattering.
fn coher(c: &mut Ctx, lat: i32, natom: i32, b: &mut RV, nbe: &mut i32, maxb: &mut i32, emax: f64) {
    let gr1 = 2.4573e-8f64;
    let gr2 = 6.700e-8f64;
    let gr3 = 12.011f64;
    let gr4 = 5.50f64;
    let be1 = 2.2856e-8f64;
    let be2 = 3.5832e-8f64;
    let be3 = 9.01f64;
    let be4 = 7.53f64;
    let beo1 = 2.695e-8f64;
    let beo2 = 4.39e-8f64;
    let beo3 = 12.5f64;
    let beo4 = 1.0f64;
    let al1 = 4.04e-8f64;
    let al3 = 26.7495f64;
    let al4 = 1.495f64;
    let pb1 = 4.94e-8f64;
    let pb3 = 207.0f64;
    let pb4 = 1.0f64;
    let fe1 = 2.86e-8f64;
    let fe3 = 55.454f64;
    let fe4 = 12.9f64;
    let twothd = 0.666666666667f64;
    let sqrt3 = 1.732050808f64;
    let toler = 1.0e-6f64;
    let eps = 0.05f64;

    //--write header
    let time = c.timer();
    c.io.so(&format!(
        "\n bragg edges for coherent elastic scattering{}{}s",
        " ".repeat(25),
        f(time, 8, 1)
    ));

    //--initialize.
    let twopis = powc(2.0 * PI, 2);
    let amne = AMASSN * AMU;
    let econ = EV * 8.0 * (amne / HBAR) / HBAR;
    let recon = 1.0 / econ;
    let tsqx = econ / 20.0;
    let (a, cc, amsc, scoh);
    if lat == 1 {
        // graphite constants.
        a = gr1;
        cc = gr2;
        amsc = gr3;
        scoh = gr4 / natom as f64;
    } else if lat == 2 {
        //  beryllium constants
        a = be1;
        cc = be2;
        amsc = be3;
        scoh = be4 / natom as f64;
    } else if lat == 3 {
        //  beryllium oxide constants
        a = beo1;
        cc = beo2;
        amsc = beo3;
        scoh = beo4 / natom as f64;
    } else if lat == 4 {
        // aluminum constants
        a = al1;
        cc = 0.0;
        amsc = al3;
        scoh = al4 / natom as f64;
    } else if lat == 5 {
        // lead constants
        a = pb1;
        cc = 0.0;
        amsc = pb3;
        scoh = pb4 / natom as f64;
    } else if lat == 6 {
        // iron constants
        a = fe1;
        cc = 0.0;
        amsc = fe3;
        scoh = fe4 / natom as f64;
    } else {
        c.error("coh", "illegal lat.", " ");
    }
    let (mut c1, mut c2, mut scon) = (0.0f64, 0.0f64, 0.0f64);
    if lat < 4 {
        c1 = 4.0 / (3.0 * a * a);
        c2 = 1.0 / (cc * cc);
        scon = scoh * powc(4.0 * PI, 2) / (2.0 * a * a * cc * sqrt3 * econ);
    } else if (4..=5).contains(&lat) {
        c1 = 3.0 / (a * a);
        scon = scoh * powc(4.0 * PI, 2) / (16.0 * a * a * a * econ);
    } else if lat == 6 {
        c1 = 2.0 / (a * a);
        scon = scoh * powc(4.0 * PI, 2) / (8.0 * a * a * a * econ);
    }
    let wint = 0.0f64;
    let t2 = HBAR / (2.0 * AMU * amsc);
    let ulim = econ * emax;
    let ifl = 1i32;
    let nw = *maxb;
    let mut k: i32;

    let tausq = |m1: i32, m2: i32, m3: i32| -> f64 {
        (c1 * (m1 * m1 + m2 * m2 + m1 * m2) as f64 + ((m3 * m3) as f64 * c2)) * twopis
    };
    let taufcc = |m1: i32, m2: i32, m3: i32| -> f64 {
        c1 * ((m1 * m1 + m2 * m2 + m3 * m3) as f64 + twothd * m1 as f64 * m2 as f64 + twothd * m1 as f64 * m3 as f64
            - twothd * m2 as f64 * m3 as f64)
            * twopis
    };
    let taubcc = |m1: i32, m2: i32, m3: i32| -> f64 {
        c1 * (m1 * m1 + m2 * m2 + m3 * m3 + m1 * m2 + m2 * m3 + m1 * m3) as f64 * twopis
    };

    if lat <= 3 {
        //--compute lattice factors for hexagonal lattices
        let phi = ulim / twopis;
        let mut i1m = int(a * sqrt(phi));
        i1m = i1m + 1;
        k = 0;
        for i1 in 1..=i1m {
            let l1 = i1 - 1;
            let mut i2m = int((l1 as f64 + sqrt(3.0 * (a * a * phi - (l1 * l1) as f64))) / 2.0);
            i2m = i2m + 1;
            for i2 in i1..=i2m {
                let l2 = i2 - 1;
                let x = phi - c1 * (l1 * l1 + l2 * l2 - l1 * l2) as f64;
                let mut i3m = 0;
                if x > ZERO {
                    i3m = int(cc * sqrt(x));
                }
                i3m = i3m + 1;
                for i3 in 1..=i3m {
                    let l3 = i3 - 1;
                    let mut w1 = 2.0f64;
                    if l1 == l2 {
                        w1 = 1.0;
                    }
                    let mut w2 = 2.0f64;
                    if l1 == 0 || l2 == 0 {
                        w2 = 1.0;
                    }
                    if l1 == 0 && l2 == 0 {
                        w2 = 1.0;
                    }
                    if l1 == 0 && l2 == 0 {
                        w2 = w2 / 2.0;
                    }
                    let mut w3 = 2.0f64;
                    if l3 == 0 {
                        w3 = 1.0;
                    }
                    for sgn in [1, -1] {
                        let tsq = tausq(l1, sgn * l2, l3);
                        if tsq > ZERO && tsq <= ulim {
                            let tau = sqrt(tsq);
                            let w = exp(-tsq * t2 * wint) * w1 * w2 * w3 / tau;
                            let ff = w * formf(lat, l1, sgn * l2, l3);
                            if k <= 0 || tsq <= tsqx {
                                k = k + 1;
                                if 2 * k > nw {
                                    c.error("coh", "storage exceeded", " ");
                                }
                                b[ifl + 2 * k - 2] = tsq;
                                b[ifl + 2 * k - 1] = ff;
                            } else {
                                let mut m = 0;
                                let mut idone = 0;
                                while m < k && idone == 0 {
                                    m = m + 1;
                                    if tsq >= b[ifl + 2 * m - 2] && tsq < (1.0 + eps) * b[ifl + 2 * m - 2] {
                                        b[ifl + 2 * m - 1] = b[ifl + 2 * m - 1] + ff;
                                        idone = 1;
                                    }
                                }
                                if idone == 0 {
                                    k = k + 1;
                                    if 2 * k > nw {
                                        c.error("coh", "storage exceeded", " ");
                                    }
                                    b[ifl + 2 * k - 2] = tsq;
                                    b[ifl + 2 * k - 1] = ff;
                                }
                            }
                        }
                    }
                }
            }
        }
    } else {
        //--compute lattice factors for fcc (lat 4,5) or bcc (lat 6) lattices
        let i1m = 15;
        k = 0;
        for i1 in -i1m..=i1m {
            for i2 in -i1m..=i1m {
                for i3 in -i1m..=i1m {
                    let tsq = if lat <= 5 { taufcc(i1, i2, i3) } else { taubcc(i1, i2, i3) };
                    if tsq > ZERO && tsq <= ulim {
                        let tau = sqrt(tsq);
                        let w = exp(-tsq * t2 * wint) / tau;
                        let ff = w * formf(lat, i1, i2, i3);
                        k = k + 1;
                        if 2 * k > nw {
                            c.error("coh", "storage exceeded", " ");
                        }
                        b[ifl + 2 * k - 2] = tsq;
                        b[ifl + 2 * k - 1] = ff;
                    }
                }
            }
        }
    }
    let imax = k - 1;

    //--sort lattice factors
    for m in 1..=imax {
        let jmin = m + 1;
        for j in jmin..=k {
            if b[ifl + 2 * j - 2] < b[ifl + 2 * m - 2] {
                let st = b[ifl + 2 * m - 2];
                let sf = b[ifl + 2 * m - 1];
                b[ifl + 2 * m - 2] = b[ifl + 2 * j - 2];
                b[ifl + 2 * m - 1] = b[ifl + 2 * j - 1];
                b[ifl + 2 * j - 2] = st;
                b[ifl + 2 * j - 1] = sf;
            }
        }
    }
    k = k + 1;
    b[ifl + 2 * k - 2] = ulim;
    b[ifl + 2 * k - 1] = b[ifl + 2 * k - 3];

    //--convert to practical units
    //--and combine duplicate bragg edges.
    let mut bel = -1.0f64;
    let mut j = 0;
    for m in 1..=k {
        let be = b[ifl + 2 * m - 2] * recon;
        let bs = b[ifl + 2 * m - 1] * scon;
        if be - bel < toler {
            b[ifl + 2 * j - 1] = b[ifl + 2 * j - 1] + bs;
        } else {
            j = j + 1;
            b[ifl + 2 * j - 2] = be;
            b[ifl + 2 * j - 1] = bs;
            bel = be;
        }
    }
    *nbe = j;
    *maxb = 2 * *nbe;
    c.io.so(&format!("\n   found{} edges below{} ev", i(*nbe, 5), f(emax, 6, 2)));
}

/// Form factors for the lattice.
fn formf(lat: i32, l1: i32, l2: i32, l3: i32) -> f64 {
    let c1 = 7.54f64;
    let c2 = 4.24f64;
    let c3 = 11.31f64;
    if lat == 1 {
        // graphite.
        let k = l3 / 2;
        if 2 * k != l3 {
            powc(sin(PI * (l1 - l2) as f64 / 3.0), 2)
        } else {
            (6.0 + 10.0 * cos(2.0 * PI * (l1 - l2) as f64 / 3.0)) / 4.0
        }
    } else if lat == 2 {
        // beryllium.
        1.0 + cos(2.0 * PI * (2 * l1 + 4 * l2 + 3 * l3) as f64 / 6.0)
    } else if lat == 3 {
        // beryllium oxide.
        (1.0 + cos(2.0 * PI * (2 * l1 + 4 * l2 + 3 * l3) as f64 / 6.0)) * (c1 + c2 + c3 * cos(3.0 * PI * l3 as f64 / 4.0))
    } else if lat == 4 || lat == 5 {
        // fcc lattices.
        let e1 = 2.0 * PI * l1 as f64;
        let e2 = 2.0 * PI * (l1 + l2) as f64;
        let e3 = 2.0 * PI * (l1 + l3) as f64;
        powc(1.0 + cos(e1) + cos(e2) + cos(e3), 2) + powc(sin(e1) + sin(e2) + sin(e3), 2)
    } else if lat == 6 {
        // bcc lattices.
        let e1 = 2.0 * PI * (l1 + l2 + l3) as f64;
        powc(1.0 + cos(e1), 2) + powc(sin(e1), 2)
    } else {
        // (formf is not set in the Fortran)
        0.0
    }
}
