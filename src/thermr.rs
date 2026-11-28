//! THERMR: thermal scattering cross sections and energy-angle
//! distributions, ported from NJOY2016 `thermr.f90` (module thermm).
//!
//! The routines keep the Fortran structure: module variables are fields of
//! `Thermr`, `save`d locals are fields too, and labels become loops or
//! `match` states where the control flow needs them. Arrays keep the
//! Fortran (1-based) indices.

#![allow(clippy::too_many_arguments)]
#![allow(clippy::collapsible_if)]
#![allow(clippy::collapsible_else_if)]
#![allow(clippy::needless_range_loop)]

use crate::endf::Ctx;
use crate::farr::RV;
use crate::ffmt;
use crate::fio::{Kind, Val};
use crate::fmath::{cos, exp, int, log, nint, sin, sqrt};
use crate::mathm::legndr;
use crate::util::physics::{AMASSN, AMU, BK, EV, HBAR, PI};
use crate::util::sigfig;

const NBUF: i32 = 1000;
const NWSCR: i32 = 500000;
/// Extra words past the Fortran array sizes, so that writes the Fortran
/// does past an array end (without checking) do not stop the port first.
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

/// THERMR module variables (and the `save`d locals of sigcoh and calcem).
pub struct Thermr {
    nendf: i32,
    nin: i32,
    nout: i32,
    nscr: i32,
    nscr2: i32,
    matde: i32,
    nbin: i32,
    iprint: i32,
    ncds: i32,
    matdp: i32,
    natom: i32,
    ntemp: i32,
    iinc: i32,
    iform: i32,
    ncdse: i32,
    za: f64,
    awr: f64,
    tol: f64,
    emax: f64,
    sb: f64,
    az: f64,
    tevz: f64,
    teff: f64,
    sb2: f64,
    az2: f64,
    teff2: f64,
    cliq: f64,
    lasym: i32,
    lat: i32,
    scr: RV,
    bufo: RV,
    bufn: RV,
    tempr: RV,
    fl: Option<RV>,
    esi: Option<RV>,
    xsi: Option<RV>,
    // sigcoh saved
    sc_k: i32,
    sc_recon: f64,
    sc_scon: f64,
    // calcem saved
    ce_nwtab: i32,
    ce_sabmin: f64,
    ce_nl: i32,
    ce_nlt: i32,
    ce_nlp: i32,
    ce_nlp1: i32,
    ce_nl1: i32,
    ce_nnl: i32,
    ce_jmax: i32,
    ce_nne: i32,
}

impl Default for Thermr {
    fn default() -> Self {
        Self::new()
    }
}

impl Thermr {
    pub fn new() -> Self {
        Thermr {
            nendf: 0,
            nin: 0,
            nout: 0,
            nscr: 0,
            nscr2: 0,
            matde: 0,
            nbin: 0,
            iprint: 0,
            ncds: 0,
            matdp: 0,
            natom: 0,
            ntemp: 0,
            iinc: 0,
            iform: 0,
            ncdse: 0,
            za: 0.0,
            awr: 0.0,
            tol: 0.0,
            emax: 0.0,
            sb: 0.0,
            az: 0.0,
            tevz: 0.0,
            teff: 0.0,
            sb2: 0.0,
            az2: 0.0,
            teff2: 0.0,
            cliq: 0.0,
            lasym: 0,
            lat: 0,
            scr: RV::new(0),
            bufo: RV::new(0),
            bufn: RV::new(0),
            tempr: RV::new(0),
            fl: None,
            esi: None,
            xsi: None,
            sc_k: 0,
            sc_recon: 0.0,
            sc_scon: 0.0,
            ce_nwtab: 0,
            ce_sabmin: 0.0,
            ce_nl: 0,
            ce_nlt: 0,
            ce_nlp: 0,
            ce_nlp1: 0,
            ce_nl1: 0,
            ce_nnl: 0,
            ce_jmax: 0,
            ce_nne: 0,
        }
    }

    /// `thermr`: the module driver.
    pub fn thermr(&mut self, c: &mut Ctx) {
        const S1099: f64 = 3.76e0;
        const A1099: f64 = 15.86e0;
        const S1095: f64 = 4.74e0;
        const A1095: f64 = 11.90e0;
        const EPS: f64 = 1.0e-4;
        const BREAK: f64 = 3000.0e0;
        const SMALL: f64 = 1.0e-10;
        const UP: f64 = 1.00001e0;
        const ZERO: f64 = 0.0;
        let (mut nb, mut nw) = (0i32, 0i32);
        let mut icoh = 0i32;
        let mut mtref = 0i32;
        let mut idis = 0i32;
        let mut ex = [0.0f64; 3];
        let (mut e, mut enext, mut s) = (0.0f64, 0.0f64, 0.0f64);

        //--start
        let time = c.timer();
        c.io.so(&format!(
            "\n thermr...compute thermal scattering cross sections and distributions{}s",
            ffmt::f(time, 8, 1)
        ));
        c.io.se(&format!("\n thermr...{}{}s", " ".repeat(59), ffmt::f(time, 8, 1)));
        let mut sz2 = 0.0f64;
        self.scr = RV::new(NWSCR as usize + SLACK);
        self.bufn = RV::new(NBUF as usize);
        self.bufo = RV::new(NBUF as usize);
        let mut indexc = -1i32;
        let mut indexi = -1i32;

        //--read user input
        let v = c.io.input.read(&[Kind::Int, Kind::Int, Kind::Int]);
        self.nendf = geti(&v[0], self.nendf);
        self.nin = geti(&v[1], self.nin);
        self.nout = geti(&v[2], self.nout);
        if self.nin == 0 {
            c.error("thermr", "nin=0.", " ");
        }
        if self.nin < 0 && self.nout > 0 {
            c.error("thermr", "mode conversion not allowed.", " ");
        }
        if self.nin > 0 && self.nout < 0 {
            c.error("thermr", "mode conversion not allowed.", " ");
        }
        c.openz(self.nendf, 0);
        c.openz(self.nin, 0);
        c.openz(self.nout, 1);
        self.iprint = 0;
        self.ntemp = 1;
        let v = c.io.input.read(&[Kind::Int; 10]);
        self.matde = geti(&v[0], self.matde);
        self.matdp = geti(&v[1], self.matdp);
        self.nbin = geti(&v[2], self.nbin);
        self.ntemp = geti(&v[3], self.ntemp);
        self.iinc = geti(&v[4], self.iinc);
        icoh = geti(&v[5], icoh);
        self.iform = geti(&v[6], self.iform);
        self.natom = geti(&v[7], self.natom);
        mtref = geti(&v[8], mtref);
        self.iprint = geti(&v[9], self.iprint);
        if !(221..=250).contains(&mtref) {
            c.error("thermr", "illegal reference mt.", " ");
        }
        let ntemp = self.ntemp;
        self.tempr = RV::new(ntemp.max(1) as usize);
        let mut eftemp = RV::new(ntemp.max(1) as usize);
        let mut eftmp2 = RV::new(ntemp.max(1) as usize);
        let kinds: Vec<Kind> = (0..ntemp).map(|_| Kind::Real).collect();
        let v = c.io.input.read(&kinds);
        for i in 1..=ntemp {
            self.tempr[i] = getr(&v[(i - 1) as usize], self.tempr[i]);
        }
        for i in 1..=ntemp {
            eftemp[i] = 0.0;
            eftmp2[i] = 0.0;
            if self.matde == 0 {
                eftemp[i] = self.tempr[i];
            }
        }
        let v = c.io.input.read(&[Kind::Real, Kind::Real]);
        self.tol = getr(&v[0], self.tol);
        self.emax = getr(&v[1], self.emax);

        //--check endf format
        c.tpidio(self.nendf, 0, 0, &mut self.scr.all(), &mut nb, &mut nw);
        c.contio(self.nendf, 0, 0, &mut self.scr.all(), &mut nb, &mut nw);
        c.contio(self.nendf, 0, 0, &mut self.scr.all(), &mut nb, &mut nw);
        if c.n1h != 0 {
            c.iverf = 4;
        } else if c.n2h == 0 {
            c.iverf = 5;
        } else {
            c.iverf = 6;
        }
        c.skiprz(self.nendf, -1);

        //--set up for sct approximation (pre endf-6 only)
        if c.iverf < 6 {
            sz2 = 0.0;
            self.az2 = 0.0;
            if self.matde == 1099 {
                sz2 = S1099;
            }
            if self.matde == 1099 {
                self.az2 = A1099;
            }
            if self.matde == 1095 {
                sz2 = S1095;
            }
            if self.matde == 1095 {
                self.az2 = A1095;
            }
            self.sb2 = sz2;
            if self.az2 != ZERO {
                let q = (self.az2 + 1.0) / self.az2;
                self.sb2 = sz2 * (q * q);
            }
            if self.matde != 0 {
                gateff(&self.tempr, &mut eftemp, ntemp, self.matde);
                if sz2 != ZERO {
                    gatef2(&self.tempr, &mut eftmp2, ntemp, self.matde);
                }
            }
        }

        //--write parameters.
        c.io.so(&format!(
            "\n unit for endf tape ................... {}\n unit for input pendf tape ............ {}\n unit for output pendf tape ........... {}",
            ffmt::i(self.nendf, 10),
            ffmt::i(self.nin, 10),
            ffmt::i(self.nout, 10)
        ));
        c.io.so(&format!(
            "\n material to be processed (endf) ...... {}\n material to be processed (pendf) ..... {}\n number of angle bins ................. {}\n number of temperatures ............... {}\n inelastic option ..................... {}\n elastic option ....................... {}\n MF6 format option .................... {}\n number of principal atoms ............ {}\n reference mt ......................... {}\n print option (0 min, 1 max) .......... {}",
            ffmt::i(self.matde, 10),
            ffmt::i(self.matdp, 10),
            ffmt::i(self.nbin, 10),
            ffmt::i(ntemp, 10),
            ffmt::i(self.iinc, 10),
            ffmt::i(icoh, 10),
            ffmt::i(self.iform, 10),
            ffmt::i(self.natom, 10),
            ffmt::i(mtref, 10),
            ffmt::i(self.iprint, 10)
        ));
        c.io.so(&format!(" temperatures (kelvin) ................ {}", ffmt::e1p(self.tempr[1], 10, 4)));
        for i in 2..=ntemp {
            c.io.so(&format!("{}{}", " ".repeat(40), ffmt::e1p(self.tempr[i], 10, 4)));
        }
        c.io.so(&format!(
            " tolerance ............................ {}\n max energy for thermal treatment ..... {}",
            ffmt::e1p(self.tol, 10, 4),
            ffmt::e1p(self.emax, 10, 4)
        ));
        if self.iinc == 2 && c.iverf < 6 {
            c.io.so("\n parameters for sct app.\n -----------------------");
            c.io.so(&format!(" effective temperatures ............... {}", ffmt::e1p(eftemp[1], 10, 4)));
            for i in 2..=ntemp {
                c.io.so(&format!("{}{}", " ".repeat(40), ffmt::e1p(eftemp[i], 10, 4)));
            }
            if sz2 != ZERO {
                c.io.so(&format!(
                    " free xsec for second atom ............ {}\n mass ratio for second atom ........... {}\n eff. temps for second atom ........... {}",
                    ffmt::e1p(sz2, 10, 4),
                    ffmt::e1p(self.az2, 10, 4),
                    ffmt::e1p(eftmp2[1], 10, 4)
                ));
                for i in 2..=ntemp {
                    c.io.so(&format!("{}{}", " ".repeat(40), ffmt::e1p(eftmp2[i], 10, 4)));
                }
            }
        }
        c.io.so(&format!("\n endf uses endf-{} format", ffmt::i(c.iverf, 1)));

        if self.nendf == 0 && self.iinc == 2 {
            c.error("thermr", "No ENDF tape given for iinc=2 (nendf is 0)", " ");
        }

        //--initialize i/o units
        let mut iold = 10i32;
        let mut inew = 11i32;
        self.nscr = 12;
        if self.nout < 0 {
            self.nscr = -self.nscr;
        }
        self.nscr2 = 13;
        if self.nout < 0 {
            self.nscr2 = -self.nscr2;
        }
        c.openz(-iold, 1);
        c.openz(-inew, 1);
        c.openz(self.nscr, 1);
        c.openz(self.nscr2, 1);
        c.repoz(self.nendf);
        c.repoz(self.nout);
        c.repoz(self.nin);
        c.repoz(self.nscr2);

        //--copy through to desired material
        c.nsh = 0;
        c.tpidio(self.nin, self.nout, 0, &mut self.scr.all(), &mut nb, &mut nw);
        loop {
            c.contio(self.nin, 0, 0, &mut self.scr.all(), &mut nb, &mut nw);
            if c.math == self.matdp {
                break;
            }
            if c.math < 0 {
                c.error("thermr", "desired material not on pendf tape.", " ");
            }
            c.contio(0, self.nout, 0, &mut self.scr.all(), &mut nb, &mut nw);
            c.tomend(self.nin, self.nout, 0, &mut self.scr.all());
        }
        c.contio(self.nin, 0, 0, &mut self.scr.at(7), &mut nb, &mut nw);
        let iverp = if c.n1h != 0 {
            4
        } else if c.n2h == 0 {
            5
        } else {
            6
        };
        c.skiprz(self.nin, -1);
        c.io.so(&format!("\n pendf uses endf-{} format", ffmt::i(iverp, 1)));

        //--loop over desired temperatures
        let mut itemp = 1i32;
        let mut icopy = 0i32;
        let mut t;
        loop {
            // label 130
            c.contio(0, 0, self.nscr2, &mut self.scr.all(), &mut nb, &mut nw);
            if c.math != self.matdp {
                c.error("thermr", "desired material not on pendf tape.", " ");
            }
            let mut nwb = 0i32;
            self.za = self.scr[1];
            self.awr = self.scr[2];
            self.az = self.scr[2];
            if iverp >= 5 {
                c.contio(self.nin, 0, self.nscr2, &mut self.scr.all(), &mut nb, &mut nw);
            }
            if iverp >= 6 {
                c.contio(self.nin, 0, self.nscr2, &mut self.scr.all(), &mut nb, &mut nw);
            }
            c.hdatio(self.nin, 0, self.nscr2, &mut self.scr.all(), &mut nb, &mut nw);
            t = self.scr[1];
            if (t - self.tempr[itemp]).abs() > EPS * self.tempr[itemp] {
                if t > self.tempr[itemp] {
                    c.error("thermr", "desired temperature not on tape.", " ");
                }
                //--count temperatures to be copied
                if itemp == 1 {
                    icopy += 1;
                } else {
                    let templ = self.tempr[itemp - 1];
                    let mut it = 0;
                    loop {
                        it += 1;
                        if it >= itemp {
                            break;
                        }
                        if (t - self.tempr[it]).abs() <= EPS * self.tempr[it] {
                            continue;
                        }
                        if (t - templ).abs() <= EPS * templ {
                            continue;
                        }
                        if t < templ {
                            continue;
                        }
                        icopy += 1;
                    }
                }
                // label 160
                if icopy == 0 {
                    c.repoz(self.nscr2);
                }
                let ntape = if icopy == 0 { 0 } else { self.nscr2 };
                c.tomend(self.nin, 0, ntape, &mut self.scr.all());
                c.contio(self.nin, 0, 0, &mut self.scr.all(), &mut nb, &mut nw);
                continue;
            }

            //--copy elastic cross section to scratch file
            // label 180
            c.findf(self.matdp, 3, 2, self.nin);
            c.contio(self.nin, 0, 0, &mut self.scr.all(), &mut nb, &mut nw);
            let mut nex = 2i32;
            let mut ne = 0i32;
            e = 0.0;
            c.gety1(e, &mut enext, &mut idis, &mut s, self.nin, &mut self.scr.all());
            let emaxs = self.emax;
            if t > BREAK {
                self.emax = self.emax * t / BREAK;
            }
            loop {
                e = enext;
                c.gety1(e, &mut enext, &mut idis, &mut s, self.nin, &mut self.scr.all());
                ex[1] = e;
                ex[2] = s;
                if e > self.emax * (1.0 + SMALL) {
                    ex[2] = 0.0;
                }
                ne += 1;
                if e > self.emax * (1.0 + SMALL) {
                    ne = -ne;
                }
                c.loada(ne, &ex[1..], nex, inew, &mut self.bufn.0, NBUF);
                if ne < 0 {
                    break;
                }
                if enext > self.emax * (1.0 + SMALL) {
                    enext = self.emax;
                }
                if (e - self.emax).abs() < SMALL * self.emax {
                    enext = UP * self.emax;
                }
            }
            let mut np = -ne;
            self.emax = emaxs;
            std::mem::swap(&mut iold, &mut inew);
            c.skiprz(self.nin, -1);

            //--read the elastic data (endf-6 format)
            let mut lthr = 0;
            if c.iverf >= 6 && self.nendf != 0 {
                c.findf(self.matde, 7, 0, self.nendf);
                c.contio(self.nendf, 0, 0, &mut self.scr.all(), &mut nb, &mut nw);
                if c.mth == 2 {
                    lthr = c.l1h;
                    icoh = 10 * lthr;
                    let temp = self.tempr[itemp];
                    self.rdelas(c, temp, lthr, self.nendf, &mut nwb, &mut indexc, &mut indexi);
                }
            }

            //--compute incoherent inelastic cross sections
            if self.iinc != 0 {
                let mut temp = self.tempr[itemp];
                self.teff = eftemp[itemp];
                self.teff2 = 0.0;
                if sz2 > ZERO {
                    self.teff2 = eftmp2[itemp];
                }
                if self.iinc == 2 {
                    temp = t;
                }
                self.calcem(c, temp, itemp, &mut iold, &mut inew, &mut np, &mut nex, mtref);
            }

            //--compute elastic cross sections
            if icoh > 0 && icoh <= 10 {
                self.coh(c, icoh, itemp, &mut iold, &mut inew, &mut np, &mut nex, mtref + 1, indexc);
            } else if icoh > 10 && icoh <= 20 {
                self.iel(c, icoh, itemp, &mut iold, &mut inew, &mut np, &mut nex, mtref + 1, indexi);
            } else if icoh > 20 {
                self.coh(c, 10, itemp, &mut iold, &mut inew, &mut np, &mut nex, mtref + 1, indexc);
                self.iel(c, 20, itemp, &mut iold, &mut inew, &mut np, &mut nex, mtref + 2, indexi);
            }

            //--write results on output tape
            self.tpend(c, iold, itemp, np, nex, icoh, &mut icopy, mtref);

            //--continue temperature loop
            self.fl = None;
            if itemp == ntemp {
                break;
            }
            itemp += 1;
            c.repoz(self.nscr2);
            c.contio(self.nin, 0, 0, &mut self.scr.all(), &mut nb, &mut nw);
        }

        //--thermr is finished.
        c.totend(self.nin, self.nout, 0, &mut self.scr.all());
        c.repoz(self.nendf);
        c.repoz(self.nin);
        c.repoz(self.nout);
        c.closz(self.nendf);
        c.closz(self.nin);
        c.closz(self.nout);
        c.closz(-iold);
        c.closz(-inew);
        c.closz(self.nscr);
        c.closz(self.nscr2);
        self.scr = RV::new(0);
        self.bufn = RV::new(0);
        self.bufo = RV::new(0);
        self.esi = None;
        self.xsi = None;
        self.fl = None;
        let time = c.timer();
        c.io.so(&format!("{}{}s\n {}*******", " ".repeat(69), ffmt::f(time, 8, 1), "**********".repeat(7)));
    }

    /// `rdelas`: read the ENDF-6 elastic data (MF7/MT2) into `fl`.
    fn rdelas(&mut self, c: &mut Ctx, temp: f64, lthr: i32, nin: i32, nwds: &mut i32, indexc: &mut i32, indexi: &mut i32) {
        const NA: i32 = 1000000;
        let mut a = RV::new(NA as usize + SLACK);
        let (mut nb, mut nw) = (0i32, 0i32);
        *indexc = -1;
        *indexi = -1;

        let mut l = 1i32;
        c.tab1io(nin, 0, 0, &mut a.at(l), &mut nb, &mut nw);
        l += nw;
        while nb != 0 {
            c.moreio(nin, 0, 0, &mut a.at(l), &mut nb, &mut nw);
            l += nw;
            if l > NA {
                c.error("rdelas", "too much elastic data", "increase na");
            }
        }
        let np = c.n2h;
        let nr = c.n1h;
        let lt = c.l1h;
        *nwds = l - 1;

        if lthr != 2 {
            *indexc = 1;
            for _j in 1..=lt {
                let mut l = *nwds + 1;
                c.listio(nin, 0, 0, &mut a.at(l), &mut nb, &mut nw);
                l += nw;
                while nb != 0 {
                    c.moreio(nin, 0, 0, &mut a.at(l), &mut nb, &mut nw);
                    l += nw;
                }
                let tnow = a[*nwds + 1];
                if (tnow - temp).abs() < temp / 1000.0 + 5.0 {
                    let mut l = *nwds + 6;
                    let mut k = 6 + 2 * nr;
                    for _i in 1..=np {
                        l += 1;
                        k += 2;
                        a[k] = a[l];
                    }
                }
            }
        }

        if lthr == 2 {
            *indexi = 1;
        } else if lthr == 3 {
            *indexi = *nwds + 1;
            let mut l = *indexi;
            c.tab1io(nin, 0, 0, &mut a.at(l), &mut nb, &mut nw);
            l += nw;
            while nb != 0 {
                c.moreio(nin, 0, 0, &mut a.at(l), &mut nb, &mut nw);
                l += nw;
                if l > NA {
                    c.error("rdelas", "too much elastic data", "increase na");
                }
            }
            *nwds = l - 1;
        }

        let mut fl = RV::new(*nwds as usize + SLACK);
        for l in 1..=*nwds {
            fl[l] = a[l];
        }
        self.fl = Some(fl);
    }
}

/// Effective temperatures for the sct approximation.
const GATEFF_TABL: [[f64; 3]; 67] = [
    [1002.0, 296.0, 1396.8],
    [1002.0, 350.0, 1411.6],
    [1002.0, 400.0, 1427.4],
    [1002.0, 450.0, 1444.9],
    [1002.0, 500.0, 1464.1],
    [1002.0, 600.0, 1506.8],
    [1002.0, 800.0, 1605.8],
    [1002.0, 1000.0, 1719.8],
    [1004.0, 296.0, 940.91],
    [1004.0, 350.0, 961.62],
    [1004.0, 400.0, 982.93],
    [1004.0, 450.0, 1006.1],
    [1004.0, 500.0, 1030.9],
    [1004.0, 600.0, 1085.1],
    [1004.0, 800.0, 1209.0],
    [1004.0, 1000.0, 1350.0],
    [1064.0, 296.0, 405.64],
    [1064.0, 400.0, 484.22],
    [1064.0, 500.0, 568.53],
    [1064.0, 600.0, 657.66],
    [1064.0, 700.0, 749.69],
    [1064.0, 800.0, 843.63],
    [1064.0, 1000.0, 1035.0],
    [1064.0, 1200.0, 1229.3],
    [1065.0, 296.0, 713.39],
    [1065.0, 400.0, 754.68],
    [1065.0, 500.0, 806.67],
    [1065.0, 600.0, 868.38],
    [1065.0, 700.0, 937.64],
    [1065.0, 800.0, 1012.7],
    [1065.0, 1000.0, 1174.9],
    [1065.0, 1200.0, 1348.2],
    [1065.0, 1600.0, 1712.9],
    [1065.0, 2000.0, 2091.0],
    [1095.0, 296.0, 1165.9],
    [1095.0, 350.0, 1177.8],
    [1095.0, 400.0, 1191.4],
    [1095.0, 450.0, 1207.7],
    [1095.0, 500.0, 1226.0],
    [1095.0, 600.0, 1268.7],
    [1095.0, 800.0, 1373.4],
    [1095.0, 1000.0, 1497.7],
    [1096.0, 296.0, 317.27],
    [1096.0, 400.0, 416.29],
    [1096.0, 500.0, 513.22],
    [1096.0, 600.0, 611.12],
    [1096.0, 700.0, 709.60],
    [1096.0, 800.0, 808.43],
    [1096.0, 1000.0, 1006.8],
    [1096.0, 1200.0, 1205.7],
    [1097.0, 296.0, 806.79],
    [1097.0, 400.0, 829.98],
    [1097.0, 500.0, 868.44],
    [1097.0, 600.0, 920.08],
    [1097.0, 700.0, 981.82],
    [1097.0, 800.0, 1051.1],
    [1097.0, 1000.0, 1205.4],
    [1097.0, 1200.0, 1373.4],
    [1099.0, 296.0, 596.4],
    [1099.0, 400.0, 643.9],
    [1099.0, 500.0, 704.6],
    [1099.0, 600.0, 775.3],
    [1099.0, 800.0, 935.4],
    [1099.0, 1000.0, 1109.8],
    [1099.0, 1200.0, 1292.3],
    [1114.0, 296.0, 1222.0],
    [1114.0, 350.0, 1239.0],
];

const GATEF2_TABL: [[f64; 3]; 16] = [
    [1095.0, 296.0, 685.54],
    [1095.0, 350.0, 712.02],
    [1095.0, 400.0, 738.97],
    [1095.0, 450.0, 768.10],
    [1095.0, 500.0, 799.22],
    [1095.0, 600.0, 866.63],
    [1095.0, 800.0, 1017.3],
    [1095.0, 1000.0, 1182.3],
    [1099.0, 296.0, 427.8],
    [1099.0, 400.0, 502.8],
    [1099.0, 500.0, 584.3],
    [1099.0, 600.0, 671.3],
    [1099.0, 700.0, 761.6],
    [1099.0, 800.0, 854.2],
    [1099.0, 1000.0, 1043.7],
    [1099.0, 1200.0, 1236.6],
];

fn gat_common(tabl: &[[f64; 3]], temp: &RV, ef: &mut RV, ntemp: i32, mat: i32) {
    for i in 1..=ntemp {
        if ef[i] == 0.0 {
            for row in tabl {
                let jmat = nint(row[0]);
                if jmat == mat {
                    let test = 5.0;
                    let diff = row[1] - temp[i];
                    if diff.abs() <= test {
                        ef[i] = row[2];
                    }
                }
            }
            if ef[i] == 0.0 {
                ef[i] = temp[i];
            }
        }
    }
}

/// `gateff`: default effective temperatures.
fn gateff(temp: &RV, eftemp: &mut RV, ntemp: i32, mat: i32) {
    gat_common(&GATEFF_TABL, temp, eftemp, ntemp, mat);
}

/// `gatef2`: default effective temperatures for the second atom.
fn gatef2(temp: &RV, eftmp2: &mut RV, ntemp: i32, mat: i32) {
    gat_common(&GATEF2_TABL, temp, eftmp2, ntemp, mat);
}

/// `terp`: Lagrangian interpolation of order il1 in a table (x, y).
/// `x` and `y` are 1-based views of length nl.
fn terp(x: &[f64], y: &[f64], nl: i32, arg: f64, il1: i32) -> f64 {
    const SMALL: f64 = 1.0e-10;
    let xx = |k: i32| x[(k - 1) as usize];
    let yy = |k: i32| y[(k - 1) as usize];
    let mut il = il1;
    let l: i32;
    'l260: {
        if nl <= il {
            if nl != il {
                il = nl;
            }
            l = 1;
            break 'l260;
        }
        let il2 = il / 2;
        let mut iadd = il % 2;
        let (ilow, ihi, iusel, iuseh, ibeg, iend, last);
        if xx(1) <= xx(nl) {
            ilow = il2 + 1;
            ihi = nl - il2 - iadd;
            iusel = 1;
            iuseh = nl - il + 1;
            ibeg = ilow + 1;
            iend = ihi - 1;
            last = iend - il2 + 1;
            iadd = 0;
        } else {
            ilow = nl - il2;
            ihi = il2 + iadd + 1;
            iusel = nl - il + 1;
            iuseh = 1;
            ibeg = ihi + 1;
            iend = ilow - 1;
            last = 2;
            iadd = 1 - iadd;
        }
        if (arg - xx(ilow)).abs() < SMALL * arg {
            return yy(ilow);
        }
        if !(arg > xx(ilow)) {
            l = iusel;
            break 'l260;
        }
        if (xx(ihi) - arg).abs() < SMALL * arg {
            return yy(ihi);
        }
        if !(xx(ihi) > arg) {
            l = iuseh;
            break 'l260;
        }
        for n in ibeg..=iend {
            let m = if iusel > 1 { nl - n + 1 } else { n };
            if (xx(m) - arg).abs() < SMALL * arg {
                return yy(m);
            }
            if xx(m) > arg {
                l = m - il2 + iadd;
                break 'l260;
            }
        }
        l = last;
    }
    let mut sum = 0.0f64;
    for i in 1..=il {
        let mut p = 1.0f64;
        let mut pk = 1.0f64;
        let inn = l + i - 1;
        for ip in 1..=il {
            if ip != i {
                let inp = l + ip - 1;
                p *= arg - xx(inp);
                pk *= xx(inn) - xx(inp);
            }
        }
        sum += p * yy(inn) / pk;
    }
    sum
}

/// `form`: structure factor for the built-in coherent elastic lattices.
fn form(lat: i32, l1: i32, l2: i32, l3: i32) -> f64 {
    const BEO1: f64 = 7.54e0;
    const BEO2: f64 = 4.24e0;
    const BEO3: f64 = 11.31e0;
    if lat == 1 {
        let i = l3 / 2;
        if 2 * i != l3 {
            let s = sin(PI * (l1 - l2) as f64 / 3.0);
            s * s
        } else {
            (6.0 + 10.0 * cos(2.0 * PI * (l1 - l2) as f64 / 3.0)) / 4.0
        }
    } else if lat == 2 {
        1.0 + cos(2.0 * PI * (2 * l1 + 4 * l2 + 3 * l3) as f64 / 6.0)
    } else if lat == 3 {
        (1.0 + cos(2.0 * PI * (2 * l1 + 4 * l2 + 3 * l3) as f64 / 6.0))
            * (BEO1 + BEO2 + BEO3 * cos(3.0 * PI * l3 as f64 / 4.0))
    } else {
        // the Fortran leaves form undefined here
        0.0
    }
}

/// `upstk`: push (e, s) onto the stack below the top.
/// stk is column-major with nx rows; columns are 1-based.
fn upstk(e: f64, s: &[f64], nl: i32, nx: i32, i: &mut i32, stk: &mut [f64]) {
    let at = |r: i32, c: i32| ((r - 1) + nx * (c - 1)) as usize;
    *i += 1;
    stk[at(1, *i)] = stk[at(1, *i - 1)];
    stk[at(1, *i - 1)] = e;
    for j in 1..=nl {
        stk[at(1 + j, *i)] = stk[at(1 + j, *i - 1)];
        stk[at(1 + j, *i - 1)] = s[j as usize];
    }
}

impl Thermr {
    /// `coh`: coherent elastic scattering cross section on a union grid.
    fn coh(
        &mut self, c: &mut Ctx, lat: i32, itemp: i32, iold: &mut i32, inew: &mut i32, ne: &mut i32, nex: &mut i32,
        mtref: i32, index: i32,
    ) {
        const NLMAX: i32 = 6;
        const IMAX: i32 = 20;
        const HALF: f64 = 0.5;
        const SMALL: f64 = 1.0e-10;
        const TOLMIN: f64 = 1.0e-6;
        const EPS: f64 = 3.0e-5;
        let (mut nb, mut nw) = (0i32, 0i32);
        let mut s = [0.0f64; NLMAX as usize + 2];
        let mut ej = [0.0f64; 21 + 8];
        let mut ex = [0.0f64; 21 + 8];
        let mut x = [0.0f64; 6 + 4];
        let mut y = [0.0f64; 6 + 4];
        let mut z = [0.0f64; 6 + 4];
        let mut b = RV::new(12 + 8);

        //--initialize
        let nl = 1i32;
        let temp = self.tempr[itemp];
        if nl > NLMAX {
            c.error("coh", "too many legendre orders.", " ");
        }
        let nx = nl + 1;
        let nj = *nex + nl;
        for i in 1..=nl {
            ex[(*nex + i) as usize] = 0.0;
        }
        let mut nlt = 5i32;
        if *ne < nlt {
            nlt = *ne;
        }
        let mut nlt1 = nlt - 1;
        let mut stk = vec![0.0f64; (nx * (IMAX + 8)) as usize];
        let sk = |r: i32, cc: i32| ((r - 1) + nx * (cc - 1)) as usize;

        //--determine the bragg edges
        let mut nbragg = nl;
        let mut nlv = nl;
        let mut e = 0.0f64;
        let mut enext = 0.0f64;
        let mut en = 0.0f64;
        let emax = self.emax;
        let natom = self.natom;
        self.sigcoh(c, e, &mut enext, &mut s, &mut nbragg, lat, temp, emax, natom, index);
        let mut ix = 1i32;
        let mut j = 0i32;
        let mut iex = 0i32;
        let mut i = 0i32;
        let mut xm = 0.0f64;
        let mut ym = 0.0f64;
        let mut st = 100;
        loop {
            match st {
                100 => {
                    iex += 1;
                    c.finda(iex, &mut ex[1..], *nex, *iold, &mut self.bufo.0, NBUF);
                    if ex[1] > enext * (1.0 + SMALL) {
                        st = 105;
                        continue;
                    }
                    x[1] = ex[1];
                    y[1] = ex[2];
                    z[1] = ex[*nex as usize];
                    j += 1;
                    c.loada(j, &ex[1..], nj, *inew, &mut self.bufn.0, NBUF);
                }
                105 => {
                    ix += 1;
                    x[ix as usize] = ex[1];
                    y[ix as usize] = ex[2];
                    z[ix as usize] = ex[*nex as usize];
                    if ix < nlt {
                        st = 100;
                        continue;
                    }
                    e = enext;
                    self.sigcoh(c, e, &mut enext, &mut s, &mut nlv, lat, temp, emax, natom, index);
                    stk[sk(1, 1)] = e;
                    for il in 1..=nlv {
                        stk[sk(1 + il, 1)] = s[il as usize];
                    }
                    i = 1;
                    st = 120;
                }
                120 => {
                    e = enext;
                    self.sigcoh(c, e, &mut enext, &mut s, &mut nlv, lat, temp, emax, natom, index);
                    upstk(e, &s, nlv, nx, &mut i, &mut stk);
                    st = 125;
                }
                125 => {
                    st = 140;
                    for k in 1..=nlt {
                        ix = k;
                        if x[ix as usize] > stk[sk(1, i)] * (1.0 + SMALL) {
                            st = 135;
                            break;
                        }
                    }
                    if st == 140 {
                        ix = nlt + 1;
                    }
                }
                135 => {
                    if !(x[ix as usize] >= stk[sk(1, i - 1)] * (1.0 - SMALL)) {
                        e = x[ix as usize];
                        self.sigcoh(c, e, &mut en, &mut s, &mut nlv, lat, temp, emax, natom, index);
                        upstk(e, &s, nlv, nx, &mut i, &mut stk);
                    }
                    st = 140;
                }
                140 => {
                    st = 160;
                    if i == IMAX {
                        continue;
                    }
                    xm = HALF * (stk[sk(1, i - 1)] + stk[sk(1, i)]);
                    xm = sigfig(xm, 7, 0);
                    if stk[sk(1, i - 1)] - stk[sk(1, i)] < EPS * xm {
                        continue;
                    }
                    self.sigcoh(c, xm, &mut en, &mut s, &mut nlv, lat, temp, emax, natom, index);
                    for il in 1..=nlv {
                        c.terp1(stk[sk(1, i)], stk[sk(1 + il, i)], stk[sk(1, i - 1)], stk[sk(1 + il, i - 1)], xm, &mut ym, 2);
                        let mut test = self.tol * s[il as usize].abs();
                        if test < TOLMIN {
                            test = TOLMIN;
                        }
                        if (s[il as usize] - ym).abs() > test {
                            st = 210;
                            break;
                        }
                    }
                }
                160 => {
                    j += 1;
                    ej[1] = stk[sk(1, i)];
                    st = 170;
                }
                170 => {
                    if ej[1] <= x[3] * (1.0 + SMALL) || iex == *ne {
                        st = 190;
                        continue;
                    }
                    for k in 1..=nlt1 as usize {
                        x[k] = x[k + 1];
                        y[k] = y[k + 1];
                        z[k] = z[k + 1];
                    }
                    iex += 1;
                    c.finda(iex, &mut ex[1..], *nex, *iold, &mut self.bufo.0, NBUF);
                    x[nlt as usize] = ex[1];
                    y[nlt as usize] = ex[2];
                    z[nlt as usize] = ex[*nex as usize];
                    if iex == *ne {
                        nlt -= 1;
                    }
                    if iex == *ne {
                        nlt1 -= 1;
                    }
                }
                190 => {
                    ej[2] = terp(&x[1..], &y[1..], nlt, ej[1], nlt1);
                    ej[*nex as usize] = terp(&x[1..], &z[1..], nlt, ej[1], nlt1);
                    let hi = stk[sk(1, i)] > self.emax * (1.0 + SMALL);
                    if hi {
                        ej[*nex as usize] = 0.0;
                    }
                    for il in 1..=nlv {
                        ej[(*nex + il) as usize] = stk[sk(1 + il, i)];
                        if hi {
                            ej[(*nex + il) as usize] = 0.0;
                        }
                    }
                    if hi {
                        st = 230;
                        continue;
                    }
                    c.loada(j, &ej[1..], nj, *inew, &mut self.bufn.0, NBUF);
                    i -= 1;
                    st = if i > 1 { 125 } else { 120 };
                }
                210 => {
                    upstk(xm, &s, nlv, nx, &mut i, &mut stk);
                    st = 125;
                }
                _ => {
                    // 230
                    *ne = j;
                    j = -j;
                    c.loada(j, &ej[1..], nj, *inew, &mut self.bufn.0, NBUF);
                    std::mem::swap(iold, inew);
                    *nex = nj;
                    let ltt = 7;
                    b[1] = self.za;
                    b[2] = self.awr;
                    b[3] = 0.0;
                    b[4] = ltt as f64;
                    b[5] = 0.0;
                    b[6] = 0.0;
                    c.math = self.matdp;
                    c.mfh = 6;
                    c.mth = mtref;
                    c.contio(0, 0, self.nscr, &mut b.all(), &mut nb, &mut nw);
                    b[1] = 1.0;
                    b[2] = 1.0;
                    b[3] = -nbragg as f64;
                    b[4] = 0.0;
                    b[5] = 1.0;
                    b[6] = 2.0;
                    b[7] = 2.0;
                    b[8] = 2.0;
                    b[9] = 1.0e-5;
                    b[10] = 1.0;
                    b[11] = self.emax;
                    b[12] = 1.0;
                    nw = 12;
                    c.tab1io(0, 0, self.nscr, &mut b.all(), &mut nb, &mut nw);
                    self.ncdse = 3;
                    c.asend(0, self.nscr);
                    return;
                }
            }
        }
    }

    /// `sigcoh`: coherent elastic cross section (built-in lattices or the
    /// ENDF-6 Bragg edges in `fl`). e=0 initializes.
    fn sigcoh(
        &mut self, c: &mut Ctx, e: f64, enext: &mut f64, s: &mut [f64], nl: &mut i32, lat: i32, temp: f64, emax: f64,
        natom: i32, index: i32,
    ) {
        const ND: i32 = 10;
        const DWF1: [f64; 10] = [2.1997e0, 2.7448e0, 3.2912e0, 3.8510e0, 4.4210e0, 4.9969e0, 6.1624e0, 7.3387e0, 9.6287e0, 11.992e0];
        const DWF2: [f64; 10] = [3.16663e0, 3.88842e0, 4.62944e0, 5.40517e0, 6.19880e0, 7.0042e0, 8.63665e0, 10.2865e0, 0.0, 0.0];
        const DWF3: [f64; 10] = [2.153e0, 2.6374e0, 3.1348e0, 3.6513e0, 4.1798e0, 4.7164e0, 5.8052e0, 6.9068e0, 0.0, 0.0];
        const TMP: [f64; 10] = [296.0, 400.0, 500.0, 600.0, 700.0, 800.0, 1000.0, 1200.0, 1600.0, 2000.0];
        const GR1: f64 = 2.4573e-8;
        const GR2: f64 = 6.700e-8;
        const GR3: f64 = 12.011e0;
        const GR4: f64 = 5.50e0;
        const BE1: f64 = 2.2856e-8;
        const BE2: f64 = 3.5832e-8;
        const BE3: f64 = 9.01e0;
        const BE4: f64 = 7.53e0;
        const BEO1: f64 = 2.695e-8;
        const BEO2: f64 = 4.39e-8;
        const BEO3: f64 = 12.5e0;
        const BEO4: f64 = 1.0e0;
        const SQRT3: f64 = 1.732050808e0;
        const CW: f64 = 0.658173e-15;
        const HALF: f64 = 0.5;
        const EPS: f64 = 0.05;
        const ZERO: f64 = 0.0;

        if e > ZERO {
            //--compute cross sections at this energy (label 210)
            let fl = self.fl.as_ref().unwrap();
            let re = 1.0 / e;
            for il in 1..=*nl as usize {
                s[il] = 0.0;
            }
            let mut last = 0;
            let mut elim = 0.0f64;
            let mut p = [0.0f64; 6 + 64];
            let k = self.sc_k;
            for i in 1..=k {
                let tsq = fl[2 * i - 1];
                elim = tsq * self.sc_recon;
                if elim >= e {
                    break;
                }
                let mut f = fl[2 * i];
                if e > emax {
                    f = 0.0;
                }
                let u = 1.0 - 2.0 * elim * re;
                let lmax = *nl - 1;
                legndr(u, &mut p, lmax);
                for il in 1..=*nl as usize {
                    s[il] += f * p[il - 1];
                }
                if i == k {
                    last = 1;
                }
            }
            for il in 1..=*nl as usize {
                s[il] = s[il] * self.sc_scon * re;
            }
            if last == 1 {
                elim = emax;
            }
            if elim > emax {
                elim = emax;
            }
            *enext = sigfig(elim, 7, -1);
            if e > sigfig(*enext, 7, -1) {
                *enext = sigfig(elim, 7, 1);
            }
            return;
        }

        //--initialize
        let twopis = (2.0 * PI) * (2.0 * PI);
        let amne = AMASSN * AMU;
        let econ = EV * 8.0 * (amne / HBAR) / HBAR;
        self.sc_recon = 1.0 / econ;
        let tsqx = econ / 20.0;
        let nord = 2;
        if lat == 10 {
            //--ENDF-6 bragg edges (label 200)
            let fl = self.fl.as_mut().unwrap();
            let nr = int(fl[index + 4]);
            let np = int(fl[index + 5]);
            self.sc_k = np;
            *nl = np;
            let mut blast = 0.0f64;
            self.sc_scon = 1.0;
            *enext = fl[index + 2 * nr + 6];
            *enext = sigfig(*enext, 7, -1);
            for i in 1..=*nl {
                let l = index + 2 * (i - 1);
                fl[l] = fl[l + 2 * nr + 6] * econ;
                fl[l + 1] = fl[l + 2 * nr + 7] - blast;
                blast = fl[l + 2 * nr + 7];
            }
            return;
        }
        let (a, cc, amsc, scoh, wal2);
        if lat == 1 {
            a = GR1;
            cc = GR2;
            amsc = GR3;
            scoh = GR4 / natom as f64;
            wal2 = terp(&TMP, &DWF1, ND, temp, nord);
        } else if lat == 2 {
            a = BE1;
            cc = BE2;
            amsc = BE3;
            scoh = BE4 / natom as f64;
            wal2 = terp(&TMP, &DWF2, ND, temp, nord);
        } else if lat == 3 {
            a = BEO1;
            cc = BEO2;
            amsc = BEO3;
            scoh = BEO4 / natom as f64;
            wal2 = terp(&TMP, &DWF3, ND, temp, nord);
        } else {
            c.error("sigcoh", "illegal lat.", " ");
        }
        let c1 = 4.0 / (3.0 * a * a);
        let c2 = 1.0 / (cc * cc);
        self.sc_scon = scoh * ((4.0 * PI) * (4.0 * PI)) / (2.0 * a * a * cc * SQRT3 * econ);
        let wint = CW * amsc * wal2;
        let t2 = HBAR / (2.0 * AMU * amsc);
        let ulim = econ * emax;
        let mut nw = 1000000i32;
        let mut wrk = RV::new(nw as usize + 4);

        //--compute and sort lattice factors
        let tausq = |m1: i32, m2: i32, m3: i32| (c1 * (m1 * m1 + m2 * m2 + m1 * m2) as f64 + ((m3 * m3) as f64 * c2)) * twopis;
        let phi = ulim / twopis;
        let mut i1m = int(a * sqrt(phi));
        i1m += 1;
        let mut k = 0i32;
        // add (tsq, f) into the list, merging close values above tsqx
        let add = |c: &mut Ctx, k: &mut i32, wrk: &mut RV, tsq: f64, f: f64| {
            if !(*k > 0 && tsq > tsqx) {
                *k += 1;
                if 2 * *k > nw {
                    c.error("sigcoh", "storage exceeded.", " ");
                }
                wrk[2 * *k - 1] = tsq;
                wrk[2 * *k] = f;
                return;
            }
            for i in 1..=*k {
                if tsq < wrk[2 * i - 1] || tsq >= (1.0 + EPS) * wrk[2 * i - 1] {
                    continue;
                }
                wrk[2 * i] += f;
                return;
            }
            *k += 1;
            if 2 * *k > nw {
                c.error("sigcoh", "storage exceeded.", " ");
            }
            wrk[2 * *k - 1] = tsq;
            wrk[2 * *k] = f;
        };
        for i1 in 1..=i1m {
            let l1 = i1 - 1;
            let mut i2m = int(HALF * (l1 as f64 + sqrt(3.0 * (a * a * phi - (l1 * l1) as f64))));
            i2m += 1;
            for i2 in i1..=i2m {
                let l2 = i2 - 1;
                let x = phi - c1 * (l1 * l1 + l2 * l2 - l1 * l2) as f64;
                let mut i3m = 0;
                if x > ZERO {
                    i3m = int(cc * sqrt(x));
                }
                i3m += 1;
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
                        w2 = HALF;
                    }
                    let mut w3 = 2.0f64;
                    if l3 == 0 {
                        w3 = 1.0;
                    }
                    let tsq = tausq(l1, l2, l3);
                    if !(tsq <= ZERO || tsq > ulim) {
                        let tau = sqrt(tsq);
                        let w = exp(-(tsq * t2 * wint)) * w1 * w2 * w3 / tau;
                        let f = w * form(lat, l1, l2, l3);
                        add(c, &mut k, &mut wrk, tsq, f);
                    }
                    let tsq = tausq(l1, -l2, l3);
                    if tsq <= ZERO || tsq > ulim {
                        continue;
                    }
                    let tau = sqrt(tsq);
                    let w = exp(-(tsq * t2 * wint)) * w1 * w2 * w3 / tau;
                    let f = w * form(lat, l1, -l2, l3);
                    add(c, &mut k, &mut wrk, tsq, f);
                }
            }
        }
        let imax = k - 1;
        for i in 1..=imax {
            let jmin = i + 1;
            for j in jmin..=k {
                if wrk[2 * j - 1] < wrk[2 * i - 1] {
                    let st = wrk[2 * i - 1];
                    let sf = wrk[2 * i];
                    wrk[2 * i - 1] = wrk[2 * j - 1];
                    wrk[2 * i] = wrk[2 * j];
                    wrk[2 * j - 1] = st;
                    wrk[2 * j] = sf;
                }
            }
        }
        k += 1;
        wrk[2 * k - 1] = ulim;
        wrk[2 * k] = wrk[2 * k - 2];
        nw = 2 * k;
        *enext = self.sc_recon * wrk[1];
        *enext = sigfig(*enext, 7, -1);
        let mut fl = RV::new(nw as usize + SLACK);
        for i in 1..=nw {
            fl[i] = wrk[i];
        }
        self.fl = Some(fl);
        self.sc_k = k;
        *nl = k;
    }

    /// `iel`: incoherent elastic scattering.
    fn iel(
        &mut self, c: &mut Ctx, mat: i32, itemp: i32, iold: &mut i32, inew: &mut i32, ne: &mut i32, nex: &mut i32,
        mtref: i32, index: i32,
    ) {
        const NUPMAX: i32 = 10;
        const TMP: [f64; 8] = [296.0, 400.0, 500.0, 600.0, 700.0, 800.0, 1000.0, 1200.0];
        const DWH: [f64; 8] = [8.4795e0, 9.0854e0, 9.8196e0, 10.676e0, 11.625e0, 12.643e0, 14.822e0, 17.125e0];
        const DWZ: [f64; 8] = [1.9957e0, 2.6546e0, 3.2946e0, 3.9380e0, 4.5835e0, 5.2302e0, 6.5260e0, 7.8236e0];
        const C11A: f64 = 162.88e0;
        const C11B: f64 = 296.0e0;
        const C11C: f64 = 34.957e0;
        const C11D: f64 = 350.0e0;
        const C11E: f64 = 40.282e0;
        const C12A: f64 = 81.44e0;
        const C13A: f64 = 6.3366e0;
        const UP: f64 = 1.1e0;
        const DN: f64 = 0.9e0;
        let (mut nb, mut nw) = (0i32, 0i32);
        let mut ex = [0.0f64; 21 + 8];
        let mut ej = [0.0f64; 21 + 8];
        let mut dwa = 0.0f64;

        if self.iprint == 2 {
            c.io.so("\n iel\n");
        }
        let temp = self.tempr[itemp];
        let nj = *nex + 1;

        //--initialize
        if mat == 11 {
            self.sb = C11A;
            dwa = C11C + (temp - C11B) * (C11E - C11C) / (C11D - C11B);
        } else if mat == 12 {
            dwa = terp(&TMP, &DWH, 8, temp, 3);
            self.sb = C12A;
        } else if mat == 13 {
            dwa = terp(&TMP, &DWZ, 8, temp, 3);
            self.sb = C13A;
        } else if mat == 20 {
            let fl = self.fl.as_mut().unwrap();
            self.sb = fl[index];
            let nr = int(fl[index + 4]);
            let np = int(fl[index + 5]);
            if np == 1 {
                let tt1 = fl[index + 2 * nr + 6];
                if (temp - tt1).abs() > temp / 10.0 {
                    c.error("iel", "bad temperature for debye-waller factor", " ");
                }
                dwa = fl[index + 2 * nr + 7];
            } else {
                let tt1 = fl[index + 2 * nr + 6];
                let ttn = fl[index + 2 * nr + 2 * np + 4];
                if temp < DN * tt1 || temp > UP * ttn {
                    c.error("iel", "bad temperature for debye-waller factor", " ");
                }
                if tt1 > temp {
                    fl[index + 2 * nr + 6] = temp;
                }
                if ttn < temp {
                    fl[index + 2 * nr + 2 * np + 4] = temp;
                }
                let mut ip = 2;
                let mut ir = 1;
                let mut tnxt = 0.0;
                let mut idis = 0;
                c.terpa(&mut dwa, temp, &mut tnxt, &mut idis, &fl.at(index), &mut ip, &mut ir);
            }
        } else {
            c.error("iel", "unknown material identifier.", " ");
        }
        let c1 = self.sb / (2 * self.natom) as f64;

        //--determine the size of the esi grid
        let esi = self.esi.as_ref().unwrap();
        let mut nne = 0i32;
        while esi[nne + 1] != 0.0 {
            nne += 1;
        }
        let mut xie = RV::new(nne.max(1) as usize + 8);

        //--write head and tab records for mf6
        let ltt = 6;
        c.math = self.matdp;
        c.mfh = 6;
        c.mth = mtref;
        let scr = &mut self.scr;
        scr[1] = self.za;
        scr[2] = self.awr;
        scr[3] = 0.0;
        scr[4] = ltt as f64;
        scr[5] = 0.0;
        scr[6] = 0.0;
        nw = 6;
        c.contio(0, 0, self.nscr, &mut scr.all(), &mut nb, &mut nw);
        scr[1] = 1.0;
        scr[2] = 1.0;
        scr[3] = -2.0;
        scr[4] = 1.0;
        scr[5] = 1.0;
        scr[6] = 2.0;
        scr[7] = 2.0;
        scr[8] = 2.0;
        scr[9] = 1.0e-5;
        scr[10] = 1.0;
        scr[11] = self.emax;
        scr[12] = 1.0;
        nw = 12;
        c.tab1io(0, 0, self.nscr, &mut scr.all(), &mut nb, &mut nw);
        scr[1] = temp;
        scr[2] = 0.0;
        scr[3] = 3.0;
        scr[5] = 1.0;
        scr[6] = nne as f64;
        scr[7] = nne as f64;
        scr[8] = 2.0;
        nw = 8;
        c.tab2io(0, 0, self.nscr, &mut scr.all(), &mut nb, &mut nw);
        let n = self.nbin;
        let mut nup = n;
        if nup > NUPMAX {
            nup = NUPMAX;
        }
        let nup1 = nup + 1;
        nw = n + 8;

        //--compute the cross section and the equally probable cosines
        for iex in 1..=nne {
            let e = esi[iex];
            let c2 = 2.0 * e * dwa;
            let mut u = -1.0f64;
            let rc2 = 1.0 / c2;
            let x1 = exp(-(2.0 * c2));
            let r1x1 = 1.0 / (1.0 - x1);
            let xsec = c1 * rc2 * (1.0 - x1);
            scr[1] = 0.0;
            scr[2] = e;
            scr[3] = 0.0;
            scr[4] = 0.0;
            scr[5] = (n + 2) as f64;
            scr[6] = (n + 2) as f64;
            scr[7] = e;
            scr[8] = 1.0;
            for iu in 1..=n {
                let x2 = exp(-(c2 * (1.0 - u)));
                let unow = 1.0 + rc2 * log((1.0 - x1) / n as f64 + x2);
                scr[8 + iu] = n as f64 * rc2 * (exp(-(c2 * (1.0 - unow))) * (c2 * unow - 1.0) - x2 * (c2 * u - 1.0)) * r1x1;
                u = unow;
            }
            xie[iex] = xsec;
            if self.iprint == 2 {
                let mut line = format!("    {}{}", ffmt::e1p(e, 12, 4), ffmt::e1p(xsec, 12, 4));
                for iu in 1..=nup {
                    line.push_str(&ffmt::f(scr[8 + iu], 9, 4));
                }
                c.io.so(&line);
            }
            if self.iprint == 2 && n > nup {
                let mut line = " ".repeat(28);
                let mut cnt = 0;
                for iu in nup1..=n {
                    if cnt == 10 {
                        c.io.so(&line);
                        line = " ".repeat(28);
                        cnt = 0;
                    }
                    line.push_str(&ffmt::f(scr[8 + iu], 9, 4));
                    cnt += 1;
                }
                c.io.so(&line);
            }
            c.listio(0, 0, self.nscr, &mut scr.all(), &mut nb, &mut nw);
        }

        //--add the cross section to the scratch tape
        for iex in 1..=*ne {
            c.finda(iex, &mut ex[1..], *nex, *iold, &mut self.bufo.0, NBUF);
            for ix in 1..=*nex as usize {
                ej[ix] = ex[ix];
            }
            ej[nj as usize] = terp(&esi.0, &xie.0, nne, ex[1], 3);
            if iex == *ne {
                ej[nj as usize] = 0.0;
            }
            let mut iet = iex;
            if iex == *ne {
                iet = -iet;
            }
            c.loada(iet, &ej[1..], nj, *inew, &mut self.bufn.0, NBUF);
        }
        std::mem::swap(iold, inew);
        *nex = nj;
        c.asend(0, self.nscr);
        self.ncdse = 5 + *ne * ((n + 11) / 6);
    }
}

/// The S(alpha,beta) table and the module values that `sig`, `sigl` and
/// `sigu` read (fixed while calcem runs over its energy grid).
#[derive(Clone)]
struct Sab {
    nalpha: i32,
    nbeta: i32,
    alpha: RV,
    beta: RV,
    /// sab(nalpha,nbeta), column-major, plus slack
    sab: Vec<f64>,
    az: f64,
    az2: f64,
    tevz: f64,
    iinc: i32,
    lat: i32,
    lasym: i32,
    cliq: f64,
    teff: f64,
    teff2: f64,
    sb: f64,
    sb2: f64,
    thr6: f64,
    /// alpha(2..nalpha) and beta(2..nbeta) are in ascending order, so that
    /// `first_above` can use a binary search
    alpha_up: bool,
    beta_up: bool,
}

/// The loop `do i=1,n1; k=i; if (x.lt.v(i+1)) exit; enddo` (k=0 if n1<1):
/// the first i with x<v(i+1), else n1. With `up` (v(2..n1+1) ascending) a
/// binary search gives the same i.
#[inline]
fn first_above(v: &RV, n1: i32, x: f64, up: bool) -> i32 {
    if n1 < 1 {
        return 0;
    }
    if up {
        let p = v.0[1..=n1 as usize].partition_point(|&w| !(x < w)) as i32;
        return if p < n1 { p + 1 } else { n1 };
    }
    let mut k = 0;
    for i in 1..=n1 {
        k = i;
        if x < v[i + 1] {
            break;
        }
    }
    k
}

/// v(2..n) is in ascending order (no NaN).
fn ascending(v: &RV, n: i32) -> bool {
    n < 2 || v.0[1..n as usize].windows(2).all(|w| w[0] <= w[1])
}

const SABFLG: f64 = -225.0e0;

impl Sab {
    #[inline]
    fn s(&self, ia: i32, ib: i32) -> f64 {
        self.sab[((ia - 1) + self.nalpha * (ib - 1)) as usize]
    }

    /// `sig`: the double-differential scattering cross section.
    fn sig(&self, e: f64, ep: f64, u: f64, tev: f64) -> f64 {
        const SIGMIN: f64 = 1.0e-10;
        const AMIN: f64 = 1.0e-6;
        const TEST1: f64 = 0.2e0;
        const TEST2: f64 = 30.0e0;
        const ZERO: f64 = 0.0;
        let (nalpha, nbeta) = (self.nalpha, self.nbeta);
        let alpha = &self.alpha;
        let beta = &self.beta;
        let az = self.az;
        let tevz = self.tevz;

        let rtev = 1.0 / tev;
        let bb = (ep - e) * rtev;
        let mut a = (e + ep - 2.0 * u * sqrt(e * ep)) / (az * tev);
        if a < AMIN {
            a = AMIN;
        }
        let sigc = sqrt(ep / e) * rtev / 2.0;
        let mut b = bb.abs();
        let c = sqrt(4.0 * PI);
        if self.iinc == 1 {
            //--free-gas scattering (label 200)
            let mut s = 0.0;
            let arg = (a + bb) * (a + bb) / (4.0 * a);
            if -arg > SABFLG {
                s = exp(-arg) / (c * sqrt(a));
            }
            let mut sig = sigc * self.sb * s;
            if sig < SIGMIN {
                sig = 0.0;
            }
            return sig;
        }
        // iinc=2 (other values are stopped by calcem before the first call)
        let mut go170 = false;
        let mut s = 0.0f64;
        'l160: {
            if self.lat == 1 {
                b = b * tev / tevz;
            }
            if self.lat == 1 {
                a = a * tev / tevz;
            }
            if a > alpha[nalpha] {
                go170 = true;
                break 'l160;
            }
            if self.lasym == 1 {
                let mut bbm = bb;
                if self.lat == 1 {
                    bbm = bb * tev / tevz;
                }
                if bbm > beta[nbeta] || bbm < beta[1] {
                    go170 = true;
                    break 'l160;
                }
            } else if b > beta[nbeta] {
                go170 = true;
                break 'l160;
            }
            let nb1 = nbeta - 1;
            let na1 = nalpha - 1;
            let mut bbb = b;
            if self.lasym == 1 && bb < ZERO {
                bbb = -b;
            }
            let mut ib = first_above(beta, nb1, bbb, self.beta_up);
            let mut ia = first_above(alpha, na1, a, self.alpha_up);
            if !(self.cliq == ZERO || a >= alpha[1]) && self.lasym != 1 && !(b > TEST1) {
                s = self.s(1, 1) + log(alpha[1] / a) / 2.0 - self.cliq * (b * b) / a;
                if s < SABFLG {
                    s = SABFLG;
                }
                break 'l160;
            }
            // label 150
            if !(a * az < TEST2 && b < TEST2) {
                if self.s(ia, ib) <= SABFLG
                    || self.s(ia + 1, ib) <= SABFLG
                    || self.s(ia, ib + 1) <= SABFLG
                    || self.s(ia + 1, ib + 1) <= SABFLG
                {
                    go170 = true;
                    break 'l160;
                }
            }
            // label 155
            if ia + 1 == nalpha {
                ia -= 1;
            }
            if ib + 1 == nbeta {
                ib -= 1;
            }
            let s1 = self.terpq(alpha[ia], self.s(ia, ib), alpha[ia + 1], self.s(ia + 1, ib), alpha[ia + 2], self.s(ia + 2, ib), a);
            let s2 = self.terpq(
                alpha[ia],
                self.s(ia, ib + 1),
                alpha[ia + 1],
                self.s(ia + 1, ib + 1),
                alpha[ia + 2],
                self.s(ia + 2, ib + 1),
                a,
            );
            let s3 = self.terpq(
                alpha[ia],
                self.s(ia, ib + 2),
                alpha[ia + 1],
                self.s(ia + 1, ib + 2),
                alpha[ia + 2],
                self.s(ia + 2, ib + 2),
                a,
            );
            s = self.terpq(beta[ib], s1, beta[ib + 1], s2, beta[ib + 2], s3, bbb);
        }
        if !go170 {
            // label 160
            let mut sig = 0.0;
            if s - bb / 2.0 > SABFLG {
                sig = exp(s - bb / 2.0);
            }
            sig = sigc * self.sb * sig;
            if sig < SIGMIN {
                sig = 0.0;
            }
            return sig;
        }
        // label 170: short collision time approximation
        if self.lat == 1 {
            b = b * tevz * rtev;
        }
        if self.lat == 1 {
            a = a * tevz * rtev;
        }
        let tfff = self.teff;
        let tfff2 = self.teff2;
        let mut s = 0.0;
        let arg = (a - b) * (a - b) * tev / (4.0 * a * tfff) + (b + bb) / 2.0;
        if -arg > SABFLG {
            s = exp(-arg) / (c * sqrt(a * tfff * rtev));
        }
        let mut sig = sigc * self.sb * s;
        if self.sb2 > ZERO {
            a = a * az / self.az2;
            let arg = (a - b) * (a - b) * tev / (4.0 * a * tfff2) + (b + bb) / 2.0;
            let mut s2 = 0.0;
            if -arg > SABFLG {
                s2 = exp(-arg) / (c * sqrt(a * tfff2 * rtev));
            }
            sig += sigc * self.sb2 * s2;
        }
        if sig < SIGMIN {
            sig = 0.0;
        }
        sig
    }

    /// `terpq`: quadratic interpolation with checks.
    fn terpq(&self, x1: f64, y1: f64, x2: f64, y2: f64, x3: f64, y3: f64, x: f64) -> f64 {
        const STEP: f64 = 2.0e0;
        let mut y = 0.0f64;
        if x < x1 {
            if y1 > y2 {
                y = y1;
            } else {
                crate::endf::terp1(x1, y1, x2, y2, x, &mut y, 3, self.thr6);
            }
        } else if x > x3 {
            if y3 > y2 {
                y = y3;
            } else {
                crate::endf::terp1(x2, y2, x3, y3, x, &mut y, 2, self.thr6);
            }
        } else if (y1 - y2).abs() > STEP || (y2 - y3).abs() > STEP {
            if x < x2 {
                crate::endf::terp1(x1, y1, x2, y2, x, &mut y, 2, self.thr6);
            } else {
                crate::endf::terp1(x2, y2, x3, y3, x, &mut y, 2, self.thr6);
            }
        } else {
            let mut b = (y2 - y1) * (x3 - x1) / ((x2 - x1) * (x3 - x2));
            b -= (y3 - y1) * (x2 - x1) / ((x3 - x1) * (x3 - x2));
            let mut c = (y3 - y1) / ((x3 - x1) * (x3 - x2));
            c -= (y2 - y1) / ((x2 - x1) * (x3 - x2));
            y = y1 + b * (x - x1) + c * (x - x1) * (x - x1);
        }
        if y < SABFLG {
            y = SABFLG;
        }
        y
    }

    /// `sigl`: the angular distribution at (e, ep) as equally probable
    /// cosines (nlin<0) or Legendre moments (nlin>0). s is 1-based.
    fn sigl(&self, c: &mut Ctx, nbin: &mut i32, e: f64, ep: f64, nlin: i32, tev: f64, s: &mut [f64], tolin: f64) {
        const IMAX: usize = 20;
        const ZERO: f64 = 0.0;
        const XTOL: f64 = 0.00001e0;
        const YTOL: f64 = 0.001e0;
        const ONE: f64 = 1.0;
        const HALF: f64 = 0.5e0;
        const THIRD: f64 = 0.333333333e0;
        const SIGMIN: f64 = 1.0e-32;
        const EPS: f64 = 1.0e-3;
        const SHADE: f64 = 0.99999999e0;
        let mut x = [0.0f64; IMAX + 2];
        let mut y = [0.0f64; IMAX + 2];
        let mut p = [0.0f64; 70];

        let mut b = (ep - e) / tev;
        let tol = HALF * tolin;
        let mut nl = nlin;
        if nl < 0 {
            nl = -nl;
        }
        b = b.abs();
        if self.lat == 1 && self.iinc == 2 {
            b = b * tev / self.tevz;
        }
        let s1bb = sqrt(1.0 + b * b);
        let mut seep = 0.0f64;
        if ep != ZERO {
            seep = 1.0 / sqrt(e * ep);
        }
        let x2_of = |seep: f64| -> f64 {
            let mut x2 = if ep == ZERO { 0.0 } else { HALF * (e + ep - (s1bb - 1.0) * self.az * tev) * seep };
            if x2.abs() > 1.0 - EPS {
                x2 = 0.99e0;
            }
            sigfig(x2, 8, 0)
        };

        //--adaptive calculation of cross section
        let mut i = 3usize;
        let mut sum = 0.0f64;
        x[3] = -1.0;
        let mut xl = x[3];
        y[3] = self.sig(e, ep, x[3], tev);
        let mut yl = y[3];
        x[2] = x2_of(seep);
        y[2] = self.sig(e, ep, x[2], tev);
        x[1] = 1.0;
        y[1] = self.sig(e, ep, x[1], tev);
        let mut ymax = y[2];
        if y[1] > ymax {
            ymax = y[1];
        }
        if y[3] > ymax {
            ymax = y[3];
        }
        if ymax < EPS {
            ymax = EPS;
        }
        let mut at110 = true;
        loop {
            if at110 {
                // label 110
                if i != IMAX {
                    let mut xm = HALF * (x[i - 1] + x[i]);
                    xm = sigfig(xm, 8, 0);
                    let ym = HALF * (y[i - 1] + y[i]);
                    let yt = self.sig(e, ep, xm, tev);
                    let test = tol * yt.abs() + tol * ymax / 50.0;
                    let test2 = ym + ymax / 100.0;
                    let conv = (yt - ym).abs() <= test && (y[i - 1] - y[i]).abs() <= test2 && (x[i - 1] - x[i]) < HALF;
                    if !conv && !(x[i - 1] - x[i] < XTOL) {
                        i += 1;
                        x[i] = x[i - 1];
                        y[i] = y[i - 1];
                        x[i - 1] = xm;
                        y[i - 1] = yt;
                        continue;
                    }
                }
            }
            // label 120
            sum += HALF * (y[i] + yl) * (x[i] - xl);
            xl = x[i];
            yl = y[i];
            i -= 1;
            if i > 1 {
                at110 = true;
                continue;
            }
            if i == 1 {
                at110 = false;
                continue;
            }
            break;
        }
        s[1] = sum;
        if !(sum > SIGMIN) {
            for il in 1..=nl as usize {
                s[il] = 0.0;
            }
            return;
        }

        //--prime the adaptive integration
        *nbin = nl - 1;
        let nbin = *nbin;
        let rnbin = ONE / nbin as f64;
        let fract = sum * rnbin;
        let rfract = 1.0 / fract;
        sum = 0.0;
        let mut gral = 0.0f64;
        for il in 2..=nl as usize {
            s[il] = 0.0;
        }
        let mut j = 0i32;
        i = 3;
        x[3] = -1.0;
        xl = x[3];
        y[3] = self.sig(e, ep, x[3], tev);
        x[2] = x2_of(seep);
        y[2] = self.sig(e, ep, x[2], tev);
        x[1] = 1.0;
        y[1] = self.sig(e, ep, x[1], tev);
        ymax = y[1];
        if y[2] > ymax {
            ymax = y[2];
        }
        if y[3] > ymax {
            ymax = y[3];
        }
        if ymax < EPS {
            ymax = EPS;
        }

        //--compute next point of the adaptive integration
        let mut xn = 0.0f64;
        let mut xil = 0.0f64;
        let mut st = 150;
        loop {
            match st {
                150 => {
                    st = 160;
                    if i == IMAX {
                        continue;
                    }
                    let mut xm = HALF * (x[i - 1] + x[i]);
                    xm = sigfig(xm, 8, 0);
                    let ym = HALF * (y[i - 1] + y[i]);
                    let yt = self.sig(e, ep, xm, tev);
                    let test = tol * yt.abs() + tol * ymax / 50.0;
                    let test2 = ym + ymax / 100.0;
                    if (yt - ym).abs() <= test && (y[i - 1] - y[i]).abs() <= test2 && (x[i - 1] - x[i]) < HALF {
                        continue;
                    }
                    if x[i - 1] - x[i] < XTOL {
                        continue;
                    }
                    i += 1;
                    x[i] = x[i - 1];
                    y[i] = y[i - 1];
                    x[i - 1] = xm;
                    y[i - 1] = yt;
                    st = 150;
                }
                160 => {
                    let add = HALF * (y[i] + yl) * (x[i] - xl);
                    if x[i] == xl {
                        st = 250;
                        continue;
                    }
                    xil = 1.0 / (x[i] - xl);
                    if i == 1 && j == nbin - 1 {
                        // label 165
                        xn = x[i];
                        j += 1;
                        st = 190;
                        continue;
                    }
                    if sum + add >= fract * SHADE && j < nbin - 1 {
                        st = 170;
                        continue;
                    }
                    sum += add;
                    gral = gral
                        + HALF * (yl * x[i] - y[i] * xl) * (x[i] + xl)
                        + THIRD * (y[i] - yl) * (x[i] * x[i] + x[i] * xl + xl * xl);
                    st = 250;
                }
                170 => {
                    j += 1;
                    st = 175;
                    if yl < SIGMIN {
                        continue;
                    }
                    let test = (fract - sum) * (y[i] - yl) / ((x[i] - xl) * (yl * yl));
                    if test.abs() > YTOL {
                        continue;
                    }
                    xn = xl + (fract - sum) / yl;
                    if xn > x[i] {
                        st = 180;
                        continue;
                    }
                    if xn >= xl && xn <= x[i] {
                        st = 190;
                        continue;
                    }
                    c.error("sigl", "no legal solution.", " ");
                }
                175 => {
                    let f = (y[i] - yl) * xil;
                    let rf = 1.0 / f;
                    let mut disc = (yl * rf) * (yl * rf) + 2.0 * (fract - sum) * rf;
                    if disc < ZERO {
                        let strng = format!("disc={}", ffmt::e1p(disc, 12, 4));
                        c.mess("sigl", &strng, "set to abs value and continue");
                        disc = disc.abs();
                    }
                    if f > ZERO {
                        xn = xl - (yl * rf) + sqrt(disc);
                    }
                    if f < ZERO {
                        xn = xl - (yl * rf) - sqrt(disc);
                    }
                    if xn > xl && xn <= x[i] {
                        st = 190;
                        continue;
                    }
                    if xn > xl && xn < (x[i] + YTOL * (x[i] - xl)) {
                        st = 180;
                        continue;
                    }
                    c.error("sigl", "no legal solution (quadratic path).", " ");
                }
                180 => {
                    xn = x[i];
                    st = 190;
                }
                190 => {
                    let yn = yl + (y[i] - yl) * (xn - xl) * xil;
                    gral += (xn - xl)
                        * (yl * HALF * (xn + xl)
                            + (y[i] - yl) * xil * (-(xl * HALF * (xn + xl)) + THIRD * (xn * xn + xn * xl + xl * xl)));
                    let xbar = gral * rfract;
                    if nlin >= 0 {
                        legndr(xbar, &mut p, nl);
                        for il in 2..=nl as usize {
                            s[il] += p[il - 1] * rnbin;
                        }
                    } else {
                        s[(j + 1) as usize] = xbar;
                    }
                    xl = xn;
                    yl = yn;
                    sum = 0.0;
                    gral = 0.0;
                    if j == nbin {
                        return;
                    }
                    st = if xl < x[i] { 160 } else { 250 };
                }
                _ => {
                    // label 250
                    xl = x[i];
                    yl = y[i];
                    i -= 1;
                    if i > 1 {
                        st = 150;
                    } else if i == 1 {
                        st = 160;
                    } else {
                        return;
                    }
                }
            }
        }
    }

    /// `sigu`: the secondary-energy distribution at (e, u); s(1)=integral,
    /// s(2)=number of points, then (ep, pdf) pairs. s is 1-based.
    fn sigu(&self, e: f64, u: f64, tev: f64, s: &mut [f64], nemax: i32, tolin: f64) {
        const IMAX: usize = 20;
        const ZERO: f64 = 0.0;
        const HALF: f64 = 0.5e0;
        const TOLMIN: f64 = 1.0e-6;
        const BMAX: f64 = 20.0;
        let mut x = [0.0f64; IMAX + 2];
        let mut y = [0.0f64; IMAX + 2];
        let nbeta = self.nbeta;
        let beta = &self.beta;
        let az = self.az;
        let tevz = self.tevz;
        let tol = tolin;
        for i in 1..=(2 * nemax) as usize {
            s[i] = 0.0;
        }
        let root1 = (u * sqrt(e) + sqrt(u * u * e + (az - 1.0) * (az + 1.0) * e)) / (az + 1.0);
        let mut sum = 0.0f64;
        x[1] = 0.0;
        y[1] = self.sig(e, x[1], u, tev);
        let mut jbeta = -nbeta;
        if self.lasym > 0 {
            jbeta = 1;
        }
        let mut j = 0i32;
        let mut xl = 0.0f64;
        let mut yl = 0.0f64;
        let mut i = 0usize;
        let mut st = 111;
        loop {
            match st {
                111 => {
                    x[2] = x[1];
                    y[2] = y[1];
                    loop {
                        // label 113
                        if jbeta == 0 {
                            jbeta = 1;
                        }
                        if jbeta <= 0 {
                            if self.lat == 1 {
                                x[1] = e - beta[-jbeta] * tevz;
                            } else {
                                x[1] = e - beta[-jbeta] * tev;
                            }
                            x[1] = sigfig(x[1], 8, 0);
                            if x[1] == e {
                                x[1] = sigfig(e, 8, -1);
                            }
                        } else {
                            if self.lat == 1 {
                                x[1] = e + beta[jbeta] * tevz;
                            } else {
                                x[1] = e + beta[jbeta] * tev;
                            }
                        }
                        if x[1] > x[2] {
                            break;
                        }
                        jbeta += 1;
                    }
                    // label 116
                    let r2 = root1 * root1;
                    if u < ZERO && r2 > (1.01f32 as f64) * x[2] && r2 < x[1] {
                        x[1] = r2;
                    }
                    x[1] = sigfig(x[1], 8, 0);
                    y[1] = self.sig(e, x[1], u, tev);
                    i = 2;
                    st = 150;
                }
                150 => {
                    st = 160;
                    if i == IMAX {
                        continue;
                    }
                    if i > 3 && HALF * (y[i - 1] + y[i]) * (x[i - 1] - x[i]) < TOLMIN {
                        continue;
                    }
                    let mut xm = HALF * (x[i - 1] + x[i]);
                    xm = sigfig(xm, 8, 0);
                    if xm <= x[i] || xm >= x[i - 1] {
                        continue;
                    }
                    let ym = HALF * (y[i - 1] + y[i]);
                    let yt = self.sig(e, xm, u, tev);
                    let test = tol * yt.abs();
                    if (yt - ym).abs() <= test {
                        continue;
                    }
                    i += 1;
                    x[i] = x[i - 1];
                    y[i] = y[i - 1];
                    x[i - 1] = xm;
                    y[i - 1] = yt;
                    st = 150;
                }
                160 => {
                    j += 1;
                    s[(2 * j + 1) as usize] = x[i];
                    s[(2 * j + 2) as usize] = y[i];
                    if j > 1 {
                        sum += (y[i] + yl) * (x[i] - xl);
                    }
                    xl = x[i];
                    yl = y[i];
                    if j >= nemax - 1 {
                        break;
                    }
                    if jbeta > 0 {
                        // beta(nbeta+1) is read past the end in the Fortran;
                        // either way the result is the same (see notes)
                        if jbeta <= nbeta && beta[jbeta] > BMAX {
                            break;
                        }
                    }
                    i -= 1;
                    if i > 1 {
                        st = 150;
                        continue;
                    }
                    jbeta += 1;
                    if jbeta <= nbeta {
                        st = 111;
                        continue;
                    }
                    if i == 1 {
                        st = 160;
                        continue;
                    }
                    break;
                }
                _ => unreachable!(),
            }
        }
        s[1] = sum;
        s[2] = j as f64;
    }
}

const NGRID: i32 = 118;
const EGRID: [f64; 118] = [
    1.0e-5, 1.78e-5, 2.5e-5, 3.5e-5, 5.0e-5, 7.0e-5, 1.0e-4, 1.26e-4, 1.6e-4, 2.0e-4, 0.000253, 0.000297, 0.000350,
    0.00042, 0.000506, 0.000615, 0.00075, 0.00087, 0.001012, 0.00123, 0.0015, 0.0018, 0.00203, 0.002277, 0.0026,
    0.003, 0.0035, 0.004048, 0.0045, 0.005, 0.0056, 0.006325, 0.0072, 0.0081, 0.009108, 0.01, 0.01063, 0.0115,
    0.012397, 0.0133, 0.01417, 0.015, 0.016192, 0.0182, 0.0199, 0.020493, 0.0215, 0.0228, 0.0253, 0.028, 0.030613,
    0.0338, 0.0365, 0.0395, 0.042757, 0.0465, 0.050, 0.056925, 0.0625, 0.069, 0.075, 0.081972, 0.09, 0.096, 0.1035,
    0.111573, 0.120, 0.128, 0.1355, 0.145728, 0.160, 0.172, 0.184437, 0.20, 0.2277, 0.2510392, 0.2705304, 0.2907501,
    0.3011332, 0.3206421, 0.3576813, 0.39, 0.4170351, 0.45, 0.5032575, 0.56, 0.625, 0.70, 0.78, 0.86, 0.95, 1.05,
    1.16, 1.28, 1.42, 1.55, 1.70, 1.855, 2.02, 2.18, 2.36, 2.59, 2.855, 3.12, 3.42, 3.75, 4.07, 4.46, 4.90, 5.35,
    5.85, 6.40, 7.00, 7.65, 8.40, 9.15, 9.85, 10.00,
];

/// Clamp a cosine to [-1,1] with the calcem message (tag "1cos=" or "2cos").
fn clamp_cos(c: &mut Ctx, v: &mut f64, tag: &str, enow: f64, ep: f64) {
    const UNITY: f64 = 1.0;
    if *v > UNITY {
        if *v > UNITY + 0.0005 {
            let strng = format!("{}{}, set to 1.  enow,e'={}{}", tag, ffmt::f(*v, 7, 4), ffmt::e1p(enow, 12, 5), ffmt::e1p(ep, 12, 5));
            c.mess("calcem", &strng, "");
        }
        *v = UNITY;
    }
    if *v < -UNITY {
        if *v < -(UNITY + 0.0005) {
            let strng = format!("{}{}, set to -1.  enow,e'={}{}", tag, ffmt::f(*v, 7, 4), ffmt::e1p(enow, 12, 5), ffmt::e1p(ep, 12, 5));
            c.mess("calcem", &strng, "");
        }
        *v = -UNITY;
    }
}


/// The inputs of the per-energy work of `calcem` that do not change in its
/// incident energy loop.
#[derive(Clone)]
struct CeIn {
    sp: Sab,
    nl: i32,
    nnl: i32,
    nbeta: i32,
    jmax: i32,
    tev: f64,
    tol: f64,
    temp: f64,
}

/// A growing 0-based array, written where `scr` would be written.
#[derive(Default)]
struct Grow(Vec<f64>);
impl std::ops::Index<i32> for Grow {
    type Output = f64;
    fn index(&self, i: i32) -> &f64 {
        &self.0[i as usize]
    }
}
impl std::ops::IndexMut<i32> for Grow {
    fn index_mut(&mut self, i: i32) -> &mut f64 {
        let i = i as usize;
        if i >= self.0.len() {
            self.0.resize(i + 1, 0.0);
        }
        &mut self.0[i]
    }
}

/// The result of one incident energy of the e-mu-ep loop of `calcem`.
#[derive(Default)]
struct EmuOut {
    /// false: not computed (a worker gave up); the main program computes it
    ok: bool,
    enow: f64,
    xsi: f64,
    ubar: f64,
    p2: f64,
    p3: f64,
    j: i32,
    /// the last value given to jnz, if any
    jnz: Option<i32>,
    /// the last value that sigl gave to nbin, if any
    nbin: Option<i32>,
    /// scr(7..): the values written to scr
    sc: Grow,
}

/// One incident energy ie of the e-mu-ep loop of `calcem` (loop 311 to
/// label 430), up to the jnz test. The energies do not depend on each
/// other except through jnz, which is returned for the caller to apply.
fn emu_energy(inp: &CeIn, c: &mut Ctx, ie: i32) -> EmuOut {
    const NLMAX: i32 = 65;
    const IMAX: i32 = 20;
    const TOLMIN: f64 = 5.0e-7;
    const HALF: f64 = 0.5e0;
    const THERM: f64 = 0.0253e0;
    const BREAK: f64 = 3000.0e0;
    const EM9: f64 = 1.0e-9;
    const ZERO: f64 = 0.0;
    const UUMIN: f64 = 0.00001e0;
    let egrid = |i: i32| EGRID[(i - 1) as usize];
    let sp = &inp.sp;
    let (nl, nnl, nbeta, jmax, tev, tol, temp) = (inp.nl, inp.nnl, inp.nbeta, inp.jmax, inp.tev, inp.tol, inp.temp);
    let mut x = [0.0f64; IMAX as usize + 2];
    let mut y = vec![0.0f64; (NLMAX * (IMAX + 2)) as usize];
    let yi = |il: i32, i: i32| ((il - 1) + NLMAX * (i - 1)) as usize;
    let mut yt = [0.0f64; NLMAX as usize + 2];
    let mut p = [0.0f64; 8];
    let mut sc = Grow::default();
    let mut nbin = i32::MIN;
    let mut jnz: Option<i32> = None;
    let mut xlast = 0.0f64;
    let mut ylast = 0.0f64;
    let mut ulast = 0.0f64;
    let mut u2last = 0.0f64;
    let mut u3last = 0.0f64;
    let mut xm = 0.0f64;
    let mut ym = 0.0f64;
    let mut xsi_e;
    let mut ubar_e;
    let mut p2_e;
    let mut p3_e;
    let mut enow;
    enow = egrid(ie);
    if temp > BREAK {
        let tone = THERM / BK;
        let elo = egrid(1);
        enow = elo * exp(log(enow / elo) * log((temp / tone) * egrid(NGRID) / elo) / log(egrid(NGRID) / elo));
    }
    enow = sigfig(enow, 8, 0);
                    xsi_e = 0.0;
    ubar_e = 0.0;
    p2_e = 0.0;
    p3_e = 0.0;
    let mut ep = 0.0f64;
    x[1] = ep;
    sp.sigl(c, &mut nbin, enow, ep, nnl, tev, &mut yt, tol);
    for il in 1..=nl {
        y[yi(il, 1)] = yt[il as usize];
    }
    let mut jbeta = -nbeta;
    if sp.lasym > 0 {
        jbeta = 1;
    }
    let mut j = 0i32;
    let mut iskip = 0;
    let mut i;
    // label 311: loop over beta values
    'beta: loop {
        x[2] = x[1];
        for il in 1..=nl {
            y[yi(il, 2)] = y[yi(il, 1)];
        }
        loop {
            // label 313
            if jbeta == 0 {
                jbeta = 1;
            }
            let bt = if sp.lat == 1 { sp.tevz } else { tev };
            if jbeta <= 0 {
                ep = enow - sp.beta[-jbeta] * bt;
                if ep == enow {
                    ep = sigfig(enow, 8, -1);
                } else {
                    ep = sigfig(ep, 8, 0);
                }
            } else {
                ep = enow + sp.beta[jbeta] * bt;
                if ep == enow {
                    ep = sigfig(enow, 8, 1);
                    iskip = 1;
                } else {
                    ep = sigfig(ep, 8, 0);
                }
            }
            if ep > x[2] {
                break;
            }
            jbeta += 1;
        }
        // label 316
        ep = sigfig(ep, 8, 0);
        x[1] = ep;
        sp.sigl(c, &mut nbin, enow, ep, nnl, tev, &mut yt, tol);
        for il in 1..=nl {
            y[yi(il, 1)] = yt[il as usize];
        }
        i = 2i32;
        loop {
            // label 330
            let iu_ = i as usize;
            let mut refine = false;
            'r360: {
                if i == IMAX {
                    break 'r360;
                }
                if iskip == 1 {
                    iskip = 0;
                    break 'r360;
                }
                if HALF * (y[yi(1, i - 1)] + y[yi(1, i)]) * (x[iu_ - 1] - x[iu_]) < TOLMIN {
                    break 'r360;
                }
                xm = HALF * (x[iu_ - 1] + x[iu_]);
                xm = sigfig(xm, 8, 0);
                if xm <= x[iu_] || xm >= x[iu_ - 1] {
                    break 'r360;
                }
                sp.sigl(c, &mut nbin, enow, xm, nnl, tev, &mut yt, tol);
                let mut uu = 0.0f64;
                let mut uum = 0.0f64;
                for k in 1..=nl {
                    c.terp1(x[iu_], y[yi(k, i)], x[iu_ - 1], y[yi(k, i - 1)], xm, &mut ym, 2);
                    if k > 1 {
                        uu += yt[k as usize];
                    }
                    if k > 1 {
                        uum += ym;
                    }
                    let test = tol * yt[k as usize].abs();
                    let mut test2 = test;
                    if k > 1 {
                        test2 = tol;
                    }
                    if (yt[k as usize] - ym).abs() > test2 {
                        refine = true;
                        break 'r360;
                    }
                }
                let test = 2.0 * tol * uu.abs() + UUMIN;
                if (uu - uum).abs() > test {
                    refine = true;
                }
            }
            if refine {
                // label 410
                i += 1;
                let iu_ = i as usize;
                x[iu_] = x[iu_ - 1];
                x[iu_ - 1] = xm;
                for il in 1..=nl {
                    y[yi(il, i)] = y[yi(il, i - 1)];
                    y[yi(il, i - 1)] = yt[il as usize];
                }
                continue;
            }
            // label 360
            j += 1;
            if j >= jmax {
                c.error("calcem", "storage exceeded.", " ");
            }
            if j > 1 {
                xsi_e += (x[iu_] - xlast) * (y[yi(1, i)] + ylast) * HALF;
            }
            if j > 1 {
                let mut uu = 0.0f64;
                let mut u2 = 0.0f64;
                let mut u3 = 0.0f64;
                let nll = 3;
                for il in 2..=nl {
                    legndr(y[yi(il, i)], &mut p, nll);
                    uu += p[1];
                    u2 += p[2];
                    u3 += p[3];
                }
                uu /= (nl - 1) as f64;
                uu *= y[yi(1, i)];
                u2 /= (nl - 1) as f64;
                u3 /= (nl - 1) as f64;
                u2 *= y[yi(1, i)];
                u3 *= y[yi(1, i)];
                ubar_e += HALF * (x[iu_] - xlast) * (uu + ulast);
                p2_e += HALF * (x[iu_] - xlast) * (u2 + u2last);
                p3_e += HALF * (x[iu_] - xlast) * (u3 + u3last);
            }
            if !(j != 3 || xsi_e >= TOLMIN) {
                j = 2;
            }
            // label 380
            let jscr = 7 + (j - 1) * (nl + 1);
            sc[jscr] = x[iu_];
            if y[yi(1, i)] >= EM9 {
                sc[1 + jscr] = sigfig(y[yi(1, i)], 9, 0);
            } else {
                sc[1 + jscr] = sigfig(y[yi(1, i)], 8, 0);
            }
            for il in 2..=nl {
                let mut v = sigfig(y[yi(il, i)], 9, 0);
                clamp_cos(c, &mut v, "1cos=", enow, sc[jscr]);
                sc[il + jscr] = v;
            }
            xlast = x[iu_];
            ylast = y[yi(1, i)];
            if ylast != ZERO {
                jnz = Some(j);
            }
            ulast = 0.0;
            u2last = 0.0;
            u3last = 0.0;
            let nll = 3;
            for il in 2..=nl {
                legndr(y[yi(il, i)], &mut p, nll);
                ulast += p[1];
                u2last += p[2];
                u3last += p[3];
            }
            ulast = ulast * y[yi(1, i)] / (nl - 1) as f64;
            u2last = u2last * y[yi(1, i)] / (nl - 1) as f64;
            u3last = u3last * y[yi(1, i)] / (nl - 1) as f64;
            i -= 1;
            if i >= 2 {
                continue;
            }
            jbeta += 1;
            if jbeta <= nbeta {
                continue 'beta;
            }
            for il in 1..=nl {
                y[yi(il, i)] = 0.0;
            }
            break 'beta;
        }
    }

    //--add the last point (label 430)
    let iu_ = i as usize;
    j += 1;
    xsi_e += (x[iu_] - xlast) * (y[yi(1, i)] + ylast) * HALF;
    let uu = 0.0f64;
    let u2 = 0.0f64;
    let u3 = 0.0f64;
    ubar_e += HALF * (x[iu_] - xlast) * (uu + ulast);
    p2_e += HALF * (x[iu_] - xlast) * (u2 + u2last);
    p3_e += HALF * (x[iu_] - xlast) * (u3 + u3last);
    xsi_e = sigfig(xsi_e, 9, 0);
    sc[7 + (nl + 1) * (j - 1)] = x[iu_];
    let jscr = 7 + (j - 1) * (nl + 1);
    if y[yi(1, i)] >= EM9 {
        sc[1 + jscr] = sigfig(y[yi(1, i)], 9, 0);
    } else {
        sc[1 + jscr] = sigfig(y[yi(1, i)], 8, 0);
    }
    for il in 2..=nl {
        let mut v = sigfig(y[yi(il, i)], 9, 0);
        clamp_cos(c, &mut v, "2cos", enow, sc[jscr]);
        sc[il + jscr] = v;
    }
    if y[yi(1, 1)] != ZERO {
        jnz = Some(j);
    }
    EmuOut {
        ok: true,
        enow,
        xsi: xsi_e,
        ubar: ubar_e,
        p2: p2_e,
        p3: p3_e,
        j,
        jnz,
        nbin: if nbin == i32::MIN { None } else { Some(nbin) },
        sc,
    }
}

/// The result of one incident energy of the e-ep-mu loop of `calcem`.
#[derive(Default)]
struct EpmOut {
    /// false: not computed (a worker gave up); the main program computes it
    ok: bool,
    enow: f64,
    xsi: f64,
    ubar: f64,
    sum: f64,
    nmu: i32,
    /// uj(1..nmu), sj(1..nmu) (index 0 unused)
    uj: Vec<f64>,
    sj: Vec<f64>,
    /// for each cosine uj(il): yu(0..2*nep+2) from sigu, nep=nint(yu(2))
    yus: Vec<Vec<f64>>,
}

/// One incident energy ie of the e-ep-mu loop of `calcem` (labels 530 to
/// 580), and the sigu result for each cosine that the loop writes later.
fn epm_energy(inp: &CeIn, c: &mut Ctx, ie: i32) -> EpmOut {
    const NEMAX: i32 = 5000;
    const MUMAX: i32 = 300;
    const IMAX: i32 = 20;
    const TOLMIN: f64 = 5.0e-7;
    const HALF: f64 = 0.5e0;
    const BREAK: f64 = 3000.0e0;
    let egrid = |i: i32| EGRID[(i - 1) as usize];
    let sp = &inp.sp;
    let (tev, tol, temp) = (inp.tev, inp.tol, inp.temp);
    let mut x = [0.0f64; IMAX as usize + 2];
    let mut yy = [0.0f64; IMAX as usize + 2];
    let mut yu = vec![0.0f64; 2 * NEMAX as usize + 8];
    let mut uj = vec![0.0f64; MUMAX as usize + 2];
    let mut sj = vec![0.0f64; MUMAX as usize + 2];
    let xsi_e;
    let mut ubar_e;
    let mut enow;
    enow = egrid(ie);
    if ie > 1 && temp > BREAK {
        enow = enow * temp / BREAK;
    }
    enow = sigfig(enow, 8, 0);
                    let mut j = 0i32;
    let mut sum = 0.0f64;
    let mut u = -1.0f64;
    x[2] = u;
    sp.sigu(enow, u, tev, &mut yu, NEMAX, tol);
    yy[2] = yu[1];
    let mut xl = x[2];
    let mut yl = yy[2];
    u = 1.0;
    x[1] = u;
    sp.sigu(enow, u, tev, &mut yu, NEMAX, tol);
    yy[1] = yu[1];
    let mut i = 2usize;
    loop {
        // label 530
        let mut split = false;
        let mut xm = 0.0;
        'l560: {
            if i == IMAX as usize {
                break 'l560;
            }
            xm = HALF * (x[i - 1] + x[i]);
            xm = sigfig(xm, 7, 0);
            if xm <= x[i] || xm >= x[i - 1] {
                break 'l560;
            }
            sp.sigu(enow, xm, tev, &mut yu, NEMAX, tol);
            if x[i - 1] - x[i] > 0.25 {
                split = true;
                break 'l560;
            }
            let ym = yy[i] + (xm - x[i]) * (yy[i - 1] - yy[i]) / (x[i - 1] - x[i]);
            if (yu[1] - ym).abs() > 2.0 * tol * ym + TOLMIN {
                split = true;
            }
        }
        if split {
            // label 575
            i += 1;
            x[i] = x[i - 1];
            x[i - 1] = xm;
            yy[i] = yy[i - 1];
            yy[i - 1] = yu[1];
            continue;
        }
        // label 560
        j += 1;
        if j > MUMAX - 1 {
            c.error("calcem", "too many angles", "see mumax");
        }
        uj[j as usize] = x[i];
        sj[j as usize] = yy[i];
        if j > 1 {
            sum += HALF * (yy[i] + yl) * (x[i] - xl);
            xl = x[i];
            yl = yy[i];
        }
        i -= 1;
        if i >= 2 {
            continue;
        }
        break;
    }
    // label 580
    j += 1;
    uj[j as usize] = x[1];
    sj[j as usize] = yy[1];
    let nmu = j;
    ubar_e = 0.0;
    sum += HALF * (yy[1] + yl) * (x[1] - xl);
    xsi_e = sum / 2.0;
    for i in 2..=nmu as usize {
        ubar_e += HALF * (uj[i] - uj[i - 1]) * (sj[i] + sj[i - 1]) * (uj[i] + uj[i - 1]);
    }
    ubar_e = HALF * ubar_e / sum;
    let mut yus = Vec::with_capacity(nmu.max(0) as usize);
    for il in 1..=nmu {
        let u = uj[il as usize];
        sp.sigu(enow, u, tev, &mut yu, NEMAX, tol);
        let nep = nint(yu[2]).max(0) as usize;
        yus.push(yu[..(2 * nep + 3).min(yu.len())].to_vec());
    }
    uj.truncate(nmu.max(0) as usize + 1);
    sj.truncate(nmu.max(0) as usize + 1);
    EpmOut { ok: true, enow, xsi: xsi_e, ubar: ubar_e, sum, nmu, uj, sj, yus }
}

/// Start workers that compute the energies 1..nne of a `calcem` loop ahead
/// of the main program with f (None with one thread). A worker uses a
/// detached Ctx: if the work tried to write anything or stopped, its result
/// is the default (not ok) and the main program computes that energy itself,
/// so the output never depends on the threads.
fn ce_workers<R: Send + Default + 'static>(
    inp: &std::sync::Arc<CeIn>, nne: i32, f: fn(&CeIn, &mut Ctx, i32) -> R,
) -> Option<crate::prefetch::RunAhead<R>> {
    let nthr = crate::util::threads();
    if nthr <= 1 || nne < 2 {
        return None;
    }
    let workers: Vec<(std::sync::Arc<CeIn>, Ctx)> =
        (0..nthr).map(|_| (inp.clone(), Ctx::new(crate::fio::Io::detached()))).collect();
    Some(crate::prefetch::RunAhead::start(workers, 1, nne as usize + 1, 4 * nthr, move |w, q| {
        w.1.io.touched = false;
        let r = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| f(&w.0, &mut w.1, q as i32)));
        match r {
            Ok(o) if !w.1.io.touched => o,
            _ => R::default(),
        }
    }))
}

impl Thermr {
    /// `calcem`: incoherent inelastic cross sections and distributions.
    fn calcem(
        &mut self, c: &mut Ctx, temp: f64, itemp: i32, iold: &mut i32, inew: &mut i32, ne: &mut i32, nex: &mut i32,
        mtref: i32,
    ) {
        const NLMAX: i32 = 65;
        const NEMAX: i32 = 5000;
        const MUMAX: i32 = 300;
        const IMAX: i32 = 20;
        const EPS: f64 = 1.0e-4;
        const THERM: f64 = 0.0253e0;
        const ZERO: f64 = 0.0;
        const TENTH: f64 = 0.1e0;
        const UP: f64 = 1.1e0;
        const DN: f64 = 0.9e0;
        const YUMIN: f64 = 2.0e-7;
        const NLPMX: i32 = 10;
        let egrid = |i: i32| EGRID[(i - 1) as usize];
        let (mut nb, mut nw) = (0i32, 0i32);
        let mut ex = [0.0f64; IMAX as usize + 2];
        let mut yu = vec![0.0f64; 2 * NEMAX as usize + 8];
        let mut ubar = [0.0f64; NGRID as usize + 2];
        let mut p2 = [0.0f64; NGRID as usize + 2];
        let mut p3 = [0.0f64; NGRID as usize + 2];
        let mut uj = [0.0f64; MUMAX as usize + 2];
        let mut sj = [0.0f64; MUMAX as usize + 2];

        //--initialize
        self.tevz = THERM;
        if itemp == 1 {
            self.ce_sabmin = exp(SABFLG);
            self.ce_nwtab = 8;
            self.ce_nl = self.nbin + 1;
            if self.ce_nl > NLMAX {
                c.error("calcem", "nl too large for binning.", " ");
            }
            self.ce_nlp = self.ce_nl + 1;
            if self.ce_nlp > NLPMX {
                self.ce_nlp = NLPMX;
            }
            self.ce_nlp1 = self.ce_nlp + 1;
            self.ce_nl1 = self.ce_nl + 1;
            self.ce_nnl = -self.ce_nl;
            self.ce_jmax = (NWSCR - 50) / self.ce_nl;
            self.ce_nne = 0;
            for i in 1..=NGRID {
                if egrid(i) > self.emax && self.ce_nne == 0 {
                    self.ce_nne = i;
                }
            }
            if self.ce_nne == 0 {
                self.ce_nne = NGRID;
            }
            self.esi = Some(RV::new(self.ce_nne as usize + 1 + 8));
            self.xsi = Some(RV::new(self.ce_nne as usize + 1 + 8));
            self.ce_nlt = 5;
            if self.ce_nne <= self.ce_nlt {
                self.ce_nlt = self.ce_nne - 1;
            }
        }
        let (nl, nlt, nlp, nlp1, nl1, nnl, jmax, nne) = (
            self.ce_nl,
            self.ce_nlt,
            self.ce_nlp,
            self.ce_nlp1,
            self.ce_nl1,
            self.ce_nnl,
            self.ce_jmax,
            self.ce_nne,
        );
        c.repoz(self.nscr);
        c.nsc = 0;
        self.ncds = 0;
        let scr = &mut self.scr;
        let nendf = self.nendf;
        let t;
        let (nalpha, nbeta);
        let mut alpha;
        let mut beta;
        let mut sab;
        if self.iinc != 1 {
            //--read in s(alpha,beta,t) data
            let matd = self.matde;
            let mut itprnt = 0;
            c.findf(matd, 7, 4, nendf);
            c.contio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
            self.lasym = c.n1h;
            if self.lasym > 1 {
                c.error("calcem", "isabt=1 pendf tape found", "thermr cannot process this format");
            }
            self.lat = c.l2h;
            c.listio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
            let ilog = c.l1h;
            let ni = c.n1h;
            let smz = scr[7] / self.natom as f64;
            self.az = scr[9];
            let q = (self.az + 1.0) / self.az;
            self.sb = smz * (q * q);
            if c.iverf >= 6 {
                self.sb2 = 0.0;
                self.az2 = 0.0;
                if ni > 6 {
                    if scr[13] == ZERO {
                        self.az2 = scr[15];
                    }
                    if self.az2 != ZERO {
                        let q = (self.az2 + 1.0) / self.az2;
                        self.sb2 = scr[14] * (q * q);
                    }
                    self.sb2 = self.sb2 / scr[18] / self.natom as f64;
                    if ni > 12 {
                        if scr[19] == ZERO {
                            c.error("calcem", "only 2 sct atoms allowed.", " ");
                        }
                    }
                }
            }
            c.tab2io(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
            nbeta = nint(scr[6]);
            let loc = 1 + nw;
            c.tab1io(nendf, 0, 0, &mut scr.at(loc), &mut nb, &mut nw);
            let t1 = c.c1h;
            let mut b = c.c2h;
            let lt = c.l1h;
            let mut it = 0;
            nalpha = c.n2h;
            let mut l = loc + nw;
            while nb != 0 {
                c.moreio(nendf, 0, 0, &mut scr.at(l), &mut nb, &mut nw);
                l += nw;
            }
            alpha = RV::new(nalpha.max(0) as usize + 4);
            l = loc + 4 + 2 * c.n1h;
            for i in 1..=nalpha {
                alpha[i] = scr[2 * i + l];
            }
            beta = RV::new(nbeta.max(0) as usize + 4);
            sab = vec![0.0f64; (nalpha.max(1) * (nbeta.max(1) + 4)) as usize + 8];
            let mut itrunc = 0;
            let mut ib = 1i32;
            beta[ib] = b;
            let mut tt;
            loop {
                // label 120
                tt = scr[loc];
                while (tt - temp).abs() > temp / 500.0 {
                    if tt > temp {
                        c.error("calcem", "desired temperature not found.", " ");
                    }
                    it += 1;
                    if it > lt {
                        c.error("calcem", "desired temperature not found.", " ");
                    }
                    c.listio(nendf, 0, 0, &mut scr.at(loc), &mut nb, &mut nw);
                    let mut l = loc + nw;
                    while nb != 0 {
                        c.moreio(nendf, 0, 0, &mut scr.at(l), &mut nb, &mut nw);
                        l += nw;
                    }
                    tt = scr[loc];
                }
                if (tt - temp).abs() > EPS && itprnt == 0 {
                    itprnt = 1;
                    let diff = (tt - temp).abs();
                    let m = format!("\n difference between temperatures desired and found is {}", ffmt::e1p(diff, 10, 2));
                    c.io.so(&m);
                    c.io.se(&m);
                }
                let mut l = loc + 6;
                if tt == t1 {
                    l = l + 2 * nint(scr[loc + 4]) + 1;
                }
                for ia in 1..=nalpha {
                    let k = ((ia - 1) + nalpha * (ib - 1)) as usize;
                    if ilog == 0 {
                        if scr[l] > self.ce_sabmin {
                            sab[k] = log(scr[l]);
                        }
                        if scr[l] <= self.ce_sabmin {
                            if scr[l] > ZERO {
                                itrunc = 1;
                            }
                            sab[k] = SABFLG;
                        }
                    } else {
                        sab[k] = scr[l];
                    }
                    if tt == t1 {
                        l += 1;
                    }
                    l += 1;
                }
                while it < lt {
                    it += 1;
                    c.listio(nendf, 0, 0, &mut scr.at(loc), &mut nb, &mut nw);
                    let mut l = loc;
                    while nb != 0 {
                        l += nw;
                        c.moreio(nendf, 0, 0, &mut scr.at(l), &mut nb, &mut nw);
                    }
                }
                if ib == nbeta {
                    break;
                }
                ib += 1;
                c.tab1io(nendf, 0, 0, &mut scr.at(loc), &mut nb, &mut nw);
                let mut l = loc;
                while nb != 0 {
                    l += nw;
                    c.moreio(nendf, 0, 0, &mut scr.at(l), &mut nb, &mut nw);
                }
                b = scr[loc + 1];
                beta[ib] = b;
                it = 0;
            }
            t = tt;
            // label 180
            if itrunc == 1 {
                let tmax = 2.0 * (-SABFLG - 5.0) * BK * temp;
                let m = format!(
                    "\n ***warning***sab contains numbers that are too small to represent on this machine.\n results may be bad for transfers larger than {} ev.",
                    ffmt::f(tmax, 8, 3)
                );
                c.io.so(&m);
                c.io.se(&m);
            }
            let mut tmax = beta[nbeta] * BK * temp;
            if self.lat == 1 {
                tmax = tmax * self.tevz / (BK * temp);
            }
            let test = 5.0;
            if tmax < test {
                let m = format!(
                    "\n ***warning***maximum value of beta limits the allowed energy transfer\n the sct approx. will be used for transfers larger than {} ev.",
                    ffmt::f(tmax, 6, 3)
                );
                c.io.so(&m);
                c.io.se(&m);
            }
            if c.iverf >= 6 {
                //--read effective temperatures
                c.tab1io(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                let nr = nint(scr[5]);
                let np = nint(scr[6]);
                let tempt = temp;
                let mut tnxt = 0.0;
                let mut idis = 0;
                if np == 1 {
                    let tt1 = scr[7 + 2 * nr];
                    if (tempt - tt1).abs() > tempt / 10.0 {
                        c.error("calcem", "bad temperature for teff", " ");
                    }
                    self.teff = scr[8 + 2 * nr];
                } else {
                    let mut ip = 2;
                    let mut ir = 1;
                    let tt1 = scr[7 + 2 * nr];
                    let ttn = scr[5 + 2 * nr + 2 * np];
                    if tempt < DN * tt1 || tempt > UP * ttn {
                        c.error("calcem", "bad temperature for teff", " ");
                    }
                    if tt1 > tempt {
                        scr[7 + 2 * nr] = tempt;
                    }
                    if ttn < tempt {
                        scr[5 + 2 * nr + 2 * np] = tempt;
                    }
                    c.terpa(&mut self.teff, tempt, &mut tnxt, &mut idis, &scr.all(), &mut ip, &mut ir);
                }
                if self.sb2 != ZERO {
                    c.tab1io(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                    let nr = nint(scr[5]);
                    let np = nint(scr[6]);
                    if np <= 1 {
                        let tt1 = scr[7 + 2 * nr];
                        if (tempt - tt1).abs() > tempt / 10.0 {
                            c.error("calcem", "bad temperature for teff2", " ");
                        }
                        self.teff2 = scr[8 + 2 * nr];
                    } else {
                        let mut ip = 2;
                        let mut ir = 1;
                        let tt1 = scr[7 + 2 * nr];
                        let ttn = scr[5 + 2 * nr + 2 * np];
                        if tempt < DN * tt1 || tempt > UP * ttn {
                            c.error("calcem", "bad temperature for teff2", " ");
                        }
                        if tt1 > tempt {
                            scr[7 + 2 * nr] = tempt;
                        }
                        if ttn < tempt {
                            scr[5 + 2 * nr + 2 * np] = tempt;
                        }
                        c.terpa(&mut self.teff2, tempt, &mut tnxt, &mut idis, &scr.all(), &mut ip, &mut ir);
                    }
                }
            }
        } else {
            //--set up for free-gas scattering (label 200)
            nalpha = 1;
            nbeta = 45;
            alpha = RV::new(nalpha as usize + 4);
            beta = RV::new(nbeta as usize + 4);
            sab = vec![0.0f64; (nalpha * (nbeta + 4)) as usize + 8];
            const BV: [f64; 45] = [
                0.0, TENTH, 2.0, 4.0, 6.0, 8.0, 10.0, 15.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 55.0, 60.0, 65.0, 70.0,
                75.0, 80.0, 100.0, 120.0, 140.0, 160.0, 180.0, 200.0, 250.0, 300.0, 350.0, 400.0, 500.0, 600.0, 700.0,
                800.0, 900.0, 1000.0, 1250.0, 1500.0, 1750.0, 2000.0, 2250.0, 2500.0, 2750.0, 3000.0, 3500.0,
            ];
            for i in 1..=45 {
                beta[i] = BV[(i - 1) as usize];
            }
            t = temp;
            let smz = 1.0;
            let q = (self.az + 1.0) / self.az;
            self.sb = smz * (q * q);
            self.lat = 0;
        }

        // label 300
        let tev = t * BK;
        let alpha_up = ascending(&alpha, nalpha);
        let beta_up = ascending(&beta, nbeta);
        let mut sp = Sab {
            alpha_up,
            beta_up,
            nalpha,
            nbeta,
            alpha,
            beta,
            sab,
            az: self.az,
            az2: self.az2,
            tevz: self.tevz,
            iinc: self.iinc,
            lat: self.lat,
            lasym: self.lasym,
            cliq: 0.0,
            teff: 0.0,
            teff2: 0.0,
            sb: self.sb,
            sb2: self.sb2,
            thr6: c.thr6,
        };
        let esi = self.esi.as_mut().unwrap();
        let xsi = self.xsi.as_mut().unwrap();
        let tol = self.tol;
        let mut enow;
        if self.iform != 1 {
            //--set up for e-mu-ep format
            let ltt = 5;
            c.math = self.matdp;
            c.mfh = 6;
            c.mth = mtref;
            self.teff *= BK;
            self.teff2 *= BK;
            scr[1] = self.za;
            scr[2] = self.awr;
            scr[3] = 0.0;
            scr[4] = ltt as f64;
            scr[5] = 1.0;
            scr[6] = 0.0;
            nw = 6;
            c.contio(0, 0, self.nscr, &mut scr.all(), &mut nb, &mut nw);
            self.ncds += 1;
            scr[1] = 1.0;
            scr[2] = 1.0;
            scr[3] = -1.0;
            scr[4] = 1.0;
            scr[5] = 1.0;
            scr[6] = 2.0;
            scr[7] = 2.0;
            scr[8] = 2.0;
            scr[9] = 1.0e-5;
            scr[10] = 1.0;
            scr[11] = self.emax;
            scr[12] = 1.0;
            nw = 12;
            c.tab1io(0, 0, self.nscr, &mut scr.all(), &mut nb, &mut nw);
            self.ncds += 2;
            scr[1] = temp;
            scr[2] = 0.0;
            scr[3] = 3.0;
            scr[4] = 1.0;
            scr[5] = 1.0;
            scr[6] = nne as f64;
            scr[7] = nne as f64;
            scr[8] = 2.0;
            nw = 8;
            c.tab2io(0, 0, self.nscr, &mut scr.all(), &mut nb, &mut self.ce_nwtab);
            self.ncds += 2;
            self.cliq = 0.0;
            if self.iinc != 1 && sp.s(1, 1) > sp.s(2, 1) {
                self.cliq = (sp.s(1, 1) - sp.s(1, 2)) * sp.alpha[1] / (sp.beta[2] * sp.beta[2]);
            }
            sp.cliq = self.cliq;
            sp.teff = self.teff;
            sp.teff2 = self.teff2;
            if self.iinc != 1 && self.iinc != 2 {
                c.error("sig", "illegal option.", " ");
            }

            //--loop over given incident energy grid.
            let inp = std::sync::Arc::new(CeIn { sp: sp.clone(), nl, nnl, nbeta, jmax, tev, tol, temp });
            let mut ra = ce_workers(&inp, nne, emu_energy);
            let mut jnz = 0i32;
            for ie in 1..=nne {
                let iu = ie as usize;
                let o = match ra.as_mut().and_then(|r| r.take(ie as usize)) {
                    Some(o) if o.ok => o,
                    _ => emu_energy(&inp, c, ie),
                };
                enow = o.enow;
                esi[ie] = enow;
                xsi[ie] = o.xsi;
                ubar[iu] = o.ubar;
                p2[iu] = o.p2;
                p3[iu] = o.p3;
                if let Some(nb) = o.nbin {
                    self.nbin = nb;
                }
                for k in 7..o.sc.0.len() {
                    scr[k as i32] = o.sc.0[k];
                }
                if let Some(z) = o.jnz {
                    jnz = z;
                }
                let mut j = o.j;
                if jnz < j {
                    j = jnz + 1;
                }
                if self.iprint == 2 {
                    ubar[iu] /= xsi[ie];
                    p2[iu] /= xsi[ie];
                    p3[iu] /= xsi[ie];
                    c.io.so(&format!(
                        "\n incident energy ={}     cross section ={}     mubar,p2,p3 ={}{}{}",
                        ffmt::e1p(enow, 13, 6),
                        ffmt::e1p(xsi[ie], 13, 6),
                        ffmt::e1p(ubar[iu], 12, 4),
                        ffmt::e1p(p2[iu], 12, 4),
                        ffmt::e1p(p3[iu], 12, 4)
                    ));
                    c.io.so("\n     exit energy           pdf       cosines");
                    c.io.so(&format!("   ---------------     -----------  {}", "-".repeat(88)));
                    let mut ll = 6;
                    for _jj in 1..=j {
                        let mut line = String::new();
                        for ii in 1..=nlp {
                            let v = scr[ll + ii];
                            match ii {
                                1 => line.push_str(&format!("  {}", ffmt::e1p(v, 15, 8))),
                                2 => line.push_str(&format!("     {}", ffmt::e1p(v, 12, 5))),
                                _ => line.push_str(&ffmt::f(v, 11, 6)),
                            }
                        }
                        c.io.so(&line);
                        if nl1 > nlp {
                            let mut line = " ".repeat(34);
                            let mut cnt = 0;
                            for ii in nlp1..=nl1 {
                                if cnt == 8 {
                                    c.io.so(&line);
                                    line = " ".repeat(34);
                                    cnt = 0;
                                }
                                line.push_str(&ffmt::f(scr[ll + ii], 11, 6));
                                cnt += 1;
                            }
                            c.io.so(&line);
                        }
                        ll += nl1;
                    }
                }
                scr[1] = 0.0;
                scr[2] = enow;
                scr[3] = 0.0;
                scr[4] = 0.0;
                scr[5] = ((nl + 1) * j) as f64;
                scr[6] = (nl + 1) as f64;
                self.ncds += 1 + (j * (nl + 1) + 5) / 6;
                c.listio(0, 0, self.nscr, &mut scr.all(), &mut nb, &mut nw);
                let mut loc = 1;
                while nb != 0 {
                    loc += nw;
                    c.moreio(0, 0, self.nscr, &mut scr.at(loc), &mut nb, &mut nw);
                }
            }
        } else {
            //--set up for e-ep-mu format (label 510)
            let ltt = 4;
            c.math = self.matdp;
            c.mfh = 6;
            c.mth = mtref;
            self.teff *= BK;
            self.teff2 *= BK;
            scr[1] = self.za;
            scr[2] = self.awr;
            scr[3] = 0.0;
            scr[4] = ltt as f64;
            scr[5] = 1.0;
            scr[6] = 0.0;
            c.contio(0, 0, self.nscr, &mut scr.all(), &mut nb, &mut nw);
            self.ncds += 1;
            scr[1] = 1.0;
            scr[2] = 1.0;
            scr[3] = 0.0;
            scr[4] = 7.0;
            scr[5] = 1.0;
            scr[6] = 2.0;
            scr[7] = 2.0;
            scr[8] = 2.0;
            scr[9] = 1.0e-5;
            scr[10] = 1.0;
            scr[11] = self.emax;
            scr[12] = 1.0;
            nw = 12;
            c.tab1io(0, 0, self.nscr, &mut scr.all(), &mut nb, &mut nw);
            self.ncds += 2;
            scr[1] = 0.0;
            scr[2] = 0.0;
            scr[3] = 0.0;
            scr[4] = 0.0;
            scr[5] = 1.0;
            scr[6] = nne as f64;
            scr[7] = nne as f64;
            scr[8] = 2.0;
            nw = 8;
            c.tab2io(0, 0, self.nscr, &mut scr.all(), &mut nb, &mut nw);
            self.ncds += 2;
            self.cliq = 0.0;
            if self.iinc != 1 && sp.s(1, 1) > sp.s(2, 1) {
                self.cliq = (sp.s(1, 1) - sp.s(1, 2)) * sp.alpha[1] / (sp.beta[2] * sp.beta[2]);
            }
            sp.cliq = self.cliq;
            sp.teff = self.teff;
            sp.teff2 = self.teff2;
            if self.iinc != 1 && self.iinc != 2 {
                c.error("sig", "illegal option.", " ");
            }
            let inp = std::sync::Arc::new(CeIn { sp: sp.clone(), nl, nnl, nbeta, jmax, tev, tol, temp });
            let mut ra = ce_workers(&inp, nne, epm_energy);
            for ie in 1..=nne {
                let iu = ie as usize;
                let o = match ra.as_mut().and_then(|r| r.take(ie as usize)) {
                    Some(o) if o.ok => o,
                    _ => epm_energy(&inp, c, ie),
                };
                enow = o.enow;
                esi[ie] = enow;
                xsi[ie] = o.xsi;
                ubar[iu] = o.ubar;
                let sum = o.sum;
                let nmu = o.nmu;
                uj[..o.uj.len()].copy_from_slice(&o.uj);
                sj[..o.sj.len()].copy_from_slice(&o.sj);
                if self.iprint == 2 {
                    c.io.so(&format!(
                        "\n{} enow {}   xsec {}   mubar  {}",
                        ffmt::i(ie, 5),
                        ffmt::e1p(enow, 13, 6),
                        ffmt::e1p(xsi[ie], 13, 6),
                        ffmt::e1p(ubar[iu], 13, 6)
                    ));
                    c.io.so(&format!("      num of mu {}", ffmt::i(nmu, 5)));
                    c.io.so("\n            mu            theta      dsigma/dmu");
                    for i in 1..=nmu {
                        let iu2 = i as usize;
                        c.io.so(&format!(
                            "{} {} {} {}",
                            ffmt::i(i, 5),
                            ffmt::f(uj[iu2], 15, 8),
                            ffmt::f(crate::fmath::acos(uj[iu2]) * 180.0 / PI, 12, 4),
                            ffmt::e1p(sj[iu2] / 2.0, 14, 7)
                        ));
                    }
                }
                c.mth = mtref;
                scr[1] = 0.0;
                scr[2] = enow;
                scr[3] = 0.0;
                scr[4] = 0.0;
                scr[5] = 1.0;
                scr[6] = nmu as f64;
                scr[7] = nmu as f64;
                scr[8] = 2.0;
                nw = 8;
                c.tab2io(0, 0, self.nscr, &mut scr.all(), &mut nb, &mut nw);
                self.ncds += 1;
                for il in 1..=nmu {
                    let u = uj[il as usize];
                    let yv = &o.yus[(il - 1) as usize];
                    yu[..yv.len()].copy_from_slice(yv);
                    let mut nep = nint(yu[2]);
                    let mut j = 0i32;
                    for i in 1..=nep {
                        j = nep - i;
                        if yu[(2 * (nep - i) + 4) as usize] / sum > YUMIN {
                            break;
                        }
                    }
                    nep = j;
                    if self.iprint == 2 {
                        c.io.so(&format!("\n mu = {}", ffmt::f(u, 15, 8)));
                        c.io.so(&format!(" (e-prime, pdf);  num of e-prime {}", ffmt::i(nep, 5)));
                        c.io.so("");
                        let mut line = String::new();
                        let mut cnt = 0;
                        for i in 1..=nep {
                            if cnt == 3 {
                                c.io.so(&line);
                                line.clear();
                                cnt = 0;
                            }
                            if cnt > 0 {
                                line.push(' ');
                            }
                            line.push_str(&format!(
                                " {} {}",
                                ffmt::e1p(yu[(2 * i + 1) as usize], 14, 7),
                                ffmt::e1p(yu[(2 * i + 2) as usize] / sum, 14, 7)
                            ));
                            cnt += 1;
                        }
                        if nep > 0 {
                            c.io.so(&line);
                        }
                    }
                    scr[1] = 0.0;
                    scr[2] = u;
                    scr[3] = 0.0;
                    scr[4] = 0.0;
                    scr[5] = 1.0;
                    scr[6] = nep as f64;
                    scr[7] = nep as f64;
                    scr[8] = 2.0;
                    let mut k = 8;
                    let mut istart = 1;
                    loop {
                        // label 595
                        let mut iend = nep;
                        if iend - istart >= c.npage / 2 {
                            iend = istart + c.npage / 2 - 1;
                        }
                        j = k - 1;
                        let mut ib = istart - 1;
                        loop {
                            j += 2;
                            ib += 1;
                            scr[j] = yu[(1 + 2 * ib) as usize];
                            scr[j + 1] = yu[(2 + 2 * ib) as usize] * 2.0 / sum;
                            if ib >= iend {
                                break;
                            }
                        }
                        nw = j + 1;
                        if k != 0 {
                            k = 0;
                            c.tab1io(0, 0, self.nscr, &mut scr.all(), &mut nb, &mut nw);
                        } else {
                            c.moreio(0, 0, self.nscr, &mut scr.all(), &mut nb, &mut nw);
                        }
                        if nb == 0 {
                            break;
                        }
                        istart = iend + 1;
                    }
                    self.ncds += 1 + (j * (nep + 1) + 5) / 6;
                }
            }
        }

        //--add the inelastic cross section to the scratch file (label 610)
        for k in 1..=20 {
            ex[k] = 0.0;
        }
        let mut nexn = *nex;
        if *nex == 2 {
            nexn = 3;
        }
        let mut ie = 0;
        let mut ien = 1;
        while ien > 0 {
            ie += 1;
            c.finda(ie, &mut ex[1..], *nex, *iold, &mut self.bufo.0, NBUF);
            if self.iinc == 1 {
                ex[3] = ex[2];
            } else {
                enow = ex[1];
                let mut xs = terp(&esi.0, &xsi.0, nne, enow, nlt);
                if ie == *ne {
                    xs = 0.0;
                }
                ex[3] = xs;
            }
            ien = ie;
            if ie == *ne {
                ien = -ien;
            }
            c.loada(ien, &ex[1..], nexn, *inew, &mut self.bufn.0, NBUF);
        }
        esi[nne + 1] = 0.0;
        std::mem::swap(iold, inew);
        *nex = nexn;
        c.asend(0, self.nscr);
    }
}

impl Thermr {
    /// `tpend`: write the new MF1 directory, MF3 and MF6 sections on the
    /// output PENDF tape for this temperature.
    fn tpend(&mut self, c: &mut Ctx, iold: i32, itemp: i32, ne: i32, nex: i32, icoh: i32, icopy: &mut i32, mtref: i32) {
        const EM9: f64 = 1.0e-9;
        const ETOP: f64 = 20.0e6;
        const ZERO: f64 = 0.0;
        let (mut nb, mut nw) = (0i32, 0i32);
        let mut ex = [0.0f64; 21 + 8];
        let mut sav = RV::new(6 + 8);
        let (nin, nout, nscr, nscr2) = (self.nin, self.nout, self.nscr, self.nscr2);
        let scr = &mut self.scr;

        let temp = self.tempr[itemp];
        let matd = self.matdp;
        c.repoz(nscr);
        c.repoz(nscr2);

        //--copy any temperatures that come before this one
        let mut jcopy = 0;
        if *icopy > 0 {
            while jcopy < *icopy {
                c.contio(nscr2, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                c.tomend(nscr2, nout, 0, &mut scr.all());
                jcopy += 1;
            }
            *icopy = 0;
        }

        //--read the old directory
        c.findf(matd, 1, 451, nin);
        let mut indx = 1i32;
        c.contio(nin, 0, 0, &mut scr.at(indx), &mut nb, &mut nw);
        let n6 = c.n2h;
        c.contio(nin, 0, 0, &mut scr.at(indx + 6), &mut nb, &mut nw);
        let mut nx = 0i32;
        let iverp;
        if c.n1h != 0 {
            iverp = 4;
            nx = n6;
        } else if c.n2h == 0 {
            iverp = 5;
        } else {
            iverp = 6;
        }
        c.skiprz(nin, -1);
        indx += nw;
        if iverp >= 5 {
            c.contio(nin, 0, 0, &mut scr.at(indx), &mut nb, &mut nw);
        }
        if iverp >= 5 {
            indx += nw;
        }
        if iverp >= 6 {
            c.contio(nin, 0, 0, &mut scr.at(indx), &mut nb, &mut nw);
        }
        if iverp >= 6 {
            indx += nw;
        }
        c.hdatio(nin, 0, 0, &mut scr.at(indx), &mut nb, &mut nw);
        scr[indx] = temp;
        if iverp >= 5 {
            nx = nint(scr[indx + 5]);
        }
        indx += nw;
        while nb != 0 {
            if indx > NWSCR {
                c.error("tpend", "storage exceeded.", " ");
            }
            c.moreio(nin, 0, 0, &mut scr.at(indx), &mut nb, &mut nw);
            indx += nw;
        }
        nw = 6 * nx;
        let mut dico = RV::new(nw.max(0) as usize + 64);
        nw = 6 * nx;
        if self.iinc > 0 {
            nw += 12;
        }
        if icoh > 0 {
            nw += 12;
        }
        if icoh > 20 {
            nw += 12;
        }
        let mut dicn = RV::new(nw.max(0) as usize + 256);
        let mut nc = nx;
        c.dictio(nin, 0, 0, &mut dico.all(), &mut nb, &mut nc);
        for i in 1..=nw {
            dicn[i] = 0.0;
        }

        //--prepare the new directory
        let mut mf2 = 0;
        let mut mfd = 3i32;
        let mut mtd = mtref;
        let lim = 6 * (nx - 1);
        let mut ix = 3i32;
        nc = 3 + (ne + 2) / 3;
        let mut i = 0i32;
        let mut j = 0i32;
        let mut mfi = 0i32;
        let mut mti = 0i32;
        let mut st = 172;
        loop {
            match st {
                172 => {
                    mti = 0;
                    if i > lim {
                        st = 176;
                        continue;
                    }
                    mfi = nint(dico[i + 3]);
                    if mfi == 2 {
                        mf2 = 1;
                    }
                    mti = nint(dico[i + 4]);
                    st = 173;
                }
                173 => {
                    st = 174;
                    if mfi > mfd || (mfi == mfd && mti >= mtd) {
                        st = 176;
                    }
                }
                174 => {
                    dicn[j + 3] = dico[i + 3];
                    dicn[j + 4] = dico[i + 4];
                    dicn[j + 5] = dico[i + 5];
                    i += 6;
                    j += 6;
                    st = 172;
                }
                176 => {
                    if mfd > 6 {
                        break;
                    }
                    dicn[j + 3] = mfd as f64;
                    dicn[j + 4] = mtd as f64;
                    dicn[j + 5] = nc as f64;
                    j += 6;
                    if mtd == mti {
                        i += 6;
                        mfi = nint(dico[i + 3]);
                        mti = nint(dico[i + 4]);
                    }
                    if ix == nex {
                        st = 180;
                        continue;
                    }
                    ix += 1;
                    mtd += 1;
                    if mfd == 6 {
                        nc = self.ncdse;
                    }
                }
                180 => {
                    if mfd == 6 || self.iinc == 0 {
                        st = 184;
                        continue;
                    }
                    nc = self.ncds;
                    ix = 3;
                    mfd = 6;
                    mtd = mtref;
                    st = if i > lim { 176 } else { 173 };
                }
                _ => {
                    // 184
                    mfd = 1000;
                    if mtd == mti {
                        i += 6;
                    }
                    st = if i <= lim { 174 } else { 176 };
                }
            }
        }

        //--write new directory (label 190)
        let nxn = j / 6;
        if iverp < 5 {
            scr[6] = nxn as f64;
        }
        c.contio(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
        indx = 1 + nw;
        if iverp >= 5 {
            c.contio(0, nout, 0, &mut scr.at(indx), &mut nb, &mut nw);
        }
        if iverp >= 5 {
            indx += nw;
        }
        if iverp >= 6 {
            c.contio(0, nout, 0, &mut scr.at(indx), &mut nb, &mut nw);
        }
        if iverp >= 6 {
            indx += nw;
        }
        if iverp >= 5 {
            scr[indx + 5] = nxn as f64;
        }
        c.hdatio(0, nout, 0, &mut scr.at(indx), &mut nb, &mut nw);
        indx += nw;
        while nb != 0 {
            if indx > 1 + NWSCR - 1 {
                c.error("tpend", "storage exceeded.", " ");
            }
            c.moreio(0, nout, 0, &mut scr.at(indx), &mut nb, &mut nw);
            indx += nw;
        }
        dicn[5] = dico[5] + nxn as f64 - nx as f64;
        let mut nxn_m = nxn;
        c.dictio(0, nout, 0, &mut dicn.all(), &mut nb, &mut nxn_m);
        c.tofend(nin, nout, 0, &mut scr.all());
        if mf2 == 1 {
            c.tofend(nin, nout, 0, &mut scr.all());
        }
        let mut mfn = 3i32;
        let mut mtn = mtref;
        let (mut mat, mut mf, mut mt);

        //--copy through to the thermal mf3 sections
        loop {
            c.contio(nin, 0, 0, &mut sav.all(), &mut nb, &mut nw);
            mat = c.math;
            mf = c.mfh;
            mt = c.mth;
            if c.math != matd {
                break;
            }
            if c.mfh == 0 {
                continue;
            }
            if c.mfh > mfn {
                break;
            }
            if c.mfh == mfn && c.mth >= mtn {
                break;
            }
            c.math = mat;
            c.mfh = mf;
            c.mth = mt;
            c.contio(0, nout, 0, &mut sav.all(), &mut nb, &mut nw);
            c.tosend(nin, nout, 0, &mut scr.all());
        }

        //--write the new mf3 sections (label 210)
        ix = 2;
        loop {
            // label 215
            ix += 1;
            scr[1] = self.za;
            scr[2] = self.awr;
            scr[3] = 0.0;
            scr[4] = 0.0;
            scr[5] = 0.0;
            scr[6] = 0.0;
            c.math = matd;
            c.mfh = 3;
            c.mth = mtn;
            c.contio(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
            scr[1] = temp;
            scr[2] = 0.0;
            scr[5] = 1.0;
            scr[6] = (ne + 1) as f64;
            scr[7] = (ne + 1) as f64;
            scr[8] = 2.0;
            let mut k = 8;
            let mut istart = 1;
            loop {
                // label 220
                let mut iend = ne + 1;
                if iend - istart >= c.npage / 2 {
                    iend = istart + c.npage / 2 - 1;
                }
                let mut j = k - 1;
                let mut ib = istart - 1;
                loop {
                    j += 2;
                    ib += 1;
                    if ib >= ne + 1 {
                        ex[1] = ETOP;
                        ex[ix as usize] = 0.0;
                    } else {
                        c.finda(ib, &mut ex[1..], nex, iold, &mut self.bufo.0, NBUF);
                    }
                    scr[j] = sigfig(ex[1], 7, 0);
                    if ex[ix as usize] >= EM9 {
                        scr[j + 1] = sigfig(ex[ix as usize], 7, 0);
                    } else {
                        scr[j + 1] = sigfig(ex[ix as usize], 6, 0);
                    }
                    if ib >= iend {
                        break;
                    }
                }
                nw = j + 1;
                if k != 0 {
                    k = 0;
                    c.tab1io(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                } else {
                    c.moreio(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                }
                if nb == 0 {
                    break;
                }
                istart = iend + 1;
            }
            // label 230
            c.asend(nout, 0);
            if mf == mfn && mt == mtn {
                c.tosend(nin, 0, 0, &mut scr.all());
                c.contio(nin, 0, 0, &mut sav.all(), &mut nb, &mut nw);
                mat = c.math;
                mf = c.mfh;
                mt = c.mth;
            }
            if ix == nex {
                break;
            }
            mtn += 1;
        }
        // label 240
        if mf != 3 {
            c.afend(nout, 0);
        } else {
            c.math = mat;
            c.mfh = mf;
            c.mth = mt;
            c.contio(0, nout, 0, &mut sav.all(), &mut nb, &mut nw);
            c.tofend(nin, nout, 0, &mut scr.all());
            c.contio(nin, 0, 0, &mut sav.all(), &mut nb, &mut nw);
            mat = c.math;
            mf = c.mfh;
            mt = c.mth;
        }

        //--copy through to the thermal mf6 sections (label 250)
        if !(self.iinc == 0 && icoh == 0) {
            let mut ifend = 0;
            mfn = 6;
            mtn = mtref;
            ix = 2;
            loop {
                // label 252
                if mf == 0 || mf > mfn || (mf == mfn && mt >= mtn) {
                    break;
                }
                c.math = mat;
                c.mfh = mf;
                c.mth = mt;
                c.contio(0, nout, 0, &mut sav.all(), &mut nb, &mut nw);
                c.tosend(nin, nout, 0, &mut scr.all());
                c.contio(nin, 0, 0, &mut sav.all(), &mut nb, &mut nw);
                ifend = 0;
                if c.mfh == 0 && c.mth == 0 {
                    ifend = 1;
                }
                mat = c.math;
                mf = c.mfh;
                mt = c.mth;
            }
            loop {
                // label 254
                ix += 1;
                c.contio(nscr, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                let ltt = nint(scr[4]);
                scr[4] = 1.0;
                scr[5] = 1.0;
                c.contio(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                if ltt == 6 {
                    c.tab1io(nscr, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                    c.tab2io(nscr, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                    let nne = c.n2h;
                    for _ie in 1..=nne {
                        c.listio(nscr, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                        scr[8] = 1.0;
                        c.listio(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                    }
                } else if ltt == 5 {
                    c.tab1io(nscr, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                    c.tab2io(nscr, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                    let nne = c.n2h;
                    let xsi = self.xsi.as_ref().unwrap();
                    for ie in 1..=nne {
                        if xsi[ie] == ZERO {
                            c.error("tpend", "cross section=0.", " ");
                        }
                        let rxsec = 1.0 / xsi[ie];
                        let mut indx = 1;
                        c.listio(nscr, 0, 0, &mut scr.at(indx), &mut nb, &mut nw);
                        indx += nw;
                        while nb != 0 {
                            if indx > NWSCR {
                                c.error("tpend", "storage exceeded.", " ");
                            }
                            c.moreio(nscr, 0, 0, &mut scr.at(indx), &mut nb, &mut nw);
                            indx += nw;
                        }
                        let nl1 = nint(scr[6]);
                        let j = nint(scr[5]) / nl1;
                        let indx = 8 - nl1;
                        for i in 1..=j {
                            let yy = scr[nl1 * i + indx] * rxsec;
                            if yy >= EM9 {
                                scr[nl1 * i + indx] = sigfig(yy, 7, 0);
                            } else {
                                scr[nl1 * i + indx] = sigfig(yy, 6, 0);
                            }
                        }
                        let mut indx = 1;
                        c.listio(0, nout, 0, &mut scr.at(indx), &mut nb, &mut nw);
                        while nb != 0 {
                            indx += nw;
                            c.moreio(0, nout, 0, &mut scr.at(indx), &mut nb, &mut nw);
                        }
                    }
                } else if ltt == 4 {
                    c.tab1io(nscr, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                    c.tab2io(nscr, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                    let nne = c.n2h;
                    for _ie in 1..=nne {
                        c.tab2io(nscr, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                        let nmu = c.n2h;
                        for _imu in 1..=nmu {
                            let mut indx = 1;
                            c.tab1io(nscr, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                            while nb != 0 {
                                indx += nw;
                                c.moreio(nscr, nout, 0, &mut scr.at(indx), &mut nb, &mut nw);
                            }
                        }
                    }
                }
                c.tosend(nscr, nout, 0, &mut scr.all());
                if mf != mfn && ix == nex && ifend == 0 {
                    c.afend(nout, 0);
                }
                if mf == mfn && mt == mtn {
                    c.tosend(nin, 0, 0, &mut scr.all());
                    c.contio(nin, 0, 0, &mut sav.all(), &mut nb, &mut nw);
                    mat = c.math;
                    mf = c.mfh;
                    mt = c.mth;
                }
                if ix == nex {
                    break;
                }
                mtn += 1;
            }
        }

        //--copy rest of material to output tape (label 295)
        c.math = mat;
        c.mfh = mf;
        c.mth = mt;
        c.contio(0, nout, 0, &mut sav.all(), &mut nb, &mut nw);
        if mat != 0 {
            c.tomend(nin, nout, 0, &mut scr.all());
        }
        let sec = c.timer();
        let m = format!("\n wrote thermal data for temp ={}{}{}s", ffmt::e1p(temp, 11, 4), " ".repeat(28), ffmt::f(sec, 8, 1));
        c.io.so(&m);
        c.io.se(&m);
    }
}
