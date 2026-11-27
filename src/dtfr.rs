//! DTFR: DTF-format tables (and plots in VIEWR format) from GROUPR or
//! GAMINR output, ported from NJOY2016 `dtfr.f90` (module dtfm).
//!
//! The routines keep the Fortran structure and 1-based indices. The
//! module variables (which keep their values from one DTFR run to the
//! next, as in the Fortran) are the fields of `Dtfr`.

#![allow(clippy::too_many_lines)]
#![allow(clippy::collapsible_if)]
#![allow(clippy::collapsible_else_if)]
#![allow(clippy::assign_op_pattern)]
#![allow(clippy::if_same_then_else)]

use crate::endf::Ctx;
use crate::farr::RV;
use crate::ffmt;
use crate::fio::{Kind, Val};
use crate::fmath::{int, log, log10, nint, pow, powi};
use crate::util::dater;

const NEDMAX: usize = 50;
const NDIM: i32 = 7000;
const NWAMAX: i32 = 40000;
const NWSMAX: i32 = 500000;
/// Extra words past the Fortran array sizes.
const SLACK: usize = 20000;
const ZERO: f64 = 0.0;
const TEN: f64 = 10.0;

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

/// A list-directed character value of length n (blank padded), or `d`.
fn getc(v: &Option<Val>, n: usize, d: &str) -> String {
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

/// `n` values per line in `1p,ew.d`, after `lead` on each line.
fn lines(vals: &[f64], n: usize, w: usize, d: usize, lead: &str) -> String {
    let mut s = String::new();
    for (k, x) in vals.iter().enumerate() {
        if k % n == 0 {
            if k > 0 {
                s.push('\n');
            }
            s.push_str(lead);
        }
        s.push_str(&ffmt::e1p(*x, w, d));
    }
    s
}

/// Module variables of `dtfm`.
pub struct Dtfr {
    nout: i32,
    npend: i32,
    nin: i32,
    nplot: i32,
    egn: RV,
    ngn: i32,
    egg: RV,
    ngp: i32,
    iprint: i32,
    ifilm: i32,
    iedit: i32,
    nlmax: i32,
    ng: i32,
    iptotl: i32,
    ipingp: i32,
    itabl: i32,
    ned: i32,
    nptabl: i32,
    jgthrm: i32,
    hednam: Vec<String>,
    jped: Vec<i32>,
    mted: Vec<i32>,
    multed: Vec<i32>,
    ids: Vec<i32>,
    hisnam: String,
    matd: i32,
    jz: i32,
    dtemp: f64,
    kpos: i32,
    ked: i32,
    lpos: i32,
    led: i32,
    iphph: i32,
    mti: i32,
    mtc: i32,
    nlc: i32,
    bot: f64,
    axl: f64,
    axd: f64,
    top: f64,
    ixo: i32,
    iyo: i32,
    nplt: i32,
    xmin: f64,
    xmax: f64,
    ymin: f64,
    ymax: f64,
    sigz: [f64; 11],
    nsigz: i32,
    x: RV,
    y: RV,
    z: RV,
    a: RV,
    sig: RV,
    /// `nedpos` of `ruin` (a local without save; kept here so that an
    /// unset value is the one from the last run, as near as can be)
    nedpos: i32,
}

impl Default for Dtfr {
    fn default() -> Self {
        let n = 3 + NEDMAX + 1;
        Dtfr {
            nout: 0,
            npend: 0,
            nin: 0,
            nplot: 0,
            egn: RV::new(2501 + SLACK),
            ngn: 0,
            egg: RV::new(400 + SLACK),
            ngp: 0,
            iprint: 0,
            ifilm: 0,
            iedit: 0,
            nlmax: 0,
            ng: 0,
            iptotl: 0,
            ipingp: 0,
            itabl: 0,
            ned: 0,
            nptabl: 0,
            jgthrm: 0,
            hednam: vec!["\0\0\0\0\0\0".to_string(); n],
            jped: vec![0; n],
            mted: vec![0; n],
            multed: vec![0; n],
            ids: vec![0; n],
            hisnam: "\0\0\0\0\0\0".to_string(),
            matd: 0,
            jz: 0,
            dtemp: 0.0,
            kpos: 0,
            ked: 0,
            lpos: 0,
            led: 0,
            iphph: 0,
            mti: 0,
            mtc: 0,
            nlc: 0,
            bot: 0.0,
            axl: 0.0,
            axd: 0.0,
            top: 0.0,
            ixo: 0,
            iyo: 0,
            nplt: 0,
            xmin: 0.0,
            xmax: 0.0,
            ymin: 0.0,
            ymax: 0.0,
            sigz: [0.0; 11],
            nsigz: 0,
            x: RV::new(NDIM as usize + SLACK),
            y: RV::new(NDIM as usize + SLACK),
            z: RV::new(NDIM as usize + SLACK),
            a: RV::new(NWAMAX as usize + SLACK),
            sig: RV::new(NWSMAX as usize + SLACK),
            nedpos: 0,
        }
    }
}

impl Dtfr {
    pub fn new() -> Self {
        Self::default()
    }

    /// Write one line of plot data (unit 0 is standard error in gfortran).
    fn pl(&self, c: &mut Ctx, s: &str) {
        if self.nplot == 0 {
            // unit 0 is nsyso (ERROR_UNIT, connected to 'output')
            c.io.so(s);
        } else {
            c.io.write_line(self.nplot, s);
        }
    }

    /// Convert output of groupr to dtf format.
    pub fn dtfr(&mut self, c: &mut Ctx) {
        let (mut nb, mut nw) = (0i32, 0i32);
        let mut spect = RV::new(2501 + SLACK);
        let mut fcap = RV::new(2501 + SLACK);
        let mut ffis = RV::new(2501 + SLACK);
        let size = 0.30f64;

        //--read  user input and write heading
        let time = c.timer();
        c.io.so(&format!(
            "\n dtfr...produce dtf format from groupr output{}{}s",
            " ".repeat(24),
            ffmt::f(time, 8, 1)
        ));
        c.io.se(&format!("\n dtfr...{}{}s", " ".repeat(61), ffmt::f(time, 8, 1)));
        self.ruin(c);
        let mut nscr = 11;
        if self.nin < 0 {
            nscr = -nscr;
        }
        c.openz(nscr, 1);

        //--initialize plotting
        let mut csz = size;
        if self.ifilm > 1 {
            csz = csz / 2.0;
        }
        self.pl(c, &format!("1 2{} /", ffmt::f(csz, 6, 3)));

        let nin = self.nin;
        //--select next material, dilution, and temperature.
        'mat: loop {
            // label 105
            let v = c.io.input.read(&[Kind::Char(6), Kind::Int, Kind::Int, Kind::Real]);
            self.hisnam = getc(&v[0], 6, "      ");
            self.matd = geti(&v[1], 0);
            let jzd = geti(&v[2], 1);
            self.dtemp = getr(&v[3], 300.0);
            let mut cnorm = 0.0f64;
            let mut cnm = 0.0f64;
            let mut dnorm = 0.0f64;
            let mut ilmax = 0i32;
            let mut ipmax = 0i32;
            self.ixo = 0;
            self.iyo = 0;
            if self.matd == 0 {
                break 'mat;
            }
            c.io.so(&format!(
                "\n\n mat={}  iso={}  sigzero no={}  temp={}",
                ffmt::i(self.matd, 4),
                ffmt::a(&self.hisnam, 6),
                ffmt::i(jzd, 2),
                ffmt::e1p(self.dtemp, 10, 3)
            ));

            //--find desired material on input tape, process header section,
            //--and copy rest of material to scratch tape
            c.repoz(nin);
            c.tpidio(nin, 0, 0, &mut self.a.all(), &mut nb, &mut nw);
            c.repoz(nscr);
            c.nsc = 0;
            let mut loc;
            loop {
                // label 130
                c.contio(nin, 0, 0, &mut self.a.all(), &mut nb, &mut nw);
                if c.math == self.matd {
                    // label 135
                    loc = 7;
                    c.listio(nin, 0, 0, &mut self.a.at(loc), &mut nb, &mut nw);
                    let temp = self.a[7];
                    let mut test = 1.0f64;
                    test = test / 1000.0;
                    if (temp - self.dtemp).abs() <= test {
                        break;
                    }
                    c.tomend(nin, 0, 0, &mut self.a.all());
                    continue;
                }
                if c.math == -1 {
                    c.io.so("\n ***material not found");
                    continue 'mat;
                }
                if c.mfh == 0 {
                    continue;
                }
                c.tomend(nin, 0, 0, &mut self.a.all());
            }
            // label 136
            while nb != 0 {
                loc += nw;
                c.moreio(nin, 0, 0, &mut self.a.at(loc), &mut nb, &mut nw);
            }
            self.nsigz = nint(self.a[4]);
            let ntw = nint(self.a[6]);
            self.ngn = nint(self.a[9]);
            if self.ngn != self.ng {
                c.error("dtfr", "number of neutron groups disagrees with number requested.", " ");
            }
            let ngg = nint(self.a[10]);
            if self.ngp > 0 && ngg != self.ngp {
                c.error("dtfr", "number of gamma groups disagrees with number requested.", " ");
            }
            let mut jbase = 12;
            jbase = jbase + ntw;
            for j in 1..=self.nsigz {
                self.sigz[j as usize] = self.a[j + jbase];
            }
            jbase = jbase + self.nsigz;
            let ngnp1 = self.ngn + 1;
            for j in 1..=ngnp1 {
                self.egn[j] = self.a[j + jbase];
            }
            jbase = jbase + ngnp1;
            let nggp1 = ngg + 1;
            for j in 1..=nggp1 {
                self.egg[j] = self.a[j + jbase];
            }
            c.tofend(nin, 0, 0, &mut self.a.all());
            c.tomend(nin, 0, nscr, &mut self.a.all());

            //--find desired material and temperature on pendf tape
            let npend = self.npend;
            if npend != 0 {
                c.repoz(npend);
                c.tpidio(npend, 0, 0, &mut self.a.all(), &mut nb, &mut nw);
                c.findf(self.matd, 1, 451, npend);
                let test = self.dtemp / 1000.0 + 1.0;
                let mut idone = 0;
                while idone == 0 {
                    c.contio(npend, 0, 0, &mut self.a.all(), &mut nb, &mut nw);
                    if c.math != self.matd {
                        c.error("dtfr", "desired temperature not on pendf", " ");
                    }
                    c.contio(npend, 0, 0, &mut self.a.all(), &mut nb, &mut nw);
                    if c.n1h != 0 {
                        c.iverf = 4;
                    } else if c.n2h == 0 {
                        c.iverf = 5;
                    } else {
                        c.iverf = 6;
                    }
                    c.skiprz(npend, -1);
                    if c.iverf >= 5 {
                        c.contio(npend, 0, 0, &mut self.a.all(), &mut nb, &mut nw);
                    }
                    if c.iverf >= 6 {
                        c.contio(npend, 0, 0, &mut self.a.all(), &mut nb, &mut nw);
                    }
                    c.hdatio(npend, 0, 0, &mut self.a.all(), &mut nb, &mut nw);
                    let t = self.a[1];
                    if (t - self.dtemp).abs() >= test {
                        c.tomend(npend, 0, 0, &mut self.a.all());
                    } else {
                        idone = 1;
                    }
                }
            }

            //--generate requested neutron and photon tables.
            // label 129
            let mut il = 1i32;
            if self.nlmax == 0 {
                il = 0;
            }
            let mut ip = 0i32;
            if self.nlmax == 0 {
                ip = 1;
            }
            let mut i618 = 0i32;
            let ng = self.ng;
            let itabl = self.itabl;
            let iptotl = self.iptotl;
            let ned = self.ned;
            loop {
                // label 110
                c.repoz(nscr);
                let mut ig = self.ngn;
                let mut nws = 0i32;
                if il > 0 {
                    nws = itabl * ng;
                }
                if ip > 0 {
                    nws = self.ngp * ng;
                }
                if nws > NWSMAX {
                    c.error("dtfr", "not enough storage for table.", " ");
                }
                let n3 = ned + 3;
                for i in 1..=NWSMAX {
                    if i <= n3 {
                        self.ids[i as usize] = 0;
                    }
                    self.sig[i] = 0.0;
                }
                let (mut nl, mut nz) = (0i32, 0i32);
                let (mut ng2, mut ig2lo, mut lz);
                let mut jg = 0i32;

                //--read next record.  process according to type.
                loop {
                    // label 150
                    if ig >= self.ngn {
                        c.contio(nscr, 0, 0, &mut self.a.all(), &mut nb, &mut nw);
                        if c.math == 0 {
                            break;
                        }
                        if c.mfh == 0 || c.mth == 0 {
                            continue;
                        }
                        if c.mth == 501 {
                            self.mted[(ned + 3) as usize] = 501;
                        }
                        if c.mfh == 23 {
                            self.iphph = 1;
                        }
                        nl = nint(self.a[3]);
                        if il > nl || ip > nl {
                            // label 155
                            c.tosend(nscr, 0, 0, &mut self.a.all());
                            continue;
                        }
                        nz = nint(self.a[4]);
                        self.jz = jzd;
                        if self.jz > nz {
                            self.jz = 1;
                        }
                    }
                    // label 160
                    c.listio(nscr, 0, 0, &mut self.a.all(), &mut nb, &mut nw);
                    ng2 = nint(self.a[3]);
                    ig2lo = nint(self.a[4]);
                    ig = nint(self.a[6]);
                    lz = 6;
                    let mut nwa = 1 + nw;
                    while nb != 0 {
                        c.moreio(nscr, 0, 0, &mut self.a.at(nwa), &mut nb, &mut nw);
                        nwa = nwa + nw;
                        if nwa > NWAMAX {
                            c.error("dtfr", "not enough storage for record.", " ");
                        }
                    }
                    let mfh = c.mfh;
                    let mth = c.mth;
                    let jz = self.jz;
                    let mut lab: i32 = if il == 1 && mth == 455 {
                        400
                    } else if il > 0 && mfh == 3 {
                        200
                    } else if il > 0 && mfh == 23 {
                        200
                    } else if il > 0 && mfh == 6 {
                        300
                    } else if il > 0 && mfh == 26 {
                        300
                    } else if ip > 0 && (mfh == 16 || mfh == 17) {
                        600
                    } else {
                        150
                    };
                    let a = &self.a;
                    let sig = &mut self.sig;
                    loop {
                        match lab {
                            150 => break,
                            //--process neutron cross sections.
                            200 => {
                                lab = if mth == 1 || mth == 501 {
                                    210
                                } else if mth == 102 || mth == 18 {
                                    220
                                } else if ned > 0 {
                                    240
                                } else {
                                    260
                                };
                            }
                            // save total.
                            210 => {
                                if il > nl {
                                    break;
                                }
                                jg = ng - ig + 1;
                                let loca = lz + il + nl * ((jz - 1) + nz);
                                let mut locs = iptotl + itabl * (jg - 1);
                                self.ids[(ned + 3) as usize] = 1;
                                sig[locs] = a[loca];
                                if il > 1 {
                                    break;
                                }
                                locs = locs - 2;
                                sig[locs] = sig[locs] + a[loca];
                                self.ids[(ned + 1) as usize] = 1;
                                if ned == 0 {
                                    break;
                                }
                                lab = 240;
                            }
                            // save capture and fission self-shielding factors
                            220 => {
                                jg = ng - ig + 1;
                                let loca = lz + il + nl * nz;
                                if mth == 18 {
                                    ffis[jg] = 1.0;
                                }
                                if mth == 18 && a[loca] != ZERO {
                                    ffis[jg] = a[loca + nl * (jz - 1)] / a[loca];
                                }
                                if mth == 102 {
                                    fcap[jg] = 1.0;
                                }
                                if mth == 102 && a[loca] != ZERO {
                                    fcap[jg] = a[loca + nl * (jz - 1)] / a[loca];
                                }
                                lab = if ned > 0 { 240 } else { 260 };
                            }
                            // save edit cross sections
                            240 => {
                                if il > 1 {
                                    break;
                                }
                                for ied in 1..=ned as usize {
                                    if mth == 1 && self.mted[ied] == 300 {
                                        let loca = lz + il + nl * (jz - 1);
                                        let jpos = self.jped[ied];
                                        let locs = jpos + itabl * (jg - 1);
                                        sig[locs] = a[loca] * self.multed[ied] as f64 + sig[locs];
                                        self.ids[ied] = ied as i32;
                                    } else if mth == self.mted[ied] {
                                        jg = ng - ig + 1;
                                        let loca = lz + il + nl * ((jz - 1) + nz);
                                        let locs = self.jped[ied] + itabl * (jg - 1);
                                        sig[locs] = sig[locs] + self.multed[ied] as f64 * a[loca];
                                        self.ids[ied] = ied as i32;
                                    }
                                }
                                lab = 260;
                            }
                            // thermal correction for total and absorption
                            260 => {
                                if il > 1 {
                                    break;
                                }
                                if mth != 2 && mth != self.mti && mth != self.mtc {
                                    break;
                                }
                                jg = ng - ig + 1;
                                if jg < self.jgthrm {
                                    break;
                                }
                                let loca = lz + il + nl * ((jz - 1) + nz);
                                let locs = iptotl + itabl * (jg - 1);
                                if mth == 2 {
                                    sig[locs] = sig[locs] - a[loca];
                                }
                                if mth == self.mti || mth == self.mtc {
                                    sig[locs] = sig[locs] + a[loca];
                                }
                                if mth == 2 {
                                    sig[locs - 2] = sig[locs - 2] - a[loca];
                                }
                                if mth == self.mti || mth == self.mtc {
                                    sig[locs - 2] = sig[locs - 2] + a[loca];
                                }
                                break;
                            }
                            //--accumulate neutron transfer matrices.
                            300 => {
                                if mth == 18 {
                                    lab = 340;
                                    continue;
                                }
                                if mth == 19 || mth == 20 || mth == 21 || mth == 38 {
                                    lab = 340;
                                    continue;
                                }
                                ilmax = il;
                                jg = ng - ig + 1;
                                if mth == 2 && jg >= self.jgthrm {
                                    break;
                                }
                                if mth < 200 || mth > 250 {
                                    lab = 305;
                                    continue;
                                }
                                if jg < self.jgthrm {
                                    break;
                                }
                                if mth == self.mti {
                                    lab = 305;
                                    continue;
                                }
                                if mth == self.mtc && mth - self.mtc + 1 == il {
                                    lab = 305;
                                    continue;
                                }
                                break;
                            }
                            305 => {
                                for k in 2..=ng2 {
                                    let mut jg2 = ng - ig2lo - k + 3;
                                    let mut jpos = self.ipingp + (jg2 - jg);
                                    if jpos > itabl {
                                        jg2 = jg2 - jpos + itabl;
                                    }
                                    if jpos > itabl {
                                        jpos = itabl;
                                    }
                                    if jpos <= iptotl {
                                        jg2 = jg2 + iptotl + 1 - jpos;
                                    }
                                    if jpos <= iptotl {
                                        jpos = iptotl + 1;
                                    }
                                    if jpos <= itabl {
                                        let loca = lz + il + nl * ((jz - 1) + nz * (k - 1));
                                        let locs = jpos + itabl * (jg2 - 1);
                                        sig[locs] = sig[locs] + a[loca];
                                        if il <= 1 {
                                            let locb = iptotl - 2 + itabl * (jg - 1);
                                            sig[locb] = sig[locb] - a[loca];
                                        }
                                    }
                                }
                                break;
                            }
                            //--special branch for fission reactions.
                            //--accumulate nu*sigf and chi.
                            340 => {
                                if il > 1 {
                                    break;
                                }
                                if mth == 18 {
                                    i618 = 1;
                                }
                                if mth == 19 && i618 > 0 {
                                    c.error("dtfr", "mt18 already processed, mt19 not allowed.", " ");
                                }
                                self.ids[(ned + 2) as usize] = 1;
                                if ig == 0 {
                                    lab = 360;
                                    continue;
                                }
                                if ig2lo == 0 {
                                    lab = 380;
                                    continue;
                                }
                                jg = ng - ig + 1;
                                for k in 2..=ng2 {
                                    let loca = lz + il + nl * ((jz - 1) + nz * (k - 1));
                                    let mut locs = iptotl - 1 + itabl * (jg - 1);
                                    let mut sss = a[loca];
                                    if mth == 18 || mth == 19 {
                                        sss = sss * ffis[jg];
                                    }
                                    sig[locs] = sig[locs] + sss;
                                    if self.kpos > 0 {
                                        let locf = lz + il + nl * (jz - 1);
                                        let jg2 = ng - ig2lo - k + 3;
                                        locs = self.kpos + itabl * (jg2 - 1);
                                        sig[locs] = sig[locs] + a[locf] * sss;
                                        cnorm = cnorm + a[locf] * sss;
                                        self.ids[self.ked as usize] = self.ked;
                                    }
                                }
                                lab = 390;
                            }
                            360 => {
                                for k in 1..=ng {
                                    let jgk = ng - k + 1;
                                    jg = jgk;
                                    spect[jgk] = 0.0;
                                    if k >= ig2lo && k <= (ig2lo + ng2 - 1) {
                                        let loca = lz + k - ig2lo + 1;
                                        spect[jgk] = a[loca];
                                    }
                                }
                                break;
                            }
                            380 => {
                                if self.kpos == 0 {
                                    break;
                                }
                                jg = ng - ig + 1;
                                let loca = lz + il + nl * (jz - 1) + jz;
                                let locf = lz + il;
                                let locs = iptotl - 1 + itabl * (jg - 1);
                                let mut sss = a[loca];
                                if mth == 18 || mth == 19 {
                                    sss = sss * ffis[jg];
                                }
                                sig[locs] = sig[locs] + sss;
                                cnm = cnm + sss * a[locf];
                                lab = 390;
                            }
                            390 => {
                                if ig >= ng && cnm != ZERO {
                                    for k in 1..=ng {
                                        let locs = self.kpos + itabl * (k - 1);
                                        sig[locs] = sig[locs] + cnm * spect[k];
                                        cnorm = cnorm + cnm * spect[k];
                                        self.ids[self.ked as usize] = self.ked;
                                    }
                                }
                                break;
                            }
                            //--special branch for delayed neutrons.
                            400 => {
                                if mfh == 5 {
                                    lab = 440;
                                    continue;
                                }
                                jg = ng - ig + 1;
                                let locs = iptotl - 1 + itabl * (jg - 1);
                                let loca = lz + 1;
                                sig[locs] = sig[locs] + a[loca + 1] * a[loca + 2] * ffis[jg];
                                self.ids[(ned + 2) as usize] = 1;
                                dnorm = dnorm + a[loca] * a[loca + 1] * a[loca + 2] * ffis[jg];
                                lab = 200;
                            }
                            440 => {
                                for kg in 2..=ng2 {
                                    jg = ng - ig2lo - kg + 3;
                                    for id in 1..=nl {
                                        let loca = lz + id + nl * (kg - 1);
                                        if self.lpos > 0 {
                                            let locs = self.lpos + itabl * (jg - 1);
                                            sig[locs] = sig[locs] + a[loca];
                                            self.ids[self.led as usize] = self.led;
                                        }
                                        if self.kpos > 0 {
                                            if dnorm == ZERO {
                                                c.error("dtfr", "delayed nubar required", "to add delayed chi to total.");
                                            }
                                            let locs = self.kpos + itabl * (jg - 1);
                                            sig[locs] = sig[locs] + dnorm * a[loca];
                                            cnorm = cnorm + dnorm * a[loca];
                                            self.ids[self.ked as usize] = self.ked;
                                        }
                                    }
                                }
                                break;
                            }
                            //--accumulate photon production matrices.
                            600 => {
                                let ngp = self.ngp;
                                jg = ng - ig + 1;
                                let mut ff = 1.0f64;
                                if mth == 102 {
                                    ff = fcap[jg];
                                }
                                if mth == 18 {
                                    ff = ffis[jg];
                                }
                                if ig == 0 {
                                    for k in 1..=ngp {
                                        let jgp = ngp - k + 1;
                                        spect[jgp] = 0.0;
                                        if k >= ig2lo && k <= (ig2lo + ng2 - 1) {
                                            let loca = lz + k - ig2lo + 1;
                                            spect[jgp] = a[loca];
                                        }
                                    }
                                } else if ig2lo == 0 {
                                    let loca = lz + ip + nl * ((jz - 1) + nz);
                                    for k in 1..=ngp {
                                        let locs = k + ngp * (jg - 1);
                                        sig[locs] = sig[locs] + a[loca] * ff * spect[k];
                                    }
                                } else {
                                    ipmax = ip;
                                    for k in 2..=ng2 {
                                        let jpos = ngp - ig2lo - k + 3;
                                        if jpos >= 1 && jpos <= ngp {
                                            let loca = lz + ip + nl * ((jz - 1) + nz * (k - 1));
                                            let locs = jpos + ngp * (jg - 1);
                                            sig[locs] = sig[locs] + a[loca] * ff;
                                        }
                                    }
                                }
                                break;
                            }
                            _ => unreachable!(),
                        }
                    }
                }

                //--accumulation complete
                //--make transport corrections if desired.
                //--normalize fission chi vector if desired.
                //--print out dtf tables if desired
                // label 800
                if il == 1 && self.kpos != 0 && cnorm != ZERO {
                    for k in 1..=ng {
                        let locs = self.kpos + itabl * (k - 1);
                        self.sig[locs] = self.sig[locs] / cnorm;
                    }
                }
                self.dtfout(c, il, ip, ilmax, ipmax);
                if ip == 0 {
                    il += 1;
                }
                if il > self.nlmax {
                    il = 0;
                }
                if il == 0 {
                    ip += 1;
                }
                if ip > self.nptabl {
                    ip = 0;
                }
                if il > 0 {
                    continue;
                }
                if ip > 0 {
                    continue;
                }
                break;
            }
        }

        //--dtfr is finished.
        // terminate plotting system
        self.pl(c, " 99/");
        c.closz(self.nplot);
        c.closz(nin);
        c.closz(self.npend);
        c.closz(self.nout);
        c.closz(nscr);
        let time = c.timer();
        c.io.so(&format!("{}{}s\n {}", " ".repeat(69), ffmt::f(time, 8, 1), "*".repeat(77)));
    }

    /// Read user input.
    fn ruin(&mut self, c: &mut Ctx) {
        const KMTED: [i32; 50] = [
            2, 4, 16, 6, 7, 8, 9, 17, 102, 107, 103, 19, 20, 21, 18, 104, 105, 106, 111, 112, 113, 114, 108, 109, 32,
            35, 34, 22, 24, 23, 25, 28, 29, 30, 33, 36, 37, 38, 470, 471, 455, 300, 301, 443, 443, 444, 452, 1, 0, 0,
        ];
        const KJPED: [i32; 50] = [
            1, 2, 3, 3, 3, 3, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
            27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 39, 40, 41, 42, 43, 44, 45,
        ];
        const KMULTD: [i32; 50] = [
            1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
            1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1,
        ];
        const HMTID: [&str; 50] = [
            "  els", "  ins", "  n2n", "  n3n", "  ngm", "  nal", "   np", " fdir", "  nnf", " n2nf", " ftot", "   nd",
            "   nt", " nhe3", "  n2p", "  npa", " nt2a", " nd2a", "  n2a", "  n3a", "  nnd", "nnd2a", "nnhe3", "  nna",
            " n2na", " nn3a", " n3na", "  nnp", " nn2a", "n2n2a", "  nnt", "nnt2a", "  n4n", " n3nf", "  chi", " chid",
            "  nud", "  phi", "theat", "kerma", "tdame", " nusf", " totl", "     ", "     ", "     ", "     ", "     ",
            "     ", "     ",
        ];

        //--read initial parameters
        self.mti = 0;
        self.mtc = 0;
        self.nlc = 0;
        self.nout = 0;
        self.npend = 0;
        self.nplot = 0;
        let v = c.io.input.read(&[Kind::Int, Kind::Int, Kind::Int, Kind::Int]);
        self.nin = geti(&v[0], self.nin);
        self.nout = geti(&v[1], self.nout);
        self.npend = geti(&v[2], self.npend);
        self.nplot = geti(&v[3], self.nplot);
        c.openz(self.nin, 0);
        c.openz(self.nout, 1);
        c.openz(self.npend, 0);
        c.openz(self.nplot, 1);
        self.ifilm = 0;
        self.iedit = 0;
        let v = c.io.input.read(&[Kind::Int, Kind::Int, Kind::Int]);
        self.iprint = geti(&v[0], self.iprint);
        self.ifilm = geti(&v[1], self.ifilm);
        self.iedit = geti(&v[2], self.iedit);
        if self.nplot == 0 {
            self.ifilm = 0;
        }
        c.io.so(&format!(
            "\n input gendf unit ..................... {}\n output unit .......................... {}\n input pendf unit ..................... {}\n output plot data ..................... {}",
            ffmt::i(self.nin, 10),
            ffmt::i(self.nout, 10),
            ffmt::i(self.npend, 10),
            ffmt::i(self.nplot, 10)
        ));
        c.io.so(&format!(
            " print option (0 min, 1 max) .......... {}\n film option (none, 1/plot, 4/plot ) .. {}\n edit option (0 table, 1 sep.)......... {}",
            ffmt::i(self.iprint, 10),
            ffmt::i(self.ifilm, 10),
            ffmt::i(self.iedit, 10)
        ));

        //--read parameter values for iedit=0
        if self.iedit == 0 {
            self.ned = 0;
            let mut nthrm = 0;
            let v = c.io.input.read(&[Kind::Int; 7]);
            self.nlmax = geti(&v[0], self.nlmax);
            self.ng = geti(&v[1], self.ng);
            self.iptotl = geti(&v[2], self.iptotl);
            self.ipingp = geti(&v[3], self.ipingp);
            self.itabl = geti(&v[4], self.itabl);
            self.ned = geti(&v[5], self.ned);
            nthrm = geti(&v[6], nthrm);
            if self.ipingp <= self.iptotl {
                c.error("dtfr", "iping.le.iptotl.", " ");
            }
            self.jgthrm = self.ng + 1 - nthrm;
            if nthrm != 0 {
                self.mtc = 0;
                self.nlc = 0;
                let v = c.io.input.read(&[Kind::Int, Kind::Int, Kind::Int]);
                self.mti = geti(&v[0], self.mti);
                self.mtc = geti(&v[1], self.mtc);
                self.nlc = geti(&v[2], self.nlc);
            }

        //--set standard values for iedit=1 (claw format)
        } else {
            self.nlmax = 5;
            self.ng = 30;
            let v = c.io.input.read(&[Kind::Int, Kind::Int]);
            self.nlmax = geti(&v[0], self.nlmax);
            self.ng = geti(&v[1], self.ng);
            self.iptotl = 43 + 3;
            self.ned = 48;
            self.jgthrm = self.ng + 1;
            self.ipingp = self.iptotl + 1;
            self.itabl = self.iptotl + self.ng;
            self.nedpos = self.iptotl - 3;
            self.nptabl = 0;
            for i in 1..=self.nedpos as usize {
                self.hednam[i] = format!("{} ", HMTID[i - 1]);
            }
            for i in 1..=self.ned as usize {
                self.jped[i] = KJPED[i - 1];
                self.mted[i] = KMTED[i - 1];
                self.multed[i] = KMULTD[i - 1];
            }
        }

        //--write out neutron table specifications
        c.io.so(&format!(
            " number of neutron tables ............. {}\n number of neutron groups ............. {}",
            ffmt::i(self.nlmax, 10),
            ffmt::i(self.ng, 10)
        ));
        if self.nlmax > 0 {
            c.io.so(&format!(
                " position of total .................... {}\n position of in-group ................. {}\n table length ......................... {}",
                ffmt::i(self.iptotl, 10),
                ffmt::i(self.ipingp, 10),
                ffmt::i(self.itabl, 10)
            ));
            if self.jgthrm <= self.ng {
                c.io.so(&format!(" thermal break group .................. {}", ffmt::i(self.jgthrm, 10)));
            }
            if self.jgthrm <= self.ng {
                c.io.so(&format!(
                    " thermal incoherent mt ................ {}\n thermal coherent mt .................. {}\n no. coherent legendre orders ......... {}",
                    ffmt::i(self.mti, 10),
                    ffmt::i(self.mtc, 10),
                    ffmt::i(self.nlc, 10)
                ));
            }

            //--read edit descriptions for iedit=0
            if self.iedit == 0 {
                if self.ned + 3 > NEDMAX as i32 {
                    c.error("dtfr", "not enough storage for edits.", " ");
                }
                if self.ned != 0 {
                    self.nedpos = self.iptotl - 3;
                    if self.nedpos != 0 {
                        let k = vec![Kind::Char(6); self.nedpos.max(0) as usize];
                        let v = c.io.input.read(&k);
                        for i in 1..=self.nedpos as usize {
                            if v[i - 1].is_some() {
                                self.hednam[i] = getc(&v[i - 1], 6, "");
                            }
                        }
                    }
                    let k = vec![Kind::Int; 3 * self.ned as usize];
                    let v = c.io.input.read(&k);
                    for j in 1..=self.ned as usize {
                        self.jped[j] = geti(&v[3 * j - 3], self.jped[j]);
                        self.mted[j] = geti(&v[3 * j - 2], self.mted[j]);
                        self.multed[j] = geti(&v[3 * j - 1], self.multed[j]);
                    }
                }
            }

            //--write out edit specifications
            if self.ned > 0 {
                c.io.so(" edit cross sections ..................  name   position reaction multiplicity");
                let mut jlast = 0;
                let mut jedit = 0usize;
                self.kpos = 0;
                self.lpos = 0;
                for i in 1..=self.ned as usize {
                    let mut hword = "      ".to_string();
                    if self.nedpos != 0 && self.jped[i] != jlast {
                        jedit += 1;
                        jlast = self.jped[i];
                        hword = self.hednam[jedit].clone();
                    }
                    if self.mted[i] == 470 {
                        self.kpos = self.jped[i];
                    }
                    if self.mted[i] == 470 {
                        self.ked = i as i32;
                    }
                    if self.mted[i] == 471 {
                        self.lpos = self.jped[i];
                    }
                    if self.mted[i] == 471 {
                        self.led = i as i32;
                    }
                    c.io.so(&format!(
                        "{}{}   {}       {}        {}",
                        " ".repeat(40),
                        ffmt::a(&hword, 6),
                        ffmt::i(self.jped[i], 2),
                        ffmt::i(self.mted[i], 3),
                        ffmt::i(self.multed[i], 2)
                    ));
                }
            }
        }

        //--add three standard edits
        let ned = self.ned as usize;
        let iptotl = self.iptotl as usize;
        self.hednam[iptotl - 1] = "nusigf".to_string();
        self.jped[ned + 2] = self.iptotl - 1;
        self.mted[ned + 2] = 1000;
        self.hednam[iptotl - 2] = "absorp".to_string();
        self.jped[ned + 1] = self.iptotl - 2;
        self.mted[ned + 1] = 1000;
        self.hednam[iptotl] = "total ".to_string();
        self.jped[ned + 3] = self.iptotl;
        self.mted[ned + 3] = 1;

        //--read and display photon table specifications
        self.nptabl = 0;
        self.ngp = 0;
        let v = c.io.input.read(&[Kind::Int, Kind::Int]);
        self.nptabl = geti(&v[0], self.nptabl);
        self.ngp = geti(&v[1], self.ngp);
        c.io.so(&format!(" number of photon tables .............. {}", ffmt::i(self.nptabl, 10)));
        if self.nptabl > 0 {
            c.io.so(&format!(" number of photon groups .............. {}", ffmt::i(self.ngp, 10)));
        }
        c.io.so("\n");
        c.repoz(self.nout);
    }

    /// Write and/or plot group constants in DTF format.
    fn dtfout(&mut self, c: &mut Ctx, il: i32, ip: i32, ilmax: i32, ipmax: i32) {
        let nout = self.nout;
        let ng = self.ng;
        let itabl = self.itabl;
        let iptotl = self.iptotl;
        let ned = self.ned;
        let mut dat = [0.0f64; 7];
        let id = "njoy    ";

        if self.iedit == 0 {
            //--special format 0 -- output dtf tables with internal edits.
            let nw;
            let head;
            if il != 0 {
                nw = itabl * ng;
                head = format!(
                    " il= {} table{} gp{} pos, mat={} iz={} temp={}",
                    ffmt::i(il, 1),
                    ffmt::i(ng, 3),
                    ffmt::i(itabl, 3),
                    ffmt::i(self.matd, 5),
                    ffmt::i(self.jz, 2),
                    ffmt::e1p(self.dtemp, 12, 5)
                );
            } else {
                nw = self.ngp * ng;
                head = format!(
                    " ip= {} table{} gp{} pos, mat={} iz={} temp={}",
                    ffmt::i(ip, 1),
                    ffmt::i(ng, 3),
                    ffmt::i(self.ngp, 3),
                    ffmt::i(self.matd, 5),
                    ffmt::i(self.jz, 2),
                    ffmt::e1p(self.dtemp, 12, 5)
                );
            }
            if nout > 0 {
                c.io.write_line(nout, "");
                c.io.write_line(nout, &head);
            }
            if self.iprint == 1 {
                c.io.so(&format!("\n{}", head));
            }
            let vals: Vec<f64> = (1..=nw).map(|j| self.sig[j]).collect();
            if nout > 0 {
                let s = lines(&vals, 6, 12, 4, "");
                for l in s.split('\n') {
                    c.io.write_line(nout, l);
                }
            }
            if self.iprint == 1 {
                c.io.so(&lines(&vals, 6, 12, 4, " "));
            }
        } else if il != 0 {
            //--special format 1 -- td6 format.  write separate edits
            if il <= 1 {
                let hdat = dater();
                let s = format!(
                    "{}     edit xsec ({}x{}) proc by {} on {}",
                    ffmt::a(&self.hisnam, 6),
                    ffmt::i(ng, 3),
                    ffmt::i(ned, 3),
                    ffmt::a(id, 6),
                    ffmt::a(&hdat, 8)
                );
                if nout > 0 {
                    c.io.write_line(nout, &s);
                }
                if self.iprint == 1 {
                    c.io.so(&format!(" {}", s));
                }
                let ifis = self.ids[(ned + 2) as usize];
                let nj = iptotl - 3;
                for j in 1..=nj {
                    let mut take = false;
                    if j <= 7 {
                        take = true;
                    } else if ifis > 0 && j <= 11 {
                        take = true;
                    } else if j >= nj - 1 {
                        take = true;
                    } else {
                        for i in 1..=ned as usize {
                            if self.jped[i] == j && self.ids[i] > 0 {
                                take = true;
                                break;
                            }
                        }
                    }
                    if !take {
                        continue;
                    }
                    // label 150
                    let mut jpos = j;
                    if j == nj - 1 && ifis == 0 {
                        continue;
                    }
                    if j == nj - 1 {
                        jpos = iptotl - 1;
                    }
                    if j == nj {
                        jpos = iptotl;
                    }
                    let hmti = self.hednam[j as usize].clone();
                    let mut iseq = 0;
                    let mut k = 0usize;
                    for ig in 1..=ng {
                        let locs = jpos + itabl * (ig - 1);
                        k += 1;
                        dat[k] = self.sig[locs];
                        if k >= 6 || ig >= ng {
                            iseq += 1;
                            if k != 6 {
                                for ii in k + 1..=6 {
                                    dat[ii] = 0.0;
                                }
                            }
                            let mut s = String::new();
                            for i in 1..=6 {
                                s.push_str(&ffmt::e1p(dat[i], 12, 5));
                            }
                            s.push_str(&ffmt::a(&self.hisnam, 2));
                            s.push_str(&ffmt::a(&hmti, 5));
                            s.push_str(&ffmt::i(iseq, 1));
                            if nout > 0 {
                                c.io.write_line(nout, &s);
                            }
                            if self.iprint == 1 {
                                c.io.so(&format!(" {}", s));
                            }
                            k = 0;
                        }
                    }
                }
            }

            //--write reduced neutron table with edits removed.
            // label 200
            let l = il - 1;
            if il <= ilmax {
                let med = iptotl - 3;
                let nws = ng * itabl;
                let mut iseq = 0;
                let ltabn = itabl - iptotl + 3;
                let s = format!(
                    "{} l={} n-n table ({}x{})",
                    ffmt::a(&self.hisnam, 6),
                    ffmt::i(l, 1),
                    ffmt::i(ltabn, 3),
                    ffmt::i(ng, 3)
                );
                if nout > 0 {
                    c.io.write_line(nout, &s);
                }
                if self.iprint == 1 {
                    c.io.so(&format!(" {}", s));
                }
                let mut k = 0usize;
                for ig in 1..=ng {
                    for jpos in 1..=itabl {
                        if jpos > med {
                            k += 1;
                            let locs = jpos + itabl * (ig - 1);
                            dat[k] = self.sig[locs];
                            if k >= 6 || locs >= nws {
                                iseq += 1;
                                if k != 6 {
                                    for i in k + 1..=6 {
                                        dat[i] = 0.0;
                                    }
                                }
                                let mut s = String::new();
                                for i in 1..=6 {
                                    s.push_str(&ffmt::e1p(dat[i], 12, 5));
                                }
                                s.push_str(&ffmt::a(&self.hisnam, 2));
                                s.push_str(&ffmt::i(l, 2));
                                s.push_str(&ffmt::i(iseq, 4));
                                if nout > 0 {
                                    c.io.write_line(nout, &s);
                                }
                                if self.iprint == 1 {
                                    c.io.so(&s);
                                }
                                k = 0;
                            }
                        }
                    }
                }
            }
        } else {
            //--write photon table
            // label 250
            let l = ip - 1;
            if ip <= ipmax {
                let ngp = self.ngp;
                let mut iseq = 0;
                if nout > 0 {
                    c.io.write_line(
                        nout,
                        &format!(
                            " {} l={}h n-p table ({}x{})",
                            ffmt::a(&self.hisnam, 6),
                            ffmt::i(l, 1),
                            ffmt::i(ngp, 3),
                            ffmt::i(ng, 3)
                        ),
                    );
                }
                if self.iprint == 1 {
                    c.io.so(&format!(
                        "{} l={} n-p table ({}x{})",
                        ffmt::a(&self.hisnam, 6),
                        ffmt::i(l, 1),
                        ffmt::i(ngp, 3),
                        ffmt::i(ng, 3)
                    ));
                }
                let mut k = 0usize;
                let nws = ng * ngp;
                for ig in 1..=ng {
                    for igp in 1..=ngp {
                        k += 1;
                        let locs = igp + ngp * (ig - 1);
                        dat[k] = self.sig[locs];
                        if k >= 6 || locs >= nws {
                            iseq += 1;
                            if k != 6 {
                                for i in k + 1..=6 {
                                    dat[i] = 0.0;
                                }
                            }
                            let mut s = String::new();
                            for i in 1..=6 {
                                s.push_str(&ffmt::e1p(dat[i], 12, 5));
                            }
                            s.push_str(&ffmt::a(&self.hisnam, 2));
                            s.push_str(&ffmt::i(l, 2));
                            s.push_str(&ffmt::i(iseq, 4));
                            if nout > 0 {
                                c.io.write_line(nout, &s);
                            }
                            if self.iprint == 1 {
                                c.io.so(&s);
                            }
                            k = 0;
                        }
                    }
                }
            }
        }

        //--plot dtf edit quantities and scattering tables
        // label 300
        if self.ifilm > 0 && il == 1 {
            self.ploted(c);
        }
        if self.ifilm > 0 && il == 1 && il <= ilmax {
            self.plotnn(c, il);
        }
        if self.ifilm > 0 && ip == 1 && ip <= ipmax {
            self.plotnp(c, ip);
        }
    }

    /// Generate plots of all edit and cross section positions.
    /// Overlay corresponding pointwise cross sections from pendf.
    fn ploted(&mut self, c: &mut Ctx) {
        let qu = "'";
        let one = 1.0f64;
        let ww = 5.0f64;
        let hh = 3.75f64;
        let ned = self.ned;
        let itabl = self.itabl;

        //--plot cross sections (non-zero reactions only)
        self.nplt = self.iptotl;
        let deltay = 16.0f64;
        for i in 1..=self.nplt {
            let nk = ned + 3;
            let mut ksave = 0i32;
            let mut hit = false;
            for k in 1..=nk {
                ksave = k;
                if self.jped[k as usize] == i && self.ids[k as usize] > 0 {
                    hit = true;
                    break;
                }
            }
            if !hit {
                continue;
            }
            let jpos = i;
            let mut mt = self.mted[ksave as usize];
            if mt == 300 {
                continue;
            }
            if mt == 1 && i != self.iptotl {
                continue;
            }
            let hedn = self.hednam[jpos as usize].clone();

            //--draw histogram.
            let labelz = format!("{}    {}", ffmt::a(&self.hisnam, 6), ffmt::a(&hedn, 6));
            for k in 1..=self.ngn {
                let locs = jpos + itabl * (self.ngn - k);
                self.y[k] = self.sig[locs];
            }
            let nh = self.histod();
            if self.ifilm == 1 {
                self.pl(c, " 1/");
            } else {
                let xll = ww * self.ixo as f64;
                let yll = hh - hh * self.iyo as f64;
                let xur = ww;
                let yur = hh;
                let s = format!(
                    "{}{}{}{}{}/",
                    ffmt::f(xll, 6, 2),
                    ffmt::f(yll, 6, 2),
                    ffmt::f(xur, 6, 2),
                    ffmt::f(yur, 6, 2),
                    ffmt::f(ZERO, 6, 2)
                );
                if self.ixo == 0 && self.iyo == 0 {
                    self.pl(c, &format!(" 1 0 1. 1.{}", s));
                } else {
                    self.pl(c, &format!(" -1 0 1. 1.{}", s));
                }
            }
            let mut xtag = 0.0f64;
            if (500..=522).contains(&mt) {
                if self.ifilm == 1 {
                    xtag = pow(TEN, self.xmax) / 10.0;
                }
                if self.ifilm == 2 {
                    xtag = 3.0 * pow(TEN, self.xmax) / 100.0;
                }
            } else {
                xtag = self.xmin + (self.xmax - self.xmin) / 50.0;
                xtag = pow(TEN, xtag);
            }
            let mut ytag = self.ymax - (self.ymax - self.ymin) / 50.0;
            ytag = pow(TEN, ytag);
            self.pl(c, &format!(" {}<{}{}/", qu, ffmt::a(&labelz, 20), qu));
            self.pl(c, " /");
            self.pl(c, &format!(" 4 0 3 1{}{}/", ffmt::e(xtag, 10, 2), ffmt::e(ytag, 10, 2)));
            self.pl(
                c,
                &format!(
                    "{}{}{}/",
                    ffmt::e1p(pow(TEN, self.xmin), 10, 2),
                    ffmt::e1p(pow(TEN, self.xmax), 10, 2),
                    ffmt::e1p(one, 10, 2)
                ),
            );
            self.pl(c, &format!(" {}<e>nergy (e<v>){}/", qu, qu));
            self.pl(
                c,
                &format!(
                    "{}{}{}/",
                    ffmt::e1p(pow(TEN, self.ymin), 10, 2),
                    ffmt::e1p(pow(TEN, self.ymax), 10, 2),
                    ffmt::e1p(one, 10, 2)
                ),
            );
            self.pl(c, &format!(" {}<c>ross <s>ection (barns){}/", qu, qu));
            self.pl(c, "/");
            self.pl(c, " 0 0 0 0 1/");
            self.pl(c, &format!(" {}<m.g.>{}/", qu, qu));
            self.pl(c, " 0");
            for ih in 1..=nh {
                let s = format!("{}{}/", ffmt::e1p(self.x[ih], 13, 4), ffmt::e1p(self.z[ih], 13, 4));
                self.pl(c, &s);
            }
            self.pl(c, " /");
            if self.ifilm != 1 {
                self.ixo += 1;
                if self.ixo > 1 {
                    self.ixo = 0;
                }
                if self.ixo == 0 {
                    self.iyo += 1;
                }
                if self.iyo > 1 {
                    self.iyo = 0;
                }
            }
            if self.iphph == 1 {
                if c.iverf >= 6 {
                    if mt < 501 || mt > 522 {
                        continue;
                    }
                } else {
                    if mt < 501 || mt > 602 {
                        continue;
                    }
                }
            } else {
                if c.iverf >= 6 {
                    if mt > 150 && mt < 301 {
                        continue;
                    }
                    if mt > 450 && mt < 600 {
                        continue;
                    }
                    if mt > 849 {
                        continue;
                    }
                } else {
                    if mt > 150 && mt < 301 {
                        continue;
                    }
                    if mt > 450 && mt < 700 {
                        continue;
                    }
                    if mt > 799 {
                        continue;
                    }
                }
            }

            //--draw smooth plot(s)
            let mut j = 0i32;
            let mut k = ksave;
            let mut iy1 = 10i32;
            loop {
                // label 140
                j += 1;
                let l1 = format!("mt{}   ", ffmt::i(mt, 3));
                let npts = self.dpend(c, mt);
                let (lin, nsym, ndash);
                if j < 5 {
                    lin = 0;
                    nsym = 0;
                    ndash = j;
                } else if j >= 5 && j < 9 {
                    lin = 1 + npts / 5;
                    nsym = 4;
                    ndash = j - 4;
                } else {
                    lin = 1 + npts / 5;
                    nsym = 3;
                    ndash = j - 8;
                }
                self.pl(c, &format!("{}/ pendf plot for mt={}", ffmt::i(j + 1, 3), ffmt::i(mt, 4)));
                self.pl(c, "/");
                self.pl(c, &format!("{}{}{} 0 1/", ffmt::i(lin, 5), ffmt::i(nsym, 3), ffmt::i(ndash, 3)));
                self.pl(c, &format!(" {}<{}{}/", qu, ffmt::a(&l1, 8), qu));
                self.pl(c, " 0");
                for ipt in 1..=npts {
                    let s = format!("{}{}/", ffmt::e1p(self.x[ipt], 13, 4), ffmt::e1p(self.y[ipt], 13, 4));
                    self.pl(c, &s);
                }
                self.pl(c, " /");
                iy1 = iy1 + nint(deltay);
                // label 150
                let mut next = false;
                loop {
                    k += 1;
                    if k > ned {
                        break;
                    }
                    if self.ids[k as usize] == 0 || self.jped[k as usize] != i {
                        continue;
                    }
                    mt = self.mted[k as usize];
                    next = true;
                    break;
                }
                if !next {
                    break;
                }
            }
            let _ = iy1;
        }
    }

    /// Prepare histogram data for plotting and choose scales for plots.
    /// Uses the module arrays x, y, z with ngn groups and egn; returns n.
    fn histod(&mut self) -> i32 {
        let two = 2.0f64;
        let five = 5.0f64;
        let ptwo = 0.2f64;
        let pfive = 0.5f64;
        let small = 1.0e-10f64;
        let big = 1.0e10f64;
        let fact = 1.7f64;
        let ng = self.ngn;
        let eg = &self.egn;
        let x = &mut self.x;
        let y = &self.y;
        let z = &mut self.z;

        //--construct histogram removing leading and trailing zeros
        let mut yh = small;
        let mut yl = big;
        let mut n = 1i32;
        z[n] = 0.0;
        for i in 1..=ng {
            if y[i] != ZERO || n != 1 {
                n += 2;
                z[n - 1] = y[i];
                z[n] = y[i];
                x[n - 2] = eg[i];
                x[n - 1] = eg[i];
                if y[i] > yh {
                    yh = y[i];
                }
                if y[i] > ZERO && y[i] < yl {
                    yl = y[i];
                }
            }
        }
        x[n] = eg[ng + 1];
        n += 1;
        x[n] = eg[ng + 1];
        z[n] = 0.0;

        //--set scales for plotting
        z[n] = 0.0;
        yh = fact * yh;
        let ayh = log10(yh);
        let mut iph = int(ayh);
        if iph >= 0 {
            iph += 1;
        }
        let ayl = log10(yl);
        let mut ipl = int(ayl - 1.0);
        if (iph - ipl) > 6 {
            ipl = iph - 6;
        }
        if iph == ipl {
            iph += 1;
        }
        self.ymin = ipl as f64;
        self.ymax = iph as f64;
        if iph - ipl <= 3 {
            let mut ynow = self.ymin;
            if yl > two * pow(TEN, ynow) {
                self.ymin = ynow + log10(two);
            }
            if yl > five * pow(TEN, ynow) {
                self.ymin = ynow + log10(five);
            }
            ynow = self.ymax;
            if yh < pfive * pow(TEN, ynow) {
                self.ymax = ynow + log10(pfive);
            }
            if yh < ptwo * pow(TEN, ynow) {
                self.ymax = ynow + log10(ptwo);
            }
        }
        self.bot = pow(TEN, self.ymin);
        self.top = pow(TEN, self.ymax);
        let axh = log10(x[n]);
        let mut iph = int(axh);
        if axh - iph as f64 != ZERO {
            iph += 1;
        }
        self.xmax = iph as f64;
        let ipl = int(log10(x[1]) + small);
        let mut ipx = -4;
        if ipl > ipx {
            ipx = -2;
        }
        if ipl > ipx {
            ipx = 4;
        }
        if ipl > ipx {
            ipx = 5;
        }
        if ipl > ipx {
            ipx = 6;
        }
        self.xmin = ipx as f64;
        if iph - ipx <= 2 {
            let mut xnow = self.xmax;
            if x[n] < pfive * pow(TEN, xnow) {
                self.xmax = xnow + log10(pfive);
            }
            if x[n] < ptwo * pow(TEN, xnow) {
                self.xmax = xnow + log10(ptwo);
            }
            xnow = self.xmin;
            if x[1] > two * pow(TEN, xnow) {
                self.xmin = xnow + log10(two);
            }
            if x[1] > five * pow(TEN, xnow) {
                self.xmin = xnow + log10(five);
            }
        }
        self.axl = self.xmin;
        self.axd = self.xmax - self.xmin;
        for i in 1..=n {
            if z[i] <= self.bot {
                z[i] = self.bot;
            }
        }
        n
    }

    /// Retrieve point cross sections from pendf file and thin for
    /// plotting if necessary (into the module arrays x, y); returns npts.
    fn dpend(&mut self, c: &mut Ctx, mt: i32) -> i32 {
        let ehigh = 5.0e7f64;
        let one = 1.0f64;
        let nstep = 50;
        let (mut nb, mut nw) = (0i32, 0i32);
        let npend = self.npend;
        let mat = self.matd;

        //--find desired reaction.
        let mut npts = 0i32;
        let mut mft = 3;
        if self.iphph == 1 {
            mft = 23;
        }
        c.findf(mat, mft, mt, npend);
        c.contio(npend, 0, 0, &mut self.a.all(), &mut nb, &mut nw);
        let mut e = 0.0f64;
        let mut enext = 0.0f64;
        let mut idis = 0i32;
        let mut s = 0.0f64;
        c.gety1(e, &mut enext, &mut idis, &mut s, npend, &mut self.a.all());
        let fact = 2000.0 / self.axd;
        let step = pow(TEN, one / nstep as f64);
        let elow = pow(TEN, self.axl);

        //--read this reaction, thinning or thickening as necessary.
        let mut ixlast = -100i32;
        let mut elast = enext;
        let mut ns = 0i32;
        loop {
            // label 210
            e = enext;
            if e > step * elast {
                e = step * elast;
            }
            elast = e;
            c.gety1(e, &mut enext, &mut idis, &mut s, npend, &mut self.a.all());
            if enext < elow {
                continue;
            }
            if s < ZERO {
                s = -s;
            }
            if s > self.top {
                s = self.top;
            }
            if s < self.bot {
                s = self.bot;
            }
            if e < elow {
                e = elow;
            }
            let ix = int((log10(e) - self.axl) * fact);
            // 0: label 230, 1: label 260, 2: label 270 (after 250 or direct)
            let lab;
            if ix != ixlast {
                ixlast = ix;
                ns = 0;
                lab = 230;
            } else if s == self.y[npts] {
                lab = 260;
            } else {
                lab = match ns {
                    1 => 230,
                    2 | 3 => 240,
                    4 => 250,
                    5 => 270,
                    _ => 240,
                };
            }
            let mut lab = lab;
            loop {
                match lab {
                    240 => {
                        if (s - self.y[npts]) * (self.y[npts] - self.y[npts - 1]) > ZERO {
                            lab = 260;
                        } else {
                            lab = 230;
                        }
                    }
                    230 => {
                        npts += 1;
                        if npts > NDIM {
                            c.error("dpend", "npts exceeds ndim.", " ");
                        }
                        ns += 1;
                        lab = 260;
                    }
                    250 => {
                        ns = 5;
                        let s2 = self.y[npts - 2];
                        let s1 = self.y[npts - 1];
                        if !(s2 > s1) {
                            self.y[npts - 2] = s1;
                            self.y[npts - 1] = s2;
                        }
                        lab = 270;
                    }
                    270 => {
                        if !(s <= self.y[npts - 2]) {
                            self.y[npts - 2] = s;
                            lab = 260;
                        } else if s >= self.y[npts - 1] {
                            lab = 260;
                        } else {
                            self.y[npts - 1] = s;
                            lab = 260;
                        }
                    }
                    _ => {
                        // label 260
                        self.x[npts] = e;
                        self.y[npts] = s;
                        break;
                    }
                }
            }
            if enext <= ehigh {
                continue;
            }
            return npts;
        }
    }

    /// Generate an isometric plot of a neutron table.
    /// This routine destroys the table.
    fn plotnn(&mut self, c: &mut Ctx, il: i32) {
        let qu = "'";
        let emin = 1.0e3f64;
        let emax = 1.0e8f64;
        let big = 1.0e10f64;

        //--remove edits and cross sections.
        let ngn = self.ng;
        let ltabn = self.itabl - self.iptotl;
        let mut k = 0i32;
        let mut top = -big;
        for ig in 1..=ngn {
            let mut j = self.itabl * (ig - 1) + self.iptotl;
            let ig2 = ngn - ig + 1;
            let dl = log(self.egn[ig2 + 1]) - log(self.egn[ig2]);
            for _ip in 1..=ltabn {
                k += 1;
                j += 1;
                let mut sigj = self.sig[j];
                if self.egn[ig2] < emin {
                    sigj = 0.0;
                }
                if self.egn[ig2 + 1] > emax {
                    sigj = 0.0;
                }
                self.sig[k] = sigj / dl;
                if self.sig[k] > top {
                    top = self.sig[k];
                }
            }
        }
        let mut i = int(log10(top));
        if i >= 0 {
            i += 1;
        }
        let bot = powi(TEN, i - 5);
        let top = powi(TEN, i);
        let lim = ngn * ltabn;
        for k in 1..=lim {
            if self.sig[k] < bot {
                self.sig[k] = bot;
            }
        }

        //--plot the modified table
        let l = il - 1;
        let mut ititle = String::new();
        if self.iphph == 0 {
            ititle = format!("{} l={} neut-neut table", ffmt::a(&self.hisnam, 6), ffmt::i(l, 1));
        }
        if self.iphph == 1 {
            ititle = format!("{} l={} phot-phot table", ffmt::a(&self.hisnam, 6), ffmt::i(l, 1));
        }
        self.xmin = self.egn[1];
        if self.xmin < emin {
            self.xmin = emin;
        }
        let mut i = int(log10(self.xmin));
        if i < 0 {
            i -= 1;
        }
        self.xmin = powi(TEN, i);
        self.xmax = self.egn[ngn + 1];
        if self.xmax > emax {
            self.xmax = emax;
        }
        let mut i = int(log10(self.xmax));
        if i >= 0 {
            i += 1;
        }
        self.xmax = powi(TEN, i);
        self.isohead(c);
        let (xmin, xmax) = (self.xmin, self.xmax);
        self.pl(c, &format!(" {}<{}{}/", qu, ffmt::a(&ititle, 26), qu));
        self.pl(c, " /");
        self.pl(c, " -4 2/");
        self.pl(c, &format!("{}{} 1./", ffmt::e1p(xmin, 12, 4), ffmt::e1p(xmax, 12, 4)));
        self.pl(c, &format!(" {}<s>ec. <e>nergy{}/", qu, qu));
        self.pl(c, &format!("{}{} 1./", ffmt::e1p(xmin, 12, 4), ffmt::e1p(xmax, 12, 4)));
        self.pl(c, &format!(" {}<e>nergy (e<v>){}/", qu, qu));
        self.pl(c, &format!("{}{} 1./", ffmt::e1p(bot, 12, 4), ffmt::e1p(top, 12, 4)));
        self.pl(c, &format!(" {}<x>sec/leth{}/", qu, qu));
        self.isotail(c);
        let mut ig2l = 0i32;
        for ig1 in 1..=ngn {
            if self.egn[ig1] >= xmin && self.egn[ig1 + 1] <= xmax {
                let i = ngn - ig1 + 1;
                self.pl(c, &format!(" {}/", ffmt::e1p(self.egn[ig1], 12, 4)));
                let mut zlast = bot;
                for ig2 in 1..=ngn {
                    if self.egn[ig2] >= xmin && self.egn[ig2 + 1] <= xmax {
                        let j = ngn - ig2 + 1;
                        let k = j - i + self.ipingp - self.iptotl;
                        if k >= 1 {
                            let mut ss = self.sig[k + ltabn * (j - 1)];
                            if k > ltabn {
                                ss = bot;
                            }
                            if ss < bot {
                                ss = bot;
                            }
                            self.pl(c, &format!(" {}{}/", ffmt::e1p(self.egn[ig2], 12, 4), ffmt::e1p(zlast, 12, 4)));
                            self.pl(c, &format!(" {}{}/", ffmt::e1p(self.egn[ig2], 12, 4), ffmt::e1p(ss, 12, 4)));
                            zlast = ss;
                            ig2l = ig2;
                        }
                    }
                }
                self.pl(c, &format!(" {}{}/", ffmt::e1p(self.egn[ig2l + 1], 12, 4), ffmt::e1p(zlast, 12, 4)));
                self.pl(c, &format!(" {}{}/", ffmt::e1p(self.egn[ig2l + 1], 12, 4), ffmt::e1p(bot, 12, 4)));
                self.pl(c, " /");
            }
        }
        self.pl(c, " /");
        if self.ifilm == 1 {
            return;
        }
        self.next_frame();
    }

    /// The frame lines that start an isometric plot, and its title.
    fn isohead(&mut self, c: &mut Ctx) {
        let ww = 5.0f64;
        let hh = 3.75f64;
        if self.ifilm == 1 {
            self.pl(c, " 1/ 3d data");
        } else {
            let xll = ww * self.ixo as f64;
            let yll = hh - hh * self.iyo as f64;
            let xur = ww;
            let yur = hh;
            let s = format!(
                "{}{}{}{}{}/ 3d data",
                ffmt::f(xll, 6, 2),
                ffmt::f(yll, 6, 2),
                ffmt::f(xur, 6, 2),
                ffmt::f(yur, 6, 2),
                ffmt::f(ZERO, 6, 2)
            );
            if self.ixo == 0 && self.iyo == 0 {
                self.pl(c, &format!(" 1 0 1. 1.{}", s));
            } else {
                self.pl(c, &format!(" -1 0 1. 1.{}", s));
            }
        }
    }

    /// The lines after the axis cards of an isometric plot.
    fn isotail(&mut self, c: &mut Ctx) {
        self.pl(c, "/");
        if self.ifilm == 1 {
            self.pl(c, " 15. -15. 15. -2.5 6.5 2.5/");
        } else {
            self.pl(c, " 12. -12. 12. -2. 5.25 2./");
        }
        self.pl(c, " 1/ 3d data");
    }

    fn next_frame(&mut self) {
        self.ixo += 1;
        if self.ixo > 1 {
            self.ixo = 0;
        }
        if self.ixo == 0 {
            self.iyo += 1;
        }
        if self.iyo > 1 {
            self.iyo = 0;
        }
    }

    /// Generate an isometric plot of a photon table.
    /// This routine destroys the table.
    fn plotnp(&mut self, c: &mut Ctx, ip: i32) {
        let qu = "'";
        let emin = 1.0e3f64;
        let emax = 1.0e8f64;
        let ngn = self.ngn;
        let ngp = self.ngp;

        let ititle = format!("{} l={} neut-phot table", ffmt::a(&self.hisnam, 6), ffmt::i(ip - 1, 1));
        self.isohead(c);
        let mut top = 0.0f64;
        for ign in 1..=ngn {
            let i = ngn - ign + 1;
            for igp in 1..=ngp {
                let j = ngp - igp + 1;
                let ss = self.sig[j + ngp * (i - 1)] / (self.egg[igp + 1] - self.egg[igp]);
                if self.egn[ign] >= emin && self.egn[ign + 1] <= emax {
                    if self.egg[igp] >= emin && self.egg[igp + 1] <= emax {
                        if ss > top {
                            top = ss;
                        }
                    }
                }
            }
        }
        let mut i = int(log10(top));
        if i >= 0 {
            i += 1;
        }
        let bot = powi(TEN, i - 5);
        let top = powi(TEN, i);
        let mut xmin = self.egn[1];
        if xmin < emin {
            xmin = emin;
        }
        let mut i = int(log10(xmin));
        if i < 0 {
            i -= 1;
        }
        xmin = powi(TEN, i);
        let mut xmax = self.egn[ngn + 1];
        if xmax > emax {
            xmax = emax;
        }
        let mut i = int(log10(xmax));
        if i >= 0 {
            i += 1;
        }
        xmax = powi(TEN, i);
        self.pl(c, &format!(" {}<{}{}/", qu, ffmt::a(&ititle, 26), qu));
        self.pl(c, " /");
        self.pl(c, " -4 2/");
        self.pl(c, &format!("{}{} 1./", ffmt::e1p(xmin, 12, 4), ffmt::e1p(xmax, 12, 4)));
        self.pl(c, &format!(" {}<g>amma <e>nergy{}/", qu, qu));
        self.pl(c, &format!("{}{} 1./", ffmt::e1p(xmin, 12, 4), ffmt::e1p(xmax, 12, 4)));
        self.pl(c, &format!(" {}<e>nergy (e<v>){}/", qu, qu));
        self.pl(c, &format!("{}{} 1./", ffmt::e1p(bot, 12, 4), ffmt::e1p(top, 12, 4)));
        self.pl(c, &format!(" {}<x>sec/e<v>{}/", qu, qu));
        self.isotail(c);
        let mut igpl = 0i32;
        let mut ss = 0.0f64;
        for ign in 1..=ngn {
            if self.egn[ign] >= xmin && self.egn[ign + 1] <= xmax {
                let i = ngn - ign + 1;
                self.pl(c, &format!(" {}/", ffmt::e1p(self.egn[ign], 12, 4)));
                let mut zlast = bot;
                for igp in 1..=ngp {
                    let j = ngp - igp + 1;
                    ss = self.sig[j + ngp * (i - 1)] / (self.egg[igp + 1] - self.egg[igp]);
                    if ss < bot {
                        ss = bot;
                    }
                    self.pl(c, &format!(" {}{}/", ffmt::e1p(self.egg[igp], 12, 4), ffmt::e1p(zlast, 12, 4)));
                    self.pl(c, &format!(" {}{}/", ffmt::e1p(self.egg[igp], 12, 4), ffmt::e1p(ss, 12, 4)));
                    zlast = ss;
                    igpl = igp;
                }
                self.pl(c, &format!(" {}{}/", ffmt::e1p(self.egg[igpl + 1], 12, 4), ffmt::e1p(zlast, 12, 4)));
                self.pl(c, &format!(" {}{}/", ffmt::e1p(self.egg[igpl + 1], 12, 4), ffmt::e1p(ss, 12, 4)));
                self.pl(c, " /");
            }
        }
        self.pl(c, " /");
        if self.ifilm == 1 {
            return;
        }
        self.next_frame();
    }
}
