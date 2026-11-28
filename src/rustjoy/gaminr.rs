//! GAMINR: multigroup photon interaction cross sections and heating
//! kerma factors, ported from NJOY2016 `gaminr.f90` (module gaminm).
//!
//! The routines keep the Fortran structure and 1-based indices. The
//! module variables and the saved locals of gtflx, gpanel and gtff
//! (which keep their values from one GAMINR run to the next, as in the
//! Fortran) are the fields of `Gaminr`.

#![allow(clippy::too_many_lines)]
#![allow(clippy::too_many_arguments)]
#![allow(clippy::collapsible_if)]
#![allow(clippy::collapsible_else_if)]
#![allow(clippy::assign_op_pattern)]

use crate::endf::Ctx;
use crate::farr::{pack_text, unpack_text, RM, R3, RV};
use crate::ffmt;
use crate::fio::{Kind, Val};
use crate::fmath::{nint, sqrt};
use crate::util::a10;
use crate::util::physics::EPAIR;

const NGMAX: usize = 400;
const IWMAX: usize = 200;
const NWPFF: i32 = 5000;
/// Extra words past the Fortran array sizes.
const SLACK: usize = 4096;
const ZERO: f64 = 0.0;
const EMAX: f64 = 1.0e12;

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

/// `n` blank-padded to 66 characters and cut there.
fn pad(s: &str, n: usize) -> String {
    let mut t: String = s.chars().take(n).collect();
    while t.len() < n {
        t.push(' ');
    }
    t
}

const QP2: [f64; 2] = [-1.0e0, 1.0e0];
const QW2: [f64; 2] = [1.0e0, 1.0e0];
const QP6: [f64; 6] = [-1.0e0, -0.76505532e0, -0.28523152e0, 0.28523152e0, 0.76505532e0, 1.0e0];
const QW6: [f64; 6] = [0.06666667e0, 0.37847496e0, 0.55485838e0, 0.55485838e0, 0.37847496e0, 0.06666667e0];
const QP10: [f64; 10] = [
    -1.0e0,
    -0.9195339082e0,
    -0.7387738651e0,
    -0.4779249498e0,
    -0.1652789577e0,
    0.1652789577e0,
    0.4779249498e0,
    0.7387738651e0,
    0.9195339082e0,
    1.0e0,
];
const QW10: [f64; 10] = [
    0.0222222222e0,
    0.1333059908e0,
    0.2248893420e0,
    0.2920426836e0,
    0.3275397612e0,
    0.3275397612e0,
    0.2920426836e0,
    0.2248893420e0,
    0.1333059908e0,
    0.0222222222e0,
];

/// Module variables of `gaminm` and the saved locals.
pub struct Gaminr {
    igg: i32,
    ngg: i32,
    egg: RV,
    iwt: i32,
    wght: RV,
    matb: i32,
    lord: i32,
    nendf: i32,
    npend: i32,
    ngam1: i32,
    ngam2: i32,
    ntw: i32,
    /// title(i)(1:4); only these characters are ever used
    title: Vec<String>,
    matd: i32,
    mfd: i32,
    mtd: i32,
    ig2pp: i32,
    iprint: i32,
    /// gaminr local `integer::nz=1` (implicit save)
    nz: i32,
    // gtflx saved
    fip: i32,
    fir: i32,
    // gpanel saved
    p_nq: i32,
    p_enext: f64,
    p_elast: f64,
    p_slst: RM,
    p_flst: RM,
    p_idisc: i32,
    p_ng1: i32,
    p_ig1: i32,
    // gtff saved
    zz: f64,
    /// gpanel local `nqp` (not saved; kept so that an unset value is
    /// the one from the last call)
    nqp: i32,
}

impl Default for Gaminr {
    fn default() -> Self {
        Gaminr {
            igg: 0,
            ngg: 0,
            egg: RV::new(NGMAX + 2),
            iwt: 0,
            wght: RV::new(IWMAX + 2),
            matb: 0,
            lord: 0,
            nendf: 0,
            npend: 0,
            ngam1: 0,
            ngam2: 0,
            ntw: 0,
            title: vec!["    ".to_string(); 18],
            matd: 0,
            mfd: 0,
            mtd: 0,
            ig2pp: 0,
            iprint: 0,
            nz: 1,
            fip: 0,
            fir: 0,
            p_nq: 0,
            p_enext: 0.0,
            p_elast: 0.0,
            p_slst: RM::new(10, 10),
            p_flst: RM::new(10, 10),
            p_idisc: 0,
            p_ng1: 0,
            p_ig1: 0,
            zz: 0.0,
            nqp: 0,
        }
    }
}

impl Gaminr {
    pub fn new() -> Self {
        Self::default()
    }

    /// Compute multigroup photon cross sections.
    pub fn gaminr(&mut self, c: &mut Ctx) {
        let (mut nb, mut nwds) = (0i32, 0i32);
        let mut nw: i32;
        let mflst: [i32; 9] = [23, 23, 26, 23, 26, 23, 26, 23, 23];
        let mtlst: [i32; 9] = [501, 502, 502, 504, 504, 516, 516, 602, 621];
        let mtlst6: [i32; 9] = [501, 502, 502, 504, 504, 516, 516, 522, 525];
        let nmlst: [&str; 9] = ["totl", "coht", "coht", "inch", "inch", "pair", "pair", "abst", "heat"];
        let ekn = 12.4e3f64;
        let ncnt = 9;

        //--write heading and read user input
        let mut time = c.timer();
        c.io.so(&format!(
            "\n gaminr...produce photon interaction cross sections{}{}s",
            " ".repeat(18),
            ffmt::f(time, 8, 1)
        ));
        c.io.se(&format!("\n gaminr...{}{}s", " ".repeat(59), ffmt::f(time, 8, 1)));
        self.ruing(c);

        //--allocate storage arrays
        nw = c.npage + 50;
        let mut scr = RV::new(nw as usize + NGMAX + SLACK);
        let mut sdat = RV::new(nw as usize + SLACK);
        let mut pff = RV::new(NWPFF as usize + SLACK);
        let ntoth = self.ngg * 2;
        let mut toth = RV::new(ntoth.max(0) as usize + 2);
        let mut zref = 101.0f64;
        let mut znow = 0.0f64;
        let ngg = self.ngg;
        let mut akn = R3::new((self.lord + 1).max(1) as usize, (ngg + 3).max(1) as usize, ngg.max(1) as usize);
        let mut ng2s = vec![0i32; NGMAX + 2];
        let mut ig2s = vec![0i32; NGMAX + 2];
        let mut z = RM::new(10, 10);
        let mut rt = [0.0f64; 17];
        let mut mtname = " ".repeat(60);
        let mut thresh = 0.0f64;
        let mut en = 0.0f64;
        let mut idis = 0i32;
        let mut ng2 = 0i32;
        let mut iglo = 0i32;
        let mut nq = 0i32;
        let mut igzero = 0i32;
        let mut za: f64;
        let mut awr: f64;
        c.repoz(self.nendf);

        //--determine endf format being used
        c.repoz(self.nendf);
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
        c.repoz(self.nendf);
        c.repoz(self.npend);
        c.tpidio(self.nendf, 0, 0, &mut scr.at(1), &mut nb, &mut nw);
        c.tpidio(self.npend, 0, 0, &mut scr.at(1), &mut nb, &mut nwds);

        //--locate position for new material on old gamout tape.
        c.repoz(self.ngam1);
        c.repoz(self.ngam2);
        c.nsh = 0;
        if self.ngam1 == 0 {
            c.tpidio(0, self.ngam2, 0, &mut scr.at(1), &mut nb, &mut nwds);
        }
        let mut go110 = !(self.ngam1 == 0 || self.ngam2 == 0);
        if go110 {
            c.tpidio(self.ngam1, self.ngam2, 0, &mut scr.at(1), &mut nb, &mut nwds);
        }

        'mat: loop {
            if go110 {
                // 110
                loop {
                    c.contio(self.ngam1, 0, 0, &mut scr.at(1), &mut nb, &mut nw);
                    if c.math == -1 {
                        c.skiprz(self.ngam1, -1);
                        break;
                    }
                    if c.math == 1 {
                        continue;
                    }
                    if c.math == self.matb {
                        c.tomend(self.ngam1, 0, 0, &mut scr.at(1));
                        break;
                    }
                    if c.math > self.matb {
                        c.skiprz(self.ngam1, -1);
                        break;
                    }
                    c.contio(0, self.ngam2, 0, &mut scr.at(1), &mut nb, &mut nw);
                    c.tomend(self.ngam1, self.ngam2, 0, &mut scr.at(1));
                }
            }

            //--search for desired material on pendf tape
            // 160
            c.findf(self.matb, 1, 451, self.npend);
            c.contio(self.npend, 0, 0, &mut scr.at(1), &mut nb, &mut nw);
            za = c.c1h;
            awr = c.c2h;
            if c.iverf >= 5 {
                c.contio(self.npend, 0, 0, &mut scr.at(1), &mut nb, &mut nw);
            }
            if c.iverf >= 6 {
                c.contio(self.npend, 0, 0, &mut scr.at(1), &mut nb, &mut nw);
            }
            c.hdatio(self.npend, 0, 0, &mut scr.at(1), &mut nb, &mut nw);
            for i in 1..=17 {
                rt[i - 1] = scr[6 + i as i32];
            }
            let tt = unpack_text(&rt);
            c.io.so(&format!(
                "\n processing mat {}\n {}\n {}",
                ffmt::i(self.matb, 4),
                "-".repeat(66),
                String::from_utf8_lossy(&tt)
            ));

            //--write head record for this material on gamout tape.
            if self.ngam2 != 0 {
                c.nsh = 1;
                c.math = self.matb;
                c.mfh = 1;
                c.mth = 451;
                scr[1] = za;
                scr[2] = awr;
                scr[3] = 0.0;
                scr[4] = self.nz as f64;
                scr[5] = -1.0;
                scr[6] = self.ntw as f64;
                c.contio(0, self.ngam2, 0, &mut scr.at(1), &mut nb, &mut nw);
                scr[1] = 0.0;
                scr[2] = 0.0;
                scr[3] = self.ngg as f64;
                scr[4] = 0.0;
                scr[5] = 0.0;
                scr[6] = 0.0;
                nw = 6;
                for i in 1..=self.ntw {
                    nw += 1;
                    // tt(i)=title(i); scr(nw)=rt(i)
                    let mut t = unpack_text(&rt);
                    let ti = self.title[i as usize].as_bytes();
                    let k = 4 * (i as usize - 1);
                    for m in 0..4 {
                        t[k + m] = ti[m];
                    }
                    pack_text(&t, &mut rt);
                    scr[nw] = rt[(i - 1) as usize];
                }
                nw += 1;
                scr[nw] = EMAX;
                let nggp1 = self.ngg + 1;
                for i in 1..=nggp1 {
                    nw += 1;
                    if self.egg[i] < EPAIR {
                        self.ig2pp = i;
                    }
                    scr[nw] = self.egg[i];
                }
                nw += 1;
                scr[nw] = 0.0;
                let np = nw - 6;
                scr[5] = np as f64;
                let mut loc = 1;
                nw = np;
                if nw > c.npage {
                    nw = c.npage;
                }
                nw += 6;
                c.listio(0, self.ngam2, 0, &mut scr.at(loc), &mut nb, &mut nw);
                let nleft = np - nw + 6;
                while nb != 0 {
                    loc += nw;
                    nw = nleft;
                    if nw > c.npage {
                        nw = c.npage;
                    }
                    c.moreio(0, self.ngam2, 0, &mut scr.at(loc), &mut nb, &mut nw);
                }
                c.afend(self.ngam2, 0);
            }

            //--main entry point for loop over reactions
            self.matd = self.matb;
            for i in 1..=ntoth {
                toth[i] = 0.0;
            }
            let mut icnt = 0usize;
            let mut mflg = 0i32;
            // 220
            loop {
                let mut auto = mflg == -1;
                if !auto {
                    let v = c.io.input.read(&[Kind::Int, Kind::Int, Kind::Char(66)]);
                    self.mfd = geti(&v[0], self.mfd);
                    self.mtd = geti(&v[1], self.mtd);
                    let text = getc(&v[2], 66, &" ".repeat(66));
                    mflg = self.mfd;
                    if mflg == -1 {
                        auto = true;
                    } else if self.mfd == 0 {
                        break;
                    } else {
                        mtname = pad(&text, 60);
                    }
                }
                if auto {
                    // 250
                    icnt += 1;
                    if icnt > ncnt {
                        break;
                    }
                    self.mfd = mflst[icnt - 1];
                    self.mtd = mtlst[icnt - 1];
                    if c.iverf >= 6 {
                        self.mtd = mtlst6[icnt - 1];
                    }
                    mtname = pad(nmlst[icnt - 1], 60);
                }
                // 270
                time = c.timer();

                // special branch for total heating assigned to mt=621.
                // for version 6 format, heating is mt=525
                if self.mtd == 621 || self.mtd == 525 {
                    // 400
                    let mut ig = 0;
                    let mut ig1 = 0;
                    let nl = 1;
                    while ig < self.ngg {
                        ig += 1;
                        ng2 = 2;
                        let ig2lo = 1;
                        let mut ans = R3::new(1, 1, 2);
                        let l = 1 + 2 * (ig - 1);
                        ans[(1, 1, 1)] = toth[l];
                        ans[(1, 1, 2)] = toth[l + 1];
                        if ig1 == 0 {
                            ig1 = 1;
                            self.group_heading(c, time, &mtname);
                            c.math = self.matb;
                            c.mfh = self.mfd;
                            c.mth = self.mtd;
                            scr[1] = za;
                            scr[2] = awr;
                            scr[3] = nl as f64;
                            scr[4] = self.nz as f64;
                            scr[5] = 0.0;
                            scr[6] = self.ngg as f64;
                            c.contio(0, self.ngam2, 0, &mut scr.at(1), &mut nb, &mut nwds);
                            self.dspla(c, -1, &mut ans, nl, self.nz, ng2, ig2lo, &mut igzero);
                        }
                        self.dspla(c, ig, &mut ans, nl, self.nz, ng2, ig2lo, &mut igzero);
                        self.write_ans(c, &mut scr, &ans, nl, ng2, ig2lo, ig, igzero);
                    }
                    continue;
                }

                // branch for normal reactions
                let mut nl = 1;
                if self.mfd == 26 {
                    nl = self.lord + 1;
                }
                if self.mtd == 516 {
                    nl = 1;
                }
                self.nz = 1;
                let nz = self.nz;
                let mut ans = R3::new(nl as usize, nz as usize, (self.ngg + 3) as usize);
                let mut ff = RM::new(nl as usize, (self.ngg + 3) as usize);
                let e = 0.0f64;
                self.gtsig(c, e, &mut thresh, &mut idis, &mut z, &mut sdat);
                self.gtflx(c, e, &mut en, &mut idis, &mut z, nl, nz);
                self.gtff(c, e, &mut en, &mut idis, &mut ff, nl, &mut ng2, &mut iglo, &mut nq, &mut pff);
                let mut ig1 = 0;
                let mut ig = 0;
                if self.mfd == 26 && self.mtd == 504 {
                    znow = nint(za / 1000.0) as f64;
                    if znow <= zref {
                        zref = znow;
                        znow = -znow;
                    }
                }

                //--loop over initial energy groups
                loop {
                    // 300
                    ig += 1;
                    let mut elo = self.egg[ig];
                    let ehi = self.egg[ig + 1];
                    let mut ig2lo;
                    if znow > ZERO && elo >= znow * ekn {
                        // 325
                        for i in 1..=3 + self.ngg {
                            for j in 1..=nl {
                                ans[(j, 1, i)] = akn[(j, i, ig)];
                                if i > 1 {
                                    ans[(j, 1, i)] = ans[(j, 1, i)] * znow;
                                }
                            }
                        }
                        ng2 = ng2s[ig as usize];
                        ig2lo = ig2s[ig as usize];
                    } else {
                        ig2lo = 0;
                        ng2 = 2;
                        if ehi <= thresh {
                            // 380
                            if ig < self.ngg {
                                continue;
                            }
                            break;
                        }
                        let mut enext = ehi;
                        for it in 1..=self.ngg + 3 {
                            for iz in 1..=nz {
                                for il in 1..=nl {
                                    ans[(il, iz, it)] = 0.0;
                                }
                            }
                        }
                        let mut idone = 0;
                        while idone == 0 {
                            self.gpanel(
                                c, &mut elo, &mut enext, &mut ans, &mut ff, nl, nz, &mut ng2, &mut ig2lo, &mut sdat,
                                &mut pff,
                            );
                            if enext == ehi {
                                idone = 1;
                            } else {
                                elo = enext;
                                enext = ehi;
                            }
                        }
                    }
                    // 330

                    //--print results for this initial energy group
                    if ig1 == 0 {
                        ig1 = 1;
                        self.group_heading(c, time, &mtname);
                        c.math = self.matb;
                        c.mfh = self.mfd;
                        c.mth = self.mtd;
                        scr[1] = za;
                        scr[2] = awr;
                        scr[3] = nl as f64;
                        scr[4] = nz as f64;
                        scr[5] = 0.0;
                        scr[6] = self.ngg as f64;
                        c.contio(0, self.ngam2, 0, &mut scr.at(1), &mut nb, &mut nwds);
                        self.dspla(c, -1, &mut ans, nl, nz, ng2, ig2lo, &mut igzero);
                    }
                    if znow < ZERO {
                        for i in 1..=3 + self.ngg {
                            for j in 1..=nl {
                                akn[(j, i, ig)] = ans[(j, 1, i)];
                                if i > 1 {
                                    akn[(j, i, ig)] = akn[(j, i, ig)] / znow.abs();
                                }
                            }
                        }
                        ng2s[ig as usize] = ng2;
                        ig2s[ig as usize] = ig2lo;
                    }
                    self.dspla(c, ig, &mut ans, nl, nz, ng2, ig2lo, &mut igzero);

                    //--accumulate total heating.
                    if ng2 != 2 && self.mtd != 502 {
                        let l = 1 + 2 * (ig - 1);
                        toth[l] = ans[(1, 1, 1)];
                        toth[l + 1] = toth[l + 1] + ans[(1, 1, ng2)] * ans[(1, 1, 1)];
                        ng2 -= 1;
                        if self.mtd == 504 {
                            ng2 -= 1;
                        }
                    }

                    //--write results on output tape
                    self.write_ans(c, &mut scr, &ans, nl, ng2, ig2lo, ig, igzero);

                    //--continue the loop over groups
                    // 380
                    if ig < self.ngg {
                        continue;
                    }
                    break;
                }
                if ig1 == 0 {
                    c.io.so(&format!(
                        "\n threshold above highest energy bound for mt={}",
                        ffmt::i(self.mtd, 3)
                    ));
                }
                znow = 0.0;
            }

            //--loop over desired materials.
            // 420
            c.amend(self.ngam2, 0);
            self.matb = 0;
            let v = c.io.input.read(&[Kind::Int]);
            self.matb = geti(&v[0], self.matb);
            if self.matb == 0 {
                break 'mat;
            }
            go110 = self.ngam1 != 0;
        }

        //--gaminr is complete.
        // 500
        if self.ngam1 == 0 {
            c.atend(self.ngam2, 0);
        }
        if self.ngam1 != 0 {
            c.totend(self.ngam1, self.ngam2, 0, &mut scr.at(1));
        }
        c.repoz(self.ngam1);
        c.repoz(self.ngam2);
        c.repoz(self.nendf);
        c.repoz(self.npend);
        c.closz(self.ngam1);
        c.closz(self.ngam2);
        c.closz(self.nendf);
        c.closz(self.npend);
        let time = c.timer();
        c.io.so(&format!("{}{}s\n {}", " ".repeat(69), ffmt::f(time, 8, 1), "*".repeat(77)));
    }

    fn group_heading(&self, c: &mut Ctx, time: f64, mtname: &str) {
        c.io.so(&format!("\n group constants{}{}s", " ".repeat(53), ffmt::f(time, 8, 1)));
        c.io.so(&format!(
            " for mf{} and mt{} {}",
            ffmt::i(self.mfd, 2),
            ffmt::i(self.mtd, 3),
            mtname
        ));
    }

    /// Write the results for one group on the output tape.
    fn write_ans(&self, c: &mut Ctx, scr: &mut RV, ans: &R3, nl: i32, ng2: i32, ig2lo: i32, ig: i32, igzero: i32) {
        let (mut nb, mut jw) = (0i32, 0i32);
        let nw = nl * self.nz * ng2;
        let lim = nw;
        if self.ngam2 != 0 && (igzero != 0 || ig == self.ngg) {
            c.math = self.matb;
            c.mfh = self.mfd;
            c.mth = self.mtd;
            scr[1] = 0.0;
            scr[2] = 0.0;
            scr[3] = ng2 as f64;
            scr[4] = ig2lo as f64;
            scr[5] = nw as f64;
            scr[6] = ig as f64;
            let mut j = 6;
            let mut ibase = 6;
            for i in 1..=lim {
                j += 1;
                let mut jj = (i - 1) / nl;
                let ii = i - nl * jj;
                jj += 1;
                scr[j] = ans[(ii, 1, jj)];
                if j >= c.npage + ibase || i == lim {
                    jw = j;
                    if ibase != 0 {
                        c.listio(0, self.ngam2, 0, &mut scr.at(1), &mut nb, &mut jw);
                        ibase = 0;
                        j = 0;
                    } else {
                        c.moreio(0, self.ngam2, 0, &mut scr.at(1), &mut nb, &mut jw);
                        j = 0;
                    }
                }
            }
            let _ = jw;
            if ig == self.ngg {
                c.asend(self.ngam2, 0);
            }
        }
    }

    /// Read user input for group averaging.
    fn ruing(&mut self, c: &mut Ctx) {
        //--read and display user input
        self.ngam1 = 0;
        self.ngam2 = 0;
        let v = c.io.input.read(&[Kind::Int, Kind::Int, Kind::Int, Kind::Int]);
        self.nendf = geti(&v[0], self.nendf);
        self.npend = geti(&v[1], self.npend);
        self.ngam1 = geti(&v[2], self.ngam1);
        self.ngam2 = geti(&v[3], self.ngam2);
        c.openz(self.nendf, 0);
        c.openz(self.npend, 0);
        c.openz(self.ngam1, 0);
        c.openz(self.ngam2, 1);
        self.iprint = 1;
        let v = c.io.input.read(&[Kind::Int, Kind::Int, Kind::Int, Kind::Int, Kind::Int]);
        self.matb = geti(&v[0], self.matb);
        self.igg = geti(&v[1], self.igg);
        self.iwt = geti(&v[2], self.iwt);
        self.lord = geti(&v[3], self.lord);
        self.iprint = geti(&v[4], self.iprint);
        let v = c.io.input.read(&[Kind::Char(66)]);
        let text = getc(&v[0], 66, &" ".repeat(66));
        // read(text,'(16a4,a2)') title
        let tb = text.as_bytes();
        for i in 1..=16 {
            self.title[i] = String::from_utf8_lossy(&tb[4 * (i - 1)..4 * i]).to_string();
        }
        self.title[17] = format!("{}  ", String::from_utf8_lossy(&tb[64..66]));
        c.io.so(&format!(
            "\n unit for endf tape ................... {}\n unit for pendf tape .................. {}\n unit for input gamout tape ........... {}\n unit for output gamout tape .......... {}",
            ffmt::i(self.nendf, 10),
            ffmt::i(self.npend, 10),
            ffmt::i(self.ngam1, 10),
            ffmt::i(self.ngam2, 10)
        ));
        c.io.so(&format!(
            " mat to be processed .................. {}\n gamma group option ................... {}\n weight function option ............... {}\n legendre order ....................... {}\n print option (0 min, 1 max) .......... {}",
            ffmt::i(self.matb, 10),
            ffmt::i(self.igg, 10),
            ffmt::i(self.iwt, 10),
            ffmt::i(self.lord, 10),
            ffmt::i(self.iprint, 10)
        ));
        let mut t = String::new();
        for i in 1..=16 {
            t.push_str(&self.title[i]);
        }
        t.push_str(&self.title[17][..2]);
        c.io.so(&format!("\n run title\n {}\n      {}", "-".repeat(38), t));
        c.io.so(" ");
        if self.ngam2 > 0 {
            self.ntw = 1;
            self.title[1] = "    ".to_string();
        }
        self.genggp(c);
        self.gnwtf(c);
    }

    /// Generate the requested gamma group structure, or read it in.
    fn genggp(&mut self, c: &mut Ctx) {
        const EG2: [f64; 95] = [
            0.005e0, 0.01e0, 0.015e0, 0.02e0, 0.03e0, 0.035e0, 0.04e0, 0.045e0, 0.055e0, 0.06e0, 0.065e0, 0.075e0,
            0.08e0, 0.09e0, 0.1e0, 0.12e0, 0.14e0, 0.15e0, 0.16e0, 0.19e0, 0.22e0, 0.26e0, 0.3e0, 0.325e0, 0.35e0,
            0.375e0, 0.4e0, 0.425e0, 0.45e0, 0.5e0, 0.525e0, 0.55e0, 0.575e0, 0.6e0, 0.625e0, 0.65e0, 0.675e0,
            0.7e0, 0.75e0, 0.8e0, 0.825e0, 0.865e0, 0.9e0, 1.0e0, 1.125e0, 1.2e0, 1.25e0, 1.33e0, 1.42e0, 1.5e0,
            1.6e0, 1.66e0, 1.75e0, 1.875e0, 2.0e0, 2.166e0, 2.333e0, 2.5e0, 2.666e0, 2.833e0, 3.0e0, 3.166e0,
            3.333e0, 3.5e0, 3.65e0, 3.8e0, 3.9e0, 4.0e0, 4.2e0, 4.4e0, 4.5e0, 4.7e0, 5.0e0, 5.2e0, 5.4e0, 5.5e0,
            5.75e0, 6.0e0, 6.25e0, 6.5e0, 6.75e0, 7.0e0, 7.25e0, 7.5e0, 7.75e0, 8.0e0, 8.5e0, 9.0e0, 9.5e0, 10.0e0,
            10.6e0, 11.0e0, 12.0e0, 14.0e0, 20.0e0,
        ];
        const EG3: [f64; 13] =
            [0.01e0, 0.10e0, 0.50e0, 1.0e0, 2.0e0, 3.0e0, 4.0e0, 5.0e0, 6.0e0, 7.0e0, 8.0e0, 9.0e0, 20.0e0];
        const EG4: [f64; 22] = [
            0.01e0, 0.1e0, 0.2e0, 0.4e0, 1.0e0, 1.5e0, 2.0e0, 2.5e0, 3.0e0, 3.5e0, 4.0e0, 4.5e0, 5.0e0, 5.5e0, 6.0e0,
            6.5e0, 7.0e0, 7.5e0, 8.0e0, 10.0e0, 12.0e0, 14.0e0,
        ];
        const EG5: [f64; 23] = [
            0.01e0, 0.03e0, 0.06e0, 0.10e0, 0.15e0, 0.30e0, 0.45e0, 0.60e0, 0.80e0, 1.0e0, 1.33e0, 1.66e0, 2.0e0,
            2.5e0, 3.0e0, 3.5e0, 4.0e0, 5.0e0, 6.0e0, 7.0e0, 8.0e0, 10.0e0, 14.0e0,
        ];
        const EG6: [f64; 49] = [
            0.001e0, 0.01e0, 0.02e0, 0.03e0, 0.045e0, 0.06e0, 0.08e0, 0.1e0, 0.15e0, 0.2e0, 0.3e0, 0.4e0, 0.45e0,
            0.5e0, 0.525e0, 0.6e0, 0.7e0, 0.8e0, 0.9e0, 1.0e0, 1.125e0, 1.2e0, 1.33e0, 1.5e0, 1.66e0, 1.875e0, 2.0e0,
            2.333e0, 2.5e0, 2.666e0, 3.0e0, 3.5e0, 4.0e0, 4.5e0, 5.0e0, 5.5e0, 6.0e0, 6.5e0, 7.0e0, 7.5e0, 8.0e0,
            9.0e0, 10.0e0, 12.0e0, 14.0e0, 17.0e0, 20.0e0, 30.0e0, 50.0e0,
        ];
        const EG7: [f64; 25] = [
            1.0e4, 3.0e4, 6.0e4, 1.0e5, 2.0e5, 3.0e5, 5.0e5, 5.25e5, 7.5e5, 1.0e6, 1.33e6, 1.66e6, 2.0e6, 2.5e6,
            3.0e6, 4.0e6, 5.0e6, 6.0e6, 7.0e6, 8.0e6, 9.0e6, 1.0e7, 1.2e7, 1.7e7, 3.0e7,
        ];
        const EG8: [f64; 39] = [
            0.01e0, 0.02e0, 0.03e0, 0.045e0, 0.06e0, 0.07e0, 0.075e0, 0.10e0, 0.15e0, 0.20e0, 0.30e0, 0.40e0, 0.45e0,
            0.510e0, 0.512e0, 0.60e0, 0.70e0, 0.80e0, 1.0e0, 1.33e0, 1.50e0, 1.66e0, 2.0e0, 2.5e0, 3.0e0, 3.5e0,
            4.0e0, 4.5e0, 5.0e0, 5.5e0, 6.0e0, 6.5e0, 7.0e0, 7.5e0, 8.0e0, 10.0e0, 12.0e0, 14.0e0, 20.0e0,
        ];
        const EG10: [f64; 43] = [
            1.0e3, 1.0e4, 2.0e4, 3.0e4, 4.5e4, 6.0e4, 7.0e4, 7.5e4, 1.0e5, 1.50e5, 2.00e5, 3.00e5, 4.00e5, 4.50e5,
            5.10e5, 5.12e5, 6.00e5, 7.00e5, 8.00e5, 1.00e6, 1.33e6, 1.34e6, 1.50e6, 1.66e6, 2.00e6, 2.50e6, 3.00e6,
            3.50e6, 4.00e6, 4.50e6, 5.00e6, 5.50e6, 6.00e6, 6.50e6, 7.00e6, 7.50e6, 8.00e6, 1.00e7, 1.20e7, 1.40e7,
            2.00e7, 3.00e7, 5.00e7,
        ];
        let mev = 1.0e6f64;

        //--select structure
        self.ngg = 0;
        self.egg[1] = 0.0;
        let igg = self.igg;
        if igg == 0 {
            return;
        }
        let set = |egg: &mut RV, t: &[f64], n: usize, scale: bool| {
            for ig in 1..=n {
                egg[ig as i32] = if scale { t[ig - 1] * mev } else { t[ig - 1] };
            }
        };

        //--group structure is read in.
        if igg == 1 {
            let v = c.io.input.read(&[Kind::Int]);
            self.ngg = geti(&v[0], self.ngg);
            let ngp = self.ngg + 1;
            if ngp > NGMAX as i32 {
                c.error("genggp", "too many groups.", " ");
            }
            let k = vec![Kind::Real; ngp.max(0) as usize];
            let v = c.io.input.read(&k);
            for i in 1..=ngp {
                self.egg[i] = getr(&v[(i - 1) as usize], self.egg[i]);
            }
        } else if igg == 2 {
            self.ngg = 94;
            set(&mut self.egg, &EG2, 95, true);
        } else if igg == 3 {
            self.ngg = 12;
            set(&mut self.egg, &EG3, 13, true);
        } else if igg == 4 {
            self.ngg = 21;
            set(&mut self.egg, &EG4, 22, true);
        } else if igg == 5 {
            self.ngg = 22;
            set(&mut self.egg, &EG5, 23, true);
        } else if igg == 6 {
            self.ngg = 48;
            set(&mut self.egg, &EG6, 49, true);
        } else if igg == 7 {
            self.ngg = 24;
            set(&mut self.egg, &EG7, 25, false);
        } else if igg == 8 || igg == 9 {
            self.ngg = 38;
            if igg == 8 {
                self.ngg = 36;
            }
            let ngm = self.ngg + 1;
            set(&mut self.egg, &EG8, ngm as usize, true);
            if igg != 9 {
                // remove group bounds eg8(7) and eg8(39) if igg=8
                for ig in 7..=ngm {
                    self.egg[ig] = EG8[ig as usize] * mev;
                }
            }
        } else if igg == 10 {
            self.ngg = 42;
            set(&mut self.egg, &EG10, 43, false);
        } else {
            c.error("genggp", "illegal group structure.", " ");
        }

        //--display group structure
        let name = match igg {
            1 => "read in",
            2 => "csewg 94 group",
            3 => "lanl 12 group",
            4 => "steiner 21-group",
            5 => "straker 22 group",
            6 => "lanl 48-group",
            7 => "lanl 24-group",
            8 => "vitamin-c 36-group",
            9 => "vitamin-e 38-group",
            _ => "vitamin-j 42-group",
        };
        c.io.so(&format!("\n gamma group structure......{}", name));
        for ig in 1..=self.ngg {
            c.io.so(&format!(
                " {}  {}  - {}",
                ffmt::i(ig, 5),
                ffmt::e1p(self.egg[ig], 12, 5),
                ffmt::e1p(self.egg[ig + 1], 12, 5)
            ));
        }
    }

    /// Set up the weight function, or read in an arbitrary TAB1 one.
    fn gnwtf(&mut self, c: &mut Ctx) {
        const WT1: [f64; 16] = [
            0.0e0, 0.0e0, 0.0e0, 0.0e0, 1.0e0, 4.0e0, 4.0e0, 5.0e0, 1.0e3, 1.0e-4, 1.0e5, 1.0e0, 1.0e7, 1.0e-2,
            3.0e7, 1.0e-4,
        ];
        if self.iwt == 1 {
            c.io.so("\n weight function......read in");
            let k = vec![Kind::Real; IWMAX];
            let v = c.io.input.read(&k);
            for i in 1..=IWMAX {
                self.wght[i as i32] = getr(&v[i - 1], self.wght[i as i32]);
            }
        } else if self.iwt == 2 {
            c.io.so("\n weight function......constant for all l");
        } else if self.iwt == 3 {
            c.io.so("\n weight function......1/e with rolloffs");
            for i in 1..=16 {
                self.wght[i as i32] = WT1[i - 1];
            }
        } else {
            c.error("gnwtf", "illegal iwt", " ");
        }
    }

    /// Retrieve or compute the weight function.
    fn gtflx(&mut self, c: &mut Ctx, e: f64, enext: &mut f64, idis: &mut i32, flux: &mut RM, nl: i32, _nz: i32) {
        let step = 1.05e0f64;
        *idis = 0;
        if e == ZERO {
            self.fip = 2;
            self.fir = 1;
            *enext = EMAX;
        } else {
            if self.iwt != 2 {
                let mut wtf = 0.0f64;
                c.terpa(&mut wtf, e, enext, idis, &self.wght.at(1), &mut self.fip, &mut self.fir);
                for il in 1..=nl {
                    flux[(1, il)] = wtf;
                }
                let enxt = step * e;
                if *enext > enxt {
                    *idis = 0;
                }
                if *enext > enxt {
                    *enext = enxt;
                }
            } else {
                for il in 1..=nl {
                    flux[(1, il)] = 1.0;
                }
                *enext = EMAX;
                *idis = 0;
            }
        }
    }

    /// Group constant integrals for one panel.
    fn gpanel(
        &mut self,
        c: &mut Ctx,
        elo: &mut f64,
        ehi: &mut f64,
        ans: &mut R3,
        ff: &mut RM,
        nl: i32,
        nz: i32,
        ng: &mut i32,
        iglo: &mut i32,
        sdat: &mut RV,
        pff: &mut RV,
    ) {
        let rndoff = 1.000002e0f64;
        let delta = 0.999995e0f64;
        let mut sig = RM::new(10, 10);
        let mut flux = RM::new(10, 10);
        let mut en = 0.0f64;
        let mut idiscf = 0i32;
        let mut eq = 0.0f64;
        let mut wq = 0.0f64;
        let ehigh;

        //--retrieve factors in integrands at lower boundry.
        if *elo > *ehi {
            c.error("gpanel", "elo gt ehi.", " ");
        }
        if *elo != self.p_elast {
            if *elo * rndoff < *ehi {
                *elo = *elo * rndoff;
            }
            self.p_elast = *elo;
            let mut slst = std::mem::take(&mut self.p_slst);
            let mut flst = std::mem::take(&mut self.p_flst);
            let (mut enx, mut idisc) = (self.p_enext, self.p_idisc);
            self.gtsig(c, *elo, &mut enx, &mut idisc, &mut slst, sdat);
            self.gtflx(c, *elo, &mut en, &mut idiscf, &mut flst, nl, nz);
            if en == enx && idiscf > idisc {
                idisc = idiscf;
            }
            if en < enx {
                idisc = idiscf;
            }
            if en < enx {
                enx = en;
            }
            let (mut ng1, mut ig1, mut nq) = (self.p_ng1, self.p_ig1, self.p_nq);
            self.gtff(c, *elo, &mut en, &mut idiscf, ff, nl, &mut ng1, &mut ig1, &mut nq, pff);
            if en == enx && idiscf > idisc {
                idisc = idiscf;
            }
            if en < enx {
                idisc = idiscf;
            }
            if en < enx {
                enx = en;
            }
            nq += 2;
            if nq > 10 {
                nq = 10;
            }
            self.p_slst = slst;
            self.p_flst = flst;
            self.p_enext = enx;
            self.p_idisc = idisc;
            self.p_ng1 = ng1;
            self.p_ig1 = ig1;
            self.p_nq = nq;
        }

        //--retrieve cross section and flux at upper boundry.
        if self.p_enext < delta * *ehi {
            *ehi = self.p_enext;
            let mut eh = *ehi;
            if self.p_idisc > 0 && *ehi * delta > *elo {
                eh = *ehi * delta;
            }
            ehigh = eh;
        } else {
            ehigh = delta * *ehi;
        }
        let (mut enx, mut idisc) = (self.p_enext, self.p_idisc);
        self.gtsig(c, ehigh, &mut enx, &mut idisc, &mut sig, sdat);
        self.gtflx(c, ehigh, &mut en, &mut idiscf, &mut flux, nl, nz);
        if en == enx && idiscf > idisc {
            idisc = idiscf;
        }
        if en < enx {
            idisc = idiscf;
        }
        if en < enx {
            enx = en;
        }

        //--compute group fluxes assuming flux is linear over panel.
        let aq = (ehigh + *elo) / 2.0;
        let bq = (ehigh - *elo) / 2.0;
        for iz in 1..=nz {
            for il in 1..=nl {
                ans[(il, iz, 1)] = ans[(il, iz, 1)] + (flux[(iz, il)] + self.p_flst[(iz, il)]) * bq;
            }
        }

        //--integrate over panel using lobatto quadrature.
        let nq = self.p_nq;
        let (mut ng1, mut ig1) = (self.p_ng1, self.p_ig1);
        let mut nqp = self.nqp;
        for iq in 1..=nq {
            let q = (iq - 1) as usize;
            if nq == 2 {
                eq = aq + bq * QP2[q];
                wq = bq * QW2[q];
            } else if nq == 6 {
                eq = aq + bq * QP6[q];
                wq = bq * QW6[q];
            } else if nq == 10 {
                eq = aq + bq * QP10[q];
                wq = bq * QW10[q];
            }
            let t1 = (eq - *elo) / (*ehi - *elo);

            //--retrieve the feed function at the quadrature point.
            //--first point was last point of previous panel.
            if iq > 1 {
                if eq > ehigh {
                    eq = ehigh;
                }
                self.gtff(c, eq, &mut en, &mut idiscf, ff, nl, &mut ng1, &mut ig1, &mut nqp, pff);
            }

            //--accumulate the ng*nl*nz integrals simultaneously.
            //--assuming that the reaction rate is linear across the panel.
            if *iglo == 0 {
                *iglo = ig1;
            }
            for iz in 1..=nz {
                for il in 1..=nl {
                    let a = sig[(iz, 1)] * flux[(iz, il)];
                    let b = self.p_slst[(iz, 1)] * self.p_flst[(iz, il)];
                    let mut rr = b + (a - b) * t1;
                    rr = rr * wq;
                    for ig in 1..=ng1 {
                        let igt = ig1 + ig - *iglo + 1;
                        if igt > 1 {
                            ans[(il, iz, igt)] = ans[(il, iz, igt)] + rr * ff[(il, ig)];
                            if igt > *ng {
                                *ng = igt;
                            }
                        }
                    }
                }
            }
        }
        self.p_ng1 = ng1;
        self.p_ig1 = ig1;
        self.nqp = nqp;

        //--save cross section and flux for next panel.
        //--determine next point and quadrature order from last ff.
        self.p_elast = ehigh;
        for iz in 1..=nz {
            self.p_slst[(iz, 1)] = sig[(iz, 1)];
            for il in 1..=nl {
                self.p_flst[(iz, il)] = flux[(iz, il)];
            }
        }
        if en == enx && idiscf > idisc {
            idisc = idiscf;
        }
        if en < enx {
            idisc = idiscf;
        }
        if en < enx {
            enx = en;
        }
        if enx <= *ehi {
            enx = rndoff * *ehi;
        }
        self.p_enext = enx;
        self.p_idisc = idisc;
        let mut nq = nqp + 2;
        if nq > 10 {
            nq = 10;
        }
        self.p_nq = nq;
    }

    /// Display the group constants for one group.
    fn dspla(&self, c: &mut Ctx, ig: i32, ans: &mut R3, nl: i32, _nz: i32, ng2: i32, ig2lo: i32, igzero: &mut i32) {
        let small = 1.0e-9f64;
        let mut result = [0.0f64; 64];
        let mut field: Vec<String> = vec![" ".repeat(10); 64];
        let mfd = self.mfd;
        let mtd = self.mtd;

        //--write appropriate heading.
        *igzero = 0;
        if ig < 0 {
            if self.iprint == 1 {
                c.io.so(" ");
                if mfd != 26 {
                    if mtd == 501 || mtd == 502 || mtd == 504 {
                        c.io.so(" gamma  sigma  \n group  (barns)");
                    }
                    if mtd == 516 {
                        c.io.so(" gamma  sigma  \n group  (barns)");
                    }
                    if mtd == 602 {
                        c.io.so(" gamma  sigma      heating\n group  (barns)    (ev-barns)");
                    }
                    if mtd == 621 {
                        c.io.so(" gamma  heating\n group  (ev-barns)");
                    }
                    if mtd == 522 {
                        c.io.so(" gamma  sigma      heating\n group  (barns)    (ev-barns)");
                    }
                    if mtd == 525 {
                        c.io.so(" gamma  heating\n group  (ev-barns)");
                    }
                } else {
                    if mtd == 502 || mtd == 504 {
                        c.io.so(" initl  final  cross sections vs\n group  group  legendre order   ");
                    }
                    if mtd == 516 {
                        c.io.so(" initl  final  cross sections vs\n group  group  legendre order   ");
                    }
                }
            }
            return;
        }

        //--write out results in appropriate format.
        //--cross sections and heating values.
        if mfd == 23 {
            for i in 2..=ng2 {
                let k = (i - 1) as usize;
                result[k] = 0.0;
                if ans[(1, 1, 1)] != ZERO && ans[(1, 1, i)].abs() >= small {
                    result[k] = ans[(1, 1, i)] / ans[(1, 1, 1)];
                    ans[(1, 1, i)] = result[k];
                    if ans[(1, 1, i)] >= small {
                        *igzero = 1;
                    }
                }
            }
            if result[1] != ZERO && self.iprint == 1 {
                for i in 2..=ng2 {
                    let k = (i - 1) as usize;
                    field[k] = a10(result[k]);
                }
                // (1x,i3,2x,1p,2a11)
                let mut s = format!(" {}  ", ffmt::i(ig, 3));
                let mut n = 0;
                for i in 2..=ng2 {
                    if n == 2 {
                        s.push_str(&format!("\n {}  ", "   "));
                        n = 0;
                    }
                    s.push(' ');
                    s.push_str(&field[(i - 1) as usize]);
                    n += 1;
                }
                c.io.so(&s);
            }

        //--scattering matrices.
        } else if mfd == 26 {
            for ig2 in 2..=ng2 {
                let igt = ig2lo + ig2 - 2;
                if mtd == 516 && ig2 > 2 {
                    result[1] = ans[(1, 1, ig2)] / ans[(1, 1, 1)];
                    ans[(1, 1, ig2)] = result[1];
                    if ans[(1, 1, ig2)] >= small {
                        *igzero = 1;
                    }
                    if self.iprint == 1 {
                        c.io.so(&format!(" {}   heat  {}", ffmt::i(ig, 3), ffmt::e1p(result[1], 11, 3)));
                    }
                } else if igt <= ig {
                    let mut max = nl;
                    if max > 6 {
                        max = 6;
                    }
                    for il in 1..=max {
                        let k = il as usize;
                        result[k] = ans[(il, 1, ig2)] / ans[(il, 1, 1)];
                        ans[(il, 1, ig2)] = result[k];
                        if ans[(il, 1, ig2)] >= small {
                            *igzero = 1;
                        }
                        field[k] = a10(result[k]);
                    }
                    if self.iprint == 1 {
                        let mut s = format!(" {}    {}  ", ffmt::i(ig, 3), ffmt::i(igt, 3));
                        for i in 1..=max {
                            s.push(' ');
                            s.push_str(&field[i as usize]);
                        }
                        c.io.so(&s);
                    }
                    if nl > 6 {
                        for il in 7..=nl {
                            let k = il as usize;
                            result[k] = ans[(il, 1, ig2)] / ans[(il, 1, 1)];
                            ans[(il, 1, ig2)] = result[k];
                            if ans[(il, 1, ig2)] >= small {
                                *igzero = 1;
                            }
                            field[k] = a10(result[k]);
                        }
                        if self.iprint == 1 {
                            // (13x,1p,6a11)
                            let mut s = " ".repeat(13);
                            let mut n = 0;
                            for i in 7..=nl {
                                if n == 6 {
                                    s.push('\n');
                                    s.push_str(&" ".repeat(13));
                                    n = 0;
                                }
                                s.push(' ');
                                s.push_str(&field[i as usize]);
                                n += 1;
                            }
                            c.io.so(&s);
                        }
                    }
                } else {
                    result[1] = ans[(1, 1, ig2)] / ans[(1, 1, 1)];
                    ans[(1, 1, ig2)] = result[1];
                    if ans[(1, 1, ig2)] >= small {
                        *igzero = 1;
                    }
                    field[1] = a10(result[1]);
                    if self.iprint == 1 {
                        if igt == ig + 1 {
                            c.io.so(&format!(" {}   xsec   {}", ffmt::i(ig, 3), field[1]));
                        }
                        if igt == ig + 2 {
                            c.io.so(&format!(" {}   heat   {}", ffmt::i(ig, 3), field[1]));
                        }
                    }
                }
            }
        }
    }

    /// Retrieve the reaction cross section defined by mfd and mtd.
    fn gtsig(&mut self, c: &mut Ctx, e: f64, enext: &mut f64, idis: &mut i32, sig: &mut RM, sdat: &mut RV) {
        let (mut nb, mut nw) = (0i32, 0i32);
        let mut s = 0.0f64;
        if e == ZERO {
            let mf = 23;
            let mt = self.mtd;
            c.findf(self.matd, mf, mt, self.npend);
            c.contio(self.npend, 0, 0, &mut sdat.at(1), &mut nb, &mut nw);
            c.gety1(e, enext, idis, &mut s, self.npend, &mut sdat.at(1));
        } else {
            c.gety1(e, enext, idis, &mut s, self.npend, &mut sdat.at(1));
            sig[(1, 1)] = s;
        }
    }

    /// Compute the feed function or yield for the desired reaction type.
    fn gtff(
        &mut self,
        c: &mut Ctx,
        e: f64,
        enext: &mut f64,
        idisc: &mut i32,
        ff: &mut RM,
        nl: i32,
        ng: &mut i32,
        iglo: &mut i32,
        nq: &mut i32,
        pff: &mut RV,
    ) {
        let c1 = 57.03156e-6f64;
        let c2 = 0.249467e0f64;
        let c3 = 1.95693e-6f64;
        let c4 = 0.0485262e0f64;
        let c5 = 20.60744e0f64;
        let toler = 1.0e-6f64;
        let rndoff = 1.0000001e0f64;
        let two = 2.0e0f64;
        let half = 0.5e0f64;
        let close = 0.99999e0f64;
        let mfd = self.mfd;
        let mtd = self.mtd;
        let (mut nb, mut nw) = (0i32, 0i32);
        let mut arg = [0.0f64; 64];
        let mut pl = [0.0f64; 64];
        let mut idis = 0i32;
        let mut ip = 0i32;
        let mut ir = 0i32;
        let mut snow = 0.0f64;
        let mut xnext = 0.0f64;
        let mut test: f64;

        //--photon interaction cross sections and heating.
        if mfd == 23 {
            *nq = 0;
            *iglo = 1;
            *ng = 0;
            if mtd == 501 || mtd == 502 || mtd == 504 {
                *ng = 1;
            }
            if mtd == 516 {
                *ng = 1;
            }
            if mtd == 602 {
                *ng = 2;
            }
            if mtd == 522 {
                *ng = 2;
            }
            if *ng == 0 {
                c.error("gtff", "illegal reaction for cross section.", " ");
            }
            *enext = EMAX;
            *idisc = 0;
            if e == ZERO {
                return;
            }
            ff[(1, 1)] = 1.0;
            if mtd == 602 {
                ff[(1, 2)] = e;
            }
            if mtd == 522 {
                ff[(1, 2)] = e;
            }

        //--photon scattering cross sections.
        //--initialize.
        } else if mfd == 26 {
            if e == ZERO || mtd != 502 {
                if e == ZERO || mtd != 504 {
                    if mtd != 516 {
                        c.findf(self.matd, 27, mtd, self.nendf);
                        let mut l = 1i32;
                        c.contio(self.nendf, 0, 0, &mut pff.at(l), &mut nb, &mut nw);
                        c.tab1io(self.nendf, 0, 0, &mut pff.at(l), &mut nb, &mut nw);
                        self.zz = pff[l + 1];
                        l += nw;
                        while nb != 0 {
                            if l > NWPFF {
                                c.error("gtff", "insufficient storage for form factor.", " ");
                            }
                            c.moreio(self.nendf, 0, 0, &mut pff.at(l), &mut nb, &mut nw);
                            l += nw;
                        }
                        *ng = 1;
                        if mtd == 504 {
                            *ng = self.ngg + 2;
                        }
                        *enext = EMAX;
                        *idisc = 0;
                        return;
                    }
                }
            }

            //--photon coherent scattering.
            //--compute all legendre components at this energy.
            if mtd == 502 {
                let mut igp = 1;
                while e >= self.egg[igp + 1] {
                    igp += 1;
                }
                for il in 1..=nl {
                    ff[(il, 1)] = 0.0;
                }
                *ng = 1;
                *iglo = igp;
                let mut xnow = 0.0f64;
                let mut unow = 1.0f64;
                ip = 2;
                ir = 1;
                c.terpa(&mut snow, xnow, &mut xnext, &mut idis, &pff.at(1), &mut ip, &mut ir);
                let fact = 2.0 * c2 / (c1 * c1 * e * e);
                let xlim = sqrt(two) * c1 * e;
                let c1e = 1.0 / ((c1 * e) * (c1 * e));
                for il in 1..=nl {
                    arg[il as usize] = 2.0 * (snow * snow);
                }
                let stest = toler * snow;

                //--loop over panels of coherent form factor.
                let mut idone = 0;
                while idone == 0 {
                    let aq = (xnext + xnow) / 2.0;
                    let bq = (xnext - xnow) / 2.0;
                    *nq = 6;
                    if nl > 4 {
                        *nq = 10;
                    }
                    for iq in 1..=*nq {
                        let q = (iq - 1) as usize;
                        let (xq, wq) = if *nq == 6 {
                            (aq + bq * QP6[q], bq * QW6[q])
                        } else {
                            (aq + bq * QP10[q], bq * QW10[q])
                        };
                        xnow = xq * rndoff;
                        if iq > 1 {
                            unow = 1.0 - c1e * xnow * xnow;
                            test = -1.0;
                            if unow < test {
                                unow = -1.0;
                            }
                            c.terpa(&mut snow, xnow, &mut xnext, &mut idis, &pff.at(1), &mut ip, &mut ir);
                            crate::mathm::legndr(unow, &mut pl[1..], nl);
                            for il in 1..=nl {
                                let k = il as usize;
                                arg[k] = (1.0 + unow * unow) * snow * snow * pl[k];
                            }
                        }
                        for il in 1..=nl {
                            ff[(il, 1)] = ff[(il, 1)] + wq * fact * xnow * arg[il as usize];
                        }
                    }
                    test = -1.0;
                    if unow <= test {
                        idone = 1;
                    } else {
                        if xnow >= xlim {
                            idone = 1;
                        } else {
                            if snow < stest {
                                idone = 1;
                            } else {
                                if xnext > xlim {
                                    xnext = xlim;
                                }
                            }
                        }
                    }
                }
                let sigcoh = ff[(1, 1)];
                for il in 1..=nl {
                    ff[(il, 1)] = ff[(il, 1)] / sigcoh;
                }
                *enext = self.egg[igp + 1];
                *idisc = 1;
                *nq = 0;

            //--photon incoherent scattering and heating.
            //--compute for all legendre orders and all sink groups.
            } else if mtd == 504 {
                let mut igp = 1;
                while e >= self.egg[igp + 1] {
                    igp += 1;
                }
                let enow = c3 * e;
                let enowi = 1.0 / enow;
                let enow2 = enow * enow;
                let mut pnow = enow;
                let mut xnow = 0.0f64;
                let xzz = c5 * sqrt(enow / 500.0);
                let qq = 2.0 * enow * (1.0 + enow) / (1.0 + 2.0 * enow);
                let q2m = qq * qq;
                let zz = self.zz;
                if xzz <= zz {
                    ip = 2;
                    ir = 1;
                    c.terpa(&mut snow, xnow, &mut xnext, &mut idis, &pff.at(1), &mut ip, &mut ir);
                } else {
                    snow = zz;
                    xnext = EMAX;
                }
                for il in 1..=nl {
                    arg[il as usize] = 0.0;
                    let lim = igp + 2;
                    for ig in 1..=lim {
                        ff[(il, ig)] = 0.0;
                    }
                }
                let ig = igp;
                let mut siginc = 0.0f64;
                let mut ebar = 0.0f64;
                if e == self.egg[igp] {
                    igp -= 1;
                }
                let mut unow = 0.0f64;
                let mut unext;
                let mut pnext;
                test = 0.0;

                //--loop over panels defined by group and form factor breaks.
                let mut ifini = 0;
                while ifini == 0 {
                    let mut idone = 0;
                    unext = 0.0;
                    while idone == 0 {
                        let t = c4 * xnext;
                        let q2 = t * t;
                        unext = -1.0;
                        if q2 > q2m {
                            idone = 1;
                        } else {
                            unext = 1.0 - ((1.0 - q2 * enowi) - sqrt(1.0 + q2)) / (q2 - enow2 - 2.0 * enow);
                            test = -1.0;
                            if unext < test {
                                unext = -1.0;
                            }
                            if unext < close {
                                idone = 1;
                            } else {
                                xnow = xnext * rndoff;
                                c.terpa(&mut snow, xnow, &mut xnext, &mut idis, &pff.at(1), &mut ip, &mut ir);
                            }
                        }
                    }
                    pnext = enow / (1.0 + 2.0 * enow);
                    if igp > 0 {
                        pnext = c3 * self.egg[igp];
                    }
                    let px = enow / (1.0 + enow * (1.0 - unext));
                    if px > pnext {
                        pnext = px;
                    }
                    if pnext > pnow / rndoff {
                        pnext = pnow / rndoff;
                    }
                    let aq = (pnext + pnow) / 2.0;
                    let bq = (pnext - pnow) / 2.0;
                    *nq = 6;
                    if nl > 6 {
                        *nq = 10;
                    }
                    for iq in 1..=*nq {
                        let q = (iq - 1) as usize;
                        let (uq, wq) = if *nq == 6 {
                            (aq + bq * QP6[q], -(c2 * bq * QW6[q]))
                        } else {
                            (aq + bq * QP10[q], -(c2 * bq * QW10[q]))
                        };
                        pnow = uq;
                        if iq != 1 {
                            let pnowi = 1.0 / pnow;
                            unow = 1.0 + enowi - pnowi;
                            test = 1.0;
                            if unow > test {
                                unow = 1.0;
                            }
                            if xzz <= zz {
                                let rm2 = (1.0 - unow) / 2.0;
                                let rm = sqrt(rm2);
                                let rt = 1.0 + 2.0 * enow * rm2;
                                xnow = c5 * 2.0 * enow * rm * sqrt(rt + enow2 * rm2) / rt;
                                xnow = xnow * rndoff;
                                c.terpa(&mut snow, xnow, &mut xnext, &mut idis, &pff.at(1), &mut ip, &mut ir);
                            }
                            crate::mathm::legndr(unow, &mut pl[1..], nl);
                            let dk = unow - 1.0;
                            let fact = snow * (enow * pnowi + pnow * enowi + dk * (2.0 + dk)) / enow2;
                            for il in 1..=nl {
                                let k = il as usize;
                                arg[k] = fact * pl[k];
                            }
                        }
                        if igp != 0 {
                            for il in 1..=nl {
                                ff[(il, igp)] = ff[(il, igp)] + wq * arg[il as usize];
                            }
                        }
                        siginc = siginc + wq * arg[1];
                        ebar = ebar + wq * arg[1] * pnow / c3;
                    }
                    if unow < -close {
                        ifini = 1;
                    } else {
                        if igp > 0 {
                            if pnext <= c3 * self.egg[igp] {
                                igp -= 1;
                            }
                        }
                    }
                }

                //--strip sink groups with zero cross section.
                if igp <= 0 {
                    igp = 1;
                }
                *iglo = igp;
                *ng = ig - igp + 1;
                for i in 1..=*ng {
                    let j = *iglo + i - 1;
                    for il in 1..=nl {
                        ff[(il, i)] = ff[(il, j)];
                    }
                }
                ebar = ebar / siginc;
                ff[(1, *ng + 1)] = siginc;
                ff[(1, *ng + 2)] = (e - ebar) * siginc;
                *ng += 2;
                for il in 1..=nl {
                    for i in 1..=*ng {
                        ff[(il, i)] = ff[(il, i)] / siginc;
                    }
                }

                //--determine next critical point for incoherent scattering.
                if c3 * e >= half {
                    *enext = EMAX;
                } else {
                    let mut idone = 0;
                    let mut i = 0;
                    while i < self.ngg && idone == 0 {
                        i += 1;
                        if c3 * self.egg[i] >= test {
                            idone = 1;
                        } else {
                            *enext = self.egg[i] / (1.0 - 2.0 * c3 * self.egg[i]);
                            if *enext > e {
                                idone = 2;
                            }
                        }
                    }
                    if idone < 2 {
                        *enext = EMAX;
                    }
                }
                if self.egg[ig + 1] < *enext {
                    *enext = self.egg[ig + 1];
                }
                *idisc = 1;
                *nq += 2;

            //--pair production matrix
            } else if mtd == 516 {
                *enext = EMAX;
                *idisc = 0;
                if e > ZERO {
                    let mut igp = 1;
                    while e >= self.egg[igp + 1] {
                        igp += 1;
                    }
                    let _ = igp;
                    for il in 1..=nl {
                        ff[(il, 1)] = 0.0;
                    }
                    *ng = 2;
                    *iglo = self.ig2pp;
                    *nq = 0;
                    ff[(1, 1)] = 2.0;
                    ff[(1, 2)] = e - 2.0 * EPAIR;
                }
            }

        //--bad file type
        } else {
            c.error("gtff", "illegal file type.", " ");
        }
        let _ = (ip, ir, idis);
    }
}
