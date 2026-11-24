//! Port of NJOY2016 `broadr.f90` (module `broadm`): Doppler broadening and
//! thinning of pointwise cross sections.
//!
//! The structure follows the Fortran routine by routine. Labels of the
//! Fortran `go to` targets are kept as comments or as states of small
//! `match` loops, so the two sources can be read side by side.

use crate::endf::Ctx;
use crate::farr::{IV, RM, RV};
use crate::ffmt;
use crate::fio::{Kind, Val};
use crate::fmath::{exp, nint, powi, sqrt};
use crate::mathm::erfc;
use crate::util::physics::{BK, PI};
use crate::util::sigfig;

const NTT: usize = 160;
const NBUF: i32 = 1000;
const NSTACK: usize = 12;

/// Module variables of `broadm`, plus the `SAVE`d locals of `broadn`.
/// They keep their values between calls and between `broadr` runs in one
/// deck, as in the Fortran.
#[derive(Clone)]
pub struct Broadr {
    nunr: i32,
    intunr: i32,
    nsig0: i32,
    ncyc: i32,
    e1: f64,
    eumax: f64,
    alpha: f64,
    tempk: f64,
    tempef: f64,
    errthn: f64,
    errmax: f64,
    errint: f64,
    thnmax: f64,
    klow: i32,
    khigh: i32,
    mlow: i32,
    mhigh: i32,
    nlow: i32,
    nhigh: i32,
    n2in: i32,
    n2out: i32,
    n2left: i32,
    mpage: i32,
    nreac: i32,
    ntx: i32,
    // funky/hunky
    f: [f64; 6],
    aa: f64,
    x: f64,
    y: f64,
    oy: f64,
    xx: f64,
    yy: f64,
    h: [f64; 6],
    alast: f64,
    s1: f64,
    s2: f64,
    qmtr: RV,
    emtr: RV,
    // broadn: save ks,es,js,ss,tt,sn,dl ; nmax=10 (implicit save)
    bn_ks: [i32; NSTACK + 1],
    bn_js: [i32; NSTACK + 1],
    bn_es: [f64; NSTACK + 1],
    bn_ss: RM,
    bn_tt: RV,
    bn_sn: RV,
    bn_dl: RV,
    // thinb's tt is a plain local that is read before it is set on every
    // page after the first (undefined in the Fortran); it is kept here.
    th_tt: RV,
    /// broadn: bsigma values computed ahead by worker threads (see BsCache)
    bpf: BpfSlot,
}

/// Holder of the broadn prefetch; a copy of the data (for a worker) never
/// carries it.
#[derive(Default)]
struct BpfSlot(Option<Box<BsCache>>);

impl Clone for BpfSlot {
    fn clone(&self) -> Self {
        BpfSlot(None)
    }
}

/// bsigma values that worker threads computed ahead for one broadn call.
///
/// bsigma's result depends only on its arguments (kn, en), the page data
/// (e, s) and the page bounds; it sets every work variable it reads. Its
/// one side effect is the "broadening truncated" message, which depends
/// only on aamin; a value taken from here repeats that message. broadn
/// itself runs unchanged and asks `bsigma_c` for each value, so a value
/// that a worker did not compute is computed by the main program as
/// before. The output is the same with any number of threads.
struct BsCache {
    /// the values of the current panel: (kn, en bits) -> nreac values, aamin
    map: std::collections::HashMap<(i32, u64), Vec<f64>>,
    /// the panel whose values are in `map`
    cur: Option<usize>,
    ra: crate::prefetch::RunAhead<Vec<(i32, u64, Vec<f64>)>>,
    nthr: usize,
    hits: u64,
    misses: u64,
    wait: f64,
}

impl Default for Broadr {
    fn default() -> Self {
        Self::new()
    }
}

impl Broadr {
    pub fn new() -> Self {
        Broadr {
            nunr: 0,
            intunr: 0,
            nsig0: 0,
            ncyc: 0,
            e1: 0.0,
            eumax: 0.0,
            alpha: 0.0,
            tempk: 0.0,
            tempef: 0.0,
            errthn: 0.0,
            errmax: 0.0,
            errint: 0.0,
            thnmax: 0.0,
            klow: 0,
            khigh: 0,
            mlow: 0,
            mhigh: 0,
            nlow: 0,
            nhigh: 0,
            n2in: 0,
            n2out: 0,
            n2left: 0,
            mpage: 0,
            nreac: 0,
            ntx: 0,
            f: [0.0; 6],
            aa: 0.0,
            x: 0.0,
            y: 0.0,
            oy: 0.0,
            xx: 0.0,
            yy: 0.0,
            h: [0.0; 6],
            alast: 0.0,
            s1: 0.0,
            s2: 0.0,
            qmtr: RV::new(NTT),
            emtr: RV::new(NTT),
            bn_ks: [0; NSTACK + 1],
            bn_js: [0; NSTACK + 1],
            bn_es: [0.0; NSTACK + 1],
            bn_ss: RM::new(NTT, NSTACK),
            bn_tt: RV::new(NTT + 1),
            bn_sn: RV::new(NTT),
            bn_dl: RV::new(NTT),
            th_tt: RV::new(NTT + 1),
            bpf: BpfSlot(None),
        }
    }

    // ================================================================ broadr
    #[allow(unused_assignments)]
    pub fn broadr(&mut self, c: &mut Ctx) {
        let therm = 0.0253e0f64;
        let fact = 0.99999e0f64;
        let e6pt5 = 6.5e6f64;
        let tbreak = 0.5e0f64;
        let big = 9.0e9f64;
        let zero = 0.0f64;

        let namax: i32 = 15000000;
        let nwscr: i32 = 2500;
        self.nunr = 0;

        //--initialize and assign i/o units.
        let mut nscr1 = 10;
        let mut iold = 12;
        let mut inew = 13;
        c.openz(-iold, 1);
        c.openz(-inew, 1);

        //--read and display user input.
        let time = c.timer();
        c.io.so(&format!(
            "\n broadr...doppler broadening of endf data{}{}s",
            " ".repeat(28),
            ffmt::f(time, 8, 1)
        ));
        c.io.se(&format!("\n broadr...{}{}s", " ".repeat(59), ffmt::f(time, 8, 1)));
        // locals that a null value leaves undefined in the Fortran are 0
        let v = c.io.input.read(&[Kind::Int, Kind::Int, Kind::Int]);
        let nendf = geti(&v[0], 0);
        let nin = geti(&v[1], 0);
        let nout = geti(&v[2], 0);
        if nin < 0 && nout > 0 {
            c.error("broadr", "nin and nout must be same mode.", " ");
        }
        if nin > 0 && nout < 0 {
            c.error("broadr", "nin and nout must be same mode.", " ");
        }
        if nin < 0 && nout <= 0 {
            nscr1 = -nscr1;
        }
        c.openz(nin, 0);
        c.openz(nendf, 0);
        c.openz(nout, 1);
        c.openz(nscr1, 1);
        let mut istart = 0;
        let mut istrap = 0;
        let mut temp1 = 0.0f64;
        let mut ntemp2 = 1;
        let v = c.io.input.read(&[Kind::Int, Kind::Int, Kind::Int, Kind::Int, Kind::Real]);
        let mut mat1 = geti(&v[0], 0);
        ntemp2 = geti(&v[1], ntemp2);
        istart = geti(&v[2], istart);
        istrap = geti(&v[3], istrap);
        temp1 = getr(&v[4], temp1);
        let mut thnmx = e6pt5;
        self.errmax = 0.0;
        self.errint = 0.0;
        let v = c.io.input.read(&[Kind::Real, Kind::Real, Kind::Real, Kind::Real]);
        self.errthn = getr(&v[0], self.errthn);
        thnmx = getr(&v[1], thnmx);
        self.errmax = getr(&v[2], self.errmax);
        self.errint = getr(&v[3], self.errint);
        if self.errmax == zero {
            self.errmax = 10.0 * self.errthn;
        }
        if self.errint == zero {
            self.errint = self.errthn / 20000.0;
        }
        let mut temp2 = RV::new(ntemp2.max(0) as usize);
        let kinds: Vec<Kind> = (0..ntemp2.max(0)).map(|_| Kind::Real).collect();
        let v = c.io.input.read(&kinds);
        for i in 1..=ntemp2 {
            temp2[i] = getr(&v[(i - 1) as usize], 0.0);
        }
        thnmx = sigfig(thnmx, 7, 0);
        c.io.so(&format!(
            "\n unit for input endf tape ............. {}\n unit for input pendf tape ............ {}\n unit for output pendf tape ........... {}\n material to be processed ............. {}\n number of final temperatures ......... {}\n restart (0 no, 1 yes) ................ {}\n bootstrap (0 no, 1 yes) .............. {}\n starting material temperature ........ {}k\n thinning tolerance ................... {}\n max. energy .......................... {}\n errmax for thinning .................. {}\n errint for thinning .................. {}",
            ffmt::i(nendf, 10),
            ffmt::i(nin, 10),
            ffmt::i(nout, 10),
            ffmt::i(mat1, 10),
            ffmt::i(ntemp2, 10),
            ffmt::i(istart, 10),
            ffmt::i(istrap, 10),
            ffmt::f(temp1, 10, 1),
            ffmt::f(self.errthn, 10, 3),
            ffmt::e1p(thnmx, 10, 3),
            ffmt::e1p(self.errmax, 10, 3),
            ffmt::e1p(self.errint, 10, 3)
        ));
        {
            // '(... ,1p,e10.3/(40x,e10.3))': one value leaves an empty record
            let mut s = format!(" final temperatures ................... {}", ffmt::e1p(temp2[1], 10, 3));
            for i in 2..=ntemp2 {
                s.push('\n');
                s.push_str(&" ".repeat(40));
                s.push_str(&ffmt::e1p(temp2[i], 10, 3));
            }
            if ntemp2 == 1 {
                s.push('\n');
            }
            c.io.so(&s);
        }

        //--allocate storage.
        let mut scr = RV::new(nwscr as usize);
        let mut bufo = vec![0.0f64; NBUF as usize];
        let mut bufn = vec![0.0f64; NBUF as usize];

        // locals of the Fortran routine
        let (mut nb, mut nw) = (0i32, 0i32);
        let mut loop_ = 0;
        let mut tt = RV::new(NTT + 1);
        let mut mtr = IV::new(NTT);
        let mut mti = IV::new(NTT);
        let mut nppmt: Option<Vec<i32>> = None;
        let mut nutot: Option<RV> = None;
        let mut unr = RV::new(1);
        let mut b = RV::new(1);
        let mut e = RV::new(1);
        let mut s = RM::new(1, 1);
        let mut lnu = 0;
        let mut lrp = 0;
        let mut lru = 0;
        let mut npp = 0;
        let mut temp;
        let mut awr = 0.0f64;
        let mut eresh = 0.0f64;
        let mut emin = 0.0f64;
        let mut nwunr = 1;
        let mut tempin = 0.0f64;
        let mut it = 0;
        let mut mt4br = 0;
        let mut mt103 = 0;
        let mut mt104 = 0;
        let mut mt105 = 0;
        let mut mt106 = 0;
        let mut mt107 = 0;
        let (mut mpmin, mut mpmax, mut mdmin, mut mdmax, mut mtmin, mut mtmax) = (0, 0, 0, 0, 0, 0);
        let (mut m3min, mut m3max, mut m4min, mut m4max) = (0, 0, 0, 0);
        // read before set on some paths in the Fortran (stack values there)
        let mut mt4i = 0i32;
        let mut ftev = 0.0f64;
        let mut enext = 0.0f64;
        let mut idis = 0i32;
        let mut sig = 0.0f64;
        let mut en = 0.0f64;
        let mut sun = 0.0f64;

        let mut lbl = 110;
        'main: loop {
            match lbl {
                110 => {
                    //--search input endf tape for some parameters ...
                    c.repoz(nendf);
                    c.tpidio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                    c.findf(mat1, 1, 0, nendf);
                    //--read input file emax and make sure thnmx is legal
                    let mut emax = 20.0e6f64;
                    c.contio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                    lrp = c.l1h;
                    c.contio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                    if c.n1h == 0 && c.n2h != 0 {
                        c.contio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                        if c.c2h > zero {
                            emax = c.c2h;
                        }
                    }
                    if thnmx.abs() > emax {
                        thnmx = 0.9999f64 * emax;
                        c.mess("broadr", "reset thnmx to just less than input file's emax", "");
                    }
                    //--jump to the end of mf1/mt451
                    c.tosend(nendf, 0, 0, &mut scr.all());
                    //--read input file nu-bar, if present
                    lnu = 0;
                    loop {
                        // 101
                        c.contio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                        if c.mfh != 1 {
                            break;
                        }
                        if c.mth != 452 {
                            // 106
                            c.tosend(nendf, 0, 0, &mut scr.all());
                            continue;
                        }
                        lnu = c.l2h;
                        if lnu != 2 {
                            c.listio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                            let nwl = 6 + c.n1h;
                            let mut nt = RV::new(nwl.max(0) as usize);
                            for i in 1..=nwl {
                                nt[i] = scr[i];
                            }
                            nutot = Some(nt);
                            break;
                        }
                        // 103
                        c.tab1io(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                        let nwt = 6 + 2 * c.n1h + 2 * c.n2h;
                        let mut nt = RV::new(nwt.max(0) as usize);
                        let mut l1 = 1;
                        loop {
                            // 104
                            for i in 1..=nw {
                                nt[l1 + i - 1] = scr[i];
                            }
                            if nb == 0 {
                                break;
                            }
                            l1 += nw;
                            c.moreio(nendf, 0, 0, &mut nt.at(l1), &mut nb, &mut nw);
                        }
                        nutot = Some(nt);
                        break;
                    }
                    // 107
                    //--read the (lru) resonance flag
                    npp = 0;
                    if lrp == 0 {
                        lru = 0;
                    } else {
                        if c.mfh <= 1 {
                            if c.mfh == 1 {
                                c.tofend(nendf, 0, 0, &mut scr.all());
                            }
                            c.contio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                        }
                        if c.mfh == 2 {
                            c.contio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                            c.contio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                            lru = c.l1h;
                            let lrf = c.l2h;
                            //--make a list of the lrf7 mt's
                            if lru == 1 && lrf == 7 {
                                c.contio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                                c.listio(nendf, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                                npp = c.l1h;
                                let mut l: Vec<i32> = (1..=npp).map(|i| nint(scr[16 + 12 * (i - 1)])).collect();
                                //--sort the nppmt list into ascending order
                                if npp >= 2 {
                                    for i in 0..(npp - 1) as usize {
                                        for j in (i + 1)..npp as usize {
                                            if l[j] < l[i] {
                                                l.swap(i, j);
                                            }
                                        }
                                    }
                                }
                                nppmt = Some(l);
                            }
                        }
                    }

                    //--search for desired mat1 at temp1 on input tape.
                    if loop_ == 0 {
                        c.nsh = 0;
                        c.repoz(nin);
                        c.repoz(nout);
                        c.tpidio(nin, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                    }
                    lbl = 115;
                }
                115 => {
                    c.contio(nin, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                    if c.math == -1 {
                        // 120
                        c.mess("broadr", "desired mat and temp not on tape", " ");
                        lbl = 400;
                        continue;
                    }
                    if c.math == mat1 {
                        lbl = 130;
                        continue;
                    }
                    c.tomend(nin, 0, 0, &mut scr.all());
                }
                130 => {
                    let mut no = 0;
                    c.nsc = 0;
                    self.thnmax = thnmx;
                    eresh = self.thnmax;
                    emin = 1.0;
                    if self.thnmax < zero {
                        emin = -self.thnmax;
                    }
                    c.repoz(nscr1);
                    if istart == 1 {
                        no = nout;
                    }
                    c.contio(0, no, nscr1, &mut scr.all(), &mut nb, &mut nw);
                    let mut nx = c.n2h;
                    c.contio(nin, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                    if c.n1h != 0 {
                        c.iverf = 4;
                    } else if c.n2h == 0 {
                        c.iverf = 5;
                    } else {
                        c.iverf = 6;
                    }
                    c.io.so(&format!("\n files are in endf-{} format", ffmt::i(c.iverf, 1)));
                    c.skiprz(nin, -1);
                    if c.iverf >= 5 {
                        c.contio(nin, no, nscr1, &mut scr.all(), &mut nb, &mut nw);
                    }
                    if c.iverf >= 6 {
                        c.contio(nin, no, nscr1, &mut scr.all(), &mut nb, &mut nw);
                    }
                    c.hdatio(nin, no, nscr1, &mut scr.all(), &mut nb, &mut nw);
                    temp = scr[1];
                    if c.iverf != 4 {
                        nx = c.n2h;
                    }
                    if no != 0 {
                        c.io.so(&format!("\n output {}k from restart", ffmt::e1p(temp, 12, 4)));
                    }
                    while nb != 0 {
                        c.moreio(nin, no, nscr1, &mut scr.all(), &mut nb, &mut nw);
                    }
                    c.dictio(nin, no, nscr1, &mut scr.all(), &mut nb, &mut nx);
                    nwunr = 1;
                    let mut i = 1;
                    while i <= nx {
                        let mf = nint(scr[i + 2]);
                        let mt = nint(scr[i + 3]);
                        let ncds = nint(scr[i + 4]);
                        if mf == 2 && mt == 152 {
                            nwunr = 6 * ncds;
                        }
                        i += 6;
                    }
                    c.tofend(nin, no, nscr1, &mut scr.all());
                    unr = RV::new(nwunr.max(1) as usize);
                    c.contio(nin, no, nscr1, &mut scr.all(), &mut nb, &mut nw);
                    if c.mfh == 2 {
                        c.contio(nin, no, nscr1, &mut scr.all(), &mut nb, &mut nw);
                        c.contio(nin, no, nscr1, &mut scr.all(), &mut nb, &mut nw);
                        eresh = c.c2h;
                        c.tosend(nin, no, nscr1, &mut scr.all());
                        c.contio(nin, no, nscr1, &mut scr.all(), &mut nb, &mut nw);
                        if c.mfh != 0 {
                            if c.mth == 152 {
                                let mut th = self.thnmax;
                                self.stounx(c, nin, no, nscr1, &mut unr, &mut th);
                                self.thnmax = th;
                            }
                            c.tofend(nin, no, nscr1, &mut scr.all());
                        }
                        c.contio(nin, no, nscr1, &mut scr.all(), &mut nb, &mut nw);
                    }
                    c.tomend(nin, no, nscr1, &mut scr.all());
                    let diff = (temp - temp1).abs();
                    let test = 1.0 + temp1 / 1000.0;
                    if diff > test {
                        lbl = 115;
                        continue;
                    }

                    if lrp == 0 {
                        c.io.so(&format!(
                            "\n non-resonance nuclide, input pendf limit{}= {} eV.",
                            " ".repeat(5),
                            ffmt::e1p(eresh, 14, 5)
                        ));
                    } else if lrp == 1 && lru == 1 {
                        c.io.so(&format!(
                            "\n resolved resonance range upper limit{}= {} eV.",
                            " ".repeat(9),
                            ffmt::e1p(eresh, 14, 5)
                        ));
                    } else if lrp == 1 && lru == 2 {
                        c.io.so(&format!(
                            "\n unresolved resonance range lower limit{}= {} eV.",
                            " ".repeat(7),
                            ffmt::e1p(eresh, 14, 5)
                        ));
                    }
                    if eresh > e6pt5 {
                        eresh = e6pt5;
                        c.io.so(&format!(" - reset to {} eV.", ffmt::e1p(e6pt5, 14, 5)));
                    }
                    if thnmx == 0.0 {
                        self.thnmax = eresh;
                    } else if thnmx > 0.0 {
                        self.thnmax = thnmx.min(eresh);
                    }
                    if lrp == 1 {
                        emin = eresh;
                    }

                    //--copy low threshold reactions from nscr1 to
                    //--scratch storage using the grid of mt1.
                    c.repoz(nscr1);
                    c.tofend(nscr1, 0, 0, &mut scr.all());
                    c.contio(nscr1, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                    c.contio(nscr1, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                    c.contio(nscr1, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                    if c.c2h < self.thnmax {
                        self.thnmax = c.c2h;
                    }
                    c.tofend(nscr1, 0, 0, &mut scr.all());
                    c.contio(nscr1, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                    awr = c.c2h;
                    let mut enow = 0.0f64;
                    c.gety1(enow, &mut enext, &mut idis, &mut sig, nscr1, &mut scr.all());
                    self.nreac = 0;
                    self.ntx = 1;
                    self.n2in = 0;
                    let mut idnx = 0;
                    mt4br = 0;
                    mt103 = 0;
                    mt104 = 0;
                    mt105 = 0;
                    mt106 = 0;
                    mt107 = 0;
                    if c.iverf >= 6 {
                        mpmin = 600;
                        mpmax = 649;
                        mdmin = 650;
                        mdmax = 699;
                        mtmin = 700;
                        mtmax = 749;
                        m3min = 750;
                        m3max = 799;
                        m4min = 800;
                        m4max = 849;
                    } else {
                        mpmin = 700;
                        mpmax = 718;
                        mdmin = 720;
                        mdmax = 738;
                        mtmin = 740;
                        mtmax = 758;
                        m3min = 760;
                        m3max = 768;
                        m4min = 780;
                        m4max = 798;
                    }
                    enow = enext;
                    while enext < big {
                        self.n2in += 1;
                        if idnx > 0 && enext * fact > enow {
                            enext = enext * fact;
                        }
                        enow = enext;
                        c.gety1(enow, &mut enext, &mut idis, &mut sig, nscr1, &mut scr.all());
                        let mth = c.mth;
                        self.getunx(c, mth, enow, &mut en, &mut sun, &unr);
                        if en > zero {
                            sig = sig - sun;
                        }
                        tt[1] = enow;
                        let mut kt = self.n2in;
                        if enext >= big {
                            kt = -kt;
                        }
                        c.loada(kt, &tt.0, self.ntx, iold, &mut bufo, NBUF);
                        if idis > 0 {
                            idnx = 1 - idnx;
                        }
                        if idis == 0 {
                            idnx = 0;
                        }
                    }
                    c.tosend(nscr1, 0, 0, &mut scr.all());
                    let mut n = 1;
                    loop {
                        // 160
                        c.contio(nscr1, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                        if c.mfh == 0 {
                            break;
                        }
                        let mth = c.mth;
                        let mut skip = false;
                        if c.iverf < 6 {
                            if mth > 150 {
                                skip = true;
                            }
                        } else {
                            if mth > 200 && mth < 600 {
                                skip = true;
                            }
                            if mth > 850 {
                                skip = true;
                            }
                        }
                        if mth == 3 || mth == 4 || (46..=49).contains(&mth) || mth == 19 {
                            skip = true;
                        }
                        let mut take = false;
                        if !skip {
                            let enow0 = 0.0f64;
                            c.gety1(enow0, &mut enext, &mut idis, &mut sig, nscr1, &mut scr.all());
                            if mth == 18 {
                                take = true;
                            } else if let Some(pm) = &nppmt {
                                if n <= npp && mth == pm[(n - 1) as usize] {
                                    n += 1;
                                    take = true;
                                }
                            }
                            if !take {
                                if enext <= emin {
                                    take = true;
                                } else if lrp == 1 && enext < eresh {
                                    take = true;
                                } else if lrp == 0 && fact * enext < self.thnmax {
                                    self.thnmax = fact * enext;
                                }
                            }
                        }
                        if !take {
                            // 165
                            c.tosend(nscr1, 0, 0, &mut scr.all());
                            continue;
                        }
                        // 170
                        self.nreac += 1;
                        if self.nreac >= NTT as i32 {
                            c.error("broadr", "too many low-threshold reactions.", " ");
                        }
                        self.ntx = self.nreac + 1;
                        let nreac = self.nreac;
                        mtr[nreac] = mth;
                        if (51..=91).contains(&mth) {
                            mt4br = 1;
                        }
                        if mth == 103 {
                            mt103 = 1;
                        }
                        if mth == 104 {
                            mt104 = 1;
                        }
                        if mth == 105 {
                            mt105 = 1;
                        }
                        if mth == 106 {
                            mt106 = 1;
                        }
                        if mth == 107 {
                            mt107 = 1;
                        }
                        mti[nreac] = 0;
                        let mut kt = 1;
                        let mut k = 0;
                        while kt > 0 {
                            k += 1;
                            c.finda(k, &mut tt.0, nreac, iold, &mut bufo, NBUF);
                            let enow = tt[1];
                            c.gety1(enow, &mut enext, &mut idis, &mut sig, nscr1, &mut scr.all());
                            self.getunx(c, mth, enow, &mut en, &mut sun, &unr);
                            if en > zero {
                                sig = sig - sun;
                            }
                            tt[1 + nreac] = sig;
                            kt = k;
                            if k == self.n2in {
                                kt = -k;
                            }
                            c.loada(kt, &tt.0, self.ntx, inew, &mut bufn, NBUF);
                        }
                        c.tosend(nscr1, 0, 0, &mut scr.all());
                        std::mem::swap(&mut iold, &mut inew);
                    }
                    // 190
                    if self.thnmax == zero && thnmx == 0.0 {
                        self.thnmax = e6pt5;
                    }
                    nppmt = None;
                    if self.thnmax < zero {
                        self.thnmax = -self.thnmax;
                    }
                    c.io.so(&format!(
                        "\n final maximum energy for broadening/thinning = {} eV",
                        ffmt::e1p(self.thnmax, 14, 5)
                    ));
                    tempin = temp1;
                    if self.nreac > 0 {
                        lbl = 200;
                        continue;
                    }
                    c.mess("broadr", "no broadenable reactions", " ");
                    it = 0;
                    lbl = 195;
                }
                195 => {
                    it += 1;
                    self.tempk = temp2[it];
                    c.io.so(&format!(
                        " copied mat{} from{} to{} k",
                        ffmt::i(mat1, 4),
                        ffmt::e1p(tempin, 12, 4),
                        ffmt::e1p(self.tempk, 12, 4)
                    ));
                    lbl = 215;
                }
                200 => {
                    //--assign storage for paging cross sections.
                    let nwb = c.npage + 50;
                    b = RV::new(nwb as usize);
                    self.mpage = namax / (4 * (self.nreac + 1));
                    let ndim = 3 * self.mpage;
                    // one slack element: bsigma looks at e(k+1) for k=khigh
                    e = RV::new(ndim as usize + 1);
                    s = RM::new(self.nreac as usize, ndim as usize + 1);

                    //--get q-values and lab system threshold energy
                    //  for the mt's to be broadened.
                    c.repoz(nscr1);
                    for i in 1..=self.nreac {
                        c.findf(mat1, 3, mtr[i], nscr1);
                        c.contio(nscr1, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                        c.contio(nscr1, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                        self.qmtr[i] = scr[2];
                        if self.qmtr[i] == 0.0 {
                            self.emtr[i] = 0.0;
                        } else {
                            self.emtr[i] = -self.qmtr[i] * (awr + 1.0) / awr;
                        }
                    }

                    //--doppler broaden and thin these reactions.
                    it = 0;
                    tempin = temp1;
                    lbl = 210;
                }
                210 => {
                    it += 1;
                    self.tempk = temp2[it];
                    self.tempef = self.tempk - tempin;
                    if self.tempef > zero {
                        self.alpha = awr / (BK * self.tempef);
                    }
                    let time = c.timer();
                    c.io.so(&format!("{}{}s", " ".repeat(69), ffmt::f(time, 8, 1)));
                    let brk = 9999.0f64;
                    if ntemp2 > 1 && self.tempk < brk {
                        c.io.se(&format!("{} deg{}{}s", ffmt::f(self.tempk, 11, 1), " ".repeat(54), ffmt::f(time, 8, 1)));
                    } else if ntemp2 > 1 {
                        c.io.se(&format!(
                            "{} deg ({} eV){}{}s",
                            ffmt::e1p(self.tempk, 11, 3),
                            ffmt::e1p(BK * self.tempk, 9, 3),
                            " ".repeat(39),
                            ffmt::f(time, 8, 1)
                        ));
                    }
                    let nreac = self.nreac;
                    self.bfile3(c, iold, inew, &mut bufo, &mut bufn, &mut e, &mut s);
                    {
                        let mut t = format!(
                            "\n broadened mat{} from{} to{} k\n         points in={}  points out={}\n         mt",
                            ffmt::i(mat1, 4),
                            ffmt::e1p(tempin, 12, 4),
                            ffmt::e1p(self.tempk, 12, 4),
                            ffmt::i(self.n2in, 7),
                            ffmt::i(self.n2out, 7)
                        );
                        // 16i4,(/11x,16i4): values 17-32 start a new record; each
                        // later group of 16 follows an empty record (format
                        // reversion ends a record, then '/' ends another);
                        // exactly 16 values leave a trailing empty record
                        for i in 1..=nreac {
                            if i > 16 && (i - 1) % 16 == 0 {
                                if i > 17 {
                                    t.push('\n');
                                }
                                t.push('\n');
                                t.push_str(&" ".repeat(11));
                            }
                            t.push_str(&ffmt::i(mtr[i], 4));
                        }
                        if nreac == 16 {
                            t.push('\n');
                        }
                        c.io.so(&t);
                    }
                    let mut nset = 0;
                    let mut j = 0;
                    while nset < nreac {
                        j += 1;
                        c.finda(j, &mut tt.0, self.ntx, inew, &mut bufn, NBUF);
                        nset = 0;
                        for k in 1..=nreac {
                            if j == 1 {
                                mti[k] = 0;
                            }
                            if mti[k] == 0 && tt[1 + k] > zero {
                                mti[k] = j;
                            }
                            if mti[k] > 0 {
                                nset += 1;
                            }
                        }
                    }
                    lbl = 215;
                }
                215 => {
                    //--correct dictionary to reflect thinning.
                    c.repoz(nscr1);
                    c.contio(nscr1, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                    let mut nxc = c.n2h;
                    if c.iverf >= 5 {
                        c.contio(nscr1, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                    }
                    if c.iverf == 6 {
                        c.contio(nscr1, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                    }
                    c.hdatio(nscr1, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                    if c.iverf != 4 {
                        nxc = c.n2h;
                    }
                    scr[1] = self.tempk;
                    c.hdatio(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                    while nb != 0 {
                        c.moreio(nscr1, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                    }
                    nw = nxc;
                    if nxc * 6 > nwscr {
                        c.error("broadr", "storage exceeded.", " ");
                    }
                    c.dictio(nscr1, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                    if self.nreac != 0 {
                        let mut i = 1;
                        while i <= nw {
                            let mf = nint(scr[i + 2]);
                            if mf == 3 {
                                let mt = nint(scr[i + 3]);
                                if mt == 1 || mt == 3 || mt == 19 {
                                    scr[i + 4] = (3 + (self.n2out + 2) / 3) as f64;
                                }
                                if mt == 4 && mt4br != 0 {
                                    for j in 1..=self.nreac {
                                        let mut jjj = mti[j];
                                        if jjj > 1 {
                                            jjj -= 1;
                                        }
                                        if mtr[j] == 51 {
                                            mt4i = jjj;
                                        }
                                    }
                                    scr[i + 4] = (3 + (self.n2out - mt4i + 3) / 3) as f64;
                                }
                                for j in 1..=self.nreac {
                                    let mut jjj = mti[j];
                                    if jjj > 1 {
                                        jjj -= 1;
                                    }
                                    if mt == mtr[j] {
                                        scr[i + 4] = (3 + (self.n2out - jjj + 3) / 3) as f64;
                                    }
                                }
                            }
                            i += 6;
                        }
                    }
                    nw = nxc;
                    c.dictio(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                    c.tofend(nscr1, nout, 0, &mut scr.all());

                    //--write out thermal quantities if t is close to 293.6 Kelvin.
                    let eone = tbreak;
                    let tev = self.tempk * BK;
                    if (tev - therm).abs() <= therm / 1000.0 {
                        self.thermal(c, tev, therm, eone, lnu, &nutot, &mut tt, &mtr, inew, &mut bufn, &mut ftev);
                    }
                    lbl = 245;
                }
                245 => {
                    //--output the broadened cross sections in endf format.
                    c.contio(nscr1, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                    if c.mfh == 0 {
                        lbl = 320;
                        continue;
                    }
                    if c.mfh == 3 {
                        lbl = 250;
                        continue;
                    }
                    if !(c.mfh == 2 && c.mth == 152) {
                        c.tosend(nscr1, nout, 0, &mut scr.all());
                        c.contio(nscr1, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                        if c.mfh == 0 {
                            continue; // 245
                        }
                    }
                    // 247: update unresolved temperature
                    c.listio(nscr1, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                    scr[1] = self.tempk;
                    c.listio(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                    c.tofend(nscr1, nout, 0, &mut scr.all());
                }
                250 => {
                    let mth = c.mth;
                    let mut ireac = 0;
                    let mut broadened = false;
                    if self.nreac != 0 {
                        if mth == 1 || mth == 3 || (mth == 4 && mt4br == 1) || (46..=49).contains(&mth) || mth == 19 {
                            broadened = true;
                        } else {
                            for i in 1..=self.nreac {
                                if mth == mtr[i] {
                                    ireac = i;
                                }
                            }
                        }
                    }
                    if ireac > 0 {
                        broadened = true;
                    }
                    if !broadened {
                        // this section wasnt broadened or thinned.
                        // copy it to output.
                        c.tab1io(nscr1, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                        let mut nbo = nb;
                        let mut nwo = nw;
                        c.tab1io(0, nout, 0, &mut scr.all(), &mut nbo, &mut nwo);
                        while nb != 0 {
                            c.moreio(nscr1, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                        }
                    } else {
                        self.write_section(
                            c, mth, ireac, nscr1, nout, &mut scr, &mut b, &mut tt, &mtr, &mti, inew, &mut bufn, &unr,
                            [mt103, mt104, mt105, mt106, mt107],
                            [mpmin, mpmax, mdmin, mdmax, mtmin, mtmax, m3min, m3max, m4min, m4max],
                        );
                    }
                    // 310
                    c.tosend(nscr1, 0, 0, &mut scr.all());
                    c.asend(nout, 0);
                    lbl = 245;
                }
                320 => {
                    c.tomend(nscr1, nout, 0, &mut scr.all());

                    //--loop over requested final temperatures.
                    if it == ntemp2 {
                        lbl = 400;
                        continue;
                    }
                    if istrap != 0 {
                        std::mem::swap(&mut iold, &mut inew);
                        self.n2in = self.n2out;
                        tempin = self.tempk;
                    }
                    lbl = if self.nreac == 0 { 195 } else { 210 };
                }
                400 => {
                    //--loop over desired materials.
                    let v = c.io.input.read(&[Kind::Int]);
                    mat1 = geti(&v[0], mat1);
                    if mat1 == 0 {
                        break 'main;
                    }
                    nutot = None;
                    let _ = lnu;
                    loop_ += 1;
                    lbl = 110;
                }
                _ => unreachable!(),
            }
        }

        //--broadr is finished
        c.atend(nout, 0);
        c.repoz(nout);
        c.repoz(nin);
        c.closz(nin);
        c.closz(nendf);
        c.closz(nout);
        c.closz(-iold);
        c.closz(-inew);
        c.closz(nscr1);
        let time = c.timer();
        c.io.so(&format!("{}{}s\n {}*******", " ".repeat(69), ffmt::f(time, 8, 1), "**********".repeat(7)));
        let _ = (e, s, b, unr, sun, en);
    }

    /// Part of `broadr` (label 220 to 240): thermal quantities near 293.6 K.
    #[allow(clippy::too_many_arguments)]
    fn thermal(
        &mut self,
        c: &mut Ctx,
        tev: f64,
        therm: f64,
        eone: f64,
        lnu: i32,
        nutot: &Option<RV>,
        tt: &mut RV,
        mtr: &IV,
        inew: i32,
        bufn: &mut [f64],
        ftev: &mut f64,
    ) {
        let zero = 0.0f64;
        let step = 1.01e0f64;
        c.io.so(&format!(
            "\n    thermal quantities at{} K ={} eV\n    -----------------------------------------",
            ffmt::f(self.tempk, 6, 1),
            ffmt::f(tev, 7, 4)
        ));
        let picon = 2.0 / sqrt(PI);
        let mut llf = 0;
        let mut llc = 0;
        let mut sc = 0.0f64;
        let mut sf = 0.0f64;
        for j in 1..=self.nreac {
            if mtr[j] == 18 {
                llf = j + 1;
            }
            if mtr[j] == 102 {
                llc = j + 1;
            }
        }
        let mut fint = 0.0f64;
        let mut cint = 0.0f64;
        let mut alint = 0.0f64;
        let mut etint = 0.0f64;
        let mut v1int = 0.0f64;
        let mut ssf = 0.0f64;
        let mut slf = 0.0f64;
        let mut ssc = 0.0f64;
        let mut slc = 0.0f64;
        let mut elast = 0.0f64;
        let mut flast = 0.0f64;
        let mut fnul = 0.0f64;
        let mut ctev;
        let mut fnu = 0.0f64;
        let mut xnext = 0.0f64;
        let mut idis = 0i32;
        let empty = RV::new(1);
        let nut = nutot.as_ref().unwrap_or(&empty);
        let nu_at = |c: &Ctx, x: f64, fnu: &mut f64, xnext: &mut f64, idis: &mut i32| {
            if lnu == 1 {
                *fnu = nut[7];
            } else {
                let mut ir = 1;
                let mut ip = 2;
                let mut v = nut.0.clone();
                c.terpa(fnu, x, xnext, idis, &crate::farr::F(&mut v[..]), &mut ip, &mut ir);
            }
        };
        let mut j = 0;
        loop {
            // 220
            j += 1;
            c.finda(j, &mut tt.0, self.ntx, inew, bufn, NBUF);
            loop {
                // 225
                let mut enow = tt[1];
                if enow > step * elast && elast > zero {
                    enow = step * elast;
                    if llf > 0 {
                        sf = slf + (tt[llf] - slf) * (enow - elast) / (tt[1] - elast);
                    }
                    if llc > 0 {
                        sc = slc + (tt[llc] - slc) * (enow - elast) / (tt[1] - elast);
                    }
                } else {
                    if llf > 0 {
                        sf = tt[llf];
                    }
                    if llc > 0 {
                        sc = tt[llc];
                    }
                }
                // thermal cross sections at 0.0253 eV
                if elast < therm && enow >= therm {
                    if llf > 0 {
                        let ss = slf + (therm - elast) * (sf - slf) / (enow - elast);
                        c.io.so(&format!("           thermal fission xsec:{}", ffmt::e1p(ss, 12, 4)));
                        nu_at(c, therm, &mut fnu, &mut xnext, &mut idis);
                        c.io.so(&format!("          thermal fission nubar:{}", ffmt::e1p(fnu, 12, 4)));
                    }
                    if llc > 0 {
                        let ss = slc + (therm - elast) * (sc - slc) / (enow - elast);
                        c.io.so(&format!("           thermal capture xsec:{}", ffmt::e1p(ss, 12, 4)));
                    }
                }
                // thermal cross sections at tev
                ctev = zero;
                if elast < tev && enow >= tev {
                    if llf > 0 {
                        let ss = slf + (tev - elast) * (sf - slf) / (enow - elast);
                        *ftev = ss;
                        nu_at(c, tev, &mut fnu, &mut xnext, &mut idis);
                    }
                    if llc > 0 {
                        let ss = slc + (tev - elast) * (sc - slc) / (enow - elast);
                        ctev = ss;
                    }
                }
                // fission nubar at this energy
                fnu = 0.0;
                if llf > 0 {
                    nu_at(c, enow, &mut fnu, &mut xnext, &mut idis);
                }
                // thermal integrals
                let brk = 1.0f64;
                if enow < brk {
                    let fnow = enow * exp(-enow / tev);
                    if llc > 0 && j > 1 {
                        cint = cint + (fnow * sc + flast * slc) * (enow - elast) / 2.0;
                    }
                    if llf > 0 && j > 1 {
                        fint = fint + (fnow * sf + flast * slf) * (enow - elast) / 2.0;
                        if sf > zero && slf > zero {
                            alint = alint + (fnow * sc / sf + flast * slc / slf) * (enow - elast) / 2.0;
                        }
                        if (slf + slc) != 0.0 && (sf + sc) != 0.0 {
                            etint = etint
                                + (fnow * fnu * sf / (sf + sc) + flast * fnul * slf / (slf + slc)) * (enow - elast) / 2.0;
                        }
                        v1int = v1int
                            + (fnow * (fnu * sf - sf - sc) + flast * (fnul * slf - slf - slc)) * (enow - elast) / 2.0;
                    }
                    flast = fnow;
                }
                // resonance integrals
                if elast < eone && enow >= eone {
                    if llf > 0 {
                        let ss = slf + (eone - elast) * (sf - slf) / (enow - elast);
                        ssf = (sf / enow + ss / elast) * (enow - eone) / 2.0;
                    }
                    if llc > 0 {
                        let ss = slc + (eone - elast) * (sc - slc) / (enow - elast);
                        ssc = (sc / enow + ss / elast) * (enow - eone) / 2.0;
                    }
                }
                if elast >= eone {
                    if llf > 0 {
                        ssf = ssf + (sf / enow + slf / elast) * (enow - elast) / 2.0;
                    }
                    if llc > 0 {
                        ssc = ssc + (sc / enow + slc / elast) * (enow - elast) / 2.0;
                    }
                }
                // continue the energy loop
                elast = enow;
                fnul = fnu;
                if llf > 0 {
                    slf = sf;
                }
                if llc > 0 {
                    slc = sc;
                }
                if enow < tt[1] {
                    continue; // 225
                }
                break;
            }
            if j >= self.n2out {
                break; // 230
            }
        }
        // 230: print out the results
        let tev2 = powi(tev, 2);
        cint = cint / tev2;
        c.io.so(&format!("       thermal capture integral:{}", ffmt::e1p(cint, 12, 4)));
        if ctev != zero {
            c.io.so(&format!("       thermal capture g-factor:{}", ffmt::e1p(picon * cint / ctev, 12, 4)));
        }
        c.io.so(&format!("     capture resonance integral:{}", ffmt::e1p(ssc, 12, 4)));
        if llf > 0 {
            fint = fint / tev2;
            alint = alint / tev2;
            etint = etint / tev2;
            v1int = v1int / tev2;
            c.io.so(&format!("       thermal fission integral:{}", ffmt::e1p(fint, 12, 4)));
            if *ftev != zero {
                c.io.so(&format!("       thermal fission g-factor:{}", ffmt::e1p(picon * fint / *ftev, 12, 4)));
            }
            c.io.so(&format!("         thermal alpha integral:{}", ffmt::e1p(alint, 12, 4)));
            c.io.so(&format!("           thermal eta integral:{}", ffmt::e1p(etint, 12, 4)));
            c.io.so(&format!("            thermal k1 integral:{}", ffmt::e1p(v1int, 12, 4)));
            c.io.so(&format!("                  equivalent k1:{}", ffmt::e1p(picon * v1int, 12, 4)));
            c.io.so(&format!("     fission resonance integral:{}", ffmt::e1p(ssf, 12, 4)));
        }
        c.io.so("    -----------------------------------------");
    }

    /// Part of `broadr` (label 275 to 300): write one broadened MF3 section.
    #[allow(clippy::too_many_arguments, unused_assignments)]
    fn write_section(
        &mut self,
        c: &mut Ctx,
        mth: i32,
        ireac: i32,
        nscr1: i32,
        nout: i32,
        scr: &mut RV,
        b: &mut RV,
        tt: &mut RV,
        mtr: &IV,
        mti: &IV,
        inew: i32,
        bufn: &mut [f64],
        unr: &RV,
        mt10x: [i32; 5],
        lim: [i32; 10],
    ) {
        let zero = 0.0f64;
        let [mt103, mt104, mt105, mt106, mt107] = mt10x;
        let [mpmin, mpmax, mdmin, mdmax, mtmin, mtmax, m3min, m3max, m4min, m4max] = lim;
        let (mut nb, mut nw) = (0i32, 0i32);
        let mut enext = 0.0f64;
        let mut idis = 0i32;
        let mut sig = 0.0f64;
        let (mut en, mut sun) = (0.0f64, 0.0f64);
        // 275
        let ee = 0.0f64;
        c.gety1(ee, &mut enext, &mut idis, &mut sig, nscr1, &mut b.all());
        let mut jjj = 1;
        for j in 1..=self.nreac {
            if mth == mtr[j] {
                jjj = mti[j];
            }
            if mth == 4 && mtr[j] == 51 {
                jjj = mti[j];
            }
            if mth == 4 && jjj == 1 && mtr[j] == 91 {
                jjj = mti[j];
            }
        }
        if jjj > 1 {
            jjj -= 1;
        }
        let nnn = self.n2out - jjj + 1;
        scr[1] = c.c1h;
        scr[2] = c.c2h;
        scr[3] = c.l1h as f64;
        scr[4] = c.l2h as f64;
        scr[5] = 1.0;
        scr[6] = nnn as f64;
        scr[7] = nnn as f64;
        scr[8] = 2.0;
        let mut k = 8;
        let j1 = c.npage / 2;
        let mut j = 0;
        loop {
            // 305
            j += 1;
            c.finda(j, &mut tt.0, self.ntx, inew, bufn, NBUF);
            if j < jjj {
                continue;
            }
            k += 1;
            scr[k] = tt[1];
            k += 1;
            if ireac != 0 {
                scr[k] = tt[1 + ireac];
                let e0 = scr[k - 1];
                self.getunx(c, mth, e0, &mut en, &mut sun, unr);
                if en > zero {
                    scr[k] = scr[k] + sun;
                }
            } else if tt[1] > self.thnmax {
                let mut y = 0.0f64;
                c.gety1(tt[1], &mut enext, &mut idis, &mut y, nscr1, &mut b.all());
                scr[k] = y;
            } else {
                scr[k] = 0.0;
                for i in 1..=self.nreac {
                    let mut iflag = 0;
                    if mth == 3 && mtr[i] == 2 {
                        iflag = 1;
                    }
                    if mth == 4 && (mtr[i] < 51 || mtr[i] > 91) {
                        iflag = 1;
                    }
                    if mth == 19 && mtr[i] != 18 {
                        iflag = 1;
                    }
                    if (46..=49).contains(&mth) && mtr[i] != mth - 40 {
                        iflag = 1;
                    }
                    if mth >= 201 {
                        iflag = 1;
                    }
                    // Don't include partial xs if its sum is already available
                    if mt103 == 1 && mtr[i] >= mpmin && mtr[i] <= mpmax {
                        iflag = 1;
                    }
                    if mt104 == 1 && mtr[i] >= mdmin && mtr[i] <= mdmax {
                        iflag = 1;
                    }
                    if mt105 == 1 && mtr[i] >= mtmin && mtr[i] <= mtmax {
                        iflag = 1;
                    }
                    if mt106 == 1 && mtr[i] >= m3min && mtr[i] <= m3max {
                        iflag = 1;
                    }
                    if mt107 == 1 && mtr[i] >= m4min && mtr[i] <= m4max {
                        iflag = 1;
                    }
                    if iflag == 0 {
                        scr[k] = scr[k] + tt[1 + i];
                    }
                }
                self.getunx(c, 1, tt[1], &mut en, &mut sun, unr);
                if en > zero {
                    scr[k] = scr[k] + sun;
                }
            }
            scr[k] = sigfig(scr[k], 7, 0);
            if j - jjj + 1 > j1 {
                // 300
                if k != c.npage && j != self.n2out {
                    continue;
                }
                nw = k;
                c.moreio(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                k = 0;
                if j < self.n2out {
                    continue;
                }
                break;
            }
            if j - jjj + 1 < j1 && j != self.n2out {
                continue;
            }
            nw = k;
            c.tab1io(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
            k = 0;
            if j < self.n2out {
                continue;
            }
            break;
        }
    }

    // ================================================================ stounx
    fn stounx(&mut self, c: &mut Ctx, nin: i32, nout: i32, nscr: i32, unr: &mut RV, thnmax: &mut f64) {
        let zero = 0.0f64;
        let (mut nb, mut nw) = (0i32, 0i32);
        //--read in the data
        self.nunr = 0;
        self.intunr = c.n2h;
        let mut l = 1;
        c.listio(nin, 0, 0, &mut unr.at(l), &mut nb, &mut nw);
        let nx = nint(unr[l + 2]);
        self.nsig0 = nint(unr[l + 3]);
        self.ncyc = 1 + nx * self.nsig0;
        self.nunr = c.n2h;
        while nb != 0 {
            l += nw;
            c.moreio(nin, 0, 0, &mut unr.at(l), &mut nb, &mut nw);
        }

        //--squeeze out sigmazero data if present
        if self.nsig0 > 1 {
            c.mess("stounx", "sigma zero data removed from unresolved", "regenerate it if needed.");
            let mut l = 8;
            let mut k = l + self.nsig0 - 1;
            for _i in 1..=self.nunr {
                unr[l] = unr[k];
                l += 1;
                k += self.nsig0;
                for _j in 1..=nx {
                    unr[l] = unr[k];
                    l += 1;
                    k += self.nsig0;
                }
            }
            self.nsig0 = 1;
            self.ncyc = 1 + nx * self.nsig0;
            let np = 6 + self.nsig0 + self.nunr * (1 + nx);
            unr[4] = self.nsig0 as f64;
            unr[6] = np as f64;
        }

        //--write mt152 to output files
        let mut l = 1;
        c.listio(0, nout, nscr, &mut unr.at(l), &mut nb, &mut nw);
        while nb != 0 {
            l += nw;
            c.moreio(0, nout, nscr, &mut unr.at(l), &mut nb, &mut nw);
        }
        let mut l = 7 + self.nsig0;
        for i in 1..=self.nunr {
            let e = unr[l];
            if i == 1 {
                self.e1 = e.abs();
            }
            if i == 1 || e < zero {
                self.eumax = e.abs();
            }
            if e < zero {
                unr[l] = -e;
            }
            l += self.ncyc;
        }
        if *thnmax > 0.0 && self.eumax < *thnmax {
            *thnmax = self.eumax;
        }
        if self.eumax > self.e1 {
            return;
        }
        //--no overlap.  unresolved not needed.
        self.nunr = 0;
    }

    // ================================================================ getunx
    fn getunx(&self, c: &Ctx, mt: i32, e: f64, enext: &mut f64, sun: &mut f64, unr: &RV) {
        let up = 1.00001e0f64;
        *enext = -1.0;
        *sun = 0.0;
        if self.nunr == 0 {
            return;
        }
        if e < self.e1 {
            return;
        }
        if e > self.eumax {
            return;
        }
        let mut ix = 0;
        if mt == 1 {
            ix = 1;
        }
        if mt == 2 {
            ix = 2;
        }
        if mt == 18 || mt == 19 {
            ix = 3;
        }
        if mt == 102 {
            ix = 4;
        }
        if ix == 0 {
            return;
        }
        let ibase = 7 + self.nsig0;
        let lim = self.nunr - 1;
        // loc1 is undefined in the Fortran if lim < 1; 0 here
        let mut loc1 = 0;
        for i in 1..=lim {
            let loc = ibase + (i - 1) * self.ncyc;
            loc1 = loc + self.ncyc;
            if e >= unr[loc] && e < unr[loc1] {
                c.terp1(unr[loc], unr[loc + ix], unr[loc1], unr[loc1 + ix], e, sun, self.intunr);
            }
        }
        if *sun == 0.0 {
            *sun = unr[loc1 + ix];
            *enext = self.eumax * up;
        } else {
            *enext = unr[loc1 + 1];
        }
    }

    // ================================================================ bfile3
    #[allow(clippy::too_many_arguments)]
    fn bfile3(&mut self, c: &mut Ctx, iold: i32, inew: i32, bufo: &mut [f64], bufn: &mut [f64], e: &mut RV, s: &mut RM) {
        let mut tt = RV::new(NTT + 1);
        //--set broadening indices to middle page
        self.n2left = self.n2in;
        self.nlow = self.mpage + 1;
        self.nhigh = 2 * self.mpage;
        if self.n2left < self.mpage {
            self.nhigh = self.nlow + self.n2left - 1;
        }
        self.klow = self.nlow;
        self.mlow = self.nlow;
        self.mhigh = 3 * self.mpage;
        let mut nj = 0i32;
        let mut l = 0;
        let nreac = self.nreac;

        let mut load = true;
        loop {
            if load {
                // 110: load desired points
                if self.n2left < self.mhigh - self.mlow + 1 {
                    self.mhigh = self.mlow + self.n2left - 1;
                }
                for k in self.mlow..=self.mhigh {
                    l += 1;
                    c.finda(l, &mut tt.0, self.ntx, iold, bufo, NBUF);
                    e[k] = tt[1];
                    for i in 1..=nreac {
                        s[(i, k)] = tt[i + 1];
                    }
                }
                self.khigh = self.mhigh;
                self.n2left -= self.mhigh - self.mlow + 1;
            }
            // 120: broaden or thin points in middle page
            if self.tempef > 0.0 {
                self.broadn(c, &mut nj, inew, bufn, e, s);
            } else {
                self.thinb(c, &mut nj, inew, bufn, e, s);
            }
            //--has last point been broadened or thinned.
            if nj < 0 {
                break;
            }
            //--shift remaining unbroadened points one page forward
            let nw = self.khigh - self.nlow + 1;
            let kk = nw;
            for k in 1..=nw {
                e[k] = e[self.nlow + k - 1];
                for i in 1..=nreac {
                    s[(i, k)] = s[(i, self.nlow + k - 1)];
                }
            }
            self.klow = 1;
            self.mlow = kk + 1;
            if self.n2left > 0 {
                load = true;
                continue;
            }
            self.khigh = kk;
            self.nhigh = self.khigh;
            self.mhigh = self.khigh;
            load = false;
        }
        // 150
        self.n2out = -nj;
    }

    // ================================================================ broadn
    #[allow(unused_assignments)]
    fn broadn(&mut self, c: &mut Ctx, j: &mut i32, inew: i32, bufn: &mut [f64], e: &mut RV, s: &mut RM) {
        let therm = 0.0253e0f64;
        let one = 1.0e0f64;
        let estp = 4.1e0f64;
        let half = 0.5e0f64;
        let small = 1.0e-9f64;
        let step = 2.01e0f64;
        let rmax = 3.0e0f64;
        let errmin = 1.0e-15f64;
        let trange = 0.4999e0f64;
        let ssmall = 1.0e-6f64;
        let tenth = 1.0e-1f64;
        let zero = 0.0f64;
        let nmax = 10;
        let nreac = self.nreac;
        let alpha = self.alpha;

        //--convert energies just loaded to velocity-like terms
        if self.mhigh >= self.mlow {
            for n in self.mlow..=self.mhigh {
                e[n] = sqrt(alpha * e[n]);
            }
        }
        let mut k = self.nlow;
        let mut klast = self.nlow - 1;
        let mut et = 0.0f64;
        let mut em = 0.0f64;
        let mut is: usize;
        let km = 0i32;

        //--set up first node in stack
        let mut tt1old = zero;
        let mut lbl;
        if self.klow == 1 {
            lbl = 110;
        } else {
            klast = k;
            self.bn_ks[2] = k;
            et = powi(e[k], 2) / alpha;
            self.bn_es[2] = sigfig(et, 7, 0);
            let xt = sqrt(alpha * self.bn_es[2]);
            for i in 1..=nreac {
                let mut dn = s[(i, k + 1)] - s[(i, k)];
                if dn.abs() < s[(i, k)].abs() / 1000.0 {
                    dn = 0.0;
                }
                self.bn_dl[i] = 1.0;
                if dn < zero {
                    self.bn_dl[i] = -1.0;
                }
            }
            let mut sn2 = RV::new(NTT);
            self.bsigma(c, k, xt, &mut sn2, e, s);
            for i in 1..=nreac {
                self.bn_ss[(i, 2)] = sn2[i];
            }
            self.bn_js[2] = 0;
            for i in 1..=nreac {
                if (self.bn_ss[(i, 2)] - s[(i, k)]).abs() > self.errthn * s[(i, k)] {
                    self.bn_js[2] = 1;
                }
            }
            k += 1;
            lbl = 120;
        }
        is = 2;
        // worker threads compute bsigma values ahead (exact; see BsCache)
        if !(lbl == 110 && self.bn_es[1] > self.thnmax) {
            self.bpf_init(e, s, k, klast);
        }
        let mut pidx = 0usize;
        loop {
            match lbl {
                110 => {
                    if self.bn_es[1] > self.thnmax {
                        lbl = 190;
                        continue;
                    }
                    lbl = 120;
                }
                120 => {
                    //--locate next node for stack
                    et = powi(e[k], 2) / alpha;
                    let mut go130 = false;
                    if et >= self.thnmax || k >= self.nhigh {
                        go130 = true;
                    } else {
                        let test = sigfig(et, 7, 0);
                        if (self.bn_es[2] - test).abs() < small * test {
                            // 126
                            k += 1;
                            continue;
                        }
                        if k > klast + nmax || et > step * self.bn_es[2] {
                            go130 = true;
                        } else {
                            let test = sigfig(et, 3, 0);
                            if (et - test).abs() < small * test || (et - therm).abs() < small * therm {
                                go130 = true;
                            } else {
                                for i in 1..=nreac {
                                    let mut dn = s[(i, k + 1)] - s[(i, k)];
                                    if dn.abs() < s[(i, k)].abs() / 1000.0 {
                                        dn = 0.0;
                                    }
                                    if dn >= zero {
                                        dn = 1.0;
                                    } else {
                                        dn = -1.0;
                                    }
                                    if dn != self.bn_dl[i] {
                                        go130 = true;
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    if !go130 {
                        // 126
                        k += 1;
                        continue;
                    }
                    // 130
                    self.bn_es[1] = sigfig(et, 7, 0);
                    let xt = sqrt(alpha * self.bn_es[1]);
                    let mut sn1 = RV::new(NTT);
                    pidx += 1;
                    self.bpf_ensure(pidx, e, s);
                    self.bsigma_c(c, k, xt, &mut sn1, e, s);
                    for i in 1..=nreac {
                        self.bn_ss[(i, 1)] = sn1[i];
                    }
                    self.bn_ks[1] = k;
                    self.bn_js[1] = 0;
                    for i in 1..=nreac {
                        if (self.bn_ss[(i, 1)] - s[(i, k)]).abs() > self.errthn * s[(i, k)] {
                            self.bn_js[1] = 1;
                        }
                    }
                    //--add points between nodes if needed
                    is = 2;
                    lbl = 140;
                }
                140 => {
                    if is >= NSTACK {
                        lbl = 150;
                        continue;
                    }
                    if self.bn_ks[is - 1] == self.bn_ks[is] + 1 && self.bn_js[is - 1] == 0 && self.bn_js[is] == 0 {
                        lbl = 150;
                        continue;
                    }
                    // compare true function to linear approximation
                    em = half * (self.bn_es[is - 1] + self.bn_es[is]);
                    let mut ndig = 9;
                    if em > tenth && em < one {
                        ndig = 8;
                    }
                    if em > sigfig(self.bn_es[is], 7, 1) {
                        if em < sigfig(self.bn_es[is - 1], 7, -1) {
                            em = sigfig(em, 7, 0);
                        }
                    } else {
                        em = sigfig(em, ndig, 0);
                    }
                    if em < sigfig(self.bn_es[is], ndig, 1) || em > sigfig(self.bn_es[is - 1], ndig, -1) {
                        lbl = 150;
                        continue;
                    }
                    let mut errt = self.errthn;
                    if self.bn_es[is - 1] < trange {
                        errt = errt / 5.0;
                    }
                    let mut errm = self.errmax;
                    if self.bn_es[is - 1] < trange {
                        errm = errm / 5.0;
                    }
                    let xm = sqrt(alpha * em);
                    let mut sn = std::mem::take(&mut self.bn_sn);
                    self.bsigma_c(c, klast, xm, &mut sn, e, s);
                    self.bn_sn = sn;
                    let dx = self.bn_es[is - 1] - self.bn_es[is];
                    let f = (em - self.bn_es[is]) / dx;
                    let test = one - one / 100.0;
                    if f > test {
                        lbl = 170;
                        continue;
                    }
                    let mut stot = zero;
                    for i in 1..=nreac {
                        stot = stot + self.bn_sn[i];
                    }
                    let mut go170 = false;
                    for i in 1..=nreac {
                        if stot < errmin {
                            continue;
                        }
                        let sni = self.bn_sn[i];
                        if (sni / stot).abs() < ssmall {
                            continue;
                        }
                        let a = self.bn_ss[(i, is as i32 - 1)];
                        let bb = self.bn_ss[(i, is as i32)];
                        if a > rmax * bb || a < bb / rmax {
                            go170 = true;
                            break;
                        }
                        let si = f * a + (1.0 - f) * bb;
                        let dy = (sni - si).abs();
                        if dy <= errt * sni.abs() + errmin {
                            continue;
                        }
                        if dy > errm * sni.abs() + errmin || dy * dx / 2.0 > self.errint * em {
                            go170 = true;
                            break;
                        }
                    }
                    if go170 {
                        lbl = 170;
                        continue;
                    }
                    // don't allow a big increase in the step size
                    let est = estp * (self.bn_es[is] - self.bn_tt[1]);
                    if *j > 3 && dx > est {
                        lbl = 170;
                        continue;
                    }
                    lbl = 150;
                }
                150 => {
                    // converged. store top point of stack.
                    self.bn_tt[1] = self.bn_es[is];
                    for i in 1..=nreac {
                        if self.bn_tt[1] > self.emtr[i] && tt1old >= self.emtr[i] {
                            self.bn_tt[1 + i] = self.bn_ss[(i, is as i32)];
                        } else {
                            self.bn_tt[1 + i] = 0.0;
                        }
                    }
                    tt1old = self.bn_tt[1];
                    *j += 1;
                    c.loada(*j, &self.bn_tt.0, self.ntx, inew, bufn, NBUF);
                    is -= 1;
                    if is > 1 {
                        lbl = 140;
                        continue;
                    }
                    // stack exhausted.  get new node.
                    self.bn_es[2] = self.bn_es[1];
                    self.bn_ks[2] = self.bn_ks[1];
                    self.bn_js[2] = self.bn_js[1];
                    for i in 1..=nreac {
                        self.bn_ss[(i, 2)] = self.bn_ss[(i, 1)];
                        let mut dn = s[(i, k + 1)] - s[(i, k)];
                        if dn.abs() < s[(i, k)].abs() / 1000.0 {
                            dn = 0.0;
                        }
                        if dn >= zero {
                            self.bn_dl[i] = 1.0;
                        } else {
                            self.bn_dl[i] = -1.0;
                        }
                    }
                    if k >= self.nhigh {
                        // 180: finished with this page.
                        self.bn_ks[2] = self.nlow - 1;
                        if self.nhigh < self.khigh {
                            self.bpf_done();
                            return;
                        }
                        lbl = 190;
                        continue;
                    }
                    if self.bn_es[1] > self.thnmax {
                        lbl = 190;
                        continue;
                    }
                    klast = k;
                    k += 1;
                    lbl = 120;
                }
                170 => {
                    // not converged.  add midpoint to stack and continue.
                    is += 1;
                    self.bn_es[is] = self.bn_es[is - 1];
                    self.bn_es[is - 1] = em;
                    self.bn_ks[is] = self.bn_ks[is - 1];
                    self.bn_ks[is - 1] = km;
                    self.bn_js[is] = self.bn_js[is - 1];
                    self.bn_js[is - 1] = 0;
                    if km == 0 {
                        self.bn_js[is - 1] = 1;
                    }
                    for i in 1..=nreac {
                        self.bn_ss[(i, is as i32)] = self.bn_ss[(i, is as i32 - 1)];
                        self.bn_ss[(i, is as i32 - 1)] = self.bn_sn[i];
                        if km != 0 && (self.bn_sn[i] - s[(i, km)]).abs() > self.errthn * s[(i, km)] + errmin {
                            self.bn_js[is - 1] = 1;
                        }
                    }
                    lbl = 140;
                }
                190 => {
                    // copy energies above thnmax to output
                    // flag last point with minus.
                    loop {
                        *j += 1;
                        et = powi(e[k], 2) / alpha;
                        self.bn_tt[1] = et;
                        for i in 1..=nreac {
                            self.bn_tt[1 + i] = s[(i, k)];
                        }
                        if k == self.nhigh && self.khigh == self.nhigh {
                            *j = -*j;
                        }
                        c.loada(*j, &self.bn_tt.0, self.ntx, inew, bufn, NBUF);
                        k += 1;
                        if k > self.nhigh {
                            break;
                        }
                    }
                    self.bn_ks[2] = self.nlow - 1;
                    self.bpf_done();
                    return;
                }
                _ => unreachable!(),
            }
        }
    }

    // ================================================================ bsigma
    fn bsigma(&mut self, c: &mut Ctx, kn: i32, en: f64, sn: &mut RV, e: &RV, s: &RM) {
        let aamin = self.bsigma_core(kn, en, sn, e, s);
        bsigma_mess(c, aamin);
    }

    /// bsigma without its message; returns aamin (the message is printed
    /// when aamin < atop).
    fn bsigma_core(&mut self, kn: i32, en: f64, sn: &mut RV, e: &RV, s: &RM) -> f64 {
        let one = 1.0e0f64;
        let two = 2.0e0f64;
        let atop = 4.0e0f64;
        let sigmin = 1.0e-15f64;
        let fzero = [0.0, 0.5e0, 1.0e0 / (2.0e0 * sqrt(PI)), 0.25e0, 1.0e0 / (2.0e0 * sqrt(PI)), 0.375e0];
        let nreac = self.nreac;
        let mut aamin = atop;
        self.aa = 0.0;
        let mut sbt = [0.0f64; NTT + 1];

        //--initialize integral
        let khm1 = self.khigh - 1;

        //--search for panel containing en
        let mut k = kn;
        if e[k + 1] >= e[k] {
            while en >= e[k + 1] {
                k += 1;
            }
        }

        //--compute desired cross sections
        self.y = en;
        self.oy = -one / self.y;
        self.yy = self.y * self.y;
        self.f[1..6].copy_from_slice(&fzero[1..6]);
        self.alast = 0.0;

        //--loop over intervals below current point
        let mut l = k + 1;
        let mut go170 = false;
        for _ll in self.klow..=k {
            l -= 1;
            self.x = e[l];
            let xp = e[l + 1];
            self.xx = xp * xp;
            if self.xx <= self.x * self.x {
                continue;
            }
            self.aa = self.y - self.x;
            self.hunky();
            let denom = one / (self.xx - self.x * self.x);
            for i in 1..=nreac {
                let slope = (s[(i, l + 1)] - s[(i, l)]) * denom;
                sbt[i as usize] = sbt[i as usize] + s[(i, l + 1)] * self.s1 + slope * self.s2;
            }
            if self.aa > atop {
                go170 = true;
                break;
            }
        }
        if !go170 {
            if self.aa < aamin && self.nlow != self.klow {
                aamin = self.aa;
            }
            // continue cross sections as 1/v to x=0.
            self.x = 0.0;
            self.xx = 0.0;
            self.aa = self.y;
            self.hunky();
            let klow = self.klow;
            for i in 1..=nreac {
                sbt[i as usize] =
                    sbt[i as usize] - s[(i, klow)] * e[klow] * (powi(self.oy, 2) * self.h[2] + self.oy * self.h[1]);
            }
        }

        //--reinitialize f functions to zero distance values
        // 170
        self.f[1..6].copy_from_slice(&fzero[1..6]);
        self.alast = 0.0;
        // reverse sign of 1/y to use same s1 and s2 functions
        self.oy = -self.oy;
        let mut go210 = false;
        // no intervals above the last point
        if k != self.khigh {
            //--loop over intervals above current one
            let mut go200 = true;
            for l in k..=khm1 {
                self.x = e[l + 1];
                let xm = e[l];
                self.xx = xm * xm;
                if self.xx >= self.x * self.x {
                    continue;
                }
                self.aa = self.x - self.y;
                self.hunky();
                let denom = one / (self.x * self.x - self.xx);
                for i in 1..=nreac {
                    let slope = (s[(i, l + 1)] - s[(i, l)]) * denom;
                    sbt[i as usize] = sbt[i as usize] + s[(i, l)] * self.s1 + slope * self.s2;
                }
                if self.aa > atop {
                    go200 = false;
                    go210 = true;
                    break;
                }
            }
            if go200 && self.aa < aamin && self.nhigh != self.khigh {
                aamin = self.aa;
            }
        }
        if !go210 {
            // 200: continue cross sections as constant to x=infinity
            let factor = (self.f[3] * self.oy + two * self.f[2]) * self.oy + self.f[1];
            let khigh = self.khigh;
            for i in 1..=nreac {
                sbt[i as usize] = sbt[i as usize] + s[(i, khigh)] * factor;
            }
        }

        //--add low energy term, if applicable
        // 210
        self.y = -self.y;
        self.aa = -self.y;
        'l230: {
            if self.aa > atop {
                break 'l230;
            }
            // extend cross section to x=0.0 as 1/v.
            self.funky();
            self.alast = self.aa;
            self.oy = -self.oy;
            let klow = self.klow;
            self.aa = e[klow] - self.y;
            self.xx = 0.0;
            self.hunky();
            for i in 1..=nreac {
                sbt[i as usize] =
                    sbt[i as usize] - s[(i, klow)] * e[klow] * (powi(self.oy, 2) * self.h[2] + self.oy * self.h[1]);
            }
            if self.aa > atop {
                break 'l230;
            }
            // loop over other intervals
            for l in klow..=khm1 {
                self.x = e[l + 1];
                self.aa = self.x - self.y;
                self.xx = powi(e[l], 2);
                if self.xx == self.x * self.x {
                    continue;
                }
                self.hunky();
                let denom = one / (self.x * self.x - self.xx);
                for i in 1..=nreac {
                    let slope = (s[(i, l + 1)] - s[(i, l)]) * denom;
                    sbt[i as usize] = sbt[i as usize] - s[(i, l)] * self.s1 - slope * self.s2;
                }
                if self.aa > atop {
                    break 'l230;
                }
            }
            // extend cross section to x=infinity as constant
            let factor = (self.f[3] * self.oy + two * self.f[2]) * self.oy + self.f[1];
            let khigh = self.khigh;
            for i in 1..=nreac {
                sbt[i as usize] = sbt[i as usize] - s[(i, khigh)] * factor;
            }
        }

        //--return broadened values
        // 230
        for i in 1..=nreac {
            if sbt[i as usize] < sigmin {
                sbt[i as usize] = 0.0;
            }
            sn[i] = sbt[i as usize];
        }

        aamin
    }

    /// bsigma, taking the value from the prefetch cache when it is there.
    fn bsigma_c(&mut self, c: &mut Ctx, kn: i32, en: f64, sn: &mut RV, e: &RV, s: &RM) {
        if let Some(pf) = self.bpf.0.as_mut() {
            if let Some(v) = pf.map.get(&(kn, en.to_bits())) {
                let nreac = self.nreac as usize;
                for i in 1..=nreac {
                    sn[i as i32] = v[i - 1];
                }
                let aamin = v[nreac];
                pf.hits += 1;
                bsigma_mess(c, aamin);
                return;
            }
            pf.misses += 1;
        }
        self.bsigma(c, kn, en, sn, e, s);
    }

    /// The stack nodes that broadn's loop will use, found by the same walk
    /// (label 120 and the end of label 150), which does not depend on any
    /// bsigma value. Node 0 is the lower node already on the stack.
    fn bpf_nodes(&self, e: &RV, s: &RM, k0: i32, klast0: i32) -> (Vec<(i32, f64)>, Vec<i32>) {
        let therm = 0.0253e0f64;
        let small = 1.0e-9f64;
        let step = 2.01e0f64;
        let zero = 0.0f64;
        let nmax = 10;
        let nreac = self.nreac;
        let alpha = self.alpha;
        let mut nodes = vec![(self.bn_ks[2], self.bn_es[2])];
        let mut klasts = vec![0];
        let mut dl = self.bn_dl.clone();
        let mut es2 = self.bn_es[2];
        let mut k = k0;
        let mut klast = klast0;
        // highest k whose e(k), s(:,k+1) exist (the walk never needs more
        // than broadn reads; this only keeps it from failing where broadn
        // would not go)
        let kmax = (e.len()).min((s.v.len() / s.n1.max(1)) as i32 - 1);
        loop {
            // 120
            let et;
            loop {
                if k < 1 || k > kmax {
                    return (Vec::new(), Vec::new());
                }
                let t = powi(e[k], 2) / alpha;
                let mut go130 = false;
                if t >= self.thnmax || k >= self.nhigh {
                    go130 = true;
                } else {
                    let test = sigfig(t, 7, 0);
                    if (es2 - test).abs() < small * test {
                        k += 1;
                        continue;
                    }
                    if k > klast + nmax || t > step * es2 {
                        go130 = true;
                    } else {
                        let test = sigfig(t, 3, 0);
                        if (t - test).abs() < small * test || (t - therm).abs() < small * therm {
                            go130 = true;
                        } else {
                            for i in 1..=nreac {
                                let mut dn = s[(i, k + 1)] - s[(i, k)];
                                if dn.abs() < s[(i, k)].abs() / 1000.0 {
                                    dn = 0.0;
                                }
                                if dn >= zero {
                                    dn = 1.0;
                                } else {
                                    dn = -1.0;
                                }
                                if dn != dl[i] {
                                    go130 = true;
                                    break;
                                }
                            }
                        }
                    }
                }
                if go130 {
                    et = t;
                    break;
                }
                k += 1;
            }
            // 130
            let es1 = sigfig(et, 7, 0);
            nodes.push((k, es1));
            klasts.push(klast);
            // 150, stack exhausted
            es2 = es1;
            for i in 1..=nreac {
                let mut dn = s[(i, k + 1)] - s[(i, k)];
                if dn.abs() < s[(i, k)].abs() / 1000.0 {
                    dn = 0.0;
                }
                dl[i] = if dn >= zero { 1.0 } else { -1.0 };
            }
            if k >= self.nhigh || es1 > self.thnmax {
                break;
            }
            klast = k;
            k += 1;
        }
        (nodes, klasts)
    }

    /// Set up the prefetch for the broadn loop that starts with node k,
    /// klast (the stack already holds the lower node). Panel p (from 1) is
    /// the interval between nodes p-1 and p.
    fn bpf_init(&mut self, e: &RV, s: &RM, k: i32, klast: i32) {
        self.bpf = BpfSlot(None);
        let nthr = crate::util::threads();
        if nthr < 2 {
            return;
        }
        let (nodes, klasts) = self.bpf_nodes(e, s, k, klast);
        if nodes.len() < 3 {
            return;
        }
        let nreac = self.nreac as usize;
        let alpha = self.alpha;
        let errthn = self.errthn;
        // node 0 is on the stack already
        let node0: (Vec<f64>, i32) = ((1..=nreac as i32).map(|i| self.bn_ss[(i, 2)]).collect(), self.bn_js[2]);
        let npan = nodes.len();
        let nodes = std::sync::Arc::new(nodes);
        let klasts = std::sync::Arc::new(klasts);
        let es = std::sync::Arc::new(e.clone());
        let ss = std::sync::Arc::new(s.clone());
        // bsigma at each node (values, aamin) and js, computed once
        type NodeVal = (Vec<f64>, f64, i32);
        let nv: std::sync::Arc<Vec<std::sync::OnceLock<NodeVal>>> =
            std::sync::Arc::new((0..npan).map(|_| std::sync::OnceLock::new()).collect());
        let mut workers = Vec::with_capacity(nthr);
        for _ in 0..nthr {
            workers.push(self.clone());
        }
        let ra = crate::prefetch::RunAhead::start(workers, 1, npan, 256 * nthr, move |w: &mut Broadr, q| {
            let (e, s) = (&*es, &*ss);
            let node = |w: &mut Broadr, idx: usize| -> NodeVal {
                let (k, ek) = nodes[idx];
                let xt = sqrt(alpha * ek);
                let mut sn = RV::new(NTT);
                let aamin = w.bsigma_core(k, xt, &mut sn, e, s);
                let v: Vec<f64> = (1..=nreac).map(|i| sn[i as i32]).collect();
                let mut js = 0;
                for i in 1..=nreac {
                    if (v[i - 1] - s[(i as i32, k)]).abs() > errthn * s[(i as i32, k)] {
                        js = 1;
                    }
                }
                (v, aamin, js)
            };
            // upper node first: the lower one is the upper node of the panel
            // before, which another worker is computing at the same time
            let up = nv[q].get_or_init(|| node(w, q)).clone();
            let lo = if q == 1 { (node0.0.clone(), 0.0, node0.1) } else { nv[q - 1].get_or_init(|| node(w, q - 1)).clone() };
            let mut out = Vec::new();
            let (k, ek) = nodes[q];
            let mut uv = up.0.clone();
            uv.push(up.1);
            out.push((k, sqrt(alpha * ek).to_bits(), uv));
            w.bpf_panel(nodes[q - 1], &(lo.0, lo.2), nodes[q], &(up.0, up.2), klasts[q], e, s, &mut out);
            out
        });
        self.bpf = BpfSlot(Some(Box::new(BsCache { map: Default::default(), cur: None, ra, nthr, hits: 0, misses: 0, wait: 0.0 })));
    }

    /// Make the values of panel p (computed by the workers) the current ones.
    fn bpf_ensure(&mut self, p: usize, _e: &RV, _s: &RM) {
        if let Some(pf) = self.bpf.0.as_mut() {
            if pf.cur == Some(p) {
                return;
            }
            let t0 = std::time::Instant::now();
            pf.map.clear();
            if let Some(v) = pf.ra.take(p) {
                for (kn, bits, vals) in v {
                    pf.map.insert((kn, bits), vals);
                }
            }
            pf.cur = Some(p);
            pf.wait += t0.elapsed().as_secs_f64();
        }
    }

    /// Refine one panel the way broadn does (labels 140, 150, 170) and keep
    /// every bsigma value computed. The last stored point before the panel
    /// (bn_tt(1)) is not known here, so the step-size test is not applied
    /// until the first point of the panel is stored, and j > 3 is assumed.
    #[allow(clippy::too_many_arguments)]
    fn bpf_panel(
        &mut self,
        lo: (i32, f64),
        lov: &(Vec<f64>, i32),
        up: (i32, f64),
        upv: &(Vec<f64>, i32),
        klast: i32,
        e: &RV,
        s: &RM,
        out: &mut Vec<(i32, u64, Vec<f64>)>,
    ) {
        let one = 1.0e0f64;
        let estp = 4.1e0f64;
        let half = 0.5e0f64;
        let rmax = 3.0e0f64;
        let errmin = 1.0e-15f64;
        let trange = 0.4999e0f64;
        let ssmall = 1.0e-6f64;
        let tenth = 1.0e-1f64;
        let zero = 0.0f64;
        let nreac = self.nreac as usize;
        let alpha = self.alpha;
        let mut es = [0.0f64; NSTACK + 1];
        let mut ks = [0i32; NSTACK + 1];
        let mut js = [0i32; NSTACK + 1];
        let mut ss = vec![0.0f64; (NSTACK + 1) * nreac];
        es[2] = lo.1;
        ks[2] = lo.0;
        js[2] = lov.1;
        ss[2 * nreac..3 * nreac].copy_from_slice(&lov.0);
        es[1] = up.1;
        ks[1] = up.0;
        js[1] = upv.1;
        ss[nreac..2 * nreac].copy_from_slice(&upv.0);
        let mut sn = RV::new(NTT);
        let mut tt1: Option<f64> = None;
        let mut is = 2usize;
        loop {
            // 140
            let mut go150 = true;
            let mut em = 0.0;
            if is < NSTACK && !(ks[is - 1] == ks[is] + 1 && js[is - 1] == 0 && js[is] == 0) {
                em = half * (es[is - 1] + es[is]);
                let mut ndig = 9;
                if em > tenth && em < one {
                    ndig = 8;
                }
                if em > sigfig(es[is], 7, 1) {
                    if em < sigfig(es[is - 1], 7, -1) {
                        em = sigfig(em, 7, 0);
                    }
                } else {
                    em = sigfig(em, ndig, 0);
                }
                if !(em < sigfig(es[is], ndig, 1) || em > sigfig(es[is - 1], ndig, -1)) {
                    let mut errt = self.errthn;
                    if es[is - 1] < trange {
                        errt = errt / 5.0;
                    }
                    let mut errm = self.errmax;
                    if es[is - 1] < trange {
                        errm = errm / 5.0;
                    }
                    let xm = sqrt(alpha * em);
                    let aamin = self.bsigma_core(klast, xm, &mut sn, e, s);
                    let mut v: Vec<f64> = (1..=nreac).map(|i| sn[i as i32]).collect();
                    v.push(aamin);
                    out.push((klast, xm.to_bits(), v));
                    let dx = es[is - 1] - es[is];
                    let f = (em - es[is]) / dx;
                    let test = one - one / 100.0;
                    let mut go170 = f > test;
                    if !go170 {
                        let mut stot = zero;
                        for i in 1..=nreac {
                            stot = stot + sn[i as i32];
                        }
                        for i in 1..=nreac {
                            if stot < errmin {
                                continue;
                            }
                            let sni = sn[i as i32];
                            if (sni / stot).abs() < ssmall {
                                continue;
                            }
                            let a = ss[(is - 1) * nreac + i - 1];
                            let bb = ss[is * nreac + i - 1];
                            if a > rmax * bb || a < bb / rmax {
                                go170 = true;
                                break;
                            }
                            let si = f * a + (1.0 - f) * bb;
                            let dy = (sni - si).abs();
                            if dy <= errt * sni.abs() + errmin {
                                continue;
                            }
                            if dy > errm * sni.abs() + errmin || dy * dx / 2.0 > self.errint * em {
                                go170 = true;
                                break;
                            }
                        }
                    }
                    if !go170 {
                        if let Some(t1) = tt1 {
                            let est = estp * (es[is] - t1);
                            if dx > est {
                                go170 = true;
                            }
                        }
                    }
                    if go170 {
                        go150 = false;
                    }
                }
            }
            if go150 {
                // 150
                tt1 = Some(es[is]);
                is -= 1;
                if is > 1 {
                    continue;
                }
                return;
            }
            // 170 (km = 0)
            is += 1;
            es[is] = es[is - 1];
            es[is - 1] = em;
            ks[is] = ks[is - 1];
            ks[is - 1] = 0;
            js[is] = js[is - 1];
            js[is - 1] = 1;
            for i in 0..nreac {
                ss[is * nreac + i] = ss[(is - 1) * nreac + i];
                ss[(is - 1) * nreac + i] = sn[(i + 1) as i32];
            }
        }
    }

    /// End of a broadn call: drop the cache.
    fn bpf_done(&mut self) {
        if let Some(pf) = self.bpf.0.take() {
            if crate::util::pf_stats() {
                eprintln!("broadr prefetch: {} hits, {} misses, {} threads, {:.2}s waiting", pf.hits, pf.misses, pf.nthr, pf.wait);
            }
        }
    }

    // ================================================================ thinb
    fn thinb(&mut self, c: &mut Ctx, j: &mut i32, inew: i32, bufn: &mut [f64], e: &RV, s: &RM) {
        let stpmax = 1.24e0f64;
        let eps = 1.0e-5f64;
        let zero = 0.0f64;
        let nreac = self.nreac;
        let nmin = self.nlow + 2;
        let mut ilow = self.nlow + 1;
        let tt = &mut self.th_tt;
        // `1.+eps`: the real(4) literal 1. is widened to real(8) first
        let one_eps = 1.0f64 + eps;

        //--keep first point
        if self.klow != 1 {
            *j += 1;
            tt[1] = e[self.nlow];
            for i in 1..=nreac {
                tt[1 + i] = s[(i, self.nlow)];
            }
            c.loada(*j, &tt.0, self.ntx, inew, bufn, NBUF);
        }

        //--loop over all other points
        // nm1 is read at label 180 even if the loop never set it (Fortran: undefined)
        let mut nm1 = 0i32;
        let mut n = nmin;
        let mut go180 = false;
        while n <= self.nhigh {
            // define index to last candidate for thinning
            nm1 = n - 1;
            let mut denom = e[n] - tt[1];
            if denom == zero {
                n += 1;
                continue;
            }
            denom = 1.0 / denom;
            let mut fail = false;
            'outer: for i in 1..=nreac {
                let slope = (s[(i, n)] - tt[1 + i]) * denom;
                // test allowable error at all intervening points
                for k in ilow..=nm1 {
                    if e[k] >= stpmax * tt[1] {
                        fail = true;
                        break 'outer;
                    }
                    let sp = tt[1 + i] + slope * (e[k] - tt[1]);
                    if (sp - s[(i, k)]).abs() > self.errthn * s[(i, k)] {
                        fail = true;
                        break 'outer;
                    }
                    if e[k] * one_eps > self.thnmax {
                        fail = true;
                        break 'outer;
                    }
                }
            }
            if fail {
                // 140: at least one point in interval fails test
                *j += 1;
                tt[1] = e[nm1];
                for i in 1..=nreac {
                    tt[1 + i] = s[(i, nm1)];
                }
                if tt[1] * one_eps >= self.thnmax {
                    go180 = true;
                    break;
                }
                c.loada(*j, &tt.0, self.ntx, inew, bufn, NBUF);
                ilow = n;
            }
            n += 1;
        }
        if !go180 {
            //--finished with this page.
            return;
        }
        //--copy energies above thnmax to scratch storage.
        loop {
            // 180
            if nm1 == self.nhigh && self.khigh == self.nhigh {
                *j = -*j;
            }
            c.loada(*j, &tt.0, self.ntx, inew, bufn, NBUF);
            nm1 += 1;
            if nm1 > self.nhigh {
                return;
            }
            // 160
            *j += 1;
            tt[1] = e[nm1];
            for i in 1..=nreac {
                tt[i + 1] = s[(i, nm1)];
            }
        }
    }

    // ================================================================ hunky
    fn hunky(&mut self) {
        let small = 1.0e-12f64;
        let toler = 1.0e-5f64;
        let two = 2.0e0f64;
        let four = 4.0e0f64;
        let six = 6.0e0f64;
        // initialize h functions to current f functions
        for k in 1..=5 {
            self.h[k] = self.f[k];
        }
        // compute new f functions
        self.funky();
        // subtract new f functions to obtain integral over interval
        for k in 1..=5 {
            self.h[k] = self.h[k] - self.f[k];
            if self.h[k].abs() <= small * self.f[k].abs() {
                self.h[k] = 0.0;
            }
            if self.h[k].abs() <= toler * self.f[k].abs() && self.aa != self.alast {
                let n = k as i32 - 1;
                self.h[k] = hnabb(n, self.alast, self.aa);
            }
        }
        self.alast = self.aa;
        // compute contributions from two integral terms
        let (h, oy, yy, xx) = (&self.h, self.oy, self.yy, self.xx);
        self.s1 = (h[3] * oy + two * h[2]) * oy + h[1];
        self.s2 = ((h[5] + (six * yy - xx) * h[3]) * oy + (four * h[4] + (four * yy - two * xx) * h[2])) * oy
            + (yy - xx) * h[1];
    }

    // ================================================================ funky
    fn funky(&mut self) {
        let half = 0.5e0f64;
        let two = 2.0e0f64;
        let three = 3.0e0f64;
        let alim = 10.0e0f64;
        let resqpi = 1.0 / sqrt(PI);
        let a = self.aa;
        let asq = a * a;
        let mut expo = 0.0f64;
        if a < alim {
            expo = exp(-asq);
        }
        self.f[1] = 0.0;
        if a < alim {
            self.f[1] = half * erfc(a);
        }
        expo = expo * resqpi;
        self.f[2] = half * expo;
        expo = expo * a;
        self.f[3] = half * (self.f[1] + expo);
        expo = expo * a;
        self.f[4] = half * (two * self.f[2] + expo);
        expo = expo * a;
        self.f[5] = half * (three * self.f[3] + expo);
    }
}

/// bsigma's truncation warning.
fn bsigma_mess(c: &mut Ctx, aamin: f64) {
    let atop = 4.0e0f64;
    if aamin < atop {
        let strng = format!("broadening truncated at a={}", ffmt::f(aamin, 3, 1));
        c.mess("bsigma", &strng, " ");
    }
}

// ==================================================================== hnabb
/// For b.ge.a, h_n(a,b) by a direct Taylor series of the defining integral.
pub fn hnabb(n: i32, aa: f64, bb: f64) -> f64 {
    let aerr = 1.0e30f64;
    let rerr = 1.0e-8f64;
    let pow2 = [0.0, 1.4142135623731e0, 2.0e0, 2.8284271247462e0, 4.0e0, 5.6568542494924e0];
    let explim = 100.0e0f64;
    let zero = 0.0f64;
    let resqpi = 1.0 / sqrt(PI);
    let mut cm = [0.0f64; 51];
    let mut cmstar = [0.0f64; 51];

    //--transform limits so that b ge a
    let mut sign = 1.0f64;
    let a;
    let b;
    if bb < aa {
        a = bb.abs();
        b = aa.abs();
        sign = -sign;
    } else {
        a = aa.abs();
        b = bb.abs();
    }
    if bb < zero && n % 2 != 0 {
        sign = -sign;
    }
    let h = (b - a) * pow2[1];
    let x = pow2[1] * a;
    let xx = x * x;
    let asq = a * a;
    let mut con = 0.0f64;
    if asq < explim {
        con = exp(-asq) * resqpi / pow2[(n + 1) as usize];
    }
    let mut mflag = 0;
    let mut k = n;
    let mut kd = 0i32;
    cm[1] = 1.0;
    let mut xk = 1.0f64;
    if k != 0 {
        xk = powi(x, k);
    }
    let mut s = h * xk;
    let mut fact = h;
    for m in 2..=50i32 {
        fact = fact * h / m as f64;
        let kstar = k;
        let kdstar = kd + 1;
        for j in 1..=kdstar as usize {
            cmstar[j] = cm[j];
        }
        k = n - m + 1;
        if k < 0 {
            let kk = k % 2;
            k = 0;
            if kk != 0 {
                k = 1;
            }
        }
        kd = (n + m - 1 - k) / 2;
        cm[(kd + 1) as usize] = -cmstar[kdstar as usize];
        let mut qmn = cm[(kd + 1) as usize];
        if kd != 0 {
            for j in 1..=kd {
                let jalpha = (2 * j + k - 1 - kstar) / 2;
                let mut beta = 0.0f64;
                if jalpha != 0 {
                    beta = cmstar[jalpha as usize];
                }
                cm[j as usize] = (2 * j + k - 1) as f64 * cmstar[(jalpha + 1) as usize] - beta;
            }
            for j in 1..=kd {
                qmn = qmn * xx + cm[(kd + 1 - j) as usize];
            }
        }
        xk = 1.0;
        if k != 0 {
            xk = powi(x, k);
        }
        let term = fact * xk * qmn;
        s = s + term;
        let xn1 = (n + 1) as f64;
        let xn2 = (h * x).abs();
        if (m as f64) < xn1.max(xn2) {
            continue;
        }
        let test = aerr + rerr * s.abs();
        if term.abs() > test {
            mflag = 0;
            continue;
        }
        if mflag == 1 {
            break;
        }
        mflag = 1;
    }
    con * s * sign
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
