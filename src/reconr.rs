//! Port of NJOY2016 `reconr.f90` (module `reconm`).
//!
//! The Fortran module variables are the fields of `Reconr`; each Fortran
//! subroutine is a method with the same name. Array indices are kept 1-based
//! (see `farr`) so that index arithmetic matches the Fortran text.

#![allow(clippy::too_many_arguments, clippy::needless_range_loop, unused_assignments)]

use crate::endf::Ctx;
use crate::farr::{IV, RM, RV};
use crate::ffmt;
use crate::fio::{Kind, Val};
use crate::fmath::{log10, nint, sqrt, pow};
use crate::samm::Samm;
use crate::util::physics::{AMASSN, AMU, EV, HBAR};
use crate::util::sigfig;

pub const NMTMAX: usize = 10;
pub const MAXUNR: i32 = 500;
pub const NODMAX: i32 = 800000;
pub const NBUFG: i32 = 2000;
pub const NBUFR: i32 = 2000;
pub const NBUF: i32 = 2000;
pub const NBUFL: i32 = 2000;

/// Extra words allocated past the Fortran array bounds. The Fortran checks
/// some bounds only after writing; this slack keeps those writes from
/// panicking here so that the same error message is produced.
const SLACK: usize = 4096;

#[derive(Default)]
pub struct Reconr {
    pub nin: i32,
    pub nout: i32,
    pub zain: f64,
    pub awin: f64,
    pub zai: f64,
    pub el: f64,
    pub eh: f64,
    pub err: f64,
    pub errmax: f64,
    pub errint: f64,
    pub za: f64,
    pub awr: f64,
    pub tempr: f64,
    pub q18: f64,
    pub tempi: f64,
    pub elis: f64,
    pub sta: f64,
    pub efmax: f64,
    pub lis: i32,
    pub lis0: i32,
    pub nfor: i32,
    pub lrel: i32,
    pub nver: i32,
    pub lfw: i32,
    pub mata: i32,
    pub itype: i32,
    pub lrp: i32,
    pub lfi: i32,
    pub lssf: i32,
    pub lrx: i32,
    pub mtr4: i32,
    pub mtr18: i32,
    pub mtr: [i32; NMTMAX + 1],
    pub mtrt: [i32; NMTMAX + 1],
    pub nmtr: i32,
    pub mt103: i32,
    pub mt104: i32,
    pub mt105: i32,
    pub mt106: i32,
    pub mt107: i32,
    pub mpmin: i32,
    pub mpmax: i32,
    pub mdmin: i32,
    pub mdmax: i32,
    pub mtmin: i32,
    pub mtmax: i32,
    pub m3min: i32,
    pub m3max: i32,
    pub m4min: i32,
    pub m4max: i32,
    pub nxc: i32,
    pub ngo: i32,
    pub mtr522: i32,
    pub ncards: i32,
    pub isot: [i32; 21],
    pub modet: [i32; 21],
    pub ibaset: [i32; 21],
    pub isect: i32,
    pub nsect: i32,
    pub abnt: [f64; 21],
    pub elt: [f64; 21],
    pub eht: [f64; 21],
    pub eresl: f64,
    pub eresr: f64,
    pub eresh: f64,
    pub eresu: f64,
    pub eresm: f64,
    pub nodes: i32,
    pub nunr: i32,
    pub maxres: i32,
    pub elst: [f64; 21],
    pub slst: RM,
    pub enxt: [f64; 21],
    pub snxt: RM,
    pub ilast: [i32; 21],
    pub mmtres: [i32; 11],
    pub mcards: [i32; 11],
    pub nmtres: i32,
    pub ncoef: i32,
    pub nsig: i32,
    pub nresp: i32,
    pub spin: f64,
    pub ascat: f64,
    pub card: Vec<[u8; 4]>,
    pub enode: RV,
    pub eunr: RV,
    pub res: RV,
    pub tr: Option<RM>,
    pub ti: Option<RM>,
    pub dict: RV,
    pub mfs: IV,
    pub mts: IV,
    pub ncs: IV,
    pub sunr: Option<RV>,
    pub samm: Samm,
    // saved locals of csunr1/csunr2/... are added where those are ported
    pub sv: Saved,
}

/// Fortran `save`d local variables and other persistent locals.
#[derive(Default)]
pub struct Saved {
    pub text66: Vec<u8>,
    /// lunion: `real(kr)::elim=.99e6_kr` (implicitly SAVEd)
    pub elim: Option<f64>,
    /// csmlbw2: `save enow`
    pub enow_mlbw2: f64,
    /// csrmat locals ri, si (left unchanged if the matrix is singular)
    pub csr_ri: [[f64; 4]; 4],
    pub csr_si: [[f64; 4]; 4],
    /// cshybr local eeff (used before set when there are no channels)
    pub hy_eeff: f64,
    pub iest_u1: i32,
    pub ifirst_u1: i32,
    pub spot_u1: f64,
    pub spot_u2: f64,
    /// csaa locals
    pub aa_delta: f64,
    pub aa_ki: i32,
    /// recout local z(17) (equivalenced with tz(17))
    pub recout_z: Vec<f64>,
}

fn txt(b: &[u8]) -> String {
    String::from_utf8_lossy(b).to_string()
}

/// Write 4-character items with format `(1x,17a4)` and format reversion.
fn lines_17a4(items: &[[u8; 4]]) -> String {
    let mut out = String::new();
    for (k, chunk) in items.chunks(17).enumerate() {
        if k > 0 {
            out.push('\n');
        }
        out.push(' ');
        for c in chunk {
            out.push_str(&txt(c));
        }
    }
    if items.is_empty() {
        out.push(' ');
    }
    out
}

/// Split a 66-character text into 17 `a4` items as `read(text,'(16a4,a2)')`.
fn text_to_a4(text: &[u8]) -> Vec<[u8; 4]> {
    let mut t = text.to_vec();
    while t.len() < 66 {
        t.push(b' ');
    }
    let mut v = Vec::with_capacity(17);
    for i in 0..16 {
        v.push([t[4 * i], t[4 * i + 1], t[4 * i + 2], t[4 * i + 3]]);
    }
    v.push([t[64], t[65], b' ', b' ']);
    v
}

fn a4_to_bytes(v: &[[u8; 4]]) -> Vec<u8> {
    let mut b = Vec::with_capacity(v.len() * 4);
    for c in v {
        b.extend_from_slice(c);
    }
    b
}

/// The first 4*n bytes of a real array viewed as character(4) items.
fn real_to_a4(r: &[f64], n: usize) -> Vec<[u8; 4]> {
    let mut bytes = Vec::new();
    for x in r {
        bytes.extend_from_slice(&x.to_le_bytes());
    }
    (0..n).map(|i| [bytes[4 * i], bytes[4 * i + 1], bytes[4 * i + 2], bytes[4 * i + 3]]).collect()
}

pub fn so_lines6e(vals: &[f64]) -> String {
    // '(3x,1p,6e12.5)' with format reversion
    let mut out = String::new();
    for (k, chunk) in vals.chunks(6).enumerate() {
        if k > 0 {
            out.push('\n');
        }
        out.push_str("   ");
        for v in chunk {
            out.push_str(&ffmt::e1p(*v, 12, 5));
        }
    }
    out
}

impl Reconr {
    pub fn new() -> Self {
        let mut r = Reconr::default();
        r.slst = RM::new(20, 4);
        r.snxt = RM::new(20, 4);
        r.sv.text66 = vec![b' '; 66];
        r.sv.recout_z = vec![0.0; 17];
        r.sv.aa_ki = 1;
        r
    }

    fn read_list(&mut self, c: &mut Ctx, kinds: &[Kind]) -> Vec<Option<Val>> {
        c.io.input.read(kinds)
    }

    // ================================================================ reconr
    pub fn reconr(&mut self, c: &mut Ctx) {
        let elow = 1.0e-5f64;
        let ehigh = 20.0e6f64;
        let elarge = 9.0e9f64;
        let (mut nb, mut nw) = (0i32, 0i32);
        let mut nrtot = 0i32;
        let mut n6 = 0i32;
        let mut nx = 0i32;
        let mut iold = 0i32;
        let mut inew = 0i32;
        let mut ngrid = 0i32;
        let mut intunr = 0i32;

        //--set samrml options
        self.samm.want_partial_derivs = false;
        self.samm.want_angular_dist = false;
        self.samm.want_samrml_bw = false;
        self.samm.want_samrml_rm = false;

        //--define local i/o units
        let mut nscr1 = -10;
        let nscr2 = 11;
        let nscr3 = 12;

        //--write heading and read user input
        let time = c.timer();
        c.io.so(&format!(
            "\n reconr...reconstruct pointwise cross sections in pendf format{}{}s",
            " ".repeat(7),
            ffmt::f(time, 8, 1)
        ));
        c.io.se(&format!("\n reconr...{}{}s", " ".repeat(59), ffmt::f(time, 8, 1)));
        let v = self.read_list(c, &[Kind::Int, Kind::Int]);
        // nendf,npend are uninitialised locals in Fortran; 0 here
        let nendf = v[0].as_ref().map(|x| x.int()).unwrap_or(0);
        let npend = v[1].as_ref().map(|x| x.int()).unwrap_or(0);
        c.openz(nendf, 0);
        c.openz(npend, 1);
        let v = self.read_list(c, &[Kind::Char(66)]);
        if let Some(Val::Char(s)) = &v[0] {
            self.sv.text66 = s.as_bytes().to_vec();
        }
        let tlabel = text_to_a4(&self.sv.text66);
        c.io.so(&format!(
            "\n unit for endf tape ................... {}\n unit for pendf tape .................. {}",
            ffmt::i(nendf, 10),
            ffmt::i(npend, 10)
        ));
        c.io.so(&format!("\n label for pendf tape\n {}\n{}", "-".repeat(66), lines_17a4(&tlabel)));

        //--read tape identification record from input tape.
        //--write new tape identification label on output tape.
        self.nin = nendf;
        self.nout = npend;
        c.repoz(self.nin);
        c.repoz(self.nout);
        let mut z = RV::new(17);
        c.tpidio(self.nin, 0, 0, &mut z.all(), &mut nb, &mut nw);
        let tz = real_to_a4(&z.0, nw.clamp(0, 17) as usize);
        c.io.so(&format!("\n tape label\n {}\n{}", "-".repeat(66), lines_17a4(&tz)));
        c.math = 1;
        c.mfh = 0;
        c.mth = 0;
        c.nsh = 0;
        nw = 17;
        let mut rlabel = RV::new(17);
        crate::farr::pack_text(&a4_to_bytes(&tlabel), &mut rlabel.0);
        c.tpidio(0, self.nout, 0, &mut rlabel.all(), &mut nb, &mut nw);
        if npend > 0 {
            nscr1 = 10;
        }
        self.mata = 1;

        //--loop through users materials
        while self.mata != 0 {
            self.enode = RV::new(NODMAX as usize + SLACK);

            //--read users input for next material.
            self.ruina(c);
            if self.mata != 0 {
                //--open files for this material
                c.openz(-nscr2, 1);
                c.openz(-nscr3, 1);
                c.openz(nscr1, 1);
                c.repoz(nscr1);
                c.repoz(nscr2);
                c.repoz(nscr3);
                let mut nscr4 = 13;
                let nwscr = c.npage + 50;
                let mut scr = RV::new(nwscr as usize + SLACK);

                //--position to required material
                c.findf(self.mata, 1, 451, self.nin);

                //--process file 1
                self.tempi = 0.0;
                c.contio(self.nin, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                self.za = c.c1h;
                self.awr = c.c2h;
                self.lrp = c.l1h;
                self.lfi = c.l2h;
                n6 = c.n2h;
                c.contio(self.nin, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                if c.n1h != 0 {
                    c.iverf = 4;
                    nx = n6;
                } else if c.n2h == 0 {
                    c.iverf = 5;
                } else {
                    c.iverf = 6;
                }
                c.skiprz(self.nin, -1);
                if c.iverf >= 5 {
                    c.contio(self.nin, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                    self.elis = c.c1h;
                    self.sta = c.c2h;
                    self.lis = c.l1h;
                    self.lis0 = c.l2h;
                    self.nfor = c.n2h;
                }
                let mut nsub = 10;
                self.zain = 1.0;
                self.awin = 1.0;
                if c.iverf >= 6 {
                    c.contio(self.nin, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                    nsub = c.n1h;
                    self.zain = (nsub / 10) as f64;
                    self.awin = c.c1h;
                    self.efmax = c.c2h;
                    self.lrel = c.l1h;
                    self.nver = c.n2h;
                }
                if 10 * nint(self.zain) != nsub && nsub != 3 {
                    c.error("reconr", "illegal nsub for reconr", " ");
                }
                c.hdatio(self.nin, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                if c.iverf >= 6 {
                    self.tempi = scr[1];
                }
                for i in 1..=17 {
                    z[i] = scr[6 + i];
                }
                let tz = real_to_a4(&z.0, 17);
                c.io.so(&format!(
                    "\n processing mat {} in endf-{} format\n {}\n{}\n",
                    ffmt::i(self.mata, 4),
                    ffmt::i(c.iverf, 1),
                    "-".repeat(66),
                    lines_17a4(&tz)
                ));
                while nb != 0 {
                    c.moreio(self.nin, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                }
                if c.iverf != 4 {
                    nx = c.n2h;
                }
                self.anlyzd(c, self.nin, nx);
                c.tofend(self.nin, 0, 0, &mut scr.all());

                //--scan file 2 to check whether sammy method is needed.
                self.lrx = 0;
                self.nmtres = 0;
                if self.lrp == 1 {
                    // isammy = 1
                    self.res = RV::new(self.maxres.max(0) as usize + SLACK);
                    let nin = self.nin;
                    let maxres = self.maxres;
                    self.samm.s2sammy(c, nin, &mut self.res, maxres, &mut self.mmtres, &mut self.nmtres);
                    if self.nmtres > 0 {
                        for i in 1..=self.nmtres as usize {
                            if self.mmtres[i] == 103 {
                                self.mmtres[i] = 600;
                            }
                            if self.mmtres[i] == 104 {
                                self.mmtres[i] = 650;
                            }
                            if self.mmtres[i] == 105 {
                                self.mmtres[i] = 700;
                            }
                            if self.mmtres[i] == 106 {
                                self.mmtres[i] = 750;
                            }
                            if self.mmtres[i] == 107 {
                                self.mmtres[i] = 800;
                            }
                        }
                    } else {
                        self.samm.desammy();
                    }
                    c.repoz(self.nin);
                    c.findf(self.mata, 2, 0, self.nin);
                    self.res = RV::default();
                }

                //--process file 2
                self.eunr = RV::new(MAXUNR as usize + SLACK);
                self.nunr = 0;
                if self.lrp != 1 {
                    c.io.so("\n mat has no resonance parameters");
                    nrtot = 0;
                    self.eresl = elow;
                    self.eresr = ehigh;
                    self.eresh = ehigh;
                    self.eresu = elarge;
                    self.eresm = elarge;
                    if self.lrp == 0 || self.lrp > 1 {
                        c.contio(self.nin, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                        c.contio(self.nin, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                        c.contio(self.nin, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                        c.contio(self.nin, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                        self.spin = c.c1h;
                        self.ascat = c.c2h;
                    } else if self.lrp < 0 {
                        self.spin = 0.0;
                        self.ascat = 0.0;
                    }
                } else {
                    nrtot = self.maxres;
                    self.res = RV::new(nrtot.max(0) as usize + SLACK);
                    self.ncoef = 1;
                    for i in 1..=10 {
                        self.mcards[i] = 0;
                    }
                    self.rdfil2(c, &mut nrtot, &mut intunr);
                    if self.nmtres == 0 {
                        self.nsig = 5;
                    } else {
                        self.nsig = 4;
                        let mut s = String::from("\n   samm resonance reactions:");
                        // '10i5' with reversion (nmtres <= 10)
                        for i in 1..=self.nmtres as usize {
                            s.push_str(&ffmt::i(self.mmtres[i], 5));
                        }
                        c.io.so(&s);
                        if self.nmtres > 2 {
                            for i in 3..=self.nmtres as usize {
                                if self.mmtres[i] != 18 {
                                    self.nsig += 1;
                                }
                            }
                        }
                    }
                    c.io.so(&format!("   samm max legendre order:{}", ffmt::i(self.ncoef - 1, 3)));
                    if self.ncoef > 1 {
                        c.io.so("   generating File 4 for resonance angular distributions");
                    }
                    if self.lrp == 1 {
                        c.io.so("\n mat has no unresolved resonance parameters");
                    }
                    if self.lrp == 3 {
                        self.genunr(c, self.nin, intunr);
                    }
                    if self.lssf > 0 {
                        self.eresh = self.eresr;
                    }
                }
                self.eunr = RV::default();

                //--process file 3 and file 13
                self.lunion(c, self.nin, nscr1, nscr2, &mut ngrid);
                self.enode = RV::default();

                //--compute resonance cross sections
                if self.lrp == 1 || self.lrp == 3 {
                    nrtot = self.resxs(c, nscr2, nscr3, nscr4, nrtot);
                    self.res = RV::default();
                }

                //--clean up after sammy calculation
                if self.nmtres > 0 {
                    self.samm.desammy();
                }

                //--do energy merge and add resonance cross sections
                if npend < 0 {
                    nscr4 = -nscr4;
                }
                self.emerge(c, nscr1, nscr2, ngrid, nscr3, nrtot, &mut iold, &mut inew, nscr4);

                //--write output tape.
                self.recout(c, iold, nscr4, nrtot);

                //--deallocate arrays for this material
                self.card.clear();
                self.mfs = IV::default();
                self.mts = IV::default();
                self.ncs = IV::default();
                self.tr = None;
                self.ti = None;
                self.sunr = None;
                // inquire(unit=16,exist=there) is always true in gfortran
                c.closz(16);
                let _ = scr;
            } else {
                self.enode = RV::default();
            }
        }

        //--reconstruction complete.
        c.atend(self.nout, 0);
        c.repoz(self.nout);
        c.closz(self.nout);
        c.closz(self.nin);
        let time = c.timer();
        c.io.so(&format!("{}{}s\n {}*******", " ".repeat(69), ffmt::f(time, 8, 1), "**********".repeat(7)));
    }

    // ================================================================= ruina
    fn ruina(&mut self, c: &mut Ctx) {
        let zero = 0.0f64;
        self.ncards = 0;
        let mut ngrid = 0;
        let v = self.read_list(c, &[Kind::Int, Kind::Int, Kind::Int]);
        if let Some(x) = &v[0] {
            self.mata = x.int();
        }
        if let Some(x) = &v[1] {
            self.ncards = x.int();
        }
        if let Some(x) = &v[2] {
            ngrid = x.int();
        }
        if self.mata == 0 {
            return;
        }
        self.tempr = 0.0;
        self.errmax = -1.0;
        self.errint = -1.0;
        let v = self.read_list(c, &[Kind::Real, Kind::Real, Kind::Real, Kind::Real]);
        if let Some(x) = &v[0] {
            self.err = x.real();
        }
        if let Some(x) = &v[1] {
            self.tempr = x.real();
        }
        if let Some(x) = &v[2] {
            self.errmax = x.real();
        }
        if let Some(x) = &v[3] {
            self.errint = x.real();
        }
        if self.errmax <= zero {
            self.errmax = 10.0 * self.err;
        }
        if self.errmax < self.err {
            self.errmax = self.err;
        }
        if self.errint <= zero {
            self.errint = self.err / 20000.0;
        }
        if self.ncards > 0 {
            self.card = Vec::with_capacity(17 * self.ncards as usize);
            for _ in 1..=self.ncards {
                let mut text = vec![b' '; 66];
                let v = self.read_list(c, &[Kind::Char(66)]);
                if let Some(Val::Char(s)) = &v[0] {
                    text = s.as_bytes().to_vec();
                }
                self.card.extend(text_to_a4(&text));
            }
        } else {
            self.ncards = 1;
            self.card = text_to_a4(&[b' '; 66]);
        }
        if ngrid > 0 {
            let kinds = vec![Kind::Real; ngrid as usize];
            let v = self.read_list(c, &kinds);
            for i in 1..=ngrid {
                if let Some(x) = &v[(i - 1) as usize] {
                    self.enode[i] = x.real();
                }
            }
            for i in 1..=ngrid {
                self.enode[i] = sigfig(self.enode[i], 7, 0);
            }
        }

        //--write out user input
        c.io.so(&format!(
            "\n material to be processed ............. {}\n reconstruction tolerance ............. {}\n reconstruction temperature ........... {}k\n resonance-integral-check tolerance ... {}\n max resonance-integral error ......... {}",
            ffmt::i(self.mata, 10),
            ffmt::f(self.err, 10, 3),
            ffmt::f(self.tempr, 10, 2),
            ffmt::f(self.errmax, 10, 3),
            ffmt::e1p(self.errint, 10, 3)
        ));
        if self.ncards > 0 {
            c.io.so(&format!("\n descriptive cards for pendf tape\n {}", "-".repeat(66)));
            c.io.so(&lines_17a4(&self.card));
        }
        c.io.so(" ");
        if ngrid > 0 {
            c.io.so(&format!(" no. users energy grid points ......... {}", ffmt::i(ngrid, 10)));
            let vals: Vec<f64> = (1..=ngrid).map(|i| self.enode[i]).collect();
            c.io.so(&so_lines6e(&vals));
        }
        self.nodes = ngrid;

        //--construct w table if needed for psi-chi.
        if self.tempr > zero {
            let mut tr = RM::new(62, 62);
            let mut ti = RM::new(62, 62);
            wtab(&mut tr, &mut ti);
            self.tr = Some(tr);
            self.ti = Some(ti);
        }
    }

    // ================================================================ anlyzd
    fn anlyzd(&mut self, c: &mut Ctx, nin: i32, nx: i32) {
        let mut nw = 6 * nx;
        self.dict = RV::new(nw.max(0) as usize + 6);
        nw = nx;
        let mut nb = 0;
        c.dictio(nin, 0, 0, &mut self.dict.all(), &mut nb, &mut nw);

        let mut igam = 0;
        let mut nxn = 0;
        self.nmtr = 0;
        self.mtr4 = 0;
        self.mtr18 = 0;
        self.mt103 = 0;
        self.mt104 = 0;
        self.mt105 = 0;
        self.mt106 = 0;
        self.mt107 = 0;
        if c.iverf >= 6 {
            self.mpmin = 600;
            self.mpmax = 649;
            self.mdmin = 650;
            self.mdmax = 699;
            self.mtmin = 700;
            self.mtmax = 749;
            self.m3min = 750;
            self.m3max = 799;
            self.m4min = 800;
            self.m4max = 849;
        } else {
            self.mpmin = 700;
            self.mpmax = 718;
            self.mdmin = 720;
            self.mdmax = 738;
            self.mtmin = 740;
            self.mtmax = 758;
            self.m3min = 760;
            self.m3max = 768;
            self.m4min = 780;
            self.m4max = 798;
        }
        self.mtr522 = 0;
        let mut i = 1;
        while i <= nw {
            if self.nmtr >= NMTMAX as i32 {
                c.error("anlyzd", "too many redundant reactions.", " ");
            }
            let mfi = nint(self.dict[i + 2]);
            if mfi == 2 {
                let mti = nint(self.dict[i + 3]);
                if mti == 151 {
                    self.maxres = 12 * nint(self.dict[i + 4]);
                }
            } else if mfi == 3 {
                nxn += 1;
                if self.nmtr == 0 {
                    self.mtr[1] = 1;
                }
                if self.nmtr == 0 {
                    self.nmtr = 1;
                }
                let mti = nint(self.dict[i + 3]);
                if mti == 19 {
                    self.nmtr += 1;
                    self.mtr[self.nmtr as usize] = 18;
                    self.mtr18 = 1;
                }
                if (51..=91).contains(&mti) && self.mtr4 == 0 {
                    self.nmtr += 1;
                    self.mtr[self.nmtr as usize] = 4;
                    self.mtr4 = 1;
                }
                if mti >= self.mpmin && mti <= self.mpmax && self.mt103 == 0 {
                    self.nmtr += 1;
                    self.mtr[self.nmtr as usize] = 103;
                    self.mt103 = 1;
                }
                if mti >= self.mdmin && mti <= self.mdmax && self.mt104 == 0 {
                    self.nmtr += 1;
                    self.mtr[self.nmtr as usize] = 104;
                    self.mt104 = 1;
                }
                if mti >= self.mtmin && mti <= self.mtmax && self.mt105 == 0 {
                    self.nmtr += 1;
                    self.mtr[self.nmtr as usize] = 105;
                    self.mt105 = 1;
                }
                if mti >= self.m3min && mti <= self.m3max && self.mt106 == 0 {
                    self.nmtr += 1;
                    self.mtr[self.nmtr as usize] = 106;
                    self.mt106 = 1;
                }
                if mti >= self.m4min && mti <= self.m4max && self.mt107 == 0 {
                    self.nmtr += 1;
                    self.mtr[self.nmtr as usize] = 107;
                    self.mt107 = 1;
                }
            } else if mfi == 10 {
                nxn += 1;
            } else if mfi == 12 {
                nxn += 1;
                let mti = nint(self.dict[i + 3]);
                if mti == 3 {
                    self.nmtr += 1;
                }
                if mti == 3 {
                    self.mtr[self.nmtr as usize] = 3;
                }
            } else if mfi == 13 {
                nxn += 1;
            } else if mfi == 23 {
                nxn += 1;
                igam = 1;
                if c.iverf < 6 {
                    self.zain = 0.0;
                    self.awin = 0.0;
                }
                if self.nmtr == 0 {
                    self.mtr[1] = 501;
                }
                if self.nmtr == 0 {
                    self.nmtr = 1;
                }
                if self.mtr522 == 0 {
                    let mti = nint(self.dict[i + 3]);
                    if (534..=572).contains(&mti) {
                        self.nmtr += 1;
                        self.mtr[self.nmtr as usize] = 522;
                        self.mtr522 = 1;
                    }
                }
            }
            i += 6;
        }
        for i in 1..=self.nmtr as usize {
            self.mtrt[i] = 1000000;
        }
        for i in 1..=self.nmtr as usize {
            if self.mtr[i] == 1 {
                self.mtrt[i] = 1;
            }
            if self.mtr[i] == 501 {
                self.mtrt[i] = 1;
            }
        }
        if self.nmtr > 2 {
            let i2 = self.nmtr - 1;
            for i in 2..=i2 as usize {
                let j1 = i + 1;
                for j in j1..=self.nmtr as usize {
                    if self.mtr[j] <= self.mtr[i] {
                        let mtsave = self.mtr[i];
                        self.mtr[i] = self.mtr[j];
                        self.mtr[j] = mtsave;
                    }
                }
            }
        }
        self.dict = RV::default();
        nxn += 4;
        self.mfs = IV::new(nxn as usize);
        self.mts = IV::new(nxn as usize);
        self.ncs = IV::new(nxn as usize);
        self.mfs[1] = 1;
        self.mts[1] = 451;
        self.ncs[1] = self.ncards + 2;
        self.nxc = 1;
        if igam == 1 {
            return;
        }
        self.mfs[2] = 2;
        self.mts[2] = 151;
        self.ncs[2] = 4;
        self.nxc = 2;
    }

    // ================================================================ rdfil2
    fn rdfil2(&mut self, c: &mut Ctx, nrtot: &mut i32, intunr: &mut i32) {
        let therm = 0.0253e0f64;
        let eps = 1.0e-10f64;
        let elow = 1.0e-5f64;
        let elarge = 9.0e9f64;
        let small = 1.0e-6f64;
        let (mut nb, mut nw) = (0, 0);
        let nin = self.nin;

        self.eresl = elarge;
        self.eresh = 0.0;
        self.eresr = 0.0;
        self.eresu = elarge;
        self.eresm = elarge;

        let mut jnow = 1;
        c.contio(nin, 0, 0, &mut self.res.at(jnow), &mut nb, &mut nw);
        if c.mfh != 2 {
            c.io.so("\n mat has no file 2.");
            *nrtot = 0;
            self.lrp = 0;
            c.tosend(nin, 0, 0, &mut self.res.all());
            return;
        }
        let nis = c.n1h;
        self.nsect = 0;

        for i in 1..=nis {
            c.contio(nin, 0, 0, &mut self.res.at(jnow), &mut nb, &mut nw);
            jnow += nw;
            if jnow > self.maxres {
                c.error("rdfil2", "res storage exceeded", " ");
            }
            self.zai = c.c1h;
            let abn = c.c2h;
            self.lfw = c.l2h;
            let ner = c.n1h;

            for j in 1..=ner {
                c.contio(nin, 0, 0, &mut self.res.at(jnow), &mut nb, &mut nw);
                self.el = c.c1h;
                self.eh = c.c2h;
                let lru = c.l1h;
                let lrf = c.l2h;

                if self.el < self.eresl {
                    self.eresl = self.el;
                }
                if self.eh > self.eresh {
                    self.eresh = self.eh;
                }
                if lru <= 1 && self.eh > self.eresr {
                    self.eresr = self.eh;
                }
                if lru != 0 {
                    if (self.el - elow).abs() <= small {
                        self.nodes += 1;
                        self.enode[self.nodes] = self.el;
                    } else {
                        self.nodes += 1;
                        self.enode[self.nodes] = sigfig(self.el, 7, -1);
                        self.nodes += 1;
                        self.enode[self.nodes] = sigfig(self.el, 7, 1);
                        if lru == 2 {
                            if self.el < self.eresu {
                                self.eresu = self.el;
                            }
                            if self.nunr >= MAXUNR {
                                c.error("rdfil2", "storage in eunr exceeded.", " ");
                            }
                            self.nunr += 1;
                            self.eunr[self.nunr] = sigfig(self.el, 7, -1);
                            self.nunr += 1;
                            self.eunr[self.nunr] = sigfig(self.el, 7, 1);
                        }
                    }
                    self.nodes += 1;
                    self.enode[self.nodes] = sigfig(self.eh, 7, -1);
                    self.nodes += 1;
                    self.enode[self.nodes] = sigfig(self.eh, 7, 1);
                    if lru == 2 {
                        if self.eh < self.eresm {
                            self.eresm = self.eh;
                        }
                        if self.nunr >= MAXUNR {
                            c.error("rdfil2", "storage in eunr exceeded.", " ");
                        }
                        self.nunr += 1;
                        self.eunr[self.nunr] = sigfig(self.eh, 7, -1);
                        self.nunr += 1;
                        self.eunr[self.nunr] = sigfig(self.eh, 7, 1);
                    }
                }

                let nro = c.n1h;
                if lru == 1 {
                    self.res[jnow + 4] = nro as f64;
                } else {
                    self.res[jnow + 4] = self.lfw as f64;
                }
                let naps = c.n2h;
                self.res[jnow + 5] = naps as f64;
                let mut mode = lrf + 10 * (lru - 1);
                if lru == 0 {
                    mode = 0;
                }
                if lru == 2 {
                    self.lrp = 3;
                }
                self.nsect += 1;
                let ns = self.nsect as usize;
                self.isot[ns] = i;
                self.abnt[ns] = abn;
                self.elt[ns] = self.el;
                self.eht[ns] = self.eh;
                self.modet[ns] = mode;
                self.ibaset[ns] = jnow;
                jnow += nw;
                if jnow > self.maxres {
                    c.error("rdfil2", "res storage exceeded", " ");
                }

                *intunr = 5;
                if mode <= 3 && self.nmtres == 0 {
                    self.rdf2bw(c, nin, &mut jnow, lrf, nro, naps, mode);
                } else if mode <= 3 && self.nmtres > 0 {
                    self.samm.rdsammy(
                        c, nin, j, &mut jnow, nro, naps, mode, self.el, self.eh, &mut self.enode,
                        &mut self.nodes, NODMAX, &mut self.spin, &mut self.ascat, &mut self.res, self.maxres,
                    );
                    self.samm.ppsammy(c, j, &mut self.ncoef, &mut self.nresp);
                } else if mode == 4 {
                    self.rdf2aa(c, nin, &mut jnow, mode);
                } else if mode == 6 {
                    self.rdf2hy(c, nin, &mut jnow);
                } else if mode == 7 && self.nmtres > 0 {
                    self.samm.rdsammy(
                        c, nin, j, &mut jnow, nro, naps, mode, self.el, self.eh, &mut self.enode,
                        &mut self.nodes, NODMAX, &mut self.spin, &mut self.ascat, &mut self.res, self.maxres,
                    );
                    self.samm.ppsammy(c, j, &mut self.ncoef, &mut self.nresp);
                } else if mode == 11 && self.lfw == 0 {
                    self.rdf2u0(c, nin, &mut jnow, nro);
                } else if mode == 11 && self.lfw == 1 {
                    self.rdf2u1(c, nin, &mut jnow, nro);
                } else if mode == 12 {
                    self.rdf2u2(c, nin, &mut jnow, nro, intunr);
                } else {
                    c.error("rdfil2", "illegal resonance mode", " ");
                }

                if self.el < therm {
                    self.nodes += 1;
                    self.enode[self.nodes] = therm;
                }
            }
        }
        if self.eresr < self.eresl {
            self.eresr = self.eresl;
        }
        if self.eresr > self.eresh {
            self.eresr = self.eresh;
        }

        c.tofend(nin, 0, 0, &mut self.res.at(jnow));
        *nrtot = jnow;
        order(&mut self.enode, &mut self.nodes);
        if self.nunr == 0 {
            return;
        }

        order(&mut self.eunr, &mut self.nunr);
        let mut j = 0;
        self.el = 0.0;
        for i in 1..=self.nunr {
            let mut e = self.eunr[i];
            if e > (1.0 - eps) * self.eresu && e < (1.0 + eps) * self.eresh && e >= self.el {
                j += 1;
                if e < self.eresr {
                    e = -e;
                }
                self.eunr[j] = e;
                self.el = sigfig(e.abs(), 7, 1);
            }
        }
        self.nunr = j;
        if self.eresu < self.eresr {
            c.io.so(&format!(
                "\n mat has resolved-unresolved overlap in range\n     {}{}",
                ffmt::e1p(self.eresu, 12, 3),
                ffmt::e1p(self.eresr, 12, 3)
            ));
        }
        if self.eresm < self.eresh {
            c.io.so(&format!(
                "\n mat has unresolved-smooth overlap in range\n     {}{}",
                ffmt::e1p(self.eresm, 12, 3),
                ffmt::e1p(self.eresh, 12, 3)
            ));
        }
    }

    // ================================================================ rdf2bw
    fn rdf2bw(&mut self, c: &mut Ctx, nin: i32, jnow: &mut i32, lrf: i32, nro: i32, naps: i32, mode: i32) {
        // facts sets nothing for l > 4: the Fortran then keeps the values of
        // the previous call, so these live for the whole routine.
        let (mut sec, mut pec) = (0.0f64, 0.0f64);
        let (mut ser, mut per) = (0.0f64, 0.0f64);
        let rc1 = 0.123e0f64;
        let rc2 = 0.08e0f64;
        let third = 0.333333333e0f64;
        let gxmin = 1.0e-5f64;
        let zero = 0.0f64;
        let cwaven = sqrt(2.0 * AMASSN * AMU * EV) * 1.0e-12 / HBAR;
        let (mut nb, mut nw) = (0, 0);

        if nro != 0 {
            if mode != 2 {
                c.error("rdf2bw", "energy-dep scattering length", "only works for mlbw");
            }
            c.tab1io(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
            *jnow += nw;
            if *jnow > self.maxres {
                c.error("rdfil2", "res storage exceeded", " ");
            }
        }

        c.contio(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
        *jnow += nw;
        if *jnow > self.maxres {
            c.error("rdf2bw", "res storage exceeded", " ");
        }
        let nls = c.n1h;
        let mut lad = 0;
        if mode == 3 {
            lad = c.l1h;
        }
        if lad != 0 {
            c.mess("rdf2bw", "calculation of angular distribution not installed.", " ");
        }
        self.lssf = 0;
        self.spin = c.c1h;
        self.ascat = c.c2h;

        for _l in 1..=nls {
            c.listio(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
            let awri = self.res[*jnow];
            let mut ascatl = self.ascat;
            if lrf == 3 {
                ascatl = self.res[*jnow + 1];
            }
            if ascatl == zero {
                ascatl = self.ascat;
            }
            let arat = awri / (awri + 1.0);
            let aw = AMASSN * awri;
            let mut ra = rc1 * pow(aw, third) + rc2;
            let mut ral = ra;
            if naps == 1 {
                ral = ascatl;
                ra = self.ascat;
            }
            let mut jj = *jnow + nw;
            while nb != 0 {
                c.moreio(nin, 0, 0, &mut self.res.at(jj), &mut nb, &mut nw);
                jj += nw;
                if jj > self.maxres {
                    c.error("rdf2bw", "res storage exceeded", " ");
                }
            }
            let ll = nint(self.res[*jnow + 2]);
            let nrs = nint(self.res[*jnow + 5]);
            if nrs == 0 {
                let strng = format!("nrs=0 for SLBW/MLBW and l={}", ffmt::i(ll, 2));
                c.mess("rdf2bw", &strng, "malformed ENDF file, check evaluation");
                *jnow += 6;
                continue;
            }
            let ncyc = nint(self.res[*jnow + 4]) / nrs;
            let qx = self.res[*jnow + 1];
            self.lrx = nint(self.res[*jnow + 3]);
            if self.lrx != 0 {
                self.nsig = 5;
            }
            *jnow += 6;

            for _n in 1..=nrs {
                let jn = *jnow;
                if self.res[jn] > self.el && self.res[jn] <= self.eh {
                    if self.nodes >= NODMAX {
                        c.error("rdf2bw", "storage in enode exceeded.", " ");
                    }
                    self.nodes += 1;
                    let mut hw = self.res[jn + 2] / 2.0;
                    if mode == 3 {
                        hw = hw + (self.res[jn + 3] + self.res[jn + 4].abs() + self.res[jn + 5].abs()) / 2.0;
                    }
                    let mut ndig = 5;
                    if self.res[jn] > zero {
                        ndig = 2 + nint(log10(self.res[jn] / (hw / 10.0)));
                    }
                    if ndig < 5 {
                        ndig = 5;
                    }
                    if ndig > 9 {
                        ndig = 9;
                    }
                    self.enode[self.nodes] = sigfig(self.res[jn], ndig, 0);
                    if self.res[jn] + hw < self.eh {
                        if self.nodes >= NODMAX {
                            c.error("rdf2bw", "storage in enode exceeded.", " ");
                        }
                        let ehalf = self.res[jn] + hw;
                        self.nodes += 1;
                        self.enode[self.nodes] = sigfig(ehalf, ndig, 0);
                    }
                    if self.res[jn] - hw > self.el {
                        if self.nodes >= NODMAX {
                            c.error("rdf2bw", "storage in enode exceeded.", " ");
                        }
                        let ehalf = self.res[jn] - hw;
                        self.nodes += 1;
                        self.enode[self.nodes] = sigfig(ehalf, ndig, 0);
                    }
                }
                let rho = cwaven * arat * sqrt(self.res[jn].abs()) * ral;
                facts(ll, rho, &mut ser, &mut per);
                self.res[jj] = ser;
                self.res[jj + 1] = per;
                self.res[jj + 2] = 0.0;
                self.res[jj + 3] = 0.0;
                self.res[jj + 4] = 0.0;
                self.res[jj + 5] = 0.0;
                if self.lrx != 0 {
                    let rhoc = cwaven * arat * sqrt((self.res[jn] + qx / arat).abs()) * ra;
                    let mut lp = ll;
                    if ll == 0 {
                        lp = 2;
                    }
                    if ll == 2 {
                        lp = 0;
                    }
                    facts(lp, rhoc, &mut sec, &mut pec);
                    self.res[jj + 2] = self.res[jn + 2] - self.res[jn + 3] - self.res[jn + 4] - self.res[jn + 5];
                    if self.res[jj + 2] < gxmin * self.res[jn + 2] {
                        self.res[jj + 2] = 0.0;
                    }
                    if self.res[jn] < -qx / arat {
                        self.res[jj + 2] = 0.0;
                    }
                    self.res[jj + 3] = pec;
                }
                jj += 6;
                if jj > self.maxres {
                    c.error("rdf2bw", "res storage exceeded", " ");
                }
                *jnow += ncyc;
            }
            *jnow = jj;
        }
    }

    // ================================================================ rdf2aa
    fn rdf2aa(&mut self, c: &mut Ctx, nin: i32, jnow: &mut i32, _mode: i32) {
        let zero = 0.0f64;
        let (mut nb, mut nw) = (0, 0);
        c.contio(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
        *jnow += nw;
        if *jnow > self.maxres {
            c.error("rdf2aa", "res storage exceeded", " ");
        }
        let nls = c.n1h;
        // note: the Fortran declares a local lssf here, so the module
        // variable is not changed
        self.spin = c.c1h;
        self.ascat = c.c2h;

        c.listio(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
        *jnow += nw;
        while nb != 0 {
            if *jnow > self.maxres {
                c.error("rdf2aa", "res storage exceeded", " ");
            }
            c.moreio(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
            *jnow += nw;
        }

        for l in 1..=nls {
            c.contio(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
            *jnow += nw;
            if *jnow > self.maxres {
                c.error("rdf2aa", "res storage exceeded", " ");
            }
            let njs = c.n1h;
            for n in 1..=njs {
                c.listio(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
                let nen = c.n2h;
                let jen = 12;
                let nloc = *jnow - 6;
                *jnow += nw;
                while nb != 0 {
                    if *jnow > self.maxres {
                        c.error("rdf2aa", "res storage exceeded", " ");
                    }
                    c.moreio(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
                    *jnow += nw;
                }
                if n == 1 && l == 1 {
                    for ien in 1..=nen {
                        let ener = self.res[nloc + ien * jen];
                        if ener > self.el && ener <= self.eh {
                            if self.nodes >= NODMAX {
                                c.error("rdf2aa", "storage in enode exceeded.", " ");
                            }
                            self.nodes += 1;
                            let hw = self.res[nloc + ien * jen + 1];
                            let mut ndig = 5;
                            if ener > zero {
                                ndig = 2 + nint(log10(ener / (hw / 10.0)));
                            }
                            if ndig < 5 {
                                ndig = 5;
                            }
                            if ndig > 9 {
                                ndig = 9;
                            }
                            self.enode[self.nodes] = sigfig(ener, ndig, 0);
                            if ener + hw <= self.eh {
                                if self.nodes >= NODMAX {
                                    c.error("rdf2aa", "storage in enode exceeded.", " ");
                                }
                                let ehalf = ener + hw;
                                self.nodes += 1;
                                self.enode[self.nodes] = sigfig(ehalf, ndig, 0);
                            }
                            if ener - hw > self.el {
                                if self.nodes >= NODMAX {
                                    c.error("rdf2aa", "storage in enode exceeded.", " ");
                                }
                                let ehalf = ener - hw;
                                self.nodes += 1;
                                self.enode[self.nodes] = sigfig(ehalf, ndig, 0);
                            }
                        }
                    }
                }
            }
        }
    }

    // ================================================================ rdf2hy
    fn rdf2hy(&mut self, c: &mut Ctx, nin: i32, jnow: &mut i32) {
        let zero = 0.0f64;
        let (mut nb, mut nw) = (0, 0);
        c.contio(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
        *jnow += nw;
        if *jnow > self.maxres {
            c.error("rdf2hy", "res storage exceeded", " ");
        }
        let nls = c.n1h;
        let lad = c.l1h;
        if lad != 0 {
            c.mess("rdf2hy", "calculation of angular distribution not installed.", " ");
        }
        self.lssf = 0;
        self.spin = c.c1h;
        self.ascat = c.c2h;

        c.contio(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
        let _nire = c.n1h;
        let ncre = c.n2h;
        if _nire + ncre > 0 {
            c.mess("rdf2hy", "hybrid competing reactions not yet added to file 3.", " ");
        }
        *jnow += nw;
        c.contio(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
        *jnow += nw;
        c.listio(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
        *jnow += nw;
        if ncre != 0 {
            for _ncr in 1..=ncre {
                for _lil in 1..=4 {
                    c.tab1io(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
                    *jnow += nw;
                }
            }
        }
        if nls == 0 {
            return;
        }
        for _l in 1..=nls {
            c.contio(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
            let nss = c.n1h;
            *jnow += nw;
            for _ns in 1..=nss {
                c.contio(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
                let njs = c.n1h;
                *jnow += nw;
                for _nj in 1..=njs {
                    c.listio(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
                    let lbk = c.l1h;
                    let lps = c.l2h;
                    let nlsj = c.n2h;
                    let mut ll = *jnow + 6;
                    for _ir in 1..=nlsj {
                        let mut hw = 0.0;
                        for ig in 1..=7 {
                            hw = hw + self.res[ll + ig];
                        }
                        hw = hw / 2.0;
                        let enow = self.res[ll];
                        let mut ndig = 5;
                        if enow > zero {
                            ndig = 2 + nint(log10(enow / (hw / 10.0)));
                        }
                        if ndig < 5 {
                            ndig = 5;
                        }
                        if ndig > 9 {
                            ndig = 9;
                        }
                        if self.nodes + 3 >= NODMAX {
                            c.error("rdf2hy", "storage in enode exceeded.", " ");
                        }
                        if enow > self.el && enow <= self.eh {
                            self.nodes += 1;
                            self.enode[self.nodes] = sigfig(enow, ndig, 0);
                            if enow - hw > self.el {
                                let ehalf = enow - hw;
                                self.nodes += 1;
                                self.enode[self.nodes] = sigfig(ehalf, ndig, 0);
                            }
                            if enow + hw < self.eh {
                                let ehalf = enow + hw;
                                self.nodes += 1;
                                self.enode[self.nodes] = sigfig(ehalf, ndig, 0);
                            }
                        }
                        ll += 12;
                    }
                    *jnow += nw;
                    if lbk != 0 {
                        c.tab1io(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
                        *jnow += nw;
                        c.tab1io(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
                        *jnow += nw;
                    }
                    if lps != 0 {
                        c.tab1io(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
                        *jnow += nw;
                        c.tab1io(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
                        *jnow += nw;
                    }
                }
            }
        }
    }

    // ================================================================ rdf2u0
    fn rdf2u0(&mut self, c: &mut Ctx, nin: i32, jnow: &mut i32, nro: i32) {
        let (mut nb, mut nw) = (0, 0);
        if nro == 1 {
            c.tab1io(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
            self.res[*jnow + 1] = -(nro as f32) as f64;
            let mut jj = *jnow + nw;
            while nb != 0 {
                c.moreio(nin, 0, 0, &mut self.res.at(jj), &mut nb, &mut nw);
                jj += nw;
                if jj > self.maxres {
                    c.error("rdf2u0", "storage in res exceeded", " ");
                }
            }
            *jnow = jj;
        }
        c.contio(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
        *jnow += nw;
        if *jnow > self.maxres {
            c.error("rdf2u0", "res storage exceeded", " ");
        }
        let nls = c.n1h;
        self.lssf = c.l1h;
        self.spin = c.c1h;
        self.ascat = c.c2h;

        for _l in 1..=nls {
            c.listio(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
            let mut jj = *jnow + nw;
            while nb != 0 {
                c.moreio(nin, 0, 0, &mut self.res.at(jj), &mut nb, &mut nw);
                jj += nw;
                if jj > self.maxres {
                    c.error("rdf2u0", "res storage exceeded", " ");
                }
            }
            *jnow = jj;

            let mut ener = self.el;
            while ener < self.eh {
                let mut enut = 0.0;
                let mut ii = 0;
                // Fortran: do while (ii.le.ngridu .and. enut.lt.ener+ener/100)
                // With ii=78 and the condition still true it reads
                // egridu(79), past the end of the array (undefined value).
                while ii <= NGRIDU && enut < ener + ener / 100.0 {
                    ii += 1;
                    enut = egridu(ii);
                }
                ener = enut;
                if ener < self.eh {
                    if self.nodes >= NODMAX {
                        c.error("rdf2u0", "storage in enode exceeded.", " ");
                    }
                    self.nodes += 1;
                    self.enode[self.nodes] = sigfig(ener, 7, 0);
                    if self.nunr >= MAXUNR {
                        c.error("rdf2u0", "storage in eunr exceeded.", " ");
                    }
                    self.nunr += 1;
                    self.eunr[self.nunr] = sigfig(ener, 7, 0);
                }
            }
        }
    }

    // ================================================================ rdf2u1
    fn rdf2u1(&mut self, c: &mut Ctx, nin: i32, jnow: &mut i32, nro: i32) {
        let wide = 1.26e0f64;
        let (mut nb, mut nw) = (0, 0);
        if nro == 1 {
            c.tab1io(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
            self.res[*jnow + 1] = -(nro as f32) as f64;
            let mut jj = *jnow + nw;
            while nb != 0 {
                c.moreio(nin, 0, 0, &mut self.res.at(jj), &mut nb, &mut nw);
                jj += nw;
                if jj > self.maxres {
                    c.error("rdf2u1", "storage in res exceeded", " ");
                }
            }
            *jnow = jj;
        }
        c.listio(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
        let nloc = *jnow + 5;
        *jnow += nw;
        if *jnow > self.maxres {
            c.error("rdf2u1", "res storage exceeded", " ");
        }
        while nb != 0 {
            c.moreio(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
        }
        self.lssf = c.l1h;
        let ne = c.n1h;
        let nls = c.n2h;
        self.spin = c.c1h;
        self.ascat = c.c2h;

        for n in 1..=ne {
            let mut ener = self.res[nloc + n];
            if ener >= self.el && ener < self.eh {
                let enex = self.res[nloc + n + 1];
                let mut ifill = 0;
                if enex > wide * ener {
                    ifill = 1;
                }
                if self.nodes >= NODMAX {
                    c.error("rdf2u1", "storage in enode exceeded.", " ");
                }
                self.nodes += 1;
                self.enode[self.nodes] = sigfig(ener, 7, 0);
                if self.nunr >= MAXUNR {
                    c.error("rdf2u1", "storage in eunr exceeded.", " ");
                }
                self.nunr += 1;
                self.eunr[self.nunr] = sigfig(ener, 7, 0);
                while ifill != 0 {
                    let mut enut = 0.0;
                    let mut ii = 0;
                    while ii < NGRIDU && enut < ener + ener / 100.0 {
                        ii += 1;
                        enut = egridu(ii);
                    }
                    ener = enut;
                    if ener >= enex {
                        ifill = 0;
                    }
                    if ifill == 1 {
                        if self.nodes >= NODMAX {
                            c.error("rdf2u1", "storage in enode exceeded.", " ");
                        }
                        self.nodes += 1;
                        self.enode[self.nodes] = sigfig(ener, 7, 0);
                        if self.nunr >= MAXUNR {
                            c.error("rdf2u1", "storage in eunr exceeded.", " ");
                        }
                        self.nunr += 1;
                        self.eunr[self.nunr] = sigfig(ener, 7, 0);
                    }
                }
            }
        }

        for _l in 1..=nls {
            c.contio(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
            *jnow += nw;
            if *jnow > self.maxres {
                c.error("rdf2u1", "res storage exceeded", " ");
            }
            let njs = c.n1h;
            for _n in 1..=njs {
                c.listio(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
                *jnow += nw;
                while nb != 0 {
                    if *jnow > self.maxres {
                        c.error("rdf2u1", "res storage exceeded", " ");
                    }
                    c.moreio(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
                    *jnow += nw;
                }
            }
        }
    }

    // ================================================================ rdf2u2
    fn rdf2u2(&mut self, c: &mut Ctx, nin: i32, jnow: &mut i32, nro: i32, intunr: &mut i32) {
        let wide = 1.26e0f64;
        let (mut nb, mut nw) = (0, 0);
        if nro == 1 {
            c.tab1io(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
            self.res[*jnow + 1] = -(nro as f32) as f64;
            let mut jj = *jnow + nw;
            while nb != 0 {
                c.moreio(nin, 0, 0, &mut self.res.at(jj), &mut nb, &mut nw);
                jj += nw;
                if jj > self.maxres {
                    c.error("rdf2u2", "storage in res exceeded", " ");
                }
            }
            *jnow = jj;
        }
        c.contio(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
        *jnow += nw;
        if *jnow > self.maxres {
            c.error("rdf2u2", "res storage exceeded", " ");
        }
        let nls = c.n1h;
        self.lssf = c.l1h;
        self.spin = c.c1h;
        self.ascat = c.c2h;

        for l in 1..=nls {
            c.contio(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
            *jnow += nw;
            if *jnow > self.maxres {
                c.error("rdf2u2", "res storage exceeded", " ");
            }
            let njs = c.n1h;
            for n in 1..=njs {
                c.listio(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
                *intunr = c.l1h;
                let nen = c.n2h;
                let jen = 6;
                let nloc = *jnow + 6;
                *jnow += nw;
                while nb != 0 {
                    if *jnow > self.maxres {
                        c.error("rdf2u2", "res storage exceeded", " ");
                    }
                    c.moreio(nin, 0, 0, &mut self.res.at(*jnow), &mut nb, &mut nw);
                    *jnow += nw;
                }
                if n == 1 && l == 1 {
                    for ien in 1..=nen {
                        let mut ener = self.res[nloc + (ien - 1) * jen];
                        if ener >= self.el && ener < self.eh {
                            let enex = self.res[nloc + (ien - 1) * jen + jen];
                            let mut ifill = 0;
                            if enex > wide * ener {
                                ifill = 1;
                            }
                            if self.nodes >= NODMAX {
                                c.error("rdf2u2", "storage in enode exceeded.", " ");
                            }
                            self.nodes += 1;
                            self.enode[self.nodes] = sigfig(ener, 7, 0);
                            if self.nunr >= MAXUNR {
                                c.error("rdf2u2", "storage in eunr exceeded.", " ");
                            }
                            self.nunr += 1;
                            self.eunr[self.nunr] = sigfig(ener, 7, 0);
                            while ifill != 0 {
                                let mut enut = 0.0;
                                let mut ii = 0;
                                while ii < NGRIDU && enut < ener + ener / 1000.0 {
                                    ii += 1;
                                    enut = egridu(ii);
                                }
                                ener = enut;
                                if ener >= enex {
                                    ifill = 0;
                                }
                                if ifill == 1 {
                                    if self.nodes >= NODMAX {
                                        c.error("rdf2u2", "storage in enode exceeded.", " ");
                                    }
                                    self.nodes += 1;
                                    self.enode[self.nodes] = sigfig(ener, 7, 0);
                                    if self.nunr >= MAXUNR {
                                        c.error("rdf2u2", "storage in eunr exceeded.", " ");
                                    }
                                    self.nunr += 1;
                                    self.eunr[self.nunr] = sigfig(ener, 7, 0);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    // ================================================================ genunr
    fn genunr(&mut self, c: &mut Ctx, nin: i32, intunr: i32) {
        let big = 1.0e10f64;
        let zero = 0.0f64;
        let mut sigu = RV::new(5);
        let nw0 = 6 * self.nunr + 13;
        let mut sunr = RV::new(nw0 as usize + SLACK);
        sunr[1] = self.za;
        sunr[2] = self.awr;
        sunr[3] = self.lssf as f64;
        sunr[4] = 0.0;
        sunr[5] = 0.0;
        sunr[6] = intunr as f64;
        sunr[7] = self.tempr;
        sunr[8] = 0.0;
        sunr[9] = 5.0;
        sunr[10] = 1.0;
        sunr[11] = (1 + 6 * self.nunr) as f64;
        sunr[12] = self.nunr as f64;
        sunr[13] = big;

        let mut e = 0.0;
        self.csunr1(c, e, &mut sigu);
        self.csunr2(c, e, &mut sigu);
        let mut l = 14;
        for i in 1..=self.nunr {
            e = self.eunr[i];
            sunr[l] = e;
            e = e.abs();
            for k in 1..=4 {
                sunr[l + k] = 0.0;
            }
            for j in 1..=self.nsect as usize {
                let mode = self.modet[j];
                if mode >= 11 && e >= self.elt[j] && e < self.eht[j] {
                    let abn = self.abnt[j];
                    self.isect = j as i32;
                    if mode == 11 {
                        self.csunr1_ref(c, &mut e, &mut sigu);
                    }
                    if mode == 12 {
                        self.csunr2_ref(c, &mut e, &mut sigu);
                    }
                    for k in 1..=4 {
                        sunr[l + k] = sunr[l + k] + abn * sigu[k];
                        sunr[l + k] = sigfig(sunr[l + k], 7, 0);
                    }
                }
            }
            sunr[l + 5] = sunr[l + 1];
            l += 6;
        }
        if self.lssf != 0 {
            self.sunr = Some(sunr);
            return;
        }

        let mut scr = RV::new((c.npage + 50) as usize + SLACK);
        let mut mt1 = 0;
        c.findf(self.mata, 3, 0, nin);
        c.mfh = 3;
        let (mut nb, mut nw) = (0, 0);
        while c.mfh != 0 {
            c.contio(nin, 0, 0, &mut scr.all(), &mut nb, &mut nw);
            if c.mfh != 0 {
                if mt1 == 0 {
                    mt1 = c.mth;
                }
                let mut ix = 0;
                if c.mth == 1 {
                    ix = 1;
                }
                if c.mth == 2 {
                    ix = 2;
                }
                if c.mth == 18 {
                    ix = 3;
                }
                if c.mth == 102 {
                    ix = 4;
                }
                if ix != 0 {
                    e = 0.0;
                    let (mut enext, mut idis, mut bkg) = (0.0, 0, 0.0);
                    c.gety1(e, &mut enext, &mut idis, &mut bkg, nin, &mut scr.all());
                    l = 14;
                    for _i in 1..=self.nunr {
                        e = sunr[l];
                        if e >= zero {
                            c.gety1(e, &mut enext, &mut idis, &mut bkg, nin, &mut scr.all());
                            sunr[l + ix] = sunr[l + ix] + bkg;
                            sunr[l + ix] = sigfig(sunr[l + ix], 7, 0);
                            if ix == 1 {
                                sunr[l + 5] = sunr[l + 5] + bkg;
                                sunr[l + 5] = sigfig(sunr[l + 5], 7, 0);
                            }
                        }
                        l += 6;
                    }
                }
                c.tosend(nin, 0, 0, &mut scr.all());
            }
        }
        c.findf(self.mata, 3, mt1, nin);
        self.sunr = Some(sunr);
    }

    // ================================================================ sigunr
    pub fn sigunr(&self, e: f64, sigu: &mut [f64; 5], thr6: f64) {
        let sunr = self.sunr.as_ref().unwrap();
        let intunr = nint(sunr[6]);
        let l = 14;
        let mut en = 0.0;
        let mut i = 1;
        let mut l2 = 0;
        while i < self.nunr && e >= en {
            i += 1;
            l2 = l + 6 * (i - 1);
            en = sunr[l2].abs();
        }
        let l1 = l2 - 6;
        let e1 = sunr[l1].abs();
        let e2 = sunr[l2].abs();
        for i in 1..=4 {
            if e >= e1 && e <= e2 {
                crate::endf::terp1(e1, sunr[l1 + i], e2, sunr[l2 + i], e, &mut sigu[i as usize], intunr, thr6);
            } else {
                sigu[i as usize] = 0.0;
            }
        }
    }
}

const NGRIDU: i32 = 78;
const EGRIDU: [f64; 78] = [
    1.0e1, 1.25e1, 1.5e1, 1.7e1, 2.0e1, 2.5e1, 3.0e1, 3.5e1, 4.0e1, 5.0e1, 6.0e1, 7.2e1, 8.5e1, 1.0e2,
    1.25e2, 1.5e2, 1.7e2, 2.0e2, 2.5e2, 3.0e2, 3.5e2, 4.0e2, 5.0e2, 6.0e2, 7.2e2, 8.5e2, 1.0e3,
    1.25e3, 1.5e3, 1.7e3, 2.0e3, 2.5e3, 3.0e3, 3.5e3, 4.0e3, 5.0e3, 6.0e3, 7.2e3, 8.5e3, 1.0e4,
    1.25e4, 1.5e4, 1.7e4, 2.0e4, 2.5e4, 3.0e4, 3.5e4, 4.0e4, 5.0e4, 6.0e4, 7.2e4, 8.5e4, 1.0e5,
    1.25e5, 1.5e5, 1.7e5, 2.0e5, 2.5e5, 3.0e5, 3.5e5, 4.0e5, 5.0e5, 6.0e5, 7.2e5, 8.5e5, 1.0e6,
    1.25e6, 1.5e6, 1.7e6, 2.0e6, 2.5e6, 3.0e6, 3.5e6, 4.0e6, 5.0e6, 6.0e6, 7.2e6, 8.5e6,
];

/// egridu(ii); index 79 is past the end of the Fortran array (rdf2u0 can
/// read it). That read is undefined in Fortran; here it stops the loop.
fn egridu(ii: i32) -> f64 {
    if ii >= 1 && ii <= NGRIDU {
        EGRIDU[(ii - 1) as usize]
    } else {
        f64::INFINITY
    }
}

// ===================================================================== order
pub fn order(x: &mut RV, n: &mut i32) {
    let small = 1.0e-10f64;
    if *n <= 2 {
        return;
    }
    let mut m = *n;
    let mut i = 0;
    while i < m - 1 {
        i += 1;
        let mut j = i;
        while j < m {
            j += 1;
            if x[j] < x[i] {
                let tsave = x[j];
                x[j] = x[i];
                x[i] = tsave;
            }
        }
        if i > 1 && (x[i] - x[i - 1]).abs() <= small * x[i] {
            m -= 1;
            if i >= m {
                *n = m;
                return;
            }
            for k in i..=m {
                x[k] = x[k + 1];
            }
            i -= 1;
        }
    }
    if (x[m] - x[m - 1]).abs() <= small * x[m] {
        m -= 1;
    }
    *n = m;
}

// ===================================================================== facts
pub fn facts(l: i32, rho: f64, se: &mut f64, pe: &mut f64) {
    let r2 = rho * rho;
    if l == 0 {
        *se = 0.0;
        *pe = rho;
    } else if l == 1 {
        let den = 1.0 + r2;
        *pe = r2 * rho / den;
        *se = -1.0 / den;
    } else if l == 2 {
        let r4 = r2 * r2;
        let den = 3.0 * r2 + r4 + 9.0;
        *pe = r4 * rho / den;
        *se = -(18.0 + 3.0 * r2) / den;
    } else if l == 3 {
        let r4 = r2 * r2;
        let r6 = r4 * r2;
        let den = 225.0 + 45.0 * r2 + 6.0 * r4 + r6;
        *pe = r6 * rho / den;
        *se = -(675.0 + 90.0 * r2 + 6.0 * r4) / den;
    } else if l == 4 {
        let r4 = r2 * r2;
        let r6 = r4 * r2;
        let r8 = r4 * r4;
        let den = 11025.0 + 1575.0 * r2 + 135.0 * r4 + 10.0 * r6 + r8;
        *pe = r8 * rho / den;
        *se = -(44100.0 + 4725.0 * r2 + 270.0 * r4 + 10.0 * r6) / den;
    }
}

// ======================================================================
// lunion, resxs, sigma
// ======================================================================

impl Reconr {
    /// Unionize the energy grids of all sections and linearize.
    fn lunion(&mut self, c: &mut Ctx, nin: i32, nout: i32, ngrid: i32, ngr: &mut i32) {
        let elow = 1.0e-5f64;
        let up = 1.001e0f64;
        let dn = 0.999e0f64;
        let therm = 0.0253e0f64;
        let small = 1.0e-9f64;
        let ten = 10.0e0f64;
        let ssmall = 1.0e-30f64;
        let stpmin = 1.001e0f64;
        let ovfact = 5.3e0f64;
        let trange = 0.4999e0f64;
        let emax = 19.0e6f64;
        let zero = 0.0f64;

        // `elim` has an initial value, so it is a SAVEd local in Fortran.
        if self.eresr < self.sv_elim() {
            self.sv.elim = Some(self.eresr);
        }
        let elim = self.sv_elim();
        let ndim: i32 = 50;
        let mut nmf12 = 0;
        let mut bufo = vec![0.0f64; NBUF as usize];
        let mut bufn = vec![0.0f64; NBUF as usize];
        let mut x = RV::new(ndim as usize + 2);
        let mut y = RV::new(ndim as usize + 2);
        let mut scr = RV::new((c.npage + 50) as usize + SLACK);
        let mut aa = [0.0f64; 1];

        let stpmax = 1.0 + sqrt(ovfact * self.err);

        let mut iold = 14;
        let mut inew = 15;
        c.openz(-iold, 1);
        c.openz(-inew, 1);
        c.repoz(-iold);
        c.repoz(-inew);
        if self.nodes > 0 {
            for i in 1..=self.nodes {
                let mut j = i;
                if i == self.nodes {
                    j = -i;
                }
                aa[0] = self.enode[i];
                c.loada(j, &aa, 1, iold, &mut bufo, NBUF);
            }
        }
        self.ngo = self.nodes;

        c.repoz(nout);
        c.nsh = 0;
        let mut mfl = 0;
        let (mut nb, mut nw) = (0i32, 0i32);
        let mut nss;
        let mut iss = 0;
        let mut awrx = 0.0f64;
        let (mut in_, mut ig, mut ir, mut jr, mut lr, mut npr, mut ibase, mut nbta, mut inta) =
            (0i32, 0i32, 0i32, 0i32, 0i32, 0i32, 0i32, 0i32, 0i32);
        let mut intn = 0i32;
        let (mut er, mut sr, mut ernext, mut srnext, mut eg) = (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
        let (mut enl, mut snl, mut erl, mut srl, mut et) = (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
        let (mut en, mut sn, mut xm, mut ym, mut yl) = (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
        let mut i = 0i32;
        let mut test;
        let mut lbl = 110;
        loop {
            match lbl {
                110 => {
                    c.contio(nin, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                    nss = 1;
                    if c.mfh == 10 {
                        nss = c.n1h;
                    }
                    iss = nss;
                    if c.math > 0 {
                        lbl = 120;
                        continue;
                    }
                    c.amend(nout, 0);
                    c.atend(nout, 0);
                    lbl = 410;
                }
                120 => {
                    if c.mfh > 0 {
                        lbl = 130;
                        continue;
                    }
                    if mfl == 3 || mfl == 10 || mfl == 13 || mfl == 23 {
                        c.afend(nout, 0);
                    }
                    if mfl == 12 && nmf12 != 0 {
                        c.afend(nout, 0);
                    }
                    mfl = 0;
                    lbl = 110;
                }
                130 => {
                    let mf = c.mfh;
                    if mf == 3 || mf == 10 || mf == 12 || mf == 13 || mf == 23 {
                        lbl = 140;
                        continue;
                    }
                    c.tofend(nin, 0, 0, &mut scr.all());
                    lbl = 110;
                }
                140 => {
                    mfl = c.mfh;
                    let mth = c.mth;
                    let mfh = c.mfh;
                    lbl = 180;
                    if mth == 501 {
                        lbl = 150;
                    } else if mth == 522 && self.mtr522 > 0 {
                        lbl = 150;
                    } else if mth == 460 {
                        lbl = 150;
                    } else if mfh == 12 && nint(scr[3]) != 1 {
                        lbl = 150;
                    } else {
                        if mfh == 12 {
                            scr[5] = 1.0;
                        }
                        if mfh == 13 {
                            scr[5] = 1.0;
                        }
                        if mfh != 3 {
                            lbl = 180;
                        } else if mth == 1 || mth == 3 {
                            lbl = 150;
                        } else if mth == 4 && self.mtr4 > 0 {
                            lbl = 150;
                        } else if mth == 103 && self.mt103 > 0 {
                            lbl = 150;
                        } else if mth == 104 && self.mt104 > 0 {
                            lbl = 150;
                        } else if mth == 105 && self.mt105 > 0 {
                            lbl = 150;
                        } else if mth == 106 && self.mt106 > 0 {
                            lbl = 150;
                        } else if mth == 107 && self.mt107 > 0 {
                            lbl = 150;
                        } else if mth == 101 || mth == 120 || mth == 151 {
                            lbl = 150;
                        } else if (251..=300).contains(&mth) && mth != 261 {
                            lbl = 150;
                        } else if mth == 18 && self.mtr18 > 0 {
                            lbl = 150;
                        } else if c.iverf < 6 {
                            if (451..700).contains(&mth) || mth >= 800 {
                                lbl = 150;
                            }
                        } else if (451..600).contains(&mth) || (mth > 850 && mth <= 870) || mth > 891 {
                            lbl = 150;
                        }
                    }
                }
                150 => {
                    c.tosend(nin, 0, 0, &mut scr.all());
                    lbl = 110;
                }
                180 => {
                    if c.mfh == 12 {
                        nmf12 += 1;
                    }
                    c.contio(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                    if c.mfh != 23 && self.awin != zero {
                        awrx = c.c2h / self.awin;
                    }
                    lbl = 181;
                }
                181 => {
                    c.tab1io(nin, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                    if c.mfh == 23 {
                        lbl = 190;
                        continue;
                    }
                    let qx = c.c2h;
                    if c.mth == 19 && self.mtr18 > 0 {
                        self.q18 = qx;
                    }
                    let mut thrx = if self.awin != 0.0 { -qx * (awrx + 1.0) / awrx } else { -qx };
                    c.thr6 = thrx;
                    if c.thr6 < zero {
                        c.thr6 = 0.0;
                    }
                    if thrx <= zero {
                        lbl = 190;
                        continue;
                    }
                    let thrxx = sigfig(thrx, 7, 1);
                    let mut l = 7 + 2 * nint(scr[5]);
                    if scr[l + 1] != zero {
                        let text = format!("xsec nonzero at threshold for mt={}", ffmt::i(c.mth, 3));
                        c.mess("lunion", &text, "adjusted using jump in xsec");
                    }
                    if scr[l] >= thrxx {
                        lbl = 190;
                        continue;
                    }
                    thrx = thrxx;
                    c.io.so(&format!(
                        "\n changed threshold from{} to{} for mt{}.",
                        ffmt::e1p(scr[l], 13, 6),
                        ffmt::e1p(thrx, 13, 6),
                        ffmt::i(c.mth, 3)
                    ));
                    scr[l] = thrx;
                    let lsave = l;
                    while scr[l + 2] <= scr[l] {
                        if l > lsave + 20 {
                            c.error("lunion", "ill-behaved threshold.", " ");
                        }
                        scr[l + 2] = sigfig(scr[l], 7, 1);
                        l += 2;
                    }
                    lbl = 190;
                }
                190 => {
                    c.tab1io(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                    in_ = 0;
                    ig = 1;
                    ir = 0;
                    jr = 1;
                    lr = 0;
                    npr = nint(scr[6]);
                    ibase = 6 + 2 * nint(scr[5]);
                    nbta = 1;
                    lbl = 205;
                }
                205 => {
                    ir += 1;
                    if lr + ir > nbta {
                        jr += 1;
                    }
                    nbta = nint(scr[5 + 2 * jr]);
                    inta = nint(scr[6 + 2 * jr]);
                    er = scr[ibase + ir * 2 - 1];
                    sr = scr[ibase + ir * 2];
                    eg = er;
                    let mth = c.mth;
                    if mth == 2 || mth == 18 || mth == 19 || mth == 102 {
                        lbl = 210;
                        continue;
                    }
                    if ibase + 2 * ir < nw || nb == 0 {
                        lbl = 207;
                        continue;
                    }
                    c.moreio(nin, nout, 0, &mut scr.at(ibase + 1), &mut nb, &mut nw);
                    lr += ir;
                    ir = 0;
                    nw = ibase + nw;
                    lbl = 207;
                }
                207 => {
                    ernext = scr[ibase + ir * 2 + 1];
                    srnext = scr[ibase + ir * 2 + 2];
                    if sr < ssmall && srnext < ssmall && ir < npr - 1 {
                        lbl = 205;
                        continue;
                    }
                    if (er - ernext).abs() > small * er {
                        lbl = 210;
                        continue;
                    }
                    er = sigfig(er, 7, 0);
                    ernext = sigfig(ernext, 7, 1);
                    scr[ibase + ir * 2 + 1] = ernext;
                    lbl = 210;
                }
                210 => {
                    if ig > self.ngo {
                        lbl = 220;
                        continue;
                    }
                    c.finda(ig, &mut aa, 1, iold, &mut bufo, NBUF);
                    eg = aa[0];
                    if eg.abs() >= er * (1.0 - small) {
                        lbl = 220;
                        continue;
                    }
                    enl = eg;
                    in_ += 1;
                    aa[0] = enl;
                    c.loada(in_, &aa, 1, inew, &mut bufn, NBUF);
                    ig += 1;
                    lbl = 210;
                }
                220 => {
                    test = elow;
                    test = sigfig(test, 7, 1);
                    if er < test {
                        lbl = 225;
                        continue;
                    }
                    if er - sigfig(eg.abs(), 7, -1) > 1.0e-8 * er {
                        // 222
                        er = eg.abs();
                        lbl = 225;
                        continue;
                    }
                    if sr == zero {
                        lbl = 225;
                        continue;
                    }
                    enl = sigfig(er, 7, 0);
                    er = sigfig(er, 7, 1);
                    in_ += 1;
                    aa[0] = enl;
                    c.loada(in_, &aa, 1, inew, &mut bufn, NBUF);
                    lbl = 225;
                }
                225 => {
                    enl = er;
                    snl = sr;
                    erl = er;
                    srl = sr;
                    if ig < self.ngo && eg.abs() <= er * (1.0 + small) {
                        ig += 1;
                        c.finda(ig, &mut aa, 1, iold, &mut bufo, NBUF);
                        eg = aa[0];
                    }
                    ir += 1;
                    er = scr[ibase + ir * 2 - 1];
                    sr = scr[ibase + ir * 2];
                    if lr + ir > nbta {
                        jr += 1;
                        nbta = nint(scr[5 + 2 * jr]);
                        inta = nint(scr[6 + 2 * jr]);
                    }
                    if ir < nbta {
                        ernext = scr[ibase + 2 * ir + 1];
                        srnext = scr[ibase + 2 * ir + 2];
                        if (er - ernext).abs() <= small * er {
                            er = sigfig(er, 7, -1);
                            ernext = sigfig(ernext, 7, 1);
                            scr[ibase + 2 * ir + 1] = ernext;
                        }
                    }
                    et = 0.0;
                    lbl = 240;
                }
                240 => {
                    if lr + ir > c.n2h {
                        lbl = 350;
                        continue;
                    }
                    if ig >= self.ngo {
                        lbl = 250;
                        continue;
                    }
                    if eg.abs() >= er * (1.0 - small) {
                        lbl = 250;
                        continue;
                    }
                    en = eg.abs();
                    c.terp1(erl.abs(), srl, er, sr, en, &mut sn, inta);
                    if eg < zero {
                        en = -en;
                    }
                    intn = inta;
                    ig += 1;
                    c.finda(ig, &mut aa, 1, iold, &mut bufo, NBUF);
                    eg = aa[0];
                    lbl = 300;
                }
                250 => {
                    if inta != 1 || lr + ir >= c.n2h || er < (1.0 + small) * enl {
                        lbl = 255;
                        continue;
                    }
                    if (et - enl).abs() < small * et {
                        // 252
                        er = sigfig(er, 7, 1);
                        lbl = 255;
                        continue;
                    }
                    en = er;
                    et = sigfig(en, 7, -1);
                    en = et;
                    sn = snl;
                    intn = inta;
                    erl = en;
                    srl = sr;
                    lbl = 300;
                }
                255 => {
                    en = er;
                    sn = sr;
                    intn = inta;
                    erl = er;
                    srl = sr;
                    lbl = 260;
                }
                260 => {
                    ir += 1;
                    if lr + ir > c.n2h {
                        lbl = 280;
                        continue;
                    }
                    er = scr[ibase + ir * 2 - 1];
                    sr = scr[ibase + ir * 2];
                    if lr + ir > nbta {
                        jr += 1;
                        nbta = nint(scr[5 + 2 * jr]);
                        inta = nint(scr[6 + 2 * jr]);
                    }
                    if ibase + 2 * ir >= nw && nb != 0 {
                        c.moreio(nin, nout, 0, &mut scr.at(ibase + 1), &mut nb, &mut nw);
                        lr += ir;
                        ir = 0;
                        nw = ibase + nw;
                    }
                    if lr + ir >= c.n2h {
                        lbl = 280;
                        continue;
                    }
                    ernext = scr[ibase + 2 * ir + 1];
                    srnext = scr[ibase + 2 * ir + 2];
                    if (ernext - er).abs() > small * er {
                        lbl = 280;
                        continue;
                    }
                    if (srnext - sr).abs() < small * sr {
                        lbl = 260;
                        continue;
                    }
                    er = sigfig(er, 7, -1);
                    ernext = sigfig(ernext, 7, 1);
                    scr[ibase + 2 * ir + 1] = ernext;
                    lbl = 280;
                }
                280 => {
                    if eg.abs() > erl * (1.0 + small) {
                        lbl = 300;
                        continue;
                    }
                    ig += 1;
                    if ig > self.ngo {
                        lbl = 300;
                        continue;
                    }
                    c.finda(ig, &mut aa, 1, iold, &mut bufo, NBUF);
                    eg = aa[0];
                    lbl = 300;
                }
                300 => {
                    if en.abs() <= enl.abs() * (1.0 + small) {
                        lbl = 240;
                        continue;
                    }
                    i = 2;
                    x[2] = enl.abs();
                    y[2] = snl;
                    x[1] = en.abs();
                    y[1] = sn;
                    lbl = 310;
                }
                310 => {
                    if x[i - 1] / x[i] < stpmin {
                        lbl = 325;
                        continue;
                    }
                    if x[i] > elim {
                        lbl = 320;
                        continue;
                    }
                    let mut found = false;
                    for ipwr in -4..=5 {
                        xm = crate::fmath::powi(ten, ipwr);
                        if x[i - 1] > up * xm && x[i] < dn * xm {
                            found = true;
                            break;
                        }
                        xm = xm / 2.0;
                        xm = sigfig(xm, 1, 0);
                        if x[i - 1] > up * xm && x[i] < dn * xm {
                            found = true;
                            break;
                        }
                        xm = 4.0 * xm / 10.0;
                        xm = sigfig(xm, 1, 0);
                        if x[i - 1] > up * xm && x[i] < dn * xm {
                            found = true;
                            break;
                        }
                    }
                    if !found {
                        xm = therm;
                        if x[i - 1] > up * xm && x[i] < dn * xm {
                            found = true;
                        }
                    }
                    if found {
                        // 315
                        c.terp1(x[i], y[i], x[i - 1], y[i - 1], xm, &mut ym, intn);
                        lbl = 330;
                    } else {
                        lbl = 320;
                    }
                }
                320 => {
                    if intn > 2 {
                        // 322
                        xm = (x[i - 1] + x[i]) / 2.0;
                        xm = sigfig(xm, 7, 0);
                        c.terp1(x[i], y[i], x[i - 1], y[i - 1], xm, &mut ym, intn);
                        c.terp1(x[i], y[i], x[i - 1], y[i - 1], xm, &mut yl, 2);
                        let mut errn = self.err;
                        if x[i - 1] < trange {
                            errn = self.err / 5.0;
                        }
                        test = errn * ym + ssmall;
                        if (ym - yl).abs() < test {
                            lbl = 325;
                            continue;
                        }
                        lbl = 330;
                        continue;
                    }
                    if x[i] > elim {
                        lbl = 325;
                        continue;
                    }
                    if x[i - 1] / x[i] <= stpmax {
                        lbl = 325;
                        continue;
                    }
                    xm = (x[i - 1] + x[i]) / 2.0;
                    xm = sigfig(xm, 7, 0);
                    c.terp1(x[i], y[i], x[i - 1], y[i - 1], xm, &mut ym, 2);
                    lbl = 330;
                }
                325 => {
                    in_ += 1;
                    let mut ent = -x[i];
                    if (enl.abs() - x[i]).abs() < small * x[i] {
                        ent = enl;
                    }
                    aa[0] = ent;
                    c.loada(in_, &aa, 1, inew, &mut bufn, NBUF);
                    i -= 1;
                    if i > 1 {
                        lbl = 310;
                        continue;
                    }
                    enl = en;
                    snl = sn;
                    lbl = 240;
                }
                330 => {
                    i += 1;
                    if i > ndim {
                        c.error("lunion", "exceeded stack.", " ");
                    }
                    // 340
                    x[i] = x[i - 1];
                    y[i] = y[i - 1];
                    x[i - 1] = xm;
                    y[i - 1] = ym;
                    lbl = 310;
                }
                350 => {
                    in_ += 1;
                    let inx = -in_;
                    aa[0] = en;
                    c.loada(inx, &aa, 1, inew, &mut bufn, NBUF);
                    let ngn = in_;
                    self.ngo = ngn;
                    std::mem::swap(&mut iold, &mut inew);
                    iss -= 1;
                    if iss > 0 {
                        lbl = 181;
                        continue;
                    }
                    c.tosend(nin, 0, 0, &mut scr.all());
                    c.asend(nout, 0);
                    lbl = 110;
                }
                _ => {
                    // 410
                    *ngr = self.ngo;
                    let mut ngneg = 0;
                    let mut jg = 0;
                    let mut egl = 0.0f64;
                    for ig in 1..=self.ngo {
                        c.finda(ig, &mut aa, 1, iold, &mut bufo, NBUF);
                        let mut eg = aa[0];
                        let mut skip = (eg.abs() - egl).abs() <= small * egl;
                        if !skip {
                            if eg < zero {
                                ngneg += 1;
                            }
                            eg = eg.abs();
                            if !(ig == 1 || ig == self.ngo || eg > emax) {
                                if (eg - self.eresl).abs() <= small * self.eresl
                                    || (eg - self.eresr).abs() <= small * self.eresr
                                    || (eg - self.eresu).abs() <= small * self.eresu
                                    || (eg - self.eresm).abs() <= small * self.eresm
                                    || (eg - self.eresh).abs() <= small * self.eresh
                                {
                                    skip = true;
                                }
                            }
                        }
                        if skip {
                            *ngr -= 1;
                            continue;
                        }
                        // 415
                        egl = eg;
                        jg += 1;
                        let mut igt = jg;
                        if ig == self.ngo {
                            igt = -igt;
                        }
                        aa[0] = eg;
                        c.loada(igt, &aa, 1, ngrid, &mut bufn, NBUF);
                    }
                    self.ngo = *ngr;
                    let ngpos = self.ngo - ngneg;
                    let time = c.timer();
                    c.io.so(&format!(
                        "\n number of user and resonance nodes           ={}\n points in initial unionized grid             ={}\n points added by linearization                ={}{}{}s",
                        ffmt::i(self.nodes, 8),
                        ffmt::i(ngpos, 8),
                        ffmt::i(ngneg, 8),
                        " ".repeat(14),
                        ffmt::f(time, 8, 1)
                    ));
                    return;
                }
            }
        }
    }

    fn sv_elim(&self) -> f64 {
        self.sv.elim.unwrap_or(0.99e6)
    }

    /// Reconstruct the resonance cross sections on the union grid.
    /// Returns the new nrtot.
    fn resxs(&mut self, c: &mut Ctx, ngrid: i32, nout: i32, nscr: i32, nrtot_in: i32) -> i32 {
        let half = 0.5e0f64;
        let estp = 4.1e0f64;
        let small = 1.0e-10f64;
        let one = 1.0f64;
        let ten = 10.0f64;
        let trange = 0.4999e0f64;
        let zero = 0.0f64;
        let tenth = 0.1e0f64;
        let ndim: i32 = 30;
        let nscrl = 16;
        let nsig = self.nsig;
        let _ = nrtot_in;

        if self.lssf == 1 && self.eresl == self.eresu {
            return 0;
        }

        c.openz(-nscr, 1);
        c.openz(-nscrl, 1);
        let mut bufr = vec![0.0f64; NBUFR as usize];
        let mut bufg = vec![0.0f64; NBUFG as usize];
        let mut bufl = vec![0.0f64; NBUFL as usize];
        let mut x = RV::new(ndim as usize + 2);
        let mut y = RV::new(ndim as usize + 2);
        let mut sigs = RM::new((nsig - 1).max(1) as usize, ndim as usize + 2);
        let mut sig = RV::new(nsig as usize + 1);
        let mut sigl = RV::new(nsig as usize + 1);
        let mut res = RV::new(nsig as usize + 2);
        let mut sigd = RM::new(10, 10);
        let mut resl = RV::new(101);
        let mut dm = RV::new(10);
        let mut aa = [0.0f64; 1];
        let mut ngneg = 0;
        let mut nmax = 0;
        let mut cmax = 0.0f64;
        let mut cint = 0.0f64;
        let mut eint = 0.0f64;
        let mut emax = 0.0f64;
        let mut fint = 0.0f64;
        let mut fmax = 0.0f64;
        for i in 1..=nsig + 1 {
            res[i] = 0.0;
        }
        let mut elo = self.eresl;
        let mut iii = crate::fmath::int(log10(elo * (1.0 + small)) + small);
        if log10(elo * (1.0 + small)) < zero {
            iii -= 1;
        }
        let mut xxx = crate::fmath::powi(ten, iii);
        let mut ehi;
        if 2.0 * xxx > elo * (1.0 + small) && elo >= one * (1.0 - small) {
            ehi = 2.0 * xxx;
        } else if 5.0 * xxx > elo * (1.0 + small) && elo >= one * (1.0 - small) {
            ehi = 5.0 * xxx;
        } else {
            ehi = 10.0 * xxx;
        }
        if ehi > self.eresr {
            ehi = self.eresr;
        }
        let e = 0.0;
        self.csunr1(c, e, &mut sig);
        self.csunr2(c, e, &mut sig);
        c.io.so(&format!(
            "\n{b15}estimated maximum error due to\n{b15}resonance integral check (errmax,errint)\n\n    upper      elastic   percent   capture   percent   fission   percent\n    energy     integral   error    integral   error    integral   error\n{e}",
            b15 = " ".repeat(15),
            e = ffmt::e1p(elo, 10, 2)
        ));

        let mut in_ = 0;
        let mut ig = 0;
        let mut egl;
        loop {
            ig += 1;
            c.finda(ig, &mut aa, 1, ngrid, &mut bufg, NBUFG);
            egl = aa[0];
            if egl >= self.eresl {
                break;
            }
        }
        self.sigma(c, egl, &mut sigl, &mut sigd);

        let mut i = 0i32;
        // xm and dx are plain locals that keep their last values between
        // passes, exactly as in the Fortran (label 150 can be reached
        // without recomputing them).
        let mut xm = 0.0f64;
        let mut dx = 0.0f64;
        let mut eg = 0.0f64;
        let mut lbl = 105;
        loop {
            match lbl {
                105 => {
                    ig += 1;
                    c.finda(ig, &mut aa, 1, ngrid, &mut bufg, NBUFG);
                    eg = aa[0];
                    self.sigma(c, eg, &mut sig, &mut sigd);
                    if eg >= self.eresh {
                        lbl = 190;
                        continue;
                    }
                    i = 2;
                    x[2] = egl;
                    y[2] = sigl[1];
                    x[1] = eg;
                    y[1] = sig[1];
                    for j in 1..=nsig - 1 {
                        sigs[(j, 2)] = sigl[j + 1];
                        sigs[(j, 1)] = sig[j + 1];
                    }
                    if (x[2] < self.eresr && x[1] > self.eresr)
                        || (x[2] < self.eresu && x[1] > self.eresu)
                        || (x[2] < self.eresm && x[1] > self.eresm)
                    {
                        lbl = 150;
                        continue;
                    }
                    lbl = 115;
                }
                115 => {
                    xm = half * (x[i] + x[i - 1]);
                    dx = x[i - 1] - x[i];
                    let mut ndig = 9;
                    if xm > tenth && xm < one {
                        ndig = 8;
                    }
                    if !(xm > sigfig(x[i], ndig, 1) && xm < sigfig(x[i - 1], ndig, -1)) {
                        lbl = 150;
                        continue;
                    }
                    // 135
                    if xm > sigfig(x[i], 7, 1) && xm < sigfig(x[i - 1], 7, -1) {
                        xm = sigfig(xm, 7, 0);
                    } else {
                        xm = sigfig(xm, ndig, 0);
                    }
                    self.sigma(c, xm, &mut sig, &mut sigd);
                    let fr2 = (xm - x[i]) / dx;
                    let fr1 = 1.0 - fr2;
                    for j in 1..=nsig - 1 {
                        let sl = fr1 * sigs[(j, i)] + fr2 * sigs[(j, i - 1)];
                        dm[j] = (sig[j + 1] - sl).abs();
                    }
                    let mut errn = self.err;
                    if x[i - 1] < trange {
                        errn = self.err / 5.0;
                    }
                    let mut errm = self.errmax;
                    if x[i - 1] < trange {
                        errm = self.errmax / 5.0;
                    }
                    let mut no = 0;
                    for j in 1..=nsig - 1 {
                        if dm[j] > errn * sig[j + 1] {
                            no = 1;
                        }
                    }
                    if no == 0 {
                        // 145
                        let est = estp * (x[i] - res[1]);
                        if in_ > 3 && dx > est && xm < self.eresr {
                            lbl = 175;
                            continue;
                        }
                        lbl = 150;
                        continue;
                    }
                    let mut go175 = false;
                    for j in 1..=nsig - 1 {
                        if dm[j] > errm * sig[j + 1] {
                            go175 = true;
                            break;
                        }
                    }
                    if !go175 {
                        let tsti = 2.0 * self.errint * xm / dx;
                        for j in 1..=nsig - 1 {
                            if dm[j] >= tsti {
                                go175 = true;
                                break;
                            }
                        }
                    }
                    if go175 {
                        lbl = 175;
                        continue;
                    }
                    let dm1 = dm[1];
                    let dm2 = dm[2];
                    let dm3 = dm[3];
                    nmax += 1;
                    cmax = cmax + dm3 * dx / (2.0 * xm);
                    emax = emax + dm1 * dx / (2.0 * xm);
                    fmax = fmax + dm2 * dx / (2.0 * xm);
                    lbl = 150;
                }
                150 => {
                    in_ += 1;
                    res[1] = x[i];
                    if (egl - x[i]).abs() <= small * egl {
                        res[1] = egl;
                    } else {
                        ngneg += 1;
                    }
                    res[2] = y[i];
                    for j in 1..=nsig - 1 {
                        res[j + 2] = sigs[(j, i)];
                    }
                    c.loada(in_, &res.0, nsig + 1, nscr, &mut bufr, NBUFR);
                    if self.ncoef > 1 {
                        let r1 = res[1];
                        self.sigma(c, r1, &mut sig, &mut sigd);
                        for j in 1..=self.nmtres {
                            if self.mcards[j as usize] == 0 && sigd[(1, j)] > zero {
                                self.mcards[j as usize] = in_;
                            }
                        }
                        let mut k = 1;
                        resl[k] = res[1];
                        for l in 1..=self.ncoef {
                            for j in 1..=self.nmtres {
                                k += 1;
                                resl[k] = sigd[(l, j)];
                            }
                        }
                        c.loada(in_, &resl.0, 1 + self.ncoef * self.nmtres, nscrl, &mut bufl, NBUFL);
                    }
                    let mut c1 = sigs[(3, i)];
                    let mut c2 = sigs[(3, i - 1)];
                    cint = cint + (c1 + c2) * dx / (2.0 * xm);
                    c1 = sigs[(1, i)];
                    c2 = sigs[(1, i - 1)];
                    eint = eint + (c1 + c2) * dx / (2.0 * xm);
                    c1 = sigs[(2, i)];
                    c2 = sigs[(2, i - 1)];
                    fint = fint + (c1 + c2) * dx / (2.0 * xm);
                    if !(x[i - 1] < ehi * (1.0 - small) || ehi > self.eresh || elo >= ehi * (1.0 - small)) {
                        if cint != zero {
                            cmax = 100.0 * cmax / cint;
                        }
                        if eint != zero {
                            emax = 100.0 * emax / eint;
                        }
                        if fint != zero {
                            fmax = 100.0 * fmax / fint;
                        }
                        c.io.so(&format!(
                            "{} {}{}{}{}{}{}",
                            ffmt::e1p(x[i - 1], 10, 2),
                            ffmt::e1p(eint, 12, 2),
                            ffmt::fp(emax, 8, 3, 1),
                            ffmt::e1p(cint, 12, 2),
                            ffmt::fp(cmax, 8, 3, 1),
                            ffmt::e1p(fint, 12, 2),
                            ffmt::fp(fmax, 8, 3, 1)
                        ));
                        elo = x[i - 1];
                        iii = crate::fmath::int(log10(elo * (1.0 + small)) + small);
                        if log10(elo * (1.0 + small)) < 0.0 {
                            iii -= 1;
                        }
                        xxx = crate::fmath::powi(ten, iii);
                        if 2.0 * xxx > elo * (1.0 + small) && elo >= one * (1.0 - small) {
                            ehi = 2.0 * xxx;
                        } else if 5.0 * xxx > elo * (1.0 + small) && elo >= one * (1.0 - small) {
                            ehi = 5.0 * xxx;
                        } else {
                            ehi = 10.0 * xxx;
                        }
                        if ehi > self.eresr {
                            ehi = self.eresr;
                        }
                        if elo >= self.eresr {
                            ehi = self.eresh;
                        }
                        cint = 0.0;
                        eint = 0.0;
                        fint = 0.0;
                        cmax = 0.0;
                        emax = 0.0;
                        fmax = 0.0;
                    }
                    // 165
                    i -= 1;
                    if i > 1 {
                        lbl = 115;
                        continue;
                    }
                    egl = eg;
                    sigl[1] = y[1];
                    for j in 1..=nsig - 1 {
                        sigl[j + 1] = sigs[(j, 1)];
                    }
                    lbl = 105;
                }
                175 => {
                    i += 1;
                    if i > ndim {
                        c.error("resxs", "stack exceeded.", " ");
                    }
                    x[i] = x[i - 1];
                    y[i] = y[i - 1];
                    x[i - 1] = xm;
                    y[i - 1] = sig[1];
                    for j in 1..=nsig - 1 {
                        sigs[(j, i)] = sigs[(j, i - 1)];
                        sigs[(j, i - 1)] = sig[j + 1];
                    }
                    lbl = 115;
                }
                _ => {
                    // 190
                    in_ += 1;
                    let mut inf = -in_;
                    res[1] = egl;
                    res[2] = sigl[1];
                    for j in 1..=nsig - 1 {
                        res[j + 2] = sigs[(j, 1)];
                    }
                    let mut nrtot = in_;
                    c.loada(inf, &res.0, nsig + 1, nscr, &mut bufr, NBUFR);
                    if self.ncoef > 1 {
                        let r1 = res[1];
                        self.sigma(c, r1, &mut sig, &mut sigd);
                        let mut k = 1;
                        resl[k] = res[1];
                        for l in 1..=self.ncoef {
                            for j in 1..=self.nmtres {
                                k += 1;
                                resl[k] = sigd[(l, j)];
                            }
                        }
                        c.loada(inf, &resl.0, 1 + self.ncoef * self.nmtres, nscrl, &mut bufl, NBUFL);
                        for ii in 1..=self.nmtres as usize {
                            if self.mcards[ii] > 0 {
                                self.mcards[ii] = 3 + 2 * (nrtot - self.mcards[ii] + 1);
                            }
                        }
                    }

                    let mut ip = 1;
                    let mut xv = vec![0.0f64; (nsig + 1) as usize];
                    c.finda(ip, &mut xv, nsig + 1, nscr, &mut bufr, NBUFR);
                    in_ = 1;
                    c.loada(in_, &xv, nsig + 1, nout, &mut bufg, NBUFG);
                    while ip < nrtot {
                        ip += 1;
                        c.finda(ip, &mut xv, nsig + 1, nscr, &mut bufr, NBUFR);
                        in_ += 1;
                        inf = in_;
                        if ip == nrtot {
                            inf = -in_;
                        }
                        c.loada(inf, &xv, nsig + 1, nout, &mut bufg, NBUFG);
                    }
                    nrtot = in_;
                    let _time = c.timer();
                    c.io.so(&format!(
                        "\n points added by resonance reconstruction     ={}\n points affected by resonance integral check  ={}\n final number of resonance points             ={}",
                        ffmt::i(ngneg, 8),
                        ffmt::i(nmax, 8),
                        ffmt::i(nrtot, 8)
                    ));
                    c.closz(-nscr);
                    return nrtot;
                }
            }
        }
    }

    /// Cross sections at energy e: sum over sections of the resonance
    /// contributions plus the unresolved table.
    fn sigma(&mut self, c: &mut Ctx, e: f64, sig: &mut RV, sigd: &mut RM) {
        let zero = 0.0f64;
        let nsig = self.nsig;
        for i in 1..=nsig {
            sig[i] = 0.0;
        }
        let ncoef = self.ncoef.max(1);
        let nmt = self.nmtres.max(1);
        let mut siga = RM::new(ncoef as usize, nmt as usize);
        if self.ncoef > 1 {
            for i in 1..=self.ncoef {
                for j in 1..=self.nmtres {
                    siga[(i, j)] = 0.0;
                }
            }
        }
        let mut sigp = RV::new(nsig.max(5) as usize + 1);
        let mut sigx = RM::new(1, 1);
        for i in 1..=self.nsect {
            self.isect = i;
            let abn = self.abnt[i as usize];
            if e >= self.elt[i as usize] && e < self.eht[i as usize] {
                let mode = self.modet[i as usize];
                if mode != 11 && mode != 12 {
                    if mode == 0 {
                        self.csnorp(e, &mut sigp);
                    } else if mode == 1 {
                        self.csslbw(c, e, &mut sigp);
                    } else if mode == 2 && self.tempr == zero {
                        self.csmlbw(c, e, &mut sigp);
                    } else if mode == 2 && self.tempr > zero {
                        self.csmlbw2(c, e, &mut sigp);
                    } else if mode == 3 && self.nmtres == 0 {
                        self.csrmat(c, e, &mut sigp);
                    } else if mode == 3 && self.nmtres > 0 {
                        self.samm.cssammy(c, e, &mut sigp, &mut siga, &mut sigx, &self.mmtres, self.nmtres, self.ncoef, self.nresp, i, &self.res);
                    } else if mode == 4 {
                        self.csaa(c, e, &mut sigp);
                    } else if mode == 5 {
                        c.error("sigma", "general r-matrix not installed.", " ");
                    } else if mode == 6 {
                        self.cshybr(c, e, &mut sigp);
                    } else if mode == 7 {
                        self.samm.cssammy(c, e, &mut sigp, &mut siga, &mut sigx, &self.mmtres, self.nmtres, self.ncoef, self.nresp, i, &self.res);
                    } else {
                        c.error("sigma", "illegal option.", " ");
                    }
                    for j in 1..=nsig {
                        if sigp[j] < zero {
                            sigp[j] = 0.0;
                        }
                        sig[j] = sig[j] + sigp[j] * abn;
                    }
                }
            }
        }
        if self.ncoef > 1 {
            for i in 1..=self.ncoef {
                for j in 1..=self.nmtres {
                    sigd[(i, j)] = siga[(i, j)];
                }
            }
        }
        if e < self.eresu {
            return;
        }
        if e >= self.eresh {
            return;
        }
        let mut sigu = [0.0f64; 5];
        self.sigunr(e, &mut sigu, c.thr6);
        for j in 1..=4 {
            sig[j] = sig[j] + sigu[j as usize];
        }
    }
}

// ======================================================================
// resonance cross section formulas
// ======================================================================

fn cwaven() -> f64 {
    sqrt(2.0 * AMASSN * AMU * EV) * 1.0e-12 / HBAR
}

#[inline]
fn sq(x: f64) -> f64 {
    x * x
}

/// Complex arithmetic as gfortran emits it (Fortran rules: naive multiply,
/// Smith's algorithm for division, libm cexp/cabs).
#[derive(Clone, Copy, Debug, Default)]
#[repr(C)]
pub struct Cx {
    pub re: f64,
    pub im: f64,
}

extern "C" {
    fn cexp(z: Cx) -> Cx;
    fn cabs(z: Cx) -> f64;
}

impl Cx {
    pub fn new(re: f64, im: f64) -> Cx {
        Cx { re, im }
    }
    pub fn r(x: f64) -> Cx {
        Cx { re: x, im: 0.0 }
    }
    pub fn add(self, b: Cx) -> Cx {
        Cx { re: self.re + b.re, im: self.im + b.im }
    }
    pub fn sub(self, b: Cx) -> Cx {
        Cx { re: self.re - b.re, im: self.im - b.im }
    }
    pub fn mul(self, b: Cx) -> Cx {
        let t1 = self.re * b.re;
        let t2 = self.im * b.im;
        let t3 = self.re * b.im;
        let t4 = self.im * b.re;
        Cx { re: t1 - t2, im: t3 + t4 }
    }
    pub fn div(self, b: Cx) -> Cx {
        let (a, bb, c, d) = (self.re, self.im, b.re, b.im);
        if c.abs() < d.abs() {
            let ratio = c / d;
            let div = (c * ratio) + d;
            let tr = (a * ratio) + bb;
            let ti = (bb * ratio) - a;
            Cx { re: tr / div, im: ti / div }
        } else {
            let ratio = d / c;
            let div = (d * ratio) + c;
            let tr = (bb * ratio) + a;
            let ti = bb - (a * ratio);
            Cx { re: tr / div, im: ti / div }
        }
    }
    pub fn exp(self) -> Cx {
        unsafe { cexp(self) }
    }
    pub fn abs(self) -> f64 {
        unsafe { cabs(self) }
    }
}

impl Reconr {
    fn csnorp(&self, _e: f64, sigp: &mut RV) {
        for i in 1..=self.nsig {
            sigp[i] = 0.0;
        }
        let inow = self.ibaset[self.isect as usize];
        let ap = self.res[inow + 7];
        sigp[2] = 4.0 * crate::util::physics::PI * ap * ap;
    }

    fn csslbw(&mut self, c: &mut Ctx, e: f64, sigp: &mut RV) {
        let pi = crate::util::physics::PI;
        let bk = crate::util::physics::BK;
        // facts sets nothing for l > 4: the Fortran then keeps the values of
        // the previous call, so these live for the whole routine.
        let (mut se, mut pe) = (0.0f64, 0.0f64);
        let rc1 = 0.123e0f64;
        let rc2 = 0.08e0f64;
        let third = 0.333333333e0f64;
        let zero = 0.0f64;
        let rpi = sqrt(pi);
        let cwaven = cwaven();
        for i in 1..=self.nsig {
            sigp[i] = 0.0;
        }
        let tbk = self.tempr * bk;
        let mut ki = 1;
        let mut inow = self.ibaset[self.isect as usize];
        let naps = nint(self.res[inow + 5]);
        let awri = self.res[inow + 12];
        let ap = self.res[inow + 7];
        let aw = AMASSN * awri;
        let mut ra = rc1 * pow(aw, third) + rc2;
        if naps == 1 {
            ra = ap;
        }
        let spi = self.res[inow + 6];
        let spifac = 1.0 / (2.0 * spi + 1.0);
        let nls = nint(self.res[inow + 10]);
        inow += 12;
        let arat = awri / (awri + 1.0);
        let k = cwaven * arat * sqrt(e.abs());
        let pifac = pi / (k * k);
        let rho = k * ra;
        let rhoc = k * ap;
        let _ = c;

        for _l in 1..=nls {
            let nrs = nint(self.res[inow + 5]);
            let ll = nint(self.res[inow + 2]);
            let qx = self.res[inow + 1];
            self.lrx = nint(self.res[inow + 3]);
            let mut in_ = inow + nrs * 6 + 6;
            facts(ll, rho, &mut se, &mut pe);
            let mut pec = 0.0;
            if self.lrx != 0 {
                let rhop = cwaven * arat * sqrt((e + qx / arat).abs()) * ra;
                let mut lp = ll;
                if ll == 0 {
                    lp = 2;
                }
                if ll == 2 {
                    lp = 0;
                }
                let mut sec = 0.0;
                facts(lp, rhop, &mut sec, &mut pec);
                if e + qx / arat < zero {
                    pec = 0.0;
                }
            }
            let phi = facphi(ll, rhoc);
            let cos2p = crate::fmath::cos(2.0 * phi);
            let sin2p = crate::fmath::sin(2.0 * phi);
            let sinsq = sq(crate::fmath::sin(phi));
            let mut spot = (4 * (2 * ll + 1)) as f64 * pifac * sinsq;

            if self.tempr > zero {
                for _i in 1..=nrs {
                    inow += 6;
                    let er = self.res[inow];
                    let aj = self.res[inow + 1];
                    let gn = self.res[inow + 3];
                    let gg = self.res[inow + 4];
                    let gf = self.res[inow + 5];
                    let ser = self.res[in_];
                    let per = self.res[in_ + 1];
                    let rper = 1.0 / per;
                    let gc = self.res[in_ + 2];
                    let pex = self.res[in_ + 3];
                    in_ += 6;
                    let gx = gg + gf;
                    let gj = (2.0 * aj + 1.0) * spifac / 2.0;
                    let erp = er + gn * (ser - se) * rper / 2.0;
                    let _edelt = e - erp;
                    let gne = gn * pe * rper;
                    let mut gtt = gne + gx;
                    if gc != zero {
                        gtt = gtt + gc * pec / pex;
                    }
                    let ex = 2.0 * (e - erp) / gtt;
                    let delta = sqrt(4.0 * tbk * e / awri);
                    let theta = gtt / delta;
                    let ax = theta * ex / 2.0;
                    let y = theta / 2.0;
                    let (mut rew, mut aimw) = (0.0, 0.0);
                    quickw(ax, y, &mut rew, &mut aimw, &mut ki, self.tr.as_ref().unwrap(), self.ti.as_ref().unwrap());
                    let psi = rpi * theta * rew / 2.0;
                    let chi = rpi * theta * aimw / 2.0;
                    let smax = 4.0 * pifac * gj * gne / sq(gtt);
                    sigp[2] = sigp[2] + smax * ((cos2p * gtt - gx) * psi + sin2p * chi * gtt);
                    sigp[3] = sigp[3] + smax * gf * psi;
                    sigp[4] = sigp[4] + smax * gg * psi;
                    if self.lrx != 0 && gc != zero {
                        sigp[5] = sigp[5] + smax * gc * pec * psi / pex;
                    }
                }
                inow = in_;
            } else {
                for _i in 1..=nrs {
                    inow += 6;
                    let er = self.res[inow];
                    let aj = self.res[inow + 1];
                    let gn = self.res[inow + 3];
                    let gg = self.res[inow + 4];
                    let gf = self.res[inow + 5];
                    let ser = self.res[in_];
                    let per = self.res[in_ + 1];
                    let rper = 1.0 / per;
                    let gc = self.res[in_ + 2];
                    let pex = self.res[in_ + 3];
                    in_ += 6;
                    let gx = gg + gf;
                    let gj = (2.0 * aj + 1.0) * spifac / 2.0;
                    let erp = er + gn * (ser - se) * rper / 2.0;
                    let edelt = e - erp;
                    let gne = gn * pe * rper;
                    let mut gtt = gne + gx;
                    if gc != zero {
                        gtt = gtt + gc * pec / pex;
                    }
                    let comfac = pifac * gj * gne / (sq(edelt) + gtt * gtt / 4.0);
                    let add = comfac * (gne * cos2p - 2.0 * gx * sinsq + 2.0 * edelt * sin2p);
                    sigp[2] = sigp[2] + add;
                    sigp[3] = sigp[3] + comfac * gf;
                    sigp[4] = sigp[4] + comfac * gg;
                    if self.lrx != 0 && gc != zero {
                        sigp[5] = sigp[5] + comfac * gc * pec / pex;
                    }
                }
                inow = in_;
            }
            sigp[2] = sigfig(sigp[2], 8, 0);
            spot = sigfig(spot, 8, 0);
            sigp[2] = sigp[2] + spot;
        }
        sigp[1] = sigp[2] + sigp[3] + sigp[4];
        if self.lrx != 0 {
            sigp[1] = sigp[1] + sigp[5];
        }
    }

    fn csmlbw(&mut self, c: &mut Ctx, e: f64, sigp: &mut RV) {
        let pi = crate::util::physics::PI;
        let rc1 = 0.123e0f64;
        let rc2 = 0.08e0f64;
        // facts sets nothing for l > 4: the Fortran then keeps the values of
        // the previous call, so these live for the whole routine.
        let (mut se, mut pe) = (0.0f64, 0.0f64);
        let third = 0.333333333e0f64;
        let half = 0.5e0f64;
        let zero = 0.0f64;
        let cwaven = cwaven();
        if self.tempr > zero {
            c.error("csmlbw", "not coded for temperature gt 0 deg k.", " ");
        }
        for i in 1..=self.nsig {
            sigp[i] = 0.0;
        }
        let mut inow = self.ibaset[self.isect as usize];
        let nro = nint(self.res[inow + 4]);
        let naps = nint(self.res[inow + 5]);
        let mut ape = 0.0;
        if nro != 0 {
            let iro = inow + 6;
            inow = inow + 6 + 2 * nint(self.res[iro + 4]) + 2 * nint(self.res[iro + 5]);
            let mut ip = 2;
            let mut ir = 1;
            let (mut enx, mut idx) = (0.0, 0);
            let sub = crate::farr::F(&mut self.res.0[(iro - 1) as usize..]);
            c.terpa(&mut ape, e, &mut enx, &mut idx, &sub, &mut ip, &mut ir);
        }
        let awri = self.res[inow + 12];
        let mut ap = self.res[inow + 7];
        let aw = AMASSN * awri;
        let mut ra = rc1 * pow(aw, third) + rc2;
        if nro == 0 {
            if naps == 1 {
                ra = ap;
            }
        } else if naps == 0 {
            ap = ape;
        } else if naps == 1 {
            ap = ape;
            ra = ape;
        } else {
            ra = ap;
            ap = ape;
        }
        let spi = self.res[inow + 6];
        let den = 4.0 * spi + 2.0;
        let nls = nint(self.res[inow + 10]);
        let arat = awri / (awri + 1.0);
        let k = cwaven * arat * sqrt(e.abs());
        let pifac = pi / (k * k);
        let rho = k * ra;
        let rhoc = k * ap;
        inow += 12;
        let mut sigj = RM::new(10, 2);
        let mut gj = RV::new(10);

        for _l in 1..=nls {
            let nrs = nint(self.res[inow + 5]);
            let ll = nint(self.res[inow + 2]);
            let qx = self.res[inow + 1];
            self.lrx = nint(self.res[inow + 3]);
            facts(ll, rho, &mut se, &mut pe);
            let mut pec = 0.0;
            if self.lrx != 0 {
                let rhop = cwaven * arat * sqrt((e + qx / arat).abs()) * ra;
                let mut lp = ll;
                if ll == 0 {
                    lp = 2;
                }
                if ll == 2 {
                    lp = 0;
                }
                let mut sec = 0.0;
                facts(lp, rhop, &mut sec, &mut pec);
                if e + qx / arat < zero {
                    pec = 0.0;
                }
            }
            let phi = facphi(ll, rhoc);
            let cos2p = 1.0 - crate::fmath::cos(2.0 * phi);
            let sin2p = crate::fmath::sin(2.0 * phi);
            let mut sum = 0.0;
            let fl = ll as f64;
            let ajmin = ((spi - fl).abs() - half).abs();
            let ajmax = spi + fl + half;
            let nj = nint(ajmax - ajmin + 1.0);
            let mut aj = ajmin;
            for i in 1..=nj {
                gj[i] = (2.0 * aj + 1.0) / den;
                aj = aj + 1.0;
                sum = sum + gj[i];
            }
            let diff = 2.0 * fl + 1.0 - sum;
            for ii in 1..=2 {
                for i in 1..=nj {
                    sigj[(i, ii)] = 0.0;
                }
            }
            inow += 6;
            let mut in_ = inow + nrs * 6;
            for _i in 1..=nrs {
                let er = self.res[inow];
                let j = nint(self.res[inow + 1] - ajmin + 1.0);
                let gn = self.res[inow + 3];
                let gg = self.res[inow + 4];
                let gf = self.res[inow + 5];
                let ser = self.res[in_];
                let per = self.res[in_ + 1];
                let rper = 1.0 / per;
                let gc = self.res[in_ + 2];
                let pex = self.res[in_ + 3];
                in_ += 6;
                inow += 6;
                let erp = er + gn * (ser - se) * rper / 2.0;
                let edelt = e - erp;
                let gne = gn * pe * rper;
                let gx = gg + gf;
                let mut gtt = gne + gx;
                if gc != zero {
                    gtt = gtt + gc * pec / pex;
                }
                let x = 2.0 * edelt / gtt;
                let mut comfac = 2.0 * gne / gtt / (1.0 + x * x);
                sigj[(j, 1)] = sigj[(j, 1)] + comfac;
                sigj[(j, 2)] = sigj[(j, 2)] + comfac * x;
                comfac = comfac * gj[j] / gtt;
                sigp[3] = sigp[3] + comfac * gf;
                sigp[4] = sigp[4] + comfac * gg;
                if self.lrx != 0 && gc != zero {
                    sigp[5] = sigp[5] + comfac * gc * pec / pex;
                }
            }
            for j in 1..=nj {
                let add = gj[j] * (sq(cos2p - sigj[(j, 1)]) + sq(sin2p + sigj[(j, 2)]));
                sigp[2] = sigp[2] + add;
            }
            sigp[2] = sigp[2] + 2.0 * diff * cos2p;
            inow = in_;
        }
        sigp[2] = sigp[2] * pifac;
        sigp[3] = sigp[3] * 2.0 * pifac;
        sigp[4] = sigp[4] * 2.0 * pifac;
        sigp[1] = sigp[2] + sigp[3] + sigp[4];
        if self.lrx != 0 {
            sigp[1] = sigp[1] + sigp[5];
        }
    }

    fn csmlbw2(&mut self, _c: &mut Ctx, e: f64, sigp: &mut RV) {
        let pi = crate::util::physics::PI;
        let bk = crate::util::physics::BK;
        let rc1 = 0.123e0f64;
        let rc2 = 0.08e0f64;
        // facts sets nothing for l > 4: the Fortran then keeps the values of
        // the previous call, so these live for the whole routine.
        let (mut se, mut pe) = (0.0f64, 0.0f64);
        let third = 0.333333333e0f64;
        let tmin = 10e0f64;
        let zero = 0.0f64;
        let rpi = sqrt(pi);
        let cwaven = cwaven();
        for i in 1..=self.nsig {
            sigp[i] = 0.0;
        }
        let tbk = self.tempr * bk;
        let mut ki = 1;
        let mut inow = self.ibaset[self.isect as usize];
        let naps = nint(self.res[inow + 5]);
        let awri = self.res[inow + 12];
        let ap = self.res[inow + 7];
        let aw = AMASSN * awri;
        let mut ra = rc1 * pow(aw, third) + rc2;
        if naps == 1 {
            ra = ap;
        }
        let spi = self.res[inow + 6];
        let spifac = 1.0 / (2.0 * spi + 1.0);
        let nls = nint(self.res[inow + 10]);
        inow += 12;
        let arat = awri / (awri + 1.0);
        let k = cwaven * arat * sqrt(e.abs());
        let pifac = pi / (k * k);
        let rho = k * ra;
        let rhoc = k * ap;

        let mut igh = 0;
        if (e - 1.0e-5).abs() < 1.0e-7 {
            self.sv.enow_mlbw2 = e;
            igh = 1;
        }
        let enow = self.sv.enow_mlbw2;
        let mut test = enow + enow / 20.0;
        if enow + 10.0 < test {
            test = enow + 10.0;
        }
        if e >= test {
            self.sv.enow_mlbw2 = e;
            igh = 1;
        }

        for _l in 1..=nls {
            let nrs = nint(self.res[inow + 5]);
            let ll = nint(self.res[inow + 2]);
            let qx = self.res[inow + 1];
            self.lrx = nint(self.res[inow + 3]);
            let in_ = inow + nrs * 6 + 6;
            facts(ll, rho, &mut se, &mut pe);
            let mut pec = 0.0;
            if self.lrx != 0 {
                let rhop = cwaven * arat * sqrt((e + qx).abs()) * ra;
                let mut sec = 0.0;
                facts(ll, rhop, &mut sec, &mut pec);
                if e + qx / arat < zero {
                    pec = 0.0;
                }
            }
            let phi = facphi(ll, rhoc);
            let cos2p = crate::fmath::cos(2.0 * phi);
            let sin2p = crate::fmath::sin(2.0 * phi);
            let sinsq = sq(crate::fmath::sin(phi));
            let spot = (4 * (2 * ll + 1)) as f64 * pifac * sinsq;
            let mut sml = 0.0;
            inow += 6;
            for i in 1..=nrs {
                let er = self.res[inow + 6 * (i - 1)];
                let aj = self.res[inow + 6 * (i - 1) + 1];
                let _gt = self.res[inow + 6 * (i - 1) + 2];
                let gn = self.res[inow + 6 * (i - 1) + 3];
                let gg = self.res[inow + 6 * (i - 1) + 4];
                let gf = self.res[inow + 6 * (i - 1) + 5];
                let ser = self.res[in_ + 6 * (i - 1)];
                let per = self.res[in_ + 6 * (i - 1) + 1];
                let rper = 1.0 / per;
                let gc = self.res[in_ + 6 * (i - 1) + 2];
                let pex = self.res[in_ + 6 * (i - 1) + 3];
                let rpex = 1.0 / pex;
                let gx = gg + gf;
                let gj = (2.0 * aj + 1.0) * spifac / 2.0;
                let erp = er + gn * (ser - se) * rper / 2.0;
                let gne = gn * pe * rper;
                let mut gtt = gne + gx;
                if gc != zero {
                    gtt = gtt + gc * pec * rpex;
                }
                if igh == 1 {
                    let mut g = 0.0;
                    let mut h = 0.0;
                    for j in 1..=nrs {
                        let ajj = self.res[inow + 6 * (j - 1) + 1];
                        if j != i && ajj == aj && (i - j).abs() < 200 {
                            let mut erj = self.res[inow + 6 * (j - 1)];
                            let sej = self.res[in_ + 6 * (j - 1)];
                            let pej = self.res[in_ + 6 * (j - 1) + 1];
                            let gcj = self.res[in_ + 6 * (j - 1) + 2];
                            let pcj = self.res[in_ + 6 * (j - 1) + 3];
                            let mut gnj = self.res[inow + 6 * (j - 1) + 3];
                            let ggj = self.res[inow + 6 * (j - 1) + 4];
                            let gfj = self.res[inow + 6 * (j - 1) + 5];
                            erj = erj + (sej - se) * gnj / (2.0 * pej);
                            gnj = gnj * pe / pej;
                            let mut gtj = gnj + ggj + gfj;
                            if gcj != zero {
                                gtj = gtj + gcj * pec / pcj;
                            }
                            g = g + gne * gnj * (gtt + gtj) / (sq(erp - erj) + sq(gtt + gtj) / 4.0);
                            h = h + gne * gnj * (erp - erj) / (sq(erp - erj) + sq(gtt + gtj) / 4.0);
                        }
                    }
                    g = g / 2.0;
                    self.res[in_ + 6 * (i - 1) + 4] = g / e;
                    self.res[in_ + 6 * (i - 1) + 5] = h / e;
                }
                let ex = 2.0 * (e - erp) / gtt;
                let psi;
                let chi;
                if self.tempr > tmin {
                    let delta = sqrt(4.0 * tbk * e / awri);
                    let theta = gtt / delta;
                    let ax = theta * ex / 2.0;
                    let y = theta / 2.0;
                    let (mut rew, mut aimw) = (0.0, 0.0);
                    quickw(ax, y, &mut rew, &mut aimw, &mut ki, self.tr.as_ref().unwrap(), self.ti.as_ref().unwrap());
                    psi = rpi * theta * rew / 2.0;
                    chi = rpi * theta * aimw / 2.0;
                } else {
                    psi = 1.0 / (1.0 + sq(ex));
                    chi = ex / (1.0 + sq(ex));
                }
                let g = self.res[in_ + 6 * (i - 1) + 4];
                let h = self.res[in_ + 6 * (i - 1) + 5];
                let smax = 4.0 * pifac * gj * gne / sq(gtt);
                sigp[2] = sigp[2] + smax * ((cos2p * gtt - gx) * psi + sin2p * gtt * chi);
                sml = sml + smax * e * (g * psi + h * chi) * gtt / gne;
                sigp[3] = sigp[3] + smax * gf * psi;
                sigp[4] = sigp[4] + smax * gg * psi;
            }
            inow += 12 * nrs;
            let ssl = sigp[2] + spot;
            sigp[2] = ssl + sml;
        }
        sigp[1] = sigp[2] + sigp[3] + sigp[4];
    }

    fn csrmat(&mut self, c: &mut Ctx, e: f64, sigp: &mut RV) {
        let pi = crate::util::physics::PI;
        let rc1 = 0.123e0f64;
        let rc2 = 0.08e0f64;
        let third = 0.333333333e0f64;
        let quar = 0.25e0f64;
        // facts sets nothing for l > 4: the Fortran then keeps the values of
        // the previous call, so these live for the whole routine.
        let (mut se, mut dum1) = (0.0f64, 0.0f64);
        let haf = 0.5e0f64;
        let uno = 1.0e0f64;
        let two = 2.0e0f64;
        let four = 4.0e0f64;
        let small = 3.0e-4f64;
        let zero = 0.0f64;
        let cwaven = cwaven();
        if self.tempr > zero {
            c.error("csrmat", "not coded for temperatures gt 0 deg k.", " ");
        }
        for i in 1..=self.nsig {
            sigp[i] = 0.0;
        }
        let mut inow = self.ibaset[self.isect as usize];
        let naps = nint(self.res[inow + 5]);
        let awri = self.res[inow + 12];
        let ap = self.res[inow + 7];
        let aw = AMASSN * awri;
        let mut ra = rc1 * pow(aw, third) + rc2;
        if naps == 1 {
            ra = ap;
        }
        let spi = self.res[inow + 6];
        let gjd = 2.0 * (2.0 * spi + 1.0);
        let nls = nint(self.res[inow + 10]);
        let arat = awri / (awri + 1.0);
        let k = cwaven * arat * sqrt(e.abs());
        let pifac = pi / (k * k);
        let mut rho;
        let mut rhoc;
        let mut gfa;
        let mut gfb;
        let mut gf = 0.0f64;
        inow += 12;
        let mut r = [[0.0f64; 4]; 4];
        let mut s = [[0.0f64; 4]; 4];
        let mut ri = std::mem::take(&mut self.sv.csr_ri);
        let mut si = std::mem::take(&mut self.sv.csr_si);
        let mut in_ = inow;

        for _l in 1..=nls {
            let inowb = inow;
            let ll = nint(self.res[inow + 2]);
            let nrs = nint(self.res[inow + 5]);
            if nrs == 0 {
                let strng = format!("nrs=0 for Reich-Moore and l={}", ffmt::i(ll, 2));
                c.mess("csrmat", &strng, "malformed ENDF file, check evaluation");
                inow += 6;
                continue;
            }
            let ncyc = nint(self.res[inow + 4]) / nrs;
            let apl = self.res[inow + 1];
            rhoc = k * ap;
            rho = k * ra;
            if apl != zero {
                rhoc = k * apl;
            }
            if apl != zero && naps == 1 {
                rho = k * apl;
            }
            facts(ll, rho, &mut se, &mut dum1);
            let pe = dum1;
            let phi = facphi(ll, rhoc);
            let phid = phi;
            let p1 = crate::fmath::cos(2.0 * phid);
            let p2 = crate::fmath::sin(2.0 * phid);
            let fl = ll as f64;
            let ajmin = ((spi - fl).abs() - haf).abs();
            let ajmax = spi + fl + haf;
            let numj = nint(ajmax - ajmin + 1.0);
            let mut ajc = ajmin - 1.0;
            let jjl = if ll != 0 && (fl > spi - haf && fl <= spi) { 0 } else { 1 };
            for jj in 1..=numj {
                inow = inowb;
                ajc = ajc + 1.0;
                let gj = (2.0 * ajc + 1.0) / gjd;
                let mut kchanl = 0;
                let idone = 0;
                while kchanl < 2 && idone == 0 {
                    kchanl += 1;
                    inow = inowb;
                    let mut kpstv = 0;
                    let mut kngtv = 0;
                    for j in 1..=3 {
                        for i in 1..=3 {
                            s[j][i] = 0.0;
                            r[j][i] = 0.0;
                        }
                    }
                    inow += 6;
                    in_ = inow + nrs * 6;
                    for _i in 1..=nrs {
                        let aj = self.res[inow + 1].abs();
                        if (aj - ajc).abs() <= quar {
                            if self.res[inow + 1] < zero {
                                kngtv += 1;
                            }
                            if self.res[inow + 1] > zero {
                                kpstv += 1;
                            }
                            let mut iskip = 0;
                            if kchanl == 1 && self.res[inow + 1] < zero {
                                iskip = 1;
                            }
                            if kchanl == 2 && self.res[inow + 1] > zero {
                                iskip = 1;
                            }
                            if iskip == 0 {
                                let er = self.res[inow];
                                let gn = self.res[inow + 2];
                                let gg = self.res[inow + 3];
                                gfa = self.res[inow + 4];
                                gfb = self.res[inow + 5];
                                let per = self.res[in_ + 1];
                                let a1 = sqrt(gn * pe / per);
                                let mut a2 = 0.0;
                                if gfa != zero {
                                    a2 = sqrt(gfa.abs());
                                }
                                if gfa < zero {
                                    a2 = -a2;
                                }
                                let mut a3 = 0.0;
                                if gfb != zero {
                                    a3 = sqrt(gfb.abs());
                                }
                                if gfb < zero {
                                    a3 = -a3;
                                }
                                let diff = er - e;
                                let den = diff * diff + quar * gg * gg;
                                if den == 0.0 {
                                    c.io.se(&format!(
                                        "er,e,diff,gg={}{}{}{}",
                                        ffmt::e1p(er, 11, 4),
                                        ffmt::e1p(e, 11, 4),
                                        ffmt::e1p(diff, 11, 4),
                                        ffmt::e1p(gg, 11, 4)
                                    ));
                                }
                                let de2 = haf * diff / den;
                                let gg4 = quar * gg / den;
                                r[1][1] = r[1][1] + gg4 * a1 * a1;
                                s[1][1] = s[1][1] - de2 * a1 * a1;
                                if gfa != zero || gfb != zero {
                                    r[1][2] = r[1][2] + gg4 * a1 * a2;
                                    s[1][2] = s[1][2] - de2 * a1 * a2;
                                    r[1][3] = r[1][3] + gg4 * a1 * a3;
                                    s[1][3] = s[1][3] - de2 * a1 * a3;
                                    r[2][2] = r[2][2] + gg4 * a2 * a2;
                                    s[2][2] = s[2][2] - de2 * a2 * a2;
                                    r[3][3] = r[3][3] + gg4 * a3 * a3;
                                    s[3][3] = s[3][3] - de2 * a3 * a3;
                                    r[2][3] = r[2][3] + gg4 * a2 * a3;
                                    s[2][3] = s[2][3] - de2 * a2 * a3;
                                    gf = 1.0;
                                }
                            }
                        }
                        inow += ncyc;
                        in_ += 6;
                    }

                    let mut kk = 0;
                    if kchanl == 1 {
                        if kpstv > 0 {
                            if kngtv == 0 {
                                kk = if jj > jjl && jj < numj { 2 } else { 1 };
                            } else if kngtv > 0 {
                                kk = 1;
                            }
                        } else if kpstv == 0 {
                            if kngtv == 0 {
                                kk = if jj > jjl && jj < numj { 2 } else { 1 };
                            } else if kngtv > 0 {
                                kk = 0;
                            }
                        }
                    } else if kchanl == 2 {
                        if kpstv > 0 {
                            if kngtv > 0 {
                                kk = 1;
                            }
                        } else if kpstv == 0 && kngtv > 0 {
                            kk = if jj > jjl && jj < numj { 2 } else { 1 };
                        }
                    }
                    if kk != 0 {
                        let mut termf;
                        let mut termt;
                        let mut termn;
                        if gf != zero {
                            r[1][1] = uno + r[1][1];
                            r[2][2] = uno + r[2][2];
                            r[3][3] = uno + r[3][3];
                            r[2][1] = r[1][2];
                            s[2][1] = s[1][2];
                            r[3][1] = r[1][3];
                            s[3][1] = s[1][3];
                            r[3][2] = r[2][3];
                            s[3][2] = s[2][3];
                            frobns(&mut r, &s, &mut ri, &mut si);
                            let t1 = ri[1][2];
                            let t2 = si[1][2];
                            let t3 = ri[1][3];
                            let t4 = si[1][3];
                            termf = four * gj * (t1 * t1 + t2 * t2 + t3 * t3 + t4 * t4);
                            let u11r = p1 * (two * ri[1][1] - uno) + two * p2 * si[1][1];
                            let u11i = p2 * (uno - two * ri[1][1]) + two * p1 * si[1][1];
                            termt = two * gj * (uno - u11r);
                            termn = gj * (sq(uno - u11r) + sq(u11i));
                        } else {
                            let dd = r[1][1];
                            let rr = uno + dd;
                            let ss = s[1][1];
                            let amag = sq(rr) + sq(ss);
                            let rri = rr / amag;
                            let ssi = -ss / amag;
                            let uur = p1 * (two * rri - uno) + two * p2 * ssi;
                            let uui = p2 * (uno - two * rri) + two * p1 * ssi;
                            if dd.abs() < small && phid.abs() < small {
                                let mut xx = 2.0 * dd;
                                xx = xx + 2.0 * (dd * dd + ss * ss + phid * phid + p2 * ss);
                                xx = xx - 2.0 * phid * phid * (dd * dd + ss * ss);
                                xx = xx / amag;
                                termt = two * gj * xx;
                                termn = gj * (sq(xx) + sq(uui));
                            } else {
                                termt = two * gj * (uno - uur);
                                termn = gj * (sq(uno - uur) + sq(uui));
                            }
                            termf = 0.0;
                        }
                        if kk == 2 {
                            termn = termn + two * gj * (1.0 - p1);
                            termt = termt + two * gj * (1.0 - p1);
                        }
                        let termg = termt - termf - termn;
                        sigp[2] = sigp[2] + termn;
                        sigp[4] = sigp[4] + termg;
                        sigp[3] = sigp[3] + termf;
                        sigp[1] = sigp[1] + termt;
                        let _ = &mut termf;
                    }
                }
            }
            inow = in_;
        }
        sigp[1] = pifac * sigp[1];
        sigp[2] = pifac * sigp[2];
        sigp[3] = pifac * sigp[3];
        sigp[4] = pifac * sigp[4];
        self.sv.csr_ri = ri;
        self.sv.csr_si = si;
    }

    fn cshybr(&mut self, c: &mut Ctx, e: f64, sigp: &mut RV) {
        let pi = crate::util::physics::PI;
        let zero = 0.0f64;
        let uno = 1.0f64;
        let eye = Cx::new(zero, uno);
        let one = Cx::new(uno, zero);
        let cwaven = cwaven();
        if self.tempr > zero {
            c.error("cshybr", "doppler broadening not provided for hybrid r-mtrix.", " ");
        }
        for i in 1..=4 {
            sigp[i] = 0.0;
        }
        let mut lc = [0i32; 5];
        let mut icrpt = [[0i32; 5]; 5];
        let mut qre = [0.0f64; 5];
        let mut gc = [0.0f64; 5];
        let mut sigc = [0.0f64; 5];
        let mut inow = self.ibaset[self.isect as usize];
        let spi = self.res[inow + 6];
        let nls = nint(self.res[inow + 10]);
        let ngre = nint(self.res[inow + 14]);
        let nfre = nint(self.res[inow + 15]);
        let nire = nint(self.res[inow + 16]);
        let ncre = nint(self.res[inow + 17]);
        let nrchan = nire + ncre;
        qre[1] = self.res[inow + 29];
        qre[2] = self.res[inow + 30];
        qre[3] = self.res[inow + 31];
        qre[4] = self.res[inow + 32];
        if nrchan > 0 {
            for n in 1..=nrchan as usize {
                sigc[n] = 0.0;
            }
        }
        inow += 34;
        if ncre != 0 {
            for n in 1..=ncre as usize {
                for l in 1..=4 {
                    icrpt[l][n] = inow;
                    let nw = 6 + 2 * nint(self.res[inow + 4]) + 2 * nint(self.res[inow + 5]);
                    inow += nw;
                }
            }
        }
        let awri = self.res[inow];
        let arat = awri / (awri + 1.0);
        let kfac = cwaven * arat;
        let k = kfac * sqrt(e.abs());
        let pifac = pi / (k * k);
        let mut eeff = self.sv.hy_eeff;
        let mut eenn = 0.0;
        let mut idis = 0;

        'l140: for _l in 1..=nls {
            let ll = nint(self.res[inow + 2]);
            let nss = nint(self.res[inow + 4]);
            inow += 6;
            for _ns in 1..=nss {
                let njs = nint(self.res[inow + 4]);
                inow += 6;
                for _nj in 1..=njs {
                    let aj = self.res[inow];
                    let gj = (2.0 * aj + 1.0) / (4.0 * spi + 2.0);
                    let ac = self.res[inow + 1];
                    let lbk = nint(self.res[inow + 2]);
                    let lps = nint(self.res[inow + 3]);
                    let nrs = nint(self.res[inow + 5]);
                    inow += 6;
                    let rho = k * ac;
                    let phir = facphi(ll, rho);
                    let phi0 = Cx::new(phir, zero);
                    let (mut se, mut pe) = (0.0, 0.0);
                    facts(ll, rho, &mut se, &mut pe);
                    let mut rfcnn = Cx::new(zero, zero);
                    if nrs == 0 {
                        continue 'l140;
                    }
                    for _i in 1..=nrs {
                        let er = self.res[inow];
                        let gn = self.res[inow + 1];
                        let gg = self.res[inow + 2];
                        let gf = self.res[inow + 3];
                        inow += 3;
                        let mut gelim = gg + gf;
                        if nrchan != 0 {
                            for ii in 1..=nrchan as usize {
                                gc[ii] = 0.0;
                                lc[ii] = nint(self.res[inow + ii as i32 + 4]);
                                eeff = e + (qre[ii] / arat);
                                if eeff > zero {
                                    let ereff = (er + (qre[ii] / arat)).abs();
                                    if ii as i32 <= nire {
                                        let phres = kfac * sqrt(ereff) * ac;
                                        let (mut sd, mut pd) = (0.0, 0.0);
                                        facts(lc[ii], phres, &mut sd, &mut pd);
                                        let phen = kfac * sqrt(eeff) * ac;
                                        let (mut sn, mut pn) = (0.0, 0.0);
                                        facts(lc[ii], phen, &mut sn, &mut pn);
                                        gc[ii] = pn * self.res[inow + ii as i32] / pd;
                                    } else if lc[ii] + 1 <= 4 {
                                        let itpt = icrpt[(lc[ii] + 1) as usize][ii];
                                        let mut ip = 2;
                                        let mut ir = 1;
                                        let (mut pd, mut pn) = (0.0, 0.0);
                                        let sub = crate::farr::F(&mut self.res.0[(itpt - 1) as usize..]);
                                        c.terpa(&mut pd, ereff, &mut eenn, &mut idis, &sub, &mut ip, &mut ir);
                                        c.terpa(&mut pn, eeff, &mut eenn, &mut idis, &sub, &mut ip, &mut ir);
                                        gc[ii] = pn * self.res[inow + ii as i32] / pd;
                                    }
                                }
                                gelim = gelim + gc[ii];
                            }
                        }
                        inow += 9;
                        let rhres = kfac * sqrt(er.abs()) * ac;
                        let (mut ser, mut per) = (0.0, 0.0);
                        facts(ll, rhres, &mut ser, &mut per);
                        let gne = gn * pe / per;
                        let gnr = gn / (2.0 * per);
                        let gtt = gne + gelim;
                        let edelt = er - e;
                        let den = sq(edelt) + gtt * gtt / 4.0;
                        let comfac = pifac * gj * gne / den;
                        // rfcnn=rfcnn+gnr/(edelt-eye*gelim/2)
                        let t = eye.mul(Cx::r(gelim)).div(Cx::r(2.0));
                        let d = Cx::r(edelt).sub(t);
                        rfcnn = rfcnn.add(Cx::r(gnr).div(d));
                        if ngre != 0 {
                            sigp[4] = sigp[4] + comfac * gg;
                        }
                        if nfre != 0 {
                            sigp[3] = sigp[3] + comfac * gf;
                        }
                        if nrchan > 0 {
                            for iii in 1..=nrchan as usize {
                                sigc[iii] = sigc[iii] + comfac * gc[iii];
                            }
                        }
                    }
                    let mut r0 = Cx::new(zero, zero);
                    if lbk != 0 {
                        let mut ip = 2;
                        let mut ir = 1;
                        let (mut rr0, mut ri0) = (0.0, 0.0);
                        {
                            let sub = crate::farr::F(&mut self.res.0[(inow - 1) as usize..]);
                            c.terpa(&mut rr0, eeff, &mut eenn, &mut idis, &sub, &mut ip, &mut ir);
                        }
                        let np = nint(self.res[inow + 5]);
                        let nr = nint(self.res[inow + 4]);
                        inow = inow + 6 + 2 * nr + 2 * np;
                        {
                            let sub = crate::farr::F(&mut self.res.0[(inow - 1) as usize..]);
                            c.terpa(&mut ri0, eeff, &mut eenn, &mut idis, &sub, &mut ip, &mut ir);
                        }
                        let np = nint(self.res[inow + 5]);
                        let nr = nint(self.res[inow + 4]);
                        inow = inow + 6 + 2 * nr + 2 * np;
                        r0 = Cx::new(rr0, ri0);
                    }
                    rfcnn = rfcnn.add(r0);
                    let mut phi = phi0;
                    if lps != 0 {
                        let mut ip = 2;
                        let mut ir = 1;
                        let (mut phro, mut phri) = (0.0, 0.0);
                        {
                            let sub = crate::farr::F(&mut self.res.0[(inow - 1) as usize..]);
                            c.terpa(&mut phro, eeff, &mut eenn, &mut idis, &sub, &mut ip, &mut ir);
                        }
                        let np = nint(self.res[inow + 5]);
                        let nr = nint(self.res[inow + 4]);
                        inow = inow + 6 + 2 * nr + 2 * np;
                        {
                            let sub = crate::farr::F(&mut self.res.0[(inow - 1) as usize..]);
                            c.terpa(&mut phri, eeff, &mut eenn, &mut idis, &sub, &mut ip, &mut ir);
                        }
                        let np = nint(self.res[inow + 5]);
                        let nr = nint(self.res[inow + 4]);
                        inow = inow + 6 + 2 * nr + 2 * np;
                        phi = Cx::new(phro, phri);
                    }
                    // exphi=exp(-2*eye*phi)
                    let exphi = Cx::r(-2.0).mul(eye).mul(phi).exp();
                    // slsjnn=(one-exphi)-eye*pe*rfcnn*(one+exphi)
                    let t = eye.mul(Cx::r(pe)).mul(rfcnn).mul(one.add(exphi));
                    let mut slsjnn = one.sub(exphi).sub(t);
                    let dd = one.sub(eye.mul(Cx::r(pe)).mul(rfcnn));
                    slsjnn = slsjnn.div(dd);
                    let rpart = slsjnn.re;
                    let sabs = slsjnn.abs();
                    sigp[1] = sigp[1] + 2.0 * pifac * gj * rpart;
                    sigp[2] = sigp[2] + pifac * gj * sabs * sabs;
                    sigp[1] = sigp[2] + sigp[3] + sigp[4];
                }
            }
        }
        self.sv.hy_eeff = eeff;
        let _ = sigc;
    }

    /// Unresolved resonance cross sections (LRU=2, LRF=1 or LFW=0).
    /// `e` is passed by reference in Fortran and is changed on some paths.
    fn csunr1(&mut self, c: &mut Ctx, e_in: f64, sigp: &mut RV) {
        let mut e = e_in;
        self.csunr1_ref(c, &mut e, sigp);
    }

    fn csunr1_ref(&mut self, c: &mut Ctx, e: &mut f64, sigp: &mut RV) {
        let pi = crate::util::physics::PI;
        let wide = 1.26e0f64;
        let rc1 = 0.123e0f64;
        let rc2 = 0.08e0f64;
        let third = 0.333333333e0f64;
        let zero = 0.0f64;
        let cwaven = cwaven();
        let ne_init = 0;
        let mut ne = ne_init;
        if *e <= zero {
            for i in 1..=20 {
                self.ilast[i] = 0;
            }
            return;
        }
        let isect = self.isect;
        let is = isect as usize;
        let ee = *e;
        let mut ifirst = self.ilast[is];
        self.sv.ifirst_u1 = ifirst;
        let mut inow = self.ibaset[is];
        let mut lfw = nint(self.res[inow + 4]);
        let mut nro = 0;
        if self.res[inow + 7] < zero {
            nro = -nint(self.res[inow + 7]);
            inow = inow + 6 + 2 * nint(self.res[inow + 10]) + 2 * nint(self.res[inow + 11]);
        }
        let mut enext;
        let mut elast;
        let mut ie;
        let mut iest = self.sv.iest_u1;
        // unfac sets nothing for l >= 3: the Fortran then reuses the vl/ps
        // of the previous (l, j) step, so they must live for the whole call.
        let (mut vl, mut ps) = (0.0f64, 0.0f64);
        let mut lbl;
        if ifirst > 0 {
            // 100
            enext = self.enxt[is];
            elast = self.elst[is];
            ie = self.ilast[is];
            if lfw == 0 {
                lbl = 130;
            } else {
                iest = inow + 12;
                if *e < enext && enext >= wide * elast {
                    lbl = 130;
                } else if *e <= enext {
                    lbl = 210;
                } else if ie == ne {
                    lbl = 120;
                } else {
                    lbl = 110;
                }
            }
        } else {
            enext = *e;
            elast = 0.0;
            self.enxt[is] = *e;
            self.elst[is] = 0.0;
            ie = 1;
            lbl = 120;
        }
        loop {
            match lbl {
                110 => {
                    elast = enext;
                    self.elst[is] = elast;
                    for i in 1..=4 {
                        self.slst[(isect, i)] = self.snxt[(isect, i)];
                    }
                    ie += 1;
                    enext = self.res[iest + ie - 1];
                    self.enxt[is] = enext;
                    self.ilast[is] = ie;
                    lbl = 120;
                }
                120 => {
                    *e = enext;
                    lbl = 130;
                }
                130 => {
                    for i in 2..=4 {
                        sigp[i] = 0.0;
                    }
                    inow = self.ibaset[is];
                    lfw = nint(self.res[inow + 4]);
                    let naps = nint(self.res[inow + 5]);
                    let ay;
                    if self.res[inow + 7] > 0.0 {
                        ay = self.res[inow + 7];
                    } else {
                        let iro = inow + 6;
                        let mut ip = 2;
                        let mut ir = 1;
                        let (mut a, mut enx, mut idx) = (0.0, 0.0, 0);
                        {
                            let sub = crate::farr::F(&mut self.res.0[(iro - 1) as usize..]);
                            c.terpa(&mut a, *e, &mut enx, &mut idx, &sub, &mut ip, &mut ir);
                        }
                        ay = a;
                        inow = inow + 6 + 2 * nint(self.res[iro + 4]) + 2 * nint(self.res[iro + 5]);
                    }
                    let spi = self.res[inow + 6];
                    let aaa = self.res[inow + 7];
                    let nls;
                    let awri;
                    if lfw != 1 {
                        ne = 0;
                        nls = nint(self.res[inow + 10]);
                        inow += 12;
                        awri = self.res[inow];
                    } else {
                        ne = nint(self.res[inow + 10]);
                        nls = nint(self.res[inow + 11]);
                        iest = inow + 12;
                        inow = inow + 12 + ne;
                        awri = self.res[inow];
                    }
                    let rat = awri / (awri + 1.0);
                    let aw = awri * AMASSN;
                    let aa;
                    if naps == 0 {
                        aa = rc1 * pow(aw, third) + rc2;
                    } else if naps == 1 {
                        aa = ay;
                    } else if naps == 2 && nro == 1 {
                        aa = aaa;
                    } else {
                        c.error("csunr1", "illegal naps", " ");
                    }
                    let cnst = (2.0 * (pi * pi)) / sq(cwaven * rat);
                    let mut spot = self.sv.spot_u1;
                    for _l in 1..=nls {
                        let mut njs = 0;
                        if lfw == 0 {
                            njs = nint(self.res[inow + 5]);
                        }
                        if lfw == 1 {
                            njs = nint(self.res[inow + 4]);
                        }
                        let ll = nint(self.res[inow + 2]);
                        inow += 6;
                        for j in 1..=njs {
                            let mut gfx = 0.0;
                            let gxx = 0.0;
                            let muf;
                            let dx;
                            let aj;
                            let amun;
                            let gnox;
                            let ggx;
                            let mut ifst = 0;
                            if lfw <= 0 {
                                muf = 0;
                                dx = self.res[inow];
                                aj = self.res[inow + 1];
                                amun = self.res[inow + 2];
                                gnox = self.res[inow + 3];
                                ggx = self.res[inow + 4];
                                inow += 6;
                            } else {
                                muf = nint(self.res[inow + 3]);
                                let nep6 = nint(self.res[inow + 4]);
                                dx = self.res[inow + 6];
                                aj = self.res[inow + 7];
                                amun = self.res[inow + 8];
                                gnox = self.res[inow + 9];
                                ggx = self.res[inow + 10];
                                ifst = inow + 12;
                                inow = inow + 6 + nep6;
                            }
                            let gj = (2.0 * aj + 1.0) / (4.0 * spi + 2.0);
                            let mu = crate::fmath::int(amun);
                            let nu = muf;
                            let lamda = 0;
                            if lfw != 0 {
                                let mut i1 = iest;
                                if ie > 1 {
                                    i1 = iest + ie - 2;
                                }
                                let mut i2 = ifst;
                                if ie > 1 {
                                    i2 = ifst + ie - 2;
                                }
                                c.terp1(self.res[i1], self.res[i2], self.res[i1 + 1], self.res[i2 + 1], *e, &mut gfx, 2);
                            }
                            let e2 = sqrt(*e);
                            let k = rat * e2 * cwaven;
                            let rho = k * aa;
                            let rhoc = k * ay;
                            unfac(ll, rho, rhoc, amun, &mut vl, &mut ps);
                            vl = vl * e2;
                            if j == 1 {
                                spot = 4.0 * pi * (2 * ll + 1) as f64 * sq(crate::fmath::sin(ps) / k);
                            }
                            let gnx = gnox * vl;
                            let diff = gxx;
                            let den = *e * dx;
                            let temp = cnst * gj * gnx / den;
                            let terg = temp * ggx;
                            let ters = temp * gnx;
                            let terf = temp * gfx;
                            let mut gs = gnrl(gnx, gfx, ggx, mu, nu, lamda, diff, 1);
                            let mut gcc = gnrl(gnx, gfx, ggx, mu, nu, lamda, diff, 2);
                            let mut gff = gnrl(gnx, gfx, ggx, mu, nu, lamda, diff, 3);
                            gcc = gcc * terg;
                            gff = gff * terf;
                            gs = gs * ters;
                            let mut add = cnst * gj * 2.0 * gnx * sq(crate::fmath::sin(ps));
                            add = add / (*e * dx);
                            gs = gs - add;
                            sigp[2] = sigp[2] + gs;
                            sigp[3] = sigp[3] + gff;
                            sigp[4] = sigp[4] + gcc;
                        }
                        sigp[2] = sigp[2] + spot;
                    }
                    self.sv.spot_u1 = spot;
                    sigp[1] = sigp[2] + sigp[3] + sigp[4];
                    if ifirst == 0 && lfw == 0 {
                        ifirst = 1;
                    }
                    if lfw == 0 {
                        break;
                    }
                    if ifirst > 0 && *e < enext && enext > wide * elast {
                        break;
                    }
                    for i in 1..=4 {
                        self.snxt[(isect, i)] = sigp[i];
                    }
                    if ifirst > 0 && *e == enext && enext > wide * elast {
                        // 200
                        *e = ee;
                        lbl = 130;
                        continue;
                    }
                    if ifirst > 0 {
                        lbl = 210;
                        continue;
                    }
                    if ifirst < 0 {
                        // 205
                        ifirst = 1;
                        lbl = 210;
                        continue;
                    }
                    ifirst = -1;
                    lbl = 110;
                }
                _ => {
                    // 210
                    *e = ee;
                    for i in 2..=4 {
                        let y1 = self.slst[(isect, i)];
                        let y2 = self.snxt[(isect, i)];
                        let mut v = 0.0;
                        c.terp1(elast, y1, enext, y2, *e, &mut v, 2);
                        sigp[i] = v;
                    }
                    sigp[1] = sigp[2] + sigp[3] + sigp[4];
                    break;
                }
            }
        }
        self.sv.iest_u1 = iest;
        self.sv.ifirst_u1 = ifirst;
    }

    /// Unresolved resonance cross sections (LRU=2, LRF=2).
    fn csunr2(&mut self, c: &mut Ctx, e_in: f64, sigp: &mut RV) {
        let mut e = e_in;
        self.csunr2_ref(c, &mut e, sigp);
    }

    fn csunr2_ref(&mut self, c: &mut Ctx, e: &mut f64, sigp: &mut RV) {
        let pi = crate::util::physics::PI;
        let wide = 1.26e0f64;
        let rc1 = 0.123e0f64;
        let rc2 = 0.08e0f64;
        let third = 0.333333333e0f64;
        let small = 1.0e-8f64;
        let zero = 0.0f64;
        let cwaven = cwaven();
        if *e <= zero {
            for i in 1..=20 {
                self.ilast[i] = 0;
            }
            return;
        }
        let isect = self.isect;
        let is = isect as usize;
        let ee = *e;
        let mut ifirst = self.ilast[is];
        let mut inow = self.ibaset[is];
        let mut nro = 0;
        let inoww;
        if self.res[inow + 7] < 0.0 {
            nro = -nint(self.res[inow + 7]);
            inoww = inow + 6 + 2 * nint(self.res[inow + 10]) + 2 * nint(self.res[inow + 11]);
        } else {
            inoww = inow;
        }
        let mut int = nint(self.res[inoww + 20]);
        let mut ne = nint(self.res[inoww + 23]);
        let iloc = inoww + 30;
        let mut enext;
        let mut elast;
        let mut ie;
        let mut lbl;
        // unfac sets nothing for l >= 3: the Fortran then reuses the vl/ps
        // of the previous (l, j) step, so they must live for the whole call.
        let (mut vl, mut ps) = (0.0f64, 0.0f64);
        if ifirst > 0 {
            enext = self.enxt[is];
            elast = self.elst[is];
            ie = self.ilast[is];
            if *e < enext && enext >= wide * elast {
                lbl = 130;
            } else if *e <= enext {
                lbl = 200;
            } else if ie == ne {
                lbl = 120;
            } else {
                lbl = 110;
            }
        } else {
            enext = *e;
            elast = 0.0;
            self.enxt[is] = *e;
            self.elst[is] = 0.0;
            ie = 1;
            lbl = 120;
        }
        loop {
            match lbl {
                110 => {
                    elast = enext;
                    self.elst[is] = elast;
                    for i in 1..=4 {
                        self.slst[(isect, i)] = self.snxt[(isect, i)];
                    }
                    ie += 1;
                    enext = self.res[iloc + (ie - 1) * 6];
                    self.enxt[is] = enext;
                    self.ilast[is] = ie;
                    lbl = 120;
                }
                120 => {
                    *e = enext;
                    lbl = 130;
                }
                130 => {
                    for i in 2..=4 {
                        sigp[i] = 0.0;
                    }
                    inow = self.ibaset[is];
                    let naps = nint(self.res[inow + 5]);
                    let ay;
                    if self.res[inow + 7] > 0.0 {
                        ay = self.res[inow + 7];
                    } else {
                        let iro = inow + 6;
                        let mut ip = 2;
                        let mut ir = 1;
                        let (mut a, mut enx, mut idx) = (0.0, 0.0, 0);
                        {
                            let sub = crate::farr::F(&mut self.res.0[(iro - 1) as usize..]);
                            c.terpa(&mut a, *e, &mut enx, &mut idx, &sub, &mut ip, &mut ir);
                        }
                        ay = a;
                        inow = inow + 6 + 2 * nint(self.res[iro + 4]) + 2 * nint(self.res[iro + 5]);
                    }
                    let spi = self.res[inow + 6];
                    let aaa = self.res[inow + 7];
                    let awri = self.res[inow + 12];
                    let nls = nint(self.res[inow + 10]);
                    inow += 12;
                    let rat = awri / (awri + 1.0);
                    let aw = awri * AMASSN;
                    let aa;
                    if naps == 0 {
                        aa = rc1 * pow(aw, third) + rc2;
                    } else if naps == 1 {
                        aa = ay;
                    } else if naps == 2 && nro == 1 {
                        aa = aaa;
                    } else {
                        c.error("csunr2", "illegal naps", " ");
                    }
                    let cnst = (2.0 * (pi * pi)) / sq(cwaven * rat);
                    let e2 = sqrt(*e);
                    let k = cwaven * rat * e2;
                    let rho = k * aa;
                    let rhoc = k * ay;
                    let mut spot = self.sv.spot_u2;
                    for _l in 1..=nls {
                        let njs = nint(self.res[inow + 4]);
                        let ll = nint(self.res[inow + 2]);
                        inow += 6;
                        for j in 1..=njs {
                            let aj = self.res[inow];
                            int = nint(self.res[inow + 2]);
                            ne = nint(self.res[inow + 5]);
                            let amux = self.res[inow + 8];
                            let amun = self.res[inow + 9];
                            let amuf = self.res[inow + 11];
                            let jnow = nint(self.res[inow + 4]) + inow + 6;
                            inow += 12;
                            let gj = (2.0 * aj + 1.0) / (4.0 * spi + 2.0);
                            let mu = nint(amun);
                            let nu = nint(amuf);
                            let lamda = nint(amux);
                            let mut i1 = inow;
                            if i1 + 6 < jnow {
                                while self.res[i1 + 6] <= *e - small {
                                    i1 += 6;
                                }
                            }
                            let i2 = i1 + 6;
                            let (mut dx, mut gxx, mut gnox, mut ggx, mut gfx) = (0.0, 0.0, 0.0, 0.0, 0.0);
                            c.terp1(self.res[i1], self.res[i1 + 1], self.res[i2], self.res[i2 + 1], *e, &mut dx, 2);
                            c.terp1(self.res[i1], self.res[i1 + 2], self.res[i2], self.res[i2 + 2], *e, &mut gxx, 2);
                            c.terp1(self.res[i1], self.res[i1 + 3], self.res[i2], self.res[i2 + 3], *e, &mut gnox, 2);
                            c.terp1(self.res[i1], self.res[i1 + 4], self.res[i2], self.res[i2 + 4], *e, &mut ggx, 2);
                            c.terp1(self.res[i1], self.res[i1 + 5], self.res[i2], self.res[i2 + 5], *e, &mut gfx, 2);
                            if gxx < small {
                                gxx = 0.0;
                            }
                            if gfx < small {
                                gfx = 0.0;
                            }
                            unfac(ll, rho, rhoc, amun, &mut vl, &mut ps);
                            vl = vl * e2;
                            if j <= 1 {
                                spot = 4.0 * pi * (2 * ll + 1) as f64 * sq(crate::fmath::sin(ps) / k);
                            }
                            let gnx = gnox * vl;
                            let gamma = ggx;
                            let galpha = gnx;
                            let gbeta = gfx;
                            let diff = gxx;
                            let den = *e * dx;
                            let temp = cnst * gj * gnx / den;
                            let terg = temp * gamma;
                            let ters = temp * gnx;
                            let terf = temp * gbeta;
                            let mut gs = gnrl(galpha, gbeta, gamma, mu, nu, lamda, diff, 1);
                            let mut gcc = gnrl(galpha, gbeta, gamma, mu, nu, lamda, diff, 2);
                            let mut gff = gnrl(galpha, gbeta, gamma, mu, nu, lamda, diff, 3);
                            gcc = gcc * terg;
                            gff = gff * terf;
                            gs = gs * ters;
                            let mut add = cnst * gj * 2.0 * gnx * sq(crate::fmath::sin(ps));
                            add = add / (*e * dx);
                            gs = gs - add;
                            sigp[2] = sigp[2] + gs;
                            sigp[3] = sigp[3] + gff;
                            sigp[4] = sigp[4] + gcc;
                            inow = jnow;
                        }
                        sigp[2] = sigp[2] + spot;
                    }
                    self.sv.spot_u2 = spot;
                    sigp[1] = sigp[2] + sigp[3] + sigp[4];
                    if ifirst > 0 && *e < enext && enext >= wide * elast {
                        break;
                    }
                    for i in 1..=4 {
                        self.snxt[(isect, i)] = sigp[i];
                    }
                    if ifirst > 0 && *e == enext && enext >= wide * elast {
                        *e = ee;
                        lbl = 130;
                        continue;
                    }
                    if ifirst > 0 {
                        lbl = 200;
                        continue;
                    }
                    if ifirst < 0 {
                        ifirst = 1;
                        lbl = 200;
                        continue;
                    }
                    ifirst = -1;
                    lbl = 110;
                }
                _ => {
                    // 200
                    *e = ee;
                    for i in 2..=4 {
                        let y1 = self.slst[(isect, i)];
                        let y2 = self.snxt[(isect, i)];
                        let mut v = 0.0;
                        c.terp1(elast, y1, enext, y2, *e, &mut v, int);
                        sigp[i] = v;
                    }
                    sigp[1] = sigp[2] + sigp[3] + sigp[4];
                    break;
                }
            }
        }
    }

    /// Adler-Adler cross sections.
    fn csaa(&mut self, c: &mut Ctx, e: f64, sigp: &mut RV) {
        let pi = crate::util::physics::PI;
        let bk = crate::util::physics::BK;
        let zero = 0.0f64;
        let cadler = 6.5099897e5f64;
        let rpi = sqrt(pi);
        let cwaven = cwaven();
        for i in 1..=4 {
            sigp[i] = 0.0;
        }
        let mut inow = self.ibaset[self.isect as usize];
        let awri = self.res[inow + 12];
        let ap = self.res[inow + 7];
        let nls = nint(self.res[inow + 10]);
        let li = nint(self.res[inow + 14]);
        inow += 18;
        let mut delta = self.sv.aa_delta;
        if self.tempr > zero {
            delta = 1.0 / sqrt(4.0 * bk * self.tempr * e / awri);
        }
        self.sv.aa_delta = delta;
        let arat = awri / (awri + 1.0);
        let eroot = sqrt(e.abs());
        let esq = e * e;
        let ecub = esq * e;
        let k = cwaven * arat * eroot;
        let cc = cadler / (arat * arat);
        let omg = 2.0 * k * ap;
        let snf = crate::fmath::sin(omg);
        let csf = crate::fmath::cos(omg);
        let mut bakt = 0.0;
        let mut bakf = 0.0;
        if li < 5 {
            c.error("csaa", "bad li value", " ");
        }
        let r = &self.res;
        if li != 6 {
            bakt = r[inow] + r[inow + 1] / e + r[inow + 2] / esq;
            bakt = bakt + r[inow + 3] / ecub + r[inow + 4] * e + r[inow + 5] * esq;
            bakt = bakt * cc / eroot;
        }
        inow += 6;
        if li != 5 {
            bakf = r[inow] + r[inow + 1] / e + r[inow + 2] / esq;
            bakf = bakf + r[inow + 3] / ecub + r[inow + 4] * e + r[inow + 5] * esq;
            bakf = bakf * cc / eroot;
        }
        inow += 6;
        let mut bakc = r[inow] + r[inow + 1] / e + r[inow + 2] / esq;
        bakc = bakc + r[inow + 3] / ecub + r[inow + 4] * e + r[inow + 5] * esq;
        bakc = bakc * cc / eroot;
        inow += 6;
        let mut ki = self.sv.aa_ki.max(1);
        for _l in 1..=nls {
            let njs = nint(self.res[inow + 4]);
            inow += 6;
            for _j in 1..=njs {
                let nlj = nint(self.res[inow + 5]);
                inow += 6;
                for _n in 1..=nlj {
                    let mut ist = inow;
                    for i in 1..=4 {
                        if i != 2 {
                            let de = self.res[ist];
                            let dw = self.res[ist + 1];
                            let gr = self.res[ist + 2];
                            let gi = self.res[ist + 3];
                            ist += 4;
                            let x = (de - e) / dw;
                            let psi;
                            let chi;
                            if self.tempr <= zero {
                                psi = 1.0 / (1.0 + x * x);
                                chi = x * psi;
                            } else {
                                let theta = 2.0 * dw * delta;
                                let ax = theta * x / 2.0;
                                let y = theta / 2.0;
                                let (mut rew, mut aimw) = (0.0, 0.0);
                                quickw(ax, y, &mut rew, &mut aimw, &mut ki, self.tr.as_ref().unwrap(), self.ti.as_ref().unwrap());
                                psi = rpi * theta * rew / 2.0;
                                chi = rpi * theta * aimw / 2.0;
                            }
                            if i == 1 {
                                sigp[i] = sigp[i] + ((gr * csf + gi * snf) * psi + (gi * csf - gr * snf) * chi) / dw;
                            }
                            if i > 1 {
                                sigp[i] = sigp[i] + (gr * psi + gi * chi) / dw;
                            }
                        }
                    }
                    inow = ist;
                }
            }
        }
        self.sv.aa_ki = ki;
        sigp[1] = (sigp[1] + bakt) * cc / eroot;
        sigp[1] = sigp[1] + 2.0 * cc * (1.0 - csf) / e;
        sigp[3] = (sigp[3] + bakf) * cc / eroot;
        sigp[4] = (sigp[4] + bakc) * cc / eroot;
        if li != 6 {
            sigp[2] = sigp[1] - (sigp[3] + sigp[4]);
        } else {
            sigp[1] = sigp[3] + sigp[4];
        }
    }
}

pub fn facphi(l: i32, rho: f64) -> f64 {
    let test = 1.0e-6f64;
    let atan = crate::fmath::atan;
    let r2 = rho * rho;
    let mut phi;
    if l == 0 {
        phi = rho;
    } else if l == 1 {
        phi = rho - atan(rho);
    } else if l == 2 {
        phi = rho - atan(3.0 * rho / (3.0 - r2));
        if phi / rho < test {
            phi = 0.0;
        }
    } else if l == 3 {
        phi = rho - atan((15.0 * rho - rho * r2) / (15.0 - 6.0 * r2));
        if phi / rho < test {
            phi = 0.0;
        }
    } else {
        let r4 = r2 * r2;
        let top = 105.0 * rho - 10.0 * r2 * rho;
        let bot = 105.0 - 45.0 * r2 + r4;
        phi = rho - atan(top / bot);
        if phi / rho < test {
            phi = 0.0;
        }
    }
    phi
}

pub fn unfac(l: i32, rho: f64, rhoc: f64, amun: f64, vl: &mut f64, ps: &mut f64) {
    let atan = crate::fmath::atan;
    let r2 = rho * rho;
    if l == 0 {
        *vl = amun;
        *ps = rhoc;
    } else if l == 1 {
        *vl = amun * r2 / (1.0 + r2);
        *ps = rhoc - atan(rhoc);
    } else if l == 2 {
        let r4 = r2 * r2;
        *vl = amun * r4 / (9.0 + 3.0 * r2 + r4);
        *ps = rhoc - atan(3.0 * rhoc / (3.0 - rhoc * rhoc));
    }
}

const QW: [f64; 40] = [
    1.1120413e-1, 2.3546798e-1, 2.8440987e-1, 2.2419127e-1, 0.10967668e0, 0.030493789e0,
    0.0042930874e0, 2.5827047e-4, 4.9031965e-6, 1.4079206e-8, 0.033773418e0, 0.079932171e0,
    0.12835937e0, 0.17652616e0, 0.21347043e0, 0.21154965e0, 0.13365186e0, 0.022630659e0,
    1.6313638e-5, 2.745383e-31, 3.3376214e-4, 0.018506108e0, 0.12309946e0, 0.29918923e0,
    0.33431475e0, 0.17766657e0, 0.042695894e0, 4.0760575e-3, 1.1766115e-4, 5.0989546e-7,
    1.7623788e-3, 0.021517749e0, 0.080979849e0, 0.18797998e0, 0.30156335e0, 0.29616091e0,
    0.10775649e0, 2.5171914e-3, 8.9630388e-10, 0.0e0,
];
const QP: [f64; 40] = [
    3.0013465e-3, 7.8592886e-2, 4.3282415e-1, 1.3345267e0, 3.0481846e0, 5.8263198e0, 9.9452656e0,
    1.5782128e1, 23.996824e0, 36.216208e0, 1.3219203e-2, 7.2349624e-2, 0.19089473e0, 0.39528842e0,
    0.74083443e0, 1.3498293e0, 2.5297983e0, 5.2384894e0, 13.821772e0, 75.647525e0, 1.0004488e-3,
    0.026197629e0, 0.14427472e0, 0.44484223e0, 1.0160615e0, 1.9421066e0, 3.3150885e0, 5.2607092e0,
    7.9989414e0, 12.072069e0, 0.013219203e0, 0.072349624e0, 0.19089473e0, 0.39528842e0,
    0.74083443e0, 1.3498293e0, 2.5297983e0, 5.2384894e0, 13.821772e0, 75.647525e0,
];

#[inline]
fn qw(j: i32, m: i32) -> f64 {
    QW[((j - 1) + 10 * (m - 1)) as usize]
}
#[inline]
fn qp(j: i32, m: i32) -> f64 {
    QP[((j - 1) + 10 * (m - 1)) as usize]
}

/// Fluctuation integrals for the unresolved range.
pub fn gnrl(galpha: f64, gbeta: f64, gamma: f64, mu: i32, nu: i32, lamda: i32, df: f64, id: i32) -> f64 {
    let zero = 0.0f64;
    let mut s = 0.0f64;
    if galpha <= zero || gamma <= zero || gbeta < zero {
        return s;
    }
    if gbeta > zero && df < zero {
        return s;
    }
    if gbeta == zero && df == zero {
        if id == 1 {
            for j in 1..=10 {
                s = s + qw(j, mu) * qp(j, mu) * qp(j, mu) / (galpha * qp(j, mu) + gamma);
            }
        } else if id == 2 {
            for j in 1..=10 {
                s = s + qw(j, mu) * qp(j, mu) / (galpha * qp(j, mu) + gamma);
            }
        }
        return s;
    }
    if gbeta == zero && df > zero {
        if id == 1 {
            for j in 1..=10 {
                let xj = qp(j, mu);
                for k in 1..=10 {
                    s = s + qw(j, mu) * qw(k, lamda) * xj * xj / (galpha * xj + gamma + df * qp(k, lamda));
                }
            }
        } else if id == 2 {
            for j in 1..=10 {
                let xj = qp(j, mu);
                for k in 1..=10 {
                    s = s + qw(j, mu) * qw(k, lamda) * xj / (galpha * xj + gamma + df * qp(k, lamda));
                }
            }
        }
        return s;
    }
    if gbeta > zero && df == 0.0 {
        if id == 1 {
            for j in 1..=10 {
                let xj = qp(j, mu);
                for k in 1..=10 {
                    s = s + qw(j, mu) * qw(k, nu) * xj * xj / (galpha * xj + gbeta * qp(k, nu) + gamma);
                }
            }
        } else if id == 2 {
            for j in 1..=10 {
                let xj = qp(j, mu);
                for k in 1..=10 {
                    s = s + qw(j, mu) * qw(k, nu) * xj / (galpha * xj + gbeta * qp(k, nu) + gamma);
                }
            }
        } else if id == 3 {
            for j in 1..=10 {
                for k in 1..=10 {
                    s = s + qw(j, mu) * qw(k, nu) * qp(j, mu) * qp(k, nu) / (galpha * qp(j, mu) + gbeta * qp(k, nu) + gamma);
                }
            }
        }
        return s;
    }
    if gbeta > zero && df > zero {
        if id == 1 {
            for j in 1..=10 {
                let xj = qp(j, mu);
                for k in 1..=10 {
                    for l in 1..=10 {
                        s = s + qw(j, mu) * qw(k, nu) * qw(l, lamda) * xj * xj
                            / (galpha * xj + gbeta * qp(k, nu) + gamma + df * qp(l, lamda));
                    }
                }
            }
        } else if id == 2 {
            for j in 1..=10 {
                for k in 1..=10 {
                    for l in 1..=10 {
                        s = s + qw(j, mu) * qw(k, nu) * qw(l, lamda) * qp(j, mu)
                            / (galpha * qp(j, mu) + gbeta * qp(k, nu) + gamma + df * qp(l, lamda));
                    }
                }
            }
        } else if id == 3 {
            for j in 1..=10 {
                for k in 1..=10 {
                    for l in 1..=10 {
                        s = s + qw(j, mu) * qw(k, nu) * qw(l, lamda) * qp(j, mu) * qp(k, nu)
                            / (galpha * qp(j, mu) + gbeta * qp(k, nu) + gamma + df * qp(l, lamda));
                    }
                }
            }
        }
        return s;
    }
    s
}

/// 3x3 complex matrix inversion (Frobenius-Schur): (a + i b)^-1 = c + i d.
pub fn frobns(a: &mut [[f64; 4]; 4], b: &[[f64; 4]; 4], c: &mut [[f64; 4]; 4], d: &mut [[f64; 4]; 4]) {
    let mut q = [[0.0f64; 4]; 4];
    for i in 1..=3 {
        for j in 1..=3 {
            c[i][j] = a[i][j];
        }
    }
    let mut ind = 0;
    thrinv(a, 3, &mut ind);
    if ind != 1 {
        abcmat(a, b, &mut q);
        abcmat(b, &q, d);
        for i in 1..=3 {
            for j in 1..=3 {
                c[i][j] = c[i][j] + d[i][j];
            }
        }
        thrinv(c, 3, &mut ind);
        abcmat(&q, c, d);
        for i in 1..=3 {
            for j in 1..=3 {
                d[i][j] = -d[i][j];
            }
        }
    }
}

pub fn thrinv(d: &mut [[f64; 4]; 4], n: usize, kimerr: &mut i32) {
    let zero = 0.0f64;
    let uno = 1.0f64;
    let mut s = [0.0f64; 4];
    *kimerr = 0;
    for j in 1..=n {
        for i in 1..=j {
            d[i][j] = -d[i][j];
            d[j][i] = d[i][j];
        }
        d[j][j] = uno + d[j][j];
    }
    for lr in 1..=n {
        let fooey = uno - d[lr][lr];
        if fooey == zero {
            *kimerr = 1;
            return;
        }
        d[lr][lr] = uno / fooey;
        for j in 1..=n {
            s[j] = d[lr][j];
            if j != lr {
                d[j][lr] = d[j][lr] * d[lr][lr];
                d[lr][j] = d[j][lr];
            }
        }
        for j in 1..=n {
            if j != lr {
                for i in 1..=j {
                    if i != lr {
                        d[i][j] = d[i][j] + d[i][lr] * s[j];
                        d[j][i] = d[i][j];
                    }
                }
            }
        }
    }
}

pub fn abcmat(a: &[[f64; 4]; 4], b: &[[f64; 4]; 4], c: &mut [[f64; 4]; 4]) {
    for i in 1..=3 {
        for j in 1..=3 {
            c[i][j] = 0.0;
            for k in 1..=3 {
                c[i][j] = c[i][j] + a[i][k] * b[k][j];
            }
        }
    }
}

// ======================================================================
// emerge, recout
// ======================================================================

impl Reconr {
    fn emerge(
        &mut self,
        c: &mut Ctx,
        nin: i32,
        ngrid: i32,
        ngp: i32,
        nres: i32,
        nrtot: i32,
        iold_out: &mut i32,
        inew_out: &mut i32,
        nscr: i32,
    ) {
        let test = 1.0e-10f64;
        let finity = 0.99e12f64;
        let small = 1.0e-8f64;
        let zero = 0.0f64;
        let one = 1.0f64;
        let nsig = self.nsig;
        let mut bufo = vec![0.0f64; NBUF as usize];
        let mut bufn = vec![0.0f64; NBUF as usize];
        let mut bufg = vec![0.0f64; NBUFG as usize];
        let mut bufr = vec![0.0f64; NBUFR as usize];
        let mut scr = RV::new((c.npage + 50) as usize + SLACK);
        let mut res = vec![0.0f64; (nsig + 1) as usize];
        let mut tot = [0.0f64; 13]; // tot(1..) at tot[0..]
        let mut aa = [0.0f64; 1];
        let mut nneg = 0;
        let ntot = self.nmtr + 1;
        for i in 1..=self.nmtr {
            tot[i as usize] = 0.0;
        }
        let mut iold = 14;
        let mut inew = 15;

        let mut in_ = 0;
        let mut ig = 1;
        c.finda(ig, &mut aa, 1, ngrid, &mut bufg, NBUFG);
        let mut eg = aa[0];
        let mut er = finity;
        let mut ir = 0;
        if nrtot > 0 {
            ir = 1;
            c.finda(ir, &mut res, nsig + 1, nres, &mut bufr, NBUFR);
            er = res[0];
        }
        // labels 130-160
        let mut en;
        'l130: loop {
            in_ += 1;
            en = eg;
            if er < eg {
                en = er;
            }
            tot[0] = en;
            let mut inn = in_;
            if ig == ngp && en == eg {
                inn = -in_;
            }
            c.loada(inn, &tot, ntot, iold, &mut bufo, NBUF);
            // 135
            loop {
                if en != er {
                    break;
                }
                if ir < nrtot {
                    ir += 1;
                    c.finda(ir, &mut res, nsig + 1, nres, &mut bufr, NBUFR);
                    er = res[0];
                    continue;
                }
                er = finity;
                break;
            }
            // 150
            loop {
                if en != eg {
                    continue 'l130;
                }
                if ig == ngp {
                    break 'l130;
                }
                ig += 1;
                c.finda(ig, &mut aa, 1, ngrid, &mut bufg, NBUFG);
                eg = aa[0];
            }
        }
        // 160
        self.ngo = in_;

        c.openz(nscr, 1);
        c.repoz(nin);
        c.repoz(nscr);
        c.nsc = 0;
        let (mut nb, mut nw) = (0i32, 0i32);
        let mut nss = 1;
        let mut iss = 0;
        let mut awrx = 0.0f64;
        let mut lbl = 210;
        loop {
            match lbl {
                210 => {
                    c.contio(nin, 0, nscr, &mut scr.all(), &mut nb, &mut nw);
                    nss = 1;
                    if c.mfh == 10 {
                        nss = c.n1h;
                    }
                    iss = nss;
                    if c.math > 0 {
                        lbl = 220;
                        continue;
                    }
                    c.atend(0, nscr);
                    break;
                }
                _ => {
                    // 220
                    if c.mth == 0 {
                        lbl = 210;
                        continue;
                    }
                    if iss == nss {
                        if c.mfh != 23 && self.awin != zero {
                            awrx = c.c2h / self.awin;
                        }
                        let qx = c.c2h;
                        if self.awin != 0.0 {
                            c.thr6 = -qx * (awrx + 1.0) / awrx;
                        } else {
                            c.thr6 = -qx;
                        }
                        if c.thr6 < zero {
                            c.thr6 = 0.0;
                        }
                    }
                    let mut e = 0.0;
                    let (mut thresh, mut idis, mut sn) = (0.0, 0, 0.0);
                    c.gety1(e, &mut thresh, &mut idis, &mut sn, nin, &mut scr.all());
                    thresh = sigfig(thresh, 7, 0);

                    self.itype = 0;
                    if nrtot != 0 && c.mfh == 3 {
                        let mth = c.mth;
                        if mth == 2 {
                            self.itype = 2;
                        }
                        if mth == 18 {
                            self.itype = 3;
                        }
                        if mth == 19 {
                            self.itype = 3;
                        }
                        if mth == 102 {
                            self.itype = 4;
                        }
                        if self.nmtres > 2 {
                            for i in 3..=self.nmtres {
                                if mth != 18 && mth != 19 && mth == self.mmtres[i as usize] {
                                    self.itype = 2 + i;
                                }
                            }
                        }
                        if self.lrx != 0 && mth == 51 {
                            self.itype = 5;
                        }
                    }
                    let itype = self.itype;
                    er = finity;
                    if itype != 0 {
                        ir = 1;
                        c.finda(ir, &mut res, nsig + 1, nres, &mut bufr, NBUFR);
                        er = res[0];
                        if er < thresh {
                            thresh = sigfig(er, 7, 0);
                        }
                    }
                    in_ = 0;
                    let mut ith = 0;
                    let mut it = 0;
                    ig = 0;
                    let mth = c.mth;
                    let mfh = c.mfh;
                    loop {
                        // 340
                        ig += 1;
                        c.finda(ig, &mut tot, ntot, iold, &mut bufo, NBUF);
                        eg = tot[0];
                        let sg = tot[1];
                        sn = 0.0;
                        if !(thresh - eg > test * thresh && itype == 0) {
                            let mut enext = 0.0;
                            c.gety1(eg, &mut enext, &mut idis, &mut sn, nin, &mut scr.all());
                            if thresh > one && (thresh - eg).abs() < test * thresh {
                                sn = 0.0;
                            }
                            if eg >= self.eresr && eg < self.eresh && itype > 0 && itype < 5 {
                                sn = 0.0;
                            }
                            // 345
                            loop {
                                if er - eg > test * eg {
                                    break;
                                }
                                if (eg - er).abs() < test * eg {
                                    // 355
                                    sn = sn + res[itype as usize];
                                    if ir >= nrtot {
                                        er = finity;
                                    } else {
                                        ir += 1;
                                        c.finda(ir, &mut res, nsig + 1, nres, &mut bufr, NBUFR);
                                        er = res[0];
                                    }
                                    break;
                                }
                                ir += 1;
                                c.finda(ir, &mut res, nsig + 1, nres, &mut bufr, NBUFR);
                                er = res[0];
                            }
                            // 360
                            in_ += 1;
                            if itype == 2 && sn <= small {
                                nneg += 1;
                                if nneg == 1 {
                                    if sn <= 0.0 {
                                        c.mess("emerge", "nonpositive elastic cross sections found.", " ");
                                    } else {
                                        c.mess("emerge", "tiny elastic cross sections replaced by the lower limit.", " ");
                                    }
                                }
                                sn = small;
                            }
                            sn = sigfig(sn, 7, 0);
                            tot[1] = sn;
                            if ith == 0 && sn > zero {
                                ith = in_;
                            }
                            let mut inn = in_;
                            if ig == self.ngo {
                                inn = -in_;
                            }
                            c.loada(inn, &tot, 2, ngrid, &mut bufg, NBUFG);
                            tot[1] = sg;
                        }
                        // 370
                        it += 1;
                        let skip440 = (thresh - eg > test * thresh && itype == 0) || mth == 10;
                        if !skip440 {
                            if mfh == 3 {
                                if !((46..=49).contains(&mth) || (mth > 200 && mth < self.mpmin)) {
                                    for imtr in 1..=self.nmtr as usize {
                                        let m = self.mtr[imtr];
                                        let add = match m {
                                            1 => true,
                                            3 => mth != 2,
                                            4 => !(mth < 51 || mth > 91) && sn != zero,
                                            103 => !(mth < self.mpmin || mth > self.mpmax) && sn != zero,
                                            104 => !(mth < self.mdmin || mth > self.mdmax) && sn != zero,
                                            105 => !(mth < self.mtmin || mth > self.mtmax) && sn != zero,
                                            106 => !(mth < self.m3min || mth > self.m3max) && sn != zero,
                                            107 => !(mth < self.m4min || mth > self.m4max) && sn != zero,
                                            18 => !(mth < 19 || (mth > 21 && mth != 38)),
                                            _ => false,
                                        };
                                        if add {
                                            // 400
                                            if it < self.mtrt[imtr] {
                                                self.mtrt[imtr] = it;
                                                if self.mtrt[imtr] > 1 {
                                                    self.mtrt[imtr] -= 1;
                                                }
                                            }
                                            tot[imtr] = tot[imtr] + sn;
                                        }
                                    }
                                }
                            }
                            // 430
                            if mfh == 23 {
                                for imtr in 1..=self.nmtr as usize {
                                    let m = self.mtr[imtr];
                                    let add = if m == 501 {
                                        !(mth == 515 || mth == 517)
                                    } else if m == 522 {
                                        !(mth < 534 || mth > 572)
                                    } else {
                                        false
                                    };
                                    if add {
                                        if it < self.mtrt[imtr] {
                                            self.mtrt[imtr] = it;
                                            if self.mtrt[imtr] > 1 {
                                                self.mtrt[imtr] -= 1;
                                            }
                                        }
                                        tot[imtr] = tot[imtr] + sn;
                                    }
                                }
                            }
                        }
                        // 440
                        let mut itt = it;
                        if ig == self.ngo {
                            itt = -it;
                        }
                        c.loada(itt, &tot, ntot, inew, &mut bufn, NBUF);
                        if ig < self.ngo {
                            continue;
                        }
                        break;
                    }
                    if ith > 1 {
                        ith -= 1;
                    }
                    if iss == nss {
                        self.nxc += 1;
                    }
                    let nxc = self.nxc;
                    self.mfs[nxc] = mfh;
                    self.mts[nxc] = mth;
                    if iss == nss {
                        self.ncs[nxc] = 1;
                    }
                    self.ncs[nxc] = self.ncs[nxc] + 2 + ((in_ - ith + 1) + 2) / 3;
                    scr[5] = 1.0;
                    scr[6] = (in_ - ith + 1) as f64;
                    scr[7] = (in_ - ith + 1) as f64;
                    scr[8] = 2.0;
                    let mut k = 8;
                    if ith > 0 {
                        for ib in 1..=ith {
                            c.finda(ib, &mut tot, 2, ngrid, &mut bufg, NBUFG);
                        }
                    }
                    let mut istart = ith;
                    nb = 1;
                    while nb != 0 {
                        let mut iend = in_;
                        if iend - istart >= c.npage / 2 {
                            iend = istart + c.npage / 2 - 1;
                        }
                        let mut j = k - 1;
                        let mut ib = istart - 1;
                        while ib < iend {
                            j += 2;
                            ib += 1;
                            c.finda(ib, &mut tot, 2, ngrid, &mut bufg, NBUFG);
                            scr[j] = tot[0];
                            scr[j + 1] = tot[1];
                        }
                        nw = j + 1;
                        if k != 0 {
                            k = 0;
                            c.tab1io(0, 0, nscr, &mut scr.all(), &mut nb, &mut nw);
                        } else {
                            c.moreio(0, 0, nscr, &mut scr.all(), &mut nb, &mut nw);
                        }
                        if nb != 0 {
                            istart = iend + 1;
                        }
                    }
                    iss -= 1;
                    if iss > 0 {
                        lbl = 220;
                        continue;
                    }
                    c.tosend(nin, 0, 0, &mut scr.all());
                    c.asend(0, nscr);
                    std::mem::swap(&mut iold, &mut inew);
                    let _ = e;
                    e = 0.0;
                    let _ = e;
                    lbl = 210;
                }
            }
        }
        // 510
        if nneg > 0 {
            c.io.so(&format!("\n number of nonpositive cross sections removed ={}", ffmt::i(nneg, 8)));
        }
        c.io.so(&format!(" number of points in final unionized grid     ={}", ffmt::i(self.ngo, 8)));
        c.closz(-inew);
        c.closz(nin);
        c.closz(-ngrid);
        c.closz(-nres);
        *iold_out = iold;
        *inew_out = inew;
    }

    fn recout(&mut self, c: &mut Ctx, iold: i32, nscr: i32, nrtot: i32) {
        let emin = 1.0e-5f64;
        let zero = 0.0f64;
        let nout = self.nout;
        let mut scr = RV::new((c.npage + 50) as usize + SLACK);
        let mut bufo = vec![0.0f64; NBUF as usize];
        let mut bufl = vec![0.0f64; NBUFL as usize];
        let nresl = (1 + self.ncoef * self.nmtres).max(1) as usize;
        let mut resl = vec![0.0f64; nresl];
        let mut tot = [0.0f64; 13];
        let (mut nb, mut nw) = (0i32, 0i32);
        let mut i152 = 0;
        if self.nunr > 0 {
            i152 = 1;
        }
        let iang = 16;

        let mut no3 = 0;
        if nint(self.zain) == 1 {
            for i in 1..=self.nxc {
                if self.mfs[i] == 3 {
                    no3 = 1;
                }
            }
            if no3 == 0 {
                self.nxc += 1;
            }
        }
        let nxc = self.nxc;
        let nmtr = self.nmtr;

        scr[1] = self.za;
        scr[2] = self.awr;
        if c.iverf == 6 {
            self.lrp = 2;
        }
        scr[3] = self.lrp as f64;
        scr[4] = self.lfi as f64;
        scr[5] = 0.0;
        scr[6] = 0.0;
        if c.iverf == 4 {
            scr[6] = (nxc + nmtr + i152) as f64;
        }
        c.math = self.mata;
        c.mfh = 1;
        c.mth = 451;
        c.nsh = 1;
        nw = 6;
        c.contio(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
        if c.iverf != 4 {
            scr[1] = self.elis;
            scr[2] = self.sta;
            scr[3] = self.lis as f64;
            scr[4] = self.lis0 as f64;
            scr[5] = 0.0;
            scr[6] = self.nfor as f64;
            if c.iverf == 6 {
                scr[6] = 6.0;
            }
            c.contio(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
        }
        if c.iverf == 6 {
            scr[1] = self.awin;
            scr[2] = self.efmax;
            scr[3] = self.lrel as f64;
            scr[4] = 0.0;
            scr[5] = crate::fmath::int(10.0 * self.zain) as f64;
            scr[6] = self.nver as f64;
            c.contio(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
        }
        scr[1] = self.tempr;
        if self.tempr == 0.0 {
            scr[1] = self.tempi;
        }
        scr[2] = self.err;
        scr[3] = 0.0;
        if c.iverf == 6 {
            scr[3] = 1.0;
        }
        scr[4] = 0.0;
        let nwd = self.ncards * 17;
        scr[5] = nwd as f64;
        scr[6] = 0.0;
        if c.iverf >= 5 {
            scr[6] = (nxc + nmtr + i152) as f64;
        }
        self.ncs[1] = self.ncs[1] + nxc + nmtr + i152;
        if self.ncoef > 1 {
            scr[6] = scr[6] + (self.nmtres - 1) as f64;
            self.ncs[1] = self.ncs[1] + self.nmtres - 1;
        }
        nw = 6;
        let mut nc = 0;
        for _i in 1..=self.ncards {
            // tz(j)=card(j+nc); scr(j+nw)=z(j): z is a local equivalenced
            // with tz, so only its first 68 bytes are set.
            let mut text = Vec::with_capacity(68);
            for j in 1..=17 {
                text.extend_from_slice(&self.card[(j + nc - 1) as usize]);
            }
            crate::farr::pack_text(&text, &mut self.sv.recout_z);
            for j in 1..=17 {
                scr[j + nw] = self.sv.recout_z[(j - 1) as usize];
            }
            nw += 17;
            nc += 17;
        }
        let mut nxcc = nxc;
        if nint(self.zain) == 1 && no3 == 0 {
            nxcc -= 1;
        }
        scr[6] = (nxcc + nmtr + i152) as f64;
        c.hdatio(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
        let mut no2 = 1;
        if nxc != 0 {
            nw = (nxc + nmtr + i152) * 6;
            if self.ncoef > 1 {
                nw += (self.nmtres - 1) * 6;
            }
            let mut dict = RV::new(nw.max(0) as usize + 64);
            let mut j = 0;
            let mut imtr = 1usize;
            let mut idone = 0;
            for i in 1..=nxcc {
                if self.mfs[i] == 2 {
                    no2 = 0;
                }
                if self.mfs[i] == 3 || self.mfs[i] == 23 {
                    while imtr as i32 <= nmtr && self.mts[i] >= self.mtr[imtr] {
                        dict[j + 1] = 0.0;
                        dict[j + 2] = 0.0;
                        dict[j + 3] = self.mfs[i] as f64;
                        if self.mfs[i] == 3 && self.mtr[imtr] == 1 && self.awin == 0.0 {
                            dict[j + 4] = 3.0;
                        } else {
                            dict[j + 4] = self.mtr[imtr] as f64;
                        }
                        let np = self.ngo - self.mtrt[imtr] + 1;
                        dict[j + 5] = (3 + (np + 2) / 3) as f64;
                        dict[j + 6] = 0.0;
                        j += 6;
                        imtr += 1;
                    }
                }
                if no3 == 0 && nint(self.zain) == 1 && self.mfs[i] > 3 {
                    no3 = -1;
                    dict[j + 1] = 0.0;
                    dict[j + 2] = 0.0;
                    dict[j + 3] = 3.0;
                    dict[j + 4] = 1.0;
                    dict[j + 5] = 4.0;
                    dict[j + 6] = 0.0;
                    j += 6;
                }
                if self.ncoef > 1 && (i == nxc || self.mfs[i] > 4) && idone == 0 {
                    idone = 1;
                    for k in 1..=self.nmtres as usize {
                        if k == 1 || k > 2 {
                            dict[j + 1] = 0.0;
                            dict[j + 2] = 0.0;
                            dict[j + 3] = 4.0;
                            if k == 1 {
                                dict[j + 4] = self.mmtres[1] as f64;
                                dict[j + 5] = self.mcards[1] as f64;
                            }
                            if k > 2 {
                                dict[j + 4] = self.mmtres[k] as f64;
                                dict[j + 5] = self.mcards[k] as f64;
                            }
                            dict[j + 6] = 0.0;
                            j += 6;
                        }
                    }
                }
                dict[j + 1] = 0.0;
                dict[j + 2] = 0.0;
                dict[j + 3] = self.mfs[i] as f64;
                dict[j + 4] = self.mts[i] as f64;
                dict[j + 5] = self.ncs[i] as f64;
                dict[j + 6] = 0.0;
                j += 6;
                if self.mfs[i] == 2 && i152 != 0 {
                    dict[j + 1] = 0.0;
                    dict[j + 2] = 0.0;
                    dict[j + 3] = 2.0;
                    dict[j + 4] = 152.0;
                    dict[j + 5] = (3 + self.nunr) as f64;
                    dict[j + 6] = 0.0;
                    j += 6;
                }
            }
            nw = nxcc + nmtr + i152;
            if self.ncoef > 1 {
                nw = nw + self.nmtres - 1;
            }
            c.dictio(0, nout, 0, &mut dict.all(), &mut nb, &mut nw);
        }
        c.asend(nout, 0);
        c.afend(nout, 0);

        if no2 <= 0 {
            scr[1] = self.za;
            scr[2] = self.awr;
            scr[3] = 0.0;
            scr[4] = 0.0;
            scr[5] = 1.0;
            scr[6] = 0.0;
            c.mfh = 2;
            c.mth = 151;
            c.contio(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
            scr[2] = 1.0;
            c.contio(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
            scr[1] = emin;
            scr[2] = self.eresh;
            if self.eresr < self.eresh {
                scr[2] = self.eresr;
            }
            scr[5] = 0.0;
            c.contio(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
            scr[1] = self.spin;
            scr[2] = self.ascat;
            c.contio(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
            c.asend(nout, 0);
            if i152 > 0 {
                c.mfh = 2;
                c.mth = 152;
                let sunr = self.sunr.as_mut().unwrap();
                c.contio(0, nout, 0, &mut sunr.all(), &mut nb, &mut nw);
                let mut l = 7;
                c.listio(0, nout, 0, &mut sunr.at(l), &mut nb, &mut nw);
                while nb != 0 {
                    l += nw;
                    c.moreio(0, nout, 0, &mut sunr.at(l), &mut nb, &mut nw);
                }
                c.asend(nout, 0);
            }
            c.afend(nout, 0);
        }

        let ntot = nmtr + 1;
        c.repoz(nscr);
        c.contio(nscr, 0, 0, &mut scr.all(), &mut nb, &mut nw);
        self.za = c.c1h;
        self.awr = c.c2h;
        let mut lis3 = c.l1h;
        let mut lfs = c.l2h;
        let n1l = c.n1h;
        let n2l = c.n2h;
        let mfl = c.mfh;
        let mut mtl = c.mth;
        let mut lbl;
        if nint(self.zain) == 1 && (mfl == 0 || mfl > 3) {
            scr[3] = 0.0;
            scr[4] = 99.0;
            scr[5] = 0.0;
            scr[6] = 0.0;
            c.mfh = 3;
            c.mth = 1;
            c.contio(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
            scr[1] = self.tempr;
            scr[2] = 0.0;
            scr[3] = 0.0;
            scr[4] = 0.0;
            scr[5] = 1.0;
            scr[6] = 2.0;
            scr[7] = 2.0;
            scr[8] = 2.0;
            scr[9] = emin;
            scr[10] = 0.0;
            scr[11] = self.eresh;
            scr[12] = 0.0;
            nw = 12;
            c.tab1io(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
            c.asend(nout, 0);
            c.afend(nout, 0);
            lbl = if mfl == 0 { 272 } else { 270 };
        } else {
            c.mth = 1;
            if self.awin == 0.0 {
                c.mth = 3;
            }
            if c.mfh == 23 {
                c.mth = 501;
            }
            scr[3] = 0.0;
            scr[4] = 99.0;
            c.contio(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
            lbl = 210;
        }
        let mut imtr = 1usize;
        let mut mtd = 0;
        loop {
            match lbl {
                210 => {
                    let mut istart = self.mtrt[imtr];
                    let np = self.ngo - istart + 1;
                    if np <= 0 {
                        let strng = format!(
                            "for mf{} mt{} np:{}-{}+1={}",
                            ffmt::i(c.mfh, 2),
                            ffmt::i(c.mth, 3),
                            ffmt::i(self.ngo, 6),
                            ffmt::i(istart, 6),
                            ffmt::i(np, 6)
                        );
                        c.error("recout", &strng, " ");
                    }
                    scr[1] = self.tempr;
                    scr[2] = 0.0;
                    if c.mth == 18 {
                        scr[2] = sigfig(self.q18, 7, 0);
                    }
                    scr[3] = 0.0;
                    scr[4] = 0.0;
                    scr[5] = 1.0;
                    scr[6] = np as f64;
                    scr[7] = np as f64;
                    scr[8] = 2.0;
                    let mut k = 8;
                    if istart > 1 {
                        let last = istart - 1;
                        let mut i = 0;
                        while i < last {
                            i += 1;
                            c.finda(i, &mut tot, ntot, iold, &mut bufo, NBUF);
                        }
                    }
                    nb = 1;
                    while nb != 0 {
                        let mut iend = self.ngo;
                        if iend - istart >= c.npage / 2 {
                            iend = istart + c.npage / 2 - 1;
                        }
                        let mut j = k - 1;
                        let mut ib = istart - 1;
                        while ib < iend {
                            j += 2;
                            ib += 1;
                            c.finda(ib, &mut tot, ntot, iold, &mut bufo, NBUF);
                            scr[j] = tot[0];
                            scr[j + 1] = sigfig(tot[imtr], 7, 0);
                        }
                        nw = j + 1;
                        if k != 0 {
                            k = 0;
                            c.tab1io(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                        } else {
                            c.moreio(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                        }
                        if nb != 0 {
                            istart = iend + 1;
                        }
                    }
                    c.asend(nout, 0);
                    imtr += 1;
                    if imtr as i32 > nmtr {
                        lbl = 270;
                        continue;
                    }
                    mtd = self.mtr[imtr];
                    lbl = 250;
                }
                250 => {
                    if mtl > mtd {
                        // 260
                        c.mth = mtd;
                        scr[1] = self.za;
                        scr[2] = self.awr;
                        scr[3] = 0.0;
                        scr[4] = 99.0;
                        scr[5] = 0.0;
                        scr[6] = 0.0;
                        c.contio(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                        lbl = 210;
                        continue;
                    }
                    c.mth = mtl;
                    scr[1] = self.za;
                    scr[2] = self.awr;
                    scr[3] = lis3 as f64;
                    scr[4] = lfs as f64;
                    scr[5] = 0.0;
                    scr[6] = 0.0;
                    c.contio(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                    c.tosend(nscr, nout, 0, &mut scr.all());
                    c.contio(nscr, 0, 0, &mut scr.all(), &mut nb, &mut nw);
                    lis3 = c.l1h;
                    lfs = c.l2h;
                    mtl = c.mth;
                    lbl = 250;
                }
                270 => {
                    scr[1] = self.za;
                    scr[2] = self.awr;
                    scr[3] = lis3 as f64;
                    scr[4] = lfs as f64;
                    scr[5] = n1l as f64;
                    scr[6] = n2l as f64;
                    c.mfh = mfl;
                    c.mth = mtl;
                    c.contio(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                    c.tofend(nscr, nout, 0, &mut scr.all());
                    lbl = 272;
                }
                _ => {
                    // 272
                    if self.ncoef > 1 {
                        let nmt = self.nmtres;
                        let na = 1 + self.ncoef * nmt;
                        for imt in 2..=nmt {
                            let mut l = 0;
                            if imt == 2 {
                                l = 1;
                            }
                            if imt > 2 {
                                l = imt;
                            }
                            let mut k = 0;
                            for i in 1..=nrtot {
                                c.finda(i, &mut resl, na, iang, &mut bufl, NBUFL);
                                if resl[(l + nmt) as usize] != zero {
                                    break;
                                }
                                k = i;
                            }
                            c.mfh = 4;
                            c.mth = self.mmtres[imt as usize];
                            if imt == 2 {
                                c.mth = self.mmtres[1];
                            }
                            scr[1] = self.za;
                            scr[2] = self.awr;
                            scr[3] = 0.0;
                            scr[4] = 1.0;
                            scr[5] = 0.0;
                            scr[6] = 0.0;
                            c.contio(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                            scr[1] = 0.0;
                            scr[2] = self.awr;
                            scr[3] = 0.0;
                            scr[4] = 2.0;
                            scr[5] = 0.0;
                            scr[6] = 0.0;
                            c.contio(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                            scr[1] = 0.0;
                            scr[2] = 0.0;
                            scr[3] = 0.0;
                            scr[4] = 0.0;
                            scr[5] = 1.0;
                            scr[6] = (nrtot - k) as f64;
                            scr[7] = (nrtot - k) as f64;
                            scr[8] = 2.0;
                            nw = 8;
                            c.tab2io(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                            if imt == 2 {
                                l = 1;
                            }
                            if imt > 2 {
                                l = imt;
                            }
                            for i in k + 1..=nrtot {
                                c.finda(i, &mut resl, na, iang, &mut bufl, NBUFL);
                                scr[1] = 0.0;
                                scr[2] = resl[0];
                                scr[3] = 0.0;
                                scr[4] = 0.0;
                                scr[5] = 0.0;
                                scr[6] = 0.0;
                                if resl[l as usize] != zero {
                                    for j in 2..=self.ncoef {
                                        let v = resl[(l + (j - 1) * nmt) as usize];
                                        if v != zero {
                                            scr[5 + j] = v;
                                            scr[5] = scr[5] + 1.0;
                                        }
                                    }
                                } else {
                                    scr[5] = 1.0;
                                    scr[7] = 0.0;
                                }
                                c.listio(0, nout, 0, &mut scr.all(), &mut nb, &mut nw);
                            }
                            c.asend(nout, 0);
                        }
                        c.afend(nout, 0);
                    }
                    c.tomend(nscr, nout, 0, &mut scr.all());
                    c.closz(-iold);
                    c.closz(nscr);
                    return;
                }
            }
        }
    }
}

// ======================================================================
// complex probability function (psi-chi broadening)
// ======================================================================

pub fn quickw(ax: f64, y: f64, rew: &mut f64, aimw: &mut f64, ki: &mut i32, tr: &RM, ti: &RM) {
    let pi = crate::util::physics::PI;
    let break1 = 36.0e0f64;
    let break2 = 144.0e0f64;
    let break3 = 10000.0e0f64;
    let c1 = 0.2752551e0f64;
    let c2 = 2.724745e0f64;
    let c3 = 0.5124242e0f64;
    let c4 = 0.05176536e0f64;
    let c5 = 1.1283792e0f64;
    let zero = 0.0f64;
    let rpi = sqrt(pi);
    let mut aki = 1.0;
    if ax < zero {
        aki = -1.0;
    }
    let x = ax.abs();
    let test = x * x + y * y;
    if test < break1 {
        let ii = crate::fmath::int(x * 10.0);
        let jj = crate::fmath::int(y * 10.0);
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
        *rew = a1 * tr[(i, n)] + a2 * tr[(i - 1, j)] + a3 * tr[(i, j)] + a4 * tr[(i + 1, j)]
            + a5 * tr[(i, j + 1)]
            + pq * tr[(i + 1, j + 1)];
        if *ki > 0 {
            *aimw = a1 * ti[(i, n)] + a2 * ti[(i - 1, j)] + a3 * ti[(i, j)] + a4 * ti[(i + 1, j)]
                + a5 * ti[(i, j + 1)]
                + pq * ti[(i + 1, j + 1)];
            *aimw = *aimw * aki;
        }
    } else if test < break2 {
        let a1 = x * x - y * y;
        let a2 = 2.0 * x * y;
        let a3 = a2 * a2;
        let a4 = a1 - c1;
        let a5 = a1 - c2;
        let d1 = c3 / (a4 * a4 + a3);
        let d2 = c4 / (a5 * a5 + a3);
        *rew = d1 * (a2 * x - a4 * y) + d2 * (a2 * x - a5 * y);
        if *ki > 0 {
            *aimw = d1 * (a4 * x + a2 * y) + d2 * (a5 * x + a2 * y);
            *aimw = *aimw * aki;
        }
    } else if test < break3 {
        let a1 = (x * x - y * y) * 2.0;
        let a2 = 4.0 * x * y;
        let a4 = a1 - 1.0;
        let d1 = c5 / (a4 * a4 + a2 * a2);
        *rew = d1 * (a2 * x - a4 * y);
        if *ki > 0 {
            *aimw = d1 * (a4 * x + a2 * y);
            *aimw = *aimw * aki;
        }
    } else {
        let a1 = 1.0 / (rpi * test);
        *rew = y * a1;
        if *ki > 0 {
            *aimw = x * a1 * aki;
        }
    }
}

pub fn wtab(rw: &mut RM, aimw: &mut RM) {
    let tenth = 0.1e0f64;
    let x0 = -tenth;
    let y0 = -tenth;
    let dx = tenth;
    let dy = tenth;
    let nx = 62;
    let ny = 62;
    let mut x = [0.0f64; 63];
    let mut y = [0.0f64; 63];
    y[1] = y0;
    for i in 1..=ny - 1 {
        y[i + 1] = i as f64 * dy + y0;
    }
    x[1] = x0;
    for i in 1..=nx - 1 {
        x[i + 1] = i as f64 * dx + x0;
    }
    for i in 1..=nx {
        for j in 1..=ny {
            let (mut rwt, mut aimt) = (0.0, 0.0);
            w(x[i], y[j], &mut rwt, &mut aimt);
            rw[(i as i32, j as i32)] = rwt;
            aimw[(i as i32, j as i32)] = aimt;
        }
    }
}

/// Complex probability function w(z) by continued fractions/asymptotics.
pub fn w(rez: f64, aim1: f64, rew: &mut f64, aimw: &mut f64) {
    let pi = crate::util::physics::PI;
    let c1 = 1.1283792e0f64;
    let c2 = 1.5e0f64;
    let up = 1.0e15f64;
    let dn = 1.0e-15f64;
    let brk1 = 1.25e0f64;
    let brk2 = 5.0e0f64;
    let brk3 = 1.863636e0f64;
    let brk4 = 4.1e0f64;
    let brk5 = 1.71e0f64;
    let brk6 = 2.89e0f64;
    let brk7 = 1.18e0f64;
    let brk8 = 5.76e0f64;
    let brk9 = 1.5e0f64;
    let eps = 1.0e-7f64;
    let zero = 0.0f64;
    let rpi = sqrt(pi);
    let exp = crate::fmath::exp;
    let cos = crate::fmath::cos;
    let sin = crate::fmath::sin;

    *rew = 0.0;
    *aimw = 0.0;
    let mut aimz = aim1.abs();
    let abrez = rez.abs();
    if abrez + aimz == zero {
        *rew = 1.0;
        *aimw = 0.0;
        return;
    }
    // locals (all uninitialised in Fortran until set)
    let (mut tempc, mut tempd) = (0.0f64, 0.0f64);
    let r2 = rez * rez;
    let ai2 = aimz * aimz;
    let kw;
    let (mut rv, mut ak, mut el, mut a, mut b, mut h, mut g, mut am, mut aak) =
        (0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64);
    let (mut c, mut d, mut tempm, mut temel) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
    let (mut temp1, mut aj, mut ajsig, mut sig2p) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
    let (mut expc, mut exps) = (0.0f64, 0.0f64);
    let (mut ajp, mut temp4) = (0.0f64, 0.0f64);
    let mut lbl;
    let go350 = abrez + brk1 * aimz - brk2 > zero || abrez + brk3 * aimz - brk4 > zero;
    let go340 = !go350 && (r2 + brk5 * ai2 - brk6 < zero || r2 + brk7 * ai2 - brk8 >= zero);
    if !go350 && !go340 && aimz - brk9 >= zero {
        // 350
        if aim1 >= zero {
            kw = 1;
            lbl = 370;
        } else {
            kw = 2;
            aimz = aim1;
            lbl = 420;
        }
    } else if go350 {
        if aim1 >= zero {
            kw = 1;
            lbl = 370;
        } else {
            kw = 2;
            aimz = aim1;
            lbl = 420;
        }
    } else {
        // 340
        kw = 2;
        aimz = aim1;
        lbl = 420;
    }
    loop {
        match lbl {
            370 => {
                rv = 2.0 * (r2 - ai2);
                ak = 4.0 * rez * aimz;
                el = ak;
                h = 0.0;
                b = 0.0;
                a = 0.0;
                tempm = 0.0;
                temel = 0.0;
                g = 1.0;
                c = -c1 * aimz;
                d = c1 * rez;
                am = rv - 1.0;
                aak = 1.0;
                lbl = 380;
            }
            380 => {
                let ajtemp = 2.0 * aak;
                temp4 = (1.0 - ajtemp) * ajtemp;
                ajp = rv - (4.0 * aak + 1.0);
                lbl = 480;
            }
            390 => {
                aak = aak + 1.0;
                let prr = *rew;
                let pii = *aimw;
                let amagn = tempm * tempm + temel * temel;
                *rew = (tempc * tempm + tempd * temel) / amagn;
                *aimw = (tempm * tempd - temel * tempc) / amagn;
                if (*rew - prr).abs() >= eps || (*aimw - pii).abs() >= eps {
                    lbl = 380;
                    continue;
                }
                return;
            }
            420 => {
                temp1 = r2 + ai2;
                let temp2 = 2.0 * temp1 * temp1;
                aj = -(r2 - ai2) / temp2;
                ak = 2.0 * rez * aimz / temp2;
                c = 0.0;
                b = 0.0;
                ajsig = 0.0;
                d = 0.0;
                g = 0.0;
                h = 0.0;
                el = 0.0;
                a = 1.0;
                am = 1.0;
                let sigp = c2;
                let expon = exp(temp2 * aj);
                expc = expon * cos(temp2 * ak);
                exps = -expon * sin(temp2 * ak);
                sig2p = 2.0 * sigp;
                lbl = 430;
            }
            430 => {
                let aj4sig = 4.0 * ajsig;
                let aj4sm1 = aj4sig - 1.0;
                let temp3 = 1.0 / (aj4sm1 * (aj4sig + 3.0));
                let tt4 = sig2p * (2.0 * ajsig - 1.0);
                temp4 = tt4 / (aj4sm1 * (aj4sig + 1.0) * (aj4sig - 3.0) * aj4sm1);
                ajp = aj + temp3;
                lbl = 480;
            }
            440 => {
                ajsig = ajsig + 1.0;
                let temp7 = rpi * (am * am + el * el);
                let reff = (aimz * (c * am + d * el) - rez * (am * d - c * el)) / temp7 / temp1;
                let aimf = (aimz * (am * d - c * el) + rez * (c * am + d * el)) / temp7 / temp1;
                let prr = *rew;
                let pii = *aimw;
                *rew = expc - reff;
                *aimw = exps - aimf;
                if (*rew - prr).abs() >= eps || (*aimw - pii).abs() >= eps {
                    // 470
                    sig2p = 2.0 * ajsig;
                    lbl = 430;
                    continue;
                }
                return;
            }
            _ => {
                // 480
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
                if (tempm.abs() + temel.abs() - up) < zero {
                    // 500
                    if (tempm.abs() + temel.abs() - dn) <= zero {
                        c = up * c;
                        d = up * d;
                        am = up * am;
                        el = up * el;
                        tempc = up * tempc;
                        tempd = up * tempd;
                        tempm = up * tempm;
                        temel = up * temel;
                    }
                } else {
                    c = dn * c;
                    d = dn * d;
                    am = dn * am;
                    el = dn * el;
                    tempc = dn * tempc;
                    tempd = dn * tempd;
                    tempm = dn * tempm;
                    temel = dn * temel;
                }
                // 520: computed goto
                lbl = match kw {
                    1 => 390,
                    2 => 440,
                    _ => return,
                };
            }
        }
    }
}
