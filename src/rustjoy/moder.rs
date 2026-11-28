//! MODER: change the mode of an ENDF-type tape (text <-> blocked binary),
//! or select materials from one or more tapes onto a new tape. Ported from
//! NJOY2016 `moder.f90` (module modem).

use crate::endf::Ctx;
use crate::farr::RV;
use crate::ffmt;
use crate::fio::{Kind, Val};
use crate::fmath::nint;

/// Size of the record buffer `a` (nbuf in the Fortran).
const NBUF: usize = 1000000;

/// MODER module variables.
#[derive(Default)]
pub struct Moder {
    ntw: i32,
    ng: i32,
}

fn geti(v: &Option<Val>, d: i32) -> i32 {
    match v {
        Some(Val::Int(i)) => *i,
        Some(Val::Real(x)) => *x as i32,
        _ => d,
    }
}

/// `read(nsysi,*) nin,nout`-style integer card.
fn read_ints(c: &mut Ctx, n: usize, defaults: &[i32]) -> Vec<i32> {
    let kinds: Vec<Kind> = (0..n).map(|_| Kind::Int).collect();
    let v = c.io.input.read(&kinds);
    (0..n).map(|i| geti(&v[i], defaults[i])).collect()
}

impl Moder {
    pub fn new() -> Self {
        Moder::default()
    }

    // ================================================================ moder
    pub fn moder(&mut self, c: &mut Ctx) {
        let (mut nb, mut nw) = (0i32, 0i32);
        let nscr = 0;
        let mut ninl = 0;
        let mut matl = -2;
        let mut a = RV::new(NBUF);
        let mut z = RV::new(20);
        let mut matd = 0;
        let mut mend = 0;
        let mut inout = 1;

        //--read user input, initialize, and write output header.
        let time = c.timer();
        c.io.so(&format!(
            "\n moder...change the mode of an endf tape or njoy output tape{}{}s",
            " ".repeat(9),
            ffmt::f(time, 8, 1)
        ));
        c.io.se(&format!("\n moder...{}{}s", " ".repeat(60), ffmt::f(time, 8, 1)));
        // nin and nout are locals of the Fortran routine: unset if missing
        let v = read_ints(c, 2, &[0, 0]);
        let mut nin = v[0];
        let nout = v[1];
        c.openz(nout, 1);
        if nin.abs() >= 20 {
            c.openz(nin, 0);
        }
        c.nsh = 0;
        c.nsp = 0;
        c.nsc = 0;
        c.repoz(nout);
        let mut no = nout;
        let mut loop_ = 0;
        if nin.abs() >= 1 && nin.abs() <= 19 {
            loop_ = 1;
        }
        let mut lbl;
        if loop_ == 1 {
            // 110
            c.io.so(&format!("\n put materials from various tapes on output tape {}", ffmt::i(nout, 3)));
            let v = c.io.input.read(&[Kind::Char(105)]);
            let strng = match &v[0] {
                Some(Val::Char(s)) => s.clone(),
                _ => " ".repeat(105),
            };
            // read(strng,'(16a4,a2)') zc(1..17): 66 characters, the a2 item
            // is padded with 2 blanks in its character(4) variable
            let mut text: Vec<u8> = strng.as_bytes().to_vec();
            text.resize(105, b' ');
            let mut t68 = text[..66].to_vec();
            t68.push(b' ');
            t68.push(b' ');
            let nz = 17;
            c.io.so(&format!(
                "\n tape id for nout\n ----------------\n {}",
                String::from_utf8_lossy(&text[..66])
            ));
            let mut zr = [0.0f64; 20];
            crate::farr::pack_text(&t68, &mut zr);
            for i in 1..=20 {
                z[i] = zr[(i - 1) as usize];
            }
            c.math = 1;
            c.mfh = 0;
            c.mth = 0;
            inout = 1;
            if nin == 2 {
                inout = 2;
            }
            if nin == 3 {
                inout = 3;
            }
            if inout == 1 {
                c.io.so("\n processing endf or pendf tape.");
                c.io.se("\n processing endf or pendf tape.");
            } else if inout == 2 {
                c.io.so("\n processing gendf tape");
                c.io.se("\n processing gendf tape");
            } else if inout == 3 {
                c.io.so("\n processing covariance tape");
                c.io.se("\n processing covariance tape");
            }
            let mut nzz = nz;
            c.tpidio(0, nout, nscr, &mut z.all(), &mut nb, &mut nzz);
            no = 0;
            c.io.so(&format!("\n{}nin{}matd\n{}---{}----", " ".repeat(10), " ".repeat(10), " ".repeat(10), " ".repeat(10)));
            lbl = 130;
        } else {
            c.io.so(&format!(
                "\n input unit (+ for coded, - for bb) ... {}\n output unit (+ for coded, - for bb) .. {}",
                ffmt::i(nin, 10),
                ffmt::i(nout, 10)
            ));
            lbl = 140;
        }

        loop {
            match lbl {
                130 => {
                    nin = 0;
                    let v = read_ints(c, 2, &[0, matd]);
                    nin = v[0];
                    matd = v[1];
                    if nin == 0 {
                        lbl = 1000;
                        continue;
                    }
                    c.io.so(&format!("{}{}{}{}", " ".repeat(10), ffmt::i(nin, 3), " ".repeat(10), ffmt::i(matd, 4)));
                    if inout == 1 && matd <= matl {
                        c.error("moder", "endf materials must be in ascending order.", " ");
                    }
                    if nin != ninl {
                        c.closz(ninl);
                    }
                    if nin != ninl {
                        c.openz(nin, 0);
                    }
                    if nin == ninl && inout >= 2 {
                        lbl = 140;
                        continue;
                    }
                    if nin == ninl && matd > matl {
                        lbl = 180;
                        continue;
                    }
                    lbl = 140;
                }
                140 => {
                    c.repoz(nin);
                    c.repoz(nscr);
                    c.tpidio(nin, no, nscr, &mut a.all(), &mut nb, &mut nw);
                    if loop_ > 0 {
                        lbl = 180;
                        continue;
                    }
                    let mut rb = [0.0f64; 17];
                    for i in 1..=17 {
                        rb[(i - 1) as usize] = a[i];
                    }
                    let hb = crate::farr::unpack_text(&rb);
                    // write(strng,'(16a4,a2)') hb: 64 + 2 characters
                    if hb[..66].iter().any(|&b| b != b' ') {
                        c.io.so(&format!(
                            "\n tape label\n {}--\n  {}",
                            "----".repeat(9),
                            String::from_utf8_lossy(&hb[..68])
                        ));
                    } else {
                        c.io.so(" tape id is blank");
                    }
                    lbl = 205;
                }
                180 => {
                    c.contio(nin, 0, 0, &mut a.all(), &mut nb, &mut nw);
                    if c.n1h == -1 {
                        lbl = 190;
                        continue;
                    }
                    c.findf(matd, 0, 0, nin);
                    ninl = nin;
                    matl = matd;
                    lbl = 205;
                }
                190 => {
                    if c.math != matd {
                        // 200
                        c.tomend(nin, 0, 0, &mut a.all());
                        c.contio(nin, 0, 0, &mut a.all(), &mut nb, &mut nw);
                        if c.math != -1 {
                            continue;
                        }
                        let strng = format!("mat={} not found on gendf tape", ffmt::i(matd, 4));
                        c.mess("moder", &strng, " ");
                        lbl = 130;
                        continue;
                    }
                    ninl = nin;
                    matl = matd;
                    lbl = 700;
                }
                205 => {
                    c.iverf = 0;
                    lbl = 210;
                }
                210 => {
                    //--process endf or pendf files
                    c.contio(nin, 0, 0, &mut a.all(), &mut nb, &mut nw);
                    if c.mfh == 1 && c.mth == 451 {
                        c.nsh = 1;
                        if c.n1h == -1 {
                            lbl = 700;
                            continue;
                        }
                        if c.iverf == 0 {
                            c.contio(nin, 0, 0, &mut a.at(7), &mut nb, &mut nw);
                            if c.n1h != 0 {
                                c.iverf = 4;
                            } else if c.n2h == 0 {
                                c.iverf = 5;
                            } else {
                                c.iverf = 6;
                            }
                            c.io.so(&format!("\n using endf-{} format", ffmt::i(c.iverf, 1)));
                            c.skiprz(nin, -1);
                        }
                    }
                    // 215
                    if !(loop_ == 0 || c.math == 0) && c.math != matd {
                        lbl = 130;
                        continue;
                    }
                    // 220
                    c.contio(0, nout, nscr, &mut a.all(), &mut nb, &mut nw);
                    if loop_ > 0 && c.math == 0 {
                        lbl = 130;
                        continue;
                    }
                    if c.math == -1 {
                        lbl = 1010;
                        continue;
                    }
                    if c.mfh > 0 {
                        lbl = 230;
                        continue;
                    }
                    if c.math > 0 {
                        continue;
                    }
                    c.nsh = 1;
                    c.nsp = 1;
                    c.nsc = 1;
                }
                230 => {
                    if c.mfh == 1 {
                        // 270
                        if c.mth == 451 && (c.n1h == -11 || c.n1h == -12 || c.n1h == -14) {
                            lbl = 800;
                            continue;
                        }
                        if loop_ != 0 {
                            if inout == 2 {
                                c.error("moder", "this material is not a gendf material.", " ");
                            }
                            if inout == 3 {
                                c.error("moder", "input is not an errorr output tape.", " ");
                            }
                        }
                        file1(c, nin, nout, nscr, &mut a);
                    } else {
                        match c.mfh {
                            2 => file2(c, nin, nout, nscr, &mut a),
                            3 | 23 | 27 => file3(c, nin, nout, nscr, &mut a),
                            4 | 24 => file4(c, nin, nout, nscr, &mut a),
                            5 | 25 => file5(c, nin, nout, nscr, &mut a),
                            6 | 26 => file6(c, nin, nout, nscr, &mut a),
                            7 => file7(c, nin, nout, nscr, &mut a),
                            8 => file8(c, nin, nout, nscr, &mut a),
                            9 | 10 => file9(c, nin, nout, nscr, &mut a),
                            12 | 13 | 30 => file1x(c, nin, nout, nscr, &mut a),
                            14 => file14(c, nin, nout, nscr, &mut a),
                            15 => file15(c, nin, nout, nscr, &mut a),
                            28 => file28(c, nin, nout, nscr, &mut a),
                            31 | 33 => file3x(c, nin, nout, nscr, &mut a),
                            32 => file32(c, nin, nout, nscr, &mut a),
                            34 => file34(c, nin, nout, nscr, &mut a),
                            35 => file35(c, nin, nout, nscr, &mut a),
                            40 => file40(c, nin, nout, nscr, &mut a),
                            _ => {
                                let strng = format!("conversion not coded for mf={}", ffmt::i(c.mfh, 3));
                                c.error("moder", &strng, " ");
                            }
                        }
                    }
                    // 410
                    c.contio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
                    if c.mth != 0 {
                        let mut strng = format!(
                            " {}{}{}{}{}",
                            ffmt::i(c.math, 6),
                            ffmt::i(c.mfh, 6),
                            ffmt::i(c.mth, 6),
                            ffmt::i(nb, 6),
                            ffmt::i(nw, 6)
                        );
                        for i in 1..=6 {
                            strng.push_str(&ffmt::e1p(a[i], 12, 5));
                        }
                        c.error("moder", "should have found send card", &strng);
                    }
                    lbl = 210;
                }
                700 => {
                    //--process gendf tapes
                    if loop_ != 0 {
                        if inout == 1 {
                            c.error("moder", "input is not an endf or pendf tape.", " ");
                        }
                        if inout == 3 {
                            c.error("moder", "input is not an errorr output tape.", " ");
                        }
                    }
                    self.ntw = c.n2h;
                    lbl = 705;
                }
                705 => {
                    if c.n1h != -1 {
                        let strng = format!("mat={} is not a gendf material", ffmt::i(c.math, 4));
                        c.error("moder", &strng, " ");
                    }
                    self.glstio(c, nin, nout, nscr, &mut a, &mut nb, &mut nw);
                    mend = 0;
                    lbl = 740;
                }
                710 => {
                    c.listio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
                    let ig = c.n2h;
                    while nb != 0 {
                        c.moreio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
                    }
                    if c.mfh != 1 && ig != self.ng {
                        continue;
                    }
                    lbl = 740;
                }
                740 => {
                    c.contio(nin, 0, 0, &mut a.all(), &mut nb, &mut nw);
                    if loop_ > 0 && c.math != matd && c.math != 0 {
                        lbl = 130;
                        continue;
                    }
                    if c.math == -1 {
                        lbl = 1000;
                        continue;
                    }
                    if c.math > 0 {
                        // 750
                        if mend > 0 {
                            lbl = 705;
                            continue;
                        }
                        c.contio(0, nout, nscr, &mut a.all(), &mut nb, &mut nw);
                        if c.math <= 0 || c.mfh == 0 || c.mth == 0 {
                            continue;
                        }
                        lbl = 710;
                        continue;
                    }
                    mend = 1;
                    c.amend(nout, nscr);
                    c.nsp = 1;
                    c.nsc = 1;
                }
                800 | 805 => {
                    //--process covariance tape
                    if lbl == 800 && loop_ != 0 {
                        if inout == 1 {
                            c.error("moder", "input is not an endf tape.", " ");
                        }
                        if inout == 2 {
                            c.error("moder", "input is not a gendf tape.", " ");
                        }
                    }
                    // 805
                    c.listio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
                    self.ng = c.l1h;
                    while nb != 0 {
                        c.moreio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
                    }
                    lbl = 820;
                }
                820 => {
                    c.contio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
                    if loop_ > 0 && c.math == 0 {
                        lbl = 130;
                        continue;
                    }
                    if c.math == -1 {
                        lbl = 1010;
                        continue;
                    }
                    if c.mth == 0 {
                        continue;
                    }
                    if c.mfh == 1 {
                        lbl = 805;
                        continue;
                    }
                    if c.mfh == 3 || c.mfh == 5 {
                        //--covariance tape mf3: no head cards; backspace
                        //--and read the list record
                        c.skiprz(nin, -1);
                        c.skiprz(nout, -1);
                        if nout > 0 {
                            c.nsh -= 1;
                        }
                        c.listio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
                        while nb != 0 {
                            c.moreio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
                        }
                    } else if c.mfh == 33 || c.mfh == 40 {
                        //--covariance tapes mf33 or mf40
                        let nk = c.n2h;
                        for _ik in 1..=nk {
                            c.contio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
                            let mut ig = 0;
                            while ig < self.ng {
                                c.listio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
                                ig = c.n2h;
                                while nb != 0 {
                                    c.moreio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
                                }
                            }
                        }
                    } else if c.mfh == 34 || c.mfh == 35 {
                        //--covariance tape mf34 or mf35
                        c.contio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
                        self.ng = c.n2h;
                        let mut ig = 0;
                        while ig < self.ng {
                            c.listio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
                            ig = c.n2h;
                            while nb != 0 {
                                c.moreio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
                            }
                        }
                    } else {
                        let strng = format!("illegal covariance mf={}", ffmt::i(c.mfh, 3));
                        c.error("moder", &strng, " ");
                    }
                }
                1000 | 1010 => {
                    //--moder complete
                    if lbl == 1000 {
                        c.atend(nout, nscr);
                    }
                    c.repoz(nout);
                    c.repoz(nin);
                    c.closz(nin);
                    c.closz(nout);
                    let time = c.timer();
                    c.io.so(&format!("{}{}s\n {}*******", " ".repeat(69), ffmt::f(time, 8, 1), "**********".repeat(7)));
                    return;
                }
                _ => unreachable!(),
            }
        }
    }

    // =============================================================== glstio
    /// Convert a GENDF File 1 MT451 list record to the opposite format
    /// (text GOUT tapes have no titles).
    fn glstio(&mut self, c: &mut Ctx, nin: i32, nout: i32, nscr: i32, a: &mut RV, nb: &mut i32, nw: &mut i32) {
        //--read in rest of header record.
        //--binary outputs can also be written.
        let mut n1 = nout;
        if n1 > 0 {
            n1 = 0;
        }
        let mut n2 = nscr;
        if n2 > 0 {
            n2 = 0;
        }
        c.contio(0, n1, n2, &mut a.at(1), nb, nw);
        self.ntw = c.n2h;
        let mut l = 7;
        c.listio(nin, n1, n2, &mut a.at(l), nb, nw);
        self.ng = c.l1h;
        while *nb != 0 {
            l += *nw;
            c.moreio(nin, n1, n2, &mut a.at(l), nb, nw);
        }
        *nw = l + *nw - 7;
        if nout <= 0 && nscr <= 0 {
            return;
        }

        //--remove title from binary input and write bcd outputs, if any.
        n1 = nout;
        if n1 < 0 {
            n1 = 0;
        }
        n2 = nscr;
        if n2 < 0 {
            n2 = 0;
        }
        a[6] = 1.0;
        let (mut nbt, mut nwt) = (0i32, 0i32);
        c.contio(0, n1, n2, &mut a.at(1), &mut nbt, &mut nwt);
        let noff = self.ntw - 1;
        *nw -= noff;
        for i in 2..=*nw {
            let l = 12 + i;
            a[l] = a[l + noff];
        }
        a[11] -= noff as f64;
        a[13] = 0.0;
        let mut l = 7;
        c.listio(0, n1, n2, &mut a.at(l), nb, nw);
        while *nb != 0 {
            l += *nw;
            c.moreio(0, n1, n2, &mut a.at(l), nb, nw);
        }
    }
}

// ---------------------------------------------------------------------
// conversion of each file (the records are read from nin and written to
// nout/nscr in their modes)
// ---------------------------------------------------------------------

/// listio followed by its moreio continuation records.
fn list(c: &mut Ctx, nin: i32, nout: i32, nscr: i32, a: &mut RV) {
    let (mut nb, mut nw) = (0i32, 0i32);
    c.listio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
    while nb != 0 {
        c.moreio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
    }
}

/// tab1io followed by its moreio continuation records.
fn tab1(c: &mut Ctx, nin: i32, nout: i32, nscr: i32, a: &mut RV) {
    let (mut nb, mut nw) = (0i32, 0i32);
    c.tab1io(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
    while nb != 0 {
        c.moreio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
    }
}

fn cont(c: &mut Ctx, nin: i32, nout: i32, nscr: i32, a: &mut RV) {
    let (mut nb, mut nw) = (0i32, 0i32);
    c.contio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
}

fn tab2(c: &mut Ctx, nin: i32, nout: i32, nscr: i32, a: &mut RV) {
    let (mut nb, mut nw) = (0i32, 0i32);
    c.tab2io(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
}

/// File 1: general information.
fn file1(c: &mut Ctx, nin: i32, nout: i32, nscr: i32, a: &mut RV) {
    let (mut nb, mut nw) = (0i32, 0i32);
    let mth = c.mth;
    if mth == 451 {
        //--hollerith descriptive data and tape dictionary.
        let mut nx = 0;
        if c.iverf <= 4 {
            nx = c.n2h;
        }
        if c.iverf >= 5 {
            cont(c, nin, nout, nscr, a);
        }
        if c.iverf >= 6 {
            cont(c, nin, nout, nscr, a);
        }
        c.hdatio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
        if c.iverf >= 5 {
            nx = c.n2h;
        }
        while nb != 0 {
            c.moreio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
        }
        nw = nx;
        c.dictio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
    } else if mth == 452 {
        //--total fission yield (nu-bar).
        if c.l2h != 2 {
            list(c, nin, nout, nscr, a);
        } else {
            tab1(c, nin, nout, nscr, a);
        }
    } else if mth == 453 {
        //--induced reaction branching ratios.
        let ns = c.n1h;
        for _i in 1..=ns {
            list(c, nin, nout, nscr, a);
            let npr = c.n2h;
            for _j in 1..=npr {
                list(c, nin, nout, nscr, a);
            }
        }
    } else if mth == 454 {
        //--fission product yield data.
        let lep1 = c.l1h;
        for _i in 1..=lep1 {
            list(c, nin, nout, nscr, a);
        }
    } else if mth == 455 {
        //--delayed fission yield (nu-bar).
        let lnd = c.l2h;
        list(c, nin, nout, nscr, a);
        if lnd != 2 {
            list(c, nin, nout, nscr, a);
        } else {
            tab1(c, nin, nout, nscr, a);
        }
    } else if mth == 456 {
        //--prompt fission yield (nu-bar).
        let lnp = c.l2h;
        if lnp != 2 {
            list(c, nin, nout, nscr, a);
        } else {
            tab1(c, nin, nout, nscr, a);
        }
    } else if mth == 457 {
        //--radioactive decay data.
        let ns457 = c.n2h;
        for _i in 1..=2 {
            list(c, nin, nout, nscr, a);
        }
        let mut ns = 0;
        while ns < ns457 {
            list(c, nin, nout, nscr, a);
            ns += 1;
        }
    } else if mth == 458 {
        //--components of energy release due to fission.
        let lfc = c.l2h;
        let nfc = c.n2h;
        c.listio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
        let nply = c.l2h;
        while nb != 0 {
            c.moreio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
        }
        if lfc == 1 {
            if nply != 0 {
                c.error("file1", "bad NPLY in mt=458.", " ");
            }
            for _i in 1..=nfc {
                tab1(c, nin, nout, nscr, a);
            }
        } else if lfc != 0 {
            c.error("file1", "bad LFC in mt=458.", " ");
        }
    } else if mth == 460 {
        //--beta-delayed photon spectra
        let lo = c.l1h;
        if lo == 1 {
            let ng460 = c.n1h;
            for _ng in 1..=ng460 {
                tab1(c, nin, nout, nscr, a);
            }
        } else if lo == 2 {
            c.listio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
        } else {
            c.error("file1", "bad LO in mt=460.", " ");
        }
    } else {
        //--illegal mt for file 1
        c.error("file1", "illegal mt.", " ");
    }
}

/// File 2: resonance parameters.
fn file2(c: &mut Ctx, nin: i32, nout: i32, nscr: i32, a: &mut RV) {
    if c.mth == 151 {
        file2a(c, nin, nout, nscr, a);
    } else if c.mth == 152 || c.mth == 153 {
        //--special point-unresolved and probability-table formats
        list(c, nin, nout, nscr, a);
    } else {
        c.error("file2", "illegal mt.", " ");
    }
}

fn file2a(c: &mut Ctx, nin: i32, nout: i32, nscr: i32, a: &mut RV) {
    let nis = c.n1h;
    for _i in 1..=nis {
        cont(c, nin, nout, nscr, a);
        let ner = c.n1h;
        let lfw = c.l2h;
        for _n in 1..=ner {
            cont(c, nin, nout, nscr, a);
            let lru = c.l1h;
            let lrf = c.l2h;
            let nro = c.n1h;
            if nro != 0 {
                let (mut nb, mut nw) = (0i32, 0i32);
                c.tab1io(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
            }
            if lru == 2 && lrf == 1 && lfw == 1 {
                list(c, nin, nout, nscr, a);
                let nls = c.n2h;
                for _l in 1..=nls {
                    cont(c, nin, nout, nscr, a);
                    let njs = c.n1h;
                    for _j in 1..=njs {
                        list(c, nin, nout, nscr, a);
                    }
                }
            } else {
                cont(c, nin, nout, nscr, a);
                if lrf != 7 {
                    let nls = c.n1h;
                    if nls != 0 {
                        for _l in 1..=nls {
                            if lru == 2 && lrf == 2 {
                                cont(c, nin, nout, nscr, a);
                                let njs = c.n1h;
                                for _j in 1..=njs {
                                    list(c, nin, nout, nscr, a);
                                }
                            } else {
                                list(c, nin, nout, nscr, a);
                                if (lru == 1 && lrf == 4) || (lru == 2 && lfw != 0) {
                                    cont(c, nin, nout, nscr, a);
                                    let njs = c.n1h;
                                    for _j in 1..=njs {
                                        list(c, nin, nout, nscr, a);
                                    }
                                }
                            }
                        }
                    }
                } else {
                    let njs = c.n1h;
                    list(c, nin, nout, nscr, a);
                    for _j in 1..=njs {
                        list(c, nin, nout, nscr, a);
                        let kbk = c.l1h;
                        let kps = c.l2h;
                        list(c, nin, nout, nscr, a);
                        if kbk > 0 {
                            for _l in 1..=kbk {
                                cont(c, nin, nout, nscr, a);
                                let lbk = c.l2h;
                                if lbk == 1 {
                                    tab1(c, nin, nout, nscr, a);
                                    tab1(c, nin, nout, nscr, a);
                                } else if lbk == 2 || lbk == 3 {
                                    list(c, nin, nout, nscr, a);
                                }
                            }
                        }
                        if kps == 1 {
                            let (mut nb, mut nw) = (0i32, 0i32);
                            c.listio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
                            let lps = c.n1h;
                            if lps == 1 {
                                tab1(c, nin, nout, nscr, a);
                            }
                        }
                    }
                }
            }
        }
    }
}

/// File 3 (and 23, 27): cross sections.
fn file3(c: &mut Ctx, nin: i32, nout: i32, nscr: i32, a: &mut RV) {
    let (mut nb, mut nw) = (0i32, 0i32);
    c.tab1io(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
    let mut lt = c.l1h;
    while nb != 0 {
        c.moreio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
    }
    while lt > 0 {
        list(c, nin, nout, nscr, a);
        lt -= 1;
    }
}

/// File 4 (and 24): angular distributions.
fn file4(c: &mut Ctx, nin: i32, nout: i32, nscr: i32, a: &mut RV) {
    let (mut nb, mut nw) = (0i32, 0i32);
    let lvt = c.l1h;
    let ltt = c.l2h;
    let mut nk = c.n1h;
    if nk == 0 {
        nk = 1;
    }
    for _ik in 1..=nk {
        if lvt == 0 {
            cont(c, nin, nout, nscr, a);
        }
        if lvt != 0 || c.l1h != 1 {
            if lvt == 1 {
                list(c, nin, nout, nscr, a);
            }
            tab2(c, nin, nout, nscr, a);
            let ne = c.n2h;
            if ltt != 2 {
                for _ie in 1..=ne {
                    c.listio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
                    let mut lt = c.l1h;
                    while nb != 0 {
                        c.moreio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
                    }
                    while lt > 0 {
                        list(c, nin, nout, nscr, a);
                        lt -= 1;
                    }
                }
            } else {
                for _ie in 1..=ne {
                    c.tab1io(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
                    let mut lt = c.l1h;
                    while nb != 0 {
                        c.moreio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
                    }
                    while lt > 0 {
                        list(c, nin, nout, nscr, a);
                        lt -= 1;
                    }
                }
            }
        }
        if ltt >= 3 {
            tab2(c, nin, nout, nscr, a);
            let ne = c.n2h;
            for _ie in 1..=ne {
                tab1(c, nin, nout, nscr, a);
            }
        }
    }
}

/// File 5 (and 25): secondary energy distributions.
fn file5(c: &mut Ctx, nin: i32, nout: i32, nscr: i32, a: &mut RV) {
    let nk = c.n1h;
    for _k in 1..=nk {
        tab1(c, nin, nout, nscr, a);
        let lf = c.l2h;
        if lf == 1 {
            tab2(c, nin, nout, nscr, a);
            let ne = c.n2h;
            for _ie in 1..=ne {
                tab1(c, nin, nout, nscr, a);
            }
        } else if lf == 5 || lf == 11 {
            tab1(c, nin, nout, nscr, a);
            tab1(c, nin, nout, nscr, a);
        } else if lf == 7 || lf == 9 || lf == 12 {
            tab1(c, nin, nout, nscr, a);
        } else if lf == 10 {
            list(c, nin, nout, nscr, a);
        } else if lf != 3 {
            c.error("file5", "illegal lf.", " ");
        }
    }
}

/// File 6 (and 26): particle energy-angle distributions.
fn file6(c: &mut Ctx, nin: i32, nout: i32, nscr: i32, a: &mut RV) {
    let (mut nb, mut nw) = (0i32, 0i32);
    let nk = c.n1h;
    for _ik in 1..=nk {
        c.tab1io(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
        let law = c.l2h;
        while nb != 0 {
            c.moreio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
        }
        if law == 6 {
            cont(c, nin, nout, nscr, a);
        } else if law == 1 || law == 2 || law == 5 {
            tab2(c, nin, nout, nscr, a);
            let ne = c.n2h;
            for _ie in 1..=ne {
                list(c, nin, nout, nscr, a);
            }
        } else if law == 7 {
            tab2(c, nin, nout, nscr, a);
            let ne = c.n2h;
            for _ie in 1..=ne {
                tab2(c, nin, nout, nscr, a);
                let nmu = c.n2h;
                for _imu in 1..=nmu {
                    tab1(c, nin, nout, nscr, a);
                }
            }
        } else if law < 0 {
            if c.mth != 18 {
                let strng = format!("illegal endf6 law for mt={}", ffmt::i(c.mth, 3));
                c.error("file6", &strng, " ");
            }
        } else if law != 0 && law != 3 && law != 4 {
            c.error("file6", "illegal endf6 law.", " ");
        }
    }
}

/// File 7: thermal neutron scattering law.
fn file7(c: &mut Ctx, nin: i32, nout: i32, nscr: i32, a: &mut RV) {
    let (mut nb, mut nw) = (0i32, 0i32);
    let zero = 0.0f64;
    if c.mth == 4 {
        //--incoherent inelastic (mt=4)
        c.listio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
        let b = a[7];
        let mut it2 = 0;
        let mut ia1 = [0i32; 4];
        if c.iverf >= 6 {
            it2 = c.n2h;
            if it2 > 0 && it2 <= 3 {
                for i in 1..=it2 {
                    let ii = 6 * i + 7;
                    ia1[i as usize] = nint(a[ii]);
                }
            } else if it2 > 3 {
                c.error("file7", "bad NS in mt=4.", " ");
            }
        }
        while nb != 0 {
            c.moreio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
        }
        if b == zero {
            return;
        }
        tab2(c, nin, nout, nscr, a);
        let n = c.n2h;
        for _i in 1..=n {
            c.tab1io(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
            let mut lt = c.l1h;
            while nb != 0 {
                c.moreio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
            }
            while lt > 0 {
                list(c, nin, nout, nscr, a);
                lt -= 1;
            }
        }

        //--effective temperatures
        if c.iverf >= 6 {
            tab1(c, nin, nout, nscr, a);
            if it2 != 0 {
                for i in 1..=it2 {
                    if ia1[i as usize] == 0 {
                        tab1(c, nin, nout, nscr, a);
                    }
                }
            }
        }
    } else if c.mth == 2 && c.iverf >= 6 {
        //--elastic (mt=2)
        let lthr = c.l1h;

        //--coherent (lthr=1) or mixed (lthr=3)
        if lthr == 1 || lthr == 3 {
            c.tab1io(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
            let mut lt = c.l1h;
            while nb != 0 {
                c.moreio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
            }
            while lt > 0 {
                list(c, nin, nout, nscr, a);
                lt -= 1;
            }
        }

        //--incoherent (lthr=2) or mixed (lthr=3)
        if lthr == 2 || lthr == 3 {
            tab1(c, nin, nout, nscr, a);
        }

        //--illegal lthr for mt=2
        if !(1..=3).contains(&lthr) {
            let strng = format!("illegal value of lthr={}", ffmt::i(lthr, 4));
            c.error("file7", &strng, " ");
        }
    } else if c.mth == 451 {
        //--general information (mt=451)
        let n = c.l1h;
        for _i in 1..=n {
            list(c, nin, nout, nscr, a);
        }
    } else {
        let strng = format!("illegal mt={}", ffmt::i(c.mth, 3));
        c.error("file7", &strng, " ");
    }
}

/// File 8: radioactive decay and fission product yield data.
fn file8(c: &mut Ctx, nin: i32, nout: i32, nscr: i32, a: &mut RV) {
    let (mut nb, mut nw) = (0i32, 0i32);
    if c.mth == 454 || c.mth == 459 {
        //--fission product yield data.
        let ne = c.l1h;
        for _ie in 1..=ne {
            list(c, nin, nout, nscr, a);
        }
    } else if c.mth == 457 {
        //--radioactive decay data.
        let ns = c.n2h;
        c.listio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
        list(c, nin, nout, nscr, a);
        if ns > 0 {
            for _is in 1..=ns {
                c.listio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
                let lcon = c.l1h;
                if lcon != 1 {
                    let ne = c.n2h;
                    for _ie in 1..=ne {
                        c.listio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
                    }
                }
                if lcon != 0 {
                    tab1(c, nin, nout, nscr, a);
                }
            }
        }
    } else {
        //--radioactive nuclide production.
        let ns = c.n1h;
        let no = c.n2h;
        if no != 1 {
            for _lfs in 1..=ns {
                list(c, nin, nout, nscr, a);
            }
        } else {
            for _lfs in 1..=ns {
                cont(c, nin, nout, nscr, a);
            }
        }
    }
}

/// File 9 or 10: yields or cross sections for radioactive nuclide
/// production.
fn file9(c: &mut Ctx, nin: i32, nout: i32, nscr: i32, a: &mut RV) {
    let ns = c.n1h;
    for _lfs in 1..=ns {
        tab1(c, nin, nout, nscr, a);
    }
}

/// File 12, 13 (and 30): photon yields or photon cross sections.
fn file1x(c: &mut Ctx, nin: i32, nout: i32, nscr: i32, a: &mut RV) {
    let lo = c.l1h;
    if c.mfh == 12 && c.mth == 460 {
        //--special path for beta delayed gammas (mt=460)
        if lo == 1 {
            let ng460 = c.n1h;
            for _ng in 1..=ng460 + 1 {
                tab1(c, nin, nout, nscr, a);
            }
        } else if lo == 2 {
            list(c, nin, nout, nscr, a);
        } else {
            c.error("file12", "bad LO in mt=460.", " ");
        }
    } else if lo != 2 {
        //--normal path for other mf1x formats
        let nk = c.n1h;
        let mut nkp = nk + 1;
        if nk == 1 {
            nkp -= 1;
        }
        for _k in 1..=nkp {
            tab1(c, nin, nout, nscr, a);
        }
    } else {
        list(c, nin, nout, nscr, a);
    }
}

/// File 14: photon angular distributions.
fn file14(c: &mut Ctx, nin: i32, nout: i32, nscr: i32, a: &mut RV) {
    let li = c.l1h;
    if li == 1 {
        return;
    }
    let ltt = c.l2h;
    let nk = c.n1h;
    let ni = c.n2h;
    for i in 1..=nk {
        if i <= ni {
            cont(c, nin, nout, nscr, a);
        } else {
            tab2(c, nin, nout, nscr, a);
            let ne = c.n2h;
            if ltt != 2 {
                for _j in 1..=ne {
                    list(c, nin, nout, nscr, a);
                }
            } else {
                for _j in 1..=ne {
                    tab1(c, nin, nout, nscr, a);
                }
            }
        }
    }
}

/// File 15: photon secondary energy distributions.
fn file15(c: &mut Ctx, nin: i32, nout: i32, nscr: i32, a: &mut RV) {
    tab1(c, nin, nout, nscr, a);
    let lf = c.l2h;
    if lf == 1 {
        //--tabulated distribution.
        tab2(c, nin, nout, nscr, a);
        let ne = c.n2h;
        for _i in 1..=ne {
            tab1(c, nin, nout, nscr, a);
        }
    } else if lf == 2 {
        //--special lanl format.
        tab2(c, nin, nout, nscr, a);
        let ne = c.l2h;
        for _i in 1..=ne {
            list(c, nin, nout, nscr, a);
        }
    } else {
        c.error("file15", "illegal lf.", " ");
    }
}

/// File 28: atomic relaxation data.
fn file28(c: &mut Ctx, nin: i32, nout: i32, nscr: i32, a: &mut RV) {
    let nss = c.n1h;
    for _i in 1..=nss {
        list(c, nin, nout, nscr, a);
    }
}

/// File 31 or 33: cross section covariances.
fn file3x(c: &mut Ctx, nin: i32, nout: i32, nscr: i32, a: &mut RV) {
    let mut nl = c.l2h;
    if c.iverf >= 5 {
        nl = c.n2h;
    }
    if nl == 0 {
        return;
    }
    for _l in 1..=nl {
        cont(c, nin, nout, nscr, a);
        let nc = c.n1h;
        let ni = c.n2h;
        if nc != 0 {
            for _n in 1..=nc {
                if c.iverf >= 5 {
                    cont(c, nin, nout, nscr, a);
                }
                list(c, nin, nout, nscr, a);
            }
        }
        if ni != 0 {
            for _i in 1..=ni {
                list(c, nin, nout, nscr, a);
            }
        }
    }
}

/// File 32: resonance-parameter covariances.
fn file32(c: &mut Ctx, nin: i32, nout: i32, nscr: i32, a: &mut RV) {
    let (mut nb, mut nw) = (0i32, 0i32);
    let nis = c.n1h;
    for _i in 1..=nis {
        cont(c, nin, nout, nscr, a);
        let ner = c.n1h;
        for _j in 1..=ner {
            cont(c, nin, nout, nscr, a);
            let lru = c.l1h;
            let lrf = c.l2h;
            let nro = c.n1h;
            if nro > 0 {
                cont(c, nin, nout, nscr, a);
                let ni = c.n2h;
                for _k in 1..=ni {
                    list(c, nin, nout, nscr, a);
                }
            }
            cont(c, nin, nout, nscr, a);
            let lcomp = c.l2h;
            let nls = c.n1h;
            let isr = c.n2h;
            if lrf != 7 {
                if isr != 0 {
                    if isr == 1 && lrf <= 2 {
                        cont(c, nin, nout, nscr, a);
                    } else if isr == 1 && lrf == 3 {
                        c.listio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
                    } else {
                        c.error("file32", "unknown isr", " ");
                    }
                }
                if lru == 2 {
                    for _l in 1..=nls {
                        list(c, nin, nout, nscr, a);
                    }
                    list(c, nin, nout, nscr, a);
                } else if lcomp == 0 {
                    for _l in 1..=nls {
                        list(c, nin, nout, nscr, a);
                    }
                } else if lcomp == 1 {
                    cont(c, nin, nout, nscr, a);
                    let nsrs = c.n1h;
                    let nlrs = c.n2h;
                    for _k in 1..=nsrs {
                        list(c, nin, nout, nscr, a);
                    }
                    if nlrs > 0 {
                        for _k in 1..=nsrs {
                            list(c, nin, nout, nscr, a);
                        }
                    }
                } else if lcomp == 2 {
                    list(c, nin, nout, nscr, a);
                    cont(c, nin, nout, nscr, a);
                    let nn = c.n1h;
                    let mut ndigit = c.l1h;
                    if !(2..=6).contains(&ndigit) {
                        if ndigit == 0 {
                            ndigit = 2;
                            c.mess("file32", "1illegal value of ndigit", "set default, ndigit=2");
                        } else {
                            c.error("file32", "illegal value of ndigit", " ");
                        }
                    }
                    for _k in 1..=nn {
                        nw = ndigit;
                        c.intgio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
                    }
                } else {
                    c.error("file32", "illegal value of lcomp", " ");
                }
            } else {
                if isr != 0 {
                    c.error("file32", "isr>0 not supported for lrf=7", " ");
                }
                if lcomp == 1 {
                    cont(c, nin, nout, nscr, a);
                    let nsg = c.n1h;
                    for _k in 1..=nsg {
                        cont(c, nin, nout, nscr, a);
                        let njsx = c.l1h;
                        for _l in 1..=njsx {
                            list(c, nin, nout, nscr, a);
                        }
                    }
                    list(c, nin, nout, nscr, a);
                } else if lcomp == 2 {
                    let njs = c.n1h;
                    list(c, nin, nout, nscr, a);
                    for _k in 1..=njs {
                        list(c, nin, nout, nscr, a);
                        list(c, nin, nout, nscr, a);
                    }
                    cont(c, nin, nout, nscr, a);
                    let nn = c.n1h;
                    let mut ndigit = c.l1h;
                    if !(2..=6).contains(&ndigit) {
                        if ndigit == 0 {
                            ndigit = 2;
                            c.mess("file32", "2illegal value of ndigit", "set default, ndigit=2");
                        } else {
                            c.error("file32", "illegal value of ndigit", " ");
                        }
                    }
                    for _k in 1..=nn {
                        nw = ndigit;
                        c.intgio(nin, nout, nscr, &mut a.all(), &mut nb, &mut nw);
                    }
                } else {
                    c.error("file32", "illegal value of lcomp", " ");
                }
            }
        }
    }
}

/// File 34: covariances for angular distributions.
fn file34(c: &mut Ctx, nin: i32, nout: i32, nscr: i32, a: &mut RV) {
    let nmt1 = c.n2h;
    let mth = c.mth;
    for _i in 1..=nmt1 {
        cont(c, nin, nout, nscr, a);
        let mt1 = c.l2h;
        let nl = c.n1h;
        let nl1 = c.n2h;
        for l in 1..=nl {
            for l1 in 1..=nl1 {
                if mt1 != mth || l1 >= l {
                    cont(c, nin, nout, nscr, a);
                    let ni = c.n2h;
                    for _j in 1..=ni {
                        list(c, nin, nout, nscr, a);
                    }
                }
            }
        }
    }
}

/// File 35: covariances for secondary energy distributions.
fn file35(c: &mut Ctx, nin: i32, nout: i32, nscr: i32, a: &mut RV) {
    let nk = c.n1h;
    for _i in 1..=nk {
        list(c, nin, nout, nscr, a);
    }
}

/// File 40.
fn file40(c: &mut Ctx, nin: i32, nout: i32, nscr: i32, a: &mut RV) {
    let ns = c.n1h;
    for _i in 1..=ns {
        cont(c, nin, nout, nscr, a);
        let nl = c.n2h;
        for _j in 1..=nl {
            cont(c, nin, nout, nscr, a);
            let nc = c.n1h;
            let ni = c.n2h;
            for _k in 1..=nc {
                cont(c, nin, nout, nscr, a);
                list(c, nin, nout, nscr, a);
            }
            if ni != 0 {
                for _l in 1..=ni {
                    list(c, nin, nout, nscr, a);
                }
            }
        }
    }
}
