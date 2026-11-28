//! Port of DICTIN (PREPRO 2026-1): make a new section index (MF=1/MT=451
//! dictionary) for each evaluation. Input `DICTIN.INP`, listing
//! `DICTIN.LST`.

use crate::endfio::Endf;
use crate::ffmt::i;
use crate::fio::{self, Rd, L};
use crate::linear::blank;

const MAXTAB: usize = 120000;

fn ctl(mat: i32, mf: i32, mt: i32) -> String {
    format!("{}{}{}", i(mat, 4), i(mf, 2), i(mt, 3))
}

struct D {
    io: Endf,
}

impl D {
    fn rd(&mut self) -> Vec<u8> {
        self.io.u.read_req(10)
    }
    /// (66X,I4,I2,I3)
    fn rd390(&mut self) -> (i32, i32, i32) {
        let l = self.rd();
        (fio::rd_i(&l, 66, 4, 10), fio::rd_i(&l, 70, 2, 10), fio::rd_i(&l, 72, 3, 10))
    }
    /// (72X,I3)
    fn rd400(&mut self) -> i32 {
        let l = self.rd();
        fio::rd_i(&l, 72, 3, 10)
    }
    fn w(&mut self, l: L) {
        self.io.u.write(11, l);
    }
}

pub fn run() -> ! {
    let mut io = Endf::new("DICTIN");
    io.timer();
    io.inp = 2;
    io.outp = 3;
    io.itape = 10;
    io.otape = 11;
    io.u.open(3, "DICTIN.LST");
    let istat1 = if io.u.open_old(2, "DICTIN.INP") { 0 } else { 1 };
    let d60 = format!(" {}", "-".repeat(60));
    io.both(format!(" Create New Section Index (DICTIN 2026-1)\n{d60}"));
    // READIN
    let mut namein = None;
    let mut nameout = None;
    if istat1 != 1 {
        if let Rd::Line(l) = io.u.read(2) {
            let mut n = fio::rd_a(&l, 0, 72);
            if blank(&n) {
                n = fio::rd_a(b"ENDFIN.IN", 0, 72);
            }
            namein = Some(n);
            if let Rd::Line(l) = io.u.read(2) {
                let mut n = fio::rd_a(&l, 0, 72);
                if blank(&n) {
                    n = fio::rd_a(b"ENDFIN.OUT", 0, 72);
                }
                nameout = Some(n);
            }
        }
    }
    let namein = namein.unwrap_or_else(|| fio::rd_a(b"ENDFB.IN", 0, 72));
    let nameout = nameout.unwrap_or_else(|| fio::rd_a(b"ENDFB.OUT", 0, 72));
    let s = L::new()
        .s(" ENDF/B Input and Output Data Filenames\n")
        .s(&d60)
        .s("\n ")
        .b(&namein)
        .s("\n ")
        .b(&nameout);
    io.both(&s);
    let no = String::from_utf8_lossy(&nameout).into_owned();
    io.u.open(11, &no);
    let ni = String::from_utf8_lossy(&namein).into_owned();
    if !io.u.open_old(10, &ni) {
        let s = L::new().s("\n\n ERROR - Opening Old file\n ").b(&namein).s("\n\n");
        io.both(&s);
        io.enderror();
    }
    let mut d = D { io };
    let mut lasmat = -99999;
    d.io.both(format!("{d60}\n First Pass. Reading ENDF/B Data\n{d60}\n  MAT   Messages\n{d60}"));
    let mut index = 0usize;
    let mut mattab = vec![0i32; MAXTAB + 2];
    let mut mftab = vec![0i32; MAXTAB + 2];
    let mut mttab = vec![0i32; MAXTAB + 2];
    let mut ncard = vec![0i32; MAXTAB + 2];
    d.rd();
    loop {
        // label 10
        let (mat, mf, mt) = d.rd390();
        if mt <= 0 {
            if mat < 0 {
                break;
            }
            continue;
        }
        if mat != lasmat {
            if mf == 1 && mt == 451 {
                if mat > lasmat {
                    d.io.both(i(mat, 5));
                } else {
                    d.io.both(format!("{}   WARNING Material is not in MAT Order", i(mat, 5)));
                }
                lasmat = mat;
            } else {
                if mat > lasmat {
                    d.io.both(format!("{}   No Section MF=1, MT=451...Material Copied", i(mat, 5)));
                } else {
                    d.io.both(format!(
                        "{}   WARNING Material is not in MAT Order\n        No Section MF=1, MT=451...Material Copied",
                        i(mat, 5)
                    ));
                }
                lasmat = mat;
                loop {
                    let (m2, _, _) = d.rd390();
                    if m2 <= 0 {
                        break;
                    }
                }
                continue;
            }
        }
        // label 90
        index += 1;
        if index > MAXTAB {
            d.io.both(format!("\n\n\nOver {} Sections---Execution Terminated", i(MAXTAB as i32, 5)));
            d.io.endit();
        }
        mattab[index] = mat;
        mftab[index] = mf;
        mttab[index] = mt;
        let mut kcard = 0;
        loop {
            let mt2 = d.rd400();
            kcard += 1;
            if mt2 <= 0 {
                break;
            }
        }
        ncard[index] = kcard;
    }
    if index == 0 {
        d.io.both(format!("{d60}\n No Section MF=1, MT=451 on Tape...Execution Terminated"));
        d.io.enderror();
    }
    // second pass
    let mut ndex = 1usize;
    d.io.u.rewind(10);
    d.io.both(format!("{d60}\n Second Pass. Creating New Index\n{d60}\n ENDF/B  MAT MF  MT  Lines  Mod\n Format"));
    let l = d.rd();
    let card1 = fio::rd_a(&l, 0, 80);
    let mut icard = 1;
    d.w(L::new().b(&card1));
    let mut ocard = 1;
    // SAVEd locals
    let mut iverse;
    let (mut mfold, mut mtold, mut modold) = (0, 0, 0);
    let fmthol = ["IV", " V", "VI"];
    loop {
        // label 160
        let l = d.rd();
        let card2 = fio::rd_a(&l, 0, 55);
        let n2 = fio::rd_i(&l, 55, 11, 10);
        let mut mat = fio::rd_i(&l, 66, 4, 10);
        let mut mf = fio::rd_i(&l, 70, 2, 10);
        let mut mt = fio::rd_i(&l, 72, 3, 10);
        icard += 1;
        if mt <= 0 {
            d.w(L::new().s(&" ".repeat(66)).s(&ctl(mat, mf, mt)));
            ocard += 1;
            if mat < 0 {
                break;
            }
            continue;
        }
        if !(mf == 1 && mt == 451) {
            // label 340
            d.w(L::new().b(&card2).s(&i(n2, 11)).s(&ctl(mat, mf, mt)));
            ocard += 1;
            loop {
                let l = d.rd();
                let card5 = fio::rd_a(&l, 0, 66);
                mat = fio::rd_i(&l, 66, 4, 10);
                mf = fio::rd_i(&l, 70, 2, 10);
                mt = fio::rd_i(&l, 72, 3, 10);
                icard += 1;
                if mt <= 0 {
                    break;
                }
                d.w(L::new().b(&card5).s(&ctl(mat, mf, mt)));
                ocard += 1;
            }
            d.w(L::new().s(&" ".repeat(66)).s(&ctl(mat, mf, mt)));
            ocard += 1;
            continue;
        }
        let mut nxcold = n2;
        let mut itop = ndex;
        while itop <= index {
            if mat != mattab[itop] {
                break;
            }
            itop += 1;
        }
        let nxc = itop as i32 - ndex as i32;
        let itop = itop - 1;
        let rd430 = |d: &mut D| {
            let l = d.rd();
            (
                fio::rd_a(&l, 0, 44),
                fio::rd_i(&l, 44, 11, 10),
                fio::rd_i(&l, 55, 11, 10),
                fio::rd_i(&l, 66, 4, 10),
                fio::rd_i(&l, 70, 2, 10),
                fio::rd_i(&l, 72, 3, 10),
            )
        };
        let w430 = |c: &[u8], a: i32, b: i32, mat: i32, mf: i32, mt: i32| L::new().b(c).s(&i(a, 11)).s(&i(b, 11)).s(&ctl(mat, mf, mt));
        let (card3, nwd0, n3, m1, f1, t1) = rd430(&mut d);
        let mut nwd = nwd0;
        mat = m1;
        mf = f1;
        mt = t1;
        icard += 1;
        if nwd > 0 && n3 == 0 {
            // label 210
            d.w(L::new().b(&card2).s(&i(nxc, 11)).s(&ctl(mat, mf, mt)));
            ocard += 1;
            d.w(w430(&card3, nwd + 1, n3, mat, mf, mt));
            ocard += 1;
            iverse = 1;
        } else if n3 > 0 {
            let (card4, a, b, m, f, t) = rd430(&mut d);
            nwd = a;
            nxcold = b;
            mat = m;
            mf = f;
            mt = t;
            icard += 1;
            d.w(L::new().b(&card2).s(&i(n2, 11)).s(&ctl(mat, mf, mt)));
            ocard += 1;
            d.w(w430(&card3, 0, n3, mat, mf, mt));
            ocard += 1;
            d.w(w430(&card4, nwd, nxcold, mat, mf, mt));
            ocard += 1;
            let (card4, a, b, m, f, t) = rd430(&mut d);
            nwd = a;
            nxcold = b;
            mat = m;
            mf = f;
            mt = t;
            icard += 1;
            d.w(w430(&card4, nwd + 1, nxc, mat, mf, mt));
            ocard += 1;
            iverse = 3;
        } else {
            // label 200
            let (card4, a, b, m, f, t) = rd430(&mut d);
            nwd = a;
            nxcold = b;
            mat = m;
            mf = f;
            mt = t;
            icard += 1;
            d.w(L::new().b(&card2).s(&i(n2, 11)).s(&ctl(mat, mf, mt)));
            ocard += 1;
            d.w(w430(&card3, 0, n3, mat, mf, mt));
            ocard += 1;
            d.w(w430(&card4, nwd + 1, nxc, mat, mf, mt));
            ocard += 1;
            iverse = 2;
        }
        // label 220
        d.io.both(&d60);
        for _ in 0..nwd.max(0) {
            let l = d.rd();
            let card5 = fio::rd_a(&l, 0, 66);
            icard += 1;
            d.w(L::new().b(&card5).s(&ctl(mat, mf, mt)));
            ocard += 1;
        }
        d.w(L::new().s(" **************** Program DICTIN (VERSION 2026-1) ****************").s(&ctl(mat, mf, mt)));
        ocard += 1;
        let mut ii = 1;
        let rd460 = |d: &mut D| {
            let l = d.rd();
            (fio::rd_i(&l, 22, 11, 10), fio::rd_i(&l, 33, 11, 10), fio::rd_i(&l, 55, 11, 10))
        };
        if nxcold > 0 {
            (mfold, mtold, modold) = rd460(&mut d);
            icard += 1;
        }
        for j in ndex..=itop {
            if mftab[j] == 1 && mttab[j] == 451 {
                ncard[j] = ncard[j] + nxc - nxcold;
            }
            let mut md = 1;
            if ii <= nxcold {
                loop {
                    // label 260
                    if mfold > mftab[j] {
                        break;
                    }
                    if mfold == mftab[j] {
                        if mtold == mttab[j] {
                            md = modold;
                            break;
                        }
                        if mtold > mttab[j] {
                            break;
                        }
                    }
                    // label 270
                    ii += 1;
                    if ii > nxcold {
                        break;
                    }
                    (mfold, mtold, modold) = rd460(&mut d);
                    icard += 1;
                }
            }
            // label 290
            if mftab[j] == 1 && mttab[j] == 451 {
                ncard[j] += 1;
                let s = format!(
                    "  {}   {}{}{}{}{}",
                    fmthol[iverse - 1],
                    i(mattab[j], 5),
                    i(mftab[j], 3),
                    i(mttab[j], 4),
                    i(ncard[j], 7),
                    i(md, 5)
                );
                d.io.both(&s);
            } else {
                let s = format!("       {}{}{}{}{}", i(mattab[j], 5), i(mftab[j], 3), i(mttab[j], 4), i(ncard[j], 7), i(md, 5));
                d.io.both(&s);
            }
            d.w(L::new()
                .s(&" ".repeat(22))
                .s(&i(mftab[j], 11))
                .s(&i(mttab[j], 11))
                .s(&i(ncard[j], 11))
                .s(&i(md, 11))
                .s(&ctl(mat, mf, mt)));
            ocard += 1;
        }
        ndex = itop + 1;
        d.w(L::new().s(&" ".repeat(66)).s(&ctl(mat, mf, 0)));
        ocard += 1;
        loop {
            let mt2 = d.rd400();
            icard += 1;
            if mt2 <= 0 {
                break;
            }
        }
    }
    d.io.both(format!("{d60}\n Original Line Count{}\n Final Line Count   {}\n{d60}", i(icard, 7), i(ocard, 7)));
    d.io.endit();
}
