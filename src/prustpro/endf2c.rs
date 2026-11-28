//! Port of ENDF2C (PREPRO 2026-1): convert ENDF data to the standard
//! FORTRAN, C and C++ number format. Reads `ENDFB.IN`, writes `ENDFB.OUT`
//! and the report `ENDF2C.LST`. No input parameters.

use crate::endfio::{ctl, dp, in9, out10, Endf};
use crate::ffmt;
use crate::fio::{self, Rd, L};

struct E2c {
    io: Endf,
    line66: [u8; 66],
    mat: i32,
    mf: i32,
    mt: i32,
    iverse: i32,
}

impl E2c {
    /// `read(10,'(66A1,i4,i2,i3)')`; None at end of file.
    fn read10(&mut self) -> Option<()> {
        match self.io.u.read(10) {
            Rd::Eof => None,
            Rd::Line(l) => {
                let l = fio::pad(l, 75);
                self.line66.copy_from_slice(&l[..66]);
                self.mat = fio::rd_i(&l, 66, 4, 10);
                self.mf = fio::rd_i(&l, 70, 2, 10);
                self.mt = fio::rd_i(&l, 72, 3, 10);
                Some(())
            }
        }
    }

    fn read10_req(&mut self) {
        if self.read10().is_none() {
            fio::runtime_error(10, "ENDFB.IN", "End of file");
        }
    }

    fn line(&self) -> L {
        L::new().b(&self.line66).s(&ctl(self.mat, self.mf, self.mt))
    }

    fn write11(&mut self) {
        let l = self.line();
        self.io.u.write(11, l);
    }

    /// `STANDARD(NFIELD)`.
    fn standard(&mut self, nfield: usize) {
        for i in 0..nfield {
            let f = &mut self.line66[11 * i..11 * i + 11];
            if !f.contains(&b'.') {
                continue;
            }
            let mut save = [0u8; 11];
            save.copy_from_slice(f);
            let zin = in9(f);
            out10(zin, f);
            let zix = in9(f);
            let diff = if zin.abs() <= 0.0 { zin - zix } else { (zin - zix) / zin };
            if diff.abs() > 1.0e-10 {
                let s = L::new()
                    .s(" ")
                    .b(&save)
                    .s(" ")
                    .b(f)
                    .s("  ")
                    .s(&dp(diff, 20, 12))
                    .s(&ffmt::i(self.mat, 5))
                    .s(&ffmt::i(self.mf, 3))
                    .s(&ffmt::i(self.mt, 4))
                    .s(&ffmt::i(self.iverse, 5));
                fio::out(&s);
                self.io.u.write(3, &s);
            }
        }
    }

    /// `FILE1`: MF/MT=1/451 with two added comment lines.
    fn file1(&mut self) {
        self.iverse = 0;
        self.standard(2);
        self.write11();
        let (mut n1a, mut n2a) = (0, 0);
        let mut nwd = 0;
        let mut found = false;
        for iline in 1..=3 {
            self.read10_req();
            self.standard(2);
            let xn1 = in9(&self.line66[44..55]);
            let xn2 = in9(&self.line66[55..66]);
            let n1in = crate::endfio::fint(xn1);
            let n2in = crate::endfio::fint(xn2);
            match iline {
                1 => {
                    n1a = n1in;
                    n2a = n2in;
                    if !(n1a <= 0 || n2a != 0) {
                        self.iverse = 4;
                        nwd = n1a;
                        found = true;
                    }
                }
                2 => {
                    if n2a <= 0 {
                        self.iverse = 5;
                        nwd = n1in;
                        found = true;
                    }
                }
                _ => {
                    self.iverse = 6;
                    nwd = n1in;
                    found = true;
                }
            }
            if found {
                break;
            }
            self.write11();
        }
        let _ = n1a;
        let nwdout = nwd + 2;
        let l = L::new()
            .b(&self.line66[..44])
            .s(&ffmt::i(nwdout, 11))
            .b(&self.line66[55..66])
            .s(&ctl(self.mat, self.mf, self.mt));
        self.io.u.write(11, l);
        for _ in 0..nwd.max(0) {
            self.read10_req();
            self.write11();
        }
        let c = ctl(self.mat, self.mf, self.mt);
        self.io.u.write(
            11,
            format!(" {}**** Program ENDF2C (Version 2026-1) ***{}{}", "*".repeat(13), "*".repeat(12), c),
        );
        self.io.u.write(
            11,
            format!("     Convert ENDF Data to Standard FORTRAN, C and C++ Format      {}", c),
        );
        // copy to SEND: format (22x,44A1,i4,i2,i3)
        loop {
            let l = fio::pad(self.io.u.read_req(10), 75);
            let mat = fio::rd_i(&l, 66, 4, 10);
            let mf = fio::rd_i(&l, 70, 2, 10);
            let mt = fio::rd_i(&l, 72, 3, 10);
            self.mat = mat;
            self.mf = mf;
            self.mt = mt;
            if mt > 0 {
                let o = L::new().s(&" ".repeat(22)).b(&l[22..66]).s(&ctl(mat, mf, mt));
                self.io.u.write(11, o);
                continue;
            }
            break;
        }
        let c = ctl(self.mat, self.mf, self.mt);
        self.io.u.write(11, format!("{}{}", " ".repeat(66), c));
    }
}

pub fn run() -> ! {
    let mut io = Endf::new("ENDF2C");
    io.inp = 2;
    io.outp = 3;
    io.itape = 10;
    io.otape = 11;
    io.timer();
    io.u.open(3, "ENDF2C.LST");
    let filein = "ENDFB.IN    ";
    let fileout = "ENDFB.OUT   ";
    let dash = format!(" {}", "-".repeat(79));
    let hdr = format!(
        "\n\n Convert ENDF Data to Standard FORTRAN, C and C++ Format\n             Program ENDF2C (Version 2026-1)\n{dash}\n ENDF Input and OutOut Filenames\n {filein}\n {fileout}\n{dash}\n ENDF Tape Label (TPID)\n{dash}"
    );
    fio::out(&hdr);
    io.u.write(3, &hdr);
    io.u.open(10, filein);
    io.u.open(11, fileout);
    let mut e = E2c { io, line66: [b' '; 66], mat: 0, mf: 0, mt: 0, iverse: 0 };
    let mut matlast = -99999;
    e.iverse = 0;
    e.read10_req();
    let l = e.line();
    fio::out(&l);
    e.io.u.write(3, &l);
    let colhdr = format!(
        " Input       Output         Fractional Diff.    MAT MF  MT  Seq# Format\n{dash}"
    );
    let mut have = false;
    if e.mf == 0 && e.mt == 0 {
        e.write11();
        let s = format!("{dash}\n{colhdr}");
        fio::out(&s);
        e.io.u.write(3, &s);
    } else {
        let s = format!(
            "{dash}\n WARNING - This DOES NOT look like a TPID. Expect MF = MT = 0 on TPID.\n           Will Assume file has NO TPID and TRY to Translate this line.\n{dash}\n{colhdr}"
        );
        fio::out(&s);
        e.io.u.write(3, &s);
        have = true;
    }
    loop {
        if !have && e.read10().is_none() {
            let s = format!("{dash}\n ERROR - ENDF file does not end with TEND Line.\n{dash}");
            fio::out(&s);
            e.io.u.write(3, &s);
            e.io.enderror();
        }
        have = false;
        if e.mat != matlast {
            matlast = e.mat;
            e.iverse = 0;
        }
        if e.mt <= 0 {
            let c = ctl(e.mat, e.mf, e.mt);
            e.io.u.write(11, format!("{}{}", " ".repeat(66), c));
            if e.mat < 0 {
                e.io.endit();
            }
            continue;
        }
        if e.mf == 1 && e.mt == 451 {
            e.file1();
            continue;
        }
        e.standard(6);
        e.write11();
    }
}
