//! Port of LINEAR (PREPRO 2026-1): linearize ENDF/B cross sections
//! (MF=1 MT=452/455/456, MF=3, 9, 10, 23, 27) to within an allowable
//! uncertainty. Input `LINEAR.INP`, listing `LINEAR.LST`.
//!
//! Arrays are 1-based as in the Fortran (element 0 is not used).

use crate::endfio::{f9s, fint, in9, incore10, out9, terpit, zaholm, Endf};
use crate::ffmt;
use crate::fio::{self, Rd, Val, L};

const MAXPAGE: usize = 120000;
const ETHERMAL: f64 = 2.53e-02;

/// One line of `LINEAR.INP` (None at end of file).
fn rd_line(io: &mut Endf) -> Option<Vec<u8>> {
    let u = io.inp;
    match io.u.read(u) {
        Rd::Line(l) => Some(l),
        Rd::Eof => None,
    }
}

/// Iw input for a READ with ERR=: None on a bad field.
pub fn rdi(l: &[u8], s: usize, w: usize) -> Option<i32> {
    let f = fio::rd_a(l, s, w);
    ffmt::read_i(&f).ok()
}

pub fn blank(b: &[u8]) -> bool {
    b.iter().all(|&c| c == b' ')
}

struct Linear {
    io: Endf,
    // /COUNTS/
    matnow: i32,
    mfnow: i32,
    matin: i32,
    matout: i32,
    n2toti: i32,
    n2toto: i32,
    // /PAGER/
    npage: usize,
    npage2: usize,
    np1p1: usize,
    np2m1: usize,
    estart: f64,
    // /SLIM/
    istart: i32,
    nothin: i32,
    ithin1: usize,
    ithin2: usize,
    ithin3: usize,
    mtend: i32,
    // /OLDTAB/
    ixylow: usize,
    ixyhi: usize,
    left: usize,
    // /NEWTAB/
    n2core: usize,
    n2scr: usize,
    n2tot: usize,
    n2p1: usize,
    n2p2: usize,
    // /OKERR3/
    ener3: [f64; 22],
    er3: [f64; 22],
    kerr3: i32,
    maxer3: usize,
    iza: i32,
    fmthol: [u8; 4],
    liso: i32,
    nocoment: i32,
    temp3: f64,
    iverse: i32,
    iscr: i32,
    monitr: i32,
    field: [[u8; 11]; 6],
    // /MATZA/
    modget: i32,
    nmatza: usize,
    matmin: [i32; 102],
    mfmin: [i32; 102],
    mtmin: [i32; 102],
    matmax: [i32; 102],
    mfmax: [i32; 102],
    mtmax: [i32; 102],
    istat1: i32,
    istat2: i32,
    namein: Vec<u8>,
    nameout: Vec<u8>,
    nbt: [i32; 101],
    int: [i32; 101],
    xin: Vec<f64>,
    yin: Vec<f64>,
    xout: Vec<f64>,
    yout: Vec<f64>,
    // SAVEd locals
    progdoc: [[u8; 66]; 8],
    nubar_ylin: f64,
    filex_yaprox: f64,
    filex_xsave: f64,
    filex_ysave: f64,
    th_imused: i32,
    th_dsign: f64,
    th_slpmax: f64,
    th_slpmin: f64,
    miner3: usize,
    nubar_hold: (f64, f64, i32, i32, i32, i32),
    lambda: Vec<[u8; 66]>,
}

impl Linear {
    fn new() -> Linear {
        let mut pd = [[b' '; 66]; 8];
        let t: [&str; 8] = [
            " ***************** Program LINEAR (VERSION 2026-1) ***********",
            " For All Data Greater than12345678901 barns in Absolute Value ",
            " Data Linearized to Within an Accuracy of12345678901 per-cent ",
            " Data Linearized Using Energy Dependent Uncertainty           ",
            "      Energy    Accuracy                                      ",
            "        (eV)  (per-cent)                                      ",
            " ----------- -----------                                      ",
            " 12345678901 12345678901                                      ",
        ];
        for (k, s) in t.iter().enumerate() {
            pd[k][..s.len()].copy_from_slice(s.as_bytes());
        }
        Linear {
            io: Endf::new("LINEAR"),
            matnow: 0,
            mfnow: 0,
            matin: 0,
            matout: 0,
            n2toti: 0,
            n2toto: 0,
            npage: 0,
            npage2: 0,
            np1p1: 0,
            np2m1: 0,
            estart: 0.0,
            istart: 0,
            nothin: 0,
            ithin1: 0,
            ithin2: 0,
            ithin3: 0,
            mtend: 0,
            ixylow: 0,
            ixyhi: 0,
            left: 0,
            n2core: 0,
            n2scr: 0,
            n2tot: 0,
            n2p1: 0,
            n2p2: 0,
            ener3: [0.0; 22],
            er3: [0.0; 22],
            kerr3: 0,
            maxer3: 0,
            iza: 0,
            fmthol: [0u8; 4],
            liso: 0,
            nocoment: 0,
            temp3: 0.0,
            iverse: 0,
            iscr: 12,
            monitr: 0,
            field: [[0u8; 11]; 6],
            modget: 0,
            nmatza: 0,
            matmin: [0; 102],
            mfmin: [0; 102],
            mtmin: [0; 102],
            matmax: [0; 102],
            mfmax: [0; 102],
            mtmax: [0; 102],
            istat1: 0,
            istat2: 0,
            namein: vec![b' '; 72],
            nameout: vec![b' '; 72],
            nbt: [0; 101],
            int: [0; 101],
            xin: vec![0.0; MAXPAGE + 2],
            yin: vec![0.0; MAXPAGE + 2],
            xout: vec![0.0; 2 * MAXPAGE + 4],
            yout: vec![0.0; 2 * MAXPAGE + 4],
            progdoc: pd,
            nubar_ylin: 0.0,
            filex_yaprox: 0.0,
            filex_xsave: 0.0,
            filex_ysave: 0.0,
            th_imused: 0,
            th_dsign: 0.0,
            th_slpmax: 0.0,
            th_slpmin: 0.0,
            miner3: 2,
            nubar_hold: (0.0, 0.0, 0, 0, 0, 0),
            lambda: vec![[b' '; 66]; 51],
        }
    }

    fn both<S: AsRef<[u8]>>(&mut self, s: S) {
        self.io.both(s);
    }

    fn dash() -> String {
        format!("  {}", "-".repeat(78))
    }

    fn zab(&self) -> [u8; 12] {
        zaholm(self.iza, self.liso)
    }

    /// The listing line start `2X,12A1,I5,I4,I4,3X,A2,3X`.
    fn id_line(&self) -> L {
        let h = &self.io.h;
        L::new()
            .s("  ")
            .b(&self.zab())
            .s(&ffmt::i(h.mat, 5))
            .s(&ffmt::i(h.mf, 4))
            .s(&ffmt::i(h.mt, 4))
            .s("   ")
            .b(&self.fmthol[..2])
            .s("   ")
    }

    // ------------------------------------------------------------------

    fn fileio(&mut self) {
        self.io.inp = 2;
        self.io.outp = 3;
        self.io.itape = 10;
        self.io.otape = 11;
        self.iscr = 12;
        self.io.u.open(3, "LINEAR.LST");
        self.io.u.open_scratch(12);
        self.istat1 = if self.io.u.open_old(2, "LINEAR.INP") { 0 } else { 1 };
    }

    fn filio2(&mut self) {
        let no = fio::L(self.nameout.clone());
        let no = String::from_utf8_lossy(&no.0).into_owned();
        self.io.u.open(11, &no);
        let ni = String::from_utf8_lossy(&self.namein).into_owned();
        self.istat2 = if self.io.u.open_old(10, &ni) { 0 } else { 1 };
    }

    fn readin(&mut self) {
        let errmin = 1.0e-03;
        let xcinp = 1.0e-30;
        let mess1 = [" MAT", "  ZA"];
        let mess2 = [" Off", "  On"];
        let mut ok = false;
        if self.istat1 != 1 {
            if let Some(l) = rd_line(&mut self.io) {
                let a = rdi(&l, 0, 11);
                let b = rdi(&l, 11, 11);
                let f = fio::rd_a(&l, 22, 11);
                let c = rdi(&l, 33, 11);
                let d = rdi(&l, 44, 11);
                // ERR= leaves the items read before the bad one
                if let Some(a) = a {
                    self.modget = a;
                    if let Some(b) = b {
                        self.monitr = b;
                        self.field[0].copy_from_slice(&f);
                        if let Some(c) = c {
                            self.nothin = c;
                            if let Some(d) = d {
                                self.nocoment = d;
                                ok = true;
                            }
                        }
                    }
                }
            }
        }
        if ok {
            self.io.err.xclow = in9(&self.field[0]);
            self.io.err.xclow = xcinp;
        } else {
            self.istat1 = 1;
            self.modget = 0;
            self.monitr = 0;
            self.io.err.xclow = xcinp;
            self.nothin = 1;
            self.nocoment = 0;
        }
        self.modget = if self.modget <= 0 { 0 } else { 1 };
        self.monitr = if self.monitr <= 0 { 0 } else { 1 };
        self.nothin = if self.nothin <= 0 { 0 } else { 1 };
        let s = format!(
            "  Retrieval Criteria-----------       {}\n  Monitor Mode-----------------       {}",
            mess1[self.modget as usize], mess2[self.monitr as usize]
        );
        self.both(&s);
        out9(self.io.err.xclow, &mut self.field[0]);
        let s = L::new()
            .s("  Minimum Cross Section--------")
            .b(&self.field[0])
            .s(" (No longer an Input Option)");
        self.both(&s);
        if self.nothin == 0 {
            self.both("  Keep Evaluated Data Points---         No");
        } else {
            self.both("  Keep Evaluated Data Points---        Yes");
        }
        if self.nocoment != 0 {
            self.both("  Add Comments to MF/MT=1/451--         No");
        }
        // filenames
        let mut got_in = false;
        let mut got_out = false;
        if self.istat1 != 1 {
            if let Some(l) = rd_line(&mut self.io) {
                self.namein = fio::rd_a(&l, 0, 72);
                if blank(&self.namein) {
                    self.namein = fio::rd_a(b"ENDFB.IN", 0, 72);
                }
                got_in = true;
                if let Some(l) = rd_line(&mut self.io) {
                    self.nameout = fio::rd_a(&l, 0, 72);
                    if blank(&self.nameout) {
                        self.nameout = fio::rd_a(b"ENDFB.OUT", 0, 72);
                    }
                    got_out = true;
                }
            }
        }
        if !got_in {
            self.namein = fio::rd_a(b"ENDFB.IN", 0, 72);
        }
        if !got_out {
            self.nameout = fio::rd_a(b"ENDFB.OUT", 0, 72);
            self.istat1 = 1;
        }
        let s = L::new()
            .s(&Self::dash())
            .s("\n  ENDF/B Input and Output Data Filenames\n  ")
            .b(&self.namein)
            .s("\n  ")
            .b(&self.nameout);
        self.both(&s);
        self.filio2();
        if self.istat2 == 1 {
            let s = L::new().s("\n\n ERROR - Opening ENDF/B data file\n ").b(&self.namein).s("\n\n");
            self.both(&s);
            self.io.enderror();
        }
        let d = Self::dash();
        if self.modget == 0 {
            self.both(format!("{d}\n  Requested Ranges\n{d}\n       Minimum       Maximum   \n     MAT MF  MT    MAT MF  MT\n{d}"));
        } else {
            self.both(format!("{d}\n  Requested Ranges\n{d}\n       Minimum       Maximum   \n      ZA MF  MT     ZA MF  MT\n{d}"));
        }
        let rng = |l: &[u8]| -> Option<[i32; 6]> {
            Some([
                rdi(l, 0, 6)?,
                rdi(l, 6, 2)?,
                rdi(l, 8, 3)?,
                rdi(l, 11, 6)?,
                rdi(l, 17, 2)?,
                rdi(l, 19, 3)?,
            ])
        };
        let mut have = false;
        if self.istat1 != 1 {
            if let Some(l) = rd_line(&mut self.io) {
                // ERR= after a partial read keeps the values read; the
                // values are then reset at label 120 anyway.
                if let Some(v) = rng(&l) {
                    self.matmin[1] = v[0];
                    self.mfmin[1] = v[1];
                    self.mtmin[1] = v[2];
                    self.matmax[1] = v[3];
                    self.mfmax[1] = v[4];
                    self.mtmax[1] = v[5];
                    have = true;
                }
            }
        }
        if !have {
            self.istat1 = 1;
            self.matmin[1] = 0;
            self.mfmin[1] = 0;
            self.mtmin[1] = 0;
            self.matmax[1] = 0;
            self.mfmax[1] = 0;
            self.mtmax[1] = 0;
        }
        if self.matmin[1] <= 0 {
            self.matmin[1] = 0;
        }
        if self.mfmin[1] <= 0 {
            self.mfmin[1] = 0;
        }
        if self.mtmin[1] <= 0 {
            self.mtmin[1] = 0;
        }
        let rline = |a: &[i32; 102], b: &[i32; 102], c: &[i32; 102], d: &[i32; 102], e: &[i32; 102], f: &[i32; 102], k: usize| {
            format!(
                "  {}{}{}{}{}{}",
                ffmt::i(a[k], 6),
                ffmt::i(b[k], 3),
                ffmt::i(c[k], 4),
                ffmt::i(d[k], 7),
                ffmt::i(e[k], 3),
                ffmt::i(f[k], 4)
            )
        };
        let mut goto170 = false;
        if !(self.matmin[1] > 0 || self.mfmin[1] > 0 || self.mtmin[1] > 0 || self.matmax[1] > 0 || self.mfmax[1] > 0 || self.mtmax[1] > 0) {
            self.matmax[1] = 9999;
            self.mfmax[1] = 99;
            self.mtmax[1] = 999;
            self.modget = 0;
            let s = format!(
                "{} (Default Option)",
                rline(&self.matmin, &self.mfmin, &self.mtmin, &self.matmax, &self.mfmax, &self.mtmax, 1)
            );
            self.both(&s);
            self.nmatza = 2;
            goto170 = true;
        } else {
            self.fix_range(1);
            let s = rline(&self.matmin, &self.mfmin, &self.mtmin, &self.matmax, &self.mfmax, &self.mtmax, 1);
            self.both(&s);
            let mut n = 2usize;
            while n <= 101 {
                self.nmatza = n;
                if self.istat1 == 1 {
                    goto170 = true;
                    break;
                }
                let v = rd_line(&mut self.io).and_then(|l| rng(&l));
                let Some(v) = v else {
                    self.istat1 = 1;
                    goto170 = true;
                    break;
                };
                self.matmin[n] = v[0];
                self.mfmin[n] = v[1];
                self.mtmin[n] = v[2];
                self.matmax[n] = v[3];
                self.mfmax[n] = v[4];
                self.mtmax[n] = v[5];
                if v.iter().all(|&x| x <= 0) {
                    goto170 = true;
                    break;
                }
                self.fix_range(n);
                let s = rline(&self.matmin, &self.mfmin, &self.mtmin, &self.matmax, &self.mfmax, &self.mtmax, n);
                self.both(&s);
                n += 1;
            }
            if !goto170 {
                self.both("  Over 100 Ranges----Execution Terminated");
                self.io.enderror();
            }
        }
        let _ = goto170;
        self.nmatza -= 1;
        // allowable uncertainty
        #[derive(PartialEq)]
        enum St {
            L190,
            L200,
        }
        let mut st;
        if self.istat1 == 1 {
            st = St::L190;
        } else {
            match rd_line(&mut self.io) {
                Some(l) => {
                    let f = fio::rd_a(&l, 0, 22);
                    self.field[0].copy_from_slice(&f[..11]);
                    self.field[1].copy_from_slice(&f[11..22]);
                    self.ener3[1] = in9(&self.field[0]);
                    self.er3[1] = in9(&self.field[1]);
                    st = St::L200;
                }
                None => st = St::L190,
            }
        }
        loop {
            if st == St::L190 {
                self.istat1 = 1;
                self.ener3[1] = 0.0;
                self.er3[1] = 0.0;
            }
            // label 200
            out9(self.ener3[1], &mut self.field[0]);
            let hdr = format!("{d}\n  Allowable Uncertainty\n{d}\n      Energy Uncertainty   per-cent\n{d}\n");
            if self.er3[1] > 0.0 {
                let percnt = 100.0 * self.er3[1];
                out9(self.er3[1], &mut self.field[1]);
                let s = L::new().s(&hdr).s(" ").b(&self.field[0]).s(" ").b(&self.field[1]).s(&ffmt::f(percnt, 11, 4));
                self.both(&s);
            } else {
                self.er3[1] = errmin;
                let percnt = 100.0 * errmin;
                out9(self.er3[1], &mut self.field[1]);
                let s = L::new()
                    .s(&hdr)
                    .s(" ")
                    .b(&self.field[0])
                    .s(" ")
                    .b(&self.field[1])
                    .s(&ffmt::f(percnt, 11, 4))
                    .s(" (Default Option)");
                self.both(&s);
                if !(self.ener3[1] > 0.0) {
                    self.maxer3 = 2;
                    break;
                }
            }
            // label 220
            let mut m = 2usize;
            let mut restart = false;
            let mut done = false;
            while m <= 21 {
                self.maxer3 = m;
                if self.istat1 == 1 {
                    done = true;
                    break;
                }
                let Some(l) = rd_line(&mut self.io) else {
                    restart = true;
                    break;
                };
                let f = fio::rd_a(&l, 0, 22);
                self.field[0].copy_from_slice(&f[..11]);
                self.field[1].copy_from_slice(&f[11..22]);
                self.ener3[m] = in9(&self.field[0]);
                self.er3[m] = in9(&self.field[1]);
                if self.ener3[m] <= 0.0 && self.er3[m] <= 0.0 {
                    done = true;
                    break;
                }
                out9(self.ener3[m], &mut self.field[0]);
                if self.er3[m] > 0.0 {
                    let percnt = 100.0 * self.er3[m];
                    out9(self.er3[m], &mut self.field[1]);
                    let s = L::new().s(" ").b(&self.field[0]).s(" ").b(&self.field[1]).s(&ffmt::f(percnt, 11, 4));
                    self.both(&s);
                } else {
                    self.er3[m] = errmin;
                    let percnt = 100.0 * errmin;
                    out9(self.er3[m], &mut self.field[1]);
                    let s = L::new()
                        .s(" ")
                        .b(&self.field[0])
                        .s(" ")
                        .b(&self.field[1])
                        .s(&ffmt::f(percnt, 11, 4))
                        .s(" (Default Option)");
                    self.both(&s);
                }
                if self.ener3[m] < self.ener3[m - 1] {
                    self.both("  Energies MUST be in Ascending Order----Execution Terminated");
                    self.io.enderror();
                }
                m += 1;
            }
            if restart {
                st = St::L190;
                continue;
            }
            if !done {
                self.both("  Over 20 Ranges----Execution Terminated");
                self.io.enderror();
            }
            break;
        }
        // label 260
        self.maxer3 -= 1;
        self.kerr3 = 0;
        if self.maxer3 > 1 {
            self.kerr3 = 1;
        }
        self.io.err.errxc3 = self.er3[self.maxer3];
        self.io.err.seterror(1);
    }

    fn fix_range(&mut self, n: usize) {
        if self.matmin[n] <= 0 {
            self.matmin[n] = 0;
        }
        if self.mfmin[n] <= 0 {
            self.mfmin[n] = 0;
        }
        if self.mtmin[n] <= 0 {
            self.mtmin[n] = 0;
        }
        if self.matmax[n] < self.matmin[n] {
            self.matmax[n] = self.matmin[n];
        }
        if self.matmax[n] <= 0 && self.modget == 0 {
            self.matmax[n] = 9999;
        }
        if self.matmax[n] <= 0 && self.modget != 0 {
            self.matmax[n] = 999999;
        }
        if self.mfmax[n] <= 0 {
            self.mfmax[n] = 99;
        }
        if self.mtmax[n] <= 0 {
            self.mtmax[n] = 999;
        }
    }

    fn mat_totals(&mut self) {
        let d = Self::dash();
        let s = format!("{d}\n{} MAT Totals{}{}\n{d}", " ".repeat(48), ffmt::i(self.matin, 9), ffmt::i(self.matout, 9));
        self.both(&s);
        self.io.maxie3(1);
        self.io.timemat();
        self.n2toti += self.matin;
        self.n2toto += self.matout;
    }

    fn nextmt(&mut self) {
        loop {
            self.io.conti();
            if self.io.h.mt <= 0 {
                if self.io.h.mat < 0 {
                    // label 90
                    if self.matnow > 0 {
                        let m = self.matnow;
                        self.io.outf(m);
                        self.io.outm();
                        if !(self.matin <= 0 && self.matout <= 0) {
                            self.mat_totals();
                        }
                    }
                    self.io.h.mat = -1;
                    self.io.h.mf = 0;
                    self.io.h.mt = 0;
                    self.matin = 0;
                    self.matout = 0;
                    return;
                }
                continue;
            }
            self.iza = fint(self.io.h.c1);
            let mut low = if self.modget != 0 { 1 } else { 0 };
            let mut iway = 1;
            let (math, mfh, mth) = (self.io.h.mat, self.io.h.mf, self.io.h.mt);
            let mut found = false;
            for k in 1..=self.nmatza {
                if self.modget == 0 {
                    if math < self.matmin[k] {
                        low = 1;
                        continue;
                    }
                    if !(math == self.matmin[k] || math <= self.matmax[k]) {
                        continue;
                    }
                    low = 1;
                } else {
                    if self.iza < self.matmin[k] {
                        continue;
                    }
                    if !(self.iza == self.matmin[k] || self.iza <= self.matmax[k]) {
                        continue;
                    }
                }
                // label 50
                if iway < 2 {
                    iway = 2;
                }
                if mfh < self.mfmin[k] {
                    continue;
                }
                if mfh != self.mfmin[k] && mfh > self.mfmax[k] {
                    continue;
                }
                iway = 3;
                if mth < self.mtmin[k] {
                    continue;
                }
                if mth == self.mtmin[k] || mth <= self.mtmax[k] {
                    found = true;
                    break;
                }
            }
            if found {
                // label 110
                if self.matnow > 0 {
                    if self.matnow == math && self.mfnow == mfh {
                        self.matnow = math;
                        self.mfnow = mfh;
                        return;
                    }
                    let m = self.matnow;
                    self.io.outf(m);
                    if self.matnow == math {
                        self.matnow = math;
                        self.mfnow = mfh;
                        return;
                    }
                    self.io.outm();
                    if !(self.matin <= 0 && self.matout <= 0) {
                        self.mat_totals();
                    }
                }
                self.matin = 0;
                self.matout = 0;
                self.fmthol = *b"6   ";
                self.iverse = 6;
                self.temp3 = 0.0;
                self.io.maxie0();
                self.matnow = math;
                self.mfnow = mfh;
                return;
            }
            if low <= 0 {
                // label 90
                if self.matnow > 0 {
                    let m = self.matnow;
                    self.io.outf(m);
                    self.io.outm();
                    if !(self.matin <= 0 && self.matout <= 0) {
                        self.mat_totals();
                    }
                }
                self.io.h.mat = -1;
                self.io.h.mf = 0;
                self.io.h.mt = 0;
                self.matin = 0;
                self.matout = 0;
                return;
            }
            match iway {
                1 => self.io.skipm(),
                2 => self.io.skipf(),
                _ => self.io.skips(),
            }
        }
    }

    fn file1(&mut self) {
        let fmttab: [[u8; 4]; 3] = [*b"4   ", *b"5   ", *b"6   "];
        for k in 62..66 {
            self.progdoc[0][k] = b'*';
        }
        let (mut c1, mut c2, mut l1, mut l2, mut n1, mut n2) = self.io.cardi_ld();
        self.iverse = 4;
        let lisox = l2;
        if !(n1 > 0 && n2 == 0) {
            let n2x = n2;
            self.io.cardo(c1, c2, l1, l2, n1, n2);
            (c1, c2, l1, l2, n1, n2) = self.io.cardi_ld();
            self.iverse = 5;
            if n2x > 0 {
                self.io.cardo(c1, c2, l1, l2, n1, n2);
                (c1, c2, l1, l2, n1, n2) = self.io.cardi_ld();
                self.iverse = 6;
                self.liso = lisox;
                self.temp3 = c1;
                l1 = 1;
            }
        }
        let _ = l1;
        self.fmthol = fmttab[(self.iverse - 4) as usize];
        let mut n1out = n1;
        if self.nocoment == 0 {
            if self.maxer3 <= 1 {
                n1out = n1 + 3;
            }
            if self.maxer3 > 1 {
                n1out = n1 + 6 + self.maxer3 as i32;
            }
        }
        self.io.cardo(c1, c2, l1, l2, n1out, n2);
        for _ in 0..n1.max(0) {
            self.io.copy1();
        }
        if self.nocoment == 0 {
            let p = self.progdoc[0];
            self.io.hollyo(&p);
            out9(self.io.err.xclow, &mut self.progdoc[1][26..37]);
            let p = self.progdoc[1];
            self.io.hollyo(&p);
            if self.maxer3 <= 1 {
                let percnt = 100.0 * self.er3[1];
                out9(percnt, &mut self.progdoc[2][41..52]);
                let p = self.progdoc[2];
                self.io.hollyo(&p);
            } else {
                for k in 3..7 {
                    let p = self.progdoc[k];
                    self.io.hollyo(&p);
                }
                for i in 1..=self.maxer3 {
                    let percnt = 100.0 * self.er3[i];
                    out9(self.ener3[i], &mut self.progdoc[7][1..12]);
                    out9(percnt, &mut self.progdoc[7][13..24]);
                    let p = self.progdoc[7];
                    self.io.hollyo(&p);
                }
            }
        }
        self.io.copys();
    }

    fn polynu(e: f64, cn: &[f64], n: usize) -> f64 {
        let mut sum = 0.0;
        let mut epower = 1.0;
        for i in 1..=n {
            sum += cn[i] * epower;
            epower *= e;
        }
        sum
    }

    fn nubar(&mut self) {
        let emin = 1.0e-05;
        let emax = 2.0e+07;
        if self.io.h.l2 == 2 {
            self.io.conto();
            if self.io.h.mt == 455 {
                let (c1, c2, l1, l2, n1, n2) = self.io.cardi_ld();
                self.io.cardo(c1, c2, l1, l2, n1, n2);
                let mut i = 1;
                while i <= n1 {
                    self.io.copy1();
                    i += 6;
                }
            }
            let nothinx = self.nothin;
            self.nothin = 1;
            self.filex();
            self.nothin = nothinx;
            return;
        }
        if self.io.h.mt == 455 {
            self.nubar_hold = self.io.cardi();
            let n1l = self.nubar_hold.4;
            if n1l > 300 {
                return self.nubar_overflow();
            }
            let mut ii = 0;
            let mut i = 1;
            while i <= n1l {
                ii += 1;
                self.lambda[ii] = self.io.hollyi();
                i += 6;
            }
        }
        let (c1, c2, l1, l2, n1, n2) = self.io.cardi_ld();
        if !(self.iverse < 6 || n1 != 1) {
            self.io.conto();
            if self.io.h.mt == 455 {
                self.put_lambda();
            }
            self.io.cardo(c1, c2, l1, l2, n1, n2);
            self.io.copy1();
            self.io.copys();
            let s = self.id_line().s(&" ".repeat(24)).s("Spontaneous");
            self.both(&s);
            return;
        }
        self.io.h.l2 = 2;
        self.io.conto();
        if self.io.h.mt == 455 {
            self.put_lambda();
        }
        let n1u = n1.max(0) as usize;
        let mut xin = std::mem::take(&mut self.xin);
        self.io.listi(&mut xin[1..], n1u);
        self.xin = xin;
        self.n2scr = 0;
        if n1 <= 2 {
            self.n2tot = 2;
            self.xout[1] = emin;
            self.xout[2] = emax;
            self.yout[1] = Self::polynu(self.xout[1], &self.xin, n1u);
            self.yout[2] = Self::polynu(self.xout[2], &self.xin, n1u);
            self.copout();
            return;
        }
        self.n2tot = 1;
        self.xout[1] = emin;
        self.yout[1] = Self::polynu(emin, &self.xin, n1u);
        let mut ebase = emin;
        'l120: loop {
            ebase *= 10.0;
            if ebase > emax {
                ebase = emax;
            }
            let ybase = Self::polynu(ebase, &self.xin, n1u);
            'l130: loop {
                if self.n2tot >= self.npage2 {
                    return self.nubar_overflow();
                }
                let n2m1 = self.n2tot;
                self.n2tot += 1;
                self.xout[self.n2tot] = ebase;
                self.yout[self.n2tot] = ybase;
                loop {
                    // label 140
                    let mut xmid = 0.5 * (self.xout[n2m1] + self.xout[self.n2tot]);
                    incore10(&mut xmid);
                    let ymid = Self::polynu(xmid, &self.xin, n1u);
                    if self.kerr3 != 0 {
                        self.errok3(xmid);
                    }
                    let (y1, y2) = (self.yout[n2m1], self.yout[self.n2tot]);
                    let mut accept = false;
                    if xmid <= self.xout[n2m1] || xmid >= self.xout[self.n2tot] {
                        accept = true;
                    } else if y1.abs() <= self.io.err.xclow && y2.abs() <= self.io.err.xclow {
                        accept = true;
                    } else if (ymid > y1 && ymid > y2) || (ymid < y1 && ymid < y2) {
                        if (ymid - self.nubar_ylin).abs() <= (ymid * self.io.err.erminmax).abs() {
                            accept = true;
                        }
                    } else {
                        self.nubar_ylin = 0.5 * (y1 + y2);
                        if (ymid - self.nubar_ylin).abs() <= (ymid * self.io.err.errxc3).abs() {
                            accept = true;
                        }
                    }
                    if accept {
                        break;
                    }
                    self.xout[self.n2tot] = xmid;
                    self.yout[self.n2tot] = ymid;
                }
                // label 160
                if self.xout[self.n2tot] < ebase {
                    continue 'l130;
                }
                if self.xout[self.n2tot] < emax {
                    continue 'l120;
                }
                break 'l120;
            }
        }
        self.copout();
    }

    fn put_lambda(&mut self) {
        let (c1l, c2l, l1l, l2l, n1l, n2l) = self.nubar_hold;
        self.io.cardo(c1l, c2l, l1l, l2l, n1l, n2l);
        let mut ii = 0;
        let mut i = 1;
        while i <= n1l {
            ii += 1;
            let p = self.lambda[ii];
            self.io.hollyo(&p);
            i += 6;
        }
    }

    fn nubar_overflow(&mut self) {
        let s = self.id_line().s(" Table Overflow - Execution Terminated");
        self.both(&s);
        self.io.enderror();
    }

    fn errok3(&mut self, e: f64) {
        let ener3 = &self.ener3;
        let er3 = &self.er3;
        let maxer3 = self.maxer3;
        if e <= ener3[1] {
            self.miner3 = 2;
            self.io.err.errxc3 = er3[1];
            return;
        }
        let mut nower3 = self.miner3;
        let mut hit = 0; // 1 = label 20, 2 = label 70
        while nower3 <= maxer3 {
            if e < ener3[nower3] {
                hit = 1;
                break;
            }
            if e == ener3[nower3] {
                hit = 2;
                break;
            }
            nower3 += 1;
        }
        if hit == 0 {
            self.miner3 = maxer3;
            self.io.err.errxc3 = er3[maxer3];
            return;
        }
        if hit == 2 {
            self.miner3 = nower3;
            self.io.err.errxc3 = er3[nower3];
            return;
        }
        let mut nm1 = nower3 - 1;
        if e == ener3[nm1] {
            self.miner3 = nm1;
            if self.miner3 <= 1 {
                self.miner3 = 2;
            }
            self.io.err.errxc3 = er3[nm1];
            return;
        }
        if !(e > ener3[nm1]) {
            nower3 = 2;
            let mut found = 0;
            while nower3 <= maxer3 {
                if e < ener3[nower3] {
                    found = 1;
                    break;
                }
                if e == ener3[nower3] {
                    found = 2;
                    break;
                }
                nower3 += 1;
            }
            if found == 0 {
                self.miner3 = maxer3;
                self.io.err.errxc3 = er3[maxer3];
                return;
            }
            if found == 2 {
                self.miner3 = nower3;
                self.io.err.errxc3 = er3[nower3];
                return;
            }
            nm1 = nower3 - 1;
        }
        // label 50
        self.miner3 = nower3;
        self.io.err.errxc3 = ((ener3[nower3] - e) * er3[nm1] + (e - ener3[nm1]) * er3[nower3]) / (ener3[nower3] - ener3[nm1]);
    }

    /// `XYPAGE`: point N of the input table.
    fn xypage(&mut self, n: usize) -> (f64, f64) {
        let n2in = self.io.ld.n2.max(0) as usize;
        if n <= self.ixylow || n > n2in {
            let s = format!(
                " Paging Index Error...Index={}MAT/MF/MT={}{}{}",
                ffmt::i(n as i32, 5),
                ffmt::i(self.io.h.mat, 5),
                ffmt::i(self.io.h.mf, 3),
                ffmt::i(self.io.h.mt, 4)
            );
            self.io.lst(&s);
            return (0.0, 0.0);
        }
        if n > self.ixyhi {
            let mut nn = self.npage;
            if self.left < nn {
                nn = self.left;
            }
            self.left -= nn;
            self.ixylow = self.ixyhi;
            self.ixyhi += nn;
            let mut xin = std::mem::take(&mut self.xin);
            let mut yin = std::mem::take(&mut self.yin);
            self.io.pointi(&mut xin[1..], &mut yin[1..], nn);
            self.io.maxie2(xin[nn], yin[nn]);
            self.xin = xin;
            self.yin = yin;
        }
        let index = n - self.ixylow;
        (self.xin[index], self.yin[index])
    }

    fn filex(&mut self) {
        let (_c1, q, _l1, _l2, n1, n2in) = self.io.cardi_ld();
        let n1u = n1.max(0) as usize;
        let mut nbt = self.nbt;
        let mut int = self.int;
        self.io.terpi(&mut nbt[1..], &mut int[1..], n1u.min(100));
        self.nbt = nbt;
        self.int = int;
        self.io.terptest(&self.nbt[1..], &self.int[1..], n1, n2in, 1);
        let mut im6 = 0;
        for ir in 1..=n1u {
            if self.int[ir] == 6 {
                im6 = 1;
            }
        }
        let (math, mfh, mth) = (self.io.h.mat, self.io.h.mf, self.io.h.mt);
        self.io.maxie1(math, mfh, mth);
        self.istart = 1;
        let mut kstart = 1usize;
        self.ixylow = 0;
        self.ixyhi = 0;
        let n2in_u = n2in.max(0) as usize;
        self.left = n2in_u;
        self.io.elast = 0.0;
        self.mtend = 0;
        self.ithin1 = 1;
        self.ithin2 = 2;
        self.ithin3 = 2;
        self.n2core = 1;
        self.n2p1 = 2;
        self.n2p2 = 3;
        self.n2scr = 0;
        (self.xout[1], self.yout[1]) = self.xypage(1);
        if !(mfh != 3 || self.nothin > 0) {
            let mut go30 = self.yout[1] != 0.0;
            if !go30 {
                loop {
                    kstart += 1;
                    if kstart > n2in_u {
                        (self.xout[1], self.yout[1]) = self.xypage(1);
                        kstart = 1;
                        break;
                    }
                    (self.xout[2], self.yout[2]) = self.xypage(kstart);
                    if self.yout[2] != 0.0 {
                        self.n2core = 2;
                        self.n2p1 = 3;
                        self.n2p2 = 4;
                        go30 = true;
                        break;
                    }
                    self.xout[1] = self.xout[2];
                    self.yout[1] = self.yout[2];
                }
            }
            if go30 && !(self.xout[1] <= self.estart || self.yout[1] == 0.0) {
                self.xout[2] = self.xout[1];
                self.yout[2] = self.yout[1];
                self.yout[1] = 0.0;
                self.n2core = 2;
                self.n2p1 = 3;
                self.n2p2 = 4;
            }
        }
        // label 40
        let mut nr2 = kstart as i64;
        for ir in 1..=n1u {
            let nr1 = nr2 + 1;
            nr2 = self.nbt[ir] as i64;
            if kstart as i64 > nr2 {
                nr2 = kstart as i64;
                continue;
            }
            let intype = self.int[ir];
            if intype == 2 {
                let mut npt = nr1;
                while npt <= nr2 {
                    let (xend, yend) = self.xypage(npt as usize);
                    if !(mfh == 9 || mfh == 10) && !(xend <= ETHERMAL || self.xout[self.n2core] >= ETHERMAL) {
                        let ythermal = terpit(ETHERMAL, xend, self.xout[self.n2core], yend, self.yout[self.n2core], intype);
                        if self.n2core >= self.npage2 {
                            self.thinit();
                        }
                        self.n2core += 1;
                        self.xout[self.n2core] = ETHERMAL;
                        self.yout[self.n2core] = ythermal;
                    }
                    if self.n2core >= self.npage2 {
                        self.thinit();
                    }
                    self.n2core += 1;
                    self.n2p1 = self.n2core + 1;
                    self.n2p2 = self.n2p1 + 1;
                    self.xout[self.n2core] = xend;
                    self.yout[self.n2core] = yend;
                    npt += 1;
                }
                continue;
            }
            if intype < 2 {
                let mut npt = nr1;
                while npt <= nr2 {
                    if self.n2core >= self.np2m1 {
                        self.thinit();
                    }
                    self.n2core += 2;
                    self.n2p1 = self.n2core + 1;
                    self.n2p2 = self.n2p1 + 1;
                    let c = self.n2core;
                    (self.xout[c], self.yout[c]) = self.xypage(npt as usize);
                    if self.xout[c] == self.xout[c - 2] {
                        self.xout[c - 1] = self.xout[c];
                        self.yout[c - 1] = self.yout[c];
                        self.n2core -= 1;
                        self.n2p1 = self.n2core + 1;
                        self.n2p2 = self.n2p1 + 1;
                    } else {
                        self.xout[c - 1] = self.xout[c];
                        self.yout[c - 1] = self.yout[c - 2];
                    }
                    npt += 1;
                }
                continue;
            }
            // intype > 2
            let mut npt = nr1;
            while npt <= nr2 {
                let mut imsaved = 0;
                let (mut xend, mut yend) = self.xypage(npt as usize);
                'l110: loop {
                    let mut xtop = xend;
                    if !(mfh == 9 || mfh == 10) && xtop > ETHERMAL && self.xout[self.n2core] < ETHERMAL {
                        xtop = ETHERMAL;
                    }
                    if xtop > 2.0 * self.xout[self.n2core] {
                        xtop = 2.0 * self.xout[self.n2core];
                    }
                    if xtop < xend {
                        let ytop = terpit(xtop, xend, self.xout[self.n2core], yend, self.yout[self.n2core], intype);
                        imsaved = 1;
                        self.filex_xsave = xend;
                        self.filex_ysave = yend;
                        xend = xtop;
                        yend = ytop;
                    }
                    loop {
                        // label 130
                        let mut xn2p2 = xend;
                        let mut yn2p2 = yend;
                        let xn2 = self.xout[self.n2core];
                        let yn2 = self.yout[self.n2core];
                        let xclow = self.io.err.xclow;
                        let mut go190 = yn2 == yn2p2 || (yn2.abs() <= xclow && yn2p2.abs() <= xclow);
                        let mut xn2p1 = 0.0;
                        let mut yn2p1 = 0.0;
                        if !go190 {
                            loop {
                                // label 140
                                xn2p1 = if intype == 4 { 0.5 * (xn2 + xn2p2) } else { (xn2 * xn2p2).sqrt() };
                                incore10(&mut xn2p1);
                                if xn2p1 <= xn2 || xn2p1 >= xn2p2 {
                                    go190 = true;
                                    break;
                                }
                                if yn2.abs() < xclow && yn2p2.abs() < xclow {
                                    go190 = true;
                                    break;
                                }
                                yn2p1 = terpit(xn2p1, xn2, xn2p2, yn2, yn2p2, intype);
                                if self.kerr3 != 0 {
                                    self.errok3(xn2p1);
                                }
                                let e3 = self.io.err.errxc3;
                                let ok;
                                if (yn2p1 > yn2 && yn2p1 > yn2p2) || (yn2p1 < yn2 && yn2p1 < yn2p2) {
                                    ok = (yn2p1 - self.filex_yaprox).abs() <= (yn2p1 * e3).abs();
                                } else {
                                    self.filex_yaprox = terpit(xn2p1, xn2, xn2p2, yn2, yn2p2, 2);
                                    ok = (yn2p1 - self.filex_yaprox).abs() <= (yn2p1 * e3).abs();
                                }
                                if ok {
                                    break;
                                }
                                xn2p2 = xn2p1;
                                yn2p2 = yn2p1;
                            }
                        }
                        if !go190 && self.nothin > 0 {
                            go190 = true;
                        }
                        if !go190 {
                            // label 180
                            if self.n2core >= self.np2m1 {
                                self.thinit();
                            }
                            let (p1, p2) = (self.n2p1, self.n2p2);
                            self.xout[p1] = xn2p1;
                            self.yout[p1] = yn2p1;
                            self.xout[p2] = xn2p2;
                            self.yout[p2] = yn2p2;
                            self.n2core = self.n2p2;
                            self.n2p1 = self.n2core + 1;
                            self.n2p2 = self.n2p1 + 1;
                        } else {
                            if self.n2core >= self.npage2 {
                                self.thinit();
                            }
                            let p1 = self.n2p1;
                            self.xout[p1] = xn2p2;
                            self.yout[p1] = yn2p2;
                            self.n2core = self.n2p1;
                            self.n2p1 = self.n2core + 1;
                            self.n2p2 = self.n2p1 + 1;
                        }
                        // label 200
                        if self.xout[self.n2core] < xend {
                            continue;
                        }
                        break;
                    }
                    if imsaved <= 0 {
                        break 'l110;
                    }
                    imsaved = 0;
                    xend = self.filex_xsave;
                    yend = self.filex_ysave;
                }
                npt += 1;
            }
        }
        self.mtend = 1;
        self.thinit();
        self.copout();
        if im6 > 0 && q < 0.0 {
            self.both(" WARNING - the Above Section Uses Interpolation Law 6\n but, Has a Negative Q-Value? Results Will be Unreliable.");
        }
    }

    fn write_page(&mut self, which: usize) {
        let base = if which == 1 { 0 } else { self.npage };
        let mut rec = Vec::with_capacity(2 * self.npage);
        for i in 1..=self.npage {
            rec.push(Val::D(self.xout[base + i]));
        }
        for i in 1..=self.npage {
            rec.push(Val::D(self.yout[base + i]));
        }
        let u = self.iscr;
        self.io.u.uwrite(u, rec);
    }

    fn read_page(&mut self, which: usize) {
        let u = self.iscr;
        let rec = self.io.u.uread_req(u);
        let base = if which == 1 { 0 } else { self.npage };
        let n = self.npage;
        for i in 1..=n {
            self.xout[base + i] = rec[i - 1].d();
            self.yout[base + i] = rec[n + i - 1].d();
        }
    }

    fn monitor_line(&mut self, n: usize, es1: f64, es2: f64) {
        out9(es1, &mut self.field[0]);
        out9(es2, &mut self.field[1]);
        let s = L::new()
            .s(&" ".repeat(34))
            .s(&ffmt::i(n as i32, 8))
            .b(&self.field[0])
            .s(" to")
            .b(&self.field[1])
            .s(" eV Finished");
        self.both(&s);
    }

    fn thinit(&mut self) {
        let mfh = self.io.h.mf;
        'l160: {
            if self.nothin != 0 {
                self.ithin1 = self.n2core;
                self.ithin3 = self.n2core + 1;
                break 'l160;
            }
            if self.ithin2 > self.n2core {
                break 'l160;
            }
            if self.istart > 0 {
                self.istart = 0;
                self.th_dsign = 1.0;
                if self.yout[1] > self.yout[2] {
                    self.th_dsign = -1.0;
                }
            }
            let mut m = self.ithin2;
            while m <= self.n2core {
                let mm1 = m - 1;
                let mut dxc = self.yout[m] - self.yout[mm1];
                #[derive(Clone, Copy, PartialEq)]
                enum G {
                    L30,
                    L40,
                    L50,
                    L80,
                    L110,
                    L120,
                    L130,
                    L140,
                    L150,
                }
                let mut g = G::L30;
                if !(dxc * self.th_dsign >= 0.0) {
                    self.th_dsign = -self.th_dsign;
                    if self.xout[m] <= self.xout[mm1] {
                        g = G::L50;
                    } else if m > self.ithin3 {
                        g = G::L110;
                    }
                }
                loop {
                    match g {
                        G::L30 => {
                            if self.yout[m] * self.yout[mm1] >= 0.0 {
                                g = G::L40;
                            } else if m <= self.ithin3 {
                                g = G::L130;
                            } else {
                                g = G::L120;
                            }
                        }
                        G::L40 => {
                            g = if self.xout[m] > self.xout[mm1] { G::L80 } else { G::L50 };
                        }
                        G::L50 => {
                            dxc = dxc.abs();
                            if m != self.ithin3 {
                                g = if dxc > 0.0 { G::L120 } else if m <= self.ithin3 { G::L130 } else { G::L120 };
                            } else if dxc <= 0.0 {
                                g = G::L140;
                            } else {
                                g = if m <= self.ithin3 { G::L130 } else { G::L120 };
                            }
                        }
                        G::L80 => {
                            if !(mfh == 9 || mfh == 10) && self.xout[m] == ETHERMAL {
                                g = if m <= self.ithin3 { G::L130 } else { G::L120 };
                                continue;
                            }
                            let dxn = self.xout[m] - self.xout[self.ithin1];
                            let slope = (self.yout[m] - self.yout[self.ithin1]) / dxn;
                            if m == self.ithin3 {
                                if self.kerr3 != 0 {
                                    let x = self.xout[m];
                                    self.errok3(x);
                                }
                                let dslope = self.io.err.errxc3 * self.yout[m] / dxn;
                                self.th_slpmax = slope + dslope;
                                self.th_slpmin = slope - dslope;
                                g = G::L150;
                                continue;
                            }
                            if slope > self.th_slpmax || slope < self.th_slpmin {
                                g = G::L110;
                                continue;
                            }
                            if self.kerr3 != 0 {
                                let x = self.xout[m];
                                self.errok3(x);
                            }
                            let dslope = self.io.err.errxc3 * self.yout[m] / dxn;
                            let slp1 = slope + dslope;
                            if slp1 < self.th_slpmax {
                                self.th_slpmax = slp1;
                            }
                            let slp2 = slope - dslope;
                            if slp2 > self.th_slpmin {
                                self.th_slpmin = slp2;
                            }
                            g = G::L150;
                        }
                        G::L110 => {
                            self.ithin1 += 1;
                            let i1 = self.ithin1;
                            self.xout[i1] = self.xout[mm1];
                            self.yout[i1] = self.yout[mm1];
                            self.ithin3 = m;
                            g = G::L80;
                        }
                        G::L120 => {
                            self.ithin1 += 1;
                            let i1 = self.ithin1;
                            self.xout[i1] = self.xout[mm1];
                            self.yout[i1] = self.yout[mm1];
                            g = G::L130;
                        }
                        G::L130 => {
                            self.ithin1 += 1;
                            let i1 = self.ithin1;
                            self.xout[i1] = self.xout[m];
                            self.yout[i1] = self.yout[m];
                            g = G::L140;
                        }
                        G::L140 => {
                            self.ithin3 = m + 1;
                            g = G::L150;
                        }
                        G::L150 => break,
                    }
                }
                m += 1;
            }
        }
        // label 160
        if self.mtend <= 0 {
            let n2old = self.n2core;
            self.n2core = self.ithin1;
            self.n2p1 = self.n2core + 1;
            self.n2p2 = self.n2p1 + 1;
            if self.ithin3 <= n2old {
                let n2p = self.ithin1 + 1;
                self.xout[n2p] = self.xout[n2old];
                self.yout[n2p] = self.yout[n2old];
                self.n2core = n2p;
                self.n2p1 = self.n2core + 1;
                self.n2p2 = self.n2p1 + 1;
            }
            if self.ithin1 > self.np1p1 {
                if self.n2scr == 0 && self.th_imused > 0 {
                    let u = self.iscr;
                    self.io.u.rewind(u);
                }
                self.write_page(1);
                self.th_imused = 1;
                self.n2scr += self.npage;
                if self.monitr > 0 {
                    let (es1, es2) = (self.xout[1], self.xout[self.npage]);
                    let n = self.n2scr;
                    self.monitor_line(n, es1, es2);
                }
                let mut ll = 0;
                for lll in self.np1p1..=self.n2core {
                    ll += 1;
                    self.xout[ll] = self.xout[lll];
                    self.yout[ll] = self.yout[lll];
                }
                self.ithin1 -= self.npage;
                self.n2core -= self.npage;
                self.n2p1 = self.n2core + 1;
                self.n2p2 = self.n2core + 2;
            }
            self.ithin2 = self.n2core + 1;
            if self.ithin3 > n2old {
                self.ithin3 = self.ithin2;
            } else {
                self.ithin3 = 0;
            }
            return;
        }
        // label 230
        if self.ithin3 <= self.n2core {
            self.ithin1 += 1;
            let (i1, c) = (self.ithin1, self.n2core);
            self.xout[i1] = self.xout[c];
            self.yout[i1] = self.yout[c];
        }
        if self.n2scr > 0 {
            self.write_page(1);
            if self.monitr > 0 {
                let npout = if self.ithin1 < self.npage { self.ithin1 } else { self.npage };
                let (es1, es2) = (self.xout[1], self.xout[npout]);
                let n = self.n2scr + npout;
                self.monitor_line(n, es1, es2);
            }
            if self.ithin1 > self.npage {
                self.write_page(2);
                if self.monitr > 0 {
                    let (es1, es2) = (self.xout[self.np1p1], self.xout[self.ithin1]);
                    let n = self.n2scr + self.ithin1;
                    self.monitor_line(n, es1, es2);
                }
            }
            self.n2scr += self.ithin1;
            self.ithin1 = 0;
            let u = self.iscr;
            self.io.u.endfile(u);
            self.io.u.rewind(u);
        }
        self.n2tot = self.n2scr + self.ithin1;
    }

    fn copout(&mut self) {
        let ld = self.io.ld;
        let (c1, q, l1, l2, n2in) = (ld.c1, ld.c2, ld.l1, ld.l2, ld.n2);
        self.io.cardo(c1, q, l1, l2, 1, self.n2tot as i32);
        self.io.terpo(&[self.n2tot as i32], &[2], 1);
        let mfh = self.io.h.mf;
        let mth = self.io.h.mt;
        if mfh != 3 {
            let s = self.id_line().s(&" ".repeat(24)).s(&ffmt::i(n2in, 9)).s(&ffmt::i(self.n2tot as i32, 9));
            self.both(&s);
        } else {
            let temp;
            if self.iverse >= 6 || l2 > 0 {
                temp = self.temp3;
            } else {
                temp = c1;
                if mth == 1 {
                    self.temp3 = temp;
                }
            }
            out9(temp, &mut self.field[0]);
            out9(q, &mut self.field[1]);
            let s = self
                .id_line()
                .s(" ")
                .b(&self.field[0])
                .s(" ")
                .b(&self.field[1])
                .s(&ffmt::i(n2in, 9))
                .s(&ffmt::i(self.n2tot as i32, 9));
            self.both(&s);
        }
        self.matin += n2in;
        self.matout += self.n2tot as i32;
        self.io.minus3 = 0;
        self.io.implus = 0;
        if self.n2scr == 0 {
            let x = std::mem::take(&mut self.xout);
            let y = std::mem::take(&mut self.yout);
            self.io.pointo(&x[1..], &y[1..], self.n2tot);
            self.xout = x;
            self.yout = y;
        } else {
            let mut ithinx = self.npage2 as i64;
            let mut lp = 1usize;
            // N2SCR goes negative in the Fortran; it is set again before
            // any later use.
            let mut n2scr = self.n2scr as i64;
            loop {
                self.read_page(1);
                if n2scr > self.npage as i64 {
                    self.read_page(2);
                }
                if ithinx > n2scr {
                    ithinx = n2scr;
                }
                n2scr -= self.npage2 as i64;
                self.n2scr = n2scr.max(0) as usize;
                let x = std::mem::take(&mut self.xout);
                let y = std::mem::take(&mut self.yout);
                self.io.pointo(&x[1..], &y[1..], ithinx.max(0) as usize);
                self.xout = x;
                self.yout = y;
                lp += self.npage2;
                if lp > self.n2tot {
                    break;
                }
            }
        }
        if mfh != 9 && mfh != 10 {
            let m = self.io.h.mat;
            self.io.outs(m, mfh);
        }
        let col = " ".repeat(21);
        if self.io.minus3 > 0 {
            let s = format!("{col}WARNING - Above Cross Section Negative at{} Energies", ffmt::i(self.io.minus3, 6));
            self.both(&s);
        }
        if self.io.implus <= 0 {
            let s = format!("{col}WARNING - Above Cross Sections NOT Positive at Any Energy");
            self.both(&s);
        }
    }
}

pub fn run() -> ! {
    let mut s = Linear::new();
    s.io.timer();
    s.fileio();
    s.npage = MAXPAGE;
    s.npage2 = 2 * MAXPAGE;
    s.np1p1 = s.npage + 1;
    s.np2m1 = s.npage2 - 1;
    s.estart = 1.0e9;
    let d = Linear::dash();
    s.both(format!("  Linearize ENDF/B Cross Sections (LINEAR 2026-1)\n{d}"));
    s.readin();
    s.io.copyl();
    let l = L::new()
        .s(&d)
        .s("\n  ENDF/B Tape Label\n")
        .s(&d)
        .s("\n  ")
        .b(&s.io.card)
        .s(&ffmt::i(s.io.mfield[0], 4))
        .s("\n")
        .s(&d)
        .s("\n    Material    MAT  MF  MT  ENDF/B      Kelvin     Q-Value   Points   Points\n")
        .s(&" ".repeat(28))
        .s(" Format")
        .s(&" ".repeat(22))
        .s("eV       In      Out\n")
        .s(&d);
    s.both(&l);
    s.n2toti = 0;
    s.n2toto = 0;
    s.matnow = 0;
    s.mfnow = 0;
    loop {
        s.nextmt();
        let (math, mfh, mth) = (s.io.h.mat, s.io.h.mf, s.io.h.mt);
        if mfh == 1 && (mth == 452 || mth == 455 || mth == 456) {
            s.nubar();
            continue;
        }
        s.io.conto();
        if math <= 0 {
            break;
        }
        if mfh == 1 && mth == 451 {
            s.file1();
            continue;
        }
        if mfh == 3 || mfh == 23 {
            s.filex();
            continue;
        }
        if mfh == 9 || mfh == 10 {
            let ns = s.io.h.n1;
            for _ in 0..ns.max(0) {
                s.filex();
            }
            s.io.outs(math, mfh);
            continue;
        }
        if mfh == 27 {
            let nothinx = s.nothin;
            if mth == 506 {
                s.nothin = 1;
            }
            s.filex();
            s.nothin = nothinx;
            continue;
        }
        s.io.copys();
    }
    if s.matnow <= 0 {
        s.both(format!(
            "  WARNING - No Data Found That Satisfied Retrieval Criteria.\n{} Therefore No Data was Linearized or Written to Output File.\n{d}",
            " ".repeat(12)
        ));
        s.io.enderror();
    }
    let t = format!("{d}\n{}Tape Totals{}{}\n{d}", " ".repeat(48), ffmt::i(s.n2toti, 9), ffmt::i(s.n2toto, 9));
    s.both(&t);
    s.io.maxie4(1);
    let _ = f9s;
    s.io.endit();
}
