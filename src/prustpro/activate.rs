//! Port of ACTIVATE (PREPRO 2026-1): make MF=10 activation cross sections
//! from MF=3 cross sections and MF=9 multiplicities. Input `ACTIVATE.INP`,
//! listing `ACTIVATE.LST`.

use crate::endfio::{ctl, dp, fint, in9, out9, terpit, zaholm, Endf};
use crate::ffmt::i;
use crate::fio::{self, Rd, L};
use crate::linear::blank;

const MAXMT: usize = 999;
const MAXPOINT: usize = 120000;
const ISCR1: i32 = 12;
const ISCR2: i32 = 14;
const ISCR3: i32 = 15;
const ISCR4: i32 = 16;

fn dash() -> String {
    format!(" {}", "-".repeat(78))
}

struct Act {
    io: Endf,
    mttab3: Vec<i32>,
    mttab9: Vec<i32>,
    mttab10: Vec<i32>,
    x9: Vec<f64>,
    y9: Vec<f64>,
    n9: usize,
    zabcd: [u8; 12],
    zabcd2: [u8; 12],
    izanow: i32,
    matnow: i32,
    n2mati: i32,
    n2mato: i32,
    n2toti: i32,
    n2toto: i32,
    liso: i32,
    lstiza: i32,
    lstmat: i32,
    lstmf: i32,
    line66: [u8; 66],
    nbt: [i32; 101],
    int: [i32; 101],
    fmthol: [u8; 4],
    field: [[u8; 11]; 6],
    mf3: i32,
    // FILE1 SAVEd IVERSE
    f1_iverse: i32,
}

impl Act {
    /// READ(u,'(66A1,I4,I2,I3,I5)') LINE66,MATH,MFH,MTH; None at EOF.
    fn rd210(&mut self, u: i32) -> Option<()> {
        match self.io.u.read(u) {
            Rd::Eof => None,
            Rd::Line(l) => {
                let l = fio::pad(l, 75);
                self.line66.copy_from_slice(&l[..66]);
                self.io.h.mat = fio::rd_i(&l, 66, 4, u);
                self.io.h.mf = fio::rd_i(&l, 70, 2, u);
                self.io.h.mt = fio::rd_i(&l, 72, 3, u);
                Some(())
            }
        }
    }

    fn rd210_req(&mut self, u: i32) {
        if self.rd210(u).is_none() {
            fio::runtime_error(u, "", "End of file");
        }
    }

    fn w210(&mut self, u: i32) {
        let l = L::new().b(&self.line66).s(&ctl(self.io.h.mat, self.io.h.mf, self.io.h.mt));
        self.io.u.write(u, l);
    }

    fn report(&mut self) {
        if self.n2mati > 0 || self.n2mato > 0 {
            let s = format!("{d}\n{} MAT Totals{}{}\n{d}", " ".repeat(45), i(self.n2mati, 7), i(self.n2mato, 7), d = dash());
            self.io.both(&s);
        } else {
            let s = L::new().s(" ").b(&self.zabcd2).s(&i(self.lstmat, 5)).s(&" ".repeat(8)).s(" WARNING - No MF=10 to create - will copy MAT");
            self.io.both(&s);
        }
        self.io.maxie3(1);
        self.io.timemat();
    }

    fn syncend(&mut self) {
        if self.io.otape <= 0 {
            return;
        }
        if self.lstmat > 0 && !(self.io.h.mat == self.lstmat && self.io.h.mf == self.lstmf) {
            let m = self.lstmat;
            self.io.outf(m);
            if self.io.h.mat != self.lstmat {
                self.report();
                self.n2mati = 0;
                self.n2mato = 0;
                self.io.outm();
            }
        }
        self.lstiza = fint(self.io.h.c1);
        self.lstmat = self.io.h.mat;
        self.lstmf = self.io.h.mf;
        self.zabcd2 = zaholm(self.lstiza, self.liso);
    }

    fn contog(&mut self) {
        self.syncend();
        self.io.conto();
    }

    fn nextmt(&mut self) {
        loop {
            self.io.conti();
            if self.io.h.mt > 0 {
                if self.matnow != self.io.h.mat {
                    self.fmthol = *b"6   ";
                }
                self.izanow = fint(self.io.h.c1);
                self.zabcd = zaholm(self.izanow, self.liso);
                self.matnow = self.io.h.mat;
                return;
            }
            if self.io.h.mat > 0 {
                continue;
            }
            if self.io.h.mat < 0 {
                self.io.h.mat = -1;
                self.io.h.mf = 0;
                self.io.h.mt = 0;
                self.n2mati = 0;
                self.n2mato = 0;
            }
            return;
        }
    }

    fn file1(&mut self) {
        let fmttab: [[u8; 4]; 3] = [*b"4   ", *b"5   ", *b"6   "];
        let texts = [
            " *************** PROGRAM ACTIVATE (VERSION 2026-1) ***********",
            " MF=10 Activation Cross Sections Defined by Combining MF=3    ",
            " Cross Section and MF=9 Multiplier. MF=9 Deleted/MF=8 Updated.",
        ];
        let mut pd = [[b' '; 66]; 3];
        for (k, s) in texts.iter().enumerate() {
            pd[k][..s.len()].copy_from_slice(s.as_bytes());
        }
        for k in 62..66 {
            pd[0][k] = b'*';
        }
        let (mut c1, mut c2, mut l1, mut l2, mut n1, mut n2) = self.io.cardi_ld();
        self.f1_iverse = 4;
        let lisox = l2;
        if !(n1 > 0 && n2 == 0) {
            let n2x = n2;
            self.io.cardo(c1, c2, l1, l2, n1, n2);
            (c1, c2, l1, l2, n1, n2) = self.io.cardi_ld();
            self.f1_iverse = 5;
            if n2x > 0 {
                self.io.cardo(c1, c2, l1, l2, n1, n2);
                (c1, c2, _, l2, n1, n2) = self.io.cardi_ld();
                self.f1_iverse = 6;
                self.liso = lisox;
                l1 = 1;
            }
        }
        self.fmthol = fmttab[(self.f1_iverse - 4) as usize];
        self.io.cardo(c1, c2, l1, l2, n1 + 3, n2);
        for _ in 0..n1.max(0) {
            self.io.copy1();
        }
        for p in &pd {
            self.io.hollyo(p);
        }
        self.io.copys();
    }

    fn file3(&mut self) {
        self.matnow = self.io.h.mat;
        self.io.u.rewind(ISCR2);
        let io1 = self.io.otape;
        self.contog();
        self.io.otape = ISCR2;
        self.io.conto();
        self.io.otape = io1;
        loop {
            self.rd210_req(ISCR1);
            if self.io.h.mf <= 0 {
                break;
            }
            if self.io.h.mt <= 0 {
                let (m, f) = (self.io.h.mat, self.io.h.mf);
                self.io.outs(m, f);
            } else {
                let o = self.io.otape;
                self.w210(o);
            }
            self.w210(ISCR2);
        }
        self.io.u.endfile(ISCR2);
        self.io.u.rewind(ISCR2);
    }

    fn copy_to(&mut self, u: i32) {
        self.io.u.rewind(u);
        let io1 = self.io.otape;
        self.io.otape = u;
        self.io.conto();
        self.io.otape = io1;
        loop {
            self.rd210_req(ISCR1);
            self.w210(u);
            if self.io.h.mf <= 0 {
                break;
            }
        }
        self.io.u.endfile(u);
        self.io.u.rewind(u);
    }

    fn intlaw(&mut self) {
        let (_c1, _c2, _l1, _l2, n1, n2) = self.io.cardi_ld();
        let n1u = n1.max(0) as usize;
        let mut nbt = self.nbt;
        let mut int = self.int;
        self.io.terpi(&mut nbt[1..], &mut int[1..], n1u.min(100));
        self.nbt = nbt;
        self.int = int;
        self.io.terptest(&self.nbt[1..], &self.int[1..], n1, n2, 3);
    }

    fn filex(&mut self) {
        if self.mf3 == 0 {
            self.mf3 = 3;
        }
        let kza = fint(self.io.h.c1);
        self.zabcd = zaholm(kza, self.liso);
        self.io.h.mf = 10;
        let mtin = self.io.h.mt;
        let n2in = self.io.ld.n2;
        let n2out = n2in;
        let (m, t) = (self.io.h.mat, self.io.h.mt);
        self.io.maxie1(m, 3, t);
        let mut i9 = 1usize;
        let mut enow = 0.0;
        let mut xcnow = 0.0;
        let n9 = self.n9;
        let mut i1 = 1;
        while i1 <= n2in {
            let mut i2 = i1 + 2;
            if i2 > n2in {
                i2 = n2in;
            }
            let k2 = (i2 - i1) + 1;
            match self.io.u.read(ISCR2) {
                Rd::Eof => {
                    self.io.both(" ERROR - END of ENDF Data while Reading MF=3 Data");
                    self.io.enderror();
                }
                Rd::Line(l) => {
                    let l = fio::pad(l, 75);
                    self.line66.copy_from_slice(&l[..66]);
                    self.io.h.mat = fio::rd_i(&l, 66, 4, ISCR2);
                    self.mf3 = fio::rd_i(&l, 70, 2, ISCR2);
                    self.io.h.mt = fio::rd_i(&l, 72, 3, ISCR2);
                }
            }
            let mut kk = 0usize;
            for _ in 0..k2 {
                enow = in9(&self.line66[kk..kk + 11]);
                xcnow = in9(&self.line66[kk + 11..kk + 22]);
                let mut k9 = i9;
                let mut hit = 0;
                while k9 <= n9 {
                    if enow < self.x9[k9] {
                        hit = 1;
                        break;
                    }
                    if enow == self.x9[k9] {
                        hit = 2;
                        break;
                    }
                    k9 += 1;
                }
                let xmult;
                if hit == 1 && k9 > 1 {
                    xmult = terpit(enow, self.x9[k9 - 1], self.x9[k9], self.y9[k9 - 1], self.y9[k9], 2);
                } else if hit == 2 {
                    xmult = self.y9[k9];
                } else {
                    let s = format!(
                        " WARNING - MF3 E={} Not in Energy Range of\n           MF9 E={}{} eV Data.\n           Will Assume MF9 Data = 0 at E={}",
                        dp(enow, 12, 5),
                        dp(self.x9[1], 12, 5),
                        dp(self.x9[n9], 12, 5),
                        dp(enow, 12, 5)
                    );
                    self.io.both(&s);
                    xmult = 0.0;
                    k9 = 1;
                }
                i9 = k9;
                xcnow *= xmult;
                out9(xcnow, &mut self.line66[kk + 11..kk + 22]);
                kk += 22;
            }
            let o = self.io.otape;
            self.w210(o);
            i1 += 3;
        }
        self.io.maxie2(enow, xcnow);
        let c2 = self.io.ld.c2;
        out9(c2, &mut self.field[0]);
        let s = L::new()
            .s(" ")
            .b(&self.zabcd)
            .s(&i(self.io.h.mat, 5))
            .s(&i(self.io.h.mf, 4))
            .s(&i(mtin, 4))
            .s("   ")
            .b(&self.fmthol[..2])
            .s(&i(self.io.ld.l1, 8))
            .s(&i(self.io.ld.l2, 6))
            .s("  ")
            .b(&self.field[0])
            .s(&i(n2in, 7))
            .s(&i(n2out, 7));
        self.io.both(&s);
        self.n2mati += n2in;
        self.n2mato += n2out;
        self.n2toti += n2in;
        self.n2toto += n2out;
    }

    fn create10(&mut self) {
        self.io.u.rewind(ISCR2);
        self.io.u.rewind(ISCR3);
        self.io.u.rewind(ISCR4);
        let isave = self.io.itape;
        self.io.itape = ISCR3;
        for imt in 1..=MAXMT {
            if self.mttab9[imt] > 0 {
                if self.mttab3[imt] <= 0 {
                    self.io.both(format!(" ERROR - MT={} MF=9 Data Required MF=3 Data", i(imt as i32, 5)));
                    self.io.enderror();
                }
                self.io.itape = ISCR3;
                loop {
                    self.io.conti();
                    if self.io.h.mt <= 0 {
                        continue;
                    }
                    if self.io.h.mt == imt as i32 {
                        break;
                    }
                    if self.io.h.mt > imt as i32 {
                        self.io.both(format!(" ERROR - Cannot Locate MF=  9, MT={} Data", i(imt as i32, 5)));
                        self.io.enderror();
                    }
                    self.io.skips();
                }
                self.io.h.mf = 10;
                self.contog();
                let ns = self.io.h.n1;
                for _ in 0..ns.max(0) {
                    self.io.itape = ISCR3;
                    self.intlaw();
                    let n9 = self.io.ld.n2;
                    if n9 < 2 || n9 > MAXPOINT as i32 {
                        self.io.both(format!(" ERROR - N2={} (MUST be 2 to {})\n Execution Terminated", i(n9, 8), i(MAXPOINT as i32, 8)));
                        self.io.enderror();
                    }
                    self.n9 = n9 as usize;
                    self.io.elast = 0.0;
                    let mut x = std::mem::take(&mut self.x9);
                    let mut y = std::mem::take(&mut self.y9);
                    self.io.pointi(&mut x[1..], &mut y[1..], self.n9);
                    self.x9 = x;
                    self.y9 = y;
                    let (m, t) = (self.io.h.mat, self.io.h.mt);
                    self.io.maxie1(m, 9, t);
                    let (xm, ym) = (self.x9[self.n9], self.y9[self.n9]);
                    self.io.maxie2(xm, ym);
                    let ld = self.io.ld;
                    let (c19, c29, l19, l29) = (ld.c1, ld.c2, ld.l1, ld.l2);
                    self.io.itape = ISCR2;
                    self.io.u.rewind(ISCR2);
                    loop {
                        self.io.conti();
                        if self.io.h.mt <= 0 {
                            continue;
                        }
                        if self.io.h.mt == imt as i32 {
                            break;
                        }
                        if self.io.h.mt > imt as i32 {
                            self.io.both(format!(" ERROR - Cannot Locate MF=  3, MT={} Data", i(imt as i32, 5)));
                            self.io.enderror();
                        }
                        self.io.skips();
                    }
                    self.intlaw();
                    self.io.ld.c1 = c19;
                    self.io.ld.c2 = c29;
                    self.io.ld.l1 = l19;
                    self.io.ld.l2 = l29;
                    self.io.h.mf = 10;
                    let ld = self.io.ld;
                    self.io.cardo(ld.c1, ld.c2, ld.l1, ld.l2, ld.n1, ld.n2);
                    let n1u = ld.n1.max(0) as usize;
                    let (nbt, int) = (self.nbt, self.int);
                    self.io.terpo(&nbt[1..], &int[1..], n1u.min(100));
                    self.filex();
                }
                let (m, f) = (self.io.h.mat, self.io.h.mf);
                self.io.outs(m, f);
                continue;
            }
            if self.mttab10[imt] <= 0 {
                continue;
            }
            self.io.h.mf = 10;
            self.syncend();
            let cant = |a: &mut Act| -> ! {
                a.io.both(format!(" ERROR - Cannot Locate MF= 10, MT={} Data", i(imt as i32, 5)));
                a.io.enderror();
            };
            loop {
                if self.rd210(ISCR4).is_none() {
                    cant(self);
                }
                if self.io.h.mt == imt as i32 {
                    break;
                }
                if self.io.h.mt > imt as i32 {
                    cant(self);
                }
                // the Fortran skips on ITAPE (ISCR3 here), not ISCR4
                self.io.skips();
            }
            loop {
                if self.io.h.mt <= 0 {
                    let (m, f) = (self.io.h.mat, self.io.h.mf);
                    self.io.outs(m, f);
                    break;
                }
                let o = self.io.otape;
                self.w210(o);
                if self.rd210(ISCR4).is_none() {
                    cant(self);
                }
            }
        }
        self.io.itape = isave;
    }
}

pub fn run() -> ! {
    let mut io = Endf::new("ACTIVATE");
    io.timer();
    io.inp = 2;
    io.outp = 3;
    io.itape = 10;
    io.otape = 11;
    io.u.open(3, "ACTIVATE.LST");
    for u in [ISCR1, ISCR2, ISCR3, ISCR4] {
        io.u.open_fscratch(u);
    }
    let istat1 = if io.u.open_old(2, "ACTIVATE.INP") { 0 } else { 1 };
    let isave = io.itape;
    io.both(format!(" ENDF Activation Cross Sections (ACTIVATE 2026-1)\n{}", dash()));
    // READIN
    let mut namein = None;
    let mut nameout = None;
    if istat1 != 1 {
        if let Rd::Line(l) = io.u.read(2) {
            let mut n = fio::rd_a(&l, 0, 72);
            if blank(&n) {
                n = fio::rd_a(b"ENDFB.IN", 0, 72);
            }
            namein = Some(n);
            if let Rd::Line(l) = io.u.read(2) {
                let mut n = fio::rd_a(&l, 0, 72);
                if blank(&n) {
                    n = fio::rd_a(b"ENDFB.OUT", 0, 72);
                }
                nameout = Some(n);
            }
        }
    }
    let namein = namein.unwrap_or_else(|| fio::rd_a(b"ENDFB.IN", 0, 72));
    let nameout = nameout.unwrap_or_else(|| fio::rd_a(b"ENDFB.OUT", 0, 72));
    let s = L::new().s(" ENDF Input and Output Data Filenames\n").s(&dash()).s("\n ").b(&namein).s("\n ").b(&nameout);
    io.both(&s);
    let ni = String::from_utf8_lossy(&namein).into_owned();
    if !io.u.open_old(10, &ni) {
        let s = L::new().s("\n\n ERROR - Opening ENDF data file\n ").b(&namein).s("\n\n");
        io.both(&s);
        io.enderror();
    }
    let no = String::from_utf8_lossy(&nameout).into_owned();
    io.u.open(11, &no);
    io.copyl();
    let s = L::new()
        .s(&dash())
        .s("\n ENDF Tape Label\n")
        .s(&dash())
        .s("\n ")
        .b(&io.card)
        .s(&i(io.mfield[0], 4))
        .s("\n")
        .s(&dash())
        .s("\n   Material    MAT  MF  MT  ENDF    ZAP Final      Q-Value Points Points\n                            Format      State           eV     In    Out\n")
        .s(&dash());
    io.both(&s);
    let mut a = Act {
        io,
        mttab3: vec![0; MAXMT + 1],
        mttab9: vec![0; MAXMT + 1],
        mttab10: vec![0; MAXMT + 1],
        x9: vec![0.0; MAXPOINT + 2],
        y9: vec![0.0; MAXPOINT + 2],
        n9: 0,
        zabcd: [0u8; 12],
        zabcd2: [0u8; 12],
        izanow: 0,
        matnow: 0,
        n2mati: 0,
        n2mato: 0,
        n2toti: 0,
        n2toto: 0,
        liso: 0,
        lstiza: 0,
        lstmat: -9999,
        lstmf: -99,
        line66: [b' '; 66],
        nbt: [0; 101],
        int: [0; 101],
        fmthol: [0u8; 4],
        field: [[0u8; 11]; 6],
        mf3: 0,
        f1_iverse: 6,
    };
    let mut imall10 = 0;
    let mt9: &[u8; 11] = b"          9";
    let mt10: &[u8; 11] = b"         10";
    'l10: loop {
        let mut imfile9 = 0;
        a.n2mati = 0;
        a.n2mato = 0;
        a.liso = 0;
        for k in 1..=MAXMT {
            a.mttab3[k] = 0;
            a.mttab9[k] = 0;
            a.mttab10[k] = 0;
        }
        a.io.u.rewind(ISCR1);
        a.io.itape = isave;
        let mut lisin = 0;
        a.rd210_req(isave);
        if a.io.h.mat < 0 {
            break 'l10;
        }
        a.matnow = a.io.h.mat;
        let zax = in9(&a.line66[..11]);
        a.izanow = fint(zax);
        a.zabcd = zaholm(a.izanow, a.liso);
        a.io.maxie0();
        loop {
            // label 20
            if a.io.h.mt > 0 {
                let (mf, mt) = (a.io.h.mf, a.io.h.mt as usize);
                if mf == 3 {
                    a.mttab3[mt.min(MAXMT)] = 1;
                }
                if mf == 9 {
                    imfile9 += 1;
                    a.mttab9[mt.min(MAXMT)] = 1;
                }
                if mf == 10 {
                    a.mttab10[mt.min(MAXMT)] = 1;
                }
                if mf == 8 && a.line66[22..33] == mt9[..] {
                    a.line66[22..33].copy_from_slice(mt10);
                }
            }
            // label 30
            a.w210(ISCR1);
            if a.io.h.mt == 451 {
                lisin += 1;
            }
            if lisin == 2 {
                let xliso = in9(&a.line66[33..44]);
                a.liso = fint(xliso);
                a.zabcd = zaholm(a.izanow, a.liso);
            }
            a.rd210_req(isave);
            if a.io.h.mat > 0 {
                continue;
            }
            break;
        }
        a.w210(ISCR1);
        a.io.u.write(ISCR1, format!("{}  -1 0  0    0", " ".repeat(66)));
        a.io.u.endfile(ISCR1);
        a.io.u.rewind(ISCR1);
        imfile9 = 0;
        let mut imhere = 0;
        for k in 1..=MAXMT {
            if a.mttab10[k] != 0 {
                let s = L::new()
                    .s(" ")
                    .b(&a.zabcd)
                    .s(&i(a.matnow, 5))
                    .s(&i(10, 4))
                    .s(&i(k as i32, 4))
                    .s(" WARNING - MF=10 data already exists - will copy MT.");
                a.io.both(&s);
                a.mttab9[k] = 0;
                a.mttab3[k] = 0;
                imhere = 1;
                imall10 = 1;
            } else if a.mttab9[k] > 0 {
                imfile9 += 1;
            }
        }
        if imhere > 0 && imfile9 > 0 {
            a.io.both(dash());
        }
        if imfile9 == 0 {
            a.io.itape = ISCR1;
            a.io.conti();
            a.contog();
            a.io.itape = isave;
            loop {
                if a.rd210(ISCR1).is_none() {
                    // label 180
                    a.report();
                    a.io.outm();
                    a.lstmat = 0;
                    a.lstmf = 0;
                    break 'l10;
                }
                if a.io.h.mat <= 0 {
                    break;
                }
                if a.io.h.mt <= 0 {
                    let (m, f) = (a.io.h.mat, a.io.h.mf);
                    a.io.outs(m, f);
                } else {
                    let o = a.io.otape;
                    a.w210(o);
                    a.lstmat = a.io.h.mat;
                    a.lstmf = a.io.h.mf;
                }
            }
            // label 170
            a.report();
            a.io.outm();
            a.lstmat = 0;
            a.lstmf = 0;
            continue 'l10;
        }
        // label 70
        a.io.itape = ISCR1;
        let mut done = false;
        loop {
            a.nextmt();
            if a.io.h.mat <= 0 {
                a.create10();
                break;
            }
            let (mf, mt) = (a.io.h.mf, a.io.h.mt);
            if mf == 1 && mt == 451 {
                a.contog();
                a.file1();
                continue;
            }
            match mf {
                3 => a.file3(),
                9 => a.copy_to(ISCR3),
                10 if imfile9 != 0 => a.copy_to(ISCR4),
                _ if mf <= 10 => {
                    a.contog();
                    a.io.copys();
                }
                _ => {
                    let hsave = a.io.h;
                    a.create10();
                    a.io.h = hsave;
                    a.contog();
                    loop {
                        a.nextmt();
                        if a.io.h.mat <= 0 {
                            let m = a.lstmat;
                            a.io.outf(m);
                            break;
                        }
                        a.contog();
                        a.io.copys();
                    }
                    done = true;
                }
            }
            if done {
                break;
            }
        }
        // label 170
        a.report();
        a.io.outm();
        a.lstmat = 0;
        a.lstmf = 0;
    }
    // label 190
    a.io.h.mat = -1;
    a.io.h.mf = 0;
    a.io.h.mt = 0;
    a.contog();
    let s = format!("{d}\n{}Tape Totals{}{}\n{d}", " ".repeat(45), i(a.n2toti, 7), i(a.n2toto, 7), d = dash());
    a.io.both(&s);
    if imall10 > 0 {
        a.io.both(format!(
            " WARNING - The purpose of this code is to create MF=10 data\n           if it does not already exist - if it does exist\n           it is copied without checking whether or not it is\n           consistent with MF=9 data. Be WARNED that for the\n           MF=10 created here MF=9 and MF=10 are consistent.\n           But for MF=10 that already existed there is no\n           check that MF=9 and MF=10 are consistent.\n{}",
            dash()
        ));
    }
    a.io.maxie4(1);
    a.io.endit();
}
