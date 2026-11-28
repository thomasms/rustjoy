//! Port of RELABEL (PREPRO 2026-1): renumber the statement labels of a
//! FORTRAN program in increasing order (10, 20, 30, ... in each routine).
//! Input `RELABEL.INP` (two filenames), listing `RELABEL.LST`, scratch file
//! `RELABEL.001`.
//!
//! CARD1 is 1-based; bytes after column 80 stand for the memory after the
//! COMMON block (the Fortran can look a few columns past 80).

use crate::endfio::Endf;
use crate::ffmt::i;
use crate::fio::{self, Rd, L};
use crate::linear::blank;

const MAXLABEL: usize = 100000;
const BLANK: u8 = b' ';
const PARENL: u8 = b'(';
const PARENR: u8 = b')';
const COMMA: u8 = b',';

const ENDITX: [&[u8; 8]; 2] = [b"END     ", b"end     "];
const FMTIT: [&[u8; 6]; 2] = [b"FORMAT", b"format"];
const NLENG: [usize; 9] = [0, 5, 4, 3, 3, 5, 5, 6, 6];
const SPECAL: [[&[u8; 6]; 8]; 2] = [
    [b"GO TO ", b"GOTO  ", b"IF(   ", b"DO    ", b"READ( ", b"OPEN( ", b"WRITE(", b"CLOSE("],
    [b"go to ", b"goto  ", b"if(   ", b"do    ", b"read( ", b"open( ", b"write(", b"close("],
];

struct Rl {
    io: Endf,
    card1: [u8; 100],
    lb1: i32,
    ltop: i32,
    lb2: i32,
    iway: i32,
    ipath: i32,
    ntab: usize,
    inend: i32,
    ic10: i32,
    ic11: i32,
    ic12a: i32,
    ic12b: i32,
    iblank: i32,
    icfmt: i32,
    oldtab: Vec<[u8; 5]>,
    newtab: Vec<[u8; 5]>,
    imused: Vec<i32>,
    iscr: i32,
    // PASS1 SAVEd locals
    ntape: i32,
    ic12x: i32,
}

impl Rl {
    #[inline]
    fn c(&self, k: i32) -> u8 {
        if (1..100).contains(&k) {
            self.card1[k as usize]
        } else {
            BLANK
        }
    }
    #[inline]
    fn set(&mut self, k: i32, v: u8) {
        if (1..100).contains(&k) {
            self.card1[k as usize] = v;
        }
    }

    fn c80(&self) -> &[u8] {
        &self.card1[1..81]
    }

    /// `READX`.
    fn readx(&mut self, unit: i32, which: u8) -> i32 {
        match self.io.u.read(unit) {
            Rd::Eof => -1,
            Rd::Line(l) => {
                let a = fio::rd_a(&l, 0, 80);
                self.card1[1..81].copy_from_slice(&a);
                match which {
                    10 => self.ic10 += 1,
                    _ => self.ic12a += 1,
                }
                self.lb1 = 6;
                self.skipl1();
                1
            }
        }
    }

    /// `WRITEO` to OTAPE (counts in IC11).
    fn writeo(&mut self) {
        let u = self.io.otape;
        if u <= 0 {
            return;
        }
        let mut n = 80;
        while n >= 1 && self.card1[n] == b' ' {
            n -= 1;
        }
        if n == 0 {
            n = 1;
        }
        let l = L::new().b(&self.card1[1..=n]);
        self.io.u.write(u, l);
        self.ic11 += 1;
    }

    fn writes(&mut self) {
        let u = self.iscr;
        let l = L::new().b(&self.card1[1..81]);
        self.io.u.write(u, l);
    }

    fn lst<S: AsRef<[u8]>>(&mut self, s: S) {
        self.io.lst(s);
    }

    fn msg(&mut self, tail: &str, ic: i32) {
        let s = L::new().s(" ").b(self.c80()).s(" (").s(&i(ic, 5)).s(tail);
        self.lst(&s);
    }

    fn msg2(&mut self, tail: &str, ic: i32) {
        let s = L::new().s(" ").b(self.c80()).s("(").s(&i(ic, 5)).s(tail);
        self.lst(&s);
    }

    /// `SKIPR`.
    fn skipr(&mut self) {
        if self.lb1 < self.ltop {
            let ll = self.lb1 + 1;
            for k in ll..=self.ltop {
                if self.c(k) != BLANK {
                    self.lb1 = k - 1;
                    return;
                }
            }
        }
        self.lb1 = self.ltop;
    }

    /// `SKIPL1`.
    fn skipl1(&mut self) {
        for k in (7..=80).rev() {
            if self.c(k) != BLANK {
                self.ltop = k;
                return;
            }
        }
        self.ltop = 6;
    }

    /// `SAMEKEY(7,CARD1,KEYWORD,NKEY,IMOK)`.
    fn samekey(&self, lb1: i32, key: [&[u8]; 2], nkey: usize) -> bool {
        let mut lb4 = lb1 - 1;
        for l in 0..nkey {
            loop {
                lb4 += 1;
                if lb4 > 72 {
                    return false;
                }
                let c = self.c(lb4);
                if c == key[0][l] || c == key[1][l] {
                    break;
                }
                if c == b' ' {
                    continue;
                }
                return false;
            }
        }
        true
    }

    /// `IMSAME(LB1,CARD1,NTRY,IMOK)`: on a match LB1 moves to its end.
    fn imsame(&mut self, ntry: usize) -> bool {
        let mut lb4 = self.lb1;
        let need = NLENG[ntry];
        for l in 0..need {
            loop {
                lb4 += 1;
                if lb4 > 72 {
                    return false;
                }
                let c = self.c(lb4);
                if c == SPECAL[0][ntry - 1][l] || c == SPECAL[1][ntry - 1][l] {
                    break;
                }
                if c == b' ' {
                    continue;
                }
                return false;
            }
        }
        self.lb1 = lb4;
        true
    }

    /// `CARDIN`.
    fn cardin(&mut self) {
        if self.inend > 0 {
            return;
        }
        loop {
            // label 20
            let t = self.io.itape;
            if self.readx(t, 10) < 0 {
                self.inend = 1;
                return;
            }
            if self.ic10 == 1 {
                let s = L::new().s("\n ").b(self.c80()).s("\n ").s(&"-".repeat(79));
                self.lst(&s);
            }
            if self.card1[1] == b'C' || self.card1[1] == b'c' {
                self.writes();
                continue;
            }
            let lab_blank = (1..=5).all(|m| self.card1[m] == BLANK);
            if !lab_blank {
                if self.card1[6] != BLANK {
                    let ic = self.ic10;
                    self.msg(" ) Line Labelled into Column 6", ic);
                    self.writes();
                    continue;
                }
                for _ in 0..4 {
                    if self.card1[5] != BLANK {
                        break;
                    }
                    let mut jj = 4;
                    for _ in 0..4 {
                        self.card1[jj + 1] = self.card1[jj];
                        jj -= 1;
                    }
                    self.card1[1] = BLANK;
                }
            }
            // label 90
            self.skipl1();
            if self.ltop <= 6 {
                let ic = self.ic10;
                self.msg(" ) Blank Line Deleted", ic);
                self.ic10 -= 1;
                self.iblank += 1;
                continue;
            }
            self.lb1 = 6;
            self.skipr();
            if self.lb1 >= self.ltop {
                continue;
            }
            if self.card1[6] == BLANK {
                let mm = 6;
                let ishs = self.lb1 - mm;
                if ishs > 0 {
                    for l in self.lb1..=self.ltop {
                        let v = self.c(l);
                        self.set(l - ishs, v);
                    }
                    let lltop = self.ltop - ishs + 1;
                    for l in lltop..=self.ltop {
                        self.set(l, BLANK);
                    }
                    self.lb1 = mm;
                    self.ltop -= ishs;
                }
            }
            self.writes();
            return;
        }
    }

    /// `INSERT`.
    fn insert(&mut self) {
        let mut iuse = [BLANK; 5];
        let mut lb3 = self.lb2 + self.lb1;
        let mut ll = 5i32;
        for _ in 0..self.lb2.max(0) {
            if ll >= 1 {
                iuse[(ll - 1) as usize] = self.c(lb3);
            }
            lb3 -= 1;
            ll -= 1;
        }
        let mut nuspot = 0usize;
        let mut found = false;
        for n in 1..=self.ntab {
            if self.oldtab[n] == iuse {
                nuspot = n;
                found = true;
                break;
            }
        }
        let no_equiv = |s: &mut Rl, iuse: &[u8; 5]| {
            let ic = s.ic12a;
            let l = L::new().s(" ").b(s.c80()).s(" (").s(&i(ic, 5)).s(" ) No Label Equivalent to ").b(iuse);
            s.lst(&l);
            s.lb1 += s.lb2;
        };
        if !found {
            no_equiv(self, &iuse);
            return;
        }
        if self.ipath != 2 {
            self.imused[nuspot] += 1;
            self.lb1 += self.lb2;
            return;
        }
        if self.imused[nuspot] <= 0 {
            no_equiv(self, &iuse);
            return;
        }
        if self.newtab[nuspot] == iuse {
            self.lb1 += self.lb2;
            return;
        }
        // label 90
        let mut nustar = 5i32;
        for k in 1..=5 {
            if self.newtab[nuspot][k - 1] != BLANK {
                nustar = k as i32;
                break;
            }
        }
        let nusize = 6 - nustar;
        let mut l3 = self.lb1 + self.lb2 + 1;
        let nmove = nusize - self.lb2;
        if nmove > 0 {
            if l3 <= self.ltop {
                let newtop = self.ltop + nmove;
                let mut do160 = true;
                if newtop > 72 {
                    let mut card2 = [BLANK; 100];
                    for k in 1..=self.ltop {
                        card2[k as usize] = self.c(k);
                    }
                    for k in self.lb1..=80 {
                        self.set(k, b' ');
                    }
                    let ic = self.ic12a;
                    self.msg(" ) Split Line to Avoid Overflow Col. 72", ic);
                    self.writeo();
                    for k in 1..=self.lb1 {
                        self.set(k, b' ');
                    }
                    self.card1[6] = b'1';
                    let mut ktop = 7;
                    for k in self.lb1..=self.ltop {
                        ktop += 1;
                        self.set(ktop, card2[k as usize]);
                    }
                    let ic = self.ic12a;
                    self.msg(" ) Overflowed into Columns 1-72", ic);
                    self.lb1 = 8;
                    self.ltop = ktop;
                    l3 = self.lb1 + self.lb2 + 1;
                    if l3 > self.ltop {
                        do160 = false;
                    }
                }
                if do160 {
                    let mut mm = self.ltop;
                    while mm >= l3 {
                        let v = self.c(mm);
                        self.set(mm + nmove, v);
                        mm -= 1;
                    }
                }
            }
        } else if nmove < 0 {
            if l3 <= self.ltop {
                for mm in l3..=self.ltop {
                    let v = self.c(mm);
                    self.set(mm + nmove, v);
                }
            }
            let ktop = self.ltop + nmove + 1;
            if ktop <= self.ltop {
                for k in ktop..=self.ltop {
                    self.set(k, BLANK);
                }
            }
        }
        // label 220
        self.ltop += nmove;
        for k in nustar..=5 {
            self.lb1 += 1;
            let v = self.newtab[nuspot][(k - 1) as usize];
            self.set(self.lb1, v);
        }
        if self.ltop > 72 {
            let ic = self.ic12a;
            self.msg(" ) Line Overflowed Columns 1-72", ic);
        }
    }

    /// `GOTOIT`.
    fn gotoit(&mut self) {
        loop {
            self.skipr();
            if self.lb1 < self.ltop {
                if self.c(self.lb1 + 1) == PARENL {
                    self.lb1 += 1;
                    self.repeatx();
                    return;
                }
                for lb2 in 1..=6 {
                    if self.c(lb2 + self.lb1) == BLANK {
                        self.lb2 = lb2 - 1;
                        self.insert();
                        self.writeo();
                        return;
                    }
                }
                self.lb2 = 7;
            }
            // label 40
            self.writeo();
            let card2: Vec<u8> = self.c80().to_vec();
            let iend = self.readx(self.iscr, 12);
            if !(iend < 0 || self.card1[6] == b' ' || self.card1[1] == b'C' || self.card1[1] == b'c') {
                continue;
            }
            let ic = self.ic12a;
            let s = L::new().s(" ").b(&card2).s(" (").s(&i(ic - 1, 5)).s(" ) Cannot Decode GO TO Label");
            self.lst(&s);
            self.msg(" ) Cannot Decode GO TO Label", ic);
            self.writeo();
            return;
        }
    }

    /// `IFIT`.
    fn ifit(&mut self) {
        let mut level = 1;
        let mut myquote = 0;
        'l10: loop {
            let ll = self.lb1 + 1;
            let mut k = ll;
            let mut closed = false;
            while k <= self.ltop {
                self.lb1 = k;
                if self.c(k) == b'\'' {
                    myquote = 1 - myquote;
                }
                if myquote != 1 {
                    if self.c(k) == PARENR {
                        level -= 1;
                        if level <= 0 {
                            closed = true;
                            break;
                        }
                    } else if self.c(k) == PARENL {
                        level += 1;
                    }
                }
                k += 1;
            }
            if !closed {
                self.lb1 = ll.max(self.ltop + 1);
                self.writeo();
                if self.readx(self.iscr, 12) < 0 {
                    self.iway = -1;
                    return;
                }
                // label 60
                if self.card1[6] == BLANK || self.card1[1] == b'C' || self.card1[1] == b'c' {
                    let ic = self.ic12a;
                    self.msg(" ) Expect IF Continuation", ic);
                    self.iway = 1;
                    return;
                }
                self.skipl1();
                if self.ltop <= 6 {
                    let ic = self.ic12a;
                    self.msg(" ) Expect IF Continuation", ic);
                    self.iway = 1;
                    return;
                }
                self.lb1 = 6;
                self.skipr();
                if self.lb1 >= self.ltop {
                    self.iway = -1;
                    return;
                }
                continue 'l10;
            }
            break;
        }
        // label 90
        self.skipr();
        if self.lb1 >= self.ltop {
            loop {
                // label 100
                self.writeo();
                if self.readx(self.iscr, 12) < 0 {
                    self.iway = -1;
                    return;
                }
                if self.card1[6] == BLANK || self.card1[1] == b'C' || self.card1[1] == b'c' {
                    let ic = self.ic12a;
                    self.msg(" ) Expect IF Continuation", ic);
                    self.iway = 1;
                    return;
                }
                self.skipl1();
                if self.ltop <= 6 {
                    continue;
                }
                self.lb1 = 6;
                self.skipr();
                break;
            }
        }
        // label 110
        if self.imsame(1) {
            self.gotoit();
            return;
        }
        if self.imsame(2) {
            if self.ltop >= 72 {
                let ic = self.ic12a;
                let s = L::new().s("  WARNING - Cannot Convert GOTO to GO TO\n ").b(self.c80()).s(" (").s(&i(ic, 5)).s(" )");
                self.lst(&s);
            } else {
                let mut ll = self.ltop;
                let lb1p3 = self.lb1 - 1;
                let ic = self.ic12a;
                let s = L::new().s("  WARNING - Converted GOTO to GO TO\n ").b(self.c80()).s(" (").s(&i(ic, 5)).s(" )");
                self.lst(&s);
                for _ in lb1p3..=self.ltop {
                    let v = self.c(ll);
                    self.set(ll + 1, v);
                    ll -= 1;
                }
                self.set(ll + 1, BLANK);
                let s = L::new().s(" ").b(self.c80());
                self.lst(&s);
                self.lb1 += 1;
            }
            self.gotoit();
            return;
        }
        if self.imsame(5) || self.imsame(7) {
            self.readit(1);
            return;
        }
        if self.imsame(6) || self.imsame(8) {
            self.readit(2);
            return;
        }
        let mut alfok = 1;
        for l in 1..=6 {
            let c = self.c(self.lb1 + l);
            if c == COMMA {
                self.repeatx();
                return;
            }
            if c.is_ascii_digit() {
                if l == 1 {
                    alfok = 0;
                }
                continue;
            }
            if alfok <= 0 {
                break;
            }
            if c.is_ascii_uppercase() || c.is_ascii_lowercase() {
                continue;
            }
            break;
        }
        self.writeo();
    }

    /// `REPEATX`.
    fn repeatx(&mut self) {
        // IEND is set to 0 only here and at label 120; READX also sets it
        let mut iend = 0;
        #[derive(PartialEq, Clone, Copy)]
        enum G {
            L10,
            L90,
            L100,
        }
        let mut g = G::L10;
        loop {
            match g {
                G::L10 => {
                    let mut found = None;
                    for l in self.lb1 + 1..=self.ltop {
                        let c = self.c(l);
                        if c != BLANK && c != COMMA && c != PARENR {
                            found = Some(l);
                            break;
                        }
                    }
                    let Some(l) = found else {
                        g = G::L100;
                        continue;
                    };
                    self.lb1 = l - 1;
                    let mut stop = 0; // 1 blank, 2 paren, 3 comma
                    let mut lb2 = 1;
                    while lb2 <= 7 {
                        let c = self.c(lb2 + self.lb1);
                        if c == BLANK {
                            stop = 1;
                            break;
                        }
                        if c == COMMA {
                            stop = 3;
                            break;
                        }
                        if c == PARENR {
                            stop = 2;
                            break;
                        }
                        lb2 += 1;
                    }
                    if stop == 0 {
                        let ic = self.ic12a;
                        self.msg(" ) Cannot Decode Repeated Labels", ic);
                        self.writeo();
                        return;
                    }
                    if stop == 1 {
                        iend = 1;
                    }
                    if stop == 2 {
                        iend = 2;
                    }
                    self.lb2 = lb2 - 1;
                    if self.lb2 > 0 {
                        self.insert();
                    }
                    if iend > 0 {
                        if iend == 2 {
                            self.writeo();
                            return;
                        }
                        if self.lb1 >= self.ltop {
                            g = G::L100;
                            continue;
                        }
                    } else {
                        self.lb1 += 1;
                    }
                    g = G::L90;
                }
                G::L90 => {
                    self.skipr();
                    g = if self.lb1 < self.ltop { G::L10 } else { G::L100 };
                }
                G::L100 => {
                    self.writeo();
                    iend = self.readx(self.iscr, 12);
                    if iend < 0 {
                        self.iway = -1;
                        return;
                    }
                    if self.card1[6] != BLANK && self.card1[1] != b'C' && self.card1[1] != b'c' {
                        self.skipl1();
                        if self.ltop <= 6 {
                            continue;
                        }
                        self.lb1 = 6;
                        iend = 0;
                        g = G::L90;
                        continue;
                    }
                    self.iway = 1;
                    return;
                }
            }
        }
    }

    /// `DOIT`.
    fn doit(&mut self) {
        self.skipr();
        if self.lb1 < self.ltop {
            for lb2 in 1..=6 {
                if self.c(lb2 + self.lb1) == BLANK {
                    self.lb2 = lb2 - 1;
                    self.insert();
                    self.writeo();
                    return;
                }
            }
            self.lb2 = 7;
        }
        let ic = self.ic12a;
        self.msg(" ) Cannot Decode DO Label", ic);
        self.writeo();
    }

    /// `READIT(MYWAY)`.
    fn readit(&mut self, myway: i32) {
        let fmteq: [&[u8; 4]; 2] = [b"FMT=", b"fmt="];
        let endeq: [&[u8; 4]; 2] = [b"END=", b"end="];
        let erreq: [&[u8; 4]; 2] = [b"ERR=", b"err="];
        let mut myparen = 0;
        let mut level = 0;
        let mut ipass = if myway == 2 { 3 } else { 1 };
        let ltopm3 = self.ltop.min(72);
        let mut llj = 0;
        if self.lb1 > ltopm3 {
            let ic = self.ic12a;
            self.msg(" ) Cannot Find Balancing Parenthesis", ic);
            self.writeo();
            return;
        }
        'l10: loop {
            let mut ll = self.lb1 - 1;
            // label 20
            let found110 = 'l20: loop {
                ll += 1;
                if ll > 99 {
                    // the Fortran would scan memory past the line
                    break 'l20 false;
                }
                if self.c(ll) == PARENL {
                    level += 1;
                }
                if self.c(ll) == PARENR {
                    level -= 1;
                    if level <= 0 {
                        self.writeo();
                        return;
                    }
                }
                if self.c(ll) == BLANK {
                    continue;
                }
                if ipass == 1 {
                    if self.c(ll) == COMMA {
                        ipass = 2;
                    }
                    continue;
                }
                let m = |s: &Rl, key: &[&[u8; 4]; 2], llj: &mut i32| -> bool {
                    for j in 1..=4 {
                        *llj = ll + j - 1;
                        let c = s.c(*llj);
                        if c != key[0][(j - 1) as usize] && c != key[1][(j - 1) as usize] {
                            return false;
                        }
                    }
                    true
                };
                if m(self, &fmteq, &mut llj) || m(self, &endeq, &mut llj) || m(self, &erreq, &mut llj) {
                    break 'l20 true;
                }
                // label 90
                if ipass == 2 {
                    llj = ll - 1;
                    break 'l20 true;
                }
                if ll < ltopm3 {
                    continue;
                }
                break 'l20 false;
            };
            if found110 {
                self.lb1 = llj;
                let mut hit = 0;
                let mut lb2 = 1;
                while lb2 <= 7 {
                    let c = self.c(self.lb1 + lb2);
                    if c == PARENR {
                        hit = 1;
                        break;
                    }
                    if c == COMMA {
                        hit = 2;
                        break;
                    }
                    lb2 += 1;
                }
                if hit > 0 {
                    if hit == 1 {
                        myparen = 1;
                    }
                    self.lb2 = lb2 - 1;
                    self.insert();
                    ipass = 3;
                    if myparen <= 0 {
                        continue 'l10;
                    }
                    self.writeo();
                    return;
                }
                self.lb2 = 8;
            }
            // label 130
            let ic = self.ic12a;
            self.msg(" ) Cannot Find Balancing Parenthesis", ic);
            self.writeo();
            return;
        }
    }

    /// `PASS1`.
    fn pass1(&mut self) {
        if self.ipath != 2 {
            self.ntape = self.io.otape;
            self.io.otape = 0;
            self.ic12x = self.ic12a;
            self.ic12a = self.ic12b;
            for j in 1..=self.ntab {
                self.imused[j] = 0;
            }
        } else if self.ntab > 0 {
            let mut j = 0;
            for jj in 1..=self.ntab {
                if self.imused[jj] <= 0 {
                    self.newtab[jj] = [BLANK; 5];
                    continue;
                }
                j += 1;
                self.newtab[jj] = [BLANK, BLANK, BLANK, BLANK, b'0'];
                let mut j2 = j;
                let mut k = 4;
                loop {
                    let i1 = j2 / 10;
                    let i2 = j2 - 10 * i1;
                    self.newtab[jj][k - 1] = b'0' + i2 as u8;
                    if i1 <= 0 {
                        break;
                    }
                    j2 = i1;
                    if k == 1 {
                        break;
                    }
                    k -= 1;
                }
            }
        }
        // label 80
        let u = self.iscr;
        self.io.u.rewind(u);
        if !(self.ntab == 0 || self.ipath == 1) && (1..=self.ntab).all(|n| self.oldtab[n] == self.newtab[n]) {
            self.ntab = 0;
        }
        #[derive(PartialEq, Clone, Copy)]
        enum G {
            L110,
            L120,
            L130,
            L350,
            L360,
        }
        let mut g = G::L110;
        loop {
            match g {
                G::L110 => {
                    let u = self.iscr;
                    g = if self.readx(u, 12) < 0 { G::L120 } else { G::L130 };
                }
                G::L120 => {
                    self.lst(" ****** Program ERROR...End of File on Scratch ******");
                    g = G::L360;
                }
                G::L130 => {
                    if self.card1[1] == b'C' || self.card1[1] == b'c' {
                        g = G::L350;
                        continue;
                    }
                    self.skipl1();
                    if self.ltop <= 6 {
                        g = G::L350;
                        continue;
                    }
                    self.lb1 = 6;
                    self.skipr();
                    if self.lb1 >= self.ltop {
                        g = G::L350;
                        continue;
                    }
                    if self.samekey(7, [ENDITX[0], ENDITX[1]], 8) {
                        self.writeo();
                        if !(self.ntab == 0 || self.ipath == 2) {
                            let mut iuse = 0;
                            for l in 1..=self.ntab {
                                if self.imused[l] > 0 {
                                    continue;
                                }
                                if iuse <= 0 {
                                    self.lst(" Statement Labels NOT Referenced (Will be Set to Blank)");
                                }
                                iuse += 1;
                                let s = L::new().s(" ").b(&self.oldtab[l].clone());
                                self.lst(&s);
                            }
                        }
                        g = G::L360;
                        continue;
                    }
                    // label 150
                    if self.ntab == 0 {
                        g = G::L350;
                        continue;
                    }
                    let mut go230 = false;
                    let mut nfound = 0;
                    if self.card1[5] == BLANK {
                        go230 = true;
                    } else {
                        for n in 1..=self.ntab {
                            if self.card1[1..6] == self.oldtab[n] {
                                nfound = n;
                                break;
                            }
                        }
                        if nfound == 0 {
                            let ic = self.ic12a;
                            self.msg(" ) No Equivalent Statement Label", ic);
                            go230 = true;
                        }
                    }
                    if !go230 {
                        let n = nfound;
                        let mut setnew = true;
                        if !(self.ipath != 2 || self.imused[n] > 0) {
                            if self.samekey(7, [FMTIT[0], FMTIT[1]], 6) {
                                let ic = self.ic12a;
                                self.msg2(" ) Unused Format Deleted", ic);
                                self.icfmt += 1;
                                // label 190
                                loop {
                                    let u = self.iscr;
                                    if self.readx(u, 12) < 0 {
                                        g = G::L120;
                                        break;
                                    }
                                    if self.card1[1] == b'C' || self.card1[1] == b'c' {
                                        g = G::L350;
                                        break;
                                    }
                                    if self.card1[6] == b' ' {
                                        g = G::L130;
                                        break;
                                    }
                                    self.icfmt += 1;
                                }
                                continue;
                            }
                            let ic = self.ic12a;
                            self.msg2(" ) WARNING - Statement Label Set Blank", ic);
                            setnew = true;
                        }
                        if setnew {
                            let nt = self.newtab[n];
                            self.card1[1..6].copy_from_slice(&nt);
                        }
                    }
                    // label 230
                    let mut ntry = 0;
                    for k in 1..=8 {
                        if self.imsame(k) {
                            ntry = k;
                            break;
                        }
                    }
                    if ntry == 0 {
                        g = G::L350;
                        continue;
                    }
                    self.iway = 0;
                    match ntry {
                        2 => {
                            if self.ltop >= 72 {
                                let ic = self.ic12a;
                                let s = L::new().s("  WARNING - Cannot Convert GOTO to GO TO\n ").b(self.c80()).s("(").s(&i(ic, 5)).s(" )");
                                self.lst(&s);
                            } else {
                                let mut ll = self.ltop;
                                let lb1m1 = self.lb1 - 1;
                                let ic = self.ic12a;
                                let s = L::new().s("  WARNING - Converted GOTO to GO TO\n ").b(self.c80()).s(" (").s(&i(ic, 5)).s(" )");
                                self.lst(&s);
                                for _ in lb1m1..=self.ltop {
                                    let v = self.c(ll);
                                    self.set(ll + 1, v);
                                    ll -= 1;
                                }
                                self.set(ll + 1, BLANK);
                                let s = L::new().s(" ").b(self.c80());
                                self.lst(&s);
                                self.lb1 += 1;
                            }
                            self.gotoit();
                        }
                        1 => self.gotoit(),
                        3 => self.ifit(),
                        4 => self.doit(),
                        5 | 7 => self.readit(1),
                        _ => self.readit(2),
                    }
                    g = if self.iway < 0 {
                        G::L120
                    } else if self.iway == 0 {
                        G::L110
                    } else {
                        G::L130
                    };
                }
                G::L350 => {
                    self.writeo();
                    g = G::L110;
                }
                G::L360 => {
                    if self.ipath != 2 {
                        self.io.otape = self.ntape;
                        self.ic12b = self.ic12a;
                        self.ic12a = self.ic12x;
                    }
                    return;
                }
            }
        }
    }
}

pub fn run() -> ! {
    let mut io = Endf::new("RELABEL");
    io.timer();
    io.inp = 2;
    io.outp = 3;
    io.itape = 10;
    io.otape = 11;
    io.u.open(3, "RELABEL.LST");
    let mut istat1 = if io.u.open_old(2, "RELABEL.INP") { 0 } else { 1 };
    let d = format!(" {}", "-".repeat(79));
    io.both(format!(" Incremental Statement Labels (RELABEL 2026-1)\n{d}"));
    // READIN
    let mut namein = None;
    let mut nameout = None;
    if istat1 != 1 {
        if let Rd::Line(l) = io.u.read(2) {
            let mut n = fio::rd_a(&l, 0, 72);
            if blank(&n) {
                n = fio::rd_a(b"RELABEL.IN", 0, 72);
            }
            namein = Some(n);
            if let Rd::Line(l) = io.u.read(2) {
                let mut n = fio::rd_a(&l, 0, 72);
                if blank(&n) {
                    n = fio::rd_a(b"RELABEL.OUT", 0, 72);
                }
                nameout = Some(n);
            }
        }
    }
    let namein = namein.unwrap_or_else(|| fio::rd_a(b"RELABEL.IN", 0, 72));
    let nameout = match nameout {
        Some(n) => n,
        None => {
            istat1 = 1;
            fio::rd_a(b"RELABEL.OUT", 0, 72)
        }
    };
    let _ = istat1;
    let s = L::new().s(" Program Input and Output Filenames\n ").b(&namein).s("\n ").b(&nameout).s("\n").s(&d);
    io.both(&s);
    let ni = String::from_utf8_lossy(&namein).into_owned();
    if !io.u.open_old(10, &ni) {
        let s = L::new().s("\n\n ERROR - Opening original file\n ").b(&namein).s("\n\n");
        io.both(&s);
        io.enderror();
    }
    io.u.open(12, "RELABEL.001");
    let no = String::from_utf8_lossy(&nameout).into_owned();
    io.u.open(11, &no);
    let mut r = Rl {
        io,
        card1: [BLANK; 100],
        lb1: 0,
        ltop: 0,
        lb2: 0,
        iway: 0,
        ipath: 0,
        ntab: 0,
        inend: 0,
        ic10: 0,
        ic11: 0,
        ic12a: 0,
        ic12b: 0,
        iblank: 0,
        icfmt: 0,
        oldtab: vec![[0u8; 5]; MAXLABEL + 1],
        // blank COMMON starts as zero bytes; PASS1 path 1 copies entries
        // that path 2 has not set yet
        newtab: vec![[0u8; 5]; MAXLABEL + 1],
        imused: vec![0; MAXLABEL + 1],
        iscr: 12,
        ntape: 0,
        ic12x: 0,
    };
    r.cardin();
    #[derive(PartialEq, Clone, Copy)]
    enum G {
        L10,
        L20,
        L60,
        L80,
    }
    let mut g = if r.inend <= 0 { G::L20 } else { G::L80 };
    loop {
        match g {
            G::L10 => {
                r.ntab = 0;
                let u = r.iscr;
                r.io.u.rewind(u);
                r.cardin();
                if r.inend > 0 {
                    g = G::L80;
                    continue;
                }
                let s = L::new().s("\n ").b(r.c80()).s("\n ").s(&"-".repeat(79));
                r.io.both(&s);
                g = G::L20;
            }
            G::L20 => {
                r.cardin();
                if r.inend > 0 {
                    g = G::L60;
                    continue;
                }
                if r.card1[5] != BLANK {
                    r.ntab += 1;
                    if r.ntab > MAXLABEL {
                        r.io.both(format!(" Over {} Statement Labels...Execution Terminated", i(MAXLABEL as i32, 6)));
                        r.io.enderror();
                    }
                    let n = r.ntab;
                    let v: [u8; 5] = r.card1[1..6].try_into().unwrap();
                    r.oldtab[n] = v;
                    continue;
                }
                if !r.samekey(7, [ENDITX[0], ENDITX[1]], 8) {
                    continue;
                }
                let u = r.iscr;
                r.io.u.endfile(u);
                r.ipath = 1;
                r.pass1();
                r.ipath = 2;
                r.pass1();
                g = G::L10;
            }
            G::L60 => {
                let u = r.iscr;
                r.io.u.endfile(u);
                r.io.u.rewind(u);
                r.io.both(" ****** End of File Before End of Routine ******");
                loop {
                    let u = r.iscr;
                    if r.readx(u, 12) < 0 {
                        break;
                    }
                    r.writeo();
                }
                g = G::L80;
            }
            G::L80 => {
                let mut lines = vec![format!(" {}\n {} Lines Translated", " ".repeat(70), i(r.ic10, 8))];
                if r.iblank != 0 {
                    lines.push(format!(" {} Blank Lines Deleted", i(r.iblank, 8)));
                }
                if r.icfmt != 0 {
                    lines.push(format!(" {} Unused Format Lines Deleted", i(r.icfmt, 8)));
                }
                if r.ic10 != r.ic11 || r.ic10 != r.ic12a || r.ic10 != r.ic12b {
                    lines.push(format!(
                        "\n\n\n Line Count NOT Conserved\n Read---------------{}\n Written------------{}\n Scratch (Pass 1)---{}\n Scratch (Pass 2)---{}",
                        i(r.ic10, 8),
                        i(r.ic11, 8),
                        i(r.ic12b, 8),
                        i(r.ic12a, 8)
                    ));
                }
                for l in &lines {
                    r.io.lst(l);
                }
                for l in &lines {
                    fio::out(l);
                }
                r.io.endit();
            }
        }
    }
}
