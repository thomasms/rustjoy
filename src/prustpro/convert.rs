//! Port of CONVERT (PREPRO 2026-1): switch the optional statements of a
//! FORTRAN program on or off (lines between `C***` keyword pairs).
//! Input `CONVERT.INP`, listing `CONVERT.LST`.

use crate::endfio::Endf;
use crate::ffmt::i;
use crate::fio::{self, Rd, L};

/// Write to CONVERT.LST and the screen.
fn both<S: AsRef<[u8]>>(io: &mut Endf, s: S) {
    io.u.write(3, s.as_ref());
    fio::out(s);
}
use crate::linear::blank;

const KEYWRD: [&[u8; 10]; 13] = [
    b"DOUBLE    ",
    b"SINGLE    ",
    b"CHARACTER ",
    b"INTEGER   ",
    b"STOP      ",
    b"EXIT      ",
    b"CIVIC     ",
    b"NAMES     ",
    b"REWIND    ",
    b"ENDFILE   ",
    b"SAVE      ",
    b"ACTION    ",
    b"CLOSE     ",
];
const PROGRM: &[u8; 10] = b"PROGRAM   ";
const NOID: &[u8; 10] = b"NOID      ";

/// `SAMKEY`: word2 at word1(nstart); on a match nstart moves past it.
/// word1 is 1-based with 72 columns.
fn samkey(word1: &[u8], nstart: &mut usize, word2: &[u8; 10]) -> bool {
    let mut inow = *nstart as i64 - 1;
    for k in 0..10 {
        inow += 1;
        if word2[k] == b' ' {
            if inow <= 72 && word1[inow as usize] != b' ' {
                return false;
            }
            *nstart = inow as usize;
            return true;
        }
        if inow > 72 {
            return false;
        }
        if word2[k] != word1[inow as usize] {
            return false;
        }
    }
    false
}

pub fn run() -> ! {
    let mut io = Endf::new("CONVERT");
    // CONVERT keeps its units in COMMON/UNITS/, so /ENDFIO/ (used by the
    // timer) stays zero: the times go to unit 0 (standard error).
    io.timer();
    io.u.open(3, "CONVERT.LST");
    let mut istat1 = if io.u.open_old(2, "CONVERT.INP") { 0 } else { 1 };
    let eq = format!("  {}", "=".repeat(78));
    both(&mut io, format!(" Convert FORTRAN Programs (CONVERT 2026-1)\n{eq}"));
    let mut ierr = 0;
    if istat1 == 1 {
        both(&mut io, "\n\n ERROR - opening CONVERT.INP input parameter file\n\n");
        io.enderror();
    }
    let mut iline = 0;
    let mut iend = 0;
    let mut imon = 0usize;
    let mut turnon = [0i32; 14];
    let mut modify = [[0i32; 2]; 14];
    let mut iprog = 0;
    let mut inoid = 0;
    let l = io.u.read_req(2);
    let mut keyin = vec![b' '; 73];
    keyin[1..].copy_from_slice(&fio::rd_a(&l, 0, 72));
    both(&mut io, format!("  Input Keyword    Comments\n{eq}"));
    let mut luse = 1usize;
    let mut lfield = 0;
    let fld = |n: i32, w: &[u8]| L::new().s("  (").s(&i(n, 2)).s(") ").b(w);
    // DO 100 LNOW=LOWL,72 with GO TO 30 restarting it at LUSE
    'l30: loop {
        let lowl = luse;
        let mut lnow = lowl;
        while lnow <= 72 {
            if keyin[lnow] == b' ' {
                lnow += 1;
                continue;
            }
            lfield += 1;
            luse = lnow;
            let mut hit = false;
            for n in 1..=13 {
                if samkey(&keyin, &mut luse, KEYWRD[n - 1]) {
                    turnon[n] = 1;
                    both(&mut io, fld(lfield, KEYWRD[n - 1]));
                    hit = true;
                    break;
                }
            }
            if !hit {
                if samkey(&keyin, &mut luse, PROGRM) {
                    both(&mut io, fld(lfield, PROGRM));
                    iprog = 1;
                    hit = true;
                } else if samkey(&keyin, &mut luse, NOID) {
                    both(&mut io, fld(lfield, NOID));
                    inoid = 1;
                    hit = true;
                }
            }
            if hit {
                if luse < 72 {
                    continue 'l30;
                }
                lnow += 1;
                continue;
            }
            // label 60
            let mut dumfld = [b' '; 10];
            let mut lfill = 0;
            luse = lnow;
            let mut stopped = false;
            while luse <= 72 {
                if keyin[luse] == b' ' {
                    stopped = true;
                    break;
                }
                lfill += 1;
                if lfill <= 10 {
                    dumfld[lfill - 1] = keyin[luse];
                }
                luse += 1;
            }
            if !stopped {
                luse = 72;
            }
            both(&mut io, fld(lfield, &dumfld).s("  Cannot Recognize...Input ERROR"));
            ierr = 1;
            luse += 1;
            if luse < 72 {
                continue 'l30;
            }
            lnow += 1;
        }
        break;
    }
    if ierr != 0 {
        let mut s = format!("{eq}\n  The Following are the ONLY Legal Input Keywords\n{eq}");
        let mut k = L::new();
        for w in KEYWRD.iter() {
            k.push("  ");
            k.pushb(&w[..]);
            k.push("\n");
        }
        k.push("  ");
        k.pushb(PROGRM);
        k.push("\n  ");
        k.pushb(NOID);
        k.push(&format!("\n{eq}\n  Execution Terminated"));
        s.push('\n');
        let all = L::new().s(&s).b(&k.0);
        both(&mut io, &all);
        io.enderror();
    }
    let mut ipass = 0;
    let hdr = format!("{eq}\n  Required Modifications\n{eq}\n  Keyword    Comments\n{eq}");
    for i1 in [1usize, 3, 5] {
        let i2 = i1 + 1;
        if turnon[i1] == 0 && turnon[i2] == 0 {
            turnon[i1] = 1;
            let s = L::new().s("  ").b(KEYWRD[i1 - 1]).s(" Activated.....Standard Option");
            if ipass == 0 {
                io.u.write(3, &hdr);
            }
            io.u.write(3, &s);
            if ipass == 0 {
                fio::out(&hdr);
            }
            fio::out(&s);
            ipass = 1;
        }
        if turnon[i1] == 1 && turnon[i2] == 1 {
            turnon[i2] = 0;
            let s = L::new().s("  ").b(KEYWRD[i2 - 1]).s(" Deactivated...Comflicts with ").b(KEYWRD[i1 - 1]);
            if ipass == 0 {
                io.u.write(3, &hdr);
            }
            io.u.write(3, &s);
            if ipass == 0 {
                fio::out(&hdr);
            }
            fio::out(&s);
            ipass = 1;
        }
    }
    if turnon[2] == 1 {
        both(&mut io, format!("{eq}\n  WARNING...You have requested the output program be setup\n            for SINGLE PRECISION. The output program should\n            only be used on long word length computers (e.g.,\n            CDC computers). Errors may occur during execution\n            if the output program is used on a short word\n            length computer (e.g., IBM computers)."));
    }
    // file names
    let mut namein = None;
    let mut nameout = None;
    if let Rd::Line(l) = io.u.read(2) {
        let mut n = fio::rd_a(&l, 0, 72);
        if blank(&n) {
            n = fio::rd_a(b"CONVERT.IN", 0, 72);
        }
        namein = Some(n);
        if let Rd::Line(l) = io.u.read(2) {
            let mut n = fio::rd_a(&l, 0, 72);
            if blank(&n) {
                n = fio::rd_a(b"CONVERT.OUT", 0, 72);
            }
            nameout = Some(n);
        }
    }
    let namein = namein.unwrap_or_else(|| fio::rd_a(b"CONVERT.IN", 0, 72));
    let nameout = match nameout {
        Some(n) => n,
        None => {
            istat1 = 1;
            fio::rd_a(b"CONVERT.OUT", 0, 72)
        }
    };
    let _ = istat1;
    let d = format!("  {}", "-".repeat(78));
    let s = L::new().s(&d).s("\n  Input and Output Data Filenames\n  ").b(&namein).s("\n  ").b(&nameout);
    both(&mut io, &s);
    let no = String::from_utf8_lossy(&nameout).into_owned();
    io.u.open(11, &no);
    let ni = String::from_utf8_lossy(&namein).into_owned();
    if !io.u.open_old(10, &ni) {
        let s = L::new().s("\n\n ERROR - opening old file\n ").b(&namein).s("\n\n");
        both(&mut io, &s);
        io.enderror();
    }
    let mut line1 = vec![b' '; 81];
    let readin = |io: &mut Endf, line1: &mut Vec<u8>, iline: &mut i32, iend: &mut i32| match io.u.read(10) {
        Rd::Line(l) => {
            line1[1..].copy_from_slice(&fio::rd_a(&l, 0, 80));
            *iline += 1;
        }
        Rd::Eof => *iend = 1,
    };
    readin(&mut io, &mut line1, &mut iline, &mut iend);
    let s = L::new()
        .s(&eq)
        .s("\n  First Line of FORTRAN Program\n")
        .s(&eq)
        .s("\n  ")
        .b(&line1[1..])
        .s("\n")
        .s(&eq)
        .s("\n   Index  Line and Comments");
    both(&mut io, &s);
    let comin = b"C**";
    let wline = |io: &mut Endf, line1: &[u8]| {
        if inoid == 0 {
            io.u.write(11, L::new().b(&line1[1..]));
        } else {
            let mut kmax = 72;
            while kmax >= 1 && line1[kmax] == b' ' {
                kmax -= 1;
            }
            if kmax == 0 {
                kmax = 1;
            }
            io.u.write(11, L::new().b(&line1[1..=kmax]));
        }
    };
    // PROGRAM line(s)
    let mut is_prog = false;
    for k in 1..=72 {
        if line1[k] == PROGRM[0] {
            let mut ii = k;
            is_prog = samkey(&line1, &mut ii, PROGRM);
            break;
        }
    }
    if is_prog {
        loop {
            line1[1] = if iprog == 0 { comin[0] } else { b' ' };
            wline(&mut io, &line1);
            readin(&mut io, &mut line1, &mut iline, &mut iend);
            if line1[6] != b' ' {
                continue;
            }
            break;
        }
    }
    // label 310 for the current line, then 300 for each next one
    let mut first = true;
    loop {
        if !first {
            readin(&mut io, &mut line1, &mut iline, &mut iend);
            if iend > 0 {
                break;
            }
        }
        // the first line enters at label 310, after the C--- test of 300
        line_300(&mut line1, !first, iline, &mut io, &mut imon, &mut turnon, &mut modify, comin);
        first = false;
        wline(&mut io, &line1);
    }
    if imon != 0 {
        let s = format!(
            "{eq}\n          ERROR in FORTRAN File...Keywords are Not Paired\n          Expected {}\n          Found End of File",
            String::from_utf8_lossy(KEYWRD[imon - 1])
        );
        both(&mut io, &s);
    }
    both(&mut io, format!("{eq}\n  Summary of Changes\n{eq}\n  Keyword     Statements Turned\n                     On     Off  Status\n{eq}"));
    for n in 1..=13 {
        let flag = if turnon[n] == 1 { " ON " } else { "OFF " };
        let s = L::new().s("  ").b(KEYWRD[n - 1]).s("   ").s(&i(modify[n][0], 8)).s(&i(modify[n][1], 8)).s("     ").s(flag);
        both(&mut io, &s);
    }
    let flag = if iprog == 1 { " ON " } else { "OFF " };
    both(&mut io, L::new().s("  ").b(PROGRM).s(&" ".repeat(24)).s(flag));
    let flag = if inoid == 1 { " ON " } else { "OFF " };
    both(&mut io, L::new().s("  ").b(NOID).s(&" ".repeat(24)).s(flag));
    both(&mut io, format!("{eq}\n  End of Run{} Lines Processed", i(iline, 7)));
    io.endit();
}

/// Labels 300-380 for one line (after the line is read): the C--- lines
/// are copied as they are, the C*** keyword lines switch a keyword on or
/// off, and the lines between them are switched.
#[allow(clippy::too_many_arguments)]
fn line_300(
    line1: &mut [u8],
    dash_test: bool,
    iline: i32,
    io: &mut Endf,
    imon: &mut usize,
    turnon: &mut [i32; 14],
    modify: &mut [[i32; 2]; 14],
    comin: &[u8; 3],
) {
    if dash_test && line1[1] == comin[0] && line1[2] == b'-' && line1[3] == b'-' {
        return;
    }
    if line1[1] == comin[0] && line1[2] == comin[1] && line1[3] == comin[2] {
        let mut start = 0;
        for k in 4..=72 {
            if line1[k] != comin[1] && line1[k] != b' ' {
                start = k;
                break;
            }
        }
        if start > 0 {
            for n in 1..=13 {
                let mut ii = start;
                if !samkey(line1, &mut ii, KEYWRD[n - 1]) {
                    continue;
                }
                if *imon > 0 {
                    if *imon == n {
                        *imon = 0;
                        return;
                    }
                    let eq = format!("  {}", "=".repeat(78));
                    let s = L::new()
                        .s(&eq)
                        .s("\n")
                        .s(&i(iline, 8))
                        .s("  ")
                        .b(&line1[1..])
                        .s("\n          ERROR in FORTRAN File...Keywords are Not Paired\n          Expected ")
                        .b(KEYWRD[*imon - 1])
                        .s("\n           Found    ")
                        .b(KEYWRD[n - 1]);
                    both(io, &s);
                }
                *imon = n;
                return;
            }
        }
    }
    // label 370
    if *imon == 0 {
        return;
    }
    if turnon[*imon] == 0 {
        line1[1] = comin[0];
        modify[*imon][1] += 1;
    } else {
        line1[1] = b' ';
        modify[*imon][0] += 1;
    }
}
