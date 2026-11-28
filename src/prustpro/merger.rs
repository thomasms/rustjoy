//! Port of MERGER (PREPRO 2026-1): merge ENDF/B data from several files
//! into one file in MAT/MF/MT order. Input `MERGER.INP`, listing
//! `MERGER.LST`.

use crate::endfio::{fint, in9, Endf};
use crate::ffmt::i;
use crate::fio::{self, Rd, L};
use crate::linear::{blank, rdi};

const MAXTAPE: usize = 1000;

fn dash() -> String {
    format!(" {}", "-".repeat(72))
}
fn star() -> String {
    format!(" {}", "*".repeat(72))
}

struct Merger {
    io: Endf,
    ntape: Vec<i32>,
    igethow: i32,
    mywish: usize,
    // IWISH(3,2,MAXTAPE) as [k][j] per request
    iwish: Vec<[[i32; 3]; 2]>,
    matgot: Vec<i32>,
    mtgot: Vec<i32>,
    icardgot: Vec<i32>,
    matnow: Vec<i32>,
    nowish: Vec<i32>,
    fileout: Vec<u8>,
    filein: Vec<Vec<u8>>,
    zatab: Vec<f64>,
    izatab: Vec<i32>,
    nunit: usize,
    istat1: i32,
    // PHILIP
    philip_name: [u8; 9],
    ix1: usize,
    ix2: usize,
}

/// A line read with format (2A4,A3,13A4,A3,I4,I2,I3[,I5]): the 66 text
/// columns and MAT, MF, MT.
fn rd40(l: &[u8], unit: i32) -> ([u8; 66], [i32; 3]) {
    let mut c = [b' '; 66];
    c.copy_from_slice(&fio::rd_a(l, 0, 66));
    (c, [fio::rd_i(l, 66, 4, unit), fio::rd_i(l, 70, 2, unit), fio::rd_i(l, 72, 3, unit)])
}

fn w40(c: &[u8], n: [i32; 3], seq: i32) -> L {
    L::new().b(c).s(&i(n[0], 4)).s(&i(n[1], 2)).s(&i(n[2], 3)).s(&i(seq, 5))
}

fn w460(n: [i32; 3], seq: i32) -> String {
    format!("{}{}{}{}{}", " ".repeat(66), i(n[0], 4), i(n[1], 2), i(n[2], 3), i(seq, 5))
}

impl Merger {
    fn rd_inp(&mut self) -> Option<Vec<u8>> {
        match self.io.u.read(2) {
            Rd::Line(l) => Some(l),
            Rd::Eof => None,
        }
    }

    fn philip(&mut self, f: &mut Vec<u8>) {
        if f.iter().take(60).any(|&c| c != b' ') {
            return;
        }
        let fillx = b"1234567890";
        self.ix1 += 1;
        self.philip_name[8] = fillx[self.ix1 - 1];
        if self.ix1 >= 10 {
            self.ix1 = 0;
            self.ix2 += 1;
            self.philip_name[7] = fillx[self.ix2 - 1];
        }
        f[..9].copy_from_slice(&self.philip_name);
    }

    fn ask4it(&mut self) {
        let upper = [9999, 99, 999];
        let mut mw = 1usize;
        let mut ended = false;
        while mw <= 100 {
            self.mywish = mw;
            self.matgot[mw] = 0;
            self.mtgot[mw] = 0;
            self.icardgot[mw] = 0;
            self.matnow[mw] = 0;
            let Some(l) = self.rd_inp() else {
                ended = true;
                break;
            };
            let v = [rdi(&l, 0, 6), rdi(&l, 6, 2), rdi(&l, 8, 3), rdi(&l, 11, 6), rdi(&l, 17, 2), rdi(&l, 19, 3), rdi(&l, 22, 11)];
            // items before a bad one are stored
            let mut bad = false;
            for (k, x) in v.iter().enumerate() {
                match x {
                    Some(x) => {
                        if k < 6 {
                            self.iwish[mw][k / 3][k % 3] = *x;
                        } else {
                            self.nowish[mw] = *x;
                        }
                    }
                    None => {
                        bad = true;
                        break;
                    }
                }
            }
            if bad {
                ended = true;
                break;
            }
            if !self.iwish[mw].iter().flatten().any(|&x| x > 0) {
                ended = true;
                break;
            }
            if mw == 1 {
                let s = if self.igethow == 0 {
                    format!("{d}\n Retrieval Criteria\n{d}\n Request *   Lower Limit *   Upper Limit  \n  Number *    MAT MF  MT *    MAT MF  MT  \n{d}", d = dash())
                } else {
                    format!("{d}\n Retrieval Criteria\n{d}\n Request *   Lower Limit *   Upper Limit  \n  Number *     ZA MF  MT *     ZA MF  MT  \n{d}", d = dash())
                };
                self.io.both(&s);
            }
            if self.iwish[mw][1][0] < self.iwish[mw][0][0] {
                self.iwish[mw][1][0] = self.iwish[mw][0][0];
            }
            for k in 0..3 {
                if self.iwish[mw][1][k] <= 0 {
                    self.iwish[mw][1][k] = upper[k];
                }
            }
            if self.nowish[mw] <= 0 {
                self.nowish[mw] = mw as i32;
            }
            let w = self.iwish[mw];
            let s = format!(
                "{} *{}{}{} *{}{}{}",
                i(self.nowish[mw], 8),
                i(w[0][0], 7),
                i(w[0][1], 3),
                i(w[0][2], 4),
                i(w[1][0], 7),
                i(w[1][1], 3),
                i(w[1][2], 4)
            );
            self.io.both(&s);
            mw += 1;
        }
        if ended {
            self.mywish -= 1;
            if self.mywish == 0 {
                self.io.both(format!("{}\n Retrieve All Data", dash()));
            }
            return;
        }
        // after 100 requests: one more line must be blank
        if let Some(l) = self.rd_inp() {
            let mut vals = Vec::new();
            for (s, w) in [(0, 6), (6, 2), (8, 3), (11, 6), (17, 2), (19, 3)] {
                match rdi(&l, s, w) {
                    Some(x) => vals.push(x),
                    None => break,
                }
            }
            if vals.len() == 6 && vals.iter().any(|&x| x > 0) {
                self.io.both(format!("{}\n Over 100 Range Requests for Retrieval\n Execution Terminated", dash()));
                self.io.enderror();
            }
        }
        self.mywish = 100;
    }

    fn filio2(&mut self) {
        let fo = String::from_utf8_lossy(&self.fileout).into_owned();
        self.io.u.open(10, &fo);
        let mut istat2 = 0;
        for n in 1..=self.nunit {
            let u = self.ntape[n];
            let f = String::from_utf8_lossy(&self.filein[n]).into_owned();
            if !self.io.u.open_old(u, &f) {
                istat2 = 1;
                let s = L::new().s(" ERROR - opening ").b(&self.filein[n]);
                self.io.both(&s);
            }
        }
        if istat2 == 1 {
            self.io.enderror();
        }
    }
}

pub fn run() -> ! {
    let mut io = Endf::new("MERGER");
    io.timer();
    io.inp = 2;
    io.outp = 3;
    io.itape = 11;
    io.otape = 10;
    io.u.open(3, "MERGER.LST");
    let istat1 = if io.u.open_old(2, "MERGER.INP") { 0 } else { 1 };
    let mut m = Merger {
        io,
        ntape: vec![0; MAXTAPE + 1],
        igethow: 0,
        mywish: 0,
        iwish: vec![[[0; 3]; 2]; MAXTAPE + 1],
        matgot: vec![0; MAXTAPE + 1],
        mtgot: vec![0; MAXTAPE + 1],
        icardgot: vec![0; MAXTAPE + 1],
        matnow: vec![0; MAXTAPE + 1],
        nowish: vec![0; MAXTAPE + 1],
        fileout: vec![b' '; 72],
        filein: vec![vec![b' '; 72]; MAXTAPE + 1],
        zatab: vec![0.0; MAXTAPE + 1],
        izatab: vec![0; MAXTAPE + 1],
        nunit: 0,
        istat1,
        philip_name: *b"ENDFB.IN1",
        ix1: 0,
        ix2: 0,
    };
    m.io.both(format!("  Merge ENDF/B Data into MAT/MF/MT Order (MERGER 2026-1)\n{}", star()));
    if m.istat1 == 1 {
        m.io.both("\n\n ERROR - Opening MERGER.INP input parameter file\n\n");
        m.io.enderror();
    }
    m.io.both(format!(" Interpretation of Input Parameters\n{}", star()));
    m.io.otape = 10;
    let l = m.io.u.read_req(2);
    m.fileout = fio::rd_a(&l, 0, 72);
    if blank(&m.fileout) {
        m.fileout = fio::rd_a(b"ENDFB.OUT", 0, 72);
    }
    let s = L::new().s(" Output Filename\n").s(&dash()).s("\n ").b(&m.fileout).s("\n").s(&dash());
    m.io.both(&s);
    let l = m.io.u.read_req(2);
    let kcard = fio::rd_a(&l, 0, 66);
    let mut ktape = fio::rd_i(&l, 66, 4, 2);
    m.igethow = fio::rd_i(&l, 70, 2, 2);
    if m.igethow != 0 {
        m.igethow = 1;
    }
    let messy = [" MAT", "  ZA"];
    let s = L::new()
        .s(" Merged File Label\n")
        .s(&dash())
        .s("\n ")
        .b(&kcard)
        .s(&i(ktape, 4))
        .s("\n")
        .s(&dash());
    m.io.lst(&s);
    m.io.lst(format!(" Retrieval Criteria-----------       {}", messy[m.igethow as usize]));
    fio::out(&s);
    fio::out(format!(" Retrieval Criteria-----------       {}", messy[m.igethow as usize]));
    m.io.both(format!("{d}\n Unit Input Filenames\n{d}", d = dash()));
    m.io.itape = 10;
    let mut nunit = 1usize;
    let mut found_end = false;
    while nunit <= MAXTAPE {
        let l = m.io.u.read_req(2);
        let mut f = fio::rd_a(&l, 0, 72);
        let t = String::from_utf8_lossy(&f).trim_end().to_string();
        if t == "END" || t == "end" {
            found_end = true;
            break;
        }
        m.io.itape += 1;
        m.ntape[nunit] = m.io.itape;
        m.philip(&mut f);
        m.filein[nunit] = f;
        let s = L::new().s(&i(m.ntape[nunit], 5)).s(" ").b(&m.filein[nunit]);
        m.io.both(&s);
        nunit += 1;
    }
    if !found_end {
        m.io.both(format!("{}\n ERROR - Over{} Units Specified for Merge\n Execution Terminated", dash(), i(MAXTAPE as i32, 5)));
        m.io.enderror();
    }
    m.nunit = nunit - 1;
    m.ask4it();
    if m.nunit < 1 {
        m.io.both(format!("{}\n Two or More Units Required for Merge\n Execution Terminated", dash()));
        m.io.enderror();
    }
    if m.nunit == 1 && m.mywish == 0 {
        m.io.both(format!("{}\n No Units Specified for Merge\n Execution Terminated", dash()));
        m.io.enderror();
    }
    m.filio2();
    m.io.both(format!(
        "{d}\n{s}\n Retrieved Data\n{s}\n     ZA  MAT MF  MT  Lines  Unit Request Message\n{}  Number\n{d}",
        " ".repeat(32),
        d = dash(),
        s = star()
    ));
    let mut tapcrd = 1;
    let mut matlas = 0;
    let mut noseq;
    let mut matout = 0;
    let mut mtout = 0;
    let mut mtdup = 0;
    let mut mtsort = 0;
    let mut kunit = 0;
    let nun = m.nunit;
    let mut icard = vec![[b' '; 66]; nun + 1];
    let mut ncard = vec![[0i32; 3]; nun + 1];
    let otape = 10;
    for iu in 1..=nun {
        let t = m.ntape[iu];
        m.io.itape = t;
        let first = match m.io.u.read(t) {
            Rd::Line(l) => rdi(&l, 66, 4).map(|v| (fio::rd_a(&l, 0, 66), v)),
            Rd::Eof => None,
        };
        let mut empty = true;
        if let Some((c, n1)) = first {
            icard[iu].copy_from_slice(&c);
            ncard[iu][0] = n1;
            empty = false;
            kunit += 1;
            if kunit <= 1 {
                if blank(&kcard) {
                    if ktape <= 0 {
                        ktape = ncard[1][0];
                    }
                    let c1 = icard[1];
                    m.io.u.write(otape, w40(&c1, [ktape, 0, 0], 0));
                } else {
                    if ktape <= 0 {
                        ktape = ncard[1][0];
                    }
                    m.io.u.write(otape, w40(&kcard, [ktape, 0, 0], 0));
                }
            }
            // label 110
            loop {
                let l = m.io.u.read_req(t);
                let (c, n) = rd40(&l, t);
                icard[iu] = c;
                ncard[iu] = n;
                m.zatab[iu] = in9(&icard[iu][..11]);
                if ncard[iu][2] > 0 {
                    break;
                }
                if ncard[iu][0] < 0 {
                    empty = true;
                    break;
                }
            }
        }
        if empty {
            ncard[iu] = [99999, 0, 0];
            m.zatab[iu] = 999999.0;
        }
        m.izatab[iu] = fint(m.zatab[iu]);
    }
    if kunit <= 0 {
        let e = format!(" {}", "=".repeat(72));
        m.io.both(format!("{e}\n ERROR - ALL data files are EMPTY = Check input and files.\n         Execution terminated\n{e}"));
        m.io.enderror();
    }
    // SAVEd locals
    let mut kwish = 0usize;
    let mut mflast = 0;
    let mut mtlast = 0;
    let mut lstape = 0;
    let mut matcrd = 0;
    let mut imadup;
    let (mut mat, mut mf, mut mt) = (0, 0, 0);
    noseq = 0;
    'l160: loop {
        let mut ismall = 1usize;
        if nun >= 2 {
            for iu in 2..=nun {
                for j in 0..3 {
                    if ncard[iu][j] < ncard[ismall][j] {
                        ismall = iu;
                        break;
                    }
                    if ncard[iu][j] > ncard[ismall][j] {
                        break;
                    }
                }
            }
        }
        let t = m.ntape[ismall];
        m.io.itape = t;
        #[derive(PartialEq)]
        enum G {
            L250,
            L310,
            L370,
            L400,
        }
        let mut g = G::L250;
        if m.mywish > 0 {
            let mut imok = m.igethow;
            let mut hit = false;
            kwish = 1;
            while kwish <= m.mywish {
                let w = m.iwish[kwish];
                let ok;
                if m.igethow != 0 {
                    ok = !(m.izatab[ismall] < w[0][0] || m.izatab[ismall] > w[1][0]);
                } else if ncard[ismall][0] > w[1][0] {
                    ok = false;
                } else {
                    imok = 1;
                    ok = ncard[ismall][0] >= w[0][0];
                }
                if ok && (1..3).all(|j| !(ncard[ismall][j] < w[0][j] || ncard[ismall][j] > w[1][j])) {
                    hit = true;
                    break;
                }
                kwish += 1;
            }
            if !hit {
                g = if imok <= 0 { G::L400 } else { G::L310 };
            }
        }
        if g == G::L250 {
            mat = ncard[ismall][0];
            mf = ncard[ismall][1];
            mt = ncard[ismall][2];
            #[derive(PartialEq)]
            enum H {
                L280,
                L300,
                L320,
                L330,
                L340,
            }
            let h;
            if matlas <= 0 {
                h = H::L320;
            } else if mat < matlas {
                h = H::L280;
            } else if mat > matlas {
                noseq = 0;
                m.io.u.write(otape, w460([matlas, 0, 0], noseq));
                m.io.u.write(otape, w460([0, 0, 0], noseq));
                matcrd += 2;
                tapcrd += matcrd;
                noseq = 0;
                matout += 1;
                m.io.both(format!("{d}\n          MAT Lines{}\n{d}", i(matcrd, 7), d = dash()));
                h = H::L320;
            } else if mf < mflast {
                h = H::L280;
            } else if mf > mflast {
                noseq = 0;
                m.io.u.write(otape, w460([matlas, 0, 0], noseq));
                matcrd += 1;
                h = H::L330;
            } else if mt == mtlast {
                h = if t != lstape { H::L300 } else { H::L340 };
            } else if mt > mtlast {
                h = H::L330;
            } else {
                h = H::L280;
            }
            let idl = |m: &Merger, ism: usize| format!("{}{}{}{}", i(m.izatab[ism], 7), i(mat, 5), i(mf, 3), i(mt, 4));
            match h {
                H::L280 => {
                    let s = format!("{}       {}{} Section is NOT in MAT/MF/MT Order (Skipped)", idl(&m, ismall), i(t, 6), i(m.nowish[kwish], 8));
                    m.io.both(&s);
                    mtsort += 1;
                    g = G::L310;
                }
                H::L300 => {
                    let s = format!("{}       {}{} Duplicate MAT/MF/MT (Skipped)", idl(&m, ismall), i(t, 6), i(m.nowish[kwish], 8));
                    m.io.both(&s);
                    mtdup += 1;
                    g = G::L310;
                }
                _ => {
                    if h == H::L320 {
                        matcrd = 0;
                    }
                    imadup = if h == H::L340 { 1 } else { 0 };
                    // label 340
                    noseq = nxtseq(noseq);
                    m.io.u.write(otape, w40(&icard[ismall], ncard[ismall], noseq));
                    let mut mtcard = 1;
                    loop {
                        let l = m.io.u.read_req(t);
                        let (c, n) = rd40(&l, t);
                        icard[ismall] = c;
                        ncard[ismall] = n;
                        noseq = nxtseq(noseq);
                        if ncard[ismall][2] <= 0 {
                            break;
                        }
                        m.io.u.write(otape, w40(&icard[ismall], ncard[ismall], noseq));
                        mtcard += 1;
                    }
                    noseq = if ncard[ismall][1] != 0 { 99999 } else { 0 };
                    m.io.u.write(otape, w460(ncard[ismall], noseq));
                    mtout += 1;
                    mtcard += 1;
                    matcrd += mtcard;
                    let s = if imadup == 0 {
                        format!("{}{}{}{}", idl(&m, ismall), i(mtcard, 7), i(t, 6), i(m.nowish[kwish], 8))
                    } else {
                        format!(
                            "{}{}{}{} WARNING - Duplicate MAT/MF/MT On Same File (Copied)",
                            idl(&m, ismall),
                            i(mtcard, 7),
                            i(t, 6),
                            i(m.nowish[kwish], 8)
                        )
                    };
                    m.io.both(&s);
                    m.mtgot[kwish] += 1;
                    m.icardgot[kwish] += mtcard;
                    if m.matnow[kwish] != mat {
                        m.matgot[kwish] += 1;
                    }
                    m.matnow[kwish] = mat;
                    g = G::L370;
                }
            }
        }
        if g == G::L310 {
            loop {
                let l = m.io.u.read_req(t);
                ncard[ismall] = [fio::rd_i(&l, 66, 4, t), fio::rd_i(&l, 70, 2, t), fio::rd_i(&l, 72, 3, t)];
                if ncard[ismall][2] <= 0 {
                    break;
                }
            }
            g = G::L370;
        }
        if g == G::L370 {
            // MAT/MF/MT are those of the last section looked at in label 250
            matlas = mat;
            mflast = mf;
            mtlast = mt;
            lstape = t;
            loop {
                let l = m.io.u.read_req(t);
                let (c, n) = rd40(&l, t);
                icard[ismall] = c;
                ncard[ismall] = n;
                m.zatab[ismall] = in9(&icard[ismall][..11]);
                if ncard[ismall][2] > 0 {
                    m.izatab[ismall] = fint(m.zatab[ismall]);
                    continue 'l160;
                }
                if ncard[ismall][0] < 0 {
                    break;
                }
            }
        }
        // label 400
        ncard[ismall] = [99999, 0, 0];
        m.zatab[ismall] = 999999.0;
        m.izatab[ismall] = fint(m.zatab[ismall]);
        if (1..=nun).any(|iu| ncard[iu][0] < 10000) {
            continue 'l160;
        }
        break;
    }
    if matlas <= 0 {
        m.io.both(format!(" WARNING - No MAT/MF/MT Combination Satisfies Retrieval Criteria\n{}", dash()));
        m.io.u.write(otape, format!("{}-1 0  0    0", " ".repeat(68)));
        m.io.enderror();
    }
    noseq = 0;
    m.io.u.write(otape, w460([matlas, 0, 0], noseq));
    m.io.u.write(otape, w460([0, 0, 0], noseq));
    matcrd += 2;
    m.io.u.write(otape, format!("{}-1 0  0    0", " ".repeat(68)));
    tapcrd += matcrd + 1;
    matout += 1;
    let s = format!(
        "{d}\n          MAT Lines{}\n{d}\n         Tape Lines{}\n{d}\n{s}\n Retrieval Statistics\n{s}\n Merged Materials---{}\n Merged Sections----{}\n Merged Lines-------{}",
        i(matcrd, 7),
        i(tapcrd, 7),
        i(matout, 11),
        i(mtout, 11),
        i(tapcrd, 11),
        d = dash(),
        s = star()
    );
    let mut lines = vec![s];
    if mtdup > 0 {
        lines.push(format!(" Duplicate Sections-{} (Skipped)", i(mtdup, 11)));
    }
    if mtsort > 0 {
        lines.push(format!(" Unsorted Sections--{} (ENDF/B Tapes were NOT in MAT/MF/MT Order)", i(mtsort, 11)));
    }
    for l in &lines {
        m.io.lst(l);
    }
    for l in &lines {
        fio::out(l);
    }
    if m.mywish > 0 {
        let hdr = if m.igethow == 0 {
            format!("{d}\n Retrieved Data for Each Request\n{d}\n Request *   Lower Limit *   Upper Limit *            Retrieved\n  Number *    MAT MF  MT *    MAT MF  MT *   Materials    Sections  Lines\n{d}", d = dash())
        } else {
            format!("{d}\n Retrieved Data for Each Request\n{d}\n Request *   Lower Limit *   Upper Limit *            Retrieved\n  Number *     ZA MF  MT *     ZA MF  MT *   Materials    Sections  Lines\n{d}", d = dash())
        };
        m.io.both(&hdr);
        for k in 1..=m.mywish {
            let w = m.iwish[k];
            let s = format!(
                "{} *{}{}{} *{}{}{} *{}{}{}",
                i(m.nowish[k], 8),
                i(w[0][0], 7),
                i(w[0][1], 3),
                i(w[0][2], 4),
                i(w[1][0], 7),
                i(w[1][1], 3),
                i(w[1][2], 4),
                i(m.matgot[k], 12),
                i(m.mtgot[k], 12),
                i(m.icardgot[k], 7)
            );
            m.io.both(&s);
        }
        let s = format!(
            "{d}\n  Totals{}{}{}{}\n{d}\n Note..The Line Count of Retrieved Data for Each Request Does\n       NOT Include FEND, MEND AND TEND Lines. Therefore the  \n       Total Number of Lines Actually Retrieved and Written  \n       on the Merged File Will be More than the Sum of the   \n       Number of Lines for All Requests.                     ",
            " ".repeat(34),
            i(matout, 12),
            i(mtout, 12),
            i(tapcrd, 7),
            d = dash()
        );
        m.io.both(&s);
    }
    m.io.endit();
}

/// `NXTSEQ`.
pub fn nxtseq(n: i32) -> i32 {
    let nn = n + 1;
    if nn == 100000 {
        1
    } else {
        nn
    }
}
