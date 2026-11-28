//! Port of MIXER (PREPRO 2026-1): mix the cross sections of up to 10
//! materials into a composite material. Input `MIXER.INP`, listing
//! `MIXER.LST`.

use crate::endfio::{fint, in9, out9, Endf};
use crate::ffmt::{f, i};
use crate::fio::{self, Rd, Val, L};
use crate::linear::{blank, rdi};

const MAXPOINT: usize = 120000;
const MAXPART: usize = 10;
const NT: usize = MAXPART + 1;

fn dash() -> String {
    format!(" {}", "-".repeat(78))
}

struct Mx {
    io: Endf,
    awre: [f64; NT + 1],
    dense: [f64; NT + 1],
    tempk: [f64; NT + 1],
    qvalue1: [f64; NT + 1],
    qvalue2: [f64; NT + 1],
    atom10: [f64; NT + 1],
    fracta: [f64; NT + 1],
    fractg: [f64; NT + 1],
    izaget: [i32; NT + 1],
    matget: [i32; NT + 1],
    mtget: [i32; NT + 1],
    iverse: [i32; NT + 1],
    nmat: usize,
    ixylow: [usize; NT + 1],
    ixyhi: [usize; NT + 1],
    iscr: [i32; NT + 1],
    npage: usize,
    npagem: usize,
    nptab: [usize; NT + 1],
    iptab: [usize; NT + 1],
    iuse: [usize; NT + 1],
    iusen: usize,
    xa: [f64; NT + 1],
    ya: [f64; NT + 1],
    yb: [f64; NT + 1],
    nbt: [i32; 101],
    int: [i32; 101],
    xtab: Vec<Vec<f64>>,
    ytab: Vec<Vec<f64>>,
    igoof: i32,
    izaout: i32,
    matout: i32,
    mfout: i32,
    mtout: i32,
    mytype: i32,
    denout: f64,
    atombc: f64,
    awrout: f64,
    title: Vec<u8>,
    field6: [[u8; 11]; 6],
    namein: Vec<u8>,
    nameout: Vec<u8>,
}

impl Mx {
    fn both<S: AsRef<[u8]>>(&mut self, s: S) {
        self.io.both(s);
    }

    fn err_eof(&mut self) -> ! {
        self.both("\n\n ERROR - End of file reading MIXER.INP\n\n");
        self.io.enderror();
    }

    fn rd(&mut self) -> Option<Vec<u8>> {
        match self.io.u.read(2) {
            Rd::Line(l) => Some(l),
            Rd::Eof => None,
        }
    }

    fn readin(&mut self) {
        self.igoof = 0;
        let Some(l1) = self.rd() else { self.err_eof() };
        let Some(l2) = self.rd() else { self.err_eof() };
        self.title[..68].copy_from_slice(&fio::rd_a(&l1, 0, 68));
        self.title[68..136].copy_from_slice(&fio::rd_a(&l2, 0, 68));
        let s = L::new()
            .s(" Problem Identification\n")
            .s(&dash())
            .s("\n ")
            .b(&self.title[..68])
            .s("\n ")
            .b(&self.title[68..136])
            .s("\n")
            .s(&dash());
        self.both(&s);
        let Some(l) = self.rd() else { self.err_eof() };
        self.namein = fio::rd_a(&l, 0, 72);
        if blank(&self.namein) {
            self.namein = fio::rd_a(b"ENDFB.IN", 0, 72);
        }
        let Some(l) = self.rd() else { self.err_eof() };
        self.nameout = fio::rd_a(&l, 0, 72);
        if blank(&self.nameout) {
            self.nameout = fio::rd_a(b"ENDFB.OUT", 0, 72);
        }
        let s = L::new()
            .s(" ENDF/B Input and Output Data Filenames\n")
            .s(&dash())
            .s("\n ")
            .b(&self.namein)
            .s("\n ")
            .b(&self.nameout);
        self.both(&s);
        self.both(format!("{d}\n Requested Composition (Input Parameters)\n{d}\n     ZA MF  MT\n{d}", d = dash()));
        let Some(l) = self.rd() else { self.err_eof() };
        let v = [rdi(&l, 0, 11), rdi(&l, 11, 6), rdi(&l, 17, 2), rdi(&l, 19, 3), rdi(&l, 22, 11)];
        if v.iter().any(|x| x.is_none()) {
            self.err_eof();
        }
        self.izaout = v[0].unwrap();
        self.matout = v[1].unwrap();
        self.mfout = v[2].unwrap();
        self.mtout = v[3].unwrap();
        self.mytype = v[4].unwrap();
        if self.izaout <= 0 {
            self.izaout = 1;
        }
        if self.matout <= 0 {
            self.matout = 1;
        }
        if self.mfout != 3 && self.mfout != 23 {
            self.mfout = 3;
        }
        if self.mtout <= 0 {
            self.mtout = 1;
        }
        let mut sumin = 0.0;
        let mut nmat = 1usize;
        let mut go100 = false;
        while nmat <= 10 {
            let rec = self.rd().and_then(|l| Some((rdi(&l, 0, 11)?, rdi(&l, 18, 4)?, fio::rd_a(&l, 22, 11))));
            let Some((za, mt, fld)) = rec else {
                go100 = true;
                break;
            };
            self.izaget[nmat] = za;
            self.mtget[nmat] = mt;
            self.field6[0].copy_from_slice(&fld);
            self.dense[nmat] = in9(&self.field6[0]);
            if za == 0 && mt == 0 && self.dense[nmat] == 0.0 {
                go100 = true;
                break;
            }
            let ok = za > 0 && mt > 0 && self.dense[nmat] > 0.0;
            out9(self.dense[nmat], &mut self.field6[0]);
            let head = format!("{}{}{}", i(za, 7), i(self.mfout, 3), i(mt, 4));
            let tail = if !ok {
                self.igoof = 1;
                " (ERROR in Data - All MUST be Positive)"
            } else if self.mytype == 0 {
                " Grams"
            } else {
                " Atoms"
            };
            let s = L::new().s(&head).b(&self.field6[0]).s(tail);
            self.both(&s);
            sumin += self.dense[nmat];
            nmat += 1;
        }
        if !go100 {
            // after 10 materials one more line must be blank
            let rec = self.rd().and_then(|l| Some((rdi(&l, 0, 11)?, rdi(&l, 18, 4)?, fio::rd_a(&l, 22, 11))));
            if let Some((za, mt, fld)) = rec {
                self.field6[0].copy_from_slice(&fld);
                let dd = in9(&self.field6[0]);
                if !(za <= 0 && mt <= 0 && dd <= 0.0) {
                    self.both("\n\n\n Over 10 Materials in Composition - Execution Terminated.");
                    self.igoof = 1;
                    return;
                }
            }
            nmat = 11;
        }
        // label 100
        if self.igoof > 0 {
            self.nmat = nmat;
            return;
        }
        nmat -= 1;
        self.nmat = nmat;
        if nmat < 1 {
            self.both(format!(" Material Count={} (MUST be at Least 1)", i(nmat as i32, 2)));
            self.igoof = 1;
            return;
        }
        out9(sumin, &mut self.field6[0]);
        let unit = if self.mytype == 0 { " Grams" } else { " Atoms" };
        let s = L::new()
            .s(&dash())
            .s("\n")
            .s(&format!("{}{}{}", i(self.izaout, 7), i(self.mfout, 3), i(self.mtout, 4)))
            .b(&self.field6[0])
            .s(unit)
            .s("\n")
            .s(&dash())
            .s("\n Processing ENDF Input Data\n")
            .s(&dash())
            .s("\n     ZA  MAT MF  MT     Kelvin     Q-Value Points AMU-Weight Atoms/b-cm *\n")
            .s(&" ".repeat(40))
            .s("eV\n")
            .s(&dash());
        self.both(&s);
    }

    /// Format 210 of GETIT.
    fn line210(&mut self, im: usize) {
        out9(self.tempk[im], &mut self.field6[0]);
        out9(self.qvalue2[im], &mut self.field6[1]);
        out9(self.awre[im], &mut self.field6[2]);
        out9(self.atom10[im], &mut self.field6[3]);
        let s = L::new()
            .s(&format!("{}{}{}{}", i(self.izaget[im], 7), i(self.matget[im], 5), i(self.mfout, 3), i(self.mtget[im], 4)))
            .b(&self.field6[0])
            .s(" ")
            .b(&self.field6[1])
            .s(&i(self.nptab[im] as i32, 7))
            .b(&self.field6[2])
            .b(&self.field6[3]);
        self.both(&s);
    }

    fn file1in(&mut self, im: usize) {
        let (_c1, _c2, _l1, _l2, n1, n2) = self.io.cardi_ld();
        self.iverse[im] = 4;
        if n1 > 0 {
            return;
        }
        let n2x = n2;
        self.io.cardi_ld();
        self.iverse[im] = 5;
        if n2x <= 0 {
            return;
        }
        let (c1, ..) = self.io.cardi_ld();
        self.iverse[im] = 6;
        self.tempk[im] = c1;
    }

    fn odump(&mut self, it: usize) {
        let n = self.npage;
        let mut rec = Vec::with_capacity(2 * n);
        for k in 1..=n {
            rec.push(Val::D(self.xtab[it][k]));
        }
        for k in 1..=n {
            rec.push(Val::D(self.ytab[it][k]));
        }
        let u = self.iscr[it];
        self.io.u.uwrite(u, rec);
    }

    fn idump(&mut self, it: usize) {
        let n = self.npage;
        let u = self.iscr[it];
        let rec = self.io.u.uread_req(u);
        for k in 1..=n {
            self.xtab[it][k] = rec[k - 1].d();
            self.ytab[it][k] = rec[n + k - 1].d();
        }
    }

    fn ipag10(&mut self, it: usize) {
        let factor = self.atom10[it];
        let nscr = self.iscr[it];
        self.io.elast = 0.0;
        let n2in = self.nptab[it];
        let mut ihigh = 0;
        let mut npt = 1;
        while npt <= n2in {
            let mut nnpt = npt + self.npagem;
            if nnpt > n2in {
                nnpt = n2in;
            }
            ihigh = nnpt - npt + 1;
            let mut x = std::mem::take(&mut self.xtab[it]);
            let mut y = std::mem::take(&mut self.ytab[it]);
            self.io.pointi(&mut x[1..], &mut y[1..], ihigh);
            for k in 1..=ihigh {
                y[k] *= factor;
            }
            self.xtab[it] = x;
            self.ytab[it] = y;
            if n2in > self.npage {
                if npt == 1 {
                    self.io.u.rewind(nscr);
                }
                self.odump(it);
            }
            npt += self.npage;
        }
        let (xm, ym) = (self.xtab[it][ihigh], self.ytab[it][ihigh]);
        self.io.maxie2(xm, ym);
        if n2in <= self.npage {
            self.ixylow[it] = 0;
            self.ixyhi[it] = n2in;
            return;
        }
        self.io.u.endfile(nscr);
        self.io.u.rewind(nscr);
        self.idump(it);
        self.ixylow[it] = 0;
        self.ixyhi[it] = self.npage;
    }

    fn opag10(&mut self, ii: usize, it: usize) -> (f64, f64) {
        if !(ii > 0 && ii <= self.nptab[it]) {
            let s = format!(" OPAG10..Index={} (MUST be 1 to{})", i(ii as i32, 6), i(self.nptab[it] as i32, 6));
            self.io.lst(&s);
            return (0.0, 0.0);
        }
        if ii <= self.ixylow[it] {
            let s = format!(" OPAG10..Index={} Less Than{} on Scratch{}", i(it as i32, 6), i(ii as i32, 6), i(self.ixylow[it] as i32, 3));
            self.io.lst(&s);
            return (0.0, 0.0);
        }
        while ii > self.ixyhi[it] {
            self.ixylow[it] = self.ixyhi[it];
            self.ixyhi[it] += self.npage;
            self.idump(it);
        }
        let ic = ii - self.ixylow[it];
        (self.xtab[it][ic], self.ytab[it][ic])
    }

    fn getit(&mut self) {
        self.io.maxie0();
        self.igoof = 0;
        for im in 1..=self.nmat {
            self.nptab[im] = 0;
            self.matget[im] = 0;
            self.iverse[im] = 6;
            self.qvalue1[im] = 0.0;
            self.qvalue2[im] = 0.0;
            self.tempk[im] = 0.0;
            self.atom10[im] = 0.0;
            self.awre[im] = 0.0;
        }
        self.atombc = 0.0;
        self.denout = 0.0;
        self.io.skip1();
        let nmat = self.nmat;
        loop {
            // label 20
            self.io.conti();
            if self.io.h.mt <= 0 {
                if self.io.h.mat < 0 {
                    break;
                }
                continue;
            }
            let iza = fint(self.io.h.c1);
            let mut imat = 1;
            while imat <= nmat {
                if iza == self.izaget[imat] {
                    break;
                }
                imat += 1;
            }
            if imat > nmat {
                self.io.skipm();
                continue;
            }
            // label 50
            let awrin = self.io.h.c2;
            for k in 1..=nmat {
                if iza == self.izaget[k] {
                    if self.matget[k] <= 0 {
                        self.awre[k] = 1.008665 * awrin;
                        if self.mytype == 0 {
                            self.atom10[k] = 0.60247 * self.dense[k] / self.awre[k];
                        } else {
                            self.atom10[k] = self.dense[k];
                            self.dense[k] = self.awre[k] * self.atom10[k] / 0.60247;
                        }
                        // the Fortran sums with IMAT (the first match), not I
                        self.atombc += self.atom10[imat];
                        self.denout += self.dense[imat];
                    }
                    self.matget[k] = self.io.h.mat;
                }
            }
            let mfh = self.io.h.mf;
            if self.mfout == 3 && mfh == 1 {
                self.file1in(imat);
                for k in 1..=nmat {
                    if iza == self.izaget[k] {
                        self.iverse[k] = self.iverse[imat];
                        self.tempk[k] = self.tempk[imat];
                    }
                }
                self.io.skipf();
                continue;
            }
            if mfh < self.mfout {
                self.io.skipf();
                continue;
            }
            if mfh > self.mfout {
                self.io.skipm();
                continue;
            }
            // label 80
            let mth = self.io.h.mt;
            let mut im = 0;
            for k in 1..=nmat {
                if self.izaget[k] != iza || self.mtget[k] != mth {
                    continue;
                }
                if self.nptab[k] == 0 {
                    im = k;
                    break;
                }
            }
            if im == 0 {
                self.io.skips();
                continue;
            }
            // label 100
            let (q1, q2, l1, l2, n1, n2) = self.io.cardi();
            self.qvalue1[im] = q1;
            self.qvalue2[im] = q2;
            self.io.ld.l1 = l1;
            self.io.ld.l2 = l2;
            self.io.ld.n1 = n1;
            self.io.ld.n2 = n2;
            let n1u = n1.max(0) as usize;
            let mut nbt = self.nbt;
            let mut int = self.int;
            self.io.terpi(&mut nbt[1..], &mut int[1..], n1u.min(100));
            self.nbt = nbt;
            self.int = int;
            self.io.terptest(&self.nbt[1..], &self.int[1..], n1, n2, 3);
            self.nptab[im] = n2.max(0) as usize;
            self.matget[im] = self.io.h.mat;
            self.line210(im);
            let (m, f, t) = (self.io.h.mat, self.io.h.mf, self.io.h.mt);
            self.io.maxie1(m, f, t);
            self.ipag10(im);
            if (1..=nmat).any(|k| self.nptab[k] == 0) {
                self.io.skips();
                continue;
            }
            if self.igoof != 0 {
                self.both(format!("{}\n Execution Terminated due to ERROR in Evaluated Data", dash()));
                self.io.enderror();
            }
            self.awrout = 0.60247 * self.denout / self.atombc;
            return;
        }
        // label 120
        self.igoof = 0;
        let mut ifound = 0;
        for im in 1..=nmat {
            if self.nptab[im] > 0 {
                ifound += 1;
                continue;
            }
            self.line210(im);
            if self.matget[im] <= 0 {
                let s = format!("{}          ERROR - ENDF Data Does not include ZA={}", i(self.izaget[im], 7), i(self.izaget[im], 7));
                self.both(&s);
                self.igoof += 1;
            } else {
                let s = format!(
                    "{}{}     WARNING - ZA={} Does not include MF/MT={}/{}",
                    i(self.izaget[im], 7),
                    i(self.matget[im], 5),
                    i(self.izaget[im], 7),
                    i(self.mfout, 2),
                    i(self.mtget[im], 3)
                );
                self.both(&s);
            }
        }
        if self.igoof > 0 {
            self.both(format!("{}\n ERROR - Missing ZA - Check ENDF Data", dash()));
            self.io.enderror();
        }
        if ifound == 0 {
            self.both(format!("{}\n ERROR - NO data found to MIX", dash()));
            self.io.enderror();
        }
        self.awrout = 0.60247 * self.denout / self.atombc;
        if ifound < nmat {
            self.both(format!("{}\n WARNING - Output includes ONLY data found", dash()));
        }
    }

    fn mixit(&mut self) {
        let mut xnow = 0.0f64;
        self.io.maxie3(1);
        let xmaximum = self.io.maxie.xmaximum;
        let nmat = self.nmat;
        let mut imfirst = 0;
        let mut ixtend = [0i32; NT + 1];
        for im in 1..=nmat {
            ixtend[im] = 0;
            if self.nptab[im] == 0 {
                continue;
            }
            if imfirst == 0 {
                imfirst = im;
            }
            if self.io.maxie.xmaxie[im] < xmaximum && self.io.maxie.ymaxie[im] > 0.0 {
                ixtend[im] = 1;
            }
        }
        self.both(format!(
            "{d}\n Processing ENDF Output MIX\n{d}\n     ZA  MAT MF  MT     Kelvin     Q-Value Points AMU-Weight Atoms/b-cm *\n{}eV",
            " ".repeat(40),
            d = dash()
        ));
        let un = self.iusen;
        self.nptab[un] = 0;
        self.ixylow[un] = 0;
        self.ixyhi[un] = 0;
        let mut iptn = 0usize;
        let nscr = self.iscr[un];
        for im in 1..=nmat {
            self.iptab[im] = 1;
            if self.nptab[im] == 0 {
                continue;
            }
            let iu = self.iuse[im];
            (self.xa[im], self.ya[im]) = self.opag10(1, iu);
        }
        loop {
            // label 40
            let xlast = xnow;
            xnow = self.xa[imfirst];
            for im in 1..=nmat {
                if self.nptab[im] == 0 {
                    continue;
                }
                if self.xa[im] < xnow {
                    xnow = self.xa[im];
                }
            }
            for im in 1..=nmat {
                if self.nptab[im] == 0 {
                    continue;
                }
                if !(self.xa[im] > xnow) {
                    self.yb[im] = self.ya[im];
                    self.iptab[im] += 1;
                    if self.iptab[im] <= self.nptab[im] {
                        let (p, iu) = (self.iptab[im], self.iuse[im]);
                        (self.xa[im], self.ya[im]) = self.opag10(p, iu);
                    } else {
                        if ixtend[im] != 0 {
                            self.xa[im] = self.io.maxie.xmaxie[im];
                            ixtend[im] = 0;
                        } else {
                            self.xa[im] = xmaximum;
                        }
                        self.ya[im] = 0.0;
                    }
                    continue;
                }
                if self.iptab[im] <= 1 {
                    self.yb[im] = 0.0;
                    continue;
                }
                if !(self.xa[im] > xlast) {
                    self.yb[im] = self.ya[im];
                    continue;
                }
                self.yb[im] = ((xnow - xlast) * self.ya[im] + (self.xa[im] - xnow) * self.yb[im]) / (self.xa[im] - xlast);
            }
            if iptn >= self.npage {
                if self.nptab[un] == 0 {
                    self.io.u.rewind(nscr);
                }
                self.odump(un);
                self.nptab[un] += iptn;
                iptn = 0;
            }
            iptn += 1;
            let mut imthru = 0;
            let mut sum = 0.0;
            for im in 1..=nmat {
                if self.nptab[im] == 0 {
                    continue;
                }
                sum += self.yb[im];
                if self.iptab[im] <= self.nptab[im] {
                    imthru = 1;
                }
            }
            self.xtab[un][iptn] = xnow;
            self.ytab[un][iptn] = sum;
            if xnow < xlast {
                let ii = iptn + self.nptab[un];
                out9(xnow, &mut self.field6[0]);
                let s = L::new()
                    .s(&" ".repeat(16))
                    .s(" Mixed Energy NOT in Ascending Energy Order")
                    .s(&i(ii as i32, 6))
                    .b(&self.field6[0]);
                self.both(&s);
            }
            if imthru > 0 {
                continue;
            }
            break;
        }
        self.nptab[un] += iptn;
        self.ixylow[un] = 0;
        if self.nptab[un] <= self.npage {
            self.ixyhi[un] = self.nptab[un];
            return;
        }
        self.odump(un);
        self.ixyhi[un] = self.npage;
        self.io.u.endfile(nscr);
        self.io.u.rewind(nscr);
        self.idump(un);
    }

    fn outit(&mut self) {
        let tlabel = format!("{:<66}", " Composite Material in ENDF/B-6 Format");
        let texts = [
            " ***************** Program MIXER (VERSION 2026-1) ************",
            " ----------------------------------------                     ",
            " Composition                                                  ",
            " ----------------------------------------                     ",
            " Isotope MF   MT  Atom-Fract   Gram-Fract                     ",
            " ----------------------------------------                     ",
            "   12345 12  123 12345678901  12345678901                     ",
            "                                1        451                  ",
            "                                                              ",
        ];
        let mut pd = [[b' '; 66]; 10];
        for (k, s) in texts.iter().enumerate() {
            pd[k + 1][..s.len()].copy_from_slice(s.as_bytes());
        }
        for k in 62..66 {
            pd[1][k] = b'*';
        }
        let nmat = self.nmat;
        let (mut q1, mut q2, mut tk) = (1.0e12f64, 1.0e12f64, 1.0e12f64);
        for k in 1..=nmat {
            if self.nptab[k] > 0 {
                if self.qvalue1[k] < q1 {
                    q1 = self.qvalue1[k];
                }
                if self.qvalue2[k] < q2 {
                    q2 = self.qvalue2[k];
                }
                if self.tempk[k] < tk {
                    tk = self.tempk[k];
                }
            }
        }
        out9(tk, &mut self.field6[0]);
        out9(q2, &mut self.field6[1]);
        out9(self.awrout, &mut self.field6[2]);
        out9(self.atombc, &mut self.field6[3]);
        let un = self.iusen;
        let s = L::new()
            .s(&dash())
            .s("\n")
            .s(&format!("{}{}{}{}", i(self.izaout, 7), i(self.matout, 5), i(self.mfout, 3), i(self.mtout, 4)))
            .b(&self.field6[0])
            .s(" ")
            .b(&self.field6[1])
            .s(&i(self.nptab[un] as i32, 7))
            .b(&self.field6[2])
            .b(&self.field6[3])
            .s("\n")
            .s(&dash())
            .s("\n * - Note b-cm = 1.0E-24 Cubic Centimeters\n")
            .s(&dash());
        self.both(&s);
        self.both(format!(
            " Atom and Gram Fractions\n{d}\n     ZA  MAT MF  MT     Atom        Gram\n{}Fraction    Fraction\n{d}",
            " ".repeat(22),
            d = dash()
        ));
        for k in 1..=nmat {
            self.fracta[k] = self.atom10[k] / self.atombc;
            self.fractg[k] = self.dense[k] / self.denout;
            let s = format!(
                "{}{}{}{}{}{}",
                i(self.izaget[k], 7),
                i(self.matget[k], 5),
                i(self.mfout, 3),
                i(self.mtget[k], 4),
                f(self.fracta[k], 12, 8),
                f(self.fractg[k], 12, 8)
            );
            self.both(&s);
        }
        let s = format!(
            "{d}\n{}{}{}{}{}{}\n{d}",
            i(self.izaout, 7),
            i(self.matout, 5),
            i(self.mfout, 3),
            i(self.mtout, 4),
            f(1.0, 12, 8),
            f(1.0, 12, 8),
            d = dash()
        );
        self.both(&s);
        let arwmix = self.awrout / 1.008665;
        let kpt = self.nptab[un];
        self.io.h.mat = self.matout;
        self.io.h.mf = 0;
        self.io.h.mt = 0;
        self.io.hollyo(tlabel.as_bytes());
        let nwd = 11 + nmat as i32;
        let zaout = self.izaout as f64;
        self.io.h.mat = self.matout;
        self.io.h.mf = 1;
        self.io.h.mt = 451;
        let nsub = if self.mfout == 23 { 0 } else { 10 };
        self.io.cardo(zaout, arwmix, 0, 0, 0, 0);
        self.io.cardo(0.0, 0.0, 0, 0, 0, 6);
        let xmass = if self.mfout == 3 { 1.0 } else { 0.0 };
        let xmaximum = self.io.maxie.xmaximum;
        self.io.cardo(xmass, xmaximum, 0, 0, nsub, 6);
        self.io.cardo(tk, 0.0, 1, 0, nwd, 2);
        self.io.hollyo(&pd[1]);
        let t1 = self.title[..66].to_vec();
        self.io.hollyo(&t1);
        let t2 = self.title[68..134].to_vec();
        self.io.hollyo(&t2);
        for k in 2..=6 {
            self.io.hollyo(&pd[k]);
        }
        for lm in 1..=nmat {
            self.fracta[lm] = self.atom10[lm] / self.atombc;
            self.fractg[lm] = self.dense[lm] / self.denout;
            intout(self.izaget[lm], &mut pd[7][3..8]);
            intout(self.mfout, &mut pd[7][9..11]);
            intout(self.mtget[lm], &mut pd[7][13..16]);
            out9(self.fracta[lm], &mut pd[7][17..28]);
            out9(self.fractg[lm], &mut pd[7][30..41]);
            self.io.hollyo(&pd[7]);
        }
        self.io.hollyo(&pd[6]);
        intout(self.izaout, &mut pd[7][3..8]);
        intout(self.mfout, &mut pd[7][9..11]);
        intout(self.mtout, &mut pd[7][13..16]);
        out9(1.0, &mut pd[7][17..28]);
        out9(1.0, &mut pd[7][30..41]);
        self.io.hollyo(&pd[7]);
        self.io.hollyo(&pd[6]);
        intout(nwd + 6, &mut pd[8][44..55]);
        intout(0, &mut pd[8][55..66]);
        self.io.hollyo(&pd[8]);
        let nline = (kpt as i32 + 11) / 3;
        intout(self.mfout, &mut pd[9][22..33]);
        intout(self.mtout, &mut pd[9][33..44]);
        intout(nline, &mut pd[9][44..55]);
        intout(0, &mut pd[9][55..66]);
        self.io.hollyo(&pd[9]);
        let mo = self.matout;
        self.io.outs(mo, 1);
        self.io.outf(mo);
        self.io.h.mf = self.mfout;
        self.io.h.mt = self.mtout;
        self.io.cardo(zaout, arwmix, 0, 0, 0, 0);
        self.io.cardo(q1, q2, 0, 0, 1, kpt as i32);
        self.io.terpo(&[kpt as i32], &[2], 1);
        let factor = 1.0 / self.atombc;
        let mut lpt = 1;
        while lpt <= kpt {
            let mut lptp2 = lpt + self.npagem;
            if lptp2 > kpt {
                lptp2 = kpt;
            }
            let mut lx = 0;
            for l3 in lpt..=lptp2 {
                lx += 1;
                let (x, y) = self.opag10(l3, un);
                self.xtab[1][lx] = x;
                self.ytab[1][lx] = factor * y;
            }
            let x = std::mem::take(&mut self.xtab[1]);
            let y = std::mem::take(&mut self.ytab[1]);
            self.io.pointo(&x[1..], &y[1..], lx);
            self.xtab[1] = x;
            self.ytab[1] = y;
            lpt += self.npage;
        }
        let mf = self.mfout;
        self.io.outs(mo, mf);
        self.io.outf(mo);
        self.io.outm();
        self.io.outt();
    }
}

/// `INTOUT`: integer right adjusted in a field of LENGTH characters.
fn intout(v: i32, field: &mut [u8]) {
    for c in field.iter_mut() {
        *c = b' ';
    }
    let mut ii = v;
    for k in (0..field.len()).rev() {
        let kk = ii / 10;
        let ll = ii - 10 * kk;
        field[k] = b'0' + ll.rem_euclid(10) as u8;
        ii = kk;
        if ii <= 0 {
            break;
        }
    }
}

pub fn run() -> ! {
    let mut io = Endf::new("MIXER");
    io.timer();
    io.inp = 2;
    io.outp = 3;
    io.itape = 10;
    io.otape = 11;
    io.u.open(3, "MIXER.LST");
    let istat1 = if io.u.open_old(2, "MIXER.INP") { 0 } else { 1 };
    let mut m = Mx {
        io,
        awre: [0.0; NT + 1],
        dense: [0.0; NT + 1],
        tempk: [0.0; NT + 1],
        qvalue1: [0.0; NT + 1],
        qvalue2: [0.0; NT + 1],
        atom10: [0.0; NT + 1],
        fracta: [0.0; NT + 1],
        fractg: [0.0; NT + 1],
        izaget: [0; NT + 1],
        matget: [0; NT + 1],
        mtget: [0; NT + 1],
        iverse: [0; NT + 1],
        nmat: 0,
        ixylow: [0; NT + 1],
        ixyhi: [0; NT + 1],
        iscr: [0; NT + 1],
        npage: MAXPOINT,
        npagem: MAXPOINT - 1,
        nptab: [0; NT + 1],
        iptab: [0; NT + 1],
        iuse: [0; NT + 1],
        iusen: 0,
        xa: [0.0; NT + 1],
        ya: [0.0; NT + 1],
        yb: [0.0; NT + 1],
        nbt: [0; 101],
        int: [0; 101],
        xtab: vec![vec![0.0; MAXPOINT + 2]; NT + 1],
        ytab: vec![vec![0.0; MAXPOINT + 2]; NT + 1],
        igoof: 0,
        izaout: 0,
        matout: 0,
        mfout: 0,
        mtout: 0,
        mytype: 0,
        denout: 0.0,
        atombc: 0.0,
        awrout: 0.0,
        title: vec![0u8; 136],
        field6: [[0u8; 11]; 6],
        namein: vec![b' '; 72],
        nameout: vec![b' '; 72],
    };
    m.both(format!(" Mix Energy Dependent ENDF/B Cross Sections (MIXER 2026-1)\n{}", dash()));
    if istat1 == 1 {
        m.both("\n\n ERROR - Opening MIXER.INP input parameter file\n\n");
        m.io.enderror();
    }
    m.readin();
    if m.igoof > 0 {
        m.io.enderror();
    }
    // FILIO2
    for k in 1..=m.nmat {
        m.iscr[k] = k as i32 + 11;
        m.iuse[k] = k;
    }
    m.iscr[11] = 22;
    m.iusen = 11;
    let no = String::from_utf8_lossy(&m.nameout).into_owned();
    m.io.u.open(11, &no);
    let ni = String::from_utf8_lossy(&m.namein).into_owned();
    let istat2 = if m.io.u.open_old(10, &ni) { 0 } else { 1 };
    for k in 1..=m.nmat {
        let u = m.iscr[k];
        m.io.u.open_scratch(u);
    }
    let u = m.iscr[11];
    m.io.u.open_scratch(u);
    if istat2 == 1 {
        let s = L::new().s("\n\n ERROR - Opening ENDF/B formatted file\n ").b(&m.namein).s("\n\n");
        m.both(&s);
        m.io.enderror();
    }
    m.getit();
    m.mixit();
    m.outit();
    m.io.maxie4(1);
    m.io.endit();
}
