//! RESXSR: a RESXS file (CCCC-like pointwise resonance cross sections:
//! elastic, fission and capture against temperature) from PENDF tapes,
//! ported from NJOY2016 `resxsr.f90` (module resxsm).
//!
//! The Fortran builds each binary record in one array seen three ways
//! (`equivalence (a(1),ia(1),ha(1))`: 4-byte reals, 4-byte integers and
//! 8-byte words). The port keeps that array as bytes (`Mem`), so every
//! record has the same bytes, including the 4-byte rounding of the reals
//! and the blanks that the A edit puts after 6 characters.

#![allow(clippy::too_many_lines)]

use crate::endf::Ctx;
use crate::farr::RV;
use crate::ffmt;
use crate::fio::{fatal_runtime, Kind, RecW, Val};

/// Extra words past the Fortran array sizes.
const SLACK: usize = 20000;

/// `real(k4)::a(10000)`, `integer(k4)::ia(10000)`, `real(k8)::ha(5000)`,
/// equivalenced.
struct Mem(Vec<u8>);

impl Mem {
    fn a_set(&mut self, i: i32, v: f64) {
        let o = 4 * (i - 1) as usize;
        self.0[o..o + 4].copy_from_slice(&(v as f32).to_le_bytes());
    }
    fn ia_set(&mut self, i: i32, v: i32) {
        let o = 4 * (i - 1) as usize;
        self.0[o..o + 4].copy_from_slice(&v.to_le_bytes());
    }
    fn ia_get(&self, i: i32) -> i32 {
        let o = 4 * (i - 1) as usize;
        i32::from_le_bytes(self.0[o..o + 4].try_into().unwrap())
    }
    fn ha_set(&mut self, i: i32, v: &[u8; 8]) {
        let o = 8 * (i - 1) as usize;
        self.0[o..o + 8].copy_from_slice(v);
    }
    /// Write the first `nbytes` bytes as one unformatted record.
    fn write(&self, c: &mut Ctx, unit: i32, nbytes: usize) {
        let mut r = RecW::new();
        r.b.extend_from_slice(&self.0[..nbytes]);
        c.io.write_rec(unit, r);
    }
    /// Read `nbytes` bytes of the next record into the start of the array.
    fn read(&mut self, c: &mut Ctx, unit: i32, nbytes: usize) {
        let m = &mut self.0;
        c.io.read_rec_with(unit, |r| {
            if r.b.len() < nbytes {
                fatal_runtime("Fortran runtime error: I/O past end of record on unformatted file");
            }
            m[..nbytes].copy_from_slice(&r.b[..nbytes]);
        });
    }
}

/// `read(s,'(a6)') x` for an 8-byte variable: 6 characters, then 2 blanks.
fn a6(s: &[u8]) -> [u8; 8] {
    let mut h = [b' '; 8];
    for (j, b) in s.iter().take(6).enumerate() {
        h[j] = *b;
    }
    h
}

/// A list-directed character value of length n (blank padded).
fn chars(v: &Option<Val>, n: usize) -> Vec<u8> {
    let mut out = vec![b' '; n];
    if let Some(Val::Char(s)) = v {
        for (j, b) in s.as_bytes().iter().take(n).enumerate() {
            out[j] = *b;
        }
    }
    out
}

fn geti(v: &Option<Val>, d: i32) -> i32 {
    match v {
        Some(Val::Int(i)) => *i,
        Some(Val::Real(x)) => *x as i32,
        _ => d,
    }
}

fn getr(v: &Option<Val>, d: f64) -> f64 {
    match v {
        Some(Val::Real(x)) => *x,
        Some(Val::Int(i)) => *i as f64,
        _ => d,
    }
}

/// The saved variable of `resxsr` (`integer::nx=31` has the save
/// attribute, and ENDF-4 input changes it).
pub struct Resxsr {
    nx: i32,
}

impl Default for Resxsr {
    fn default() -> Self {
        Resxsr { nx: 31 }
    }
}

impl Resxsr {
    pub fn new() -> Self {
        Self::default()
    }

    /// Construct a RESXS resonance cross section file.
    pub fn resxsr(&mut self, c: &mut Ctx) {
        let (mut nb, mut nw) = (0i32, 0i32);
        let mult = 2i32;
        let nbuf = 5000i32;
        let nblok = 5000i32;
        let hfile = a6(b"resxsr");
        let hblank = [b' '; 8];
        let mut mem = Mem(vec![0u8; 40000 + 8 * SLACK]);
        let mut edat = RV::new(5000 + SLACK);
        let nsig = 31 + 2000;
        let mut sigs = vec![0.0f64; nsig];
        let mut sigi = vec![0.0f64; nsig];
        let mut sigl = vec![0.0f64; nsig];
        // stack(50,31), column major, kept flat
        let mut stack = vec![0.0f64; 50 * nsig + SLACK];
        let st = |i: i32, j: i32| ((i - 1) + 50 * (j - 1)) as usize;
        let mut temp = vec![0.0f64; 10 + 1 + 100];
        let mut hmat: Vec<Vec<u8>> = vec![vec![b' '; 8]; 21];
        let mut mat = [0i32; 21];
        let mut ninm = [0i32; 21];
        let mut ntemp = [0i32; 21];
        let mut locm = [0i32; 21];
        let mut awrm = [0.0f64; 21];

        //--assign storage pointers
        let iscr = 1i32;
        let iscrh = 1i32;
        let mut bufo = vec![0.0f64; nbuf as usize];
        let mut bufn = vec![0.0f64; nbuf as usize];

        //--assign i/o units
        let nscr = 10;
        let mut inew = 11;
        let mut iold = 12;
        c.openz(-nscr, 0);
        c.openz(-inew, 0);
        c.openz(-iold, 0);
        let mut irec = 0i32;

        //--read user input
        let v = c.io.input.read(&[Kind::Int]);
        let nout = geti(&v[0], 0);
        let v = c.io.input.read(&[Kind::Int, Kind::Int, Kind::Int, Kind::Real, Kind::Real, Kind::Real]);
        let nmat = geti(&v[0], 0);
        let maxt = geti(&v[1], 0);
        let nholl = geti(&v[2], 0);
        let efirst = getr(&v[3], 0.0);
        let elast = getr(&v[4], 0.0);
        let eps = getr(&v[5], 0.0);
        let v = c.io.input.read(&[Kind::Char(12), Kind::Int]);
        let uid = chars(&v[0], 12);
        let ivers = geti(&v[1], 0);
        let huse1 = a6(&uid[0..6]);
        let huse2 = a6(&uid[6..12]);
        for _ in 1..=nholl {
            // the descriptive text is read but not written to the file
            let _ = c.io.input.read(&[Kind::Char(72)]);
        }
        for i in 1..=nmat as usize {
            let v = c.io.input.read(&[Kind::Char(8), Kind::Int, Kind::Int]);
            if v[0].is_some() {
                hmat[i] = chars(&v[0], 8);
            }
            mat[i] = geti(&v[1], mat[i]);
            ninm[i] = geti(&v[2], ninm[i]);
        }

        //--loop over desired materials
        let mut je = 0i32;
        let mut jx;
        for im in 1..=nmat as usize {
            let matd = mat[im];
            locm[im] = irec;
            let nin = ninm[im];
            c.openz(nin, 0);
            c.tpidio(nin, 0, 0, &mut edat.all(), &mut nb, &mut nw);

            //--initialize loada/finda file
            for ix in 0..self.nx as usize {
                sigs[ix] = 0.0;
            }
            sigs[0] = efirst;
            c.loada(1, &sigs, self.nx, inew, &mut bufn, nbuf);
            sigs[0] = elast;
            c.loada(-2, &sigs, self.nx, inew, &mut bufn, nbuf);

            //--loop over temperatures
            let mut itemp = 0i32;
            jx = 0i32;
            loop {
                // label 110
                c.contio(nin, 0, 0, &mut edat.all(), &mut nb, &mut nw);
                if c.math != matd {
                    break;
                }
                if itemp == maxt {
                    break;
                }
                let awr = c.c2h;
                awrm[im] = awr;
                c.contio(nin, 0, 0, &mut edat.all(), &mut nb, &mut nw);
                if c.n1h != 0 {
                    c.iverf = 4;
                } else if c.n2h == 0 {
                    c.iverf = 5;
                } else {
                    c.iverf = 6;
                }
                c.skiprz(nin, -1);
                if c.iverf == 4 {
                    self.nx = c.n2h;
                }
                if c.iverf >= 5 {
                    c.contio(nin, 0, 0, &mut edat.all(), &mut nb, &mut nw);
                }
                if c.iverf >= 6 {
                    c.contio(nin, 0, 0, &mut edat.all(), &mut nb, &mut nw);
                }
                c.hdatio(nin, 0, 0, &mut edat.all(), &mut nb, &mut nw);
                itemp += 1;
                ntemp[im] = itemp;
                temp[itemp as usize] = edat[1];
                let hm = String::from_utf8_lossy(&hmat[im]).to_string();
                c.io.so(&format!(
                    " found: {}{}{}",
                    ffmt::a(&hm, 8),
                    ffmt::i(c.math, 5),
                    ffmt::e1p(temp[itemp as usize], 12, 4)
                ));

                //--loop over resonance reactions
                c.findf(matd, 3, 0, nin);
                loop {
                    // label 120
                    c.contio(nin, 0, 0, &mut edat.all(), &mut nb, &mut nw);
                    if c.mfh == 0 {
                        break;
                    }
                    if !(c.mth == 2 || c.mth == 18 || c.mth == 102) {
                        c.tosend(nin, 0, 0, &mut edat.all());
                        continue;
                    }
                    jx += 1;

                    //--loop through the energy grid for this reaction
                    //--and add the cross sections to the union set
                    let nx = self.nx;
                    let mut e = 0.0f64;
                    let mut enext = 0.0f64;
                    let mut idis = 0i32;
                    let mut s = 0.0f64;
                    c.gety1(e, &mut enext, &mut idis, &mut s, nin, &mut edat.all());
                    std::mem::swap(&mut iold, &mut inew);
                    let mut ie = 1i32;
                    c.finda(ie, &mut sigs, nx, iold, &mut bufo, nbuf);
                    e = sigs[0];
                    c.gety1(e, &mut enext, &mut idis, &mut s, nin, &mut edat.all());
                    sigs[jx as usize] = s;
                    je = 1;
                    c.loada(je, &sigs, nx, inew, &mut bufn, nbuf);
                    'l150: loop {
                        ie += 1;
                        for ix in 0..nx as usize {
                            sigl[ix] = sigs[ix];
                        }
                        c.finda(ie, &mut sigi, nx, iold, &mut bufo, nbuf);
                        let ei = sigi[0];
                        loop {
                            // label 160
                            e = ei;
                            if e > enext {
                                e = enext;
                            }
                            c.gety1(e, &mut enext, &mut idis, &mut s, nin, &mut edat.all());
                            if jx != 1 {
                                for ix in 2..=jx as usize {
                                    let mut y = sigs[ix - 1];
                                    c.terp1(sigl[0], sigl[ix - 1], sigi[0], sigi[ix - 1], e, &mut y, 2);
                                    sigs[ix - 1] = y;
                                }
                            }
                            sigs[0] = e;
                            sigs[jx as usize] = s;
                            je += 1;
                            if e > elast - elast / 200000.0 {
                                break 'l150;
                            }
                            c.loada(je, &sigs, nx, inew, &mut bufn, nbuf);
                            if e == ei {
                                continue 'l150;
                            }
                        }
                    }
                    // label 190
                    c.loada(-je, &sigs, nx, inew, &mut bufn, nbuf);
                    c.io.so(&format!("  mt={}    ne={}", ffmt::i(c.mth, 4), ffmt::i(je, 6)));
                    c.tosend(nin, 0, 0, &mut edat.all());
                }

                //--continue temperature loop for this material
                c.tomend(nin, 0, 0, &mut edat.all());
            }

            //--thin the cross sections for this material
            let nx = self.nx;
            let ne = je;
            std::mem::swap(&mut inew, &mut iold);
            let mut ie = 0i32;
            je = 0;
            let mut n = 0i32;
            loop {
                ie += 1;
                c.finda(ie, &mut sigs, nx, iold, &mut bufo, nbuf);
                if ie == ne {
                    break;
                }
                n += 1;
                for j in 1..=nx {
                    stack[st(n, j)] = sigs[(j - 1) as usize];
                }
                if n < 3 {
                    continue;
                }
                let lim = n - 1;
                let mut fail = false;
                'o: for i in 2..=lim {
                    for j in 2..=nx {
                        let mut test = 0.0f64;
                        c.terp1(
                            stack[st(1, 1)],
                            stack[st(1, j)],
                            stack[st(n, 1)],
                            stack[st(n, j)],
                            stack[st(i, 1)],
                            &mut test,
                            2,
                        );
                        if (test - stack[st(i, j)]).abs() > eps * stack[st(i, j)] {
                            fail = true;
                            break 'o;
                        }
                    }
                }
                if !fail {
                    continue;
                }
                for j in 1..=nx {
                    sigs[(j - 1) as usize] = stack[st(1, j)];
                }
                je += 1;
                c.loada(je, &sigs, nx, inew, &mut bufn, nbuf);
                for i in 1..=2 {
                    for j in 1..=nx {
                        stack[st(i, j)] = stack[st(lim + i - 1, j)];
                    }
                }
                n = 2;
            }
            je += 1;
            c.loada(-je, &sigs, nx, inew, &mut bufn, nbuf);
            c.io.so(&format!(" after thinning, ne={}", ffmt::i(je, 6)));

            //--write the resxs material control record
            let n = ntemp[im];
            let nwds = mult + n + 3;
            irec += 1;
            mem.ha_set(iscrh, &a6(&hmat[im]));
            mem.a_set(iscr + mult, awrm[im]);
            for i in 1..=n {
                mem.a_set(iscr + mult + i, temp[i as usize]);
            }
            mem.ia_set(iscr + mult + n + 1, jx / n);
            mem.ia_set(iscr + mult + n + 2, je);
            mem.write(c, nscr, 4 * nwds as usize);

            //--write the thinned cross-section blocks
            let nn = 1 + jx;
            let mut nwds = 0i32;
            for i in 1..=je {
                c.finda(i, &mut sigs, nx, inew, &mut bufn, nbuf);
                if nwds + nn > nblok {
                    irec += 1;
                    mem.write(c, nscr, 4 * nwds as usize);
                    nwds = 0;
                }
                for j in 1..=nn {
                    mem.a_set(iscr - 1 + nwds + j, sigs[(j - 1) as usize]);
                }
                nwds += nn;
            }
            if nwds != 0 {
                irec += 1;
                mem.write(c, nscr, 4 * nwds as usize);
            }
        }

        //--write the output in resxs format
        //--start with the file identification record
        c.openz(-nout, 1);
        let nwds = 3 * mult + 1;
        mem.ha_set(iscrh, &hfile);
        mem.ha_set(iscrh + 1, &huse1);
        mem.ha_set(iscrh + 2, &huse2);
        mem.ia_set(iscr + 3 * mult, ivers);
        c.repoz(-nout);
        mem.write(c, nout, 4 * nwds as usize);

        //--file control record
        mem.a_set(iscr, efirst);
        mem.a_set(iscr + 1, elast);
        mem.ia_set(iscr + 2, nholl);
        mem.ia_set(iscr + 3, nmat);
        mem.ia_set(iscr + 4, nblok);
        mem.write(c, nout, 20);

        //--set hollerith identification
        for i in 1..=mult * nholl {
            mem.ha_set(i, &hblank);
        }
        mem.write(c, nout, 8 * nholl.max(0) as usize);

        //--file data
        for i in 1..=nmat {
            mem.ha_set(iscrh - 1 + i, &a6(&hmat[i as usize]));
            mem.ia_set(iscr - 1 + nmat * mult + i, ntemp[i as usize]);
            mem.ia_set(iscr - 1 + nmat * (mult + 1) + i, locm[i as usize]);
        }
        mem.write(c, nout, 16 * nmat.max(0) as usize);

        //--copy the material control and cross-section blocks
        //--from the scratch file to nout
        c.repoz(-nscr);
        for i in 1..=nmat as usize {
            let nt = ntemp[i];
            let nwds = mult + nt + 3;
            mem.read(c, nscr, 4 * nwds as usize);
            mem.write(c, nout, 4 * nwds as usize);
            let nreac = mem.ia_get(iscr + mult + nt + 1);
            let nener = mem.ia_get(iscr + mult + nt + 2);
            let nn = 1 + nreac * nt;
            let mut nw = nener * nn;
            let mut nwds = (nblok / nn) * nn;
            let mut nb = nw / nwds;
            if nw > nb * nwds {
                nb += 1;
            }
            for _ in 1..=nb {
                if nw < nwds {
                    nwds = nw;
                }
                mem.read(c, nscr, 4 * nwds as usize);
                mem.write(c, nout, 4 * nwds as usize);
                nw -= nwds;
            }
        }
    }
}
