//! POWR lib=2: the EPRI-CELL thermal library (LIBRAR input).

use super::{a4s, geti, getr, getc, pe, rows, rowsi, rowsj, Powr, SLACK};
use crate::endf::Ctx;
use crate::farr::RV;
use crate::ffmt;
use crate::fio::Kind;
use crate::fmath::nint;

const ZERO: f64 = 0.0;

impl Powr {
    /// Read data for adding nuclides to thermal library from a gout tape.
    pub(super) fn therm(&mut self, c: &mut Ctx, iprint: i32, igprnt: &mut i32) {
        let watt = 1.0e13f64;
        let ng = self.ngendf;
        let nout = self.nout;
        let (mut nb, mut nw) = (0i32, 0i32);
        self.ngnd = 35;
        let ngnd = self.ngnd;
        c.io.so("\n ***epri-cell thermal library***/");

        //--allocate storage.
        let mut s0 = RV::new(ngnd as usize + SLACK);
        let mut s1 = RV::new(ngnd as usize + SLACK);
        let mut xa = RV::new(ngnd as usize + SLACK);
        let mut xf = RV::new(ngnd as usize + SLACK);
        let ngndsq = ngnd * ngnd;
        let mut b = RV::new(ngndsq as usize + SLACK);
        let mut p = RV::new(ngndsq as usize + SLACK);
        let ngnd1 = ngnd + 1;
        let mut grp = RV::new(ngnd1 as usize + SLACK);
        let il = 1;
        let lz = 6;
        let mut igi = [0i32; 36 + 1];
        let mut ngi = [0i32; 36 + 1];
        let mut z = [0.0f64; 11];
        let mut aname: Vec<String> = vec!["    ".to_string(); 3];
        let mut id4 = 0;
        let mut id5 = 0;

        //--read and write user input.
        let mut n = 0;
        let mut matd = 0i32;
        let (mut idtemp, mut itrc, mut mti, mut mtc) = (0i32, 0i32, 0i32, 0i32);
        let (mut i318, mut mf3mt2, mut iflx0, mut iflx1, mut ifound) = (0, 0, 0, 0, 0);
        let (mut nl, mut nz) = (0i32, 0i32);
        let (mut idm, mut idma, mut idb, mut ida, mut idf, mut id1, mut id2, mut id3) = (0, 0, 0, 0, 0, 0, 0, 0);
        let mut nlm = 0;
        let (mut nuclid, mut nsigz, mut ntw) = (0i32, 0i32, 0i32);
        let mut awr = 0.0f64;
        let (mut ig, mut igc) = (0i32, 0i32);
        let mut nbm = 0i32;
        let mut lbl = 100;
        loop {
            match lbl {
                100 => {
                    idtemp = 300;
                    let v = c.io.input.read(&[Kind::Int, Kind::Int, Kind::Char(12)]);
                    matd = geti(&v[0], matd);
                    idtemp = geti(&v[1], idtemp);
                    let idname = getc(&v[2], 12, &" ".repeat(12));
                    if matd == 0 {
                        lbl = 300;
                        continue;
                    }
                    aname = a4s(&idname, 3);
                    itrc = 0;
                    mti = 0;
                    mtc = 0;
                    let v = c.io.input.read(&[Kind::Int, Kind::Int, Kind::Int]);
                    itrc = geti(&v[0], itrc);
                    mti = geti(&v[1], mti);
                    mtc = geti(&v[2], mtc);
                    for i in 1..=8 {
                        z[i] = 0.0;
                    }
                    let v = c.io.input.read(&[Kind::Real; 8]);
                    for i in 1..=8 {
                        z[i] = getr(&v[i - 1], z[i]);
                    }
                    self.valu[2] = z[1];
                    self.valu[3] = z[2];
                    self.valu[4] = z[3];
                    self.valu[5] = z[4];
                    self.valu[6] = z[5] * watt;
                    self.valu[7] = z[6] * watt;
                    self.valu[8] = 0.0;
                    self.valu[9] = z[7];
                    self.valu[10] = z[8];
                    n += 1;
                    c.io.so(&format!(
                        "\n\n mat to be processed .................. {}\n \
                         temperature id ....................... {}\n \
                         isotope name ......................... {}{}{}",
                        ffmt::i(matd, 10),
                        ffmt::i(idtemp, 10),
                        aname[0],
                        aname[1],
                        &aname[2][..2]
                    ));
                    c.io.so(&format!(
                        " transport corr. option (0 no, 1 yes) . {}\n \
                         mt inelastic ......................... {}\n \
                         mt elastic ........................... {}",
                        ffmt::i(itrc, 10),
                        ffmt::i(mti, 10),
                        ffmt::i(mtc, 10)
                    ));
                    c.io.so(&format!(
                        " xi ................................... {}\n \
                         alpha ................................ {}\n \
                         mubar ................................ {}\n \
                         nu ................................... {}\n \
                         kappa fission ........................ {}\n \
                         kappa capture ........................ {}\n \
                         lambda ............................... {}\n \
                         sigma s .............................. {}",
                        pe(z[1], 10, 3),
                        pe(z[2], 10, 3),
                        pe(z[3], 10, 3),
                        pe(z[4], 10, 3),
                        pe(z[5], 10, 3),
                        pe(z[6], 10, 3),
                        pe(z[7], 10, 3),
                        pe(z[8], 10, 3)
                    ));

                    //--read data for desired nuclide.
                    i318 = 0;
                    mf3mt2 = 0;
                    for i in 1..=ngndsq {
                        b[i] = 0.0;
                        p[i] = 0.0;
                    }
                    iflx0 = 0;
                    iflx1 = 0;
                    c.repoz(ng);
                    c.contio(ng, 0, 0, &mut self.scr.at(1), &mut nb, &mut nw);
                    ifound = 0;
                    lbl = 115;
                }
                115 => {
                    c.contio(ng, 0, 0, &mut self.scr.at(1), &mut nb, &mut nw);
                    if c.math < 0 {
                        lbl = 220;
                        continue;
                    }
                    if c.math == 0 || c.mfh == 0 || c.mth == 0 {
                        lbl = 115;
                        continue;
                    }
                    if c.math == matd {
                        lbl = 120;
                        continue;
                    }
                    if ifound > 0 {
                        lbl = 220;
                        continue;
                    }
                    c.tomend(ng, 0, 0, &mut self.scr.at(1));
                    lbl = 115;
                }
                120 => {
                    nl = c.l1h;
                    nz = c.l2h;
                    if c.mfh == 3 {
                        lbl = 152;
                        continue;
                    }
                    if c.mfh == 6 && c.mth == mti {
                        lbl = 152;
                        continue;
                    }
                    if c.mfh == 6 && c.mth == mtc && mtc != 0 {
                        lbl = 152;
                        continue;
                    }
                    if c.mfh == 1 {
                        lbl = 130;
                        continue;
                    }
                    c.tosend(ng, 0, 0, &mut self.scr.at(1));
                    lbl = 115;
                }
                130 => {
                    let iza = nint(c.c1h);
                    nuclid = iza * 10;
                    awr = c.c2h;
                    nsigz = c.l2h;
                    ntw = c.n2h;
                    c.listio(ng, 0, 0, &mut self.scr.at(1), &mut nb, &mut nw);
                    let itemp = nint(c.c1h);
                    if idtemp == itemp {
                        lbl = 132;
                        continue;
                    }
                    if ifound > 0 {
                        lbl = 220;
                        continue;
                    }
                    c.tomend(ng, 0, 0, &mut self.scr.at(1));
                    lbl = 115;
                }
                132 => {
                    ifound = 1;
                    idm = nuclid;
                    idma = idtemp;
                    idb = 0;
                    ida = 0;
                    idf = 0;
                    id1 = 0;
                    id2 = 0;
                    id3 = 0;
                    id4 = 0;
                    id5 = 0;
                    nlm = 0;
                    let mut ia = 1;
                    while nb != 0 {
                        ia += nw;
                        c.moreio(ng, 0, 0, &mut self.scr.at(ia), &mut nb, &mut nw);
                        if ia + nw - 1 > self.nwscr {
                            c.error("therm", "storage exceeded.", " ");
                        }
                    }
                    self.valu[1] = awr;
                    c.io.so(&format!(" atomic weight ........................ {}", pe(self.valu[1], 10, 3)));
                    self.ngn = c.l1h;
                    if self.ngn != ngnd && self.ngn != self.ngnf {
                        c.error("therm", "ngn is not equal to ngnd.", " ");
                    }
                    if *igprnt <= 0 {
                        let ia = 1 + lz + ntw + nsigz;
                        let lim = self.ngn + 1;
                        if self.ngn == ngnd {
                            for i in 1..=lim {
                                grp[i] = self.scr[i - 1 + ia];
                            }
                        } else {
                            for i in 1..=lim {
                                let ig = self.icgrp(i);
                                if ig != 0 {
                                    grp[ig] = self.scr[ia + i - 1];
                                }
                            }
                        }
                        c.io.so(&format!(
                            "\n24h neutron group structure\n -----------------------\n{}",
                            rowsj(&grp.0[..ngnd1 as usize], 6, |x| pe(x, 12, 4))
                        ));
                        *igprnt = 1;
                    }
                    for i in 1..=ngnd {
                        self.cflux[i] = 0.0;
                        self.cflux[i + ngnd] = 0.0;
                        s0[i] = 0.0;
                        s1[i] = 0.0;
                        xa[i] = 0.0;
                        xf[i] = 0.0;
                    }
                    lbl = 115;
                }
                152 => {
                    c.listio(ng, 0, 0, &mut self.scr.at(1), &mut nb, &mut nw);
                    let mut ia = 1;
                    while nb != 0 {
                        ia += nw;
                        c.moreio(ng, 0, 0, &mut self.scr.at(ia), &mut nb, &mut nw);
                        if ia + nw - 1 > self.nwscr {
                            c.error("therm", "storage exceeded.", " ");
                        }
                    }
                    ig = c.n2h;
                    igc = self.icgrp(ig);
                    if !(mf3mt2 > 0 || iflx0 == 0 || (c.mfh == 3 && c.mth == 1)) {
                        mf3mt2 = 1;
                        for i in 1..=ngnd {
                            if self.cflux[i] != ZERO {
                                self.cflux[i] = 1.0 / self.cflux[i];
                            }
                            if self.cflux[ngnd + i] != ZERO {
                                self.cflux[ngnd + i] = 1.0 / self.cflux[ngnd + i];
                            }
                        }
                    }
                    // 160
                    if c.mfh == 6 {
                        lbl = 185;
                        continue;
                    }
                    if igc == 0 || igc > ngnd {
                        lbl = 166;
                        continue;
                    }
                    if c.mth == 1 {
                        // coarse group flux
                        // p0
                        self.cflux[igc] = self.cflux[igc] + self.scr[1 + lz];
                        iflx0 = 1;
                        if nl > 1 {
                            // p1
                            self.cflux[ngnd + igc] = self.cflux[ngnd + igc] + self.scr[lz + 2];
                            iflx1 = 1;
                        }
                        lbl = 166;
                        continue;
                    }
                    if c.mth > 1 && iflx0 == 0 {
                        c.error("therm", "p0 flux not found.", " ");
                    }
                    if c.mth == 18 || (c.mth == 19 && i318 == 0) {
                        // fission
                        idf = 1;
                        if c.mth == 18 {
                            i318 = 1;
                        }
                        xf[igc] = xf[igc] + self.scr[1 + lz + nl * nz] * self.scr[1 + lz] * self.cflux[igc];
                        lbl = 166;
                        continue;
                    }
                    if c.mth >= 102 && c.mth <= 150 {
                        // capture
                        xa[igc] = xa[igc] + self.scr[1 + lz + nl * nz] * self.scr[1 + lz] * self.cflux[igc];
                        lbl = 166;
                        continue;
                    }
                    c.tosend(ng, 0, 0, &mut self.scr.at(1));
                    lbl = 115;
                }
                166 => {
                    lbl = if ig < self.ngn { 152 } else { 115 };
                }
                185 => {
                    // scattering matrix
                    if ida <= 0 {
                        c.io.so(
                            "\n scattering cross sections and mubar\n        sig0        sig1       mubar\n -----------------------------------",
                        );
                        ida = 1;
                    }
                    let ig2lo = c.l2h;
                    let ng2 = c.l1h;
                    if nl > nlm {
                        nlm = nl;
                    }
                    ig = c.n2h;
                    igc = self.icgrp(ig);
                    if igc != 0 && igc <= ngnd {
                        for i in 2..=ng2 {
                            let loca = 1 + (il - 1) + lz + nl * nz * (i - 1);
                            let ig2 = ig2lo + i - 2;
                            let ig2c = self.icgrp(ig2);
                            if ig2c != 0 && ig2c <= ngnd {
                                let loc = ig2c + ngnd * (igc - 1);
                                p[loc] = p[loc] + self.scr[loca] * self.scr[1 + lz] * self.cflux[igc];
                                s0[igc] = s0[igc] + self.scr[loca] * self.scr[1 + lz] * self.cflux[igc];
                                if nl != 1 {
                                    let mut ndexc = ngnd + igc;
                                    if iflx1 == 0 {
                                        ndexc = igc;
                                    }
                                    let mut ndexf = 1 + lz + 1;
                                    if iflx1 == 0 {
                                        ndexf = 1 + lz;
                                    }
                                    s1[igc] = s1[igc] + self.scr[loca + 1] * self.scr[ndexf] * self.cflux[ndexc];
                                }
                            }
                        }
                    }
                    lbl = if ig < self.ngn { 152 } else { 115 };
                }
                //--prepare final matrix
                220 => {
                    if ifound == 0 {
                        let strng = format!("mat={} temp {} not found", ffmt::i(matd, 4), ffmt::i(idtemp, 2));
                        c.error("therm", &strng, " ");
                    }
                    if ida != 0 {
                        for ig in 1..=ngnd {
                            let mut amubar = 0.0;
                            if s0[ig] == ZERO {
                                amubar = 0.0;
                            }
                            if s0[ig] != ZERO {
                                amubar = s1[ig] / s0[ig];
                            }
                            c.io.so(&format!("{}{}{}", pe(s0[ig], 12, 4), pe(s1[ig], 12, 4), pe(amubar, 12, 4)));
                        }
                        // convert to xsec per unit energy.
                        for ig in 1..=ngnd {
                            for ig2 in 1..=ngnd {
                                p[ig2 + ngnd * (ig - 1)] = p[ig2 + ngnd * (ig - 1)] * (grp[1 + ig] - grp[ig])
                                    / (grp[1 + ig2] - grp[ig2]);
                            }
                        }
                        // pack scattering matrix
                        packa(&p, ngnd, &mut b, &mut nbm, &mut igi, &mut ngi);
                        id3 = nbm;
                    }

                    //--compute absorption and transport cross sections
                    for i in 1..=ngnd {
                        if idf != 0 {
                            xa[i] = xa[i] + xf[i];
                        }
                        if ida != 0 {
                            s1[i] = xa[i] + s0[i] - s1[i];
                        }
                    }
                    if self.valu[10] == 0.0 {
                        self.valu[10] = s0[ngnd];
                    }

                    //--write out data for this nuclide.
                    if itrc == 0 {
                        self.valu[4] = 1.0 - self.valu[4];
                    }
                    if itrc == 1 {
                        self.valu[4] = -1.0;
                    }
                    let mut s = format!("{}{}", ffmt::i(idm, 10), ffmt::i(idma, 10));
                    for x in [idb, ida, idf, id1, id2, id3, id4, id5] {
                        s.push_str(&ffmt::i(x, 5));
                    }
                    c.io.write_line(nout, &s);
                    let mut s = format!("{}{}{}", aname[0], aname[1], &aname[2][..2]);
                    for i in 1..=6 {
                        s.push_str(&ffmt::f(self.valu[i], 10, 4));
                    }
                    c.io.write_line(nout, &s);
                    let mut s = String::new();
                    for i in 7..=13 {
                        s.push_str(&ffmt::f(self.valu[i], 10, 4));
                    }
                    c.io.write_line(nout, &s);
                    let mut s = String::new();
                    for i in 14..=20 {
                        s.push_str(&ffmt::f(self.valu[i], 10, 4));
                    }
                    c.io.write_line(nout, &s);
                    let g = ngnd as usize;
                    for s in rows(&xa.0[..g], 6, |x| pe(x, 12, 5)) {
                        c.io.write_line(nout, &s);
                    }
                    if idf == 1 {
                        for s in rows(&xf.0[..g], 6, |x| pe(x, 12, 5)) {
                            c.io.write_line(nout, &s);
                        }
                    }
                    if itrc == 1 {
                        for s in rows(&s1.0[..g], 6, |x| pe(x, 12, 5)) {
                            c.io.write_line(nout, &s);
                        }
                    }
                    if ida == 1 {
                        for s in rowsi(&igi[1..=g], 12, 6) {
                            c.io.write_line(nout, &s);
                        }
                        for s in rowsi(&ngi[1..=g], 12, 6) {
                            c.io.write_line(nout, &s);
                        }
                        for s in rows(&b.0[..nbm.max(0) as usize], 6, |x| pe(x, 12, 5)) {
                            c.io.write_line(nout, &s);
                        }
                    }
                    c.io.so(&format!(
                        "\n storage pointers\n ----------------\n      idm        {}\n      idma       {}\n      \
                         idb        {}\n      ida        {}\n      idf        {}\n      id1        {}\n      \
                         id2        {}\n      id3        {}",
                        ffmt::i(idm, 9),
                        ffmt::i(idma, 9),
                        ffmt::i(idb, 9),
                        ffmt::i(ida, 9),
                        ffmt::i(idf, 9),
                        ffmt::i(id1, 9),
                        ffmt::i(id2, 9),
                        ffmt::i(id3, 9)
                    ));
                    if iprint != 0 {
                        c.io.so("\n absorption cross section\n ------------------------");
                        c.io.so(&rowsj(&xa.0[..g], 6, |x| pe(x, 12, 4)));
                        if idf != 0 {
                            c.io.so("\n fission cross section\n ---------------------");
                            c.io.so(&rowsj(&xf.0[..g], 6, |x| pe(x, 12, 4)));
                        }
                        if ida != 0 {
                            if itrc == 1 {
                                c.io.so("\n\n transport cross section\n -----------------------");
                            }
                            if itrc == 1 {
                                c.io.so(&rowsj(&s1.0[..g], 6, |x| pe(x, 12, 4)));
                            }
                            c.io.so("\n\n scattering matrix\n -----------------");
                            tprnt(c, &b, &igi, &ngi, ngnd);
                        }
                    }
                    lbl = 100;
                }
                //--therm is finished.
                _ => {
                    c.repoz(nout);
                    let _ = (n, nlm, id4, id5);
                    return;
                }
            }
        }
    }
}

/// Matrix printing routine.
fn tprnt(c: &mut Ctx, a: &RV, igi: &[i32], ngi: &[i32], nix: i32) {
    let ncol = 5;
    let ncol1 = ncol - 1;
    let mut h = String::from("\n46h final initl    cross section vs initial group\n group group    +0");
    for i in 1..=ncol1 {
        h.push_str(&format!("{}+{}", " ".repeat(10), ffmt::i(i, 1)));
    }
    h.push_str(&format!("\n {}--------\n", "----------".repeat(7)));
    c.io.so(&h);
    let mut k = 0;
    for i in 1..=nix {
        let mut jlo = igi[i as usize];
        let mut ng = ngi[i as usize];
        let mut lim = ng;
        if lim > ncol {
            lim = ncol;
        }
        let mut s = format!("{}{}  ", ffmt::i(i, 5), ffmt::i(jlo, 6));
        for j in 1..=lim {
            s.push_str(&pe(a[k + j], 12, 4));
        }
        c.io.so(&s);
        loop {
            // 110
            k += lim;
            if lim == ng {
                break;
            }
            jlo += lim;
            ng -= lim;
            lim = ng;
            if lim > ncol {
                lim = ncol;
            }
            let mut s = format!("     {}  ", ffmt::i(jlo, 6));
            for j in 1..=lim {
                s.push_str(&pe(a[k + j], 12, 4));
            }
            c.io.so(&s);
        }
    }
}

/// Pack a scattering matrix into banded form by removing unneeded
/// zeroes.
fn packa(tpa: &RV, igroup: i32, a: &mut RV, k: &mut i32, igi: &mut [i32], ng: &mut [i32]) {
    let small = 1.0e-6f64;
    let t = |i: i32, j: i32| tpa[i + igroup * (j - 1)];
    *k = 0;
    for i in 1..=igroup {
        let iu = i as usize;
        igi[iu] = 0;
        ng[iu] = 0;
        let mut idone = 0;
        let mut j = 0;
        let mut jj = 0;
        while j < igroup && idone == 0 {
            j += 1;
            jj = igroup - j + 1;
            if t(i, jj).abs() >= small {
                idone = 1;
            }
        }
        let jhi = jj;
        for j in 1..=jhi {
            if t(i, j).abs() >= small || igi[iu] != 0 {
                if t(i, j).abs() >= small && igi[iu] == 0 {
                    igi[iu] = j;
                }
                *k += 1;
                a[*k] = t(i, j);
            }
        }
        ng[iu] = jhi - igi[iu] + 1;
    }
}
