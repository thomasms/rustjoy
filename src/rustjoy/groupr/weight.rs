//! GROUPR group structures, weight functions and the flux calculator
//! (gengpn, gengpg, genwtf, getwtf, genflx, getfwt).

use super::*;

/// Read the next list-directed values for n reals.
fn read_reals(c: &mut Ctx, n: i32, a: &mut RV) {
    let v = c.io.input.read(&vec![Kind::Real; n.max(0) as usize]);
    for i in 1..=n {
        a[i] = getr(&v[(i - 1) as usize], a[i]);
    }
}

/// Generate the requested neutron group structure, or read it in.
pub fn gengpn(c: &mut Ctx, ign: i32, ngn: &mut i32, egn: &mut RV) {
    let ezero = 1.0e7f64;
    let tenth = 0.10e0f64;
    let eighth = 0.125e0f64;
    let quart = 0.25e0f64;
    let bgam2 = 27.631021e0f64;
    let tgam2 = -0.53062825e0f64;
    let u187a = -15.5e0f64;
    let u187b = -14.125e0f64;
    let u187c = -5.875e0f64;
    let e187d = 2.6058e4f64;
    let e187e = 6.868e0f64;
    let sanda = 1.0e-4f64;
    let sandb = 1.0e-6f64;
    let sandc = 2.8e-4f64;
    let sandd = 1.0e6f64;
    let sande = 1.0e5f64;
    let uu80 = 0.6931472e0f64;
    let e175 = 1.284e7f64;

    let mut lflag = 0;
    let ngp: i32;
    let copy = |egn: &mut RV, t: &[f64], n: usize| {
        for ig in 1..=n {
            egn[ig as i32] = t[ig - 1];
        }
    };
    if ign.abs() == 1 {
        let v = c.io.input.read(&[Kind::Int]);
        *ngn = geti(&v[0], *ngn);
        ngp = *ngn + 1;
        *egn = alloc(ngp);
        read_reals(c, ngp, egn);
        for ig in 1..=*ngn {
            if egn[ig] >= egn[ig + 1] {
                c.error("gengpn", "read-in group structure is out of order.", " ");
            }
        }
    } else if ign == 2 {
        *ngn = 240;
        ngp = *ngn + 1;
        *egn = alloc(ngp);
        copy(egn, &N_GL2, ngp as usize);
        lflag = 1;
    } else if ign == 3 {
        *ngn = 30;
        ngp = *ngn + 1;
        *egn = alloc(ngp);
        copy(egn, &N_EG3, ngp as usize);
    } else if ign == 4 {
        *ngn = 27;
        ngp = *ngn + 1;
        *egn = alloc(ngp);
        copy(egn, &N_GL4, ngp as usize);
        lflag = 1;
    } else if ign == 5 {
        *ngn = 50;
        ngp = *ngn + 1;
        *egn = alloc(ngp);
        copy(egn, &N_GL5, ngp as usize);
        lflag = 1;
    } else if ign == 6 {
        *ngn = 68;
        let mut u = -quart;
        let du = quart;
        ngp = *ngn + 1;
        *egn = alloc(ngp);
        for ig in 1..=ngp {
            u = u + du;
            egn[70 - ig] = u;
        }
        lflag = 1;
    } else if ign == 7 {
        *ngn = 100;
        ngp = *ngn + 1;
        *egn = alloc(ngp);
        let mut u = -(4.0 * tenth);
        let mut du = tenth;
        for ig in 1..=99 {
            u = u + du;
            egn[101 - ig] = u;
            if ig == 49 {
                du = quart;
            }
        }
        egn[1] = bgam2;
        egn[101] = tgam2;
        lflag = 1;
    } else if ign == 8 {
        *ngn = 35;
        ngp = *ngn + 1;
        *egn = alloc(ngp);
        copy(egn, &N_EG6, ngp as usize);
    } else if ign == 9 {
        *ngn = 69;
        ngp = *ngn + 1;
        *egn = alloc(ngp);
        copy(egn, &N_EG9, ngp as usize);
    } else if ign == 10 {
        *ngn = 187;
        ngp = *ngn + 1;
        *egn = alloc(ngp);
        for ig in 1..=48 {
            egn[ig] = N_EG10A[(ig - 1) as usize];
        }
        let mut u = u187a;
        for ig in 49..=59 {
            egn[ig] = ezero * exp(u);
            u = u + eighth;
        }
        egn[60] = e187e;
        u = u187b;
        for ig in 61..=126 {
            egn[ig] = ezero * exp(u);
            u = u + eighth;
        }
        egn[127] = e187d;
        u = u187c;
        for ig in 128..=175 {
            egn[ig] = ezero * exp(u);
            u = u + eighth;
        }
        for ig in 176..=188 {
            egn[ig] = N_EG10B[(ig - 176) as usize];
        }
    } else if ign == 11 {
        *ngn = 70;
        ngp = *ngn + 1;
        *egn = alloc(ngp);
        copy(egn, &N_EG11, ngp as usize);
    } else if ign == 12 || ign == 15 {
        *ngn = 620;
        if ign == 15 {
            *ngn = 640;
        }
        ngp = *ngn + 1;
        *egn = alloc(ngp);
        egn[1] = sanda;
        for ig in 1..=8 {
            let delta = N_DELTL[ig - 1] * sandb;
            let n1 = N_NDELTA[ig - 1];
            let n2 = N_NDELTA[ig] - 1;
            for n in n1..=n2 {
                egn[n] = egn[n - 1] + delta;
            }
        }
        egn[21] = sandc;
        for ig in 46..=450 {
            egn[ig] = egn[ig - 45] * 10.0;
        }
        egn[451] = sandd;
        for ig in 452..=ngp {
            egn[ig] = egn[ig - 1] + sande;
        }
    } else if ign == 13 {
        *ngn = 80;
        let mut u = uu80;
        ngp = *ngn + 1;
        *egn = alloc(ngp);
        for ig in 1..=ngp {
            egn[82 - ig] = ezero * exp(u);
            if ig <= *ngn {
                u = u - N_U80[(ig - 1) as usize];
            }
        }
        egn[81] = 2.0 * ezero;
    } else if ign == 14 {
        *ngn = 100;
        ngp = *ngn + 1;
        *egn = alloc(ngp);
        egn[101] = -(4.0 * tenth);
        let mut ic = 0usize;
        for ig in 2..=101 {
            if ig == N_IG14[ic] {
                ic += 1;
            }
            egn[102 - ig] = egn[103 - ig] + N_GL14[ic - 1];
        }
        lflag = 1;
    } else if ign == 16 || ign == 17 {
        if ign == 16 {
            *ngn = 174;
        }
        if ign == 17 {
            *ngn = 175;
        }
        ngp = *ngn + 1;
        *egn = alloc(ngp.max(176));
        for ig in 1..=84 {
            egn[ig] = N_EG15A[(ig - 1) as usize];
        }
        for ig in 85..=175 {
            egn[ig] = N_EG15B[(ig - 85) as usize];
        }
        if ign != 16 {
            egn[166] = e175;
            for ig in 167..=176 {
                egn[ig] = N_EG15B[(ig - 86) as usize];
            }
        }
    } else if ign == 18 {
        *ngn = 172;
        ngp = *ngn + 1;
        *egn = alloc(ngp);
        for ig in 1..=ngp {
            egn[ig] = N_EG18[(174 - ig - 1) as usize];
        }
    } else if ign == 19 {
        *ngn = 33;
        ngp = *ngn + 1;
        *egn = alloc(ngp);
        copy(egn, &N_EG19, 34);
    } else if ign == 20 {
        *ngn = 1968;
        ngp = *ngn + 1;
        *egn = alloc(ngp);
        let parts: [&[f64]; 5] = [&N_EG20A, &N_EG20B, &N_EG20C, &N_EG20D, &N_EG20E];
        for (k, p) in parts.iter().enumerate() {
            for ig in 1..=392 {
                egn[ig + 392 * k as i32] = p[(ig - 1) as usize];
            }
        }
        for ig in 1..=9 {
            egn[ig + 1960] = N_EG20F[(ig - 1) as usize];
        }
    } else if ign == 21 {
        *ngn = 315;
        ngp = *ngn + 1;
        *egn = alloc(ngp);
        copy(egn, &N_EG21, 316);
    } else if ign == 22 {
        *ngn = 172;
        ngp = *ngn + 1;
        *egn = alloc(ngp);
        copy(egn, &N_EG22, 173);
    } else if ign == 23 {
        *ngn = 175;
        ngp = *ngn + 1;
        *egn = alloc(ngp);
        copy(egn, &N_EG23, 176);
    } else if ign == 24 {
        *ngn = 281;
        ngp = *ngn + 1;
        *egn = alloc(ngp);
        copy(egn, &N_EG24, 282);
    } else if ign == 25 {
        *ngn = 295;
        ngp = *ngn + 1;
        *egn = alloc(ngp);
        copy(egn, &N_EG25, 296);
    } else if ign == 26 {
        *ngn = 361;
        ngp = *ngn + 1;
        *egn = alloc(ngp);
        copy(egn, &N_EG26, 362);
    } else if ign == 27 {
        *ngn = 315;
        ngp = *ngn + 1;
        *egn = alloc(ngp);
        copy(egn, &N_EG27, 316);
    } else if ign == 28 {
        *ngn = 89;
        ngp = *ngn + 1;
        *egn = alloc(ngp);
        copy(egn, &N_EG28, 90);
    } else if ign == 29 {
        *ngn = 660;
        ngp = *ngn + 1;
        *egn = alloc(ngp);
        copy(egn, &N_EG29, ngp as usize);
    } else if ign == 30 {
        *ngn = 1025;
        ngp = *ngn + 1;
        *egn = alloc(ngp);
        copy(egn, &N_EG30, ngp as usize);
    } else if ign == 31 {
        *ngn = 1067;
        ngp = *ngn + 1;
        *egn = alloc(ngp);
        copy(egn, &N_EG31, ngp as usize);
    } else if ign == 32 {
        *ngn = 1102;
        ngp = *ngn + 1;
        *egn = alloc(ngp);
        copy(egn, &N_EG32, ngp as usize);
    } else if ign == 33 {
        *ngn = 142;
        ngp = *ngn + 1;
        *egn = alloc(ngp);
        copy(egn, &N_EG33, ngp as usize);
    } else if ign == 34 {
        *ngn = 618;
        ngp = *ngn + 1;
        *egn = alloc(ngp);
        copy(egn, &N_EG618, ngp as usize);
    } else if ign == 35 {
        *ngn = 99;
        ngp = *ngn + 1;
        *egn = alloc(ngp);
        copy(egn, &N_EG35, ngp as usize);
    } else if ign == 36 {
        *ngn = 1962;
        ngp = *ngn + 1;
        *egn = alloc(ngp);
        for ig in 1..=392 {
            egn[ig] = N_EG20A[(ig - 1) as usize];
        }
        for ig in 1..=392 {
            egn[ig + 392] = N_EG20B[(ig - 1) as usize];
        }
        for ig in 785..=ngp {
            egn[ig] = egn[ig - 1] * exp(1.0f64 / 120.0f64);
        }
    } else {
        c.error("gengpn", "illegal group structure.", " ");
    }

    //--convert lethargy grid to energies
    if lflag == 1 {
        let ngp = *ngn + 1;
        for ig in 1..=ngp {
            egn[ig] = sigfig(ezero * exp(-egn[ig]), 7, 0);
        }
    }

    //--display group structure - except for ign=-1
    let names: [&str; 37] = [
        "",
        "read in",
        "csewg 240 group",
        "lanl 30 group",
        "anl 27 group",
        "rrd 50 group",
        "gam-i 68 group",
        "gam-ii 100 group",
        "laser-thermos 35 group",
        "epri-cpm 69 group",
        "lanl 187-group",
        "lanl 70-group",
        "sand-ii 620 group",
        "lanl 80-group",
        "eurlib 100-group",
        "sand-iia 640-group",
        "vitamin-e 174-group",
        "vitamin-j 175-group",
        "xmas 172-group",
        "ecco  33-group",
        "ecco  1968-group",
        "tripoli 315-group",
        "xmas lwpc 172-group",
        "vit-j lwpc 175-group",
        "shem cea 281-group",
        "shem epm 295-group",
        "shem cea/epm 361-group",
        "shem epm 315-group",
        "rahab aecl 89-group",
        "ccfe 660-group",
        "ukaea 1025-group",
        "ukaea 1067-group",
        "ukaea 1102-group",
        "ukaea 142-group",
        "lanl 618-group",
        "apollo 99-group",
        "ecco  1962-group",
    ];
    if ign >= 1 && ign <= 36 {
        c.io.so(&format!("\n neutron group structure......{}", names[ign as usize]));
    }
    if ign != -1 {
        for ig in 1..=*ngn {
            c.io.so(&format!(
                " {}  {}  - {}",
                ffmt::i(ig, 5),
                ffmt::e1p(egn[ig], 12, 5),
                ffmt::e1p(egn[ig + 1], 12, 5)
            ));
        }
    }
}

/// Generate the requested gamma group structure, or read it in.
pub fn gengpg(c: &mut Ctx, igg: i32, ngg: &mut i32, egg: &mut RV) {
    let emev = 1.0e6f64;
    let ngp: i32;
    if igg == 0 {
        *ngg = 0;
        *egg = alloc(1);
        egg[1] = 0.0;
        return;
    } else if igg == 1 {
        let v = c.io.input.read(&[Kind::Int]);
        *ngg = geti(&v[0], *ngg);
        ngp = *ngg + 1;
        *egg = alloc(ngp);
        read_reals(c, ngp, egg);
    } else if igg >= 2 && igg <= 6 {
        let (n, t): (i32, &[f64]) = match igg {
            2 => (94, &G_EG2),
            3 => (12, &G_EG3),
            4 => (21, &G_EG4),
            5 => (22, &G_EG5),
            _ => (48, &G_EG6),
        };
        *ngg = n;
        ngp = n + 1;
        *egg = alloc(ngp);
        for ig in 1..=ngp {
            egg[ig] = t[(ig - 1) as usize] * emev;
        }
    } else if igg == 7 {
        *ngg = 24;
        ngp = *ngg + 1;
        *egg = alloc(ngp);
        for ig in 1..=ngp {
            egg[ig] = G_EG7[(ig - 1) as usize];
        }
    } else if igg == 8 || igg == 9 {
        *ngg = 38;
        if igg == 8 {
            *ngg = 36;
        }
        ngp = *ngg + 1;
        *egg = alloc(ngp);
        for ig in 1..=ngp {
            egg[ig] = G_EG8[(ig - 1) as usize] * emev;
        }
        if igg != 9 {
            // remove group bounds eg8(7) and eg8(39) if igg=8
            for ig in 7..=ngp {
                egg[ig] = G_EG8[ig as usize] * emev;
            }
        }
    } else if igg == 10 {
        *ngg = 42;
        ngp = *ngg + 1;
        *egg = alloc(ngp);
        for ig in 1..=ngp {
            egg[ig] = G_EG10[(ig - 1) as usize];
        }
    } else {
        c.error("gengpg", "illegal group structure.", " ");
    }

    //--display group structure
    let names: [&str; 11] = [
        "",
        "read in",
        "csewg 94 group",
        "lanl 12 group",
        "steiner 21-group",
        "straker 22 group",
        "lanl 48-group",
        "lanl 24-group",
        "vitamin-c 36-group",
        "vitamin-e 38-group",
        "vitamin-j 42-group",
    ];
    c.io.so(&format!("\n gamma group structure......{}", names[igg as usize]));
    for ig in 1..=*ngg {
        c.io.so(&format!(
            " {}  {}  - {}",
            ffmt::i(ig, 5),
            ffmt::e1p(egg[ig], 12, 5),
            ffmt::e1p(egg[ig + 1], 12, 5)
        ));
    }
}

impl Groupr {
    /// Set up the weight function, or read in an arbitrary one, or the
    /// parameters of an analytic one, or the flux calculator input.
    pub(super) fn genwtf(&mut self, c: &mut Ctx) {
        let small = 1.0e-10f64;
        let zero = 0.0f64;
        let onep5 = 1.5e0f64;

        //--read flux calculator input, if any.
        let iwtt = self.iwt.abs();
        self.nflmax = 0;
        self.sigpot = 0.0;
        if self.iwt < 0 {
            self.jsigz = 1;
            self.alpha2 = 0.0;
            self.sam = 0.0;
            self.beta = 0.0;
            self.alpha3 = 0.0;
            self.gamma = 0.0;
            let k = [
                Kind::Real,
                Kind::Real,
                Kind::Int,
                Kind::Int,
                Kind::Int,
                Kind::Real,
                Kind::Real,
                Kind::Real,
                Kind::Real,
                Kind::Real,
            ];
            let v = c.io.input.read(&k);
            self.fehi = getr(&v[0], self.fehi);
            self.sigpot = getr(&v[1], self.sigpot);
            self.nflmax = geti(&v[2], self.nflmax);
            self.ninwt = geti(&v[3], self.ninwt);
            self.jsigz = geti(&v[4], self.jsigz);
            self.alpha2 = getr(&v[5], self.alpha2);
            self.sam = getr(&v[6], self.sam);
            self.beta = getr(&v[7], self.beta);
            self.alpha3 = getr(&v[8], self.alpha3);
            self.gamma = getr(&v[9], self.gamma);
            self.ninwt = self.ninwt.abs();
            c.openz(-self.ninwt, 1);
            c.io.so(&format!(
                "\n compute flux...fehi, sigpot, nflmax ={}{}{}\n    ninwt, jsigz ={}{}",
                ffmt::f(self.fehi, 9, 1),
                ffmt::f(self.sigpot, 9, 2),
                ffmt::i(self.nflmax, 8),
                ffmt::i(-self.ninwt, 3),
                ffmt::i(self.jsigz, 4)
            ));
            if self.alpha3 != zero && self.alpha2 == zero {
                self.alpha2 = small;
            }
            if self.alpha2 > zero {
                c.io.so(&format!(
                    "       alpha2, sam ={}{}\n       beta, alpha3, gamma ={}{}{}",
                    ffmt::e1p(self.alpha2, 12, 4),
                    ffmt::f(self.sam, 8, 3),
                    ffmt::fp(self.beta, 6, 3, 0),
                    ffmt::e1p(self.alpha3, 12, 4),
                    ffmt::f(self.gamma, 8, 4)
                ));
            }
        }

        //--compute flux from data on ninwt
        if iwtt == 0 {
            c.io.so("\n compute flux......from data on ninwt");
            let v = c.io.input.read(&[Kind::Int]);
            self.ninwt = geti(&v[0], self.ninwt);
            self.ninwt = self.ninwt.abs();
            c.openz(-self.ninwt, 0);
            c.io.so(&format!("\n ninwt......{}", ffmt::i(-self.ninwt, 4)));
            return;
        }

        if iwtt == 1 {
            //--arbitrary
            c.io.so("\n weight function......read in");
            let ntmp = 10000i32;
            let mut tmp = RV::new(ntmp as usize);
            read_reals(c, ntmp, &mut tmp);
            let nr = nint(tmp[5]);
            let np = nint(tmp[6]);
            let iw = 6 + 2 * nr + 2 * np;
            if iw > ntmp {
                c.error("genwtf", "exceeded storage reading user weight function", " ");
            }
            self.wght = alloc(iw);
            for i in 1..=iw {
                self.wght[i] = tmp[i];
            }
        } else if iwtt == 2 {
            c.io.so("\n weight function......constant for all l");
        } else if iwtt == 3 {
            c.io.so("\n weight function......1/e for all l");
        } else if iwtt == 4 {
            //--1/e+fission+thermal
            let v = c.io.input.read(&[Kind::Real; 4]);
            let eb = getr(&v[0], 0.0);
            let tb = getr(&v[1], 0.0);
            let ec = getr(&v[2], 0.0);
            let tc = getr(&v[3], 0.0);
            let (ab, ac);
            if eb > 50.0 * tb {
                ab = 1.0;
                ac = 0.0;
            } else {
                ab = 1.0 / (exp(-eb / tb) * powc(eb, 2));
                ac = 1.0 / (exp(-ec / tc) * pow(ec, onep5));
            }
            self.wght = alloc(6);
            self.wght[1] = eb;
            self.wght[2] = tb;
            self.wght[3] = ab;
            self.wght[4] = ec;
            self.wght[5] = tc;
            self.wght[6] = ac;
            c.io.so(&format!(
                "\n weight function......thermal + 1/e + fission\n thermal breakpoint and temperature  {}{}\n fission breakpoint and temperature  {}{}",
                ffmt::e1p(eb, 12, 4),
                ffmt::e1p(tb, 12, 4),
                ffmt::e1p(ec, 12, 4),
                ffmt::e1p(tc, 12, 4)
            ));
        } else if iwtt == 5 {
            c.io.so("\n weight function......epri-cell lwr");
            self.wght = alloc(184);
            for i in 1..=92 {
                self.wght[i] = W_W1[(i - 1) as usize];
            }
            for i in 1..=92 {
                self.wght[i + 92] = W_W2[(i - 1) as usize];
            }
        } else if iwtt == 6 || iwtt == 7 {
            c.io.so("\n weight function......(thermal) -- (1/e) -- (fission + fusion)");
            if iwtt > 6 {
                c.io.so(&format!("{}temperature dependent", " ".repeat(22)));
            }
        } else if iwtt == 8 {
            c.io.so("\n weight function...thermal--1/e--fast reactor--fission + fusion");
            self.wght = alloc(66);
            for i in 1..=66 {
                self.wght[i] = W_W8[(i - 1) as usize];
            }
        } else if iwtt == 9 || iwtt == 10 {
            c.io.so("\n weight function......claw weight function");
            if iwtt > 9 {
                c.io.so(&format!("{}temperature dependent", " ".repeat(22)));
            }
            self.wght = alloc(106);
            for i in 1..=106 {
                self.wght[i] = W_W9[(i - 1) as usize];
            }
        } else if iwtt == 11 || iwtt == 12 {
            c.io.so("\n weight function......vitamin-e");
            if iwtt > 11 {
                c.io.so(&format!("{}temperature dependent", " ".repeat(22)));
            }
        } else {
            c.error("genwtf", "illegal weight function.", " ");
        }
    }

    /// Retrieve or compute the weight function constructed or read in
    /// by genwtf.
    pub(super) fn getwtf(&mut self, e: f64, enext: &mut f64, idis: &mut i32, _lord: i32, wtf: &mut f64, c: &Ctx) {
        let con1 = 7.45824e+07f64;
        let con2 = 1.0e0f64;
        let con3 = 1.44934e-09f64;
        let con4 = 3.90797e-02f64;
        let con5 = 2.64052e-05f64;
        let con6 = 6.76517e-02f64;
        let en1 = 0.414e0f64;
        let en2 = 2.12e6f64;
        let en3 = 1.0e7f64;
        let en4 = 1.252e7f64;
        let en5 = 1.568e7f64;
        let therm = 0.0253e0f64;
        let theta = 1.415e6f64;
        let fusion = 2.5e4f64;
        let ep = 1.407e7f64;
        let emax = 1.0e10f64;
        let s110 = 1.10e0f64;
        let s101 = 1.01e0f64;
        let s1002 = 1.002e0f64;
        let s1005 = 1.005e0f64;
        let s10001 = 1.0001e0f64;
        let tenth = 0.1e0f64;
        let half = 0.5e0f64;
        let veb = 5.0e5f64;
        let wt6a = 0.054e0f64;
        let wt6b = 1.578551e-3f64;
        let wt6c = 2.1e6f64;
        let wt6d = 2.32472e-12f64;
        let wt6e = 1.4e6f64;
        let wt6f = 2.5e4f64;
        let wt6g = 1.407e7f64;
        let wt6h = 2.51697e-11f64;
        let wt6i = 1.6e6f64;
        let wt6j = 3.3e5f64;
        let wt10a = 0.15e0f64;
        let wt10b = 300.0e0f64;
        let wt10c = 1.15e6f64;
        let exmin = -89.0e0f64;
        let two = 2.0f64;
        let zero = 0.0f64;

        //--initialize
        let mut iwtt = self.iwt.abs();
        *idis = 0;
        if e == 0.0 {
            self.gw_ip = 2;
            self.gw_ir = 1;
            self.gw_ipl = 0;
            *enext = emax;
            self.gw_step = s110;
            return;
        }
        iwtt = iwtt.abs();

        if iwtt == 1 || iwtt == 5 || iwtt == 8 || iwtt == 9 {
            //--tabulated
            c.terpa(wtf, e, enext, idis, &self.wght.at(1), &mut self.gw_ip, &mut self.gw_ir);
            if *wtf != zero {
                if self.gw_ip != self.gw_ipl {
                    self.gw_step = s10001 * pow(*enext / e, tenth);
                    self.gw_ipl = self.gw_ip;
                }
                let mut enxt = self.gw_step * e;
                if enxt > s101 * e {
                    enxt = s101 * e;
                }
                if *enext > enxt {
                    *idis = 0;
                }
                if *enext > enxt {
                    *enext = enxt;
                }
            }
        } else if iwtt == 2 {
            *wtf = 1.0;
            *enext = s101 * e;
        } else if iwtt == 3 {
            *wtf = 1.0 / e;
            *enext = s101 * e;
        } else if iwtt == 4 {
            let w = &self.wght;
            if e <= w[1] {
                *wtf = w[3] * e * exp(-e / w[2]);
                *enext = s101 * e;
                if e < w[1] && *enext > w[1] {
                    *enext = w[1];
                }
            } else if e <= w[4] {
                *wtf = 1.0 / e;
                *enext = s101 * e;
                if e < w[4] && *enext > w[4] {
                    *enext = w[4];
                }
            } else {
                *wtf = w[6] * sqrt(e) * exp(-e / w[5]);
                *enext = s101 * e;
            }
        } else if iwtt == 6 || iwtt == 7 {
            let mut tt = wt6a;
            if iwtt > 6 {
                tt = self.temp[self.jtemp] * BK;
            }
            let bb = 2.0 * tt;
            let mut cc = 1.0;
            if iwtt > 6 {
                cc = wt6b * exp(two) / powc(bb, 2);
            }
            if e <= bb {
                *wtf = cc * e * exp(-e / tt);
                *enext = s101 * e;
            } else if e <= wt6c {
                *wtf = wt6b / e;
                *enext = s101 * e;
            } else {
                *wtf = wt6d * sqrt(e) * exp(-e / wt6e);
                let pw = -powc(sqrt(e / wt6f) - sqrt(wt6g / wt6f), 2) / 2.0;
                if pw > exmin {
                    *wtf = *wtf + wt6h * exp(pw);
                }
                *enext = s101 * e;
                let mut test = wt6i;
                if (e - wt6g).abs() <= test {
                    *enext = s1005 * e;
                }
                test = wt6j;
                if (e - wt6g).abs() <= test {
                    *enext = s1002 * e;
                }
            }
        } else if iwtt == 10 {
            let ea = BK * self.temp[self.jtemp];
            let eb = wt10a * self.temp[self.jtemp] / wt10b;
            if e < eb {
                *wtf = wt10c * (e / powc(eb, 2)) * exp(-(e - eb) / ea);
                *enext = s101 * e;
                if *enext > eb {
                    *enext = eb;
                }
            } else {
                c.terpa(wtf, e, enext, idis, &self.wght.at(1), &mut self.gw_ip, &mut self.gw_ir);
                if *wtf == zero {
                    *enext = emax;
                } else {
                    if self.gw_ip != self.gw_ipl {
                        self.gw_step = s10001 * pow(*enext / e, tenth);
                        self.gw_ipl = self.gw_ip;
                    }
                    let mut enxt = self.gw_step * e;
                    if enxt > s101 * e {
                        enxt = s101 * e;
                    }
                    if *enext > enxt {
                        *idis = 0;
                    }
                    if *enext > enxt {
                        *enext = enxt;
                    }
                }
            }
        } else if iwtt == 11 || iwtt == 12 {
            *enext = s101 * e;
            if e < en1 {
                let mut tt = therm;
                if iwtt > 11 {
                    tt = self.temp[self.jtemp] * BK;
                }
                let mut cc = con1;
                if iwtt > 11 {
                    cc = con2 * exp(en1 / tt) / powc(en1, 2);
                }
                *wtf = cc * e * exp(-e / tt);
                if *enext > en1 {
                    *enext = en1;
                }
            } else if e < en2 {
                *wtf = con2 / e;
                if *enext > en2 {
                    *enext = en2;
                }
            } else if e < en3 {
                *wtf = con3 * pow(e, half) * exp(-e / theta);
                if *enext > en3 {
                    *enext = en3;
                }
            } else if e < en4 {
                *wtf = con4 / e;
                if *enext > en4 {
                    *enext = en4;
                }
            } else if e < en5 {
                *wtf = con5 * exp(-5.0 * powc(pow(e, half) - pow(ep, half), 2) / fusion);
                if (e - ep).abs() <= veb {
                    *enext = s1002 * e;
                }
                if *enext > en5 {
                    *enext = en5;
                }
            } else {
                *wtf = con6 / e;
            }
        }

        //--return enext on an even grid
        *enext = sigfig(*enext, 6, 1);
    }

    /// Flux calculator (Bondarenko model, or the slowing-down equation
    /// for an infinite mixture of absorber and moderator).
    pub(super) fn genflx(&mut self, c: &mut Ctx, b: &mut RV) {
        let up = 1.001e0f64;
        let dn = 0.999e0f64;
        let tenth = 0.1e0f64;
        let sigzmx = 1.0e10f64;
        let etop = 1.0e10f64;
        let zero = 0.0f64;
        let mut nb = 0i32;
        let mut nw: i32;
        let nsigz = self.nsigz;

        //--allocate temporary arrays
        let mut denom = RV::new(nsigz as usize + 1);
        let mut add = RV::new(nsigz as usize + 1);
        let mut factor = RV::new(nsigz as usize + 1);
        let mut tot = RV::new(nsigz as usize + 1);

        let mut nl = self.lord + 1;
        if nl == 1 && nsigz > 1 {
            nl = 2;
        }
        self.nfv = 1 + nl * nsigz + 1;
        let nfv = self.nfv;
        let mut fout = RV::new(nfv as usize + 1);
        let mut e: f64;
        let mut enext = 0.0f64;
        let mut en = 0.0f64;
        let mut idis = 0i32;
        let mut t = 0.0f64;
        let mut wtf = 0.0f64;
        let mut xc = 0.0f64;
        let mut ie: i32;
        let mut net: i32;
        let net2: i32;
        let nscr: i32;
        let mut scr1: RV;
        let mut tt = 0.0f64;
        let mut ss = 0.0f64;

        if self.nflmax == 0 {
            //--bondarenko flux model.
            nw = c.npage + 50;
            scr1 = RV::new(nw as usize + SLACK);
            nscr = self.npend;
            c.findf(self.matb, 3, 1, nscr);
            c.contio(nscr, 0, 0, &mut scr1.at(1), &mut nb, &mut nw);
            e = 0.0;
            c.gety1(e, &mut enext, &mut idis, &mut t, nscr, &mut scr1.at(1));
            if enext > up * self.egn[1] {
                c.error("genflx", "total not defined over energy range.", " ");
            }
            if self.iwt != 0 {
                self.getwtf(e, &mut en, &mut idis, 0, &mut wtf, c);
            }
            if self.iwt == 0 {
                self.getfwt(c, e, &mut en, &mut wtf, &mut xc);
            }
            ie = 0;
            net2 = 0;
            net = 0;
        } else {
            //--compute flux from slowing down equation.
            let time = c.timer();
            c.io.so(&format!("{}{}s", " ".repeat(68), ffmt::f(time, 9, 1)));
            let nin = self.npend;
            nscr = self.nfscr;
            let mut alpha = [0.0f64; 4];
            alpha[2] = self.alpha2;
            alpha[3] = self.alpha3;
            let mut nalph = 1;
            if self.alpha2 != zero {
                nalph = 2;
            }
            if self.alpha3 != zero {
                nalph = 3;
            }
            let nemax = self.nflmax;
            let nx = 3 + nsigz;
            self.felo = self.egn[1];
            let ebot = self.felo;
            if self.felo < tenth {
                self.felo = tenth;
            }
            nw = c.npage + 50;
            scr1 = RV::new(nw as usize + SLACK);
            let mut scr2 = RV::new(nw as usize + SLACK);
            e = 0.0;
            self.getwtf(e, &mut en, &mut idis, 0, &mut wtf, c);

            //--read required cross sections.
            c.repoz(nscr);
            c.findf(self.matb, 3, 1, nin);
            c.nsh = 1;
            c.tosend(nin, nscr, 0, &mut scr1.at(1));
            c.repoz(nscr);
            c.contio(nscr, 0, 0, &mut scr1.at(1), &mut nb, &mut nw);
            self.awr = scr1[2];
            alpha[1] = powc((self.awr - 1.0) / (self.awr + 1.0), 2);
            e = 0.0;
            c.gety1(e, &mut enext, &mut idis, &mut tt, nscr, &mut scr1.at(1));
            c.findf(self.matb, 3, 2, nin);
            c.contio(nin, 0, 0, &mut scr2.at(1), &mut nb, &mut nw);
            c.gety2(e, &mut en, &mut idis, &mut ss, nin, &mut scr2.at(1));
            enext = self.felo;
            let mut ne = 0;
            let mut li = 0;
            while e < self.fehi && ne < nemax {
                ne += 1;
                e = enext;
                if e > self.fehi {
                    e = self.fehi;
                }
                c.gety1(e, &mut enext, &mut idis, &mut tt, nscr, &mut scr1.at(1));
                c.gety2(e, &mut en, &mut idis, &mut ss, nin, &mut scr2.at(1));
                self.getwtf(e, &mut en, &mut idis, 0, &mut wtf, c);
                li = nx * (ne - 1);
                b[li + 1] = e;
                b[li + 2] = tt;
                b[li + 3] = ss;
                for iz in 1..=nsigz {
                    b[iz + 3 + li] = (self.sigz[iz] - self.sam) * wtf * (1.0 - self.beta);
                }
            }
            self.fehi = e;

            //--assume nr flux at fehi and above
            //--calculate contributions below fehi for nr flux
            for iz in 1..=nsigz {
                b[iz + 3 + li] = (self.sigz[iz] + self.sigpot) * wtf / (self.sigz[iz] + tt);
            }
            let lim = ne - 1;
            let (sam, beta, gamma, sigpot) = (self.sam, self.beta, self.gamma, self.sigpot);
            for k in 1..=nalph {
                for je in 1..=lim {
                    let lj = nx * (je - 1);
                    e = b[lj + 1];
                    if e >= alpha[k as usize] * self.fehi {
                        let f1 = (1.0 - alpha[k as usize] * self.fehi / e) * wtf / (1.0 - alpha[k as usize]);
                        for iz in 1..=nsigz {
                            let sz = self.sigz[iz];
                            if k == 1 {
                                b[lj + 3 + iz] = b[lj + 3 + iz] + f1 * sigpot;
                            }
                            if k == 2 {
                                b[lj + 3 + iz] = b[lj + 3 + iz] + f1 * (sam + beta * gamma * (sz - sam));
                            }
                            if k == 3 {
                                b[lj + 3 + iz] = b[lj + 3 + iz] + f1 * beta * (1.0 - gamma) * (sz - sam);
                            }
                        }
                    }
                }
            }

            //--solve slowing down equation from high to low energies
            let (mut s1, mut s2) = (0.0f64, 0.0f64);
            let mut je = ne - 1;
            while je > 0 {
                let lj = nx * (je - 1);
                let ej = b[lj + 1];
                let ejp = b[lj + 1 + nx];
                let f1 = ejp / (ejp - ej);
                let f2 = ej / (ejp - ej);
                let f3 = f1 * log(ejp / ej) - 1.0;
                let f4 = 1.0 - f2 * log(ejp / ej);
                for iz in 1..=nsigz {
                    denom[iz] = self.sigz[iz] + b[lj + 2];
                }
                for k in 1..=nalph {
                    let ak = alpha[k as usize];
                    let mut elim = ej / ak;
                    if elim > ejp {
                        elim = ejp;
                    }
                    for iz in 1..=nsigz {
                        let sz = self.sigz[iz];
                        if k == 1 {
                            s1 = b[lj + 3];
                        }
                        if k == 1 {
                            s2 = b[lj + 3 + nx];
                        }
                        if k == 2 {
                            s1 = sam + beta * gamma * (sz - sam);
                        }
                        if k == 2 {
                            s2 = s1;
                        }
                        if k == 3 {
                            s1 = beta * (1.0 - gamma) * (sz - sam);
                        }
                        if k == 3 {
                            s2 = s1;
                        }
                        let g3 = f1 * log(elim / ej) - (elim - ej) / (ejp - ej);
                        let g4 = (elim - ej) / (ejp - ej) - f2 * log(elim / ej);
                        b[lj + 3 + iz] = b[lj + 3 + iz] + g4 * s2 * b[lj + 3 + iz + nx] / (1.0 - ak);
                        denom[iz] = denom[iz] - g3 * s1 / (1.0 - ak);
                    }
                }
                for iz in 1..=nsigz {
                    b[iz + 3 + lj] = b[iz + 3 + lj] / denom[iz];
                }
                for k in 1..=nalph {
                    let ak = alpha[k as usize];
                    for iz in 1..=nsigz {
                        let sz = self.sigz[iz];
                        if k == 1 {
                            s1 = b[lj + 3];
                        }
                        if k == 1 {
                            s2 = b[lj + 3 + nx];
                        }
                        if k == 2 {
                            s1 = sam + beta * gamma * (sz - sam);
                        }
                        if k == 2 {
                            s2 = s1;
                        }
                        if k == 3 {
                            s1 = beta * (1.0 - gamma) * (sz - sam);
                        }
                        if k == 3 {
                            s2 = s1;
                        }
                        add[iz] = (f3 * s1 * b[lj + 3 + iz] + f4 * s2 * b[lj + 3 + iz + nx]) / (1.0 - ak);
                    }
                    let mut elim = 2.0 * ej;
                    let mut ie2 = je;
                    while ie2 > 1 && elim > ej {
                        ie2 -= 1;
                        let li = nx * (ie2 - 1);
                        let ei = b[li + 1];
                        elim = ei / ak;
                        if elim > ej {
                            if elim >= ejp {
                                for iz in 1..=nsigz {
                                    b[iz + 3 + li] = b[iz + 3 + li] + add[iz];
                                }
                            } else {
                                for iz in 1..=nsigz {
                                    let sz = self.sigz[iz];
                                    if k == 1 {
                                        s1 = b[lj + 3];
                                    }
                                    if k == 1 {
                                        s2 = b[lj + 3 + nx];
                                    }
                                    if k == 2 {
                                        s1 = sam + beta * gamma * (sz - sam);
                                    }
                                    if k == 2 {
                                        s2 = s1;
                                    }
                                    if k == 3 {
                                        s1 = beta * (1.0 - gamma) * (sz - sam);
                                    }
                                    if k == 3 {
                                        s2 = s1;
                                    }
                                    let g3 = f1 * log(elim / ej) - (elim - ej) / (ejp - ej);
                                    let g4 = (elim - ej) / (ejp - ej) - f2 * log(elim / ej);
                                    b[li + 3 + iz] = b[li + 3 + iz]
                                        + (g3 * s1 * b[lj + 3 + iz] + g4 * s2 * b[lj + 3 + iz + nx]) / (1.0 - ak);
                                }
                            }
                        }
                    }
                }
                je -= 1;
            }

            //--flux is converged
            //--use specified weight function below elow
            let mut enxt = 0.0f64;
            e = 0.0;
            self.getwtf(e, &mut enxt, &mut idis, 0, &mut wtf, c);
            e = b[1];
            self.getwtf(e, &mut enxt, &mut idis, 0, &mut wtf, c);
            for iz in 1..=nsigz {
                factor[iz] = b[iz + 3] / wtf;
            }
            e = 0.0;
            self.getwtf(e, &mut enxt, &mut idis, 0, &mut wtf, c);
            net = 0;
            enxt = ebot;
            while enxt < b[1] {
                net += 1;
                e = enxt;
                self.getwtf(e, &mut enxt, &mut idis, 0, &mut wtf, c);
                fout[1] = e;
                for iz in 1..=nsigz {
                    let l = (iz - 1) * nl + 1;
                    for il in 1..=nl {
                        fout[il + l] = factor[iz] * wtf;
                    }
                }
                c.loada(net, &fout.0[..], nfv, self.nflx, &mut self.wtbuf.0[..], NBUF);
            }
            c.io.so(&format!(
                " flux calculator used weight function from 1e-5 to{} ev\n {} points",
                ffmt::e1p(fout[1], 12, 4),
                ffmt::i(net, 6)
            ));
            let net1 = net;

            //--output converged flux
            let mut sigzj = sigzmx;
            let mut indx = 1;
            if self.ninwt != 0 && self.iwt != 0 {
                if self.jtemp == 1 {
                    c.repoz(self.ninwt);
                }
                let mut r = RecW::new();
                r.r(self.temp[self.jtemp]);
                c.io.write_rec(self.ninwt, r);
                sigzj = self.sigz[self.jsigz];
                indx = 2 + nl * (self.jsigz - 1);
                c.repoz(nscr);
                c.contio(nscr, 0, 0, &mut scr1.at(1), &mut nb, &mut nw);
                e = 0.0;
                c.gety1(e, &mut enext, &mut idis, &mut tt, nscr, &mut scr1.at(1));
            }
            for ie2 in 1..=ne {
                net += 1;
                let li = nx * (ie2 - 1);
                fout[1] = b[li + 1];
                for iz in 1..=nsigz {
                    let fac = (sigpot + self.sigz[iz]) / (b[li + 2] + self.sigz[iz]);
                    for il in 1..=nl {
                        let l = 1 + il + nl * (iz - 1);
                        if il == 1 {
                            fout[l] = b[li + 3 + iz];
                        }
                        if il > 1 {
                            fout[l] = fout[l - 1] * fac;
                        }
                    }
                }
                c.loada(net, &fout.0[..], nfv, self.nflx, &mut self.wtbuf.0[..], NBUF);
                if self.ninwt != 0 && self.iwt != 0 {
                    e = fout[1];
                    c.gety1(e, &mut enext, &mut idis, &mut tt, nscr, &mut scr1.at(1));
                    xc = 1.0 + tt / sigzj;
                    let mut r = RecW::new();
                    r.r(e).r(fout[indx]).r(xc);
                    c.io.write_rec(self.ninwt, r);
                }
            }
            ie = net;
            c.io.so(&format!(
                " computed flux from{} to {} ev\n {} points",
                ffmt::e1p(self.felo, 12, 4),
                ffmt::e1p(fout[1], 12, 4),
                ffmt::i(net - net1, 6)
            ));
            net2 = net;
        }

        //--complete flux using nr approximation
        let mut nett = 1;
        while nett > 0 {
            ie += 1;
            e = enext;
            if self.iwt != 0 {
                self.getwtf(e, &mut enext, &mut idis, 0, &mut wtf, c);
            }
            if self.iwt == 0 {
                self.getfwt(c, e, &mut enext, &mut wtf, &mut xc);
            }
            c.gety1(e, &mut en, &mut idis, &mut tt, nscr, &mut scr1.at(1));
            if en < enext {
                enext = en;
            }
            if en > dn * etop {
                enext = etop;
            }
            fout[1] = e;
            let ttt = tt;
            for il in 1..=nl {
                let mut mt = 1;
                if il > 1 {
                    mt = 261;
                }
                for iz in 1..=nsigz {
                    tot[iz] = ttt;
                }
                self.xtot = 0.0;
                self.getunr(c, mt, e, &mut en, &mut tot.at(1));
                if en >= zero {
                    if en < enext {
                        idis = 0;
                    }
                    if en < enext {
                        enext = en;
                    }
                }
                for iz in 1..=nsigz {
                    tt = tot[iz];
                    let mut fac = 0.0;
                    if self.iwt != 0 {
                        fac = (self.sigpot + self.sigz[iz]) / (tt + self.sigz[iz]);
                    }
                    if self.iwt == 0 {
                        fac = self.sigz[iz] * xc / (self.sigz[iz] * xc + tt);
                    }
                    let l = 1 + il + nl * (iz - 1);
                    if il == 1 {
                        fout[l] = wtf * fac;
                    }
                    if il > 1 {
                        fout[l] = fout[l - 1] * fac;
                    }
                }
            }
            fout[nfv] = self.xtot;
            net = ie;
            nett = net;
            if enext > dn * etop {
                nett = -nett;
            }
            c.loada(nett, &fout.0[..], nfv, self.nflx, &mut self.wtbuf.0[..], NBUF);
            if self.ninwt != 0 && self.iwt != 0 {
                let mut r = RecW::new();
                r.r(e).r(wtf).r(xc);
                c.io.write_rec(self.ninwt, r);
            }
        }
        self.nfp = net;
        c.io.so(&format!(
            " finished with narrow-resonance flux to{} ev\n {} points",
            ffmt::e1p(fout[1], 12, 4),
            ffmt::i(net - net2, 6)
        ));
        if self.ninwt == 0 || self.iwt == 0 {
            return;
        }
        let flag = -1.0f64;
        let mut r = RecW::new();
        r.r(flag).r(flag).r(flag);
        c.io.write_rec(self.ninwt, r);
    }

    /// Read ninwt for the weight function and xc value at e.
    pub(super) fn getfwt(&mut self, c: &mut Ctx, e: f64, enext: &mut f64, w: &mut f64, xc: &mut f64) {
        let eps = 1.0e-6f64;
        let emax = 1.0e10f64;
        let small = 1.0e-10f64;
        let zero = 0.0f64;
        let ninwt = self.ninwt;
        let rd1 = |c: &mut Ctx| -> Option<f64> {
            if c.io.rec_eof(ninwt) {
                return None;
            }
            Some(c.io.read_rec_with(ninwt, |r| r.r()))
        };
        let rd3 = |c: &mut Ctx| -> Option<(f64, f64, f64)> {
            if c.io.rec_eof(ninwt) {
                return None;
            }
            Some(c.io.read_rec_with(ninwt, |r| (r.r(), r.r(), r.r())))
        };

        //--initialize when e=0
        if !(e > zero) {
            if self.jtemp == 1 || self.temp[self.jtemp] < self.fw_templ {
                c.repoz(ninwt);
            }
            self.fw_el = 0.0;
            self.fw_wl = 0.0;
            self.fw_xl = 0.0;
            self.fw_en = 0.0;
            self.fw_wn = 0.0;
            self.fw_xn = 0.0;
            let tj = self.temp[self.jtemp];
            // search tape for desired temperature
            let mut lbl = 110;
            let mut tempn = 0.0f64;
            if self.jtemp != 1 {
                let flag = match rd1(c) {
                    Some(x) => x,
                    None => fatal_eof(ninwt),
                };
                if flag < zero {
                    lbl = 110;
                } else if (tj - flag).abs() > eps {
                    lbl = 130;
                } else {
                    tempn = flag;
                    lbl = 150;
                }
            }
            loop {
                match lbl {
                    110 => match rd1(c) {
                        None => lbl = 145,
                        Some(x) => {
                            tempn = x;
                            lbl = if (tempn - tj).abs() <= eps { 150 } else { 130 };
                        }
                    },
                    130 => match rd1(c) {
                        None => lbl = 145,
                        Some(flag) => {
                            lbl = if flag >= zero { 130 } else { 110 };
                        }
                    },
                    145 => {
                        let s = format!("temperature {}", ffmt::e1p(tj, 11, 4));
                        c.error("getfwt", &s, " ");
                    }
                    _ => break,
                }
            }
            self.fw_templ = tempn;
            *enext = emax;
            return;
        }

        //--search for e on ninwt
        if (e - self.fw_el).abs() < self.fw_el * small && self.fw_en < 0.0 {
            *w = self.fw_wl;
            *xc = self.fw_xl;
            *enext = self.fw_en;
            return;
        }
        let mut lbl;
        if (e - self.fw_en).abs() < self.fw_en * small {
            lbl = 240;
        } else if e > self.fw_el * (1.0 - small) && e < self.fw_en * (1.0 + small) {
            lbl = 230;
        } else {
            if e < self.fw_el * (1.0 - small) {
                c.error("getfwt", "requested e is out of order", " ");
            }
            lbl = 210;
        }
        loop {
            match lbl {
                210 => {
                    self.fw_el = self.fw_en;
                    self.fw_wl = self.fw_wn;
                    self.fw_xl = self.fw_xn;
                    match rd3(c) {
                        None => c.error("getfwt", "e outside range of data", " "),
                        Some((a, b2, d)) => {
                            self.fw_en = a;
                            self.fw_wn = b2;
                            self.fw_xn = d;
                        }
                    }
                    if (e - self.fw_en).abs() < self.fw_en * small {
                        lbl = 240;
                    } else if e > self.fw_el * (1.0 - small) && e < self.fw_en * (1.0 + small) {
                        lbl = 230;
                    }
                }
                230 => {
                    c.terp1(self.fw_el, self.fw_wl, self.fw_en, self.fw_wn, e, w, 2);
                    c.terp1(self.fw_el, self.fw_xl, self.fw_en, self.fw_xn, e, xc, 2);
                    *enext = self.fw_en;
                    return;
                }
                _ => {
                    // 240
                    *w = self.fw_wn;
                    *xc = self.fw_xn;
                    self.fw_el = self.fw_en;
                    self.fw_wl = self.fw_wn;
                    self.fw_xl = self.fw_xn;
                    let (a, b2, d) = match rd3(c) {
                        Some(x) => x,
                        None => fatal_eof(ninwt),
                    };
                    self.fw_en = a;
                    self.fw_wn = b2;
                    self.fw_xn = d;
                    *enext = self.fw_en;
                    *enext = emax;
                    return;
                }
            }
        }
    }
}

/// End of file on an unformatted read without END=.
fn fatal_eof(nun: i32) -> ! {
    crate::fio::fatal_runtime(&format!("Fortran runtime error: End of file (unit {})", nun))
}
