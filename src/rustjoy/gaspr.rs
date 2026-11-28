//! GASPR: gas production cross sections (MT203-207) on the PENDF tape,
//! ported from NJOY2016 `gaspr.f90` (module gasprm).
//!
//! The routine keeps the Fortran structure: labels become labeled loops
//! and blocks, and arrays keep the Fortran (1-based) indices. The long
//! reaction tables (which MT gives which gas yield) are written out in the
//! same order as the Fortran `if` chains.

#![allow(clippy::too_many_lines)]
#![allow(clippy::collapsible_if)]
#![allow(clippy::collapsible_else_if)]
#![allow(clippy::assign_op_pattern)]

use crate::endf::Ctx;
use crate::farr::{RM, RV};
use crate::ffmt;
use crate::fio::{Kind, Val};
use crate::fmath::nint;
use crate::util::sigfig;

const EMAX: f64 = 1.0e10;
const ZERO: f64 = 0.0;
const HUND: f64 = 100.0;
/// Extra words past the Fortran array sizes, so that writes the Fortran
/// does past an array end (without checking) do not stop the port first.
const SLACK: usize = 20000;

fn geti(v: &Option<Val>, d: i32) -> i32 {
    match v {
        Some(Val::Int(i)) => *i,
        Some(Val::Real(x)) => *x as i32,
        _ => d,
    }
}

/// Add gas production reactions (MT203-207) to the PENDF tape.
pub fn gaspr(c: &mut Ctx) {
    let (mut nb, mut nw) = (0i32, 0i32);
    let mut a = RV::new(10000 + SLACK);
    let mut b = RV::new(6 + 100);
    let mut six = RV::new(5000 + SLACK);
    let mut egas: RV;
    let mut sgas: RM;
    let mut nscr1 = 10;
    // set iprint nonzero for gas production printout on listing
    let iprint = 0;
    let mut nx = 0i32;
    let mut en: f64;
    let (mut enext, mut y) = (0.0f64, 0.0f64);
    let mut idis = 0i32;
    let (mut yyy, mut xnext) = (0.0f64, 0.0f64);
    let (mut l203, mut l204, mut l205, mut l206, mut l207) = (0i32, 0i32, 0i32, 0i32, 0i32);
    let mut ngas: i32;
    let (mut mfb, mut mtb) = (0i32, 0i32);

    //--read user input
    let time = c.timer();
    c.io.so(&format!(
        "\n gaspr...add gas production cross sections{}{}s",
        " ".repeat(27),
        ffmt::f(time, 8, 1)
    ));
    c.io.se(&format!(
        "\n gaspr...{}{}s",
        " ".repeat(60),
        ffmt::f(time, 8, 1)
    ));
    let v = c.io.input.read(&[Kind::Int, Kind::Int, Kind::Int]);
    let nendf = geti(&v[0], 0);
    let npend = geti(&v[1], 0);
    let noutp = geti(&v[2], 0);
    if npend < 0 {
        nscr1 = -nscr1;
    }
    if npend * noutp < 0 {
        c.error(
            "gaspr",
            "npend and noutp must both be binary",
            "or both be coded",
        );
    }
    c.io.so(&format!(
        "\n units:{}{}{}",
        ffmt::i(nendf, 6),
        ffmt::i(npend, 6),
        ffmt::i(noutp, 6)
    ));
    c.openz(nendf, 0);
    c.openz(npend, 0);
    c.openz(noutp, 1);
    c.openz(nscr1, 1);

    //--check endf tape for mf6,mt5.
    //--also set flags for the absence or presence of charged-particle
    //--reactions and the endf version dependent partial cross section
    //--mt ranges.
    c.repoz(nendf);
    c.tpidio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
    c.contio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
    let za = c.c1h;
    let awr = c.c2h;
    let n6 = c.n2h;
    c.contio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
    if c.n1h != 0 {
        c.iverf = 4;
        nx = n6;
    } else if c.n2h == 0 {
        c.iverf = 5;
    } else {
        c.iverf = 6;
    }
    c.skiprz(nendf, -1);
    if c.iverf >= 5 {
        c.contio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
    }
    let mut zain = 1.0f64;
    if c.iverf >= 6 {
        c.contio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
        let nsub = c.n1h;
        zain = (nsub / 10) as f64;
    }
    c.hdatio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
    let matd = c.math;
    while nb != 0 {
        c.moreio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
    }
    if c.iverf != 4 {
        nx = c.n2h;
    }
    nw = nx;
    c.dictio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
    let mut mf6mt5 = 0;
    let mut i = 1i32;
    let mut idone = 0;
    let (mpmin, mpmax, mdmin, mdmax, mtmin, mtmax, m3min, m3max, m4min, m4max) = if c.iverf >= 6 {
        (600, 649, 650, 699, 700, 749, 750, 799, 800, 849)
    } else {
        (700, 718, 720, 738, 740, 758, 760, 768, 780, 798)
    };
    while i <= nw && idone == 0 {
        let mfi = nint(a[1 + i + 1]);
        let mti = nint(a[1 + i + 2]);
        if mfi > 6 {
            idone = 1;
        } else {
            if mfi == 6 && mti == 5 {
                mf6mt5 = 1;
            }
            i += 6;
        }
    }
    c.repoz(nendf);
    if mf6mt5 != 0 {
        c.io.so("\n mf6,mt5 found");
        let mfd = 6;
        let mtd = 5;
        c.findf(matd, mfd, mtd, nendf);
        c.contio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
        let nk = c.n1h;
        let mut lsix = 1i32;
        l203 = 0;
        l204 = 0;
        l205 = 0;
        l206 = 0;
        l207 = 0;
        for _ik in 1..=nk {
            c.tab1io(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
            let izap = nint(c.c1h);
            let law = c.l2h;
            let mut ll = 1 + nw;
            while nb != 0 {
                c.moreio(nendf, 0, 0, &mut a.at(ll), &mut nb, &mut nw);
                ll += nw;
            }
            let lp = match izap {
                1001 => Some(&mut l203),
                1002 => Some(&mut l204),
                1003 => Some(&mut l205),
                2003 => Some(&mut l206),
                2004 => Some(&mut l207),
                _ => None,
            };
            if let Some(lp) = lp {
                *lp = lsix;
                for i in 1..=nw {
                    six[lsix + i - 1] = a[i];
                }
                lsix += nw;
            }
            if law == 6 {
                c.contio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
            } else if law == 1 || law == 2 || law == 5 {
                c.tab2io(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
                let ne = c.n2h;
                for _ie in 1..=ne {
                    c.listio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
                    while nb != 0 {
                        c.moreio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
                    }
                }
            } else if law == 7 {
                c.tab2io(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
                let ne = c.n2h;
                for _ie in 1..=ne {
                    c.tab2io(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
                    let nmu = c.n2h;
                    for _imu in 1..=nmu {
                        c.tab1io(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
                        while nb != 0 {
                            c.moreio(nendf, 0, 0, &mut a.all(), &mut nb, &mut nw);
                        }
                    }
                }
            }
        }
    }

    //--loop over all temperatures on the pendf tape
    c.repoz(npend);
    c.repoz(noutp);
    c.tpidio(npend, noutp, 0, &mut a.all(), &mut nb, &mut nw);
    let mut itemp = 0;
    loop {
        // 225
        itemp += 1;
        c.contio(npend, 0, 0, &mut b.all(), &mut nb, &mut nw);
        if c.mfh == 0 {
            break;
        }

        //--copy data through file3, include any partial
        //-charged-particle cross sections.
        c.repoz(nscr1);
        c.contio(0, nscr1, 0, &mut b.all(), &mut nb, &mut nw);
        c.tofend(npend, nscr1, 0, &mut a.all());
        c.tofend(npend, nscr1, 0, &mut a.all());
        idone = 0;
        while idone == 0 {
            c.contio(npend, 0, 0, &mut b.all(), &mut nb, &mut nw);
            mfb = c.mfh;
            mtb = c.mth;
            if c.mth > m4max || c.mth == 0 {
                idone = 1;
            } else {
                c.contio(0, nscr1, 0, &mut b.all(), &mut nb, &mut nw);
                c.tosend(npend, nscr1, 0, &mut a.all());
            }
        }
        c.contio(0, nscr1, 0, &mut b.all(), &mut nb, &mut nw);

        //--reposition npend to the location where gas production
        //--cross sections would go.
        c.repoz(npend);
        c.tpidio(npend, 0, 0, &mut a.all(), &mut nb, &mut nw);
        if itemp > 1 {
            for _i in 2..=itemp {
                c.tomend(npend, 0, 0, &mut a.all());
            }
        }
        c.tofend(npend, 0, 0, &mut a.all());
        c.tofend(npend, 0, 0, &mut a.all());
        idone = 0;
        while idone == 0 {
            c.contio(npend, 0, 0, &mut b.all(), &mut nb, &mut nw);
            if c.mth > 200 || c.mth == 0 {
                idone = 1;
            } else {
                c.tosend(npend, 0, 0, &mut a.all());
            }
        }

        //--check the reactions in file 3 for gas threshold
        c.repoz(nscr1);
        c.tofend(nscr1, 0, 0, &mut a.all());
        c.tofend(nscr1, 0, 0, &mut a.all());
        let mut thrg = EMAX;
        'l240: loop {
            c.contio(nscr1, 0, 0, &mut a.all(), &mut nb, &mut nw);
            'b245: {
                if c.mth > 200 && c.mth < mpmin {
                    break 'b245;
                }
                if c.mth > m4max || c.mth == 0 {
                    break 'l240;
                }
                let mut izg = 0;
                let mut izr = nint(za + zain);
                en = 0.0;
                c.gety1(en, &mut enext, &mut idis, &mut y, nscr1, &mut a.all());
                let lr = c.l2h;
                if c.mth >= 1 && c.mth <= 4 {
                    break 'b245;
                }
                if c.mth == 5 && mf6mt5 == 1 {
                    izg = 1;
                }
                if c.mth >= 6 && c.mth <= 10 {
                    break 'b245;
                }
                if c.mth == 11 {
                    izg = 1;
                }
                if c.mth == 11 {
                    izr = izr - 1004;
                }
                if c.mth >= 12 && c.mth <= 15 {
                    break 'b245;
                }
                if c.mth == 16 {
                    izr = izr - 2;
                }
                if c.mth == 17 {
                    izr = izr - 3;
                }
                if c.mth >= 18 && c.mth <= 21 {
                    break 'b245;
                }
                if c.mth >= 22 && c.mth <= 25 {
                    izg = 1;
                }
                if c.mth == 22 {
                    izr = izr - 2005;
                }
                if c.mth == 23 {
                    izr = izr - 6013;
                }
                if c.mth == 24 {
                    izr = izr - 2006;
                }
                if c.mth == 25 {
                    izr = izr - 2007;
                }
                if c.mth >= 26 && c.mth <= 27 {
                    break 'b245;
                }
                if c.mth >= 28 && c.mth <= 32 {
                    izg = 1;
                }
                if c.mth == 28 {
                    izr = izr - 1002;
                }
                if c.mth == 29 {
                    izr = izr - 4009;
                }
                if c.mth == 30 {
                    izr = izr - 4010;
                }
                if c.mth == 31 {
                    break 'b245;
                }
                if c.mth >= 32 && c.mth <= 36 {
                    izg = 1;
                }
                if c.mth == 32 {
                    izr = izr - 1003;
                }
                if c.mth == 33 {
                    izr = izr - 1004;
                }
                if c.mth == 34 {
                    izr = izr - 2004;
                }
                if c.mth == 35 {
                    izr = izr - 5011;
                }
                if c.mth == 36 {
                    izr = izr - 5012;
                }
                if c.mth == 37 {
                    izr = izr - 4;
                }
                if c.mth >= 38 && c.mth <= 40 {
                    break 'b245;
                }
                if c.mth >= 41 && c.mth <= 42 {
                    izg = 1;
                }
                if c.mth == 41 {
                    izr = izr - 1003;
                }
                if c.mth == 42 {
                    izr = izr - 1004;
                }
                if c.mth == 43 {
                    break 'b245;
                }
                if c.mth >= 44 && c.mth <= 45 {
                    izg = 1;
                }
                if c.mth == 44 {
                    izr = izr - 2003;
                }
                if c.mth == 45 {
                    izr = izr - 3006;
                }
                if c.mth >= 46 && c.mth <= 50 {
                    break 'b245;
                }
                if c.mth >= 51 && c.mth <= 91 {
                    izr = izr - 1;
                    if lr >= 22 && lr <= 25 {
                        izg = 1;
                    }
                    if lr == 22 {
                        izr = izr - 2004;
                    }
                    if lr == 23 {
                        izr = izr - 6012;
                    }
                    if lr == 24 {
                        izr = izr - 2005;
                    }
                    if lr == 25 {
                        izr = izr - 2006;
                    }
                    if lr >= 28 && lr <= 30 {
                        izg = 1;
                    }
                    if lr == 28 {
                        izr = izr - 1001;
                    }
                    if lr == 29 {
                        izr = izr - 4008;
                    }
                    if lr == 30 {
                        izr = izr - 4009;
                    }
                    if lr >= 32 && lr <= 36 {
                        izg = 1;
                    }
                    if lr == 32 {
                        izr = izr - 1002;
                    }
                    if lr == 33 {
                        izr = izr - 1003;
                    }
                    if lr == 34 {
                        izr = izr - 2003;
                    }
                    if lr == 35 {
                        izr = izr - 5010;
                    }
                    if lr == 36 {
                        izr = izr - 5011;
                    }
                }
                if c.mth >= 92 && c.mth <= 101 {
                    break 'b245;
                }
                if c.mth >= 103 && c.mth <= 117 {
                    izg = 1;
                }
                if c.mth == 103 {
                    izr = izr - 1001;
                }
                if c.mth == 104 {
                    izr = izr - 1002;
                }
                if c.mth == 105 {
                    izr = izr - 1003;
                }
                if c.mth == 106 {
                    izr = izr - 2003;
                }
                if c.mth == 107 {
                    izr = izr - 2004;
                }
                //--always skip levels for mt103-107
                if c.mth >= mpmin && c.mth <= mpmax {
                    break 'b245;
                }
                if c.mth >= mdmin && c.mth <= mdmax {
                    break 'b245;
                }
                if c.mth >= mtmin && c.mth <= mtmax {
                    break 'b245;
                }
                if c.mth >= m3min && c.mth <= m3max {
                    break 'b245;
                }
                if c.mth >= m4min && c.mth <= m4max {
                    break 'b245;
                }
                if c.mth == 108 {
                    izr = izr - 4008;
                }
                if c.mth == 109 {
                    izr = izr - 6012;
                }
                if c.mth == 111 {
                    izr = izr - 2002;
                }
                if c.mth == 112 {
                    izr = izr - 3005;
                }
                if c.mth == 113 {
                    izr = izr - 5011;
                }
                if c.mth == 114 {
                    izr = izr - 5010;
                }
                if c.mth == 115 {
                    izr = izr - 2003;
                }
                if c.mth == 116 {
                    izr = izr - 2004;
                }
                if c.mth == 117 {
                    izr = izr - 3006;
                }
                if c.mth >= 154 && c.mth <= 159 {
                    izg = 1;
                }
                if c.mth >= 162 && c.mth <= 200 {
                    izg = 1;
                }
                if c.mth == 152 {
                    izr = izr - 5;
                }
                if c.mth == 153 {
                    izr = izr - 6;
                }
                if c.mth == 154 {
                    izr = izr - 1005;
                }
                if c.mth == 155 {
                    izr = izr - 3007;
                }
                if c.mth == 156 {
                    izr = izr - 1005;
                }
                if c.mth == 157 {
                    izr = izr - 1005;
                }
                if c.mth == 158 {
                    izr = izr - 3007;
                }
                if c.mth == 159 {
                    izr = izr - 3007;
                }
                if c.mth == 160 {
                    izr = izr - 7;
                }
                if c.mth == 161 {
                    izr = izr - 8;
                }
                if c.mth == 162 {
                    izr = izr - 1006;
                }
                if c.mth == 163 {
                    izr = izr - 1007;
                }
                if c.mth == 164 {
                    izr = izr - 1008;
                }
                if c.mth == 165 {
                    izr = izr - 2008;
                }
                if c.mth == 166 {
                    izr = izr - 2009;
                }
                if c.mth == 167 {
                    izr = izr - 2010;
                }
                if c.mth == 168 {
                    izr = izr - 2011;
                }
                if c.mth == 169 {
                    izr = izr - 1006;
                }
                if c.mth == 170 {
                    izr = izr - 1007;
                }
                if c.mth == 171 {
                    izr = izr - 1008;
                }
                if c.mth == 172 {
                    izr = izr - 1006;
                }
                if c.mth == 173 {
                    izr = izr - 1007;
                }
                if c.mth == 174 {
                    izr = izr - 1008;
                }
                if c.mth == 175 {
                    izr = izr - 1009;
                }
                if c.mth == 176 {
                    izr = izr - 2005;
                }
                if c.mth == 177 {
                    izr = izr - 2006;
                }
                if c.mth == 178 {
                    izr = izr - 2007;
                }
                if c.mth == 179 {
                    izr = izr - 2005;
                }
                if c.mth == 180 {
                    izr = izr - 4011;
                }
                if c.mth == 181 {
                    izr = izr - 3008;
                }
                if c.mth == 182 {
                    izr = izr - 2005;
                }
                if c.mth == 183 {
                    izr = izr - 2004;
                }
                if c.mth == 184 {
                    izr = izr - 2005;
                }
                if c.mth == 185 {
                    izr = izr - 2006;
                }
                if c.mth == 186 {
                    izr = izr - 3005;
                }
                if c.mth == 187 {
                    izr = izr - 3006;
                }
                if c.mth == 188 {
                    izr = izr - 3007;
                }
                if c.mth == 189 {
                    izr = izr - 3008;
                }
                if c.mth == 190 {
                    izr = izr - 2004;
                }
                if c.mth == 191 {
                    izr = izr - 3004;
                }
                if c.mth == 192 {
                    izr = izr - 3005;
                }
                if c.mth == 193 {
                    izr = izr - 4007;
                }
                if c.mth == 194 {
                    izr = izr - 2006;
                }
                if c.mth == 195 {
                    izr = izr - 4012;
                }
                if c.mth == 196 {
                    izr = izr - 3009;
                }
                if c.mth == 197 {
                    izr = izr - 3003;
                }
                if c.mth == 198 {
                    izr = izr - 3004;
                }
                if c.mth == 199 {
                    izr = izr - 4009;
                }
                if c.mth == 200 {
                    izr = izr - 2007;
                }
                if izr == 4008 {
                    izg = 1;
                }
                if izg == 0 && (izr > 2004 || izr <= 0) {
                    break 'b245;
                }
                if enext < thrg {
                    thrg = enext;
                }
            }
            // 245
            c.tosend(nscr1, 0, 0, &mut a.all());
        }
        // 250
        if itemp == 1 {
            c.io.so(&format!(
                "\n the gas production threshold is{} ev",
                ffmt::e1p(thrg, 12, 4)
            ));
        }

        //--read through the scratch tape
        //--collecting gas production in the process
        //--use the energy grid of mt1 starting at thrg
        c.repoz(nscr1);
        c.tofend(nscr1, 0, 0, &mut a.all());
        c.tofend(nscr1, 0, 0, &mut a.all());
        c.contio(nscr1, 0, 0, &mut a.all(), &mut nb, &mut nw);
        en = 0.0;
        c.gety1(en, &mut enext, &mut idis, &mut y, nscr1, &mut a.all());
        let maxg = nint(a[6]);
        egas = RV::new(maxg.max(0) as usize);
        sgas = RM::new(5, maxg.max(0) as usize);
        enext = thrg;
        i = 0;
        while enext < EMAX {
            i += 1;
            if i > maxg {
                c.error("gaspr", "too many gas production energy points", " ");
            }
            en = enext;
            c.gety1(en, &mut enext, &mut idis, &mut y, nscr1, &mut a.all());
            egas[i] = en;
        }
        ngas = i;
        if itemp == 1 {
            c.io.so(&format!("\n found{} points", ffmt::i(ngas, 6)));
        }
        c.tosend(nscr1, 0, 0, &mut a.all());
        if itemp == 1 {
            c.io.so("\n pendf mt  mt203  mt204  mt205  mt206  mt207\n ________  _____  _____  _____  _____  _____");
        }
        for i in 1..=ngas {
            for j in 1..=5 {
                sgas[(j, i)] = 0.0;
            }
        }

        //--loop over other reactions and
        //--sum up gas production values
        'l270: loop {
            c.contio(nscr1, 0, 0, &mut a.all(), &mut nb, &mut nw);
            'b310: {
                if c.mth > 200 && c.mth < mpmin {
                    break 'b310;
                }
                if c.mth > m4max || c.mth == 0 {
                    break 'l270;
                }
                if c.mth <= 4 {
                    break 'b310;
                }
                if c.mth >= 6 && c.mth <= 10 {
                    break 'b310;
                }
                if c.mth >= 12 && c.mth <= 15 {
                    break 'b310;
                }
                if c.mth >= 18 && c.mth <= 21 {
                    break 'b310;
                }
                if c.mth >= 38 && c.mth <= 40 {
                    break 'b310;
                }
                if c.mth == 43 {
                    break 'b310;
                }
                if c.mth >= 46 && c.mth <= 50 {
                    break 'b310;
                }
                if c.mth >= 92 && c.mth <= 101 {
                    break 'b310;
                }
                //--always skip levels for mt103-107
                if c.mth >= mpmin && c.mth <= mpmax {
                    break 'b310;
                }
                if c.mth >= mdmin && c.mth <= mdmax {
                    break 'b310;
                }
                if c.mth >= mtmin && c.mth <= mtmax {
                    break 'b310;
                }
                if c.mth >= m3min && c.mth <= m3max {
                    break 'b310;
                }
                if c.mth >= m4min && c.mth <= m4max {
                    break 'b310;
                }
                if c.mth == 152 || c.mth == 153 || c.mth == 160 || c.mth == 161 {
                    break 'b310;
                }
                en = 0.0;
                c.gety1(en, &mut enext, &mut idis, &mut y, nscr1, &mut a.all());
                let lr = c.l2h;
                let mut izr = nint(za + zain);
                let (mut y203, mut y204, mut y205, mut y206, mut y207);
                y203 = 0.0;
                y204 = 0.0;
                y205 = 0.0;
                y206 = 0.0;
                y207 = 0.0;
                if c.mth == 5 && mf6mt5 == 1 {
                    if l203 > 0 {
                        y203 = 111.0;
                    }
                    if l204 > 0 {
                        y204 = 111.0;
                    }
                    if l205 > 0 {
                        y205 = 111.0;
                    }
                    if l206 > 0 {
                        y206 = 111.0;
                    }
                    if l207 > 0 {
                        y207 = 111.0;
                    }
                } else if c.mth == 11 {
                    izr = izr - 1004;
                    y204 = 1.0;
                } else if c.mth == 16 {
                    izr = izr - 2;
                } else if c.mth == 17 {
                    izr = izr - 3;
                } else if c.mth == 22 {
                    izr = izr - 2005;
                    y207 = 1.0;
                } else if c.mth == 23 {
                    izr = izr - 6013;
                    y207 = 3.0;
                } else if c.mth == 24 {
                    izr = izr - 2006;
                    y207 = 1.0;
                } else if c.mth == 25 {
                    izr = izr - 2007;
                    y207 = 1.0;
                } else if c.mth == 28 {
                    izr = izr - 1002;
                    y203 = 1.0;
                } else if c.mth == 29 {
                    izr = izr - 4009;
                    y207 = 2.0;
                } else if c.mth == 30 {
                    izr = izr - 4010;
                    y207 = 2.0;
                } else if c.mth == 32 {
                    izr = izr - 1003;
                    y204 = 1.0;
                } else if c.mth == 33 {
                    izr = izr - 1004;
                    y205 = 1.0;
                } else if c.mth == 34 {
                    izr = izr - 2004;
                    y206 = 1.0;
                } else if c.mth == 35 {
                    izr = izr - 5011;
                    y204 = 1.0;
                    y207 = 2.0;
                } else if c.mth == 36 {
                    izr = izr - 5012;
                    y205 = 1.0;
                    y207 = 2.0;
                } else if c.mth == 37 {
                    izr = izr - 4;
                } else if c.mth == 41 {
                    izr = izr - 1003;
                    y203 = 1.0;
                } else if c.mth == 42 {
                    izr = izr - 1004;
                    y203 = 1.0;
                } else if c.mth == 44 {
                    izr = izr - 2003;
                    y203 = 2.0;
                } else if c.mth == 45 {
                    izr = izr - 3006;
                    y203 = 1.0;
                    y207 = 1.0;
                } else if c.mth >= 51 && c.mth <= 91 {
                    izr = izr - 1;
                    if lr == 22 {
                        izr = izr - 2004;
                        y207 = 1.0;
                    } else if lr == 23 {
                        izr = izr - 6012;
                        y207 = 3.0;
                    } else if lr == 24 {
                        izr = izr - 2005;
                        y207 = 1.0;
                    } else if lr == 25 {
                        izr = izr - 2006;
                        y207 = 1.0;
                    } else if lr == 28 {
                        izr = izr - 1001;
                        y203 = 1.0;
                    } else if lr == 29 {
                        izr = izr - 4008;
                        y207 = 2.0;
                    } else if lr == 30 {
                        izr = izr - 4009;
                        y207 = 2.0;
                    } else if lr == 32 {
                        izr = izr - 1002;
                        y204 = 1.0;
                    } else if lr == 33 {
                        izr = izr - 1003;
                        y205 = 1.0;
                    } else if lr == 34 {
                        izr = izr - 2003;
                        y206 = 1.0;
                    } else if lr == 35 {
                        izr = izr - 5010;
                        y204 = 1.0;
                        y207 = 2.0;
                    } else if lr == 36 {
                        izr = izr - 5011;
                        y205 = 1.0;
                        y207 = 2.0;
                    } else if lr == 39 {
                        {}
                    } else if lr == 40 {
                        {}
                    }
                } else if c.mth == 103 {
                    izr = izr - 1001;
                    y203 = 1.0;
                } else if c.mth == 104 {
                    izr = izr - 1002;
                    y204 = 1.0;
                } else if c.mth == 105 {
                    izr = izr - 1003;
                    y205 = 1.0;
                } else if c.mth == 106 {
                    izr = izr - 2003;
                    y206 = 1.0;
                } else if c.mth == 107 {
                    izr = izr - 2004;
                    y207 = 1.0;
                } else if c.mth == 108 {
                    izr = izr - 4008;
                    y207 = 2.0;
                } else if c.mth == 109 {
                    izr = izr - 6012;
                    y207 = 3.0;
                } else if c.mth == 111 {
                    izr = izr - 2002;
                    y203 = 2.0;
                } else if c.mth == 112 {
                    izr = izr - 3005;
                    y203 = 1.0;
                    y207 = 1.0;
                } else if c.mth == 113 {
                    izr = izr - 5011;
                    y205 = 1.0;
                    y207 = 2.0;
                } else if c.mth == 114 {
                    izr = izr - 5010;
                    y204 = 1.0;
                    y207 = 2.0;
                } else if c.mth == 115 {
                    izr = izr - 2003;
                    y203 = 1.0;
                    y204 = 1.0;
                } else if c.mth == 116 {
                    izr = izr - 2004;
                    y203 = 1.0;
                    y205 = 1.0;
                } else if c.mth == 117 {
                    izr = izr - 3006;
                    y204 = 1.0;
                    y207 = 1.0;
                } else if c.mth == 152 {
                    izr = izr - 5;
                } else if c.mth == 153 {
                    izr = izr - 6;
                } else if c.mth == 154 {
                    izr = izr - 1005;
                    y205 = 1.0;
                } else if c.mth == 155 {
                    izr = izr - 3007;
                    y205 = 1.0;
                    y207 = 1.0;
                } else if c.mth == 156 {
                    izr = izr - 1005;
                    y203 = 1.0;
                } else if c.mth == 157 {
                    izr = izr - 1005;
                    y204 = 1.0;
                } else if c.mth == 158 {
                    izr = izr - 3007;
                    y204 = 1.0;
                    y207 = 1.0;
                } else if c.mth == 159 {
                    izr = izr - 3007;
                    y203 = 1.0;
                    y207 = 1.0;
                } else if c.mth == 160 {
                    izr = izr - 7;
                } else if c.mth == 161 {
                    izr = izr - 8;
                } else if c.mth == 162 {
                    izr = izr - 1006;
                    y203 = 1.0;
                } else if c.mth == 163 {
                    izr = izr - 1007;
                    y203 = 1.0;
                } else if c.mth == 164 {
                    izr = izr - 1008;
                    y203 = 1.0;
                } else if c.mth == 165 {
                    izr = izr - 2008;
                    y207 = 1.0;
                } else if c.mth == 166 {
                    izr = izr - 2009;
                    y207 = 1.0;
                } else if c.mth == 167 {
                    izr = izr - 2010;
                    y207 = 1.0;
                } else if c.mth == 168 {
                    izr = izr - 2011;
                    y207 = 1.0;
                } else if c.mth == 169 {
                    izr = izr - 1006;
                    y204 = 1.0;
                } else if c.mth == 170 {
                    izr = izr - 1007;
                    y204 = 1.0;
                } else if c.mth == 171 {
                    izr = izr - 1008;
                    y204 = 1.0;
                } else if c.mth == 172 {
                    izr = izr - 1006;
                    y205 = 1.0;
                } else if c.mth == 173 {
                    izr = izr - 1007;
                    y205 = 1.0;
                } else if c.mth == 174 {
                    izr = izr - 1008;
                    y205 = 1.0;
                } else if c.mth == 175 {
                    izr = izr - 1009;
                    y205 = 1.0;
                } else if c.mth == 176 {
                    izr = izr - 2005;
                    y206 = 1.0;
                } else if c.mth == 177 {
                    izr = izr - 2006;
                    y206 = 1.0;
                } else if c.mth == 178 {
                    izr = izr - 2007;
                    y206 = 1.0;
                } else if c.mth == 179 {
                    izr = izr - 2005;
                    y203 = 2.0;
                } else if c.mth == 180 {
                    izr = izr - 4011;
                    y207 = 2.0;
                } else if c.mth == 181 {
                    izr = izr - 3008;
                    y203 = 1.0;
                    y207 = 1.0;
                } else if c.mth == 182 {
                    izr = izr - 2005;
                    y204 = 1.0;
                    y205 = 1.0;
                } else if c.mth == 183 {
                    izr = izr - 2004;
                    y203 = 1.0;
                    y204 = 1.0;
                } else if c.mth == 184 {
                    izr = izr - 2005;
                    y203 = 1.0;
                    y205 = 1.0;
                } else if c.mth == 185 {
                    izr = izr - 2006;
                    y204 = 1.0;
                    y205 = 1.0;
                } else if c.mth == 186 {
                    izr = izr - 3005;
                    y203 = 1.0;
                    y206 = 1.0;
                } else if c.mth == 187 {
                    izr = izr - 3006;
                    y204 = 1.0;
                    y206 = 1.0;
                } else if c.mth == 188 {
                    izr = izr - 3007;
                    y205 = 1.0;
                    y206 = 1.0;
                } else if c.mth == 189 {
                    izr = izr - 3008;
                    y205 = 1.0;
                    y207 = 1.0;
                } else if c.mth == 190 {
                    izr = izr - 2004;
                    y203 = 2.0;
                } else if c.mth == 191 {
                    izr = izr - 3004;
                    y203 = 1.0;
                    y206 = 1.0;
                } else if c.mth == 192 {
                    izr = izr - 3005;
                    y204 = 1.0;
                    y206 = 1.0;
                } else if c.mth == 193 {
                    izr = izr - 4007;
                    y206 = 1.0;
                    y207 = 1.0;
                } else if c.mth == 194 {
                    izr = izr - 2006;
                    y203 = 2.0;
                } else if c.mth == 195 {
                    izr = izr - 4012;
                    y207 = 2.0;
                } else if c.mth == 196 {
                    izr = izr - 3009;
                    y203 = 1.0;
                    y207 = 1.0;
                } else if c.mth == 197 {
                    izr = izr - 3003;
                    y203 = 3.0;
                } else if c.mth == 198 {
                    izr = izr - 3004;
                    y203 = 3.0;
                } else if c.mth == 199 {
                    izr = izr - 4009;
                    y203 = 2.0;
                    y207 = 1.0;
                } else if c.mth == 200 {
                    izr = izr - 2007;
                    y203 = 2.0;
                }
                if izr == 1001 {
                    y203 = y203 + 1.0;
                }
                if izr == 1002 {
                    y204 = y204 + 1.0;
                }
                if izr == 1003 {
                    y205 = y205 + 1.0;
                }
                if izr == 2003 {
                    y206 = y206 + 1.0;
                }
                if izr == 2004 {
                    y207 = y207 + 1.0;
                }
                if izr == 4008 {
                    y207 = y207 + 2.0;
                }
                if y203 == ZERO && y204 == ZERO && y205 == ZERO && y206 == ZERO && y207 == ZERO {
                    break 'b310;
                }
                enext = thrg;
                if itemp == 1 {
                    let mut s = ffmt::i(c.mth, 8);
                    for yv in [y203, y204, y205, y206, y207] {
                        s.push_str("    ");
                        s.push_str(&ffmt::f(yv, 3, 1));
                    }
                    c.io.so(&s);
                }
                let mut i = 0;
                while i < ngas {
                    i += 1;
                    en = egas[i];
                    c.gety1(en, &mut enext, &mut idis, &mut y, nscr1, &mut a.all());
                    let mut ip = 2;
                    let mut ir = 1;
                    let ys = [y203, y204, y205, y206, y207];
                    let ls = [l203, l204, l205, l206, l207];
                    for g in 0..5 {
                        if ys[g] > HUND {
                            let f = six.at(ls[g]);
                            c.terpa(
                                &mut yyy, egas[i], &mut xnext, &mut idis, &f, &mut ip, &mut ir,
                            );
                            sgas[(g as i32 + 1, i)] = sgas[(g as i32 + 1, i)] + yyy * y;
                        } else {
                            sgas[(g as i32 + 1, i)] = sgas[(g as i32 + 1, i)] + ys[g] * y;
                        }
                    }
                }
            }
            // 310
            c.tosend(nscr1, 0, 0, &mut a.all());
        }
        // 330
        if itemp == 1 {
            c.io.so("\n *** means that the yield is energy dependent");
        }

        //--print out the gas-production results
        if iprint != 0 {
            if itemp == 1 {
                c.io.so("\n gas production versus energy");
                for i in 1..=ngas {
                    let mut s = ffmt::e1p(egas[i], 12, 4);
                    for j in 1..=5 {
                        s.push_str(&ffmt::e1p(sgas[(j, i)], 12, 4));
                    }
                    c.io.so(&s);
                }
            }
        }

        //--update the directory in file 1
        c.repoz(nscr1);
        let mut nsec = 0;
        let mut nn = [0i32; 6];
        for jg in 1..=5 {
            idone = 0;
            let mut i = 0;
            let mut ii = 0;
            while i < ngas && idone == 0 {
                i += 1;
                ii = i;
                if sgas[(jg, i)] != ZERO {
                    idone = 1;
                }
            }
            if idone == 1 {
                ii -= 1;
                nsec += 1;
                nn[jg as usize] = ngas - ii + 1;
            }
        }
        c.contio(nscr1, 0, 0, &mut a.all(), &mut nb, &mut nw);
        let mut k = 1 + nw;
        if c.iverf == 4 {
            nx = c.n2h;
        }
        if c.iverf >= 5 {
            c.contio(nscr1, 0, 0, &mut a.at(k), &mut nb, &mut nw);
            k += nw;
        }
        if c.iverf >= 6 {
            c.contio(nscr1, 0, 0, &mut a.at(k), &mut nb, &mut nw);
            k += nw;
        }
        c.hdatio(nscr1, 0, 0, &mut a.at(k), &mut nb, &mut nw);
        k += nw;
        if c.iverf != 4 {
            nx = c.n2h;
        }
        while nb != 0 {
            c.moreio(nscr1, 0, 0, &mut a.at(k), &mut nb, &mut nw);
            k += nw;
        }
        nw = nx;
        c.dictio(nscr1, 0, 0, &mut a.at(k), &mut nb, &mut nw);
        let mut nold = 0;
        for i in 1..=nx {
            let j = k - 1 + 6 * (i - 1);
            let mfi = nint(a[j + 3]);
            let mti = nint(a[j + 4]);
            if mfi == 3 && mti >= 203 && mti <= 207 {
                nold += 1;
            }
        }
        if nold > 0 {
            c.io.so("\n gas data on input pendf tape deleted");
        }
        let mut i = 0;
        let mut j = 0;
        idone = 0;
        while i < nx && idone == 0 {
            i += 1;
            j = k - 1 + 6 * (i - 1);
            if nint(a[j + 3]) > 3 {
                idone = 1;
            }
            if nint(a[j + 3]) == 3 && nint(a[j + 4]) > 200 {
                idone = 1;
            }
        }
        if idone == 0 {
            j = k - 1 + 6 * nx;
        } else {
            let ni = 6 * (nsec - nold);
            let nm = 6 * nx - j + k + 5;
            for i in 1..=nm {
                a[6 * nx + ni + k - i] = a[6 * nx + k - i];
            }
        }
        for (g, mt) in (1..=5).zip(203..=207) {
            let n = nn[g];
            if n != 0 {
                a[j + 1] = 0.0;
                a[j + 2] = 0.0;
                a[j + 3] = 3.0;
                a[j + 4] = mt as f64;
                a[j + 5] = ((n + 2) / 3) as f64;
                a[j + 6] = 1.0;
                j += 6;
            }
        }
        nw = 6;
        k = 1;
        if c.iverf == 4 {
            a[6] = (nx + nsec - nold) as f64;
        }
        c.contio(0, noutp, 0, &mut a.at(k), &mut nb, &mut nw);
        k += nw;
        if c.iverf >= 5 {
            c.contio(0, noutp, 0, &mut a.at(k), &mut nb, &mut nw);
            k += nw;
        }
        if c.iverf >= 6 {
            c.contio(0, noutp, 0, &mut a.at(k), &mut nb, &mut nw);
            k += nw;
        }
        if c.iverf != 4 {
            a[k + 5] = (nx + nsec - nold) as f64;
        }
        c.hdatio(0, noutp, 0, &mut a.at(k), &mut nb, &mut nw);
        k += nw;
        while nb != 0 {
            c.moreio(0, noutp, 0, &mut a.at(k), &mut nb, &mut nw);
            k += nw;
        }
        nw = nx + nsec - nold;
        c.dictio(0, noutp, 0, &mut a.at(k), &mut nb, &mut nw);
        c.tofend(nscr1, noutp, 0, &mut a.all());
        c.tofend(nscr1, noutp, 0, &mut a.all());

        //--copy file 3 down to the gas production area
        idone = 0;
        while idone == 0 {
            c.contio(nscr1, 0, 0, &mut b.all(), &mut nb, &mut nw);
            mfb = c.mfh;
            mtb = c.mth;
            if c.mth > 200 || c.mth == 0 {
                idone = 1;
            } else {
                c.contio(0, noutp, 0, &mut b.all(), &mut nb, &mut nw);
                c.tosend(nscr1, noutp, 0, &mut a.all());
            }
        }

        //--write the gas production sections in file 3
        for jg in 1..=5 {
            let mut i = 0;
            let mut ii = 0;
            let mut ifini = 0;
            while i < ngas && ifini == 0 {
                i += 1;
                ii = i;
                if sgas[(jg, i)] != ZERO {
                    ifini = 1;
                }
            }
            if ifini == 1 {
                i = ii - 1;
                if i == 0 {
                    i = 1;
                }
                c.mth = 202 + jg;
                c.math = matd;
                c.mfh = 3;
                a[1] = za;
                a[2] = awr;
                a[3] = 0.0;
                a[4] = 0.0;
                a[5] = 0.0;
                a[6] = 0.0;
                nw = 6;
                c.contio(0, noutp, 0, &mut a.all(), &mut nb, &mut nw);
                let np = ngas - i + 1;
                a[1] = 0.0;
                a[2] = 0.0;
                a[3] = 0.0;
                a[4] = 0.0;
                a[5] = 1.0;
                a[6] = np as f64;
                a[7] = np as f64;
                a[8] = 2.0;
                let mut k = 8;
                let mut istart = i;
                idone = 0;
                while idone == 0 {
                    let mut iend = ngas;
                    if (iend - istart) >= c.npage / 2 {
                        iend = istart + c.npage / 2 - 1;
                    }
                    let mut j = k - 1;
                    let mut ib = istart - 1;
                    while ib < iend {
                        j += 2;
                        ib += 1;
                        a[j] = egas[ib];
                        a[j + 1] = sigfig(sgas[(jg, ib)], 7, 0);
                    }
                    nw = j + 1;
                    if k != 0 {
                        k = 0;
                        c.tab1io(0, noutp, 0, &mut a.all(), &mut nb, &mut nw);
                    } else {
                        c.moreio(0, noutp, 0, &mut a.all(), &mut nb, &mut nw);
                    }
                    if nb == 0 {
                        idone = 1;
                    } else {
                        istart = iend + 1;
                    }
                }
                c.asend(noutp, 0);
            }
        }

        //--copy rest of this temperature to output file
        //--delete any gas sections on the input pend file
        c.mfh = mfb;
        c.mth = mtb;
        nw = 6;
        loop {
            // 520
            if c.mfh == 0 {
                break;
            }
            if c.mth >= 203 && c.mth <= 207 {
                // 525
                c.tosend(npend, 0, 0, &mut a.all());
                c.contio(npend, 0, 0, &mut b.all(), &mut nb, &mut nw);
                continue;
            }
            c.contio(0, noutp, 0, &mut b.all(), &mut nb, &mut nw);
            c.tosend(npend, noutp, 0, &mut a.all());
            c.contio(npend, 0, 0, &mut b.all(), &mut nb, &mut nw);
        }
        // 530
        c.contio(0, noutp, 0, &mut b.all(), &mut nb, &mut nw);
        c.tomend(npend, noutp, 0, &mut a.all());
    }
    // 550
    c.contio(0, noutp, 0, &mut b.all(), &mut nb, &mut nw);

    //--finished
    c.closz(nendf);
    c.closz(npend);
    c.closz(noutp);
    c.closz(nscr1);
    c.io.so(&format!("\n found{} temperatures", ffmt::i(itemp - 1, 2)));
    let time = c.timer();
    c.io.so(&format!(
        "{}{}s\n {}*******",
        " ".repeat(69),
        ffmt::f(time, 8, 1),
        "**********".repeat(7)
    ));
}
