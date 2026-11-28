//! GROUPR conver: copy the material to the scratch tapes, convert photon
//! transition probability arrays (lo=2) to yields, add mt456 if needed,
//! and build the lists used for automatic reaction processing.

use super::*;

/// Add mt to a run-length list (`mfX(i)` lists of conver): a value
/// -mt closes a run of consecutive mt numbers.
fn add_run(list: &mut [i32], idx: &mut i32, mth: i32) {
    let mut itest = 0;
    if *idx != 1 {
        if mth <= list[(*idx - 1) as usize].abs() + 1 {
            itest = 1;
        }
    }
    if itest == 1 {
        if list[(*idx - 1) as usize] < 0 {
            *idx -= 1;
        }
        list[*idx as usize] = -mth;
        *idx += 1;
    } else {
        list[*idx as usize] = mth;
        *idx += 1;
    }
}

/// Same for row k of a mf4r(6,*) list.
fn add_run2(list: &mut I2, k: i32, idx: &mut i32, mth: i32) {
    let mut itest = 0;
    if *idx != 1 {
        if mth <= list.g(k, *idx - 1).abs() + 1 {
            itest = 1;
        }
    }
    if itest == 1 {
        if list.g(k, *idx - 1) < 0 {
            *idx -= 1;
        }
        list.s(k, *idx, -mth);
        *idx += 1;
    } else {
        list.s(k, *idx, mth);
        *idx += 1;
    }
}

/// The File 6 list rule (labels 721-756): equal mt is skipped.
fn add_run6(list: &mut [i32], idx: &mut i32, mth: i32) {
    if *idx == 1 {
        list[*idx as usize] = mth;
        *idx += 1;
        return;
    }
    let prev = list[(*idx - 1) as usize];
    if mth == prev.abs() {
        return;
    }
    if mth > prev.abs() + 1 {
        list[*idx as usize] = mth;
        *idx += 1;
        return;
    }
    if prev < 0 {
        *idx -= 1;
    }
    list[*idx as usize] = -mth;
    *idx += 1;
}

fn add_run6b(list: &mut I2, k: i32, idx: &mut i32, mth: i32, skip18: bool) {
    if *idx == 1 {
        list.s(k, *idx, mth);
        *idx += 1;
        return;
    }
    if skip18 && mth == 18 {
        return;
    }
    let prev = list.g(k, *idx - 1);
    if mth == prev.abs() {
        return;
    }
    if mth > prev.abs() + 1 {
        list.s(k, *idx, mth);
        *idx += 1;
        return;
    }
    if prev < 0 {
        *idx -= 1;
    }
    list.s(k, *idx, -mth);
    *idx += 1;
}

impl Groupr {
    pub(super) fn conver(&mut self, c: &mut Ctx, nin: i32, nout: i32, nscr: i32) {
        let zero = 0.0f64;
        let mxnnth = 450;
        let (mut nb, mut nw) = (0i32, 0i32);
        let mut l2flg = 0;
        self.ndelg = 0;
        let imax: i32 = 50;
        let lmax: i32 = 500;
        let mut ngam = vec![0i32; 441];
        let mut istor = 0;
        let mut scr = RV::new((c.npage + 50) as usize + SLACK);
        let mut flsx = RV::new((c.npage + 50) as usize + SLACK);
        c.nsh = 0;
        c.nsc = 0;
        let mut imf4 = 1;
        let mut imf6 = 1;
        let mut imf10 = 1;
        let mut imf12 = 1;
        let mut imf13 = 1;
        let mut imf18 = 1;
        let mut imf2x = [1i32; 7];
        let mut irr2x = [1i32; 7];
        let mut mt0 = 0;
        let mut mt0old = 0;
        let mut mtnow = -1;
        let mut nnth = 0;
        self.izatest = 0;
        self.lastza = 0;
        let mut mtth = vec![0i32; mxnnth + 1];
        let mut eeth = vec![0.0f64; mxnnth + 1];
        let mut e = RV::new(0);
        let mut eg = RV::new(0);
        let mut es = RV::new(0);
        let mut y = RV::new(0);
        let mut aa = RV::new(0);
        let mut r = RV::new(0);
        let mut sig = RM::new(1, 1);
        let mut fl = [0.0f64; 66];

        let mut lbl = 110;
        loop {
            match lbl {
                110 => {
                    //--loop over all sections of this material
                    c.contio(nin, 0, 0, &mut scr.at(1), &mut nb, &mut nw);
                    if c.mfh == 0 && c.math != 0 {
                        lbl = 130;
                        continue;
                    }
                    if c.mfh == 1 {
                        lbl = 550;
                        continue;
                    }

                    //--get thresholds vs mt number
                    if c.mfh == 3 && c.mth != 0 {
                        let ei = 0.0f64;
                        let (mut enxt, mut jdis, mut x) = (0.0f64, 0i32, 0.0f64);
                        c.gety1(ei, &mut enxt, &mut jdis, &mut x, nin, &mut scr.at(1));
                        nnth += 1;
                        if nnth > mxnnth {
                            c.error("conver", "nnth too large", " ");
                        }
                        eeth[nnth] = enxt;
                        mtth[nnth] = c.mth;
                        c.contio(0, nout, nscr, &mut scr.at(1), &mut nb, &mut nw);
                        c.tosend(nin, nout, nscr, &mut scr.at(1));
                        lbl = 110;
                        continue;
                    }
                    if c.mfh == 12 && c.mth == 460 {
                        c.tosend(nin, 0, 0, &mut scr.at(1));
                        lbl = 110;
                        continue;
                    }
                    if c.mfh == 6 && c.mth == 18 && c.l1h != 0 {
                        c.mess("conver", "skipping new mf6/mt18 multiplicity section", "");
                        c.tosend(nin, 0, 0, &mut scr.at(1));
                        lbl = 110;
                        continue;
                    }
                    let mth = c.mth;
                    if c.mfh == 4 {
                        // mf4 list
                        let mut to211 = imf4 == 1 || mth == 18;
                        let mut to311 = false;
                        if !to211 {
                            if mth == 19 && self.mf4[(imf4 - 1) as usize] == 18 {
                                to311 = true;
                            } else if mth > self.mf4[(imf4 - 1) as usize].abs() + 1 {
                                to211 = true;
                            }
                        }
                        if to311 {
                            imf4 -= 1;
                            to211 = true;
                        }
                        if to211 {
                            self.mf4[imf4 as usize] = mth;
                            imf4 += 1;
                        } else {
                            if self.mf4[(imf4 - 1) as usize] < 0 {
                                imf4 -= 1;
                            }
                            self.mf4[imf4 as usize] = -mth;
                            imf4 += 1;
                        }
                        // 312
                        let iverf = c.iverf;
                        let mut za2 = -1.0f64;
                        if mth == 2 || (mth >= 50 && mth < 91) {
                            za2 = 1.0;
                        }
                        if iverf < 6 {
                            if mth >= 700 && mth < 719 {
                                za2 = 1001.0;
                            }
                            if mth >= 720 && mth < 739 {
                                za2 = 1002.0;
                            }
                            if mth >= 740 && mth < 759 {
                                za2 = 1003.0;
                            }
                            if mth >= 760 && mth < 779 {
                                za2 = 2003.0;
                            }
                            if mth >= 780 && mth < 799 {
                                za2 = 2004.0;
                            }
                        } else {
                            if mth >= 600 && mth < 649 {
                                za2 = 1001.0;
                            }
                            if mth >= 650 && mth < 699 {
                                za2 = 1002.0;
                            }
                            if mth >= 700 && mth < 749 {
                                za2 = 1003.0;
                            }
                            if mth >= 750 && mth < 799 {
                                za2 = 2003.0;
                            }
                            if mth >= 800 && mth < 849 {
                                za2 = 2004.0;
                            }
                            if mth >= 875 && mth < 891 {
                                za2 = 1.0;
                            }
                        }
                        if za2 < 0.0 && mth != 18 && mth != 19 && mth != 20 && mth != 21 && mth != 38 {
                            let s = format!("cannot do complete particle production for mt={}", ffmt::i(mth, 3));
                            c.mess("conver", &s, "only mf4/mf5 provided");
                            lbl = 119;
                            continue;
                        }
                        self.izat = nint(c.c1h);
                        let jza2 = nint(za2);
                        let jzar = self.izat + self.izap - jza2;
                        let k = match jza2 {
                            1001 => 1,
                            1002 => 2,
                            1003 => 3,
                            2003 => 4,
                            2004 => 5,
                            _ => 0,
                        };
                        if k > 0 {
                            add_run2(&mut self.mf4r, k, &mut irr2x[k as usize], mth);
                        }
                        if jzar > 2004 {
                            let idx = &mut irr2x[6];
                            let mut itest = 0;
                            if *idx != 1 {
                                if mth <= self.mf4r.g(6, *idx - 1).abs() + 1 {
                                    itest = 1;
                                }
                            }
                            if mth != 18 {
                                // exclude fission for residual production
                                if itest == 1 {
                                    if self.mf4r.g(6, *idx - 1) < 0 {
                                        *idx -= 1;
                                    }
                                    self.mf4r.s(6, *idx, -mth);
                                    *idx += 1;
                                } else {
                                    self.mf4r.s(6, *idx, mth);
                                    *idx += 1;
                                }
                            }
                        }
                        lbl = 119;
                        continue;
                    }
                    // 111
                    if c.mfh == 6 {
                        lbl = 620;
                        continue;
                    }
                    if c.mfh == 8 {
                        lbl = 820;
                        continue;
                    }
                    if c.mfh == 12 {
                        add_run(&mut self.mf12, &mut imf12, mth);
                    } else if c.mfh == 13 {
                        add_run(&mut self.mf13, &mut imf13, mth);
                    }
                    lbl = 119;
                }
                119 => {
                    if c.mfh == 12 {
                        // 120
                        if c.l1h == 2 {
                            lbl = 150;
                            continue;
                        }
                        c.contio(0, nout, nscr, &mut scr.at(1), &mut nb, &mut nw);
                        c.tosend(nin, nout, nscr, &mut scr.at(1));
                        lbl = 110;
                        continue;
                    }
                    let (math, mfh, mth, iverf) = (c.math, c.mfh, c.mth, c.iverf);
                    if (math == 1149 && mfh == 13 && mth >= 51 && iverf < 5)
                        || (math == 1150 && mfh == 13 && mth <= 54 && iverf < 5)
                    {
                        // 140
                        c.tosend(nin, 0, 0, &mut scr.at(1));
                        let s = format!("gamma production patch made for {}", ffmt::i(mth, 3));
                        c.mess("conver", &s, " ");
                        lbl = 110;
                        continue;
                    }
                    if mfh == 14 && l2flg == 1 {
                        lbl = 500;
                        continue;
                    }
                    c.contio(0, nout, nscr, &mut scr.at(1), &mut nb, &mut nw);
                    if c.math == 0 {
                        lbl = 600;
                        continue;
                    }
                    c.tosend(nin, nout, nscr, &mut scr.at(1));
                    lbl = 110;
                }
                130 => {
                    c.contio(0, nout, nscr, &mut scr.at(1), &mut nb, &mut nw);
                    lbl = 110;
                }
                150 => {
                    //--allocate storage
                    if istor == 0 {
                        istor = 1;
                        let imax2 = imax * imax;
                        e = RV::new(imax as usize + 1);
                        eg = RV::new(lmax as usize + 1);
                        es = RV::new(lmax as usize + 1);
                        y = RV::new(lmax as usize + 1);
                        aa = RV::new(imax2 as usize + 1);
                        r = RV::new(imax2 as usize + 1);
                        e[1] = 0.0;
                        for i in 1..=imax2 {
                            r[i] = 0.0;
                            aa[i] = 0.0;
                        }
                        for i in 1..=imax {
                            r[imax * (i - 1) + i] = 1.0;
                        }
                    }
                    let mth = c.mth;
                    let iverf = c.iverf;

                    //--check for missing mf12 mt's ... advise user if found
                    let v6 = iverf >= 6;
                    let v5 = iverf <= 5;
                    if mth == 51
                        || (v6 && (mth == 601 || mth == 651 || mth == 701 || mth == 751 || mth == 801 || mth == 876))
                        || (v5 && (mth == 701 || mth == 721 || mth == 741 || mth == 761 || mth == 781))
                    {
                        mtnow = mth;
                    } else if (mth > 51 && mth < 91)
                        || (v6 && mth > 601 && mth < 648)
                        || (v6 && mth > 651 && mth < 698)
                        || (v6 && mth > 701 && mth < 748)
                        || (v6 && mth > 751 && mth < 798)
                        || (v6 && mth > 801 && mth < 848)
                        || (v6 && mth > 876 && mth < 890)
                        || (v5 && mth > 701 && mth < 718)
                        || (v5 && mth > 721 && mth < 738)
                        || (v5 && mth > 741 && mth < 758)
                        || (v5 && mth > 761 && mth < 778)
                        || (v5 && mth > 781 && mth < 798)
                    {
                        if mtnow + 1 != mth {
                            let s = format!("mf12, mt{} may be missing", ffmt::i(mth - 1, 2));
                            c.mess("conver", &s, "discrete photon data may be incomplete");
                        }
                        mtnow = mth;
                    }

                    //--convert transition probability array to yields
                    let za = c.c1h;
                    self.awr = c.c2h;
                    let awr = self.awr;
                    let lg = c.l2h;
                    l2flg = 1;
                    let ll = 1;
                    c.listio(nin, 0, 0, &mut scr.at(ll), &mut nb, &mut nw);
                    let mut ll = ll;
                    while nb != 0 {
                        ll += nw;
                        c.moreio(nin, 0, 0, &mut flsx.at(ll), &mut nb, &mut nw);
                    }

                    //--make sure mt0 is correct for this range of mt's.
                    if mth >= 51 && mth <= 90 && mt0 != 49 {
                        mt0 = 49;
                    }
                    if v6 {
                        if mth >= 601 && mth <= 649 && mt0 != 599 {
                            mt0 = 599;
                        }
                        if mth >= 651 && mth <= 699 && mt0 != 649 {
                            mt0 = 649;
                        }
                        if mth >= 701 && mth <= 749 && mt0 != 699 {
                            mt0 = 699;
                        }
                        if mth >= 751 && mth <= 799 && mt0 != 749 {
                            mt0 = 749;
                        }
                        if mth >= 801 && mth <= 849 && mt0 != 799 {
                            mt0 = 799;
                        }
                        if mth >= 876 && mth <= 891 && mt0 != 874 {
                            mt0 = 874;
                        }
                    } else if v5 {
                        if mth >= 701 && mth <= 719 && mt0 != 699 {
                            mt0 = 699;
                        }
                        if mth >= 721 && mth <= 739 && mt0 != 719 {
                            mt0 = 719;
                        }
                        if mth >= 741 && mth <= 759 && mt0 != 739 {
                            mt0 = 739;
                        }
                        if mth >= 761 && mth <= 779 && mt0 != 759 {
                            mt0 = 759;
                        }
                        if mth >= 781 && mth <= 799 && mt0 != 779 {
                            mt0 = 779;
                        }
                    }

                    //--load the e(_) array with mf3 -q.
                    if mt0 != mt0old {
                        for i in 1..=imax {
                            e[i] = 0.0;
                        }
                        mt0old = mt0;
                        let m1 = mt0 + 2;
                        let mut m2 = 0;
                        if mt0 == 49 {
                            m2 = 91;
                        } else if v6 && mt0 != 874 {
                            m2 = m1 + 48;
                        } else if v6 && mt0 == 874 {
                            m2 = m1 + 15;
                        } else if v5 {
                            m2 = m1 + 17;
                        }
                        let mut nn = 1;
                        let mut mttst = -1;
                        while mttst < m2 && nn <= nnth {
                            mttst = mtth[nn as usize];
                            if mttst >= m1 && mttst <= m2 {
                                e[mttst - mt0] = eeth[nn as usize] * awr / (awr + 1.0);
                            }
                            nn += 1;
                        }
                    }

                    let j = mth - mt0;
                    e[j] = scr[1];
                    let n = c.n2h;
                    let jm1 = j - 1;
                    for kk in 1..=jm1 {
                        let k = jm1 - kk + 1;
                        let mut ii = 0;
                        let mut idone = 0;
                        let mut i = 0;
                        while ii < n && idone == 0 {
                            ii += 1;
                            i = ii;
                            let ei = scr[6 + (lg + 1) * i - lg];
                            if ei == zero && e[k] == zero {
                                idone = 1;
                            }
                            if ei != zero {
                                if (ei - e[k]).abs() / ei < 0.0001 {
                                    idone = 1;
                                }
                            }
                        }
                        if idone == 0 {
                            aa[(k - 1) * imax + j] = 0.0;
                            r[(k - 1) * imax + j] = 0.0;
                        } else {
                            let p = scr[7 + (lg + 1) * i - lg];
                            let mut g = 1.0;
                            if lg == 2 {
                                g = scr[6 + (lg + 1) * i];
                            }
                            aa[(k - 1) * imax + j] = p * g;
                            r[(k - 1) * imax + j] = p;
                        }
                        if k != jm1 {
                            let kp1 = k + 1;
                            for ii in kp1..=jm1 {
                                r[(k - 1) * imax + j] = r[(k - 1) * imax + j] + r[(ii - 1) * imax + j] * aa[(k - 1) * imax + ii];
                            }
                        }
                    }
                    let mut l = 0;
                    let mut ysum = 0.0f64;
                    for i in 2..=j {
                        let ja = j + 2 - i;
                        let eja = e[ja];
                        for ii in 1..=jm1 {
                            let jb = j - ii;
                            let ejb = e[jb];
                            let yy = aa[(jb - 1) * imax + ja] * r[(ja - 1) * imax + j];
                            if yy != zero {
                                l += 1;
                                if l > lmax {
                                    c.error("conver", "too many lo=2 gammas", " ");
                                }
                                eg[l] = eja - ejb;
                                es[l] = eja;
                                y[l] = yy;
                                ysum = ysum + y[l];
                            }
                        }
                    }
                    e[1] = 0.0;

                    //--arrange gamma energies in descending order
                    if l > 1 {
                        let lm1 = l - 1;
                        for i in 1..=lm1 {
                            let ip1 = i + 1;
                            for ii in ip1..=l {
                                if eg[i] < eg[ii] {
                                    let t = eg[i];
                                    eg[i] = eg[ii];
                                    eg[ii] = t;
                                    let t = y[i];
                                    y[i] = y[ii];
                                    y[ii] = t;
                                    let t = es[i];
                                    es[i] = es[ii];
                                    es[ii] = t;
                                }
                            }
                        }
                    }

                    //--output computed yields in endf lo=1 format
                    scr[1] = za;
                    scr[2] = awr;
                    scr[3] = 1.0;
                    scr[4] = 0.0;
                    scr[5] = l as f64;
                    scr[6] = 0.0;
                    let mut mtl = 49;
                    if v5 && mth >= 700 {
                        mtl = 649;
                    }
                    if v6 && mth >= 600 {
                        mtl = 549;
                    }
                    ngam[(mth - mtl) as usize] = l;
                    c.contio(0, nout, nscr, &mut scr.at(1), &mut nb, &mut nw);
                    let etop = self.emaxx;
                    let mut elow = self.ebeg;
                    for i in 1..=nnth {
                        if mtth[i as usize] == mth {
                            elow = eeth[i as usize];
                        }
                    }
                    // output tab1 sum record, if needed
                    if l > 1 {
                        scr[1] = 0.0;
                        scr[2] = 0.0;
                        scr[3] = 0.0;
                        scr[5] = 1.0;
                        scr[6] = 2.0;
                        scr[7] = 2.0;
                        scr[8] = 2.0;
                        scr[9] = elow;
                        scr[10] = ysum;
                        scr[11] = etop;
                        scr[12] = ysum;
                        nw = 12;
                        c.tab1io(0, nout, nscr, &mut scr.at(1), &mut nb, &mut nw);
                    }
                    // output rest of tab1 records
                    for i in 1..=l {
                        scr[1] = eg[i];
                        scr[2] = es[i];
                        scr[3] = 0.0;
                        scr[4] = 2.0;
                        scr[5] = 1.0;
                        scr[6] = 2.0;
                        scr[7] = 2.0;
                        scr[8] = 2.0;
                        scr[9] = elow;
                        scr[10] = y[i];
                        scr[11] = etop;
                        scr[12] = y[i];
                        nw = 12;
                        c.tab1io(0, nout, nscr, &mut scr.at(1), &mut nb, &mut nw);
                    }

                    //--copy send record and loop over remaining sections.
                    c.contio(nin, nout, nscr, &mut scr.at(1), &mut nb, &mut nw);
                    lbl = 110;
                }
                500 => {
                    //--angular distributions
                    let mth = c.mth;
                    let iverf = c.iverf;
                    let v5 = iverf <= 5;
                    let v6 = iverf >= 6;
                    let mut mt0x = 0;
                    if mth >= 50 && mth <= 90 {
                        mt0x = 50;
                    }
                    if v5 && mth >= 700 && mth <= 717 {
                        mt0x = 700;
                    }
                    if v5 && mth >= 720 && mth <= 737 {
                        mt0x = 720;
                    }
                    if v5 && mth >= 740 && mth <= 757 {
                        mt0x = 740;
                    }
                    if v5 && mth >= 760 && mth <= 777 {
                        mt0x = 760;
                    }
                    if v5 && mth >= 780 && mth <= 797 {
                        mt0x = 780;
                    }
                    if v6 && mth >= 600 && mth <= 648 {
                        mt0x = 600;
                    }
                    if v6 && mth >= 650 && mth <= 698 {
                        mt0x = 650;
                    }
                    if v6 && mth >= 700 && mth <= 748 {
                        mt0x = 700;
                    }
                    if v6 && mth >= 750 && mth <= 798 {
                        mt0x = 750;
                    }
                    if v6 && mth >= 800 && mth <= 848 {
                        mt0x = 800;
                    }
                    if v6 && mth >= 875 && mth <= 890 {
                        mt0x = 875;
                    }
                    mt0 = mt0x;
                    if mt0 != 0 {
                        let mut mtl = 49;
                        if v5 && mth >= 700 {
                            mtl = 649;
                        }
                        if v6 && mth >= 600 {
                            mtl = 549;
                        }
                        scr[3] = 1.0;
                        scr[4] = 0.0;
                        scr[6] = 0.0;
                        scr[5] = ngam[(mth - mtl) as usize] as f64;
                        c.contio(0, nout, nscr, &mut scr.at(1), &mut nb, &mut nw);
                        c.tosend(nin, 0, 0, &mut scr.at(1));
                        c.asend(nout, nscr);
                    } else {
                        c.contio(0, nout, nscr, &mut scr.at(1), &mut nb, &mut nw);
                        c.tosend(nin, nout, nscr, &mut scr.at(1));
                    }
                    lbl = 110;
                }
                550 => {
                    //--file1.  add 456=452 if necessary.
                    c.contio(0, nout, nscr, &mut scr.at(1), &mut nb, &mut nw);
                    self.izap = 1;
                    self.awrp = 1.0;
                    if c.iverf >= 6 {
                        c.contio(nin, nout, nscr, &mut scr.at(1), &mut nb, &mut nw);
                        c.contio(nin, nout, nscr, &mut scr.at(1), &mut nb, &mut nw);
                        self.izap = c.n1h / 10;
                        self.awrp = c.c1h;
                    }
                    c.tosend(nin, nout, nscr, &mut scr.at(1));
                    let mut no455 = 0;
                    c.contio(nin, nout, nscr, &mut scr.at(1), &mut nb, &mut nw);
                    if c.mfh == 0 {
                        lbl = 110;
                        continue;
                    }
                    if c.mth != 452 {
                        // 595
                        c.tofend(nin, nout, nscr, &mut scr.at(1));
                        lbl = 110;
                        continue;
                    }
                    let nnu = 20000;
                    let mut nu = RV::new(nnu as usize + SLACK);
                    let mut l = 1;
                    let lnu = c.l2h;
                    if lnu == 1 {
                        c.listio(nin, nout, nscr, &mut nu.at(l), &mut nb, &mut nw);
                    }
                    if lnu == 2 {
                        c.tab1io(nin, nout, nscr, &mut nu.at(l), &mut nb, &mut nw);
                    }
                    while nb != 0 {
                        if l + nw > nnu {
                            c.error(
                                "conver",
                                "storage for fission nu exceeded (increase the size of the nu array)",
                                " ",
                            );
                        }
                        l += nw;
                        c.moreio(nin, nout, nscr, &mut nu.at(l), &mut nb, &mut nw);
                    }
                    c.tosend(nin, nout, nscr, &mut scr.at(1));
                    // 565
                    let mut go590 = false;
                    loop {
                        c.contio(nin, 0, 0, &mut scr.at(1), &mut nb, &mut nw);
                        if c.mth == 456 {
                            go590 = true;
                            break;
                        }
                        if c.mth == 0 && no455 == 0 {
                            break;
                        }
                        if c.mfh == 0 {
                            go590 = true;
                            break;
                        }
                        if c.mth == 455 {
                            no455 = 1;
                        }
                        if c.mth > 456 && no455 == 0 {
                            break;
                        }
                        c.contio(0, nout, nscr, &mut scr.at(1), &mut nb, &mut nw);
                        if c.mth == 455 {
                            c.listio(nin, nout, nscr, &mut scr.at(1), &mut nb, &mut nw);
                            self.ndelg = c.n1h;
                        }
                        c.tosend(nin, nout, nscr, &mut scr.at(1));
                    }
                    if !go590 {
                        // 570
                        if no455 != 1 {
                            let mf = c.mfh;
                            let mt = c.mth;
                            let l1 = c.l1h;
                            let l2 = c.l2h;
                            let n1 = c.n1h;
                            let n2 = c.n2h;
                            c.mfh = 1;
                            c.mth = 456;
                            scr[3] = 0.0;
                            scr[4] = lnu as f64;
                            scr[5] = 0.0;
                            scr[6] = 0.0;
                            l = 1;
                            c.contio(0, nout, nscr, &mut scr.at(1), &mut nb, &mut nw);
                            if lnu == 1 {
                                c.listio(0, nout, nscr, &mut nu.at(l), &mut nb, &mut nw);
                            }
                            if lnu == 2 {
                                c.tab1io(0, nout, nscr, &mut nu.at(l), &mut nb, &mut nw);
                            }
                            while nb != 0 {
                                l += nw;
                                c.moreio(0, nout, nscr, &mut nu.at(l), &mut nb, &mut nw);
                            }
                            c.asend(nout, nscr);
                            c.mfh = mf;
                            c.mth = mt;
                            scr[3] = l1 as f64;
                            scr[4] = l2 as f64;
                            scr[5] = n1 as f64;
                            scr[6] = n2 as f64;
                        }
                    }
                    // 590
                    c.contio(0, nout, nscr, &mut scr.at(1), &mut nb, &mut nw);
                    if c.mfh == 0 {
                        lbl = 110;
                        continue;
                    }
                    // 595
                    c.tofend(nin, nout, nscr, &mut scr.at(1));
                    lbl = 110;
                }
                620 => {
                    //--examine contents of file 6.
                    c.contio(0, nout, nscr, &mut scr.at(1), &mut nb, &mut nw);
                    let nk = c.n1h;
                    let mut ik = 0;
                    loop {
                        // 625
                        ik += 1;
                        c.tab1io(nin, nout, nscr, &mut scr.at(1), &mut nb, &mut nw);
                        while nb != 0 {
                            c.moreio(nin, nout, nscr, &mut scr.at(1), &mut nb, &mut nw);
                        }
                        let jzap = nint(c.c1h);
                        self.law = c.l2h;

                        if self.law == 5 {
                            //--convert cp elastic ni format
                            self.izap = -self.izap;
                            self.matd = c.math;
                            self.mfd = c.mfh;
                            self.mtd = c.mth;
                            let mut e1 = 0.0f64;
                            let (mut enext, mut idis) = (0.0f64, 0i32);
                            self.getsig(c, e1, &mut enext, &mut idis, &mut sig, 1, 1);
                            c.tab2io(nin, nout, nscr, &mut scr.at(1), &mut nb, &mut nw);
                            self.spi = c.c1h;
                            let ne = c.n2h;
                            for _ie in 1..=ne {
                                let mut ll = 1;
                                c.listio(nin, 0, 0, &mut scr.at(ll), &mut nb, &mut nw);
                                let ltp = c.l1h;
                                while nb != 0 {
                                    ll += nw;
                                    c.moreio(nin, 0, 0, &mut scr.at(ll), &mut nb, &mut nw);
                                }
                                if ltp <= 2 {
                                    ll = 1;
                                    c.listio(0, nout, nscr, &mut scr.at(ll), &mut nb, &mut nw);
                                    while nb != 0 {
                                        ll += nw;
                                        c.moreio(0, nout, nscr, &mut scr.at(ll), &mut nb, &mut nw);
                                    }
                                } else {
                                    e1 = c.c2h;
                                    self.getsig(c, e1, &mut enext, &mut idis, &mut sig, 1, 1);
                                    let nm = c.n2h;
                                    for im in 1..=nm {
                                        let w = scr[5 + 2 * im];
                                        scr[6 + 2 * im] = sig[(1, 1)] * (1.0 - w) * scr[6 + 2 * im];
                                    }
                                    let ltt = 2;
                                    let lcd = 2;
                                    self.lct = 2;
                                    let mut nl = 21;
                                    let cp = dist::CoPar { awr: self.awr, awrp: self.awrp, izap: self.izap, izat: self.izat };
                                    dist::getco(c, &cp, &mut fl[1..], &mut nl, lcd, &mut scr.0[..], self.lct, ltt, &mut idis);
                                    scr[3] = 2.0;
                                    scr[5] = nl as f64;
                                    scr[6] = nl as f64;
                                    for il in 1..=nl {
                                        scr[6 + il] = fl[il as usize];
                                    }
                                    ll = 1;
                                    c.listio(0, nout, nscr, &mut scr.at(ll), &mut nb, &mut nw);
                                    while nb != 0 {
                                        ll += nw;
                                        c.moreio(0, nout, nscr, &mut scr.at(ll), &mut nb, &mut nw);
                                    }
                                }
                            }
                            self.izap = -self.izap;
                        } else {
                            //--just copy this subsection
                            c.skip6(nin, nout, nscr, &mut scr.at(1), self.law);
                        }

                        //--record reaction information
                        let mth = c.mth;
                        match jzap {
                            0 => add_run6(&mut self.mf18, &mut imf18, mth),
                            1 => add_run6(&mut self.mf6, &mut imf6, mth),
                            1001 => add_run6b(&mut self.mf6p, 1, &mut imf2x[1], mth, false),
                            1002 => add_run6b(&mut self.mf6p, 2, &mut imf2x[2], mth, false),
                            1003 => add_run6b(&mut self.mf6p, 3, &mut imf2x[3], mth, false),
                            2003 => add_run6b(&mut self.mf6p, 4, &mut imf2x[4], mth, false),
                            2004 => add_run6b(&mut self.mf6p, 5, &mut imf2x[5], mth, false),
                            _ => add_run6b(&mut self.mf6p, 6, &mut imf2x[6], mth, true),
                        }
                        // 790
                        if ik < nk {
                            continue;
                        }
                        c.tosend(nin, nout, nscr, &mut scr.at(1));
                        break;
                    }
                    lbl = 110;
                }
                820 => {
                    //--examine contents of file 8
                    let nk = c.n1h;
                    let mth = c.mth;
                    if mth <= 200 || mth >= 600 {
                        let mut ik = 0;
                        while ik < nk {
                            ik += 1;
                            let mut ll = 1;
                            c.listio(nin, 0, 0, &mut scr.at(ll), &mut nb, &mut nw);
                            while nb != 0 {
                                ll += nw;
                                c.moreio(nin, 0, 0, &mut scr.at(ll), &mut nb, &mut nw);
                            }
                            let izan = nint(c.c1h);
                            let imf = c.l1h;
                            let iis = c.l2h;
                            let k = imf10 as usize;
                            self.mf10f[k] = imf;
                            self.mf10s[k] = mth;
                            self.mf10i[k] = 10 * izan;
                            if izan == -1 && mth == 18 {
                                self.mf10i[k] = 0;
                            }
                            self.lfs8[k] = iis;
                            if izan != self.lastza {
                                if iis == 0 {
                                    self.mlfs8[k] = 0;
                                } else {
                                    self.mlfs8[k] = 1;
                                }
                            } else {
                                self.mlfs8[k] = self.mlfs8[k - 1] + 1;
                            }
                            self.lastza = izan;
                            imf10 += 1;
                        }
                    }
                    c.tosend(nin, 0, 0, &mut scr.at(1));
                    lbl = 110;
                }
                _ => {
                    // 600: conver is finished.
                    let maxr1 = MAXR1 as i32;
                    let chk: [(i32, &str); 20] = [
                        (imf4, "imf4 too big"),
                        (irr2x[1], "irr21 too big"),
                        (irr2x[2], "irr22 too big"),
                        (irr2x[3], "irr23 too big"),
                        (irr2x[4], "irr24 too big"),
                        (irr2x[4], "irr24 too big"),
                        (irr2x[5], "irr25 too big"),
                        (irr2x[6], "irr26 too big"),
                        (imf12, "imf12 too big"),
                        (imf13, "imf13 too big"),
                        (imf18, "imf18 too big"),
                        (imf6, "imf6 too big"),
                        (imf2x[1], "imf21 too big"),
                        (imf2x[2], "imf22 too big"),
                        (imf2x[3], "imf23 too big"),
                        (imf2x[4], "imf24 too big"),
                        (imf2x[5], "imf25 too big"),
                        (imf2x[6], "imf26 too big"),
                        (imf2x[6], "imf26 too big"),
                        (0, ""),
                    ];
                    for (v, m) in chk.iter() {
                        if *v > maxr1 {
                            c.error("conver", m, " ");
                        }
                    }
                    if imf10 > MAXR2 as i32 {
                        c.error("conver", "imf10 too big", " ");
                    }
                    self.mf4[imf4 as usize] = 0;
                    self.mf6[imf6 as usize] = 0;
                    let k = imf10 as usize;
                    self.mf10f[k] = 0;
                    self.mf10s[k] = 0;
                    self.mf10i[k] = 0;
                    self.lfs8[k] = 0;
                    self.mf12[imf12 as usize] = 0;
                    self.mf13[imf13 as usize] = 0;
                    self.mf18[imf18 as usize] = 0;
                    for kk in 1..=6 {
                        self.mf6p.s(kk, imf2x[kk as usize], 0);
                        self.mf4r.s(kk, irr2x[kk as usize], 0);
                    }
                    return;
                }
            }
        }
    }
}
