//! Compare the allocation-free Fortran editing routines in `ffmt` with the
//! reference versions they replaced (which were checked against gfortran).
//! Prints the number of cases and exits with status 1 on the first mismatch.
//! Usage: ffmt_check [N_RANDOM]
use rustjoy::ffmt;

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        // xorshift64*
        self.0 ^= self.0 >> 12;
        self.0 ^= self.0 << 25;
        self.0 ^= self.0 >> 27;
        self.0.wrapping_mul(0x2545F4914F6CDD1D)
    }
    fn below(&mut self, n: u64) -> u64 {
        self.next() % n
    }
}

fn value(r: &mut Rng) -> f64 {
    match r.below(6) {
        // any bit pattern (NaN, inf, subnormals included)
        0 => f64::from_bits(r.next()),
        // short decimals, many exact ties at the rounding digit
        1 => {
            let m = (r.below(2_000_001) as i64 - 1_000_000) as f64;
            let e = r.below(25) as i32 - 12;
            m * 10f64.powi(e) / 8.0
        }
        // ordinary magnitudes
        2 => {
            let m = (r.next() >> 11) as f64 / (1u64 << 53) as f64;
            let e = r.below(40) as i32 - 20;
            let s = if r.below(2) == 0 { 1.0 } else { -1.0 };
            s * m * 10f64.powi(e)
        }
        // values that round across a power of ten
        3 => {
            let e = r.below(30) as i32 - 15;
            let p = 10f64.powi(e);
            let t = (r.below(2001) as f64 - 1000.0) * 1e-12;
            (1.0 - 0.5e-9 + t) * p
        }
        4 => [0.0, -0.0, 0.5, -0.5, 9.5, 99.95, 0.05, 1e-300, 1e300, f64::MAX, f64::MIN_POSITIVE][r.below(11) as usize],
        _ => (r.below(200001) as f64 - 100000.0) / 1000.0,
    }
}

fn main() {
    let n: u64 = std::env::args().nth(1).map(|s| s.parse().unwrap()).unwrap_or(5_000_000);
    let mut r = Rng(0x9E3779B97F4A7C15);
    let mut count = 0u64;
    let fail = |what: &str, a: &str, b: &str| -> ! {
        println!("MISMATCH {}: new '{}' ref '{}'", what, a, b);
        std::process::exit(1);
    };
    for _ in 0..n {
        let x = value(&mut r);
        let w = 1 + r.below(30) as usize;
        let d = r.below(16) as usize;
        let k = r.below(9) as i32 - 3;
        let a = ffmt::fp(x, w, d, k);
        let b = ffmt::fp_ref(x, w, d, k);
        if a != b {
            fail(&format!("fp({:e},{},{},{})", x, w, d, k), &a, &b);
        }
        let ke = r.below(4) as i32;
        if ke as usize <= d + 1 && x.is_finite() {
            let a = ffmt::ep(x, w, d, ke);
            let b = ffmt::ep_ref(x, w, d, ke);
            if a != b {
                fail(&format!("ep({:e},{},{},{})", x, w, d, ke), &a, &b);
            }
        }
        let v = match r.below(3) {
            0 => r.next() as i32,
            1 => r.below(200001) as i32 - 100000,
            _ => [0, -1, i32::MIN, i32::MAX, 9, 10, -9, -10][r.below(8) as usize],
        };
        let wi = 1 + r.below(14) as usize;
        let a = ffmt::i(v, wi);
        let b = ffmt::i_ref(v, wi);
        if a != b {
            fail(&format!("i({},{})", v, wi), &a, &b);
        }
        count += 3;
    }
    // the formats NJOY uses, over a dense set of values
    for j in 0..2_000_000u64 {
        let x = value(&mut r);
        for &(w, d, k) in &[(11usize, 8usize, 0i32), (11, 7, 1), (11, 6, 2), (11, 5, 3), (11, 4, 4), (11, 3, 5), (8, 3, 1), (10, 3, 0), (8, 1, 0), (12, 5, 0), (10, 4, 0)] {
            let a = ffmt::fp(x, w, d, k);
            let b = ffmt::fp_ref(x, w, d, k);
            if a != b {
                fail(&format!("fp({:e},{},{},{})", x, w, d, k), &a, &b);
            }
            count += 1;
        }
        for &(w, d) in &[(11usize, 4usize), (12, 5), (10, 2), (12, 2), (13, 6), (11, 3), (10, 3)] {
            if !x.is_finite() {
                continue;
            }
            let a = ffmt::e1p(x, w, d);
            let b = ffmt::e1p_ref(x, w, d);
            if a != b {
                fail(&format!("e1p({:e},{},{})", x, w, d), &a, &b);
            }
            count += 1;
        }
        let _ = j;
    }
    // ENDF 11-column reals (a11), from a blank field and from a filled one
    for _ in 0..20_000_000u64 {
        let x = value(&mut r);
        if !x.is_finite() {
            continue;
        }
        let start = if r.below(2) == 0 { [b' '; 11] } else { *b"abcdefghijk" };
        let (mut a, mut b) = (start, start);
        rustjoy::endf::a11(x, &mut a);
        rustjoy::endf::a11_ref(x, &mut b);
        if a != b {
            fail(&format!("a11({:e})", x), &String::from_utf8_lossy(&a), &String::from_utf8_lossy(&b));
        }
        count += 1;
    }
    // input: random fields over a number alphabet, and ENDF-style numbers
    let alpha = b" 0123456789.+-eEdDqQx\t";
    let same = |a: &Result<f64, ffmt::ReadError>, b: &Result<f64, ffmt::ReadError>| match (a, b) {
        (Ok(x), Ok(y)) => x.to_bits() == y.to_bits(),
        (Err(x), Err(y)) => x.0 == y.0,
        _ => false,
    };
    for j in 0..20_000_000u64 {
        let mut f = Vec::new();
        if j % 2 == 0 {
            let n = 1 + r.below(14) as usize;
            for _ in 0..n {
                let c = if r.below(3) == 0 { alpha[r.below(alpha.len() as u64) as usize] } else { b"0123456789"[r.below(10) as usize] };
                f.push(c);
            }
        } else {
            let x = value(&mut r);
            let t = if x.is_finite() { rustjoy::endf::a11s(x) } else { format!("{}", x) };
            f.extend_from_slice(t.as_bytes());
        }
        let d = if r.below(4) == 0 { r.below(6) as usize } else { 0 };
        let a = ffmt::read_f(&f, d);
        let b = ffmt::read_f_ref(&f, d);
        if !same(&a, &b) {
            fail(&format!("read_f('{}',{})", String::from_utf8_lossy(&f), d), &format!("{:?}", a.map_err(|e| e.0)), &format!("{:?}", b.map_err(|e| e.0)));
        }
        let a = ffmt::read_i(&f).map_err(|e| e.0);
        let b = ffmt::read_i_ref(&f).map_err(|e| e.0);
        if a != b {
            fail(&format!("read_i('{}')", String::from_utf8_lossy(&f)), &format!("{:?}", a), &format!("{:?}", b));
        }
        let st = r.below(16) as usize;
        let w = 1 + r.below(12) as usize;
        if *ffmt::fieldx(&f, st, w) != ffmt::field(&f, st, w)[..] {
            fail("fieldx", "", "");
        }
        count += 3;
    }
    // sigfig (table of powers of ten) against sigfig_ref (log10 and ten**n)
    let sf = |x: f64, nd: i32, id: i32| {
        let a = rustjoy::util::sigfig(x, nd, id);
        let b = rustjoy::util::sigfig_ref(x, nd, id);
        if a.to_bits() != b.to_bits() {
            fail(&format!("sigfig({:e},{},{})", x, nd, id), &format!("{:e}", a), &format!("{:e}", b));
        }
    };
    for _ in 0..20_000_000u64 {
        let x = match r.below(3) {
            0 => value(&mut r),
            // energies in eV, as the modules use them
            1 => 10f64.powf(r.below(1_000_000_000) as f64 * 1e-8 - 5.0),
            // any positive finite value
            _ => f64::from_bits(r.next() >> 2),
        };
        let nd = match r.below(4) {
            0 => 7,
            1 => 8,
            2 => 9,
            _ => 1 + r.below(15) as i32,
        };
        sf(x, nd, r.below(3) as i32 - 1);
        count += 1;
    }
    // nint without libm against x.round()
    for _ in 0..20_000_000u64 {
        let x = match r.below(4) {
            0 => value(&mut r),
            1 => (r.below(2_000_001) as f64 - 1_000_000.0) / 2.0,
            2 => f64::from_bits(r.next()),
            _ => {
                let v = (r.next() >> 11) as f64 * 2f64.powi(r.below(64) as i32 - 52);
                let v = f64::from_bits((v.to_bits() as i64 + r.below(5) as i64 - 2) as u64);
                if r.below(2) == 0 { v } else { -v }
            }
        };
        let (a, b) = (rustjoy::fmath::nint(x), rustjoy::fmath::nint_ref(x));
        if a != b {
            fail(&format!("nint({:e})", x), &a.to_string(), &b.to_string());
        }
        count += 1;
    }
    // the values next to each power of ten, where log10 can round to it
    for k in -310..=310 {
        let p: f64 = format!("1e{}", k).parse().unwrap();
        for d in -3000i64..=3000 {
            let x = f64::from_bits((p.to_bits() as i64 + d) as u64);
            for nd in [7, 8, 9] {
                for id in [-1, 0, 1] {
                    sf(x, nd, id);
                    sf(-x, nd, id);
                    count += 2;
                }
            }
        }
    }
    println!("ffmt_check: {} cases identical", count);
}
