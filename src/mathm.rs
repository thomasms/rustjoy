//! The parts of NJOY2016 `mathm.f90` that the ported modules use: `erfc`
//! (SLATEC derfc) with its helpers `csevl` and `initds`.

use crate::fmath::{exp, log, sqrt};
use std::sync::OnceLock;

include!("mathm_tables.rs");

struct ErfcInit {
    nterf: i32,
    nterfc: i32,
    nterc2: i32,
    xsml: f64,
    xmax: f64,
    sqeps: f64,
}

/// SLATEC `initds`. The Fortran adds `abs(real(os(i)))`: `real()` gives a
/// default (single precision) real, so each term is rounded to f32 first.
/// (Its message for a too-short series is never reached with these tables.)
fn initds(os: &[f64], nos: i32, eta: f64) -> i32 {
    let mut err = 0.0f64;
    let mut ii = 0;
    let mut i = 0;
    while err <= eta && ii < nos {
        ii += 1;
        i = nos + 1 - ii;
        err = err + ((os[(i - 1) as usize] as f32).abs() as f64);
    }
    i
}

/// SLATEC `dcsevl` as in mathm.f90 (its range checks never fail here).
fn csevl(x: f64, cs: &[f64], n: i32) -> f64 {
    let mut b1 = 0.0f64;
    let mut b0 = 0.0f64;
    let mut b2 = 0.0f64;
    let twox = 2.0 * x;
    for i in 1..=n {
        b2 = b1;
        b1 = b0;
        let ni = n + 1 - i;
        b0 = twox * b1 - b2 + cs[(ni - 1) as usize];
    }
    (b0 - b2) / 2.0
}

fn init() -> &'static ErfcInit {
    static I: OnceLock<ErfcInit> = OnceLock::new();
    I.get_or_init(|| {
        let one = 1.0f64;
        let d1mach1 = f64::MIN_POSITIVE; // tiny(one)
        let d1mach3 = 2.0f64.powi(-53); // two**(-digits(one)), exact
        let eta = d1mach3 / 10.0;
        let sqrtpi = SQRTPI;
        let txmax = sqrt(-log(sqrtpi * d1mach1));
        ErfcInit {
            nterf: initds(&ERFCS, 21, eta),
            nterfc: initds(&ERFCCS, 59, eta),
            nterc2: initds(&ERC2CS, 49, eta),
            xsml: -sqrt(-log(sqrtpi * d1mach3)),
            xmax: txmax - log(txmax) / (2.0 * txmax) - one / 10.0,
            sqeps: sqrt(2.0 * d1mach3),
        }
    })
}

const SQRTPI: f64 = 1.77245385090551602729816748334115e0;

/// Complementary error function, as `erfc` in NJOY2016 mathm.f90.
pub fn erfc(x: f64) -> f64 {
    let c = init();
    let one = 1.0f64;
    let zero = 0.0f64;
    let mut erfc = 0.0f64;
    if x <= c.xsml {
        erfc = 2.0;
    } else if x <= c.xmax {
        let mut y = x.abs();
        if y <= one {
            if y < c.sqeps {
                erfc = 1.0 - 2.0 * x / SQRTPI;
            }
            if y >= c.sqeps {
                erfc = 1.0 - x * (1.0 + csevl(2.0 * x * x - 1.0, &ERFCS, c.nterf));
            }
        } else {
            y = y * y;
            if y <= 4.0 * one {
                let xx = (8.0 / y - 5.0) / 3.0;
                erfc = exp(-y) / x.abs() * (one / 2.0 + csevl(xx, &ERC2CS, c.nterc2));
            }
            if y > 4.0 * one {
                let xx = 8.0 / y - 1.0;
                erfc = exp(-y) / x.abs() * (one / 2.0 + csevl(xx, &ERFCCS, c.nterfc));
            }
            if x < zero {
                erfc = 2.0 - erfc;
            }
        }
    } else {
        erfc = 0.0;
    }
    erfc
}
