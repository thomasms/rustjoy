//! Exact replicas of the arithmetic that gfortran (at -O0, x86-64) emits for
//! Fortran intrinsics and operators.
//!
//! * `x**n` with a *variable* integer `n` calls libgcc's `__powidf2`; `powi`
//!   below is that algorithm, step for step.
//! * `x**k` with a *constant* integer `k` is expanded by the gfortran front end
//!   into a multiplication chain (`gfc_conv_cst_int_power`); `powc` below
//!   replicates that chain.
//! * Transcendental functions go straight to the system libm (glibc), exactly
//!   like gfortran does, so we declare them here and call them directly.

#![allow(dead_code)]

mod c {
    extern "C" {
        pub fn exp(x: f64) -> f64;
        pub fn log(x: f64) -> f64;
        pub fn log10(x: f64) -> f64;
        pub fn pow(x: f64, y: f64) -> f64;
        pub fn sin(x: f64) -> f64;
        pub fn cos(x: f64) -> f64;
        pub fn tan(x: f64) -> f64;
        pub fn atan(x: f64) -> f64;
        pub fn atan2(y: f64, x: f64) -> f64;
        pub fn sqrt(x: f64) -> f64;
        pub fn fmod(x: f64, y: f64) -> f64;
        pub fn sinh(x: f64) -> f64;
        pub fn cosh(x: f64) -> f64;
        pub fn tanh(x: f64) -> f64;
    }
}

#[inline(never)]
pub fn exp(x: f64) -> f64 { unsafe { c::exp(x) } }
#[inline(never)]
pub fn log(x: f64) -> f64 { unsafe { c::log(x) } }
#[inline(never)]
pub fn log10(x: f64) -> f64 { unsafe { c::log10(x) } }
/// `x**y` with a real exponent.
#[inline(never)]
pub fn pow(x: f64, y: f64) -> f64 { unsafe { c::pow(x, y) } }
#[inline(never)]
pub fn sin(x: f64) -> f64 { unsafe { c::sin(x) } }
#[inline(never)]
pub fn cos(x: f64) -> f64 { unsafe { c::cos(x) } }
#[inline(never)]
pub fn tan(x: f64) -> f64 { unsafe { c::tan(x) } }
#[inline(never)]
pub fn atan(x: f64) -> f64 { unsafe { c::atan(x) } }
#[inline(never)]
pub fn atan2(y: f64, x: f64) -> f64 { unsafe { c::atan2(y, x) } }
#[inline(never)]
pub fn sinh(x: f64) -> f64 { unsafe { c::sinh(x) } }
#[inline(never)]
pub fn cosh(x: f64) -> f64 { unsafe { c::cosh(x) } }
#[inline(never)]
pub fn tanh(x: f64) -> f64 { unsafe { c::tanh(x) } }
/// sqrt is correctly rounded in IEEE arithmetic; gfortran emits `sqrtsd`.
#[inline]
pub fn sqrt(x: f64) -> f64 { x.sqrt() }
/// Fortran `mod(a,p)` for reals: gfortran uses `fmod`.
#[inline(never)]
pub fn fmod(x: f64, y: f64) -> f64 { unsafe { c::fmod(x, y) } }

/// libgcc `__powidf2`: `x**m` for a run-time integer `m`.
pub fn powi(x: f64, m: i32) -> f64 {
    let mut n: u32 = if m < 0 { (m as u32).wrapping_neg() } else { m as u32 };
    let mut x = x;
    let mut y = if n % 2 != 0 { x } else { 1.0 };
    loop {
        n >>= 1;
        if n == 0 {
            break;
        }
        x = x * x;
        if n % 2 != 0 {
            y = y * x;
        }
    }
    if m < 0 { 1.0 / y } else { y }
}

/// `x**k` for a compile-time constant integer `k`. At -O0 gfortran does not
/// expand these into multiplications: it calls `__powidf2`, exactly as for a
/// run-time exponent (checked in the generated assembly).
pub fn powc(x: f64, k: i32) -> f64 {
    powi(x, k)
}

/// Fortran `nint` for real(8) -> default integer (gfortran: `lround`).
#[inline]
pub fn nint(x: f64) -> i32 {
    x.round() as i64 as i32
}

/// Fortran `int` (truncation toward zero).
#[inline]
pub fn int(x: f64) -> i32 {
    x.trunc() as i64 as i32
}

/// Integer `i**j` (gfortran `_gfortran_pow_i4_i4`).
pub fn ipow(a: i32, b: i32) -> i32 {
    if b < 0 {
        return match a {
            1 => 1,
            -1 => if b % 2 == 0 { 1 } else { -1 },
            _ => 0,
        };
    }
    let mut pow: i32 = 1;
    let mut x = a;
    let mut u = b as u32;
    while u != 0 {
        if u & 1 != 0 {
            pow = pow.wrapping_mul(x);
        }
        u >>= 1;
        if u != 0 {
            x = x.wrapping_mul(x);
        }
    }
    pow
}

/// Fortran `sign(a,b)` for reals.
#[inline]
pub fn sign(a: f64, b: f64) -> f64 {
    if b.is_sign_negative() { -a.abs() } else { a.abs() }
}

/// Fortran `aint`.
#[inline]
pub fn aint(x: f64) -> f64 { x.trunc() }

/// Fortran `anint` (round half away from zero, as `round`).
#[inline]
pub fn anint(x: f64) -> f64 { x.round() }
