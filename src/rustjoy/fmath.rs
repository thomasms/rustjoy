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
        pub fn acos(x: f64) -> f64;
    }
}

/// `exp`. On glibc, the exported `exp` is a wrapper that calls
/// `__exp_finite` and returns its value unchanged (it only sets errno on
/// overflow and underflow), so the port calls `__exp_finite` directly, which
/// saves one call level in the hot loops. The function is looked up once;
/// it is used only if it gives the same bits as `exp` on a set of test
/// values, else (or when it is missing) the port calls `exp`.
#[cfg(all(target_os = "linux", target_env = "gnu"))]
#[inline]
pub fn exp(x: f64) -> f64 {
    let f = EXP_FN.load(std::sync::atomic::Ordering::Relaxed);
    let f: unsafe extern "C" fn(f64) -> f64 = if f.is_null() { exp_init() } else { unsafe { std::mem::transmute(f) } };
    unsafe { f(x) }
}

/// The libm `exp` itself (the reference for `exp`).
#[inline(never)]
pub fn exp_libm(x: f64) -> f64 { unsafe { c::exp(x) } }

#[cfg(all(target_os = "linux", target_env = "gnu"))]
static EXP_FN: std::sync::atomic::AtomicPtr<()> = std::sync::atomic::AtomicPtr::new(std::ptr::null_mut());

#[cfg(all(target_os = "linux", target_env = "gnu"))]
#[cold]
fn exp_init() -> unsafe extern "C" fn(f64) -> f64 {
    extern "C" {
        fn dlvsym(handle: *mut std::ffi::c_void, symbol: *const std::ffi::c_char, version: *const std::ffi::c_char) -> *mut std::ffi::c_void;
    }
    unsafe extern "C" fn libm_exp(x: f64) -> f64 { c::exp(x) }
    let mut f: unsafe extern "C" fn(f64) -> f64 = libm_exp;
    // RTLD_DEFAULT is a null handle on glibc
    let p = unsafe { dlvsym(std::ptr::null_mut(), c"__exp_finite".as_ptr(), c"GLIBC_2.15".as_ptr()) };
    if !p.is_null() {
        let g: unsafe extern "C" fn(f64) -> f64 = unsafe { std::mem::transmute(p) };
        let mut ok = true;
        let mut t = 0x9e37_79b9_7f4a_7c15u64;
        for i in 0..20_000u32 {
            t ^= t << 13;
            t ^= t >> 7;
            t ^= t << 17;
            let x = match i % 4 {
                0 => (t >> 11) as f64 / (1u64 << 53) as f64 * 1500.0 - 760.0,
                1 => ((t >> 11) as f64 / (1u64 << 53) as f64 - 0.5) * 2.0,
                2 => f64::from_bits(t),
                _ => [0.0, -0.0, 709.78, 709.79, -745.13, -745.14, -708.4, f64::INFINITY, -f64::INFINITY, 1e-300][(i / 4 % 10) as usize],
            };
            let (a, b) = (unsafe { g(x) }, unsafe { c::exp(x) });
            if a.to_bits() != b.to_bits() && !(a.is_nan() && b.is_nan()) {
                ok = false;
                break;
            }
        }
        if ok {
            f = g;
        }
    }
    EXP_FN.store(f as *mut (), std::sync::atomic::Ordering::Relaxed);
    f
}

/// `exp` on other systems (macOS, musl, ...): the libm `exp`, as gfortran
/// calls it there.
#[cfg(not(all(target_os = "linux", target_env = "gnu")))]
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
#[inline]
pub fn acos(x: f64) -> f64 { unsafe { c::acos(x) } }
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
    // round half away from zero, as x.round(), without the libm call: for
    // |x| < 2^52, t = trunc(x) and x - t are exact
    if x.abs() < 4_503_599_627_370_496.0 {
        let t = x as i64;
        let f = x - t as f64;
        (t + (f >= 0.5) as i64 - (f <= -0.5) as i64) as i32
    } else {
        x.round() as i64 as i32
    }
}

/// The original `nint` (`x.round()`), the reference for `nint`.
pub fn nint_ref(x: f64) -> i32 {
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
