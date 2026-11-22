//! 1-based array views, so that ported code can keep the Fortran index
//! expressions unchanged: `a[k]` here is `a(k)` in the Fortran source.
//!
//! A Fortran actual argument `res(jnow)` passed to a dummy `a(*)` becomes
//! `res.at(jnow)`, a view whose element 1 is `res(jnow)`.

#![allow(dead_code)]

use std::ops::{Index, IndexMut};

/// Owning 1-based real array (`real(kr),dimension(:),allocatable`).
#[derive(Clone, Default, Debug)]
pub struct RV(pub Vec<f64>);

impl RV {
    pub fn new(n: usize) -> Self {
        RV(vec![0.0; n])
    }
    pub fn len(&self) -> i32 {
        self.0.len() as i32
    }
    /// View starting at Fortran index j.
    pub fn at(&mut self, j: i32) -> F<'_> {
        F(&mut self.0[(j - 1) as usize..])
    }
    pub fn all(&mut self) -> F<'_> {
        F(&mut self.0[..])
    }
}

impl Index<i32> for RV {
    type Output = f64;
    #[inline]
    fn index(&self, i: i32) -> &f64 {
        &self.0[(i - 1) as usize]
    }
}
impl IndexMut<i32> for RV {
    #[inline]
    fn index_mut(&mut self, i: i32) -> &mut f64 {
        &mut self.0[(i - 1) as usize]
    }
}

/// Owning 1-based integer array.
#[derive(Clone, Default, Debug)]
pub struct IV(pub Vec<i32>);

impl IV {
    pub fn new(n: usize) -> Self {
        IV(vec![0; n])
    }
}
impl Index<i32> for IV {
    type Output = i32;
    #[inline]
    fn index(&self, i: i32) -> &i32 {
        &self.0[(i - 1) as usize]
    }
}
impl IndexMut<i32> for IV {
    #[inline]
    fn index_mut(&mut self, i: i32) -> &mut i32 {
        &mut self.0[(i - 1) as usize]
    }
}

/// Borrowed 1-based view (`real(kr)::a(*)`).
pub struct F<'a>(pub &'a mut [f64]);

impl<'a> F<'a> {
    pub fn at(&mut self, j: i32) -> F<'_> {
        F(&mut self.0[(j - 1) as usize..])
    }
    pub fn len(&self) -> i32 {
        self.0.len() as i32
    }
}

impl<'a> Index<i32> for F<'a> {
    type Output = f64;
    #[inline]
    fn index(&self, i: i32) -> &f64 {
        &self.0[(i - 1) as usize]
    }
}
impl<'a> IndexMut<i32> for F<'a> {
    #[inline]
    fn index_mut(&mut self, i: i32) -> &mut f64 {
        &mut self.0[(i - 1) as usize]
    }
}

/// Column-major 2-D real array with 1-based indices (`a(n1,n2)`).
#[derive(Clone, Default, Debug)]
pub struct RM {
    pub n1: usize,
    pub v: Vec<f64>,
}

impl RM {
    pub fn new(n1: usize, n2: usize) -> Self {
        RM { n1, v: vec![0.0; n1 * n2] }
    }
}
impl Index<(i32, i32)> for RM {
    type Output = f64;
    #[inline]
    fn index(&self, (i, j): (i32, i32)) -> &f64 {
        &self.v[(i - 1) as usize + (j - 1) as usize * self.n1]
    }
}
impl IndexMut<(i32, i32)> for RM {
    #[inline]
    fn index_mut(&mut self, (i, j): (i32, i32)) -> &mut f64 {
        &mut self.v[(i - 1) as usize + (j - 1) as usize * self.n1]
    }
}

/// Pack 4-character Hollerith text into the bytes of a real array, as a
/// Fortran `equivalence` of `character(4)::t(17)` with `real(kr)::r(17)` does:
/// the 68 characters occupy the first 68 bytes of the 136-byte array.
pub fn pack_text(text68: &[u8], r: &mut [f64]) {
    let mut bytes: Vec<u8> = Vec::with_capacity(r.len() * 8);
    for x in r.iter() {
        bytes.extend_from_slice(&x.to_le_bytes());
    }
    for (i, b) in text68.iter().enumerate() {
        bytes[i] = *b;
    }
    for (i, x) in r.iter_mut().enumerate() {
        *x = f64::from_le_bytes(bytes[8 * i..8 * i + 8].try_into().unwrap());
    }
}

/// The first 68 bytes of a 17-word real array, as text.
pub fn unpack_text(r: &[f64]) -> Vec<u8> {
    let mut bytes: Vec<u8> = Vec::with_capacity(136);
    for x in r.iter().take(17) {
        bytes.extend_from_slice(&x.to_le_bytes());
    }
    bytes.truncate(68);
    bytes
}
