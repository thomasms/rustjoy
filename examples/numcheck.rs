//! Prints the results of the basic numeric routines for every input value,
//! in the same layout as tools/numcheck.f90, for byte comparison.
use reconr_rs::{endf, ffmt, fmath, util};
use std::io::{Read, Write};

fn main() {
    let mut b = Vec::new();
    std::io::stdin().read_to_end(&mut b).unwrap();
    let out = std::io::stdout();
    let mut o = std::io::BufWriter::new(out.lock());
    for ch in b.chunks_exact(12) {
        let x = f64::from_le_bytes(ch[..8].try_into().unwrap());
        let n = i32::from_le_bytes(ch[8..12].try_into().unwrap());
        let a = endf::a11s(x);
        let t = util::a10(x);
        let s7 = util::sigfig(x, 7, 0).to_bits();
        let s9 = util::sigfig(x, 9, -1).to_bits();
        let pw = fmath::powi(x, n).to_bits();
        let p2 = fmath::powc(x, 2).to_bits();
        let p3 = fmath::powc(x, 3).to_bits();
        let p5 = fmath::powc(x, 5).to_bits();
        let p13 = fmath::powc(x, 13).to_bits();
        let r = endf::rd_e(a.as_bytes()).to_bits();
        let e = ffmt::e1p(x, 12, 5);
        let f = ffmt::f(x, 10, 3);
        writeln!(o, "{}|{}|{:016X}|{:016X}|{:016X}|{:016X}|{:016X}|{:016X}|{:016X}|{:016X}|{}|{}", a, t, s7, s9, pw, p2, p3, p5, p13, r, e, f).unwrap();
    }
}
