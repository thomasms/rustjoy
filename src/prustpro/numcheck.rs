//! Hidden subcommand `prustpro _numcheck`: the Rust side of
//! `tools/prepro/numcheck.f` (same input, same output), to check the
//! number editing of `endfio` against the Fortran.

use crate::endfio::{in9, incore10, incore7, incore9, out10, out7, out9, out9g};
use std::io::BufRead;

pub fn run() -> ! {
    let stdin = std::io::stdin();
    for line in stdin.lock().lines() {
        let line = line.unwrap();
        let b = line.as_bytes();
        if b.first() == Some(&b'S') {
            let mut f = [b' '; 11];
            for k in 0..11 {
                if 2 + k < b.len() {
                    f[k] = b[2 + k];
                }
            }
            let x = in9(&f);
            crate::fio::out(format!("{:016X}", x.to_bits()));
            continue;
        }
        let bits = u64::from_str_radix(line[2..18].trim(), 16).unwrap();
        let x = f64::from_bits(bits);
        let (mut f, mut g, mut h, mut p) = ([b'#'; 11], [b'#'; 11], [b'#'; 11], [b'#'; 11]);
        out9(x, &mut f);
        out9g(x, &mut g);
        out10(x, &mut h);
        out7(x, &mut p);
        let (mut x1, mut x2, mut x3) = (x, x, x);
        incore10(&mut x1);
        incore9(&mut x2);
        incore7(&mut x3);
        let s = |a: &[u8; 11]| String::from_utf8_lossy(a).into_owned();
        crate::fio::out(format!(
            "{}|{}|{}|{}|{:016X} {:016X} {:016X}",
            s(&f),
            s(&g),
            s(&h),
            s(&p),
            x1.to_bits(),
            x2.to_bits(),
            x3.to_bits()
        ));
    }
    crate::fio::flush_stdout();
    std::process::exit(0);
}
