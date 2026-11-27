//! Prints erfc(x) and hnabb(n,a,b), n=0..4, in hex for each (x,a,b) triple
//! on standard input, in the layout of tools/numcheck_broadr.f90.
use rustjoy::{broadr::hnabb, mathm::erfc};
use std::io::{Read, Write};

fn main() {
    let mut b = Vec::new();
    std::io::stdin().read_to_end(&mut b).unwrap();
    let out = std::io::stdout();
    let mut o = std::io::BufWriter::new(out.lock());
    for r in b.chunks_exact(24) {
        let g = |i: usize| f64::from_le_bytes(r[i..i + 8].try_into().unwrap());
        let (x, a, bb) = (g(0), g(8), g(16));
        write!(o, "{:016X}", erfc(x).to_bits()).unwrap();
        for n in 0..5 {
            write!(o, " {:016X}", hnabb(n, a, bb).to_bits()).unwrap();
        }
        writeln!(o).unwrap();
    }
}
