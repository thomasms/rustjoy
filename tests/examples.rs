//! Compare rustjoy with NJOY2016 on the input decks in `examples/`.
//!
//! Each example is a directory in `examples/` with an NJOY input deck
//! (`input`) and its tapes (`tapeNN`, links to the ENDF files in
//! `examples/endf/`). A test copies the directory to a new work directory,
//! runs the whole deck with rustjoy and compares the results with NJOY2016
//! byte for byte: the exit status, standard output, the listing `output` and
//! every tape that the deck writes.
//!
//! There are two ways to get the NJOY2016 results:
//!
//! * Stored: `examples/NAME/reference.txt` holds the exit status and a
//!   SHA-256 checksum of each NJOY2016 output file. It is made from
//!   NJOY2016 commit ab27c64 built with its own CMake defaults (gfortran
//!   -O0) on Linux x86-64. `cargo test` uses it and needs no Fortran.
//! * Live: with `NJOY=/path/to/njoy`, the test also runs NJOY2016 on the
//!   same deck and compares the files directly. A difference is then shown
//!   with its first byte offset and, for text files, the line.
//!
//! With `NJOY` set and `RUSTJOY_UPDATE_REFERENCES=1`, the tests write
//! `reference.txt` again from the NJOY2016 results (never from rustjoy).
//!
//! Before the comparison, only bytes that differ between two runs of
//! NJOY2016 itself are masked (the same rules as `tools/compare.py`):
//! CPU times and the date and time in text files, and bytes 68..135 of the
//! tape label and of the MF1/MT451 Hollerith cards in binary tapes (the
//! Fortran fills them from uninitialised memory).
//!
//! Build with `--release` (`cargo test --release`): the debug build is many
//! times slower.

use std::collections::BTreeMap;
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

/// One test for each example directory.
macro_rules! examples {
    ($($test:ident => $dir:literal,)*) => {
        $(
            #[test]
            fn $test() {
                check_example($dir);
            }
        )*

        const EXAMPLES: &[&str] = &[$($dir),*];
    };
}

examples! {
    b10_gaspr => "b10-gaspr",
    fe_groupr => "fe-groupr",
    h2_reconr => "h2-reconr",
    mixr_h2_b10 => "mixr-h2-b10",
    parah_leapr => "parah-leapr",
    photon_gaminr => "photon-gaminr",
    pu241_unresr => "pu241-unresr",
    sn119_heatr => "sn119-heatr",
    sn119_purr => "sn119-purr",
    sn119_resxsr => "sn119-resxsr",
    zrh_thermr => "zrh-thermr",
}

/// A new example directory must also get a line in `examples!` above.
#[test]
fn every_example_has_a_test() {
    let mut found: Vec<String> = fs::read_dir(examples_dir())
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.join("input").is_file())
        .map(|p| p.file_name().unwrap().to_string_lossy().into_owned())
        .collect();
    found.sort();
    let mut listed: Vec<String> = EXAMPLES.iter().map(|s| s.to_string()).collect();
    listed.sort();
    assert_eq!(found, listed, "the examples! list in tests/examples.rs does not match examples/");
}

/// The results of one run: exit status and masked output files by name.
/// Standard output is stored under the name `stdout`.
struct RunResult {
    status: i32,
    files: BTreeMap<String, Vec<u8>>,
}

fn check_example(name: &str) {
    let src = examples_dir().join(name);
    let work = Path::new(env!("CARGO_TARGET_TMPDIR")).join("examples").join(name);
    if work.exists() {
        // Only the work directory of this example from an earlier run.
        fs::remove_dir_all(&work).unwrap();
    }

    let rust = run_deck(Path::new(env!("CARGO_BIN_EXE_rustjoy")), &src, &work.join("rustjoy"));
    let mut errors = Vec::new();

    if let Some(njoy) = std::env::var_os("NJOY") {
        let njoy = run_deck(Path::new(&njoy), &src, &work.join("njoy"));
        errors.extend(compare_runs(&njoy, &rust));
        if std::env::var_os("RUSTJOY_UPDATE_REFERENCES").is_some() {
            fs::write(src.join("reference.txt"), reference_text(&njoy)).unwrap();
        } else if let Some(reference) = read_reference(&src) {
            let e = compare_with_reference(&reference, &njoy);
            if !e.is_empty() {
                errors.push(format!(
                    "this NJOY2016 build does not give the stored results in {}/reference.txt \
                     (another NJOY2016 version, compiler or platform?):\n    {}",
                    name,
                    e.join("\n    ")
                ));
            }
        }
    }

    match read_reference(&src) {
        Some(reference) => errors.extend(compare_with_reference(&reference, &rust)),
        None => errors.push(format!("examples/{}/reference.txt is missing", name)),
    }

    if errors.is_empty() {
        // Keep the outputs only when something is wrong; some are large.
        fs::remove_dir_all(&work).unwrap();
    } else {
        panic!(
            "example {}: rustjoy differs from NJOY2016\n  {}\nThe output files are in {}",
            name,
            errors.join("\n  "),
            work.display()
        );
    }
}

fn examples_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("examples")
}

/// Copy the example to `dir` (links are copied as files), run the deck with
/// `program` there, and return its masked results. Input files that the deck
/// does not change are left out.
fn run_deck(program: &Path, src: &Path, dir: &Path) -> RunResult {
    fs::create_dir_all(dir).unwrap();
    let mut inputs = BTreeMap::new();
    for entry in fs::read_dir(src).unwrap() {
        let path = entry.unwrap().path();
        let name = path.file_name().unwrap().to_string_lossy().into_owned();
        if name == "reference.txt" || name.starts_with("README") {
            continue;
        }
        let bytes = fs::read(&path).unwrap();
        fs::write(dir.join(&name), &bytes).unwrap();
        inputs.insert(name, bytes);
    }
    let out = Command::new(program)
        .current_dir(dir)
        .stdin(fs::File::open(dir.join("input")).unwrap())
        .stdout(Stdio::piped())
        .stderr(Stdio::inherit())
        .output()
        .unwrap_or_else(|e| panic!("cannot run {}: {}", program.display(), e));
    let mut files = BTreeMap::new();
    files.insert("stdout".to_string(), mask_text(&out.stdout));
    for entry in fs::read_dir(dir).unwrap() {
        let path = entry.unwrap().path();
        let name = path.file_name().unwrap().to_string_lossy().into_owned();
        let bytes = fs::read(&path).unwrap();
        if inputs.get(&name) == Some(&bytes) {
            continue;
        }
        let masked = if name == "output" { mask_text(&bytes) } else { mask_tape(&bytes) };
        files.insert(name, masked);
    }
    RunResult { status: out.status.code().unwrap_or(-1), files }
}

fn compare_runs(njoy: &RunResult, rust: &RunResult) -> Vec<String> {
    let mut errors = Vec::new();
    if njoy.status != rust.status {
        errors.push(format!("exit status: NJOY2016 {}, rustjoy {}", njoy.status, rust.status));
    }
    let names: std::collections::BTreeSet<&String> = njoy.files.keys().chain(rust.files.keys()).collect();
    for name in names {
        match (njoy.files.get(name), rust.files.get(name)) {
            (Some(a), Some(b)) if a == b => {}
            (Some(a), Some(b)) => errors.push(format!("{}: {}", name, first_difference(a, b))),
            (Some(_), None) => errors.push(format!("{}: written by NJOY2016 but not by rustjoy", name)),
            (None, _) => errors.push(format!("{}: written by rustjoy but not by NJOY2016", name)),
        }
    }
    errors
}

/// The stored results: exit status and (length, SHA-256) by file name.
struct Reference {
    status: i32,
    files: BTreeMap<String, (usize, String)>,
}

fn read_reference(src: &Path) -> Option<Reference> {
    let text = fs::read_to_string(src.join("reference.txt")).ok()?;
    let mut status = None;
    let mut files = BTreeMap::new();
    for line in text.lines() {
        let words: Vec<&str> = line.split_whitespace().collect();
        match words.as_slice() {
            [] => {}
            [w, ..] if w.starts_with('#') => {}
            ["exit", s] => status = Some(s.parse().expect("exit status in reference.txt")),
            [hash, len, name] => {
                files.insert(name.to_string(), (len.parse().expect("length in reference.txt"), hash.to_string()));
            }
            _ => panic!("bad line in {}/reference.txt: {}", src.display(), line),
        }
    }
    Some(Reference { status: status.expect("no exit line in reference.txt"), files })
}

fn reference_text(njoy: &RunResult) -> String {
    let mut s = String::from(
        "# NJOY2016 results for this example, checked by tests/examples.rs.\n\
         # Made with NJOY2016 commit ab27c64 (gfortran, CMake defaults, Linux x86-64).\n\
         # Each file line: SHA-256 of the masked bytes, length in bytes, file name.\n",
    );
    s += &format!("exit {}\n", njoy.status);
    for (name, bytes) in &njoy.files {
        s += &format!("{} {} {}\n", sha256_hex(bytes), bytes.len(), name);
    }
    s
}

fn compare_with_reference(reference: &Reference, run: &RunResult) -> Vec<String> {
    let mut errors = Vec::new();
    if reference.status != run.status {
        errors.push(format!("exit status: reference {}, this run {}", reference.status, run.status));
    }
    for (name, (len, hash)) in &reference.files {
        match run.files.get(name) {
            None => errors.push(format!("{}: in reference.txt but not written", name)),
            Some(b) if b.len() != *len => {
                errors.push(format!("{}: {} bytes, reference has {} bytes", name, b.len(), len))
            }
            Some(b) if sha256_hex(b) != *hash => errors.push(format!(
                "{}: checksum differs from reference.txt (run with NJOY=/path/to/njoy to see the first difference)",
                name
            )),
            Some(_) => {}
        }
    }
    for name in run.files.keys() {
        if !reference.files.contains_key(name) {
            errors.push(format!("{}: written but not in reference.txt", name));
        }
    }
    errors
}

/// Where two masked files first differ, with the line for text files.
fn first_difference(a: &[u8], b: &[u8]) -> String {
    let off = a.iter().zip(b).position(|(x, y)| x != y).unwrap_or(a.len().min(b.len()));
    let mut s = format!(
        "differs at byte {} (NJOY2016 {} bytes, rustjoy {} bytes)",
        off,
        a.len(),
        b.len()
    );
    if !is_unformatted(a) {
        let line = a[..off].iter().filter(|&&c| c == b'\n').count();
        let get = |x: &[u8]| {
            x.split(|&c| c == b'\n').nth(line).map(|l| String::from_utf8_lossy(l).into_owned()).unwrap_or_default()
        };
        s += &format!("\n      line {}\n      NJOY2016: {}\n      rustjoy:  {}", line + 1, get(a), get(b));
    }
    s
}

// ---------------------------------------------------------------------------
// Masks: the same rules as mask_text and mask_tape in tools/compare.py.

/// Mask CPU times at line ends (`[ 0-9.-]{8}s` that reads as `f8.1`), and
/// `date: 99/99/99`, `time: 99:99:99` and `99/99/99 99:99:99` (9: a digit).
fn mask_text(b: &[u8]) -> Vec<u8> {
    let mut lines: Vec<Vec<u8>> = b.split(|&c| c == b'\n').map(|l| l.to_vec()).collect();
    for line in &mut lines {
        let n = line.len();
        if n >= 9 && line[n - 1] == b's' {
            let f = &line[n - 9..n - 1];
            if f.iter().all(|c| b" 0123456789.-".contains(c)) && is_time_field(f) {
                line.splice(n - 9.., b"#time#s".iter().copied());
            }
        }
    }
    let mut out = lines.join(&b'\n');
    replace_pattern(&mut out, b"date: 99/99/99", b"date: ##/##/##");
    replace_pattern(&mut out, b"time: 99:99:99", b"time: ##:##:##");
    replace_pattern(&mut out, b"99/99/99 99:99:99", b"##/##/## ##:##:##");
    out
}

/// `^ *-?\d+\.\d$`
fn is_time_field(f: &[u8]) -> bool {
    let s = f.iter().position(|&c| c != b' ').unwrap_or(f.len());
    let mut t = &f[s..];
    if t.first() == Some(&b'-') {
        t = &t[1..];
    }
    let n = t.len();
    n >= 3
        && t[n - 2] == b'.'
        && t[n - 1].is_ascii_digit()
        && t[..n - 2].iter().all(u8::is_ascii_digit)
}

/// Replace each match of `pat` (`9` is any digit, other bytes match
/// themselves) from left to right, without overlaps, with `with`.
fn replace_pattern(b: &mut Vec<u8>, pat: &[u8], with: &[u8]) {
    let n = pat.len();
    let mut i = 0;
    while i + n <= b.len() {
        let hit = pat.iter().zip(&b[i..i + n]).all(|(&p, &c)| if p == b'9' { c.is_ascii_digit() } else { p == c });
        if hit {
            b.splice(i..i + n, with.iter().copied());
            i += with.len();
        } else {
            i += 1;
        }
    }
}

fn int_at(b: &[u8], i: usize) -> i32 {
    i32::from_le_bytes([b[i], b[i + 1], b[i + 2], b[i + 3]])
}

/// A gfortran unformatted sequential file: 4-byte record lengths before and
/// after each record.
fn is_unformatted(b: &[u8]) -> bool {
    if b.len() < 8 {
        return false;
    }
    let n = int_at(b, 0);
    n > 0 && (n as usize) + 8 <= b.len() && b[4 + n as usize..8 + n as usize] == b[..4]
}

/// Mask bytes 68..135 of the tape-label record and of each 136-byte
/// Hollerith card in binary MF1/MT451 records.
fn mask_tape(b: &[u8]) -> Vec<u8> {
    let mut b = b.to_vec();
    if !is_unformatted(&b) {
        return b;
    }
    let mut off = 0;
    let mut first = true;
    let mut in_hol = false;
    while off + 4 <= b.len() {
        let n = int_at(&b, off) as usize;
        let d0 = off + 4;
        if d0 + 20 > b.len() {
            break;
        }
        let (mf, mt, nb, nw) = (int_at(&b, d0 + 4), int_at(&b, d0 + 8), int_at(&b, d0 + 12), int_at(&b, d0 + 16));
        let data = d0 + 20;
        let mut zero = |b: &mut Vec<u8>, from: usize| {
            for x in &mut b[from + 68..from + 136] {
                *x = 0;
            }
        };
        if first {
            if mf == 0 && mt == 0 && nw == 17 && n == 20 + 136 {
                zero(&mut b, data);
            }
            first = false;
        } else if mf == 1 && mt == 451 {
            let mut blocks = None;
            if !in_hol && nw > 6 && (nw - 6) % 17 == 0 && n == 20 + 8 * nw as usize {
                blocks = Some((data + 48, (nw - 6) / 17));
                in_hol = nb > 0;
            } else if in_hol && nw % 17 == 0 {
                blocks = Some((data, nw / 17));
                in_hol = nb > 0;
            }
            if let Some((at, ncard)) = blocks {
                for k in 0..ncard as usize {
                    zero(&mut b, at + 136 * k);
                }
            }
        } else {
            in_hol = false;
        }
        off += 8 + n;
    }
    b
}

// ---------------------------------------------------------------------------
// SHA-256 (FIPS 180-4), so that the test needs no dependencies.

fn sha256_hex(data: &[u8]) -> String {
    const K: [u32; 64] = [
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
        0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
        0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
        0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7, 0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
        0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
        0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
        0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2,
    ];
    let mut h: [u32; 8] =
        [0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab, 0x5be0cd19];
    let mut msg = data.to_vec();
    msg.push(0x80);
    while msg.len() % 64 != 56 {
        msg.push(0);
    }
    msg.extend_from_slice(&((data.len() as u64) * 8).to_be_bytes());
    for chunk in msg.chunks(64) {
        let mut w = [0u32; 64];
        for i in 0..16 {
            w[i] = u32::from_be_bytes([chunk[4 * i], chunk[4 * i + 1], chunk[4 * i + 2], chunk[4 * i + 3]]);
        }
        for i in 16..64 {
            let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
            let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16].wrapping_add(s0).wrapping_add(w[i - 7]).wrapping_add(s1);
        }
        let mut v = h;
        for i in 0..64 {
            let s1 = v[4].rotate_right(6) ^ v[4].rotate_right(11) ^ v[4].rotate_right(25);
            let ch = (v[4] & v[5]) ^ (!v[4] & v[6]);
            let t1 = v[7].wrapping_add(s1).wrapping_add(ch).wrapping_add(K[i]).wrapping_add(w[i]);
            let s0 = v[0].rotate_right(2) ^ v[0].rotate_right(13) ^ v[0].rotate_right(22);
            let maj = (v[0] & v[1]) ^ (v[0] & v[2]) ^ (v[1] & v[2]);
            let t2 = s0.wrapping_add(maj);
            v = [t1.wrapping_add(t2), v[0], v[1], v[2], v[3].wrapping_add(t1), v[4], v[5], v[6]];
        }
        for i in 0..8 {
            h[i] = h[i].wrapping_add(v[i]);
        }
    }
    h.iter().map(|x| format!("{:08x}", x)).collect()
}

#[test]
fn sha256_known_values() {
    assert_eq!(sha256_hex(b""), "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
    assert_eq!(sha256_hex(b"abc"), "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
    let long = vec![b'a'; 1_000_000];
    assert_eq!(sha256_hex(&long), "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0");
}

#[test]
fn masks_match_compare_py() {
    let t = b" reconr...                                                   0.1s\n\
              date: 11/27/26  time: 20:54:01\n 11/27/26 20:54:01\n";
    let m = mask_text(t);
    let s = String::from_utf8(m).unwrap();
    assert!(s.contains("#time#s\n"), "{}", s);
    assert!(s.contains("date: ##/##/##  time: ##:##:##"), "{}", s);
    assert!(s.contains("##/##/## ##:##:##"), "{}", s);
    // a line that ends in 's' but not with an f8.1 time is not changed
    assert_eq!(mask_text(b"   1234 cards"), b"   1234 cards".to_vec());
}
