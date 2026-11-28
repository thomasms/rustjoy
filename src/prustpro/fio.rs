//! Fortran sequential I/O as the PREPRO codes use it (gfortran runtime).
//!
//! * Formatted units: one line per record. A read returns the line without
//!   its line end (a CR before the LF is removed, as libgfortran does); the
//!   callers pad it with blanks (PAD='YES'). A unit opened on a file that
//!   does not exist is created on the first write, as gfortran does.
//! * Unformatted scratch units: records are kept in memory.
//! * Unit `*` (standard output) is buffered and flushed at the end.

#![allow(dead_code)]

use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::sync::Mutex;

/// One value in an unformatted record.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Val {
    I(i32),
    D(f64),
    /// CHARACTER*1 value.
    C(u8),
}

impl Val {
    pub fn d(&self) -> f64 {
        match *self {
            Val::D(x) => x,
            Val::I(i) => i as f64,
            Val::C(c) => c as f64,
        }
    }
    pub fn i(&self) -> i32 {
        match *self {
            Val::I(i) => i,
            Val::D(x) => x as i32,
            Val::C(c) => c as i32,
        }
    }
}

enum Mode {
    /// Opened, nothing read or written yet.
    Fresh,
    Read(BufReader<File>),
    Write(BufWriter<File>),
    /// End of file was read (a further read is an error in gfortran;
    /// the PREPRO codes always stop at END=).
    Eof,
}

struct FmtUnit {
    pub(crate) path: String,
    mode: Mode,
    /// Bytes read so far: a WRITE after reads starts here (the rest of
    /// the file is cut off, as for any sequential write).
    nread: u64,
}

struct Scratch {
    recs: Vec<Vec<Val>>,
    pos: usize,
}

/// A formatted scratch file (OPEN(u,STATUS='SCRATCH')), kept in memory.
struct FmtScratch {
    lines: Vec<Vec<u8>>,
    pos: usize,
}

enum Unit {
    Fmt(FmtUnit),
    Scr(Scratch),
    FScr(FmtScratch),
}

/// Outcome of reading one record.
pub enum Rd {
    Line(Vec<u8>),
    Eof,
}

pub struct Units {
    units: HashMap<i32, Unit>,
}

/// Program name printed by gfortran runtime errors.
pub static PROGRAM: Mutex<String> = Mutex::new(String::new());

/// Standard output (Fortran unit *).
pub static STDOUT: Mutex<Option<BufWriter<std::io::Stdout>>> = Mutex::new(None);

/// Write one or more lines to standard output. `s` holds whole lines
/// separated by '\n' and no final '\n'.
pub fn out<S: AsRef<[u8]>>(s: S) {
    let mut g = STDOUT.lock().unwrap();
    if g.is_none() {
        *g = Some(BufWriter::with_capacity(1 << 16, std::io::stdout()));
    }
    let w = g.as_mut().unwrap();
    let _ = w.write_all(s.as_ref());
    let _ = w.write_all(b"\n");
}

pub fn flush_stdout() {
    let mut g = STDOUT.lock().unwrap();
    if let Some(w) = g.as_mut() {
        let _ = w.flush();
    }
}

/// A gfortran runtime error, raised as a panic so that the open units are
/// flushed while the stack unwinds (gfortran flushes them before it
/// stops); `main` prints the message and exits with status 2.
pub struct RuntimeError(pub String);

/// gfortran runtime error: message on standard error and exit status 2.
pub fn runtime_error(unit: i32, path: &str, msg: &str) -> ! {
    let m = if path.is_empty() {
        format!("At line 0 of file {} (unit = {})\nFortran runtime error: {}", PROGRAM.lock().unwrap(), unit, msg)
    } else {
        format!(
            "At line 0 of file {} (unit = {}, file = '{}')\nFortran runtime error: {}",
            PROGRAM.lock().unwrap(),
            unit,
            path,
            msg
        )
    };
    std::panic::panic_any(RuntimeError(m));
}

/// Fortran STOP: flush everything and exit with status 0.
pub fn stop(units: &mut Units) -> ! {
    units.close_all();
    flush_stdout();
    std::process::exit(0);
}

impl Default for Units {
    fn default() -> Self {
        Self::new()
    }
}

impl Drop for Units {
    fn drop(&mut self) {
        self.close_all();
    }
}

impl Units {
    pub fn new() -> Units {
        Units { units: HashMap::new() }
    }

    /// OPEN(u, FILE=path) with STATUS='UNKNOWN'. The name is trimmed of
    /// trailing blanks as gfortran does.
    pub fn open(&mut self, u: i32, path: &str) {
        self.close(u);
        let p = path.trim_end().to_string();
        // gfortran creates a file that does not exist at OPEN; an existing
        // file is left as it is until the first write.
        if !std::path::Path::new(&p).exists() {
            let _ = OpenOptions::new().write(true).create(true).truncate(false).open(&p);
        }
        self.units.insert(
            u,
            Unit::Fmt(FmtUnit { path: p, mode: Mode::Fresh, nread: 0 }),
        );
    }

    /// OPEN(u, FILE=path, STATUS='OLD', ERR=...): false if the file cannot
    /// be opened for reading.
    pub fn open_old(&mut self, u: i32, path: &str) -> bool {
        let p = path.trim_end().to_string();
        match File::open(&p) {
            Ok(f) => {
                if f.metadata().map(|m| m.is_dir()).unwrap_or(false) {
                    return false;
                }
                self.close(u);
                self.units.insert(
                    u,
                    Unit::Fmt(FmtUnit {
                        path: p,
                        mode: Mode::Read(BufReader::with_capacity(1 << 16, f)),
                        nread: 0,
                    }),
                );
                true
            }
            Err(_) => false,
        }
    }

    /// OPEN(u, STATUS='SCRATCH') (formatted).
    pub fn open_fscratch(&mut self, u: i32) {
        self.close(u);
        self.units.insert(u, Unit::FScr(FmtScratch { lines: Vec::new(), pos: 0 }));
    }

    /// OPEN(u, STATUS='SCRATCH', FORM='UNFORMATTED').
    pub fn open_scratch(&mut self, u: i32) {
        self.close(u);
        self.units.insert(u, Unit::Scr(Scratch { recs: Vec::new(), pos: 0 }));
    }

    pub fn is_open(&self, u: i32) -> bool {
        self.units.contains_key(&u)
    }

    pub fn close(&mut self, u: i32) {
        if let Some(Unit::Fmt(mut f)) = self.units.remove(&u) {
            if let Mode::Write(w) = &mut f.mode {
                let _ = w.flush();
            }
        }
    }

    pub fn close_all(&mut self) {
        let keys: Vec<i32> = self.units.keys().copied().collect();
        for k in keys {
            self.close(k);
        }
    }

    fn fmt(&mut self, u: i32) -> &mut FmtUnit {
        if !self.units.contains_key(&u) {
            // An unopened unit: gfortran uses the file fort.N.
            self.open(u, &format!("fort.{}", u));
        }
        match self.units.get_mut(&u) {
            Some(Unit::Fmt(f)) => f,
            _ => runtime_error(u, "", "Formatted I/O on an unformatted unit"),
        }
    }

    /// Read one formatted record.
    pub fn read(&mut self, u: i32) -> Rd {
        if let Some(Unit::FScr(s)) = self.units.get_mut(&u) {
            if s.pos >= s.lines.len() {
                return Rd::Eof;
            }
            s.pos += 1;
            return Rd::Line(s.lines[s.pos - 1].clone());
        }
        let f = self.fmt(u);
        if let Mode::Fresh = f.mode {
            match File::open(&f.path) {
                Ok(h) => f.mode = Mode::Read(BufReader::with_capacity(1 << 16, h)),
                Err(_) => {
                    // gfortran creates the (empty) file: the read is at EOF.
                    let _ = OpenOptions::new().write(true).create(true).truncate(false).open(&f.path);
                    f.mode = Mode::Eof;
                    return Rd::Eof;
                }
            }
        }
        match &mut f.mode {
            Mode::Read(r) => {
                let mut buf = Vec::with_capacity(96);
                match r.read_until(b'\n', &mut buf) {
                    Ok(0) => {
                        f.mode = Mode::Eof;
                        Rd::Eof
                    }
                    Ok(nb) => {
                        f.nread += nb as u64;
                        if buf.last() == Some(&b'\n') {
                            buf.pop();
                            if buf.last() == Some(&b'\r') {
                                buf.pop();
                            }
                        }
                        Rd::Line(buf)
                    }
                    Err(_) => {
                        let p = f.path.clone();
                        runtime_error(u, &p, "Read error")
                    }
                }
            }
            Mode::Eof => {
                let p = f.path.clone();
                runtime_error(u, &p, "Sequential READ or WRITE not allowed after EOF marker, possibly use REWIND or BACKSPACE")
            }
            Mode::Write(_) => Rd::Eof,
            Mode::Fresh => unreachable!(),
        }
    }

    /// Read one record that must exist (a READ without END=): at end of
    /// file this is the gfortran runtime error "End of file".
    pub fn read_req(&mut self, u: i32) -> Vec<u8> {
        match self.read(u) {
            Rd::Line(l) => l,
            Rd::Eof => {
                let p = match self.units.get(&u) {
                    Some(Unit::Fmt(f)) => f.path.clone(),
                    _ => String::new(),
                };
                runtime_error(u, &p, "End of file")
            }
        }
    }

    /// Write one or more lines (separated by '\n', no final '\n').
    pub fn write<S: AsRef<[u8]>>(&mut self, u: i32, s: S) {
        // preconnected units: 0 = standard error, 6 = standard output
        if u == 0 && !self.units.contains_key(&0) {
            flush_stdout();
            let mut e = std::io::stderr();
            let _ = e.write_all(s.as_ref());
            let _ = e.write_all(b"\n");
            return;
        }
        if u == 6 && !self.units.contains_key(&6) {
            out(s);
            return;
        }
        if let Some(Unit::FScr(sc)) = self.units.get_mut(&u) {
            sc.lines.truncate(sc.pos);
            // a write of several records ('/' in the format)
            for l in s.as_ref().split(|&c| c == b'\n') {
                sc.lines.push(l.to_vec());
            }
            sc.pos = sc.lines.len();
            return;
        }
        let f = self.fmt(u);
        if !matches!(f.mode, Mode::Write(_)) {
            // position of the write: 0 for a unit just opened or rewound,
            // else after the records read
            let at = match f.mode {
                Mode::Fresh => 0,
                _ => f.nread,
            };
            let h = match OpenOptions::new().write(true).create(true).truncate(at == 0).open(&f.path) {
                Ok(h) => h,
                Err(_) => {
                    let p = f.path.clone();
                    runtime_error(u, &p, "Cannot open file")
                }
            };
            if at > 0 {
                use std::io::{Seek, SeekFrom};
                let _ = h.set_len(at);
                let mut h2 = h;
                let _ = h2.seek(SeekFrom::Start(at));
                f.mode = Mode::Write(BufWriter::with_capacity(1 << 16, h2));
            } else {
                f.mode = Mode::Write(BufWriter::with_capacity(1 << 16, h));
            }
        }
        if let Mode::Write(w) = &mut f.mode {
            let _ = w.write_all(s.as_ref());
            let _ = w.write_all(b"\n");
        }
    }

    /// REWIND of a formatted unit (reopens it for reading).
    pub fn rewind(&mut self, u: i32) {
        match self.units.get_mut(&u) {
            Some(Unit::Fmt(f)) => {
                if let Mode::Write(w) = &mut f.mode {
                    let _ = w.flush();
                }
                f.mode = Mode::Fresh;
                f.nread = 0;
            }
            Some(Unit::Scr(s)) => s.pos = 0,
            Some(Unit::FScr(s)) => s.pos = 0,
            None => {}
        }
    }

    // ---------------- unformatted scratch units ----------------

    fn scr(&mut self, u: i32) -> &mut Scratch {
        if !self.units.contains_key(&u) {
            self.open_scratch(u);
        }
        match self.units.get_mut(&u) {
            Some(Unit::Scr(s)) => s,
            _ => runtime_error(u, "", "Unformatted I/O on a formatted unit"),
        }
    }

    /// Unformatted WRITE: the record replaces everything after the
    /// current position (sequential files are truncated after a write).
    pub fn uwrite(&mut self, u: i32, rec: Vec<Val>) {
        let s = self.scr(u);
        s.recs.truncate(s.pos);
        s.recs.push(rec);
        s.pos += 1;
    }

    /// Unformatted READ; None at end of file.
    pub fn uread(&mut self, u: i32) -> Option<&[Val]> {
        let s = self.scr(u);
        if s.pos >= s.recs.len() {
            return None;
        }
        s.pos += 1;
        Some(&s.recs[s.pos - 1])
    }

    /// Unformatted READ that must succeed.
    pub fn uread_req(&mut self, u: i32) -> Vec<Val> {
        match self.uread(u) {
            Some(r) => r.to_vec(),
            None => runtime_error(u, "", "End of file"),
        }
    }

    /// Unformatted READ of an empty record list: skip one record.
    pub fn uskip(&mut self, u: i32) -> bool {
        self.uread(u).is_some()
    }

    /// END FILE: the file ends at the current position.
    pub fn endfile(&mut self, u: i32) {
        if let Some(Unit::FScr(s)) = self.units.get_mut(&u) {
            s.lines.truncate(s.pos);
            return;
        }
        if let Some(Unit::Fmt(f)) = self.units.get_mut(&u) {
            match &mut f.mode {
                Mode::Write(w) => {
                    let _ = w.flush();
                }
                Mode::Read(_) | Mode::Eof => {
                    if let Ok(h) = OpenOptions::new().write(true).open(&f.path) {
                        let _ = h.set_len(f.nread);
                    }
                    f.mode = Mode::Eof;
                }
                Mode::Fresh => {
                    let _ = OpenOptions::new().write(true).create(true).truncate(true).open(&f.path);
                    f.mode = Mode::Eof;
                }
            }
            return;
        }
        let s = self.scr(u);
        s.recs.truncate(s.pos);
    }

    pub fn backspace(&mut self, u: i32) {
        if let Some(Unit::Scr(s)) = self.units.get_mut(&u) {
            if s.pos > 0 {
                s.pos -= 1;
            }
        }
    }
}

/// A formatted record padded with blanks to at least n bytes.
pub fn pad(mut l: Vec<u8>, n: usize) -> Vec<u8> {
    if l.len() < n {
        l.resize(n, b' ');
    }
    l
}

/// Fortran `Iw` input of bytes [s, s+w) of a padded record; a bad field is
/// the gfortran runtime error.
pub fn rd_i(rec: &[u8], s: usize, w: usize, unit: i32) -> i32 {
    let mut f = [b' '; 64];
    for k in 0..w {
        if s + k < rec.len() {
            f[k] = rec[s + k];
        }
    }
    match crate::ffmt::read_i(&f[..w]) {
        Ok(v) => v,
        Err(_) => runtime_error(unit, "", "Bad value during integer read"),
    }
}

/// Fortran `Fw.d`/`Ew.d` input of bytes [s, s+w) of a padded record.
pub fn rd_f(rec: &[u8], s: usize, w: usize, d: usize, unit: i32) -> f64 {
    let mut f = [b' '; 64];
    for k in 0..w {
        if s + k < rec.len() {
            f[k] = rec[s + k];
        }
    }
    match crate::ffmt::read_f(&f[..w], d) {
        Ok(v) => v,
        Err(_) => runtime_error(unit, "", "Bad value during floating point read"),
    }
}

/// Bytes [s, s+w) of a record as text (blank padded).
pub fn rd_a(rec: &[u8], s: usize, w: usize) -> Vec<u8> {
    let mut v = vec![b' '; w];
    for k in 0..w {
        if s + k < rec.len() {
            v[k] = rec[s + k];
        }
    }
    v
}

/// A line built from text and raw record bytes.
#[derive(Default, Clone)]
pub struct L(pub Vec<u8>);

impl L {
    pub fn new() -> L {
        L(Vec::with_capacity(96))
    }
    pub fn s(mut self, t: &str) -> L {
        self.0.extend_from_slice(t.as_bytes());
        self
    }
    pub fn b(mut self, t: &[u8]) -> L {
        self.0.extend_from_slice(t);
        self
    }
    pub fn push(&mut self, t: &str) {
        self.0.extend_from_slice(t.as_bytes());
    }
    pub fn pushb(&mut self, t: &[u8]) {
        self.0.extend_from_slice(t);
    }
}

impl AsRef<[u8]> for L {
    fn as_ref(&self) -> &[u8] {
        &self.0
    }
}
