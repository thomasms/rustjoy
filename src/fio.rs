//! Emulation of Fortran sequential I/O units as gfortran implements them.
//!
//! Every unit is held in memory as a list of records plus a position, so the
//! Fortran semantics (rewind, backspace, "a write truncates the file after the
//! current record") are exact. Disk-backed units (`tapeNN`) are loaded when
//! opened and written back when closed or when the program ends. Unformatted
//! files use gfortran's on-disk layout: a 4-byte little-endian length marker
//! before and after each record.

#![allow(dead_code)]

use std::collections::HashMap;
use std::fs;
use std::io::Write;

use crate::ffmt;

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Status {
    Old,
    New,
    Scratch,
}

pub struct Unit {
    pub formatted: bool,
    pub recs: Vec<Vec<u8>>,
    pub pos: usize,
    pub path: Option<String>,
    pub dirty: bool,
}

impl Unit {
    fn load(path: &str, formatted: bool) -> std::io::Result<Vec<Vec<u8>>> {
        let bytes = fs::read(path)?;
        let mut recs = Vec::new();
        if formatted {
            let mut start = 0;
            for (i, b) in bytes.iter().enumerate() {
                if *b == b'\n' {
                    // libgfortran treats CR-LF as a record terminator
                    let end = if i > start && bytes[i - 1] == b'\r' { i - 1 } else { i };
                    recs.push(bytes[start..end].to_vec());
                    start = i + 1;
                }
            }
            if start < bytes.len() {
                recs.push(bytes[start..].to_vec());
            }
        } else {
            let mut off = 0;
            while off + 4 <= bytes.len() {
                let n = i32::from_le_bytes(bytes[off..off + 4].try_into().unwrap());
                let n = n.unsigned_abs() as usize;
                let rec = bytes[off + 4..off + 4 + n].to_vec();
                recs.push(rec);
                off += 8 + n;
            }
        }
        Ok(recs)
    }

    fn save(&self) {
        if let Some(p) = &self.path {
            let mut out = Vec::new();
            for r in &self.recs {
                if self.formatted {
                    out.extend_from_slice(r);
                    out.push(b'\n');
                } else {
                    let n = (r.len() as i32).to_le_bytes();
                    out.extend_from_slice(&n);
                    out.extend_from_slice(r);
                    out.extend_from_slice(&n);
                }
            }
            if let Err(e) = fs::write(p, out) {
                eprintln!("cannot write {}: {}", p, e);
            }
        }
    }
}

/// Builder for one unformatted record (default integers are 4 bytes).
#[derive(Default)]
pub struct RecW {
    pub b: Vec<u8>,
}

impl RecW {
    pub fn new() -> Self {
        RecW { b: Vec::new() }
    }
    pub fn i(&mut self, v: i32) -> &mut Self {
        self.b.extend_from_slice(&v.to_le_bytes());
        self
    }
    pub fn r(&mut self, v: f64) -> &mut Self {
        self.b.extend_from_slice(&v.to_le_bytes());
        self
    }
    pub fn rs(&mut self, v: &[f64]) -> &mut Self {
        for x in v {
            self.b.extend_from_slice(&x.to_le_bytes());
        }
        self
    }
}

/// Reader over one unformatted record.
pub struct RecR<'a> {
    pub b: &'a [u8],
    pub off: usize,
}

impl<'a> RecR<'a> {
    pub fn i(&mut self) -> i32 {
        if self.off + 4 > self.b.len() {
            fatal_runtime("Fortran runtime error: I/O past end of record on unformatted file");
        }
        let v = i32::from_le_bytes(self.b[self.off..self.off + 4].try_into().unwrap());
        self.off += 4;
        v
    }
    pub fn r(&mut self) -> f64 {
        if self.off + 8 > self.b.len() {
            fatal_runtime("Fortran runtime error: I/O past end of record on unformatted file");
        }
        let v = f64::from_le_bytes(self.b[self.off..self.off + 8].try_into().unwrap());
        self.off += 8;
        v
    }
}

/// A run-time error of the Fortran library: gfortran prints a message and
/// exits with status 2.
pub fn fatal_runtime(msg: &str) -> ! {
    crate::fio::flush_all_global();
    eprintln!("{}", msg);
    std::process::exit(2);
}

// The I/O system is a process-wide singleton, like the Fortran runtime.
static mut IO_PTR: *mut Io = std::ptr::null_mut();

pub fn flush_all_global() {
    unsafe {
        if !IO_PTR.is_null() {
            (*IO_PTR).flush_all();
        }
    }
}

pub struct Io {
    pub units: HashMap<i32, Unit>,
    pub out: Vec<u8>,
    pub out_path: String,
    pub input: ListInput,
    /// worker mode (see `detached`): nothing is written anywhere
    pub quiet: bool,
    /// set in worker mode when anything tried to write
    pub touched: bool,
}

/// Payload of the unwind that `Ctx::error` starts in worker mode.
pub struct WorkerAbort;

impl Io {
    pub fn new(input_text: Vec<u8>, out_path: &str) -> Box<Io> {
        let mut io = Box::new(Io {
            units: HashMap::new(),
            out: Vec::new(),
            out_path: out_path.to_string(),
            input: ListInput::new(input_text),
            quiet: false,
            touched: false,
        });
        unsafe {
            IO_PTR = &mut *io as *mut Io;
        }
        io
    }

    /// An `Io` for a worker thread that computes values ahead of the main
    /// program: it has no units, writes nothing (listing, terminal, files)
    /// and does not register itself for the global flush. Any attempt to
    /// write sets `touched`, so the worker knows its result had a side
    /// effect and must not be used in place of the real call.
    pub fn detached() -> Box<Io> {
        Box::new(Io {
            units: HashMap::new(),
            out: Vec::new(),
            out_path: String::new(),
            input: ListInput::new(Vec::new()),
            quiet: true,
            touched: false,
        })
    }

    /// Write one or more records to the listing (nsyso). `s` may contain '\n'
    /// for Fortran '/' record breaks; a final newline is added.
    pub fn so(&mut self, s: &str) {
        if self.quiet {
            self.touched = true;
            return;
        }
        self.out.extend_from_slice(s.as_bytes());
        self.out.push(b'\n');
    }

    /// Write to the terminal (nsyse).
    pub fn se(&mut self, s: &str) {
        if self.quiet {
            self.touched = true;
            return;
        }
        let stdout = std::io::stdout();
        let mut h = stdout.lock();
        let _ = h.write_all(s.as_bytes());
        let _ = h.write_all(b"\n");
    }

    pub fn flush_all(&mut self) {
        if self.quiet {
            return;
        }
        let _ = fs::write(&self.out_path, &self.out);
        for (_, u) in self.units.iter_mut() {
            if u.dirty {
                u.save();
                u.dirty = false;
            }
        }
        let _ = std::io::stdout().flush();
    }

    pub fn open(&mut self, nun: i32, formatted: bool, status: Status, path: Option<&str>) {
        if self.units.contains_key(&nun) {
            // Re-opening a connected unit without FILE= refers to the same
            // file: gfortran leaves the connection as it is.
            if path.is_none() || self.units[&nun].path.as_deref() == path {
                return;
            }
            self.close(nun, false);
        }
        let mut recs = Vec::new();
        let mut dirty = false;
        match status {
            Status::Old => {
                let p = path.unwrap();
                match Unit::load(p, formatted) {
                    Ok(r) => recs = r,
                    Err(_) => fatal_runtime(&format!(
                        "Fortran runtime error: Cannot open file '{}': No such file or directory",
                        p
                    )),
                }
            }
            Status::New => {
                dirty = true;
            }
            Status::Scratch => {}
        }
        let path = if status == Status::Scratch { None } else { path.map(|s| s.to_string()) };
        self.units.insert(nun, Unit { formatted, recs, pos: 0, path, dirty });
    }

    fn unit(&mut self, nun: i32, formatted: bool) -> &mut Unit {
        if !self.units.contains_key(&nun) {
            // implicit open, as gfortran does with fort.N
            let path = format!("fort.{}", nun);
            let recs = Unit::load(&path, formatted).unwrap_or_default();
            self.units.insert(nun, Unit { formatted, recs, pos: 0, path: Some(path), dirty: false });
        }
        let u = self.units.get_mut(&nun).unwrap();
        if u.formatted != formatted {
            fatal_runtime(&format!(
                "Fortran runtime error: {} I/O on {} unit {}",
                if formatted { "Formatted" } else { "Unformatted" },
                if u.formatted { "formatted" } else { "unformatted" },
                nun
            ));
        }
        u
    }

    pub fn rewind(&mut self, nun: i32) {
        if let Some(u) = self.units.get_mut(&nun) {
            u.pos = 0;
        }
    }

    pub fn backspace(&mut self, nun: i32) {
        if let Some(u) = self.units.get_mut(&nun) {
            if u.pos > 0 {
                u.pos -= 1;
            }
        }
    }

    pub fn close(&mut self, nun: i32, _delete: bool) {
        if let Some(u) = self.units.remove(&nun) {
            if u.dirty {
                u.save();
            }
        }
    }

    pub fn is_open(&self, nun: i32) -> bool {
        self.units.contains_key(&nun)
    }

    /// Read one formatted record.
    pub fn read_line(&mut self, nun: i32) -> Vec<u8> {
        let u = self.unit(nun, true);
        if u.pos >= u.recs.len() {
            fatal_runtime(&format!("Fortran runtime error: End of file (unit {})", nun));
        }
        let r = u.recs[u.pos].clone();
        u.pos += 1;
        r
    }

    /// Write one formatted record.
    pub fn write_line(&mut self, nun: i32, s: &str) {
        let u = self.unit(nun, true);
        u.recs.truncate(u.pos);
        u.recs.push(s.as_bytes().to_vec());
        u.pos += 1;
        u.dirty = true;
    }

    /// Read one unformatted record; the caller decodes it.
    pub fn read_rec(&mut self, nun: i32) -> Vec<u8> {
        let u = self.unit(nun, false);
        if u.pos >= u.recs.len() {
            fatal_runtime(&format!("Fortran runtime error: End of file (unit {})", nun));
        }
        let r = u.recs[u.pos].clone();
        u.pos += 1;
        r
    }

    /// Like read_rec but returns a borrowed record (avoids a copy).
    pub fn read_rec_with<T>(&mut self, nun: i32, f: impl FnOnce(&mut RecR) -> T) -> T {
        let u = self.unit(nun, false);
        if u.pos >= u.recs.len() {
            fatal_runtime(&format!("Fortran runtime error: End of file (unit {})", nun));
        }
        let p = u.pos;
        u.pos += 1;
        let mut rr = RecR { b: &u.recs[p], off: 0 };
        f(&mut rr)
    }

    pub fn write_rec(&mut self, nun: i32, rec: RecW) {
        let u = self.unit(nun, false);
        u.recs.truncate(u.pos);
        u.recs.push(rec.b);
        u.pos += 1;
        u.dirty = true;
    }
}

// --------------------------------------------------- list-directed input --

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Kind {
    Int,
    Real,
    Char(usize),
}

#[derive(Clone, Debug)]
pub enum Val {
    Int(i32),
    Real(f64),
    Char(String),
}

impl Val {
    pub fn int(&self) -> i32 {
        match self {
            Val::Int(i) => *i,
            _ => unreachable!(),
        }
    }
    pub fn real(&self) -> f64 {
        match self {
            Val::Real(x) => *x,
            _ => unreachable!(),
        }
    }
    pub fn chars(&self) -> String {
        match self {
            Val::Char(s) => s.clone(),
            _ => unreachable!(),
        }
    }
}

/// gfortran list-directed (`read(nsysi,*)`) input.
pub struct ListInput {
    lines: Vec<Vec<u8>>,
    rec: usize,
}

enum Tok {
    Value(Vec<u8>, bool), // text, quoted
    Null,
    Slash,
}

impl ListInput {
    pub fn new(text: Vec<u8>) -> Self {
        let mut lines = Vec::new();
        let mut start = 0;
        for (i, b) in text.iter().enumerate() {
            if *b == b'\n' {
                lines.push(text[start..i].to_vec());
                start = i + 1;
            }
        }
        if start < text.len() {
            lines.push(text[start..].to_vec());
        }
        ListInput { lines, rec: 0 }
    }

    /// Perform one list-directed READ for the given item kinds. Items that
    /// receive no value (null values, or a terminating slash) are `None` and
    /// must be left unchanged by the caller.
    pub fn read(&mut self, kinds: &[Kind]) -> Vec<Option<Val>> {
        let mut out: Vec<Option<Val>> = Vec::with_capacity(kinds.len());
        if self.rec >= self.lines.len() {
            fatal_runtime("Fortran runtime error: End of file (standard input)");
        }
        let mut line = self.rec;
        let mut p = 0usize;
        let mut first = true; // at start of the READ (a leading comma is a null)
        let mut after_comma = false;
        let mut pending: Option<(usize, Option<(Vec<u8>, bool)>)> = None; // repeat
        let mut done = false;
        while out.len() < kinds.len() {
            if let Some((n, v)) = pending.take() {
                let k = kinds[out.len()];
                out.push(match &v {
                    None => None,
                    Some((t, q)) => Some(Self::convert(t, *q, k)),
                });
                if n > 1 {
                    pending = Some((n - 1, v));
                }
                continue;
            }
            let kind = kinds[out.len()];
            // get next token
            let tok = self.next_token(&mut line, &mut p, &mut first, &mut after_comma, kind);
            match tok {
                Tok::Slash => {
                    done = true;
                    break;
                }
                Tok::Null => out.push(None),
                Tok::Value(t, q) => {
                    // repeat count?
                    if !q {
                        if let Some(star) = t.iter().position(|c| *c == b'*') {
                            let (r, rest) = (&t[..star], &t[star + 1..]);
                            if !r.is_empty() && r.iter().all(|c| c.is_ascii_digit()) {
                                let n: usize = String::from_utf8_lossy(r).parse().unwrap();
                                let v = if rest.is_empty() {
                                    None
                                } else {
                                    let (rt, rq) = if rest[0] == b'\'' || rest[0] == b'"' {
                                        (rest[1..rest.len() - 1].to_vec(), true)
                                    } else {
                                        (rest.to_vec(), false)
                                    };
                                    Some((rt, rq))
                                };
                                pending = Some((n, v));
                                continue;
                            }
                        }
                    }
                    out.push(Some(Self::convert(&t, q, kind)));
                }
            }
        }
        let _ = done;
        while out.len() < kinds.len() {
            out.push(None);
        }
        // the rest of the current record is skipped
        self.rec = line + 1;
        out
    }

    fn convert(t: &[u8], quoted: bool, k: Kind) -> Val {
        let s = String::from_utf8_lossy(t).to_string();
        match k {
            Kind::Int => {
                if quoted {
                    fatal_runtime("Fortran runtime error: Bad integer for item in list input");
                }
                match ffmt::read_i(t) {
                    Ok(v) if !s.contains('.') => Val::Int(v),
                    _ => fatal_runtime(&format!(
                        "Fortran runtime error: Bad integer for item in list input: '{}'",
                        s
                    )),
                }
            }
            Kind::Real => {
                if quoted {
                    fatal_runtime("Fortran runtime error: Bad real number in item of list input");
                }
                match ffmt::read_f(t, 0) {
                    Ok(v) => Val::Real(v),
                    Err(_) => fatal_runtime(&format!(
                        "Fortran runtime error: Bad real number in item of list input: '{}'",
                        s
                    )),
                }
            }
            Kind::Char(n) => {
                let mut v: Vec<u8> = t.to_vec();
                v.truncate(n);
                while v.len() < n {
                    v.push(b' ');
                }
                Val::Char(String::from_utf8_lossy(&v).to_string())
            }
        }
    }

    fn next_token(
        &mut self,
        line: &mut usize,
        p: &mut usize,
        first: &mut bool,
        after_comma: &mut bool,
        _kind: Kind,
    ) -> Tok {
        loop {
            if *line >= self.lines.len() {
                fatal_runtime("Fortran runtime error: End of file (standard input)");
            }
            let l = &self.lines[*line];
            // skip blanks
            while *p < l.len() && (l[*p] == b' ' || l[*p] == b'\t' || l[*p] == b'\r') {
                *p += 1;
            }
            if *p >= l.len() {
                // end of record acts as a blank
                *line += 1;
                *p = 0;
                continue;
            }
            let c = l[*p];
            if c == b',' {
                *p += 1;
                if *first || *after_comma {
                    *first = false;
                    *after_comma = true;
                    return Tok::Null;
                }
                *after_comma = true;
                continue;
            }
            if c == b'/' {
                *p += 1;
                return Tok::Slash;
            }
            *first = false;
            *after_comma = false;
            if c == b'\'' || c == b'"' {
                // quoted string; may continue on following records
                let q = c;
                *p += 1;
                let mut s = Vec::new();
                loop {
                    let l = &self.lines[*line];
                    if *p >= l.len() {
                        *line += 1;
                        *p = 0;
                        if *line >= self.lines.len() {
                            fatal_runtime("Fortran runtime error: End of file (standard input)");
                        }
                        continue;
                    }
                    let ch = l[*p];
                    if ch == q {
                        if *p + 1 < l.len() && l[*p + 1] == q {
                            s.push(q);
                            *p += 2;
                            continue;
                        }
                        *p += 1;
                        break;
                    }
                    s.push(ch);
                    *p += 1;
                }
                self.after_value(line, p, after_comma);
                return Tok::Value(s, true);
            }
            // unquoted value: up to blank, comma, slash or end of record
            let l = &self.lines[*line];
            let start = *p;
            while *p < l.len() {
                let ch = l[*p];
                if ch == b' ' || ch == b'\t' || ch == b',' || ch == b'/' || ch == b'\r' {
                    break;
                }
                // a repeat count followed by a quoted string, e.g. 3*'abc'
                if (ch == b'\'' || ch == b'"') && *p > start && l[*p - 1] == b'*' {
                    let q = ch;
                    *p += 1;
                    while *p < l.len() && l[*p] != q {
                        *p += 1;
                    }
                    if *p < l.len() {
                        *p += 1;
                    }
                    break;
                }
                *p += 1;
            }
            let t = l[start..*p].to_vec();
            self.after_value(line, p, after_comma);
            return Tok::Value(t, false);
        }
    }

    /// Consume a value separator after a value: blanks and at most one comma.
    fn after_value(&mut self, line: &mut usize, p: &mut usize, after_comma: &mut bool) {
        let l = &self.lines[*line];
        while *p < l.len() && (l[*p] == b' ' || l[*p] == b'\t' || l[*p] == b'\r') {
            *p += 1;
        }
        if *p < l.len() && l[*p] == b',' {
            *p += 1;
            *after_comma = true;
        } else {
            *after_comma = false;
        }
    }
}
