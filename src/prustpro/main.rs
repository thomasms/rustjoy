//! prustpro: a Rust port of the PREPRO 2026-1 ENDF/B pre-processing codes
//! (IAEA Nuclear Data Section, D.E. Cullen).
//!
//! Each PREPRO program is a subcommand: `prustpro linear` does what the
//! Fortran program `linear` does, with the same input files, output files,
//! screen output and exit status. If the executable itself is called by a
//! program name (for example through a link named `linear`), it runs that
//! program, so the links can replace the Fortran executables.

#[path = "../common/ffmt.rs"]
mod ffmt;

mod endf2c;
mod endfio;
mod fio;
mod linear;
mod numcheck;

/// The programs, in the order of the PREPRO verification run.
const PROGRAMS: &[(&str, fn() -> !)] = &[("endf2c", endf2c::run), ("linear", linear::run)];

fn usage() -> ! {
    eprintln!("prustpro {} - Rust port of PREPRO 2026-1", env!("CARGO_PKG_VERSION"));
    eprintln!();
    eprintln!("usage: prustpro PROGRAM");
    eprintln!("   or: call the executable by a program name (a link named PROGRAM)");
    eprintln!();
    eprintln!("Each program reads and writes the same files in the current");
    eprintln!("directory as the PREPRO program of the same name.");
    eprintln!();
    eprintln!("programs:");
    for (n, _) in PROGRAMS {
        eprintln!("  {}", n);
    }
    std::process::exit(1);
}

fn find(name: &str) -> Option<fn() -> !> {
    let n = name.to_ascii_lowercase();
    let n = n.strip_suffix(".exe").unwrap_or(&n);
    PROGRAMS.iter().find(|(p, _)| *p == n).map(|(_, f)| *f)
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let exe = args
        .first()
        .map(|a| {
            std::path::Path::new(a)
                .file_name()
                .map(|f| f.to_string_lossy().into_owned())
                .unwrap_or_default()
        })
        .unwrap_or_default();
    let (name, prog) = if let Some(f) = find(&exe) {
        (exe.clone(), f)
    } else {
        match args.get(1) {
            Some(a) if a == "-h" || a == "--help" => usage(),
            Some(a) if a == "_numcheck" => numcheck::run(),
            Some(a) => match find(a) {
                Some(f) => (a.clone(), f),
                None => {
                    eprintln!("prustpro: unknown program '{}'", a);
                    usage()
                }
            },
            None => usage(),
        }
    };
    *fio::PROGRAM.lock().unwrap() = format!("{}.f", name.to_ascii_lowercase());
    std::panic::set_hook(Box::new(|info| {
        if info.payload().downcast_ref::<fio::RuntimeError>().is_none() {
            eprintln!("prustpro: internal error: {}", info);
        }
    }));
    match std::panic::catch_unwind(|| prog()) {
        Ok(never) => never,
        Err(e) => {
            fio::flush_stdout();
            if let Some(fio::RuntimeError(m)) = e.downcast_ref::<fio::RuntimeError>() {
                eprintln!("{}\n\nError termination.", m);
                std::process::exit(2);
            }
            std::process::exit(101);
        }
    }
}
