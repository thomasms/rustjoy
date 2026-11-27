//! `njoy-reconr`: a Rust port of the MODER, RECONR, BROADR, HEATR, THERMR,
//! UNRESR, GASPR, PURR, LEAPR, MIXR, RESXSR and DTFR modules of NJOY2016.
//!
//! Usage is the same as `njoy`: the input deck is read from standard input,
//! the listing is written to the file `output`, terminal messages go to
//! standard output, and ENDF/PENDF tapes are `tapeNN` files in the current
//! directory. Only the `moder`, `reconr`, `broadr`, `heatr`, `thermr`,
//! `unresr`, `gaspr`, `purr`, `leapr`, `mixr`, `resxsr`,
//! `dtfr` and `gaminr` modules (and `--` comment cards)
//! are supported.

use reconr_rs::{broadr, dtfr, endf, ffmt, fio, gaminr, gaspr, heatr, leapr, mixr, moder, purr, reconr, resxsr, thermr, unresr, util};

use std::io::Read;

const VERS: &str = "2016.80 ";
const VDAY: &str = "08Oct26 ";
const LAB: &str = "lanl t2 ";
const MX: &str = "        ";

fn main() {
    let mut input = Vec::new();
    std::io::stdin().read_to_end(&mut input).expect("reading standard input");
    let io = fio::Io::new(input, "output");
    let mut ctx = endf::Ctx::new(io);
    ctx.npage = (ctx.npage / 102) * 102;

    let date = util::dater();
    let time = util::wclock();
    let star77 = "*".repeat(77);
    let hdr = format!(
        "1{s}\n \
         *{b41}*{b14}*{b18}*\n \
         *   $$   $$       $$   $$$$$   $$    $$   *{b14}*  vers: {v}  *\n \
         *   $$$  $$       $$  $$$$$$$   $$  $$    *  nuclear     *  vday: {vd}  *\n \
         *   $$$$ $$       $$  $$   $$    $$$$     *  data        *  site: {l}  *\n \
         *   $$ $$$$  $$   $$  $$   $$     $$      *  processing  *  mach: {m}  *\n \
         *   $$  $$$  $$$$$$$  $$$$$$$     $$      *  system      *  date: {d}  *\n \
         *   $$   $$   $$$$$    $$$$$      $$      *{b14}*  time: {t}  *\n \
         *{b41}*{b14}*{b18}*\n \
         {s}",
        s = star77,
        b41 = " ".repeat(41),
        b14 = " ".repeat(14),
        b18 = " ".repeat(18),
        v = VERS,
        vd = VDAY,
        l = LAB,
        m = MX,
        d = date,
        t = time
    );
    ctx.io.so(&hdr);
    ctx.io.se(&format!("\n njoy {} {}{}{} {}\n {}", VERS, VDAY, " ".repeat(38), date, time, star77));

    let mut rc = reconr::Reconr::new();
    let mut br = broadr::Broadr::new();
    let mut ht = heatr::Heatr::new();
    let mut md = moder::Moder::new();
    let mut tm = thermr::Thermr::new();
    let mut ur = unresr::Unresr::new();
    let mut pu = purr::Purr::new();
    let mut lp = leapr::Leapr::new();
    let mut mx = mixr::Mixr::new();
    let mut rx = resxsr::Resxsr::new();
    let mut dt = dtfr::Dtfr::new();
    let mut gm = gaminr::Gaminr::new();
    loop {
        let v = ctx.io.input.read(&[fio::Kind::Char(6)]);
        // a null value leaves the variable unchanged (blank on first read)
        let module = match &v[0] {
            Some(x) => x.chars(),
            None => "      ".to_string(),
        };
        let m = module.trim_end();
        if m == "stop" {
            break;
        }
        match m {
            "reconr" => {
                rc.reconr(&mut ctx);
            }
            "broadr" => {
                br.broadr(&mut ctx);
            }
            "heatr" => {
                ht.heatr(&mut ctx);
            }
            "moder" => {
                md.moder(&mut ctx);
            }
            "thermr" => {
                tm.thermr(&mut ctx);
            }
            "unresr" => {
                ur.unresr(&mut ctx);
            }
            "gaspr" => {
                gaspr::gaspr(&mut ctx);
            }
            "purr" => {
                pu.purr(&mut ctx);
            }
            "leapr" => {
                lp.leapr(&mut ctx);
            }
            "mixr" => {
                mx.mixr(&mut ctx);
            }
            "resxsr" => {
                rx.resxsr(&mut ctx);
            }
            "dtfr" => {
                dt.dtfr(&mut ctx);
            }
            "gaminr" => {
                gm.gaminr(&mut ctx);
            }
            "--" => {}
            _ => {
                let strng = format!("illegal module name---{}", ffmt::a(&module, 6));
                ctx.error("njoy", &strng, " ");
            }
        }
    }
    let secs = ctx.timer();
    ctx.io.se(&format!("{}{}s\n {}", " ".repeat(69), ffmt::f(secs, 8, 1), star77));
    ctx.io.flush_all();
}
