# Examples

Each directory here is an NJOY input deck (`input`) with the tapes that it
reads (`tapeNN`). The tapes are links to the ENDF files in [`endf/`](endf).
Together, the examples use all 14 modules that rustjoy ports.

| Example | Modules | Tapes | From |
| --- | --- | --- | --- |
| `h2-reconr` | reconr | tape20: H-2 (ENDF/B-VIII.0) | NJOY2016 test 84 |
| `sn119-heatr` | moder, reconr, broadr, heatr | tape20: Sn-119 (ENDF/B-VIII.0) | NJOY2016 test 79 |
| `sn119-purr` | moder, reconr, broadr, purr | tape20: Sn-119 (ENDF/B-VIII.0) | NJOY2016 test 36, without acer and viewr |
| `sn119-resxsr` | reconr, broadr, resxsr | tape20: Sn-119 (ENDF/B-VIII.0) | new (no NJOY2016 test calls resxsr) |
| `pu241-unresr` | moder, reconr, broadr, unresr, heatr, purr | tape20: Pu-241 (ENDF/B-VI) | NJOY2016 test 19, without acer |
| `fe-groupr` | moder, reconr, broadr, groupr | tape20: Fe (IRDFF-II) | NJOY2016 test 60 |
| `b10-gaspr` | moder, reconr, broadr, gaspr | tape20: B-10 (ENDF/B-VIII.0) | NJOY2016 test 45 |
| `mixr-h2-b10` | reconr, mixr | tape20: H-2, tape30: B-10 (ENDF/B-VIII.0) | new (no NJOY2016 test calls mixr) |
| `zrh-thermr` | moder, reconr, broadr, thermr | tape20: Zr-90, tape30: Zr in ZrH (ENDF/B-VIII.0) | NJOY2016 test 32 |
| `parah-leapr` | leapr | none | NJOY2016 test 22 |
| `photon-gaminr` | reconr, gaminr, dtfr | tape30, tape32: photon data (DLC-7E) | NJOY2016 test 03, without matxsr and viewr |

The ENDF files are copies of files in the NJOY2016 test resources
(`tests/resources`, commit `ab27c64`).

## Run an example

rustjoy reads the deck from standard input and writes the listing to
`output` and the tapes to `tapeNN` files in the current directory. To keep
this directory clean, run an example in a new directory (`cp -L` copies
the files that the links point to):

```bash
cargo build --release
mkdir /tmp/sn119-heatr
cp -L examples/sn119-heatr/* /tmp/sn119-heatr/
cd /tmp/sn119-heatr && /path/to/rustjoy/target/release/rustjoy < input
```

Run NJOY2016 the same way (`njoy < input`) in another directory to get its
results.

## Comparison with NJOY2016

`tests/examples.rs` runs every example with rustjoy and compares the
results with NJOY2016 byte for byte: exit status, standard output, the
listing and every tape that the deck writes. See the
[main README](../README.md#tests) for how to run the tests.

`reference.txt` in each directory holds the NJOY2016 results: the exit
status, and the length and SHA-256 checksum of each output file after the
masks are applied. Only bytes that differ between two runs of NJOY2016 are
masked: CPU times, dates and times, and some bytes of binary tapes that
NJOY2016 fills from uninitialised memory. The checksums replace the
outputs themselves, which are about 180 MB in total.

To add an example, make a new directory with `input` and links to the
tapes, add a line to the `examples!` list in `tests/examples.rs`, and make
its `reference.txt` with NJOY2016:

```bash
NJOY=/path/to/njoy RUSTJOY_UPDATE_REFERENCES=1 cargo test --release --test examples
```

Only the modules that rustjoy ports can be in the deck.
