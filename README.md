# sciexwiff

Pure-Rust tooling for SCIEX ZenoTOF data, reverse-engineered from the data files
alone (no Clearcore2 SDK), in two parts:

1. **`.wiff` method reader** — opens the `.wiff` OLE2 / Compound File (via the
   [`cfb`](https://crates.io/crates/cfb) crate) and decodes the acquisition
   *design* a simulator needs to replicate a run:
   - **SWATH isolation scheme** — the variable isolation windows
     (`/MethodSubtree/Method1/DeviceMethod0/SWATHMethod`).
   - **TOF → m/z calibration** coefficients
     (`/SampleSubtree/Sample1/TOFCalibrationData`).
2. **`.wiff.scan` spectra codec** (`wiffscan` module) — decode **and** encode the
   packed peak token stream that holds the actual spectra.

## The `.wiff.scan` spectra codec

An earlier version of this note said the spectra were out of reach. That is no
longer true — the `.wiff.scan` packing is now fully worked out and this crate
encodes as well as decodes it. From inspection of real ZenoTOF 7600 data:

- `.wiff2` is **encrypted / high-entropy** — not open SQLite. It is *not* required
  to read the spectra (ProteoWizard reads a `.wiff` + `.wiff.scan` with no
  `.wiff2`), so it is not needed here.
- `.wiff.scan` is **not encrypted**. It is a sequential stream of units, each
  `[protobuf metadata][ffffffff  u32 hdr  00  token-stream]`. Peaks are stored
  m/z-ascending as `(n, intensity)` with `n` an integer TOF sample index and
  `m/z = (a/5·n + b)²` (`a`,`b` = the two per-scan calibration doubles in the
  metadata). The full token grammar (gap/consecutive deltas + intensity escapes)
  is documented at the top of `src/wiffscan.rs`.
- The codec is **round-trip byte-identical** against real blocks — see the
  `parity_real_wiff_scan` test (gated on `TIMSIM_SCIEX_WIFF_SCAN=<path>`), which
  decodes→re-encodes real scan blocks and asserts the bytes match exactly.

**Oracle / legitimacy.** Everything here was derived by inspecting the *data
file* (the same legal footing as reading any file format), using **ProteoWizard
as a read-only oracle** (`msconvert` under Wine, used as intended) to check
decoded peaks against the vendor's own interpretation. No Clearcore2 binary is
decompiled, disassembled, or linked; the SDK EULA is not touched.

## Reference data

The decoding above was validated against:

> **PRIDE `PXD036786`** — *"High-throughput proteomics of nanogram-scale
> samples with Zeno SWATH DIA"* (ZenoTOF 7600, pub. 2022/12).

No vendor data is committed here; re-download from PRIDE to run. FTP archive:
`https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/12/PXD036786/`. The file
used for development (K562 0.98 ng, standard SWATH, ~53 MB `.wiff2.Zip` =
`.wiff` OLE2 11 MB + `.wiff2` 29 MB + `.wiff.scan` 32 MB):

```
20211020_Z1_ZW_001_30-0043_K562_1ul_K562_0.98ng_1.wiff2.Zip
```

(There is also a `..._zeno_0.98ng_1.wiff2.Zip`, ~96 MB, for Zeno-mode SWATH.)

Local dev layout this repo was validated against (not portable; adjust to your
machine): the bundle is unzipped under `~/thermo-raw-spike/data/sciex/` and run
as `cargo run --release -- <unzipped>/...K562_0.98ng_1.wiff`, reproducing 60
SWATH windows (399.5–899.9 m/z) and TOF calibration (4.898e-4 / -12.90).

## Usage

```
cargo run --release -- path/to/sample.wiff
```

## License

MIT OR Apache-2.0.
