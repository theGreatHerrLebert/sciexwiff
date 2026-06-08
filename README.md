# sciexwiff

Minimal pure-Rust reader for the SCIEX `.wiff` acquisition **method**, built to
extract the acquisition *design* a simulator needs to replicate a real SCIEX
run — **not** a spectra reader.

It opens the `.wiff` OLE2 / Compound File (via the [`cfb`](https://crates.io/crates/cfb)
crate) and decodes:

- **SWATH isolation scheme** — the variable isolation windows
  (`/MethodSubtree/Method1/DeviceMethod0/SWATHMethod`).
- **TOF → m/z calibration** coefficients
  (`/SampleSubtree/Sample1/TOFCalibrationData`).

## Why only the method (and not the spectra)

Established by inspection of real ZenoTOF 7600 data:

- `.wiff2` is **encrypted / high-entropy** — not open SQLite (the common web
  claim is wrong for current ZenoTOF data).
- `.wiff.scan` holds the spectra in a **proprietary packed binary** (no known
  codec; no Linux oracle — Clearcore2/WiffReader is Windows-bound and fails on
  Mono/Wine). SCIEX's EULA forbids decompiling Clearcore2.
- `.wiff` **is** open OLE2 and exposes the method (Period → Experiment cycle,
  61-window SWATH, Zeno attributes, TOF calibration).

**Conclusion: SCIEX *output* from the simulator goes via mzML** (or the pwiz
container, used as intended, for real→mzML ground truth). This crate supplies
the windows + calibration that let the simulator reproduce a real run's
acquisition layout. Reading real spectra is out of scope (and not needed for
simulated output).

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
