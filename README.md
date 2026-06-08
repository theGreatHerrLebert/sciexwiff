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

> **PRIDE `PXD036786`** — ZenoTOF 7600, K562 0.98 ng, ~53 MB `.wiff2.Zip`
> (`.wiff` OLE2 ~11 MB + `.wiff2` ~29 MB + `.wiff.scan` ~32 MB).

No vendor data is committed here; re-download from PRIDE to run.

## Usage

```
cargo run --release -- path/to/sample.wiff
```

## License

MIT OR Apache-2.0.
