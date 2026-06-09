//! Minimal pure-Rust reader for the SCIEX `.wiff` acquisition **method**
//! (SWATH isolation windows + TOF calibration). See `README.md` for the scope
//! and legitimacy notes — this reads the open OLE2 method, not the proprietary
//! `.wiff.scan` spectra.

use std::io;
use std::path::Path;

/// A SWATH isolation window (m/z bounds).
#[derive(Clone, Copy, Debug)]
pub struct SwathWindow {
    pub lower_mz: f64,
    pub upper_mz: f64,
}

impl SwathWindow {
    pub fn center_mz(&self) -> f64 {
        (self.lower_mz + self.upper_mz) / 2.0
    }
    pub fn width_mz(&self) -> f64 {
        self.upper_mz - self.lower_mz
    }
}

/// TOF → m/z calibration, of the form `m/z = (coef1·tof + coef2)²`.
#[derive(Clone, Copy, Debug)]
pub struct TofCalibration {
    pub coef1: f64,
    pub coef2: f64,
}

/// The acquisition method decoded from a `.wiff`.
#[derive(Clone, Debug)]
pub struct WiffMethod {
    pub swath_windows: Vec<SwathWindow>,
    pub tof_calibration: Option<TofCalibration>,
}

fn f64le(b: &[u8], o: usize) -> Option<f64> {
    b.get(o..o + 8)
        .map(|s| f64::from_le_bytes(s.try_into().unwrap()))
}

/// Read a stream's bytes. `Ok(None)` means the stream is genuinely absent;
/// `Err` means it exists but could not be read (corruption / short read).
fn read_stream(
    comp: &mut cfb::CompoundFile<std::fs::File>,
    name: &str,
) -> io::Result<Option<Vec<u8>>> {
    use std::io::Read;
    let mut s = match comp.open_stream(name) {
        Ok(s) => s,
        Err(e) if e.kind() == io::ErrorKind::NotFound => return Ok(None),
        Err(e) => return Err(e),
    };
    let mut v = Vec::new();
    s.read_to_end(&mut v)?;
    Ok(Some(v))
}

/// Open a `.wiff` (OLE2 compound file) and decode its SWATH method + TOF
/// calibration. The SWATH windows live in
/// `/MethodSubtree/Method1/DeviceMethod0/SWATHMethod` as 20-byte records
/// `{ f64 lower, f64 upper, u32 }` starting at offset 40; the calibration in
/// `/SampleSubtree/Sample1/TOFCalibrationData` (24-byte records from offset 32).
pub fn read_method<P: AsRef<Path>>(path: P) -> io::Result<WiffMethod> {
    let mut comp = cfb::open(path.as_ref())?;

    // SWATHMethod is required and must be well-formed: a 40-byte preamble + an
    // exact number of 20-byte `{ f64 lower, f64 upper, u32 }` records.
    let sw = read_stream(&mut comp, "/MethodSubtree/Method1/DeviceMethod0/SWATHMethod")?
        .ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "no SWATHMethod stream (not a SWATH .wiff?)")
        })?;
    const BASE: usize = 40;
    const STRIDE: usize = 20;
    if sw.len() < BASE || (sw.len() - BASE) % STRIDE != 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("SWATHMethod stream has a truncated record table ({} bytes)", sw.len()),
        ));
    }
    let n = (sw.len() - BASE) / STRIDE;
    let mut swath_windows = Vec::with_capacity(n);
    for k in 0..n {
        let o = BASE + k * STRIDE;
        let (lo, hi) = (f64le(&sw, o), f64le(&sw, o + 8));
        match (lo, hi) {
            (Some(lo), Some(hi))
                if lo.is_finite() && hi.is_finite() && hi > lo && lo.abs() < 1e7 && hi.abs() < 1e7 =>
            {
                swath_windows.push(SwathWindow { lower_mz: lo, upper_mz: hi });
            }
            _ => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("SWATHMethod record {k} is malformed (lower={lo:?}, upper={hi:?})"),
                ))
            }
        }
    }

    // TOF calibration is optional; absent or unreadably-short -> None.
    let tof_calibration = read_stream(&mut comp, "/SampleSubtree/Sample1/TOFCalibrationData")?
        .and_then(|cal| {
            let c1 = f64le(&cal, 32)?;
            let c2 = f64le(&cal, 40)?;
            (c1.is_finite() && c2.is_finite()).then_some(TofCalibration { coef1: c1, coef2: c2 })
        });

    Ok(WiffMethod {
        swath_windows,
        tof_calibration,
    })
}
