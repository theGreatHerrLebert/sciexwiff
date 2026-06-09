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

fn read_stream(comp: &mut cfb::CompoundFile<std::fs::File>, name: &str) -> Option<Vec<u8>> {
    use std::io::Read;
    let mut s = comp.open_stream(name).ok()?;
    let mut v = Vec::new();
    s.read_to_end(&mut v).ok()?;
    Some(v)
}

/// Open a `.wiff` (OLE2 compound file) and decode its SWATH method + TOF
/// calibration. The SWATH windows live in
/// `/MethodSubtree/Method1/DeviceMethod0/SWATHMethod` as 20-byte records
/// `{ f64 lower, f64 upper, u32 }` starting at offset 40; the calibration in
/// `/SampleSubtree/Sample1/TOFCalibrationData` (24-byte records from offset 32).
pub fn read_method<P: AsRef<Path>>(path: P) -> io::Result<WiffMethod> {
    let mut comp = cfb::open(path.as_ref())?;

    let mut swath_windows = Vec::new();
    if let Some(sw) = read_stream(&mut comp, "/MethodSubtree/Method1/DeviceMethod0/SWATHMethod") {
        let (base, stride) = (40usize, 20usize);
        if sw.len() > base {
            let n = (sw.len() - base) / stride;
            for k in 0..n {
                let o = base + k * stride;
                if let (Some(lo), Some(hi)) = (f64le(&sw, o), f64le(&sw, o + 8)) {
                    if lo.is_finite() && hi.is_finite() && hi > lo {
                        swath_windows.push(SwathWindow {
                            lower_mz: lo,
                            upper_mz: hi,
                        });
                    }
                }
            }
        }
    }

    let tof_calibration = read_stream(&mut comp, "/SampleSubtree/Sample1/TOFCalibrationData")
        .and_then(|cal| {
            let c1 = f64le(&cal, 32)?;
            let c2 = f64le(&cal, 40)?;
            (c1.is_finite() && c2.is_finite()).then_some(TofCalibration {
                coef1: c1,
                coef2: c2,
            })
        });

    Ok(WiffMethod {
        swath_windows,
        tof_calibration,
    })
}
