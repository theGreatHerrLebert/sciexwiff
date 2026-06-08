use std::io::Read;

fn read_stream(comp: &mut cfb::CompoundFile<std::fs::File>, name: &str) -> Option<Vec<u8>> {
    let mut s = comp.open_stream(name).ok()?;
    let mut v = Vec::new();
    s.read_to_end(&mut v).ok()?;
    Some(v)
}
fn f64le(b: &[u8], o: usize) -> f64 { f64::from_le_bytes(b[o..o + 8].try_into().unwrap()) }

fn main() {
    let path = std::env::args().nth(1).expect("usage: sciexwiff <file.wiff>");
    let mut comp = cfb::open(&path).expect("open .wiff (OLE2 compound file)");

    // --- SWATH isolation-window scheme: 60 records of {f64 lower, f64 upper, u32} from off 40 ---
    if let Some(sw) = read_stream(&mut comp, "/MethodSubtree/Method1/DeviceMethod0/SWATHMethod") {
        let (base, stride) = (40usize, 20usize);
        let n = sw.len().saturating_sub(base) / stride;
        println!("=== SWATH: {} variable isolation windows [lower, upper] ===", n);
        for k in 0..n {
            let o = base + k * stride;
            if o + 16 > sw.len() { break; }
            let lo = f64le(&sw, o);
            let hi = f64le(&sw, o + 8);
            if k < 4 || k + 2 >= n {
                println!("  win[{k:2}] = [{:.3}, {:.3}]  width {:.2}  center {:.3}", lo, hi, hi - lo, (lo + hi) / 2.0);
            } else if k == 4 {
                println!("   ...");
            }
        }
        if n > 0 {
            let lo0 = f64le(&sw, base);
            let hin = f64le(&sw, base + (n - 1) * stride + 8);
            println!("  coverage {:.1} .. {:.1} m/z", lo0, hin);
        }
    }

    // --- TOF -> m/z calibration coefficients (per-spectrum, 24-byte records) ---
    if let Some(cal) = read_stream(&mut comp, "/SampleSubtree/Sample1/TOFCalibrationData") {
        let nrec = (cal.len().saturating_sub(32)) / 24;
        println!("\n=== TOFCalibration: {} records (24B each) ===", nrec);
        for k in 0..3 {
            let o = 32 + k * 24;
            println!("  rec[{k}]: coef1={:.6e}  coef2={:.5}", f64le(&cal, o), f64le(&cal, o + 8));
        }
        println!("  (m/z = (coef1*tof + coef2)^2 form)");
    }
}
