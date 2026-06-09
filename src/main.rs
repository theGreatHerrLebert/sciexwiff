//! `sciexwiff <file.wiff>` — dump the SWATH window scheme + TOF calibration.

fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("usage: sciexwiff <file.wiff>");
    let m = sciexwiff::read_method(&path).expect("read .wiff method");

    let n = m.swath_windows.len();
    println!("=== SWATH: {n} variable isolation windows [lower, upper] ===");
    for (k, w) in m.swath_windows.iter().enumerate() {
        if k < 4 || k + 2 >= n {
            println!(
                "  win[{k:2}] = [{:.3}, {:.3}]  width {:.2}  center {:.3}",
                w.lower_mz,
                w.upper_mz,
                w.width_mz(),
                w.center_mz()
            );
        } else if k == 4 {
            println!("   ...");
        }
    }
    if let (Some(first), Some(last)) = (m.swath_windows.first(), m.swath_windows.last()) {
        println!("  coverage {:.1} .. {:.1} m/z", first.lower_mz, last.upper_mz);
    }
    if let Some(c) = m.tof_calibration {
        println!(
            "\n=== TOFCalibration: coef1={:.6e}  coef2={:.5}  (m/z = (coef1*tof + coef2)^2) ===",
            c.coef1, c.coef2
        );
    }
}
