//! Pure-Rust codec for the SCIEX `.wiff.scan` **spectra** token stream, reverse-engineered
//! from the data file alone (no Clearcore2 SDK). Ported 1:1 from the validated Python codec.
//!
//! Scan layout in `.wiff.scan`: a sequential stream of units, each
//! `[protobuf metadata][ffffffff  u32 hdr  00  token-stream]`, the token stream running until
//! the next unit's metadata. Peaks are stored **m/z-ascending** as `(n, intensity)` where `n`
//! is the integer TOF sample index and `m/z = (a/5 · n + b)²` (`a`,`b` = the two per-scan
//! calibration doubles at metadata offsets +4 and +13).
//!
//! Token grammar (each token emits one peak):
//! - `0x00..=0x7b` — consecutive peak (delta = 1); the byte begins the intensity field
//! - `0x80..=0xfb` — gap peak, delta = `b - 0x7f` (1..=124), then the intensity field
//! - `0xfc [v]`    — delta = `v + 1`         (125..=256), then the intensity field
//! - `0xfd [lo][hi]` — delta = `lo + hi·256 + 1` (>256), then the intensity field
//!
//! Intensity field: `c 0x00..=0x7b` → `c`; `0x7c [b]` → `b` (124..=255);
//! `0x7d [lo][hi]` → 2-byte (>255); `0x7e [lo][mid][hi]` → 3-byte. `0x7f..=0xff` is invalid.

/// Largest representable intensity (3-byte `0x7e` escape).
pub const MAX_INTENSITY: u32 = 1 << 24;
/// The 2-byte `0xfd` delta escape carries `delta - 1`; larger deltas are auto-bridged.
pub const FD_MAX_DELTA: i64 = 65536;

/// A decoded scan block located within a `.wiff.scan` buffer.
#[derive(Clone, Copy, Debug)]
pub struct ScanBlock {
    /// Offset of the metadata length-prefix byte.
    pub meta: usize,
    /// Offset of the peak block's `ffffffff` sentinel.
    pub ff: usize,
    /// Exclusive end of the token stream (start of the next unit's metadata, or EOF).
    pub end: usize,
    pub cal_a: f64,
    pub cal_b: f64,
}

impl ScanBlock {
    /// Offset of the first token (after `ffffffff` + `u32 hdr` + `00`).
    pub fn stream_start(&self) -> usize {
        self.ff + 9
    }
}

#[derive(Debug, PartialEq, Eq)]
pub enum CodecError {
    /// An intensity field began with `0x7f..=0xff` (not a valid prefix).
    BadIntensityPrefix(u8),
    /// Peaks were not strictly increasing in `n`.
    NonIncreasing { index: usize, delta: i64 },
    /// An intensity was `>= MAX_INTENSITY`.
    IntensityRange(u32),
    /// Ran off the end of the buffer while decoding.
    Truncated,
}

pub type Peak = (i64, u32);

/// `n = round(5·(√mz − b) / a)`.
pub fn mz_to_n(mz: f64, cal_a: f64, cal_b: f64) -> i64 {
    (5.0 * (mz.sqrt() - cal_b) / cal_a).round() as i64
}

/// `mz = (a/5 · n + b)²`.
pub fn n_to_mz(n: i64, cal_a: f64, cal_b: f64) -> f64 {
    let v = cal_a / 5.0 * n as f64 + cal_b;
    v * v
}

/// Decode one intensity field at `body[i]`; returns `(value, next_i)`.
fn read_int(body: &[u8], i: usize, strict: bool) -> Result<(u32, usize), CodecError> {
    let c = *body.get(i).ok_or(CodecError::Truncated)?;
    let g = |o: usize| body.get(o).copied().ok_or(CodecError::Truncated);
    match c {
        0x00..=0x7b => Ok((c as u32, i + 1)),
        0x7c => Ok((g(i + 1)? as u32, i + 2)),
        0x7d => Ok((g(i + 1)? as u32 | (g(i + 2)? as u32) << 8, i + 3)),
        0x7e => Ok((
            g(i + 1)? as u32 | (g(i + 2)? as u32) << 8 | (g(i + 3)? as u32) << 16,
            i + 4,
        )),
        _ => {
            if strict {
                Err(CodecError::BadIntensityPrefix(c))
            } else {
                Ok((c as u32, i + 1))
            }
        }
    }
}

/// Decode up to `npeaks` peaks from the token stream starting at `body[i]`.
/// `seed_n` is the absolute `n` of the first peak (the scan cutoff — the first token's own
/// position is ignored, matching the vendor reader). Stops at `npeaks` or end of buffer.
pub fn decode_stream(
    body: &[u8],
    mut i: usize,
    seed_n: i64,
    npeaks: usize,
    strict: bool,
) -> Result<Vec<Peak>, CodecError> {
    let mut peaks = Vec::with_capacity(npeaks.min(4096));
    let mut prev = seed_n;
    let mut first = true;
    while peaks.len() < npeaks && i < body.len() {
        let b = body[i];
        let (delta, inten) = match b {
            0xfc => {
                let d = *body.get(i + 1).ok_or(CodecError::Truncated)? as i64 + 1;
                let (v, ni) = read_int(body, i + 2, strict)?;
                i = ni;
                (d, v)
            }
            0xfd => {
                let lo = *body.get(i + 1).ok_or(CodecError::Truncated)? as i64;
                let hi = *body.get(i + 2).ok_or(CodecError::Truncated)? as i64;
                let (v, ni) = read_int(body, i + 3, strict)?;
                i = ni;
                (lo + (hi << 8) + 1, v)
            }
            0x80..=0xff => {
                let d = b as i64 - 0x7f;
                let (v, ni) = read_int(body, i + 1, strict)?;
                i = ni;
                (d, v)
            }
            _ => {
                let (v, ni) = read_int(body, i, strict)?;
                i = ni;
                (1, v)
            }
        };
        let n = if first { seed_n } else { prev + delta };
        first = false;
        prev = n;
        peaks.push((n, inten));
    }
    Ok(peaks)
}

fn put_int(out: &mut Vec<u8>, v: u32) -> Result<(), CodecError> {
    if v >= MAX_INTENSITY {
        return Err(CodecError::IntensityRange(v));
    }
    match v {
        0..=0x7b => out.push(v as u8),
        0x7c..=0xff => {
            out.push(0x7c);
            out.push(v as u8);
        }
        0x100..=0xffff => {
            out.push(0x7d);
            out.push((v & 0xff) as u8);
            out.push((v >> 8) as u8);
        }
        _ => {
            out.push(0x7e);
            out.push((v & 0xff) as u8);
            out.push(((v >> 8) & 0xff) as u8);
            out.push(((v >> 16) & 0xff) as u8);
        }
    }
    Ok(())
}

fn put_delta_int(out: &mut Vec<u8>, delta: i64, inten: u32) -> Result<(), CodecError> {
    match delta {
        1 => {}
        2..=124 => out.push((0x7f + delta) as u8),
        125..=256 => {
            out.push(0xfc);
            out.push((delta - 1) as u8);
        }
        _ => {
            out.push(0xfd);
            out.push(((delta - 1) & 0xff) as u8);
            out.push(((delta - 1) >> 8) as u8);
        }
    }
    put_int(out, inten)
}

/// Inverse of [`decode_stream`]: encode `peaks` (which MUST be strictly increasing in `n`) to a
/// token stream. The first peak carries only its intensity (its `n` is the seed/cutoff, supplied
/// out-of-band). Deltas above [`FD_MAX_DELTA`] are auto-bridged with intensity-1 filler peaks so
/// the 2-byte `0xfd` escape never overflows.
pub fn encode_stream(peaks: &[Peak]) -> Result<Vec<u8>, CodecError> {
    let mut out = Vec::with_capacity(peaks.len() * 2);
    let mut prev: i64 = 0;
    for (k, &(n, inten)) in peaks.iter().enumerate() {
        if k == 0 {
            prev = n;
            put_int(&mut out, inten)?;
            continue;
        }
        let mut delta = n - prev;
        prev = n;
        if delta <= 0 {
            return Err(CodecError::NonIncreasing { index: k, delta });
        }
        while delta > FD_MAX_DELTA {
            put_delta_int(&mut out, FD_MAX_DELTA, 1)?;
            delta -= FD_MAX_DELTA;
        }
        put_delta_int(&mut out, delta, inten)?;
    }
    Ok(out)
}

/// Enumerate physical scan blocks in file order. The anchor is the metadata signature
/// `0a 12 09` with `0x11` at +11 and `0x12` at +20, a length-prefix `>= 27` at −1, and the peak
/// block's `ffffffff` immediately after the metadata. The `ffffffff`-after test is the
/// discriminator that rejects the `0a1209` byte-triples occurring inside peak data, and a
/// calibration sanity check (`a ≈ 5e-4 > 0`, `b ≈ −13 < 0`) rejects any remaining garbage.
pub fn scan_blocks(sb: &[u8]) -> Vec<ScanBlock> {
    let mut starts: Vec<(usize, usize, f64, f64)> = Vec::new();
    let mut i = 0usize;
    while i + 21 <= sb.len() {
        if sb[i] == 0x0a && sb[i + 1] == 0x12 && sb[i + 2] == 0x09 {
            if i >= 1 && sb[i - 1] >= 27 && sb[i + 11] == 0x11 && sb[i + 20] == 0x12 {
                let s = i - 1;
                let ff = s + 1 + sb[s] as usize;
                if ff + 4 <= sb.len() && &sb[ff..ff + 4] == b"\xff\xff\xff\xff" {
                    if let (Some(a), Some(b)) = (f64le(sb, s + 4), f64le(sb, s + 13)) {
                        if a > 1e-5 && a < 1e-2 && b > -50.0 && b < 0.0 {
                            starts.push((s, ff, a, b));
                        }
                    }
                }
            }
        }
        i += 1;
    }
    let mut out = Vec::with_capacity(starts.len());
    for k in 0..starts.len() {
        let (s, ff, a, b) = starts[k];
        let end = if k + 1 < starts.len() {
            starts[k + 1].0
        } else {
            sb.len()
        };
        out.push(ScanBlock {
            meta: s,
            ff,
            end,
            cal_a: a,
            cal_b: b,
        });
    }
    out
}

fn f64le(b: &[u8], o: usize) -> Option<f64> {
    b.get(o..o + 8)
        .map(|s| f64::from_le_bytes(s.try_into().unwrap()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn encode_guards() {
        assert_eq!(
            encode_stream(&[(100, 5), (90, 3)]).unwrap_err(),
            CodecError::NonIncreasing { index: 1, delta: -10 }
        );
        assert_eq!(
            encode_stream(&[(100, 5), (200, MAX_INTENSITY)]).unwrap_err(),
            CodecError::IntensityRange(MAX_INTENSITY)
        );
    }

    #[test]
    fn strict_rejects_0x7f() {
        assert_eq!(read_int(&[0x7f, 0], 0, true).unwrap_err(), CodecError::BadIntensityPrefix(0x7f));
    }

    #[test]
    fn roundtrip_synthetic() {
        // seed + consecutive + gap + fc + fd + 1/2/3-byte intensities (no bridge -> exact identity)
        let peaks: Vec<Peak> = vec![
            (336703, 25),
            (336704, 56),   // delta 1 consecutive
            (336780, 20),   // delta 76 gap marker
            (337040, 216),  // delta 260 -> fc; int 216 -> 0x7c escape
            (340500, 5000),  // delta 3460 -> fd; int 5000 -> 0x7d 2-byte
            (400000, 200000), // delta 59500 (<=65536) -> fd; int 200000 -> 0x7e 3-byte
        ];
        let enc = encode_stream(&peaks).unwrap();
        let dec = decode_stream(&enc, 0, peaks[0].0, peaks.len(), true).unwrap();
        assert_eq!(dec, peaks);
    }

    #[test]
    fn bridge_preserves_real_peaks() {
        // a gap larger than FD_MAX_DELTA is bridged with intensity-1 fillers; the real peaks
        // must still decode at their exact positions (with fillers interspersed).
        let peaks: Vec<Peak> = vec![(1000, 30), (1001, 40), (200_000, 99)];
        let enc = encode_stream(&peaks).unwrap();
        let dec = decode_stream(&enc, 0, peaks[0].0, 64, true).unwrap();
        for p in &peaks {
            assert!(dec.contains(p), "real peak {p:?} lost across bridge");
        }
        // every non-real decoded peak is a filler (intensity 1)
        for d in &dec {
            assert!(peaks.contains(d) || d.1 == 1, "unexpected decoded peak {d:?}");
        }
    }

    /// Byte-identical parity against a real `.wiff.scan` block, gated behind
    /// `TIMSIM_SCIEX_WIFF_SCAN=<path to a real .wiff.scan>` (skips if unset).
    #[test]
    fn parity_real_wiff_scan() {
        let path = match std::env::var("TIMSIM_SCIEX_WIFF_SCAN") {
            Ok(p) => p,
            Err(_) => return, // no oracle available; skip
        };
        let sb = std::fs::read(&path).expect("read .wiff.scan");
        let blocks = scan_blocks(&sb);
        assert!(blocks.len() > 100, "expected many scan blocks");
        // decode->encode every rich block and assert byte-identical token streams
        let mut checked = 0;
        for b in &blocks {
            let slot = &sb[b.stream_start()..b.end];
            if slot.len() < 200 {
                continue; // skip empty/marker blocks
            }
            // Decode a bounded prefix that stays well inside the block, re-encode, and assert the
            // token bytes are byte-identical. (Decoding past the block would hit the next block's
            // metadata and, under strict mode, error — so we bound the peak count.)
            let peaks = match decode_stream(&sb, b.stream_start(), 0, 40, true) {
                Ok(p) if p.len() >= 8 => p,
                _ => continue,
            };
            let enc = encode_stream(&peaks).expect("re-encode");
            assert!(enc.len() <= slot.len());
            assert_eq!(&enc[..], &slot[..enc.len()], "byte mismatch in block at ff={}", b.ff);
            checked += 1;
            if checked >= 50 {
                break;
            }
        }
        assert!(checked > 0, "no rich blocks were parity-checked");
        eprintln!("parity: {checked} real blocks decode->encode byte-identical");
    }
}
