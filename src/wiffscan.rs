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
//! `0x7d [lo][hi]` → 2-byte (>255); `0x7e [lo][mid][hi]` → 3-byte; `0x80..=0xff` → raw 128..=255.
//!
//! The vendor encoder is NON-CANONICAL: some blocks encode delta 127/128 as `0xfe`/`0xff` and some
//! store intensity 128..=255 raw vs `0x7c`-escaped. The DECODER accepts all forms; the ENCODER
//! emits one canonical form (which the reader accepts). The writer therefore copies unedited blocks
//! verbatim and re-encodes only the blocks it authors — so `decode → encode` is peaks-stable but
//! not universally byte-identical (see the `rebuild` writer and the `parity_real_wiff_scan` test).

/// Exclusive upper bound on intensity (3-byte `0x7e` escape holds up to `MAX_INTENSITY - 1`).
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
        0x00..=0x7b => Ok((c as u32, i + 1)),                          // raw 0..123
        0x7c => Ok((g(i + 1)? as u32, i + 2)),                         // 1-byte escape 124..255
        0x7d => Ok((g(i + 1)? as u32 | (g(i + 2)? as u32) << 8, i + 3)),               // 2-byte
        0x7e => Ok((
            g(i + 1)? as u32 | (g(i + 2)? as u32) << 8 | (g(i + 3)? as u32) << 16,     // 3-byte
            i + 4,
        )),
        0x7f => {
            if strict { Err(CodecError::BadIntensityPrefix(c)) } else { Ok((c as u32, i + 1)) }
        }
        _ => Ok((c as u32, i + 1)), // 0x80..=0xff : raw intensity 128..255 (block-dependent)
    }
}

/// Decode up to `npeaks` peaks from the token stream starting at `body[i]`.
/// `seed_n` is the absolute `n` of the first peak (the scan cutoff — the first token's own
/// position is ignored, matching the vendor reader). Stops at `npeaks` or end of buffer.
pub fn decode_stream(
    body: &[u8],
    i: usize,
    seed_n: i64,
    npeaks: usize,
    strict: bool,
) -> Result<Vec<Peak>, CodecError> {
    Ok(decode_tracked(body, i, seed_n, npeaks, strict).0)
}

/// Like `decode_stream` but also returns, for each peak, its token byte-span [start,end) and the
/// final cursor. A trailing PARTIAL token stops decoding cleanly (matches the reader tolerating a
/// truncated tail). This is the primitive the writer uses to split payload from terminator.
pub fn decode_tracked(
    body: &[u8],
    mut i: usize,
    seed_n: i64,
    npeaks: usize,
    strict: bool,
) -> (Vec<Peak>, Vec<(usize, usize)>, usize) {
    let mut peaks = Vec::new();
    let mut spans = Vec::new();
    let mut prev = seed_n;
    let mut first = true;
    let n = body.len();
    while peaks.len() < npeaks && i < n {
        let start = i;
        let b = body[i];
        let res: Option<(i64, u32, usize)> = match b {
            0xfc => body.get(i + 1).and_then(|&d1| {
                read_int(body, i + 2, strict).ok().map(|(v, ni)| (d1 as i64 + 1, v, ni))
            }),
            0xfd => match (body.get(i + 1), body.get(i + 2)) {
                (Some(&lo), Some(&hi)) => read_int(body, i + 3, strict)
                    .ok()
                    .map(|(v, ni)| (lo as i64 + ((hi as i64) << 8) + 1, v, ni)),
                _ => None,
            },
            0x80..=0xff => read_int(body, i + 1, strict).ok().map(|(v, ni)| (b as i64 - 0x7f, v, ni)),
            _ => read_int(body, i, strict).ok().map(|(v, ni)| (1i64, v, ni)),
        };
        let (delta, inten, ni) = match res {
            Some(t) => t,
            None => break, // partial trailing token
        };
        i = ni;
        let nval = if first { seed_n } else { prev + delta };
        first = false;
        prev = nval;
        peaks.push((nval, inten));
        spans.push((start, i));
    }
    (peaks, spans, i)
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
    // Encode canonically as gap markers (0x80..=0xfb -> 1..=124) then the 0xfc/0xfd escapes.
    // NOTE: the vendor is non-canonical — some blocks encode delta 127/128 as 0xfe/0xff and some
    // intensities 128..255 raw; the DECODER accepts all of those. We emit one canonical form (the
    // reader accepts it); unedited blocks are copied verbatim so their exact bytes are preserved.
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
        let mut delta = n.checked_sub(prev).ok_or(CodecError::NonIncreasing { index: k, delta: 0 })?;
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
        if sb[i] == 0x0a && sb[i + 1] == 0x12 && sb[i + 2] == 0x09
            && i >= 1 && sb[i - 1] >= 27 && sb[i + 11] == 0x11 && sb[i + 20] == 0x12 {
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

// ============================================================================================
// Writer: rebuild a .wiff.scan (author arbitrary peaks/scan) + recompute the Idx directory.
//
// The Idx stream (108 B/record) indexes a 2-block SEGMENT: +32 = block1 ffffffff offset (rel byte
// 44), +86 = block2 ffffffff offset, +36 = block2_meta_start - block1_ff. The writer streams every
// physical block into a new buffer (edited payloads replaced, the rest verbatim) and recomputes
// +32/+86/+36 from the NEW positions — so there is no in-place offset patching and cumulative
// edits are handled correctly. PyO3-friendly: all inputs/outputs are byte slices / Vec<u8>.
// ============================================================================================

const IDX_RS: usize = 108;

/// Errors from the writer path. Kept separate from [`CodecError`] so the FFI boundary can map a
/// structural/bounds failure to a Python exception instead of panicking across it.
#[derive(Debug, PartialEq, Eq)]
pub enum WriterError {
    /// A buffer (scan or Idx) was shorter than the structure required.
    Truncated,
    /// A physical block's offsets were inconsistent (meta..ff..end out of order or out of range).
    BadBlock(usize),
    /// A segment references block indices out of range, or is structurally impossible (e.g. a
    /// block2 with no block1 in a 2-block segment), or produced an offset that underflows.
    BadSegment(usize),
}

fn u32le(b: &[u8], o: usize) -> Result<u32, WriterError> {
    let s = b.get(o..o + 4).ok_or(WriterError::Truncated)?;
    Ok(u32::from_le_bytes(s.try_into().expect("length checked")))
}
fn put_u32le(b: &mut [u8], o: usize, v: u32) -> Result<(), WriterError> {
    let s = b.get_mut(o..o + 4).ok_or(WriterError::Truncated)?;
    s.copy_from_slice(&v.to_le_bytes());
    Ok(())
}
fn find_sub(hay: &[u8], needle: &[u8]) -> Option<usize> {
    if needle.is_empty() || needle.len() > hay.len() { return None; }
    hay.windows(needle.len()).position(|w| w == needle)
}

/// One Idx record mapped onto physical block indices (block1 always present when indexed).
#[derive(Clone, Copy, Debug)]
pub struct Segment {
    pub rec: usize,
    pub b1: Option<usize>,
    pub b2: Option<usize>,
}

/// Validate that every block's `meta <= ff`, `ff + 9 <= end`, `end <= scan.len()`.
fn validate_blocks(scan: &[u8], blocks: &[ScanBlock]) -> Result<(), WriterError> {
    if scan.len() < 44 { return Err(WriterError::Truncated); }
    for (i, b) in blocks.iter().enumerate() {
        if !(b.meta <= b.ff && b.ff + 9 <= b.end && b.end <= scan.len()) {
            return Err(WriterError::BadBlock(i));
        }
    }
    Ok(())
}

/// Enumerate physical blocks and map every Idx record onto them by offset.
pub fn map_segments(scan: &[u8], idx: &[u8]) -> Result<(Vec<ScanBlock>, Vec<Segment>), WriterError> {
    let blocks = scan_blocks(scan);
    validate_blocks(scan, &blocks)?;
    let mut ff_to_i = std::collections::HashMap::with_capacity(blocks.len());
    for (i, b) in blocks.iter().enumerate() {
        ff_to_i.insert(b.ff, i);
    }
    let mut segs = Vec::new();
    for r in 0..idx.len() / IDX_RS {
        let o32 = u32le(idx, r * IDX_RS + 32)? as usize;
        let o86 = u32le(idx, r * IDX_RS + 86)? as usize;
        let b1 = if o32 > 0 { ff_to_i.get(&(44 + o32)).copied() } else { None };
        let b2 = if o86 > 0 && o86 != o32 { ff_to_i.get(&(44 + o86)).copied() } else { None };
        if b1.is_some() || b2.is_some() {
            segs.push(Segment { rec: r, b1, b2 });
        }
    }
    Ok((blocks, segs))
}

/// Decode an editable block's payload peaks and its terminator byte-length. Trailing tokens whose
/// bytes are ALL 0xff are the block terminator (2 or 4 bytes) and are dropped. Caller must have
/// validated `b` (see [`validate_blocks`]).
pub fn block_payload(scan: &[u8], b: &ScanBlock) -> Result<(Vec<Peak>, usize), WriterError> {
    let body = scan.get(b.ff + 9..b.end).ok_or(WriterError::BadBlock(0))?;
    let (mut peaks, mut spans, _) = decode_tracked(body, 0, 0, usize::MAX, false);
    let mut term = 0usize;
    while let Some(&(s, e)) = spans.last() {
        if body[s..e].iter().all(|&x| x == 0xff) {
            term += e - s;
            peaks.pop();
            spans.pop();
        } else {
            break;
        }
    }
    Ok((peaks, term))
}

/// Output of [`rebuild`]: (new `.wiff.scan` bytes, new `ffffffff` position per block, metadata
/// length per block). The latter two feed [`recompute_idx`].
pub type RebuildOutput = (Vec<u8>, Vec<usize>, Vec<usize>);

/// Rebuild the `.wiff.scan`. `edits` maps a physical block index -> new payload token bytes (as
/// from `encode_stream`). Returns (new_scan, new_ff positions, metadata lengths).
pub fn rebuild(
    scan: &[u8],
    blocks: &[ScanBlock],
    edits: &std::collections::HashMap<usize, Vec<u8>>,
) -> Result<RebuildOutput, WriterError> {
    validate_blocks(scan, blocks)?;
    let mut out = Vec::with_capacity(scan.len() + 4096);
    out.extend_from_slice(&scan[..44]);
    let mut new_ff = vec![0usize; blocks.len()];
    let mut metalen = vec![0usize; blocks.len()];
    for (i, b) in blocks.iter().enumerate() {
        metalen[i] = b.ff - b.meta;
        out.extend_from_slice(&scan[b.meta..b.ff]); // metadata
        new_ff[i] = out.len(); // this block's ffffffff
        match edits.get(&i) {
            Some(payload) => {
                out.extend_from_slice(&scan[b.ff..b.ff + 9]); // ffffffff + u32hdr + 00
                out.extend_from_slice(payload);
                let (_, mut term) = block_payload(scan, b)?; // preserve original terminator length
                if term == 0 { term = 4; }
                out.extend(std::iter::repeat_n(0xffu8, term));
            }
            None => out.extend_from_slice(&scan[b.ff..b.end]), // verbatim
        }
    }
    Ok((out, new_ff, metalen))
}

/// Recompute +32/+86/+36 for every segment from the new block positions. Returns a new Idx buffer
/// (same length as the input). A segment with a block2 but no block1 is rejected (a 2-block
/// segment cannot be described without its block1 anchor).
pub fn recompute_idx(
    idx: &[u8],
    segs: &[Segment],
    new_ff: &[usize],
    metalen: &[usize],
) -> Result<Vec<u8>, WriterError> {
    let mut idx = idx.to_vec();
    let err = |rec: usize| WriterError::BadSegment(rec);
    let off = |ff: usize, rec: usize| ff.checked_sub(44).ok_or(err(rec))
        .and_then(|o| u32::try_from(o).map_err(|_| err(rec)));
    for s in segs {
        // +36 (block2_meta relative to block1) can only be written when BOTH blocks are present.
        // A segment with only block2 (empty block1 slot, o32==0) keeps its original +36.
        if let Some(i1) = s.b1 {
            let ff1 = *new_ff.get(i1).ok_or(err(s.rec))?;
            put_u32le(&mut idx, s.rec * IDX_RS + 32, off(ff1, s.rec)?)?;
            if let Some(i2) = s.b2 {
                let ff2 = *new_ff.get(i2).ok_or(err(s.rec))?;
                let ml2 = *metalen.get(i2).ok_or(err(s.rec))?;
                let block2_meta = ff2.checked_sub(ml2).ok_or(err(s.rec))?;
                let o36 = block2_meta.checked_sub(ff1).ok_or(err(s.rec))?;
                put_u32le(&mut idx, s.rec * IDX_RS + 86, off(ff2, s.rec)?)?;
                put_u32le(&mut idx, s.rec * IDX_RS + 36, u32::try_from(o36).map_err(|_| err(s.rec))?)?;
            }
        } else if let Some(i2) = s.b2 {
            let ff2 = *new_ff.get(i2).ok_or(err(s.rec))?;
            put_u32le(&mut idx, s.rec * IDX_RS + 86, off(ff2, s.rec)?)?; // +36 unchanged (no block1)
        }
    }
    Ok(idx)
}

// ============================================================================================
// GROW writer: author arbitrary token lengths (peaks/scan) and retranslate EVERY Idx offset.
//
// The [`rebuild`] path only recomputes the *enumerated* segments, so it corrupts the ~4% of Idx
// records that reference blocks `scan_blocks` does not enumerate (empty / embedded-`ffffffff`
// blocks) once byte positions shift. The grow path avoids that by building an old→new byte
// offset MAP (from where each grown block adds/removes bytes) and translating every Idx offset
// field through it — so mapped and unmapped records alike stay valid. Only blocks with a clean
// tail (`tokens + terminator`, no embedded `ffffffff`) are grown; the caller leaves the rare
// messy blocks out of `edits` (copied verbatim).
// ============================================================================================

/// A grown block edit: replace the block's whole body `[ff+9 .. end]` with `tokens` (the
/// caller supplies a complete token stream — the reader bounds the scan by the recomputed Idx
/// span, so no terminator is required). `block` indexes into the `blocks` slice. The block must
/// have a clean tail (no embedded `ffffffff`), else replacing the body would drop a mini-block.
pub struct GrowEdit {
    pub block: usize,
    pub tokens: Vec<u8>,
}

/// Rebuild the `.wiff.scan`, replacing each edited block's body with an arbitrary-length token
/// stream. Returns the new bytes and a sorted `(old_offset, cumulative_delta)` breakpoint list
/// for [`translate_offset`] / [`retranslate_idx`]. Edits are applied in ascending block
/// position; everything between blocks (metadata, mini-blocks, trailing bytes) is copied
/// verbatim and simply shifts. Using the whole `[ff+9 .. end]` body as the edit span avoids any
/// token/terminator-boundary guess (which is fragile for blocks whose last peak encodes 0xff).
pub fn rebuild_grow(
    scan: &[u8],
    blocks: &[ScanBlock],
    edits: &[GrowEdit],
) -> Result<(Vec<u8>, Vec<(usize, i64)>), WriterError> {
    validate_blocks(scan, blocks)?;
    // Sort edits by the block body start; ensure they are disjoint and in range.
    let mut ordered: Vec<(usize, usize, &[u8], usize)> = Vec::with_capacity(edits.len());
    for e in edits {
        let b = blocks.get(e.block).ok_or(WriterError::BadBlock(e.block))?;
        let tok_start = b.ff + 9;
        let tok_end = b.end; // replace the whole body
        ordered.push((tok_start, tok_end, &e.tokens, e.block));
    }
    ordered.sort_by_key(|&(s, _, _, _)| s);
    let mut out = Vec::with_capacity(scan.len() + (1 << 16));
    let mut breakpoints: Vec<(usize, i64)> = Vec::with_capacity(ordered.len());
    let mut cum: i64 = 0;
    let mut cursor = 0usize;
    for (tok_start, tok_end, tokens, bi) in ordered {
        if tok_start < cursor || tok_end < tok_start || tok_end > scan.len() {
            return Err(WriterError::BadBlock(bi)); // overlapping / duplicate / malformed edit
        }
        out.extend_from_slice(&scan[cursor..tok_start]); // verbatim up to the token stream
        out.extend_from_slice(tokens); // the grown tokens
        cum += tokens.len() as i64 - (tok_end - tok_start) as i64;
        breakpoints.push((tok_start, cum)); // any old offset >= tok_end shifts by `cum`
        cursor = tok_end; // skip the old tokens; the terminator (from tok_end) copies next
    }
    out.extend_from_slice(&scan[cursor..]); // rest verbatim
    Ok((out, breakpoints))
}

/// Translate an old byte offset to its new position given the `rebuild_grow` breakpoints.
/// Idx offsets point at `ffffffff` block starts, which never fall inside a grown token region,
/// so the mapping is unambiguous. Returns `None` if the (possibly negative) cumulative delta
/// underflows the offset — a malformed edit set rather than a valid position.
pub fn translate_offset(old: usize, breakpoints: &[(usize, i64)]) -> Option<usize> {
    // Largest breakpoint whose token-start position is < `old` applies (a block that starts at
    // or after `old` does not shift it). Breakpoints are sorted by position.
    let i = breakpoints.partition_point(|&(pos, _)| pos < old);
    let d = if i == 0 { 0 } else { breakpoints[i - 1].1 };
    (old as i64).checked_add(d).and_then(|v| usize::try_from(v).ok())
}

/// Retranslate every Idx record's `+32` / `+86` / `+36` offsets through the `rebuild_grow`
/// breakpoint map. Handles ALL records (enumerated or not); the Idx length is unchanged.
pub fn retranslate_idx(idx: &[u8], breakpoints: &[(usize, i64)]) -> Result<Vec<u8>, WriterError> {
    // The Idx stream is `n` full 108-byte records + a short trailing tail (a real property of
    // the format); only the full records are retranslated, the tail is copied verbatim.
    let mut out = idx.to_vec();
    let tr = |abs: usize| translate_offset(abs, breakpoints).ok_or(WriterError::BadSegment(0));
    let to_rel = |abs: usize| -> Result<u32, WriterError> {
        abs.checked_sub(44)
            .ok_or(WriterError::Truncated)
            .and_then(|o| u32::try_from(o).map_err(|_| WriterError::Truncated))
    };
    for r in 0..idx.len() / IDX_RS {
        let base = r * IDX_RS;
        let o32 = u32le(idx, base + 32)? as usize;
        let o86 = u32le(idx, base + 86)? as usize;
        let o36 = u32le(idx, base + 36)? as usize;
        if o32 > 0 {
            let ff1_old = 44usize.checked_add(o32).ok_or(WriterError::Truncated)?;
            let ff1_new = tr(ff1_old)?;
            put_u32le(&mut out, base + 32, to_rel(ff1_new)?)?;
            // +36 is block2_meta relative to block1; recompute from the translated positions.
            let block2_meta_old = ff1_old.checked_add(o36).ok_or(WriterError::Truncated)?;
            let new_o36 = tr(block2_meta_old)?.checked_sub(ff1_new).ok_or(WriterError::BadSegment(r))?;
            put_u32le(&mut out, base + 36, u32::try_from(new_o36).map_err(|_| WriterError::Truncated)?)?;
        }
        if o86 > 0 {
            let ff2_new = tr(44usize.checked_add(o86).ok_or(WriterError::Truncated)?)?;
            put_u32le(&mut out, base + 86, to_rel(ff2_new)?)?;
        }
    }
    Ok(out)
}

/// Locate the `Idx` stream inside a `.wiff` CFB by an unambiguous 4-record needle (the Idx stream
/// is unchanged in length, so it can be patched in place). Returns the stream's byte offset.
pub fn locate_idx_in_wiff(raw: &[u8], idx: &[u8]) -> Option<usize> {
    let anchor = 5000 * IDX_RS;
    let needle = idx.get(anchor..anchor + IDX_RS * 4)?;
    let pos = find_sub(raw, needle)?;
    if find_sub(&raw[pos + 1..], needle).is_some() {
        return None; // not unique
    }
    let ib = pos.checked_sub(anchor)?;
    for k in [0usize, 1000, 20000, 51000] {
        if raw.get(ib + k * IDX_RS..ib + (k + 1) * IDX_RS)? != &idx[k * IDX_RS..(k + 1) * IDX_RS] {
            return None;
        }
    }
    Some(ib)
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

    /// Writer round-trip against a real file (gated on `TIMSIM_SCIEX_WIFF_SCAN`, needs the sibling
    /// `.wiff`): a no-edit rebuild must be byte-identical, and recomputing the Idx from the
    /// unchanged positions must reproduce the original Idx (proving the +32/+86/+36 formulas).
    #[test]
    fn rebuild_roundtrip_real() {
        let scan_path = match std::env::var("TIMSIM_SCIEX_WIFF_SCAN") {
            Ok(p) => p,
            Err(_) => return,
        };
        let wiff_path = scan_path.strip_suffix(".scan").expect("path ends with .scan").to_string();
        let scan = std::fs::read(&scan_path).expect("read .wiff.scan");
        let mut comp = cfb::open(&wiff_path).expect("open .wiff CFB");
        let mut idx = Vec::new();
        {
            use std::io::Read;
            comp.open_stream("SampleSubtree/Sample1/Idx")
                .expect("Idx stream")
                .read_to_end(&mut idx)
                .unwrap();
        }
        let (blocks, segs) = map_segments(&scan, &idx).expect("map_segments");
        assert!(blocks.len() > 1000 && segs.len() > 1000);
        let (new_scan, new_ff, metalen) = rebuild(&scan, &blocks, &std::collections::HashMap::new()).expect("rebuild");
        assert_eq!(new_scan, scan, "no-edit rebuild not byte-identical");
        let new_idx = recompute_idx(&idx, &segs, &new_ff, &metalen).expect("recompute_idx");
        assert_eq!(new_idx, idx, "no-edit Idx recompute != original");
        eprintln!("rebuild: {} blocks, {} indexed segments; no-edit round-trip OK", blocks.len(), segs.len());
    }

    /// Grow-writer identity: rebuilding with every block's tokens replaced by *its own* bytes
    /// must reproduce the file byte-for-byte, all deltas zero, and retranslating the Idx must
    /// leave it unchanged — proving the offset-map machinery is a correct no-op at rest. Gated
    /// on `TIMSIM_SCIEX_WIFF_SCAN`.
    #[test]
    fn rebuild_grow_identity_real() {
        let scan_path = match std::env::var("TIMSIM_SCIEX_WIFF_SCAN") {
            Ok(p) => p,
            Err(_) => return,
        };
        let wiff_path = scan_path.strip_suffix(".scan").expect("ends with .scan").to_string();
        let scan = std::fs::read(&scan_path).expect("read .wiff.scan");
        let mut comp = cfb::open(&wiff_path).expect("open .wiff CFB");
        let mut idx = Vec::new();
        {
            use std::io::Read;
            comp.open_stream("SampleSubtree/Sample1/Idx").unwrap().read_to_end(&mut idx).unwrap();
        }
        let blocks = scan_blocks(&scan);
        // Identity edits: each block's tokens = exactly its own whole body [ff+9 .. end].
        let mut edits = Vec::new();
        for (i, b) in blocks.iter().enumerate() {
            edits.push(GrowEdit { block: i, tokens: scan[b.ff + 9..b.end].to_vec() });
        }
        let (new_scan, bp) = rebuild_grow(&scan, &blocks, &edits).expect("rebuild_grow");
        assert_eq!(new_scan, scan, "identity grow not byte-identical");
        assert!(bp.iter().all(|&(_, c)| c == 0), "identity deltas must be zero");
        let new_idx = retranslate_idx(&idx, &bp).expect("retranslate_idx");
        assert_eq!(new_idx, idx, "identity retranslate != original Idx");
        eprintln!("rebuild_grow identity OK: {} blocks, Idx {} records", blocks.len(), idx.len() / IDX_RS);
    }

    #[test]
    fn translate_offset_shifts_after_grown_block() {
        // A grown block at old token-start 100 that adds 5 bytes shifts everything at/after it.
        let bp = vec![(100usize, 5i64)];
        assert_eq!(translate_offset(50, &bp), Some(50), "before the grow: unchanged");
        assert_eq!(translate_offset(100, &bp), Some(100), "at the grow start: unchanged");
        assert_eq!(translate_offset(200, &bp), Some(205), "after the grow: +5");
        // A net-negative delta that underflows an offset is rejected, not wrapped.
        assert_eq!(translate_offset(3, &[(0usize, -10i64)]), None);
    }

    // Build a synthetic .wiff.scan: 44-byte header + `n` blocks, each `[meta len-prefix +
    // 0a1209 + cal_a f64 + 11 + cal_b f64 + 12 06 ..pad..][ffffffff][hdr u32][00][tokens]`.
    // Returns (bytes, blocks). cal chosen in the sane range so `scan_blocks` enumerates them.
    fn synth_scan(token_lens: &[usize]) -> (Vec<u8>, Vec<ScanBlock>) {
        let mut b = vec![0u8; 44];
        for &tl in token_lens {
            // metadata: 0x1d len, then 0a1209 <a> 11 <b> 12 06 08 00 10 08  (28 bytes body)
            let mut meta = vec![0x0a, 0x12, 0x09];
            meta.extend_from_slice(&0.000489823_f64.to_le_bytes());
            meta.push(0x11);
            meta.extend_from_slice(&(-12.9765_f64).to_le_bytes());
            // pad the field-2 submessage so the metadata is >= 27 bytes (scan_blocks requires it)
            meta.extend_from_slice(&[0x12, 0x08, 0x08, 0x00, 0x10, 0x08, 0x00, 0x00]);
            b.push(meta.len() as u8); // len prefix (>= 27)
            b.extend_from_slice(&meta);
            b.extend_from_slice(&[0xff, 0xff, 0xff, 0xff]); // sentinel
            b.extend_from_slice(&0u32.to_le_bytes()); // hdr
            b.push(0x00);
            b.extend(std::iter::repeat_n(0x01u8, tl)); // tokens (0x01 = a peak)
        }
        let blocks = scan_blocks(&b);
        (b, blocks)
    }

    #[test]
    fn rebuild_grow_shifts_later_blocks() {
        // Three blocks; grow the first, shrink the second — the third's Idx offset must track.
        let (scan, blocks) = synth_scan(&[10, 10, 10]);
        assert_eq!(blocks.len(), 3, "three enumerated blocks");
        // Build a minimal 1-record Idx: o32 -> block0.ff, o86 -> block2.ff, o36 = block1.meta-block0.ff.
        let mut idx = vec![0u8; IDX_RS];
        let put = |idx: &mut [u8], off: usize, v: u32| idx[off..off + 4].copy_from_slice(&v.to_le_bytes());
        put(&mut idx, 32, (blocks[0].ff - 44) as u32);
        put(&mut idx, 86, (blocks[2].ff - 44) as u32);
        put(&mut idx, 36, (blocks[1].meta - blocks[0].ff) as u32);
        // Grow block0 (10 -> 20 tokens), shrink block1 (10 -> 4). Body = ff+9..end.
        let g0 = vec![0x01u8; blocks[0].end - (blocks[0].ff + 9) + 10];
        let g1 = vec![0x01u8; (blocks[1].end - (blocks[1].ff + 9)).saturating_sub(6)];
        let edits = vec![
            GrowEdit { block: 0, tokens: g0 },
            GrowEdit { block: 1, tokens: g1 },
        ];
        let (new_scan, bp) = rebuild_grow(&scan, &blocks, &edits).expect("grow");
        let new_idx = retranslate_idx(&idx, &bp).expect("retranslate");
        // The translated o32/o86/o36 must land on real block starts / boundaries in new_scan.
        let no32 = u32::from_le_bytes(new_idx[32..36].try_into().unwrap()) as usize;
        let no86 = u32::from_le_bytes(new_idx[86..90].try_into().unwrap()) as usize;
        assert_eq!(&new_scan[44 + no32..44 + no32 + 4], b"\xff\xff\xff\xff", "block0 sentinel");
        assert_eq!(&new_scan[44 + no86..44 + no86 + 4], b"\xff\xff\xff\xff", "block2 sentinel");
        // block2 shifted by (+10 from grow) + (-6 from shrink) = +4.
        assert_eq!(44 + no86, blocks[2].ff + 4, "block2 shifted by net delta");
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
        // Universal invariant: decode -> encode -> decode reproduces the SAME peaks. (Byte-identity
        // is NOT universal because the vendor encoder is non-canonical: some blocks use 0xfe/0xff
        // for delta 127/128, some store intensity 128..255 raw. We also track how many happen to be
        // byte-identical, which should be the majority — the MS1-canonical form.)
        let mut checked = 0;
        let mut byte_identical = 0;
        for b in &blocks {
            let slot = &sb[b.stream_start()..b.end];
            if slot.len() < 200 {
                continue;
            }
            let peaks = match decode_stream(&sb, b.stream_start(), 0, 40, false) {
                Ok(p) if p.len() >= 8 => p,
                _ => continue,
            };
            let enc = encode_stream(&peaks).expect("re-encode");
            // peaks round-trip: decoding our re-encoded stream must give the same peaks back
            let redec = decode_stream(&enc, 0, peaks[0].0, peaks.len(), false).expect("re-decode");
            assert_eq!(redec, peaks, "peaks round-trip mismatch in block at ff={}", b.ff);
            if enc.len() <= slot.len() && enc[..] == slot[..enc.len()] {
                byte_identical += 1;
            }
            checked += 1;
            if checked >= 50 {
                break;
            }
        }
        assert!(checked > 0, "no rich blocks were parity-checked");
        eprintln!("parity: {checked} blocks peaks-round-trip OK ({byte_identical} byte-identical)");
    }

    /// Idx-validity probe: for a real `.wiff`, check that every Idx record's block offsets (+32/+86,
    /// relative to byte 44) resolve to a real block `ffffffff` sentinel. If a minority DON'T, the Idx
    /// retranslation is corrupt for those blocks and pwiz would seek to the wrong byte. Gated on
    /// `TIMSIM_SCIEX_WIFF_SCAN`.
    #[test]
    fn probe_idx_validity() {
        let scan_path = match std::env::var("TIMSIM_SCIEX_WIFF_SCAN") {
            Ok(p) => p,
            Err(_) => return,
        };
        let wiff_path = scan_path.strip_suffix(".scan").expect("ends .scan").to_string();
        let scan = std::fs::read(&scan_path).expect("read .wiff.scan");
        let mut comp = cfb::open(&wiff_path).expect("open .wiff CFB");
        let mut idx = Vec::new();
        {
            use std::io::Read;
            comp.open_stream("SampleSubtree/Sample1/Idx")
                .expect("Idx stream")
                .read_to_end(&mut idx)
                .unwrap();
        }
        let blocks = scan_blocks(&scan);
        let ff_set: std::collections::HashSet<usize> = blocks.iter().map(|b| b.ff).collect();
        // Bodies of enumerated blocks: (ff+9 .. end). An Idx offset landing here (not at a known ff)
        // points at a mini/empty block that sits INSIDE an enumerated block's [ff+9..end] region —
        // which the grow rebuild overwrites and mis-translates.
        let mut bodies: Vec<(usize, usize)> = blocks.iter().map(|b| (b.ff + 9, b.end)).collect();
        bodies.sort();
        let inside_body = |pos: usize| bodies.iter().any(|&(s, e)| pos >= s && pos < e);
        let (mut offs, mut exact_ff, mut in_body, mut is_ffffffff, mut elsewhere) = (0, 0, 0, 0, 0);
        let mut samples: Vec<(usize, usize)> = Vec::new();
        for r in 0..idx.len() / IDX_RS {
            for field in [32usize, 86] {
                let off = u32le(&idx, r * IDX_RS + field).unwrap_or(0) as usize;
                if off == 0 {
                    continue;
                }
                offs += 1;
                let pos = 44 + off;
                if ff_set.contains(&pos) {
                    exact_ff += 1;
                } else if inside_body(pos) {
                    in_body += 1;
                    if scan.get(pos..pos + 4) == Some(&b"\xff\xff\xff\xff"[..]) {
                        is_ffffffff += 1; // a real (mini-block) ff sentinel, just not enumerated
                    }
                    if samples.len() < 6 {
                        samples.push((r, pos));
                    }
                } else {
                    elsewhere += 1;
                }
            }
        }
        eprintln!("IDX: {} blocks, {} nonzero offsets", blocks.len(), offs);
        eprintln!(
            "  exact block ff:      {} ({:.1}%)",
            exact_ff, 100.0 * exact_ff as f64 / offs as f64
        );
        eprintln!(
            "  INSIDE a block body: {} ({:.1}%)  <- points at a swallowed mini/empty block; {} are real ffffffff",
            in_body, 100.0 * in_body as f64 / offs as f64, is_ffffffff
        );
        eprintln!("  elsewhere:           {}", elsewhere);
        for (r, pos) in &samples {
            let ctx: Vec<String> = scan[*pos..(*pos + 8).min(scan.len())].iter().map(|b| format!("{b:02x}")).collect();
            eprintln!("  rec {} -> byte {} bytes=[{}]", r, pos, ctx.join(" "));
        }
    }
}
