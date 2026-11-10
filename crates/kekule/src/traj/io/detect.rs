use std::io::{Read, Seek, SeekFrom};
use std::path::Path;

use crate::traj::{
    TrajectoryFormat, TrajectoryIoError, TrajectoryIoErrorKind, TrajectoryIoOperation,
};

use super::{
    codec_context, io_context, FormatDetectionEvidence, TrajectoryFormatHint, TrajectoryIoLimits,
};

pub(crate) struct DetectionResult {
    pub(crate) format: TrajectoryFormat,
    pub(crate) evidence: Vec<FormatDetectionEvidence>,
}

/// Selects the format of the stream at `path`. A `growing` file too short to
/// show the signature its extension names is incomplete rather than of an
/// unknown format.
pub(crate) fn select_format<R: Read + Seek>(
    reader: &mut R,
    path: &Path,
    hint: TrajectoryFormatHint,
    limits: &TrajectoryIoLimits,
    source_label: &str,
    growing: bool,
) -> Result<DetectionResult, TrajectoryIoError> {
    let start = reader
        .stream_position()
        .map_err(|error| io_context(TrajectoryIoOperation::Detect, None, source_label, error))?;
    let max = limits.max_detection_bytes();
    if max == 0 {
        return Err(codec_context(
            TrajectoryIoErrorKind::LimitExceeded,
            TrajectoryIoOperation::Detect,
            None,
            source_label,
            "The detection limit is 0 bytes. Detection reads at least 1 byte.",
        ));
    }
    if max > limits.max_scratch_bytes() {
        return Err(codec_context(
            TrajectoryIoErrorKind::LimitExceeded,
            TrajectoryIoOperation::Detect,
            None,
            source_label,
            format!(
                "The detection limit of {max} bytes exceeds the buffer limit of {} bytes.",
                limits.max_scratch_bytes()
            ),
        ));
    }
    let mut prefix = Vec::new();
    prefix.try_reserve_exact(max).map_err(|_| {
        codec_context(
            TrajectoryIoErrorKind::LimitExceeded,
            TrajectoryIoOperation::Detect,
            None,
            source_label,
            "Cannot allocate the detection buffer.",
        )
    })?;
    let read_result = Read::by_ref(reader)
        .take(max as u64)
        .read_to_end(&mut prefix);
    let restore_result = reader
        .seek(SeekFrom::Start(start))
        .map_err(|error| io_context(TrajectoryIoOperation::Detect, None, source_label, error));
    restore_result?;
    read_result
        .map_err(|error| io_context(TrajectoryIoOperation::Detect, None, source_label, error))?;

    if is_compressed_wrapper(&prefix) {
        return Err(codec_context(
            TrajectoryIoErrorKind::UnsupportedVariant,
            TrajectoryIoOperation::Detect,
            None,
            source_label,
            "The file is compressed. Decompress it first.",
        ));
    }

    let signature = signature_format(&prefix);
    let extension = extension_format(path);
    match hint {
        TrajectoryFormatHint::Explicit(format) => Ok(DetectionResult {
            format,
            evidence: vec![FormatDetectionEvidence::ExplicitHint],
        }),
        TrajectoryFormatHint::Auto => match (signature, extension) {
            (Some(signature), Some(extension)) if signature == extension => Ok(DetectionResult {
                format: signature,
                evidence: vec![
                    FormatDetectionEvidence::Signature,
                    FormatDetectionEvidence::Extension,
                    FormatDetectionEvidence::ExtensionSignatureAgreement,
                ],
            }),
            (Some(signature), Some(extension)) => Err(codec_context(
                TrajectoryIoErrorKind::FormatMismatch,
                TrajectoryIoOperation::Detect,
                Some(signature),
                source_label,
                format!("The content is {signature} but the extension names {extension}."),
            )),
            (Some(signature), None) => Ok(DetectionResult {
                format: signature,
                evidence: vec![
                    FormatDetectionEvidence::Signature,
                    FormatDetectionEvidence::MissingExtension,
                ],
            }),
            (None, Some(extension)) if growing && prefix.len() < max => Err(codec_context(
                TrajectoryIoErrorKind::Incomplete,
                TrajectoryIoOperation::Detect,
                Some(extension),
                source_label,
                "The file does not yet hold its format's signature.",
            )),
            (None, Some(extension)) => Err(codec_context(
                TrajectoryIoErrorKind::UnknownFormat,
                TrajectoryIoOperation::Detect,
                Some(extension),
                source_label,
                "The content does not match the format its extension names.",
            )),
            (None, None) => Err(codec_context(
                TrajectoryIoErrorKind::UnknownFormat,
                TrajectoryIoOperation::Detect,
                None,
                source_label,
                "The content matches no supported trajectory format.",
            )),
        },
    }
}

fn is_compressed_wrapper(prefix: &[u8]) -> bool {
    prefix.starts_with(&[0x1f, 0x8b])
        || prefix.starts_with(&[0xfd, b'7', b'z', b'X', b'Z', 0x00])
        || prefix.starts_with(b"BZh")
}

fn signature_format(prefix: &[u8]) -> Option<TrajectoryFormat> {
    if prefix.len() >= 8 {
        let marker = &prefix[..4];
        if &prefix[4..8] == b"CORD"
            && (marker == 84u32.to_le_bytes() || marker == 84u32.to_be_bytes())
        {
            return Some(TrajectoryFormat::Dcd);
        }
    }
    if prefix.len() >= 4 {
        let magic = i32::from_be_bytes(prefix[..4].try_into().ok()?);
        if matches!(magic, 1995 | 2023) {
            return Some(TrajectoryFormat::Xtc);
        }
        if magic == 1993 {
            return Some(TrajectoryFormat::Trr);
        }
    }
    xyz_signature(prefix).then_some(TrajectoryFormat::Xyz)
}

fn xyz_signature(prefix: &[u8]) -> bool {
    let Ok(text) = std::str::from_utf8(prefix) else {
        return false;
    };
    let mut lines = text.lines();
    let Some(count) = lines.next() else {
        return false;
    };
    let Ok(count) = count.trim().parse::<usize>() else {
        return false;
    };
    if count == 0 || lines.next().is_none() {
        return false;
    }
    let Some(atom) = lines.next() else {
        return false;
    };
    let mut fields = atom.split_whitespace();
    let (Some(symbol), Some(x), Some(y), Some(z), None) = (
        fields.next(),
        fields.next(),
        fields.next(),
        fields.next(),
        fields.next(),
    ) else {
        return false;
    };
    crate::molecule::Element::from_symbol(symbol).is_some()
        && [x, y, z]
            .into_iter()
            .all(|value| value.parse::<f64>().is_ok_and(f64::is_finite))
}

pub(super) fn extension_format(path: &Path) -> Option<TrajectoryFormat> {
    let extension = path.extension()?.to_str()?.to_ascii_lowercase();
    match extension.as_str() {
        "xyz" => Some(TrajectoryFormat::Xyz),
        "dcd" => Some(TrajectoryFormat::Dcd),
        "xtc" => Some(TrajectoryFormat::Xtc),
        "trr" => Some(TrajectoryFormat::Trr),
        _ => None,
    }
}
