//! What went wrong in a [`TrajectoryIoError`](super::TrajectoryIoError),
//! beyond its kind: the typed facts its message is written from.

use std::fmt;

use crate::diagnostics::Capitalized;
use crate::molecule::Element;

use super::TrajectoryFormat;

/// A part of a trajectory file, frame, or codec buffer that a [`Reason`]
/// names. It displays as a noun phrase with its article.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Part {
    AtomCount,
    Block,
    BlockSize,
    /// The bytes of a TRR data block.
    BlockOf(Datum),
    BufferOf(Datum),
    Buffers,
    CompressedByteCount,
    CompressedCoordinateHeader,
    CompressedCoordinates,
    CompressedRunAtomCount,
    CompressedSize,
    Coordinate,
    CoordinateBuffer,
    CoordinateCount,
    CoordinateRecord,
    CoordinateRecordSize,
    Coordinates,
    DecodedCoordinate,
    DetectionBuffer,
    DetectionLimit,
    EnergyCount,
    Frame,
    FrameBuffers,
    FrameCount,
    FrameEnd,
    FrameHeader,
    FrameHeaderSize,
    FrameIndex,
    FrameIndexTimes,
    FrameNumber,
    FrameSize,
    FrameStep,
    FrameTime,
    FrameValues,
    FrameVectors,
    FreeAtomIndex,
    HeaderOffset,
    Index,
    Lambda,
    Line,
    MagicNumber,
    Nset,
    Padding,
    PositionBuffer,
    Positions,
    Record,
    RepeatedAtomCount,
    RecordLength,
    RecordMarker,
    RecordPayload,
    Step,
    Time,
    TitleCount,
    TrailingRecordMarker,
    Trajectory,
    /// A single value of a quantity, such as a position value.
    ValueOf(Datum),
    VectorBlock,
    VectorBuffers,
    Vectors,
    VersionLength,
    VersionString,
    VersionStringLength,
    VersionStringPadding,
    WriterBuffers,
    Cell,
}

impl Part {
    /// Whether the phrase takes a plural verb.
    const fn plural(self) -> bool {
        matches!(
            self,
            Self::Buffers
                | Self::CompressedCoordinates
                | Self::Coordinates
                | Self::FrameBuffers
                | Self::FrameIndexTimes
                | Self::FrameValues
                | Self::FrameVectors
                | Self::Positions
                | Self::VectorBuffers
                | Self::Vectors
                | Self::WriterBuffers
        )
    }

    /// `singular` or `plural` to agree with the phrase.
    pub(crate) const fn verb(self, singular: &'static str, plural: &'static str) -> &'static str {
        if self.plural() {
            plural
        } else {
            singular
        }
    }
}

impl fmt::Display for Part {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let phrase = match self {
            Self::BlockOf(datum) => return write!(formatter, "the {datum} block"),
            Self::BufferOf(datum) => return write!(formatter, "the {datum} buffer"),
            Self::ValueOf(datum) => return write!(formatter, "a {datum} value"),
            Self::AtomCount => "the atom count",
            Self::Block => "a block",
            Self::BlockSize => "a block size",
            Self::Buffers => "the buffers",
            Self::CompressedByteCount => "the compressed-byte count",
            Self::CompressedCoordinateHeader => "the compressed-coordinate header",
            Self::CompressedCoordinates => "the compressed coordinates",
            Self::CompressedRunAtomCount => "a compressed run's atom count",
            Self::CompressedSize => "the compressed coordinates' size",
            Self::Coordinate => "a coordinate",
            Self::CoordinateBuffer => "the coordinate buffer",
            Self::CoordinateCount => "the coordinate count",
            Self::CoordinateRecord => "a coordinate record",
            Self::CoordinateRecordSize => "a coordinate record's size",
            Self::Coordinates => "the coordinates",
            Self::DecodedCoordinate => "a decoded coordinate",
            Self::DetectionBuffer => "the detection buffer",
            Self::DetectionLimit => "the detection limit",
            Self::EnergyCount => "the energy count nre",
            Self::Frame => "the frame",
            Self::FrameBuffers => "the frame's buffers",
            Self::FrameCount => "the frame count",
            Self::FrameEnd => "the frame's end",
            Self::FrameHeader => "the frame header",
            Self::FrameHeaderSize => "the frame header's size",
            Self::FrameIndex => "the frame index",
            Self::FrameIndexTimes => "the frame index's times or steps",
            Self::FrameNumber => "the frame number",
            Self::FrameSize => "the frame's size",
            Self::FrameStep => "the frame's step",
            Self::FrameTime => "the frame's time",
            Self::FrameValues => "the frame's values",
            Self::FrameVectors => "the frame's vectors",
            Self::FreeAtomIndex => "a free-atom index",
            Self::HeaderOffset => "the header's byte offset",
            Self::Index => "the index",
            Self::Lambda => "the lambda value",
            Self::Line => "a line",
            Self::MagicNumber => "the magic number",
            Self::Nset => "the DCD NSET",
            Self::Padding => "the padding",
            Self::PositionBuffer => "the position buffer",
            Self::Positions => "the positions",
            Self::Record => "a record",
            Self::RepeatedAtomCount => "the repeated atom count",
            Self::RecordLength => "a record length",
            Self::RecordMarker => "the DCD record marker",
            Self::RecordPayload => "the DCD record payload",
            Self::Step => "the step",
            Self::Time => "the time",
            Self::TitleCount => "the title count NTITLE",
            Self::TrailingRecordMarker => "the DCD trailing record marker",
            Self::Trajectory => "the trajectory",
            Self::VectorBlock => "a vector block",
            Self::VectorBuffers => "the vector buffers",
            Self::Vectors => "the vectors",
            Self::VersionLength => "the version length",
            Self::VersionString => "the version string",
            Self::VersionStringLength => "the version string's length",
            Self::VersionStringPadding => "the version string's padding",
            Self::WriterBuffers => "the writer's buffers",
            Self::Cell => "the cell",
        };
        formatter.write_str(phrase)
    }
}

/// A quantity a codec reads or writes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Datum {
    Box,
    Cell,
    Position,
    Velocity,
    Force,
    Time,
    Step,
    Lambda,
}

impl Datum {
    /// The plural noun, for the values of every atom.
    const fn plural(self) -> &'static str {
        match self {
            Self::Box => "boxes",
            Self::Cell => "cell",
            Self::Position => "positions",
            Self::Velocity => "velocities",
            Self::Force => "forces",
            Self::Time => "times",
            Self::Step => "steps",
            Self::Lambda => "lambda values",
        }
    }
}

impl fmt::Display for Datum {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Box => "box",
            Self::Cell => "cell",
            Self::Position => "position",
            Self::Velocity => "velocity",
            Self::Force => "force",
            Self::Time => "time",
            Self::Step => "step",
            Self::Lambda => "lambda",
        })
    }
}

/// The 32-bit field a value does not fit.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Field32 {
    RecordLength,
    BlockSize,
    AtomCountField,
    StepField,
    LengthField,
}

impl fmt::Display for Field32 {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::RecordLength => "record length",
            Self::BlockSize => "block size",
            Self::AtomCountField => "atom-count field",
            Self::StepField => "step field",
            Self::LengthField => "length field",
        })
    }
}

/// A DCD header control field.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum DcdField {
    Nset,
    Istart,
    Nsavc,
    Namnf,
    Natom,
    Ntitle,
}

impl fmt::Display for DcdField {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::Nset => "NSET",
            Self::Istart => "ISTART",
            Self::Nsavc => "NSAVC",
            Self::Namnf => "NAMNF",
            Self::Natom => "NATOM",
            Self::Ntitle => "NTITLE",
        })
    }
}

/// A DCD header record of fixed size.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum DcdRecord {
    AtomCount,
    FreeAtomIndex,
    Control,
}

impl fmt::Display for DcdRecord {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::AtomCount => "atom-count record",
            Self::FreeAtomIndex => "free-atom index record",
            Self::Control => "control block",
        })
    }
}

/// The quantity whose unit a codec cannot convert.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum UnitOf {
    /// The time unit of a DCD time policy.
    TimePolicy,
    Of(Datum),
    /// The values of a quantity for every atom.
    AllOf(Datum),
}

impl fmt::Display for UnitOf {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::TimePolicy => formatter.write_str("time-policy"),
            Self::Of(datum) => write!(formatter, "{datum}"),
            Self::AllOf(datum) => formatter.write_str(datum.plural()),
        }
    }
}

/// What a frame lacks a datum for.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Need {
    /// Every frame of the format.
    Format(TrajectoryFormat),
    /// The time policy of a DCD header.
    TimePolicy,
    /// The unit-cell flag of a DCD header.
    CellFlag,
}

/// The atom counts a format or a header holds.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum AtomRange {
    DcdHeader,
    Trr,
    Xtc,
    Xtc1995,
}

/// An encoded frame index or file fingerprint.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Encoding {
    FileFingerprint,
    FrameIndex,
}

impl Encoding {
    const fn magic(self) -> &'static str {
        match self {
            Self::FileFingerprint => "KEKTRFPR",
            Self::FrameIndex => "KEKTRIDX",
        }
    }
}

impl fmt::Display for Encoding {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::FileFingerprint => "file fingerprint",
            Self::FrameIndex => "frame index",
        })
    }
}

/// A corruption of XTC compressed coordinates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Compressed {
    BoundsReversed,
    SmallIndexRange,
    SmallIndexOutsideTable,
    NoBytes,
    EndsBeforeLastAtom,
    SpanTooWide,
    Axis32Bit,
    OutsideDeclaredBounds,
    ExceedsDeclaredBounds,
    RunPastAtomCount,
    RunExceedsBounds,
    SmallIndexBelowTable,
    SmallIndexAboveTable,
    TrailingBytes,
    TrailingBits,
    ValueOver32Bits,
    AtomCountDiffers,
    TripleExceedsBounds,
    XExceedsBounds,
    OutsideBounds,
}

impl fmt::Display for Compressed {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        formatter.write_str(match self {
            Self::BoundsReversed => {
                "The compressed coordinates' lower bound exceeds the upper bound."
            }
            Self::SmallIndexRange => {
                "The compressed coordinates' small index is negative or too large."
            }
            Self::SmallIndexOutsideTable => {
                "The compressed coordinates' small index is outside the table."
            }
            Self::NoBytes => "The compressed coordinates have no bytes.",
            Self::EndsBeforeLastAtom => "The compressed coordinates end before the last atom.",
            Self::SpanTooWide => "The compressed coordinates span more than the decoder supports.",
            Self::Axis32Bit => {
                "A compressed axis needs 32-bit values. The decoder does not support them."
            }
            Self::OutsideDeclaredBounds => {
                "A compressed coordinate is outside the declared bounds."
            }
            Self::ExceedsDeclaredBounds => "A compressed coordinate exceeds the declared bounds.",
            Self::RunPastAtomCount => "A compressed run extends past the declared atom count.",
            Self::RunExceedsBounds => "A compressed run's coordinate exceeds the declared bounds.",
            Self::SmallIndexBelowTable => {
                "The compressed coordinates' small index drops below the table."
            }
            Self::SmallIndexAboveTable => {
                "The compressed coordinates' small index rises above the table."
            }
            Self::TrailingBytes => "The compressed coordinates have unused trailing bytes.",
            Self::TrailingBits => "The compressed coordinates have nonzero trailing bits.",
            Self::ValueOver32Bits => "A compressed value needs more than 32 bits.",
            Self::AtomCountDiffers => {
                "The compressed coordinates decode to fewer or more atoms than declared."
            }
            Self::TripleExceedsBounds => "A compressed coordinate triple exceeds its bounds.",
            Self::XExceedsBounds => "A compressed x coordinate exceeds its bounds.",
            Self::OutsideBounds => "A compressed coordinate is outside its bounds.",
        })
    }
}

/// What went wrong, in the terms its message states.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum Reason {
    // Families over a part.
    Truncated(Part),
    Missing(Part),
    Incomplete(Part),
    Negative(Part),
    NotFinite(Part),
    NotZero(Part),
    AddressRange(Part),
    Range64(Part),
    Exceeds32Bit(Part, Field32),
    Allocation(Part),
    /// The error's [`LimitExceeded`](crate::diagnostics::LimitExceeded)
    /// names the limit.
    Exceeds(Part),
    DoesNotFitSingle(Part),
    DoesNotFitSelected(Part),
    IncompatibleUnit(Option<TrajectoryFormat>, UnitOf),
    XyzLengthUnitNotLength,
    FrameLacks(Datum, Need),
    NotFullyPeriodic(TrajectoryFormat),
    AtomCountOutOfRange(usize, AtomRange),
    Compressed(Compressed),

    // Detection, opening, and writing files.
    DetectionLimitZero,
    CompressedFile,
    ContentExtensionMismatch(TrajectoryFormat, TrajectoryFormat),
    NoSignatureYet,
    NotExtensionFormat,
    NoCompleteFrameYet,
    EarlierFrameFailed,
    StreamNeedsFormat,
    UnknownExtension,
    NoHardLinks,
    NoTemporaryName,
    NoFrames(TrajectoryFormat),
    PartialWrite(TrajectoryFormat),
    FrameIntervalInvalid(TrajectoryFormat),
    ColumnMismatch,
    TopologyHasNoAtoms,
    EmptyStream,
    FileHoldsNoFrame,

    // Appending and cutting.
    XyzFramesHaveNoStep,
    FileChangedBeforeCut,
    KeepFrames {
        requested: usize,
        available: usize,
    },
    AppendHeaderDiffers,
    AppendMixedTrrWidths,
    AppendTrrPrecision {
        file: &'static str,
        options: &'static str,
    },
    AppendXtcMagic {
        file: i32,
        options: i32,
    },
    AppendXtcPrecision {
        file: f32,
        options: f32,
    },

    // Frame indexes and fingerprints.
    FrameGone {
        frame: usize,
        length: u64,
        end: u64,
    },
    IndexedBytesChanged,
    IndexForOtherFormat(TrajectoryFormat),
    PrefixMismatch,
    ChangedWhileIndexed,
    TrailingBytes(Encoding),
    ReservedByte(Encoding),
    MalformedModificationTime(Encoding),
    FingerprintTextLength,
    FingerprintTextCase,
    StampLayout(TrajectoryFormat),
    StampMismatch,
    FrameEndMismatch {
        actual: u64,
        expected: u64,
    },
    IndexStartMismatch,
    FileShorterThanIndex {
        length: u64,
        frames: usize,
        end: u64,
    },
    UnknownFormatCode(u8),
    InvalidFlags {
        flags: u8,
        format: TrajectoryFormat,
    },
    FlagsDisagree {
        flags: u8,
        frames: usize,
        format: TrajectoryFormat,
    },
    IndexLength,
    OffsetOutOfOrder(usize),
    IndexTimeNotFinite(usize),
    EncodingTruncated(Encoding),
    NoMagic(Encoding),
    UnsupportedVersion(Encoding, u32),
    ChecksumMismatch(Encoding),
    IndexFull(usize),
    BuffersWouldExceed(usize),
    AllocatedBuffersExceed(usize),

    // Header and frame contents.
    NoAtomsDeclared,
    FixedAtomsExceed {
        fixed: usize,
        atoms: usize,
    },
    HeaderFieldNegative(DcdField, i32),
    HeaderFieldNotPositive(DcdField, i32),
    RecordSize {
        record: DcdRecord,
        bytes: usize,
        expected: usize,
    },
    TitleRecordShort(usize),
    TitleRecordSize {
        bytes: usize,
        titles: usize,
        expected: usize,
    },
    FreeAtomIndexNotPositive(i32),
    FreeAtomIndexInvalid {
        index: i32,
        atoms: usize,
    },
    FirstMarkerNot84,
    NotCordHeader,
    FourDimensional,
    UnsupportedDcdLayout,
    UnitCellValues,
    CoordinateRecordSize,
    FrameEndsBeforeStart,
    FrameEndsAfterFile,
    StepRange64,
    StepNotInSequence(u64),
    TimeNotFromDelta,
    StepIntervalZero,
    StepSequenceRange,
    DeltaNotTime,
    InvalidDelta(f32),
    FrameCountField,
    CellNotFinite,
    BoxNotFinite,
    BoxZero,
    CellLengths {
        a: f64,
        b: f64,
        c: f64,
    },
    DegenerateCellAngles,
    DegenerateCellVectors,
    InvalidDcdCell,
    InvalidPeriodicCell,
    CellNotStandard,
    CellNoParameters,
    CellDegenerateParameters,
    FrameChanged,
    DecodedFrameDiffers,
    MagicOrPrecisionDiffers,
    IndexedFrameGone,
    RepeatedAtomCountDiffers,
    XtcMagic,
    TrrMagic,
    PrecisionNotInverseLength,
    PrecisionNotPositive,
    CoordinateOutOfIntegerRange,
    SpanOutOfIntegerRange,
    Xtc1995ByteCount,
    TrrVersionLength,
    TrrVersionString,
    VersionStringTooLong(usize),
    UnsupportedTrrBlocks,
    NoPositions,
    BlockSizesMismatch,
    EnergyTerms,
    TimeOrLambdaNotFinite {
        time: f64,
        lambda: f64,
    },
    NoPrecisionBlock,
    NoScalarWidth,
    ScalarWidth(usize),
    ScalarWidthDiffers,
    LambdaNotProperty,
    LambdaNotZero,
    LambdaPolicyKey,
    LambdaPropertyValue,
    AtomCountLine,
    ZeroAtoms,
    NoCommentLine,
    EndsBeforeAtomLine(usize),
    AtomLineFields(usize),
    ExtraColumns(usize),
    ContentAfterBlankLine,
    UnknownElement {
        line: usize,
        symbol: Box<str>,
    },
    ElementMismatch {
        atom: usize,
        found: Element,
        expected: Element,
    },
    CoordinateNotNumber(usize),
    CoordinateNotFiniteOnLine(usize),
    LineBufferGrowth(Box<Reason>),
    NoLineBreak,
    NotFiniteInOutputUnit,
}

impl fmt::Display for Reason {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        let text = match self {
            Self::Truncated(part) => {
                let verb = part.verb("is", "are");
                return write!(formatter, "{} {verb} truncated.", Capitalized(part));
            }
            Self::Missing(part) => return write!(formatter, "{} is missing.", Capitalized(part)),
            Self::Incomplete(part) => {
                return write!(formatter, "{} is incomplete.", Capitalized(part))
            }
            Self::Negative(part) => return write!(formatter, "{} is negative.", Capitalized(part)),
            Self::NotFinite(part) => {
                let verb = part.verb("is", "are");
                return write!(formatter, "{} {verb} not finite.", Capitalized(part));
            }
            Self::NotZero(part) => return write!(formatter, "{} is not zero.", Capitalized(part)),
            Self::AddressRange(part) => {
                let verb = part.verb("exceeds", "exceed");
                return write!(
                    formatter,
                    "{} {verb} this platform's address range.",
                    Capitalized(part)
                );
            }
            Self::Range64(part) => {
                return write!(formatter, "{} exceeds the 64-bit range.", Capitalized(part))
            }
            Self::Exceeds32Bit(part, field) => {
                let verb = part.verb("exceeds", "exceed");
                return write!(
                    formatter,
                    "{} {verb} the 32-bit {field}.",
                    Capitalized(part)
                );
            }
            Self::Allocation(part) => return write!(formatter, "Cannot allocate {part}."),
            Self::Exceeds(part) => {
                let verb = part.verb("exceeds", "exceed");
                return write!(formatter, "{} {verb} its limit.", Capitalized(part));
            }
            Self::DoesNotFitSingle(part) => {
                return write!(
                    formatter,
                    "{} does not fit single precision.",
                    Capitalized(part)
                )
            }
            Self::DoesNotFitSelected(part) => {
                return write!(
                    formatter,
                    "{} does not fit the selected precision.",
                    Capitalized(part)
                )
            }
            Self::IncompatibleUnit(Some(format), of) => {
                return write!(formatter, "The {format} {of} unit is incompatible.")
            }
            Self::IncompatibleUnit(None, of) => {
                return write!(formatter, "The {of} unit is incompatible.")
            }
            Self::XyzLengthUnitNotLength => "The XYZ length unit is not a length.",
            Self::FrameLacks(datum, need) => {
                let datum = match datum {
                    Datum::Cell | Datum::Box => "periodic cell".into(),
                    datum => datum.to_string(),
                };
                return match need {
                    Need::Format(format) => write!(
                        formatter,
                        "The frame has no {datum}. {format} frames need one."
                    ),
                    Need::TimePolicy => write!(
                        formatter,
                        "The frame has no {datum}. The header's time policy needs one."
                    ),
                    Need::CellFlag => write!(
                        formatter,
                        "The frame has no {datum}. The header's unit-cell flag needs one."
                    ),
                };
            }
            Self::NotFullyPeriodic(format) => {
                return write!(
                    formatter,
                    "The cell is not periodic along all three axes. {format} stores only fully \
                     periodic cells."
                )
            }
            Self::AtomCountOutOfRange(atoms, range) => {
                return match range {
                    AtomRange::DcdHeader => write!(
                        formatter,
                        "The trajectory has {atoms} atoms. The header holds 1 to 2,147,483,647."
                    ),
                    AtomRange::Trr => write!(
                        formatter,
                        "The topology has {atoms} atoms. TRR holds 1 to 715,827,882."
                    ),
                    AtomRange::Xtc => write!(
                        formatter,
                        "The topology has {atoms} atoms. XTC holds 1 to 2,147,483,647."
                    ),
                    AtomRange::Xtc1995 => write!(
                        formatter,
                        "The topology has {atoms} atoms. XTC 1995 holds at most 298,261,617. \
                         Use XTC 2023."
                    ),
                };
            }
            Self::Compressed(fault) => return fault.fmt(formatter),

            Self::DetectionLimitZero => {
                "The detection limit is 0 bytes. Detection reads at least 1 byte."
            }
            Self::CompressedFile => "The file is compressed. Decompress it first.",
            Self::ContentExtensionMismatch(content, extension) => {
                return write!(
                    formatter,
                    "The content is {content} but the extension names {extension}."
                )
            }
            Self::NoSignatureYet => "The file does not yet hold its format's signature.",
            Self::NotExtensionFormat => {
                "The content does not match the format its extension names."
            }
            Self::NoCompleteFrameYet => {
                "The file does not yet hold a complete frame. Open it again once more is written."
            }
            Self::EarlierFrameFailed => {
                "An earlier frame failed to write. The trajectory is incomplete."
            }
            Self::StreamNeedsFormat => "A stream has no extension. Choose the format explicitly.",
            Self::UnknownExtension => {
                "The destination's extension names no supported format. Choose the format \
                 explicitly."
            }
            Self::NoHardLinks => {
                "The destination's file system has no hard links. Writing without overwriting \
                 needs them."
            }
            Self::NoTemporaryName => "Cannot find a free name for the temporary file.",
            Self::NoFrames(format) => {
                return write!(
                    formatter,
                    "The {format} trajectory has no frame. It needs at least one."
                )
            }
            Self::PartialWrite(format) => {
                return write!(
                    formatter,
                    "An earlier {format} write failed partway. The stream is incomplete."
                )
            }
            Self::FrameIntervalInvalid(format) => {
                return write!(
                    formatter,
                    "The {format} frame interval is not a finite, nonnegative time."
                )
            }
            Self::ColumnMismatch => {
                "The frame has a time or step the first frame lacks, or lacks one it has."
            }
            Self::TopologyHasNoAtoms => "The topology has no atoms.",
            Self::EmptyStream => "The stream is empty.",
            Self::FileHoldsNoFrame => "The file holds no frame.",

            Self::XyzFramesHaveNoStep => "XYZ frames have no step. Keep frames by count.",
            Self::FileChangedBeforeCut => "The file changed before it was cut.",
            Self::KeepFrames {
                requested,
                available,
            } => {
                return write!(
                    formatter,
                    "Cannot keep {requested} frames. The file holds {available} complete frames."
                )
            }
            Self::AppendHeaderDiffers => {
                "The file's header differs from the header of the write options. Append with the \
                 options the file was written with."
            }
            Self::AppendMixedTrrWidths => {
                "The file's TRR frames mix scalar widths. One writer cannot continue them."
            }
            Self::AppendTrrPrecision { file, options } => {
                return write!(
                    formatter,
                    "The file's TRR frames are in {file}, the write options in {options}. Append \
                     with the options the file was written with."
                )
            }
            Self::AppendXtcMagic { file, options } => {
                return write!(
                    formatter,
                    "The file's XTC frames use magic {file} and the write options {options}. \
                     Append with the options the file was written with."
                )
            }
            Self::AppendXtcPrecision { file, options } => {
                return write!(
                    formatter,
                    "The file's XTC frames have precision {file} per nm and the write options \
                     {options} per nm. Append with the options the file was written with."
                )
            }

            Self::FrameGone { frame, length, end } => {
                return write!(
                    formatter,
                    "The file no longer holds frame {frame}: it is {length} bytes, and the frame \
                     ends at byte {end}. Open the file again."
                )
            }
            Self::IndexedBytesChanged => {
                "The file's indexed bytes changed since they were indexed. Open the file again."
            }
            Self::IndexForOtherFormat(format) => {
                return write!(formatter, "The index was built for a {format} file.")
            }
            Self::PrefixMismatch => "The file does not begin with the bytes the index covers.",
            Self::ChangedWhileIndexed => "The file changed while it was indexed.",
            Self::TrailingBytes(encoding) => {
                return write!(formatter, "The {encoding} has trailing bytes.")
            }
            Self::ReservedByte(encoding) => {
                return write!(formatter, "A reserved byte of the {encoding} is not zero.")
            }
            Self::MalformedModificationTime(encoding) => {
                return write!(formatter, "The {encoding} modification time is malformed.")
            }
            Self::FingerprintTextLength => {
                "The file fingerprint text is not 112 hexadecimal digits."
            }
            Self::FingerprintTextCase => "The file fingerprint text is not lowercase hexadecimal.",
            Self::StampLayout(format) => {
                return write!(
                    formatter,
                    "The {format} frame records time and step differently from frame 0."
                )
            }
            Self::StampMismatch => "The frame's time or step differs from the index's.",
            Self::FrameEndMismatch { actual, expected } => {
                return write!(
                    formatter,
                    "The frame ends at byte {actual}, but the index places the next frame at \
                     byte {expected}."
                )
            }
            Self::IndexStartMismatch => "The index does not start at the file's first frame.",
            Self::FileShorterThanIndex {
                length,
                frames,
                end,
            } => {
                return write!(
                    formatter,
                    "The file is {length} bytes, and its first {frames} frames end at byte {end}."
                )
            }
            Self::UnknownFormatCode(code) => {
                return write!(
                    formatter,
                    "The frame index names the unknown format {code}."
                )
            }
            Self::InvalidFlags { flags, format } => {
                return write!(
                    formatter,
                    "The frame index flags {flags:#04x} are invalid for {format}."
                )
            }
            Self::FlagsDisagree {
                flags,
                frames,
                format,
            } => {
                return write!(
                    formatter,
                    "The frame index flags {flags:#04x} disagree with the steps and times its \
                     {frames} {format} frames record."
                )
            }
            Self::IndexLength => "The frame index's length does not match its frame count.",
            Self::OffsetOutOfOrder(entry) => {
                return write!(
                    formatter,
                    "The frame index offset {entry} is out of order or outside the file."
                )
            }
            Self::IndexTimeNotFinite(entry) => {
                return write!(formatter, "The frame index time {entry} is not finite.")
            }
            Self::EncodingTruncated(encoding) => {
                return write!(formatter, "The {encoding} is truncated.")
            }
            Self::NoMagic(encoding) => {
                return write!(
                    formatter,
                    "The {encoding} has no {} magic.",
                    encoding.magic()
                )
            }
            Self::UnsupportedVersion(encoding, version) => {
                return write!(
                    formatter,
                    "The {encoding} version {version} is not supported."
                )
            }
            Self::ChecksumMismatch(encoding) => {
                return write!(formatter, "The {encoding} checksum does not match.")
            }
            Self::IndexFull(frames) => {
                return write!(formatter, "The frame index is full at {frames} frames.")
            }
            Self::BuffersWouldExceed(limit) => {
                return write!(
                    formatter,
                    "The buffers would exceed their limit of {limit} bytes."
                )
            }
            Self::AllocatedBuffersExceed(limit) => {
                return write!(
                    formatter,
                    "The allocated buffers exceed their limit of {limit} bytes."
                )
            }

            Self::NoAtomsDeclared => "The header declares no atoms.",
            Self::FixedAtomsExceed { fixed, atoms } => {
                return write!(
                    formatter,
                    "The header declares {fixed} fixed atoms but only {atoms} atoms."
                )
            }
            Self::HeaderFieldNegative(field, value) => {
                return write!(
                    formatter,
                    "The header field {field} is {value}. It is negative."
                )
            }
            Self::HeaderFieldNotPositive(field, value) => {
                return write!(
                    formatter,
                    "The header field {field} is {value}. It is not positive."
                )
            }
            Self::RecordSize {
                record,
                bytes,
                expected,
            } => return write!(formatter, "The {record} has {bytes} bytes, not {expected}."),
            Self::TitleRecordShort(bytes) => {
                return write!(
                    formatter,
                    "The title record has {bytes} bytes, fewer than 4."
                )
            }
            Self::TitleRecordSize {
                bytes,
                titles,
                expected,
            } => {
                return write!(
                    formatter,
                    "The title record has {bytes} bytes. Its {titles} titles need {expected}."
                )
            }
            Self::FreeAtomIndexNotPositive(index) => {
                return write!(formatter, "The free-atom index {index} is not positive.")
            }
            Self::FreeAtomIndexInvalid { index, atoms } => {
                return write!(
                    formatter,
                    "The free-atom index {index} is outside the {atoms} atoms or repeated."
                )
            }
            Self::FirstMarkerNot84 => "The first record marker is not 84 in either byte order.",
            Self::NotCordHeader => "The header is not an 84-byte CORD record.",
            Self::FourDimensional => {
                "Four-dimensional and charge coordinate records are not supported."
            }
            Self::UnsupportedDcdLayout => {
                "The header is not in the CHARMM, NAMD or OpenMM layout. Other DCD layouts are \
                 not supported."
            }
            Self::UnitCellValues => {
                "The unit-cell record does not hold six double-precision values."
            }
            Self::CoordinateRecordSize => "A coordinate record has the wrong size.",
            Self::FrameEndsBeforeStart => "The frame ends before it starts.",
            Self::FrameEndsAfterFile => "The frame ends after the end of the file.",
            Self::StepRange64 => "The frame's step exceeds the 64-bit step range.",
            Self::StepNotInSequence(step) => {
                return write!(
                    formatter,
                    "The frame's step is not {step}, the next step of the header's sequence."
                )
            }
            Self::TimeNotFromDelta => {
                "The frame's time differs from the time the step and DELTA give."
            }
            Self::StepIntervalZero => "The step interval is zero.",
            Self::StepSequenceRange => {
                "The first step or the step interval exceeds the 32-bit header field."
            }
            Self::DeltaNotTime => "The DCD DELTA is not a time.",
            Self::InvalidDelta(delta) => {
                return write!(
                    formatter,
                    "The time step DELTA is {delta}. It is not a finite, nonnegative number."
                )
            }
            Self::FrameCountField => {
                "The trajectory has more frames than the 32-bit frame count of the header holds."
            }
            Self::CellNotFinite => "The cell holds a value that is not finite.",
            Self::BoxNotFinite => "The box holds a value that is not finite.",
            Self::BoxZero => "The box is a zero matrix.",
            Self::CellLengths { a, b, c } => {
                return write!(
                    formatter,
                    "The cell lengths {a} Å, {b} Å and {c} Å are not all positive."
                )
            }
            Self::DegenerateCellAngles => "The cell angles are degenerate.",
            Self::DegenerateCellVectors => "The cell vectors are degenerate.",
            Self::InvalidDcdCell => "The DCD cell is invalid.",
            Self::InvalidPeriodicCell => "The periodic cell is invalid.",
            Self::CellNotStandard => {
                "The cell is not in standard orientation. Cell vector a lies along +x, b in the \
                 xy plane with +y, c with +z."
            }
            Self::CellNoParameters => "The cell has no finite lengths and angles.",
            Self::CellDegenerateParameters => {
                "The cell is degenerate when written as lengths and angles."
            }
            Self::FrameChanged => "The frame changed between its first and second reading.",
            Self::DecodedFrameDiffers => "The decoded frame differs from its header.",
            Self::MagicOrPrecisionDiffers => {
                "The frame's magic number or precision differs from the first frame's."
            }
            Self::IndexedFrameGone => "An indexed frame is no longer in the file.",
            Self::RepeatedAtomCountDiffers => "The repeated atom count differs from the header's.",
            Self::XtcMagic => "The magic number is not 1995 or 2023.",
            Self::TrrMagic => "The magic number is not 1993.",
            Self::PrecisionNotInverseLength => "The XTC precision is not an inverse length.",
            Self::PrecisionNotPositive => "The precision is not a finite, positive number.",
            Self::CoordinateOutOfIntegerRange => {
                "A coordinate exceeds the integer range of the selected precision."
            }
            Self::SpanOutOfIntegerRange => {
                "The coordinates span more along one axis than the integer range of the selected \
                 precision holds."
            }
            Self::Xtc1995ByteCount => {
                "The compressed coordinates exceed the 32-bit byte count of XTC 1995. Use XTC \
                 2023."
            }
            Self::TrrVersionLength => "The version length is not that of GMX_trn_file.",
            Self::TrrVersionString => "The version string is not GMX_trn_file.",
            Self::VersionStringTooLong(length) => {
                return write!(
                    formatter,
                    "The version string has {length} bytes, more than 256 or the line limit."
                )
            }
            Self::UnsupportedTrrBlocks => {
                "The frame has input record, energy, virial, pressure, topology or symbol blocks. \
                 Kekule cannot keep them."
            }
            Self::NoPositions => "The frame has no positions.",
            Self::BlockSizesMismatch => {
                "The block sizes do not match the atom count and the scalar precision."
            }
            Self::EnergyTerms => "The frame has energy terms. Kekule cannot keep them.",
            Self::TimeOrLambdaNotFinite { time, lambda } => {
                return write!(
                    formatter,
                    "The time {time} or the lambda value {lambda} is not finite."
                )
            }
            Self::NoPrecisionBlock => "The frame has no block that shows its precision.",
            Self::NoScalarWidth => "The block sizes give no whole scalar width.",
            Self::ScalarWidth(width) => {
                return write!(formatter, "The scalar width is {width} bytes, not 4 or 8.")
            }
            Self::ScalarWidthDiffers => "The frame's scalar width differs from the first frame's.",
            Self::LambdaNotProperty => "The lambda value is not a valid frame property.",
            Self::LambdaNotZero => "The lambda value is not zero.",
            Self::LambdaPolicyKey => {
                "The lambda policy names a property other than gromacs.trr.lambda."
            }
            Self::LambdaPropertyValue => {
                "The lambda property is not one finite, dimensionless number."
            }
            Self::AtomCountLine => "The atom-count line is not one nonnegative whole number.",
            Self::ZeroAtoms => "The frame has zero atoms.",
            Self::NoCommentLine => "The frame has no comment line.",
            Self::EndsBeforeAtomLine(line) => {
                return write!(formatter, "The frame ends before atom line {line}.")
            }
            Self::AtomLineFields(line) => {
                return write!(
                    formatter,
                    "Atom line {line} lacks an element or one of three coordinates."
                )
            }
            Self::ExtraColumns(line) => {
                return write!(
                    formatter,
                    "Atom line {line} has columns after its coordinates."
                )
            }
            Self::ContentAfterBlankLine => "Content follows a blank line after the last frame.",
            Self::UnknownElement { line, symbol } => {
                return write!(
                    formatter,
                    "Atom line {line} names {symbol:?}. It is neither an element symbol nor an \
                     atomic number."
                )
            }
            Self::ElementMismatch {
                atom,
                found,
                expected,
            } => {
                return write!(
                    formatter,
                    "Atom {atom} is {}. The topology has {} there.",
                    found.symbol(),
                    expected.symbol()
                )
            }
            Self::CoordinateNotNumber(line) => {
                return write!(
                    formatter,
                    "Atom line {line} has a coordinate that is not a number."
                )
            }
            Self::CoordinateNotFiniteOnLine(line) => {
                return write!(
                    formatter,
                    "Atom line {line} has a coordinate that is not finite."
                )
            }
            Self::LineBufferGrowth(reason) => {
                return write!(formatter, "The line buffer cannot grow. {reason}")
            }
            Self::NoLineBreak => "The last line has no line break yet.",
            Self::NotFiniteInOutputUnit => "A coordinate is not finite in the output unit.",
        };
        formatter.write_str(text)
    }
}
