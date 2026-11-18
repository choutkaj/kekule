# Frame-table format

A frame table holds a run's per-step records, such as energies,
temperature, volume, and bias, in a text file beside its trajectory. Rows
are keyed by MD step, so a trajectory's frames, usually a subset of the rows,
join it by step. The API lives in
[`kekule::traj::frames`](../crates/kekule/src/traj/frames/mod.rs).

The format must be appendable while a run grows, survive a crash, cut at a
checkpoint, and keep every bit of a real number. One row is one line, so all
four follow from complete lines. numpy and pandas read it with
`sep='\t', comment='#'` when no text cell holds `#`: that option cuts a line
at any `#`.

## Layout

UTF-8 text, LF line ends, tab-separated fields. The header comes first, and
the line naming the columns ends it. Readers also take CRLF line ends: they
strip one `\r` before each line feed, which no cell ends with, since the
writer escapes `\r` in text.

```text
#kekule-frames	1
#column	step	int
#column	time	real	ps
#column	stage	text
#column	potential_energy	real	kJ/mol
#column	bias.potential	real	kJ/mol
#column	bias.dihedral	real	kJ/mol
#time	time
#state	temperature	300	K
#state	pressure	1	bar
#sampling	from_step	2600000
#sampling	potential_energy	potential_energy
#sampling	method	GaMD dual boost
#bias	bias.potential
#bias	bias.dihedral
step	time	stage	potential_energy	bias.potential	bias.dihedral
500	1	nvt	-30521.4410981		
1000	2	nvt	-30530.12	1.25	0.5
```

## Header records

| Record | Meaning |
| --- | --- |
| `#kekule-frames 1` | The first line: the magic and the format version. Another version is unsupported. |
| `#column key type [unit]` | One column, in file order. The first is always `step int`. `type` is `int`, `real`, `text`, or `bool`; a `real` column names its unit. |
| `#time key` | The real column of each row's time. |
| `#state temperature value K` | The sampled temperature. |
| `#state pressure value bar` | The sampled pressure, in the isothermal-isobaric ensemble. |
| `#sampling from_step n` | The first MD step whose row is a sample of the state. |
| `#sampling potential_energy key` | The real column of the unbiased potential energy. |
| `#sampling method text` | A label of the method, such as `GaMD dual boost`. |
| `#bias key` | One bias column; each row's bias is the sum of these columns. |
| `##…` | A comment, ignored. Comments appear in the header only. |

The `#state`, `#sampling`, and `#bias` records are the run's sampling record
(`kekule::structure::Sampling`). A table with one has the temperature and
`from_step`; the potential energy and bias records name real columns of
energies. Any other record starting with `#` is unsupported, so a newer field
is never silently ignored. Records after the magic may come in any order;
`#column` and `#bias` records keep theirs. A first line that cannot become
`#kekule-frames` and a tab is no frame table, even before its line feed.

Units are kekule's unit-constant symbols, such as `ps`, `kJ/mol`, `K`, `bar`,
`nm^3`, and `1` for a pure number. Keys are kekule property keys. The line
naming the columns lists `step` and every `#column` key in order.

## Rows

One line per row, one field per column. Steps are non-negative decimal
integers and strictly increase.

| Type | Cell | Missing |
| --- | --- | --- |
| `int` | `i64` in decimal | empty |
| `real` | a finite number | empty |
| `bool` | `true` or `false` | empty |
| `text` | text with `\t`, `\n`, `\r`, and `\\` escaped | `\N` |

A real number is written in Rust's shortest round-trip form: plain for zero
and for magnitudes from 1e-4 to 1e16, in scientific notation otherwise
(`5e-324`, `1e300`). Reading it gives the same bits, `-0` included. Readers
accept Rust's `f64` grammar for finite numbers. The writer writes a value
that is not finite as missing, so a run never stops on a diverging energy.
An empty text is an empty field; a text `\N` is written `\\N`.

## Crashes and growing files

The writer writes and flushes one complete line per row. A final line
without a line feed, and NUL bytes after the last one (blocks lost at power
loss), are the pending tail: readers leave it out and report its length. A
complete line that does not parse is corrupt. A file whose header is not
complete yet is pending; a follower reads it again later.

A growing reader parses only the complete rows appended since its last read.
First it checks that the header and the last row it read are unchanged and
that the file is not shorter; otherwise, as when a run was cut and
continued, it fails and keeps the rows it had. It does not detect an edit in
place of an earlier row: writers only append, and a cut changes the last
row read or shortens the file.

## Continuing at a checkpoint

`append_frame_table` continues a table after a restart:

1. It parses the header and requires it to equal the run's: the same records
   with the same parsed values, in any order, whatever the comments.
   Otherwise the file is unchanged. A file without a complete header starts
   again.
2. It keeps the rows in order up to the first row after the checkpoint's
   step or the first line that is not a valid row, such as one a crash
   damaged, and cuts the rest and the pending tail. It syncs the cut to the
   disk. Keeping every row instead requires every complete line to be valid.
3. It returns a writer after the kept rows and reports the last kept step,
   which the run compares with its checkpoint: an earlier step means rows
   the checkpoint counts are missing.

With the same rows, the file then equals one written straight through, byte
for byte. Sync the file before a checkpoint counts its rows.

## Joining

`StructureSet::join_frames` puts a table's real and integer columns into a
set's columns by step: each geometry takes the row of its order's step, or of
its integer column `step`. The step and time columns and boolean and text
columns are not joined. A joined column replaces the set's column of the
same key, so joining a grown table again fills the new cells; a joined
column whose cells are all missing is absent from the set. The set takes the
table's sampling record and refuses a different one; it keeps its own when
the table has none.

More columns arrive as another table: `step` first, then its columns,
usually without a sampling record. A file beside a trajectory is named
`<trajectory stem>.frames.tsv`.
