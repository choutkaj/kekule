use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

use crate::bio::*;
use crate::core::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MmcifParseOptions {
    pub strict: bool,
    pub max_input_bytes: usize,
    pub max_tokens: usize,
    pub max_token_bytes: usize,
    pub max_atom_site_rows: usize,
}

impl Default for MmcifParseOptions {
    fn default() -> Self {
        Self {
            strict: true,
            max_input_bytes: 256 * 1024 * 1024,
            max_tokens: 10_000_000,
            max_token_bytes: 16 * 1024 * 1024,
            max_atom_site_rows: 5_000_000,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MmcifParseError {
    pub line: usize,
    pub message: String,
}

impl MmcifParseError {
    fn new(line: usize, message: impl Into<String>) -> Self {
        Self {
            line,
            message: message.into(),
        }
    }
}

impl fmt::Display for MmcifParseError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "mmCIF parse error at line {}: {}",
            self.line, self.message
        )
    }
}

impl std::error::Error for MmcifParseError {}

pub fn read_mmcif_str(
    input: &str,
    options: MmcifParseOptions,
) -> std::result::Result<MacroMolecule, MmcifParseError> {
    let tokens = tokenize_mmcif(input, options)?;
    let atom_site_loop = find_atom_site_loop(&tokens)
        .ok_or_else(|| MmcifParseError::new(1, "missing _atom_site loop"))?;
    build_macro_molecule_from_atom_site_loop(atom_site_loop, options)
}

#[derive(Debug, Clone)]
struct MmcifToken {
    text: String,
    line: usize,
}

#[derive(Debug, Clone)]
struct MmcifLoop<'a> {
    tags: Vec<&'a MmcifToken>,
    values: Vec<&'a MmcifToken>,
}

fn tokenize_mmcif(
    input: &str,
    options: MmcifParseOptions,
) -> std::result::Result<Vec<MmcifToken>, MmcifParseError> {
    if input.len() > options.max_input_bytes {
        return Err(MmcifParseError::new(
            1,
            "input exceeds configured byte limit",
        ));
    }
    let normalized = input.replace("\r\n", "\n").replace('\r', "\n");
    let lines = normalized.lines().collect::<Vec<_>>();
    let mut tokens = Vec::new();
    let mut line_index = 0usize;

    while line_index < lines.len() {
        let line_number = line_index + 1;
        let line = lines[line_index];
        if line.starts_with(';') {
            let mut text = String::new();
            line_index += 1;
            while line_index < lines.len() && !lines[line_index].starts_with(';') {
                if !text.is_empty() {
                    text.push('\n');
                }
                let next_len = text
                    .len()
                    .checked_add(lines[line_index].len())
                    .ok_or_else(|| MmcifParseError::new(line_number, "token length overflow"))?;
                if next_len > options.max_token_bytes {
                    return Err(MmcifParseError::new(
                        line_number,
                        "semicolon value exceeds configured token limit",
                    ));
                }
                text.push_str(lines[line_index]);
                line_index += 1;
            }
            if line_index == lines.len() {
                return Err(MmcifParseError::new(
                    line_number,
                    "unterminated semicolon text",
                ));
            }
            push_mmcif_token(&mut tokens, text, line_number, options)?;
            line_index += 1;
            continue;
        }

        let bytes = line.as_bytes();
        let mut column = 0usize;
        while column < bytes.len() {
            while column < bytes.len() && bytes[column].is_ascii_whitespace() {
                column += 1;
            }
            if column == bytes.len() || bytes[column] == b'#' {
                break;
            }
            let start = column;
            let text = if bytes[column] == b'\'' || bytes[column] == b'"' {
                let quote = bytes[column];
                column += 1;
                let value_start = column;
                while column < bytes.len() && bytes[column] != quote {
                    column += 1;
                }
                if column == bytes.len() {
                    return Err(MmcifParseError::new(
                        line_number,
                        "unterminated quoted value",
                    ));
                }
                let value = &line[value_start..column];
                column += 1;
                value.to_owned()
            } else {
                while column < bytes.len()
                    && !bytes[column].is_ascii_whitespace()
                    && bytes[column] != b'#'
                {
                    column += 1;
                }
                line[start..column].to_owned()
            };
            push_mmcif_token(&mut tokens, text, line_number, options)?;
        }
        line_index += 1;
    }

    Ok(tokens)
}

fn push_mmcif_token(
    tokens: &mut Vec<MmcifToken>,
    text: String,
    line: usize,
    options: MmcifParseOptions,
) -> std::result::Result<(), MmcifParseError> {
    if text.len() > options.max_token_bytes {
        return Err(MmcifParseError::new(
            line,
            "value exceeds configured token limit",
        ));
    }
    if tokens.len() >= options.max_tokens {
        return Err(MmcifParseError::new(
            line,
            "input exceeds configured token count",
        ));
    }
    tokens.push(MmcifToken { text, line });
    Ok(())
}

fn find_atom_site_loop(tokens: &[MmcifToken]) -> Option<MmcifLoop<'_>> {
    let mut index = 0usize;
    while index < tokens.len() {
        if tokens[index].text != "loop_" {
            index += 1;
            continue;
        }
        index += 1;
        let tag_start = index;
        while index < tokens.len() && tokens[index].text.starts_with('_') {
            index += 1;
        }
        let tags = tokens[tag_start..index].iter().collect::<Vec<_>>();
        let value_start = index;
        while index < tokens.len()
            && tokens[index].text != "loop_"
            && !tokens[index].text.starts_with("data_")
            && !tokens[index].text.starts_with('_')
        {
            index += 1;
        }
        if tags.iter().any(|tag| tag.text.starts_with("_atom_site.")) {
            return Some(MmcifLoop {
                tags,
                values: tokens[value_start..index].iter().collect(),
            });
        }
    }
    None
}

fn build_macro_molecule_from_atom_site_loop(
    atom_loop: MmcifLoop<'_>,
    options: MmcifParseOptions,
) -> std::result::Result<MacroMolecule, MmcifParseError> {
    let width = atom_loop.tags.len();
    if width == 0 || atom_loop.values.len() % width != 0 {
        let line = atom_loop.values.last().map(|token| token.line).unwrap_or(1);
        return Err(MmcifParseError::new(line, "atom-site loop has ragged rows"));
    }
    let row_count = atom_loop.values.len() / width;
    if row_count > options.max_atom_site_rows {
        let line = atom_loop
            .values
            .first()
            .map(|token| token.line)
            .unwrap_or(1);
        return Err(MmcifParseError::new(
            line,
            "atom-site loop exceeds configured row limit",
        ));
    }

    let tag_index = atom_loop.tags.iter().enumerate().try_fold(
        BTreeMap::new(),
        |mut tags, (index, token)| {
            if tags.insert(token.text.as_str(), index).is_some() {
                return Err(MmcifParseError::new(
                    token.line,
                    format!("duplicate atom-site loop tag {}", token.text),
                ));
            }
            Ok(tags)
        },
    )?;
    let mut macro_mol = MacroMolecule::default();
    let mut models = BTreeMap::<String, ModelId>::new();
    let mut chains = BTreeMap::<(String, String), ChainId>::new();
    let mut residues = BTreeMap::<MmcifResidueKey, ResidueId>::new();
    let mut ambiguous_residue = None::<MmcifAmbiguousResidueState>;
    let mut next_ambiguous_occurrence = 0usize;
    let mut conformer = Conformer::new();
    let mut saw_coordinates = false;

    for (row_index, row) in atom_loop.values.chunks(width).enumerate() {
        let line = row.first().map(|token| token.line).unwrap_or(1);
        let type_symbol = required_mmcif_value(row, &tag_index, "_atom_site.type_symbol", line)?;
        let element_symbol = canonical_mmcif_element_symbol(type_symbol);
        let element = Element::from_symbol(&element_symbol).ok_or_else(|| {
            MmcifParseError::new(line, format!("unknown atom-site element `{type_symbol}`"))
        })?;
        let group_pdb = optional_mmcif_value(row, &tag_index, "_atom_site.group_PDB");
        let atom_site_id = optional_mmcif_value(row, &tag_index, "_atom_site.id");
        let label_atom_id = optional_mmcif_value(row, &tag_index, "_atom_site.label_atom_id");
        let auth_atom_id = optional_mmcif_value(row, &tag_index, "_atom_site.auth_atom_id");
        let label_asym_id = optional_mmcif_value(row, &tag_index, "_atom_site.label_asym_id");
        let auth_asym_id = optional_mmcif_value(row, &tag_index, "_atom_site.auth_asym_id");
        let label_comp_id = optional_mmcif_value(row, &tag_index, "_atom_site.label_comp_id");
        let auth_comp_id = optional_mmcif_value(row, &tag_index, "_atom_site.auth_comp_id");
        let label_seq_id =
            optional_i32_mmcif_value(row, &tag_index, "_atom_site.label_seq_id", line)?;
        let auth_seq_id = optional_mmcif_value(row, &tag_index, "_atom_site.auth_seq_id");
        let insertion_code = optional_mmcif_value(row, &tag_index, "_atom_site.pdbx_PDB_ins_code");
        let model_key = optional_mmcif_value(row, &tag_index, "_atom_site.pdbx_PDB_model_num")
            .unwrap_or("1")
            .to_owned();
        let label_alt_id = optional_mmcif_value(row, &tag_index, "_atom_site.label_alt_id");
        let occupancy_raw = optional_mmcif_value(row, &tag_index, "_atom_site.occupancy");
        let b_factor_raw = optional_mmcif_value(row, &tag_index, "_atom_site.B_iso_or_equiv");
        let cartn_x_raw = optional_mmcif_value(row, &tag_index, "_atom_site.Cartn_x");
        let cartn_y_raw = optional_mmcif_value(row, &tag_index, "_atom_site.Cartn_y");
        let cartn_z_raw = optional_mmcif_value(row, &tag_index, "_atom_site.Cartn_z");
        let occupancy = optional_f64_mmcif_value(row, &tag_index, "_atom_site.occupancy", line)?;
        let b_factor =
            optional_f64_mmcif_value(row, &tag_index, "_atom_site.B_iso_or_equiv", line)?;
        let coordinates = optional_mmcif_point(row, &tag_index, line, options.strict)?;

        let chain_label = label_asym_id
            .or(auth_asym_id)
            .ok_or_else(|| MmcifParseError::new(line, "missing atom-site chain identifier"))?;
        let hierarchy_chain_id = auth_asym_id
            .or(label_asym_id)
            .ok_or_else(|| MmcifParseError::new(line, "missing atom-site chain identifier"))?;
        let residue_name = label_comp_id
            .or(auth_comp_id)
            .ok_or_else(|| MmcifParseError::new(line, "missing atom-site residue name"))?;
        if options.strict && label_atom_id.is_none() {
            return Err(MmcifParseError::new(
                line,
                "missing atom-site label atom id",
            ));
        }

        let model = *models
            .entry(model_key.clone())
            .or_insert_with(|| macro_mol.hierarchy.add_model(model_key.clone()));
        let chain_key = (model_key.clone(), hierarchy_chain_id.to_owned());
        let chain = if let Some(chain) = chains.get(&chain_key) {
            *chain
        } else {
            let chain = macro_mol
                .hierarchy
                .add_chain(
                    model,
                    chain_label.to_owned(),
                    auth_asym_id.map(str::to_owned),
                )
                .map_err(|error| MmcifParseError::new(line, error.to_string()))?;
            chains.insert(chain_key, chain);
            chain
        };
        let ambiguous_occurrence = if label_seq_id.is_none() && auth_seq_id.is_none() {
            Some(next_mmcif_ambiguous_occurrence(
                &mut ambiguous_residue,
                &mut next_ambiguous_occurrence,
                model_key.as_str(),
                chain_label,
                residue_name,
                insertion_code,
                label_atom_id.or(auth_atom_id),
                label_alt_id,
            ))
        } else {
            ambiguous_residue = None;
            None
        };
        let residue_key = MmcifResidueKey::from_row(
            line,
            options.strict,
            model_key.as_str(),
            chain_label,
            auth_asym_id,
            residue_name,
            auth_comp_id,
            label_seq_id,
            auth_seq_id,
            insertion_code,
            ambiguous_occurrence.unwrap_or(row_index),
        )?;
        let residue = if let Some(residue) = residues.get(&residue_key) {
            *residue
        } else {
            let residue = macro_mol
                .hierarchy
                .add_residue(
                    chain,
                    residue_name.to_owned(),
                    label_seq_id,
                    auth_seq_id.map(str::to_owned),
                    insertion_code.map(str::to_owned),
                )
                .map_err(|error| MmcifParseError::new(line, error.to_string()))?;
            let residue_record = &mut macro_mol.hierarchy.residues[residue.index()];
            residue_record.label_comp_id = label_comp_id.map(str::to_owned);
            residue_record.author_comp_id = auth_comp_id.map(str::to_owned);
            residues.insert(residue_key, residue);
            residue
        };

        let atom = macro_mol.mol.add_atom(Atom::new(element));
        if let Some(point) = coordinates {
            conformer.set_position(atom, point);
            saw_coordinates = true;
        }
        macro_mol
            .add_atom_site(
                residue,
                atom,
                AtomSiteMetadata {
                    group_pdb: group_pdb.map(str::to_owned),
                    atom_site_id: atom_site_id.map(str::to_owned),
                    type_symbol: Some(type_symbol.to_owned()),
                    label_asym_id: label_asym_id.map(str::to_owned),
                    auth_asym_id: auth_asym_id.map(str::to_owned),
                    label_atom_id: label_atom_id.map(str::to_owned),
                    auth_atom_id: auth_atom_id.map(str::to_owned),
                    label_alt_id: label_alt_id.map(str::to_owned),
                    occupancy,
                    occupancy_raw: occupancy_raw.map(str::to_owned),
                    b_factor,
                    b_factor_raw: b_factor_raw.map(str::to_owned),
                    cartn_x_raw: cartn_x_raw.map(str::to_owned),
                    cartn_y_raw: cartn_y_raw.map(str::to_owned),
                    cartn_z_raw: cartn_z_raw.map(str::to_owned),
                },
            )
            .map_err(|error| MmcifParseError::new(line, error.to_string()))?;
    }
    if saw_coordinates {
        macro_mol.mol.add_conformer(conformer);
    }

    Ok(macro_mol)
}

#[derive(Debug, Clone, PartialEq, Eq)]
struct MmcifAmbiguousResidueState {
    model_id: String,
    chain_id: String,
    component_id: String,
    insertion_code: Option<String>,
    occurrence: usize,
    atom_sites: BTreeSet<(String, Option<String>)>,
}

#[allow(clippy::too_many_arguments)]
fn next_mmcif_ambiguous_occurrence(
    current: &mut Option<MmcifAmbiguousResidueState>,
    next_occurrence: &mut usize,
    model_id: &str,
    chain_id: &str,
    component_id: &str,
    insertion_code: Option<&str>,
    atom_id: Option<&str>,
    alt_id: Option<&str>,
) -> usize {
    let insertion_code = insertion_code.map(str::to_owned);
    let atom_site = atom_id.map(|atom_id| (atom_id.to_owned(), alt_id.map(str::to_owned)));
    let same_context = current.as_ref().is_some_and(|state| {
        state.model_id == model_id
            && state.chain_id == chain_id
            && state.component_id == component_id
            && state.insertion_code == insertion_code
    });
    let repeats_atom_site = match (current.as_ref(), atom_site.as_ref()) {
        (_, None) => true,
        (Some(state), Some((atom_id, alt_id))) => {
            let prior_alt_ids = state
                .atom_sites
                .iter()
                .filter_map(|(prior_atom_id, prior_alt_id)| {
                    (prior_atom_id == atom_id).then_some(prior_alt_id)
                })
                .collect::<Vec<_>>();
            !prior_alt_ids.is_empty()
                && (alt_id.is_none()
                    || prior_alt_ids
                        .iter()
                        .any(|prior_alt_id| prior_alt_id.is_none() || *prior_alt_id == alt_id))
        }
        (None, Some(_)) => false,
    };
    if !same_context || repeats_atom_site {
        let occurrence = *next_occurrence;
        *next_occurrence = next_occurrence.saturating_add(1);
        *current = Some(MmcifAmbiguousResidueState {
            model_id: model_id.to_owned(),
            chain_id: chain_id.to_owned(),
            component_id: component_id.to_owned(),
            insertion_code,
            occurrence,
            atom_sites: BTreeSet::new(),
        });
    }
    if let (Some(state), Some(atom_site)) = (current.as_mut(), atom_site) {
        state.atom_sites.insert(atom_site);
    }
    current
        .as_ref()
        .map(|state| state.occurrence)
        .unwrap_or_default()
}

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
enum MmcifResidueKey {
    Label {
        model_id: String,
        label_chain_id: String,
        component_id: String,
        label_seq_id: i32,
        insertion_code: Option<String>,
    },
    Author {
        model_id: String,
        chain_id: String,
        component_id: String,
        auth_seq_id: String,
        insertion_code: Option<String>,
    },
    Occurrence {
        model_id: String,
        chain_id: String,
        component_id: String,
        insertion_code: Option<String>,
        occurrence: usize,
    },
}

impl MmcifResidueKey {
    #[allow(clippy::too_many_arguments)]
    fn from_row(
        line: usize,
        strict: bool,
        model_id: &str,
        chain_label: &str,
        auth_asym_id: Option<&str>,
        residue_name: &str,
        auth_comp_id: Option<&str>,
        label_seq_id: Option<i32>,
        auth_seq_id: Option<&str>,
        insertion_code: Option<&str>,
        occurrence: usize,
    ) -> std::result::Result<Self, MmcifParseError> {
        let insertion_code = insertion_code.map(str::to_owned);
        if let Some(label_seq_id) = label_seq_id {
            return Ok(Self::Label {
                model_id: model_id.to_owned(),
                label_chain_id: chain_label.to_owned(),
                component_id: residue_name.to_owned(),
                label_seq_id,
                insertion_code,
            });
        }
        if let Some(auth_seq_id) = auth_seq_id {
            return Ok(Self::Author {
                model_id: model_id.to_owned(),
                chain_id: auth_asym_id.unwrap_or(chain_label).to_owned(),
                component_id: auth_comp_id.unwrap_or(residue_name).to_owned(),
                auth_seq_id: auth_seq_id.to_owned(),
                insertion_code,
            });
        }
        if strict {
            return Err(MmcifParseError::new(
                line,
                "missing residue sequence identifier",
            ));
        }
        Ok(Self::Occurrence {
            model_id: model_id.to_owned(),
            chain_id: chain_label.to_owned(),
            component_id: residue_name.to_owned(),
            insertion_code,
            occurrence,
        })
    }
}

fn required_mmcif_value<'a>(
    row: &'a [&'a MmcifToken],
    tag_index: &BTreeMap<&str, usize>,
    tag: &str,
    line: usize,
) -> std::result::Result<&'a str, MmcifParseError> {
    optional_mmcif_value(row, tag_index, tag)
        .ok_or_else(|| MmcifParseError::new(line, format!("missing required {tag}")))
}

fn optional_mmcif_value<'a>(
    row: &'a [&'a MmcifToken],
    tag_index: &BTreeMap<&str, usize>,
    tag: &str,
) -> Option<&'a str> {
    let value = row.get(*tag_index.get(tag)?)?.text.as_str();
    (!matches!(value, "." | "?")).then_some(value)
}

fn canonical_mmcif_element_symbol(symbol: &str) -> String {
    let mut chars = symbol.chars();
    let Some(first) = chars.next() else {
        return String::new();
    };
    let mut canonical = first.to_ascii_uppercase().to_string();
    canonical.extend(chars.flat_map(char::to_lowercase));
    canonical
}

fn optional_i32_mmcif_value(
    row: &[&MmcifToken],
    tag_index: &BTreeMap<&str, usize>,
    tag: &str,
    line: usize,
) -> std::result::Result<Option<i32>, MmcifParseError> {
    optional_mmcif_value(row, tag_index, tag)
        .map(|value| {
            value
                .parse()
                .map_err(|_| MmcifParseError::new(line, format!("invalid integer {tag}")))
        })
        .transpose()
}

fn optional_f64_mmcif_value(
    row: &[&MmcifToken],
    tag_index: &BTreeMap<&str, usize>,
    tag: &str,
    line: usize,
) -> std::result::Result<Option<f64>, MmcifParseError> {
    optional_mmcif_value(row, tag_index, tag)
        .map(|value| {
            let parsed = value
                .parse()
                .map_err(|_| MmcifParseError::new(line, format!("invalid float {tag}")))?;
            if !f64::is_finite(parsed) {
                return Err(MmcifParseError::new(
                    line,
                    format!("non-finite float {tag}"),
                ));
            }
            Ok(parsed)
        })
        .transpose()
}

fn optional_mmcif_point(
    row: &[&MmcifToken],
    tag_index: &BTreeMap<&str, usize>,
    line: usize,
    strict: bool,
) -> std::result::Result<Option<Point3>, MmcifParseError> {
    let x = optional_f64_mmcif_value(row, tag_index, "_atom_site.Cartn_x", line)?;
    let y = optional_f64_mmcif_value(row, tag_index, "_atom_site.Cartn_y", line)?;
    let z = optional_f64_mmcif_value(row, tag_index, "_atom_site.Cartn_z", line)?;
    match (x, y, z) {
        (Some(x), Some(y), Some(z)) => Ok(Some(Point3::new(x, y, z))),
        (None, None, None) => Ok(None),
        _ if strict => Err(MmcifParseError::new(
            line,
            "partial atom-site coordinate triplet",
        )),
        _ => Ok(None),
    }
}
