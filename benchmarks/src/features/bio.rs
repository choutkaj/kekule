use crate::*;

pub(super) fn mmcif_document_json(input: &Input) -> Result<Value, Box<dyn Error>> {
    let document = mmcif::parse_str(&input.text)?;
    let blocks = document
        .blocks()
        .iter()
        .map(|block| {
            let mut values = BTreeMap::new();
            for entry in block.entries() {
                match entry {
                    mmcif::MmcifEntry::Item(item) => {
                        if values
                            .insert(
                                item.tag().to_ascii_lowercase(),
                                vec![item.value().text().to_owned()],
                            )
                            .is_some()
                        {
                            return Err(boxed_error("duplicate mmCIF tag"));
                        }
                    }
                    mmcif::MmcifEntry::Loop(table) => {
                        for tag in table.tags() {
                            let column = (0..table.row_count())
                                .map(|row| {
                                    table
                                        .value(row, tag)
                                        .expect("parsed column")
                                        .text()
                                        .to_owned()
                                })
                                .collect::<Vec<_>>();
                            if values.insert(tag.to_ascii_lowercase(), column).is_some() {
                                return Err(boxed_error("duplicate mmCIF tag"));
                            }
                        }
                    }
                }
            }
            Ok(json!({"name":block.name(),"values":values}))
        })
        .collect::<Result<Vec<_>, Box<dyn Error>>>()?;
    Ok(json!({"blocks":blocks}))
}
pub(super) fn dssp_record_json(fixture_path: &Input) -> Result<Value, Box<dyn Error>> {
    let input = fixture_path.text.clone();
    let document = mmcif::parse_str(&input)?;
    let interpretation = mmcif::interpret(
        &document,
        MmcifInterpretOptions {
            model_selection: MmcifModelSelection::First,
            ..MmcifInterpretOptions::default()
        },
    )?;
    let result = match dssp::assign(interpretation.model().view(), dssp::DsspOptions::default()) {
        Ok(result) => result,
        Err(dssp::DsspError::NoAnalyzableProteinResidues) => {
            return Ok(json!({
                "status": "no_analyzable_residues",
                "residues": [],
            }));
        }
        Err(error) => return Err(Box::new(error)),
    };
    let residues = result.residues().collect::<Vec<_>>();
    let residues_by_key = residues
        .iter()
        .map(|residue| (residue.key(), *residue))
        .collect::<BTreeMap<_, _>>();
    let records = residues
        .iter()
        .map(|residue| {
            let source = residue.source();
            let sequence_id = source
                .author_sequence_id
                .as_deref()
                .and_then(|value| value.parse::<i32>().ok())
                .or(source.label_sequence_id);
            json!({
                "chain_id": source.chain_author_id.as_ref().unwrap_or(&source.chain_label_id),
                "sequence_id": sequence_id,
                "insertion_code": source.insertion_code,
                "label_chain_id": source.chain_label_id,
                "label_sequence_id": source.label_sequence_id,
                "residue_name": source.residue_name,
                "residue_one_letter": dssp_residue_letter(&source.residue_name),
                "secondary_structure": residue.secondary_structure().code().to_string(),
                "chain_break": dssp_chain_break_json(residue.chain_break()),
                "phi_degrees": residue.phi_degrees(),
                "psi_degrees": residue.psi_degrees(),
                "tco": residue.tco(),
                "kappa_degrees": residue.kappa_degrees(),
                "alpha_degrees": residue.alpha_degrees(),
                "helix_positions": residue.helix_positions().map(dssp_helix_position_json),
                "sheet": residue.sheet(),
                "strand": residue.strand(),
                "ladders": residue.beta_partners().map(|partner| partner.map(|partner| partner.ladder)),
                "beta_parallel": residue.beta_partners().map(|partner| partner.map(|partner| partner.parallel)),
                "acceptors": residue.acceptors().map(|bond| dssp_bond_json(bond, &residues_by_key)),
                "donors": residue.donors().map(|bond| dssp_bond_json(bond, &residues_by_key)),
            })
        })
        .collect::<Vec<_>>();
    Ok(json!({
        "status": "ok",
        "residues": records,
    }))
}

fn dssp_chain_break_json(chain_break: dssp::DsspChainBreak) -> &'static str {
    match chain_break {
        dssp::DsspChainBreak::None => "none",
        dssp::DsspChainBreak::NewChain => "new_chain",
        dssp::DsspChainBreak::Gap => "gap",
        _ => "unknown",
    }
}

fn dssp_helix_position_json(position: dssp::DsspHelixPosition) -> &'static str {
    match position {
        dssp::DsspHelixPosition::None => "none",
        dssp::DsspHelixPosition::Start => "start",
        dssp::DsspHelixPosition::End => "end",
        dssp::DsspHelixPosition::StartAndEnd => "start_and_end",
        dssp::DsspHelixPosition::Middle => "middle",
        _ => "unknown",
    }
}

fn dssp_bond_json(
    bond: Option<dssp::DsspHydrogenBond>,
    residues: &BTreeMap<kekule::topology::ResidueId, &dssp::DsspResidue>,
) -> Value {
    let Some(bond) = bond else {
        return Value::Null;
    };
    let source = residues[&bond.partner].source();
    json!({
        "partner_chain_id": source.chain_author_id.as_ref().unwrap_or(&source.chain_label_id),
        "partner_sequence_id": dssp_sequence_id(source),
        "partner_insertion_code": source.insertion_code,
        "energy_kcal_per_mol": bond.energy_kcal_per_mol,
    })
}

fn dssp_sequence_id(source: &dssp::DsspResidueSource) -> Option<i32> {
    source
        .author_sequence_id
        .as_deref()
        .and_then(|value| value.parse::<i32>().ok())
        .or(source.label_sequence_id)
}

fn dssp_residue_letter(name: &str) -> char {
    match name.to_ascii_uppercase().as_str() {
        "ALA" => 'A',
        "ARG" => 'R',
        "ASN" => 'N',
        "ASP" => 'D',
        "CYS" => 'C',
        "GLN" => 'Q',
        "GLU" => 'E',
        "GLY" => 'G',
        "HIS" => 'H',
        "ILE" => 'I',
        "LEU" => 'L',
        "LYS" => 'K',
        "MET" => 'M',
        "PHE" => 'F',
        "PRO" => 'P',
        "SER" => 'S',
        "THR" => 'T',
        "TRP" => 'W',
        "TYR" => 'Y',
        "VAL" => 'V',
        _ => 'X',
    }
}
