use crate::*;

mod bio;
mod chemistry;
mod descriptors;
mod io;
mod smiles;
mod stereo;

use bio::{dssp_record_json, mmcif_document_json};
use descriptors::{
    aromaticity_record_json, canonical_ranking_record_json, default_perception_atom_record_json,
    hydrogen_transform_record_json, mol_parse_record_json, molecular_descriptor_record_json,
    ring_membership_record_json, ring_set_record_json, rotatable_bond_record_json,
    rotatable_bond_smiles_record_json, valence_record_json,
};
use io::{
    interpret_molfile, interpret_sdf, mol_record_json, read_small_records_by_suffix,
    read_stereo_perception_records_by_suffix, read_stereo_records_by_suffix, sdf_record_json,
    small_record, smarts_query_records_json, substructure_record_json,
};
use smiles::{
    isomeric_smiles_record_is_stereo_bearing, isomeric_smiles_record_json, smiles_parse_record_json,
};
use stereo::{stereo_cip_record_json, stereo_perception_group_record_json, stereo_record_json};

#[cfg(test)]
pub(crate) use chemistry::explicit_valence_json;
#[cfg(test)]
pub(crate) use io::{
    read_canonical_smiles_records, read_nonisomeric_smiles_records, read_smiles_records,
    smiles_unsupported_subset_reason, IndexedSmallRecord, IndexedSmilesRecord,
};
#[cfg(test)]
pub(crate) use smiles::{
    canonical_smiles_record_json, smiles_components_perceived_semantic_json,
    smiles_perceived_atoms_json, smiles_perceived_bonds_json, smiles_perceived_semantic_json,
    smiles_write_record_json,
};
#[cfg(test)]
pub(crate) use stereo::stereo_perception_record_json;

#[cfg(not(test))]
use io::{
    read_canonical_smiles_records, read_nonisomeric_smiles_records, read_smiles_records,
    IndexedSmallRecord,
};
#[cfg(not(test))]
use smiles::{canonical_smiles_record_json, smiles_write_record_json};

pub(crate) fn implementation_expected(
    benchmark: &str,
    _corpus: &str,
    fixture_path: &Path,
) -> Result<Value, Box<dyn Error>> {
    match benchmark {
        "bio.secondary-structure.dssp" => dssp_record_json(fixture_path),
        "io.mmcif.parse" => mmcif_document_json(fixture_path),
        "io.sdf.v2000.parse" => {
            let records = read_small_records_by_suffix(fixture_path)?;
            Ok(json!({ "records": records.iter().map(sdf_record_json).collect::<Vec<_>>() }))
        }
        "io.sdf.v2000.write" => {
            let records = read_small_records_by_suffix(fixture_path)?;
            let records = records
                .into_iter()
                .map(|record| {
                    let fields = record
                        .sdf_fields
                        .into_iter()
                        .map(|(name, value)| SdfDataField::new(name, value))
                        .collect();
                    SdfRecordInterpretation::new(record.title, vec![record.molecule], fields)
                })
                .collect::<Vec<_>>();
            let written = sdf::write_v2000(&records)?;
            let records = interpret_sdf(&written)?
                .into_iter()
                .enumerate()
                .map(|(index, record)| small_record(index, record))
                .collect::<Vec<_>>();
            Ok(json!({ "records": records.iter().map(sdf_record_json).collect::<Vec<_>>() }))
        }
        "io.mol.v2000.parse" => {
            let records = read_small_records_by_suffix(fixture_path)?;
            Ok(json!({ "records": records.iter().map(mol_parse_record_json).collect::<Vec<_>>() }))
        }
        "io.mol.v3000.parse" => {
            let records = read_small_records_by_suffix(fixture_path)?;
            let records = records
                .into_iter()
                .enumerate()
                .map(|(index, record)| {
                    let title = record.title;
                    let written = molfile::write_v3000(&record.molecule)?;
                    let molecule = interpret_molfile(&written)?;
                    Ok(IndexedSmallRecord {
                        record_index: index,
                        title,
                        molecule,
                        sdf_fields: BTreeMap::new(),
                    })
                })
                .collect::<Result<Vec<_>, Box<dyn Error>>>()?;
            Ok(json!({ "records": records.iter().map(mol_parse_record_json).collect::<Vec<_>>() }))
        }
        "descriptor.molecular" => {
            let mut records = read_small_records_by_suffix(fixture_path)?;
            Ok(json!({
                "records": records
                    .iter_mut()
                    .map(molecular_descriptor_record_json)
                    .collect::<Vec<_>>()
            }))
        }
        "descriptor.rotatable-bonds.rdkit-strict" => {
            if matches!(
                fixture_path
                    .extension()
                    .and_then(|extension| extension.to_str()),
                Some("smi" | "smiles" | "txt")
            ) {
                let records = read_smiles_records(fixture_path)?;
                return Ok(json!({
                    "records": records
                        .iter()
                        .map(rotatable_bond_smiles_record_json)
                        .collect::<Vec<_>>()
                }));
            }
            let records = read_small_records_by_suffix(fixture_path)?;
            Ok(json!({
                "records": records
                    .iter()
                    .map(rotatable_bond_record_json)
                    .collect::<Vec<_>>()
            }))
        }
        "io.mol.v2000.write" => {
            let records = read_small_records_by_suffix(fixture_path)?;
            let records = records
                .into_iter()
                .enumerate()
                .map(|(index, record)| {
                    let title = record.title;
                    let written = molfile::write_v2000(&record.molecule)?;
                    let molecule = interpret_molfile(&written)?;
                    Ok(IndexedSmallRecord {
                        record_index: index,
                        title,
                        molecule,
                        sdf_fields: BTreeMap::new(),
                    })
                })
                .collect::<Result<Vec<_>, Box<dyn Error>>>()?;
            Ok(json!({ "records": records.iter().map(mol_record_json).collect::<Vec<_>>() }))
        }
        "io.mol.v3000.write" => {
            let records = read_small_records_by_suffix(fixture_path)?;
            let records = records
                .into_iter()
                .enumerate()
                .map(|(index, record)| {
                    let title = record.title;
                    let written = molfile::write_v3000(&record.molecule)?;
                    let molecule = interpret_molfile(&written)?;
                    Ok(IndexedSmallRecord {
                        record_index: index,
                        title,
                        molecule,
                        sdf_fields: BTreeMap::new(),
                    })
                })
                .collect::<Result<Vec<_>, Box<dyn Error>>>()?;
            Ok(json!({ "records": records.iter().map(mol_record_json).collect::<Vec<_>>() }))
        }
        "io.smiles.parse" => {
            let records = read_nonisomeric_smiles_records(fixture_path)?;
            Ok(
                json!({ "records": records.iter().map(smiles_parse_record_json).collect::<Vec<_>>() }),
            )
        }
        "io.smiles.write" => {
            let records = read_nonisomeric_smiles_records(fixture_path)?;
            Ok(json!({
                "records": records
                    .iter()
                    .map(smiles_write_record_json)
                    .collect::<Result<Vec<_>, Box<dyn Error>>>()?
            }))
        }
        "io.smiles.canonical" => {
            let records = read_canonical_smiles_records(fixture_path)?;
            let exact_smiles = false;
            Ok(json!({
                "records": records
                    .iter()
                    .map(|record| canonical_smiles_record_json(record, exact_smiles))
                    .collect::<Result<Vec<_>, Box<dyn Error>>>()?
            }))
        }
        "io.smiles.isomeric" => {
            let records = read_canonical_smiles_records(fixture_path)?;
            let stereo_only = true;
            Ok(json!({
                "records": records
                    .iter()
                    .filter(|record| {
                        !stereo_only || isomeric_smiles_record_is_stereo_bearing(record)
                    })
                    .map(isomeric_smiles_record_json)
                    .collect::<Result<Vec<_>, Box<dyn Error>>>()?
            }))
        }
        "query.smarts" => Ok(json!({
            "records": smarts_query_records_json(fixture_path)?
        })),
        "algo.substructure.vf2" => {
            let mut records = read_small_records_by_suffix(fixture_path)?;
            Ok(json!({
                "records": records
                    .iter_mut()
                    .map(substructure_record_json)
                    .collect::<Vec<_>>()
            }))
        }
        "algo.rings.fast" => {
            let mut records = read_small_records_by_suffix(fixture_path)?;
            Ok(
                json!({ "records": records.iter_mut().map(ring_membership_record_json).collect::<Vec<_>>() }),
            )
        }
        "algo.rings.sssr" => {
            let mut records = read_small_records_by_suffix(fixture_path)?;
            Ok(
                json!({ "records": records.iter_mut().map(ring_set_record_json).collect::<Vec<_>>() }),
            )
        }
        "algo.valence.rdkit-like" => {
            let mut records = read_small_records_by_suffix(fixture_path)?;
            Ok(
                json!({ "records": records.iter_mut().map(valence_record_json).collect::<Vec<_>>() }),
            )
        }
        "chem.hydrogen-transforms" => {
            let mut records = read_small_records_by_suffix(fixture_path)?;
            Ok(json!({
                "records": records
                    .iter_mut()
                    .map(hydrogen_transform_record_json)
                    .collect::<Vec<_>>()
            }))
        }
        "chem.perception.default" => {
            let mut records = read_small_records_by_suffix(fixture_path)?;
            Ok(
                json!({ "records": records.iter_mut().map(default_perception_atom_record_json).collect::<Vec<_>>() }),
            )
        }
        "algo.aromaticity.rdkit-like" => {
            let mut records = read_small_records_by_suffix(fixture_path)?;
            Ok(
                json!({ "records": records.iter_mut().map(aromaticity_record_json).collect::<Vec<_>>() }),
            )
        }
        "algo.canonical-ranking" => {
            let mut records = read_small_records_by_suffix(fixture_path)?;
            Ok(json!({
                "records": records
                    .iter_mut()
                    .map(canonical_ranking_record_json)
                    .collect::<Vec<_>>()
            }))
        }
        "stereo.representation" => {
            let records = read_stereo_records_by_suffix(fixture_path)?;
            Ok(json!({ "records": records.iter().map(stereo_record_json).collect::<Vec<_>>() }))
        }
        "stereo.perception" => {
            let mut records = read_stereo_perception_records_by_suffix(fixture_path)?;
            Ok(json!({
                "records": records
                    .iter_mut()
                    .map(stereo_perception_group_record_json)
                    .collect::<Vec<_>>()
            }))
        }
        "stereo.cip" => {
            let mut records = read_stereo_records_by_suffix(fixture_path)?;
            let remove_plain_hydrogens = matches!(
                fixture_path.extension().and_then(|ext| ext.to_str()),
                Some("txt" | "smi" | "smiles")
            );
            Ok(json!({
                "records": records
                    .iter_mut()
                    .filter_map(|record| stereo_cip_record_json(record, remove_plain_hydrogens))
                    .collect::<Vec<_>>()
            }))
        }
        _ => Err(boxed_error(format!(
            "no implementation comparison configured for benchmark `{benchmark}`"
        ))),
    }
}
