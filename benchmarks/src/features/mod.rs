use crate::{boxed_error, dataset::Input};
use serde_json::{json, Value};
use std::error::Error;
mod bio;
mod chemistry;
mod descriptors;
mod io;
mod stereo;
mod strict;
use descriptors::*;

pub(crate) fn is_writer(feature: &str) -> bool {
    feature.ends_with(".write") || matches!(feature, "io.smiles.canonical" | "io.smiles.isomeric")
}

pub(crate) fn evaluate(feature: &str, input: &Input) -> Result<Value, Box<dyn Error>> {
    if is_writer(feature) {
        return strict::write(feature, input);
    }
    if feature.starts_with("io.smiles.")
        || feature.starts_with("io.mol.")
        || feature.starts_with("io.sdf.")
        || matches!(feature, "stereo.representation" | "stereo.perception")
    {
        return strict::molecular(feature, input);
    }
    match feature {
        "io.mmcif.parse" => return bio::mmcif_document_json(input),
        "bio.secondary-structure.dssp" => return bio::dssp_record_json(input),
        "query.smarts" => return Ok(json!({"records":io::smarts_query_records_json(input)?})),
        "stereo.cip" => {
            let mut records = io::read_stereo_perception_records_by_suffix(input)?;
            let values = records
                .iter_mut()
                .map(stereo::stereo_cip_record_json)
                .collect::<Result<Vec<_>, _>>()?;
            return Ok(json!({"records":values}));
        }
        _ => (),
    }
    let mut records = io::read_small_records_by_suffix(input)?;
    let values = records
        .iter_mut()
        .map(|record| {
            Ok(match feature {
                "descriptor.molecular" => molecular_descriptor_record_json(record)?,
                "descriptor.rotatable-bonds.rdkit-strict" => rotatable_bond_record_json(record),
                "algo.rings.fast" => ring_membership_record_json(record),
                "algo.rings.sssr" => ring_set_record_json(record)?,
                "algo.valence.rdkit-like" => valence_record_json(record)?,
                "algo.aromaticity.rdkit-like" => aromaticity_record_json(record)?,
                "algo.canonical-ranking" => canonical_ranking_record_json(record)?,
                "algo.substructure.vf2" => io::substructure_record_json(record)?,
                "chem.perception.default" => default_perception_atom_record_json(record)?,
                "chem.hydrogen-transforms" => hydrogen_transform_record_json(record)?,
                _ => return Err(boxed_error(format!("no Kekule adapter for {feature}"))),
            })
        })
        .collect::<Result<Vec<_>, Box<dyn Error>>>()?;
    Ok(json!({"records":values}))
}
