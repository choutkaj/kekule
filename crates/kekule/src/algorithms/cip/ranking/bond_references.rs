//! Read-only CIP reference selection; atom IDs resolve only complete ties.

use super::*;

/// The caller supplies a live bond. Missing coordinate-bearing neighbors are
/// normal absence; implicit hydrogens and phantom atoms are never references.
pub(crate) fn bond_reference_atoms(
    mol: &Molecule,
    bond: BondId,
    options: CipAssignmentOptions,
) -> RankingResult<Option<[AtomId; 4]>> {
    let (mut b, mut c) = mol.bond(bond).expect("validated bond").endpoints();
    if c < b {
        std::mem::swap(&mut b, &mut c);
    }
    let neighbors = |root, other| {
        let mut atoms = mol
            .neighbors(root)
            .expect("live endpoint")
            .filter(|atom| *atom != other)
            .collect::<Vec<_>>();
        // Rule 1a at the attached atom can discard lower-priority candidates
        // without expanding their ligands or requiring inferred hydrogen counts.
        let number = |atom| {
            mol.atom(atom)
                .expect("live neighbor")
                .element
                .atomic_number()
        };
        if let Some(highest) = atoms.iter().map(|&atom| number(atom)).max() {
            atoms.retain(|&atom| number(atom) == highest);
        }
        atoms
    };
    let left = neighbors(b, c);
    let right = neighbors(c, b);
    if left.is_empty() || right.is_empty() {
        return Ok(None);
    }
    if left.len() == 1 && right.len() == 1 {
        return Ok((left[0] != right[0]).then_some([left[0], b, c, right[0]]));
    }
    if let Err(error) = validate_stereo(mol) {
        return Err(CipRankingError::InvalidStereo {
            issue: error.issues.into_iter().next().expect("validation issue"),
        });
    }
    // A single available reference needs no ranking. Otherwise ligand expansion
    // must know the implicit hydrogens: missing perception is not zero H.
    for atom in mol.atom_ids() {
        if mol.implicit_hydrogens(atom).expect("live atom").is_none() {
            return Err(CipRankingError::UnknownHydrogenCount { atom });
        }
    }
    let fractions = cip_atomic_number_fractions(mol);
    let a = reference(mol, b, &left, options, &fractions)?;
    let d = reference(mol, c, &right, options, &fractions)?;
    // A closed three-atom walk is not a four-atom bond dihedral. Do not select
    // a lower-priority reference merely to produce a value.
    Ok((a != d).then_some([a, b, c, d]))
}

fn reference(
    mol: &Molecule,
    root: AtomId,
    candidates: &[AtomId],
    options: CipAssignmentOptions,
    fractions: &[AtomicNumberFraction],
) -> RankingResult<AtomId> {
    if candidates.len() == 1 {
        return Ok(candidates[0]);
    }
    let descriptors = DescriptorContext::default();
    let context = LigandBuildContext {
        mol,
        descriptor_context: &descriptors,
        options,
        atomic_number_fractions: fractions,
        // Retain multiple-bond duplicates at the root of ordinary bond ligands.
        atropisomer_mode: true,
    };
    let mut signatures = signatures(&context, root, candidates)?;
    if signatures.len() > 1
        && mol
            .stereo_elements()
            .any(|(_, stereo)| stereo.is_specified())
    {
        let graph = build_auxiliary_graph(mol, root, options, fractions, true)?;
        let mut descriptors = DescriptorContext::default();
        precompute_auxiliary_descriptors(mol, &mut descriptors, &graph, options, fractions, true);
        // Descriptor assignment can deliberately leave an occurrence unlabelled.
        // Such incomplete stereochemical ranking must not become an ID tie.
        if collect_auxiliary_occurrences_from_molecule(mol, &descriptors, &graph)
            .iter()
            .any(|occurrence| {
                occurrence.distance > 0
                    && mol
                        .stereo_element(occurrence.key.element)
                        .expect("live stereo")
                        .is_specified()
                    && !descriptors.aux_labels.contains_key(&occurrence.key)
            })
        {
            return Err(CipRankingError::UnresolvedPriority);
        }
        let candidates = signatures
            .iter()
            .map(|(carrier, _)| match carrier {
                StereoCarrier::Atom(atom) => *atom,
                _ => unreachable!("bond references are explicit atoms"),
            })
            .collect::<Vec<_>>();
        signatures = self::signatures(
            &LigandBuildContext {
                descriptor_context: &descriptors,
                ..context
            },
            root,
            &candidates,
        )?;
    }
    // Expansion returns either a unique maximum or completed, tied maxima.
    // Preserve carrier identities when lower-priority candidates are removed.
    Ok(signatures
        .iter()
        .map(|(carrier, _)| match carrier {
            StereoCarrier::Atom(atom) => *atom,
            _ => unreachable!("bond references are explicit atoms"),
        })
        .min()
        .expect("nonempty reference candidates"))
}

fn signatures(
    context: &LigandBuildContext<'_>,
    root: AtomId,
    candidates: &[AtomId],
) -> RankingResult<Vec<(StereoCarrier, LigandSignature)>> {
    expansion::maximal_carrier_signatures(
        context,
        candidates.iter().map(|&atom| {
            let carrier = StereoCarrier::Atom(atom);
            (carrier, carrier_node(carrier, root))
        }),
        |node| {
            let priority = node.priority(context);
            let mut children = Vec::new();
            node.extend(
                context.mol,
                context.atomic_number_fractions,
                context.atropisomer_mode,
                &mut children,
            );
            (priority, children)
        },
    )
}
