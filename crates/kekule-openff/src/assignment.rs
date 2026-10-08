use crate::{
    error, explicit,
    offxml::{ForceField, Rule},
    parameters::*,
    ChargeAssignment, ChargeSource, NaglModel, Result,
};
use kekule::{
    core::{AromaticityModel, AtomId, Molecule},
    perception::aromaticity,
    substructure::{PreparedTarget, SubstructureMatchOptions, TaggedQuery},
    topology::{InstanceAtomId, Topology},
    units::{Quantity, ELEMENTARY_CHARGE},
};
use std::{
    collections::{BTreeMap, BTreeSet, VecDeque},
    sync::Arc,
};

pub(crate) fn options() -> SubstructureMatchOptions {
    SubstructureMatchOptions {
        max_matches: 1_000_000,
        max_search_states: 10_000_000,
        max_candidate_pairs: 16_000_000,
        uniquify: false,
        ..Default::default()
    }
}
fn ordered<const N: usize>(mut a: [AtomId; N]) -> [AtomId; N] {
    let mut reverse = a;
    reverse.reverse();
    if reverse < a {
        a = reverse;
    }
    a
}
fn assignments<const N: usize, P>(
    target: &PreparedTarget<'_>,
    rules: &[Rule<P>],
    improper: bool,
) -> Result<BTreeMap<[AtomId; N], usize>> {
    let mut map = BTreeMap::new();
    for (i, rule) in rules.iter().enumerate() {
        let query = TaggedQuery::new(&rule.query).map_err(error)?;
        for tuple in query.find_matches(target, options(), true).map_err(error)? {
            let mut atoms: [AtomId; N] = tuple
                .try_into()
                .map_err(|_| error("incorrect parameter tag count"))?;
            if improper {
                let mut outer = [atoms[0], atoms[2], atoms[3]];
                outer.sort();
                atoms[0] = outer[0];
                atoms[2] = outer[1];
                atoms[3] = outer[2];
            } else {
                atoms = ordered(atoms);
            }
            map.insert(atoms, i);
        }
    }
    Ok(map)
}
struct Assigned {
    molecule: Molecule,
    bonds: BTreeMap<[AtomId; 2], usize>,
    angles: BTreeMap<[AtomId; 3], usize>,
    propers: BTreeMap<[AtomId; 4], usize>,
    impropers: BTreeMap<[AtomId; 4], usize>,
    constraints: BTreeMap<[AtomId; 2], usize>,
    vdw: BTreeMap<[AtomId; 1], usize>,
}
impl ForceField {
    fn assign(&self, input: &Molecule) -> Result<Assigned> {
        let mut molecule = explicit(input)?;
        aromaticity::perceive_aromaticity(&mut molecule, AromaticityModel::Mdl).map_err(error)?;
        let target = PreparedTarget::new(&molecule);
        let assigned = Assigned {
            bonds: assignments(&target, &self.bonds, false)?,
            angles: assignments(&target, &self.angles, false)?,
            propers: assignments(&target, &self.propers, false)?,
            impropers: assignments(&target, &self.impropers, true)?,
            constraints: assignments(&target, &self.constraints, false)?,
            vdw: assignments(&target, &self.vdw, false)?,
            molecule,
        };
        let m = &assigned.molecule;
        let bonds = m
            .bonds()
            .map(|(_, b)| ordered([b.a(), b.b()]))
            .collect::<BTreeSet<_>>();
        let mut angles = BTreeSet::new();
        let mut propers = BTreeSet::new();
        for center in m.atom_ids() {
            let neighbors = m.neighbors(center).map_err(error)?.collect::<Vec<_>>();
            for i in 0..neighbors.len() {
                for j in i + 1..neighbors.len() {
                    angles.insert(ordered([neighbors[i], center, neighbors[j]]));
                }
            }
        }
        for (_, b) in m.bonds() {
            for a in m.neighbors(b.a()).map_err(error)?.filter(|a| *a != b.b()) {
                for d in m
                    .neighbors(b.b())
                    .map_err(error)?
                    .filter(|d| *d != b.a() && *d != a)
                {
                    propers.insert(ordered([a, b.a(), b.b(), d]));
                }
            }
        }
        complete("Bonds", &bonds, &assigned.bonds)?;
        complete("Angles", &angles, &assigned.angles)?;
        complete("ProperTorsions", &propers, &assigned.propers)?;
        complete("vdW", &m.atom_ids().map(|a| [a]).collect(), &assigned.vdw)?;
        for atoms in assigned.impropers.keys() {
            for i in [0, 2, 3] {
                if m.bond_between(atoms[1], atoms[i]).map_err(error)?.is_none() {
                    return Err(error("improper parameter has invalid connectivity"));
                }
            }
        }
        for (pair, &rule) in &assigned.constraints {
            if !bonds.contains(pair) && self.constraints[rule].parameter.1.is_none() {
                return Err(error(
                    "a nonbonded constraint requires an explicit distance",
                ));
            }
        }
        Ok(assigned)
    }
    /// Assign ordered SMIRKS rules using last-match-wins precedence, without loading NAGL.
    pub fn label_molecule(&self, molecule: &Molecule) -> Result<MoleculeLabels> {
        let a = self.assign(molecule)?;
        Ok(MoleculeLabels {
            bonds: a
                .bonds
                .iter()
                .map(|(k, i)| (*k, self.bonds[*i].parameter.source.clone()))
                .collect(),
            angles: a
                .angles
                .iter()
                .map(|(k, i)| (*k, self.angles[*i].parameter.source.clone()))
                .collect(),
            proper_torsions: a
                .propers
                .iter()
                .map(|(k, i)| (*k, self.propers[*i].parameter.source.clone()))
                .collect(),
            improper_torsions: a
                .impropers
                .iter()
                .map(|(k, i)| (*k, self.impropers[*i].parameter.source.clone()))
                .collect(),
            constraints: a
                .constraints
                .iter()
                .map(|(k, i)| (*k, self.constraints[*i].parameter.0.clone()))
                .collect(),
            vdw: a
                .vdw
                .iter()
                .map(|(k, i)| (*k, self.vdw[*i].parameter.source.clone()))
                .collect(),
        })
    }
    /// Parameterize each reusable definition once and instantiate its parameters.
    /// Errors publish no partial result. Hydrogens must already be explicit.
    pub fn parameterize(
        &self,
        topology: Arc<Topology>,
        model: &NaglModel,
    ) -> Result<ParameterizedTopology> {
        if self
            .charge_model()
            .is_some_and(|required| required != model.identity())
        {
            return Err(error(format!(
                "NAGL model identity mismatch: force field requires {:?}, supplied {:?}",
                self.charge_model(),
                model.identity()
            )));
        }
        self.parameterize_with(topology, |m| {
            self.charges(m, |m| {
                if self.charge_model().is_none() {
                    return Err(error(
                        "incomplete library charges and no NAGLCharges handler",
                    ));
                }
                model.assign_charges(m)
            })
        })
    }

    /// Parameterize using only LibraryCharges, without loading a neural model.
    ///
    /// Complete library coverage and the correct formal charge are required for
    /// every molecule. A declared NAGL handler is not used by this explicit path.
    /// Errors publish no partial result; the exact topology binding is retained.
    pub fn parameterize_without_nagl(
        &self,
        topology: Arc<Topology>,
    ) -> Result<ParameterizedTopology> {
        self.parameterize_with(topology, |m| {
            self.charges(m, |_| {
                Err(error(
                    "incomplete library charges; a compatible NAGL model is required",
                ))
            })
        })
    }
    fn parameterize_with(
        &self,
        topology: Arc<Topology>,
        mut charges_for: impl FnMut(&Molecule) -> Result<ChargeAssignment>,
    ) -> Result<ParameterizedTopology> {
        let mut definitions = BTreeMap::new();
        for (id, definition) in topology.definitions() {
            let assigned = self.assign(definition.molecule())?;
            let charges = charges_for(&assigned.molecule)?;
            definitions.insert(id, (assigned, charges));
        }
        let mut result = ParameterizedTopology {
            topology: topology.clone(),
            bonds: vec![],
            angles: vec![],
            propers: vec![],
            impropers: vec![],
            constraints: vec![],
            vdw: vec![],
            charges: Quantity::new(vec![], ELEMENTARY_CHARGE),
            charge_sources: vec![],
            exceptions: vec![],
            settings: self.settings.clone(),
        };
        // One shared allocation per rule (or per derived parameter variant).
        let bonds: Vec<Arc<BondParameter>> = shared(&self.bonds);
        let angles: Vec<Arc<AngleParameter>> = shared(&self.angles);
        let impropers: Vec<Arc<TorsionParameter>> = shared(&self.impropers);
        let vdw: Vec<Arc<VdwParameter>> = shared(&self.vdw);
        let mut propers = BTreeMap::<(usize, usize), Arc<TorsionParameter>>::new();
        let mut constraints = BTreeMap::<(usize, Option<usize>), Arc<ConstraintParameter>>::new();
        for instance in topology.molecules() {
            let (a, charges) = &definitions[&instance.definition_id()];
            let qualify = |atom| InstanceAtomId::new(instance.id(), atom);
            for (atoms, &i) in &a.bonds {
                result.bonds.push(Interaction {
                    atoms: atoms.map(qualify),
                    parameter: Arc::clone(&bonds[i]),
                });
            }
            for (atoms, &i) in &a.angles {
                result.angles.push(Interaction {
                    atoms: atoms.map(qualify),
                    parameter: Arc::clone(&angles[i]),
                });
            }
            let mut around = BTreeMap::<[AtomId; 2], usize>::new();
            for t in a.propers.keys() {
                *around.entry(ordered([t[1], t[2]])).or_default() += 1;
            }
            for (atoms, &i) in &a.propers {
                // Automatic idivf divides by the torsions around the central bond.
                let count = around[&ordered([atoms[1], atoms[2]])];
                let parameter = propers.entry((i, count)).or_insert_with(|| {
                    let mut p = self.propers[i].parameter.clone();
                    for term in &mut p.terms {
                        if term.idivf == 0.0 {
                            term.idivf = count as f64;
                        }
                    }
                    Arc::new(p)
                });
                result.propers.push(Interaction {
                    atoms: atoms.map(qualify),
                    parameter: Arc::clone(parameter),
                });
            }
            for (atoms, &i) in &a.impropers {
                for [x, y, z] in [[0, 2, 3], [2, 3, 0], [3, 0, 2]] {
                    result.impropers.push(Interaction {
                        atoms: [atoms[1], atoms[x], atoms[y], atoms[z]].map(qualify),
                        parameter: Arc::clone(&impropers[i]),
                    });
                }
            }
            for (atoms, &i) in &a.constraints {
                let (source, distance) = &self.constraints[i].parameter;
                let bond = distance.is_none().then(|| a.bonds[atoms]);
                let parameter = constraints.entry((i, bond)).or_insert_with(|| {
                    Arc::new(ConstraintParameter {
                        source: source.clone(),
                        distance: (*distance)
                            .unwrap_or_else(|| self.bonds[a.bonds[atoms]].parameter.length),
                    })
                });
                result.constraints.push(Interaction {
                    atoms: atoms.map(qualify),
                    parameter: Arc::clone(parameter),
                });
            }
            for atom in a.molecule.atom_ids() {
                result.vdw.push(Arc::clone(&vdw[a.vdw[&[atom]]]));
            }
            result
                .charges
                .value_mut()
                .extend(charges.charges.value().iter().copied());
            result.charge_sources.push(charges.source.clone());
            for start in a.molecule.atom_ids() {
                let mut queue = VecDeque::from([(start, 0)]);
                let mut visited = BTreeSet::from([start]);
                while let Some((at, depth)) = queue.pop_front() {
                    if depth == 3 {
                        continue;
                    }
                    for next in a.molecule.neighbors(at).map_err(error)? {
                        if visited.insert(next) {
                            queue.push_back((next, depth + 1));
                            if start < next {
                                result.exceptions.push(PairException {
                                    atoms: [qualify(start), qualify(next)],
                                    vdw_scale: self.settings.vdw_scales[depth],
                                    electrostatics_scale: self.settings.electrostatics_scales
                                        [depth],
                                });
                            }
                        }
                    }
                }
            }
        }
        Ok(result)
    }
    pub fn parameterize_molecule(
        &self,
        molecule: &Molecule,
        model: &NaglModel,
    ) -> Result<ParameterizedTopology> {
        self.parameterize(
            Arc::new(Topology::from_molecule(molecule).map_err(error)?),
            model,
        )
    }
    fn charges(
        &self,
        molecule: &Molecule,
        fallback: impl FnOnce(&Molecule) -> Result<ChargeAssignment>,
    ) -> Result<ChargeAssignment> {
        let target = PreparedTarget::new(molecule);
        let mut charges = BTreeMap::new();
        let mut sources = BTreeMap::new();
        for rule in &self.library {
            for tuple in TaggedQuery::new(&rule.query)
                .map_err(error)?
                .find_matches(&target, options(), true)
                .map_err(error)?
            {
                for (atom, q) in tuple.into_iter().zip(&rule.parameter.charges) {
                    charges.insert(atom, *q);
                    sources.insert(atom, rule.parameter.source.id.clone());
                }
            }
        }
        if charges.len() == molecule.atom_count() {
            let values = molecule.atom_ids().map(|a| charges[&a]).collect::<Vec<_>>();
            if (values.iter().sum::<f64>() - molecule.formal_charge() as f64).abs() > 1e-6 {
                return Err(error(
                    "library charges do not sum to the molecular formal charge",
                ));
            }
            return Ok(ChargeAssignment {
                charges: Quantity::new(values, ELEMENTARY_CHARGE),
                source: ChargeSource::Library {
                    parameter_ids: sources
                        .into_values()
                        .collect::<BTreeSet<_>>()
                        .into_iter()
                        .collect(),
                },
            });
        }
        fallback(molecule)
    }
}

fn shared<P: Clone>(rules: &[Rule<P>]) -> Vec<Arc<P>> {
    rules
        .iter()
        .map(|rule| Arc::new(rule.parameter.clone()))
        .collect()
}

fn complete<const N: usize>(
    handler: &str,
    expected: &BTreeSet<[AtomId; N]>,
    actual: &BTreeMap<[AtomId; N], usize>,
) -> Result<()> {
    for atoms in expected {
        if !actual.contains_key(atoms) {
            return Err(error(format!("unassigned {handler} interaction {atoms:?}")));
        }
    }
    for atoms in actual.keys() {
        if !expected.contains(atoms) {
            return Err(error(format!(
                "{handler} pattern matches invalid connectivity {atoms:?}"
            )));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use kekule::{smiles, topology::TopologyBuilder};

    #[test]
    fn fourteen_scaling_uses_shortest_bond_distance() {
        let mut molecule = smiles::to_molecules("CCCC").unwrap().remove(0);
        molecule.perceive().unwrap();
        molecule.add_hydrogens().unwrap();
        let carbons = molecule
            .atoms()
            .filter(|(_, a)| a.element.atomic_number() == 6)
            .map(|(a, _)| a)
            .collect::<Vec<_>>();
        let topology = Arc::new(Topology::from_molecule(&molecule).unwrap());
        let ff = ForceField::rosemary().unwrap();
        let p = ff
            .parameterize_with(topology, |m| {
                Ok(ChargeAssignment {
                    charges: Quantity::new(vec![0.0; m.atom_count()], ELEMENTARY_CHARGE),
                    source: ChargeSource::Library {
                        parameter_ids: vec![],
                    },
                })
            })
            .unwrap();
        let pair = p
            .pair_exceptions()
            .iter()
            .find(|p| p.atoms.map(|a| a.atom()) == [carbons[0], carbons[3]])
            .unwrap();
        assert_eq!(pair.vdw_scale, 0.5);
        assert_eq!(pair.electrostatics_scale, 0.8333333333);
    }

    #[test]
    fn reused_definitions_snapshot_trefoils_constraints_and_exceptions() {
        let molecule = smiles::to_molecules("[H]C([H])=O").unwrap().remove(0);
        let mut builder = TopologyBuilder::new();
        let definition = builder.add_molecule_definition(&molecule).unwrap();
        builder.add_instance(definition).unwrap();
        builder.add_instance(definition).unwrap();
        let topology = Arc::new(builder.build().unwrap());
        let ff = ForceField::rosemary().unwrap();
        let mut calls = 0;
        let p = ff
            .parameterize_with(topology.clone(), |m| {
                calls += 1;
                Ok(ChargeAssignment {
                    charges: Quantity::new(vec![0.0; m.atom_count()], ELEMENTARY_CHARGE),
                    source: ChargeSource::Library {
                        parameter_ids: vec![],
                    },
                })
            })
            .unwrap();
        assert_eq!(calls, 1);
        assert!(Arc::ptr_eq(p.topology(), &topology));
        assert_eq!(p.charges().value(), &[0.0; 8]);
        assert_eq!(p.bonds().len(), 6);
        assert_eq!(p.angles().len(), 6);
        assert_eq!(p.improper_torsions().len(), 6);
        for t in p.improper_torsions() {
            let center = topology.atom(t.atoms[0]).unwrap();
            assert_eq!(center.element.atomic_number(), 6);
            assert!(t.parameter.terms.iter().all(|term| term.idivf == 3.0));
        }
        assert_eq!(p.constraints().len(), 4);
        for c in p.constraints() {
            let b = p.bonds().iter().find(|b| b.atoms == c.atoms).unwrap();
            assert_eq!(c.parameter.distance, b.parameter.length);
        }
        // All six pairs within each four-atom star are at distance one or two;
        // instance boundaries must never introduce exceptions.
        assert_eq!(p.pair_exceptions().len(), 12);
        assert!(p
            .pair_exceptions()
            .iter()
            .all(|e| e.vdw_scale == 0.0 && e.electrostatics_scale == 0.0));
    }

    #[test]
    fn complete_library_charges_precede_inference_partial_coverage_falls_back() {
        let ff = ForceField::rosemary().unwrap();
        let ion = crate::explicit(&smiles::to_molecules("[Na+]").unwrap().remove(0)).unwrap();
        let charges = ff
            .charges(&ion, |_| panic!("library ion must not invoke NAGL"))
            .unwrap();
        assert_eq!(charges.charges.value(), &[1.0]);
        assert!(matches!(charges.source, ChargeSource::Library { .. }));
        let methanol = crate::explicit(
            &smiles::to_molecules("[H]C([H])([H])O[H]")
                .unwrap()
                .remove(0),
        )
        .unwrap();
        let mut modified = ff;
        let h = crate::offxml::LibraryCharge {
            source: ParameterIdentity {
                id: "partial".into(),
                smirks: "[#1:1]".into(),
            },
            charges: vec![0.0],
        };
        modified.library.push(Rule {
            query: kekule::query::parse_smarts(&h.source.smirks).unwrap(),
            parameter: h,
        });
        assert_eq!(
            modified
                .charges(&methanol, |_| Err(error("fallback reached")))
                .unwrap_err()
                .to_string(),
            "fallback reached"
        );
    }
}
