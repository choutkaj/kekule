use super::{default_bond_expression, BondExpression, SmartsParseError};
use crate::core::{DoubleBondOrientation, TetrahedralOrientation};
use crate::query::{QueryAtomId, QueryBondId, QueryGraphBuilder, QueryStereoConstraint};
use std::ops::Range;

pub(super) struct BondSyntax {
    pub expression: BondExpression,
    /// Direction from the first written atom toward the second.
    pub direction: Option<bool>,
    pub span: Range<usize>,
}

impl BondSyntax {
    pub fn default_at(span: Range<usize>) -> Result<Self, SmartsParseError> {
        Ok(Self {
            expression: default_bond_expression()?,
            direction: None,
            span,
        })
    }

    pub fn ring(
        opening: Option<Self>,
        closing: Option<Self>,
        span: Range<usize>,
    ) -> Result<Self, SmartsParseError> {
        // Closing notation is written in the opposite direction to the edge.
        let closing = closing.map(|mut bond| {
            bond.direction = bond.direction.map(|up| !up);
            bond
        });
        match (opening, closing) {
            (None, None) => Self::default_at(span),
            (Some(bond), None) | (None, Some(bond)) => Ok(bond),
            (Some(mut left), Some(right)) => {
                if left.expression != right.expression
                    || matches!((left.direction, right.direction), (Some(a), Some(b)) if a != b)
                {
                    return Err(SmartsParseError::syntax(
                        left.span.start..right.span.end,
                        "conflicting bond expressions on ring closure",
                    ));
                }
                left.direction = left.direction.or(right.direction);
                Ok(left)
            }
        }
    }
}

pub(super) struct AtomStereo {
    orientation: TetrahedralOrientation,
    inline_hydrogen: bool,
    span: Range<usize>,
}

impl AtomStereo {
    pub fn parse(
        source: &str,
        start: usize,
        base: usize,
    ) -> Result<(Self, usize), SmartsParseError> {
        let bytes = source.as_bytes();
        let mut end = start + 1;
        let inverted = if bytes.get(end) == Some(&b'@') {
            end += 1;
            true
        } else if source[end..].starts_with("TH1") || source[end..].starts_with("TH2") {
            end += 3;
            bytes[end - 1] == b'2'
        } else {
            false
        };
        if bytes
            .get(end)
            .is_some_and(|byte| *byte == b'?' || byte.is_ascii_digit())
            || ["TH", "AL", "SP", "TB", "OH"]
                .iter()
                .any(|class| source[end..].starts_with(class))
        {
            return Err(SmartsParseError::unsupported(
                base + start..base + bytes.len(),
                "only specified tetrahedral stereochemical atom queries are supported",
            ));
        }
        // Only the H token immediately following chirality occupies a source
        // carrier slot. A separate ;H1 predicate merely tests hydrogen count.
        let inline_hydrogen =
            if bytes.get(end) == Some(&b'H') {
                let digits: String = source[end + 1..]
                    .chars()
                    .take_while(char::is_ascii_digit)
                    .collect();
                if digits.is_empty() {
                    true
                } else {
                    match digits.parse::<u8>() {
                        Ok(0) => false,
                        Ok(1) => true,
                        _ => return Err(SmartsParseError::unsupported(
                            base + start..base + end + 1 + digits.len(),
                            "tetrahedral inline hydrogen counts greater than one are unsupported",
                        )),
                    }
                }
            } else {
                false
            };
        Ok((
            Self {
                orientation: if inverted {
                    TetrahedralOrientation::CounterClockwise
                } else {
                    TetrahedralOrientation::Clockwise
                },
                inline_hydrogen,
                span: base + start..base + end,
            },
            end,
        ))
    }
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum SourceNeighbor {
    Atom(QueryAtomId),
    Ring(u16),
}

struct SourceAtom {
    neighbors: Vec<SourceNeighbor>,
    stereo: Option<AtomStereo>,
    root: bool,
}

struct DirectedEdge {
    a: QueryAtomId,
    b: QueryAtomId,
    up: bool,
    span: Range<usize>,
}

#[derive(Default)]
pub(super) struct StereoSyntax {
    atoms: Vec<SourceAtom>,
    doubles: Vec<(QueryBondId, QueryAtomId, QueryAtomId)>,
    directions: Vec<DirectedEdge>,
}

impl StereoSyntax {
    pub fn add_atom(
        &mut self,
        atom: QueryAtomId,
        previous: Option<QueryAtomId>,
        stereo: Option<AtomStereo>,
    ) {
        if let Some(previous) = previous {
            self.atoms[previous.index()]
                .neighbors
                .push(SourceNeighbor::Atom(atom));
        }
        self.atoms.push(SourceAtom {
            neighbors: previous.into_iter().map(SourceNeighbor::Atom).collect(),
            stereo,
            root: previous.is_none(),
        });
    }

    pub fn open_ring(&mut self, atom: QueryAtomId, label: u16) {
        self.atoms[atom.index()]
            .neighbors
            .push(SourceNeighbor::Ring(label));
    }

    pub fn close_ring(&mut self, opening: QueryAtomId, closing: QueryAtomId, label: u16) {
        let neighbor = self.atoms[opening.index()]
            .neighbors
            .iter_mut()
            .find(|neighbor| **neighbor == SourceNeighbor::Ring(label))
            .expect("open ring has a source-order placeholder");
        *neighbor = SourceNeighbor::Atom(closing);
        self.atoms[closing.index()]
            .neighbors
            .push(SourceNeighbor::Atom(opening));
    }

    pub fn add_bond(
        &mut self,
        id: QueryBondId,
        a: QueryAtomId,
        b: QueryAtomId,
        is_double: bool,
        direction: Option<bool>,
        span: Range<usize>,
    ) {
        if is_double {
            self.doubles.push((id, a, b));
        }
        if let Some(up) = direction {
            self.directions.push(DirectedEdge { a, b, up, span });
        }
    }

    pub fn install(self, builder: &mut QueryGraphBuilder) -> Result<(), SmartsParseError> {
        for (index, atom) in self.atoms.iter().enumerate() {
            let Some(stereo) = &atom.stereo else {
                continue;
            };
            let carriers: Vec<_> = atom
                .neighbors
                .iter()
                .map(|neighbor| match neighbor {
                    SourceNeighbor::Atom(atom) => *atom,
                    SourceNeighbor::Ring(_) => unreachable!("ring syntax validated before stereo"),
                })
                .collect();
            // The canonical query puts the omitted carrier last. At a root,
            // inline H is first in SMILES order: moving it past three neighbors
            // reverses parity. Elsewhere it follows the incoming neighbor.
            // With four graph neighbors, H is only a total-H predicate: it can
            // be satisfied by an explicit query hydrogen among those neighbors.
            let orientation = if stereo.inline_hydrogen && atom.root && carriers.len() == 3 {
                stereo.orientation.inverted()
            } else {
                stereo.orientation
            };
            builder
                .add_stereo_constraint(QueryStereoConstraint::Tetrahedral {
                    center: QueryAtomId::new(index as u32),
                    carriers,
                    orientation,
                })
                .map_err(|error| {
                    SmartsParseError::syntax(stereo.span.clone(), error.to_string())
                })?;
        }
        for (bond, a, b) in &self.doubles {
            let left = self.directed_neighbor(*a, *b)?;
            let right = self.directed_neighbor(*b, *a)?;
            if let (Some((left_carrier, left_up, span)), Some((right_carrier, right_up, _))) =
                (left, right)
            {
                builder
                    .add_stereo_constraint(QueryStereoConstraint::DoubleBond {
                        bond: *bond,
                        left_carrier,
                        right_carrier,
                        orientation: if left_up == right_up {
                            DoubleBondOrientation::Together
                        } else {
                            DoubleBondOrientation::Opposite
                        },
                    })
                    .map_err(|error| SmartsParseError::syntax(span, error.to_string()))?;
            }
        }
        // A lone direction carries no cis/trans relationship. As in RDKit,
        // it contributes only the single-or-aromatic predicate in that case.
        Ok(())
    }

    fn directed_neighbor(
        &self,
        endpoint: QueryAtomId,
        other: QueryAtomId,
    ) -> Result<Option<(QueryAtomId, bool, Range<usize>)>, SmartsParseError> {
        let mut selected: Option<(QueryAtomId, bool, Range<usize>)> = None;
        let mut seen = [false; 2];
        for edge in &self.directions {
            let (carrier, up) = if edge.a == endpoint {
                (edge.b, edge.up)
            } else if edge.b == endpoint {
                (edge.a, !edge.up)
            } else {
                continue;
            };
            if carrier == other {
                continue;
            }
            if let Some((_, _, span)) = &selected {
                if seen[usize::from(up)] {
                    return Err(SmartsParseError::syntax(
                        span.start.min(edge.span.start)..span.end.max(edge.span.end),
                        "contradictory directions at double-bond endpoint",
                    ));
                }
            } else {
                selected = Some((carrier, up, edge.span.clone()));
            }
            seen[usize::from(up)] = true;
        }
        Ok(selected)
    }
}
