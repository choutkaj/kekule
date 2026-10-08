use crate::{error, parameters::*, ModelIdentity, Result};
use kekule::{
    query::{parse_smarts, QueryGraph},
    substructure::TaggedQuery,
    units::*,
};
use roxmltree::Node;
use std::{collections::BTreeSet, path::Path};

mod quantity;
pub(crate) use quantity::parse_quantity;
#[cfg(test)]
mod reference_tests;

#[derive(Debug, Clone)]
pub(crate) struct Rule<P> {
    pub query: QueryGraph,
    pub parameter: P,
}
#[derive(Debug, Clone)]
pub(crate) struct LibraryCharge {
    pub source: ParameterIdentity,
    pub charges: Vec<f64>,
}
/// A compiled SMIRNOFF force field supporting the functional forms in Rosemary.
///
/// Unsupported handlers, units, interpolation, or potentials fail explicitly.
/// The packaged Rosemary release is `openff_no_water-3.0.0-alpha2b` (2026-09-11).
/// It has no water model or virtual sites. NAGL weights are loaded separately.
#[derive(Debug, Clone)]
pub struct ForceField {
    pub(crate) bonds: Vec<Rule<BondParameter>>,
    pub(crate) angles: Vec<Rule<AngleParameter>>,
    pub(crate) propers: Vec<Rule<TorsionParameter>>,
    pub(crate) impropers: Vec<Rule<TorsionParameter>>,
    pub(crate) constraints: Vec<Rule<(ParameterIdentity, Option<Quantity<f64>>)>>,
    pub(crate) vdw: Vec<Rule<VdwParameter>>,
    pub(crate) library: Vec<Rule<LibraryCharge>>,
    pub(crate) settings: NonbondedSettings,
    pub(crate) charge_model: Option<ModelIdentity>,
}
impl ForceField {
    /// Optional model identity declared by the NAGLCharges handler.
    pub fn charge_model(&self) -> Option<&ModelIdentity> {
        self.charge_model.as_ref()
    }

    /// Append compiled rules with last-match-wins precedence, preserving their IDs.
    ///
    /// Both fields must use compatible nonbonded settings. Only length conversion
    /// roundoff (8 machine epsilons relative) is allowed; the first values remain.
    /// Two declared NAGL models must match; an absent declaration does not erase
    /// an existing one. Failure leaves this force field unchanged. This composes
    /// complete supported documents, not arbitrary partial handler fragments.
    pub fn append(&mut self, other: &Self) -> Result<()> {
        if !compatible_settings(&self.settings, &other.settings) {
            return Err(error(
                "cannot combine force fields with different nonbonded settings",
            ));
        }
        if let (Some(a), Some(b)) = (self.charge_model(), other.charge_model()) {
            if a != b {
                return Err(error(
                    "cannot combine force fields with different NAGL models",
                ));
            }
        }
        self.bonds.extend_from_slice(&other.bonds);
        self.angles.extend_from_slice(&other.angles);
        self.propers.extend_from_slice(&other.propers);
        self.impropers.extend_from_slice(&other.impropers);
        self.constraints.extend_from_slice(&other.constraints);
        self.vdw.extend_from_slice(&other.vdw);
        self.library.extend_from_slice(&other.library);
        if self.charge_model.is_none() {
            self.charge_model.clone_from(&other.charge_model);
        }
        Ok(())
    }

    /// Compile the bundled, version-pinned Rosemary preset.
    pub fn rosemary() -> Result<Self> {
        Self::from_offxml(include_str!("../data/rosemary.offxml"))
    }
    /// Read one UTF-8 OFFXML file and compile its supported rules.
    ///
    /// File paths are used literally, without registry lookup or network access.
    /// Neural charges require a compatible bundle matching `NAGLCharges`;
    /// complete library charges can be assigned without a bundle.
    ///
    /// ```no_run
    /// let force_field = kekule_openff::ForceField::from_file("custom.offxml")?;
    /// # Ok::<(), Box<dyn std::error::Error>>(())
    /// ```
    pub fn from_file(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let load = || Self::from_offxml(&std::fs::read_to_string(path).map_err(error)?);
        load().map_err(|e| error(format!("OFFXML {}: {e}", path.display())))
    }
    /// Compile one SMIRNOFF 0.3 document with MDL aromaticity.
    ///
    /// Constraints, ImproperTorsions, LibraryCharges and NAGLCharges may be absent.
    /// Bonds, Angles, ProperTorsions, vdW and Electrostatics remain required.
    /// Without NAGL, every molecule must have complete library charges.
    /// Supported optional attributes use specification defaults.
    /// Unknown physics, attributes and section versions fail explicitly.
    pub fn from_offxml(xml: &str) -> Result<Self> {
        let document = roxmltree::Document::parse(xml).map_err(error)?;
        let root = document.root_element();
        if root.tag_name().name() != "SMIRNOFF" || root.tag_name().namespace().is_some() {
            return Err(error("expected SMIRNOFF root"));
        }
        require(root, "version", "0.3")?;
        require(root, "aromaticity_model", "OEAroModel_MDL")?;
        allowed_attributes(root, &["version", "aromaticity_model"])?;
        let mut seen = BTreeSet::new();
        for node in root.children().filter(Node::is_element) {
            let tag = node.tag_name().name();
            if !matches!(
                tag,
                "Author"
                    | "Date"
                    | "Constraints"
                    | "Bonds"
                    | "Angles"
                    | "ProperTorsions"
                    | "ImproperTorsions"
                    | "vdW"
                    | "Electrostatics"
                    | "LibraryCharges"
                    | "NAGLCharges"
            ) {
                return Err(error(format!("unsupported SMIRNOFF handler {tag}")));
            }
            if !seen.insert(tag) {
                return Err(error(format!("duplicate section {tag}")));
            }
            let extra: &[&str] = match tag {
                "Author" | "Date" => &[],
                "Bonds" => &[
                    "potential",
                    "fractional_bondorder_method",
                    "fractional_bondorder_interpolation",
                ],
                "Angles" => &["potential"],
                "ProperTorsions" => &[
                    "potential",
                    "default_idivf",
                    "fractional_bondorder_method",
                    "fractional_bondorder_interpolation",
                ],
                "ImproperTorsions" => &["potential", "default_idivf"],
                "vdW" => &[
                    "method",
                    "potential",
                    "combining_rules",
                    "scale12",
                    "scale13",
                    "scale14",
                    "scale15",
                    "cutoff",
                    "switch_width",
                    "periodic_method",
                    "nonperiodic_method",
                ],
                "Electrostatics" => &[
                    "method",
                    "scale12",
                    "scale13",
                    "scale14",
                    "scale15",
                    "cutoff",
                    "switch_width",
                    "periodic_potential",
                    "nonperiodic_potential",
                    "exception_potential",
                ],
                "NAGLCharges" => &["model_file", "model_file_hash"],
                _ => &[],
            };
            let mut attrs = extra.to_vec();
            if !matches!(tag, "Author" | "Date") {
                attrs.push("version");
            }
            allowed_attributes(node, &attrs)?;
            if matches!(tag, "Electrostatics" | "NAGLCharges" | "Author" | "Date")
                && node.children().any(|n| n.is_element())
            {
                return Err(error(format!("unexpected child in {tag}")));
            }
        }
        let optional_section = |name| root.children().find(|n| n.has_tag_name(name));
        let section =
            |name| optional_section(name).ok_or_else(|| error(format!("missing handler {name}")));
        let bonds = section("Bonds")?;
        let angles = section("Angles")?;
        let propers = section("ProperTorsions")?;
        let impropers = optional_section("ImproperTorsions");
        let vdw = section("vdW")?;
        let electrostatics = section("Electrostatics")?;
        let charge_model = optional_section("NAGLCharges")
            .map(|nagl| {
                require(nagl, "version", "0.3")?;
                ModelIdentity::new(
                    attr(nagl, "model_file")?.to_owned(),
                    attr(nagl, "model_file_hash")?.to_owned(),
                )
            })
            .transpose()?;
        choice(
            bonds,
            "potential",
            "harmonic",
            &["harmonic", "(k/2)*(r-length)^2"],
        )?;
        choice(angles, "potential", "harmonic", &["harmonic"])?;
        for n in [Some(propers), impropers].into_iter().flatten() {
            choice(
                n,
                "potential",
                "k*(1+cos(periodicity*theta-phase))",
                &["k*(1+cos(periodicity*theta-phase))"],
            )?;
            let divisor = n.attribute("default_idivf").unwrap_or("auto");
            if divisor != "auto" && number(divisor)? <= 0.0 {
                return Err(error("default_idivf must be positive or auto"));
            }
        }
        choice(
            vdw,
            "potential",
            "Lennard-Jones-12-6",
            &["Lennard-Jones-12-6"],
        )?;
        choice(
            vdw,
            "combining_rules",
            "Lorentz-Berthelot",
            &["Lorentz-Berthelot"],
        )?;
        for (name, versions) in [
            ("Bonds", &["0.3", "0.4"][..]),
            ("Angles", &["0.3"][..]),
            ("ProperTorsions", &["0.3", "0.4"][..]),
            ("ImproperTorsions", &["0.3"][..]),
            ("Constraints", &["0.3"][..]),
            ("vdW", &["0.3", "0.4"][..]),
            ("Electrostatics", &["0.3", "0.4"][..]),
            ("LibraryCharges", &["0.3"][..]),
        ] {
            if let Some(n) = optional_section(name) {
                let version = attr(n, "version")?;
                choice(n, "version", version, versions)?;
            }
        }
        // These declarations are inert for fixed bond/torsion parameters. Actual
        // bond-order indexed attributes remain rejected by children().
        for n in [bonds, propers] {
            let default_method = if n.has_tag_name("Bonds") && attr(n, "version")? == "0.3" {
                "none"
            } else {
                "AM1-Wiberg"
            };
            choice(
                n,
                "fractional_bondorder_method",
                default_method,
                &["none", "AM1-Wiberg"],
            )?;
            choice(
                n,
                "fractional_bondorder_interpolation",
                "linear",
                &["linear"],
            )?;
        }
        let (vdw_periodic_method, vdw_nonperiodic_method) = vdw_methods(vdw)?;
        let electrostatics_periodic_method = electrostatics_method(electrostatics)?;
        let settings = NonbondedSettings {
            vdw_cutoff: quantity_default(vdw, "cutoff", "9*angstrom", NANOMETER)?,
            vdw_switch_width: quantity_default(vdw, "switch_width", "1*angstrom", NANOMETER)?,
            electrostatics_cutoff: quantity_default(
                electrostatics,
                "cutoff",
                "9*angstrom",
                NANOMETER,
            )?,
            electrostatics_switch_width: quantity_default(
                electrostatics,
                "switch_width",
                "0*angstrom",
                NANOMETER,
            )?,
            vdw_scales: scales(vdw, "0.5")?,
            electrostatics_scales: scales(electrostatics, "0.833333")?,
            vdw_periodic_method,
            vdw_nonperiodic_method,
            electrostatics_periodic_method,
            electrostatics_nonperiodic_method: ElectrostaticsMethod::Coulomb,
        };
        if settings.vdw_scales[3] != 1.0 || settings.electrostatics_scales[3] != 1.0 {
            return Err(error("only scale15=1 is supported"));
        }
        for (cutoff, width) in [
            (&settings.vdw_cutoff, &settings.vdw_switch_width),
            (
                &settings.electrostatics_cutoff,
                &settings.electrostatics_switch_width,
            ),
        ] {
            if *cutoff.value() <= 0.0 || *width.value() < 0.0 || width.value() > cutoff.value() {
                return Err(error("invalid cutoff/switch width"));
            }
        }
        let mut ff = Self {
            bonds: vec![],
            angles: vec![],
            propers: vec![],
            impropers: vec![],
            constraints: vec![],
            vdw: vec![],
            library: vec![],
            settings,
            charge_model,
        };
        for n in children(bonds, "Bond")? {
            let p = BondParameter {
                source: identity(n)?,
                length: quantity(n, "length", NANOMETER)?,
                k: quantity(n, "k", CANONICAL_FORCE_CONSTANT_UNIT)?,
            };
            if *p.length.value() <= 0.0 || *p.k.value() < 0.0 {
                return Err(error("invalid harmonic bond parameter"));
            }
            ff.bonds.push(rule(n, 2, p)?);
        }
        let angle_k = KILOJOULE_PER_MOLE
            .try_div(RADIAN.try_powi(2).map_err(error)?)
            .map_err(error)?;
        for n in children(angles, "Angle")? {
            let p = AngleParameter {
                source: identity(n)?,
                angle: quantity(n, "angle", RADIAN)?,
                k: quantity(n, "k", angle_k)?,
            };
            if !(0.0..=std::f64::consts::PI).contains(p.angle.value()) || *p.k.value() < 0.0 {
                return Err(error("invalid harmonic angle parameter"));
            }
            ff.angles.push(rule(n, 3, p)?);
        }
        for (section, improper) in [(Some(propers), false), (impropers, true)] {
            let Some(section) = section else { continue };
            let default = section.attribute("default_idivf").unwrap_or("auto");
            for n in children(section, if improper { "Improper" } else { "Proper" })? {
                let mut terms = Vec::new();
                for i in 1..=64 {
                    let name = format!("periodicity{i}");
                    if n.attribute(name.as_str()).is_none() {
                        break;
                    }
                    let periodicity = attr(n, &name)?.parse::<u32>().map_err(error)?;
                    let divisor = n.attribute(format!("idivf{i}").as_str()).unwrap_or(default);
                    let idivf = if divisor == "auto" {
                        if improper {
                            3.0
                        } else {
                            0.0
                        }
                    } else {
                        number(divisor)?
                    };
                    if periodicity == 0 || idivf < 0.0 || (idivf == 0.0 && divisor != "auto") {
                        return Err(error("invalid torsion periodicity or divisor"));
                    }
                    let phase = quantity(n, &format!("phase{i}"), RADIAN)?;
                    // The outer-atom orientation is immaterial only for phases
                    // that are integer multiples of pi. Rosemary uses pi.
                    let turns = *phase.value() / std::f64::consts::PI;
                    if improper && (turns - turns.round()).abs() > 1e-12 {
                        return Err(error("improper phase must be an integer multiple of pi"));
                    }
                    terms.push(TorsionTerm {
                        periodicity,
                        idivf,
                        phase,
                        k: quantity(n, &format!("k{i}"), KILOJOULE_PER_MOLE)?,
                    });
                }
                if terms.is_empty() {
                    return Err(error("torsion has no Fourier terms"));
                }
                let count = n
                    .attributes()
                    .filter(|a| a.name().starts_with("periodicity"))
                    .count();
                if count != terms.len() {
                    return Err(error("torsion term indices must be contiguous"));
                }
                for a in n.attributes() {
                    for prefix in ["phase", "k", "idivf"] {
                        if let Some(index) = a
                            .name()
                            .strip_prefix(prefix)
                            .and_then(|s| s.parse::<usize>().ok())
                        {
                            if index == 0 || index > terms.len() {
                                return Err(error("orphaned torsion attribute"));
                            }
                        }
                    }
                }
                let p = rule(
                    n,
                    4,
                    TorsionParameter {
                        source: identity(n)?,
                        terms,
                    },
                )?;
                if improper {
                    ff.impropers.push(p);
                } else {
                    ff.propers.push(p);
                }
            }
        }
        for n in optional_children(optional_section("Constraints"), "Constraint")? {
            let distance = n
                .attribute("distance")
                .map(|_| quantity(n, "distance", NANOMETER))
                .transpose()?;
            if distance.as_ref().is_some_and(|x| *x.value() <= 0.0) {
                return Err(error("constraint distance must be positive"));
            }
            ff.constraints.push(rule(n, 2, (identity(n)?, distance))?);
        }
        for n in children(vdw, "Atom")? {
            let sigma = match (n.attribute("sigma"), n.attribute("rmin_half")) {
                (Some(_), None) => quantity(n, "sigma", NANOMETER)?,
                (None, Some(_)) => Quantity::new(
                    *quantity(n, "rmin_half", NANOMETER)?.value() * 2.0 / 2f64.powf(1.0 / 6.0),
                    NANOMETER,
                ),
                _ => return Err(error("vdW requires exactly one of sigma and rmin_half")),
            };
            let epsilon = quantity(n, "epsilon", KILOJOULE_PER_MOLE)?;
            if *sigma.value() <= 0.0 || *epsilon.value() < 0.0 {
                return Err(error("invalid vdW parameter"));
            }
            ff.vdw.push(rule(
                n,
                1,
                VdwParameter {
                    source: identity(n)?,
                    sigma,
                    epsilon,
                },
            )?);
        }
        for n in optional_children(optional_section("LibraryCharges"), "LibraryCharge")? {
            let query = parse_smarts(attr(n, "smirks")?).map_err(error)?;
            let size = TaggedQuery::new(&query).map_err(error)?.tags().len();
            if size == 0
                || n.attributes()
                    .filter(|a| a.name().starts_with("charge"))
                    .count()
                    != size
            {
                return Err(error("library charge count must equal tag count"));
            }
            let charges = (1..=size)
                .map(|i| {
                    quantity(n, &format!("charge{i}"), ELEMENTARY_CHARGE).map(Quantity::into_value)
                })
                .collect::<Result<Vec<_>>>()?;
            ff.library.push(rule(
                n,
                size,
                LibraryCharge {
                    source: identity(n)?,
                    charges,
                },
            )?);
        }
        Ok(ff)
    }
    pub fn nonbonded_settings(&self) -> &NonbondedSettings {
        &self.settings
    }
}

fn compatible_settings(a: &NonbondedSettings, b: &NonbondedSettings) -> bool {
    let length = |x: &Quantity<f64>, y: &Quantity<f64>| {
        let x = *x.value();
        let y = *y.value();
        (x - y).abs() <= 8.0 * f64::EPSILON * x.abs().max(y.abs())
    };
    length(&a.vdw_cutoff, &b.vdw_cutoff)
        && length(&a.vdw_switch_width, &b.vdw_switch_width)
        && length(&a.electrostatics_cutoff, &b.electrostatics_cutoff)
        && length(
            &a.electrostatics_switch_width,
            &b.electrostatics_switch_width,
        )
        && a.vdw_scales == b.vdw_scales
        && a.electrostatics_scales == b.electrostatics_scales
        && a.vdw_periodic_method == b.vdw_periodic_method
        && a.vdw_nonperiodic_method == b.vdw_nonperiodic_method
        && a.electrostatics_periodic_method == b.electrostatics_periodic_method
        && a.electrostatics_nonperiodic_method == b.electrostatics_nonperiodic_method
}

fn optional_children<'a, 'input>(
    node: Option<Node<'a, 'input>>,
    name: &str,
) -> Result<Vec<Node<'a, 'input>>> {
    node.map(|n| children(n, name))
        .transpose()
        .map(Option::unwrap_or_default)
}

fn vdw_methods(node: Node<'_, '_>) -> Result<(VdwMethod, VdwMethod)> {
    let method = |name| match name {
        "cutoff" => VdwMethod::Cutoff,
        _ => VdwMethod::NoCutoff,
    };
    if attr(node, "version")? == "0.3" {
        absent(node, &["periodic_method", "nonperiodic_method"])?;
        choice(node, "method", "cutoff", &["cutoff"])?;
        Ok((VdwMethod::Cutoff, VdwMethod::NoCutoff))
    } else {
        absent(node, &["method"])?;
        Ok((
            method(choice(
                node,
                "periodic_method",
                "cutoff",
                &["cutoff", "no-cutoff"],
            )?),
            method(choice(
                node,
                "nonperiodic_method",
                "no-cutoff",
                &["cutoff", "no-cutoff"],
            )?),
        ))
    }
}

fn electrostatics_method(node: Node<'_, '_>) -> Result<ElectrostaticsMethod> {
    if attr(node, "version")? == "0.3" {
        absent(
            node,
            &[
                "periodic_potential",
                "nonperiodic_potential",
                "exception_potential",
            ],
        )?;
        Ok(match choice(node, "method", "PME", &["PME", "Coulomb"])? {
            "PME" => ElectrostaticsMethod::Ewald3DConductingBoundary,
            _ => ElectrostaticsMethod::Coulomb,
        })
    } else {
        absent(node, &["method"])?;
        choice(node, "nonperiodic_potential", "Coulomb", &["Coulomb"])?;
        choice(node, "exception_potential", "Coulomb", &["Coulomb"])?;
        Ok(
            match choice(
                node,
                "periodic_potential",
                "Ewald3D-ConductingBoundary",
                &["Ewald3D-ConductingBoundary", "PME", "Coulomb"],
            )? {
                "Coulomb" => ElectrostaticsMethod::Coulomb,
                _ => ElectrostaticsMethod::Ewald3DConductingBoundary,
            },
        )
    }
}

fn absent(node: Node<'_, '_>, names: &[&str]) -> Result<()> {
    for name in names {
        if node.attribute(*name).is_some() {
            return Err(error(format!(
                "{name} is not valid on {} version {}",
                node.tag_name().name(),
                attr(node, "version")?
            )));
        }
    }
    Ok(())
}

fn choice<'a>(
    node: Node<'a, '_>,
    name: &str,
    default: &'a str,
    supported: &[&str],
) -> Result<&'a str> {
    let value = node.attribute(name).unwrap_or(default);
    if supported.contains(&value) {
        Ok(value)
    } else {
        Err(error(format!(
            "unsupported {name}={value:?} on {}; supported: {}",
            node.tag_name().name(),
            supported.join(", ")
        )))
    }
}

fn children<'a, 'input>(node: Node<'a, 'input>, name: &str) -> Result<Vec<Node<'a, 'input>>> {
    node.children()
        .filter(Node::is_element)
        .map(|n| {
            if !n.has_tag_name(name) || n.tag_name().namespace().is_some() {
                return Err(error(format!(
                    "unexpected {} in {}",
                    n.tag_name().name(),
                    node.tag_name().name()
                )));
            }
            if n.children().any(|child| child.is_element()) {
                return Err(error(format!("unexpected child element in {name}")));
            }
            if n.attributes().any(|a| a.name().contains("bondorder")) {
                return Err(error(
                    "fractional-bond-order interpolation is not supported",
                ));
            }
            for a in n.attributes() {
                let key = a.name();
                let common = matches!(key, "id" | "smirks" | "parent_id" | "description" | "name");
                let indexed = |prefix: &str| {
                    key.strip_prefix(prefix)
                        .is_some_and(|s| !s.is_empty() && s.bytes().all(|c| c.is_ascii_digit()))
                };
                let supported = match name {
                    "Bond" => matches!(key, "length" | "k"),
                    "Angle" => matches!(key, "angle" | "k"),
                    "Proper" | "Improper" => ["periodicity", "phase", "k", "idivf"]
                        .iter()
                        .any(|p| indexed(p)),
                    "Atom" => matches!(key, "sigma" | "rmin_half" | "epsilon"),
                    "Constraint" => key == "distance",
                    "LibraryCharge" => indexed("charge"),
                    _ => false,
                };
                if a.namespace().is_some() || (!common && !supported) {
                    return Err(error(format!("unsupported {name} attribute {key}")));
                }
            }
            Ok(n)
        })
        .collect()
}
fn attr<'a>(n: Node<'a, '_>, key: &str) -> Result<&'a str> {
    n.attribute(key)
        .ok_or_else(|| error(format!("missing {key} on {}", n.tag_name().name())))
}
fn allowed_attributes(node: Node<'_, '_>, allowed: &[&str]) -> Result<()> {
    if node.tag_name().namespace().is_some() {
        return Err(error("XML namespaces are not supported in OFFXML"));
    }
    for a in node.attributes() {
        if a.namespace().is_some() || !allowed.contains(&a.name()) {
            return Err(error(format!(
                "unsupported {} attribute {}",
                node.tag_name().name(),
                a.name()
            )));
        }
    }
    Ok(())
}
fn require(n: Node<'_, '_>, key: &str, expected: &str) -> Result<()> {
    if attr(n, key)? == expected {
        Ok(())
    } else {
        Err(error(format!(
            "unsupported {key}={:?} on {}; expected {expected:?}",
            attr(n, key)?,
            n.tag_name().name()
        )))
    }
}
fn identity(n: Node<'_, '_>) -> Result<ParameterIdentity> {
    Ok(ParameterIdentity {
        id: n.attribute("id").unwrap_or("").into(),
        smirks: attr(n, "smirks")?.into(),
    })
}
fn rule<P>(n: Node<'_, '_>, size: usize, parameter: P) -> Result<Rule<P>> {
    let query = parse_smarts(attr(n, "smirks")?).map_err(error)?;
    let tagged = TaggedQuery::new(&query).map_err(error)?;
    if tagged.tags().iter().map(|(t, _)| *t).ne(1..=size as u32) {
        return Err(error("SMIRKS tags must be exactly 1..N"));
    }
    Ok(Rule { query, parameter })
}
fn number(s: &str) -> Result<f64> {
    let x = s.parse::<f64>().map_err(error)?;
    if x.is_finite() {
        Ok(x)
    } else {
        Err(error("nonfinite numeric parameter"))
    }
}
fn scales(n: Node<'_, '_>, scale14: &str) -> Result<[f64; 4]> {
    let mut values = [0.0; 4];
    for (i, x) in values.iter_mut().enumerate() {
        *x = number(
            n.attribute(format!("scale1{}", i + 2).as_str())
                .unwrap_or(["0", "0", scale14, "1"][i]),
        )?;
        if !(0.0..=1.0).contains(x) {
            return Err(error("invalid nonbonded scale"));
        }
    }
    Ok(values)
}
fn quantity(n: Node<'_, '_>, key: &str, target: Unit) -> Result<Quantity<f64>> {
    quantity_default(n, key, attr(n, key)?, target)
}
fn quantity_default(
    n: Node<'_, '_>,
    key: &str,
    default: &str,
    target: Unit,
) -> Result<Quantity<f64>> {
    parse_quantity(n.attribute(key).unwrap_or(default), target).map_err(|e| {
        error(format!(
            "{} id={:?} attribute {key}: {e}",
            n.tag_name().name(),
            n.attribute("id").unwrap_or("")
        ))
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn unit_conversion_and_dimensional_rejection() {
        let bond = parse_quantity(
            "2 * kilocalorie_per_mole * angstrom ** -2",
            CANONICAL_FORCE_CONSTANT_UNIT,
        )
        .unwrap();
        assert!((*bond.value() - 836.8).abs() < 1e-10);
        assert!(
            (*parse_quantity("180 * degree", RADIAN).unwrap().value() - std::f64::consts::PI).abs()
                < 1e-14
        );
        assert!(parse_quantity("1 * angstrom", KILOJOULE_PER_MOLE).is_err());
        assert!(parse_quantity("NaN * angstrom", NANOMETER).is_err());
    }
}
