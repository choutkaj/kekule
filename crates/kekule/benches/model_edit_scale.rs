//! Informational cost of one edit of a ligand in a large model: deleting one
//! atom of a ligand, attaching one, or changing one's element or formal
//! charge, then publishing the edited draft and perceiving the result, while
//! an earlier snapshot stays alive.
//!
//! The model is 66 merged copies of an externally supplied mmCIF entry, each
//! appended as separate chains and molecules. SeekDelta's request K56 used
//! PDB 9UOH (KRas G12C with GDP, Mg2+ and an inhibitor, 1,526 atoms), which
//! makes 100,716 atoms:
//!
//! ```sh
//! curl -sSfLO https://files.rcsb.org/download/9UOH.cif
//! cargo bench -p kekule --bench model_edit_scale -- 9UOH.cif A1L9E
//! ```
//!
//! The arguments are the mmCIF path and the ligand's component ID; the
//! `KEKULE_EDIT_SCALE_CIF` and `KEKULE_EDIT_SCALE_LIGAND` environment variables
//! may supply them instead. Without a path the benchmark reports nothing.

use std::error::Error;
use std::hint::black_box;
use std::time::{Duration, Instant};

use kekule::geometry::Point3;
use kekule::mmcif;
use kekule::molecule::{Atom, BondOrder, Element};
use kekule::structure::ModelEditor;
use kekule::topology::{InstanceAtomId, TopologyAtomIndex};
use kekule::units::{Quantity, ANGSTROM};

const COPIES: usize = 66;
const ROUNDS: usize = 15;

fn main() -> Result<(), Box<dyn Error>> {
    let mut args = std::env::args()
        .skip(1)
        .filter(|argument| !argument.starts_with("--"));
    let Some(path) = args
        .next()
        .or_else(|| std::env::var("KEKULE_EDIT_SCALE_CIF").ok())
    else {
        println!("model_edit_scale: no mmCIF path given; nothing measured");
        return Ok(());
    };
    let ligand = args
        .next()
        .or_else(|| std::env::var("KEKULE_EDIT_SCALE_LIGAND").ok())
        .unwrap_or_else(|| "A1L9E".to_owned());

    let text = std::fs::read_to_string(&path)?;
    let copy = mmcif::parse_str(&text, Default::default())?
        .interpret(Default::default())?
        .into_model();
    let mut editor = ModelEditor::new();
    for _ in 0..COPIES {
        editor.append_model(&copy)?;
    }
    let mut model = editor.finish()?;
    model.perceive()?;
    let topology = model.topology();

    // An atom of the ligand in the middle copy, bonded to exactly one other
    // atom, so deleting it splits nothing.
    let per_copy = copy.atom_count();
    let middle = (COPIES / 2) * per_copy..(COPIES / 2 + 1) * per_copy;
    let terminal = topology
        .atoms()
        .filter(|atom| middle.contains(&atom.index().index()))
        .find(|atom| {
            atom.residue().is_some_and(|r| r.name() == ligand) && atom.neighbors().count() == 1
        })
        .ok_or("no terminal ligand atom in the middle copy")?;
    let target = terminal.id();
    let target_index = terminal.index();
    let ligand_atoms = terminal.instance().molecule().atom_count();

    println!(
        "informational only; target={}-{}; atoms={}; definitions={}; instances={}; \
         ligand {ligand} molecule atoms={ligand_atoms}; minimum and median of {ROUNDS} rounds",
        std::env::consts::OS,
        std::env::consts::ARCH,
        model.atom_count(),
        topology.definition_count(),
        topology.instance_count(),
    );
    println!("| edit | phase | min ms | median ms |");
    println!("|---|---|---:|---:|");

    let edits: [(&str, Edit); 4] = [
        ("delete one ligand atom", delete_atom),
        ("attach one ligand atom", add_atom),
        ("set one ligand atom's element", set_element),
        ("set one ligand atom's charge", set_charge),
    ];
    for (name, edit) in edits {
        let mut phases = [const { Vec::new() }; 6];
        for _ in 0..ROUNDS {
            let earlier = model.clone();
            let start = Instant::now();
            let mut editor = model.to_editor();
            let built = start.elapsed();
            let start = Instant::now();
            edit(&mut editor, target)?;
            let edited = start.elapsed();
            // A host's remap of one per-atom column from old to new dense
            // indices goes through the draft handles of every source atom.
            let start = Instant::now();
            let handles = earlier
                .topology()
                .atom_ids()
                .iter()
                .map(|&atom| editor.topology().atom_handle(atom))
                .collect::<Vec<_>>();
            let mut remapped = start.elapsed();
            let start = Instant::now();
            let (mut published, correspondence) = editor.finish_with_correspondence()?;
            let finished = start.elapsed();
            let start = Instant::now();
            published.perceive()?;
            let perceived = start.elapsed();
            let start = Instant::now();
            let map = handles
                .iter()
                .map(|handle| {
                    handle
                        .and_then(|handle| correspondence.atom(handle))
                        .map(|(_, index)| index)
                })
                .collect::<Vec<Option<TopologyAtomIndex>>>();
            remapped += start.elapsed();
            black_box((&published, &map, target_index));
            let cycle = built + edited + finished + perceived;
            for (phase, time) in phases
                .iter_mut()
                .zip([built, edited, finished, perceived, remapped, cycle])
            {
                phase.push(time);
            }
        }
        for (phase, times) in [
            "build the editor",
            "edit",
            "publish (finish_with_correspondence)",
            "perceive",
            "remap through handles",
            "whole cycle (build, edit, publish, perceive)",
        ]
        .into_iter()
        .zip(&mut phases)
        {
            times.sort_unstable();
            println!(
                "| {name} | {phase} | {:.2} | {:.2} |",
                millis(times[0]),
                millis(times[times.len() / 2])
            );
        }
    }
    Ok(())
}

type Edit = fn(&mut ModelEditor, InstanceAtomId) -> Result<(), Box<dyn Error>>;

fn delete_atom(editor: &mut ModelEditor, atom: InstanceAtomId) -> Result<(), Box<dyn Error>> {
    let handle = editor
        .topology()
        .atom_handle(atom)
        .ok_or("no source atom")?;
    editor.delete_atom(handle)?;
    Ok(())
}

/// Adds a carbon 1.5 Å from `atom`, bonded to it.
fn add_atom(editor: &mut ModelEditor, atom: InstanceAtomId) -> Result<(), Box<dyn Error>> {
    let handle = editor
        .topology()
        .atom_handle(atom)
        .ok_or("no source atom")?;
    let near = editor
        .position(handle)
        .ok_or("no source position")?
        .value_in(ANGSTROM)?;
    let added = editor.add_atom(
        Atom::new(Element::from_symbol("C").ok_or("no carbon")?),
        Quantity::new(Point3::new(near.x + 1.5, near.y, near.z), ANGSTROM),
    )?;
    editor.add_bond(handle, added, BondOrder::Single)?;
    Ok(())
}

/// Changes `atom`'s element: carbon to nitrogen, anything else to carbon.
fn set_element(editor: &mut ModelEditor, atom: InstanceAtomId) -> Result<(), Box<dyn Error>> {
    let handle = editor
        .topology()
        .atom_handle(atom)
        .ok_or("no source atom")?;
    let current = editor.topology().atom(handle).ok_or("no atom")?.clone();
    let symbol = if current.element().symbol() == "C" {
        "N"
    } else {
        "C"
    };
    let element = Element::from_symbol(symbol).ok_or("no element")?;
    editor.replace_atom(handle, current.with_element(element))?;
    Ok(())
}

/// Raises `atom`'s formal charge by one.
fn set_charge(editor: &mut ModelEditor, atom: InstanceAtomId) -> Result<(), Box<dyn Error>> {
    let handle = editor
        .topology()
        .atom_handle(atom)
        .ok_or("no source atom")?;
    let current = editor.topology().atom(handle).ok_or("no atom")?.clone();
    let charge = current.formal_charge() + 1;
    editor.replace_atom(handle, current.with_formal_charge(charge))?;
    Ok(())
}

fn millis(time: Duration) -> f64 {
    time.as_secs_f64() * 1e3
}
