"""The task table: the only place that routes dataset files to observers.

Each task is observed by Kekule (`benchmarks/observer`) and by one reference
toolkit, which emit the same facts. A task with `requires` depends on the
represented chemistry `small.parse` compares; when parse already differs in
those fact kinds, the task's case is blocked rather than counted again.
"""
from __future__ import annotations

from dataclasses import dataclass, field

SMALL_FORMATS = ("smiles", "mol", "sdf")
COORDINATE_FORMATS = ("mol", "sdf")
# Represented chemistry other tasks build on; stereo matters only where named.
GRAPH = frozenset({"atom", "bond", "component"})
STEREO = GRAPH | {"stereo", "group"}
# Fixed SMARTS target panel: curated records of these strata, first file each.
SMARTS_TARGET_STRATA = ("chembl-clinical", "chembl-stereo-rich", "ccd-components", "rdkit-structures")


@dataclass(frozen=True)
class Task:
    id: str
    dataset: str
    formats: tuple[str, ...]
    observer: str
    reference: str
    options: dict = field(default_factory=dict)
    requires: frozenset = frozenset()
    description: str = ""
    # Strata the task applies to (all when empty).
    strata: tuple[str, ...] = ()
    # Reference-side preparation both observers then read, e.g. one conformer.
    prepare: str | None = None
    # Options of the stored reference observation, when they differ from `options`.
    reference_options: dict | None = None

    @property
    def writes(self) -> str | None:
        """Output format of a writer task, read back by the reference."""
        return self.options.get("format") if self.observer == "small.write" else None


TASKS = (
    Task("small.parse", "small", SMALL_FORMATS, "small.parse", "rdkit:parse",
         description="Represented chemistry: elements, isotopes, charges, radicals, total H, maps, "
                     "bonds, components, cleaned stereo and enhanced stereo groups"),
    Task("small.rings", "small", SMALL_FORMATS, "small.rings", "rdkit:rings", requires=GRAPH,
         description="Ring membership and the size multiset of a smallest set of smallest rings"),
    Task("small.aromaticity.rdkit", "small", SMALL_FORMATS, "small.aromaticity", "rdkit:aromaticity",
         {"model": "rdkit"}, GRAPH, "Atom and bond aromaticity, RDKit model"),
    Task("small.aromaticity.mdl", "small", SMALL_FORMATS, "small.aromaticity", "rdkit:aromaticity",
         {"model": "mdl"}, GRAPH, "Atom and bond aromaticity, MDL model"),
    Task("small.conjugation", "small", SMALL_FORMATS, "small.conjugation", "rdkit:conjugation",
         requires=GRAPH, description="Bond conjugation"),
    Task("small.resonance", "small", SMALL_FORMATS, "small.resonance", "rdkit:resonance",
         requires=GRAPH, description="Conjugated groups and default-option resonance contributors"),
    Task("small.stereo.candidates", "small", SMALL_FORMATS, "small.stereo.candidates",
         "rdkit:stereo.candidates", requires=GRAPH, description="Potential stereo centres and bonds"),
    Task("small.stereo.cip", "small", SMALL_FORMATS, "small.stereo.cip", "rdkit:stereo.cip",
         requires=STEREO, description="CIP descriptors (RDKit's new labeler)"),
    Task("small.symmetry", "small", SMALL_FORMATS, "small.symmetry", "rdkit:symmetry",
         requires=GRAPH, description="Atom equivalence classes per component"),
    Task("small.formula", "small", SMALL_FORMATS, "small.formula", "rdkit:formula",
         requires=GRAPH, description="Molecular formula terms and net charge per component"),
    Task("small.rotatable", "small", SMALL_FORMATS, "small.rotatable", "rdkit:rotatable",
         requires=GRAPH, description="Strict rotatable bonds, heavy-atom graph"),
    Task("small.smarts", "small", SMALL_FORMATS, "small.smarts", "rdkit:smarts", requires=STEREO,
         description="Every RDKit query table row matched against a fixed target panel"),
    # Writers keep atom order, so Kekule's reading of its own output is compared
    # with the reference's reading of the original. SMILES writers reorder
    # atoms; their round trips are pinned by the invariant suite instead.
    Task("small.write.molfile-v2000", "small", COORDINATE_FORMATS, "small.write", "rdkit:parse",
         {"format": "molfile-v2000"}, STEREO, "V2000 Molfile: same chemistry and coordinates",
         reference_options={"coordinates": True}),
    Task("small.write.molfile-v3000", "small", COORDINATE_FORMATS, "small.write", "rdkit:parse",
         {"format": "molfile-v3000"}, STEREO, "V3000 Molfile: same chemistry and coordinates",
         reference_options={"coordinates": True}),
    Task("small.write.sdf", "small", COORDINATE_FORMATS, "small.write", "rdkit:parse",
         {"format": "sdf"}, STEREO, "SDF record: same chemistry and coordinates",
         reference_options={"coordinates": True}),
    Task("bio.cif.syntax", "bio", ("mmcif",), "bio.cif.syntax", "bio:cif.syntax",
         description="Every item value and loop column of the CIF document, nulls kept distinct"),
    Task("bio.hierarchy", "bio", ("mmcif",), "bio.hierarchy", "bio:hierarchy",
         description="Each kept atom site of the first model: identity, element, occupancy, B and coordinates"),
    Task("bio.connectivity", "bio", ("mmcif",), "bio.connectivity", "bio:connectivity",
         description="Covalent bonds and orders between kept sites, and the molecule partition"),
    Task("bio.classification", "bio", ("mmcif",), "bio.classification", "bio:classification",
         description="Residue classes: amino acid, DNA, RNA, carbohydrate, water, ion"),
    Task("bio.dssp", "bio", ("mmcif",), "bio.dssp", "bio:dssp", prepare="single-conformer",
         description="DSSP per residue on the same single-model, single-conformer file"),
    Task("bio.superposition", "bio", ("mmcif",), "bio.superposition", "bio:superposition", strata=("nmr",),
         description="Cα RMSD of every NMR model fitted onto the first"),
)

WRITTEN_FORMATS = {
    "molfile-v2000": "mol",
    "molfile-v3000": "mol",
    "sdf": "sdf",
}


def task(task_id: str) -> Task:
    for candidate in TASKS:
        if candidate.id == task_id:
            return candidate
    raise KeyError(task_id)


def select(dataset: str, patterns: list[str] | None) -> list[Task]:
    from fnmatch import fnmatch

    chosen = [t for t in TASKS if t.dataset == dataset and (not patterns or any(fnmatch(t.id, p) for p in patterns))]
    if not chosen:
        raise ValueError(f"no {dataset} task matches {patterns}")
    return chosen
