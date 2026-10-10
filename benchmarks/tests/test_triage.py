"""Comparison, signatures, known-difference triage and error kinds (stdlib only)."""
import gzip
import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from kbench import compare, errors, known, run  # noqa: E402
from kbench.datasets import DatasetError  # noqa: E402
from kbench.references import observers  # noqa: E402
from kbench.datasets import Record, SourceFile  # noqa: E402
from kbench.tasks import TASKS, WRITTEN_FORMATS, task  # noqa: E402

CONTEXT = {
    "atoms": {"0": "C|ar", "1": "N|ar", "2": "O|-1"},
    "bonds": {"0-1": "C-N|ar", "1-2": "N-O|1"},
    "source_aromatic": [[0, 1]],
    "normalized": False,
}


def signatures(task_id, kekule, reference, options=None):
    return [d.signature for d in compare.diff(task_id, kekule, reference, CONTEXT, options or {})]


class Compare(unittest.TestCase):
    def test_equal_facts_agree_regardless_of_order(self):
        facts = [["atom", 0, "element", "C"], ["bond", 0, 1, "order", 2]]
        self.assertEqual(signatures("small.parse", facts, list(reversed(facts))), [])

    def test_source_aromatic_bonds_compare_presence_not_kekule_form(self):
        self.assertEqual(signatures("small.parse", [["bond", 0, 1, "order", 2]], [["bond", 0, 1, "order", "aromatic"]]), [])
        self.assertEqual(
            signatures("small.parse", [["bond", 0, 1, "order", 3]], [["bond", 0, 1, "order", "aromatic"]]),
            ["bond.order|C-N|ar|aromatic→3"],
        )
        self.assertEqual(
            signatures("small.parse", [["bond", 1, 2, "order", 2]], [["bond", 1, 2, "order", 1]]),
            ["bond.order|N-O|1|1→2"],
        )

    def test_hydrogens_compare_only_where_both_are_defined(self):
        self.assertEqual(signatures("small.parse", [["atom", 0, "hydrogens", None]], [["atom", 0, "hydrogens", 1]]), [])
        self.assertEqual(signatures("small.parse", [["atom", 0, "hydrogens", 2]], [["atom", 0, "hydrogens", None]]), [])
        self.assertEqual(
            signatures("small.parse", [["atom", 1, "hydrogens", 0]], [["atom", 1, "hydrogens", 1]]),
            ["atom.hydrogens|N|ar|1→0"],
        )

    def test_resonance_contributors_compare_modulo_kekule_form(self):
        def facts(*contributors):
            return [["resonance", "contributors", 0, "count", len(contributors)],
                    ["resonance", "contributors", 0, "set", list(contributors)]]

        # Benzene's two Kekulé forms are one structure on either side.
        first, second = "|0-1:2,2-3:2,4-5:2", "|0-5:2,1-2:2,3-4:2"
        self.assertEqual(signatures("small.resonance", facts(first), facts(second)), [])
        self.assertEqual(signatures("small.resonance", facts(first), facts(first, second)), [])
        # Moving a charge is a different contributor: the allyl cation.
        self.assertEqual(
            signatures("small.resonance", facts("2:1|0-1:2"), facts("2:1|0-1:2", "0:1|1-2:2")),
            ["resonance.contributors.count||fewer", "resonance.contributors.set||≠"],
        )

    def test_transitions_describe_presence_counts_lists_and_types(self):
        self.assertEqual(
            signatures("small.stereo.candidates", [], [["candidate", "axis", [0, 1], "present", True]]),
            ["candidate.axis.present|C-N|ar|true→missing"],
        )
        self.assertEqual(
            signatures("small.smarts", [["query", 0, "count", 3]], [["query", 0, "count", 4]], {"queries": ["[#6]"]}),
            ["query.count|query=[#6]|fewer"],
        )
        self.assertEqual(signatures("small.symmetry", [["class", 0, "atoms", [0, 2]]], [["class", 0, "atoms", [0]]]),
                         ["class.atoms||≠"])
        self.assertEqual(signatures("small.rings", [["rings", 6, "count", 1]], []), ["rings.count|size=6|missing→1"])
        self.assertEqual(signatures("small.parse", [["atom", 2, "charge", True]], [["atom", 2, "charge", 1]]),
                         ["atom.charge|O|-1|1→true"])

    def test_coordinates_agree_within_half_a_printed_step(self):
        reference = [["atom", 0, "xyz", [1.0, 2.0, 3.0]]]
        self.assertEqual(signatures("small.write.molfile-v2000", [["atom", 0, "xyz", [1.00005, 2.0, 3.0]]], reference), [])
        self.assertEqual(signatures("small.write.molfile-v2000", [["atom", 0, "xyz", [1.0002, 2.0, 3.0]]], reference),
                         ["atom.xyz|C|ar|≠"])
        difference = compare.diff("small.write.sdf", [["atom", 0, "xyz", [1.0002, 2.0, 3.0]]], reference, CONTEXT, {})[0]
        self.assertAlmostEqual(difference.delta(), 0.0002)

    def test_dssp_partners_compare_by_kind_and_energies_only_for_the_same_partner(self):
        names = 14
        reference = ["H", 1.0, 2.0, None, None, 0.1, "A:5", -2.0] + [None] * (names - 8)
        kekule = ["H", 1.0, 2.0, None, None, 0.1, "B:?", -1.0] + [None] * (names - 8)
        found = signatures("bio.dssp", [["dssp", "A", 3, kekule]], [["dssp", "A", 3, reference]])
        self.assertEqual(found, ["dssp.acceptor_1||residue→unnumbered"])
        self.assertEqual(signatures("bio.dssp", [["dssp", "B", None, kekule]], []), ["dssp||missing→present"])

    def test_connectivity_compares_bonds_between_sites_both_sides_kept(self):
        kekule = [["bond", "1", "2", "order", 1], ["bond", "2", "9", "order", 1], ["molecule", "1", "sites", ["1", "2", "9"]]]
        reference = [["bond", "1", "2", "order", 2], ["molecule", "1", "sites", ["1", "2"]], ["molecule", "7", "sites", ["7"]]]
        self.assertEqual(signatures("bio.connectivity", kekule, reference), ["bond.order||2→1"])


KNOWN_TOML = """
[[known]]
id = "axis-candidates"
task = "small.stereo.candidates"
match = ["candidate.axis.present|*|true→missing"]
verdict = "kekule-gap"
reason = "No axis candidates yet"
examples = ["rdkit:a.sdf"]
max_cases = { curated = 3, full = 10 }

[[known]]
id = "xyz"
task = "small.write.*"
match = ["atom.xyz|*|≠"]
verdict = "kekule-bug"
reason = "Coordinates drift"
examples = ["rdkit:b.sdf"]
max_cases = { curated = 5 }
max_abs_delta = 0.001
"""


def case(record_id, task_id, kekule, reference):
    file = SourceFile(f"{record_id}.sdf", "sdf", "0" * 64, "https://x")
    value = run.Case(f"{record_id}|{task_id}", task(task_id), Record(record_id, "s", "t", record_id, (file,)), file,
                     "x", {})
    value.kekule = {"status": "ok", "value": kekule}
    value.reference = {"status": "ok", "facts": reference}
    return value


class Triage(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.path = Path(self.temporary.name) / "known.toml"
        self.path.write_text(KNOWN_TOML, encoding="utf-8")

    def tearDown(self):
        self.temporary.cleanup()

    def classify(self, cases):
        run.compare_cases(cases, {c.file.path: {"status": "ok", "context": CONTEXT} for c in cases})
        return run.triage(cases, known.load(self.path), "curated")

    def test_known_new_stale_and_bounds(self):
        axis = [["candidate", "axis", [0, 1], "present", True]]
        cases = [
            case("rdkit:a.sdf", "small.stereo.candidates", [], axis),
            case("rdkit:b.sdf", "small.write.sdf", [["atom", 0, "xyz", [1.0, 0.0, 0.0]]], [["atom", 0, "xyz", [1.0, 0.0, 0.0]]]),
            case("rdkit:c.sdf", "small.stereo.candidates", [["candidate", "tetrahedral", [2], "present", True]], []),
        ]
        summary = self.classify(cases)
        self.assertEqual([c.outcome for c in cases], ["known", "agree", "new"])
        self.assertEqual(summary["new"][0]["signature"], "candidate.tetrahedral.present|O|-1|missing→true")
        self.assertEqual(summary["stale"], ["xyz: example rdkit:b.sdf no longer shows the difference"])
        self.assertEqual(summary["known"]["axis-candidates"]["cases"], 1)
        self.assertTrue(summary["failed"])

    def test_bounds_and_numeric_limits_are_one_way(self):
        drift = [["atom", 0, "xyz", [1.01, 0.0, 0.0]]]
        cases = [case(f"rdkit:{n}.sdf", "small.write.sdf", drift, [["atom", 0, "xyz", [1.0, 0.0, 0.0]]])
                 for n in "bdefgh"]
        summary = self.classify(cases)
        self.assertIn("xyz: 6 cases exceed max_cases.curated = 5", summary["violations"])
        self.assertTrue(any("exceeds 0.001" in v for v in summary["violations"]))
        changed = known.lower_bounds(self.path, "curated", {"axis-candidates": 1, "xyz": 9})
        self.assertEqual(changed, ["axis-candidates"])
        text = self.path.read_text(encoding="utf-8")
        self.assertIn("max_cases = { curated = 1, full = 10 }", text)
        self.assertIn("max_cases = { curated = 5 }", text)

    def test_entries_cover_only_runs_of_matching_tasks(self):
        entry = known.Entry("bio-x", "bio.*", ["*"], "kekule-gap", "", [], {"full": 3})
        self.assertTrue(entry.covers(["bio.hierarchy", "bio.dssp"]))
        self.assertFalse(entry.covers(["small.parse", "small.rings"]))

    def test_one_sided_errors_are_differences_and_both_errors_agree_on_rejection(self):
        failed = case("rdkit:a.sdf", "small.parse", [], [])
        failed.kekule = {"status": "error", "error": {"kind": "parse", "message": "unsupported element symbol `CL`"}}
        both = case("rdkit:c.sdf", "small.parse", [], [])
        both.kekule = {"status": "error", "error": {"kind": "parse", "message": "x"}}
        both.reference = {"status": "error", "error": {"kind": "parse", "message": "y"}}
        summary = self.classify([failed, both])
        self.assertEqual(summary["new"][0]["signature"], "error|parse:uppercase-element|ok→parse:uppercase-element")
        self.assertEqual(both.outcome, "both-error")

    def test_parse_differences_block_dependent_tasks(self):
        parse = case("rdkit:a.sdf", "small.parse", [["bond", 1, 2, "order", 2]], [["bond", 1, 2, "order", 1]])
        rings = case("rdkit:a.sdf", "small.rings", [], [])
        self.classify([rings, parse])
        self.assertEqual((rings.outcome, rings.reason), ("blocked", "parse differs"))


class StoredReferencesTest(unittest.TestCase):
    def write(self, path, observers_fingerprint):
        lines = [
            json.dumps({"dataset": "small", "version": 1, "observers": observers_fingerprint, "tools": {}}),
            json.dumps({"prepared": "single-conformer|pdb/x.cif", "text": "data_x\n"}),
            json.dumps({"key": "rdkit:parse|a.smi|0", "value": {"status": "ok", "facts": []}}),
        ]
        path.write_bytes(gzip.compress(("\n".join(lines) + "\n").encode()))

    def test_stored_references_load_and_report_missing_entries(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "references.jsonl.gz"
            self.write(path, observers.fingerprint("small"))
            stored = run.StoredReferences(path, "small", 1)
            self.assertEqual(stored.get("rdkit:parse|a.smi|0"), {"status": "ok", "facts": []})
            self.assertEqual(stored.get("rdkit:rings|a.smi|0")["error"]["kind"], "no-reference")
            self.assertEqual(stored.prepared, {"single-conformer|pdb/x.cif": "data_x\n"})

    def test_references_from_other_observer_code_are_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "references.jsonl.gz"
            self.write(path, {"rdkit": "0" * 64})
            with self.assertRaises(DatasetError):
                run.StoredReferences(path, "small", 1)

    def test_references_for_another_dataset_revision_are_refused(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "references.jsonl.gz"
            self.write(path, observers.fingerprint("small"))
            for dataset, version in (("small", 2), ("bio", 1)):
                with self.subTest(dataset=dataset, version=version), self.assertRaises(DatasetError):
                    run.StoredReferences(path, dataset, version)

    def test_fingerprints_cover_only_the_observers_a_dataset_uses(self):
        self.assertEqual(sorted(observers.fingerprint("small")), ["rdkit"])
        self.assertEqual(sorted(observers.fingerprint("bio")), ["bio"])


class TaskTable(unittest.TestCase):
    def test_ids_are_unique_and_writers_have_readable_outputs(self):
        ids = [t.id for t in TASKS]
        self.assertEqual(len(ids), len(set(ids)))
        for t in TASKS:
            if t.writes:
                self.assertIn(t.writes, WRITTEN_FORMATS)
                self.assertEqual((t.reference, t.reference_options), ("rdkit:parse", {"coordinates": True}))


class ErrorKinds(unittest.TestCase):
    def test_messages_refine_stage_kinds(self):
        self.assertEqual(errors.refine("kekule", "parse", "R-group query atom"), "parse:query-atom")
        self.assertEqual(errors.refine("reference", "sanitize", "Explicit valence for atom # 1 N, 4"), "sanitize:valence")
        self.assertEqual(errors.refine("kekule", "parse", "unclosed ring closure"), "parse")
        self.assertEqual(errors.refine("kekule", "interpret", "unknown atom-site element `D`"), "interpret:deuterium")
        self.assertEqual(errors.refine("kekule", "parse", "unsupported element symbol `CL`"), "parse:uppercase-element")
        self.assertEqual(errors.refine("kekule", "parse", "unsupported element symbol `R1`"), "parse:query-atom")
        self.assertEqual(errors.refine("kekule", "parse", "unsupported element symbol `R`"), "parse:query-atom")
        self.assertEqual(errors.refine("kekule", "parse", "unsupported element symbol `X`"), "parse:query-atom")
        self.assertEqual(
            errors.refine("kekule", "parse", "source-stereo canonicalization reported 1 issue(s): "
                          "[UnassembledTetrahedralBondMark { bond: BondId(5), kind: WedgeDown }]"),
            "parse:unassembled-wedge",
        )
        self.assertEqual(errors.refine("kekule", "resonance", "resonance work limit exceeded (100000000)"),
                         "resonance:resource-limit")
        self.assertEqual(
            errors.refine("reference", "resonance", "enumeration work limit exceeded (5000 progress steps)"),
            "resonance:resource-limit",
        )
        self.assertEqual(
            errors.refine("kekule", "write", "V2000 cannot encode a coordinate in its fixed-width atom field"),
            "write:field-width",
        )


if __name__ == "__main__":
    unittest.main()
