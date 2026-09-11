"""Failure isolation and evidence preservation for the external probe transport."""

from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch

import probe_support


class ProbeTransport(unittest.TestCase):
    def test_process_and_protocol_failures_retain_records(self):
        lines = ["CCO", "C("]
        for process, status in (
            (subprocess.CompletedProcess([], 0, "{}\n[]\n", ""), "protocol_error"),
            (subprocess.CompletedProcess([], 1, "", "failed"), "process_error"),
        ):
            with patch.object(probe_support.subprocess, "run", return_value=process):
                values = probe_support.run_probe(Path("probe"), lines, 1)
            self.assertEqual([value["status"] for value in values], [status] * 2)

    def test_launch_failure_retains_every_input(self):
        with patch.object(probe_support.subprocess, "run", side_effect=PermissionError("cannot execute")):
            values = probe_support.run_probe(Path("probe"), ["CCO", "C("], 1)
        self.assertEqual([value["status"] for value in values], ["process_error"] * 2)
        self.assertTrue(all("cannot execute" in value["message"] for value in values))

    def test_failed_batch_isolates_inputs_in_order(self):
        timeout = subprocess.TimeoutExpired("probe", 1)
        truncated = subprocess.CompletedProcess([], 0, '{"status":"ok"}\n', "")
        for batch_failure in (timeout, truncated):
            with patch.object(probe_support.subprocess, "run", side_effect=[
                batch_failure, timeout,
                subprocess.CompletedProcess([], 0, '{"status":"ok","id":2}\n', ""),
            ]) as run:
                values = probe_support.run_probe(Path("probe"), ["first", "second"], 1)
            self.assertEqual(values, [{"status": "timeout", "seconds": 1}, {"status": "ok", "id": 2}])
            self.assertEqual([call.kwargs["input"] for call in run.call_args_list],
                             ["first\nsecond\n", "first\n", "second\n"])

    def test_report_refuses_to_replace_existing_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "report.json"
            probe_support.write_report(path, {"first": True})
            original = path.read_bytes()
            with self.assertRaises(FileExistsError):
                probe_support.write_report(path, {"replacement": True})
            self.assertEqual(path.read_bytes(), original)


if __name__ == "__main__":
    unittest.main()
