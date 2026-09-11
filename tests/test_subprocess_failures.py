import pathlib
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

from hpv16genotyper import appFunctions


class _PoolSpy:
    def __init__(self, *args):
        self.closed = False
        self.terminated = False
        self.joined = False

    def close(self):
        self.closed = True

    def terminate(self):
        self.terminated = True

    def join(self):
        self.joined = True


class ExternalToolFailureTests(unittest.TestCase):
    def test_managed_thread_pool_closes_and_joins_on_success(self):
        pool = _PoolSpy()
        with patch.object(appFunctions, "ThreadPool", return_value=pool):
            with appFunctions._managed_thread_pool(2) as managed:
                self.assertIs(managed, pool)

        self.assertTrue(pool.closed)
        self.assertFalse(pool.terminated)
        self.assertTrue(pool.joined)

    def test_active_timer_stops_when_scope_raises(self):
        events = []

        class TimerSpy:
            def startTimer(self):
                events.append("start")

            def stopTimer(self):
                events.append("stop")

        with self.assertRaises(RuntimeError):
            with appFunctions._active_timer(TimerSpy()):
                events.append("body")
                raise RuntimeError("failed")

        self.assertEqual(events, ["start", "body", "stop"])

    def test_run_external_includes_exit_code_and_stderr_summary(self):
        completed = type("Completed", (), {"returncode": 17, "stderr": "bad\ninput"})()
        with patch.object(appFunctions.subprocess, "run", return_value=completed):
            with self.assertRaises(appFunctions.ExternalToolError) as raised:
                appFunctions._run_external(["tool", "--input", "file"], "MUSCLE")

        self.assertIn("MUSCLE", str(raised.exception))
        self.assertIn("exit code 17", str(raised.exception))
        self.assertIn("bad input", str(raised.exception))

    def test_blastn_failure_stops_timer_and_joins_pool(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            query = root / "query.fa"
            query.write_text(">sample\nAAAA\n")
            (root / ".tmp").mkdir()
            (root / "BlastResults").mkdir()
            params = pd.DataFrame(
                {
                    "Value": [
                        1,
                        root,
                        root / "reference.fa",
                        1e-5,
                        7,
                    ]
                },
                index=[
                    "num_threads",
                    "out",
                    "HPV_filter_db",
                    "Gene_identification_Evalue",
                    "Gene_identification_Word_size",
                ],
            )
            functions = appFunctions.MainFunctions()
            stopped = []
            functions.stopTimerSignal.connect(stopped.append)
            pool = _PoolSpy()

            def run(command, tool, **kwargs):
                if tool == "blastn":
                    raise appFunctions.ExternalToolError(tool, command, 3, "query failed")
                return None

            with patch.object(appFunctions, "ThreadPool", return_value=pool):
                with patch.object(appFunctions, "_run_external", side_effect=run):
                    with self.assertRaises(appFunctions.ExternalToolError):
                        functions.blastn_search(
                            params,
                            query,
                            "HPV16 filter",
                            "makeblastdb",
                            "blastn",
                            total=1,
                        )

            self.assertTrue(pool.terminated)
            self.assertTrue(pool.joined)
            self.assertEqual(stopped, [1])


if __name__ == "__main__":
    unittest.main()
