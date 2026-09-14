import importlib
import io
import unittest
from unittest.mock import patch

simplot = importlib.import_module("hpv16genotyper.simplot")


class CalculateSimilaritiesTests(unittest.TestCase):
    def test_uses_full_windows_and_zero_based_midpoints(self):
        alignment = io.StringIO(
            ">query\nAAAAAACCCC\n>database\nAAAAAACTCC\n"
        )

        positions, results, alignment_length, query_length = simplot.calculate_similarities(
            "query", alignment, step=4, window=4
        )

        self.assertEqual(positions, [1.5, 5.5])
        self.assertEqual(results, {"database": [100.0, 75.0]})
        self.assertEqual(alignment_length, 10)
        self.assertEqual(query_length, 10)

    def test_preserves_exact_p_distance_before_identity_conversion(self):
        alignment = io.StringIO(
            ">query\nAAAAAA\n>database\nAAAATA\n"
        )

        positions, results, _, _ = simplot.calculate_similarities(
            "query", alignment, step=1, window=6
        )

        self.assertEqual(results["database"][0], (1 - (1 / 6)) * 100)
        self.assertEqual(positions, [2.5])


class AlignTests(unittest.TestCase):
    def test_runs_muscle_profile_alignment_with_existing_output_path(self):
        seqs = {"query": type("Record", (), {"seq": "acgt"})()}

        with patch.object(
            simplot, "_isolate_sequence", return_value="query.fa"
        ) as isolate:
            with patch.object(simplot, "_run_external") as run_external:
                result = simplot.align(
                    "muscle",
                    "profile.fa",
                    seqs,
                    "query",
                    "tmp",
                )

        self.assertEqual(result, "query.fa")
        isolate.assert_called_once_with(seqs, "query", "tmp")
        run_external.assert_called_once_with(
            [
                "muscle",
                "-profile",
                "-in1",
                "query.fa",
                "-in2",
                "profile.fa",
                "-out",
                "query.fa",
            ],
            "MUSCLE",
        )


if __name__ == "__main__":
    unittest.main()
