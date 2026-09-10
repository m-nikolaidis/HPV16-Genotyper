import importlib
import io
import pathlib
import sys
import types
import unittest


# Biopython installations without Bio.Align.Applications should still allow
# this module's calculation function to be tested.  The production import is
# intentionally unchanged.
applications = types.ModuleType("Bio.Align.Applications")
applications.MuscleCommandline = object
sys.modules.setdefault("Bio.Align.Applications", applications)
sys.path.insert(0, str(pathlib.Path(__file__).parents[1]))
simplot = importlib.import_module("simplot")


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


if __name__ == "__main__":
    unittest.main()
