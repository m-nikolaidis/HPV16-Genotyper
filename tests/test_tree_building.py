import pathlib
import tempfile
import textwrap
import unittest
from unittest.mock import patch

from Bio import Phylo

from hpv16genotyper import appFunctions


class BioNJTreeBuildingTests(unittest.TestCase):
    def _write_fake_fastme(self, root: pathlib.Path, body: str) -> pathlib.Path:
        executable = root / "fastme"
        executable.write_text(
            "#!/usr/bin/env python3\n" + textwrap.dedent(body), encoding="utf-8"
        )
        executable.chmod(0o755)
        return executable

    def test_build_trees_runs_bionj_and_preserves_original_sequence_ids(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            trees_dir = root / "Phylogenetic_Trees"
            trees_dir.mkdir()
            alignment = root / "sample_E1_aln.fa"
            alignment.write_text(
                ">sample_accession_with_a_long_name\nACGT\n"
                ">A1_reference_with_a_long_name\nACGA\n",
                encoding="utf-8",
            )
            fastme = self._write_fake_fastme(
                root,
                """
                import pathlib
                import sys

                arguments = sys.argv[1:]
                required = {
                    "--method=I",
                    "--dna=4",
                    "--branch_length=n",
                    "--nb_threads=1",
                }
                if not required.issubset(arguments):
                    sys.stderr.write("missing BioNJ options")
                    raise SystemExit(2)

                input_path = pathlib.Path(
                    next(value.split("=", 1)[1] for value in arguments
                         if value.startswith("--input_data="))
                )
                output_path = pathlib.Path(
                    next(value.split("=", 1)[1] for value in arguments
                         if value.startswith("--output_tree="))
                )
                lines = input_path.read_text(encoding="utf-8").splitlines()
                if lines[0].split() != ["2", "4"]:
                    sys.stderr.write("input was not PHYLIP")
                    raise SystemExit(3)
                identifiers = [line.split()[0] for line in lines[1:] if line.strip()]
                if len(identifiers) != 2 or any(len(name) > 10 for name in identifiers):
                    sys.stderr.write("unsafe PHYLIP identifiers")
                    raise SystemExit(4)
                output_path.write_text(
                    f"({identifiers[0]}:0.1,{identifiers[1]}:0.2);\\n",
                    encoding="utf-8",
                )
                """,
            )

            functions = appFunctions.MainFunctions()
            result = functions.build_trees(root, [alignment], str(fastme), threads=1)

            tree_path = trees_dir / "sample_E1_aln_NJ_tree.nwk"
            self.assertEqual(result, trees_dir)
            tree = Phylo.read(tree_path, "newick")
            self.assertEqual(
                {terminal.name for terminal in tree.get_terminals()},
                {
                    "sample_accession_with_a_long_name",
                    "A1_reference_with_a_long_name",
                },
            )

    def test_build_trees_reports_fastme_failures(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            (root / "Phylogenetic_Trees").mkdir()
            alignment = root / "sample_E1_aln.fa"
            alignment.write_text(">sample\nACGT\n>A1_E1\nACGA\n", encoding="utf-8")
            fastme = self._write_fake_fastme(
                root,
                """
                import sys
                sys.stderr.write("cannot calculate distances")
                raise SystemExit(9)
                """,
            )

            with self.assertRaises(appFunctions.ExternalToolError) as raised:
                appFunctions.MainFunctions().build_trees(
                    root, [alignment], str(fastme), threads=1
                )

            self.assertIn("FastME", str(raised.exception))
            self.assertIn("cannot calculate distances", str(raised.exception))

    def test_build_trees_reports_when_fastme_produces_no_tree(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            (root / "Phylogenetic_Trees").mkdir()
            alignment = root / "sample_E1_aln.fa"
            alignment.write_text(">sample\nACGT\n>A1_E1\nACGA\n", encoding="utf-8")
            fastme = self._write_fake_fastme(root, "raise SystemExit(0)\n")

            with self.assertRaises(appFunctions.ExternalToolError) as raised:
                appFunctions.MainFunctions().build_trees(
                    root, [alignment], str(fastme), threads=1
                )

            self.assertIn("FastME", str(raised.exception))
            self.assertIn("did not produce a tree", str(raised.exception))

    def test_binary_discovery_requires_fastme(self):
        with patch.object(
            appFunctions.shutil, "which", side_effect=lambda name: f"/tools/{name}"
        ):
            binaries = appFunctions._init_binaries("linux")

        self.assertEqual(
            binaries,
            [
                "/tools/makeblastdb",
                "/tools/blastn",
                "/tools/muscle",
                "/tools/fastme",
            ],
        )


if __name__ == "__main__":
    unittest.main()
