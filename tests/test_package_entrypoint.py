import importlib
from pathlib import Path
import unittest


class PackageEntrypointTests(unittest.TestCase):
    def test_package_exposes_module_entrypoint(self):
        package = importlib.import_module("hpv16genotyper")
        entrypoint = Path(package.__file__).with_name("__main__.py")

        self.assertTrue(entrypoint.is_file())
        self.assertIn("from .app import main", entrypoint.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
