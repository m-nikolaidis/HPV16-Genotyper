import hashlib
import os
import pathlib
import subprocess
import tarfile
import tempfile
import unittest


PROJECT_ROOT = pathlib.Path(__file__).resolve().parents[1]
DOWNLOAD_SCRIPT = (
    PROJECT_ROOT / "hpv16genotyper" / "misc" / "download_windows_binaries.sh"
)


def _sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


class WindowsBinaryDownloadTests(unittest.TestCase):
    def test_prepares_packaging_ready_windows_tool_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            fixtures = root / "fixtures"
            fixtures.mkdir()

            blast_tree = root / "blast-tree" / "ncbi-blast-test" / "bin"
            blast_tree.mkdir(parents=True)
            (blast_tree / "blastn.exe").write_bytes(b"BLASTN")
            (blast_tree / "makeblastdb.exe").write_bytes(b"MAKEBLASTDB")
            (blast_tree / "ncbi-runtime.dll").write_bytes(b"BLAST-DLL")
            blast_archive = fixtures / "blast.tar.gz"
            with tarfile.open(blast_archive, "w:gz") as archive:
                archive.add(blast_tree.parent, arcname="ncbi-blast-test")

            muscle_binary = fixtures / "muscle-upstream.exe"
            muscle_binary.write_bytes(b"MUSCLE")

            fastme_tree = root / "fastme-tree" / "FastME-test"
            fastme_tree.mkdir(parents=True)
            (fastme_tree / "marker").write_bytes(b"FASTME")
            configure = fastme_tree / "configure"
            configure.write_text(
                """#!/usr/bin/env bash
set -euo pipefail
[[ -f install-sh ]]
[[ " $* " == *" --host=x86_64-w64-mingw32 "* ]]
[[ " $* " != *" --disable-OpenMP "* ]]
[[ "${CC##*/}" == "x86_64-w64-mingw32-gcc" ]]
[[ "${LDFLAGS:-}" == *"-static"* ]]
mkdir -p src
printf '%s\n' 'all:' $'\tcp marker src/fastme.exe' > Makefile
""",
                encoding="utf-8",
            )
            configure.chmod(0o755)
            (fastme_tree / "install-sh").symlink_to("/missing/automake/install-sh")
            fastme_archive = fixtures / "fastme.tar.gz"
            with tarfile.open(fastme_archive, "w:gz") as archive:
                archive.add(fastme_tree, arcname="FastME-test")

            fixture_bin = root / "fixture-bin"
            fixture_bin.mkdir()
            autoreconf = fixture_bin / "autoreconf"
            autoreconf.write_text(
                """#!/usr/bin/env bash
set -euo pipefail
rm install-sh
printf '#!/usr/bin/env sh\n' > install-sh
chmod +x install-sh
""",
                encoding="utf-8",
            )
            autoreconf.chmod(0o755)
            objdump = fixture_bin / "x86_64-w64-mingw32-objdump"
            objdump.write_text(
                """#!/usr/bin/env bash
set -euo pipefail
printf '%s\n' "$*" >> "$HPV16_OBJDUMP_MARKER"
if [[ "$1" == "-f" ]]; then
    printf '%s\n' 'architecture: i386:x86-64, flags 0x0000012f:'
elif [[ "$1" == "-p" ]]; then
    printf '%s\n' 'DLL Name: KERNEL32.dll'
else
    exit 2
fi
""",
                encoding="utf-8",
            )
            objdump.chmod(0o755)

            output = root / "windows-tools"
            existing_bin = output / "bin"
            existing_bin.mkdir(parents=True)
            (existing_bin / "stale-runtime.dll").write_bytes(b"STALE")
            environment = os.environ.copy()
            environment.update(
                {
                    "HPV16_BLAST_URL": blast_archive.as_uri(),
                    "HPV16_BLAST_SHA256": _sha256(blast_archive),
                    "HPV16_MUSCLE_URL": muscle_binary.as_uri(),
                    "HPV16_MUSCLE_SHA256": _sha256(muscle_binary),
                    "HPV16_FASTME_URL": fastme_archive.as_uri(),
                    "HPV16_FASTME_SHA256": _sha256(fastme_archive),
                    "HPV16_OBJDUMP_MARKER": str(root / "objdump-called"),
                    "PATH": f"{fixture_bin}{os.pathsep}{environment['PATH']}",
                }
            )

            if not DOWNLOAD_SCRIPT.is_file():
                self.fail(f"download script is missing: {DOWNLOAD_SCRIPT}")
            completed = subprocess.run(
                ["bash", str(DOWNLOAD_SCRIPT), str(output)],
                cwd=PROJECT_ROOT,
                env=environment,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)

            binaries = output / "bin"
            self.assertEqual((binaries / "blastn.exe").read_bytes(), b"BLASTN")
            self.assertEqual(
                (binaries / "makeblastdb.exe").read_bytes(), b"MAKEBLASTDB"
            )
            self.assertEqual(
                (binaries / "ncbi-runtime.dll").read_bytes(), b"BLAST-DLL"
            )
            self.assertEqual((binaries / "muscle.exe").read_bytes(), b"MUSCLE")
            self.assertEqual((binaries / "fastme.exe").read_bytes(), b"FASTME")
            self.assertFalse((binaries / "stale-runtime.dll").exists())

            manifest = (binaries / "SHA256SUMS").read_text(encoding="utf-8")
            self.assertIn("./fastme.exe", manifest)
            verified = subprocess.run(
                ["sha256sum", "--check", "SHA256SUMS"],
                cwd=binaries,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertEqual(verified.returncode, 0, verified.stderr)
            objdump_calls = (root / "objdump-called").read_text(encoding="utf-8")
            self.assertIn("-f", objdump_calls)
            self.assertIn("-p", objdump_calls)

            previous_manifest = (binaries / "SHA256SUMS").read_bytes()
            environment["HPV16_FASTME_SHA256"] = "0" * 64
            failed = subprocess.run(
                ["bash", str(DOWNLOAD_SCRIPT), str(output)],
                cwd=PROJECT_ROOT,
                env=environment,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertNotEqual(failed.returncode, 0)
            self.assertIn("Checksum verification failed", failed.stderr)
            self.assertEqual(
                (binaries / "SHA256SUMS").read_bytes(), previous_manifest
            )
            self.assertEqual((binaries / "fastme.exe").read_bytes(), b"FASTME")

            interrupt_bin = root / "interrupt-bin"
            interrupt_bin.mkdir()
            interrupting_mv = interrupt_bin / "mv"
            interrupting_mv.write_text(
                """#!/usr/bin/env bash
set -euo pipefail
if [[ "$1" == "$HPV16_BIN_TO_INTERRUPT" ]]; then
    /usr/bin/mv "$@"
    kill -TERM "$PPID"
    sleep 1
    exit 143
fi
exec /usr/bin/mv "$@"
""",
                encoding="utf-8",
            )
            interrupting_mv.chmod(0o755)
            environment["HPV16_FASTME_SHA256"] = _sha256(fastme_archive)
            environment["HPV16_BIN_TO_INTERRUPT"] = str(binaries)
            environment["PATH"] = (
                f"{interrupt_bin}{os.pathsep}{environment['PATH']}"
            )
            interrupted = subprocess.run(
                ["bash", str(DOWNLOAD_SCRIPT), str(output)],
                cwd=PROJECT_ROOT,
                env=environment,
                text=True,
                capture_output=True,
                check=False,
            )
            self.assertNotEqual(interrupted.returncode, 0)
            self.assertEqual(
                (binaries / "SHA256SUMS").read_bytes(), previous_manifest
            )
            self.assertEqual((binaries / "fastme.exe").read_bytes(), b"FASTME")


if __name__ == "__main__":
    unittest.main()
