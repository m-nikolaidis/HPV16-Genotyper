#!/usr/bin/env bash

set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
project_root="$(cd "$script_dir/../.." && pwd -P)"
output_dir="${1:-$project_root/downloads/windows}"
cache_dir="$output_dir/cache"
bin_dir="$output_dir/bin"

blast_version="2.17.0"
fastme_revision="08e2150495640a45b89fb816edd894bec8219186"

blast_url="${HPV16_BLAST_URL:-https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/$blast_version/ncbi-blast-$blast_version+-x64-win64.tar.gz}"
blast_sha256="${HPV16_BLAST_SHA256:-ccde8788641e8f4137536aaadedfeac2f3599dbbc6166e701b5d89d19fa79038}"
muscle_url="${HPV16_MUSCLE_URL:-https://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86win32.exe}"
muscle_sha256="${HPV16_MUSCLE_SHA256:-e6ba0fe75a4b12f401be632653314fcffc5aca3dec98c9cdd20cd1ed1ce1dc6a}"
fastme_url="${HPV16_FASTME_URL:-https://gite.lirmm.fr/atgc/FastME/-/archive/$fastme_revision/FastME-$fastme_revision.tar.gz}"
fastme_sha256="${HPV16_FASTME_SHA256:-26b7ded89ff6943134ccc9306088d7d8c0d8bb7212e60b11630d73f38f453e34}"

mingw_host="${HPV16_MINGW_HOST:-x86_64-w64-mingw32}"
mingw_cc="${HPV16_MINGW_CC:-$mingw_host-gcc}"
mingw_objdump="${HPV16_MINGW_OBJDUMP:-$mingw_host-objdump}"
jobs="${HPV16_BUILD_JOBS:-2}"

require_command() {
    if ! command -v "$1" >/dev/null 2>&1; then
        printf 'Required command not found: %s\n' "$1" >&2
        exit 1
    fi
}

download() {
    local url="$1"
    local destination="$2"
    local expected_sha256="$3"

    if [[ -f "$destination" ]] && printf '%s  %s\n' "$expected_sha256" "$destination" | sha256sum --check --status; then
        printf 'Using cached %s\n' "$(basename "$destination")"
        return
    fi

    local partial="$destination.part"
    curl --fail --location --retry 3 --output "$partial" "$url"
    if ! printf '%s  %s\n' "$expected_sha256" "$partial" | sha256sum --check --status; then
        printf 'Checksum verification failed for %s\n' "$url" >&2
        rm -f "$partial"
        exit 1
    fi
    mv "$partial" "$destination"
}

for command in curl tar sha256sum autoreconf find grep make install "$mingw_cc" "$mingw_objdump"; do
    require_command "$command"
done

mkdir -p "$cache_dir"
build_dir="$(mktemp -d "$output_dir/.build.XXXXXX")"
previous_root="$(mktemp -d "$output_dir/.previous-bin.XXXXXX")"
had_previous=0
published=0

cleanup() {
    local status=$?
    trap - EXIT HUP INT TERM
    set +e
    if (( ! published && had_previous )) \
        && [[ ! -e "$bin_dir" && -e "$previous_root/bin" ]]; then
        mv "$previous_root/bin" "$bin_dir"
        if (( $? != 0 && status == 0 )); then
            status=1
        fi
    fi
    rm -rf "$build_dir" "$previous_root"
    exit "$status"
}

trap cleanup EXIT
trap 'exit 129' HUP
trap 'exit 130' INT
trap 'exit 143' TERM
staged_bin="$build_dir/bin"
mkdir -p "$staged_bin"

blast_archive="$cache_dir/ncbi-blast-$blast_version+-x64-win64.tar.gz"
muscle_download="$cache_dir/muscle3.8.31_i86win32.exe"
fastme_archive="$cache_dir/FastME-$fastme_revision.tar.gz"

download "$blast_url" "$blast_archive" "$blast_sha256"
download "$muscle_url" "$muscle_download" "$muscle_sha256"
download "$fastme_url" "$fastme_archive" "$fastme_sha256"

mkdir -p "$build_dir/blast" "$build_dir/fastme"
tar -xzf "$blast_archive" -C "$build_dir/blast" --strip-components=1
install -m 0755 "$build_dir/blast/bin/blastn.exe" "$staged_bin/blastn.exe"
install -m 0755 "$build_dir/blast/bin/makeblastdb.exe" "$staged_bin/makeblastdb.exe"
while IFS= read -r -d '' dll; do
    install -m 0644 "$dll" "$staged_bin/$(basename "$dll")"
done < <(find "$build_dir/blast/bin" -maxdepth 1 -type f -iname '*.dll' -print0)

install -m 0755 "$muscle_download" "$staged_bin/muscle.exe"

tar -xzf "$fastme_archive" -C "$build_dir/fastme"
fastme_source="$(find "$build_dir/fastme" -mindepth 1 -maxdepth 1 -type d -print -quit)"
if [[ -z "$fastme_source" || ! -x "$fastme_source/configure" ]]; then
    printf 'FastME source archive does not contain an executable configure script\n' >&2
    exit 1
fi

(
    cd "$fastme_source"
    # GitLab archives preserve maintainer-local Automake symlinks, so replace
    # them with portable helper files before configuring the cross-build.
    autoreconf --force --install
    CC="$mingw_cc" LDFLAGS="-static" ./configure \
        --host="$mingw_host"
    make -j "$jobs"
)
install -m 0755 "$fastme_source/src/fastme.exe" "$staged_bin/fastme.exe"
if ! "$mingw_objdump" -f "$staged_bin/fastme.exe" | grep --quiet 'architecture: i386:x86-64'; then
    printf 'FastME build is not an x86-64 Windows executable\n' >&2
    exit 1
fi
if "$mingw_objdump" -p "$staged_bin/fastme.exe" \
    | grep --extended-regexp --ignore-case --quiet \
        'DLL Name: (libgomp|libgcc|libwinpthread|libssp|libstdc\+\+)'; then
    printf 'FastME build has an unexpected MinGW runtime DLL dependency\n' >&2
    exit 1
fi

(
    cd "$staged_bin"
    sha256sum ./*.exe ./*.dll > SHA256SUMS
)

if [[ -e "$bin_dir" ]]; then
    had_previous=1
    mv "$bin_dir" "$previous_root/bin"
fi
if ! mv "$staged_bin" "$bin_dir"; then
    exit 1
fi
published=1

printf 'Windows tools are ready in %s\n' "$bin_dir"
