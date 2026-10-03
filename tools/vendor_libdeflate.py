#!/usr/bin/env python3
"""Copy the pinned libdeflate release into `third_party/`, and only that release.

The core takes one third-party dependency -- libdeflate, the inflate behind the
gz path -- and for a long time it took it from the host: `find_library(deflate)`
plus an optional static archive that Debian does not build PIC.  That makes the
wheel's portability a property of the machine that built it, which is the one
failure mode a binary distribution cannot have.  So the release is vendored
instead, and this script is what puts it there.

Vendoring by hand is how a tree ends up holding eleven files from one version
and a header from another, or holding a source file that CMake does not
compile, so this script checks four things and **writes nothing unless all four
hold**:

1. The tarball's SHA-256 equals the pin below.  GitHub regenerates its
   `/archive/` tarballs from time to time and the bytes then differ for the same
   tag; a mismatch stops here rather than vendoring a different set of bytes
   than the ones this pin was reviewed against.  See `--tarball` for the offline
   path and `third_party/libdeflate/README.md` for what to do when it fires.
2. The version in the extracted `libdeflate.h` equals the pinned one, so a
   mislabelled download cannot be vendored under the wrong name.
3. The list of C sources is read out of *upstream's own* `CMakeLists.txt`, not
   transcribed here, and every one of those files is present in the extracted
   tree.  The list is written to `libdeflate_sources.cmake`, which the build
   includes, so bumping the version is one script run and not an edit to two
   files that can drift apart.
4. After writing, the destination tree is walked again and every file's digest
   compared with the source's, and the file set is compared with the one
   intended -- no stale file from an older version survives, and no file
   outside the recorded set is left behind.

The result is reproducible from this command:

    python3 tools/vendor_libdeflate.py             # fetch, verify, vendor
    python3 tools/vendor_libdeflate.py --check     # verify the tree, write nothing

The second form is what a test can run: it re-downloads nothing and compares the
vendored tree against what the pin says it should be.

This is a build-time tool.  It uses the standard library plus `tarfile` and has
no runtime counterpart; nothing here ships.
"""

from __future__ import annotations

import argparse
import hashlib
import re
import shutil
import sys
import tarfile
import tempfile
import urllib.request
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DESTINATION = ROOT / "third_party" / "libdeflate"

# Disposable, and under build/ for the same reason the bench cache is: it is
# derived from a URL and a digest, so losing it costs one download.
CACHE = ROOT / "build" / "vendor-cache"

NAME = "libdeflate"
TAG = "v1.23"
VERSION = "1.23"
URL = f"https://github.com/ebiggers/libdeflate/archive/refs/tags/{TAG}.tar.gz"
SHA256 = "1ab18349b9fb0ce8a0ca4116bded725be7dcbfa709e19f6f983d99df1fb8b25f"
LICENSE = "MIT"

# What travels and what does not.  `lib/**` is the library entire -- the eleven
# C files and the headers they include -- because upstream's sources include
# their neighbours by name and a hand-picked subset is a build that breaks when
# a version adds one.  `common_defs.h` is included from `lib/lib_common.h` as
# `"../common_defs.h"`, so it belongs at the root of the vendored tree with
# `libdeflate.h`, which is the public header and the only file a consumer
# includes.  Everything else in the release -- the CMake build, the gzip
# program, the benchmark and test scripts, the CI configuration -- is not
# vendored: we build the library with our own CMakeLists, so carrying a second
# build system would be carrying a second thing to keep working.
TAKEN_FILES = ("libdeflate.h", "common_defs.h", "COPYING")
TAKEN_DIRECTORIES = ("lib",)

# Written by this script next to the sources.  Named here so that the
# destination walk in check 4 knows what it put there itself.
GENERATED = ("README.md", "libdeflate_sources.cmake")

#: Files a consumer links against nothing of, but that upstream ships in `lib/`
#: and a `.c` includes.  Recorded rather than filtered, so the count in the
#: provenance note is a count of what is on disk.
SOURCE_SUFFIX = ".c"


def sha256_of(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while chunk := handle.read(1 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def fetch(tarball: Path | None) -> Path:
    """Return a tarball whose digest is the pin, downloading it if needed."""
    if tarball is not None:
        tarball = Path(tarball)
        got = sha256_of(tarball)
        if got != SHA256:
            raise SystemExit(
                f"{tarball} has SHA-256 {got}, not the pinned {SHA256}.\n"
                "Refusing to vendor bytes that are not the reviewed ones."
            )
        return tarball

    CACHE.mkdir(parents=True, exist_ok=True)
    cached = CACHE / f"{NAME}-{VERSION}.tar.gz"
    if cached.exists() and sha256_of(cached) == SHA256:
        return cached

    print(f"downloading {URL}")
    partial = cached.with_suffix(cached.suffix + ".part")
    with urllib.request.urlopen(URL) as response, open(partial, "wb") as out:
        shutil.copyfileobj(response, out)
    got = sha256_of(partial)
    if got != SHA256:
        partial.unlink()
        raise SystemExit(
            f"{URL} has SHA-256 {got}, not the pinned {SHA256}.\n"
            "The tag may have been re-cut, or GitHub regenerated the archive.\n"
            "Do not simply update this pin: download the tarball, diff it\n"
            "against the vendored tree, and change the pin deliberately."
        )
    partial.replace(cached)
    return cached


def extract(tarball: Path, into: Path) -> Path:
    """Extract `tarball` under `into` and return the release's root."""
    with tarfile.open(tarball, "r:gz") as archive:
        for member in archive.getmembers():
            # A tarball is untrusted input even when its digest matches: a
            # member named `../../x` or a symlink out of the tree is a write
            # outside the build directory.  Both are one line to refuse.
            if not member.isreg() and not member.isdir():
                raise SystemExit(f"{member.name}: not a regular file or directory")
            target = (into / member.name).resolve()
            if not target.is_relative_to(into.resolve()):
                raise SystemExit(f"{member.name}: escapes the extraction directory")
        # 3.12 added the extraction filters and 3.14 will default to `data`,
        # which is the checks above plus permissions and ownership.  Passing it
        # where it exists is the same decision applied twice; the project floor
        # is 3.10, so the manual checks are what covers the versions before it.
        filters = {"filter": "data"} if sys.version_info >= (3, 12) else {}
        archive.extractall(into, **filters)

    roots = [p for p in into.iterdir() if p.is_dir()]
    if len(roots) != 1:
        raise SystemExit(f"expected one directory in the archive, found {roots}")
    return roots[0]


def upstream_sources(release: Path) -> list[str]:
    """The C files upstream's CMakeLists compiles, in its own order.

    Reading the list out of upstream rather than writing it here is the point:
    a hand-copied list is the thing that goes stale, and the failure it causes
    is silent -- a new source file exists in the tree and is never compiled, so
    the library builds and behaves like the version before it.
    """
    text = (release / "CMakeLists.txt").read_text(encoding="utf-8")
    names: list[str] = []
    for block in re.findall(r"(?:set\(|list\(APPEND\s+)LIB_SOURCES(.*?)\)", text, re.S):
        names += re.findall(r"([\w./-]+\.c)\b", block)
    if not names:
        raise SystemExit("no LIB_SOURCES entries found in upstream's CMakeLists")
    return names


def vendored_layout(release: Path) -> dict[str, Path]:
    """Every file this script intends to write, as `relative path -> source`."""
    layout: dict[str, Path] = {}
    for name in TAKEN_FILES:
        source = release / name
        if not source.is_file():
            raise SystemExit(f"{name} is missing from the release")
        layout[name] = source
    for name in TAKEN_DIRECTORIES:
        directory = release / name
        if not directory.is_dir():
            raise SystemExit(f"{name}/ is missing from the release")
        for source in sorted(directory.rglob("*")):
            if source.is_file():
                layout[str(source.relative_to(release))] = source
    return layout


def sources_cmake(sources: list[str]) -> str:
    lines = [
        "# Generated by tools/vendor_libdeflate.py -- do not edit by hand.",
        "#",
        "# Upstream's own LIB_SOURCES list for this release, verbatim in its order.",
        "# The build includes this file rather than repeating the list, so a",
        "# version bump is one script run and not an edit in two places.",
        "",
        "set(BIOFASTING_LIBDEFLATE_SOURCES",
    ]
    lines += [f"  {name}" for name in sources]
    lines.append(")")
    lines.append("")
    return "\n".join(lines)


def provenance(release: Path, sources: list[str], layout: dict[str, Path]) -> str:
    # The release's own entries, split into the ones that travel (as a
    # top-level name: `lib/arm/x.c` means `lib` travels) and the ones that do
    # not.  Computed from the layout rather than listed here, so a release that
    # adds a directory cannot be described wrongly.
    roots = {name.split("/")[0] for name in layout}
    dropped = "\n".join(
        f"- `{path.name}`"
        for path in sorted(release.iterdir())
        if path.name not in roots
    )
    taken = "\n".join(f"- `{name}`" for name in sorted(layout) if "/" not in name)
    return f"""# libdeflate {VERSION}, vendored

Not written by hand: `python3 tools/vendor_libdeflate.py` put these files here
and re-running it is how they are updated.  Read this before editing anything
in the tree.

## Where it came from

| | |
|---|---|
| Upstream | <https://github.com/ebiggers/libdeflate> |
| Release | tag `{TAG}` |
| Tarball | <{URL}> |
| SHA-256 | `{SHA256}` |
| Version in `libdeflate.h` | `{VERSION}` |
| Licence | {LICENSE}, in [COPYING](COPYING) |

The digest is of the tarball, checked before anything is extracted, and the
version in `libdeflate.h` is checked against the tag afterwards.  GitHub
regenerates its `/archive/` tarballs occasionally, and then the same tag
downloads as different bytes; the script stops there instead of vendoring
something unreviewed.  The fix is not to update the pin -- it is to download the
tarball, diff it against this tree, and change the pin deliberately.

## Why it is vendored at all

The core needs one inflate, and taking it from the host made the wheel's
portability a property of the machine that built it: a shared library links fine
on the builder and fails to import on a user's machine, and the static archives
Debian ships are not built with `-fPIC`, so they cannot be linked into a shared
module at all.  Vendoring the release and compiling it here is what makes the
wheel self-contained, and `build_info()["libdeflate"]` says `{VERSION} (vendored,
static)` so a gz number can be traced to the binary that produced it.

## What is here

Taken, the library whole:

{taken}
- `lib/**` -- {len([n for n in layout if n.startswith('lib/')])} files: the {len(sources)} C sources upstream's own `CMakeLists.txt` lists
  in `LIB_SOURCES`, plus every header they include.  The list travels in
  [libdeflate_sources.cmake](libdeflate_sources.cmake), which the build
  includes, so the two cannot drift apart.

Left behind, everything else in the release:

{dropped}

None of it is needed: the library is built by this project's `CMakeLists.txt`
with upstream's own compile settings (`-O2 -DNDEBUG`, C99), so carrying a second
build system would only be a second thing to keep working.  The gzip and
benchmark programs, the test scripts and the CI configuration are not part of
the library.

## Updating to a new release

1. Change `TAG`, `VERSION` and `SHA256` in `tools/vendor_libdeflate.py`.
2. `python3 tools/vendor_libdeflate.py` -- it refuses to write until the digest,
   the version and the source list all agree.
3. Re-measure the gz rows in `bench/targets.md`.  A different libdeflate is a
   different inflate, and a number quoted without it is a number about an
   unknown program.

## Compile settings worth knowing

Upstream pins `-O2 -DNDEBUG` for Release and C99, and that is what the vendored
target uses rather than this project's `-O3`: the point of vendoring is the same
code compiled the same way, and an unmeasured `-O3` here would make every gz
benchmark row incomparable with the ones measured before it.

Upstream also probes the assembler for AVX512VNNI, VPCLMULQDQ, AVX-VNNI, DOTPROD
and SHA3 support and defines `LIBDEFLATE_ASSEMBLER_DOES_NOT_SUPPORT_*` when a
compiler is paired with a binutils older than itself.  That probe is carried
into this project's `CMakeLists.txt` unchanged, so the vendored build fails
where upstream's fails instead of introducing a portability regression.
"""


def write_tree(layout: dict[str, Path], sources: list[str], release: Path) -> None:
    """Replace the destination with exactly `layout`, plus the two generated files."""
    if DESTINATION.exists():
        # Replaced rather than merged.  A merge leaves the previous release's
        # headers next to this one's, which builds and is wrong.
        shutil.rmtree(DESTINATION)
    for relative in layout:
        target = DESTINATION / relative
        target.parent.mkdir(parents=True, exist_ok=True)
    for relative, source in layout.items():
        shutil.copyfile(source, DESTINATION / relative)
    (DESTINATION / "README.md").write_text(
        provenance(release, sources, layout), encoding="utf-8"
    )
    (DESTINATION / "libdeflate_sources.cmake").write_text(
        sources_cmake(sources), encoding="utf-8"
    )


def verify(layout: dict[str, Path], sources: list[str]) -> list[str]:
    """Check the destination against `layout`.  Returns the problems found.

    Run after writing and, with `--check`, instead of writing -- so the same
    code answers "did the write work" and "is the tree still what the pin says",
    and a tree edited by hand is caught by the same comparison as a bad copy.
    """
    problems: list[str] = []
    expected = dict(layout)
    expected["README.md"] = None
    expected["libdeflate_sources.cmake"] = None

    if not DESTINATION.is_dir():
        return [f"{DESTINATION} does not exist"]

    present = {
        str(path.relative_to(DESTINATION))
        for path in DESTINATION.rglob("*")
        if path.is_file()
    }
    for missing in sorted(set(expected) - present):
        problems.append(f"missing: {missing}")
    for extra in sorted(present - set(expected)):
        problems.append(f"not part of release {TAG}: {extra}")

    for relative, source in layout.items():
        target = DESTINATION / relative
        if target.is_file() and sha256_of(target) != sha256_of(source):
            problems.append(f"differs from the tarball: {relative}")

    listed = re.findall(
        r"^\s+([\w./-]+\.c)$",
        (DESTINATION / "libdeflate_sources.cmake").read_text(encoding="utf-8"),
        re.M,
    )
    if listed != sources:
        problems.append(
            "libdeflate_sources.cmake does not match upstream's LIB_SOURCES "
            f"({len(listed)} entries against {len(sources)})"
        )
    on_disk = sorted(
        str(path.relative_to(DESTINATION))
        for path in DESTINATION.rglob(f"*{SOURCE_SUFFIX}")
    )
    for orphan in sorted(set(on_disk) - set(sources)):
        problems.append(f"vendored but never compiled: {orphan}")
    return problems


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--tarball",
        type=Path,
        help="a local copy of the release tarball, checked against the pin",
    )
    parser.add_argument(
        "--check",
        action="store_true",
        help="verify the vendored tree against the pin without writing it",
    )
    args = parser.parse_args(argv)

    tarball = fetch(args.tarball)
    with tempfile.TemporaryDirectory() as scratch:
        release = extract(tarball, Path(scratch))

        header = (release / "libdeflate.h").read_text(encoding="utf-8")
        found = re.search(r'#define LIBDEFLATE_VERSION_STRING\s+"([^"]+)"', header)
        if not found or found.group(1) != VERSION:
            raise SystemExit(
                f"the tarball says version {found and found.group(1)}, "
                f"the pin says {VERSION}"
            )

        sources = upstream_sources(release)
        layout = vendored_layout(release)
        for relative in sources:
            if relative not in layout:
                raise SystemExit(
                    f"upstream compiles {relative}, which is not in the vendored set"
                )

        if args.check:
            problems = verify(layout, sources)
            if problems:
                print(f"{DESTINATION} is not release {TAG}:")
                for problem in problems:
                    print(f"  {problem}")
                print("\nRun `python3 tools/vendor_libdeflate.py` to restore it.")
                return 1
            print(
                f"{DESTINATION} is libdeflate {VERSION}, "
                f"{len(layout)} files, {len(sources)} C sources"
            )
            return 0

        write_tree(layout, sources, release)
        problems = verify(layout, sources)
        if problems:
            for problem in problems:
                print(f"  {problem}", file=sys.stderr)
            raise SystemExit("the vendored tree does not match what was written")

    print(
        f"vendored libdeflate {VERSION} into {DESTINATION.relative_to(ROOT)}: "
        f"{len(layout)} files, {len(sources)} C sources"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
