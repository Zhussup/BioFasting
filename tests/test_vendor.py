"""The wheel carries its own inflate, and this is the test that says so.

The core has one third-party dependency, libdeflate, and for a while it took it
from the host.  That made the wheel's portability a property of the machine that
built it: an extension linked against a shared `libdeflate.so` imports on the
builder and fails on a user's machine, and it does so at *import* time, before
any of the code that would report the problem.  CI could not see it either,
because CI installed libdeflate too.  The fix was to vendor the release
(`third_party/libdeflate`, put there by `tools/vendor_libdeflate.py`) and compile
it into the wheel.

So the claim under test is not "the gz path is fast".  It is:

1. the extension has no dynamic dependency on libdeflate -- asserted on the
   built artifact rather than inferred from the build system, because the two
   can disagree when a linker flag is dropped;
2. the vendored tree is the version the build reports, so a number that names an
   inflate names the inflate that produced it;
3. every `.c` in the vendored tree is compiled and every compiled file is there,
   so a release bump cannot add a source that nothing builds;
4. the tree is what the pin says it is -- the same check
   `tools/vendor_libdeflate.py --check` runs, so a hand edit to a third-party
   file fails here instead of in a benchmark nobody re-runs.

A build configured with `-DBIOFASTING_VENDORED_LIBDEFLATE=OFF` links the host's
libdeflate on purpose.  Tests 1 and 2 skip there and say why -- asking that
build to prove the wheel is self-contained would be asking it to fail; tests 3,
4 and 5 hold either way, because they are about the tree and not about the
linkage.
"""

import re
import subprocess
import sys
from pathlib import Path

import pytest

import biofasting

VENDORED = Path(__file__).resolve().parent.parent / "third_party" / "libdeflate"

#: What a failure to *fetch* the pinned tarball looks like, as opposed to a
#: failure of the check itself.  A verification failure is printed to stdout by
#: the script and is never skipped; these strings only ever appear in the
#: traceback of a network or filesystem error, which says nothing about this
#: tree and must not be reported as one.
UNREACHABLE = (
    "urlopen",
    "URLError",
    "Temporary failure",
    "timed out",
    "Connection",
    "SSL",
)


def read_version(header: Path) -> str | None:
    """The `LIBDEFLATE_VERSION_STRING` a header declares, or None."""
    if not header.is_file():
        return None
    found = re.search(
        r'#define LIBDEFLATE_VERSION_STRING\s+"([^"]+)"',
        header.read_text(encoding="utf-8"),
    )
    return found.group(1) if found else None


def dynamic_dependencies(module: Path) -> list[str] | None:
    """The shared libraries `module` names, or None where no tool can say."""
    if sys.platform.startswith("linux"):
        command = ["readelf", "-d", str(module)]
        pattern = r"\(NEEDED\)\s+Shared library: \[([^\]]+)\]"
    elif sys.platform == "darwin":
        command = ["otool", "-L", str(module)]
        pattern = r"^\s+(\S+\.dylib)"
    else:
        return None
    try:
        output = subprocess.run(
            command, capture_output=True, text=True, check=True
        ).stdout
    except (OSError, subprocess.CalledProcessError):
        return None
    return re.findall(pattern, output, re.M)


def vendored_by_configuration() -> bool:
    """Whether the built extension says it was built from `third_party/`.

    The two artifact tests below turn on this, and the reason is not to make
    them pass: `-DBIOFASTING_VENDORED_LIBDEFLATE=OFF` is a supported way to build
    against a system libdeflate, and asking that build to prove the wheel is
    self-contained would be asking it to fail.  What the gate makes checkable is
    the claim itself -- *if* the build says vendored, then the file that got
    installed had better have no libdeflate in its dynamic section.  A build that
    reports `vendored, static` and links `libdeflate.so.0` anyway is the defect,
    and it is not a hypothetical: it is what this package shipped before.
    """
    return "vendored" in biofasting.build_info()["libdeflate"]


def test_the_extension_does_not_depend_on_a_shared_libdeflate():
    """The whole point of vendoring, asserted on the file that gets installed.

    A `libdeflate.so.1` in this list means the wheel imports only where
    libdeflate is installed.  Read from the artifact rather than inferred from
    the build system, because the two disagree exactly when a linker flag is
    dropped, and that failure arrives at *import* time on the user's machine.
    """
    recorded = biofasting.build_info()["libdeflate"]
    if not vendored_by_configuration():
        pytest.skip(f"built against the host's libdeflate by configuration: {recorded}")
    module = Path(biofasting._core.__file__).resolve()
    needs = dynamic_dependencies(module)
    if needs is None:
        pytest.skip(f"no dynamic-section reader available on {sys.platform}")
    offenders = [name for name in needs if "libdeflate" in name]
    assert not offenders, (
        f"the build reports libdeflate {recorded!r} but {module.name} links "
        f"{offenders}, so the wheel is only importable where libdeflate is "
        "installed."
    )


def test_the_vendored_tree_is_the_version_the_build_reports():
    """Two places name the inflate; a bug report depends on them agreeing.

    Not a tautology, and not a check that "the string is non-empty": the version
    comes from the header the build actually compiled against, and this reads the
    header in the source tree.  A tree restored from a stale tarball while the
    build used a system header is exactly what this catches, and so is the
    reverse -- a vendored tree sitting there unused while the build reports an
    inflate nobody can point at.
    """
    recorded = biofasting.build_info()["libdeflate"]
    if not vendored_by_configuration():
        pytest.skip(f"built against the host's libdeflate by configuration: {recorded}")
    version = read_version(VENDORED / "libdeflate.h")
    assert version is not None, "the vendored tree has no libdeflate.h"
    assert recorded == f"{version} (vendored, static)", (
        f"build_info() reports libdeflate {recorded!r}, the vendored tree is "
        f"{version}"
    )


def test_the_vendored_sources_are_the_ones_the_build_compiles():
    """Every `.c` in the tree is compiled, and every compiled file is there.

    Both directions, because the silent one is the first: a source file that is
    vendored and never compiled leaves the library building, behaving exactly
    like the release before it, with nothing anywhere saying so.
    """
    if not VENDORED.is_dir():
        pytest.skip("third_party/libdeflate is absent")
    declared = re.findall(
        r"^\s+([\w./-]+\.c)$",
        (VENDORED / "libdeflate_sources.cmake").read_text(encoding="utf-8"),
        re.M,
    )
    present = sorted(str(path.relative_to(VENDORED)) for path in VENDORED.rglob("*.c"))
    assert declared, "libdeflate_sources.cmake lists no sources"
    assert sorted(declared) == present
    for relative in declared:
        assert (VENDORED / relative).is_file()


def test_the_vendored_tree_still_matches_the_release_it_came_from():
    """Run the vendoring script's own check, so an edit cannot go unnoticed.

    This is the only check that compares the tree against the release rather
    than against itself, and it needs the release to do it -- so it is skipped,
    with the reason, on a machine that cannot reach the pinned tarball.  What it
    catches when it does run is a hand edit to a third-party file, which is
    otherwise found by a benchmark nobody re-runs.
    """
    script = VENDORED.parent.parent / "tools" / "vendor_libdeflate.py"
    if not script.is_file():
        pytest.skip("tools/vendor_libdeflate.py is absent")
    result = subprocess.run(
        [sys.executable, str(script), "--check"], capture_output=True, text=True
    )
    if result.returncode != 0 and not result.stdout.strip():
        # The script prints a verification failure to stdout and returns 1; a
        # fetch failure is a traceback on stderr instead.  Only the second is
        # somebody else's problem.
        if not any(word in result.stderr for word in UNREACHABLE):
            pytest.fail(
                "the vendoring script failed in a way that is not a fetch "
                f"error:\n{result.stderr}"
            )
        pytest.skip("the pinned tarball could not be fetched from here")
    assert result.returncode == 0, result.stdout + result.stderr


def test_the_vendored_licence_travels_with_the_source():
    """MIT requires the notice to travel, and a wheel is a redistribution."""
    if not VENDORED.is_dir():
        pytest.skip("third_party/libdeflate is absent")
    notice = (VENDORED / "COPYING").read_text(encoding="utf-8")
    assert "Permission is hereby granted" in notice
    assert "Eric Biggers" in notice
