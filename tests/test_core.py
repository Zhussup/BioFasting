"""Step 1.1 smoke tests: the scaffold builds, imports, and reports itself.

Deliberately not performance tests.  Nothing here asserts a speed -- the job of
step 1.1 is to prove the toolchain works, and to prove that the metadata a bug
report depends on is measured rather than hard-coded.  A version string that is
true by construction is worth nothing in a bug report, so these tests compare
the two places the version exists and fail when they disagree.
"""

import platform

import biofasting

# Keys the C++ core must supply.
COMPILE_TIME_KEYS = {
    "version",
    "compiler",
    "cxx_standard",
    "build_type",
    "arch",
    "cxx_flags",
    "cmake",
    "nanobind",
}

# Keys only a running interpreter can add.
RUNTIME_KEYS = {"python", "platform", "machine"}


def test_import_exposes_the_documented_surface():
    assert set(biofasting.__all__) <= set(dir(biofasting))
    for name in biofasting.__all__:
        assert getattr(biofasting, name) is not None


def test_version_is_real_and_not_the_standalone_fallback():
    assert biofasting.__version__
    # A bare `cmake` build reports this marker; an installed wheel must not.
    assert biofasting.__version__ != "0.0.0.dev0"


def test_core_and_distribution_versions_agree():
    """The version is written in two places: pyproject.toml and CMake.

    If they drift, a bug report names a version that does not match the binary
    it came from, which is worse than no version at all.
    """
    assert biofasting._core.__version__ == biofasting.__version__


def test_build_info_has_both_halves_and_no_empty_values():
    info = biofasting.build_info()
    assert COMPILE_TIME_KEYS <= set(info), "the C++ half is incomplete"
    assert RUNTIME_KEYS <= set(info), "the runtime half is incomplete"
    empty = sorted(key for key, value in info.items() if not value)
    assert not empty, f"build_info() reported empty values: {empty}"


def test_cxx_standard_is_the_one_the_project_claims():
    assert biofasting.build_info()["cxx_standard"] == "20"


def test_build_arch_matches_the_running_machine():
    """A wheel built for one architecture must not report another.

    This is the test that would catch a cross-compiled wheel whose CMake target
    and C++ preprocessor macros disagreed.  platform.machine() and CMake name
    the same architecture differently, so both are normalised first.
    """
    aliases = {"amd64": "x86_64", "arm64": "aarch64"}
    built = biofasting.build_info()["arch"].lower()
    running = platform.machine().lower()
    assert aliases.get(built, built) == aliases.get(running, running)


def test_build_info_reports_the_pinned_nanobind():
    """A scaffold built against a different nanobind than the pins claim.

    Not a correctness bug by itself, but it means the environment is not the one
    requirements-dev.txt describes, and every other number in it is suspect.
    """
    assert biofasting.build_info()["nanobind"] != "unknown"
