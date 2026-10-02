"""Step 1.1: the runtime ISA dispatch interface exists and is self-consistent.

No kernel uses it yet -- step 1.2 is the first -- so there is no speed to check
here.  What is checked is the invariant every future kernel will rest on: the
selected level is one the CPU actually supports, and the ladder always has a
portable floor.  A dispatch that picked a rung the CPU lacks would be a SIGILL
in a user's pipeline, and it would not fail on the machine that built it.
"""

import biofasting


def test_baseline_is_always_the_floor():
    """Index 0 is the promise, not a feature: see src/core/cpu_features.hpp."""
    assert biofasting.cpu_levels()[0] == "baseline"


def test_ladder_has_no_duplicates():
    levels = biofasting.cpu_levels()
    assert len(set(levels)) == len(levels)


def test_supported_levels_are_a_subset_of_the_ladder():
    assert set(biofasting.supported_cpu_levels()) <= set(biofasting.cpu_levels())


def test_baseline_is_always_supported():
    """If the extension imported at all, its floor is available by definition."""
    assert "baseline" in biofasting.supported_cpu_levels()


def test_supported_levels_are_ascending():
    """detect() builds the list by appending rungs upward; nothing reorders it.

    best_cpu_level() relies on that ordering, so the ordering is asserted here
    rather than trusted.
    """
    ladder = biofasting.cpu_levels()
    positions = [ladder.index(level) for level in biofasting.supported_cpu_levels()]
    assert positions == sorted(positions)


def test_best_level_is_the_highest_supported():
    ladder = biofasting.cpu_levels()
    supported = biofasting.supported_cpu_levels()
    best = biofasting.best_cpu_level()
    assert best in supported
    assert best == max(supported, key=ladder.index)


def test_detection_is_stable_across_calls():
    """The answer is cached in a function-local static; it must not drift."""
    assert biofasting.supported_cpu_levels() == biofasting.supported_cpu_levels()
