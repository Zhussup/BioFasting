"""BioFasting: a C/C++-backed core for sequence bioinformatics.

Nothing here is a public API yet.  Step 1.1 built the scaffold, and the first
kernels arrive in step 1.2; what this module provides today is the ability to
state exactly which build is running, which is the precondition for every
performance claim the project intends to make.  A number without a build
identity behind it cannot be reproduced, and an unreproducible number is not
evidence.

The public surface is re-exported from the compiled core rather than
reimplemented, so there is one implementation of each answer and not two.
"""

import platform
from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _distribution_version

from . import _core

__all__ = [
    "__version__",
    "best_cpu_level",
    "build_info",
    "cpu_levels",
    "supported_cpu_levels",
]

try:
    __version__ = _distribution_version("biofasting")
except PackageNotFoundError:  # an extension copied out of its wheel
    __version__ = _core.__version__


def build_info() -> dict[str, str]:
    """Return everything needed to reproduce a bug report or a benchmark run.

    The compile-time half (version, compiler, flags, target architecture) comes
    from the C++ core.  The interpreter and platform half can only be known at
    runtime, so it is added here -- a binary cannot recover from the CPU it is
    executing on which compiler produced it.
    """
    info = dict(_core.build_info())
    info["python"] = platform.python_version()
    info["platform"] = platform.platform()
    info["machine"] = platform.machine()
    return info


def cpu_levels() -> list[str]:
    """Return the ISA ladder this build defines, ascending, baseline first."""
    return _core.cpu_levels()


def supported_cpu_levels() -> list[str]:
    """Return the subset of :func:`cpu_levels` the running CPU supports."""
    return _core.supported_cpu_levels()


def best_cpu_level() -> str:
    """Return the highest supported level: the rung a kernel dispatches to."""
    return _core.best_cpu_level()
