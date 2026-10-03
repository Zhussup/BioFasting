"""BioFasting: a C/C++-backed core for sequence bioinformatics.

Step 1.1 built the scaffold and step 1.2 has begun: the FASTQ kernel is here
(`open_fastq`), so is the FASTA index (`open_fasta`), and so are the sequence
operations (`reverse_complement`, `gc_fraction`, `count_kmers`) and the
translation of a reading frame (`translate`) -- the first kernels in the project
that dispatch on the CPU at runtime.  `biofasting.sequtils` adds the other half
of `Bio.SeqUtils`: the functions that measure a sequence rather than change it,
from `GC123` to `molecular_weight`, plus the four checksums in
`biofasting.checksum`.  What the scaffold
provides besides them is the ability to
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
from .annotations import AnnotationReference, read_annotations, to_reference
from .checksum import crc32, crc64, gcg, seguid
from .fasta import FastaGrid, open_fasta, open_fasta_grid, read_fasta
from .fastq import FastqGrid, open_fastq, open_fastq_grid, open_fastq_index, scan_fastq
from .genbank import (
    FLAT_FILE_FORMATS,
    FlatFileRecord,
    open_genbank,
    read_genbank,
)
from .interop import (
    fasta_seqrecords,
    fastq_seqrecords,
    from_seqrecord,
    to_seqrecord,
    write_fasta,
    write_fastq,
    write_qual,
)
from .seqops import BiopythonWarning, count_kmers, gc_fraction, reverse_complement, translate
from .seqfeature import (
    Feature,
    Location,
    Part,
    Position,
    Qualifier,
    read_features,
    to_location,
    to_seqfeature,
)
from .isoelectric_point import IsoelectricPoint
from .protparam import ProteinAnalysis
from .sequtils import (
    GC123,
    GC_skew,
    CodonAdaptationIndex,
    molecular_weight,
    nt_search,
    seq1,
    seq3,
    six_frame_translations,
    xGC_skew,
)

__all__ = [
    "__version__",
    "AnnotationReference",
    "BiopythonWarning",
    "CodonAdaptationIndex",
    "FastaGrid",
    "FastqGrid",
    "GC123",
    "GC_skew",
    "IsoelectricPoint",
    "ProteinAnalysis",
    "FLAT_FILE_FORMATS",
    "Feature",
    "FlatFileRecord",
    "Location",
    "Part",
    "Position",
    "Qualifier",
    "best_cpu_level",
    "build_info",
    "count_kmers",
    "cpu_levels",
    "crc32",
    "crc64",
    "fasta_seqrecords",
    "fastq_seqrecords",
    "from_seqrecord",
    "gc_fraction",
    "gcg",
    "molecular_weight",
    "nt_search",
    "open_fasta",
    "open_fasta_grid",
    "open_fastq",
    "open_fastq_grid",
    "open_fastq_index",
    "open_genbank",
    "read_annotations",
    "read_fasta",
    "read_features",
    "read_genbank",
    "reverse_complement",
    "scan_fastq",
    "seq1",
    "seq3",
    "seguid",
    "seqops_level",
    "six_frame_translations",
    "supported_cpu_levels",
    "to_location",
    "to_seqfeature",
    "to_seqrecord",
    "to_reference",
    "translate",
    "write_fasta",
    "write_fastq",
    "write_qual",
    "xGC_skew",
]


def seqops_level() -> str:
    """Return the rung the sequence kernels dispatch to.

    Not the same question as :func:`best_cpu_level`: a CPU can support AVX-512
    while no kernel here uses it, and this is the answer a benchmark has to
    report, because it is the code that actually ran.  Set
    ``BIOFASTING_SEQOPS_LEVEL=baseline`` before the first call to pin it to the
    scalar path -- it only ever lowers the rung, never raises it.
    """
    return _core.seqops_level()

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
