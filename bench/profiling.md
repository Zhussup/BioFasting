# Profiling notes (Phase 0, step 0.4)

The runner answers *how much* slower Biopython is.  This answers *why*, and
therefore which parts of it a C core can actually delete.

Reproduce:

```sh
.venv/bin/python bench/profile_paths.py --all --repeat 1 --top 8
.venv/bin/py-spy record -o /tmp/prof.svg --format speedscope \
    -- .venv/bin/python bench/profile_paths.py --workload seqio_fastq
```

`profile_paths.py` reuses the runner's own workload factories, so the code being
profiled is the code that was measured — the two sets of numbers describe the
same work.

## How to read these numbers

**Shares, not seconds.**  `cProfile` instruments every Python call.  On
call-heavy code that inflates wall time roughly 3x (`seqio_fastq`: 9.9 s
profiled against 3.18 s in `bench/results/latest.json` for the same pass).  A
profiled second is not a real second; a profiled *percentage* is still a fair
statement of where the real time goes.

**Two Python numbers, not one.**  The harness frame — the loop in
`profile_paths.py` that consumes records — is counted separately from library
Python.  It exists to keep the comparison honest: `cProfile` charges time spent
inside a C extension to whichever Python frame called it, so for a C-backed
implementation the harness loop silently absorbs the C library's time.  Counting
it separately means "library Python" measures the same thing in every row.

**Where cProfile is the wrong instrument.**  For an implementation whose work
happens in C, `cProfile` cannot see inside it; it reports the harness frame as
~90% and that tells you nothing.  Those rows are read with `py-spy`, which
samples the real process without instrumentation and so attributes time to the
native frames.

## What the shares say

The decisive number is **library Python share**: Python bytecode inside
Biopython (or inside the harness's own parse loop) that a C implementation
deletes outright, as opposed to time already spent in C that a rewrite would
merely have to match.

<!-- filled from the profile run -->
