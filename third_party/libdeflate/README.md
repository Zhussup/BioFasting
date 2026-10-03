# libdeflate 1.23, vendored

Not written by hand: `python3 tools/vendor_libdeflate.py` put these files here
and re-running it is how they are updated.  Read this before editing anything
in the tree.

## Where it came from

| | |
|---|---|
| Upstream | <https://github.com/ebiggers/libdeflate> |
| Release | tag `v1.23` |
| Tarball | <https://github.com/ebiggers/libdeflate/archive/refs/tags/v1.23.tar.gz> |
| SHA-256 | `1ab18349b9fb0ce8a0ca4116bded725be7dcbfa709e19f6f983d99df1fb8b25f` |
| Version in `libdeflate.h` | `1.23` |
| Licence | MIT, in [COPYING](COPYING) |

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
wheel self-contained, and `build_info()["libdeflate"]` says `1.23 (vendored,
static)` so a gz number can be traced to the binary that produced it.

## What is here

Taken, the library whole:

- `COPYING`
- `common_defs.h`
- `libdeflate.h`
- `lib/**` -- 38 files: the 11 C sources upstream's own `CMakeLists.txt` lists
  in `LIB_SOURCES`, plus every header they include.  The list travels in
  [libdeflate_sources.cmake](libdeflate_sources.cmake), which the build
  includes, so the two cannot drift apart.

Left behind, everything else in the release:

- `.cirrus.yml`
- `.github`
- `.gitignore`
- `CMakeLists.txt`
- `NEWS.md`
- `README.md`
- `libdeflate-config.cmake.in`
- `libdeflate.pc.in`
- `programs`
- `scripts`

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
