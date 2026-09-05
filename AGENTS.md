# AGENTS.md

Guidance for coding agents working in this repository (the TypedLinearAlgebra
library, `git@github.com:FrancoisCarouge/TypedLinearAlgebra.git`). This file is
the repository root; all paths below are relative to it.

## Build & test

The repository root is the CMake source directory. The canonical loop — the same
one CI runs (`.github/workflows/pipeline.yml`) — is:

```shell
cmake -S . -B build -G Ninja
cmake --build build --parallel
ctest --test-dir build --parallel --verbose
```

- Requires CMake ≥ 4.3 (`cmake_minimum_required(VERSION "4.3")`) and a C++26
  compiler. Local builds use whatever `cc`/`c++` resolve to on `$PATH` unless you
  pass `-DCMAKE_CXX_COMPILER=`. CI exercises `clang++-20`, `clang++-21` (Ninja)
  and `g++-14`, `g++-15` (Unix Makefiles) on `ubuntu-26.04`, plus MSVC `cl`
  (generator `Visual Studio 18 2026`) on `windows-2025`, in both Debug and
  Release.
- MSVC can't take `CMAKE_CXX_STANDARD 26` (CMake doesn't record C++26 support for
  it, so requesting it hard-fails feature checks). `support/support.cmake`
  `unset()`s `CMAKE_CXX_STANDARD` under `if(MSVC)`, and `support/CMakeLists.txt`'s
  MSVC branch restores the dialect via `/std:c++latest` directly — the root
  `CMakeLists.txt` stays compiler-agnostic. Replicate that split rather than
  adding an `if(compiler)` branch at the root.
- To run a single test, use CTest's `-R`/`--tests-regex` against the generated
  test name (see the naming convention below), e.g. `ctest --test-dir build -R
  typed_linear_algebra_eigen_addition_1x2_eigen_pass --output-on-failure`.
- Formatting (enforced by `.github/workflows/format.yml`): `clang-format-22
  --Werror -i -style=file` on `.hpp`/`.tpp`/`.cpp`, and `cmake-format -i` on
  `CMakeLists.txt`/`*.cmake`.
- `.clang-tidy` runs `Checks: '*'` with `WarningsAsErrors: '*'` (minus a short
  denylist) — a very strict baseline; don't casually suppress warnings.
  `.github/workflows/clang_tidy.yml` runs it with `clang++-20` against
  `compile_commands.json`.
- Pre-commit hooks (`.pre-commit-config.yaml`): gitleaks, shellcheck, cpplint,
  end-of-file-fixer, trailing-whitespace.
- Install: `sudo cmake --install build`.

## Architecture

### Core header

`include/fcarouge/typed_linear_algebra.hpp` is the single public entry point
(`fcarouge::typed_matrix`). It declares concepts, the `typed_matrix<Matrix,
RowIndexes, ColumnIndexes>` class template, operators, and the
`element_caster`/`cast` customization point, then pulls in implementation from
`include/fcarouge/typed_linear_algebra_internal/` via `.tpp` includes at the
bottom of the header (`algorithm/{add,divide,equal_to,magnitude,matrix_product,
matrix_vector_product,minus,product,scale,substract,transposed}.tpp`, `cast.tpp`,
`chrono.tpp`, `common_type.tpp`, `format.tpp`, `tuple.tpp`,
`typed_linear_algebra.tpp`, plus the non-`.tpp` `utility.hpp`). `typed_matrix`
**composes** (does not inherit) an underlying backend `Matrix`; each element's
type comes from the `RowIndexes`/`ColumnIndexes` tuples, not from storage.

Key design tradeoff baked into the API (documented in-header): there is no
runtime-indexed safe element access — `at<i,j>()` is necessarily a compile-time
template API, because converting a runtime index to a dependent template type
isn't possible in C++. Lvalue reference element assignment is likewise not
offered; see the "Lessons Learned" section of `README.md` for the full rationale
before trying to "fix" either of these.

### Backend plug-in pattern

The library is backend-agnostic: `support/<backend>/fcarouge/linalg.hpp` is what
test/sample code actually includes (`#include "fcarouge/linalg.hpp"`), and
*which* backend's header is picked up is determined purely by which CMake target
you link against — not by any `#ifdef` in the core library. Backends live under
`support/` (wired in by `support/CMakeLists.txt`, gated on `BUILD_TESTING`):

- `eigen`, `eigexed`, `nested_typed_eigen` — Eigen-backed variants.
- `kokkos` — Kokkos/mdspan-backed.
- `mp_units`, `mp_units_eigen`, `mp_units_std` — mp-units (quantities/units)
  integration; also defines `fcarouge/mp_units.hpp`.
- `au`, `au_eigen`, `au_std` — Au (quantities/units) integration; also defines
  `fcarouge/au.hpp`. Unlike mp-units, Au has no vector/matrix quantity
  representation, so `fcarouge/au.hpp` only needs `element_caster`
  specializations — no tuple-decomposition/`disable_representation` machinery.
  Two rough edges to know before touching Au-typed code: an unqualified `matrix *
  bare_quantity` expression can resolve to Au's own hidden-friend `operator*` via
  ADL instead of this library's, requiring an explicit `fcarouge::operator*(...)`
  call (see `test/multiplication/scalar_au_eigen.cpp`); and Au's
  `Quantity::operator=` is lvalue-only, which is incompatible with the
  `add()`/`substract()` std::linalg-style free functions' element-compatibility
  check (no `au_std` `add()`/`substract()` tests exist for this reason — use
  `operator+`/`operator-` instead).
- `nholthaus`, `nholthaus_eigen`, `nholthaus_std` — nholthaus/units integration
  (`FetchContent` of `github.com/nholthaus/units`, or a system `units` package);
  also defines `fcarouge/nholthaus.hpp`.
- `chrono_eigen`, `chrono_std` — test/sample glue for the `std::chrono::duration`
  element-type integration wired into the core header via `chrono.tpp`; each only
  defines `fcarouge/linalg.hpp`, there is no bare `chrono` target.
- `main` — shared CTest `main()` driver linked into `pass`/`fail` test
  executables.
- Built-in scalar types and `std::linalg` need no plug-in at all.

Adding a new backend means adding a `support/<name>/` directory with its own
`fcarouge/linalg.hpp`, a `typed_linear_algebra_<name>` CMake target, and an
`add_subdirectory("<name>")` line in `support/CMakeLists.txt` (see
`support/support.cmake` and existing backends for the pattern).

### Test/benchmark generation (`support/support.cmake`)

Tests and benchmarks are not hand-declared per backend; four CMake functions
generate them from a test name + a `BACKENDS` list:

- `pass(NAME BACKENDS ...)` — compiles `NAME.cpp` once per backend and registers
  it as a CTest test (assertions run at static-init via `main` from
  `support/main`).
- `build(NAME BACKENDS ...)` — compile-only (`OBJECT` library,
  `EXCLUDE_FROM_ALL`, no `main` linked); the CTest test just builds the target.
  Use it for `static_assert`-only tests such as the concept checks in
  `same_as_typed_matrix/`, `distinct/`, `nested/`, `underlying/`.
- `fail(NAME BACKENDS ...)` — same shape as `pass`, but expects a **compile
  failure** (`EXCLUDE_FROM_ALL` target invoked via `cmake --build --target`,
  `WILL_FAIL TRUE`).
- `bench(NAME SIZE BACKENDS ...)` — configures `NAME.cpp` → `NAME_SIZE.cpp` and
  registers a nanobench-based benchmark.

Generated target/test names follow
`typed_linear_algebra_<backend>_<dir>_<name>_pass` (`pass` and `build` both use
the `_pass` suffix; `..._<name>` for `fail`, `..._<name>_<size>_bench`), where
`<dir>` is the calling directory's name. `test/` is organized one directory per
operation/feature/concept — `addition/`, `assign/`, `at/`, `common_with/`,
`constructor/`, `copy/`, `distinct/`, `division/`, `element/`, `equal_to/`,
`format/`, `interconvertible/`, `magnitude/`, `matrix_product/`,
`matrix_vector_product/`, `minus/`, `mp_units/`, `multiplication/`, `nested/`,
`operator/`, `row_typed_matrix/`, `same_as_typed_matrix/`, `scale/`,
`structured_bindings/`, `substraction/`, `transposed/`, `underlying/` — each with
its own `CMakeLists.txt` listing `pass(...)`/`build(...)`/`fail(...)` calls per
test file and backend combination (some dirs, e.g. `element/` and
`row_typed_matrix/`, mix `build()` concept checks with `fail()` cases);
`test/CMakeLists.txt` itself only `add_subdirectory`s them. Follow that pattern
when adding a new test case rather than writing a bespoke `add_test`.

Test files themselves are minimal: `#include "fcarouge/linalg.hpp"`, then a
single `[[maybe_unused]] const auto test{[] { ...; assert(...); return 0; }()};`
block inside `namespace fcarouge::test { namespace { ... } }` — no test
framework, just `<cassert>` run at static-init time via `main` from
`support/main`.

### Other directories

- `sample/` — usage examples (e.g. `au_eigen.cpp`, `au_std.cpp`,
  `chrono_eigen.cpp`, `chrono_std.cpp`, `mp_units_eigen.cpp`, `mp_units_std.cpp`,
  `nholthaus_eigen.cpp`, `nholthaus_std.cpp`), built the same `pass`-style way
  via `sample/CMakeLists.txt`.
- `benchmark/` — perf comparisons across raw Eigen, mdspan, and typed wrappers;
  results plotted via `plot.cpp` to `plot.png`.
- `documentation/` — Doxygen config/theme.
- `cmake/`, `pkgconfig/` — install/find_package export files
  (`fcarouge-typed-linear-algebra-config.cmake.in`, `.pc.in`).
- The root `CMakeLists.txt` only `add_subdirectory`s `benchmark`, `pkgconfig`,
  `sample`, `support`, `test` when `PROJECT_IS_TOP_LEVEL` — so consumers using
  `FetchContent`/`find_package` only pull in `cmake/` + `include/`.

## Conventions

- Header-only library: keep the public API under `include/fcarouge/`;
  implementation details belong in `typed_linear_algebra_internal/` and are not
  part of the public surface.
- Every source/CMake file carries the Unlicense SPDX header block — copy it
  verbatim (version/URL match the root `CMakeLists.txt`) into any new file.
- The public CMake target consumers link against is
  `fcarouge-typed-linear-algebra::tlinalg`; the package name is
  `fcarouge-typed-linear-algebra`.
- The author writes precise, terminology-careful `@note`/`@todo`/`@warning`
  doxygen comments explaining design rationale directly in headers — match that
  register when editing docs/comments rather than simplifying.
