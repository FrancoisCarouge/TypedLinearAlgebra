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
  `CMakeLists.txt`/`*.cmake`. There is **no `.clang-format` file** in the repo,
  so `-style=file` falls back to clang-format's built-in LLVM style (80
  columns, 2-space indent, template declarations always broken). Hand-written
  layout will not match it — always finish an edit by running
  `clang-format-22 -i -style=file` (that exact version; v18/20/21 disagree on
  wrapping) and `cmake-format -i` on every file you touched. The CI check is
  `find . -not -path './.git/*' \( -iname '*.hpp' -o -iname '*.cpp' -o -iname
  '*.tpp' \) | xargs clang-format-22 --Werror --dry-run -style=file`.
- `.clang-tidy` runs `Checks: '*'` with `WarningsAsErrors: '*'` (minus a short
  denylist) — a very strict baseline; don't casually suppress warnings.
  `.github/workflows/clang_tidy.yml` configures the build with `clang++-20`
  (`-DCMAKE_EXPORT_COMPILE_COMMANDS=ON`), then `run-clang-tidy -p build` over
  every `benchmark|sample|support|test/*.cpp` except `*_fail.cpp`, sharded 6
  ways, using the runner's default `clang-tidy` (v21 as of writing, unpinned)
  and `compile_commands.json`. The codebase has **zero `NOLINT` comments** —
  keep it that way; rework the code instead. Gotchas that bite test files:
  - `misc-include-cleaner`: every `std::` symbol needs its own direct `#include`
    even if transitively available (`std::milli` → `<ratio>`, `std::identity` →
    `<functional>`, …). `.clang-tidy`'s `IgnoreHeaders` whitelists only the
    Eigen/mp-units/au/`fcarouge/*` facade headers, never the standard library.
  - `google-explicit-constructor` / `hicpp-explicit-conversions` reject any
    non-`explicit` converting constructor or conversion operator — you can't
    write an implicit-conversion fixture; test a shared base class instead.
  - Run it before finishing: `clang-tidy-21 <file> -- -std=c++26 -Iinclude
    -Isupport/<backend> -Isupport/eigen <dep -isystem flags>`.
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
- `armadillo`, `armadilloxed` — Armadillo-backed, mirroring `eigen`/`eigexed`.
  `armadillo.hpp` fetches Armadillo (`gitlab.com/conradsnicta/armadillo-code`,
  `SOURCE_SUBDIR` set to a `CMakeLists.txt`-less dir so its own build never
  runs) for headers only and compiles it with `ARMA_DONT_USE_{LAPACK,BLAS,
  WRAPPER,ARPACK,SUPERLU,HDF5,FFTW3}` — no runtime library is linked. Because
  `solve()`/`inv()` need LAPACK, **matrix division (`operator/`) is unsupported
  with any Armadillo backend**; `test/division/` has no Armadillo entries and
  the `mp_units_armadillo` sample notes the gap (compare the `au_std`
  `add()`/`substract()` omission). Facade specifics: `is_armadillo` keys off
  Armadillo's `arma::is_arma_type` trait (a `fixed<R,C>`'s CRTP base is
  `Mat<eT>`, not itself); a constrained `std::formatter` partial specialization
  plus a `std::format_kind = disabled` opt-out (Armadillo matrices are
  ranges); ADL `operator==`/`!=` returning `bool` in `namespace arma` (raw
  Armadillo `==` yields an element-wise expression). The bare `armadillo`
  backend is dropped from a few test lines where raw-Armadillo semantics differ
  (`{{1},{4},{7}}` column-literal ctor ambiguity); `armadilloxed` covers those.
- `au_armadillo`, `mp_units_armadillo`, `nholthaus_armadillo`, `chrono_armadillo`
  — the unit-library combos over Armadillo, mirroring the `*_eigen` combos.
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

### Algorithms (`include/fcarouge/typed_linear_algebra_internal/algorithm/`)

One `.tpp` per algorithm (`add`, `divide`, `equal_to`, `magnitude`,
`matrix_product`, `matrix_vector_product`, `minus`, `product`, `scale`,
`substract`, `transposed`), each `#include`d near the bottom of
`typed_linear_algebra.hpp` (alphabetically, ~lines 527-537) and **re-declared for
documentation** in the `//! @name Algorithms` block just below those includes
(redeclaring an already-defined function — that block is the curated public API
listing doxygen renders). All live in `namespace fcarouge`, guarded by
`FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_ALGORITHM_<NAME>_TPP`. `namespace tla =
typed_linear_algebra_internal` is available (aliased in the core header; some
`.tpp`s re-alias it).

Core principle: an algorithm never re-implements representation arithmetic. It
(1) computes the *result's* `row_indexes`/`column_indexes` tuples, (2) statically
checks element-type compatibility, (3) delegates the numerics to the backend via
`.data()` (the underlying backend matrix). Only `magnitude` is implemented atop
the strong types (`cast` round-trips + `using std::sqrt`) — flagged in-file as a
tradeoff to revisit.

Two deliberately mirrored API families:

- **Expression / operator style** — mirrors Eigen and math notation: `operator+`,
  `operator-` (binary and unary negation), `operator*`, `operator/`,
  `operator==`, plus `magnitude`, `transposed`. `[[nodiscard]] constexpr auto`,
  return by value, result built with `make_typed_matrix<row_indexes,
  column_indexes>(<backend expression on .data()>)`. Works on every backend.
  `make_typed_matrix` composes an element-accessible backend expression (Eigen)
  as-is, but materializes one that is not (Armadillo `Glue`/`Op`) through its
  `.eval()` member first, so the resulting typed matrix stays usable.
- **`std::linalg` free-function style** — mirrors the `std::linalg` names and
  out-parameter signatures: `add`, `matrix_product`, `matrix_vector_product`,
  `scale`. `constexpr void f(inputs..., result&)`, caller pre-allocates
  `result`, body is `using std::linalg::f; f(lhs.data(), ..., result.data());`.
  Whole file wrapped in `#ifdef __cpp_lib_linalg` / `#include <linalg>` (`@todo`
  remove for native C++26), so only the `*_std` backends exercise them.
  `matrix_vector_product` also reshapes the vector's rank-two n-by-1 storage into
  the rank-one span `std::linalg` wants via a local `as_vector_span` mdspan
  helper. `transposed` straddles both families: `if constexpr (requires {
  value.data().transpose(); })` (Eigen member), else `requires { value.data().t();
  }` (Armadillo member), else `std::linalg::transposed`.

Rank dispatch: overloads constrained on `rank_typed_matrix<0|1|2> auto` (0 =
singleton, 1 = row/column vector, 2 = matrix; `rank` computed in `utility.hpp`).
Singletons usually get a dedicated simpler overload plus scalar-interop overloads
against the `other` concept (any non-typed-matrix operand), which route through
`element_caster` / `cast<underlying, T>`. Constraint vocabulary lives in
`utility.hpp`: `same_as_typed_matrix`, `same_shape`, `uniform_typed_matrix`,
`row_typed_matrix`, `column_typed_matrix`, plus `multipliable*` / `scalable_by`
in `product.tpp`.

Type-level index math: `tla::product` / `tla::quotient` (the `multiplies` /
`divides` metafunction objects in `utility.hpp`, specialized over `std::tuple`
and `std::identity`; `std::identity` is the dimensionless "1" index). Element
checks run as `tla::for_constexpr<N>` loops of `static_assert(requires {
std::declval<lhs_element>() OP std::declval<rhs_element>(); }, "<Operation>
requires compatible element types.")`; `add()` additionally probes assignability
into the result element. `matrix_product`, `matrix_vector_product`, and the
`operator/` overloads still carry `@todo`s for element verification — don't
assume it is there.

Checks route through the *typed* element access, never the raw storage: storage
holds bare representation oblivious to the indexes, so verifying it directly
would silently accept e.g. a length matrix vs. an area matrix whose reps happen
to match (see `equal_to.tpp` `@details`); it also sidesteps backends (mdspan)
whose storage type has no `operator==`.

`scalable_by` (in `product.tpp`) exists so the deduced-return scalar overloads
stay SFINAE-friendly: an unconstrained deduced return forces the compiler to
instantiate — and hard-error in — the body when third-party code (mp-units)
merely `requires`-probes multiplication by an element type. Constrain any new
deduced-return scalar overload the same way.

Known Au friction (see also the `au` backend note above): Au's `Quantity`
hidden-friend `operator*` can win ADL over this library's, needing an explicit
`fcarouge::operator*(...)`; and Au's lvalue-only `Quantity::operator=` fails the
`add()`/`substract()` result-assignability probe, so those free functions have no
`au_std` tests — use `operator+` / `operator-` instead. The codebase spells it
"substract" / "substraction" throughout (sic) — match it.

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
`format/`, `magnitude/`, `matrix_product/`,
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

Concept checks that inspect element types (`uniform_typed_matrix`,
`distinct_typed_matrix`, `element/`) have traps — `distinct_typed_matrix`
especially, because it is the only check that runs `std::is_convertible` /
`std::common_type` between *different* element types:

- In a two-index `matrix<Rep, RowIndexes, ColumnIndexes>` an element's type is
  `RowType * ColType` (`typed_linear_algebra_internal::product`). That product
  must be well-formed for *every* row/col index pair, so keep `std::identity`
  (via `std::tuple<std::identity>`) on one axis — `product<duration, duration>`,
  `product<year, month>`, etc. are ill-formed and hard-error.
- Au: checking `is_convertible` between two same-dimension quantities of
  different *representation* (`QuantityD` vs `QuantityI`, i.e. a float→int rep
  conversion) instantiates Au's `overflow_boundary.hh`, whose
  `ValueOfMaxFloatNotExceedingMaxInt::max_mantissa` has a local `ONE` that
  shadows `au::ONE` — MSVC C4459, fatal under `/WX`. Compare quantities of the
  same representation (different unit, e.g. `Meters` vs `Kilo<Meters>`, is
  fine); if a rep-conversion test is unavoidable, add `/wd4459` to
  `support/au/CMakeLists.txt` next to the existing `/wd4244`.
- Prefer distinctness assertions that rest on standard-mandated behaviour over
  QoI: `duration` pairs always share a `std::common_type`; a `duration` vs a
  `time_point` reliably does not (`time_point(duration)` is mandated
  `explicit`). Two `std::chrono` calendar types (`year`/`month`/`day`) having
  no `common_type` is *not* clearly mandated — don't lean on it.

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

## Recipe: adding an algorithm

1. **Header.** Create
   `include/fcarouge/typed_linear_algebra_internal/algorithm/<name>.tpp`: verbatim
   Unlicense SPDX block, include guard
   `FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_ALGORITHM_<NAME>_TPP`, `namespace
   fcarouge`. Choose the API family — operator/expression (works on every
   backend) or `std::linalg` free function (wrap the whole file in `#ifdef
   __cpp_lib_linalg` + `#include <linalg>`, use the `void f(inputs..., result&)`
   out-parameter signature). When `std::linalg` has the operation, mirror its
   name and semantics and add `//! @see std::linalg::<name>`; otherwise mirror
   Eigen/common nomenclature (README already tracks the mapping).
2. **Overloads per rank.** Provide `rank_typed_matrix<2>`, `<1>`, `<0>` overloads
   as the operation admits; add `other`-operand overloads for scalar interop
   (route values through `cast<underlying, T>`). Constrain any deduced-return
   scalar overload on `scalable_by` or an equivalent `requires` so it stays
   SFINAE-friendly.
3. **Result type.** Compute `row_indexes` / `column_indexes` with `tla::product`
   / `tla::quotient` over the operand index tuples, then `return
   make_typed_matrix<row_indexes, column_indexes>(<expression on .data()>)` — or
   write into `result.data()` for the free-function family.
4. **Element checks.** `tla::for_constexpr` over rows/columns (or `rows *
   columns` for rank 1) with `static_assert(requires { std::declval<lhs_element>()
   OP std::declval<rhs_element>(); }, "<Operation> requires compatible element
   types.")`; also assert result-element assignability for out-parameter
   algorithms. Add `same_shape` / size `static_assert`s as the math requires.
5. **Wire into `typed_linear_algebra.hpp`.** Add the
   `#include "typed_linear_algebra_internal/algorithm/<name>.tpp"` line
   (alphabetical, with the rest, ~line 527) **and** a matching forward
   re-declaration with full `@brief`/`@details`/`@see`/`@todo` doxygen in the
   `//! @name Algorithms` block below the includes.
6. **README.** Add a row to the "Operations" table.
7. **Tests.** New `test/<name>/` directory (or reuse an existing operation dir)
   with its own `CMakeLists.txt` (SPDX comment block, then `fail(...)` lines then
   `pass(...)` lines), and an alphabetical `add_subdirectory("<name>")` in
   `test/CMakeLists.txt`. `test/CMakeLists.txt` only aggregates — never a bespoke
   `add_test`.
   - **File naming:** `<shape>[_<backend>][_fail].cpp` — `1x1` singleton,
     `1x2`/`2x1` vectors, `2x2`/`2x3`/`3x2` matrices, `RxC_R'xC'` for products,
     plus semantic names (`row`, `column`, `scalar`, `matrix_rhs`, `vector_lhs`,
     `rank_mismatch`). Generated test name:
     `typed_linear_algebra_<backend>_<name>_<file>_pass`.
   - **Backends & types to cover:** expression-style → `eigen`, `eigexed`,
     `nested_typed_eigen` with plain `double`, plus `au_eigen`,
     `mp_units_eigen`, `chrono_eigen`, `nholthaus_eigen` for typed elements;
     `std::linalg`-style → `au_std`, `mp_units_std`, `chrono_std`,
     `nholthaus_std` only (they need `<linalg>`/`<mdspan>` storage). Cover at
     least one plain and two unit backends, and one mixed-type vector (e.g.
     `position, velocity`) to prove per-element typing. `chrono` also exercises
     semantic rejection (time_point + time_point must not compile).
   - **`pass` body:** the standard static-init form —
     `[[maybe_unused]] const auto test{[] -> int { ...; assert(...); return 0;
     }()};` inside `namespace fcarouge::test { namespace { ... } }`, `#include
     "fcarouge/linalg.hpp"` + `<cassert>` + backend unit headers, with a `//!
     @test` one-liner on top. Assert the numeric result **and** the result type
     (`static_assert(std::same_as<decltype(r)::row_indexes, ...>)`), across the
     accessor forms (`.at<i>()`, `.at<i_i>()`, `[i_i]`, `(i_i)`) where relevant.
   - **`fail` test per rejected shape/rank/type mismatch:** keep the correct line
     as an `// Intended:` comment directly above the broken one; `fail(...)` sets
     `WILL_FAIL`.
   - Use `build(...)` instead of `pass(...)` for `static_assert`-only tests.
8. **Sample** (optional): add a `sample/<name>.cpp` if the algorithm is
   user-facing and illustrative, wired via `sample/CMakeLists.txt` the same
   `pass`-style way.
9. **Verify.** `cmake --build build --parallel && ctest --test-dir build
   --parallel`; single test: `ctest --test-dir build -R
   typed_linear_algebra_<backend>_<name>_<file>_pass --output-on-failure`. Keep
   `clang-format-22 --Werror` / `clang-tidy '*'` clean and doxygen
   warning-free.

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
