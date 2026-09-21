# Changelog

All notable changes to TypedLinearAlgebra are documented here. The format
follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and versions
follow [Semantic Versioning](https://semver.org/).

## [Unreleased]

## [0.3.0] - 2026-09-12

### Added

- Operations: unary `minus` (with rank-0 support), `magnitude`,
  a generalized `transposed`, `equal_to`, and `matrix_vector_product`.
- Backend plug-ins: Au (`au`, `au_eigen`, `au_std`) and nholthaus/units
  (`nholthaus`, `nholthaus_eigen`, `nholthaus_std`), plus `std::chrono`
  duration integration (`chrono_eigen`, `chrono_std`).
- `distinct_typed_matrix` and `uniform_typed_matrix` concepts for
  element-type verification.
- Roadmap section and additional reference/use-case documentation, including
  AI agent guidance (`AGENTS.md`).

### Changed

- Matrix-matrix product constraints tightened; more tests for `transposed`,
  `scale`, and `matrix_product`.
- Division disambiguated from the other operators, with added quantity tests.

### Compiler & Build

- CI now covers Ubuntu 26.04 and Clang 19-22, alongside continued MSVC
  C++26 support and interprocedural-optimization/strict-aliasing hardening.
- Sanitizer coverage expanded with memory, control-flow integrity, and
  integer sanitizers; clang-tidy jobs sharded for turnaround.

### Fixed

- A typo in the `magnitude` assertion message.

### Dependencies

- Routine `dependabot` bumps to `codeql-action`, `harden-runner`,
  `checkout`, `cff-validator`, and `ossf/scorecard-action`.

## [0.2.0] - 2026-06-30

### Added

- Safe write interface and `at<i,j>()` compile-time element accessors.
- Index-literal operators, with expanded documentation and tests.
- Compile-fail CMake test support (the `fail()` test generator).
- Initial `benchmark/` suite, including matrix-product benchmarks.

### Changed

- Element and index-count checks now consistently require the number of
  indexes to match rank across `at`, construction, addition, substraction,
  multiplication, transpose, and operator tests.
- Documentation clarified around member and `explicit` semantics; test
  suite reorganized.

### Breaking Changes & Build Requirements

- **C++ Standard:** now targets C++26.
- **Compiler:** added Visual Studio 2026 / MSVC and GCC 15 pipelines,
  with MSVC warnings-as-errors enforced.
- Hardened structured-binding access was removed in favor of the safer
  write interface above; downstream code relying on the old accessor
  shape needs to migrate to `at<i,j>()`.

### Fixed

- Compile-time accessor and access-constraint bugs.
- gsl-lite `C4875` warning workaround on MSVC.
- A temporary mp-units regression workaround on MSVC.

### Dependencies

- Routine `dependabot` bumps to `codeql-action`, `harden-runner`,
  `checkout`, `add-pr-comment`, `actions-gh-pages`, `cff-validator`, and
  `dependency-review-action`.

## [0.1.0] - 2026-03-16

### Added

- Initial release: the `fcarouge::typed_matrix` strongly-typed facade over
  a linear algebra backend, composing (not inheriting) the underlying
  matrix representation.
- Core algorithms and operators: `add`, `divide`, `magnitude`,
  `matrix_product`, `product`, `scale`, `substract`, plus structured
  bindings and compile-time subscripting.
- Backend plug-ins: `std::linalg`, Eigen, Kokkos, and mp-units.
- `std::formatter` specialization for typed matrices.
- CMake install/export support (`fcarouge-typed-linear-algebra::tlinalg`),
  Doxygen documentation, and initial samples.

### Compiler & Build

- CI established across multiple compilers and Windows, with CppCheck and
  CodeQL integrated.

[Unreleased]: https://github.com/FrancoisCarouge/TypedLinearAlgebra/compare/0.3.0...HEAD
[0.3.0]: https://github.com/FrancoisCarouge/TypedLinearAlgebra/compare/0.2.0...0.3.0
[0.2.0]: https://github.com/FrancoisCarouge/TypedLinearAlgebra/compare/0.1.0...0.2.0
[0.1.0]: https://github.com/FrancoisCarouge/TypedLinearAlgebra/releases/tag/0.1.0
