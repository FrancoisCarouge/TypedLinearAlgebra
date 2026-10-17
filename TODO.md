# TODO

## Modules

* Use `import std;` in place of the standard library includes.
* Provide the library as a named module, alongside the headers.
* Provide the backend plug-ins as modules.
* Install, test, and exercise the modules in CI.

## Future Operations

The following top-10 operations and their variations are to be implemented. They are sorted by most commonly used, needed by end-users, in the industry. They are a partial superset of all supported backends. Prioritizing std::linalg and other standard library support.

| Priority | Operation Description | Prototype |
| --- | --- | --- |
| 1 | **Norms**: L1 sum of magnitudes and index of the largest magnitude of a vector, complementing the vector L2 `magnitude` and the `matrix_frob_norm`, `matrix_one_norm`, and `matrix_inf_norm` matrix norms. `std::linalg` `vector_abs_sum`, `vector_idx_abs_max`; BLAS ASUM, IAMAX; Eigen `.lpNorm<1>()`; NumPy `linalg.norm`. | `[[nodiscard]] constexpr auto vector_abs_sum(const rank_typed_matrix<1> auto &v);` |
| 2 | **Outer product**: symmetric, Hermitian, and rank-2 variants of the `matrix_rank_1_update` overwriting `A = x yᵀ` and updating `A = E + x yᵀ` rank-1 updates. `std::linalg` `symmetric_matrix_rank_1_update`, `hermitian_matrix_rank_1_update`, `symmetric_matrix_rank_2_update`, `hermitian_matrix_rank_2_update`; BLAS SYR, HER, SYR2, HER2. | `constexpr void symmetric_matrix_rank_1_update(const rank_typed_matrix<1> auto &x, rank_typed_matrix<2> auto &A, auto triangle);` |
| 3 | **Triangular solve**: solve `Ax = b` or `AX = B` for a lower or upper triangular `A`, the companion of the Cholesky and LU factorizations. `std::linalg` `triangular_matrix_vector_solve`, `triangular_matrix_matrix_left_solve`, `triangular_matrix_matrix_right_solve`; BLAS TRSV, TRSM; Eigen `.triangularView<Mode>().solve()`; SciPy `linalg.solve_triangular`. | `constexpr void triangular_matrix_vector_solve(const rank_typed_matrix<2> auto &A, auto triangle, auto diagonal, const column_typed_matrix auto &b, column_typed_matrix auto &x);` |
| 4 | **Inverse**: the matrix inverse, whose indexes are the reciprocal, transposed indexes of the operand. Eigen `.inverse()`; Armadillo `inv()`; LAPACK GETRF, GETRI; NumPy `linalg.inv`. | `[[nodiscard]] constexpr auto inverse(const rank_typed_matrix<2> auto &A);` |
| 5 | **Element-wise product**: Hadamard product, and quotient, of two matrices of identical shapes. Eigen `.cwiseProduct()`, `.cwiseQuotient()`; Armadillo `%`, `/`; NumPy `multiply`. | `[[nodiscard]] constexpr auto hadamard_product(const same_as_typed_matrix auto &lhs, const same_as_typed_matrix auto &rhs);` |
| 6 | **Eigendecomposition**: eigenvalues and eigenvectors of a square matrix, and of a self-adjoint matrix. Eigen `EigenSolver`, `SelfAdjointEigenSolver`; Armadillo `eig_gen()`, `eig_sym()`; LAPACK GEEV, SYEV; NumPy `linalg.eig`, `linalg.eigh`. | `[[nodiscard]] constexpr auto eigen_decomposition(const rank_typed_matrix<2> auto &A);` |
| 7 | **Cross product**: vector product of two three-element vectors. Eigen `.cross()`; Armadillo `cross()`; NumPy `cross`. | `[[nodiscard]] constexpr auto cross(const rank_typed_matrix<1> auto &lhs, const rank_typed_matrix<1> auto &rhs);` |
| 8 | **Singular value decomposition**: factors `A = U Σ Vᵀ`, the basis of the pseudo-inverse and of least-squares solutions. Eigen `JacobiSVD`, `BDCSVD`; Armadillo `svd()`; LAPACK GESVD, GESDD; NumPy `linalg.svd`. | `[[nodiscard]] constexpr auto singular_value_decomposition(const rank_typed_matrix<2> auto &A);` |
| 9 | **Trace**: sum of the diagonal elements of a square matrix whose diagonal elements are of compatible types. Eigen `.trace()`; Armadillo `trace()`; NumPy `trace`. | `[[nodiscard]] constexpr auto trace(const rank_typed_matrix<2> auto &A);` |
| 10 | **Cholesky decomposition**: lower triangular factor `L` of a symmetric positive-definite `A = L Lᵀ`, and its `L D Lᵀ` variant. Eigen `LLT`, `LDLT`; Armadillo `chol()`; LAPACK POTRF; NumPy `linalg.cholesky`; SciPy `linalg.cholesky`, `linalg.cho_factor`. | `[[nodiscard]] constexpr auto cholesky(const rank_typed_matrix<2> auto &A);` |

The `std::linalg` operations come first, then the others by decreasing usage, measured as call counts in the non-test code of about 4,000 scientific Python repositories surveyed by [python-api-inspect](https://github.com/Quansight-Labs/python-api-inspect): norm 6,690; inverse 2,088; outer product 1,990; element-wise product 1,911; eigendecomposition 1,798; cross product 1,506; singular value decomposition 1,190; trace 761; Cholesky 751; determinant 663; least-squares 568; triangular solve 153.

## Named values

Elevate the identity and zero named values from Kalman into this library, for reuse by any typed linear algebra client. Kalman pull request [#1026](https://github.com/FrancoisCarouge/Kalman/pull/1026) gives them their constant-initialized form (finding F41 in `BUG.md`):

- `one<Type>` and `zero<Type>`, `inline constexpr` variables, defined only when the backend's value is a constant expression, so that no client copies them before their dynamic initialization.
- The `ones<Type>` and `zeros<Type>` factories, the customization points a backend specializes, building the value at the point of use when it is not a constant expression, with `make_one<Type>()` and `make_zero<Type>()` choosing between the two.
- The `constant_evaluable` concept they rely on.
- Their specializations for `typed_matrix` (today in Kalman's `support/eigen_typed` and `support/quantity`) and for mp-units quantities (today in Kalman's `support/mp_units`), next to this library's `element_caster` specializations.

Today each client carries its own copy: Kalman's two facades, and ArduPilot's `AP_LinearAlgebra_Kalman.h`, which re-specializes the four for its backend. Kalman would then use this library's named values; its scalar `one<double>` and `zero<double>` specializations could stay in Kalman or move too.

## Proposals from the ArduPilot integration

Options for the owner to decide, not implemented. Finding numbers refer to `BUG.md`.

**In-place evaluation (F28).** Every typed operation returns a new matrix, so `P = P - K * H * P` on a 24-state covariance holds several 24x24 temporaries: 19 kB of stack in double precision on Cortex-M7, where ArduPilot allows 2.6 kB per function. ArduPilot works around it by checking the expression's type with `decltype` and writing the backend storage by hand.

- A. Compound assignment, `P -= K * HP;`, defined when the right-hand side's type converts to the left-hand side's, delegating to the backend's `-=` (Eigen evaluates `K * HP` into the destination without a temporary of P's size). Smallest change; the backend decides the temporaries.
- B. An explicit `assign(destination, expression)` evaluating a typed expression into existing storage, with the same type check. More general (also `=`), and makes the in-place intent visible at the call site.
- C. Rank-one updates without `<linalg>`: the existing `matrix_rank_1_update` signature, implemented by a loop over the backend storage when `std::linalg` is unavailable (libstdc++ 14 and 15, Arm GCC 15.2). Covers `P -= K (H P)` exactly, the EKF's most common step.

**Blocks and selections (F26, F39).** Grouped states and sparse Jacobians need typed sub-matrices.

- A. `block<Row, Column, Rows, Columns>(matrix)` returning a typed copy whose element types are the selected ones.
- B. The same as a typed view over the backend storage, writable, which would also give in-place updates of a block.
- C. `select<indexes...>(vector)` for non-contiguous selections, such as the five non-zero elements of an EKF3 observation Jacobian.
