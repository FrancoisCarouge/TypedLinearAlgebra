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
