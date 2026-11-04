# Sparse Matrices

A typed matrix storage whose sparsity pattern is known at compile time. The pattern is part of the type: only the stored elements occupy memory and cost arithmetic operations, every other element is a structural zero. There is no memory allocation, no exception, and no input-output, and every operation is `constexpr`.

```cpp
constexpr auto pattern{sparse::make_pattern<2, 2>({{0, 0}, {1, 0}})};

using matrix_type = sparse_matrix<double, pattern>;

using sparse::structural_zero;

// sizeof(value) == 2 * sizeof(double)
matrix_type value{sparse::make_matrix<matrix_type>(1., structural_zero,
                                                   2., structural_zero)};

value.at<1, 0>(4.);         // A stored element.
value.at<1, 1>();           // A structural zero reads as a zero.
// value.at<1, 1>(4.);      // Does not compile.
auto square{value * value}; // Stores the structural product only.
```

- [Reference](#reference)
  - [CMake](#cmake)
  - [Include](#include)
  - [Namespace](#namespace)
  - [Plug-ins](#plug-ins)
  - [Samples](#samples)
  - [Aliases](#aliases)
  - [Class Matrix](#class-matrix)
  - [Class Pattern](#class-pattern)
  - [Patterns](#patterns)
  - [Pattern Algebra](#pattern-algebra)
  - [Element Construction](#element-construction)
  - [Element Access](#element-access)
  - [Operations](#operations)
  - [Concepts](#concepts)
  - [Format](#format)
- [Notes](#notes)
- [Example](#example)
- [See Also](#see-also)

## Reference

### CMake

Link against the target of the plug-in of your element types, listed below. The plug-ins are defined in-tree when `BUILD_TESTING` is enabled and are not installed: copy or adapt them, as for any other [plug-in](https://github.com/FrancoisCarouge/TypedLinearAlgebra/blob/master/README.md#plug-ins).

### Include

Include the plug-in's facade, which provides the `matrix`, `column_vector`, and `row_vector` aliases of the [dense](#patterns) pattern:

```cpp
#include "fcarouge/linalg.hpp"
```

For the storage and the `sparse_matrix` aliases only, include instead:

```cpp
#include "fcarouge/sparse.hpp"
```

### Namespace

The `sparse_matrix`, `sparse_column_vector`, and `sparse_row_vector` aliases live in the `fcarouge` namespace. The storage, the patterns, and `make_matrix` live in the `fcarouge::sparse` namespace.

### Plug-ins

| Element Types | CMake Target | Facade |
| --- | --- | --- |
| Built-in Types | `typed_linear_algebra_sparse` | `support/sparse/fcarouge/linalg.hpp` |
| Au | `typed_linear_algebra_au_sparse` | `support/au_sparse/fcarouge/linalg.hpp` |
| mp-units | `typed_linear_algebra_mp_units_sparse` | `support/mp_units_sparse/fcarouge/linalg.hpp` |
| nholthaus/units | `typed_linear_algebra_nholthaus_sparse` | `support/nholthaus_sparse/fcarouge/linalg.hpp` |
| std::chrono | `typed_linear_algebra_chrono_sparse` | `support/chrono_sparse/fcarouge/linalg.hpp` |

### Samples

| Interoperability | Sample |
| --- | --- |
| mp-units + Sparse | [Usage example: `sample/mp_units_sparse.cpp`](https://github.com/FrancoisCarouge/TypedLinearAlgebra/blob/master/sample/mp_units_sparse.cpp). |

### Aliases

Strongly typed matrix and vectors of a compile-time sparsity pattern. Each one is a [`typed_matrix`](https://github.com/FrancoisCarouge/TypedLinearAlgebra/blob/master/README.md#class-typed-matrix) composing a [`sparse::matrix`](#class-matrix), with all of its member types, member functions, and operations, subject to [Element Access](#element-access).

```cpp
template <typename Representation, sparse::sparsity_pattern auto Pattern,
          typename RowIndexes = /* Pattern.rows of Representation */,
          typename ColumnIndexes = /* Pattern.columns of Representation */>
using sparse_matrix;

template <typename Representation, sparse::sparsity_pattern auto Pattern,
          typename... Types>
using sparse_column_vector;

template <typename Representation, sparse::sparsity_pattern auto Pattern,
          typename... Types>
using sparse_row_vector;
```

| Template Parameter | Definition |
| --- | --- |
| `Representation` | The type of the stored elements. |
| `Pattern` | The sparsity pattern, of the shape of the indexes. |
| `RowIndexes` | The tuple type of the row indexes. |
| `ColumnIndexes` | The tuple type of the column indexes. |
| `Types` | The types of the elements of the vector, one for each of its rows or columns. |

A pattern of another shape than the indexes does not satisfy the constraints of the alias.

### Class Matrix

Compile-time sparsity pattern matrix storage. Stores the elements of the pattern in place, in row-major order.

Also documented in the [fcarouge/sparse.hpp](https://github.com/FrancoisCarouge/TypedLinearAlgebra/blob/master/support/sparse/fcarouge/sparse.hpp) header.

#### Declaration

```cpp
namespace sparse {
template <typename Type, sparsity_pattern auto Pattern>
struct matrix;
}
```

#### Template Parameters

| Template Parameter | Definition |
| --- | --- |
| `Type` | The type of the stored elements. |
| `Pattern` | The sparsity pattern of the matrix. |

#### Member Variables

| Member Variable | Definition |
| --- | --- |
| `structure` | Static. The sparsity pattern of the matrix. |
| `values` | The `std::array` of the stored elements, zero-initialized, in row-major order. |

#### Member Functions

The [operations](#operations) on the composing typed matrix are the preferred interface. The vector, norm, and sum members below serve them.

| Member Function | Definition |
| --- | --- |
| `(default constructor)` | Construct a zero matrix. |
| `(explicit constructor)` | Construct from every stored element, in row-major order. Each one is exactly of `Type`, and there is at least one. |
| `(conversion constructor)` | Implicitly convert from a matrix of the same element type and shape, whose stored positions this pattern all stores. |
| `operator[row, column] const` | Read the specified element by value. A structural zero reads as `Type{}`. |
| `operator[row, column]` | Access the specified element by reference. Only without a structural zero. |
| `transpose()` | The transposed matrix, of the transposed pattern. |
| `dot(other)` | The sum of the products of the elements stored by both vectors, of the same size and of any orientation. |
| `norm()` | The Frobenius norm. Only for a floating-point `Type`. |
| `cwiseAbs()` | The element-wise absolute values, of the same pattern. |
| `colwise()` | The sums of the columns. |
| `rowwise()` | The sums of the rows. |

### Class Pattern

The positions of the stored elements of a matrix. A structural type, usable as a template argument, of strictly ascending row-major positions: two equal patterns are one and the same type of matrix. Prefer the [named patterns and factories](#patterns), which establish that order.

#### Declaration

```cpp
namespace sparse {
struct position;

inline constexpr std::size_t absent;

template <std::size_t Rows, std::size_t Columns, std::size_t Count>
struct pattern;
}
```

#### Template Parameters

| Template Parameter | Definition |
| --- | --- |
| `Rows` | The count of rows of the matrix. |
| `Columns` | The count of columns of the matrix. |
| `Count` | The count of stored elements of the matrix. |

#### Member Variables

| Member Variable | Definition |
| --- | --- |
| `rows` | The count of rows. |
| `columns` | The count of columns. |
| `count` | The count of stored elements. |
| `positions` | The `std::array` of the stored positions, in row-major order. A `position` has a `row` and a `column`. |

#### Member Functions

| Member Function | Definition |
| --- | --- |
| `lower_bound(position)` | The index of the first stored position not before the given one. |
| `index_of(row, column)` | The index of the stored element at the given position, or `absent` for a structural zero. |
| `stores(row, column)` | Whether the given position is a stored element. |
| `operator==` | Whether two patterns, of any sizes, store the same positions of the same shape. |

### Patterns

```cpp
namespace sparse {
template <std::size_t Rows, std::size_t Columns>
consteval auto make_pattern();

template <std::size_t Rows, std::size_t Columns, std::size_t Count>
consteval auto make_pattern(const position (&positions)[Count]);

template <std::size_t Rows, std::size_t Columns, auto Predicate>
consteval auto generate();
}
```

| Pattern | Definition |
| --- | --- |
| `make_pattern<R, C>()` | No position. |
| `make_pattern<R, C>({{i, j}, ...})` | The listed positions, in any order. A duplicate or out-of-range position does not compile. |
| `generate<R, C, predicate>()` | The positions satisfying a stateless predicate of a row and a column. |
| `dense<R, C>` | Every position: a matrix without structural zeros. |
| `zero<R, C>` | No position: the zero matrix. |
| `banded<R, C, L, U>` | The `L` sub-diagonals, the main diagonal, and the `U` super-diagonals. |
| `diagonal<N>` | The main diagonal of a square matrix. |
| `tridiagonal<N>` | The three central diagonals of a square matrix. |
| `lower_triangular<N>` | The main diagonal and below of a square matrix. |
| `upper_triangular<N>` | The main diagonal and above of a square matrix. |

### Pattern Algebra

The compile-time computations the operations derive their results with, each a variable template of the `fcarouge::sparse` namespace.

| Variable | Definition |
| --- | --- |
| `pattern_union<Lhs, Rhs>` | The positions stored by either pattern of the same shape: the pattern of a sum or a difference. |
| `pattern_intersection<Lhs, Rhs>` | The positions stored by both patterns of the same shape. |
| `pattern_transpose<Pattern>` | The transposed pattern. |
| `pattern_product<Lhs, Rhs>` | The structural pattern of a matrix product: a position is stored whenever a term reaches it. |
| `product_term_count<Lhs, Rhs>` | The count of elementary products of a matrix product. |
| `is_subset_of<Inner, Outer>` | Whether a pattern stores every position of another, of the same shape. |

### Element Construction

#### `sparse::make_matrix`

```cpp
namespace sparse {
struct structural_zero_t;

inline constexpr structural_zero_t structural_zero;

template <typename TypedMatrix, typename... Values>
[[nodiscard]] constexpr auto make_matrix(const Values &...values)
    -> TypedMatrix;
}
```

Construct a typed matrix from a value at every position, in row-major order. Does not participate in overload resolution unless every position has a value, every stored element's value is implicitly convertible to its element's type, and every structural zero's value is `structural_zero`: a `0.` at a structural zero does not compile, nor does `structural_zero` at a stored element.

| Parameter | Definition |
| --- | --- |
| `TypedMatrix` | The typed matrix over a sparse storage to construct. |
| `values` | The value of each position, in row-major order. |

### Element Access

A structural zero has nothing to refer to. On a pattern with a structural zero, a read is by value, and a write that would be lost does not compile, whatever the value. The position of `at<i, j>(value)` is known at compile time, so it writes a stored element. A pattern without a structural zero is an ordinary dense storage.

| Typed Matrix Access | Without Structural Zero | With Structural Zero |
| --- | --- | --- |
| Read: `at<i, j>()`, `at<Type>()`, `operator[i, j]`, `operator(i, j)`, structured bindings | Yes. | Yes, by value. A structural zero reads as a zero of the element's type. |
| Write: `at<i, j>(value)` | Yes. | Stored positions only. |
| Construction from a list, a tuple-like vector, or a singleton value | Yes. | Does not compile. |
| Construction from a compatible typed matrix | Yes. | Yes, when the pattern stores each of its stored positions, and a conversion between units keeps a zero. |
| Construction from a `sparse::matrix`, in stored order | Yes. | Yes. |
| `sparse::make_matrix<Matrix>(values...)` | Yes. | Yes, with `structural_zero` at each structural zero. |

The rejections of the member write and of the constructions are `static_assert` failures, which a `requires` expression cannot detect. To ask whether a position can be written, use `decltype(matrix)::matrix::structure.stores(i, j)`. The rejections of `sparse::make_matrix` are constraints, which it can detect.

### Operations

The [operations](https://github.com/FrancoisCarouge/TypedLinearAlgebra/blob/master/README.md#operations) of the typed matrix below are supported. `K` is the count of stored elements, `R` and `C` the counts of rows and columns, and `T` the [`product_term_count`](#pattern-algebra). Each result pattern is computed at compile time.

| Operation | Result Pattern | Run-Time Complexity |
| --- | --- | --- |
| `+`, binary `-` | The union of the patterns. | `O(K)` of the result. |
| Unary `-` | The same pattern. | `O(K)`. |
| `*`, of matrices and vectors | The structural product of the patterns. | `O(K + T)` of the result. |
| `*`, `/`, by a scalar | The same pattern. | `O(K)`. |
| `/`, of singletons | A singleton. | `O(1)`. |
| `==` | Not applicable. Compares every position: a structural zero equals a stored zero. | `O(R * C * log K)`. |
| `dot` | Not applicable. Returns an element. | `O(K)`, of the elements stored by both vectors. |
| `magnitude`, `matrix_frob_norm` | Not applicable. Returns an element. Floating-point elements only. | `O(K)`. |
| `matrix_one_norm` | Not applicable. Returns an element. | `O(K + C)`. |
| `matrix_inf_norm` | Not applicable. Returns an element. | `O(K + R)`. |
| `transposed` | The transposed pattern. | `O(K)`. |
| `at<i, j>(value)` | Not applicable. Writes one element. | `O(1)`. |
| `sparse::make_matrix` | The pattern of the matrix to construct. | `O(K)`. |
| Read of an element | Not applicable. Returns an element. | `O(log K)`. |

The following are not supported yet:

- Matrix division, other than of singletons, which needs a solver, a factorization.
- The `std::linalg` out-parameter algorithms: `add`, `matrix_product`, `matrix_rank_1_update`, `matrix_vector_product`, and `scale`, which need `std::mdspan` storage.

### Concepts

| Concept | Definition |
| --- | --- |
| `sparse::sparsity_pattern` | Concept of a compile-time sparsity pattern. |

### Format

The typed matrix formats as for any other backend, structural zeros included, so the format does not depend on the pattern. Hosted implementation only.

```cpp
std::println("{}", value); // [[1, 0], [4, 0]]
```

## Notes

- A pattern is structural: elements cancelling out at run time stay stored.
- A structural zero stays an exact zero, for which no term exists. Unlike a dense computation, a non-finite operand does not propagate through it: `0 * infinity` is never computed, and is not a not-a-number value.
- A matrix converts implicitly to a pattern storing all of its positions, and never to a narrower one, which could lose an element.
- A conversion between units keeps each element's value, a structural zero included, when the conversion keeps a zero: seconds from minutes, meters from kilometers. A conversion with an offset, as kelvin from degrees Celsius, does not compile into a pattern with a structural zero, which cannot hold the converted zero, 273.15 kelvin.
- `sizeof` of a matrix with stored elements is that of its `std::array` of stored elements. A matrix of an empty pattern still occupies the storage of an empty `std::array`, which is not guaranteed to be nothing.
- The storage is freestanding. On GCC and Clang, the `typed_linear_algebra_sparse_freestanding` target verifies it with `-ffreestanding -fno-exceptions -fno-rtti`.
- At compile time, and on a freestanding implementation, the Frobenius norm and `magnitude` use their own square root, whose result is not guaranteed to be identical to the last digit to that of `std::sqrt`.
- The compile-time cost suits tens to low hundreds of stored elements. The pattern operations are near-linear in the stored elements and verified at 100 x 100, but the typed `*` also verifies each of its R x C x K element combinations: verified at 32 x 32, while 64 x 64 exceeds Clang's default constant evaluation limit.
- With Au, an unqualified `matrix * matrix` may resolve to the hidden-friend `operator*` of Au's quantity. Call `fcarouge::operator*(lhs, rhs)` explicitly.

## Example

A constant-velocity state transition, of strongly typed elements with mp-units, linked against `typed_linear_algebra_mp_units_sparse`.

```cpp
#include "fcarouge/linalg.hpp"

#include <print>
#include <tuple>

namespace {
using namespace fcarouge;
using namespace mp_units;
using namespace mp_units::si::unit_symbols;

using position = quantity<isq::length[m], double>;
using velocity = quantity<isq::speed[m / s], double>;
using inverse_position = quantity<one / isq::length[m], double>;
using inverse_velocity = quantity<one / isq::speed[m / s], double>;

// The bottom left element is a structural zero, without storage.
using transition =
    sparse_matrix<double, sparse::upper_triangular<2>,
                  std::tuple<position, velocity>,
                  std::tuple<inverse_position, inverse_velocity>>;

using state =
    sparse_column_vector<double, sparse::dense<2, 1>, position, velocity>;

static_assert(sizeof(transition) == 3 * sizeof(double));
} // namespace

auto main() -> int {
  // Each stored element is written with its own type: a ratio, a duration.
  transition model;

  model.at<0, 0>(1. * one);
  model.at<0, 1>(2. * s);
  model.at<1, 1>(1. * one);

  // Does not compile, the value would be lost:
  // model.at<1, 0>(1. / s);

  // A structural zero reads with its own type, an inverse duration.
  std::println("Structural zero: {}", model.at<1, 0>());

  // The product only computes the terms the patterns admit.
  const state initial{1. * m, 3. * m / s};
  const state predicted{model * initial};

  std::println("Position: {}", predicted.at<0>());
  std::println("Velocity: {}", predicted.at<1>());
}
```

Output:

```text
Structural zero: 0 1/s
Position: 7 m
Velocity: 3 m/s
```

## See Also

| Reference | Description |
| --- | --- |
| [Typed Linear Algebra](https://github.com/FrancoisCarouge/TypedLinearAlgebra/blob/master/README.md) | The typed matrix, its operations, and the other plug-ins. |
| [`support/freestanding`](https://github.com/FrancoisCarouge/TypedLinearAlgebra/tree/master/support/freestanding) | A fixed-size dense matrix over a `std::array`, relying only on the freestanding standard library. |
| [`test/sparse`](https://github.com/FrancoisCarouge/TypedLinearAlgebra/tree/master/test/sparse) | The verification of this plug-in, including the rejections above. |
