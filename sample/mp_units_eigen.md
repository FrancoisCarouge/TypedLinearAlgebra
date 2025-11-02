# Quantity-Safe Linear Algebra Use Case: Eigen + mp-units

François Carouge (francois.carouge@gmail.com)

October 11, 2025

## Abstract

This work presents a practical approach to achieving quantity-safe linear algebra in C++ through the composition of the [Eigen](https://eigen.tuxfamily.org/) and [mp-units](https://mpusz.github.io/mp-units) libraries, enabled by the [FrancoisCarouge/TypedLinearAlgebra](https://github.com/FrancoisCarouge/TypedLinearAlgebra) library. The proposed method integrates dimensional analysis directly into linear algebra computations, ensuring compile-time enforcement of unit correctness while preserving the efficiency and flexibility of established numerical backends. By separating type semantics from data storage, the design maintains compatibility with Eigen’s expression templates while extending type safety to operations on physical quantities. The paper documents the library’s design principles, implementation strategies, and accumulated experience, including trade-offs in type representation and performance considerations. Lessons learned highlight the feasibility and challenges of embedding strong typing into numerical computation, while preliminary extensions demonstrate the applicability of the approach beyond physical units. This contribution offers both a demonstration of safe, dimensionally consistent linear algebra and a foundation for extending type-aware computation into broader C++ scientific and engineering domains.

## Introduction

Scientific and engineering applications frequently rely on linear algebra to model, simulate, and analyze physical systems. In these domains, correctness does not stop at numerical accuracy: operations must also respect the physical dimensions of the underlying quantities. While dimensional analysis provides a systematic way to ensure that equations remain physically valid, traditional linear algebra libraries in C++ treat all scalar values as interchangeable numeric types. As a result, mistakes such as adding a velocity to a length or multiplying incompatible quantities can go undetected.

C++ offers strong typing, but its native type system does not directly enforce unit correctness in mathematical operations. To address this gap, libraries such as [mp-units](https://mpusz.github.io/mp-units) have emerged, enabling compile-time dimensional analysis, unit- and character-safe quantities. Separately, numerical libraries such as [Eigen](https://eigen.tuxfamily.org/) provide efficient and expressive support for linear algebra through expression templates and optimized storage backends. However, these two capabilities—unit safety and linear algebra—have traditionally been developed in isolation, leaving developers without a unified approach to both.

This paper explores how the [FrancoisCarouge/TypedLinearAlgebra](https://github.com/FrancoisCarouge/TypedLinearAlgebra) library can bridge that gap by composing Eigen and mp-units behind a facade for quantity-safe linear algebra. The approach separates type semantics from storage representation: physical units and dimensional correctness are captured at the type level, while Eigen retains responsibility for numerical efficiency and data layout. In this way, dimensionally safe operations can be expressed naturally while leveraging the performance optimizations of existing linear algebra infrastructure.

We describe the design principles of the library, illustrate the use of type-safe matrices in practice, and discuss the challenges encountered during implementation. In particular, we highlight trade-offs in type representation, template complexity, and expression template interactions. Finally, we present lessons learned from applying quantity-safe linear algebra in real-world scenarios and outline opportunities for extending this work to broader semantic domains beyond physical units.

The remainder of the paper is organized as follows. Section 2 provides background on type safety, linear algebra, and dimensional analysis in C++. Section 3 compares three approaches to typed linear algebra. Section 4 presents case studies that demonstrate the framework in practice. Section 5 documnts the library. Section 6 outlines opportunities for extending type safety beyond physical quantities and standard linear algebra. Section 7 concludes with reflections on the contribution and potential directions for future work.

## 2 Type safety, Linear Algebra, and Dimensional Analysis

Type safety is a concept in programming that aims to prevent type errors during program compilation or execution. Type safety ensures that operations are performed on compatible data types, thus avoiding unexpected behavior and improving code reliability. A strongly typed programming language is one in which types are enforced strictly, meaning the language does not allow operations or conversions between mismatched types without explicit handling by the programmer. C++ is mostly [strongly typed](https://en.cppreference.com/w/cpp/language/type-id.html). The strong typing in C++ promotes correctness, safety, readability, maintainability, and performance.

Linear algebra is a branch of mathematics that studies vectors, matrices, and linear transformations. Linear algebra is a foundational tool in many fields including physics, engineering, computer graphics, and machine learning. Linear algebra is used to solve systems of linear equations, analyze data, and understand geometric concepts like rotations and projections. Standard C++ offers [linear algebra algorithms](https://en.cppreference.com/w/cpp/numeric/linalg.html) `<linalg>` since C++26 from the [Basic Linear Algebra Subprograms (BLAS)](https://www.netlib.org/blas/). There also exists other linear algebra packages: [Armadillo](https://arma.sourceforge.net), [OpenGL Mathematics (GLM)](https://github.com/g-truc/glm), [Kokkos Kernels](https://github.com/kokkos/kokkos-kernels), BLAS/LAPACK, [Portable, Extensible Toolkit for Scientific Computation (PETSc)](https://petsc.org/), etc. [Eigen](https://eigen.tuxfamily.org/) is a C++ library for linear algebra: matrices, vectors, numerical solvers, and related algorithms. Eigen uses the [expression template](https://en.wikibooks.org/wiki/More_C%2B%2B_Idioms/Expression-template) technique. 

[Dimensional analysis](https://en.wikipedia.org/wiki/Dimensional_analysis) is a method used in physics, engineering, and other sciences to understand the relationships between different physical quantities by analyzing their dimensions (such as length [L], mass [M], time [T], etc.). It involves expressing physical quantities in terms of their basic dimensions and using these to check equations, derive relationships, or convert units. Dimensional analysis in C++ ensures physical correctness by enforcing unit safety at compile time. Dimensional analysis makes C++ code more reliable, readable, and maintainable when working with physical quantities or scientific computations. C++ offers some [time quantities](https://en.cppreference.com/w/cpp/chrono.html) with `<chrono>` since C++11. There also exists other units packages: [nholthaus/units](https://nholthaus.github.io/units/), [aurora-opensource/au](https://github.com/aurora-opensource/au), [Boost.Units](https://www.boost.org/library/latest/units), etc. [Mp-units](https://mpusz.github.io/mp-units) provides compile-time dimensional analysis and unit/quantity manipulation.

Linear algebra, dimensional analysis, and C++ intersect in real-world domains: robotics, autonomous vehicles, aerospace, structural engineering, medical engineering, scientific simulation, computational fluid dynamics, etc. The C++ implementers expect linear algebra and dimensional analysis to be compatible and also expect type safeties to propagate. C++ implementers seeks linear algebra and unit support [[1](https://stackoverflow.com/questions/8120126/combining-a-linear-algebra-library-with-boostunits)]:

> Is there any matrix library that supports units? *-2011 Stack Overflow user*

The [Buckingham π theorem](https://en.wikipedia.org/wiki/Buckingham_%CF%80_theorem) states the equation can be rewritten in terms of a set of dimensionless parameters. The nondimensionalization produces natural scales. The theorem can be used to systematically remove quantities from computations. Quantities still have to wrap around the computations. The Buckingham π theorem is not always applied. Strongly typed linbear algebra supports implementers.

## 3 Typed Linear Algebra

An approach in designing typed linear algebra consists of the injection of strong types into existing linear algebra libraries [[2](https://eigen.tuxfamily.org/dox/TopicCustomizing_CustomScalar.html)][[3](https://github.com/blitzpp/blitz)][[4](https://stackoverflow.com/questions/62040237/eigen-with-custom-scalar-types-matrix-multiplication-with-custom-type-fails-wit)]. Non-uniformly typed vectors or matrices are challenging: each element may have a different type and the libraries not typically equiped with a supporting type machinery or its requirements. Then the results of operations require specialization, customization points, so that multiplying two elements in a unit, computes to the square unit, while still represented by an underlying floating point type: $length * length = L^2 = area$ while for double precision floating point data type $double * double = double$.

Yet another approach in designing typed linear algebra consists of tracking every elements types. Non-uniformly typed matrices may present with up to `m x n` types. Consequently, large matrices are challenging: the number of template parameters grows quadratically.

In this library approach in designing typed linear algebra consists of not tracking every single element types. Hall explains Hart's conjecture that all `m x n` usable matrices can encode all their element dimensions in just two m-, n-vectors [5], that is `m + n` types. Daniel Withopf [[6](https://www.youtube.com/watch?v=4LmMwhM8ODI)] and Chip Hogg [[7](https://www.youtube.com/watch?v=5dhFtSu3wCo)] used the conjecture to present a typed matrix composing a traditional linear algebra matrix with strong types. The typed matrix aggregates three template parameters: the representation matrix, the row indexes, and the column indexes. Note the matrix template parameter may support the expression template facilities of the underlying representation matrix with its compiler performance optmizations and abstract expression temporary pitfalls. There exists multiple valid pairs of row index and column index to obtain one typed matrix element. Hart's conjecture is useful to reduce the quantity of stored indexes, however operations still need to verify the combination of the row and column indexes are compatible.

## 4 Case Study

```cpp
// 1-D vehicle location Kalman estimation.
using state = column_vector<position, velocity, acceleration>;
state x{0. * m, 0. * m / s, 0. * m / s2};
std::println("X: {}", x);
// X: [[0 m],
//     [0 m/s],
//     [0 m/s²]]

using estimate_uncertainty =
    matrix<std::tuple<position, velocity, acceleration>,
            std::tuple<position, velocity, acceleration>>;
estimate_uncertainty p;
p.at<0, 0>() = 500. * m2;
p.at<1, 1>() = 500. * m2 / s2;
p.at<2, 2>() = 500. * m2 / s4;
std::println("P: {}", p);
// P: [[500 m²,     0 m²/s,    0 m²/s²],
//     [  0 m²/s, 500 m²/s²,   0 m²/s³],
//     [  0 m²/s²,  0 m²/s³, 500 m²/s⁴]]

  using process_uncertainty = estimate_uncertainty;
  process_uncertainty q;
  q.at<0, 0>() = 0.01 * m2;
  q.at<0, 1>() = 0.02 * m2 / s;
  q.at<0, 2>() = 0.02 * m2 / s2;
  q.at<1, 0>() = 0.02 * m2 / s;
  q.at<1, 1>() = 0.04 * m2 / s2;
  q.at<1, 2>() = 0.04 * m2 / s2;
  q.at<2, 0>() = 0.02 * m2 / s2;
  q.at<2, 1>() = 0.04 * m2 / s3;
  q.at<2, 2>() = 0.04 * m2 / s4;
  std::println("Q: {}", q);
  // Q: [[0.01 m²,    0.02 m²/s,  0.02 m²/s²],
  //     [0.02 m²/s,  0.04 m²/s², 0.04 m²/s³],
  //     [0.02 m²/s², 0.04 m²/s³, 0.04 m²/s⁴]]

using output_uncertainty = quantity<m2>;
output_uncertainty r{9. * m2};
std::println("R: {}", r);
// R: 9 m²

using output_model = row_vector<quantity<one>, quantity<s>, quantity<s2>>;
// output_model h{1., 0., 0.}; // WHY NOT FAILING COMPIL?
output_model h{output_model::matrix::Identity()};
std::println("H: {}", h);
// H: [1, 0 s, 0 s²]

using state_transition =
    matrix<std::tuple<position, velocity, acceleration>,
            std::tuple<quantity<one / m>, quantity<s / m>, quantity<s2 / m>>>;
state_transition f{state_transition::matrix::Identity()};
f.at<0, 1>() = 1. * s;
f.at<0, 2>() = 0.5 * s2;
f.at<1, 2>() = 1. * s;
std::println("F: {}", f);
// F: [[1, 1 s, 0.5 s²],
//     [0 1/s, 1, 1 s],
//     [0 1/s², 0 1/s, 1]]

// Predict
x = f * x;
p = f * p * transposed(f) + q;

// Update
using output = position;
output z{-393.66 * m};

using innovation_uncertainty = output_uncertainty;
innovation_uncertainty si{h * p * transposed(h) + r};

using unevaluated_gain =
    decltype(std::declval<state>() / std::declval<output>());
using gain =
    matrix<unevaluated_gain::row_indexes, unevaluated_gain::column_indexes>;
gain k{p * transposed(h) / si};

using innovation = output;
innovation y{z - h * x};
x = x + k * y;

std::println("X: {}", x);
// X: [[-390.53 m],
//     [-260.36 m/s],
//     [ -86.79 m/s²]]

using unevaluated_kh =
    decltype(std::declval<gain>() * std::declval<output_model>());
using kh =
    matrix<unevaluated_kh::row_indexes, unevaluated_kh::column_indexes>;
kh i{state_transition::matrix::Identity()};
p = (i - k * h) * p * transposed(i - k * h) + k * r * transposed(k);
std::println("P: {}", p);
// P: [[8.92 m²,      5.95 m²/s,    1.98 m²/s²],
//     [5.95 m²/s,  503.98 m²/s², 334.67 m²/s³],
//     [1.98 m²/s², 334.67 m²/s³, 444.91 m²/s⁴]]
```

## 5 Library Reference

### Class Typed Matrix

Strongly typed matrix. Compose a linear algebra backend matrix into a typed matrix. Row and column indexes provide each element's index type.

Also documented in the [fcarouge/typed_linear_algebra.hpp](https://github.com/FrancoisCarouge/TypedLinearAlgebra/blob/master/include/fcarouge/typed_linear_algebra.hpp) header.

#### Declaration

```cpp
template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
class typed_matrix
```

#### Template Parameters

| Template Parameter | Definition |
| --- | --- |
| `Matrix` | The underlying linear algebra matrix. |
| `RowIndexes` | The tuple type of the row indexes. |
| `ColumnIndexes` | The tuple type of the row indexes. |

#### Member Types

| Member Type | Definition |
| --- | --- |
| `matrix` | The type of the composed matrix. |
| `underlying` | The type of the element's underlying storage. |
| `row_indexes` | The tuple with the row components of the indexes. |
| `column_indexes` | The tuple with the column components of the indexes. |
| `element<i, j>` | The type of the element at the given matrix indexes position. |

#### Member Variables

| Member Variable | Definition |
| --- | --- |
| `rows` | The count of rows. |
| `columns` | The count of columns. |

#### Member Functions

| Member Function | Definition |
| --- | --- |
| `(default constructor)` | Construct a default typed matrix. |
| `(default copy constructor)` | Copy construct the typed matrix. |
| `(default copy assignment operator)` | Copy assign a typed matrix. |
| `(default move constructor)` | Move construct a typed matrix. |
| `(default move assignment operator)` | Move construct a typed matrix. |
| `(conversion copy constructor)` | Copy construct generalization of a compatible typed matrix. |
| `(conversion copy assignment operator)` | Copy assign generalization of a compatible typed matrix. |
| `(conversion move constructor)` | Move construct generalization of a compatible typed matrix. |
| `(conversion move assignment operator)` | Move assign generalization of a compatible typed matrix. |
| `(conversion copy constructor)` | Convert construct a typed matrix from an underlying matrix. |
| `(conversion copy constructor)` | Convert construct a one-dimension uniformly typed matrix from array. |
| `(conversion copy constructor)` | Convert construct a uniformly typed matrix from list-initializers. |
| `(conversion copy constructor)` | Convert construct a row or column typed vector from elements. |
| `(conversion copy constructor)` | Convert construct a singleton typed matrix from a single value. |
| `operator[i, j]` | Access the specified element. |
| `operator(i, j)` | Access the specified element. |
| `at<i, j>()` | Access the specified element. |
| `(conversion operator)` | Access the singleton typed matrix element. |
| `(destructor)` | Destruct a default typed matrix. |

### Operations

The following useful operations are supported. This library attempts to align its nomenclature aligned with that of the primitives provided by `std::linalg`. This library attempts some compatibility with other C++ standard library primitives, ranges, and iterators.

| Operation | Definition |
| --- | --- |
| `+` | Addition where the terms are of identical shapes and addable types. |
| `-` | Substraction where the terms are of identical shapes and substractable types. |
| `*` | Multiplication where the factors are of multipliable shapes and multipliable types. |
| `/` | Solution, if there exists one, to the inverse multiplication, where the factor are of compatible shapes and types. |
| `==` | Direct, strict equality comparison, with traditional floating-point comparison pitfalls. |
| `transposed` | Transpose the input matrix. |

### Aliases

```cpp
template <typename Matrix, typename... ColumnIndexes>
typed_row_vector;

template <typename Matrix, typename... RowIndexes>
typed_column_vector;
```

### Format

A specialization of the standard formatter is provided for the typed matrix. Use `std::format` to store a formatted representation of the matrix. Standard format parameters to be supported.

### Concepts

| Concept | Definition |
| --- | --- |
| `same_as_typed_matrix` | Concept of a typed matrix type. |
| `singleton_typed_matrix` | Concept of a singleton, one-element typed matrix type. |
| `uniform_typed_matrix` | Concept of a typed matrix in which all element types are the same. |
| `one_dimension_typed_matrix` | Concept of a typed matrix with only one dimension, row, or column. |
| `row_typed_matrix` | Concept of a row typed matrix, vector. |
| `column_typed_matrix` | Concept of a column typed matrix, vector. |
| `same_shape` | Concept of typed matrices of the same shape, that is they have the same number of rows and columns. |

### Structure Element Caster

Typed matrix element conversions customization point. Specialize this template to allow conversion of element's type and underlying type.

```cpp
template <typename To, typename From>
struct element_caster;
```

## 6 Open Questions

**Index-Type Safety Required:** Type safety cannot be guaranteed at compilation time without index-type safety. The indexes can either be non-type template parameters or strong types overloadings. Converting a runtime index to a dependent template type is not possible in C++23. A proxy reference could be used to allow traditional assignment syntax but the runtime check and extra indirection are not interesting tradeoffs. A template call operator can be used for getting a type safe value but impractical syntax to set values. Without index safety, the accepted tradeoff is a templated index `at<i, j>()` method over the unchecked `operator[i, j]` and `operator(i, j)` accessors. Daniel Withopf [[6](https://www.youtube.com/watch?v=4LmMwhM8ODI)] uses strongly typed indexed to provide index safety complementing the quantity-safe linear algebra approach. A strongly-typed linear algebra must include index-type safety. Other types of safeties should be considered: reference frames, coordinate systems, or taxonomy of the matrix.

**Zero-Cost Abstraction Required:** The performance of the typed linear algebra should be identical to the performance of the underlying linear algebra backend, considering the types. The composition is intended to be a zero-cost abstraction. Measuring and comparing the library to their equivalent counterpart is left as a future exercise.

**Type Conversions:** The conversion from the underlying matrix element type references to the quantity type references is achieved through a `reinterpret_cast` conversion which works by reinterpreting the underlying bit pattern forgoing aliasing and alignement rules. The safety of the `element_caster` implementation is dubious and safer alternatives to be identified. The quantity type library only provide partial type conversion to and from standard types. It may be appropriate for the type library to extend its safe conversion support: to/from, by-value/reference, not/constant, and l/rvalue. A single template specialization could then suffice to customize the type conversions.

**No Unsafe API:** Three methods of the matrix class present risks: the default constructor, the conversion constructor for underlying matrix, and the underlying data access method:
* The default constructor could construct the underlying matrix according to its default linear algebra backend initialization, or lack thereof. Backends without guaranteed initialization have shown to result in uninitialized memory defects. The default constructor should be safe by default under all corner cases, or be gone. Adding an explicitely uninitialized constructor could be considered.
* One of the conversion constructor permits the construction of a typed matrix from its untyped, underlying matrix data. Although practical for the user and operation implementations, the constructor however defeats the purpose of the type safety and re-introduces the risk the typed matrix set out to resolve. Additional contruction strategies should help the user in constructing type safe matrices without using an unsafe construction methodology.
* The underlying data access method permits accessing the untyped underlying matrix data of a typed matrix. Although practical for the user and operation implementations, the method however defeats the purpose of the type safety and re-introduces the risk the typed matrix set out to resolve. Privatization, friendship of operations, ADL may prove useful in resolving this issue.

## 7 Beyond Unit Safety

There are more operations used in linear algebra than the ones presented in this document. More or all operations from the Eigen or the `std::linalg` would form a practical operation catalog. Identifying the most commonly used operations are implementing their facades are left as future exercises. Similarly for the exploration and integration with compatible standard library support such as algorithms, ranges, or mdspan.


An interesting research venue would be to elevate, generalize the dimension of the matrices beyond the rank-2 matrices present in common linear algebra. Replacing the row- and column-index tuples by a pack of index tuples is left for a future exploration.

## 7 Conclusion

This paper has shown that combining Eigen with mp-units through the TypedLinearAlgebra library enables quantity-safe linear algebra in C++ without sacrificing usability. By enforcing dimensional correctness at compile time while relying on Eigen for numerical efficiency, the approach makes physically consistent computation practical in real-world applications.

The results confirm both the feasibility and the limitations of this design, particularly regarding template complexity, index safety, and the scope of supported operations. Future work may extend these ideas beyond units to encompass coordinate systems, reference frames, and higher-rank tensors, broadening the role of type-aware computation in scientific and engineering software.
