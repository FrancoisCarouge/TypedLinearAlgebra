# Quantity-Safe Linear Algebra Use Case: Eigen + mp-units

This article showcases the use of the [FrancoisCarouge/TypedLinearAlgebra](https://github.com/FrancoisCarouge/TypedLinearAlgebra) to compose [Eigen](https://eigen.tuxfamily.org/) and [mp-units](https://mpusz.github.io/mp-units) in support of quantity-safe linear algebra.

## Type safety, Linear Algebra, and Dimensional Analysis

Type safety is a concept in programming that aims to prevent type errors during program compilation or execution. Type safety ensures that operations are performed on compatible data types, thus avoiding unexpected behavior and improving code reliability. A strongly typed programming language is one in which types are enforced strictly, meaning the language does not allow operations or conversions between mismatched types without explicit handling by the programmer. C++ is mostly [strongly typed](https://en.cppreference.com/w/cpp/language/type-id.html). The strong typing in C++ promotes correctness, safety, readability, maintainability, and performance.

Linear algebra is a branch of mathematics that studies vectors, matrices, and linear transformations. Linear algebra is a foundational tool in many fields, including physics, engineering, computer graphics, and machine learning. Linear algebra is used to solve systems of linear equations, analyze data, and understand geometric concepts like rotations and projections. Standard C++ offers basic [linear algebra algorithms](https://en.cppreference.com/w/cpp/numeric/linalg.html) `<linalg>` since C++26 from the [Basic Linear Algebra Subprograms (BLAS)](https://www.netlib.org/blas/). There also exists other linear algebra packages: [Armadillo](https://arma.sourceforge.net), [OpenGL Mathematics (GLM)](https://github.com/g-truc/glm), BLAS/LAPACK, [Portable, Extensible Toolkit for Scientific Computation (PETSc)](https://petsc.org/), etc. [Eigen](https://eigen.tuxfamily.org/) is a C++ template library for linear algebra: matrices, vectors, numerical solvers, and related algorithms.

Dimensional analysis is a method used in physics, engineering, and other sciences to understand the relationships between different physical quantities by analyzing their dimensions (such as length [L], mass [M], time [T], etc.). It involves expressing physical quantities in terms of their basic dimensions and using these to check equations, derive relationships, or convert units. Dimensional analysis in C++ ensures physical correctness by enforcing unit safety at compile time. Dimensional analysis makes C++ code more reliable, readable, and maintainable when working with physical quantities or scientific computations. C++ offers some [time quantities](https://en.cppreference.com/w/cpp/chrono.html) with `<chrono>` since C++11. There also exists other units packages: [nholthaus/units](https://nholthaus.github.io/units/), [aurora-opensource/au](https://github.com/aurora-opensource/au), [Boost.Units](https://www.boost.org/library/latest/units), etc. [mp-units](https://mpusz.github.io/mp-units) is a compile-time enabled feature-rich Modern C++ modular/header-only library that provides compile-time dimensional analysis and unit/quantity manipulation.

Linear algebra, dimensional analysis, and C++ intersect in real-world domains: robotics, autonomous vehicles, aerospace, structural engineering, medical engineering, scientific simulation, computational fluid dynamics, etc. The [Buckingham π theorem](https://en.wikipedia.org/wiki/Buckingham_%CF%80_theorem) states the equation can be rewritten in terms of a set of dimensionless parameters. The nondimensionalization produces natural scales. The Buckingham π theorem is not always applied.

The C++ implementers expect linear algebra and dimensional analysis to be compatible and also expect type safeties to propagate. C++ implementers seeks linear algebra and unit support [[1](https://stackoverflow.com/questions/8120126/combining-a-linear-algebra-library-with-boostunits)]:

> Is there any matrix library that supports units? *-2011 Stackoverflow user*

## Typed Linear Algebra

A variation of typed linear algebra consists of the injection of strong types into existing linear algebra libraries [[2](https://eigen.tuxfamily.org/dox/TopicCustomizing_CustomScalar.html)][[3](https://github.com/blitzpp/blitz)][[4](https://stackoverflow.com/questions/62040237/eigen-with-custom-scalar-types-matrix-multiplication-with-custom-type-fails-wit)]. Non-uniformly typed vectors or matrices are challenging: each element may have a different type and the library not equiped with the type machinery. Then the results of operations require specialization, customization points, so that multiplying two elements in meters units, computes to square meters.

Non-uniformly typed matrices may present with up to `m x n` types. Not every single type needs to be tracked. Hall explains Hart's conjecture that all `m x n` usable matrices can encode all their element dimensions in just two m-, n-vectors [5]. Daniel Withopf [[6](https://www.youtube.com/watch?v=4LmMwhM8ODI)] and Chip Hogg [[7](https://www.youtube.com/watch?v=5dhFtSu3wCo)] uses the conjecture to present a typed matrix composing a traditional linear algebra matrix with strong types.

## Lessons Learned

## Future

* Ensure the library benefits from expression templates.
* Should higher dimensions, beyond the two for linear algebra, be supported? 

## References

[[1](https://stackoverflow.com/questions/8120126/combining-a-linear-algebra-library-with-boostunits)]
[[2](https://eigen.tuxfamily.org/dox/TopicCustomizing_CustomScalar.html)]
[[3](https://github.com/blitzpp/blitz)]
[[4](https://stackoverflow.com/questions/62040237/eigen-with-custom-scalar-types-matrix-multiplication-with-custom-type-fails-wit)]
[5] Hall, Blair. (2002). Software support for physical quantities. 10.5281/zenodo.5950895. 
[[6] Daniel Withopf - Physical Units for Matrices. How hard can it be? - Meeting C++ 2021](https://www.youtube.com/watch?v=4LmMwhM8ODI)
[[7] Units Libraries and Autonomous Vehicles: Lessons from the Trenches - Chip Hogg - CppCon 2021](https://www.youtube.com/watch?v=5dhFtSu3wCo)

FUTURE

Read https://www.met.reading.ac.uk/clouds/cpp_arrays/cpp_arrays.pdf
Explore https://www.boost.org/doc/libs/1_84_0/libs/qvm/doc/html/index.html
Read Veldhuizen, T., 1995: Expression templates. C++ Report, 7, 2631-.
Read Iglberger, K., G. Hager, J. Treibig and U. R ¨ude, 2012: Expression templates revisited: A performance analysis of current methodologies. SIAM J. Sci. Comp., 34, C42–C69.