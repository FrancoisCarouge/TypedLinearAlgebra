# Typed Linear Algebra

This library provides a C++ statically type-safe linear algebra facade to the Eigen linear algebra library.

# Installation

Example of installation commands in Shell:

```shell
git clone --depth 1 "https://github.com/FrancoisCarouge/TypedLinearAlgebra"
cmake -S "TypedLinearAlgebra" -B "build"
cmake --build "build" --parallel
sudo cmake --install "build"
```

Another variation for your CMake infrastructure via fetch content:

```cmake
include(FetchContent)

FetchContent_Declare(
  fcarouge-typed-linear-algebra
  GIT_REPOSITORY "https://github.com/FrancoisCarouge/TypedLinearAlgebra"
  FIND_PACKAGE_ARGS NAMES fcarouge-typed-linear-algebra)
FetchContent_MakeAvailable(fcarouge-typed-linear-algebra)

target_link_libraries(your_target PRIVATE fcarouge-typed-linear-algebra::linalg)
```

[For more, see installation instructions](https://github.com/FrancoisCarouge/TypedLinearAlgebra/tree/master/INSTALL.md).

# Resources

## Third Party Acknowledgement

The library is designed, developed, and tested with the help of third-party tools and services acknowledged and thanked here:

- [actions-gh-pages](https://github.com/peaceiris/actions-gh-pages) to upload the documentation to GitHub pages.
- [Clang](https://clang.llvm.org) for compilation and code sanitizers.
- [CMake](https://cmake.org) for build automation.
- [cmakelang](https://pypi.org/project/cmakelang) for pretty CMake list files.
- [cppcheck](https://cppcheck.sourceforge.io) for static analysis.
- [Doxygen](https://doxygen.nl) for documentation generation.
- [Doxygen Awesome](https://github.com/jothepro/doxygen-awesome-css) for pretty documentation.
- [Eigen](https://eigen.tuxfamily.org/) for linear algebra.
- [GCC](https://gcc.gnu.org) for compilation and code sanitizers.
- [Google Benchmark](https://github.com/google/benchmark) to implement the benchmarks.
- [lcov](http://ltp.sourceforge.net/coverage/lcov.php) to process coverage information.
- [mp-units](https://github.com/mpusz/mp-units) the quantities and units library for C++.
- [MSVC](https://docs.microsoft.com/en-US/cpp/windows/latest-supported-vc-redist) for compilation and code sanitizers.
- [Valgrind](https://valgrind.org) to check for correct memory management.

## Sponsors

Become a sponsor today! Support this project with coffee and infrastructure!

[![Sponsor](https://img.shields.io/badge/Support-Sponsor-brightgreen)](http://paypal.me/francoiscarouge)

### Corporations & Institutions

*Your group logo and link here!*

### Individuals

*Your name and link here!*

Thanks everyone!

# Continuous Integration & Deployment Actions

[![Code Repository](https://img.shields.io/badge/Repository-GitHub%20%F0%9F%94%97-brightgreen)](https://github.com/FrancoisCarouge/TypedLinearAlgebra)
<br>
<br>
[![Pipeline](https://github.com/FrancoisCarouge/TypedLinearAlgebra/actions/workflows/pipeline.yml/badge.svg)](https://github.com/FrancoisCarouge/TypedLinearAlgebra/actions/workflows/pipeline.yml)
<br>
<br>
[![Sanitizer](https://github.com/FrancoisCarouge/TypedLinearAlgebra/actions/workflows/sanitizer.yml/badge.svg)](https://github.com/FrancoisCarouge/TypedLinearAlgebra/actions/workflows/sanitizer.yml)
<br>
[![Format](https://github.com/FrancoisCarouge/TypedLinearAlgebra/actions/workflows/format.yml/badge.svg)](https://github.com/FrancoisCarouge/TypedLinearAlgebra/actions/workflows/format.yml)
<br>
[![ClangTidy](https://github.com/FrancoisCarouge/TypedLinearAlgebra/actions/workflows/clang_tidy.yml/badge.svg)](https://github.com/FrancoisCarouge/TypedLinearAlgebra/actions/workflows/clang_tidy.yml)
<br>
[![CppCheck](https://github.com/FrancoisCarouge/TypedLinearAlgebra/actions/workflows/cppcheck.yml/badge.svg)](https://github.com/FrancoisCarouge/TypedLinearAlgebra/actions/workflows/cppcheck.yml)
<br>
[![Doxygen](https://github.com/FrancoisCarouge/TypedLinearAlgebra/actions/workflows/doxygen.yml/badge.svg)](https://github.com/FrancoisCarouge/TypedLinearAlgebra/actions/workflows/doxygen.yml)
<br>
[![Valgrind](https://github.com/FrancoisCarouge/TypedLinearAlgebra/actions/workflows/memory_valgrind.yml/badge.svg)](https://github.com/FrancoisCarouge/TypedLinearAlgebra/actions/workflows/memory_valgrind.yml)
<br>
<br>
[![Public Domain](https://img.shields.io/badge/License-Public%20Domain%20%F0%9F%94%97-brightgreen)](https://raw.githubusercontent.com/francoiscarouge/TypedLinearAlgebra/master/LICENSE.txt)
<br>
[![License Scan](https://app.fossa.com/api/projects/git%2Bgithub.com%2FFrancoisCarouge%2FTypedLinearAlgebra.svg?type=shield)](https://app.fossa.com/projects/git%2Bgithub.com%2FFrancoisCarouge%2FTypedLinearAlgebra?ref=badge_shield)
<br>
[![OpenSSF Best Practices](https://www.bestpractices.dev/projects/8933/badge)](https://www.bestpractices.dev/projects/8933)
<br>
<br>
[![Deploy Unit Test Code Coverage](https://github.com/FrancoisCarouge/TypedLinearAlgebra/actions/workflows/deploy_coverage.yml/badge.svg)](https://francoiscarouge.github.io/TypedLinearAlgebra/unit_test_coverage.xhtml)
<br>
[![Deploy Doxygen](https://github.com/FrancoisCarouge/TypedLinearAlgebra/actions/workflows/deploy_doxygen.yml/badge.svg)](https://francoiscarouge.github.io/TypedLinearAlgebra/index.xhtml)
<br>
<br>
[![Sponsor](https://img.shields.io/badge/Support-Sponsor%20%F0%9F%94%97-brightgreen)](http://paypal.me/francoiscarouge)

# License

<img align="right" src="http://opensource.org/trademarks/opensource/OSI-Approved-License-100x137.png">

TypedLinearAlgebra is public domain:

This is free and unencumbered software released into the public domain.

Anyone is free to copy, modify, publish, use, compile, sell, or
distribute this software, either in source code form or as a compiled
binary, for any purpose, commercial or non-commercial, and by any
means.

In jurisdictions that recognize copyright laws, the author or authors
of this software dedicate any and all copyright interest in the
software to the public domain. We make this dedication for the benefit
of the public at large and to the detriment of our heirs and
successors. We intend this dedication to be an overt act of
relinquishment in perpetuity of all present and future rights to this
software under copyright law.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR
OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
OTHER DEALINGS IN THE SOFTWARE.

For more information, please refer to <https://unlicense.org>
