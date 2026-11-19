# Installation

Requires CMake ≥ 4.3 and a C++26 compiler: GCC 14+, Clang 20+ (libstdc++ 14+ or libc++), MSVC `/std:c++latest` (Visual Studio 2026), or Apple Clang (Xcode 26).

Download and install the [latest release package](https://github.com/FrancoisCarouge/TypedLinearAlgebra/releases). Alternatively, you may install and use the library in your projects by cloning the repository, configuring, and installing the project:

```shell
git clone --depth 1 "https://github.com/FrancoisCarouge/TypedLinearAlgebra"
cmake -S "TypedLinearAlgebra" -B "build"
cmake --build "build" --parallel
sudo cmake --install "build"
```

Only the `include/` folder is the library: `cmake --install` and the release packages deliver it, with its CMake package configuration, and nothing else. The `support/`, `sample/`, `test/`, and `benchmark/` folders are development examples and tooling, not part of the delivered library.

The standard shared CMake configuration file provides the library target to use in your own target:

```cmake
find_package(fcarouge-typed-linear-algebra)
target_link_libraries(your_target PRIVATE fcarouge-typed-linear-algebra::tlinalg)
```

Alternatively, fetch the library directly from your project's CMake configuration:

```cmake
include(FetchContent)

FetchContent_Declare(
  fcarouge-typed-linear-algebra
  GIT_REPOSITORY "https://github.com/FrancoisCarouge/TypedLinearAlgebra"
  FIND_PACKAGE_ARGS NAMES fcarouge-typed-linear-algebra)
FetchContent_MakeAvailable(fcarouge-typed-linear-algebra)

target_link_libraries(your_target PRIVATE fcarouge-typed-linear-algebra::tlinalg)
```

A single, amalgamated `typed_linear_algebra.h` header is also [published](https://francoiscarouge.github.io/TypedLinearAlgebra/amalgamate/fcarouge/typed_linear_algebra.h), or generated locally with `cmake --build "build" --target "amalgamate"` into `build/amalgamate/fcarouge/typed_linear_algebra.h`.

In your sources, include the library header and add a plug-in for each strong type and linear algebra backend you use. See [the samples](https://github.com/FrancoisCarouge/TypedLinearAlgebra/tree/master/sample) for more.

```cpp
#include "fcarouge/typed_linear_algebra.hpp"
```

# Development Build & Run

## Tests & Samples

Build and run the tests and samples:

```shell
git clone --depth 1 "https://github.com/FrancoisCarouge/TypedLinearAlgebra"
cmake -S "TypedLinearAlgebra" -B "build"
cmake --build "build" --config "Debug" --parallel
ctest --test-dir "build" --build-config "Debug" --output-on-failure --parallel
```

## Benchmarks

See the [Benchmark](https://github.com/FrancoisCarouge/TypedLinearAlgebra/tree/master/benchmark) section.

## Installation Packages

### Linux

```shell
git clone --depth 1 "https://github.com/FrancoisCarouge/TypedLinearAlgebra"
cmake -S "TypedLinearAlgebra" -B "build"
cmake --build "build" --target "package" --parallel
cmake --build "build" --target "package_source" --parallel
```

### Windows

```shell
git clone --depth 1 "https://github.com/FrancoisCarouge/TypedLinearAlgebra"
cmake -S "TypedLinearAlgebra" -B "build"
cmake --build "build" --target "package" --parallel --config "Release"
cmake --build "build" --target "package_source" --parallel --config "Release"
```
