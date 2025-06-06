# Installation

Download and install the [latest release package](https://github.com/FrancoisCarouge/TypedLinearAlgebra/releases). Alternatively, you may install and use the library in your projects by cloning the repository, configuring, and installing the project:

```shell
git clone --depth 1 "https://github.com/FrancoisCarouge/TypedLinearAlgebra"
cmake -S "TypedLinearAlgebra" -B "build"
cmake --build "build" --parallel
sudo cmake --install "build"
```

The standard shared CMake configuration file provides the library target to use in your own target:

```cmake
find_package(fcarouge-typed-linear-algebra)
target_link_libraries(your_target PRIVATE fcarouge-typed-linear-algebra::linalg)
```

# Development Build & Run

## Tests

Build and run the tests:

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