# Benchmarks

Build and run the benchmarks on Windows platform:

```shell
git clone --depth 1 "https://github.com/FrancoisCarouge/TypedLinearAlgebra" "linalg"
Remove-Item -Path build -Force -Recurse
cmake -S "linalg" -B "build" -G "Visual Studio 18 2026"
cmake --build "build" --config "Release" --parallel
ctest --test-dir "build" --build-config "Release" --tests-regex "benchmark" --parallel 1
type build/benchmark/results.txt
```

# Results

Disclaimer: naive benchmark results for illustration purposes only.

| Title | Size | Median Elapsed Time (s) | Median Absolute Error (%) |
| --- | --- | --- | --- |
| Eigen Matrix-Matrix Product |     1x1     | 3.89693150707904e-10 | 0.00649927099426216 |
| Eigen Matrix-Matrix Product |     2x2     | 1.9543461792134e-09 | 0.00396471231832233 |
| Eigen Matrix-Matrix Product |     4x4     | 5.80756561197589e-09 | 0.00451574571252226 |
| Eigen Matrix-Matrix Product |     8x8     | 8.76893203883495e-08 | 0.00250655516972582 |
| Eigen Matrix-Matrix Product |    16x16    | 4.78992983904251e-07 | 0.00591234324719731 |
| Eigen Matrix-Matrix Product |    32x32    | 3.20365853658537e-06 | 0.0023648375755153 |
| Eigen Matrix-Matrix Product |    64x64    | 2.47113636363636e-05 | 0.00246791254965643 |
| Eigen Matrix-Matrix Product |   128x128   | 0.0001862 | 0.00359324469996406 |
| Typed Eigen Matrix-Matrix Product |     1x1     | 5.81106608661144e-10 | 0.00165715285785163 |
| Typed Eigen Matrix-Matrix Product |     2x2     | 1.94718378048236e-09 | 0.00337954547229306 |
| Typed Eigen Matrix-Matrix Product |     4x4     | 5.80331713854917e-09 | 0.00312892572857658 |
| Typed Eigen Matrix-Matrix Product |     8x8     | 8.76508107277338e-08 | 0.00531984945419938 |
| Typed Eigen Matrix-Matrix Product |    16x16    | 4.74411231884058e-07 | 0.00513657374075729 |
| Typed Eigen Matrix-Matrix Product |    32x32    | 3.23988764044944e-06 | 0.00237182653569531 |
