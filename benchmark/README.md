# Benchmarks

Build and run the benchmarks on all platforms:

```shell
git clone --depth 1 "https://github.com/FrancoisCarouge/TypedLinearAlgebra" "linalg"
Remove-Item -Path build -Force -Recurse
cmake -S "linalg" -B "build" -G "Visual Studio 18 2026"
cmake --build "build" --config "Release" --parallel 20
ctest --test-dir "build" --build-config "Release" --tests-regex "benchmark" --verbose --parallel 1
```

# Results

Disclaimer: naive benchmark results for illustration purposes only.

| Title | Size | Median Elapsed Time (s) | Median Absolute Error (%) |
| --- | --- | --- | --- |
| Eigen Matrix-Matrix Product |     1x1     | 3.88117474999264e-10 | 0.00488067571184614 |
| Eigen Matrix-Matrix Product |     2x2     | 2.1369873786026e-09 | 0.00189067028605782 |
| Eigen Matrix-Matrix Product |     4x4     | 5.91550106487078e-09 | 0.000924268548847455 |
| Eigen Matrix-Matrix Product |     8x8     | 8.99013245621833e-08 | 0.00694897213796576 |
| Eigen Matrix-Matrix Product |    16x16    | 4.85027185278126e-07 | 0.00175669845586784 |
| Eigen Matrix-Matrix Product |    32x32    | 3.24067796610169e-06 | 0.00395230900414392 |
| Eigen Matrix-Matrix Product |    64x64    | 2.46234042553191e-05 | 0.00135515970684615 |
| Eigen Matrix-Matrix Product |   128x128   | 0.00018736 | 0.00246120920278212 |
