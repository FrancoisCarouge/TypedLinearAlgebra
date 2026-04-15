# Benchmarks

Build and run the benchmarks on Windows platform:

```shell
git clone --depth 1 "https://github.com/FrancoisCarouge/TypedLinearAlgebra" "linalg"
Remove-Item -Path build -Force -Recurse
cmake -S "linalg" -B "build" -G "Visual Studio 18 2026" -DBUILD_BENCHMARKING=ON
cmake --build "build" --config "Release" --parallel
ctest --test-dir "build" --build-config "Release" --tests-regex "benchmark" --parallel 1
type build/benchmark/results.txt
```

# Results

Disclaimer: naive benchmark results for illustration purposes only.

MSVC 19.50.35728.0

| Title | Size | Median Elapsed Time (s) | Median Absolute Error (%) |
| --- | --- | --- | --- |
| Eigen Matrix-Matrix Product |     1x1     | 5.79392722831808e-10 | 0.00144623780651739 |
| Eigen Matrix-Matrix Product |     2x2     | 1.95589001618936e-09 | 0.00269177459439564 |
| Eigen Matrix-Matrix Product |     4x4     | 5.81199851000211e-09 | 0.00518763823073691 |
| Eigen Matrix-Matrix Product |     8x8     | 8.80302913221544e-08 | 0.00416339388806324 |
| Eigen Matrix-Matrix Product |    16x16    | 4.75733774286895e-07 | 0.00480277383354175 |
| Eigen Matrix-Matrix Product |    32x32    | 3.20241286863271e-06 | 0.000558270738808668 |
| Eigen Matrix-Matrix Product |    64x64    | 2.44744680851064e-05 | 0.00139287890659009 |
| Eigen Matrix-Matrix Product |   128x128   | 0.000187566666666667 | 0.00860369241799612 |
| Typed Eigen Matrix-Matrix Product |     1x1     | 3.8802152457711e-10 | 0.00600626775460243 |
| Typed Eigen Matrix-Matrix Product |     2x2     | 2.3254118233246e-09 | 0.00473804888396924 |
| Typed Eigen Matrix-Matrix Product |     4x4     | 5.79325954965105e-09 | 0.00240710034370305 |
| Typed Eigen Matrix-Matrix Product |     8x8     | 8.81667375584858e-08 | 0.00815057380692842 |
| Typed Eigen Matrix-Matrix Product |    16x16    | 4.77938432835821e-07 | 0.00472518473777041 |
| Typed Eigen Matrix-Matrix Product |    32x32    | 3.18621700879765e-06 | 0.00283102217796413 |
| std::mdspan Matrix-Matrix Product |     1x1     | 4.1379977341181e-10 | 0.0692817846067371 |
| std::mdspan Matrix-Matrix Product |     2x2     | 2.16724909570229e-09 | 0.0122684458015465 |
| std::mdspan Matrix-Matrix Product |     4x4     | 1.16626640363668e-08 | 0.00660626950214691 |
| std::mdspan Matrix-Matrix Product |     8x8     | 7.46279564177518e-08 | 0.00247378901236277 |
| std::mdspan Matrix-Matrix Product |    16x16    | 6.0378408458542e-07 | 0.00360272025138946 |
| std::mdspan Matrix-Matrix Product |    32x32    | 5.41064814814815e-06 | 0.00574930294847302 |
| std::mdspan Matrix-Matrix Product |    64x64    | 5.47904761904762e-05 | 0.00229536614792269 |
| std::mdspan Matrix-Matrix Product |   128x128   | 0.000741 | 0.00134770889487874 |
