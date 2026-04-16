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
| Eigen Matrix-Matrix Product |     1x1     | 5.79280640443907e-10 | 0.00149738964412317 |
| Eigen Matrix-Matrix Product |     2x2     | 2.35519984347582e-09 | 0.0150145464563847 |
| Eigen Matrix-Matrix Product |     4x4     | 5.80984338866267e-09 | 0.00604399015185017 |
| Eigen Matrix-Matrix Product |     8x8     | 8.74154030669065e-08 | 0.0046562539557634 |
| Eigen Matrix-Matrix Product |    16x16    | 4.79342934293429e-07 | 0.00755936512081264 |
| Eigen Matrix-Matrix Product |    32x32    | 3.20851063829787e-06 | 0.00475280948785159 |
| Eigen Matrix-Matrix Product |    64x64    | 2.47531914893617e-05 | 0.003131611321132 |
| Eigen Matrix-Matrix Product |   128x128   | 0.000185733333333333 | 0.00336266723903551 |
| Typed Eigen Matrix-Matrix Product |     1x1     | 5.81548538614975e-10 | 0.00376989768395853 |
| Typed Eigen Matrix-Matrix Product |     2x2     | 2.13931264173522e-09 | 0.00450181246383349 |
| Typed Eigen Matrix-Matrix Product |     4x4     | 5.80778608488658e-09 | 0.00407707131372669 |
| Typed Eigen Matrix-Matrix Product |     8x8     | 8.91929099982791e-08 | 0.010134591604394 |
| Typed Eigen Matrix-Matrix Product |    16x16    | 4.80356327089995e-07 | 0.00283070127249361 |
| Typed Eigen Matrix-Matrix Product |    32x32    | 3.18421052631579e-06 | 0.00250651024946968 |
| std::mdspan Matrix-Matrix Product |     1x1     | 3.89847953144399e-10 | 0.0109093295398069 |
| std::mdspan Matrix-Matrix Product |     2x2     | 2.17775642474226e-09 | 0.00893159004026262 |
| std::mdspan Matrix-Matrix Product |     4x4     | 1.16537275337658e-08 | 0.0093594516510802 |
| std::mdspan Matrix-Matrix Product |     8x8     | 7.45213751458176e-08 | 0.00164401529064526 |
| std::mdspan Matrix-Matrix Product |    16x16    | 6.03575547866205e-07 | 0.00323663159095816 |
| std::mdspan Matrix-Matrix Product |    32x32    | 5.33222748815166e-06 | 0.00638206037613696 |
| std::mdspan Matrix-Matrix Product |    64x64    | 5.50105263157895e-05 | 0.00384172109104879 |
| std::mdspan Matrix-Matrix Product |   128x128   | 0.0007491 | 0.00252996005326231 |
| std::mdspan/tuple Matrix-Matrix Product |     1x1     | 2.32929721125891e-09 | 0.00272458857567422 |
| std::mdspan/tuple Matrix-Matrix Product |     2x2     | 2.20989384370261e-08 | 0.001428002001545 |
| std::mdspan/tuple Matrix-Matrix Product |     4x4     | 2.04769874476987e-07 | 0.0191008736903228 |
| std::mdspan/tuple Matrix-Matrix Product |     8x8     | 2.94362244897959e-06 | 0.00302726988184359 |
| std::mdspan/tuple Matrix-Matrix Product |    16x16    | 4.91695652173913e-05 | 0.00375315404415025 |
