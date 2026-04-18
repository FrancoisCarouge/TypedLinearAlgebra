# Benchmarks

Build and run the benchmarks on Windows platform:

```shell
git clone --depth 1 "https://github.com/FrancoisCarouge/TypedLinearAlgebra" "linalg"
Remove-Item -Path build -Force -Recurse
cmake -S "linalg" -B "build" -G "Visual Studio 18 2026" -DBUILD_BENCHMARKING=ON
cmake --build "build" --config "Release" --parallel
ctest --test-dir "build" --build-config "Release" --tests-regex "bench" --parallel 1
```

# Results

Disclaimer: naive benchmark results for illustration purposes only.

MSVC 19.50.35728.0

<img src="plot.png">

| Title | Size | Median Elapsed Time (s) | Median Absolute Error (%) |
| --- | --- | --- | --- |
| Eigen Matrix-Matrix Product |     1x1     | 3.86926308934042e-10 | 0.00362170772453916 |
| Eigen Matrix-Matrix Product |     2x2     | 2.1422574761891e-09 | 0.0047132763352075 |
| Eigen Matrix-Matrix Product |     4x4     | 5.79035495617265e-09 | 0.00263924084699713 |
| Eigen Matrix-Matrix Product |     8x8     | 8.82471891935614e-08 | 0.00197539274178802 |
| Eigen Matrix-Matrix Product |    16x16    | 4.81097046413502e-07 | 0.00607296084493375 |
| Eigen Matrix-Matrix Product |    32x32    | 3.1941348973607e-06 | 0.0031633792869026 |
| Eigen Matrix-Matrix Product |    64x64    | 2.441875e-05 | 0.000866107153313199 |
| Eigen Matrix-Matrix Product |   128x128   | 0.00018515 | 0.00144235103218254 |
| Typed Eigen Matrix-Matrix Product |     1x1     | 5.79379930032255e-10 | 0.000702695077781357 |
| Typed Eigen Matrix-Matrix Product |     2x2     | 2.14140880945798e-09 | 0.00355914946730552 |
| Typed Eigen Matrix-Matrix Product |     4x4     | 5.79967653021709e-09 | 0.00210146423312762 |
| Typed Eigen Matrix-Matrix Product |     8x8     | 8.78154348713761e-08 | 0.00488743624686649 |
| Typed Eigen Matrix-Matrix Product |    16x16    | 4.79429987608426e-07 | 0.00162202324815421 |
| Typed Eigen Matrix-Matrix Product |    32x32    | 3.185e-06 | 0.00201527791737873 |
| std::mdspan Matrix-Matrix Product |     1x1     | 3.88496058575182e-10 | 0.00308631590095693 |
| std::mdspan Matrix-Matrix Product |     2x2     | 2.17730207449334e-09 | 0.0107109395221211 |
| std::mdspan Matrix-Matrix Product |     4x4     | 1.1645345319388e-08 | 0.0012460454495222 |
| std::mdspan Matrix-Matrix Product |     8x8     | 7.47712762456394e-08 | 0.00265945900784587 |
| std::mdspan Matrix-Matrix Product |    16x16    | 6.01623740201568e-07 | 0.00110916614981304 |
| std::mdspan Matrix-Matrix Product |    32x32    | 5.42979797979798e-06 | 0.00140469106799371 |
| std::mdspan Matrix-Matrix Product |    64x64    | 5.4655e-05 | 0.00161881490531467 |
| std::mdspan Matrix-Matrix Product |   128x128   | 0.0007057 | 0.00498433494730835 |
| std::mdspan/tuple Matrix-Matrix Product |     1x1     | 2.32316768458827e-09 | 0.000836487956535696 |
| std::mdspan/tuple Matrix-Matrix Product |     2x2     | 2.21534194128103e-08 | 0.00171143975547006 |
| std::mdspan/tuple Matrix-Matrix Product |     4x4     | 2.02108662613982e-07 | 0.00449747941789824 |
| std::mdspan/tuple Matrix-Matrix Product |     8x8     | 2.95777777777778e-06 | 0.0104089219330855 |
| std::mdspan/tuple Matrix-Matrix Product |    16x16    | 4.90681818181818e-05 | 0.00208765050727594 |
| Typed std::mdspan/tuple Matrix-Matrix Product |     1x1     | 2.32535507064431e-09 | 0.0015576140277047 |
| Typed std::mdspan/tuple Matrix-Matrix Product |     2x2     | 2.26439157394784e-08 | 0.00122223329852446 |
| Typed std::mdspan/tuple Matrix-Matrix Product |     4x4     | 2.00973053892216e-07 | 0.000494009371200135 |
| Typed std::mdspan/tuple Matrix-Matrix Product |     8x8     | 3.25163043478261e-06 | 0.0031593568021656 |
| Typed std::mdspan/tuple Matrix-Matrix Product |    16x16    | 4.91428571428571e-05 | 0.00590143016640947 |
