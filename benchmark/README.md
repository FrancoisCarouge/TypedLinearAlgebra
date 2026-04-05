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

MSVC 19.50.35728.0


compilers maximum length for a template parameter list

Also recursive depth!! Try to flatten?



|  N |               ns/op |                op/s |    err% |     total | benchmark
|----|--------------------:|--------------------:|--------:|----------:|:----------
|  1 |                0.33 |    3,044,970,754.34 |    0.1% |      0.01 | `matrix_product_mdspan`
|  1 |                0.71 |    1,415,882,301.04 |    1.2% |      0.01 | `matrix_product_tuple`
|  1 |                0.39 |    2,540,143,262.94 |    0.7% |      0.01 | `eigen_product`
|  2 |                1.57 |      638,224,754.11 |    0.3% |      0.01 | `matrix_product_mdspan`
|  2 |               22.08 |       45,299,654.58 |    0.8% |      0.01 | `matrix_product_tuple`
|  2 |                2.14 |      466,760,494.79 |    0.3% |      0.01 | `eigen_product`
|  4 |               11.87 |       84,217,318.12 |    0.6% |      0.01 | `matrix_product_mdspan`
|  4 |              480.31 |        2,081,974.07 |    0.2% |      0.01 | `matrix_product_tuple`
|  4 |                6.22 |      160,786,788.86 |    0.4% |      0.01 | `eigen_product`
|  8 |               74.79 |       13,370,211.97 |    0.6% |      0.01 | `matrix_product_mdspan`
|  8 |               91.16 |       10,969,258.90 |    0.3% |      0.01 | `eigen_product`
|  8 |           13,743.42 |           72,762.09 |    1.9% |      0.01 | `matrix_product_tuple`
| 16 |              602.66 |        1,659,305.99 |    0.1% |      0.01 | `matrix_product_mdspan`
| 16 |          463,400.00 |            2,157.96 |    0.7% |      0.01 | `matrix_product_tuple`
| 16 |              486.18 |        2,056,851.93 |    0.7% |      0.01 | `eigen_product`

| 32 | error C2999: maximum template instantiation depth of 1000 exceeded
> Hahaha! Can't do recursive obviously...
