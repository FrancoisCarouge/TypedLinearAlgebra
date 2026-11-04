#[[ Typed Linear Algebra
Version 0.4.0
https://github.com/FrancoisCarouge/TypedLinearAlgebra

SPDX-License-Identifier: Unlicense

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

For more information, please refer to <https://unlicense.org> ]]

fail("add_chrono_std_fail" BACKENDS "chrono_std")
fail("add_mp_units_std_fail" BACKENDS "mp_units_std")
fail("add_nholthaus_std_fail" BACKENDS "nholthaus_std")
fail("addition_2x2_au_fail" BACKENDS "au_armadillo" "au_eigen" "au_sparse")
fail(
  "addition_2x2_chrono_fail"
  BACKENDS "chrono_armadillo" "chrono_eigen" "chrono_sparse")
fail(
  "addition_2x2_mp_units_fail"
  BACKENDS "mp_units_armadillo" "mp_units_eigen" "mp_units_sparse")
fail(
  "addition_2x2_nholthaus_fail"
  BACKENDS "nholthaus_armadillo" "nholthaus_eigen" "nholthaus_sparse")
fail("addition_au_fail" BACKENDS "au_armadillo" "au_eigen" "au_sparse")
fail(
  "addition_chrono_fail"
  BACKENDS "chrono_armadillo" "chrono_eigen" "chrono_sparse")
fail(
  "addition_mp_units_fail"
  BACKENDS "mp_units_armadillo" "mp_units_eigen" "mp_units_sparse")
fail(
  "addition_nholthaus_fail"
  BACKENDS "nholthaus_armadillo" "nholthaus_eigen" "nholthaus_sparse")
fail("division_au_fail" BACKENDS "au_eigen" "au_sparse")
fail("division_chrono_fail" BACKENDS "chrono_eigen" "chrono_sparse")
fail("division_mp_units_fail" BACKENDS "mp_units_eigen" "mp_units_sparse")
fail("division_nholthaus_fail" BACKENDS "nholthaus_eigen" "nholthaus_sparse")
fail("dot_au_fail" BACKENDS "au_armadillo" "au_eigen" "au_sparse")
fail("dot_au_std_fail" BACKENDS "au_std")
fail(
  "dot_chrono_fail"
  BACKENDS "chrono_armadillo" "chrono_eigen" "chrono_sparse")
fail("dot_chrono_std_fail" BACKENDS "chrono_std")
fail(
  "dot_mp_units_fail"
  BACKENDS "mp_units_armadillo" "mp_units_eigen" "mp_units_sparse")
fail("dot_mp_units_std_fail" BACKENDS "mp_units_std")
fail(
  "dot_nholthaus_fail"
  BACKENDS "nholthaus_armadillo" "nholthaus_eigen" "nholthaus_sparse")
fail("dot_nholthaus_std_fail" BACKENDS "nholthaus_std")
fail("matrix_product_au_std_fail" BACKENDS "au_std")
fail("matrix_product_chrono_std_fail" BACKENDS "chrono_std")
fail("matrix_product_mp_units_std_fail" BACKENDS "mp_units_std")
fail("matrix_product_nholthaus_std_fail" BACKENDS "nholthaus_std")
fail("matrix_rank_1_update_au_std_fail" BACKENDS "au_std")
fail("matrix_rank_1_update_chrono_std_fail" BACKENDS "chrono_std")
fail("matrix_rank_1_update_mp_units_std_fail" BACKENDS "mp_units_std")
fail("matrix_rank_1_update_nholthaus_std_fail" BACKENDS "nholthaus_std")
fail("matrix_vector_product_au_std_fail" BACKENDS "au_std")
fail("matrix_vector_product_chrono_std_fail" BACKENDS "chrono_std")
fail("matrix_vector_product_mp_units_std_fail" BACKENDS "mp_units_std")
fail("matrix_vector_product_nholthaus_std_fail" BACKENDS "nholthaus_std")
fail("subtraction_2x2_au_fail" BACKENDS "au_armadillo" "au_eigen" "au_sparse")
fail(
  "subtraction_2x2_chrono_fail"
  BACKENDS "chrono_armadillo" "chrono_eigen" "chrono_sparse")
fail(
  "subtraction_2x2_mp_units_fail"
  BACKENDS "mp_units_armadillo" "mp_units_eigen" "mp_units_sparse")
fail(
  "subtraction_2x2_nholthaus_fail"
  BACKENDS "nholthaus_armadillo" "nholthaus_eigen" "nholthaus_sparse")
fail("subtraction_au_fail" BACKENDS "au_armadillo" "au_eigen" "au_sparse")
fail(
  "subtraction_chrono_fail"
  BACKENDS "chrono_armadillo" "chrono_eigen" "chrono_sparse")
fail(
  "subtraction_mp_units_fail"
  BACKENDS "mp_units_armadillo" "mp_units_eigen" "mp_units_sparse")
fail(
  "subtraction_nholthaus_fail"
  BACKENDS "nholthaus_armadillo" "nholthaus_eigen" "nholthaus_sparse")

build("conversion_au_std" BACKENDS "au_std")
build("conversion_chrono_std" BACKENDS "chrono_std")
build("conversion_mp_units_std" BACKENDS "mp_units_std")
build("conversion_nholthaus_std" BACKENDS "nholthaus_std")
build("product_au" BACKENDS "au_armadillo" "au_eigen" "au_sparse")
build(
  "product_chrono"
  BACKENDS "chrono_armadillo" "chrono_eigen" "chrono_sparse")
build(
  "product_mp_units"
  BACKENDS "mp_units_armadillo" "mp_units_eigen" "mp_units_sparse")
build(
  "product_nholthaus"
  BACKENDS "nholthaus_armadillo" "nholthaus_eigen" "nholthaus_sparse")

pass("cancelling_units_au" BACKENDS "au_armadillo" "au_eigen" "au_sparse")
pass(
  "cancelling_units_mp_units"
  BACKENDS "mp_units_armadillo" "mp_units_eigen" "mp_units_sparse")
pass("control_au" BACKENDS "au_armadillo" "au_eigen" "au_sparse")
pass(
  "control_chrono"
  BACKENDS "chrono_armadillo" "chrono_eigen" "chrono_sparse")
pass(
  "control_mp_units"
  BACKENDS "mp_units_armadillo" "mp_units_eigen" "mp_units_sparse")
pass(
  "control_nholthaus"
  BACKENDS "nholthaus_armadillo" "nholthaus_eigen" "nholthaus_sparse")
pass("conversion_au" BACKENDS "au_armadillo" "au_eigen" "au_sparse")
pass(
  "conversion_chrono"
  BACKENDS "chrono_armadillo" "chrono_eigen" "chrono_sparse")
pass(
  "conversion_mp_units"
  BACKENDS "mp_units_armadillo" "mp_units_eigen" "mp_units_sparse")
pass(
  "conversion_nholthaus"
  BACKENDS "nholthaus_armadillo" "nholthaus_eigen" "nholthaus_sparse")
pass("conversion_same_scale_mp_units_std" BACKENDS "mp_units_std")
pass("division_au" BACKENDS "au_eigen")
pass("division_chrono" BACKENDS "chrono_eigen")
pass("division_mp_units" BACKENDS "mp_units_eigen")
pass("division_nholthaus" BACKENDS "nholthaus_eigen")
pass("explicit_conversion_au" BACKENDS "au_armadillo" "au_eigen" "au_sparse")
pass(
  "explicit_conversion_chrono"
  BACKENDS "chrono_armadillo" "chrono_eigen" "chrono_sparse")
pass(
  "explicit_conversion_mp_units"
  BACKENDS "mp_units_armadillo" "mp_units_eigen" "mp_units_sparse")
pass(
  "explicit_conversion_nholthaus"
  BACKENDS "nholthaus_armadillo" "nholthaus_eigen" "nholthaus_sparse")
