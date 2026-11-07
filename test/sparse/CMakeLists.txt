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

build("concepts" BACKENDS "sparse")
build("constexpr_arithmetic" BACKENDS "sparse")
build("freestanding" BACKENDS "sparse_freestanding")
build("operation_count" BACKENDS "sparse")
build("pattern_algebra" BACKENDS "sparse")
build("pattern_large" BACKENDS "sparse")
build("pattern_make" BACKENDS "sparse")
build("pattern_named" BACKENDS "sparse")
build("result_patterns" BACKENDS "sparse")
build("scale_offset_nholthaus" BACKENDS "nholthaus_sparse")
build("storage_layout" BACKENDS "sparse")
fail("addition_au_fail" BACKENDS "au_sparse")
fail("addition_chrono_fail" BACKENDS "chrono_sparse")
fail("addition_mp_units_fail" BACKENDS "mp_units_sparse")
fail("addition_nholthaus_fail" BACKENDS "nholthaus_sparse")
fail("addition_shape_fail" BACKENDS "sparse")
fail("assign_narrowing_fail" BACKENDS "sparse")
fail("at_assign_chrono_fail" BACKENDS "chrono_sparse")
fail("at_assign_fail" BACKENDS "sparse")
fail("at_assign_mp_units_fail" BACKENDS "mp_units_sparse")
fail("at_runtime_assign_chrono_fail" BACKENDS "chrono_sparse")
fail("at_structural_zero_au_fail" BACKENDS "au_sparse")
fail("at_structural_zero_chrono_fail" BACKENDS "chrono_sparse")
fail("at_structural_zero_fail" BACKENDS "sparse")
fail("at_structural_zero_mp_units_fail" BACKENDS "mp_units_sparse")
fail("at_structural_zero_nholthaus_fail" BACKENDS "nholthaus_sparse")
fail("at_structural_zero_singleton_fail" BACKENDS "sparse")
fail("at_structural_zero_vector_fail" BACKENDS "sparse")
fail("list_fail" BACKENDS "sparse")
fail("make_matrix_count_fail" BACKENDS "sparse")
fail("make_matrix_structural_zero_fail" BACKENDS "sparse")
fail("make_matrix_type_fail" BACKENDS "au_sparse")
fail("make_matrix_value_fail" BACKENDS "sparse")
fail("convert_narrowing_fail" BACKENDS "sparse")
fail("division_fail" BACKENDS "sparse")
fail("multiplication_shape_fail" BACKENDS "sparse")
fail("norm_integer_fail" BACKENDS "sparse")
fail("pattern_duplicate_fail" BACKENDS "sparse")
fail("pattern_out_of_range_fail" BACKENDS "sparse")
fail("pattern_shape_fail" BACKENDS "sparse")
fail("representation_fail" BACKENDS "sparse")
fail("scale_offset_nholthaus_fail" BACKENDS "nholthaus_sparse")
pass("addition" BACKENDS "sparse")
pass("allocation" BACKENDS "sparse")
pass("at" BACKENDS "sparse")
pass("constructor" BACKENDS "sparse")
pass("differential" BACKENDS "sparse")
pass("dot" BACKENDS "sparse")
pass("equal_to" BACKENDS "sparse")
pass("format" BACKENDS "sparse")
pass("format_chrono" BACKENDS "chrono_sparse")
pass("format_mp_units" BACKENDS "mp_units_sparse")
pass("format_nholthaus" BACKENDS "nholthaus_sparse")
pass("large" BACKENDS "sparse")
pass("multiplication" BACKENDS "sparse")
pass("non_finite" BACKENDS "sparse")
pass("norm" BACKENDS "sparse")
pass("scale" BACKENDS "sparse")
pass("scale_au" BACKENDS "au_sparse")
pass("scale_chrono" BACKENDS "chrono_sparse")
pass("scale_mp_units" BACKENDS "mp_units_sparse")
pass("scale_nholthaus" BACKENDS "nholthaus_sparse")
pass("singleton" BACKENDS "sparse")
pass("subtraction" BACKENDS "sparse")
pass("transposed" BACKENDS "sparse")
pass("unit_au" BACKENDS "au_sparse")
pass("unit_chrono" BACKENDS "chrono_sparse")
pass("unit_mp_units" BACKENDS "mp_units_sparse")
pass("unit_nholthaus" BACKENDS "nholthaus_sparse")
pass("vector" BACKENDS "sparse")
