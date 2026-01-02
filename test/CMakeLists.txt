#[[ Typed Linear Algebra
Version 0.1.0
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

if(NOT BUILD_TESTING)
  return()
endif()

test("addition" BACKENDS "eigen" "eigexed" "nested_typed_eigen")
test("assign" BACKENDS "eigen" "eigexed" "nested_typed_eigen")
test("at_1x1_const" BACKENDS "eigexed" "nested_typed_eigen")
test("at_1x1" BACKENDS "eigexed" "nested_typed_eigen")
test("at_1xn_const" BACKENDS "eigexed" "nested_typed_eigen")
test("at_1xn" BACKENDS "eigexed" "nested_typed_eigen")
test("at_mx1_const" BACKENDS "eigexed" "nested_typed_eigen")
test("at_mx1" BACKENDS "eigexed" "nested_typed_eigen")
test("at_mxn_const" BACKENDS "eigexed" "nested_typed_eigen")
test("at_mxn" BACKENDS "eigexed" "nested_typed_eigen")
test("constructor_1x1_array" BACKENDS "eigen" "eigexed" "nested_typed_eigen")
test("constructor_1x1_unit" BACKENDS "quantity")
test("constructor_1x1" BACKENDS "eigen" "eigexed" "nested_typed_eigen")
test("constructor_1xn_array" BACKENDS "eigen" "eigexed" "nested_typed_eigen")
test("constructor_1xn" BACKENDS "eigen" "eigexed" "nested_typed_eigen")
test("constructor_copy" BACKENDS "eigexed" "nested_typed_eigen")
test("constructor_initializer_lists" BACKENDS "eigen" "eigexed"
     "nested_typed_eigen")
test("constructor_nx1_array" BACKENDS "eigen" "eigexed" "nested_typed_eigen")
test("constructor_nx1" BACKENDS "eigen" "eigexed" "nested_typed_eigen")
test("copy" BACKENDS "eigen" "eigexed" "nested_typed_eigen")
test("element_1x1" BACKENDS "quantity")
test("element_mxn" BACKENDS "quantity")
test("format_1x1" BACKENDS "eigen" "eigexed" "nested_typed_eigen")
test("format_1xn" BACKENDS "eigen" "eigexed" "nested_typed_eigen")
test("format_mx1" BACKENDS "eigen" "eigexed" "nested_typed_eigen")
test("format_mxn" BACKENDS "eigen" "eigexed" "nested_typed_eigen")
test("identity" BACKENDS "eigen" "eigexed" "nested_typed_eigen")
test("multiplication_arithmetic" BACKENDS "eigen" "eigexed"
     "nested_typed_eigen")
test("multiplication_rxc" BACKENDS "eigen" "eigexed" "nested_typed_eigen")
test("multiplication_sxc" BACKENDS "eigen" "eigexed" "nested_typed_eigen")
test("nested" BACKENDS "nested_typed_eigen")
test("operator_equality" BACKENDS "eigen" "eigexed" "nested_typed_eigen")
test("operator_call_1x1_const" BACKENDS "eigexed" "nested_typed_eigen")
test("operator_call_1x1" BACKENDS "eigexed" "nested_typed_eigen")
test("operator_call_1xn_const" BACKENDS "eigexed" "nested_typed_eigen")
test("operator_call_1xn" BACKENDS "eigexed" "nested_typed_eigen")
test("operator_call_mx1_const" BACKENDS "eigexed" "nested_typed_eigen")
test("operator_call_mx1" BACKENDS "eigexed" "nested_typed_eigen")
test("operator_call_mxn_const" BACKENDS "eigexed" "nested_typed_eigen")
test("operator_call_mxn" BACKENDS "eigexed" "nested_typed_eigen")
test("operator_subscript_1x1_const" BACKENDS "eigexed" "nested_typed_eigen")
test("operator_subscript_1x1" BACKENDS "eigexed" "nested_typed_eigen")
test("operator_subscript_1xn_const" BACKENDS "eigexed" "nested_typed_eigen")
test("operator_subscript_1xn" BACKENDS "eigexed" "nested_typed_eigen")
test("operator_subscript_mx1_const" BACKENDS "eigexed" "nested_typed_eigen")
test("operator_subscript_mx1" BACKENDS "eigexed" "nested_typed_eigen")
test("operator_subscript_mxn_const" BACKENDS "eigexed" "nested_typed_eigen")
test("operator_subscript_mxn" BACKENDS "eigexed" "nested_typed_eigen")
test("structured_bindings_1x1_const" BACKENDS "eigen" "eigexed")
test("structured_bindings_1x1" BACKENDS "eigen" "eigexed")
test("structured_bindings_1xn_const" BACKENDS "eigen" "eigexed")
test("structured_bindings_1xn" BACKENDS "eigen" "eigexed")
test("structured_bindings_mx1_const" BACKENDS "eigen" "eigexed")
test("structured_bindings_mx1" BACKENDS "eigen" "eigexed")
test("underlying" BACKENDS "eigexed" "nested_typed_eigen")
test("zero" BACKENDS "eigen" "eigexed" "nested_typed_eigen")
