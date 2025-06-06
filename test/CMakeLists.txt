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

test("addition" BACKENDS "eigen" "eigexed")
test("assign" BACKENDS "eigen" "eigexed")
test("at" BACKENDS "eigexed")
test("constructor_1x1_array" BACKENDS "eigen" "eigexed")
test("constructor_1x1" BACKENDS "eigen" "eigexed")
test("constructor_1xn_array" BACKENDS "eigen" "eigexed")
test("constructor_1xn" BACKENDS "eigen" "eigexed")
test("constructor_initializer_lists" BACKENDS "eigen" "eigexed")
test("constructor_nx1_array" BACKENDS "eigen" "eigexed")
test("constructor_nx1" BACKENDS "eigen" "eigexed")
test("copy" BACKENDS "eigen" "eigexed")
test("format_1x1" BACKENDS "eigen" "eigexed")
test("format_1xn" BACKENDS "eigen" "eigexed")
test("format_mx1" BACKENDS "eigen" "eigexed")
test("format_mxn" BACKENDS "eigen" "eigexed")
test("identity" BACKENDS "eigen" "eigexed")
test("multiplication_arithmetic" BACKENDS "eigen" "eigexed")
test("multiplication_rxc" BACKENDS "eigen" "eigexed")
test("multiplication_sxc" BACKENDS "eigen" "eigexed")
test("operator_bracket" BACKENDS "eigen" "eigexed")
test("operator_equality" BACKENDS "eigen" "eigexed")
test("zero" BACKENDS "eigen" "eigexed")
