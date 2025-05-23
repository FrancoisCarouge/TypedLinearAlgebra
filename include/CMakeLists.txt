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

add_library(linalg INTERFACE)
target_sources(
  linalg
  INTERFACE FILE_SET
            "typed_linear_algebra_headers"
            TYPE
            "HEADERS"
            FILES
            "fcarouge/typed_linear_algebra_internal/format.hpp"
            "fcarouge/typed_linear_algebra_internal/utility.hpp"
            "fcarouge/typed_linear_algebra.hpp")
install(
  TARGETS linalg
  EXPORT "fcarouge-typed-linear-algebra-target"
  FILE_SET "typed_linear_algebra_headers")

# Conditionally provide the namespace alias target which may be an imported
# target from a package, or an aliased target if built as part of the same
# buildsystem.
if(NOT TARGET fcarouge-typed-linear-algebra::linalg)
  add_library(fcarouge-typed-linear-algebra::linalg ALIAS linalg)
endif()
