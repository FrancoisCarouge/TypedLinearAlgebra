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

if(NOT BUILD_TESTING)
  return()
endif()

# These drivers deliberately do not link `tlinalg` or the backend glue
# libraries: they compile solely against the generated, amalgamated headers
# under `${CMAKE_BINARY_DIR}/amalgamate/fcarouge`, to prove those single
# files remain functional, self-sufficient replacements for the multi-file
# layout.

add_executable(typed_linear_algebra_amalgamate_typed_linear_algebra_driver
               "typed_linear_algebra.cpp")
target_include_directories(
  typed_linear_algebra_amalgamate_typed_linear_algebra_driver
  PRIVATE "${CMAKE_BINARY_DIR}/amalgamate/fcarouge")
target_link_libraries(
  typed_linear_algebra_amalgamate_typed_linear_algebra_driver
  PRIVATE Eigen3::Eigen typed_linear_algebra_main typed_linear_algebra_options)
add_dependencies(typed_linear_algebra_amalgamate_typed_linear_algebra_driver
                 amalgamate_typed_linear_algebra)
add_test(NAME typed_linear_algebra_amalgamate_typed_linear_algebra
         COMMAND typed_linear_algebra_amalgamate_typed_linear_algebra_driver)

add_executable(
  typed_linear_algebra_amalgamate_typed_linear_algebra_mp_units_eigen_driver
  "typed_linear_algebra_mp_units_eigen.cpp")
target_include_directories(
  typed_linear_algebra_amalgamate_typed_linear_algebra_mp_units_eigen_driver
  PRIVATE "${CMAKE_BINARY_DIR}/amalgamate/fcarouge")
target_link_libraries(
  typed_linear_algebra_amalgamate_typed_linear_algebra_mp_units_eigen_driver
  PRIVATE Eigen3::Eigen mp-units::mp-units mp-units::integrations
          typed_linear_algebra_main typed_linear_algebra_options)
add_dependencies(
  typed_linear_algebra_amalgamate_typed_linear_algebra_mp_units_eigen_driver
  amalgamate_typed_linear_algebra_mp_units_eigen)
add_test(
  NAME typed_linear_algebra_amalgamate_typed_linear_algebra_mp_units_eigen
  COMMAND
    typed_linear_algebra_amalgamate_typed_linear_algebra_mp_units_eigen_driver
)
