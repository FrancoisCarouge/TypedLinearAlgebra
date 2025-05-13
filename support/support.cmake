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

# Add a given test.
#
# * NAME The name of the test file without extension. Add a given test.
#
# * NAME The name of the test file without extension.
# * BACKENDS Optional list of backends to use against the test.
function(test TEST_NAME)
  set(multiValueArgs BACKENDS)
  cmake_parse_arguments(PARSE_ARGV 0 TEST "" "${oneValueArgs}"
                        "${multiValueArgs}")

  foreach(BACKEND IN ITEMS ${TEST_BACKENDS})
    if((CMAKE_CXX_COMPILER_ID STREQUAL "MSVC") AND (BACKEND STREQUAL "quantity"
                                                   ))
      message(STATUS "${TEST_NAME} not yet compatible with MSVC/mp-units.")
      continue()
    endif()

    add_executable(typed_linear_algebra_test_${BACKEND}_${TEST_NAME}_driver
                   "${TEST_NAME}.cpp")
    target_link_libraries(
      typed_linear_algebra_test_${BACKEND}_${TEST_NAME}_driver
      PRIVATE typed_linear_algebra_main typed_linear_algebra_${BACKEND}
              typed_linear_algebra_options)
    separate_arguments(TEST_COMMAND UNIX_COMMAND $ENV{COMMAND})
    add_test(
      NAME typed_linear_algebra_test_${BACKEND}_${TEST_NAME}
      COMMAND
        ${TEST_COMMAND}
        $<TARGET_FILE:typed_linear_algebra_test_${BACKEND}_${TEST_NAME}_driver>)
  endforeach()
endfunction(test)
