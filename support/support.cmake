#[[ Typed Linear Algebra
Version 0.2.0
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
# * NAME The name of the test file without extension.
#
# * NAME The name of the test file without extension.
# * BACKENDS Optional list of backends to use against the test.
function(pass TEST_NAME)
  set(multiValueArgs BACKENDS)
  cmake_parse_arguments(PARSE_ARGV 0 TEST "" "${oneValueArgs}"
                        "${multiValueArgs}")

  get_filename_component(CALLER "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

  foreach(BACKEND IN ITEMS ${TEST_BACKENDS})
    add_executable(
      typed_linear_algebra_test_${BACKEND}_${CALLER}_${TEST_NAME}_driver
      "${TEST_NAME}.cpp")
    target_link_libraries(
      typed_linear_algebra_test_${BACKEND}_${CALLER}_${TEST_NAME}_driver
      PRIVATE typed_linear_algebra_options typed_linear_algebra_main
              typed_linear_algebra_${BACKEND})
    separate_arguments(TEST_COMMAND UNIX_COMMAND $ENV{COMMAND})
    add_test(
      NAME typed_linear_algebra_test_${BACKEND}_${CALLER}_${TEST_NAME}
      COMMAND
        ${TEST_COMMAND}
        $<TARGET_FILE:typed_linear_algebra_test_${BACKEND}_${CALLER}_${TEST_NAME}_driver>
    )
  endforeach()
endfunction(pass)

# Add a given compile-fail test.
#
# * NAME The name of the test file without extension.
#
# * NAME The name of the test file without extension.
# * BACKENDS Optional list of backends to use against the test.
function(fail TEST_NAME)
  set(multiValueArgs BACKENDS)
  cmake_parse_arguments(PARSE_ARGV 0 TEST "" "${oneValueArgs}"
                        "${multiValueArgs}")

  get_filename_component(CALLER "${CMAKE_CURRENT_SOURCE_DIR}" NAME)

  foreach(BACKEND IN ITEMS ${TEST_BACKENDS})
    add_executable(
      typed_linear_algebra_test_${BACKEND}_${CALLER}_${TEST_NAME}_driver
      "${TEST_NAME}.cpp")
    target_link_libraries(
      typed_linear_algebra_test_${BACKEND}_${CALLER}_${TEST_NAME}_driver
      PRIVATE typed_linear_algebra_options typed_linear_algebra_main
              typed_linear_algebra_${BACKEND})
    set_target_properties(
      typed_linear_algebra_test_${BACKEND}_${CALLER}_${TEST_NAME}_driver
      PROPERTIES EXCLUDE_FROM_ALL TRUE)
    separate_arguments(TEST_COMMAND UNIX_COMMAND $ENV{COMMAND})
    add_test(
      NAME typed_linear_algebra_test_${BACKEND}_${CALLER}_${TEST_NAME}
      COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target
              typed_linear_algebra_test_${BACKEND}_${CALLER}_${TEST_NAME}_driver
    )
    set_tests_properties(
      typed_linear_algebra_test_${BACKEND}_${CALLER}_${TEST_NAME}
      PROPERTIES WILL_FAIL TRUE)
  endforeach()
endfunction(fail)
