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

FetchContent_Declare(
  boost-multi
  GIT_REPOSITORY "https://gitlab.com/correaa/boost-multi.git"
  GIT_TAG "develop"
  GIT_SHALLOW TRUE
  SOURCE_SUBDIR "include"
  SYSTEM
  FIND_PACKAGE_ARGS NAMES multi)
FetchContent_MakeAvailable(boost-multi)

add_library(typed_linear_algebra_boost_multi INTERFACE)
target_include_directories(
  typed_linear_algebra_boost_multi
  SYSTEM
  INTERFACE "${boost-multi_SOURCE_DIR}/include")

foreach(TYPES IN ITEMS "au" "chrono" "mp_units" "nholthaus")
  add_library(typed_linear_algebra_${TYPES}_std_boost_multi INTERFACE)
  target_link_libraries(
    typed_linear_algebra_${TYPES}_std_boost_multi
    INTERFACE typed_linear_algebra_${TYPES}_std typed_linear_algebra_boost_multi)
endforeach()
