#[[ Typed Linear Algebra
Version 0.3.0
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

# Writes documentation/test_matrix.html: the machine-owned body of the Test
# Matrix page -- a self-contained HTML fragment (its own <style>, then the
# operation x backend support grid built from the pass()/fail() lines under
# test/). documentation/test_matrix.md owns the title, prose, and legend and
# pulls this in verbatim with `@htmlinclude test_matrix.html`, so a regenerated
# grid never churns the prose.
#
# Run with `cmake -P`; CTest runs it after every test session
# (CTEST_CUSTOM_POST_TEST, configured from CTestCustom.cmake.in).

cmake_minimum_required(VERSION 3.20)

set(_out "${CMAKE_CURRENT_LIST_DIR}/test_matrix.html")

# operation|backend|reason a combination cannot exist. Reason: no ';', '|', '"'.
set(_na
    "addition|au_std|Au's lvalue-only Quantity::operator= is incompatible with the std::linalg-style add() (use operator+)"
    "substraction|au_std|Au's lvalue-only Quantity::operator= is incompatible with the std::linalg-style substract() (use operator-)"
)

file(GLOB _lists "${CMAKE_CURRENT_LIST_DIR}/../test/*/CMakeLists.txt")

set(_ops "")
set(_backends "")
foreach(_list IN LISTS _lists)
  get_filename_component(_op "${_list}" DIRECTORY)
  get_filename_component(_op "${_op}" NAME)

  file(STRINGS "${_list}" _decls REGEX "^[ \t]*(pass|fail)\\(\"")
  foreach(_decl IN LISTS _decls)
    string(REGEX MATCH "^[ \t]*(pass|fail)" _ "${_decl}")
    set(_kind "${CMAKE_MATCH_1}")
    string(REGEX MATCHALL "\"[^\"]+\"" _quoted "${_decl}")
    string(REPLACE "\"" "" _quoted "${_quoted}")
    list(POP_FRONT _quoted _name)

    # <rows>x<cols> shape token(s) in the file name; "any" if none
    string(REGEX MATCHALL "[0-9]+x[0-9]+" _shapes "${_name}")
    if(_shapes)
      list(REMOVE_DUPLICATES _shapes)
      list(JOIN _shapes "/" _shape)
      string(REPLACE "x" "×" _shape "${_shape}")
    else()
      set(_shape "any")
    endif()
    if(_kind STREQUAL "fail")
      set(_shape "✗${_shape}")
    endif()

    foreach(_backend IN LISTS _quoted)
      list(APPEND _ops "${_op}")
      list(APPEND _backends "${_backend}")
      list(APPEND "_cell_${_op}_${_backend}" "${_shape}")
    endforeach()
  endforeach()
endforeach()

list(REMOVE_DUPLICATES _ops)
list(SORT _ops)
list(REMOVE_DUPLICATES _backends)
list(SORT _backends) # alphabetical keeps each unit system's backends adjacent

# Doxygen splices this in verbatim (@htmlinclude), so it carries its own style.
set(_html
    [[<style>
#tm{border-collapse:collapse;font:12px/1.3 system-ui,sans-serif}
#tm th,#tm td{border:1px solid rgba(127,127,127,.35);padding:2px 6px}
#tm thead th:not(:first-child){writing-mode:vertical-rl;transform:rotate(180deg);
  padding:6px 3px;white-space:nowrap;vertical-align:bottom;font-weight:600}
#tm tbody th{text-align:left;white-space:nowrap;font-weight:600}
#tm td{text-align:center}
#tm .y{color:#2e9a3e;font-weight:700}
#tm .n{color:rgba(127,127,127,.6)}
</style>
<div style="overflow-x:auto"><table id="tm"><thead><tr><th>operation</th>]])
foreach(_backend IN LISTS _backends)
  string(APPEND _html "<th>${_backend}</th>")
endforeach()
string(APPEND _html "</tr></thead><tbody>\n")

foreach(_op IN LISTS _ops)
  string(APPEND _html "<tr><th>${_op}</th>")
  foreach(_backend IN LISTS _backends)
    if(DEFINED "_cell_${_op}_${_backend}")
      set(_shapes "${_cell_${_op}_${_backend}}")
      list(REMOVE_DUPLICATES _shapes)
      list(SORT _shapes) # pass shapes sort before the ✗ (compile-fail) ones
      list(JOIN _shapes " " _hint)
      string(APPEND _html "<td class=\"y\" title=\"${_hint}\">&#9679;</td>")
    else()
      set(_reason "")
      foreach(_entry IN LISTS _na)
        if(_entry MATCHES "^${_op}\\|${_backend}\\|(.*)$")
          set(_reason "${CMAKE_MATCH_1}")
        endif()
      endforeach()
      if(_reason)
        string(REPLACE "&" "&amp;" _reason "${_reason}")
        string(REPLACE "<" "&lt;" _reason "${_reason}")
        string(REPLACE ">" "&gt;" _reason "${_reason}")
        string(APPEND _html "<td class=\"n\" title=\"${_reason}\">&#8709;</td>")
      else()
        string(APPEND _html "<td></td>")
      endif()
    endif()
  endforeach()
  string(APPEND _html "</tr>\n")
endforeach()
string(APPEND _html "</tbody></table></div>\n")

# Only touch the file when the grid actually changed, to keep mtimes (and hence
# doc rebuilds) stable across no-op `ctest` runs.
set(_current "")
if(EXISTS "${_out}")
  file(READ "${_out}" _current)
endif()
if(_html STREQUAL _current)
  message(STATUS "test_matrix.cmake: ${_out} already current")
else()
  file(WRITE "${_out}" "${_html}")
  message(STATUS "test_matrix.cmake: wrote ${_out}")
endif()
