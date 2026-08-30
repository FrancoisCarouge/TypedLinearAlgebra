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

# Writes documentation/test_matrix.md: a Doxygen page whose body is an HTML
# support grid of operation x backend, from the pass()/fail() lines under test/.
# Run with `cmake -P`; CTest runs it after every test session
# (CTEST_CUSTOM_POST_TEST, wired in test/CMakeLists.txt).

cmake_minimum_required(VERSION 3.20)

get_filename_component(_root "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
set(_out "${CMAKE_CURRENT_LIST_DIR}/test_matrix.md")

# operation|backend|reason a combination cannot exist. Reason: no ';', '|', '"'.
set(_na
    "addition|au_std|Au's lvalue-only Quantity::operator= is incompatible with the std::linalg-style add() (use operator+)"
    "substraction|au_std|Au's lvalue-only Quantity::operator= is incompatible with the std::linalg-style substract() (use operator-)"
)

file(GLOB _lists "${_root}/test/*/CMakeLists.txt")

set(_ops "")
set(_backends "")
set(_files "")
set(_cases 0)
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

    list(APPEND _files "${_op}/${_name}")
    foreach(_backend IN LISTS _quoted)
      list(APPEND _ops "${_op}")
      list(APPEND _backends "${_backend}")
      list(APPEND "_cell_${_op}_${_backend}" "${_shape}")
      math(EXPR _cases "${_cases} + 1")
    endforeach()
  endforeach()
endforeach()

list(REMOVE_DUPLICATES _ops)
list(SORT _ops)
list(REMOVE_DUPLICATES _backends)
list(SORT _backends) # alphabetical keeps each unit system's backends adjacent
list(REMOVE_DUPLICATES _files)
list(LENGTH _ops _n_ops)
list(LENGTH _backends _n_backends)
list(LENGTH _files _n_files)

set(_md "# Test Matrix {#test_matrix}\n\n")
string(
  APPEND
  _md
  "Which operations are exercised against which backends, from the "
  "`pass()` and `fail()` declarations under `test/`. Regenerated after "
  "every `ctest` run. Hover a cell for the shapes covered (`✗` = a "
  "compile-`fail` test).\n\n")

# The grid is raw HTML so Doxygen renders it verbatim in both colour themes.
string(APPEND _md "@htmlonly\n")
string(
  APPEND
  _md
  [[<style>
#tm{border-collapse:collapse;font:12px/1.3 system-ui,sans-serif}
#tm th,#tm td{border:1px solid rgba(127,127,127,.35);padding:2px 6px}
#tm thead th:not(:first-child){writing-mode:vertical-rl;transform:rotate(180deg);
  vertical-align:bottom;font-weight:600}
#tm tbody th{text-align:left;white-space:nowrap;font-weight:600}
#tm td{text-align:center}
#tm .y{color:#2e9a3e;font-weight:700}
#tm .n{color:rgba(127,127,127,.6)}
</style>
<p><span class="y">&#9679;</span> tested&nbsp;&nbsp;
<span class="n">&#8709;</span> not applicable&nbsp;&nbsp;
blank = untested</p>
<div style="overflow-x:auto"><table id="tm"><thead><tr><th>operation</th>]])
foreach(_backend IN LISTS _backends)
  string(APPEND _md "<th>${_backend}</th>")
endforeach()
string(APPEND _md "</tr></thead><tbody>\n")

foreach(_op IN LISTS _ops)
  string(APPEND _md "<tr><th>${_op}</th>")
  foreach(_backend IN LISTS _backends)
    if(DEFINED "_cell_${_op}_${_backend}")
      set(_shapes "${_cell_${_op}_${_backend}}")
      list(REMOVE_DUPLICATES _shapes)
      list(SORT _shapes) # pass shapes sort before the ✗ (compile-fail) ones
      list(JOIN _shapes " " _hint)
      string(APPEND _md "<td class=\"y\" title=\"${_hint}\">&#9679;</td>")
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
        string(APPEND _md "<td class=\"n\" title=\"${_reason}\">&#8709;</td>")
      else()
        string(APPEND _md "<td></td>")
      endif()
    endif()
  endforeach()
  string(APPEND _md "</tr>\n")
endforeach()
string(APPEND _md "</tbody></table></div>\n@endhtmlonly\n\n")

string(APPEND _md "## Not applicable {#test_matrix_na}\n\n")
foreach(_entry IN LISTS _na)
  string(REPLACE "|" ";" _entry "${_entry}")
  list(GET _entry 0 _op)
  list(GET _entry 1 _backend)
  list(GET _entry 2 _reason)
  string(APPEND _md "- `${_op}` × `${_backend}` — ${_reason}\n")
endforeach()

string(APPEND _md
       "\n<em>${_n_files} test files, ${_cases} compiled cases, ${_n_ops} "
       "operations × ${_n_backends} backends.</em>\n")

file(WRITE "${_out}" "${_md}")
message(STATUS "test_matrix.cmake: wrote ${_out}")
