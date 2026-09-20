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

# Writes documentation/test_matrix.html: the machine-owned body of the Test
# Matrix page -- a self-contained HTML fragment (its own <style>, then the
# operation x backend support grid built from the pass() lines under test/).
# build() and fail() declarations are not support signals (the former is a
# static_assert-only compile check, the latter an intentional compile failure),
# so neither is reported. documentation/test_matrix.md owns the title, prose,
# and legend and pulls this in verbatim with `@htmlinclude test_matrix.html`, so
# a regenerated grid never churns the prose.
#
# Run with `cmake -P`; CTest runs it after every test session
# (CTEST_CUSTOM_POST_TEST, configured from CTestCustom.cmake.in).

cmake_minimum_required(VERSION 3.20)

set(_out "${CMAKE_CURRENT_LIST_DIR}/test_matrix.html")

# A backend name is <linear algebra library> or <strong type>_<linear algebra
# library> (e.g. "au_eigen" is the Au strong type over the Eigen backend). These
# backends have no strong type component, just the library itself.
set(_type_less_backends "eigen" "eigexed" "nested_typed_eigen" "armadillo"
                        "armadilloxed")

file(GLOB _lists "${CMAKE_CURRENT_LIST_DIR}/../test/*/CMakeLists.txt")

set(_ops "")
set(_backends "")
foreach(_list IN LISTS _lists)
  get_filename_component(_op "${_list}" DIRECTORY)
  get_filename_component(_op "${_op}" NAME)

  file(STRINGS "${_list}" _decls REGEX "^[ \t]*pass\\(\"")
  foreach(_decl IN LISTS _decls)
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

# Decompose each backend into its type ("au", "mp_units", ...; empty when the
# backend is type-less) and linear algebra library ("eigen", "std", ...), then
# sort by (library, type) so the library footer row below can span contiguous
# same-library columns.
set(_entries "")
foreach(_backend IN LISTS _backends)
  set(_type "")
  set(_library "${_backend}")
  if(NOT _backend IN_LIST _type_less_backends
     AND _backend MATCHES "^(.+)_(eigen|std|armadillo)$")
    set(_type "${CMAKE_MATCH_1}")
    set(_library "${CMAKE_MATCH_2}")
  endif()
  list(APPEND _entries "${_library}|${_type}|${_backend}")
endforeach()
list(SORT _entries)

set(_backends "")
set(_types "") # joined with trailing ";" per element: list(APPEND) on an empty
set(_libraries "") # list can't distinguish "no elements" from "one empty
foreach(_entry IN LISTS _entries) # element", silently dropping a leading "".
  string(REPLACE "|" ";" _fields "${_entry}")
  list(GET _fields 0 _library)
  list(GET _fields 1 _type)
  list(GET _fields 2 _backend)
  string(APPEND _backends "${_backend};")
  string(APPEND _types "${_type};")
  string(APPEND _libraries "${_library};")
endforeach()
string(REGEX REPLACE ";$" "" _backends "${_backends}")
string(REGEX REPLACE ";$" "" _types "${_types}")
string(REGEX REPLACE ";$" "" _libraries "${_libraries}")

# Doxygen splices this in verbatim (@htmlinclude), so it carries its own style.
set(_html
    [[<style>
#tm{border-collapse:collapse;font:12px/1.3 system-ui,sans-serif}
#tm th,#tm td{border:1px solid rgba(127,127,127,.35);padding:2px 6px}
#tm thead th:not(:first-child),#tm tfoot th:not(:first-child){
  writing-mode:vertical-rl;transform:rotate(180deg);padding:6px 3px;
  white-space:nowrap;vertical-align:middle;font-weight:600}
#tm thead th:not(:first-child){text-align:left}
#tm tfoot th:not(:first-child){text-align:right}
#tm tbody th{text-align:right;vertical-align:middle;white-space:nowrap;
  font-weight:600}
#tm td{text-align:center}
.y{background:#2e9a3e;color:#fff;font-weight:700;font-size:9px}
</style>
<div style="overflow-x:auto"><table id="tm"><thead><tr><th></th>]])

# Header: the strong type, one <th> per column.
foreach(_type IN LISTS _types)
  string(APPEND _html "<th>${_type}</th>")
endforeach()
string(APPEND _html "</tr></thead><tbody>\n")

foreach(_op IN LISTS _ops)
  string(APPEND _html "<tr><th>${_op}</th>")
  foreach(_backend IN LISTS _backends)
    if(DEFINED "_cell_${_op}_${_backend}")
      set(_all "${_cell_${_op}_${_backend}}")
      list(LENGTH _all _count)
      set(_shapes "${_all}")
      list(REMOVE_DUPLICATES _shapes)
      list(SORT _shapes)
      list(JOIN _shapes " " _hint)
      if(_count GREATER 1)
        string(APPEND _html "<td class=\"y\" title=\"${_hint}\">${_count}</td>")
      else()
        string(APPEND _html "<td class=\"y\" title=\"${_hint}\"></td>")
      endif()
    else()
      string(APPEND _html "<td></td>")
    endif()
  endforeach()
  string(APPEND _html "</tr>\n")
endforeach()
string(APPEND _html "</tbody><tfoot><tr><th></th>")

# Footer: the linear algebra library, below the data, one <th> spanning each
# contiguous run of columns sharing it.
list(LENGTH _libraries _n)
set(_i 0)
while(_i LESS _n)
  list(GET _libraries ${_i} _library)
  set(_span 1)
  math(EXPR _j "${_i}+1")
  while(_j LESS _n)
    list(GET _libraries ${_j} _next)
    if(NOT _next STREQUAL _library)
      break()
    endif()
    math(EXPR _span "${_span}+1")
    math(EXPR _j "${_j}+1")
  endwhile()
  if(_span GREATER 1)
    string(APPEND _html "<th colspan=\"${_span}\">${_library}</th>")
  else()
    string(APPEND _html "<th>${_library}</th>")
  endif()
  set(_i ${_j})
endwhile()
string(APPEND _html "</tr></tfoot></table></div>\n")

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
