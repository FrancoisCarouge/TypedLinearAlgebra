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

# Amalgamate the typed linear algebra headers into a single, dependency-free
# header, suitable for pasting into or `#include`-by-URL from tools such as
# Compiler Explorer that cannot resolve this project's multi-file layout.
#
# Recursively inlines every quoted `#include "..."` starting from INPUT,
# resolving each one first relative to the including file's own directory and
# then against each root in INCLUDE_DIRS, in order, and writes the result to
# OUTPUT. Angle-bracket `#include <...>` directives (standard library and
# third-party headers, such as Eigen or mp-units) are left untouched.
#
# Files are read and rewritten as whole strings, not split into per-line
# lists: `file(STRINGS)` has been observed to corrupt long, deeply nested
# statements (for example the constexpr `for_constexpr` lambdas in
# `format.tpp`), silently fusing several source lines into one.
#
# Run in script mode:
#   cmake -D "INPUT=<entry header>" -D "INCLUDE_DIRS=<root>;<root>;..." \
#         -D "OUTPUT=<generated header>" -P "amalgamate.cmake"

if(NOT INPUT OR NOT INCLUDE_DIRS OR NOT OUTPUT)
  message(FATAL_ERROR "amalgamate.cmake requires INPUT, INCLUDE_DIRS, and OUTPUT")
endif()

get_filename_component(TOP_FILE "${INPUT}" ABSOLUTE)

set_property(GLOBAL PROPERTY AMALGAMATE_SEEN "")

# Compute the fully inlined content of FILE_PATH and store it in OUT_VAR. A
# file already inlined elsewhere yields an empty string so it is not
# duplicated.
function(amalgamate_get_inlined FILE_PATH OUT_VAR)
  get_filename_component(ABSOLUTE_PATH "${FILE_PATH}" ABSOLUTE)

  get_property(SEEN GLOBAL PROPERTY AMALGAMATE_SEEN)
  list(FIND SEEN "${ABSOLUTE_PATH}" SEEN_INDEX)
  if(NOT SEEN_INDEX EQUAL -1)
    set("${OUT_VAR}" "" PARENT_SCOPE)
    return()
  endif()
  set_property(GLOBAL APPEND PROPERTY AMALGAMATE_SEEN "${ABSOLUTE_PATH}")

  if(NOT EXISTS "${ABSOLUTE_PATH}")
    message(FATAL_ERROR "amalgamate: cannot find header '${ABSOLUTE_PATH}'")
  endif()

  file(READ "${ABSOLUTE_PATH}" CONTENT)
  get_filename_component(FILE_DIR "${ABSOLUTE_PATH}" DIRECTORY)

  # Every header carries the same leading Unlicense block comment. Keep it
  # once, on the top-level file, and drop it from every other inlined file.
  if(NOT ABSOLUTE_PATH STREQUAL TOP_FILE AND CONTENT MATCHES "^/\\*")
    string(FIND "${CONTENT}" "*/" COMMENT_END)
    math(EXPR COMMENT_END "${COMMENT_END} + 2")
    string(SUBSTRING "${CONTENT}" ${COMMENT_END} -1 CONTENT)
    string(REGEX REPLACE "^[\r\n]+" "" CONTENT "${CONTENT}")
  endif()

  # Recursively inline every local #include, one at a time: each iteration
  # re-scans the (partially expanded) content for the next remaining quoted
  # include directive belonging to the *original* file, since freshly
  # inlined content never itself contains a quoted include (it was fully
  # expanded before being substituted in).
  while(CONTENT MATCHES "\n[ \t]*#[ \t]*include[ \t]*\"([^\"]+)\"[^\n]*")
    set(DIRECTIVE "${CMAKE_MATCH_0}")
    set(QUOTED "${CMAKE_MATCH_1}")

    if(EXISTS "${FILE_DIR}/${QUOTED}")
      amalgamate_get_inlined("${FILE_DIR}/${QUOTED}" INLINED)
    else()
      set(RESOLVED "")
      foreach(ROOT IN LISTS INCLUDE_DIRS)
        if(EXISTS "${ROOT}/${QUOTED}")
          set(RESOLVED "${ROOT}/${QUOTED}")
          break()
        endif()
      endforeach()
      if(NOT RESOLVED)
        message(
          FATAL_ERROR
          "amalgamate: cannot resolve #include \"${QUOTED}\" from '${ABSOLUTE_PATH}'"
        )
      endif()
      amalgamate_get_inlined("${RESOLVED}" INLINED)
    endif()

    string(REPLACE "${DIRECTIVE}" "\n${INLINED}" CONTENT "${CONTENT}")
  endwhile()

  set("${OUT_VAR}" "${CONTENT}" PARENT_SCOPE)
endfunction()

amalgamate_get_inlined("${TOP_FILE}" RESULT)

get_filename_component(OUTPUT_DIR "${OUTPUT}" DIRECTORY)
get_filename_component(OUTPUT_NAME "${OUTPUT}" NAME)
file(MAKE_DIRECTORY "${OUTPUT_DIR}")
file(
  WRITE "${OUTPUT}"
  "// Generated file: single-header amalgamation `${OUTPUT_NAME}`.\n"
  "// Do not edit directly: edit the sources under `include/` and\n"
  "// `support/` and regenerate with the `amalgamate` CMake target.\n"
  "// https://github.com/FrancoisCarouge/TypedLinearAlgebra\n\n"
  "${RESULT}")
