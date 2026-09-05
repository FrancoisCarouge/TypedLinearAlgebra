/* Typed Linear Algebra
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

For more information, please refer to <https://unlicense.org> */

#include "fcarouge/linalg.hpp"

#include <string>
#include <tuple>

namespace fcarouge::test {
namespace {
//! @test Verifies the same_as_typed_matrix concept.

// Positive: every rank of typed matrix satisfies the concept.
using singleton = matrix<double, 1, 1>;
static_assert(same_as_typed_matrix<singleton>);
using row = matrix<double, 1, 3>;
static_assert(same_as_typed_matrix<row>);
using column = matrix<double, 3, 1>;
static_assert(same_as_typed_matrix<column>);
using full = matrix<double, 2, 3>;
static_assert(same_as_typed_matrix<full>);

// Positive, special case: cv-qualification and references are stripped
// before the check, `std::remove_cvref_t` is applied first.
using m = matrix<double, 2, 2>;
static_assert(same_as_typed_matrix<m>);
static_assert(same_as_typed_matrix<const m>);
static_assert(same_as_typed_matrix<volatile m>);
static_assert(same_as_typed_matrix<const volatile m>);
static_assert(same_as_typed_matrix<m &>);
static_assert(same_as_typed_matrix<const m &>);
static_assert(same_as_typed_matrix<m &&>);
static_assert(same_as_typed_matrix<const m &&>);

// Negative, special case: fundamental and unrelated types have no nested
// `matrix`, `row_indexes`, `column_indexes` member types for the concept to
// probe. That probe is in the immediate context of the constraint, so the
// substitution failure makes the concept simply unsatisfied; it does not
// make the program ill-formed. `void` is included since it cannot have
// members at all.
static_assert(not same_as_typed_matrix<void>);
static_assert(not same_as_typed_matrix<int>);
static_assert(not same_as_typed_matrix<double>);
static_assert(not same_as_typed_matrix<std::string>);
static_assert(not same_as_typed_matrix<std::tuple<double>>);

// Negative, special case: an incomplete type also has no accessible nested
// types yet, and probing one still fails safely instead of hard-erroring on
// "invalid use of incomplete type".
struct incomplete;
static_assert(not same_as_typed_matrix<incomplete>);

// Negative: a pointer to a typed matrix is not itself a typed matrix.
// `remove_cvref_t` strips references and top-level cv-qualification, not
// indirection.
using ptr_target = matrix<double, 1, 1>;
static_assert(not same_as_typed_matrix<ptr_target *>);
static_assert(not same_as_typed_matrix<const ptr_target *>);

// Negative: the underlying backend matrix, composed into a typed matrix but
// not itself one, does not satisfy the concept.
static_assert(not same_as_typed_matrix<eigen::matrix<double, 2, 2>>);

// Negative, special case: the concept is not structural, duck typing. A type
// exposing the same `matrix`, `row_indexes`, and `column_indexes` member
// aliases as a typed matrix, but that is not literally that `typed_matrix`
// class template instantiation, does not satisfy the concept:
// `std::same_as` requires nominal type identity.
struct fake {
  using matrix = eigen::matrix<double, 1, 1>;
  using row_indexes = std::tuple<double>;
  using column_indexes = std::tuple<double>;
};
static_assert(not same_as_typed_matrix<fake>);

// Negative, special case: a type publicly inheriting from a typed matrix
// inherits its `matrix`, `row_indexes`, `column_indexes` member aliases, so
// the concept reconstructs the *base* typed matrix from them, but
// `std::same_as` then compares the derived type itself against that
// reconstructed base. They are different types, so a derived type never
// satisfies `same_as_typed_matrix`, even though it "is a" typed matrix in the
// inheritance sense and adds no members of its own. One consequence: such a
// derived type instead satisfies `other`, so operators and constructors
// overloaded on `same_as_typed_matrix` versus `other` would treat it as
// "anything else", not as a typed matrix.
using base = matrix<double, 1, 1>;
struct derived : base {};
static_assert(same_as_typed_matrix<base>);
static_assert(not same_as_typed_matrix<derived>);
static_assert(other<derived>);

// Positive: `other` is exactly the negation of `same_as_typed_matrix`,
// practical for disambiguating overloads between a typed matrix and anything
// else.
static_assert(other<int>);
static_assert(not other<matrix<double, 1, 1>>);
} // namespace
} // namespace fcarouge::test
