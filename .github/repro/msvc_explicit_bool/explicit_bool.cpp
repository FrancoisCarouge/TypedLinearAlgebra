/* Typed Linear Algebra
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

For more information, please refer to <https://unlicense.org> */

// MSVC drops the value-dependent `explicit(bool)` specifier of a member
// template once its out-of-class definition is seen: the definition, which
// does not and cannot repeat the specifier, makes the member implicit.
// See [class.conv.ctor], [class.conv.fct], and [dcl.fct.spec].

#include <type_traits>

// A constructor template, implicit from `int` only.
struct constructor {
  template <typename Type>
  explicit(not std::is_same_v<Type, int>) constructor(const Type &value);
};

static_assert(not std::is_convertible_v<double, constructor>);

template <typename Type>
constructor::constructor([[maybe_unused]] const Type &value) {}

// Fails on MSVC: the conversion became implicit.
static_assert(not std::is_convertible_v<float, constructor>);
static_assert(std::is_convertible_v<int, constructor>);

// A conversion function template, implicit to `int` only.
struct conversion {
  template <typename Type>
  explicit(not std::is_same_v<Type, int>) operator Type() const;
};

template <typename Type> conversion::operator Type() const { return Type{}; }

// Fails on MSVC: the conversion became implicit.
static_assert(not std::is_convertible_v<conversion, double>);
static_assert(std::is_convertible_v<conversion, int>);

// The same constructor template, defined in the class: correct on MSVC.
struct in_class {
  template <typename Type>
  explicit(not std::is_same_v<Type, int>)
      in_class([[maybe_unused]] const Type &value) {}
};

static_assert(not std::is_convertible_v<double, in_class>);
static_assert(std::is_convertible_v<int, in_class>);
