/* Typed Linear Algebra
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

For more information, please refer to <https://unlicense.org> */

#ifndef FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_CAST_TPP
#define FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_CAST_TPP

namespace fcarouge {
namespace tla = typed_linear_algebra_internal;

//! @todo Deduplicate, generalize the built-in casts.
template <typename To, typename From>
[[nodiscard]] constexpr To
element_caster<To, From>::operator()(From value) const {
  return value;
}

template <typename To, typename From> struct element_caster<To &, From &> {
  [[nodiscard]] constexpr To &operator()(From &value) const { return value; }
  [[nodiscard]] constexpr To &&operator()(From &&value) const {
    return std::move(value);
  }
};

template <typename To, typename From> struct element_caster<To &&, From &&> {
  [[nodiscard]] constexpr To &&operator()(From &&value) const {
    return std::move(value);
  }
};

template <typename To, typename From>
  requires requires(const From &value) { value(0, 0); }
struct element_caster<To, From> {
  [[nodiscard]] constexpr const To &operator()(const From &value) const {
    return value(0, 0);
  }
};
} // namespace fcarouge

#endif // FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_CAST_TPP
