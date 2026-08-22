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

#ifndef FCAROUGE_AU_HPP
#define FCAROUGE_AU_HPP

//! @file
//! @brief Quantities and units facade for the Au third party implementation.
//!
//! @details Supporting the conversion of typed matrix elements to and from
//! Au's quantity type.
//!
//! @note Unlike mp-units, Au has no vector or matrix quantity representation
//! (tracked upstream as
//! [aurora-opensource/au#70](https://github.com/aurora-opensource/au/issues/70)):
//! an Au `Quantity<Unit, EigenMatrixType>` means one shared unit for an
//! entire matrix, not an independently typed element per cell. This library's
//! typed matrix model only ever needs Au's plain, scalar `Quantity<Unit,
//! Rep>`, converted one element at a time through `element_caster`: no
//! tuple-like decomposition or representation customization is needed, unlike
//! `fcarouge/mp_units.hpp`.

#include <au/au.hh>

#include <concepts>
#include <type_traits>

namespace fcarouge {
// True for an Au quantity, regardless of its unit or representation.
template <typename T>
concept au_quantity =
    requires {
      typename std::remove_cvref_t<T>::Unit;
      typename std::remove_cvref_t<T>::Rep;
    } and std::same_as<std::remove_cvref_t<T>,
                       au::Quantity<typename std::remove_cvref_t<T>::Unit,
                                    typename std::remove_cvref_t<T>::Rep>>;

// Teach the typed linear algebra library how to convert underlying scalar
// types to and from Au's quantity type.
template <typename To, au_quantity From> struct element_caster<To, From> {
  [[nodiscard]] static constexpr auto operator()(From value) -> To {
    static_assert(std::same_as<typename From::Rep, std::remove_cvref_t<To>>,
                  "The underlying storage type must be identical to the "
                  "quantity representation type to guarantee the conversion "
                  "is explicitely decided by the end-user.");

    return value.in(typename From::Unit{});
  }
};

template <au_quantity To, typename From> struct element_caster<To, From> {
  [[nodiscard]] static constexpr auto operator()(From value) -> To {
    static_assert(std::same_as<typename To::Rep, std::remove_cvref_t<From>>,
                  "The underlying storage type must be identical to the "
                  "quantity representation type to guarantee the conversion "
                  "is explicitely decided by the end-user.");

    return au::make_quantity<typename To::Unit>(value);
  }
};

template <au_quantity To, typename From> struct element_caster<To &, From &> {
  // A quantity reference cannot be safely materialized out of a
  // representation reference. It would be undefined behavior even if the
  // size, padding, alignment, aliasing are controlled. Therefore the best we
  // can do is to return a constant quantity value to inform the end-user
  // lvalue reference assignment cannot be supported.
  [[nodiscard]] static constexpr auto operator()(From value) -> const To {
    static_assert(std::same_as<typename To::Rep, std::remove_cvref_t<From>>,
                  "The underlying storage type must be identical to the "
                  "quantity representation type to guarantee the conversion "
                  "is explicitely decided by the end-user.");

    return au::make_quantity<typename To::Unit>(value);
  }
};
} // namespace fcarouge

#endif // FCAROUGE_AU_HPP
