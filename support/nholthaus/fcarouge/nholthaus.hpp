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

#ifndef FCAROUGE_NHOLTHAUS_HPP
#define FCAROUGE_NHOLTHAUS_HPP

//! @file
//! @brief Quantities and units facade for the nholthaus/units implementation.
//!
//! @details Supporting the conversion of typed matrix elements to and from
//! nholthaus/units' `units::unit` container type.

#include "fcarouge/typed_linear_algebra_forward.hpp"

#include <units/core.h>

#include <concepts>
#include <type_traits>

namespace fcarouge {
// True for a nholthaus/units quantity, regardless of its conversion factor,
// representation, or numerical scale. `is_unit_v` already excludes plain
// arithmetic types.
template <typename T>
concept nholthaus_quantity = units::traits::is_unit_v<std::remove_cvref_t<T>>;

// The arithmetic storage type backing a nholthaus/units quantity (e.g.
// `double` for `units::length::meters<double>`).
template <nholthaus_quantity Quantity>
using nholthaus_representation = typename units::traits::unit_traits<
    std::remove_cvref_t<Quantity>>::underlying_type;

// Teach the typed linear algebra library how to convert underlying scalar
// types to and from nholthaus/units' quantity type.
template <typename To, nholthaus_quantity From>
struct element_caster<To, From> {
  [[nodiscard]] static constexpr auto operator()(From value) -> To {
    static_assert(
        std::same_as<nholthaus_representation<From>, std::remove_cvref_t<To>>,
        "The underlying storage type must be identical to the "
        "quantity representation type to guarantee the conversion "
        "is explicitly decided by the end-user.");

    return value.value();
  }
};

template <nholthaus_quantity To, typename From>
struct element_caster<To, From> {
  [[nodiscard]] static constexpr auto operator()(From value) -> To {
    static_assert(
        std::same_as<nholthaus_representation<To>, std::remove_cvref_t<From>>,
        "The underlying storage type must be identical to the "
        "quantity representation type to guarantee the conversion "
        "is explicitly decided by the end-user.");

    return To(value);
  }
};

template <nholthaus_quantity To, typename From>
struct element_caster<To &, From &> {
  // A quantity reference cannot be safely materialized out of a
  // representation reference. It would be undefined behavior even if the
  // size, padding, alignment, aliasing are controlled. Therefore the best we
  // can do is to return a constant quantity value to inform the end-user
  // lvalue reference assignment cannot be supported.
  [[nodiscard]] static constexpr auto operator()(From value) -> const To {
    static_assert(
        std::same_as<nholthaus_representation<To>, std::remove_cvref_t<From>>,
        "The underlying storage type must be identical to the "
        "quantity representation type to guarantee the conversion "
        "is explicitly decided by the end-user.");

    return To(value);
  }
};
} // namespace fcarouge

#endif // FCAROUGE_NHOLTHAUS_HPP
