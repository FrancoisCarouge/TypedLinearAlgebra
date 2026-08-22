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

#ifndef FCAROUGE_UNIT_HPP
#define FCAROUGE_UNIT_HPP

//! @file
//! @brief Quantities and units facade for mp-units third party implementation.
//!
//! @details Supporting quantities, values, and functions.

#include <concepts>
#include <cstddef>
#include <tuple>

#include <mp-units/framework/quantity.h>
#include <mp-units/framework/quantity_cast.h>
#include <mp-units/framework/quantity_point.h>
#include <mp-units/framework/vector_components.h>
#include <mp-units/math.h>
#include <mp-units/systems/isq/thermodynamics.h>
#include <mp-units/systems/si.h>

namespace fcarouge {
using mp_units::si::unit_symbols::A;
using mp_units::si::unit_symbols::h;
using mp_units::si::unit_symbols::km;
using mp_units::si::unit_symbols::m;
using mp_units::si::unit_symbols::m2;
using mp_units::si::unit_symbols::m3;
using mp_units::si::unit_symbols::m4;
using mp_units::si::unit_symbols::mol;
using mp_units::si::unit_symbols::N;
using mp_units::si::unit_symbols::s;
using mp_units::si::unit_symbols::s2;
using mp_units::si::unit_symbols::s3;

inline constexpr auto s4{pow<4>(s)};

using height = mp_units::quantity<mp_units::isq::height[m]>;
using position = mp_units::quantity<mp_units::isq::length[m]>;
using velocity = mp_units::quantity<mp_units::isq::velocity[m / s]>;
using acceleration = mp_units::quantity<mp_units::isq::acceleration[m / s2]>;

// A vector quantity whose specification has no named-axis decomposition (no
// `vector_components` specialization) but whose own representation is itself
// tuple-like (e.g. a fcarouge typed vector). `not requires { ...::size; }`
// mirrors mp-units's own `Decomposable` prerequisite check and is therefore
// guaranteed mutually exclusive with it: mp-units's `std::tuple_size` and
// `get<Idx>` specializations for `mp_units::quantity` require `Decomposable`,
// which itself requires `vector_components<Spec>::size` to exist.
template <auto Reference, typename Representation>
concept undecomposed_tuple_like_quantity =
    (not requires {
      mp_units::vector_components<get_quantity_spec(Reference)>::size;
    }) and requires(const Representation &representation) {
      {
        std::tuple_size<Representation>::value
      } -> std::convertible_to<std::size_t>;
      get<0>(representation);
    };
} // namespace fcarouge

// Extend the standard tuple-like protocol (`std::tuple_size`,
// `std::tuple_element`, and the ADL-found `get<Idx>` below) to vector
// quantities left undecomposed by mp-units itself, forwarding to their own
// tuple-like representation. This is what lets `fcarouge::typed_matrix`'s
// core, vanilla-C++ `other_tuple_like_vector` converting constructor accept
// such a quantity directly: `get<Idx>` is found purely through
// argument-dependent lookup on `Representation` (which names `fcarouge`, and
// the underlying matrix backend), so the core library never names any
// mp-units-specific member.
template <auto Reference, typename Representation>
  requires fcarouge::undecomposed_tuple_like_quantity<Reference, Representation>
struct std::tuple_size<mp_units::quantity<Reference, Representation>>
    : std::tuple_size<Representation> {};

template <std::size_t Index, auto Reference, typename Representation>
  requires fcarouge::undecomposed_tuple_like_quantity<Reference,
                                                      Representation> and
           (Index < std::tuple_size_v<Representation>)
struct std::tuple_element<Index,
                          mp_units::quantity<Reference, Representation>> {
  using type = mp_units::quantity<Reference,
                                  std::tuple_element_t<Index, Representation>>;
};

namespace fcarouge {
template <std::size_t Index, auto Reference, typename Representation>
  requires undecomposed_tuple_like_quantity<Reference, Representation> and
           (Index < std::tuple_size_v<Representation>)
[[nodiscard]] constexpr auto
get(const mp_units::quantity<Reference, Representation> &value) {
  return get<Index>(value.numerical_value_in(value.unit)) * value.unit;
}

// Teach the typed linear algebra library how to convert underlying scalar types
// to and from mp-units' types. Support castable quantity specifications.
template <typename To, mp_units::Quantity From>
struct element_caster<To, From> {
  [[nodiscard]] static constexpr auto operator()(mp_units::Quantity auto value)
      -> To {
    using representation = typename std::remove_cvref_t<From>::rep;

    static_assert(std::same_as<representation, std::remove_cvref_t<To>>,
                  "The underlying storage type must be identical to the "
                  "quantity representation type to guarantee the conversion is "
                  "explicitely decided by the end-user.");

    constexpr auto reference{std::remove_cvref_t<From>::reference};
    constexpr auto specification{get_quantity_spec(reference)};
    auto specified{mp_units::quantity_cast<specification>(value)};

    return specified.numerical_value_in(specified.unit);
  }
};

template <typename To, mp_units::UnitMagnitude Magnitude>
struct element_caster<To, Magnitude> {
  [[nodiscard]] static constexpr auto operator()(Magnitude value) -> To {
    return get_value<To>(value);
  }
};

template <mp_units::Quantity To, typename From>
struct element_caster<To, From> {
  [[nodiscard]] static constexpr auto operator()(From value) -> To {
    using representation = typename std::remove_cvref_t<To>::rep;

    static_assert(std::same_as<representation, std::remove_cvref_t<From>>,
                  "The underlying storage type must be identical to the "
                  "quantity representation type to guarantee the conversion is "
                  "explicitely decided by the end-user.");

    return value * To::reference;
  }
};

template <mp_units::Quantity To, typename From>
struct element_caster<To &, From &> {
  // A quantity reference cannot be safely materialized out of a representation
  // reference. It would be undefined behavior even if the size, padding,
  // alignment, aliasing are controlled. Therefore the best we can do is to
  // return a constant quantity value to inform the end-user lvalue reference
  // assignment cannot be supported.
  [[nodiscard]] static constexpr auto operator()(From value) -> const To {
    using representation = typename std::remove_cvref_t<To>::rep;

    static_assert(std::same_as<representation, std::remove_cvref_t<From>>,
                  "The underlying storage type must be identical to the "
                  "quantity representation type to guarantee the conversion is "
                  "explicitely decided by the end-user.");

    return value * To::reference;
  }
};

template <typename To, mp_units::QuantityPoint From>
struct element_caster<To, From> {
  [[nodiscard]] static constexpr auto operator()(From value) -> To {
    using representation = typename std::remove_cvref_t<From>::rep;

    static_assert(std::same_as<representation, std::remove_cvref_t<To>>,
                  "The underlying storage type must be identical to the "
                  "quantity representation type to guarantee the conversion is "
                  "explicitely decided by the end-user.");

    return value.quantity_from_zero().numerical_value_in(value.unit);
  }
};

template <mp_units::QuantityPoint To, typename From>
struct element_caster<To, From> {
  [[nodiscard]] static constexpr auto operator()(From value) -> To {
    using representation = typename std::remove_cvref_t<To>::rep;

    static_assert(std::same_as<representation, std::remove_cvref_t<From>>,
                  "The underlying storage type must be identical to the "
                  "quantity representation type to guarantee the conversion is "
                  "explicitely decided by the end-user.");

    return {value * To::unit, mp_units::default_point_origin(To::unit)};
  }
};

template <mp_units::QuantityPoint To, typename From>
struct element_caster<To &, From &> {
  // A quantity point reference cannot be safely materialized out of a
  // representation reference. It would be undefined behavior even if the size,
  // padding, alignment, aliasing are controlled. Therefore the best we can do
  // is to return a constant quantity value to inform the end-user lvalue
  // reference assignment cannot be supported.
  [[nodiscard]] static constexpr auto operator()(From value) -> const To {
    using representation = typename std::remove_cvref_t<To>::rep;

    static_assert(std::same_as<representation, std::remove_cvref_t<From>>,
                  "The underlying storage type must be identical to the "
                  "quantity representation type to guarantee the conversion is "
                  "explicitely decided by the end-user.");

    return {value * To::unit, mp_units::default_point_origin(To::unit)};
  }
};

template <typename To, mp_units::Reference From>
struct element_caster<To, From> {
  [[nodiscard]] static constexpr auto operator()([[maybe_unused]] From value)
      -> To {
    return 1.;
  }
};

template <typename Type, mp_units::UnitMagnitude Magnitude>
struct multiplies<Type, Magnitude> {
  // A type that provides its own magnitude-aware operator*(T, UnitMagnitude)
  // customization point. scale() will prefer this path when available, and the
  // return type may differ from the input (e.g. a wrapper with scaled bounds).
  [[nodiscard]] static constexpr auto operator()(const Type &lhs,
                                                 const Magnitude &rhs) -> Type;
};
} // namespace fcarouge

#endif // FCAROUGE_UNIT_HPP
