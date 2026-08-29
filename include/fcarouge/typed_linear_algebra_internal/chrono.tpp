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

#ifndef FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_CHRONO_TPP
#define FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_CHRONO_TPP

//! @file
//! @brief Time quantities facade for the standard `<chrono>` library.
//!
//! @details Supporting the conversion of typed matrix elements to and from
//! `std::chrono` durations and time points. A duration or time point carries
//! its representation and its period, so a typed vector element keeps its unit
//! (seconds, minutes, ...) while the backend stores the bare representation.
//!
//! @note The standard `<chrono>` library ships with every conforming C++
//! implementation, so this facade is part of the core header rather than an
//! external backend plug-in under `support/`.

#include "fcarouge/typed_linear_algebra_forward.hpp"

#include <chrono>
#include <concepts>
#include <type_traits>

namespace fcarouge {
// True for a `std::chrono::duration`, regardless of its period or
// representation.
template <typename T>
concept chrono_duration = std::same_as<
    std::remove_cvref_t<T>,
    std::chrono::duration<typename std::remove_cvref_t<T>::rep,
                          typename std::remove_cvref_t<T>::period>>;

// True for a `std::chrono::time_point`, regardless of its clock or duration.
template <typename T>
concept chrono_time_point = std::same_as<
    std::remove_cvref_t<T>,
    std::chrono::time_point<typename std::remove_cvref_t<T>::clock,
                            typename std::remove_cvref_t<T>::duration>>;

// Teach the typed linear algebra library how to convert underlying scalar
// types to and from `std::chrono` durations.
template <typename To, chrono_duration From> struct element_caster<To, From> {
  [[nodiscard]] static constexpr auto operator()(From value) -> To {
    static_assert(std::same_as<typename From::rep, std::remove_cvref_t<To>>,
                  "The underlying storage type must be identical to the "
                  "duration representation type to guarantee the conversion is "
                  "explicitely decided by the end-user.");

    return value.count();
  }
};

template <chrono_duration To, typename From> struct element_caster<To, From> {
  [[nodiscard]] static constexpr auto operator()(From value) -> To {
    static_assert(std::same_as<typename To::rep, std::remove_cvref_t<From>>,
                  "The underlying storage type must be identical to the "
                  "duration representation type to guarantee the conversion is "
                  "explicitely decided by the end-user.");

    return To{value};
  }
};

template <chrono_duration To, typename From>
struct element_caster<To &, From &> {
  // A duration reference cannot be safely materialized out of a representation
  // reference. The best we can do is a constant duration value to inform the
  // end-user lvalue reference assignment cannot be supported.
  [[nodiscard]] static constexpr auto operator()(From value) -> const To {
    static_assert(std::same_as<typename To::rep, std::remove_cvref_t<From>>,
                  "The underlying storage type must be identical to the "
                  "duration representation type to guarantee the conversion is "
                  "explicitely decided by the end-user.");

    return To{value};
  }
};

// Teach the typed linear algebra library how to convert underlying scalar
// types to and from `std::chrono` time points, measured from their epoch.
template <typename To, chrono_time_point From> struct element_caster<To, From> {
  [[nodiscard]] static constexpr auto operator()(From value) -> To {
    static_assert(
        std::same_as<typename From::rep, std::remove_cvref_t<To>>,
        "The underlying storage type must be identical to the time point "
        "representation type to guarantee the conversion is explicitely "
        "decided by the end-user.");

    return value.time_since_epoch().count();
  }
};

template <chrono_time_point To, typename From> struct element_caster<To, From> {
  [[nodiscard]] static constexpr auto operator()(From value) -> To {
    static_assert(
        std::same_as<typename To::rep, std::remove_cvref_t<From>>,
        "The underlying storage type must be identical to the time point "
        "representation type to guarantee the conversion is explicitely "
        "decided by the end-user.");

    return To{typename To::duration{value}};
  }
};

template <chrono_time_point To, typename From>
struct element_caster<To &, From &> {
  // See the duration reference caster: a time point reference cannot be safely
  // materialized out of a representation reference.
  [[nodiscard]] static constexpr auto operator()(From value) -> const To {
    static_assert(
        std::same_as<typename To::rep, std::remove_cvref_t<From>>,
        "The underlying storage type must be identical to the time point "
        "representation type to guarantee the conversion is explicitely "
        "decided by the end-user.");

    return To{typename To::duration{value}};
  }
};
} // namespace fcarouge

#endif // FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_CHRONO_TPP
