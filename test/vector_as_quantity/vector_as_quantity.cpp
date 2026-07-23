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

#include <cassert>
#include <format>
#include <print>
#include <utility>

namespace mp_units {
// Just like Eigen, the fcarouge::typed_matrix
// arithmetic operators return lazy expression templates; store their evaluated
// concrete type (`PlainObject`) in a quantity instead. Concrete
// matrices/vectors map to themselves.
//
// The `typename T::PlainObject` requirement is checked first and short-circuits
// the rest of the constraint: `representation_canonical_type` is instantiated
// for every representation type (including `int`, `double`, ...), and
// `Eigen::EigenBase<T>` is ill-formed for a non-Eigen `T`, so it must not be
// instantiated unless `T` already looks like an Eigen type.
template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
  requires requires { typename Matrix::PlainObject; } &&
           std::derived_from<Matrix, Eigen::EigenBase<Matrix>>
// SIMLPIFY WITH CONCEPTS
struct representation_canonical_type<
    fcarouge::typed_matrix<Matrix, RowIndexes, ColumnIndexes>> {
  using type =
      fcarouge::typed_matrix<std::remove_cvref_t<typename Matrix::PlainObject>,
                             RowIndexes, ColumnIndexes>;
};

// Just like Blaze, the fcarouge::typed_matrix
// does not expose the `value_type`/`element_type` names the library detects
// automatically, so map its `ElementType` explicitly.
template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
  requires requires { typename Matrix::PlainObject; } &&
           std::derived_from<Matrix, Eigen::EigenBase<Matrix>>
struct representation_underlying_type<
    fcarouge::typed_matrix<Matrix, RowIndexes, ColumnIndexes>> {
  using type =
      fcarouge::typed_matrix<Matrix, RowIndexes, ColumnIndexes>::underlying;
};

// DOC HERE
// NOT QUITE... TRY AGAIN
template <auto QuantitySpecification>
  requires // mp_units::Quantity<Type> &&
    fcarouge::rank_typed_matrix<QuantitySpecification::rep, 1>
  struct mp_units::vector_components<QuantitySpecification> {
  static constexpr std::size_t size = 3;
  // TypedVector::rows * TypedVector::columns;
};

} // namespace mp_units

using namespace mp_units; // REMOVE?

// WIP SUPPORT - END ///////////////////////////////////////////////////////////

namespace fcarouge::test {
namespace {

// Inject the magnitude customization-point object (CPO, a variable with a call
// operator) so it competes with the library's magnitude free function.
using mp_units::magnitude;

// COMBINE THE MP-UNITS REFERENCE SAMPLES INTO ONE FOR EVERYONE.

//! @brief The mp-units reference for Eigen.
[[maybe_unused]] const auto reference{[] {
  quantity velocity = isq::velocity(Eigen::Vector3d{30, 40, 0} * km / h);

  // From there the framework treats the vector exactly like any other
  // representation:
  // * Unit conversion scales every component:
  quantity velocity_mps = velocity.in(m / s);
  // *vector quantity × scalar quantity:
  quantity displacement = velocity * (2 * h);
  // * Euclidean magnitude → scalar quantity:
  quantity speed = magnitude(velocity);

  assert(std::format("{}", velocity) == "[[30], [40], [0]] km/h");
  assert(std::format("{}", velocity_mps) == "[[8.33333], [11.1111], [0]] m/s");
  assert(std::format("{}", displacement) == "[[60], [80], [0]] km");
  assert(std::format("{}", speed) == "50 km/h");

  // OK: [[maybe_unused]] auto [x, y, z] = Eigen::Vector3d{30, 40, 0};

  // NOK? Why not?
  // [[maybe_unused]] auto [x, y, z] = velocity;

  // static_assert(detail::Decomposable<flight_velocity, vec3>);
  // const quantity forward = get<0>(velocity);

  return 0;
}()};

// MOVE STATIC ASSERT TO A REFERENCE.CPP TEST FILE.

// Compile-time guarantee tests from:
// mp-unit's test/runtime/linear_algebra_test.cpp
using vec3 =
    typed_column_vector<Eigen::Vector<double, 3>, double, double, double>;

// The library vector type is accepted as a vector representation.
static_assert(RepresentationOf<vec3, quantity_tensor_order::vector>);

// A library scalar multiplication may yield a lazy expression template; the
// representation machinery must canonicalize it back to the concrete vector
// type rather than store the proxy (which would dangle). `decltype(vec3{}
// * 2.0)` is exactly such a proxy for Eigen/Blaze.
static_assert(
    std::same_as<
        representation_canonical_type_t<decltype(std::declval<vec3>() * 2.0)>,
        vec3>);

// Consequently, arithmetic on vector quantities stores the concrete vector
// type, not a proxy.
static_assert(std::same_as<decltype(vec3(1, 2, 3) * isq::velocity[m / s] *
                                    (2. * isq::duration[s]))::rep,
                           vec3>);

// A vector quantity is NOT a representation type: a quantity can never be
// nested as another quantity's representation (`value_type_t<quantity>` is
// the quantity itself, which `disable_representation` rejects).
static_assert(!detail::VectorRepresentation<decltype(vec3(0, 0, 0) *
                                                     isq::velocity[m / s])>);
static_assert(!RepresentationOf<decltype(vec3(0, 0, 0) * isq::velocity[m / s]),
                                quantity_tensor_order::vector>);
// `magnitude()` of a vector quantity is a scalar quantity in the same unit.
static_assert(
    QuantityOf<decltype(magnitude(vec3(3, 4, 0) * isq::velocity[m / s])),
               isq::speed>);

//! @test Verifies typed linear algebra's mp-units compatibility plugin.
[[maybe_unused]] const auto test{[] {
  using vector3d =
      typed_column_vector<Eigen::Vector<double, 3>, double, double, double>;

  quantity velocity = isq::velocity(vector3d{30, 40, 0} * km / h);
  // From there the framework treats the vector exactly like any other
  // representation:
  // * Unit conversion scales every component:
  quantity velocity_mps = velocity.in(m / s);
  // *vector quantity × scalar quantity:
  quantity displacement = velocity * (2 * h);
  // * Euclidean magnitude → scalar quantity:
  quantity speed = magnitude(velocity);

  assert(std::format("{}", velocity) == "[[30], [40], [0]] km/h");
  assert(std::format("{}", velocity_mps) ==
         "[[8.333333333333334], [11.11111111111111], [0]] m/s");
  assert(std::format("{}", displacement) == "[[60], [80], [0]] km");
  assert(std::format("{}", speed) == "50 km/h");

  // Vector of quantities:
  // using vector3v = typed_column_vector<vector3d, quantity<isq::velocity[m /
  // s]>, quantity<isq::velocity[m / s]>, quantity<isq::velocity[m / s]>>;

  // Conversion from a vector as a quantity.
  // vector3v velocities{velocity};
  // NOT COMPILING BECAUSE MISSING CONSTRUCTOR
  // CONSTRUCTOR CAN USE get<N>() API ON QUANTITY?
  // MP_UNITS IMPLEMENTS GET ON QUANTITY AS A FRIEND IMPLEMENTATION OF
  // QUANTITY
  // WHY IS THIS ONE NOT WORKING?
  auto v0 = get<0>(velocity);

  // error: no matching constructor for initialization of
  // vector3v
  // aka

  // typed_matrix<
  //  fcarouge::typed_matrix<
  //    Eigen::Matrix<
  //      double, 3, 1, 0, 3, 1>,
  //    std::tuple<
  //      double, double, double>,
  //    std::tuple<
  //      std::identity>>,
  //  tuple<
  //    mp_units::quantity<
  //      reference<
  //        velocity, derived_unit<
  //          metre, per<second>>>{}, double>,
  //    mp_units::quantity<
  //      reference<
  //        velocity, derived_unit<
  //          metre, per<second>>>{}, double>,
  //    mp_units::quantity<
  //      reference<
  //        velocity, derived_unit<
  //          metre, per<second>>>{}, double>>,
  //  tuple<std::identity>>')

  // from

  // quantity<
  //  mp_units::reference<
  //    mp_units::isq::velocity,
  //    mp_units::derived_unit<
  //      mp_units::si::kilo_<
  //        mp_units::si::metre>,
  //      mp_units::per<
  //        mp_units::non_si::hour>>>{},
  //  typed_matrix<
  //    Eigen::Matrix<
  //      double, 3, 1, 0, 3, 1>,
  //    std::tuple<
  //      double, double, double>,
  //    std::tuple<
  //      std::identity>>>')

  // assert(std::format("{}", velocity) == "[[30 km/h], [40 km/h], [0 km/h]]");

  return 0;
}()};
} // namespace
} // namespace fcarouge::test

// ADD SAMPLE FOR NAMED AXIS? INDEX ACCESS OF QUANTITIES
