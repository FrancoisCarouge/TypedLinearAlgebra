/* Typed Linear Algebra
Version 0.2.0
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

//! @file
//! @brief Unit safe linear algebra with mp-units and Eigen.
//!
//! @details Demonstrate a variety of linear algebra operations with mp-units
//! and Eigen. This library composes Eigen as the matrix' linear algebra backend
//! with index typed as mp-units types. This sample explicitly uses double
//! precision floating point numbers. This sample uses Eigen linear algebra as
//! the linear algebra backend. This sample uses mp-units types for the strongly
//! typed units.

#include "fcarouge/linalg.hpp"
#include "fcarouge/typed_linear_algebra.hpp"

#include <cstddef>
#include <format>
#include <linalg>
#include <mdspan>
#include <print>
#include <tuple>
#include <type_traits>
#include <vector>

namespace fcarouge::sample {
namespace {
// Set up heterogenously unit typed linear algebra types.
using representation = double;

template <auto QuantityReference>
using quantity = mp_units::quantity<QuantityReference, representation>;

template <mp_units::Reference auto QuantityReference>
using quantity_point =
    mp_units::quantity_point<QuantityReference,
                             mp_units::default_point_origin(QuantityReference),
                             representation>;

template <typename RowIndexes, typename ColumnIndexes>
using matrix = typed_matrix<
    std::mdspan<representation,
                std::extents<std::size_t, std::tuple_size_v<RowIndexes>,
                             std::tuple_size_v<ColumnIndexes>>>,
    RowIndexes, ColumnIndexes>;

template <typename... Types>
using column_vector = typed_column_vector<
    std::mdspan<representation, std::extents<std::size_t, sizeof...(Types), 1>>,
    Types...>;

template <typename... Types>
using row_vector = typed_row_vector<
    std::mdspan<representation, std::extents<std::size_t, 1, sizeof...(Types)>>,
    Types...>;

template <std::size_t Rows>
using column_extents = std::extents<std::size_t, Rows, 1>;

template <std::size_t Columns>
using row_extents = std::extents<std::size_t, 1, Columns>;

template <typename Extents>
constexpr std::size_t extents_size{[] {
  std::size_t size{1};
  Extents extents{};
  for (std::size_t i{0}; i < Extents::rank(); ++i) {
    size *= extents.extent(i);
  }
  return size;
}()};

// Expose a few mp-units types and unit symbols.
using mp_units::one;
using mp_units::si::unit_symbols::A;
using mp_units::si::unit_symbols::m;
using mp_units::si::unit_symbols::m2;
using mp_units::si::unit_symbols::mol;
using mp_units::si::unit_symbols::s;
using mp_units::si::unit_symbols::s2;
using mp_units::si::unit_symbols::s3;

// Shorten some mp-units quantities.
using position = quantity<mp_units::isq::length[m]>;
using velocity = quantity<mp_units::isq::velocity[m / s]>;
using acceleration = quantity<mp_units::isq::acceleration[m / s2]>;

// Set up a heterogenous column vector type.
using state = column_vector<position, velocity, acceleration>;

// TYPE CARTESIAN PRODUCT
// 1. Crosses a single tuple of accumulated types `Ts...` with a list of new
// types `Us...`. Example: tuple<A, B> x tuple<C, D> -> tuple<tuple<A,B,C>,
// tuple<A,B,D>>
template <typename Tuple, typename Types> struct cross_one_tuple;

template <typename... Ts, typename... Us>
struct cross_one_tuple<std::tuple<Ts...>, std::tuple<Us...>> {
  using type = std::tuple<std::tuple<Ts..., Us>...>;
};

// 2. Crosses a list of tuples with a list of types.
// We use decltype(std::tuple_cat(...)) to cleanly and efficiently flatten the
// nested tuples.
template <typename Tuples, typename Types> struct cross_one;

template <typename... Tuples, typename Types>
struct cross_one<std::tuple<Tuples...>, Types> {
  using type = decltype(std::tuple_cat(
      std::declval<typename cross_one_tuple<Tuples, Types>::type>()...));
};

// 3. Recursively applies `cross_one` over all provided lists.
template <typename Tuples, typename... Lists> struct cross_all {
  using type = Tuples; // Base case: no more lists to cross
};

template <typename Tuples, typename L1, typename... Ls>
struct cross_all<Tuples, L1, Ls...> {
  using type =
      typename cross_all<typename cross_one<Tuples, L1>::type, Ls...>::type;
};

// 4. Applies a metafunction `F` to the expanded elements of a tuple.
template <template <typename...> class F, typename Tuple>
struct apply_f_to_tuple;

template <template <typename...> class F, typename... Ts>
struct apply_f_to_tuple<F, std::tuple<Ts...>> {
  using type = F<Ts...>;
};

// 5. Maps the metafunction `F` over a list of tuples.
template <template <typename...> class F, typename TupleOfTuples>
struct map_apply;

template <template <typename...> class F, typename... Tuples>
struct map_apply<F, std::tuple<Tuples...>> {
  using type = std::tuple<typename apply_f_to_tuple<F, Tuples>::type...>;
};

// The main algorithm
// Computes the Cartesian product of the provided lists and applies F to each
// combination.
template <template <typename...> class F, typename... Lists>
using mp_product = typename map_apply<
    F, typename cross_all<std::tuple<std::tuple<>>, Lists...>::type>::type;

template <typename RowIndexes, typename ColumnIndexes>
using tuple = mp_product<fcarouge::tla::product, RowIndexes, ColumnIndexes>;

template <typename RowIndexes, typename ColumnIndexes>
using typed_matrix = fcarouge::typed_matrix<tuple<RowIndexes, ColumnIndexes>,
                                            RowIndexes, ColumnIndexes>;

//! @brief Strongly typed linear algebra samples.
[[maybe_unused]] const auto sample{[] {
  [[maybe_unused]] typed_matrix<std::tuple<position, velocity>,
                                std::tuple<position, velocity>> a;

  a.at<0, 0>() = 33. * m2;

  // WIP...
  assert((33. * m2 == a.at<0, 0>()));
  // operations...

  return 0;
}()};
} // namespace
} // namespace fcarouge::sample
