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

#include "fcarouge/linalg.hpp"

#include <cassert>
#include <cstddef>
#include <format>
#include <iostream>
#include <mdspan>
#include <print>
#include <tuple>
#include <type_traits>
#include <utility>

namespace mp {

namespace detail {

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

} // namespace detail

// The main algorithm
// Computes the Cartesian product of the provided lists and applies F to each
// combination.
template <template <typename...> class F, typename... Lists>
using mp_product = typename detail::map_apply<
    F,
    typename detail::cross_all<std::tuple<std::tuple<>>, Lists...>::type>::type;

} // namespace mp

// A sample metafunction to apply to our combinations
template <typename T1, typename T2, typename T3> struct MyType {};

// Type lists
using Modifiers1 = std::tuple<const int, volatile int>;
using Modifiers2 = std::tuple<char *, char &>;
using BaseTypes = std::tuple<float, double>;

int main() {
  // Generate the Cartesian product: Modifiers1 x Modifiers2 x BaseTypes
  using ProductList = mp::mp_product<MyType, Modifiers1, Modifiers2, BaseTypes>;

  // The expected result of the Cartesian product
  using ExpectedList = std::tuple<
      MyType<const int, char *, float>, MyType<const int, char *, double>,
      MyType<const int, char &, float>, MyType<const int, char &, double>,
      MyType<volatile int, char *, float>, MyType<volatile int, char *, double>,
      MyType<volatile int, char &, float>,
      MyType<volatile int, char &, double>>;

  // Verify at compile time
  static_assert(std::is_same_v<ProductList, ExpectedList>,
                "The Cartesian product did not match the expected type list!");

  std::cout << "Cartesian product successfully generated and verified at "
               "compile time.\n";
  return 0;
}

namespace fcarouge::test {

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
class tt_matrix {
public:
  using matrix = Matrix;
  using row_indexes = RowIndexes;
  using column_indexes = ColumnIndexes;
  using storage_t =
      mp::mp_product<fcarouge::tla::product, std::tuple<position, velocity>,
                     std::tuple<position, velocity>>;

  static inline constexpr auto rows{std::tuple_size_v<row_indexes>};
  static inline constexpr auto columns{std::tuple_size_v<column_indexes>};
  static inline constexpr auto rank{tla::rank<rows, columns>};

  // constexpr explicit tt_matrix(const Matrix &other);

  // Benefit: this is type safe now.
  constexpr explicit tt_matrix(const storage_t &other) : storage{other} {}

  template <auto... Indexes>
  [[nodiscard]] constexpr decltype(auto) at(this auto &&self)
    requires(sizeof...(Indexes) >= rank)
  {
    return std::get<std::get<0>(std::tuple{Indexes...}) * columns +
                    std::get<1>(std::tuple{Indexes...})>(self.storage);
  }

private:
  storage_t storage;
};

// Convert tuple to variant
template <typename T> struct tuple_to_variant;

template <typename... Ts> struct tuple_to_variant<std::tuple<Ts...>> {
  using type = std::variant<Ts...>;
};

// Deduplicate tuple types for variant well formed
// 1. Helper to check if a type T is in a list of types Us
template <typename T, typename... Us>
struct is_contained : std::disjunction<std::is_same<T, Us>...> {};

// 2. Base case: Empty tuple
template <typename... Ts> struct unique_tuple_impl {
  using type = std::tuple<>;
};

// 3. Recursive case: Process first type T, then remaining Ts
template <typename T, typename... Ts> struct unique_tuple_impl<T, Ts...> {
  using rest_unique = typename unique_tuple_impl<Ts...>::type;

  // If T is already in the unique list of the 'rest', don't add it
  using type = std::conditional_t<is_contained<T, Ts...>::value, rest_unique,
                                  decltype(std::tuple_cat(std::tuple<T>{},
                                                          rest_unique{}))>;
};

// 4. Alias for easier use
template <typename T> struct unique_tuple;

template <typename... Ts> struct unique_tuple<std::tuple<Ts...>> {
  using type = typename unique_tuple_impl<Ts...>::type;
};

template <typename T> using unique_tuple_t = typename unique_tuple<T>::type;

// mdpsan accessor
template <class T> struct accessor {
  using element_type = T;
  using reference = double &;
  using data_handle_type = T *;
  using offset_policy = accessor;

  constexpr accessor() = default;

  constexpr reference access(data_handle_type p, std::size_t i) const noexcept {
    return std::visit(
        [](auto &q) -> reference { return q.numerical_value_ref_in(q.unit); },
        p[i]);
  }

  // constexpr element_type &typed_access(data_handle_type p,
  //                                      std::size_t i) noexcept {
  //   return p[i];
  // }

  // Required: offset
  constexpr data_handle_type offset(data_handle_type p,
                                    std::size_t i) const noexcept {
    return p + i;
  }
};

// Example
using literals::operator""_i;
using representation = double;

template <auto QuantityReference>
using quantity = mp_units::quantity<QuantityReference, representation>;

namespace {
//! @test Verifies the addition operator with non-trivial types.
[[maybe_unused]] const auto test{[] {
  using position = quantity<mp_units::isq::length[m]>;
  using velocity = quantity<mp_units::isq::velocity[m / s]>;

  using all_types =
      mp::mp_product<fcarouge::tla::product, std::tuple<position, velocity>,
                     std::tuple<position, velocity>>;

  // all_types storage_a{1. * m2, 2. * m2 / s, 11. * m2 / s, 22. * m2 / s2};

  using unique_types = unique_tuple_t<all_types>;

  using common_type = typename tuple_to_variant<unique_types>::type;

  // accessor<common_type> a;

  using span = std::mdspan<common_type, std::extents<std::size_t, 2, 2>,
                           Kokkos::layout_right, accessor<common_type>>;

  common_type storage_a[]{1. * m2, 2. * m2 / s, 11. * m2 / s, 22. * m2 / s2};
  span span_a{&storage_a[0], std::extents<std::size_t, 2, 2>()};

  common_type storage_b[]{2. * m2, 2. * m2 / s, 11. * m2 / s, 22. * m2 / s2};
  span span_b{&storage_b[0], std::extents<std::size_t, 2, 2>()};

  common_type storage_r[]{0. * m2, 0. * m2 / s, 0. * m2 / s, 0. * m2 / s2};
  span span_r{&storage_r[0], std::extents<std::size_t, 2, 2>()};

  typed_matrix<decltype(span_a), std::tuple<position, velocity>,
               std::tuple<position, velocity>>
      a{span_a};
  typed_matrix<decltype(span_b), std::tuple<position, velocity>,
               std::tuple<position, velocity>>
      b{span_b};
  typed_matrix<decltype(span_r), std::tuple<position, velocity>,
               std::tuple<position, velocity>>
      r{span_r};

  r.at<0, 0>() = 33. * m2;
  assert((33. * m2 == r.at<0, 0>()));
  std::println("{}", r.at<0, 0>());

  add(a, b, r);

  std::println("{}", r.at<0, 0>());
  assert((3. * m2 == r.at<0, 0>()));

  // conlusion:
  // Possible, not interesting for size, and dispatch costs.


  // r.at<0, 0>() = 33. * m2;
  // std::println("{}", r.at<0, 0>());
  // assert((33. * m2 == r.at<0, 0>()));

  // assert((4. * m2 / s == r.at<0, 1>()));
  // assert((22. * m2 / s == r.at<1, 0>()));
  // assert((44. * m2 / s2 == r.at<1, 1>()));

  /////////////////////////////////////////////////////////////

  // std::mdspan<double, std::extents<std::size_t, 2, 2>> span_a{
  //     reinterpret_cast<double *>(&storage_a),
  //     std::extents<std::size_t, 2, 2>()};

  //////////////////////////////////////////////////////////

  // using unique_types = unique_tuple_t<all_types>;

  // using common_type = typename tuple_to_variant<unique_types>::type;

  // common_type storage_a[]{1. * m2, 2. * m2 / s, 11. * m2 / s, 22. * m2 / s2};
  // common_type storage_b[]{1. * m2, 2. * m2 / s, 11. * m2 / s, 22. * m2 / s2};
  // common_type storage_r[]{1. * m2, 2. * m2 / s, 11. * m2 / s, 22. * m2 / s2};

  // std::mdspan span_a{&storage_a[0], std::extents<std::size_t, 2, 2>{}};
  // std::mdspan span_b{&storage_b[0], std::extents<std::size_t, 2, 2>{}};
  // std::mdspan span_r{&storage_r[0], std::extents<std::size_t, 2, 2>{}};

  // typed_matrix<decltype(span_a), std::tuple<position, velocity>,
  //              std::tuple<position, velocity>>
  //     a{span_a};
  // typed_matrix<decltype(span_b), std::tuple<position, velocity>,
  //              std::tuple<position, velocity>>
  //     b{span_b};
  // typed_matrix<decltype(span_r), std::tuple<position, velocity>,
  //              std::tuple<position, velocity>>
  //     r{span_r};

  /////////////////////////////////////////////////////////////

  // using algebra = std::mdspan<double, std::extents<std::size_t, 2, 2>>;

  // using mat = tt_matrix<algebra, std::tuple<position, velocity>,
  //                       std::tuple<position, velocity>>;

  // mat::storage_t store{1. * m2, 2. * m2 / s, 11. * m2 / s, 22. * m2 / s2};

  // mat r{store};

  // assert((1. * m2 == r.at<0, 0>()));
  // assert((2. * m2 / s == r.at<0, 1>()));
  // assert((11. * m2 / s == r.at<1, 0>()));
  // assert((22. * m2 / s2 == r.at<1, 1>()));

  // r.at<0, 0>() = 3. * m2;
  // r.at<0, 1>() = 4. * m2 / s;
  // r.at<1, 0>() = 5. * m2 / s;
  // r.at<1, 1>() = 6. * m2 / s2;

  // assert((3. * m2 == r.at<0, 0>()));
  // assert((4. * m2 / s == r.at<0, 1>()));
  // assert((5. * m2 / s == r.at<1, 0>()));
  // assert((6. * m2 / s2 == r.at<1, 1>()));

  /////////////////////////////////////////////////////////

  // double storage_a[]{0., 0.};
  // double storage_b[]{0., 0.};
  // double storage_r[]{0., 0.};

  // std::mdspan span_a{&storage_a[0], std::extents<std::size_t, 1, 2>{}};
  // std::mdspan span_b{&storage_b[0], std::extents<std::size_t, 1, 2>{}};
  // std::mdspan span_r{&storage_r[0], std::extents<std::size_t, 1, 2>{}};

  // row_vector<representation, position, velocity> a{span_a};
  // row_vector<representation, position, velocity> b{span_b};
  // row_vector<representation, position, velocity> r{span_r};

  // a[0_i] = 1. * m;
  // a[1_i] = 2. * m / s;
  // b[0_i] = 3. * m;
  // b[1_i] = 4. * m / s;

  // add(a, b, r);

  // assert((4. * m == r.at<0, 0>()));
  // assert((4. * m == r(0_i, 0_i)));
  // assert((4. * m == r[0_i, 0_i]));
  // assert((4. * m == r.at<0_i, 0_i>()));
  // assert((4. * m == r.at<0>()));
  // assert((4. * m == r(0_i)));
  // assert((4. * m == r[0_i]));
  // assert((4. * m == r.at<0_i>()));

  // assert((6. * m / s == r.at<0, 1>()));
  // assert((6. * m / s == r(0_i, 1_i)));
  // assert((6. * m / s == r[0_i, 1_i]));
  // assert((6. * m / s == r.at<0_i, 1_i>()));
  // assert((6. * m / s == r.at<1>()));
  // assert((6. * m / s == r(1_i)));
  // assert((6. * m / s == r[1_i]));
  // assert((6. * m / s == r.at<1_i>()));

  return 0;
}()};
} // namespace
} // namespace fcarouge::test
