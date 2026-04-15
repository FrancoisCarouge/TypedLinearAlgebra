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

#include "fcarouge/typed_linear_algebra.hpp"

#include <nanobench.h>

#include <cstddef>
#include <format>
#include <fstream>
#include <linalg>
#include <mdspan>
#include <random>
#include <string>
#include <vector>

namespace fcarouge::benchmark {
namespace {
template <auto Size>
const std::string csv{std::format(
    "{{{{#result}}}}| {{{{title}}}} | {:5d}x{:<5d} | {{{{median(elapsed)}}}} | "
    "{{{{medianAbsolutePercentError(elapsed)}}}} |{{{{/result}}}}\n",
    Size, Size)};

template <typename Tuple,
          typename Indices =
              std::make_index_sequence<std::tuple_size<Tuple>::value>>
struct lookup;

template <typename Tuple, std::size_t... Indices>
struct lookup<Tuple, std::index_sequence<Indices...>> {
  using result = typename std::tuple_element<0, Tuple>::type &;
  using pointer = result (*)(Tuple &) noexcept;

  static constexpr pointer table[std::tuple_size<Tuple>::value] = {
      &std::get<Indices>...};
};

template <typename Tuple, std::size_t... Indices>
constexpr
    typename lookup<Tuple, std::index_sequence<Indices...>>::pointer lookup<
        Tuple,
        std::index_sequence<Indices...>>::table[std::tuple_size<Tuple>::value];

template <typename Tuple>
constexpr typename std::tuple_element<
    0, typename std::remove_reference<Tuple>::type>::type &
get(std::size_t index, Tuple &&tuple) {
  using tuple_type = typename std::remove_reference<Tuple>::type;

  if (index >= std::tuple_size<tuple_type>::value) {
    throw std::runtime_error("Out of range");
  }

  return lookup<tuple_type>::table[index](tuple);
}

template <typename Tuple, typename Type> struct accessor {
  using element_type = Type;
  using reference = Type &;
  using data_handle_type = Tuple *;
  using offset_policy = accessor;

  static inline constexpr auto size{std::tuple_size_v<Tuple>};

  [[nodiscard]] static constexpr reference access(data_handle_type data,
                                                  std::size_t index) {
    return get(index, *data);
  }

  [[nodiscard]] static constexpr data_handle_type
  offset(data_handle_type data, [[maybe_unused]] std::size_t index) noexcept {
    return data;
  }
};

//! @benchmark `std::mdspan` square matrix-matrix product.
template <auto Size> void bench() {
  using tuple =
      typed_linear_algebra_internal::tuple_n_type<double, Size * Size>;
  using mdspan = std::mdspan<double, std::extents<std::size_t, Size, Size>,
                             Kokkos::layout_right, accessor<tuple, double>>;
  using matrix =
      typed_matrix<mdspan,
                   typed_linear_algebra_internal::tuple_n_type<double, Size>,
                   typed_linear_algebra_internal::tuple_n_type<double, Size>>;

  tuple storage_a;
  tuple storage_b;
  tuple storage_r;
  mdspan mdspan_a{&storage_a};
  mdspan mdspan_b{&storage_b};
  mdspan mdspan_r{&storage_r};
  matrix a{mdspan_a};
  matrix b{mdspan_b};
  matrix r{mdspan_r};
  std::random_device device;
  std::mt19937 generator{device()};
  std::uniform_real_distribution<> distribution{0., 1.};

  for (std::size_t i{0}; i < Size; ++i) {
    for (std::size_t j{0}; j < Size; ++j) {
      a(i, j) = distribution(generator);
      b(i, j) = distribution(generator);
    }
  }

  std::ofstream results{"results.txt", std::ios::app};
  ankerl::nanobench::Bench()
      .output(nullptr)
      .title("typed matrix from std::mdspan from std::tuple")
      .run([&]() {
        matrix_product(a, b, r);
        ankerl::nanobench::doNotOptimizeAway(r);
      })
      .render(csv<Size>.c_str(), results);
}
} // namespace
} // namespace fcarouge::benchmark

int main() { fcarouge::benchmark::bench<${SIZE}>(); }
