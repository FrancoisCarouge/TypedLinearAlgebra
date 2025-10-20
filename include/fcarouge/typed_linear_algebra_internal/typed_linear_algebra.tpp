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

#ifndef FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_TYPED_LINEAR_ALGEBRA_TPP
#define FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_TYPED_LINEAR_ALGEBRA_TPP

namespace fcarouge {
namespace tla = typed_linear_algebra_internal;

//! @todo Verify types and storage (?) compatibility.
template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::typed_matrix(
    const is_typed_matrix auto &other)
    : storage{other.data()} {}

//! @todo Verify types and storage (?) compatibility.
template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes> &
typed_matrix<Matrix, RowIndexes, ColumnIndexes>::operator=(
    const is_typed_matrix auto &other) {
  storage = other.data();
  return *this;
}

//! @todo Verify types and storage (?) compatibility.
template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::typed_matrix(
    const is_typed_matrix auto &&other)
    : storage{std::forward<decltype(other)>(other).data()} {}

//! @todo Verify types and storage (?) compatibility.
template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes> &
typed_matrix<Matrix, RowIndexes, ColumnIndexes>::operator=(
    const is_typed_matrix auto &&other) {
  storage = std::forward<decltype(other)>(other).data();
  return *this;
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::typed_matrix(
    const Matrix &other)
    : storage{other} {}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::typed_matrix(
    const element<0, 0> (&elements)[typed_matrix::rows * typed_matrix::columns])
  requires is_uniform_typed_matrix<typed_matrix> and
           is_one_dimension_typed_matrix<typed_matrix>
    : storage{elements} {}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::typed_matrix(
    const auto &value)
  requires is_singleton_typed_matrix<typed_matrix>
{
  if constexpr (requires { value(std::size_t{0}, std::size_t{0}); }) {
    storage(std::size_t{0}, std::size_t{0}) =
        underlying{value(std::size_t{0}, std::size_t{0})};
  } else {
    using type = std::remove_cvref_t<decltype(value)>;
    storage(std::size_t{0}, std::size_t{0}) = cast<underlying, type>(value);
  }
}

//! @todo Verify the list sizes at runtime? Deprecate?
//! @todo Verify `Type` is `element<0, 0>`-compatible-safe?
template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
template <typename Type>
constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::typed_matrix(
    std::initializer_list<std::initializer_list<Type>> row_list)
  requires is_uniform_typed_matrix<typed_matrix>
{
  for (std::size_t i{0}; const auto &row : row_list) {
    for (std::size_t j{0}; const auto &value : row) {
      storage(std::size_t{i}, std::size_t{j}) = cast<underlying, Type>(value);
      ++j;
    }
    ++i;
  }
}

//! @todo Verify if the types are the same, or assignable, for nicer error?
//! @todo Rewrite with a fold expression over the pack?
template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::typed_matrix(
    const auto &first_value, const auto &second_value, const auto &...values)
  requires is_one_dimension_typed_matrix<typed_matrix>
{
  //! @todo Move the assert as a require clause when the compilers support it.
  static_assert(columns * rows == 2 + sizeof...(values),
                "The count of parameters must match the size of the vector.");
  std::tuple value_pack{first_value, second_value, values...};
  tla::for_constexpr<0, typed_matrix::rows * typed_matrix::columns, 1>(
      [this, &value_pack](auto position) {
        auto value{std::get<position>(value_pack)};
        using type = std::remove_cvref_t<decltype(value)>;
        storage(std::size_t{position}) = cast<underlying, type>(value);
      });
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
[[nodiscard]] constexpr typed_matrix<
    Matrix, RowIndexes, ColumnIndexes>::operator element<0, 0>(this auto &&self)
  requires is_singleton_typed_matrix<typed_matrix>
{
  return cast<element<0, 0>, underlying>(
      self.storage(std::size_t{0}, std::size_t{0}));
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
[[nodiscard]] constexpr auto &&
typed_matrix<Matrix, RowIndexes, ColumnIndexes>::operator[](this auto &&self,
                                                            std::size_t index)
  requires(is_uniform_typed_matrix<typed_matrix> and
           is_one_dimension_typed_matrix<typed_matrix>)
{
  return std::forward<decltype(self)>(self).storage(std::size_t{index});
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
[[nodiscard]] constexpr auto &&
typed_matrix<Matrix, RowIndexes, ColumnIndexes>::operator[](this auto &&self,
                                                            std::size_t row,
                                                            std::size_t column)
  requires is_uniform_typed_matrix<typed_matrix>
{
  return std::forward<decltype(self)>(self).storage(std::size_t{row},
                                                    std::size_t{column});
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
[[nodiscard]] constexpr auto &&
typed_matrix<Matrix, RowIndexes, ColumnIndexes>::operator()(this auto &&self,
                                                            std::size_t index)
  requires is_uniform_typed_matrix<typed_matrix> and
           is_one_dimension_typed_matrix<typed_matrix>
{
  return std::forward<decltype(self)>(self).storage(std::size_t{index});
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
[[nodiscard]] constexpr auto &&
typed_matrix<Matrix, RowIndexes, ColumnIndexes>::operator()(this auto &&self,
                                                            std::size_t row,
                                                            std::size_t column)
  requires is_uniform_typed_matrix<typed_matrix>
{
  return std::forward<decltype(self)>(self).storage(std::size_t{row},
                                                    std::size_t{column});
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
template <std::size_t Row, std::size_t Column>
[[nodiscard]] constexpr auto
typed_matrix<Matrix, RowIndexes, ColumnIndexes>::at() -> element<Row, Column> &
  requires(Row < rows) and (Column < columns)
{
  return cast<element<Row, Column> &, underlying &>(
      storage(std::size_t{Row}, std::size_t{Column}));
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
template <std::size_t Row, std::size_t Column>
[[nodiscard]] constexpr auto
typed_matrix<Matrix, RowIndexes, ColumnIndexes>::at() const
    -> element<Row, Column>
  requires(Row < rows) and (Column < columns)
{
  return cast<element<Row, Column>, underlying>(
      storage(std::size_t{Row}, std::size_t{Column}));
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
template <std::size_t Index>
[[nodiscard]] constexpr auto
typed_matrix<Matrix, RowIndexes, ColumnIndexes>::at() -> element<Index, 0> &
  requires is_column_typed_matrix<typed_matrix> and (Index < rows)
{
  return cast<element<Index, 0> &, underlying &>(storage(std::size_t{Index}));
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
template <std::size_t Index>
[[nodiscard]] constexpr auto
typed_matrix<Matrix, RowIndexes, ColumnIndexes>::at() const -> element<Index, 0>
  requires is_column_typed_matrix<typed_matrix> and (Index < rows)
{
  return cast<element<Index, 0>, underlying>(storage(std::size_t{Index}));
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
[[nodiscard]] constexpr auto &&
typed_matrix<Matrix, RowIndexes, ColumnIndexes>::data(this auto &&self) {
  return std::forward<decltype(self)>(self).storage;
}
} // namespace fcarouge

#endif // FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_TYPED_LINEAR_ALGEBRA_TPP
