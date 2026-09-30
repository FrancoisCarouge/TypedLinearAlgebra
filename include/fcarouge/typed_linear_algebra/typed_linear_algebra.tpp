/* Typed Linear Algebra
Version 0.4.0
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

#ifndef FCAROUGE_TYPED_LINEAR_ALGEBRA_TYPED_LINEAR_ALGEBRA_TPP
#define FCAROUGE_TYPED_LINEAR_ALGEBRA_TYPED_LINEAR_ALGEBRA_TPP

namespace fcarouge {
namespace tla = typed_linear_algebra_internal;

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::typed_matrix()
  requires std::default_initializable<Matrix>
    : storage{} {
  if constexpr (requires { Matrix::Zero(); }) {
    storage = Matrix::Zero();
  }
}

//! @todo Verify types and storage (?) compatibility.
template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::typed_matrix(
    const same_as_typed_matrix auto &other)
    : storage{other.data()} {}

//! @todo Verify types and storage (?) compatibility.
template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
constexpr auto typed_matrix<Matrix, RowIndexes, ColumnIndexes>::operator=(
    const same_as_typed_matrix auto &other)
    -> typed_matrix<Matrix, RowIndexes, ColumnIndexes> & {
  storage = other.data();
  return *this;
}

//! @todo Verify types and storage (?) compatibility.
template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::typed_matrix(
    same_as_typed_matrix auto &&other)
    : storage{std::forward<decltype(other)>(other).data()} {}

//! @todo Verify types and storage (?) compatibility.
template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
constexpr auto typed_matrix<Matrix, RowIndexes, ColumnIndexes>::operator=(
    same_as_typed_matrix auto &&other)
    -> typed_matrix<Matrix, RowIndexes, ColumnIndexes> & {
  storage = std::forward<decltype(other)>(other).data();
  return *this;
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::typed_matrix(
    Matrix other)
    : storage{std::move(other)} {}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::typed_matrix(
    const element<> (&elements)[typed_matrix::rows * typed_matrix::columns])
  requires rank_typed_matrix<typed_matrix, 0>
{
  if constexpr (requires { storage = elements; }) {
    storage = elements;
  } else {
    tla::for_constexpr<typed_matrix::rows * typed_matrix::columns>(
        [this, &elements](auto position) {
          tla::write_element<columns>(
              storage, cast<underlying, element<>>(elements[position]),
              position);
        });
  }
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
constexpr auto typed_matrix<Matrix, RowIndexes, ColumnIndexes>::operator=(
    const element<> (&elements)[typed_matrix::rows * typed_matrix::columns])
    -> typed_matrix<Matrix, RowIndexes, ColumnIndexes> &
  requires rank_typed_matrix<typed_matrix, 0>
{
  using element_type = element<>;

  if constexpr (requires { storage = elements; }) {
    storage = elements;
  } else {
    tla::for_constexpr<typed_matrix::rows * typed_matrix::columns>(
        [this, &elements](auto position) {
          tla::write_element<columns>(
              storage, cast<underlying, element_type>(elements[position]),
              position);
        });
  }
  return *this;
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::typed_matrix(
    const element<0> (&elements)[typed_matrix::rows * typed_matrix::columns])
  requires rank_typed_matrix<typed_matrix, 1> and
           uniform_typed_matrix<typed_matrix>
{
  using element_type = element<0>;

  if constexpr (requires { storage = elements; }) {
    storage = elements;
  } else {
    tla::for_constexpr<typed_matrix::rows * typed_matrix::columns>(
        [this, &elements](auto position) {
          tla::write_element<columns>(
              storage, cast<underlying, element_type>(elements[position]),
              position);
        });
  }
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
constexpr auto typed_matrix<Matrix, RowIndexes, ColumnIndexes>::operator=(
    const element<0> (&elements)[typed_matrix::rows * typed_matrix::columns])
    -> typed_matrix<Matrix, RowIndexes, ColumnIndexes> &
  requires rank_typed_matrix<typed_matrix, 1> and
           uniform_typed_matrix<typed_matrix>
{
  using element_type = element<0>;

  if constexpr (requires { storage = elements; }) {
    storage = elements;
  } else {
    tla::for_constexpr<typed_matrix::rows * typed_matrix::columns>(
        [this, &elements](auto position) {
          tla::write_element<columns>(
              storage, cast<underlying, element_type>(elements[position]),
              position);
        });
  }
  return *this;
}

//! @todo How to handle all combinations of storage and value types?
template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::typed_matrix(
    const element<> &value)
  requires rank_typed_matrix<typed_matrix, 0>
{
  tla::write_element<columns>(storage, cast<underlying, element<>>(value));
}

//! @todo How to handle all combinations of storage and value types?
template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
constexpr auto typed_matrix<Matrix, RowIndexes, ColumnIndexes>::operator=(
    const element<> &value) -> typed_matrix<Matrix, RowIndexes, ColumnIndexes> &
  requires rank_typed_matrix<typed_matrix, 0>
{
  tla::write_element<columns>(storage, cast<underlying, element<>>(value));

  return *this;
}

//! @todo Verify the list sizes at runtime? Deprecate?
//! @todo Verify `Type` is `element<0, 0>`-compatible-safe?
template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
template <typename Type>
constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::typed_matrix(
    std::initializer_list<std::initializer_list<Type>> row_list)
  requires uniform_typed_matrix<typed_matrix>
{
  for (std::size_t i{0}; const auto &row : row_list) {
    for (std::size_t j{0}; const auto &value : row) {
      tla::write_element<columns>(storage, cast<underlying, Type>(value), i, j);
      ++j;
    }
    ++i;
  }
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::typed_matrix(
    const auto &first_value, const auto &second_value, const auto &...values)
  requires rank_typed_matrix<typed_matrix, 1>
{
  //! @todo Move the assert as a require clause when the compilers support it.
  static_assert(columns * rows == 2 + sizeof...(values),
                "The count of parameters must match the size of the vector.");
  std::tuple value_pack{first_value, second_value, values...};
  tla::for_constexpr<typed_matrix::rows * typed_matrix::columns>(
      [this, &value_pack](auto position) {
        auto value{std::get<position>(value_pack)};
        using type = std::remove_cvref_t<decltype(value)>;
        static_assert(
            std::is_assignable_v<element<std::size_t{position}> &, type>,
            "The parameter type is not compatible with the element type.");
        tla::write_element<columns>(storage, cast<underlying, type>(value),
                                    position);
      });
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::typed_matrix(
    const other_tuple_like_vector<typed_matrix> auto &value) {
  tla::for_constexpr<rows * columns>([&](auto i) {
    tla::write_element<columns>(storage,
                                cast<underlying, element<i>>(get<i>(value)), i);
  });
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
[[nodiscard]] constexpr typed_matrix<Matrix, RowIndexes, ColumnIndexes>::
operator element<>(this auto &&self)
  requires rank_typed_matrix<typed_matrix, 0>
{
  return self.at();
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
template <typename... Indexes>
[[nodiscard]] constexpr auto
typed_matrix<Matrix, RowIndexes, ColumnIndexes>::operator[](this auto &&self,
                                                            Indexes... indexes)
    -> decltype(auto)
  requires(sizeof...(Indexes) == rank) and
          ((index<Indexes> && ...) or uniform_typed_matrix<typed_matrix>)
{
  return self.operator()(indexes...);
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
template <typename... Indexes>
[[nodiscard]] constexpr auto
typed_matrix<Matrix, RowIndexes, ColumnIndexes>::operator()(this auto &&self,
                                                            Indexes... indexes)
    -> decltype(auto)
  requires(sizeof...(Indexes) == rank) and
          ((index<Indexes> && ...) or uniform_typed_matrix<typed_matrix>)
{
  if constexpr ((index<Indexes> && ...)) {
    return self.template at<indexes...>();
  } else {
    // A uniform typed matrix: every element is of the first element's type.
    using self_t = std::remove_reference_t<decltype(self)>;
    using storage_access =
        decltype(tla::read_element<columns>(self.storage, indexes...));

    // An unevaluated expression template has nothing to reference, so read it
    // by value even through a non-const `self`.
    static constexpr bool by_value{std::is_const_v<self_t> ||
                                   not std::is_reference_v<storage_access>};

    using qualified_underlying =
        std::conditional_t<by_value, underlying, underlying &>;
    using qualified_element =
        std::conditional_t<by_value, tla::element_at<typed_matrix, 0, 0>,
                           tla::element_at<typed_matrix, 0, 0> &>;

    return cast<qualified_element, qualified_underlying>(
        tla::read_element<columns>(self.storage, indexes...));
  }
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
template <auto... Indexes>
[[nodiscard]] constexpr auto
typed_matrix<Matrix, RowIndexes, ColumnIndexes>::at(this auto &&self)
    -> decltype(auto)
  requires(sizeof...(Indexes) == rank)
{
  using self_t = std::remove_reference_t<decltype(self)>;
  using storage_access =
      decltype(tla::read_element<columns>(self.storage, Indexes...));

  // An unevaluated expression template has nothing to reference, so read it
  // by value even through a non-const `self`.
  static constexpr bool by_value{std::is_const_v<self_t> ||
                                 not std::is_reference_v<storage_access>};

  using qualified_underlying =
      std::conditional_t<by_value, underlying, underlying &>;
  using qualified_element =
      std::conditional_t<by_value, element<Indexes...>, element<Indexes...> &>;

  return cast<qualified_element, qualified_underlying>(
      tla::read_element<columns>(self.storage, Indexes...));
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
template <typename Type>
[[nodiscard]] constexpr decltype(auto)
typed_matrix<Matrix, RowIndexes, ColumnIndexes>::at(this auto &&self)
  requires distinct_typed_matrix<typed_matrix>
{
  using tuple = tla::tuple_typed_matrix<typed_matrix>;

  constexpr std::size_t matches{tla::count_convertible_indexes<Type, tuple>()};

  static_assert(matches != 0,
                "No element type of the typed matrix is implicitly convertible "
                "to the requested type.");
  static_assert(matches <= 1,
                "The requested type is an implicit conversion target of more "
                "than one element type: the by-type lookup is ambiguous.");

  if constexpr (matches == 1) {
    constexpr std::size_t index{
        tla::find_first_convertible_index<Type, tuple>()};

    if constexpr (rank_typed_matrix<typed_matrix, 2>) {
      return self.template at<index / columns, index % columns>();
    } else if constexpr (rank_typed_matrix<typed_matrix, 1>) {
      return self.template at<index>();
    } else {
      return self.template at<>();
    }
  }
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
template <auto... Indexes>
constexpr void typed_matrix<Matrix, RowIndexes, ColumnIndexes>::at(
    this auto &&self, const element<Indexes...> &value)
  requires(sizeof...(Indexes) == rank)
{
  // Forward the write operation to the storage through the shared element
  // access support, which selects the access pattern from the linear algebra
  // backend's available API. The value is converted to the storage underlying
  // type here, where both the element and underlying types are known.
  tla::write_element<columns>(
      self.storage, cast<underlying, element<Indexes...>>(value), Indexes...);
}

template <typename Matrix, typename RowIndexes, typename ColumnIndexes>
[[nodiscard]] constexpr auto
typed_matrix<Matrix, RowIndexes, ColumnIndexes>::data(this auto &&self)
    -> decltype(auto) {
  return std::forward_like<decltype(self)>(self.storage);
}

template <typename RowIndexes, typename ColumnIndexes>
[[nodiscard]] constexpr auto make_typed_matrix(auto &&value) {
  using type = decltype(value);
  using matrix = std::remove_cvref_t<type>;

  return typed_matrix<matrix, RowIndexes, ColumnIndexes>{
      std::forward<type>(value)};
}

template <int Index>
auto get(rank_typed_matrix<1> auto &&value) -> decltype(auto) {
  return (value.template at<Index>());
}

template <int Index>
auto get(rank_typed_matrix<0> auto &&value) -> decltype(auto) {
  static_assert(Index == 0, "The index must be zero for a rank 0 matrix.");

  return (value.at());
}
} // namespace fcarouge

#endif // FCAROUGE_TYPED_LINEAR_ALGEBRA_TYPED_LINEAR_ALGEBRA_TPP
