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

#ifndef FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_TUPLE_TPP
#define FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_TUPLE_TPP

//! @brief Tuple size specialization in support of structured bindings.
//!
//! @details Specialization for one-dimension matrices. The tuple size is the
//! number of elements in the matrix.
template <fcarouge::rank_typed_matrix<1> Type>
struct std::tuple_size<Type>
    : std::integral_constant<std::size_t, Type::rows * Type::columns> {};

//! @brief Tuple element specialization in support of structured bindings.
//!
//! @details Specialization for one-dimension matrices. The tuple element is the
//! matrix element at the given index.
template <std::size_t Index, fcarouge::rank_typed_matrix<1> Type>
  requires(Index < std::tuple_size_v<Type>)
struct std::tuple_element<Index, Type> {
  using type = typename Type::template element<Index>;
};

//! @brief Tuple size specialization in support of structured bindings.
//!
//! @details Specialization for zero-dimension matrices. The tuple size is one.
template <fcarouge::rank_typed_matrix<0> Type>
struct std::tuple_size<Type> : std::integral_constant<std::size_t, 1> {};

//! @brief Tuple element specialization in support of structured bindings.
//!
//! @details Specialization for zero-dimension matrices. The tuple element is
//! the single matrix element.
template <std::size_t Index, fcarouge::rank_typed_matrix<0> Type>
  requires(Index == 0)
struct std::tuple_element<Index, Type> {
  using type = typename Type::template element<>;
};

#endif // FCAROUGE_TYPED_LINEAR_ALGEBRA_INTERNAL_TUPLE_TPP
