/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Adrian Kummerlaender
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

#ifndef UTILITIES_META_H
#define UTILITIES_META_H

#include <type_traits>
#include <tuple>
#include <utility>

// Forward-definition of ADf type marker
struct AD;

namespace olb {

namespace utilities {

namespace meta {

/// Checks whether T can be used as a scalar arithmetic type
/**
 * Base checking for AD types to be replaced by e.g. concepts in C++20
 **/
template <typename T>
using is_arithmetic = typename std::integral_constant<bool,
     std::is_base_of<AD,T>::type::value
  || std::is_arithmetic<T>::type::value
>;

template <bool B, class T = void>
using disable_if_t = typename std::enable_if<!B,T>::type;

template <typename T, typename U = void>
using enable_if_arithmetic_t = std::enable_if_t<is_arithmetic<T>::type::value, U>;

/// Check whether a given type list contains WANTED
template <
  typename WANTED,
  typename HEAD = void, // Default argument in case the list is empty
  typename... TAIL
>
struct list_contains_item {
  using type = std::conditional_t<
    std::is_same<WANTED, HEAD>::value,
    std::true_type,
    typename list_contains_item<WANTED, TAIL...>::type
  >;
};

template <typename WANTED, typename HEAD>
struct list_contains_item<WANTED, HEAD> {
  using type = typename std::is_same<WANTED, HEAD>::type;
};


/// Get first type based on BASE contained in a given type list
/**
 * If no such list item exists, type is void.
 **/
template <
  typename BASE,
  typename HEAD = void, // Default argument in case the list is empty
  typename... TAIL
>
struct list_item_with_base {
  using type = std::conditional_t<
    std::is_base_of<BASE, HEAD>::value,
    HEAD,
    typename list_item_with_base<BASE, TAIL...>::type
  >;
};

template <typename BASE, typename HEAD>
struct list_item_with_base<BASE, HEAD> {
  using type = std::conditional_t<
    std::is_base_of<BASE, HEAD>::value,
    HEAD,
    void
  >;
};

/// Plain wrapper for type list of FIELDS
template <typename... FIELDS>
struct field_list {
  template <template<typename...> class COLLECTION>
  using decomposeInto = COLLECTION<FIELDS...>;
};

/// Cons a head FIELD onto an existing list
template <typename HEAD, typename FIELD_LIST>
struct field_list_cons {
  template <typename... FIELDS>
  using cons_helper = field_list<HEAD,FIELDS...>;

  using type = typename FIELD_LIST::template decomposeInto<cons_helper>;
};

/// Return type list of all FIELDS meeting COND
template <template <typename> class COND, typename FIELD=void, typename... FIELDS>
struct field_list_filter {
  using type = std::conditional_t<
    COND<FIELD>::value && !std::is_void<FIELD>::value,
    typename field_list_cons<FIELD, typename field_list_filter<COND, FIELDS...>::type>::type,
    typename field_list_filter<COND, FIELDS...>::type
  >;
};

/// Return either nil type list or type list containing FIELD depending on COND
template <template <typename> class COND, typename FIELD>
struct field_list_filter<COND,FIELD> {
  using type = std::conditional_t<
    COND<FIELD>::value && !std::is_void<FIELD>::value,
    field_list<FIELD>,
    field_list<>
  >;
};

template <template <typename> class COND, typename... FIELDS>
using field_list_filter_t = typename field_list_filter<COND,FIELDS...>::type;

/// Return void for any given type
/**
 * Useful for _swallowing_ variadic arguments and replacing them with arbitrary values.
 * See e.g. `src/core/fieldArrayD.h`. To be replaced by `std::void_t` when we start to
 * use C++17.
 **/
template <typename...>
using void_t = void;

/// Identity type to pass non-constructible types as value
/**
 * Aid for using fields in generic lambdas
 **/
template <typename TYPE>
struct id {
  using type = TYPE;
};

/// Function equivalent of void_t, swallows any argument
template <typename... ARGS>
void swallow(ARGS&&...) { }

/// Apply F to each element of TUPLE listed in INDICES
template <typename TUPLE, typename F, std::size_t... INDICES>
void tuple_for_each_index(TUPLE& tuple, F&& f, std::index_sequence<INDICES...>) {
  swallow((f(std::get<INDICES>(tuple)), INDICES)...);
}

/// Apply F to each element of TUPLE
template <typename TUPLE, typename F>
void tuple_for_each(TUPLE& tuple, F&& f) {
  tuple_for_each_index(tuple, std::forward<F>(f), std::make_index_sequence<std::tuple_size<TUPLE>::value>{});
}

/// Return std::array<T,D> where T is initialized with a common value (helper)
template <typename T, unsigned D, typename U, std::size_t... INDICES>
std::array<T,D> make_array(U&& u, std::index_sequence<INDICES...>) {
  return std::array<T,D>{(swallow(INDICES), u)...};
}

/// Return std::array<T,D> where T is initialized with a common value
template <typename T, unsigned D, typename U>
std::array<T,D> make_array(U&& u) {
  return make_array<T,D,U>(std::forward<U>(u), std::make_index_sequence<D>{});
}

/// Return std::array<T,D> where T is initialized using a iDim-dependent function (helper)
template <typename T, unsigned D, typename F, std::size_t... INDICES>
std::array<T,D> make_array_f(F&& f, std::index_sequence<INDICES...>) {
  return std::array<T,D>{f(INDICES)...};
}

/// Return std::array<T,D> where T is initialized using a iDim-dependent function
template <typename T, unsigned D, typename F>
std::array<T,D> make_array_f(F&& f) {
  return make_array_f<T,D,F>(std::forward<F&&>(f), std::make_index_sequence<D>{});
}

}

}

}

#endif
