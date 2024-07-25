/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Adrian Kummerlaender
 *                2021 Adrian Kummerlaender, Nicolas Hafen, Mathias J. Krause
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

#ifndef CORE_META_H
#define CORE_META_H

#include <type_traits>
#include <typeindex>
#include <tuple>
#include <set>
#include <utility>
#include <string>

// Forward-declaration of ADf type marker
struct AD;

// Define operator+ to mitigate buggy unqualified name lookup in fold expressions for Clang Versions < 12
// see https://bugs.llvm.org/show_bug.cgi?id=30738
#if defined __clang__ && __clang_major__ < 12

namespace std {

/// Concatenate two index sequences
template <size_t... Is, size_t... Js>
constexpr auto operator+(index_sequence<Is...>, index_sequence<Js...>)
{
  return index_sequence<Is..., Js...>();
}

}

#endif

namespace olb {

struct ExprBase { };
struct SimdBase { };

namespace meta {

// *INDENT-OFF*

template <typename... TYPES> struct list;

/// Returns true iff address ptr is aligned w.r.t. TYPE
template <typename TYPE>
bool is_aligned(const void* ptr) {
  auto addr = reinterpret_cast<std::uintptr_t>(ptr);
  return !(addr % alignof(TYPE));
}

/// Identity type to pass non-constructible types as value
/**
 * Aid for using fields in generic lambdas
 **/
template <typename TYPE>
struct id {
  using type = TYPE;

  operator std::type_index() const {
    return typeid(TYPE);
  }

  TYPE get() const {
    return TYPE();
  }
};

template <typename TYPE>
using id_t = typename id<TYPE>::type;

/// Returns distinct name on GCC, Clang and ICC but may return arbitrary garbage as per the standard
template <typename T>
std::string name() {
  return typeid(T).name();
}

/// Function equivalent of void_t, swallows any argument
template <typename... ARGS>
void swallow(ARGS&&...) { }

/// Checks whether T can be used as a scalar arithmetic type
/**
 * Base checking for AD types to be replaced by e.g. concepts in C++20
 **/
template <typename T>
using is_arithmetic = typename std::integral_constant<bool,
     std::is_base_of_v<AD,T>
  || std::is_base_of_v<SimdBase,T>
  || std::is_base_of_v<ExprBase,T>
  || std::is_arithmetic_v<T>
  // Pointers and enums are enabled for convenience
  // TODO: Add specific arithmetic-support checking to arithmetic Vector operators
  || std::is_pointer_v<T>
  || std::is_enum_v<T>
>;

template <typename T, typename U = void>
using enable_if_arithmetic_t = std::enable_if_t<is_arithmetic<T>::type::value, U>;

/// Returns true iff a given type list contains WANTED
template <typename WANTED, typename... TYPES>
constexpr bool contains() {
  return (std::is_same_v<WANTED, TYPES> || ... || false);
}

/// Returns true iff a given type list contains a type index
template <typename... TYPES>
bool contains(std::type_index field) {
  return ((field == typeid(TYPES)) || ... || false);
}

/// Get first type based on BASE contained in a given type list
/**
 * If no such list item exists, type is BASE.
 **/
//TODO: check, whether this can be generalized with "list_item_with_base"
template <
  typename BASE,
  typename HEAD = void, // Default argument in case the list is empty
  typename... TAIL
  >
struct list_item_with_base_default_base {
  using type = std::conditional_t<
    std::is_base_of<BASE, HEAD>::value,
    HEAD,
    typename list_item_with_base_default_base<BASE, TAIL...>::type
  >;
};

template <typename BASE, typename HEAD>
struct list_item_with_base_default_base<BASE, HEAD> {
  using type = std::conditional_t<
    std::is_base_of<BASE, HEAD>::value,
    HEAD,
    BASE
  >;
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
struct first_type_with_base {
  using type = std::conditional_t<
    std::is_base_of<BASE, HEAD>::value,
    HEAD,
    typename first_type_with_base<BASE, TAIL...>::type
  >;
};

template <typename BASE, typename HEAD>
struct first_type_with_base<BASE, HEAD> {
  using type = std::conditional_t<
    std::is_base_of<BASE, HEAD>::value,
    HEAD,
    void
  >;
};

template <typename BASE, typename... TYPES>
using first_type_with_base_t = typename first_type_with_base<BASE, TYPES...>::type;

/// Helper for computing indices in type lists
template <template<typename> typename COND, typename... TYPES>
struct index_of_first_matching;

template <template<typename> typename COND, typename HEAD, typename... TAIL>
struct index_of_first_matching<COND,HEAD,TAIL...> {
  static constexpr unsigned value = std::conditional_t<
    COND<HEAD>::value,
    std::integral_constant<unsigned, 0>,
    std::integral_constant<unsigned, 1 + index_of_first_matching<COND,TAIL...>::value>
  >::value;
};

template <template<typename> typename COND>
struct index_of_first_matching<COND> {
  static constexpr unsigned value = 0; // This way access into a list of this size fails (hacky)
};

/// Evaluates to true iff T is in TYPES
template <typename... TYPES>
struct eq {
  template <typename T>
  using type = std::integral_constant<bool, (std::is_same_v<TYPES,T> || ...)>;
};

/// Evaluates to true iff T is not in TYPES
template <typename... TYPES>
struct neq {
  template <typename T>
  using type = std::integral_constant<bool, (!std::is_same_v<TYPES,T> && ...)>;
};

/// Return type list of all FIELDS meeting COND
template <template <typename> class COND, typename HEAD=void, typename... TAIL>
struct filter {
  using type = std::conditional_t<
    COND<HEAD>::value && !std::is_void_v<HEAD>,
    typename filter<COND, TAIL...>::type::template push<HEAD>,
    typename filter<COND, TAIL...>::type
  >;
};

/// Return either nil type list or type list containing FIELD depending on COND
template <template <typename> class COND, typename TYPE>
struct filter<COND,TYPE> {
  using type = std::conditional_t<
    COND<TYPE>::value && !std::is_void<TYPE>::value,
    list<TYPE>,
    list<>
  >;
};

/// meta::list of TYPES meeting COND
template <template <typename> class COND, typename... TYPES>
using filter_t = typename filter<COND,TYPES...>::type;


/// Return type list of all FIELDS in reversed order
template <typename HEAD=void, typename... TAIL>
struct reverse {
  using type = std::conditional_t<!std::is_void_v<HEAD>,
    typename reverse<TAIL...>::type::template append<HEAD>,
    typename reverse<TAIL...>::type
  >;
};

/// Return either nil type list or type list containing FIELD in reversed order
template <typename TYPE>
struct reverse<TYPE> {
  using type = std::conditional_t<
    !std::is_void<TYPE>::value,
    list<TYPE>,
    list<>
  >;
};

/// meta::list of TYPES in reversed order
template <typename... TYPES>
using reverse_t = typename reverse<TYPES...>::type;

/// Plain wrapper for list of types
template <typename... TYPES>
struct list {
  static constexpr unsigned size = sizeof...(TYPES);

  /// Returns INDEXth type of TYPES
  template <unsigned INDEX>
  using get = typename std::tuple_element<INDEX, std::tuple<id<TYPES>...>>::type::type;

  /// Export TYPES into arbitrary variadic template COLLECTION
  /**
   * Commonly used to instantiate data classes for a meta::list of fields
   **/
  template <template<typename...> class COLLECTION>
  using decompose_into = COLLECTION<TYPES...>;

  template <typename TYPE>
  using push = list<TYPE, TYPES...>;

  template <typename... UYPES>
  using append = list<TYPES..., UYPES...>;

  /// Merge TYPES and UYPES into new list
  template <typename... UYPES>
  using include = typename filter_t<neq<TYPES...>::template type, UYPES...>::template decompose_into<
    list<TYPES...>::append
  >;

  /// Returns first type of TYPES that is derived from BASE
  template <typename BASE>
  using first_with_base = first_type_with_base_t<BASE, TYPES...>;

  /// Returns first type of TYPES that is derived from BASE
  /**
   * If such a type doesn't exist, FALLBACK is _returned_.
   **/
  template <typename BASE, typename FALLBACK>
  using first_with_base_or_fallback = std::conditional_t<
    std::is_void_v<first_with_base<BASE>>,
    FALLBACK,
    first_with_base<BASE>
  >;

  /// Index of first instance of TYPE in TYPES
  template <typename TYPE>
  static constexpr unsigned index() {
    return index_of_first_matching<eq<TYPE>::template type, TYPES...>::value;
  }

  /// Calls f for each type of TYPES by-value (in reversed order!)
  template <typename F>
  static constexpr void for_each(F f) {
    swallow((f(id<TYPES>()), 0)...);
  }

  template <typename TYPE>
  static constexpr bool contains() {
    return olb::meta::contains<TYPE,TYPES...>();
  }

};

/// Apply F to each element of TUPLE listed in INDICES
template <typename TUPLE, typename F, std::size_t... INDICES>
void tuple_for_each_index(TUPLE& tuple, F&& f, std::index_sequence<INDICES...>)
{
  if constexpr (std::is_invocable_v<F, decltype(std::get<0>(tuple)), unsigned>) {
    swallow((f(std::get<INDICES>(tuple), INDICES), INDICES)...);
  } else {
    swallow((f(std::get<INDICES>(tuple)), INDICES)...);
  }
}

/// Apply F to each element of TUPLE
template <typename TUPLE, typename F>
void tuple_for_each(TUPLE& tuple, F&& f)
{
  tuple_for_each_index(tuple, std::forward<F>(f), std::make_index_sequence<std::tuple_size<TUPLE>::value> {});
}

/// Return std::array<T,D> where T is initialized with a common value (helper)
template <typename T, unsigned D, typename U, std::size_t... INDICES>
std::array<T,D> make_array(U&& u, std::index_sequence<INDICES...>)
{
  return std::array<T,D> {(swallow(INDICES), u)...};
}

/// Return std::array<T,D> where T is initialized with a common value
template <typename T, unsigned D, typename U>
std::array<T,D> make_array(U&& u)
{
  return make_array<T,D,U>(std::forward<U>(u), std::make_index_sequence<D> {});
}

/// Return std::array<T,D> where T is initialized using a iDim-dependent function (helper)
template <typename T, unsigned D, typename F, std::size_t... INDICES>
std::array<T,D> make_array_f(F&& f, std::index_sequence<INDICES...>)
{
  return std::array<T,D> {f(INDICES)...};
}

/// Return std::array<T,D> where T is initialized using a iDim-dependent function
template <typename T, unsigned D, typename F>
std::array<T,D> make_array_f(F&& f)
{
  return make_array_f<T,D,F>(std::forward<F&&>(f), std::make_index_sequence<D> {});
}

/// Returns true iff at least one VALUE satisfies COND
template <auto... VALUE, typename COND>
constexpr bool indexed_pack_contains(COND cond) {
  std::size_t i = 0;
  return (cond(i++, VALUE) || ...);
}

#if !defined __clang__ || __clang_major__ >= 12

/// Concatenate two index sequences
template <std::size_t... Is, std::size_t... Js>
constexpr auto operator+(std::index_sequence<Is...>, std::index_sequence<Js...>)
{
  return std::index_sequence<Is..., Js...>();
}

#endif

/// Convert index sequence into an array of its values
template <std::size_t... Is>
constexpr auto array_from_index_sequence(std::index_sequence<Is...>)
{
  return std::array<std::size_t,sizeof...(Is)>{Is...};
}

/// Return index sequence of Is matching PREDICATE
template <typename PREDICATE, std::size_t... Is>
constexpr auto filter_index_sequence(PREDICATE predicate, std::index_sequence<Is...>) {
  return (
      std::conditional_t<predicate(Is),
                         std::index_sequence<Is>,
                         std::index_sequence<>>()
    + ...
    + std::index_sequence<>()
  );
}

template <typename MAP, std::size_t... Is>
constexpr auto map_index_sequence(MAP map, std::index_sequence<Is...>) {
  return std::index_sequence<map(Is)...>();
}

template <std::size_t N, std::size_t I, std::size_t... Is>
constexpr auto take_n_sequence(std::index_sequence<I,Is...>) {
  if constexpr (N > 0) {
    return std::index_sequence<I>() + take_n_sequence<N-1>(std::index_sequence<Is...>());
  } else {
    return std::index_sequence<>();
  }
  __builtin_unreachable();
}

template <std::size_t N, std::size_t I, std::size_t... Is>
constexpr auto drop_n_sequence(std::index_sequence<I,Is...>) {
  if constexpr (N > 0) {
    return drop_n_sequence<N-1>(std::index_sequence<Is...>());
  } else {
    return std::index_sequence<I,Is...>();
  }
  __builtin_unreachable();
}

template <std::size_t N>
constexpr auto zero_sequence() {
  return map_index_sequence([](std::size_t) constexpr {
    return 0;
  }, std::make_index_sequence<N>());
}

template <std::size_t... Is>
constexpr bool is_zero_sequence(std::index_sequence<Is...>) {
  return ((Is == 0) && ... && true);
}

/// Call F for each index (exlicitly unrolled loop)
template <typename F, std::size_t... INDICES>
void call_n_times(F&& f, std::index_sequence<INDICES...>) {
  (f(INDICES), ...);
}

/// Call F for each i in 0..N-1 (exlicitly unrolled loop)
template <unsigned N, typename F>
void call_n_times(F&& f) {
  return call_n_times<F>(std::forward<F&&>(f), std::make_index_sequence<N>{});
}

template <
  typename BASE,
  typename HEAD = void, //Default argument in case the list is empty
  typename ...TAIL
>
struct derived_type_in_nested {
  auto constexpr static evaluate_type(){
    if constexpr (!sizeof...(TAIL)) {
      if constexpr (std::is_same<HEAD,void>::value) {
        return id<BASE>{};  //Fallback
      } else {
        if constexpr(std::is_same<BASE,void>::value){
          return id<void>{};
        } else {
          using HEAD_EVAL = typename BASE::template derivedField<HEAD>;
          return id<HEAD_EVAL>{};
        }
      }
    } else {
      if constexpr(std::is_same<BASE,void>::value){
         return id<void>{};
      } else {
        using HEAD_EVAL = typename BASE::template derivedField<HEAD>;
        return id<typename derived_type_in_nested<HEAD_EVAL,TAIL...>::type>{};
      }
    }
    __builtin_unreachable();
  }

  using type = typename decltype(evaluate_type())::type;

  static constexpr bool contains() {
    return !std::is_same_v<type, void>;
  }
};


// *INDENT-ON*

}

}

#endif
