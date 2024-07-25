/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Adrian Kummerlaender, Mathias J. Krause
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

#ifndef DESCRIPTOR_BASE_H
#define DESCRIPTOR_BASE_H

#include "descriptorTag.h"
#include "descriptorField.h"
#include "core/meta.h"

namespace olb {

namespace descriptors {

template <typename FIELD>
using is_data_field = std::integral_constant<bool, !is_tag_field<FIELD>::value>;

/// \defgroup descriptor
//@{

/// Identificator types for descriptors
struct TYPE_LATTICE;
struct TYPE_PARTICLE;


/// Tuple of abstract field declarations
template <typename... FIELDS>
struct FIELD_TUPLE : public meta::list<FIELDS...> {
  /// Deleted constructor to enforce pure usage as type and prevent implicit narrowing conversions
  FIELD_TUPLE() = delete;

  /// Fieldlist, which can be passed if desired
  using fieldList = meta::list<FIELDS...>;

  template <template<typename...> class COLLECTION>
  using decompose_into = COLLECTION<FIELDS...>;

  /// Number of fields
  static constexpr std::size_t field_count = sizeof...(FIELDS);

  /// Returns FIELD_TUPLE with subset of FIELDS meeting COND
  template <template<typename> class COND>
  using filter = typename meta::filter_t<COND, FIELDS...>
    ::template decompose_into<FIELD_TUPLE>;

  /// Returns whether WANTED_FIELD is contained in FIELDS
  template <typename WANTED_FIELD>
  static constexpr bool provides()
  {
    return meta::contains<WANTED_FIELD,FIELDS...>();
  }

  /// Returns whether field is contained in FIELDS
  static bool provides(std::type_index field)
  {
    return meta::contains<FIELDS...>(field);
  }

};

/// Tuple of parameters to concretize field declarations
template <unsigned... PARAMS>
struct PARAMETER_TUPLE {
  template <typename FIELD>
  static constexpr std::size_t size() {
    return FIELD::template size<PARAMS...>();
  };
};

/// Tuple of concretized field declarations
template <typename PARAMETERS, typename... FIELDS>
struct CONCRETE_FIELD_TUPLE : public FIELD_TUPLE<FIELDS...> {
  using parameter_tuple_t = PARAMETERS;

  /// Returns FIELD by provided BASE
  /**
   * returns void, if FIELD_TUPE does not provide BASE field
   **/
  template<typename BASE>
  using derivedField = std::conditional_t<
    meta::contains<
      typename meta::list_item_with_base_default_base<BASE,FIELDS...>::type,
      FIELDS...
    >(),
    typename meta::list_item_with_base_default_base<BASE,FIELDS...>::type,
    void
  >;

  /// Returns whether WANTED_FIELD (last in NESTED_FIELDS) is provided by CONCRETE_FIELD_TUPLE
  /**
   * Offers the same funcionality as it overloaded pendant in in FIELD_TUPLE, however allows beeing called
   * with multiple BASE types (as currently used in the particle framework)
   **/
  template <typename ...NESTED_FIELDS>
  static constexpr bool providesNested()
  {
    return meta::derived_type_in_nested<CONCRETE_FIELD_TUPLE, NESTED_FIELDS...>::contains();
  }

  /// Returns size of FIELD if given or of whole tuple if not
  /**
   * Works for all fields, even if they are not provided by this instantiation of CONCRETE_FIELD_TUPLE.
   **/
  template <typename FIELD = void>
  static constexpr std::size_t size()
  {
    if constexpr (std::is_void_v<FIELD>) {
      return (parameter_tuple_t::template size<FIELDS>() + ... + 0);
    } else {
      return parameter_tuple_t::template size<FIELD>();
    }
    __builtin_unreachable();
  }
};

/// Base descriptor of a d-dimensional system
template <unsigned D, typename... FIELDS>
struct SPATIAL_DESCRIPTOR : public CONCRETE_FIELD_TUPLE<PARAMETER_TUPLE<D>, FIELDS...> {
  /// Number of dimensions
  static constexpr int d = D;
};

/// Base descriptor of a D-dimensional lattice with Q directions and a list of additional fields
template <unsigned D, unsigned Q, typename... FIELDS>
struct LATTICE_DESCRIPTOR : public CONCRETE_FIELD_TUPLE<PARAMETER_TUPLE<D,Q>, FIELDS...> {
  /// Type identifier
  using type = TYPE_LATTICE;
  /// Number of dimensions
  static constexpr int d = D;
  /// Number of velocities
  static constexpr int q = Q;

  using all_t = meta::list<FIELDS...>;
  using fields_t = meta::filter_t<is_data_field, FIELDS...>;
  using tags_t   = meta::filter_t<is_tag_field, FIELDS...>;

  /// Tag that describes the category of the descriptor
  /**
   * e.g. tag::RTLBM for RTLBM descriptors or tag::DEFAULT if no category tag-field is found.
   *
   * This is needed to enable tag-dispatching of descriptor function overloads.
   **/
  using category_tag = typename tags_t::template first_with_base_or_fallback<tag::CATEGORY, tag::DEFAULT>;

};

/// Base descriptor of a particle system
template <unsigned D, typename... FIELDS>
struct PARTICLE_DESCRIPTOR : public CONCRETE_FIELD_TUPLE<PARAMETER_TUPLE<D>, FIELDS...> {
  /// Type identifier
  using type = TYPE_PARTICLE;
  /// Number of dimensions
  static constexpr int d = D;




};

//@}

}

}

#endif
