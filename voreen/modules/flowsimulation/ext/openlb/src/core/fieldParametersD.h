/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Adrian Kummerlaender
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

#ifndef FIELD_PARAMETERS_D_H
#define FIELD_PARAMETERS_D_H

#include "fieldArrayD.h"

#include <stdexcept>

namespace olb {

/// Storage of a single FIELD-valued parameter
template <typename T, typename DESCRIPTOR, typename FIELD>
class ParameterD {
private:
  std::conditional_t<
    DESCRIPTOR::template size<FIELD>() == 1,
    typename FIELD::template value_type<T>,
    FieldD<T,DESCRIPTOR,FIELD>
  > _data;

public:
  using data_t = decltype(_data);

  ParameterD() = default;

  template <typename V>
  ParameterD(const ParameterD<V,DESCRIPTOR,FIELD>& rhs) any_platform:
    _data(rhs.read()) { }

  const auto& read() const any_platform {
    return _data;
  }

  auto& read() any_platform {
    return _data;
  }

  void write(FieldD<T,DESCRIPTOR,FIELD>&& value) any_platform {
    if constexpr (DESCRIPTOR::template size<FIELD>() == 1) {
      _data = value[0];
    } else {
      _data = value;
    }
  }

  void write(const FieldD<T,DESCRIPTOR,FIELD>& value) any_platform {
    if constexpr (DESCRIPTOR::template size<FIELD>() == 1) {
      _data = value[0];
    } else {
      _data = value;
    }
  }
};

/// Dynamic access interface for FIELD-valued parameters
template <typename T, typename DESCRIPTOR>
struct AbstractParameters {
  virtual ~AbstractParameters() = default;

  template <typename FIELD>
  bool provides() const {
    if (auto concrete = dynamic_cast<const ParameterD<T,DESCRIPTOR,FIELD>*>(this)) {
      return true;
    } else {
      return false;
    }
  };

  template <typename FIELD>
  auto get() const {
    if (auto concrete = dynamic_cast<const ParameterD<T,DESCRIPTOR,FIELD>*>(this)) {
      return concrete->read();
    } else {
      throw std::invalid_argument("FIELD not provided by this parameter set");
    }
  };

  template <typename FIELD>
  auto getOrFallback(typename ParameterD<T,DESCRIPTOR,FIELD>::data_t&& fallback) const {
    if (auto concrete = dynamic_cast<const ParameterD<T,DESCRIPTOR,FIELD>*>(this)) {
      return concrete->read();
    } else {
      return fallback;
    }
  }

  template <typename FIELD>
  void set(FieldD<T,DESCRIPTOR,FIELD>&& value) {
    if (auto concrete = dynamic_cast<ParameterD<T,DESCRIPTOR,FIELD>*>(this)) {
      return concrete->write(std::forward<decltype(value)>(value));
    } else {
      throw std::invalid_argument("FIELD not provided by this parameter set");
    }
  };

  template <typename FIELD>
  void set(const FieldD<T,DESCRIPTOR,FIELD>& value) {
    if (auto concrete = dynamic_cast<ParameterD<T,DESCRIPTOR,FIELD>*>(this)) {
      return concrete->write(value);
    } else {
      throw std::invalid_argument("FIELD not provided by this parameter set");
    }
  };

};

/// Set of FIELD-valued parameters
template <typename T, typename DESCRIPTOR, typename... FIELDS>
struct ParametersD final : public AbstractParameters<T,DESCRIPTOR>
                         , public ParameterD<T,DESCRIPTOR,FIELDS>...
{
  using fields_t = meta::list<FIELDS...>;

  template <typename... Fs>
  using include_fields = ParametersD<T,DESCRIPTOR,FIELDS...,Fs...>;

  ParametersD() = default;

  template <typename V>
  ParametersD(const ParametersD<V,DESCRIPTOR,FIELDS...>& rhs) any_platform :
    ParameterD<T,DESCRIPTOR,FIELDS>(static_cast<const ParameterD<T,DESCRIPTOR,FIELDS>&>(rhs))...
  { }

  template <typename FIELD>
  bool provides() const any_platform {
    return fields_t::template contains<FIELD>();
  }

  template <typename FIELD>
  const auto& get() const any_platform {
    return static_cast<const ParameterD<T,DESCRIPTOR,FIELD>*>(this)->read();
  }

  template <typename V>
  ParametersD<V,DESCRIPTOR,FIELDS...> copyAs() const any_platform {
    return ParametersD<V,DESCRIPTOR,FIELDS...>(*this);
  }

  template <typename FIELD>
  auto& get() any_platform {
    return static_cast<ParameterD<T,DESCRIPTOR,FIELD>*>(this)->read();
  }

  template <typename FIELD>
  void set(FieldD<T,DESCRIPTOR,FIELD>&& value) any_platform {
    return static_cast<ParameterD<T,DESCRIPTOR,FIELD>*>(this)->write(
      std::forward<decltype(value)>(value));
  };

  template <typename FIELD>
  void set(const FieldD<T,DESCRIPTOR,FIELD>& value) any_platform {
    return static_cast<ParameterD<T,DESCRIPTOR,FIELD>*>(this)->write(value);
  };

};

}

#endif
