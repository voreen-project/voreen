/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2017 Adrian Kummerländer
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

#ifndef FUNCTOR_PTR_H
#define FUNCTOR_PTR_H

#include <memory>
#include <type_traits>

namespace olb {

/// Smart pointer for managing the various ways of passing functors around
/**
 * There is a rich set of functors that compose other functors by accepting
 * them as constructor arguments. e.g. SuperLpNorm3D, SuperPlaneIntegralF3D
 *
 * Previously all of these functors only accepted other functors by reference
 * which prevented e.g. the use case of accepting a std::unique_ptr to a
 * indicator freshly created by SuperGeometry3D::getMaterialIndicator.
 *
 * This class hides the memory management required to support accepting either
 * functor references, non-owning pointers or owning unique pointers in a single
 * argument.
 **/
template<typename F>
class FunctorPtr {
private:
  /// Optional functor owner
  std::unique_ptr<F> _ownF;
  /// Pointer to the exposed functor
  F* const _f;

public:
  /// Constructor for transparently accepting a functor reference
  FunctorPtr(F& f);
  /// Constructor for transparently accepting a non-owning functor pointer
  FunctorPtr(F* f);
  /// Constructor for transparently accepting a owning functor unique pointer
  FunctorPtr(std::unique_ptr<F>&& f);

  /// Copy construction is disabled as it is incompatible with unique ownership
  FunctorPtr(FunctorPtr&) = delete;
  /// Move constructor
  FunctorPtr(FunctorPtr&&) = default;

  /// Perfect forwarding functor operator
  template<typename... Args>
  bool operator()(Args... args);

  /// \return reference to the exposed functor
  typename std::add_lvalue_reference<F>::type operator*() const;
  /// Enable pointer-like access to the exposed functor's members
  typename std::add_pointer<F>::type operator->() const noexcept;

  /// \return true iff a functor is exposed
  /**
   * This is useful for supporting optional functors in constructor interfaces
   * by constructing FunctorPtr using a nullptr
   **/
  operator bool() const;
};

}

#endif
