/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2016 Robin Trunk
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

#ifndef ADVECTION_DIFFUSION_FORCES_H
#define ADVECTION_DIFFUSION_FORCES_H

#include "core/unitConverter.h"

namespace olb {

template<typename T, template<typename U> class Lattice,
template<typename U> class ADLattice=descriptors::ParticleAdvectionDiffusionD3Q7Descriptor>
class advectionDiffusionForce3D {
public:
  advectionDiffusionForce3D()
  {
    initArg = 0;
  };
  virtual ~advectionDiffusionForce3D() {};
  virtual void applyForce(T force[], Cell<T,Lattice> *nsCell, Cell<T,ADLattice> *adCell, T vel[], int latticeR[])=0;
  int getInitArg()
  {
    return initArg;
  }
private:
  int initArg;
};

template<typename T, template<typename U> class Lattice,
template<typename U> class ADLattice=descriptors::ParticleAdvectionDiffusionD3Q7Descriptor>
class advDiffDragForce3D : public advectionDiffusionForce3D<T,Lattice,ADLattice> {
public:
  advDiffDragForce3D(UnitConverter<T,Lattice> const& converter_, T St_);
  advDiffDragForce3D(UnitConverter<T,Lattice> const& converter_, T pRadius_, T pRho_);
  ~advDiffDragForce3D() override {};
  void applyForce(T force[], Cell<T,Lattice> *nsCell, Cell<T,ADLattice> *adCell, T vel[], int latticeR[]) override;

private:
  int initArg;
  T dragCoeff;
};

template<typename T, template<typename U> class Lattice,
template<typename U> class ADLattice=descriptors::ParticleAdvectionDiffusionD3Q7Descriptor>
class advDiffRotatingForce3D : public advectionDiffusionForce3D<T,Lattice,ADLattice> {
public:
  advDiffRotatingForce3D(SuperGeometry3D<T>& superGeometry_,
                         const UnitConverter<T,Lattice>& converter_,
                         std::vector<T> axisPoint_,
                         std::vector<T> axisDirection_,
                         T w_, T* frac_,
                         bool centrifugeForceOn_ = true,
                         bool coriolisForceOn_ = true);
  advDiffRotatingForce3D(UnitConverter<T,Lattice> const& converter_, T pRadius_, T pRho_);
  virtual ~advDiffRotatingForce3D() {};
  void applyForce(T force[], Cell<T,Lattice> *nsCell, Cell<T,ADLattice> *adCell, T vel[], int latticeR[]);

protected:
  SuperGeometry3D<T>& sg;
  std::vector<T> axisPoint;
  std::vector<T> axisDirection;
  T invMassLessForce;
  T w;
  T* frac;
  bool centrifugeForceOn;
  bool coriolisForceOn;

};

}

#endif
