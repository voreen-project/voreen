/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Davide Dapelo
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

/** \file
 * Postprocessor to store the particle density, converted to Euler
 *  -- header file
 */
#ifndef EUL2LAGR_POST_PROCESSOR_H
#define EUL2LAGR_POST_PROCESSOR_H

#include "particles/subgrid3DLegacyFramework/particleSystem3D.h"

namespace olb {

////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Virtual base class for Eul2LagrOperator3D.
 * Its raison d'etre lies on the fact that postprocessor classes cannot be templetized in PARTICLETYPE,
 * and hence, a reference to Eul2LagrOperatorBase3D is passed.
 */
template<typename T, typename DESCRIPTOR>
class Eul2LagrOperatorBase3D {
public:
  virtual bool operator()(BlockLattice<T,DESCRIPTOR>& blockLattice) =0;
protected:
  Eul2LagrOperatorBase3D() {};
};

/*
 * Stores the particle density, converted to Euler
 */
template<typename T, typename DESCRIPTOR, template<typename U> class PARTICLETYPE>
class Eul2LagrOperator3D final : public Eul2LagrOperatorBase3D<T,DESCRIPTOR> {
public:
  Eul2LagrOperator3D(ParticleSystem3D<T,PARTICLETYPE>& pSystem, SuperGeometry<T,3>& superGeometry);
  bool operator()(BlockLattice<T,DESCRIPTOR>& blockLattice) override;
private:
  SuperGeometry<T,3>& _superGeometry;
  ParticleSystem3D<T,PARTICLETYPE>& _pSystem;
};


////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Postprocessor to store the particle density, converted to Euler.
 * Due to the impossibility of templating the class on PARTICLETYPE, the real operations
 * are outsourced in Eul2LagrOperator3D.
 */
template<typename T, typename DESCRIPTOR>
class Eul2LagrPostProcessor3D final : public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  Eul2LagrPostProcessor3D ( int x0, int x1, int y0, int y1, int z0, int z1,
                            std::shared_ptr<Eul2LagrOperatorBase3D<T,DESCRIPTOR>> eul2LagrOperator );
  Eul2LagrPostProcessor3D(std::shared_ptr<Eul2LagrOperatorBase3D<T,DESCRIPTOR>> eul2LagrOperator);
  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    return 1;
  }
  void process(BlockLattice<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain ( BlockLattice<T,DESCRIPTOR>& blockLattice,
                          int x0, int x1, int y0, int y1, int z0, int z1 ) override;
private:
  int _x0, _x1, _y0, _y1, _z0, _z1;
  std::shared_ptr<Eul2LagrOperatorBase3D<T,DESCRIPTOR>> _eul2LagrOperator;
};

/*
 * Generator of Eul2LagrPostProcessor3D
 */
template<typename T, typename DESCRIPTOR, template<typename U> class PARTICLETYPE>
class Eul2LagrPostProcessorGenerator3D final : public PostProcessorGenerator3D<T,DESCRIPTOR> {
public:
  Eul2LagrPostProcessorGenerator3D ( int x0_, int x1_, int y0_, int y1_, int z0_, int z1_,
                                     SuperParticleSystem3D<T,PARTICLETYPE>& spSys, SuperGeometry<T,3>& superGeometry );
  Eul2LagrPostProcessorGenerator3D(SuperParticleSystem3D<T,PARTICLETYPE>& spSys, SuperGeometry<T,3>& superGeometry);
  PostProcessor3D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator3D<T,DESCRIPTOR>* clone() const override;
private:
  std::shared_ptr<Eul2LagrOperator3D<T,DESCRIPTOR,PARTICLETYPE>> _eul2LagrOperator;
};


////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Postprocessor test AnalyticalRandomNormal functional
 * through storing a Normal density to eul2Lagr external field.
 */
template<typename T, typename DESCRIPTOR>
class Eul2LagrNormDistrPostProcessor3D final : public LocalPostProcessor3D<T,DESCRIPTOR> {
public:
  Eul2LagrNormDistrPostProcessor3D(int x0, int x1, int y0, int y1, int z0, int z1, T mean, T stdDev, SuperGeometry<T,3>& superGeometry);
  Eul2LagrNormDistrPostProcessor3D(T mean, T stdDev, SuperGeometry<T,3>& superGeometry);
  int extent() const override
  {
    return 1;
  }
  int extent(int whichDirection) const override
  {
    return 1;
  }
  void process(BlockLattice<T,DESCRIPTOR>& blockLattice) override;
  void processSubDomain ( BlockLattice<T,DESCRIPTOR>& blockLattice,
                          int x0, int x1, int y0, int y1, int z0, int z1 ) override;
private:
  int _x0, _x1, _y0, _y1, _z0, _z1;
  AnalyticalRandomNormal<3,T,T> _randomNormal;
  SuperGeometry<T,3>& _superGeometry;
};

/*
 * Generator of Eul2LagrNormDistrPostProcessor3D
 */
template<typename T, typename DESCRIPTOR>
class Eul2LagrNormDistrPostProcessorGenerator3D final : public PostProcessorGenerator3D<T,DESCRIPTOR> {
public:
  Eul2LagrNormDistrPostProcessorGenerator3D(int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, T mean, T stdDev, SuperGeometry<T,3>& superGeometry);
  Eul2LagrNormDistrPostProcessorGenerator3D(T mean, T stdDev, SuperGeometry<T,3>& superGeometry);
  PostProcessor3D<T,DESCRIPTOR>* generate() const override;
  PostProcessorGenerator3D<T,DESCRIPTOR>* clone() const override;
private:
  T _mean;
  T _stdDev;
  SuperGeometry<T,3>& _superGeometry;
};



} // olb

#endif
