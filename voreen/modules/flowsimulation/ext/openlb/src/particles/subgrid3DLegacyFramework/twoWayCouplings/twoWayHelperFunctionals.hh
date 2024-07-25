/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2019 Davide Dapelo
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

/* Helper functionals for Lagrangian two-way coupling methods -- generic implementation.
 */

#ifndef LB_TWO_WAY_HELPER_FUNCTIONALS_HH
#define LB_TWO_WAY_HELPER_FUNCTIONALS_HH

namespace olb {

////////////////////// Class TwoWayHelperFunctional ////////////////////////

template<typename T, typename Lattice>
TwoWayHelperFunctional<T, Lattice>::TwoWayHelperFunctional (
  UnitConverter<T, Lattice>& converter,
  SuperLattice<T, Lattice>& sLattice )
  : _converter(converter),
    _sLattice(sLattice)
{}

template<typename T, typename Lattice>
TwoWayHelperFunctional<T, Lattice>::~TwoWayHelperFunctional()
{
  _interpLatticeDensity.reset();
  _interpLatticeVelocity.reset();
}


////////////////////// Class NaiveMomentumExchange ////////////////////////

template<typename T, typename Lattice>
NaiveMomentumExchange<T, Lattice>::NaiveMomentumExchange (
  UnitConverter<T, Lattice>& converter,
  SuperLattice<T, Lattice>& sLattice,
  std::shared_ptr<SuperLatticeInterpDensity3Degree3D<T, Lattice> > interpLatticeDensity )
  : TwoWayHelperFunctional<T, Lattice>(converter, sLattice)
{
  this->_interpLatticeDensity = interpLatticeDensity;
}

template<typename T, typename Lattice>
bool NaiveMomentumExchange<T, Lattice>::operator() ( T gF[], T latticeVelF[], T latticeVelP[],
    T physPosP[], int latticeRoundedP[],
    int globic )
{
  T magU = util::sqrt( util::pow(latticeVelF[0] - latticeVelP[0],2) +
                       util::pow(latticeVelF[1] - latticeVelP[1],2) +
                       util::pow(latticeVelF[2] - latticeVelP[2],2) );

  int locIC = this->_sLattice.getLoadBalancer().loc(globic);
  T rhoL = this->_sLattice.getBlock(locIC).get(
             latticeRoundedP[0],
             latticeRoundedP[1],
             latticeRoundedP[2]).computeRho();

  gF[0] = rhoL * magU;
  gF[1] = rhoL * magU;
  gF[2] = rhoL * magU;

  return true;
}


////////////////////// Class LaddMomentumExchange ////////////////////////

template<typename T, typename Lattice>
LaddMomentumExchange<T, Lattice>::LaddMomentumExchange (
  UnitConverter<T, Lattice>& converter,
  SuperLattice<T, Lattice>& sLattice,
  std::shared_ptr<SuperLatticeInterpDensity3Degree3D<T, Lattice> > interpLatticeDensity,
  std::shared_ptr<SuperLatticeInterpPhysVelocity3D<T, Lattice> > interpLatticeVelocity )
  : TwoWayHelperFunctional<T, Lattice>(converter, sLattice)
{
  this->_interpLatticeDensity = interpLatticeDensity;
  this->_interpLatticeVelocity = interpLatticeVelocity;
}

template<typename T, typename Lattice>
bool LaddMomentumExchange<T, Lattice>::operator() ( T gF[], T latticeVelF[], T latticeVelP[],
    T physPosP[], int latticeRoundedP[],
    int globic )
{
  T physLatticeL = this->_converter.getConversionFactorLength();

  // force density gF
  gF[0] = T();
  gF[1] = T();
  gF[2] = T();

  T fiPop = T();
  T sp = T(); // dot product for particle velocity
  T faPos[3] = {T(), T(), T()}; // fAlphaPosition = particle position
  T fbPos[3] = {T(), T(), T()}; // fBetaPosition = neighbor position to particle position in direction iPop

  T fa[Lattice::q] = { T() }; // fAlpha = interpolated density to fAlphaPosition
  T fb[Lattice::q] = { T() }; // fBeta = interpolated density to fBetaPosition
  T lFU[3] = {T(), T(), T()};

  // runs through all q discrete velocity directions
  for (unsigned iPop = 0; iPop < Lattice::q; ++iPop) {
    // physical position on neighbor to get pre-streaming collision part
    faPos[0] = physPosP[0] + physLatticeL * descriptors::c<Lattice>(iPop,0);
    faPos[1] = physPosP[1] + physLatticeL * descriptors::c<Lattice>(iPop,1);
    faPos[2] = physPosP[2] + physLatticeL * descriptors::c<Lattice>(iPop,2);
    // Lagrange interpolated polynomial to get density on particle position
    this->_interpLatticeDensity->operator() (fa, faPos, globic);

    // physical position on neighbor to get pre-streaming collision part
    fbPos[0] = physPosP[0] - physLatticeL * descriptors::c<Lattice>(iPop,0);
    fbPos[1] = physPosP[1] - physLatticeL * descriptors::c<Lattice>(iPop,1);
    fbPos[2] = physPosP[2] - physLatticeL * descriptors::c<Lattice>(iPop,2);
    // Lagrange interpolated polynomial to get density on particle position
    this->_interpLatticeDensity->operator() (fb, fbPos, globic);

    // fiPop = density on fBetaPosition in direction iPop
    fiPop = fb[util::opposite<Lattice >(iPop)];
    // Get f_l of the boundary cell
    // add density fAlphaL of opposite direction to iPop
    fiPop -= fa[iPop];

    // physical velocity
    lFU[0] = -descriptors::c<Lattice>(iPop,0) * fiPop;
    lFU[1] = -descriptors::c<Lattice>(iPop,1) * fiPop;
    lFU[2] = -descriptors::c<Lattice>(iPop,2) * fiPop;

    // point product
    sp = descriptors::c<Lattice>(iPop,0) * latticeVelP[0] + descriptors::c<Lattice>(iPop,1) * latticeVelP[1]
         + descriptors::c<Lattice>(iPop,2) * latticeVelP[2];
    sp *= 2. * descriptors::invCs2<T,Lattice>()
          * descriptors::t<T,Lattice>(iPop);

    // external force density that acts on particles
    gF[0] += (lFU[0] - descriptors::c<Lattice>(iPop,0) * (sp));
    gF[1] += (lFU[1] - descriptors::c<Lattice>(iPop,1) * (sp));
    gF[2] += (lFU[2] - descriptors::c<Lattice>(iPop,2) * (sp));
  }
  gF[0] = util::fabs(gF[0]);
  gF[1] = util::fabs(gF[1]);
  gF[2] = util::fabs(gF[2]);

  return true;
}


}

#endif
