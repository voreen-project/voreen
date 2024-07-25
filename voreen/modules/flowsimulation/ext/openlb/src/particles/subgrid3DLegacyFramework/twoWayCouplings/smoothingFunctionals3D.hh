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

/* Smoothings functionals for Lagrangian two-way coupling methods -- generic implementation.
 */

#ifndef LB_SMOOTHING_FUNCTIONALS_3D_HH
#define LB_SMOOTHING_FUNCTIONALS_3D_HH

namespace olb {

////////////////////// Class SmoothingFunctional ////////////////////////

template<typename T, typename Lattice>
SmoothingFunctional<T, Lattice>::SmoothingFunctional (
  T kernelLength, UnitConverter<T, Lattice>& converter, SuperLattice<T, Lattice>& sLattice, int nVoxelInterpPoints )
  : _kernelLength(kernelLength),
    _converter(converter),
    _sLattice(sLattice),
    _nVoxelInterpPoints(nVoxelInterpPoints)
{
  if (2.0*_kernelLength < converter.getPhysDeltaX()) {
    throw std::out_of_range("2 * SmoothingFunctional::_kernelLength must be >= converter::getPhysDeltaX().");
  }
  if (_nVoxelInterpPoints < 2) {
    throw std::out_of_range("SmoothingFunctional::_nVoxelInterpPoints must be  >=2.");
  }
}

template<typename T, typename Lattice>
bool SmoothingFunctional<T, Lattice>::update(T physPosP[], int globic)
{
  // Bottom-left corner of a cube centered at the particle, with side 2*_kernelLength
  T physPosMin[3] = {T(), T(), T()};
  physPosMin[0] = physPosP[0] - _kernelLength;
  physPosMin[1] = physPosP[1] - _kernelLength;
  physPosMin[2] = physPosP[2] - _kernelLength;
  int latticePosMin[3] = {0, 0, 0};
  this->_sLattice.getCuboidGeometry().get(globic).getLatticeR (
    latticePosMin, physPosMin );

  // Top-right corner of a cube centered at the particle, with side 2*_kernelLength
  T physPosMax[3] = {T(), T(), T()};
  physPosMax[0] = physPosP[0] + _kernelLength;
  physPosMax[1] = physPosP[1] + _kernelLength;
  physPosMax[2] = physPosP[2] + _kernelLength;
  int latticePosMax[3] = {0, 0, 0};
  this->_sLattice.getCuboidGeometry().get(globic).getLatticeR (
    latticePosMax, physPosMax );

  // Clearing the _latticePosAndWeight list
  _latticePosAndWeight.clear();

  T normalizer = T();
  int iLatticePos[3] = {0, 0, 0};
  // Cycling all the cells on a cube containing a sphee centered in bubble's position and with kernel smoothing length as radius
  for (iLatticePos[0]=latticePosMin[0]; iLatticePos[0]<=latticePosMax[0]; iLatticePos[0]++) {
    for (iLatticePos[1]=latticePosMin[1]; iLatticePos[1]<=latticePosMax[1]; iLatticePos[1]++) {
      for (iLatticePos[2]=latticePosMin[2]; iLatticePos[2]<=latticePosMax[2]; iLatticePos[2]++) {

        T iPhysPos[3] = {T(), T(), T()};
        this->_sLattice.getCuboidGeometry().get(globic).getPhysR (
          iPhysPos, iLatticePos );

        // Is the voxel within a smooting kernel length from the bubble's position?
        if ( util::pow(physPosP[0] - iPhysPos[0], 2) +
             util::pow(physPosP[1] - iPhysPos[1], 2) +
             util::pow(physPosP[2] - iPhysPos[2], 2) < util::pow(_kernelLength, 2) ) {

          // Adding the voxel's position (and relative weight) to the _latticePosAndWeight list
          LatticePosAndWeight<T> item;
          item.latticePos[0] = iLatticePos[0];
          item.latticePos[1] = iLatticePos[1];
          item.latticePos[2] = iLatticePos[2];
          item.weight = this->compute(physPosP, iPhysPos);

          normalizer += item.weight;
          _latticePosAndWeight.push_back(item);
        }
      }
    }
  }

  // If normalizer is zero, then no voxels are within a kernel smoothing length from the bubble's location.
  // And it is a problem.
  if (normalizer == T()) {
    std::cout << "ERROR: SmoothingFunctional::update(...):" << std::endl
              << "[smoothingFunctional] physPosP: "
              << physPosP[0] << " "
              << physPosP[1] << " "
              << physPosP[2] << std::endl
              << "[smoothingFunctional] physPosMin: "
              << physPosMin[0] << " "
              << physPosMin[1] << " "
              << physPosMin[2] << std::endl
              << "[smoothingFunctional] physPosMax: "
              << latticePosMax[0] << " "
              << latticePosMax[1] << " "
              << latticePosMax[2] << std::endl
              << "[smoothingFunctional] normalizer: "
              << normalizer << std::endl;
    return false;
  }

  // Normalizing to one
  for (auto&& i : _latticePosAndWeight) {
    i.weight /= normalizer;
  }
  return true;
}

////////////////////// Class LinearAveragingSmoothingFunctional ////////////////////////

template<typename T, typename Lattice>
LinearAveragingSmoothingFunctional<T, Lattice>::LinearAveragingSmoothingFunctional (
  T kernelLength, UnitConverter<T, Lattice>& converter, SuperLattice<T, Lattice>& sLattice, int nVoxelInterpPoints )
  : SmoothingFunctional<T, Lattice>(kernelLength, converter, sLattice, nVoxelInterpPoints)
{}

template<typename T, typename Lattice>
T LinearAveragingSmoothingFunctional<T, Lattice>::compute(T physPosP[], T physPosL[])
{
  return this->smoothingFunction ( util::sqrt (
                                     util::pow(physPosP[0] - physPosL[0], 2) +
                                     util::pow(physPosP[1] - physPosL[1], 2) +
                                     util::pow(physPosP[2] - physPosL[2], 2) ) );
}


////////////////////// Class VolumeAveragingSmoothingFunctional ////////////////////////

template<typename T, typename Lattice>
VolumeAveragingSmoothingFunctional<T, Lattice>::VolumeAveragingSmoothingFunctional (
  T kernelLength, UnitConverter<T, Lattice>& converter, SuperLattice<T, Lattice>& sLattice, int nVoxelInterpPoints )
  : SmoothingFunctional<T, Lattice>(kernelLength, converter, sLattice, nVoxelInterpPoints)
{}

template<typename T, typename Lattice>
T VolumeAveragingSmoothingFunctional<T, Lattice>::compute(T physPosP[], T physPosL[])
{
  return this->smoothingFunction(physPosP[0] - physPosL[0])
         * this->smoothingFunction(physPosP[1] - physPosL[1])
         * this->smoothingFunction(physPosP[2] - physPosL[2]);
}


////////////////////// Class DeenSmoothingFunctional ////////////////////////

template<typename T, typename Lattice>
DeenSmoothingFunctional<T, Lattice>::DeenSmoothingFunctional (
  T kernelLength, UnitConverter<T, Lattice>& converter, SuperLattice<T, Lattice>& sLattice )
  : LinearAveragingSmoothingFunctional<T, Lattice>(kernelLength, converter, sLattice)
{}

template<typename T, typename Lattice>
T DeenSmoothingFunctional<T, Lattice>::smoothingFunction(T delta)
{
  return ( util::pow(delta, 4)/util::pow(this->_kernelLength, 5)
           - 2.*util::pow(delta, 2)/util::pow(this->_kernelLength, 3)
           + 1./this->_kernelLength
         );
}


////////////////////// Class vanWachemSmoothingFunctional ////////////////////////

template<typename T, typename Lattice>
vanWachemSmoothingFunctional<T, Lattice>::vanWachemSmoothingFunctional (
  T kernelLength, UnitConverter<T, Lattice>& converter, SuperLattice<T, Lattice>& sLattice, T radius, int nVoxelInterpPoints )
  : LinearAveragingSmoothingFunctional<T, Lattice>(kernelLength, converter, sLattice, nVoxelInterpPoints),
    _radius(radius)
{}

template<typename T, typename Lattice>
bool vanWachemSmoothingFunctional<T, Lattice>::update(T physPosP[], int globic)
{
  if ( ! SmoothingFunctional<T, Lattice>::update(physPosP, globic) ) {
    return false;
  }

  updateContinuousPhaseFraction(physPosP, globic);

  return true;
}

template<typename T, typename Lattice>
T vanWachemSmoothingFunctional<T, Lattice>::smoothingFunction(T delta)
{
  T x = delta / this->_kernelLength;
  return (4.0 * x + 1.0) * util::pow(1.0 - x, 4);
}

template<typename T, typename Lattice>
void vanWachemSmoothingFunctional<T, Lattice>::updateContinuousPhaseFraction(T physPosP[], int globic)
{
  for (auto&& i : this->_latticePosAndWeight) {
    T iPhysPos[3] = {T(), T(), T()};
    this->_sLattice.getCuboidGeometry().get(globic).getPhysR (
      iPhysPos, i.latticePos );

    T iExtremaMin[3] = {T(), T(), T()};
    iExtremaMin[0] = util::max( physPosP[0] - this->_kernelLength, iPhysPos[0] - 0.5 *this->_converter.getPhysDeltaX() );
    iExtremaMin[1] = util::max( physPosP[1] - this->_kernelLength, iPhysPos[1] - 0.5 *this->_converter.getPhysDeltaX() );
    iExtremaMin[2] = util::max( physPosP[2] - this->_kernelLength, iPhysPos[2] - 0.5 *this->_converter.getPhysDeltaX() );

    T iExtremaMax[3] = {T(), T(), T()};
    iExtremaMax[0] = util::min( physPosP[0] + this->_kernelLength, iPhysPos[0] + 0.5 *this->_converter.getPhysDeltaX() );
    iExtremaMax[1] = util::min( physPosP[1] + this->_kernelLength, iPhysPos[1] + 0.5 *this->_converter.getPhysDeltaX() );
    iExtremaMax[2] = util::min( physPosP[2] + this->_kernelLength, iPhysPos[2] + 0.5 *this->_converter.getPhysDeltaX() );

    T discretePhaseFraction = T();
    for (int nx=1; nx<=this->_nVoxelInterpPoints; nx++) {
      for (int ny=1; ny<=this->_nVoxelInterpPoints; ny++) {
        for (int nz=1; nz<=this->_nVoxelInterpPoints; nz++) {
          T physPosNm[3] = {T(), T(), T()};
          physPosNm[0] = iExtremaMin[0] + nx * (iExtremaMax[0] - iExtremaMin[0]) / ((T) this->_nVoxelInterpPoints + 1.0);
          physPosNm[1] = iExtremaMin[1] + ny * (iExtremaMax[1] - iExtremaMin[1]) / ((T) this->_nVoxelInterpPoints + 1.0);
          physPosNm[2] = iExtremaMin[2] + nz * (iExtremaMax[2] - iExtremaMin[2]) / ((T) this->_nVoxelInterpPoints + 1.0);

          discretePhaseFraction += (    util::pow(physPosNm[0] - physPosP[0], 2)
                                        +  util::pow(physPosNm[1] - physPosP[1], 2)
                                        +  util::pow(physPosNm[2] - physPosP[2], 2)
                                        <= util::pow(this->_radius, 2)
                                   ) ?  0.5 : 0.0;

          T physPosNp[3] = {T(), T(), T()};
          physPosNp[0] = iExtremaMin[0] + (nx - 1.0) * (iExtremaMax[0] - iExtremaMin[0]) / ((T) this->_nVoxelInterpPoints - 1.0);
          physPosNp[1] = iExtremaMin[1] + (ny - 1.0) * (iExtremaMax[1] - iExtremaMin[1]) / ((T) this->_nVoxelInterpPoints - 1.0);
          physPosNp[2] = iExtremaMin[2] + (nz - 1.0) * (iExtremaMax[2] - iExtremaMin[2]) / ((T) this->_nVoxelInterpPoints - 1.0);

          discretePhaseFraction += (    util::pow(physPosNp[0] - physPosP[0], 2)
                                        +  util::pow(physPosNp[1] - physPosP[1], 2)
                                        +  util::pow(physPosNp[2] - physPosP[2], 2)
                                        <= util::pow(this->_radius, 2)
                                   ) ?  0.5 : 0.0;
        }
      }
    }

    discretePhaseFraction *= (iExtremaMax[0]-iExtremaMin[0]) * (iExtremaMax[1]-iExtremaMin[1]) * (iExtremaMax[2]-iExtremaMin[2])
                             / util::pow(this->_nVoxelInterpPoints * this->_converter.getPhysDeltaX(), 3);

    i.continuousPhaseFraction = 1.0 - discretePhaseFraction;
  }
}


////////////////////// Class StepSmoothingFunctional ////////////////////////

template<typename T, typename Lattice>
StepSmoothingFunctional<T, Lattice>::StepSmoothingFunctional (
  T kernelLength, UnitConverter<T, Lattice>& converter, SuperLattice<T, Lattice>& sLattice )
  : VolumeAveragingSmoothingFunctional<T, Lattice>(kernelLength, converter, sLattice)
{}

template<typename T, typename Lattice>
T StepSmoothingFunctional<T, Lattice>::smoothingFunction(T delta)
{
  return 1.;
}


}

#endif
