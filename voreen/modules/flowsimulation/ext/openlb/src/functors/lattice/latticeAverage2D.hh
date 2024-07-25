/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Albert Mink, Mathias J. Krause, Lukas Baron
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

#ifndef LATTICE_AVERAGE_2D_HH
#define LATTICE_AVERAGE_2D_HH

#include <vector>
#include "utilities/omath.h"
#include <limits>

#include "latticeAverage2D.h"
#include "dynamics/lbm.h"  // for computation of lattice rho and velocity
#include "geometry/superGeometry.h"
#include "indicator/superIndicatorF2D.h"
#include "blockBaseF2D.h"
#include "functors/genericF.h"
#include "functors/analytical/analyticalF.h"
#include "functors/analytical/indicator/indicatorF2D.h"
#include "communication/mpiManager.h"


namespace olb {

template <typename T, typename DESCRIPTOR>
BlockLatticeAverage2D<T,DESCRIPTOR>::BlockLatticeAverage2D
(BlockLatticeF2D<T,DESCRIPTOR>& f, BlockGeometry<T,2>& blockGeometry,
 int material, T radius)
  : BlockLatticeF2D<T,DESCRIPTOR>(f.getBlock(), f.getTargetDim()),
    _f(f), _blockGeometry(blockGeometry), _material(material), _radius(radius)
{
  this->getName() = "Average("+f.getName()+")";
}


template <typename T, typename DESCRIPTOR>
bool BlockLatticeAverage2D<T,DESCRIPTOR>::operator() (T output[], const int input[])
{
  //  CuboidGeometry2D<T>& cGeometry = f.getBlockLattice2D().get_cGeometry();
  //  loadBalancer& load = f.getBlockLattice2D().get_load();

  //  //create boolean indicator functor isInSphere
  //  std::vector<T> center(3,0);
  //  center[0] = (int)cGeometry.get(load.glob(input[0])).get_globPosX() + input[1];
  //  center[1] = (int)cGeometry.get(load.glob(input[0])).get_globPosY() + input[2];
  //  center[2] = (int)cGeometry.get(load.glob(input[0])).get_globPosZ() + input[3];
  //  SphereAnalyticalF2D<bool,T> isInSphere(center,radius);

  // iterate over all cuboids & points and test for material && isInSphere
  //  std::vector<T> tmp( this->_n, T() );
  //  int numVoxels(0);
  //  if (this->_blockGeometry.getMaterial({center[0],center[1],center[2]}) == material) {
  //    for (int iC=0; iC<load.size(); iC++) {
  //      int nX = cGeometry.get(load.glob(iC)).getNx();
  //      int nY = cGeometry.get(load.glob(iC)).getNy();
  //      int nZ = cGeometry.get(load.glob(iC)).getNz();
  //      for (int iX=0; iX<nX; ++iX) {
  //        for (int iY=0; iY<nY; ++iY) {
  //          for (int iZ=0; iZ<nZ; iZ++) {
  //            std::vector<T> glob(3,0);
  //            glob[0] = (int)cGeometry.get(load.glob(iC)).get_globPosX() + iX;
  //            glob[1] = (int)cGeometry.get(load.glob(iC)).get_globPosY() + iY;
  //            glob[2] = (int)cGeometry.get(load.glob(iC)).get_globPosZ() + iZ;
  //            if (this->_blockGeometry.getMaterial({glob[0],glob[1],glob[2]}) == material
  //                && isInSphere(glob)[0]==true) {
  //              for (unsigned iD=0; iD<f(load.glob(0),0,0,0).size(); iD++) {
  //                tmp[iD]+=f(load.glob(iC),iX,iY,iZ)[iD];
  //              }
  //              numVoxels++;
  //            }
  //          }
  //        }
  //      }
  //    }

  //#ifdef PARALLEL_MODE_MPI
  //    singleton::mpi().reduceAndBcast(numVoxels, MPI_SUM);
  //#endif
  //    for (int iD=0; iD<f.getTargetDim(); iD++) {
  //#ifdef PARALLEL_MODE_MPI
  //      singleton::mpi().reduceAndBcast(tmp[iD], MPI_SUM);
  //#endif
  //      if (numVoxels>0) {
  //        tmp[iD] /= numVoxels;
  //      }
  //    }
  //  }
  //  return tmp;

  return false;
}

}
#endif
