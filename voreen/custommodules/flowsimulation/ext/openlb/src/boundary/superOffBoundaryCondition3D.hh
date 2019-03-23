/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012 Jonas Kratzke, Mathias J. Krause
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
 * A helper for initialising 3D boundaries -- generic implementation.
 */


#ifndef SUPER_OFF_BOUNDARY_CONDITION_3D_HH
#define SUPER_OFF_BOUNDARY_CONDITION_3D_HH

#include <vector>
#include <list>
#include "offBoundaryCondition3D.h"
#include "superOffBoundaryCondition3D.h"
#include "core/superLattice3D.h"
#include "core/util.h"
#include "functors/analytical/analyticalF.h"

namespace olb {

///////// class superOffBoundaryCondition3D ///////////////////////////////

template<typename T, template<typename U> class Lattice>
sOffLatticeBoundaryCondition3D<T,Lattice>::
sOffLatticeBoundaryCondition3D (SuperLattice3D<T,Lattice>& sLattice, T epsFraction )
  : clout(std::cout,"sOffLatticeBoundaryCondition3D"),
    _sLattice(sLattice),
    _epsFraction(epsFraction),
    _output(false)
{}

template<typename T, template<typename U> class Lattice>
sOffLatticeBoundaryCondition3D<T,Lattice>::
sOffLatticeBoundaryCondition3D(sOffLatticeBoundaryCondition3D<T,Lattice> const& rhs)
  : clout(std::cout,"sOffLatticeBoundaryCondition3D"),
    _sLattice(rhs._sLattice),
    _epsFraction(rhs._epsFraction),
    _output(false)
{
  _blockBCs = rhs._blockBCs;
  _overlap = rhs._overlap;
}

template<typename T, template<typename U> class Lattice>
sOffLatticeBoundaryCondition3D<T,Lattice> sOffLatticeBoundaryCondition3D<T,Lattice>::operator=(
  sOffLatticeBoundaryCondition3D<T,Lattice> rhs)
{

  sOffLatticeBoundaryCondition3D<T,Lattice> tmp(rhs);
  return tmp;
}

template<typename T, template<typename U> class Lattice>
sOffLatticeBoundaryCondition3D<T,Lattice>::
~sOffLatticeBoundaryCondition3D()
{

  for (unsigned iC=0; iC<_blockBCs.size(); iC++) {
    delete _blockBCs[iC];
  }
}

template<typename T, template<typename U> class Lattice>
void sOffLatticeBoundaryCondition3D<T,Lattice>::addZeroVelocityBoundary(SuperGeometry3D<T>& superGeometry, int material,
    IndicatorF3D<T>& indicator, std::list<int> bulkMaterials)
{
  if (_output) {
    clout << "epsFraction=" << _epsFraction << std::endl;
    clout.setMultiOutput(true);
  }
  int nCloc = _sLattice.getLoadBalancer().size();
  for (int iCloc = 0; iCloc < nCloc; iCloc++) {
    if (_output) {
      clout << "Cuboid globiC " << _sLattice.getLoadBalancer().glob(iCloc)
        << " starts to read distances for ZeroVelocity Boundary..." << std::endl;
    }
    _blockBCs[iCloc]->addZeroVelocityBoundary(superGeometry.getExtendedBlockGeometry(iCloc), material, indicator, bulkMaterials);
    if (_output) {
      clout << "Cuboid globiC " << _sLattice.getLoadBalancer().glob(iCloc)
        << " finished reading distances for ZeroVelocity Boundary." << std::endl;
    }
  }
  if (_output)
    clout.setMultiOutput(false);
  addPoints2CommBC(superGeometry, material);
}

template<typename T, template<typename U> class Lattice>
void sOffLatticeBoundaryCondition3D<T,Lattice>::
addVelocityBoundary(SuperGeometry3D<T>& superGeometry, int material, IndicatorF3D<T>& indicator, std::list<int> bulkMaterials)
{
  if (_output) {
    clout << "epsFraction=" << _epsFraction << std::endl;
    clout.setMultiOutput(true);
  }
  int nC = _sLattice.getLoadBalancer().size();
  for (int iCloc = 0; iCloc < nC; iCloc++) {
    if (_output) {
      clout << "Cuboid globiC " << _sLattice.getLoadBalancer().glob(iCloc)
        << " starts to read distances for Velocity Boundary..." << std::endl;
    }
    _blockBCs[iCloc]->addVelocityBoundary(superGeometry.getExtendedBlockGeometry(iCloc), material, indicator, bulkMaterials);
    if (_output) {
      clout << "Cuboid globiC " << _sLattice.getLoadBalancer().glob(iCloc)
        << " finished reading distances for Velocity Boundary." << std::endl;
    }
  }
  if (_output)
    clout.setMultiOutput(false);
  addPoints2CommBC(superGeometry, material);
}


template<typename T, template<typename U> class Lattice>
void sOffLatticeBoundaryCondition3D<T,Lattice>::
defineU(SuperGeometry3D<T>& superGeometry, int material, AnalyticalF3D<T,T>& u, std::list<int> bulkMaterials )
{

  int nC = _sLattice.getLoadBalancer().size();
  for (int iCloc = 0; iCloc < nC; iCloc++) {
    _blockBCs[iCloc]->defineU(superGeometry.getExtendedBlockGeometry(iCloc), material, u, bulkMaterials );
  }
}


template<typename T, template<typename U> class Lattice>
void sOffLatticeBoundaryCondition3D<T,Lattice>::
addPoints2CommBC(SuperGeometry3D<T>& superGeometry, int material)
{

  if (_overlap != 0) {
    int nC = _sLattice.getLoadBalancer().size();
    for (int iCloc = 0; iCloc < nC; iCloc++) {
      int nX = superGeometry.getBlockGeometry(iCloc).getNx();
      int nY = superGeometry.getBlockGeometry(iCloc).getNy();
      int nZ = superGeometry.getBlockGeometry(iCloc).getNz();

      for (int iX=-_overlap; iX<nX+_overlap; iX++) {
        for (int iY=-_overlap; iY<nY+_overlap; iY++) {
          for (int iZ=-_overlap; iZ<nZ+_overlap; iZ++) {
            if (iX < 0 || iX > nX - 1 ||
                iY < 0 || iY > nY - 1 ||
                iZ < 0 || iZ > nZ - 1 ) {
              int found = false;
              if (superGeometry.getBlockGeometry(iCloc).getMaterial(iX,iY,iZ)!=0) {
                for (int iXo=-_overlap; iXo<=_overlap; iXo++) {
                  for (int iYo=-_overlap; iYo<=_overlap; iYo++) {
                    for (int iZo=-_overlap; iZo<=_overlap; iZo++) {
                      int nextX = iXo + iX;
                      int nextY = iYo + iY;
                      int nextZ = iZo + iZ;
                      if (superGeometry.getBlockGeometry(iCloc).getMaterial(nextX,nextY,nextZ)==material) {
                        _sLattice.get_commBC().add_cell(iCloc, iX, iY, iZ);
                        //std::cout << "found:" <<iX<<"/"<<iY<<"/"<<iZ<<std::endl;
                        found = true;
                      }
                      if (found) {
                        break;
                      }
                    }
                    if (found) {
                      break;
                    }
                  }
                  if (found) {
                    break;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}

template<typename T, template<typename U> class Lattice>
SuperLattice3D<T,Lattice>& sOffLatticeBoundaryCondition3D<T,Lattice>::getSuperLattice()
{
  return _sLattice;
}

template<typename T, template<typename U> class Lattice>
std::vector<OffLatticeBoundaryCondition3D<T,Lattice>* >& sOffLatticeBoundaryCondition3D<T,Lattice>::getBlockBCs()
{
  return _blockBCs;
}

template<typename T, template<typename U> class Lattice>
int sOffLatticeBoundaryCondition3D<T,Lattice>::getOverlap()
{
  return _overlap;
}

template<typename T, template<typename U> class Lattice>
void sOffLatticeBoundaryCondition3D<T,Lattice>::setOverlap(int overlap)
{
  _overlap = overlap;
}

template<typename T, template<typename U> class Lattice>
void sOffLatticeBoundaryCondition3D<T,Lattice>::outputOn()
{
  _output = true;
  int nC = _sLattice.getLoadBalancer().size();
  for (int iCloc = 0; iCloc < nC; iCloc++) {
    _blockBCs[iCloc]->outputOn();
  }
}

template<typename T, template<typename U> class Lattice>
void sOffLatticeBoundaryCondition3D<T,Lattice>::outputOff()
{
  _output = false;
  int nC = _sLattice.getLoadBalancer().size();
  for (int iCloc = 0; iCloc < nC; iCloc++) {
    _blockBCs[iCloc]->outputOff();
  }
}

////////////////// Factory functions //////////////////////////////////

template<typename T, template<typename U> class Lattice, typename MixinDynamics>
void createBouzidiBoundaryCondition3D(sOffLatticeBoundaryCondition3D<T,Lattice>& sBC)
{

  int nC = sBC.getSuperLattice().getLoadBalancer().size();
  sBC.setOverlap(1);
  for (int iC=0; iC<nC; iC++) {
    OffLatticeBoundaryCondition3D<T,Lattice>* blockBC
      = createBouzidiBoundaryCondition3D<T,Lattice,MixinDynamics>(sBC.getSuperLattice().getExtendedBlockLattice(iC));
    sBC.getBlockBCs().push_back(blockBC);
  }
}

}  // namespace olb

#endif

