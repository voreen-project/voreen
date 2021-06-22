/*  This file is part of the OpenLB library
 *  Copyright (C) 2007 Mathias J. Krause, Jonathan Jeppener-Haltenhoff
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
 * The description of a 3D super lattice -- generic implementation.
 */

#ifndef SUPER_LATTICE_3D_HH
#define SUPER_LATTICE_3D_HH

#include <limits>
#include <numeric>

#include "communication/mpiManager.h"
#include "cell.h"
#include "superLattice3D.h"
#include "io/base64.h"
#include "functors/analytical/analyticalF.h"
#include "functors/lattice/superBaseF3D.h"
#include "functors/lattice/indicator/superIndicatorF3D.h"
#include "io/serializerIO.h"
#include "geometry/superGeometry3D.h"
#include "communication/loadBalancer.h"
#include "geometry/cuboidGeometry3D.h"

namespace olb {

////////////////////// Class SuperLattice3D /////////////////////////

template<typename T, typename DESCRIPTOR>
SuperLattice3D<T, DESCRIPTOR>::SuperLattice3D(SuperGeometry3D<T>& superGeometry)
  : SuperStructure3D<T>(superGeometry.getCuboidGeometry(), superGeometry.getLoadBalancer()),
    _commStream(*this),
    _commBC(*this),
    _statistics()
{
  auto& load = this->getLoadBalancer();

  int overlapBC = superGeometry.getOverlap();
  if (overlapBC >= 1) {
    _commBC_on = true;
    this->_overlap = overlapBC;
  } else {
    _commBC_on = false;
    this->_overlap = 1;
  }

  this->_communicator.init_nh();
  this->_communicator.add_cells(this->_overlap);
  this->_communicator.init();

  _extendedBlockLattices.reserve(load.size());

  for (int iC = 0; iC < load.size(); ++iC) {
    int nX = this->_cuboidGeometry.get(load.glob(iC)).getNx() + 2*this->_overlap;
    int nY = this->_cuboidGeometry.get(load.glob(iC)).getNy() + 2*this->_overlap;
    int nZ = this->_cuboidGeometry.get(load.glob(iC)).getNz() + 2*this->_overlap;
    _extendedBlockLattices.emplace_back(nX, nY, nZ);
  }

  for (int iC = 0; iC < load.size(); ++iC) {
    _blockLattices.emplace_back(
      _extendedBlockLattices[iC],
      this->_overlap, _extendedBlockLattices[iC].getNx() - this->_overlap - 1,
      this->_overlap, _extendedBlockLattices[iC].getNy() - this->_overlap - 1,
      this->_overlap, _extendedBlockLattices[iC].getNz() - this->_overlap - 1
    );
  }

#ifndef NEW_INTERIM_BLOCK_PROPAGATION
  _commStream.init_nh();
  _commStream.add_cells(1);
#endif
  _commStream.init();

  _statistics_on = true;

  if (_commBC_on) {
    _commBC.init_nh();
  }

  this->_communicationNeeded=true;
}

template<typename T, typename DESCRIPTOR>
Cell<T,DESCRIPTOR> SuperLattice3D<T,DESCRIPTOR>::get(int iC, int iX, int iY, int iZ)
{
#ifdef PARALLEL_MODE_MPI
  if (this->_loadBalancer.isLocal(iC)) {
    return _blockLattices[this->_loadBalancer.loc(iC)].get(iX, iY, iZ);
  } else {
    throw std::domain_error("iC must be local");
  }
#else
  return _blockLattices[this->_loadBalancer.loc(iC)].get(iX, iY, iZ);
#endif
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T,DESCRIPTOR>::setCopy(int iC, int iX, int iY, int iZ, ConstCell<T,DESCRIPTOR>& cell)
{
#ifdef PARALLEL_MODE_MPI
  if (this->_loadBalancer.isLocal(iC)) {
    _blockLattices[this->_loadBalancer.loc(iC)].get(iX, iY, iZ) = cell;
  }
#else
  _blockLattices[this->_loadBalancer.loc(iC)].get(iX, iY, iZ) = cell;
#endif
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T,DESCRIPTOR>::setCopy(Vector<int,4> pos, ConstCell<T,DESCRIPTOR>& cell)
{
  setCopy(pos[0], pos[1], pos[2], pos[3], cell);
}

template<typename T, typename DESCRIPTOR>
bool SuperLattice3D<T,DESCRIPTOR>::setCopy(T locX, T locY, T locZ, ConstCell<T,DESCRIPTOR>& cell)
{
  bool found = false;
  int iX, iY, iZ;
  for (int iC=0; iC<this->_loadBalancer.size(); ++iC) {
    if (this->_cuboidGeometry.get(this->_loadBalancer.glob(iC)).checkPoint(locX, locY, locZ, iX, iY, iZ, this->_overlap)) {
      _extendedBlockLattices[iC].get(iX,iY,iZ) = cell;
      found = true;
    }
  }
  return found;
}

template<typename T, typename DESCRIPTOR>
bool SuperLattice3D<T,DESCRIPTOR>::setCopy(Vector<T,3> pos, ConstCell<T,DESCRIPTOR>& cell)
{
  return setCopy(pos[0], pos[1], pos[2], cell);
}

template<typename T, typename DESCRIPTOR>
CellD<T,DESCRIPTOR> SuperLattice3D<T,DESCRIPTOR>::getCopy(int iC, int iX, int iY, int iZ) const
{
#ifdef PARALLEL_MODE_MPI
  const auto iCrank = this->_loadBalancer.rank(iC);
  if (this->_loadBalancer.isLocal(iC)) {
    auto cell  = _blockLattices[this->_loadBalancer.loc(iC)].get(iX, iY, iZ);

    int sizeOfCell = cell.getSerializedSize();
    singleton::mpi().bCast(&sizeOfCell, 1, iCrank);

    std::vector<T> cellData(sizeOfCell);
    _blockLattices[this->_loadBalancer.loc(iC)].get(iX, iY, iZ).serialize(cellData.data());

    singleton::mpi().bCast(cellData.data(), sizeOfCell, iCrank);

    return cell;
  } else {
    int sizeOfCell;
    singleton::mpi().bCast(&sizeOfCell, 1, iCrank);

    std::vector<T> cellData(sizeOfCell);

    singleton::mpi().bCast(cellData.data(), sizeOfCell, iCrank);

    CellD<T,DESCRIPTOR> cell;
    cell.unSerialize(cellData.data());
    return cell;
  }
#else
  auto cell = _blockLattices[this->_loadBalancer.loc(iC)].get(iX, iY, iZ);
  return CellD<T,DESCRIPTOR>(cell);
#endif
}

template<typename T, typename DESCRIPTOR>
CellD<T,DESCRIPTOR> SuperLattice3D<T,DESCRIPTOR>::getCopy(Vector<int,4> pos) const
{
  return getCopy(pos[0], pos[1], pos[2], pos[3]);
}

template<typename T, typename DESCRIPTOR>
CellD<T,DESCRIPTOR> SuperLattice3D<T,DESCRIPTOR>::getCopy(T locX, T locY, T locZ) const
{
  int iX, iY, iZ;
  //bool found = false;
  int iCfound = 0;

  for (int iC=0; iC<this->_cuboidGeometry.getNc(); ++iC) {
    if (this->_cuboidGeometry.get(iC).checkPoint(locX, locY, locZ, iX, iY, iZ)) {
      //found = true;
      iCfound = iC;
      break;
    }
  }

  return getCopy(iCfound, iX, iY, iZ);
}

template<typename T, typename DESCRIPTOR>
CellD<T,DESCRIPTOR> SuperLattice3D<T,DESCRIPTOR>::getCopy(Vector<T,3> pos) const
{
  return getCopy(pos[0], pos[1], pos[2]);
}


template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T, DESCRIPTOR>::initialize()
{
  if (_commBC_on) {
    _commBC.init();
  }

  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _blockLattices[iC].initialize();
  }

  this->_communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T,DESCRIPTOR>::defineDynamics(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicator, Dynamics<T, DESCRIPTOR>* dynamics)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].defineDynamics(
      indicator->getExtendedBlockIndicatorF(iC), dynamics);
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T, DESCRIPTOR>::defineDynamics(
  SuperGeometry3D<T>& sGeometry, int material, Dynamics<T, DESCRIPTOR>* dynamics)
{
  defineDynamics(sGeometry.getMaterialIndicator(material), dynamics);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T, DESCRIPTOR>::defineRho(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicator, AnalyticalF<3,T,T>& rho)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].defineRho(indicator->getExtendedBlockIndicatorF(iC),
                                         rho);
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T, DESCRIPTOR>::defineRho(
  SuperGeometry3D<T>& sGeometry, int material, AnalyticalF<3,T,T>& rho)
{
  defineRho(sGeometry.getMaterialIndicator(material), rho);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T, DESCRIPTOR>::defineU(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicator, AnalyticalF<3,T,T>& u)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].defineU(indicator->getExtendedBlockIndicatorF(iC),
                                       u);
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T, DESCRIPTOR>::defineU(
  SuperGeometry3D<T>& sGeometry, int material, AnalyticalF<3,T,T>& u)
{
  defineU(sGeometry.getMaterialIndicator(material), u);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T, DESCRIPTOR>::defineRhoU(FunctorPtr<SuperIndicatorF3D<T>>&& indicator,
    AnalyticalF<3,T,T>& rho, AnalyticalF<3,T,T>& u)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].defineRhoU(indicator->getExtendedBlockIndicatorF(iC),
                                          rho, u);
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T, DESCRIPTOR>::defineRhoU(SuperGeometry3D<T>& sGeometry, int material,
    AnalyticalF<3,T,T>& rho, AnalyticalF<3,T,T>& u)
{
  defineRhoU(sGeometry.getMaterialIndicator(material), rho, u);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T,DESCRIPTOR>::definePopulations(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicator, AnalyticalF<3,T,T>& Pop)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].definePopulations(
      indicator->getExtendedBlockIndicatorF(iC), Pop);
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T,DESCRIPTOR>::definePopulations(
  SuperGeometry3D<T>& sGeometry, int material, AnalyticalF<3,T,T>& Pop)
{
  definePopulations(sGeometry.getMaterialIndicator(material), Pop);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T,DESCRIPTOR>::definePopulations(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicator, SuperF3D<T,T>& Pop)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].definePopulations(
      indicator->getExtendedBlockIndicatorF(iC), Pop.getBlockF(iC));
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T,DESCRIPTOR>::definePopulations(
  SuperGeometry3D<T>& sGeometry, int material, SuperF3D<T,T>& Pop)
{
  definePopulations(sGeometry.getMaterialIndicator(material), Pop);
}


template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice3D<T,DESCRIPTOR>::defineField(
  SuperGeometry3D<T>& sGeometry, int material, SuperF3D<T,T>& field)
{
  if (sGeometry.getStatistics().getNvoxel(material)!=0) {
    int overlap = sGeometry.getOverlap();
    for (int iC=0; iC < this->_loadBalancer.size(); ++iC) {
      if (sGeometry.getExtendedBlockGeometry(iC).getStatistics().getNvoxel(material)!=0) {
        const int x0 = sGeometry.getExtendedBlockGeometry(iC).getStatistics().getMinLatticeR(material)[0];
        const int y0 = sGeometry.getExtendedBlockGeometry(iC).getStatistics().getMinLatticeR(material)[1];
        const int z0 = sGeometry.getExtendedBlockGeometry(iC).getStatistics().getMinLatticeR(material)[2];
        const int x1 = sGeometry.getExtendedBlockGeometry(iC).getStatistics().getMaxLatticeR(material)[0];
        const int y1 = sGeometry.getExtendedBlockGeometry(iC).getStatistics().getMaxLatticeR(material)[1];
        const int z1 = sGeometry.getExtendedBlockGeometry(iC).getStatistics().getMaxLatticeR(material)[2];
        for (int iX=x0; iX<=x1; ++iX) {
          for (int iY=y0; iY<=y1; ++iY) {
            for (int iZ=z0; iZ<=z1; ++iZ) {
              if (sGeometry.getExtendedBlockGeometry(iC).getMaterial(iX,iY,iZ) == material) {
                FieldD<T,DESCRIPTOR,FIELD> fieldTmp;
                int inputTmp[4]= {this->_loadBalancer.glob(iC),iX-overlap,iY-overlap,iZ-overlap};
                field(fieldTmp.data(), inputTmp);
                _extendedBlockLattices[iC].get(iX,iY,iZ).template setField<FIELD>(fieldTmp);
              }
            }
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice3D<T,DESCRIPTOR>::defineField(
  SuperGeometry3D<T>& sGeometry, IndicatorF3D<T>& indicator, int material, SuperF3D<T,T>& field)
{
  if (sGeometry.getStatistics().getNvoxel(material)!=0) {
    int overlap = sGeometry.getOverlap();
    for (int iC=0; iC < this->_loadBalancer.size(); ++iC) {
      if (sGeometry.getExtendedBlockGeometry(iC).getStatistics().getNvoxel(material)!=0) {
        const int x0 = sGeometry.getExtendedBlockGeometry(iC).getStatistics().getMinLatticeR(material)[0];
        const int y0 = sGeometry.getExtendedBlockGeometry(iC).getStatistics().getMinLatticeR(material)[1];
        const int z0 = sGeometry.getExtendedBlockGeometry(iC).getStatistics().getMinLatticeR(material)[2];
        const int x1 = sGeometry.getExtendedBlockGeometry(iC).getStatistics().getMaxLatticeR(material)[0];
        const int y1 = sGeometry.getExtendedBlockGeometry(iC).getStatistics().getMaxLatticeR(material)[1];
        const int z1 = sGeometry.getExtendedBlockGeometry(iC).getStatistics().getMaxLatticeR(material)[2];
        for (int iX=x0; iX<=x1; ++iX) {
          for (int iY=y0; iY<=y1; ++iY) {
            for (int iZ=z0; iZ<=z1; ++iZ) {
              if (sGeometry.getExtendedBlockGeometry(iC).getMaterial(
                    iX, iY, iZ) == material) {

                T physR[3];
                int inputTmp[4] = {this->_loadBalancer.glob(iC), iX - overlap,
                                   iY - overlap, iZ - overlap
                                  };

                sGeometry.getPhysR(physR, inputTmp);
                bool inside[1];
                indicator(inside, physR);

                if (inside[0]) {
                  FieldD<T,DESCRIPTOR,FIELD> fieldTmp;
                  field(fieldTmp.data(), inputTmp);
                  _extendedBlockLattices[iC].get(iX, iY, iZ).template setField<FIELD>(fieldTmp);
                }
              }
            }
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice3D<T,DESCRIPTOR>::defineField(SuperGeometry3D<T>& sGeometry, IndicatorF3D<T>& indicator,
    AnalyticalF<3,T,T>& field)
{
  SuperIndicatorFfromIndicatorF3D<T> indicatorF(indicator, sGeometry);
  defineField<FIELD>(indicatorF, field);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T,DESCRIPTOR>::iniEquilibrium(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicator,
  AnalyticalF<3,T,T>& rho, AnalyticalF<3,T,T>& u)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].iniEquilibrium(
      indicator->getExtendedBlockIndicatorF(iC), rho, u);
  }
  this->_communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T,DESCRIPTOR>::iniEquilibrium(
  SuperGeometry3D<T>& sGeometry, int material,
  AnalyticalF<3,T,T>& rho, AnalyticalF<3,T,T>& u)
{
  iniEquilibrium(sGeometry.getMaterialIndicator(material), rho, u);
}


template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T,DESCRIPTOR>::iniRegularized(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicator,
  AnalyticalF<3,T,T>& rho, AnalyticalF<3,T,T>& u,
  AnalyticalF<3,T,T>& pi)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].iniRegularized(
      indicator->getExtendedBlockIndicatorF(iC), rho, u, pi);
  }
  this->_communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T,DESCRIPTOR>::iniRegularized(
  SuperGeometry3D<T>& sGeometry, int material,
  AnalyticalF<3,T,T>& rho, AnalyticalF<3,T,T>& u,
  AnalyticalF<3,T,T>& pi)
{
  iniRegularized(sGeometry.getMaterialIndicator(material), rho, u, pi);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T, DESCRIPTOR>::collideAndStream()
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _blockLattices[iC].collide();
  }

#ifdef NEW_INTERIM_BLOCK_PROPAGATION
  _commStream.propagate();
#else
  _commStream.send();
  _commStream.receive();
  _commStream.write();
#endif

  if (_commBC_on) {
    _commBC.send();
    _commBC.receive();
    _commBC.write();
  }

  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].stream();
  }

  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _blockLattices[iC].postProcess();
  }

  if (_statistics_on) {
    reset_statistics();
  }

  this->_communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T,DESCRIPTOR>::stripeOffDensityOffset ( int x0, int x1, int y0,
    int y1, int z0, int z1, T offset )
{

  int locX0, locX1, locY0, locY1, locZ0, locZ1;
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    if (this->_cuboidGeometry.get(this->_loadBalancer.glob(iC)).checkInters(x0, x1, y0, y1,
        z0, z1, locX0, locX1, locY0, locY1, locZ0, locZ1, this->_overlap)) {
      _extendedBlockLattices[iC].stripeOffDensityOffset(locX0, locX1, locY0, locY1,
          locZ0, locZ1, offset);
    }
  }
  this->_communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T,DESCRIPTOR>::stripeOffDensityOffset(T offset)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].stripeOffDensityOffset(offset);
  }
  this->_communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
LatticeStatistics<T>& SuperLattice3D<T, DESCRIPTOR>::getStatistics()
{
  return _statistics;
}

template<typename T, typename DESCRIPTOR>
LatticeStatistics<T> const& SuperLattice3D<T, DESCRIPTOR>::getStatistics() const
{
  return _statistics;
}


template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T, DESCRIPTOR>::reset_statistics()
{
  T weight;
  T sum_weight = 0;
  T average_rho = 0;
  T average_energy = 0;
  T maxU = 0;
  T delta = 0;

  getStatistics().reset();

  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    delta = this->_cuboidGeometry.get(this->_loadBalancer.glob(iC)).getDeltaR();
    weight = _extendedBlockLattices[iC].getStatistics().getNumCells() * delta
             * delta * delta;
    sum_weight += weight;
    average_rho += _extendedBlockLattices[iC].getStatistics().getAverageRho()
                   * weight;
    average_energy += _extendedBlockLattices[iC].getStatistics().getAverageEnergy()
                      * weight;
    if (maxU < _extendedBlockLattices[iC].getStatistics().getMaxU()) {
      maxU = _extendedBlockLattices[iC].getStatistics().getMaxU();
    }
  }
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().reduceAndBcast(sum_weight, MPI_SUM);
  singleton::mpi().reduceAndBcast(average_rho, MPI_SUM);
  singleton::mpi().reduceAndBcast(average_energy, MPI_SUM);
  singleton::mpi().reduceAndBcast(maxU, MPI_MAX);
#endif

  average_rho = average_rho / sum_weight;
  average_energy = average_energy / sum_weight;

  getStatistics().reset(average_rho, average_energy, maxU, (int) sum_weight);
  getStatistics().incrementTime();
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    delta = this->_cuboidGeometry.get(this->_loadBalancer.glob(iC)).getDeltaR();
    _extendedBlockLattices[iC].getStatistics().reset(average_rho, average_energy,
        maxU, (int) sum_weight);
    _extendedBlockLattices[iC].getStatistics().incrementTime();
  }
}


template<typename T, typename DESCRIPTOR>
template<typename Slattice>
void SuperLattice3D<T, DESCRIPTOR>::addLatticeCoupling( LatticeCouplingGenerator3D<T, DESCRIPTOR> const& lcGen, SuperLattice3D<T,Slattice>& partnerLattice)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    std::vector< SpatiallyExtendedObject3D* > partnerOne;
    partnerOne.push_back(&partnerLattice.getExtendedBlockLattice(iC));

    int nx = getExtendedBlockLattice(iC).getNx();
    int ny = getExtendedBlockLattice(iC).getNy();
    int nz = getExtendedBlockLattice(iC).getNz();

    LatticeCouplingGenerator3D<T, DESCRIPTOR> *extractedLcGen = lcGen.clone();
    extractedLcGen->reset(1, nx-2, 1, ny-2, 1, nz-2);
    getExtendedBlockLattice(iC).addLatticeCoupling(*extractedLcGen, partnerOne);

    delete extractedLcGen;
  }
  return;
}


template<typename T, typename DESCRIPTOR>
template<typename Slattice>
void SuperLattice3D<T, DESCRIPTOR>::addLatticeCoupling(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicator,
  LatticeCouplingGenerator3D<T, DESCRIPTOR> const& lcGen,
  SuperLattice3D<T,Slattice>& partnerLattice)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    std::vector<SpatiallyExtendedObject3D*> partnerOne;
    partnerOne.push_back(&partnerLattice.getExtendedBlockLattice(iC));

    for (int iX = 1; iX < _extendedBlockLattices[iC].getNx()-1; ++iX) {
      for (int iY = 1; iY < _extendedBlockLattices[iC].getNy()-1; ++iY) {
        for (int iZ = 1; iZ < _extendedBlockLattices[iC].getNz()-1; ++iZ) {
          std::unique_ptr<LatticeCouplingGenerator3D<T, DESCRIPTOR>> extractedLcGen{
            lcGen.clone() };
          //TODO done quick and dirty
          if (extractedLcGen->extract(0, 0, 0, 0, 0, 0) ) {
            if (indicator->getExtendedBlockIndicatorF(iC)(iX, iY, iZ)) {
              extractedLcGen->shift(iX, iY, iZ, this->_loadBalancer.glob(iC));
              _extendedBlockLattices[iC].addLatticeCoupling(*extractedLcGen, partnerOne);
            }
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
template<typename Slattice>
void SuperLattice3D<T, DESCRIPTOR>::addLatticeCoupling(
  SuperGeometry3D<T>& sGeometry, int material,
  LatticeCouplingGenerator3D<T, DESCRIPTOR> const& lcGen,
  SuperLattice3D<T,Slattice>& partnerLattice)
{
  addLatticeCoupling(sGeometry.getMaterialIndicator(material),
                     lcGen, partnerLattice);
}

template<typename T, typename DESCRIPTOR>
template<typename Slattice>
void SuperLattice3D<T, DESCRIPTOR>::addLatticeCoupling(
  LatticeCouplingGenerator3D<T, DESCRIPTOR> const& lcGen,
  std::vector<SuperLattice3D<T,Slattice>*> partnerLattices)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    std::vector< SpatiallyExtendedObject3D* > partners;
    for (auto partnerLattice: partnerLattices) {
      partners.push_back(&partnerLattice->getExtendedBlockLattice(iC));
    }

    int nx = getExtendedBlockLattice(iC).getNx();
    int ny = getExtendedBlockLattice(iC).getNy();
    int nz = getExtendedBlockLattice(iC).getNz();

    LatticeCouplingGenerator3D<T, DESCRIPTOR> *extractedLcGen = lcGen.clone();
    extractedLcGen->reset(1, nx-2, 1, ny-2, 1, nz-2);
    getExtendedBlockLattice(iC).addLatticeCoupling(*extractedLcGen, partners);

    delete extractedLcGen;
  }
  return;
}


template<typename T, typename DESCRIPTOR>
template<typename Slattice>
void SuperLattice3D<T, DESCRIPTOR>::addLatticeCoupling(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicator,
  LatticeCouplingGenerator3D<T, DESCRIPTOR> const& lcGen,
  std::vector<SuperLattice3D<T,Slattice>*> partnerLattices)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    std::vector< SpatiallyExtendedObject3D* > partners;
    for (auto partnerLattice: partnerLattices) {
      partners.push_back(&partnerLattice->getExtendedBlockLattice(iC));
    }

    for (int iX = 1; iX < _extendedBlockLattices[iC].getNx()-1; ++iX) {
      for (int iY = 1; iY < _extendedBlockLattices[iC].getNy()-1; ++iY) {
        for (int iZ = 1; iZ < _extendedBlockLattices[iC].getNz()-1; ++iZ) {
          std::unique_ptr<LatticeCouplingGenerator3D<T, DESCRIPTOR>> extractedLcGen{
            lcGen.clone() };
          //TODO done quick and dirty
          if (extractedLcGen->extract(0, 0, 0, 0, 0, 0) ) {
            if (indicator->getExtendedBlockIndicatorF(iC)(iX, iY, iZ)) {
              extractedLcGen->shift(iX, iY, iZ, this->_loadBalancer.glob(iC));
              _extendedBlockLattices[iC].addLatticeCoupling(*extractedLcGen, partners);
            }
          }
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR>
template<typename Slattice>
void SuperLattice3D<T, DESCRIPTOR>::addLatticeCoupling(
  SuperGeometry3D<T>& sGeometry, int material,
  LatticeCouplingGenerator3D<T, DESCRIPTOR> const& lcGen,
  std::vector<SuperLattice3D<T,Slattice>*> partnerLattices)
{
  addLatticeCoupling(sGeometry.getMaterialIndicator(material),
                     lcGen, partnerLattices);
}


template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T, DESCRIPTOR>::addPostProcessor(
  PostProcessorGenerator3D<T, DESCRIPTOR> const& ppGen)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    int nx = getExtendedBlockLattice(iC).getNx();
    int ny = getExtendedBlockLattice(iC).getNy();
    int nz = getExtendedBlockLattice(iC).getNz();

    PostProcessorGenerator3D<T, DESCRIPTOR> *extractedPpGen = ppGen.clone();
    extractedPpGen->reset(1, nx-2, 1, ny-2, 1, nz-2);
    getExtendedBlockLattice(iC).addPostProcessor(*extractedPpGen);

    delete extractedPpGen;
  }
}


template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T, DESCRIPTOR>::addPostProcessor(
  FunctorPtr<SuperIndicatorF3D<T>>&& indicator,
  PostProcessorGenerator3D<T, DESCRIPTOR> const& ppGen)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    getExtendedBlockLattice(iC).addPostProcessor(ppGen);

    for (int iX = 1; iX < _extendedBlockLattices[iC].getNx()-1; ++iX) {
      for (int iY = 1; iY < _extendedBlockLattices[iC].getNy()-1; ++iY) {
        for (int iZ = 1; iZ < _extendedBlockLattices[iC].getNz()-1; ++iZ) {
          std::unique_ptr<PostProcessorGenerator3D<T, DESCRIPTOR>> extractedPpGen{
            ppGen.clone() };
          //TODO done quick and dirty
          if (extractedPpGen->extract(0, 0, 0, 0, 0, 0) ) {
            if (indicator->getExtendedBlockIndicatorF(iC)(iX, iY, iZ)) {
              extractedPpGen->shift(iX, iY, iZ, this->_loadBalancer.glob(iC));
              _extendedBlockLattices[iC].addPostProcessor(*extractedPpGen);
            }
          }
        }
      }
    }
  }

  this->_communicationNeeded = true;
}


template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T, DESCRIPTOR>::addPostProcessor(
  SuperGeometry3D<T>& sGeometry, int material,
  PostProcessorGenerator3D<T, DESCRIPTOR> const& ppGen)
{
  addPostProcessor(sGeometry.getMaterialIndicator(material), ppGen);
}


template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T, DESCRIPTOR>::executeCoupling()
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].executeCoupling();
  }
  this->_communicationNeeded = true;
  return;
}


template<typename T, typename DESCRIPTOR>
std::size_t SuperLattice3D<T,DESCRIPTOR>::getNblock() const
{
  return std::accumulate(_extendedBlockLattices.begin(), _extendedBlockLattices.end(), size_t(0),
                         Serializable::sumNblock());
}


template<typename T, typename DESCRIPTOR>
std::size_t SuperLattice3D<T,DESCRIPTOR>::getSerializableSize() const
{
  return std::accumulate(_extendedBlockLattices.begin(), _extendedBlockLattices.end(), size_t(0),
                         Serializable::sumSerializableSize());
}

template<typename T, typename DESCRIPTOR>
bool* SuperLattice3D<T,DESCRIPTOR>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  // NOTE: _extendedBlockLattices is correctly sized after constructing SuperLattice, so no resize should be done!
  for ( auto& bLattice : _extendedBlockLattices ) {
    registerSerializableOfConstSize(iBlock, sizeBlock, currentBlock, dataPtr, bLattice, loadingMode);
  }

  return dataPtr;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice3D<T,DESCRIPTOR>::postLoad()
{
  for (int iC=0; iC < this->_loadBalancer.size(); ++iC) {
    _extendedBlockLattices[iC].getStaticPopulationD().refreshControlStructure();
  }
}


////////// FREE FUNCTIONS //////////

template<typename T, typename DESCRIPTOR>
void setSuperExternalParticleField( SuperGeometry3D<T>& sGeometry, AnalyticalF<3,T,T>& velocity,
                                    SmoothIndicatorF3D<T,T,true>& sIndicator,
                                    SuperLattice3D<T,DESCRIPTOR>& sLattice)
{
  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    setBlockExternalParticleField( sGeometry.getExtendedBlockGeometry(iC), velocity, sIndicator,
                                   sLattice.getExtendedBlockLattice(iC));
  }
}

template<typename T, typename DESCRIPTOR>
void setSuperExternalParticleField( SuperGeometry3D<T>& sGeometry, AnalyticalF<3,T,T>& velocity,
                                    SmoothIndicatorF3D<T,T,true>& sIndicator,
                                    SuperLattice3D<T,DESCRIPTOR>& sLattice,
                                    Vector<bool,3> periodicity)
{
  Vector<T,3> min = sGeometry.getStatistics().getMinPhysR( 1 );
  Vector<T,3> max = sGeometry.getStatistics().getMaxPhysR( 1 );
  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    setBlockExternalParticleField( sGeometry.getExtendedBlockGeometry(iC), velocity, sIndicator,
                                   sLattice.getExtendedBlockLattice(iC), min, max, periodicity);
  }
}


} // namespace olb

#endif
