/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2007 Mathias J. Krause
 *                2021 Adrian Kummerlaender
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

#ifndef SUPER_LATTICE_HH
#define SUPER_LATTICE_HH

#include "superLattice.h"

#include "communication/mpiManager.h"
#include "cell.h"
#include "io/base64.h"

#include "functors/analytical/analyticalF.h"

#include "functors/lattice/superBaseF2D.h"
#include "functors/lattice/superBaseF3D.h"
#include "functors/lattice/indicator/superIndicatorF2D.h"
#include "functors/lattice/indicator/superIndicatorF3D.h"

#include "io/serializerIO.h"

#include "geometry/cuboidGeometry2D.h"
#include "geometry/cuboidGeometry3D.h"

#include "geometry/superGeometry.hh"

#include "communication/loadBalancer.h"

namespace olb {


template<typename T, typename DESCRIPTOR>
SuperLattice<T,DESCRIPTOR>::SuperLattice(SuperGeometry<T,DESCRIPTOR::d>& superGeometry)
  : SuperStructure<T,DESCRIPTOR::d>(superGeometry.getCuboidGeometry(),
                                    superGeometry.getLoadBalancer(),
                                    superGeometry.getOverlap()),
    _communicator{meta::make_array_f<std::unique_ptr<SuperCommunicator<T,SuperLattice>>,
                                     communication_stages::size>([&](unsigned) {
      return std::make_unique<SuperCommunicator<T,SuperLattice>>(*this);
    })},
    _statistics()
{
  auto& load = this->getLoadBalancer();

  for (int iC = 0; iC < load.size(); ++iC) {
    auto& cuboid = this->_cuboidGeometry.get(load.glob(iC));
    _block.emplace_back(constructUsingConcretePlatform<ConcretizableBlockLattice<T,DESCRIPTOR>>(
      load.platform(iC), cuboid.getExtent(), this->getOverlap()));
  }

  {
    auto& communicator = getCommunicator(PostCollide());
    communicator.template requestField<descriptors::POPULATION>();
    communicator.requestOverlap(1); // Required for inter-block propagation
    communicator.exchangeRequests();
  }

  {
    auto& communicator = getCommunicator(Full());
    DESCRIPTOR::fields_t::for_each([&](auto field) {
      communicator.template requestField<typename decltype(field)::type>();
    });
    // VTK output includes overlap of 1, some implicit dependencies for overlap of 2 exist
    communicator.requestOverlap(std::min(2, this->getOverlap()));
    communicator.exchangeRequests();
  }

  _statisticsEnabled = true;
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
Cell<T,DESCRIPTOR> SuperLattice<T,DESCRIPTOR>::get(LatticeR<DESCRIPTOR::d+1> latticeR)
{
#ifdef PARALLEL_MODE_MPI
  if (this->_loadBalancer.isLocal(latticeR[0])) {
    if constexpr (DESCRIPTOR::d == 3) {
      return _block[this->_loadBalancer.loc(latticeR[0])]->get(latticeR[1], latticeR[2], latticeR[3]);
    } else {
      return _block[this->_loadBalancer.loc(latticeR[0])]->get(latticeR[1], latticeR[2]);
    }
  }
  else {
    throw std::domain_error("Cuboid iC must be locally available");
  }
#else
  if constexpr (DESCRIPTOR::d == 3) {
    return _block[this->_loadBalancer.loc(latticeR[0])]->get(latticeR[1], latticeR[2], latticeR[3]);
  } else {
    return _block[this->_loadBalancer.loc(latticeR[0])]->get(latticeR[1], latticeR[2]);
  }
#endif
}

template<typename T, typename DESCRIPTOR>
template<typename... R>
std::enable_if_t<sizeof...(R) == DESCRIPTOR::d+1, Cell<T,DESCRIPTOR>>
SuperLattice<T,DESCRIPTOR>::get(R... latticeR)
{
  return get({latticeR...});
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::initialize()
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->initialize();
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
template<typename DYNAMICS>
void SuperLattice<T,DESCRIPTOR>::defineDynamics(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->template defineDynamics<DYNAMICS>(indicator->getBlockIndicatorF(iC));
  }
}

template<typename T, typename DESCRIPTOR>
template<typename DYNAMICS>
void SuperLattice<T,DESCRIPTOR>::defineDynamics(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material)
{
  defineDynamics<DYNAMICS>(sGeometry.getMaterialIndicator(material));
}

template<typename T, typename DESCRIPTOR>
template<template<typename...> typename DYNAMICS>
void SuperLattice<T,DESCRIPTOR>::defineDynamics(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->template defineDynamics<DYNAMICS<T,DESCRIPTOR>>(indicator->getBlockIndicatorF(iC));
  }
}

template<typename T, typename DESCRIPTOR>
template<template<typename...> typename DYNAMICS>
void SuperLattice<T,DESCRIPTOR>::defineDynamics(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material)
{
  defineDynamics<DYNAMICS>(sGeometry.getMaterialIndicator(material));
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::defineDynamics(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                                Dynamics<T,DESCRIPTOR>* dynamics)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->defineDynamics(indicator->getBlockIndicatorF(iC), dynamics);
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::defineDynamics(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                                                Dynamics<T,DESCRIPTOR>* dynamics)
{
  defineDynamics(sGeometry.getMaterialIndicator(material), dynamics);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::defineRho(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                           AnalyticalF<DESCRIPTOR::d,T,T>& rho)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->defineRho(indicator->getBlockIndicatorF(iC),
                         rho);
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::defineRho(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                                           AnalyticalF<DESCRIPTOR::d,T,T>& rho)
{
  defineRho(sGeometry.getMaterialIndicator(material), rho);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::defineU(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                         AnalyticalF<DESCRIPTOR::d,T,T>& u)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->defineU(indicator->getBlockIndicatorF(iC),
                        u);
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::defineU(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                                         AnalyticalF<DESCRIPTOR::d,T,T>& u)
{
  defineU(sGeometry.getMaterialIndicator(material), u);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::defineRhoU(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                            AnalyticalF<DESCRIPTOR::d,T,T>& rho,
                                            AnalyticalF<DESCRIPTOR::d,T,T>& u)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->defineRhoU(indicator->getBlockIndicatorF(iC),
                          rho, u);
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::defineRhoU(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                                            AnalyticalF<DESCRIPTOR::d,T,T>& rho,
                                            AnalyticalF<DESCRIPTOR::d,T,T>& u)
{
  defineRhoU(sGeometry.getMaterialIndicator(material), rho, u);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::definePopulations(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                                   AnalyticalF<DESCRIPTOR::d,T,T>& Pop)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->definePopulations(indicator->getBlockIndicatorF(iC), Pop);
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::definePopulations(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                                                   AnalyticalF<DESCRIPTOR::d,T,T>& Pop)
{
  definePopulations(sGeometry.getMaterialIndicator(material), Pop);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::definePopulations(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                                   SuperF<DESCRIPTOR::d,T,T>& Pop)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->definePopulations(indicator->getBlockIndicatorF(iC), Pop.getBlockF(iC));
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::definePopulations(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                                                   SuperF<DESCRIPTOR::d,T,T>& Pop)
{
  definePopulations(sGeometry.getMaterialIndicator(material), Pop);
}


template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice<T,DESCRIPTOR>::defineField(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                             FunctorPtr<SuperF<DESCRIPTOR::d,T,T>>&& field)
{
  if (field->getBlockFSize() == this->_loadBalancer.size()) {
    for (int iC=0; iC < this->_loadBalancer.size(); ++iC) {
      _block[iC]->template defineField<FIELD>(indicator->getBlockIndicatorF(iC),
                                              field->getBlockF(iC));
    }
  }
  else {
    FieldD<T,DESCRIPTOR,FIELD> fieldTmp;
    int coords[DESCRIPTOR::d+1];
    for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
      _block[iC]->forSpatialLocations([&](LatticeR<DESCRIPTOR::d> loc) {
        coords[0] = iC;
        for (int j = 0; j < DESCRIPTOR::d; ++j) {
          coords[j+1] = loc[j];
        }
        if (indicator(coords)) {
          field(fieldTmp.data(), coords);
          _block[iC]->get(loc).template setField<FIELD>(fieldTmp);
        }
      });
    }
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice<T,DESCRIPTOR>::defineField(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                             AnalyticalF<DESCRIPTOR::d,T,T>& field)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->template defineField<FIELD>(indicator->getBlockIndicatorF(iC),
                                            field);
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice<T,DESCRIPTOR>::defineField(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                                             FunctorPtr<SuperF<DESCRIPTOR::d,T,T>>&& field)
{
  defineField<FIELD>(sGeometry.getMaterialIndicator(material), std::forward<decltype(field)>(field));
}

template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice<T,DESCRIPTOR>::defineField(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                                             AnalyticalF<DESCRIPTOR::d,T,T>& field)

{
  defineField<FIELD>(sGeometry.getMaterialIndicator(material), field);
}

template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice<T,DESCRIPTOR>::defineField(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, IndicatorF3D<T>& indicator,
                                             AnalyticalF<DESCRIPTOR::d,T,T>& field)
{
  SuperIndicatorFfromIndicatorF3D<T> indicatorF(indicator, sGeometry);
  defineField<FIELD>(indicatorF, field);
}

template<typename T, typename DESCRIPTOR>
template <typename FIELD>
void SuperLattice<T,DESCRIPTOR>::setParameter(FieldD<T,DESCRIPTOR,FIELD>&& field)
{
  for (int iC=0; iC < this->getLoadBalancer().size(); ++iC) {
    _block[iC]->template setParameter<FIELD>(std::forward<decltype(field)>(field));
  }
}

template<typename T, typename DESCRIPTOR>
template <typename FIELD, typename DYNAMICS>
void SuperLattice<T,DESCRIPTOR>::setParameterOfDynamics(FieldD<T,DESCRIPTOR,FIELD>&& field)
{
  for (int iC=0; iC < this->getLoadBalancer().size(); ++iC) {
    _block[iC]->template getData<DynamicsParameters<DYNAMICS>>().template set<FIELD>(
       std::forward<decltype(field)>(field));
  }
}

template<typename T, typename DESCRIPTOR>
template <typename FIELD, template<typename...> typename DYNAMICS>
void SuperLattice<T,DESCRIPTOR>::setParameterOfDynamics(FieldD<T,DESCRIPTOR,FIELD>&& field)
{
  setParameterOfDynamics<FIELD,DYNAMICS<T,DESCRIPTOR>>(std::forward<decltype(field)>(field));
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::iniEquilibrium(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& rho,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& u)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->iniEquilibrium(indicator->getBlockIndicatorF(iC), rho, u);
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::iniEquilibrium(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& rho,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& u)
{
  iniEquilibrium(sGeometry.getMaterialIndicator(material), rho, u);
}


template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::iniRegularized(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& rho,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& u,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& pi)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->iniRegularized(indicator->getBlockIndicatorF(iC), rho, u, pi);
  }
  this->_communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::iniRegularized(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& rho,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& u,
                                                AnalyticalF<DESCRIPTOR::d,T,T>& pi)
{
  iniRegularized(sGeometry.getMaterialIndicator(material), rho, u, pi);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::collideAndStream()
{
  if (_statisticsEnabled) {
    setParameter<statistics::AVERAGE_RHO>(getStatistics().getAverageRho());
  }

  getCommunicator(PreCollide()).communicate();

  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->collide();
  }

  getCommunicator(PostCollide()).communicate();

  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->stream();
  }

  getCommunicator(PostStream()).communicate();

  executePostProcessors(PostStream());

  for (auto& f : _customPostProcessing) {
    f(*this);
  }

  getCommunicator(PostPostProcess()).communicate();

  if (_statisticsEnabled) {
    resetStatistics();
  }
  this->_communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
template<typename STAGE>
void SuperLattice<T,DESCRIPTOR>::executePostProcessors(STAGE stage)
{
  getCommunicator(stage).communicate();

  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->template postProcess<STAGE>();
  }
}

template<typename T, typename DESCRIPTOR>
template<typename STAGE>
SuperCommunicator<T,SuperLattice<T,DESCRIPTOR>>& SuperLattice<T,DESCRIPTOR>::getCommunicator(STAGE stage)
{
  if (communication_stages::template contains<STAGE>()) {
    return *_communicator[communication_stages::template index<STAGE>()];
  } else {
    auto iter = _customCommunicator.find(typeid(STAGE));
    if (iter == _customCommunicator.end()) {
      iter = std::get<0>(_customCommunicator.emplace(typeid(STAGE),
                                                     std::make_unique<SuperCommunicator<T,SuperLattice>>(*this)));
    }
    return *std::get<1>(*iter);
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::communicate()
{
  if (_communicationNeeded) {
    getCommunicator(Full()).communicate();
    _communicationNeeded = false;
  }
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::stripeOffDensityOffset(T offset)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->stripeOffDensityOffset(offset);
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
LatticeStatistics<T>& SuperLattice<T,DESCRIPTOR>::getStatistics()
{
  return _statistics;
}

template<typename T, typename DESCRIPTOR>
LatticeStatistics<T> const& SuperLattice<T,DESCRIPTOR>::getStatistics() const
{
  return _statistics;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::resetStatistics()
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
    weight = _block[iC]->getStatistics().getNumCells() * delta
             * delta * delta;
    sum_weight += weight;
    average_rho += _block[iC]->getStatistics().getAverageRho()
                   * weight;
    average_energy += _block[iC]->getStatistics().getAverageEnergy()
                      * weight;
    if (maxU < _block[iC]->getStatistics().getMaxU()) {
      maxU = _block[iC]->getStatistics().getMaxU();
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
    _block[iC]->getStatistics().reset(average_rho, average_energy,
        maxU, (int) sum_weight);
    _block[iC]->getStatistics().incrementTime();
  }
}

template<typename T, typename DESCRIPTOR>
template<typename STAGE>
void SuperLattice<T,DESCRIPTOR>::addPostProcessor(
  PostProcessorGenerator<T,DESCRIPTOR> const& ppGen)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    int overlap = getBlock(iC).getPadding();
    PostProcessorGenerator<T,DESCRIPTOR> *extractedPpGen = ppGen.clone();
    extractedPpGen->reset(-overlap+1, getBlock(iC).getExtent()+overlap-2);
    getBlock(iC).template addPostProcessor<STAGE>(*extractedPpGen);
    delete extractedPpGen;
  }
}

template<typename T, typename DESCRIPTOR>
template<typename STAGE>
void SuperLattice<T,DESCRIPTOR>::addPostProcessor(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                                  PostProcessorGenerator<T,DESCRIPTOR> const& ppGen)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    getBlock(iC).template addPostProcessor<STAGE>(indicator->getBlockIndicatorF(iC),
                                                  ppGen);
  }
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
template<typename STAGE>
void SuperLattice<T,DESCRIPTOR>::addPostProcessor(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                                  PostProcessorPromise<T,DESCRIPTOR>&& promise)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    getBlock(iC).addPostProcessor(typeid(STAGE),
                                  indicator->getBlockIndicatorF(iC),
                                  std::forward<decltype(promise)>(promise));
  }
}

template<typename T, typename DESCRIPTOR>
template<typename STAGE>
void SuperLattice<T,DESCRIPTOR>::addPostProcessor(PostProcessorPromise<T,DESCRIPTOR>&& promise)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    getBlock(iC).addPostProcessor(typeid(STAGE),
                                  std::forward<decltype(promise)>(promise));
  }
}

template<typename T, typename DESCRIPTOR>
template<typename STAGE>
void SuperLattice<T,DESCRIPTOR>::addPostProcessor(
  SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
  PostProcessorGenerator<T,DESCRIPTOR> const& ppGen)
{
  addPostProcessor<STAGE>(sGeometry.getMaterialIndicator(material), ppGen);
}

template<typename T, typename DESCRIPTOR>
template<typename PARTNER_DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::addLatticeCoupling(
  LatticeCouplingGenerator<T,DESCRIPTOR> const& lcGen,
  std::vector<SuperLattice<T,PARTNER_DESCRIPTOR>*> partnerLattices)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    std::vector<BlockStructureD<DESCRIPTOR::d>*> partners;
    for (auto partnerLattice: partnerLattices) {
      partners.push_back(&partnerLattice->getBlock(iC));
    }
    int overlap = getBlock(iC).getPadding();
    LatticeCouplingGenerator<T,DESCRIPTOR> *extractedLcGen = lcGen.clone();
    extractedLcGen->reset(-overlap+1, getBlock(iC).getExtent()+overlap-2);
    getBlock(iC).addLatticeCoupling(*extractedLcGen, partners);
    delete extractedLcGen;
  }
  return;
}

template<typename T, typename DESCRIPTOR>
template<typename PARTNER_DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::addLatticeCoupling(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                                    LatticeCouplingGenerator<T,DESCRIPTOR> const& lcGen,
                                                    std::vector<SuperLattice<T,PARTNER_DESCRIPTOR>*> partnerLattices)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    std::vector<BlockStructureD<DESCRIPTOR::d>*> partners;
    for (auto partnerLattice: partnerLattices) {
      partners.push_back(&partnerLattice->getBlock(iC));
    }
    getBlock(iC).addLatticeCoupling(indicator->getBlockIndicatorF(iC),
                                           lcGen,
                                           partners);
  }
}

template<typename T, typename DESCRIPTOR>
template<typename PARTNER_DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::addLatticeCoupling(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                                                    LatticeCouplingGenerator<T,DESCRIPTOR> const& lcGen,
                                                    std::vector<SuperLattice<T,PARTNER_DESCRIPTOR>*> partnerLattices)
{
  addLatticeCoupling(sGeometry.getMaterialIndicator(material),
                     lcGen, partnerLattices);
}

template<typename T, typename DESCRIPTOR>
template<typename... PARTNER_DESCRIPTORS>
void SuperLattice<T,DESCRIPTOR>::addLatticeCoupling(
  LatticeCouplingGenerator<T,DESCRIPTOR> const& lcGen,
  PARTNER_DESCRIPTORS&... partnerLattices)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    std::vector<BlockStructureD<DESCRIPTOR::d>*> partners{&partnerLattices.getBlock(iC)...};
    int overlap = getBlock(iC).getPadding();
    LatticeCouplingGenerator<T,DESCRIPTOR> *extractedLcGen = lcGen.clone();
    extractedLcGen->reset(-overlap+1, getBlock(iC).getExtent()+overlap-2);
    getBlock(iC).addLatticeCoupling(*extractedLcGen, partners);
    delete extractedLcGen;
  }
  return;
}

template<typename T, typename DESCRIPTOR>
template<typename... PARTNER_DESCRIPTORS>
void SuperLattice<T,DESCRIPTOR>::addLatticeCoupling(FunctorPtr<SuperIndicatorF<T,DESCRIPTOR::d>>&& indicator,
                                                    LatticeCouplingGenerator<T,DESCRIPTOR> const& lcGen,
                                                    PARTNER_DESCRIPTORS&... partnerLattices)
{
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    std::vector<BlockStructureD<DESCRIPTOR::d>*> partners{&partnerLattices.getBlock(iC)...};
    getBlock(iC).addLatticeCoupling(indicator->getBlockIndicatorF(iC),
                                    lcGen,
                                    partners);
  }
}

template<typename T, typename DESCRIPTOR>
template<typename... PARTNER_DESCRIPTORS>
void SuperLattice<T,DESCRIPTOR>::addLatticeCoupling(SuperGeometry<T,DESCRIPTOR::d>& sGeometry, int material,
                                                    LatticeCouplingGenerator<T,DESCRIPTOR> const& lcGen,
                                                    PARTNER_DESCRIPTORS&... partnerLattices)
{
  addLatticeCoupling(sGeometry.getMaterialIndicator(material),
                     lcGen, partnerLattices...);
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::executeCoupling()
{
  getCommunicator(PreCoupling()).communicate();
  for (int iC = 0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->executeCoupling();
  }
  getCommunicator(PostCoupling()).communicate();
  _communicationNeeded = true;
}

template<typename T, typename DESCRIPTOR>
std::size_t SuperLattice<T,DESCRIPTOR>::getNblock() const
{
  return std::accumulate(_block.begin(), _block.end(), size_t(0), [](std::size_t sum, auto& b) -> std::size_t {
    return sum + b->getNblock();
  });
}


template<typename T, typename DESCRIPTOR>
std::size_t SuperLattice<T,DESCRIPTOR>::getSerializableSize() const
{
  return std::accumulate(_block.begin(), _block.end(), size_t(0), [](std::size_t sum, auto& b) -> std::size_t {
    return sum + b->getSerializableSize();
  });
}

template<typename T, typename DESCRIPTOR>
bool* SuperLattice<T,DESCRIPTOR>::getBlock(std::size_t iBlock, std::size_t& sizeBlock, bool loadingMode)
{
  std::size_t currentBlock = 0;
  bool* dataPtr = nullptr;

  for (std::size_t iC=0; iC < _block.size(); ++iC) {
    registerSerializableOfConstSize(iBlock, sizeBlock, currentBlock, dataPtr, getBlock(iC), loadingMode);
  }

  return dataPtr;
}

template<typename T, typename DESCRIPTOR>
void SuperLattice<T,DESCRIPTOR>::postLoad()
{
  for (int iC=0; iC < this->_loadBalancer.size(); ++iC) {
    _block[iC]->postLoad();
  }
}


}

#endif
