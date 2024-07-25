/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2021 Mathias J. Krause, Benjamin Förster,
 *  Julius Jessberger
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


#ifndef LBSOLVER_HH
#define LBSOLVER_HH


#include "LBSolver.h"


namespace olb {

template<typename T,  typename PARAMETERS, typename LATTICES>
void LBSolver<T,PARAMETERS,LATTICES>::initialize()
{
  if (this->parameters(names::Output()).verbose) {
    clout << "Start initialization" << std::endl;
  }

  this->prepareGeometry();
  // Removes all not needed boundary voxels outside the surface
  //this->geometry().clean(this->parameters(names::Output()).verbose);
  // Removes all not needed boundary voxels inside the surface
  this->geometry().innerClean(this->parameters(names::Output()).verbose);
  geometry().checkForErrors(this->parameters(names::Output()).verbose);
  if (this->parameters(names::Output()).verbose) {
    // TODO: enable different output of the geometry
    // this->geometry().print();
    geometry().getStatistics().print();
  }

  this->_isInitialized = true;

  if (this->parameters(names::Output()).verbose) {
    clout << "Finished initialization" << std::endl;
  }
}

template<typename T,  typename PARAMETERS, typename LATTICES>
void LBSolver<T,PARAMETERS,LATTICES>::buildAndReturn()
{
  if (! this->_isInitialized) {
    this->initialize();
  }
  build();
  computeResults();
}

template<typename T,  typename PARAMETERS, typename LATTICES>
void LBSolver<T,PARAMETERS,LATTICES>::build()
{
  if (this->parameters(names::Output()).printLogConverter) {
    writeLogConverter();
  }

  _timer = std::make_unique<util::Timer<BaseType<T>>>(
    converter().getLatticeTime(this->parameters(names::Simulation()).maxTime),
    _sGeometry->getStatistics().getNvoxel());

  if constexpr (isStationary){
    // for (auto& tracer : _convergenceCheck) {
    for (unsigned i = 0; i < getNumberStationaryLattices(); ++i) {
      // tracer = std::make_unique<util::ValueTracer<T>>(
      _convergenceCheck[i] = std::make_unique<util::ValueTracer<T>>(
        converter().getLatticeTime(this->parameters(names::Stationarity()).physInterval[i]),
        this->parameters(names::Stationarity()).epsilon[i]);
    }
  }

  renewLattices();
}

template<typename T,  typename PARAMETERS, typename LATTICES>
void LBSolver<T,PARAMETERS,LATTICES>::prepareSimulation()
{
  if (this->parameters(names::Output()).verbose) {
    clout << "Start prepareSimulation" << std::endl;
  }

  this->_iT = 0;

  build();

  setInitialValues();

  meta::tuple_for_each(_sLattices, [](auto& lattice){
    lattice->initialize();
    lattice->getStatistics().reset();
    lattice->getStatistics().initialize();
  });

  if (this->parameters(names::Output()).vtkOutput) {
    prepareVTK();
  }

  if (this->parameters(names::Output()).verbose) {
    clout << "Finished prepareSimulation" << std::endl;
    clout << "Start collide-and-stream loop" << std::endl;
  }
  _timer->start();
}

template<typename T,  typename PARAMETERS, typename LATTICES>
void LBSolver<T,PARAMETERS,LATTICES>::timeStep(std::size_t iT)
{
  const std::size_t itBoundaryUpdate = util::max(
    converter().getLatticeTime(this->parameters(names::Simulation()).physBoundaryValueUpdateTime),
    (unsigned long) {1});
  if (iT % itBoundaryUpdate == 0) {
    setBoundaryValues(iT);
  }

  meta::tuple_for_each(_sLattices, [](auto& lattice) {
    // TODO: dispose also communication of fields...
    lattice->communicate();
  });

  if ( this->parameters(names::Output()).vtkOutput
    && iT % converter().getLatticeTime(this->parameters(names::Output()).vtkSaveT) == 0 ) {
    writeVTK(iT);
  }
  if ( this->parameters(names::Output()).imageOutput
    && iT % converter().getLatticeTime(this->parameters(names::Output()).imageSaveT) == 0 ) {
    writeImages(iT);
  }
  if ( this->parameters(names::Output()).verbose
    && iT % converter().getLatticeTime(this->parameters(names::Output()).logT) == 0) {
    printLog(iT);
  }
  computeResults(iT);
  getResults(iT);

  meta::tuple_for_each(_sLattices, [](auto& lattice){
    lattice->collideAndStream();
    lattice->executeCoupling();
      });


  if (this->parameters(names::Simulation()).pressureFilter) {
    // TODO: How to choose the correct lattice here?
    std::get<0>(_sLattices)->stripeOffDensityOffset(std::get<0>(_sLattices)->getStatistics().getAverageRho() - T(1));
  }

  if constexpr (isStationary){
    unsigned counter = 0;
    using parameters_t = typename PARAMETERS::template value<names::Stationarity>;
    parameters_t::stat_lattices::for_each( [&](auto Id) {
      // using Name = typename decltype(Id)::type;
      if (this->parameters(names::Stationarity()).convergenceType[counter] == parameters_t::MaxLatticeVelocity){
        _convergenceCheck[counter]->takeValue( this->lattice(Id).getStatistics().getMaxU());
      }
      else if (this->parameters(names::Stationarity()).convergenceType[counter] == parameters_t::AverageEnergy){
        _convergenceCheck[counter]->takeValue( this->lattice(Id).getStatistics().getAverageEnergy() );
      }
      else if (this->parameters(names::Stationarity()).convergenceType[counter] == parameters_t::AverageRho){
        _convergenceCheck[counter]->takeValue( this->lattice(Id).getStatistics().getAverageRho());
      }
      else {
        throw std::invalid_argument("Convergence type is not supported.\n");
      }
      counter++;
    });
  }

  checkStability(iT);
}

template<typename T,  typename PARAMETERS, typename LATTICES>
void LBSolver<T,PARAMETERS,LATTICES>::postprocessing()
{
  meta::tuple_for_each(_sLattices, [](auto& lattice) {
    // TODO: dispose also communication of fields...
    lattice->communicate();
  });

  if ( this->parameters(names::Output()).verbose ) {
    printLog(this->_iT);
  }
  if ( this->parameters(names::Output()).vtkOutput ) {
    writeVTK(this->_iT);
  }
  if ( this->parameters(names::Output()).imageOutput ) {
    writeImages(this->_iT);
  }
  getResults(this->_iT);

  _timer->stop();
  if (this->parameters(names::Output()).verbose) {
    clout << "Finished collide-and-stream loop" << std::endl;
    clout << "Start postprocessing" << std::endl;
  }
  computeResults();

  _timer->printShortSummary();

  if ( this->parameters(names::Output()).gnuplotOutput ) {
    writeGnuplot();
  }

  if (this->parameters(names::Output()).verbose) {
    clout << "Finished postprocessing" << std::endl;
  }
}

template<typename T,  typename PARAMETERS, typename LATTICES>
bool LBSolver<T,PARAMETERS,LATTICES>::exitCondition(std::size_t iT) const
{
  if (iT > this->converter().getLatticeTime(this->parameters(names::Simulation()).maxTime)){
    return true;
  }

  if constexpr (isStationary){
    if (std::all_of(_convergenceCheck.cbegin(), _convergenceCheck.cend(), [](auto& c){
      return c->hasConverged();
      } )) {
      clout << "Simulation converged." << std::endl;
      return true;
    }
  }

  return false;
}

template<typename T,  typename PARAMETERS, typename LATTICES>
bool LBSolver<T,PARAMETERS,LATTICES>::checkStability(std::size_t iT)
{
  bool result = true;
  meta::tuple_for_each(_sLattices, [&](auto& lattice){
    if (lattice->getStatistics().getMaxU() > _boundMaxU) {
      clout << "PROBLEM uMax=" << lattice->getStatistics().getMaxU() << std::endl;
      lattice->getStatistics().print(iT, converter().getPhysTime(iT));
      if (_exitMaxU) {
        exit(1);
      }
      result = false;
    }
  });
  return result;
}

template<typename T,  typename PARAMETERS, typename LATTICES>
void LBSolver<T,PARAMETERS,LATTICES>::renewLattices()
{
  meta::tuple_for_each(_sLattices, [&](auto& lattice){
    using lattice_type = typename std::remove_reference_t<decltype(lattice)>::element_type;
    lattice = std::make_shared<lattice_type>(geometry());
  });

  prepareLattices();
}

template<typename T,  typename PARAMETERS, typename LATTICES>
void LBSolver<T,PARAMETERS,LATTICES>::printLog(std::size_t iT) const
{
  using OutputParameters_t = typename PARAMETERS::template value<names::Output>;
  OutputParameters_t::printLatticeStatistics::for_each( [&](auto type){
    lattice(type.get()).getStatistics().print(iT, converter().getPhysTime(iT));
  } );

  _timer->print(iT, this->parameters(names::Output()).timerPrintMode);
}

template<typename T,  typename PARAMETERS, typename LATTICES>
void LBSolver<T,PARAMETERS,LATTICES>::writeLogConverter() const
{
  converter().print();
  converter().write(this->parameters(names::Output()).name);
}

template<typename T,  typename PARAMETERS, typename LATTICES>
void LBSolver<T,PARAMETERS,LATTICES>::prepareVTK() const
{
  // write the geometric information. Since this does not depend on the
  // lattice, we may w.l.o.g. work with the first lattice
  auto& lattice = *std::get<0>(_sLattices);
  using DESCRIPTOR = typename LATTICES::values_t::template get<0>;

  SuperVTMwriter<T,dim> vtmWriter(this->parameters(names::Output()).vtkFilename);
  /// Writes the geometry, cuboid no. and rank no. as vti file for visualization
  SuperLatticeGeometry<T,DESCRIPTOR> geometry(lattice, this->geometry());
  SuperLatticeCuboid<T,DESCRIPTOR> cuboid(lattice);
  SuperLatticeRank<T,DESCRIPTOR> rank(lattice);
  vtmWriter.write( cuboid );
  vtmWriter.write( geometry );
  vtmWriter.write( rank );
  vtmWriter.createMasterFile();
}

} // namespace olb

#endif
