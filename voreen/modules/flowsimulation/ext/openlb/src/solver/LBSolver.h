/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2012-2021 Mathias J. Krause, Benjamin Förster, 
 *  Julius Jessberger, Adrian Kummerlaender
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


#ifndef LBSOLVER_H
#define LBSOLVER_H

#include "SolverParameters.h"
#include "utilities/timer.h"
#include "utilities/typeIndexedContainers.h"
#include "utilities/typeMap.h"


namespace olb {

/// BaseSolver implements the solving process of an instationary simulation,
/// consisting of preSimulationTasks, time-stepping and postprocessing.
template <typename T, typename PARAMETERS>
class BaseSolver {

private:
  mutable OstreamManager                 clout {std::cout, "BaseSolver"};

public:
  using t = T;
  using Parameters_t = PARAMETERS;
  utilities::TypeIndexedTuple<PARAMETERS>           _parameters;

protected:
  bool                                   _isInitialized {false};
  std::size_t                            _iT {0};
  bool                                   _finishedTimeLoop {false};

public:
  BaseSolver(utilities::TypeIndexedTuple<PARAMETERS> params) : _parameters(params)
  { }


  void solve()
  {
    if (! _isInitialized) {
      initialize();
    }

    prepareSimulation();

    do
    {
      timeStep(_iT);
      ++_iT;
    } while (! exitCondition(_iT));
    
    _finishedTimeLoop = true;
    postprocessing();
  }

  /// Actions that shall be executed once after construction
  virtual void initialize() = 0;

protected:
  /// Actions that shall be executed before the time-stepping
  virtual void prepareSimulation() = 0;

  /// Defines what to do in a single time step
  virtual void timeStep(std::size_t iT) = 0;

  /// Condition, when to exit the time-stepping loop
  /// Returns true if the loop shall be continued
  virtual bool exitCondition(std::size_t iT) const = 0;

  /// Actions that shall be executed after the time-stepping
  virtual void postprocessing() = 0;

public:
  /// Access to parameter structs as parameters(KEY())
  template <typename KEY>
  auto& parameters(KEY = KEY()) {
    return *_parameters.template get<KEY>();
  }
  template <typename KEY>
  auto& parameters(KEY = KEY()) const {
    return *_parameters.template get<KEY>();
  }
  template <typename KEY>
  auto& parameters(meta::id<KEY>) {
    return *_parameters.template get<KEY>();
  }
  template <typename KEY>
  auto& parameters(meta::id<KEY>) const {
    return *_parameters.template get<KEY>();
  }
};


/// LBSolver is a generic solver for Lattice-Boltzmann problems. It holds
/// geometry and lattices and defines the abstract methods of BaseSolver.
/// Every simulation app should be able to inherit from this class.
template<
  typename T,
  typename PARAMETERS,
  typename LATTICES
>
class LBSolver : public BaseSolver<T,PARAMETERS> {

private:
  mutable OstreamManager                             clout {std::cout, "LBSolver"};

protected:
  static constexpr unsigned dim = LATTICES::values_t::template get<0>::d;

  template<typename... DESCRIPTORS>
  using SuperLattices = std::tuple<std::shared_ptr<SuperLattice<T,DESCRIPTORS>>...>;

  std::shared_ptr<SuperGeometry<T,dim>>              _sGeometry;
  std::shared_ptr<CuboidGeometry<T,dim>>             _cGeometry;
  std::shared_ptr<LoadBalancer<T>>                   _loadBalancer;

  typename LATTICES::values_t::template decompose_into<SuperLattices>  _sLattices;

  std::unique_ptr<util::Timer<BaseType<T>>>          _timer;

  static constexpr bool isStationary = PARAMETERS::keys_t::template contains<names::Stationarity>();
  static constexpr unsigned getNumberStationaryLattices() {
    if constexpr (isStationary) {
      return PARAMETERS::template value<names::Stationarity>::numberOfStationaryLattices;
    } else {
      return 0;
    }
    __builtin_unreachable();
  }
  std::array<std::unique_ptr<util::ValueTracer<T>>,
    getNumberStationaryLattices()> _convergenceCheck;

  bool                                               _exitMaxU {false};
  BaseType<T>                                        _boundMaxU {1.0};

public:
  LBSolver(utilities::TypeIndexedTuple<PARAMETERS> params) : LBSolver::BaseSolver(params)
  { }

  /// Build geometry, lattice and call computeResults
  // Allows to return the geometric/ lattice data without solving
  void buildAndReturn();

  /// Set up geometry
  void initialize() override;

protected:
  /// Set up lattice and initialize fields
  void prepareSimulation() override;

  /// Collide-and-stream + additional computations
  void timeStep(std::size_t iT) override;

  /// Evaluate results
  void postprocessing() override;


  /// Define the geometry
  // Inheritant has to provide this method
  virtual void prepareGeometry() = 0;

  /// Choose dynamics and boundary conditions
  // Inheritant has to provide this method
  virtual void prepareLattices() = 0;

  /// Define fields and initialize lattice populations
  // Inheritant has to provide this method
  // TODO: One may split this into fields and initial values/ population
  virtual void setInitialValues() = 0;

  /// Update fields and boundary values
  // Inheritant has to provide this method
  // TODO: One may split this into fields and boundary values (and 
  // postprocessors...)
  virtual void setBoundaryValues(std::size_t iT) = 0;

  /// Computation of results and output with full flexibility
  // Inheritant may provide this method
  virtual void getResults(std::size_t iT) { };

  /// Perform further computations (compute errors etc.)
  // Inheritant may provide these method
  // This function is called at every time step
  virtual void computeResults(std::size_t iT) { };
  // This function is called after the simulation
  virtual void computeResults() { };

  virtual bool exitCondition(std::size_t iT) const override;

  /// check stability: maxU should be <= _boundMaxU for a stable simulation
  /// Returns true if this fulfilled
  virtual bool checkStability(std::size_t iT);

private:
  /// Set up geometry and lattices
  void build();
  /// Construct and set up a new lattice
  void renewLattices();

  // ------------ Output generation -------------------------------------------
protected:
  virtual void printLog(std::size_t iT) const;

  virtual void writeLogConverter() const;

  /** Write geometric information for vtk output
   * The default version writes geometry, cuboid, rank and works for several
   * lattices. This method may be overridden in case e.g. discrete normals
   * shall be visualized (cf. example cavity2d).
   */
  virtual void prepareVTK() const;

  // Inheritant may override this method for vtk output
  virtual void writeVTK(std::size_t iT) const { };

  // Inheritant may override this method for image output
  virtual void writeImages(std::size_t iT) const { };

  // Inheritant may override this method for gnuplot output
  virtual void writeGnuplot() const { };

  // ------------ Access to converters and lattices ---------------------------
  template<typename... ARGS>
  auto& converter(ARGS&&... args) {
    using SimulationParameters_t = typename PARAMETERS::template value<names::Simulation>;
    if constexpr(std::is_invocable_v<decltype(&SimulationParameters_t::converter),ARGS...>) {
      return this->parameters(names::Simulation()).converter(args...);
    }
    else {
      return *this->parameters(names::Simulation()).converter;
    }
    __builtin_unreachable();
  }
  template<typename... ARGS>
  auto& converter(ARGS&&... args) const {
    using SimulationParameters_t = typename PARAMETERS::template value<names::Simulation>;
    if constexpr(std::is_invocable_v<decltype(&SimulationParameters_t::converter),ARGS...>) {
      return this->parameters(names::Simulation()).converter(args...);
    }
    else {
      return *this->parameters(names::Simulation()).converter;
    }
    __builtin_unreachable();
  }
  template <typename KEY>
  auto& lattice(KEY = KEY()) {
    return *std::get<(LATTICES::keys_t::template index<KEY>())>(_sLattices);
  }
  template <typename KEY>
  auto& lattice(KEY = KEY()) const {
    return *std::get<(LATTICES::keys_t::template index<KEY>())>(_sLattices);
  }
  template <typename KEY>
  auto& lattice(meta::id<KEY>) {
    return *std::get<(LATTICES::keys_t::template index<KEY>())>(_sLattices);
  }
  template <typename KEY>
  auto& lattice(meta::id<KEY>) const {
    return *std::get<(LATTICES::keys_t::template index<KEY>())>(_sLattices);
  }
  auto& lattice() {
    static_assert((LATTICES::size == 1), "Lattice name must be provided");
    return *std::get<0>(_sLattices);
  }
  auto& lattice() const {
    static_assert((LATTICES::size == 1), "Lattice name must be provided");
    return *std::get<0>(_sLattices);
  }
  auto& geometry() {
    return *_sGeometry;
  }
  auto& geometry() const {
    return *_sGeometry;
  }
  /*auto& results() {
    return *(this->_simulationResults);
  }
  auto& results() const {
    return *(this->_simulationResults);
  }*/
};


/// Returns a function that encapsulates the solving process.
template<typename T, typename SOLVER>
std::function<T (const std::vector<T>&)> getCallable(std::shared_ptr<SOLVER> solver){
  return [=](const std::vector<T>& control) -> T {
    solver->parameters(names::Opti()).applyControl(control);
    solver->solve();
    return solver->parameters(names::Results()).objective;
  };
}


// ------------------ Creator functions for XML interface ---------------------

template <class SOLVER>
std::shared_ptr<SOLVER> createLBSolver(XMLreader const& xml)
{
  utilities::TypeIndexedTuple<typename SOLVER::BaseSolver::Parameters_t> paramsTuple;
  meta::tuple_for_each(paramsTuple.tuple, [&](auto& params){ 
    using parameter_type = typename std::remove_reference_t<decltype(params)>::element_type;
    params = createParameters<parameter_type>(xml);
  });

  auto result = std::make_shared<SOLVER>(paramsTuple);
  return result;
}


} // namespace olb

#endif
