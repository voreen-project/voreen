/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2008 Orestis Malaspinas
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


/* porousPlate3dSolver:
 * This example is identical to porousPlate3d, but it is implemented on the
 * base of the Solver-class which holds the standard variables (like the
 * super geometry) and executes the standard steps that are necessary for
 * simulation.
 */


#include "olb3D.h"
#include "olb3D.hh"

using namespace olb;
using namespace olb::names;

using NSDESCRIPTOR = descriptors::D3Q19<descriptors::FORCE>;
using TDESCRIPTOR = descriptors::D3Q7<descriptors::VELOCITY>;
// map lattice names to the descriptors
using LATTICES = meta::map<
  NavierStokes, NSDESCRIPTOR,
  Temperature, TDESCRIPTOR
>;

/// Struct that keeps all simulation-specific parameters
template<typename T>
struct PorousPlate3dSimulationParameters : public parameters::SimulationBase<T>
{
  const T lx = 1.0;      // length of the channel
  const T ly = 1.0;      // height of the channel
  T lz = 1.0;
  const int N = 20;      // resolution of the model
  const T tau = 1.;      // relaxation time
  const T Re = 5.;       // Reynolds number
  const T Ra = 100.;     // Rayleigh number
  const T Pr = 0.71;     // Prandtl number

  const T Tcold = 273.15;
  const T Thot = 274.15;

  std::shared_ptr<ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const> converter;

  PorousPlate3dSimulationParameters()
  {
    this->converter = std::make_shared<ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const> (
      (T) 1.0 / N, // physDeltaX
      (T) 1.0 / N * 1.0 / 1e-3 * (tau - 0.5) / 3 / N, // physDeltaT
      (T) 1.0, // charPhysLength
      (T) util::sqrt( 9.81 * Ra * 1e-3 * 1e-3 / Pr / 9.81 / (Thot - Tcold) / util::pow(1.0, 3) * (Thot - Tcold) * 1.0 ), // charPhysVelocity
      (T) 1e-3, // physViscosity
      (T) 1.0, // physDensity
      (T) 0.03, // physThermalConductivity
      (T) Pr * 0.03 / 1e-3 / 1.0, // physSpecificHeatCapacity
      (T) Ra * 1e-3 * 1e-3 / Pr / 9.81 / (Thot - Tcold) / util::pow(1.0, 3), // physThermalExpansionCoefficient
      (T) Tcold, // charPhysLowTemperature
      (T) Thot // charPhysHighTemperature
    );

    lz = converter->getPhysDeltaX() * 3.;
    this->maxTime = 1e4;
  }
};

/// Map names to parameter structs.
// Enables access to instances of parameter structs by names.
// Default versions are used for Stationarity and Output.
template<typename T>
using Parameters = meta::map<
  Simulation, PorousPlate3dSimulationParameters<T>,
  Stationarity, parameters::Stationarity<T,Temperature>,
  Output, parameters::OutputBase<T>
>;

/// Implement analytical solution (velocity) as a functor
template <typename T, typename S>
class AnalyticalVelocityPorousPlate3D : public AnalyticalF3D<T, S> {
private:
  T _Re;
  T _u0;
  T _v0;
  T _ly;
public:
  AnalyticalVelocityPorousPlate3D(T Re, T u0, T v0, T ly) : AnalyticalF3D<T, S>(3),
    _Re(Re), _u0(u0), _v0(v0), _ly(ly) {
    this->getName() = "AnalyticalVelocityPorousPlate3D";
  }

  bool operator()(T output[3], const S x[3]) override {
    output[0] = _u0*((util::exp(_Re* x[1] / _ly) - 1) / (util::exp(_Re) - 1));
    output[1] = _v0;
    output[2] = 0.0;
    return true;
  }
};

/// Implement analytical solution (temperature) as a functor
template <typename T, typename S>
class AnalyticalTemperaturePorousPlate3D : public AnalyticalF3D<T, S> {
private:
  T _Re;
  T _Pr;
  T _ly;
  T _T0;
  T _deltaT;
public:
  AnalyticalTemperaturePorousPlate3D(T Re, T Pr, T ly, T T0, T deltaT) : AnalyticalF3D<T, S>(1),
    _Re(Re), _Pr(Pr), _ly(ly), _T0(T0), _deltaT(deltaT) {
    this->getName() = "AnalyticalTemperaturePorousPlate3D";
  }

  bool operator()(T output[1], const S x[3]) override {
    output[0] = _T0 + _deltaT*((util::exp(_Pr*_Re*x[1] / _ly) - 1) / (util::exp(_Pr*_Re) - 1));
    return true;
  }
};

/// Implement analytical solution (heat flux) as a functor
template <typename T, typename S>
class AnalyticalHeatFluxPorousPlate3D : public AnalyticalF3D<T, S> {
private:
  T _Re;
  T _Pr;
  T _deltaT;
  T _ly;
  T _lambda;
public:
  AnalyticalHeatFluxPorousPlate3D(T Re, T Pr, T deltaT, T ly,T lambda) : AnalyticalF3D<T, S>(3),
    _Re(Re), _Pr(Pr), _deltaT(deltaT), _ly(ly), _lambda(lambda) {
    this->getName() = "AnalyticalHeatFluxPorousPlate3D";
  }

  bool operator()(T output[3], const S x[3]) override {
    output[0] = 0;
    output[1] = - _lambda * _Re * _Pr * _deltaT / _ly * (util::exp(_Pr * _Re * x[1] / _ly))/(util::exp(_Pr * _Re) - 1);
    output[2] = 0;
    return true;
  }
};


/// Solver class: implements all simulation-specific functions and routines
template<typename T>
class PorousPlate3dSolver : public LBSolver<T,Parameters<T>,LATTICES> {
private:
  mutable OstreamManager            clout {std::cout, "PorousPlate3dSolver"};

public:
  PorousPlate3dSolver(utilities::TypeIndexedTuple<Parameters<T>> params)
   : PorousPlate3dSolver::LBSolver(params)
  { }

protected:
  void prepareGeometry() override
  {
    std::vector<T> extend(3,T());
    extend[0] = this->parameters(Simulation()).lx;
    extend[1] = this->parameters(Simulation()).ly;
    extend[2] = this->parameters(Simulation()).lz;
    const std::vector<T> origin(3,T());
    IndicatorCuboid3D<T> cuboid(extend, origin);

    /// Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
    const unsigned noOfCuboids = this->parameters(Simulation()).noC * singleton::mpi().getSize();
#else
    const unsigned noOfCuboids = 7;
#endif
    this->_cGeometry = std::make_shared<CuboidGeometry3D<T>>(
      cuboid,
      this->converter().getPhysDeltaX(),
      noOfCuboids);
    this->_cGeometry->setPeriodicity(true,false, true);

    /// Instantiation of a loadBalancer
    this->_loadBalancer = std::make_shared<HeuristicLoadBalancer<T>> (*this->_cGeometry);

    /// Instantiation of a superGeometry
    this->_sGeometry = std::make_shared<SuperGeometry<T,3>> (
      *this->_cGeometry,
      *this->_loadBalancer,
      this->parameters(Simulation()).overlap);

    this->geometry().rename(0,2);
    this->geometry().rename(2,1,{0,1,0});

    std::vector<T> extendBottom( 3, T(0) );
    extendBottom[0] = this->parameters(Simulation()).lx;
    extendBottom[1] = this->converter().getPhysLength(1);
    extendBottom[2] = this->parameters(Simulation()).lz;
    IndicatorCuboid3D<T> bottom(extendBottom, origin);
    /// Set material number for bottom
    this->geometry().rename(2,3,1,bottom);
  }

  void prepareLattices() override
  {
    std::vector<T> dir{0.0, 1.0, 0.0};

    T boussinesqForcePrefactor = 9.81
      / this->converter().getConversionFactorVelocity()
      * this->converter().getConversionFactorTime()
      * this->converter().getCharPhysTemperatureDifference()
      * this->converter().getPhysThermalExpansionCoefficient();

    NavierStokesAdvectionDiffusionCouplingGenerator3D<T,NSDESCRIPTOR>
    coupling(
      0, this->converter().getLatticeLength(this->parameters(Simulation()).lx),
      0, this->converter().getLatticeLength(this->parameters(Simulation()).ly),
      0, this->converter().getLatticeLength(this->parameters(Simulation()).lz),
      boussinesqForcePrefactor,
      this->converter().getLatticeTemperature(this->parameters(Simulation()).Tcold),
      1., dir);

    this->lattice(NavierStokes()).addLatticeCoupling(coupling, this->lattice(Temperature()));

    const T Tomega  = this->converter().getLatticeThermalRelaxationFrequency();
    const T NSomega = this->converter().getLatticeRelaxationFrequency();

    this->lattice(Temperature()).template defineDynamics<NoDynamics>(this->geometry(), 0);
    this->lattice(NavierStokes()).template defineDynamics<NoDynamics>(this->geometry(), 0);

    this->lattice(Temperature()).template defineDynamics<AdvectionDiffusionBGKdynamics>(
      this->geometry().getMaterialIndicator({1, 2, 3}));
    this->lattice(NavierStokes()).template defineDynamics<ForcedBGKdynamics>(
      this->geometry().getMaterialIndicator({1, 2, 3}));

    this->lattice(Temperature()).template setParameter<descriptors::OMEGA>(Tomega);
    this->lattice(NavierStokes()).template setParameter<descriptors::OMEGA>(NSomega);

    /// sets boundary
    setLocalVelocityBoundary<T,NSDESCRIPTOR>(
      this->lattice(NavierStokes()),
      NSomega,
      this->geometry().getMaterialIndicator({2, 3}));
    setAdvectionDiffusionTemperatureBoundary<T,TDESCRIPTOR>(
      this->lattice(Temperature()),
      Tomega,
      this->geometry().getMaterialIndicator({2, 3}));
  }

  void setInitialValues() override
  {
    /// for each material set the defineRhoU and the Equilibrium
    std::vector<T> zero(3,T());
    AnalyticalConst3D<T,T> u(zero);
    AnalyticalConst3D<T,T> rho(1.);
    AnalyticalConst3D<T,T> force(zero);

    T u_Re = this->converter().getLatticeVelocity(
      this->parameters(Simulation()).Re
      * this->converter().getPhysViscosity()
      / this->converter().getCharPhysLength() );
    AnalyticalConst3D<T,T> u_top(this->converter().getCharLatticeVelocity(), u_Re, 0.0);
    AnalyticalConst3D<T,T> u_bot(0.0, u_Re, 0.0);

    auto& NSlattice = this->lattice(NavierStokes());
    auto& Tlattice = this->lattice(Temperature());
    const auto& everywhere = this->geometry().getMaterialIndicator({1,2,3});

    NSlattice.defineRhoU(this->geometry(), 1, rho, u);
    NSlattice.iniEquilibrium(this->geometry(), 1, rho, u);
    NSlattice.defineRhoU(this->geometry(), 2, rho, u_top);
    NSlattice.iniEquilibrium(this->geometry(), 2, rho, u_top);
    NSlattice.defineRhoU(this->geometry(), 3, rho, u_bot);
    NSlattice.iniEquilibrium(this->geometry(), 3, rho, u_bot);
    NSlattice.template defineField<descriptors::FORCE>(everywhere, force);

    AnalyticalConst3D<T,T> Cold(this->converter().getLatticeTemperature(this->parameters(Simulation()).Tcold));
    AnalyticalConst3D<T,T> Hot(this->converter().getLatticeTemperature(this->parameters(Simulation()).Thot));

    Tlattice.defineRho(this->geometry(), 1, Cold);
    Tlattice.iniEquilibrium(this->geometry(), 1, Cold, u);
    Tlattice.defineRho(this->geometry(), 2, Hot);
    Tlattice.iniEquilibrium(this->geometry(), 2, Hot, u);
    Tlattice.defineRho(this->geometry(), 3, Cold);
    Tlattice.iniEquilibrium(this->geometry(), 3, Cold, u);
    Tlattice.template defineField<descriptors::VELOCITY>(everywhere, u);
  }

  void setBoundaryValues(std::size_t iT) override { }

  void writeImages(std::size_t iT) const override {
    SuperLatticePhysTemperature3D<T,NSDESCRIPTOR,TDESCRIPTOR> temperature(this->lattice(Temperature()), this->converter());

    BlockReduction3D2D<T> planeReduction( temperature, Vector<T,3>({0,0,1}), 600, BlockDataSyncMode::ReduceOnly );
    // write output as JPEG
    heatmap::plotParam<T> jpeg_Param;
    jpeg_Param.maxValue = this->parameters(Simulation()).Thot;
    jpeg_Param.minValue = this->parameters(Simulation()).Tcold;
    heatmap::write(planeReduction, iT, jpeg_Param);
  }

  void writeVTK(std::size_t iT) const override
  {
    SuperVTMwriter3D<T> vtkWriter(this->parameters(Output()).vtkFilename);
    auto& NSlattice = this->lattice(NavierStokes());
    auto& Tlattice = this->lattice(Temperature());

    SuperLatticePhysVelocity3D velocity(NSlattice, this->converter());
    SuperLatticePhysPressure3D pressure(NSlattice, this->converter());
    SuperLatticePhysTemperature3D<T,NSDESCRIPTOR,TDESCRIPTOR> temperature(Tlattice, this->converter());

    AnalyticalHeatFluxPorousPlate3D<T,T> HeatFluxSol(
      this->parameters(Simulation()).Re, 
      this->parameters(Simulation()).Pr, 
      this->converter().getCharPhysTemperatureDifference(), 
      this->converter().getCharPhysLength(), 
      this->converter().getThermalConductivity());
    SuperLatticePhysHeatFlux3D<T,NSDESCRIPTOR,TDESCRIPTOR> HeatFlux1(
      Tlattice,
      this->converter());
    SuperLatticeFfromAnalyticalF3D<T,TDESCRIPTOR> HeatFluxSolLattice(
      HeatFluxSol,
      Tlattice);
    vtkWriter.addFunctor( pressure );
    vtkWriter.addFunctor( velocity );
    vtkWriter.addFunctor( temperature );
    vtkWriter.addFunctor( HeatFlux1 );
    vtkWriter.addFunctor( HeatFluxSolLattice );

    vtkWriter.write(iT);
  }

  void computeResults(std::size_t iT) override {
    if (iT % this->converter().getLatticeTime(this->parameters(Output()).logT) == 0) {
      error();
    }
  }

  void computeResults() override {
    error();
  }

private:
  void error() const {
    /// Compute errors w.r.t. analytical reference
    OstreamManager clout( std::cout, "error" );
    int input[1] = { };
    T result[1]  = { };

    auto indicatorF = this->geometry().getMaterialIndicator(1);

    T u_Re = this->parameters(Simulation()).Re * this->converter().getPhysViscosity() / this->converter().getCharPhysLength();
    AnalyticalVelocityPorousPlate3D<T,T> uSol(this->parameters(Simulation()).Re,this->converter().getCharPhysVelocity(), u_Re, this->converter().getCharPhysLength());
    SuperLatticePhysVelocity3D<T,NSDESCRIPTOR> u(this->lattice(NavierStokes()),this->converter());

    SuperAbsoluteErrorL2Norm3D<T> absVelocityErrorNormL2(u, uSol, indicatorF);
    absVelocityErrorNormL2(result, input);
    clout << "velocity-L2-error(abs)=" << result[0];
    SuperRelativeErrorL2Norm3D<T> relVelocityErrorNormL2(u, uSol, indicatorF);
    relVelocityErrorNormL2(result, input);
    clout << "; velocity-L2-error(rel)=" << result[0] << std::endl;

    AnalyticalTemperaturePorousPlate3D<T,T> tSol(this->parameters(Simulation()).Re, this->parameters(Simulation()).Pr,  this->converter().getCharPhysLength(),  this->converter().getCharPhysLowTemperature(), this->converter().getCharPhysTemperatureDifference());
    SuperLatticePhysTemperature3D<T, NSDESCRIPTOR, TDESCRIPTOR> t(this->lattice(Temperature()), this->converter());

    SuperAbsoluteErrorL2Norm3D<T> absTemperatureErrorNormL2(t, tSol, indicatorF);
    absTemperatureErrorNormL2(result, input);
    clout << "temperature-L2-error(abs)=" << result[0];
    SuperRelativeErrorL2Norm3D<T> relTemperatureErrorNormL2(t, tSol, indicatorF);
    relTemperatureErrorNormL2(result, input);
    clout << "; temperature-L2-error(rel)=" << result[0] << std::endl;

    AnalyticalHeatFluxPorousPlate3D<T,T> HeatFluxSol(this->parameters(Simulation()).Re, this->parameters(Simulation()).Pr, this->converter().getCharPhysTemperatureDifference(), this->converter().getCharPhysLength(), this->converter().getThermalConductivity());
    SuperLatticePhysHeatFlux3D<T,NSDESCRIPTOR,TDESCRIPTOR> HeatFlux(this->lattice(Temperature()),this->converter());

    SuperAbsoluteErrorL2Norm3D<T> absHeatFluxErrorNormL2(HeatFlux, HeatFluxSol, indicatorF);
    absHeatFluxErrorNormL2(result, input);
    clout << "heatFlux-L2-error(abs)=" << result[0];
    SuperRelativeErrorL2Norm3D<T> relHeatFluxErrorNormL2(HeatFlux, HeatFluxSol, indicatorF);
    relHeatFluxErrorNormL2(result, input);
    clout << "; heatFlux-L2-error(rel)=" << result[0] << std::endl;
  }
};


int main(int argc, char *argv[])
{
  using T = double;

  olbInit(&argc, &argv);

  // create instances of parameter structs
  utilities::TypeIndexedTuple<Parameters<T>> params;
  params.template get<Simulation>() = std::make_shared<PorousPlate3dSimulationParameters<T>>();
  params.template get<Stationarity>() = std::make_shared<parameters::Stationarity<T,Temperature>>(
    parameters::Stationarity<T,Temperature>::AverageEnergy, 1.0, 1.e-7);
  params.template get<Output>() = std::make_shared<parameters::OutputBase<T>>(
    "thermalPorousPlate3d", "../../../", "./tmp/",
    true, true, 100., 0,
    false, "", 0.,
    true, "thermalPorousPlate3d", 100.,
    true, "thermalPorousPlate3d", 100.);

  // create instance of solver class
  PorousPlate3dSolver<T> porousPlate(params);
  // run simulation
  porousPlate.solve();
}
