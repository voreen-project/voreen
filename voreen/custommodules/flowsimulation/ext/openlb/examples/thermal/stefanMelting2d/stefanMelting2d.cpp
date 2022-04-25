/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2020 Maximilian Gaedke, Larissa Dietz
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

/* The solution for the melting problem (solid-liquid phase change)
is found using the lattice Boltzmann method after Rongzong Huang
and Huiying Wu (2015)[1]. The equilibrium distribution
function for the temperature is modified in order to deal with the
latent-heat source term. That way, iteration steps or solving a group
of linear equations is avoided, which results in enhanced efficiency.
The phase interface is located by the current total enthalpy, and
its movement is considered by the immersed moving boundary scheme after
Noble and Torczynski (1998)[2]. This method was validated by the
problem of conduction-induced melting in a semi-infinite space,
comparing its results to analytical solutions.

[1] Rongzong Huang, Huiying Wu, Phase interface effects in the total enthalpy-based lattice
Boltzmann model for solid–liquid phase change, Journal of Computational Physics 294 (2015) 346–362.

[2] D. Noble, J. Torczynski, A lattice-Boltzmann method for partially saturated
computational cells, Int. J. Modern Phys. C 9 (8) (1998) 1189–1202.
 */

#include "olb2D.h"
#include "olb2D.hh"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = double;

using NSDESCRIPTOR = D2Q9<POROSITY,VELOCITY_SOLID,FORCE>;
using TDESCRIPTOR  = D2Q5<VELOCITY,TEMPERATURE>;

using TotalEnthalpyAdvectionDiffusionDynamics = TotalEnthalpyAdvectionDiffusionBGKdynamics<T,TDESCRIPTOR>;

// Parameters for the simulation setup
const T lx  = 1.0;      // length of the channel
T ly  = lx / 8.;        // height of the channel
int     N = 128;        // resolution of the model
T       tau = 1.;       // relaxation time

const T Re = 20.;       // Reynolds number
const T Pr = 0.71;      // Prandtl number
T maxPhysT;             // simulation time in s (calculated later as a function of k)

const T Tcold = 0.5;    // cold temperature
const T Thot = 1.5;     // hot temperature

const T lambda = 1./6.; // W / m K
const T cp_s = 1.0;     // J / kg K
const T cp_l = 1.0;     // J / kg K
const T cp_ref = 2.0 * cp_s * cp_l / (cp_s + cp_l); // J /kg K
const T density = 1.0;  // kg / m^3

const T Ste = 0.01;     // Stephan number
T k;

const T L = cp_l * (Thot - Tcold) / Ste; // J / kg

T lattice_Hcold, lattice_Hhot;
T physDeltaX, physDeltaT;

size_t iT = 0;

T sum_error_melt = 0.0;      //cumulative error melting (L2-norm)
T sum_error_temp = 0.0;      //cumulative error temperature (L2-norm)
T sum_error_melt_inf = 0.0;  //cumulative error melting (Linf-norm))
T sum_error_temp_inf = 0.0;  //cumulative error temperature inf(Linf-norm)
int num_error_samples = 25;  //number of samples for error averaging
int num_errors = 0;          //number of current samples in the sums

T func ( T k )
{
  return Ste / util::exp(k*k) / erf(k) - k * util::sqrt(M_PI);
}

T func_deriv (T k)
{
  return -(2 * Ste * util::exp(-k*k)*k/erf(k))-(2 * Ste * util::exp(-2*k*k)/(util::sqrt(M_PI)*(erf(k)*erf(k))))-util::sqrt(M_PI);
}

template <typename T, typename S>
class AnalyticalMeltFraction2D : public AnalyticalF2D<T, S> {
public:
  AnalyticalMeltFraction2D() : AnalyticalF2D<T, S>(1)
  {
    this->getName() = "AnalyticalMeltFraction2D";
  };

  bool operator()(T output[1], const S x[2])
  {
    T X = 2.0 * k * util::sqrt(lambda / density / cp_l * iT * physDeltaT);
    if (util::nearZero(x[0]-X)) {
      output[0] = 0.5;
    }
    else if (x[0] < X) {
      output[0] = 1.0;
    }
    else if (x[0] > X) {
      output[0] = 0.0;
    }

    return true;
  };
};

template <typename T, typename S>
class AnalyticalTemperature2D : public AnalyticalF2D<T, S> {
public:
  AnalyticalTemperature2D() : AnalyticalF2D<T, S>(1)
  {
    this->getName() = "AnalyticalTemperature2D";
  };

  bool operator()(T output[1], const S x[2])
  {
    //calculates maximum distance of melt front compared to position of starting point
    T X = 2.0 * k * util::sqrt(lambda / density / cp_l * iT * physDeltaT);
    if (x[0] < X) {
      output[0] = Thot - (Thot - Tcold) / erf(k) * erf(x[0] * k / X);
    }
    if (x[0] >= X) {
      output[0] = Tcold;
    }

    return true;
  };
};

//finds root of func(T x),
T func_newton(T startValue, int maxIterations, T maxError)
{
  int i = 0;
  T x = startValue;
  T frac = func(x)/func_deriv(x);

  while (i < maxIterations && util::abs(frac) > maxError) {
    frac = func(x)/func_deriv(x);
    x = x - frac;
    i++;
  }
  return x;
}

void error( SuperGeometry<T,2>& superGeometry,
            SuperLattice<T, NSDESCRIPTOR>& NSlattice,
            SuperLattice<T, TDESCRIPTOR>& ADlattice,
            ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter,
            T Re )
{
  OstreamManager clout( std::cout, "error" );

  int input[1] = { };
  T result[1]  = { };

  auto indicatorF = superGeometry.getMaterialIndicator({1});

  SuperLatticeField2D<T, NSDESCRIPTOR, POROSITY> liquid_frac(NSlattice);
  liquid_frac.getName() = "liquid fraction";
  AnalyticalMeltFraction2D<T,T> melt_solution;

  SuperLatticeFfromAnalyticalF2D<T,TDESCRIPTOR> melt_solution_lattice(melt_solution, ADlattice);

  SuperAbsoluteErrorL2Norm2D<T> absMeltFractionErrorNormL2(liquid_frac, melt_solution, indicatorF);
  absMeltFractionErrorNormL2(result, input);
  sum_error_melt += result[0];
  clout << "melt-fraction-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm2D<T> relMeltFractionErrorNormL2(liquid_frac, melt_solution, indicatorF);
  relMeltFractionErrorNormL2(result, input);
  clout << "; melt-fraction-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm2D<T> absMeltFractionErrorNormLinf(liquid_frac, melt_solution, indicatorF);
  absMeltFractionErrorNormLinf(result, input);
  sum_error_melt_inf += result[0];
  clout << "melt-fraction-Linf-error(abs)=" << result[0];
  SuperRelativeErrorLinfNorm2D<T> relMeltFractionErrorNormLinf(liquid_frac, melt_solution, indicatorF);
  relMeltFractionErrorNormLinf(result, input);
  clout << "; melt-fraction-Linf-error(rel)=" << result[0] << std::endl;

  SuperLatticeField2D<T, TDESCRIPTOR, TEMPERATURE> temperature(ADlattice);
  temperature.getName() = "temperature";
  AnalyticalTemperature2D<T,T> temp_solution;
  SuperLatticeFfromAnalyticalF2D<T,TDESCRIPTOR> temp_solution_lattice(temp_solution, ADlattice);

  SuperAbsoluteErrorL2Norm2D<T> absTemperatureErrorNormL2(temperature, temp_solution, indicatorF);
  absTemperatureErrorNormL2(result, input);
  sum_error_temp += result[0];
  clout << "temperature-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm2D<T> relTemperatureErrorNormL2(temperature, temp_solution, indicatorF);
  relTemperatureErrorNormL2(result, input);
  clout << "; temperature-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm2D<T> absTemperatureErrorNormLinf(temperature, temp_solution, indicatorF);
  absTemperatureErrorNormLinf(result, input);
  sum_error_temp_inf += result[0];
  clout << "temperature-Linf-error(abs)=" << result[0];
  SuperRelativeErrorLinfNorm2D<T> relTemperatureErrorNormLinf(temperature, temp_solution, indicatorF);
  relTemperatureErrorNormLinf(result, input);
  clout << "; temperature-Linf-error(rel)=" << result[0] << std::endl;

  num_errors++;

}

/// Stores geometry information in form of material numbers
void prepareGeometry(SuperGeometry<T,2>& superGeometry,
                     ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter)
{

  OstreamManager clout(std::cout,"prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0,2);
  superGeometry.rename(2,1,{1,0});
  std::vector<T> extend( 2, T(0) );
  extend[0] = converter.getPhysLength(1);
  extend[1] = ly + 2.*converter.getPhysLength(1);
  std::vector<T> origin( 2, T(0) );
  origin[0] = -converter.getPhysLength(1);
  IndicatorCuboid2D<T> left(extend, origin);
  superGeometry.rename(2,3,1,left);

  /// Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  /// Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice( ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter,
                     SuperLattice<T, NSDESCRIPTOR>& NSlattice,
                     SuperLattice<T, TDESCRIPTOR>& ADlattice,
                     SuperGeometry<T,2>& superGeometry )
{
  T Tomega  = converter.getLatticeThermalRelaxationFrequency();
  T NSomega = converter.getLatticeRelaxationFrequency();

  ADlattice.defineDynamics<NoDynamics>(superGeometry, 0);
  NSlattice.defineDynamics<NoDynamics>(superGeometry, 0);

  ADlattice.defineDynamics<TotalEnthalpyAdvectionDiffusionDynamics>(superGeometry, 1);
  ADlattice.defineDynamics<BounceBack>(superGeometry, 2);
  ADlattice.defineDynamics<TotalEnthalpyAdvectionDiffusionDynamics>(superGeometry, 3);
  NSlattice.defineDynamics<ForcedPSMBGKdynamics>(superGeometry, 1);
  NSlattice.defineDynamics<ForcedPSMBGKdynamics>(superGeometry, 2);
  NSlattice.defineDynamics<ForcedPSMBGKdynamics>(superGeometry, 3);

  setLocalVelocityBoundary(NSlattice, NSomega, superGeometry, 2);
  setLocalVelocityBoundary(NSlattice, NSomega, superGeometry, 3);

  setRegularizedTemperatureBoundary(ADlattice, Tomega, superGeometry, 3);

  NSlattice.setParameter<descriptors::OMEGA>(NSomega);

  ADlattice.setParameter<descriptors::OMEGA>(Tomega);
  ADlattice.setParameter<TotalEnthalpy::T_S>(Tcold);
  ADlattice.setParameter<TotalEnthalpy::T_L>(Tcold);
  ADlattice.setParameter<TotalEnthalpy::CP_S>(cp_s);
  ADlattice.setParameter<TotalEnthalpy::CP_L>(cp_l);
  ADlattice.setParameter<TotalEnthalpy::LAMBDA_S>(cp_ref / descriptors::invCs2<T,TDESCRIPTOR>() * (tau - 0.5));
  ADlattice.setParameter<TotalEnthalpy::LAMBDA_L>(cp_ref / descriptors::invCs2<T,TDESCRIPTOR>() * (tau - 0.5));
  ADlattice.setParameter<TotalEnthalpy::L>(L);

}

void setBoundaryValues(ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter,
                       SuperLattice<T, NSDESCRIPTOR>& NSlattice,
                       SuperLattice<T, TDESCRIPTOR>& ADlattice,
                       int iT, SuperGeometry<T,2>& superGeometry)
{

  if (iT == 0) {

    /// for each material set the defineRhoU and the Equilibrium
    std::vector<T> zero(2,T());
    AnalyticalConst2D<T,T> u(zero);
    AnalyticalConst2D<T,T> rho(1.);
    AnalyticalConst2D<T,T> force(zero);

    T u_Re = converter.getLatticeVelocity( Re * converter.getPhysViscosity() / converter.getCharPhysLength() );

    AnalyticalConst2D<T,T> u_top(converter.getCharLatticeVelocity(), u_Re);
    AnalyticalConst2D<T,T> u_bot(0.0, u_Re);

    NSlattice.defineRhoU(superGeometry.getMaterialIndicator({1,2,3}), rho, u);
    NSlattice.iniEquilibrium(superGeometry.getMaterialIndicator({1,2,3}), rho, u);

    AnalyticalConst2D<T,T> Cold(lattice_Hcold);
    AnalyticalConst2D<T,T> Hot(lattice_Hhot);

    ADlattice.defineField<descriptors::VELOCITY>(superGeometry.getMaterialIndicator({1, 2, 3}), u);
    ADlattice.defineRho(superGeometry.getMaterialIndicator({1,2}), Cold);
    ADlattice.iniEquilibrium(superGeometry.getMaterialIndicator({1,2}), Cold, u);
    ADlattice.defineRho(superGeometry, 3, Hot);
    ADlattice.iniEquilibrium(superGeometry, 3, Hot, u);

    /// Make the lattice ready for simulation
    NSlattice.initialize();
    ADlattice.initialize();
  }
}

void getResults(ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter,
                SuperLattice<T, NSDESCRIPTOR>& NSlattice,
                SuperLattice<T, TDESCRIPTOR>& ADlattice, int iT,
                SuperGeometry<T,2>& superGeometry,
                util::Timer<T>& timer,
                bool converged)
{
  OstreamManager clout(std::cout,"getResults");

  SuperVTMwriter2D<T> vtkWriter("stefanMelting1d");
  SuperLatticeGeometry2D<T, NSDESCRIPTOR> geometry2(NSlattice, superGeometry);
  SuperLatticePhysVelocity2D<T, NSDESCRIPTOR> velocity(NSlattice, converter);
  SuperLatticePhysPressure2D<T, NSDESCRIPTOR> pressure(NSlattice, converter);
  SuperLatticePhysHeatFlux2D<T, NSDESCRIPTOR, TDESCRIPTOR> heatflux(ADlattice, converter);
  SuperLatticeDensity2D<T, TDESCRIPTOR> enthalpy(ADlattice);
  enthalpy.getName() = "enthalpy";
  SuperLatticeField2D<T, NSDESCRIPTOR, POROSITY> liquid_frac(NSlattice);
  liquid_frac.getName() = "liquid fraction";
  SuperLatticeField2D<T, TDESCRIPTOR, TEMPERATURE> temperature(ADlattice);
  temperature.getName() = "temperature";

  AnalyticalMeltFraction2D<T,T> melt_solution;
  SuperLatticeFfromAnalyticalF2D<T,TDESCRIPTOR> melt_solution_lattice(melt_solution, ADlattice);
  AnalyticalTemperature2D<T,T> temp_solution;
  SuperLatticeFfromAnalyticalF2D<T,TDESCRIPTOR> temp_solution_lattice(temp_solution, ADlattice);

  vtkWriter.addFunctor( geometry2 );
  vtkWriter.addFunctor( pressure );
  vtkWriter.addFunctor( velocity );
  vtkWriter.addFunctor( enthalpy );
  vtkWriter.addFunctor( temperature );
  vtkWriter.addFunctor( heatflux );
  vtkWriter.addFunctor( liquid_frac );
  vtkWriter.addFunctor( melt_solution_lattice );
  vtkWriter.addFunctor( temp_solution_lattice );

  const int vtkIter = converter.getLatticeTime(maxPhysT/T(num_error_samples));

  if (iT == 0) {
    /// Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry2D<T, NSDESCRIPTOR> geometry(NSlattice, superGeometry);
    SuperLatticeCuboid2D<T, NSDESCRIPTOR> cuboid(NSlattice);
    SuperLatticeRank2D<T, NSDESCRIPTOR> rank(NSlattice);
    vtkWriter.write(geometry);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);

    vtkWriter.createMasterFile();
  }

  /// Writes the VTK files
  if (iT%vtkIter == 0 || converged) {
    ADlattice.getStatistics().print(iT,converter.getPhysTime(iT));
    timer.print(iT);
    error(superGeometry, NSlattice, ADlattice, converter, Re);

    vtkWriter.write(iT);
  }
}

int main(int argc, char *argv[])
{

  /// === 1st Step: Initialization ===
  OstreamManager clout(std::cout,"main");
  olbInit(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");

  k = func_newton (1e-5, 100, 2e-7);

  // run simulation until one forth of the domain is molten
  maxPhysT = util::pow(0.25 / 2.0 / k, 2) * density * cp_l / lambda;

  if (argc >= 2) {
    N = atoi(argv[1]);
  }
  if (argc >= 3) {
    tau = atof(argv[2]);
  }

  physDeltaX = lx / N;
  physDeltaT = density * cp_ref / lambda / descriptors::invCs2<T,NSDESCRIPTOR>() * (tau - 0.5) * physDeltaX * physDeltaX;

  lattice_Hcold = cp_s * Tcold;
  lattice_Hhot = cp_l * Thot;

  ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const converter(
    (T) physDeltaX, // physDeltaX
    (T) physDeltaT, // physDeltaT
    (T) lx, // charPhysLength
    (T) 1.0, // charPhysVelocity
    (T) lambda / density / cp_l, // physViscosity
    (T) density, // physDensity
    (T) lambda, // physThermalConductivity
    (T) cp_l, // physSpecificHeatCapacity
    (T) 1.0, // physThermalExpansionCoefficient
    (T) Tcold, // charPhysLowTemperature
    (T) Thot // charPhysHighTemperature
  );
  converter.print();

  clout << "H_cold: " << lattice_Hcold << std::endl;
  clout << "H_hot: " << lattice_Hhot << std::endl;
  clout << "k: " << std::setprecision(17) << k << std::endl;
  clout << "Ste: " << Ste << std::endl;
  clout << "lattice cp: " << converter.getLatticeSpecificHeatCapacity(cp_l) << std::endl;

  std::vector<T> extend(2,T());
  extend[0] = lx;
  extend[1] = ly;
  std::vector<T> origin(2,T());
  IndicatorCuboid2D<T> cuboid(extend, origin);


  /// Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif
//the cuboids not needed are removed and too big ones are shrinked
  CuboidGeometry2D<T> cuboidGeometry(cuboid, converter.getPhysDeltaX(), noOfCuboids);
  cuboidGeometry.setPeriodicity(false, true);

  /// Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

  /// Instantiation of a superGeometry
  SuperGeometry<T,2> superGeometry(cuboidGeometry, loadBalancer);

  prepareGeometry(superGeometry, converter);

  /// === 3rd Step: Prepare Lattice ===

  SuperLattice<T, TDESCRIPTOR> ADlattice(superGeometry);
  SuperLattice<T, NSDESCRIPTOR> NSlattice(superGeometry);



  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
  // This coupling must be necessarily be put on the Navier-Stokes lattice!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//

  std::vector<T> dir{0.0, 1.0};

  T boussinesqForcePrefactor = 9.81 / converter.getConversionFactorVelocity() * converter.getConversionFactorTime() *
                               converter.getCharPhysTemperatureDifference() * converter.getPhysThermalExpansionCoefficient();

  TotalEnthalpyPhaseChangeCouplingGenerator2D<T,NSDESCRIPTOR,TotalEnthalpyAdvectionDiffusionDynamics>
  coupling(0, converter.getLatticeLength(lx), 0, converter.getLatticeLength(ly),
           boussinesqForcePrefactor, converter.getLatticeTemperature(Tcold), 1., dir);

  NSlattice.addLatticeCoupling(superGeometry, 1, coupling, ADlattice);

  //prepareLattice and setBoundaryConditions
  prepareLattice(converter, NSlattice, ADlattice, superGeometry);

  /// === 4th Step: Main Loop with Timer ===
  util::Timer<T> timer(converter.getLatticeTime(maxPhysT), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for (iT = 0; iT < converter.getLatticeTime(maxPhysT)+1; ++iT) {

    /// === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues(converter, NSlattice, ADlattice, iT, superGeometry);

    /// === 6th Step: Collide and Stream Execution ===
    NSlattice.executeCoupling();
    std::vector<T> zero(2,T());
    AnalyticalConst2D<T,T> u(zero);
    ADlattice.defineField<descriptors::VELOCITY>(superGeometry, 1, u);
    /// === 7th Step: Computation and Output of the Results ===
    getResults(converter, NSlattice, ADlattice, iT, superGeometry, timer, false);

    // NSlattice.collideAndStream();
    ADlattice.collideAndStream();

  }
  clout << "liquid_fraction-L2-error(abs,mean)=" << sum_error_melt / T(num_errors) << std::endl;
  clout << "temperature-L2-error(abs,mean)=" << sum_error_temp / T(num_errors) << std::endl;
  clout << "liquid_fraction-Linf-error(abs,mean)=" << sum_error_melt_inf / T(num_errors) << std::endl;
  clout << "temperature-Linf-error(abs,mean)=" << sum_error_temp_inf / T(num_errors) << std::endl;

  timer.stop();
  timer.printSummary();
}
