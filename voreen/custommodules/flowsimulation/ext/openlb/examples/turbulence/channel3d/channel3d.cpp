/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2020 Jonathan Jeppener-Haltenhoff, Marc Haußmann, Mathias J. Krause
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

/* channel3d.cpp:
 * This example examines a wall-bounded flow at high Reynolds numbers.
 * The dynamics follow the BGK collision operator and Smagorinsky-Lilly turbulence model.
 * The near wall region is modelled by a wall function.
 *
 * The example shows the usage of turbulent wall models.
 *
 * The results are published in
 * Haußmann, M. et al. 2019: Large-eddy simulation coupled with wall models for turbulent
 * channel flows at high Reynolds numbers with a lattice Boltzmann method — Application
 * to Coriolis mass flowmeter. In: Computers & Mathematics with Applications. 78 (10). 3285-3302
 */

#include "olb3D.h"
#include "olb3D.hh"   // include full template code

using namespace olb;
using namespace olb::descriptors;
using namespace olb::util;

typedef double T;
typedef WallFunctionForcedD3Q19Descriptor DESCRIPTOR;

// Mathmatical constants
const T pi = util::acos(-1);

// Parameters for the simulation setup
const int N = 18;
const T physRefL = 1.0;     // half channel height in meters
const T lx = 2. * pi * physRefL;  // streamwise length in meters
const T ly = 2. * pi * physRefL;  // spanwise length in meters
const T lz = 2. * physRefL;       // wall-normal length in meters

// Choose friction reynolds number ReTau
#define Case_ReTau_1000
//#define Case_ReTau_2000

// Wallfunction parameters
const T latticeWallDistance = 0.5;  // lattice distance to boundary
const int rhoMethod = 0;    // method for density reconstruction
// 0: Zou-He
// 1: extrapolation
// 2: constant
const int fneqMethod = 1;   // method for fneq reconstruction
// 0: regularized NEBB (Latt)
// 1: extrapolation NEQ (Guo Zhaoli)
// 2: regularized second order finite Differnce
// 3: equilibrium scheme
const int wallProfile = 0;    // wallfunction profile
// 0: Musker profile
// 1: power law profile

// Reynolds number based on the friction velocity
#if defined (Case_ReTau_1000)
T ReTau = 1000.512;
#elif defined (Case_ReTau_2000)
T ReTau = 1999.756;
#endif
// Characteristic physical kinematic viscosity from DNS Data
// http://turbulence.ices.utexas.edu/channel2015/data/LM_Channel_1000_mean_prof.dat
#if defined (Case_ReTau_1000)
T charPhysNu = 5./100000.;
#elif defined (Case_ReTau_2000)
T charPhysNu = 2.3/100000.;
#endif

// number of forcing updates over simulation time
#if defined (Case_ReTau_1000)
T fluxUpdates = 2000;
#elif defined (Case_ReTau_2000)
T fluxUpdates = 4000;
#endif

// physical simulated length adapted for lattice distance to boundary in meters
const T adaptedPhysSimulatedLength = 2 * physRefL - 2 * ((2. / T(N + 2 * latticeWallDistance)) * latticeWallDistance);
// Characteristic physical mean bulk velocity from Dean correlations in meters - Malaspinas and Sagaut (2014)
const T charPhysU = ( util::pow((8.0/0.073), (4.0/7.0)) * util::pow((T)ReTau, (8.0/7.0)) ) * charPhysNu / (2. * physRefL);
// Time of the simulation in seconds
const T charPhysT = physRefL / (ReTau * charPhysNu / physRefL);

const T physConvergeTime = 120. * charPhysT;  // time until until statistics sampling in seconds
const T physStatisticsTime = 40. * charPhysT ;  // statistics sampling time in seconds
const T maxPhysT = physConvergeTime + physStatisticsTime; // max. simulation time in seconds
const T statisticsSave = 1./25.;    // time between statistics samples in seconds

// seed the rng with time if SEED_WITH_TIME is set, otherwise just use a fixed seed.
#if defined SEED_WITH_TIME
#include <chrono>
auto seed = std::chrono::system_clock().now().time_since_epoch().count();
std::default_random_engine generator(seed);
#else
std::default_random_engine generator(0x1337533DAAAAAAAA);
#endif

// Compute mean lattice velocity from musker wallfunction
T computeLatticeVelocity()
{
  T Ma_max = 0.1;
  T c_s = 1/util::sqrt(3.0);
  T latticeUMax = Ma_max * c_s;
  Musker<T,T> musker_tmp(charPhysNu, physRefL, 1.);
  T charPhysU_tau = ReTau * charPhysNu / physRefL;
  T tau_w[1];
  tau_w[0] = util::pow(charPhysU_tau,2.);
  T charPhysUMax[1];
  musker_tmp(charPhysUMax,tau_w);
  T latticeU = charPhysU * latticeUMax / charPhysUMax[0];
  return latticeU;
}

template <typename T, typename S>
class Channel3D : public AnalyticalF3D<T,S> {

protected:
  T turbulenceIntensity;
  T maxVelocity;
  T distanceToWall;
  T obst_z;
  T obst_r;
  T a;
  T b;

public:
  Channel3D(UnitConverter<T,DESCRIPTOR> const& converter, T frac) : AnalyticalF3D<T,S>(3)
  {
    turbulenceIntensity = 0.05;
    distanceToWall = -converter.getPhysDeltaX()/2.;
    maxVelocity = converter.getLatticeVelocity(converter.getCharPhysVelocity()*(8./7.)); // Centerline Velocity
    obst_z = physRefL + distanceToWall;
    obst_r = physRefL;
    a = -1.;
    b = 1.;
  };

  bool operator()(T output[], const S input[])
  {
    std::uniform_real_distribution<T> distribution(a, b);
    T nRandom1 = distribution(generator);
    T nRandom2 = distribution(generator);
    T nRandom3 = distribution(generator);

    T u_calc = maxVelocity*util::pow(((obst_r-util::abs(input[2] - obst_z))/obst_r), 1./7.);

    output[0] = turbulenceIntensity*nRandom1*maxVelocity + u_calc;
    output[1] = turbulenceIntensity*nRandom2*maxVelocity;
    output[2] = turbulenceIntensity*nRandom3*maxVelocity;

    return true;
  };
};

template <typename T, typename S>
class TrackedForcing3D : public AnalyticalF3D<T,S> {

protected:
  T um;
  T utau;
  T h2;
  T aveVelocity;

public:
  TrackedForcing3D(UnitConverter<T,DESCRIPTOR> const& converter, int ReTau) : AnalyticalF3D<T, S>(3)
  {
    um = converter.getCharPhysVelocity();
    utau = ReTau * converter.getPhysViscosity()/(converter.getCharPhysLength()/2.);
    h2 = converter.getCharPhysLength()/2.;
    aveVelocity = um;
  };

  void updateAveVelocity(T newVel)
  {
    aveVelocity = newVel;
  }

  bool operator()(T output[], const S input[])
  {
    output[0] = util::pow(utau, 2)/h2 + (um - aveVelocity)*um/h2;
    output[1] = 0;
    output[2] = 0;

    return true;
  };
};

void prepareGeometry(SuperGeometry<T,3>& superGeometry, IndicatorF3D<T>& indicator,
                     UnitConverter<T,DESCRIPTOR> const& converter)
{
  OstreamManager clout(std::cout,"prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0,2,indicator);
  superGeometry.rename(2,1,{0,0,1});
  superGeometry.clean();
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  olb::Vector<T, 3> PhyMax = superGeometry.getStatistics().getMaxPhysR(2);
  olb::Vector<T, 3> PhyMin = superGeometry.getStatistics().getMinPhysR(2);
  clout << "Dimension of the channel in meters: x = " << PhyMax[0] - PhyMin[0];
  clout << " ; y = " << PhyMax[1] - PhyMin[1];
  clout << " ; z = " << PhyMax[2] - PhyMin[2] << std::endl;

  clout << "Prepare Geometry ... OK" << std::endl;
}

// set up initial conditions
void setInitialConditions(SuperLattice<T, DESCRIPTOR>& sLattice,
                          UnitConverter<T,DESCRIPTOR> const& converter,
                          SuperGeometry<T,3>& superGeometry,
                          AnalyticalScaled3D<T, T>& forceSolScaled,
                          TrackedForcing3D<T, T>& forceSol)
{

  OstreamManager clout(std::cout, "setInitialConditions");
  clout << "Set initial conditions ..." << std::endl;

  AnalyticalConst3D<T, T> rho(1.);
  Channel3D<T, T> uSol(converter, 1.);

  sLattice.defineRhoU(superGeometry, 1, rho, uSol);
  sLattice.iniEquilibrium(superGeometry, 1, rho, uSol);

  sLattice.defineRhoU(superGeometry, 2, rho, uSol);
  sLattice.iniEquilibrium(superGeometry, 2, rho, uSol);

  AnalyticalConst3D<T,T> TauEff(1./converter.getLatticeRelaxationFrequency());

  sLattice.defineField<TAU_EFF>(superGeometry, 1, TauEff);
  sLattice.defineField<TAU_EFF>(superGeometry, 2, TauEff);

  // Force Initialization
  forceSol.updateAveVelocity(converter.getCharPhysVelocity()); // New average velocity

  // Initialize force
  sLattice.defineField<FORCE>(superGeometry, 1, forceSolScaled);
  sLattice.defineField<FORCE>(superGeometry, 2, forceSolScaled);
  // Tau_w Initialization
  T tau_w_guess = 0.0; // Wall shear stress in phys units
  AnalyticalConst3D<T, T> tau_w_ini(tau_w_guess);
  AnalyticalScaled3D<T, T> tau_w_ini_scaled(tau_w_ini, 1. / ( converter.getConversionFactorForce() * util::pow(converter.getConversionFactorLength(),2.) ) );

  sLattice.defineField<TAU_W>(superGeometry, 1, tau_w_ini_scaled);
  sLattice.defineField<TAU_W>(superGeometry, 2, tau_w_ini_scaled);

  clout << "Set initial conditions ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice(SuperLattice<T,DESCRIPTOR>& sLattice,
                    UnitConverter<T,DESCRIPTOR> const& converter,
                    SuperGeometry<T,3>& superGeometry,
                    AnalyticalScaled3D<T, T>& forceSolScaled,
                    TrackedForcing3D<T, T>& forceSol,
                    wallFunctionParam<T> const& wallFunctionParam
                   )
{

  OstreamManager clout(std::cout,"prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  /// Material=0 -->do nothing
  sLattice.defineDynamics<NoDynamics>(superGeometry, 0);
  /// Material=1 -->bulk dynamics
  sLattice.defineDynamics<SmagorinskyForcedBGKdynamics>(superGeometry, 1);
  /// Material = 2 --> boundary node + wallfunction
  sLattice.defineDynamics<ExternalTauEffLESForcedBGKdynamics>(superGeometry, 2);
  setWallFunctionBoundary<T,DESCRIPTOR>(sLattice, superGeometry, 2, converter, wallFunctionParam);

  /// === Set Initial Conditions == ///
  setInitialConditions(sLattice, converter, superGeometry, forceSolScaled, forceSol);

  sLattice.setParameter<descriptors::OMEGA>( converter.getLatticeRelaxationFrequency() );
  sLattice.setParameter<collision::LES::Smagorinsky>(0.12);

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

/// Computes the pressure drop between the voxels before and after the cylinder
void getResults(SuperLattice<T, DESCRIPTOR>& sLattice,
                UnitConverter<T,DESCRIPTOR> const& converter, size_t iT,
                SuperGeometry<T,3>& superGeometry, Timer<double>& timer,
                SuperLatticeTimeAveragedF3D<T>& sAveragedVel)
{
  OstreamManager clout(std::cout, "getResults");

  const T checkstatistics = (T)maxPhysT/200.;

  SuperVTMwriter3D<T> vtmWriter("channel3d");
  SuperLatticeGeometry3D<T, DESCRIPTOR> geometry(sLattice, superGeometry);
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(sLattice, converter);
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure(sLattice, converter);
  vtmWriter.addFunctor(geometry);
  vtmWriter.addFunctor(velocity);
  vtmWriter.addFunctor(pressure);
  vtmWriter.addFunctor(sAveragedVel);

  if (iT == 0) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry3D<T, DESCRIPTOR> geometry(sLattice, superGeometry);
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid(sLattice);
    SuperLatticeRank3D<T, DESCRIPTOR> rank(sLattice);

    vtmWriter.write(geometry);
    vtmWriter.write(cuboid);
    vtmWriter.write(rank);

    vtmWriter.createMasterFile();
  }

  // Writes output on the console
  if (iT % converter.getLatticeTime(checkstatistics) == 0) {
    // Timer console output
    timer.update(iT);
    timer.printStep(2);
    // Lattice statistics console output
    sLattice.getStatistics().print(iT, iT * converter.getPhysDeltaT());
    clout << "Max. physical velocity(m/s): " << converter.getPhysVelocity(sLattice.getStatistics().getMaxU()) << std::endl;
    clout << "Max u+:" << converter.getPhysVelocity(sLattice.getStatistics().getMaxU())
          / (ReTau * charPhysNu / (converter.getCharPhysLength()/2.)) << std::endl;
  }

  if (iT%converter.getLatticeTime(statisticsSave) == 0 && iT > converter.getLatticeTime(physConvergeTime)) {
    // Add ensemble to temporal averaged velocity
    sLattice.communicate();
    sAveragedVel.addEnsemble();
  }

  if (iT%converter.getLatticeTime(maxPhysT/20)==0 || iT==converter.getLatticeTime(maxPhysT)-1) {
    // Writes the vtk files
    vtmWriter.write(iT);
  }
}

int main(int argc, char* argv[])
{

  /// === 1st Step: Initialization ===
  olbInit(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");
  OstreamManager clout(std::cout, "main");
  // display messages from every single mpi process
  //clout.setMultiOutput(true);

  UnitConverterFromResolutionAndLatticeVelocity<T,DESCRIPTOR> const converter(
    int {N},                  // resolution: number of voxels per charPhysL
    (T)   computeLatticeVelocity(),     // latticeU : mean lattice velocity
    (T)   adaptedPhysSimulatedLength, // charPhysLength: reference length of simulation geometry
    (T)   1.0,                // charPhysVelocity: mean bulk velocity in __m / s__
    (T)   charPhysNu,           // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.0                 // physDensity: physical density in __kg / m^3__
  );
  converter.print();

  clout << "----------------------------------------------------------------------" << std::endl;
  clout << "Converge time(s): " << physConvergeTime << std::endl;
  clout << "Lattice converge time: " << converter.getLatticeTime(physConvergeTime) << std::endl;
  clout << "Max. Phys. simulation time(s): " << maxPhysT << std::endl;
  clout << "Max. Lattice simulation time: " << converter.getLatticeTime(maxPhysT) << std::endl;
  clout << "Frequency Statistics Save(Hz): " << 1./statisticsSave << std::endl;
  clout << "Statistics save period(s): " << statisticsSave << std::endl;
  clout << "Lattice statistics save period: " << converter.getLatticeTime(statisticsSave) << std::endl;
  clout << "----------------------------------------------------------------------" << std::endl;
  clout << "Channel height(m): " << adaptedPhysSimulatedLength << std::endl;
  clout << "y+ value: " << (ReTau * converter.getPhysViscosity() / (physRefL)) * ((2. / T(N + 2 * latticeWallDistance)) * latticeWallDistance)/ converter.getPhysViscosity() << std::endl;
  clout << "y+ value spacing: " << (ReTau * converter.getPhysViscosity() / (physRefL)) * (converter.getPhysDeltaX()) / converter.getPhysViscosity() << std::endl;
  clout << "----------------------------------------------------------------------" << std::endl;

#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif

  Vector<T,3> extend(lx, ly, adaptedPhysSimulatedLength);
  extend[2] += (1./8.)*converter.getPhysDeltaX();

  Vector<T,3> origin(0., 0., 0.);
  IndicatorCuboid3D<T> cuboid(extend, origin );

  CuboidGeometry3D<T> cuboidGeometry( cuboid, converter.getPhysDeltaX(), noOfCuboids );

  cuboidGeometry.setPeriodicity(true, true, false);

  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  SuperGeometry<T,3> superGeometry(cuboidGeometry, loadBalancer);

  prepareGeometry(superGeometry, cuboid, converter);

  /// === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice(superGeometry);

  //forcing of the channel
  TrackedForcing3D<T,T> forceSol(converter, ReTau);
  AnalyticalScaled3D<T,T> forceSolScaled(forceSol, 1./(converter.getConversionFactorForce()/converter.getConversionFactorMass()));

  int input[3];
  T output[5];

  Vector<T,3> normal( 1, 0, 0 );
  std::vector<int> normalvec{ 1, 0, 0 };
  Vector<T,3> center;
  for (int i = 0 ; i <3; ++i) {
    center[i] = origin[i] + extend[i]/2;
  }
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity(sLattice, converter);
  std::vector<int> materials;
  materials.push_back( 1 );
  materials.push_back( 2 );

  std::list<int> materialslist;
  materialslist.push_back( 1 );
  materialslist.push_back( 2 );

  wallFunctionParam<T> wallFunctionParam;
  wallFunctionParam.bodyForce = true;
  wallFunctionParam.wallProfile = wallProfile;
  wallFunctionParam.rhoMethod = rhoMethod;
  wallFunctionParam.fneqMethod = fneqMethod;
  wallFunctionParam.latticeWalldistance = latticeWallDistance;
  wallFunctionParam.vonKarman = 0.375;
  wallFunctionParam.curved = false;

  prepareLattice(sLattice, converter, superGeometry,
                 forceSolScaled, forceSol, wallFunctionParam);

  SuperPlaneIntegralFluxVelocity3D<T> velFlux(sLattice, converter, superGeometry, center, normal, materials,
      BlockDataReductionMode::Discrete);

  /// === 5th Step: Definition of turbulent Statistics Objects ===

  SuperLatticePhysVelocity3D<T, DESCRIPTOR> sVel(sLattice, converter);
  SuperLatticeTimeAveragedF3D<T> sAveragedVel(sVel);

  /// === 4th Step: Main Loop with Timer ===
  Timer<double> timer(converter.getLatticeTime(maxPhysT), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for (size_t iT=0; iT<converter.getLatticeTime(maxPhysT); ++iT) {
    /// === 6th Step: Computation and Output of the Results ===
    getResults(sLattice, converter, iT, superGeometry, timer, sAveragedVel);

    if ( iT%converter.getLatticeTime(maxPhysT/fluxUpdates)==0 || iT == 0 ) {

      velFlux(output, input);
      T flux = output[0];
      T area = output[1];
      forceSol.updateAveVelocity(flux/area);
      sLattice.defineField<FORCE>(superGeometry, 1, forceSolScaled);
      sLattice.defineField<FORCE>(superGeometry, 2, forceSolScaled);
    }

    /// === 7th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

  }
  timer.stop();
  timer.printSummary();

  return 0;
}

