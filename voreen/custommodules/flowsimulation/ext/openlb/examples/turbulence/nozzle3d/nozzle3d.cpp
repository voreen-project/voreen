/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2015 Mathias J. Krause, Patrick Nathan
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

/* nozzle3d.cpp:
 * This example examines a turbulent flow in a nozzle injection tube. At the
 * main inlet, either a block profile or a power 1/7 profile is imposed as a
 * Dirchlet velocity boundary condition, whereas at the outlet a
 * Dirichlet pressure condition is set by p=0 (i.e. rho=1).
 *
 * The example shows the usage of turbulent models.
 *
 * The results of a simular simulation setup are publish in TODO
 */

#include "olb3D.h"
#include "olb3D.hh"

#include <chrono>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = float;

// Choose your turbulent model of choice
//#define USE_RLB
#define USE_SMAGORINSKY //default
//#define USE_CONSISTENT_STRAIN_SMAGORINSKY
//#define USE_SHEAR_SMAGORINSKY
//#define USE_KRAUSE

#ifdef USE_SHEAR_SMAGORINSKY
using DESCRIPTOR = D3Q19<AV_SHEAR>;
#else
using DESCRIPTOR = D3Q19<>;
#endif

// seed the rng with time if SEED_WITH_TIME is set, otherwise just use a fixed seed.
#if defined SEED_WITH_TIME
const auto seed = std::chrono::system_clock().now().time_since_epoch().count();
#else
const auto seed = 0x1337533DAAAAAAAA;
#endif

// Parameters for the simulation setup
const int N = 3;                 // resolution of the model, for RLB N>=5, others N>2, but uneven N>=5 recommended
const int inflowProfileMode = 0; // block profile (mode=0), power profile (mode=1)
const T maxPhysT = 200.;         // max. simulation time in s, SI unit

template <typename T, typename _DESCRIPTOR>
class TurbulentVelocity3D : public AnalyticalF3D<T,T> {
private:
  // block profile (mode=0), power profile (mode=1)
  int _mode;
  T _u0;
  T _nu;
  T _charL;
  T _dx;

  std::default_random_engine _generator;

public:
  TurbulentVelocity3D(UnitConverter<T,_DESCRIPTOR> const& converter, int mode)
    : AnalyticalF3D<T,T>(3)
    , _mode(mode)
    , _u0(converter.getCharLatticeVelocity())
    , _nu(converter.getPhysViscosity())
    , _charL(converter.getCharPhysLength())
    , _dx(converter.getConversionFactorLength())
    , _generator(seed)
  {
    this->getName() = "turbulentVelocity3d";
  };

  bool operator()( T output[], const BaseType<T> input[] ) override
  {
    T y = input[1];
    T z = input[2];

    T a = -1., b = 1.;
    std::uniform_real_distribution<T> distribution(a, b);
    T nRandom = distribution(_generator);

    // power profile inititalization
    if ( _mode==1 ) {
      T obst_y = 5.5+_dx;
      T obst_z = 5.5+_dx;
      T obst_r = 0.5;

      T B      = 5.5;
      T kappa  = 0.4;
      T ReTau  = 183.6;

      T u_calc = _u0/7.*( 2.*_nu*ReTau/( _charL*kappa )*util::log( util::fabs( 2.*ReTau/_charL*( obst_r - util::sqrt( util::pow( y - obst_y, 2. )
               + util::pow( z - obst_z, 2. ) ) )*1.5*( 1 + util::sqrt( util::pow( y - obst_y, 2. )
               + util::pow( z - obst_z, 2. ) )/obst_r )/( 1 + 2.*util::pow( util::sqrt( util::pow( y - obst_y, 2. )
               + util::pow( z - obst_z, 2. ) )/obst_r, 2. ) ) ) + B ) );

      output[0] = u_calc + 0.15*_u0*nRandom;
      output[1] = 0.15*_u0*nRandom;
      output[2] = 0.15*_u0*nRandom;
    // block profile inititalization
    } else {
      output[0] = _u0 + 0.05*_u0*nRandom;
      output[1] = 0.05*_u0*nRandom;
      output[2] = 0.05*_u0*nRandom;
    }

    return true;
  };
};


void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter, IndicatorF3D<T>& indicator, SuperGeometry<T,3>& superGeometry )
{

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  // Sets material number for fluid and boundary
  superGeometry.rename( 0,2,indicator );

  Vector<T,3> origin( T(),
                      5.5*converter.getCharPhysLength()+converter.getConversionFactorLength(),
                      5.5*converter.getCharPhysLength()+converter.getConversionFactorLength() );

  Vector<T,3> extend( 4.*converter.getCharPhysLength()+5*converter.getConversionFactorLength(),
                      5.5*converter.getCharPhysLength()+converter.getConversionFactorLength(),
                      5.5*converter.getCharPhysLength()+converter.getConversionFactorLength() );

  IndicatorCylinder3D<T> inletCylinder( extend, origin, converter.getCharPhysLength() );
  superGeometry.rename( 2,1,inletCylinder );

  origin[0]=4.*converter.getCharPhysLength();
  origin[1]=5.5*converter.getCharPhysLength()+converter.getConversionFactorLength();
  origin[2]=5.5*converter.getCharPhysLength()+converter.getConversionFactorLength();

  extend[0]=80.*converter.getCharPhysLength();
  extend[1]=5.5*converter.getCharPhysLength()+converter.getConversionFactorLength();
  extend[2]=5.5*converter.getCharPhysLength()+converter.getConversionFactorLength();

  IndicatorCylinder3D<T> injectionTube( extend, origin, 5.5*converter.getCharPhysLength() );
  superGeometry.rename( 2,1,injectionTube );

  origin[0]=converter.getConversionFactorLength();
  origin[1]=5.5*converter.getCharPhysLength()+converter.getConversionFactorLength();
  origin[2]=5.5*converter.getCharPhysLength()+converter.getConversionFactorLength();

  extend[0]=T();
  extend[1]=5.5*converter.getCharPhysLength()+converter.getConversionFactorLength();
  extend[2]=5.5*converter.getCharPhysLength()+converter.getConversionFactorLength();

  IndicatorCylinder3D<T> cylinderIN( extend, origin, converter.getCharPhysLength() );
  superGeometry.rename( 1,3,cylinderIN );


  origin[0]=80.*converter.getCharPhysLength()-converter.getConversionFactorLength();
  origin[1]=5.5*converter.getCharPhysLength()+converter.getConversionFactorLength();
  origin[2]=5.5*converter.getCharPhysLength()+converter.getConversionFactorLength();

  extend[0]=80.*converter.getCharPhysLength();
  extend[1]=5.5*converter.getCharPhysLength()+converter.getConversionFactorLength();
  extend[2]=5.5*converter.getCharPhysLength()+converter.getConversionFactorLength();

  IndicatorCylinder3D<T> cylinderOUT( extend, origin, 5.5*converter.getCharPhysLength() );
  superGeometry.rename( 1,4,cylinderOUT );

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice( SuperLattice<T,DESCRIPTOR>& sLattice,
                     UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperGeometry<T,3>& superGeometry )
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=0 -->do nothing
  sLattice.defineDynamics<NoDynamics>(superGeometry, 0);

  // Material=1 -->bulk dynamics
  // Material=3 -->bulk dynamics (inflow)
  // Material=4 -->bulk dynamics (outflow)
  auto bulkIndicator = superGeometry.getMaterialIndicator({1, 3, 4});
#if defined(USE_RLB)
  using BulkDynamics = RLBdynamics<T,DESCRIPTOR>;
  sLattice.defineDynamics<BulkDynamics>(bulkIndicator);
#elif defined(USE_SMAGORINSKY)
  using BulkDynamics = SmagorinskyBGKdynamics<T,DESCRIPTOR>;
  sLattice.defineDynamics<BulkDynamics>(bulkIndicator);
  sLattice.setParameter<collision::LES::Smagorinsky>(0.15);
#elif defined(USE_SHEAR_SMAGORINSKY)
  using BulkDynamics = ShearSmagorinskyBGKdynamics<T,DESCRIPTOR>;
  sLattice.defineDynamics<BulkDynamics>(bulkIndicator);
  sLattice.setParameter<collision::LES::Smagorinsky>(0.15);
#elif defined(USE_KRAUSE)
  using BulkDynamics = KrauseBGKdynamics<T,DESCRIPTOR>;
  sLattice.defineDynamics<BulkDynamics>(bulkIndicator);
  sLattice.setParameter<collision::LES::Smagorinsky>(0.15);
#else // USE_CONSISTENT_STRAIN_SMAGORINSKY
  using BulkDynamics = ConStrainSmagorinskyBGKdynamics<T,DESCRIPTOR>;
  sLattice.defineDynamics<BulkDynamics>(bulkIndicator);
  sLattice.setParameter<collision::LES::Smagorinsky>(0.15);
#endif

  // Material=2 -->bounce back
  sLattice.defineDynamics<BounceBack>(superGeometry, 2);

  setInterpolatedVelocityBoundary(sLattice, omega, superGeometry, 3);
  setInterpolatedPressureBoundary(sLattice, omega, superGeometry, 4);

  sLattice.setParameter<descriptors::OMEGA>(omega);

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues(UnitConverter<T,DESCRIPTOR> const&converter,
                       SuperLattice<T,DESCRIPTOR>& lattice,
                       SuperGeometry<T,3>& superGeometry,
                       std::size_t iT)
{
  if (iT == 0) {
    AnalyticalConst3D<T,T> rhoF(1);
    Vector<T,3> velocity{};
    AnalyticalConst3D<T,T> uF(velocity);

    lattice.defineRhoU(superGeometry.getMaterialIndicator({1, 2, 3, 4}), rhoF, uF);
    lattice.iniEquilibrium(superGeometry.getMaterialIndicator({1, 2, 3, 4}), rhoF, uF);

    lattice.initialize();
  }

  const auto maxStartT =  converter.getLatticeTime(5);
  const auto startIterT = converter.getLatticeTime(0.05);

  if (iT < maxStartT && iT % startIterT == 0) {
    auto uSol = std::shared_ptr<AnalyticalF3D<T,T>>(new TurbulentVelocity3D<T,DESCRIPTOR>(converter, inflowProfileMode));

    PolynomialStartScale<T,std::size_t> scale(maxStartT, 1);
    T frac{};
    scale(&frac, &iT);

    lattice.defineU(superGeometry, 3, *(frac * uSol));

    lattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
}

void getResults( SuperLattice<T, DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR> const& converter, std::size_t iT,
                 SuperGeometry<T,3>& superGeometry, util::Timer<T>& timer )
{
  OstreamManager clout( std::cout,"getResults" );
  SuperVTMwriter3D<T> vtmWriter( "nozzle3d" );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry3D<T, DESCRIPTOR> geometry( sLattice, superGeometry );
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( sLattice );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }

  if (iT % converter.getLatticeTime(0.5) == 0) {
    timer.update(iT);
    timer.printStep();
    sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
    if (std::isnan(sLattice.getStatistics().getAverageRho())) {
      std::exit(-1);
    }
  }

  // Writes the vtk files
  if (iT % converter.getLatticeTime(2) == 0) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    // Create the data-reading functors...
    SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
    SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );
    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( pressure );
    vtmWriter.write( iT );

    SuperEuklidNorm3D<T, DESCRIPTOR> normVel( velocity );
    BlockReduction3D2D<T> planeReduction( normVel, {0, 1, 0} );
    // write output as JPEG
    heatmap::write(planeReduction, iT);
  }
}



int main( int argc, char* argv[] )
{

  // === 1st Step: Initialization ===

  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );
  // display messages from every single mpi process
  // clout.setMultiOutput(true);

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
    int {N},        // resolution: number of voxels per charPhysL
    (T)   0.500018, // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   1,        // charPhysLength: reference length of simulation geometry
    (T)   1,        // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   0.0002,   // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.0       // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("nozzle3d");

  Vector<T,3> origin;
  Vector<T,3> extend( 80.*converter.getCharPhysLength(), 11.*converter.getCharPhysLength()+2.*converter.getConversionFactorLength(), 11.*converter.getCharPhysLength()+2.*converter.getConversionFactorLength() );

  IndicatorCuboid3D<T> cuboid( extend,origin );

  CuboidGeometry3D<T> cuboidGeometry( cuboid, converter.getConversionFactorLength(), singleton::mpi().getSize() );
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // === 2nd Step: Prepare Geometry ===

  SuperGeometry<T,3> superGeometry( cuboidGeometry, loadBalancer );
  prepareGeometry( converter, cuboid, superGeometry );

  // === 3rd Step: Prepare Lattice ===

  SuperLattice<T, DESCRIPTOR> sLattice( superGeometry );

  prepareLattice( sLattice, converter, superGeometry );

  // === 4th Step: Main Loop with Timer ===

  setBoundaryValues( converter, sLattice, superGeometry, 0 );

  util::Timer<T> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for (std::size_t iT = 0; iT <= converter.getLatticeTime(maxPhysT); ++iT) {
    // === 5ath Step: Apply filter
#ifdef ADM
    SuperLatticeADM3D<T, DESCRIPTOR> admF(sLattice, 0.01, 2);
    admF.execute( superGeometry, 1 );
#endif
#if defined(USE_SHEAR_SMAGORINSKY)
    sLattice.setParameter<descriptors::LATTICE_TIME>(iT);
#endif

    setBoundaryValues(converter, sLattice, superGeometry, iT);

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults(sLattice, converter, iT, superGeometry, timer);
  }

  timer.stop();
  timer.printSummary();
}
