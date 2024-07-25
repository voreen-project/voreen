/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2019 Fabian Klemens, Davide Dapelo, Mathias J. Krause
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

/* porousPoiseuille2d.cpp:
 * This example examines a 2D Poseuille flow with porous media.
 * Two porous media LB methods can be used here:
 * Spaid and Phelan (doi:10.1063/1.869392), or
 * Guo and Zhao (doi:10.1103/PhysRevE.66.036304)
 */


#include "olb2D.h"
#include "olb2D.hh"

using namespace olb;
using namespace olb::descriptors;

typedef double T;

//#define BGK
#define SPAID_PHELAN
//#define GUO_ZHAO


T Kin = 1e-2;               // Permeability
T epsilon = 1.;             // Porosity (Spaid and Phelan can only handle epsilon=1)

T p0;                       // Initial pressure at inlet
T dp;                       // Pressure gradient
T mu;                       // Dynamic viscosity

#ifdef BGK
typedef D2Q9<> DESCRIPTOR;
#define DYNAMICS BGKdynamics
#endif
// Porous media
#ifdef SPAID_PHELAN
typedef D2Q9<POROSITY> DESCRIPTOR;
#define DYNAMICS PorousBGKdynamics
#elif defined GUO_ZHAO
typedef D2Q9<FORCE,EPSILON,K,NU,BODY_FORCE> DESCRIPTOR;
#define DYNAMICS GuoZhaoBGKdynamics
#endif



// Functor to convert physical to lattice velocity
template <typename T>
class PhysicalToLatticeVelocityF2D: public AnalyticalF2D<T,T> {
protected:
  AnalyticalF2D<T,T>* f;
  UnitConverter<T,DESCRIPTOR> converter;

public:
  PhysicalToLatticeVelocityF2D(AnalyticalF2D<T,T>* f_, UnitConverter<T,DESCRIPTOR> const& converter_)
    : AnalyticalF2D<T,T>(2), f(f_), converter(converter_) {};

  bool operator()(T output[], const T x[]) override
  {
    (*f)(output, x);
    for (int i=0; i<2; ++i) {
      output[i] = converter.getLatticeVelocity( output[i] );
    }
    return true;
  };
};

/// Velocity profile in a pipe filled with isotropic porous media
template <typename T>
class PorousPoiseuille2D : public AnalyticalF2D<T,T> {
protected:
  T K, dp, mu, radius, wallOffset;

public:
  PorousPoiseuille2D(T K_, T dp_, T mu_, T radius_, T wallOffset_=0.)
    : AnalyticalF2D<T,T>(2), K(K_), dp(dp_), mu(mu_), radius(radius_),
      wallOffset(wallOffset_)
  { }

  bool operator()(T output[], const T x[]) override
  {
    T r = util::sqrt(epsilon/K);
    output[0] = dp / mu * K
     * (1. - util::cosh(r*(x[1] - radius)) / util::cosh(r*(radius-wallOffset)));
    output[1] = 0.;

    if ( x[1] < wallOffset || x[1] > 2.*radius - wallOffset ) {
      output[0] = 0.;
    }

    return true;
  }
};


typedef enum {bounceBack, local, interpolated} BoundaryType;

// Parameters for the simulation setup
BoundaryType boundaryType = interpolated;
const T lx  = 2.;             // length of the channel
const T ly  = 1.;             // height of the channel
int N = 50;                   // resolution of the model
const T Re = 1.;              // Reynolds number
const T maxPhysT = 20.;       // max. simulation time in s, SI unit
const T physInterval = 0.25;  // interval for the convergence check in s
const T residuum = 1e-5;      // residuum for the convergence check


// Stores geometry information in form of material numbers
void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter,
                      SuperGeometry<T,2>& superGeometry )
{

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( 0,2 );

  superGeometry.rename( 2,1,{1,1} );

  Vector<T,2> extend;
  Vector<T,2> origin;
  T physSpacing = converter.getPhysDeltaX();

  // Set material number for inflow
  extend[1] = ly;
  extend[0] = physSpacing / 2;
  origin[0] -= physSpacing / 4;
  IndicatorCuboid2D<T> inflow( extend, origin );
  superGeometry.rename( 2,3,1,inflow );

  // Set material number for outflow
  origin[0] = lx - physSpacing / 4;
  IndicatorCuboid2D<T> outflow( extend, origin );
  superGeometry.rename( 2,4,1,outflow );

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice( UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperLattice<T, DESCRIPTOR>& sLattice,
                     SuperGeometry<T,2>& superGeometry )
{

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=0 -->do nothing
  sLattice.defineDynamics<NoDynamics>(superGeometry, 0);

  // Material=1 -->bulk dynamics
  #ifdef GUO_ZHAO
  Dynamics<T,DESCRIPTOR>* bulkDynamics = new DYNAMICS<T,DESCRIPTOR>(omega);
  sLattice.defineDynamics(superGeometry, 1, bulkDynamics);
  #else
  sLattice.defineDynamics<DYNAMICS>(superGeometry, 1);
  #endif

  if (boundaryType == bounceBack) {
    sLattice.defineDynamics<BounceBack>(superGeometry, 2);
  }
  else {
    if (boundaryType == local) {
      setLocalVelocityBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 2);
    }
    else {
      setInterpolatedVelocityBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 2);
    }
  }

  // Material=3 -->bulk dynamics
  // Material=4 -->bulk dynamics
  #ifdef GUO_ZHAO
  sLattice.defineDynamics(superGeometry.getMaterialIndicator({3,4}), bulkDynamics);
  #else
  sLattice.defineDynamics<DYNAMICS>(superGeometry.getMaterialIndicator({3,4}));
  #endif

  if (boundaryType == local) {
    setLocalVelocityBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 3);
    setLocalPressureBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 4);
  }
  else {
    setInterpolatedVelocityBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 3);
    setInterpolatedPressureBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 4);
  }

#ifdef SPAID_PHELAN
  T tau = converter.getLatticeRelaxationTime();
  T nu = converter.getLatticeViscosity();
  T h = converter.getPhysDeltaX();
  T d = 1. - (h*h*nu*tau/Kin);
  clout << "Lattice Porosity: " << d << std::endl;
  clout << "Kmin: " << h*h*nu*tau << std::endl;
  if (Kin < h*h*nu*tau) {
    clout << "WARNING: Chosen K is too small!" << std::endl;
    exit(1);
  }
  AnalyticalConst2D<T,T> porosity( d );
  for (int i: {
         0,1,2,3,4
       }) {
    sLattice.defineField<POROSITY>(superGeometry, i, porosity);
  }
#elif defined GUO_ZHAO
  AnalyticalConst2D<T,T> eps( epsilon );
  AnalyticalConst2D<T,T> Nu( converter.getLatticeViscosity() );
  AnalyticalConst2D<T,T> k( Kin/util::pow(converter.getPhysDeltaX(), 2.) );
  for (int i: {
         0,1,2,3,4
       }) {
    sLattice.defineField<EPSILON>(superGeometry, i, eps);
    sLattice.defineField<NU>(superGeometry, i, Nu);
    sLattice.defineField<K>(superGeometry, i, k);
  }
#endif

  // Initial conditions
  // Pressure for Poiseuille flow with maximum velocity of charU at K->infty
  p0 = 8.*converter.getPhysViscosity()*converter.getCharPhysVelocity()*lx/( ly*ly );
  // Pressure for PorousPoiseuille with maximum velocity of charU for every permeability K
  //p0 = converter.getCharPhysVelocity() * converter.getPhysViscosity() * converter.getPhysDensity() / Kin * lx / (1. - 1./util::cosh(util::sqrt(1./Kin)*ly/2.));

  T p0L = converter.getLatticePressure(p0);
  AnalyticalLinear2D<T,T> rho( -p0L/lx*invCs2<T,DESCRIPTOR>(), 0, p0L*invCs2<T,DESCRIPTOR>()+1 );

  dp = p0/lx;
  mu = converter.getPhysViscosity()*converter.getPhysDensity();

  const T wallOffset = (boundaryType == bounceBack) ? 0.5 * converter.getPhysDeltaX() : 0.;
  PorousPoiseuille2D<T> uSol( Kin, dp, mu, ly/2., wallOffset );
  PhysicalToLatticeVelocityF2D<T> u( &uSol, converter );

  // Initialize all values of distribution functions to their local equilibrium
  for (int i: {
         0,1,2,3,4
       }) {
    sLattice.defineRhoU( superGeometry, i, rho, u );
    sLattice.iniEquilibrium( superGeometry, i, rho, u );
  }

  sLattice.setParameter<descriptors::OMEGA>(omega);

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Compute error norms
void error( SuperGeometry<T,2>& superGeometry,
            SuperLattice<T, DESCRIPTOR>& sLattice,
            UnitConverter<T,DESCRIPTOR> const& converter )
{

  OstreamManager clout( std::cout,"error" );

  int tmp[]= { };
  T result[2]= { };

  const T wallOffset = (boundaryType == bounceBack) ? 0.5 * converter.getPhysDeltaX() : 0.;
  PorousPoiseuille2D<T> uSol( Kin, dp, mu, ly/2., wallOffset );
  SuperLatticePhysVelocity2D<T,DESCRIPTOR> u( sLattice,converter );
  auto indicatorF = superGeometry.getMaterialIndicator(1);

  SuperAbsoluteErrorL1Norm2D<T> absVelocityErrorNormL1(u, uSol, indicatorF);
  absVelocityErrorNormL1(result, tmp);
  clout << "velocity-L1-error(abs)=" << result[0];
  SuperRelativeErrorL1Norm2D<T> relVelocityErrorNormL1(u, uSol, indicatorF);
  relVelocityErrorNormL1(result, tmp);
  clout << "; velocity-L1-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorL2Norm2D<T> absVelocityErrorNormL2(u, uSol, indicatorF);
  absVelocityErrorNormL2(result, tmp);
  clout << "velocity-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm2D<T> relVelocityErrorNormL2(u, uSol, indicatorF);
  relVelocityErrorNormL2(result, tmp);
  clout << "; velocity-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm2D<T> absVelocityErrorNormLinf(u, uSol, indicatorF);
  absVelocityErrorNormLinf(result, tmp);
  clout << "velocity-Linf-error(abs)=" << result[0];
  SuperRelativeErrorLinfNorm2D<T> relVelocityErrorNormLinf(u, uSol, indicatorF);
  relVelocityErrorNormLinf(result, tmp);
  clout << "; velocity-Linf-error(rel)=" << result[0] << std::endl;

  // pressure error
  T p0L = converter.getLatticePressure( p0 );
  AnalyticalLinear2D<T,T> pressureSol( -converter.getPhysPressure( p0L )/lx, 0, converter.getPhysPressure( p0L ) );
  SuperLatticePhysPressure2D<T,DESCRIPTOR> pressure( sLattice,converter );

  SuperAbsoluteErrorL1Norm2D<T> absPressureErrorNormL1(pressure, pressureSol, indicatorF);
  absPressureErrorNormL1(result, tmp);
  clout << "pressure-L1-error(abs)=" << result[0];
  SuperRelativeErrorL1Norm2D<T> relPressureErrorNormL1(pressure, pressureSol, indicatorF);
  relPressureErrorNormL1(result, tmp);
  clout << "; pressure-L1-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorL2Norm2D<T> absPressureErrorNormL2(pressure, pressureSol, indicatorF);
  absPressureErrorNormL2(result, tmp);
  clout << "pressure-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm2D<T> relPressureErrorNormL2(pressure, pressureSol, indicatorF);
  relPressureErrorNormL2(result, tmp);
  clout << "; pressure-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm2D<T> absPressureErrorNormLinf(pressure, pressureSol, indicatorF);
  absPressureErrorNormLinf(result, tmp);
  clout << "pressure-Linf-error(abs)=" << result[0];
  SuperRelativeErrorLinfNorm2D<T> relPressureErrorNormLinf(pressure, pressureSol, indicatorF);
  relPressureErrorNormLinf(result, tmp);
  clout << "; pressure-Linf-error(rel)=" << result[0] << std::endl;
}

// Output to console and files
void getResults( SuperLattice<T,DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                 SuperGeometry<T,2>& superGeometry, util::Timer<T>& timer, bool hasConverged )
{

  OstreamManager clout( std::cout,"getResults" );

  SuperVTMwriter2D<T> vtmWriter( "porousPoiseuille2d" );
  SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure2D<T, DESCRIPTOR> pressure( sLattice, converter );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );

  const T wallOffset = (boundaryType == bounceBack) ? 0.5 * converter.getPhysDeltaX() : 0.;
  PorousPoiseuille2D<T> uSol( Kin, dp, mu, ly/2., wallOffset );
  SuperLatticeFfromAnalyticalF2D<T,DESCRIPTOR> uSolF( uSol, sLattice);
  vtmWriter.addFunctor( uSolF );

  const int vtmIter  = converter.getLatticeTime( maxPhysT/20. );
  const int statIter = converter.getLatticeTime( maxPhysT/20. );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry2D<T, DESCRIPTOR> geometry( sLattice, superGeometry );
    SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank2D<T, DESCRIPTOR> rank( sLattice );
    superGeometry.rename( 0,2 );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );

    vtmWriter.createMasterFile();
  }

  // Writes the vtm files and profile text file
  if ( iT%vtmIter==0 || hasConverged ) {
    vtmWriter.write( iT );

    SuperEuklidNorm2D<T, DESCRIPTOR> normVel( velocity );
    BlockReduction2D2D<T> planeReduction( normVel, 600, BlockDataSyncMode::ReduceOnly );
    // write output as JPEG
    heatmap::write(planeReduction, iT);
  }

  if ( hasConverged ) {
    Gnuplot<T> gplot( "centerVelocity" );
    T Ly = converter.getLatticeLength( ly );
    for ( int iY=0; iY<=Ly; ++iY ) {
      T dx = 1. / T(converter.getResolution());
      T point[2]= {T(),T()};
      point[0] = lx/2.;
      point[1] = ( T )iY/Ly;

      const T wallOffset = (boundaryType == bounceBack) ? 0.5 * converter.getPhysDeltaX() : 0.;
      PorousPoiseuille2D<T> uSol( Kin, dp, mu, ly/2., wallOffset );
      T analytical[2] = {T(),T()};
      uSol( analytical,point );
      SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity( sLattice, converter );
      AnalyticalFfromSuperF2D<T> intpolateVelocity( velocity, true );
      T numerical[2] = {T(),T()};
      intpolateVelocity( numerical,point );
      gplot.setData( iY*dx, {analytical[0],numerical[0]}, {"analytical","numerical"} );
    }
    // Create PNG file
    gplot.writePNG();
  }

  // Writes output on the console
  if ( iT%statIter==0 || hasConverged ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );

    // Error norms
    error( superGeometry, sLattice, converter );
  }
}

int main( int argc, char* argv[] )
{

  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );

  if (argc > 1) {
    if (argv[1][0]=='-'&&argv[1][1]=='h') {
      OstreamManager clout( std::cout,"help" );
      clout<<"Usage: program [Resolution] [Permeability] [BoundaryType]"<<std::endl;
      clout<<"BoundaryType: 0=bounceBack, 1=local, 2=interpolated"<<std::endl;
      clout<<"Default: Resolution=50, Permeability=1e-2, BoundaryTpe=interpolated"<<std::endl;
      return 0;
    }
  }

  if (argc > 1) {
    N = atoi(argv[1]);
    if (N < 1) {
      std::cerr << "Fluid domain is too small" << std::endl;
      return 1;
    }
  }

  if (argc > 2) {
    Kin = atof(argv[2]);
    if (Kin < 0) {
      std::cerr << "Permeability must be non-negative" << std::endl;
      return 2;
    }
  }

  if (argc > 3) {
    int boundaryTypeNumber = atoi(argv[3]);
    if (boundaryTypeNumber < 0 || boundaryTypeNumber > (int) interpolated) {
      std::cerr << "Unknown boundary type" << std::endl;
      return 3;
    }
    boundaryType = (BoundaryType) boundaryTypeNumber;
  }

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
    int {N},     // resolution: number of voxels per charPhysL
    (T)   0.8,   // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   1,     // charPhysLength: reference length of simulation geometry
    (T)   1,     // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   1./Re, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.     // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("porousPoiseuille2d");


  // === 2nd Step: Prepare Geometry ===
  Vector<T,2> extend( lx, ly );
  Vector<T,2> origin;
  IndicatorCuboid2D<T> cuboid( extend, origin );

  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif
  CuboidGeometry2D<T> cuboidGeometry( cuboid, converter.getConversionFactorLength(), noOfCuboids );

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  const int overlap = (boundaryType == bounceBack) ? 2 : 3;
  SuperGeometry<T,2> superGeometry( cuboidGeometry, loadBalancer, overlap );

  prepareGeometry( converter, superGeometry );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice( superGeometry );

  //prepareLattice and setBoundaryConditions
  prepareLattice( converter, sLattice, superGeometry );

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  util::ValueTracer<T> converge( converter.getLatticeTime( physInterval ), residuum );
  timer.start();

  for ( std::size_t iT = 0; iT < converter.getLatticeTime( maxPhysT ); ++iT ) {
    if ( converge.hasConverged() ) {
      clout << "Simulation converged." << std::endl;
      getResults( sLattice, converter, iT, superGeometry, timer, converge.hasConverged() );

      break;
    }

    // === 5th Step: Definition of Initial and Boundary Conditions ===
    // in this application no boundary conditions have to be adjusted

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, iT, superGeometry, timer, converge.hasConverged()  );
    converge.takeValue( sLattice.getStatistics().getAverageEnergy(), true );
  }

  timer.stop();
  timer.printSummary();
}
