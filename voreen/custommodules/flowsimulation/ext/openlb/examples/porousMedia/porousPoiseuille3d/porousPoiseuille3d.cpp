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

/* porousPoiseuille3d.cpp:
 * This example examines a 3D Poseuille flow with porous media.
 * Two porous media LB methods can be used here:
 * Spaid and Phelan (doi:10.1063/1.869392), or
 * Guo and Zhao (doi:10.1103/PhysRevE.66.036304)
 */


#include "olb3D.h"
#include "olb3D.hh"

using namespace olb;
using namespace olb::descriptors;

typedef double T;

#define SPAID_PHELAN

T Kin = 1e-2;               // Permeability
T epsilon = 1.;             // Porosity (Spaid and Phelan can only handle epsilon=1)

T p0;                       // Initial pressure at inlet
T dp;                       // Pressure gradient
T mu;                       // Dynamic viscosity

#ifdef BGK
typedef D3Q19<> DESCRIPTOR;
#define DYNAMICS BGKdynamics
#endif
// Porous media
#ifdef SPAID_PHELAN
typedef D3Q19<POROSITY> DESCRIPTOR;
#define DYNAMICS PorousBGKdynamics
#elif defined GUO_ZHAO
typedef D3Q19<FORCE,EPSILON,K,NU,BODY_FORCE> DESCRIPTOR;
#define DYNAMICS GuoZhaoBGKdynamics
#endif

// Parameters for the simulation setup
const T length  = 2.;         // length of the pie
const T diameter  = 1.;       // diameter of the pipe
int N = 21;                   // resolution of the model
const T physU = 1.;           // physical velocity
const T Re = 1.;              // Reynolds number
const T physRho = 1.;         // physical density
const T tau = 0.8;            // lattice relaxation time
const T maxPhysT = 20.;       // max. simulation time in s, SI unit
const T residuum = 1e-5;      // residuum for the convergence check

// Scaled Parameters
const T radius  = diameter/2.;            // radius of the pipe
const T physInterval = 0.0125*maxPhysT;   // interval for the convergence check in s


template <typename T>
class PhysicalToLatticeVelocityF3D: public AnalyticalF3D<T,T> {
protected:
  AnalyticalF3D<T,T>* f;
  UnitConverter<T,DESCRIPTOR> converter;

public:
  PhysicalToLatticeVelocityF3D(AnalyticalF3D<T,T>* f_, UnitConverter<T,DESCRIPTOR> const& converter_)
    : AnalyticalF3D<T,T>(3), f(f_), converter(converter_) {};

  bool operator()(T output[], const T x[]) override
  {
    (*f)(output, x);
    for (int i=0; i<3; ++i) {
      output[i] = converter.getLatticeVelocity( output[i] );
    }
    return true;
  };
};


// Approximation of the modified Bessel function (doi:10.1088/1742-6596/1043/1/012003)
T besselApprox( T x )
{
  return util::cosh(x) / util::pow( 1 + 0.25*util::pow(x,2), 0.25 ) * ( 1 + 0.24273*util::pow(x,2) )/( 1 + 0.43023*util::pow(x,2) );
}

/// Functional to calculate velocity profile on pipe with porous media.
template <typename T>
class PorousPoiseuille3D : public AnalyticalF3D<T,T> {
protected:
  T K, mu, dp, radius;
  bool trunc;

public:
  PorousPoiseuille3D(T K_, T mu_, T dp_, T radius_)
    : AnalyticalF3D<T,T>(3), K(K_), mu(mu_), dp(dp_), radius(radius_) {};
  bool operator()(T output[], const T x[]) override
  {
    T r = util::sqrt(epsilon/K);
    T dist = util::sqrt( util::pow(x[1]-radius, 2.) + util::pow(x[2]-radius, 2.) );
    output[0] = K / mu * dp * ( 1. - besselApprox(r*dist) / besselApprox(r*radius)  );
    output[1] = 0.;
    output[2] = 0.;

    if ( dist > radius ) {
      output[0] = 0.;
    }

    return true;
  };
};

// Stores geometry information in form of material numbers
void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter,
                      SuperGeometry<T,3>& superGeometry )
{

  OstreamManager clout(std::cout, "prepareGeometry");

  clout << "Prepare Geometry ..." << std::endl;

  Vector<T, 3> center0(-converter.getPhysDeltaX() * 0.2, radius, radius);
  Vector<T, 3> center1(length, radius, radius);
  IndicatorCylinder3D<T> pipe(center0, center1, radius);

  superGeometry.rename(0, 2);

  superGeometry.rename(2, 1, pipe);

  Vector<T, 3> origin(0, radius, radius);
  Vector<T, 3> extend = origin;

  // Set material number for inflow
  origin[0] = -converter.getPhysDeltaX() * 2;
  extend[0] = converter.getPhysDeltaX() * 2;
  IndicatorCylinder3D<T> inflow(origin, extend, radius);
  superGeometry.rename(2, 3, 1, inflow);

  // Set material number for outflow
  origin[0] = length - 2 * converter.getPhysDeltaX();
  extend[0] = length + 2 * converter.getPhysDeltaX();
  IndicatorCylinder3D<T> outflow(extend, origin, radius);
  superGeometry.rename(2, 4, 1, outflow);

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice(SuperLattice<T, DESCRIPTOR>& sLattice,
                    UnitConverter<T, DESCRIPTOR>const& converter,
                    SuperGeometry<T,3>& superGeometry)
{

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=0 -->do nothing
  sLattice.defineDynamics<NoDynamics>(superGeometry, 0);

  // Material=1 -->bulk dynamics
  #ifdef GUO_ZHAO
  Dynamics<T,DESCRIPTOR>* dynamics = new DYNAMICS<T,DESCRIPTOR>(omega);
  sLattice.defineDynamics(superGeometry, 1, dynamics);
  #else
  sLattice.defineDynamics<DYNAMICS>(superGeometry, 1);
  #endif

  Vector<T, 3> center0(0, radius, radius);
  Vector<T, 3> center1(length, radius, radius);

  AnalyticalConst3D<T,T> zero(0.);
  AnalyticalConst3D<T,T> one(1.);

  T nu = (tau-0.5)/3.;
  T h = converter.getPhysDeltaX();
#ifdef SPAID_PHELAN
  T d = 1. - (h*h*nu*tau/Kin);
  clout << "Lattice Porosity: " << d << std::endl;
  clout << "Kmin: " << h*h*nu*tau << std::endl;
  if (Kin < h*h*nu*tau) {
    clout << "WARNING: Chosen K is too small!" << std::endl;
    exit(1);
  }
  AnalyticalConst3D<T,T> porosity(d);
  for (int i: {
         0,1,2,3,4
       }) {
    sLattice.defineField<POROSITY>(superGeometry, i, porosity);
  }
#elif defined GUO_ZHAO
  AnalyticalConst3D<T,T> Nu(nu);
  AnalyticalConst3D<T,T> k(Kin/(h*h));
  AnalyticalConst3D<T,T> eps(epsilon);
  for (int i: {
         0,1,2,3,4
       }) {
    sLattice.defineField<EPSILON>(superGeometry, i, eps);
    sLattice.defineField<NU>(superGeometry, i, Nu);
    sLattice.defineField<K>(superGeometry, i, k);
  }
#endif

  // Bouzidi
  sLattice.defineDynamics<NoDynamics>(superGeometry, 2);
  center0[0] -= 0.5*converter.getPhysDeltaX();
  center1[0] += 0.5*converter.getPhysDeltaX();
  IndicatorCylinder3D<T> pipe(center0, center1, radius);
  setBouzidiZeroVelocityBoundary<T,DESCRIPTOR>(sLattice, superGeometry, 2, pipe);
  // Interp
  //sLattice.defineDynamics<DYNAMICS>(superGeometry, 2);
  //setInterpolatedVelocityBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 2);


  // Material=3 --> bulk dynamics
  #ifdef GUO_ZHAO
  sLattice.defineDynamics(superGeometry, 3, dynamics);
  #else
  sLattice.defineDynamics<DYNAMICS>(superGeometry, 3);
  #endif
  setInterpolatedVelocityBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 3);

  // Material=4 --> bulk dynamics
  #ifdef GUO_ZHAO
  sLattice.defineDynamics(superGeometry, 4, dynamics);
  #else
  sLattice.defineDynamics<DYNAMICS>(superGeometry, 4);
  #endif
  setInterpolatedPressureBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 4);

  // Initial conditions
  // Pressure for Poiseuille flow with maximum velocity of charU at K->infty
  p0 = 4. * converter.getPhysViscosity() * converter.getCharPhysVelocity() * length / (radius * radius);
  T p0L = converter.getLatticePressure(p0);
  AnalyticalLinear3D<T, T> rho(-p0L / length * invCs2<T,DESCRIPTOR>(), 0, 0, p0L * invCs2<T,DESCRIPTOR>() + 1);

  dp = p0/length;
  mu = converter.getPhysViscosity()*converter.getPhysDensity();

  //CirclePoiseuille3D<T> uSol( {0., radius, radius}, {1, 0, 0}, converter.getCharPhysVelocity(), radius );
  PorousPoiseuille3D<T> uSol( Kin, mu, dp, radius );
  PhysicalToLatticeVelocityF3D<T> u(&uSol, converter);

  // Initialize all values of distribution functions to their local equilibrium
  for (int i: {
         0,1,2,3,4
       }) {
    sLattice.defineRhoU(superGeometry, i, rho, u);
    sLattice.iniEquilibrium(superGeometry, i, rho, u);
  }

  sLattice.setParameter<descriptors::OMEGA>(omega);

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

/// Compute error norms
void error( SuperGeometry<T,3>& superGeometry,
            SuperLattice<T, DESCRIPTOR>& sLattice,
            UnitConverter<T,DESCRIPTOR> const& converter,
            AnalyticalF3D<T,T>& uSol)
{

  OstreamManager clout( std::cout,"error" );


  int tmp[]= { };
  T result[2]= { };
  auto indicatorF = superGeometry.getMaterialIndicator(1);
  SuperLatticePhysVelocity3D<T,DESCRIPTOR> u( sLattice, converter );

  // Velocity error
  SuperAbsoluteErrorL1Norm3D<T> absVelocityErrorNormL1(u, uSol, indicatorF);
  absVelocityErrorNormL1(result, tmp);
  clout << "velocity-L1-error(abs)=" << result[0];
  SuperRelativeErrorL1Norm3D<T> relVelocityErrorNormL1(u, uSol, indicatorF);
  relVelocityErrorNormL1(result, tmp);
  clout << "; velocity-L1-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorL2Norm3D<T> absVelocityErrorNormL2(u, uSol, indicatorF);
  absVelocityErrorNormL2(result, tmp);
  clout << "velocity-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm3D<T> relVelocityErrorNormL2(u, uSol, indicatorF);
  relVelocityErrorNormL2(result, tmp);
  clout << "; velocity-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm3D<T> absVelocityErrorNormLinf(u, uSol, indicatorF);
  absVelocityErrorNormLinf(result, tmp);
  clout << "velocity-Linf-error(abs)=" << result[0];
  SuperRelativeErrorLinfNorm3D<T> relVelocityErrorNormLinf(u, uSol, indicatorF);
  relVelocityErrorNormLinf(result, tmp);
  clout << "; velocity-Linf-error(rel)=" << result[0] << std::endl;

  // Pressure error
  T p0L = converter.getLatticePressure(p0);
  AnalyticalLinear3D<T,T> pressureSol( -converter.getPhysPressure( p0L )/length, 0, 0, converter.getPhysPressure( p0L ) );
  SuperLatticePhysPressure3D<T,DESCRIPTOR> pressure( sLattice,converter );

  SuperAbsoluteErrorL1Norm3D<T> absPressureErrorNormL1(pressure, pressureSol, indicatorF);
  absPressureErrorNormL1(result, tmp);
  clout << "pressure-L1-error(abs)=" << result[0];
  SuperRelativeErrorL1Norm3D<T> relPressureErrorNormL1(pressure, pressureSol, indicatorF);
  relPressureErrorNormL1(result, tmp);
  clout << "; pressure-L1-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorL2Norm3D<T> absPressureErrorNormL2(pressure, pressureSol, indicatorF);
  absPressureErrorNormL2(result, tmp);
  clout << "pressure-L2-error(abs)=" << result[0];
  SuperRelativeErrorL2Norm3D<T> relPressureErrorNormL2(pressure, pressureSol, indicatorF);
  relPressureErrorNormL2(result, tmp);
  clout << "; pressure-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm3D<T> absPressureErrorNormLinf(pressure, pressureSol, indicatorF);
  absPressureErrorNormLinf(result, tmp);
  clout << "pressure-Linf-error(abs)=" << result[0];
  SuperRelativeErrorLinfNorm3D<T> relPressureErrorNormLinf(pressure, pressureSol, indicatorF);
  relPressureErrorNormLinf(result, tmp);
  clout << "; pressure-Linf-error(rel)=" << result[0] << std::endl;
}


// Output to console and files
void getResults( SuperLattice<T,DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                 SuperGeometry<T,3>& superGeometry, util::Timer<T>& timer, bool hasConverged )
{

  OstreamManager clout( std::cout,"getResults" );

  SuperVTMwriter3D<T> vtmWriter( "porousPoiseuille3d" );
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );

  //CirclePoiseuille3D<T> uSol( {0., radius, radius}, {1, 0, 0}, converter.getCharPhysVelocity(), radius );
  PorousPoiseuille3D<T> uSol( Kin, mu, dp, radius );
  SuperLatticeFfromAnalyticalF3D<T,DESCRIPTOR> uSolF( uSol, sLattice );
  vtmWriter.addFunctor( uSolF );

  const int vtmIter  = converter.getLatticeTime( maxPhysT/20. );
  const int statIter = converter.getLatticeTime( maxPhysT/20. );

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

  // Writes the vtm files and profile text file
  if ( iT%vtmIter==0 || hasConverged ) {
    vtmWriter.write( iT );
  }


  // Writes output on the console
  if ( iT%statIter==0 || hasConverged ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );

    // Calculate inflow and outflow flux
    std::vector<int> materials = { 1, 3, 4 };
    Vector<T,3> normal( 1, 0, 0 );
    auto mode = BlockDataReductionMode::Discrete;
    Vector<T,3> posInflow = superGeometry.getStatistics().getMinPhysR( 1 );
    Vector<T,3> posOutflow = superGeometry.getStatistics().getMaxPhysR( 1 );

    SuperPlaneIntegralFluxVelocity3D<T> vFluxIn( sLattice, converter,
        superGeometry, posInflow, normal, materials, mode );
    SuperPlaneIntegralFluxPressure3D<T> pFluxIn( sLattice, converter,
        superGeometry, posInflow, normal, materials, mode );
    SuperPlaneIntegralFluxVelocity3D<T> vFluxOut( sLattice, converter,
        superGeometry, posOutflow, normal, materials, mode );
    SuperPlaneIntegralFluxPressure3D<T> pFluxOut( sLattice, converter,
        superGeometry, posOutflow, normal, materials, mode );

    vFluxIn.print( "Inflow" );
    pFluxIn.print( "Inflow" );
    vFluxOut.print( "Outflow" );
    pFluxOut.print( "Outflow" );

    // Error norms
    AnalyticalFfromSuperF3D<T> intpolatePressure( pressure, true );

    T point1[3] = {0, radius, radius};
    T point2[3] = {0, radius, radius};

    point1[0] = length*0.5 - length*0.01;
    point2[0] = length*0.5 + length*0.01;

    T p1, p2;
    intpolatePressure( &p1,point1 );
    intpolatePressure( &p2,point2 );

    clout << "pressure1=" << p1;
    clout << "; pressure2=" << p2;

    T pressureDrop = p1-p2;
    clout << "; pressureDrop=" << pressureDrop << std::endl;

    error(superGeometry, sLattice, converter, uSol);

    // Gnuplot
    Gnuplot<T> gplot( "velocityProfile" );
    T uAnalytical[3] = {};
    T uNumerical[3] = {};
    AnalyticalFfromSuperF3D<T> intpolateVelocity( velocity, true );
    for (int i=0; i<101; i++) {
      T yInput = diameter*i/100.;
      T input[3] = {length*0.5, yInput, radius};
      uSol(uAnalytical, input);
      intpolateVelocity(uNumerical, input);
      gplot.setData( yInput, {uAnalytical[0], uNumerical[0]}, {"analytical","numerical"} );
    }

    // Create PNG file
    gplot.writePNG();
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
      clout<<"Usage: program [Resolution] [Permeability]" <<std::endl;
      clout<<"Default: Resolution=21, Permeability=1e-2" <<std::endl;
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
      std::cerr << "Permeabilty must be greater than 0" << std::endl;
      return 2;
    }
  }

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
    int {N},                  // resolution: number of voxels per charPhysL
    (T)   tau,                // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   diameter,           // charPhysLength: reference length of simulation geometry
    (T)   physU,              // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   diameter*physU/Re,  // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   physRho             // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  //converter.print();
  // Writes the converter log in a file
  converter.write("porousPoiseuille3d");


  // === 2nd Step: Prepare Geometry ===

  Vector<T, 3> center0(0, radius, radius);
  Vector<T, 3> center1(length + 0.5 * converter.getPhysDeltaX(), radius, radius);
  IndicatorCylinder3D<T> pipe(center0, center1, radius);
  IndicatorLayer3D<T> extendedDomain(pipe, converter.getPhysDeltaX());

  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = 2*singleton::mpi().getSize();
#else // ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = 6;
#endif // ifdef PARALLEL_MODE_MPI
  CuboidGeometry3D<T> cuboidGeometry(extendedDomain, converter.getPhysDeltaX(), noOfCuboids);

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

  // Instantiation of a superGeometry
  SuperGeometry<T,3> superGeometry(cuboidGeometry, loadBalancer);

  prepareGeometry(converter, superGeometry);

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice( superGeometry );

  //prepareLattice and setBoundaryConditions
  prepareLattice(sLattice, converter, superGeometry);

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
