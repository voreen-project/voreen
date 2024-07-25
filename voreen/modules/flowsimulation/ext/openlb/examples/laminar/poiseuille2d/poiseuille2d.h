/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2007, 2012 Jonas Latt, Mathias J. Krause
 *  Vojtech Cvrcek, Peter Weisbrod
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

#include "olb2D.h"
#include "olb2D.hh"


using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = double;

//#define ENABLE_MRT
#ifdef ENABLE_MRT
using DESCRIPTOR = D2Q9<tag::MRT,FORCE>;
using BulkDynamics       = MRTdynamics<T,DESCRIPTOR>;
using ForcedBulkDynamics = ForcedMRTdynamics<T,DESCRIPTOR>;
#else
using DESCRIPTOR = D2Q9<FORCE>;
using BulkDynamics       = BGKdynamics<T,DESCRIPTOR>;
using ForcedBulkDynamics = ForcedBGKdynamics<T,DESCRIPTOR>;
#endif

typedef enum {forced, nonForced} FlowType;

typedef enum {bounceBack, local, interpolated, freeSlip, partialSlip} BoundaryType;


// Parameters for the simulation setup
FlowType flowType = forced;
BoundaryType boundaryType = interpolated;
const T lx  = 2.;             // length of the channel
const T ly  = 1.;             // height of the channel
const T Re = 10.;             // Reynolds number
const T maxPhysT = 30.;       // max. simulation time in s, SI unit
const T physInterval = 0.25;  // interval for the convergence check in s
const T residuum = 1e-9;      // residuum for the convergence check
const T tuner = 0.99;         // for partialSlip only: 0->bounceBack, 1->freeSlip


// variables for eoc analysis
T velocityL1AbsError = 0;
T velocityL2AbsError = 0;
T strainRateL1AbsError = 0;
T strainRateL2AbsError = 0;
T pressureL1AbsError = 0;
T pressureL2AbsError = 0;
T velocityLinfAbsError = 0;
T pressureLinfAbsError = 0;
T strainRateLinfAbsError = 0;

// Stores geometry information in form of material numbers
void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter,
                      SuperGeometry<T,2>& superGeometry )
{

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( 0,2 );
  superGeometry.rename( 2,1,{1,1} );

  if (flowType == nonForced) {
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
  }

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

  if (flowType == forced) {
    sLattice.defineDynamics<ForcedBulkDynamics>(superGeometry, 1);
  }
  else {
    sLattice.defineDynamics<BulkDynamics>(superGeometry, 1);
  }

  // Material=0 -->do nothing
  sLattice.defineDynamics<NoDynamics>(superGeometry, 0);

  if (boundaryType == bounceBack) {
    sLattice.defineDynamics<BounceBack>(superGeometry, 2);
  }
  else if (boundaryType == freeSlip) {
    sLattice.defineDynamics<NoDynamics>(superGeometry, 2);
    setSlipBoundary(sLattice, superGeometry, 2);
  }
  else if (boundaryType == partialSlip) {
    sLattice.defineDynamics<NoDynamics>(superGeometry, 2);
    setPartialSlipBoundary(sLattice, tuner, superGeometry, 2);
  }
  else {
    if (boundaryType == local) {
      setLocalVelocityBoundary(sLattice, omega, superGeometry, 2);
    }
    else {
      setInterpolatedVelocityBoundary(sLattice, omega, superGeometry, 2);
    }
  }

  if (flowType == nonForced) {
    if (boundaryType == local) {
      setLocalVelocityBoundary(sLattice, omega, superGeometry, 3);
      setLocalPressureBoundary(sLattice, omega, superGeometry, 4);
    }
    else {
      setInterpolatedVelocityBoundary(sLattice, omega, superGeometry, 3);
      setInterpolatedPressureBoundary(sLattice, omega, superGeometry, 4);
    }
  }

  // Initial conditions
  T Lx = converter.getLatticeLength( lx );
  T Ly = converter.getLatticeLength( ly );
  if (boundaryType == bounceBack) {
    Lx -= T(1);
    Ly -= T(1);
  }

  if (flowType == forced) {
    std::vector<T> poiseuilleForce( 2,T() );
    poiseuilleForce[0]
     = 8.*converter.getLatticeViscosity()
      * converter.getCharLatticeVelocity() / ( Ly*Ly );
    AnalyticalConst2D<T,T> force( poiseuilleForce );

    // Initialize force
    sLattice.defineField<FORCE>(superGeometry, 1, force);
    sLattice.defineField<FORCE>(superGeometry, 2, force);
  }
  else {
    T p0 = 8. * converter.getLatticeViscosity()
      * converter.getCharLatticeVelocity()*Lx/( Ly*Ly );
    AnalyticalLinear2D<T,T> rho(
      -p0/lx*invCs2<T,DESCRIPTOR>(), 0, p0*invCs2<T,DESCRIPTOR>()+1 );

    const T maxVelocity = converter.getCharLatticeVelocity();
    const T radius = (boundaryType == bounceBack) ?
      T(0.5) * (ly - converter.getPhysDeltaX()) : T(0.5) * ly;
    std::vector<T> axisPoint( 2,T() );
    axisPoint[0] = lx/2.;
    axisPoint[1] = ly/2.;
    std::vector<T> axisDirection( 2,T() );
    axisDirection[0] = 1;
    axisDirection[1] = 0;
    Poiseuille2D<T> u( axisPoint, axisDirection, maxVelocity, radius );

    std::vector<T> zero(2, T());
    AnalyticalConst2D<T, T> u0(zero);

    // Initialize all values of distribution functions to their local equilibrium
    sLattice.defineRhoU(superGeometry, 0, rho, u0);
    sLattice.iniEquilibrium(superGeometry, 0, rho, u0);
    sLattice.defineRhoU( superGeometry, 1, rho, u );
    sLattice.iniEquilibrium( superGeometry, 1, rho, u );
    sLattice.defineRhoU( superGeometry, 2, rho, u );
    sLattice.iniEquilibrium( superGeometry, 2, rho, u );
    sLattice.defineRhoU( superGeometry, 3, rho, u );
    sLattice.iniEquilibrium( superGeometry, 3, rho, u );
    sLattice.defineRhoU( superGeometry, 4, rho, u );
    sLattice.iniEquilibrium( superGeometry, 4, rho, u );
  }

  sLattice.template setParameter<descriptors::OMEGA>(omega);

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

  // velocity error
  const T maxVelocity = converter.getCharPhysVelocity();
  const T radius = (boundaryType == bounceBack) ?
    T(0.5) * (ly - converter.getPhysDeltaX()) : T(0.5) * ly;
  std::vector<T> axisPoint( 2,T() );
  axisPoint[0] = lx/2.;
  axisPoint[1] = ly/2.;
  std::vector<T> axisDirection( 2,T() );
  axisDirection[0] = 1;
  axisDirection[1] = 0;
  Poiseuille2D<T> uSol( axisPoint, axisDirection, maxVelocity, radius );
  SuperLatticePhysVelocity2D<T,DESCRIPTOR> u( sLattice,converter );
  auto indicatorF = superGeometry.getMaterialIndicator(1);

  SuperAbsoluteErrorL1Norm2D<T> absVelocityErrorNormL1(u, uSol, indicatorF);
  absVelocityErrorNormL1(result, tmp);
  clout << "velocity-L1-error(abs)=" << result[0];
  velocityL1AbsError = result[0];
  SuperRelativeErrorL1Norm2D<T> relVelocityErrorNormL1(u, uSol, indicatorF);
  relVelocityErrorNormL1(result, tmp);
  clout << "; velocity-L1-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorL2Norm2D<T> absVelocityErrorNormL2(u, uSol, indicatorF);
  absVelocityErrorNormL2(result, tmp);
  clout << "velocity-L2-error(abs)=" << result[0];
  velocityL2AbsError = result[0];
  SuperRelativeErrorL2Norm2D<T> relVelocityErrorNormL2(u, uSol, indicatorF);
  relVelocityErrorNormL2(result, tmp);
  clout << "; velocity-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm2D<T> absVelocityErrorNormLinf(u, uSol, indicatorF);
  absVelocityErrorNormLinf(result, tmp);
  clout << "velocity-Linf-error(abs)=" << result[0];
  velocityLinfAbsError = result[0];
  SuperRelativeErrorLinfNorm2D<T> relVelocityErrorNormLinf(u, uSol, indicatorF);
  relVelocityErrorNormLinf(result, tmp);
  clout << "; velocity-Linf-error(rel)=" << result[0] << std::endl;

  // strainRate error
  PoiseuilleStrainRate2D<T,T,DESCRIPTOR> sSol( converter, T(2) * radius );
  SuperLatticePhysStrainRate2D<T,DESCRIPTOR> s( sLattice,converter );

  SuperAbsoluteErrorL1Norm2D<T> absStrainRateErrorNormL1(s, sSol, indicatorF);
  absStrainRateErrorNormL1(result, tmp);
  clout << "strainRate-L1-error(abs)=" << result[0];
  strainRateL1AbsError = result[0];
  SuperRelativeErrorL1Norm2D<T> relStrainRateErrorNormL1(s, sSol, indicatorF);
  relStrainRateErrorNormL1(result, tmp);
  clout << "; strainRate-L1-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorL2Norm2D<T> absStrainRateErrorNormL2(s, sSol, indicatorF);
  absStrainRateErrorNormL2(result, tmp);
  clout << "strainRate-L2-error(abs)=" << result[0];
  strainRateL2AbsError = result[0];
  SuperRelativeErrorL2Norm2D<T> relStrainRateErrorNormL2(s, sSol, indicatorF);
  relStrainRateErrorNormL2(result, tmp);
  clout << "; strainRate-L2-error(rel)=" << result[0] << std::endl;

  SuperAbsoluteErrorLinfNorm2D<T> absStrainRateErrorNormLinf(s, sSol, indicatorF);
  absStrainRateErrorNormLinf(result, tmp);
  clout << "strainRate-Linf-error(abs)=" << result[0];
  strainRateLinfAbsError = result[0];
  SuperRelativeErrorLinfNorm2D<T> relStrainRateErrorNormLinf(s, sSol, indicatorF);
  relStrainRateErrorNormLinf(result, tmp);
  clout << "; strainRate-Linf-error(rel)=" << result[0] << std::endl;

  if (flowType == nonForced) {
    // pressure error
    int Lx = converter.getLatticeLength( lx );
    int Ly = converter.getLatticeLength( ly );
    T p0 = 8.*converter.getLatticeViscosity()*converter.getCharLatticeVelocity()*Lx/T( Ly*Ly );
    AnalyticalLinear2D<T,T> pressureSol( -converter.getPhysPressure( p0 )/lx, 0, converter.getPhysPressure( p0 ) );
    SuperLatticePhysPressure2D<T,DESCRIPTOR> pressure( sLattice,converter );

    SuperAbsoluteErrorL1Norm2D<T> absPressureErrorNormL1(pressure, pressureSol, indicatorF);
    absPressureErrorNormL1(result, tmp);
    clout << "pressure-L1-error(abs)=" << result[0];
    pressureL1AbsError = result[0];
    SuperRelativeErrorL1Norm2D<T> relPressureErrorNormL1(pressure, pressureSol, indicatorF);
    relPressureErrorNormL1(result, tmp);
    clout << "; pressure-L1-error(rel)=" << result[0] << std::endl;

    SuperAbsoluteErrorL2Norm2D<T> absPressureErrorNormL2(pressure, pressureSol, indicatorF);
    absPressureErrorNormL2(result, tmp);
    clout << "pressure-L2-error(abs)=" << result[0];
    pressureL2AbsError = result[0];
    SuperRelativeErrorL2Norm2D<T> relPressureErrorNormL2(pressure, pressureSol, indicatorF);
    relPressureErrorNormL2(result, tmp);
    clout << "; pressure-L2-error(rel)=" << result[0] << std::endl;

    SuperAbsoluteErrorLinfNorm2D<T> absPressureErrorNormLinf(pressure, pressureSol, indicatorF);
    absPressureErrorNormLinf(result, tmp);
    clout << "pressure-Linf-error(abs)=" << result[0];
    pressureLinfAbsError = result[0];
    SuperRelativeErrorLinfNorm2D<T> relPressureErrorNormLinf(pressure, pressureSol, indicatorF);
    relPressureErrorNormLinf(result, tmp);
    clout << "; pressure-Linf-error(rel)=" << result[0] << std::endl;
  }
}

// Output to console and files
void getResults( SuperLattice<T,DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR> const& converter, std::size_t iT,
                 SuperGeometry<T,2>& superGeometry, util::Timer<T>& timer, bool hasConverged,
                 Gnuplot<T>& gplot,
                 bool eoc)
{
  OstreamManager clout( std::cout,"getResults" );

  //VTK Data only if no EOC analysis
  if(eoc == false){

    SuperVTMwriter2D<T> vtmWriter( "poiseuille2d" );
    const int vtmIter  = converter.getLatticeTime( maxPhysT/20. );

    SuperLatticePhysVelocity2D<T,DESCRIPTOR> velocity( sLattice, converter );
    SuperLatticeGeometry2D<T,DESCRIPTOR> materials( sLattice, superGeometry );
    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor(materials);

    SuperLatticePhysPressure2D<T,DESCRIPTOR> pressure( sLattice, converter );
    if (flowType == nonForced) {
      vtmWriter.addFunctor( pressure );
    }

    const T maxVelocity = converter.getCharPhysVelocity();
    const T radius = (boundaryType == bounceBack) ?
      T(0.5) * (ly - converter.getPhysDeltaX()) : T(0.5) * ly;
    std::vector<T> axisPoint( 2,T() );
    axisPoint[0] = lx/2.;
    axisPoint[1] = ly/2.;
    std::vector<T> axisDirection( 2,T() );
    axisDirection[0] = 1;
    axisDirection[1] = 0;
    Poiseuille2D<T> analyticalVelocity( axisPoint, axisDirection, maxVelocity, radius );
    SuperLatticeFfromAnalyticalF2D<T,DESCRIPTOR> analyticalVelocityLattice(
      analyticalVelocity, sLattice);
    analyticalVelocityLattice.getName() = "analytical solution";
    if ((boundaryType != freeSlip) && (boundaryType != partialSlip)) {
      vtmWriter.addFunctor(analyticalVelocityLattice);
    }

    if ( iT==0 ) {
      // Writes the geometry, cuboid no. and rank no. as vti file for visualization
      SuperLatticeGeometry2D<T, DESCRIPTOR> geometry( sLattice, superGeometry );
      SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid( sLattice );
      SuperLatticeRank2D<T, DESCRIPTOR> rank( sLattice );
      SuperLatticeDiscreteNormal2D<T, DESCRIPTOR> discreteNormal(
        sLattice, superGeometry, superGeometry.getMaterialIndicator({2, 3}) );
      SuperLatticeDiscreteNormalType2D<T, DESCRIPTOR> discreteNormalType(
        sLattice, superGeometry, superGeometry.getMaterialIndicator({2, 3, 4, 5}) );
      vtmWriter.write( geometry );
      vtmWriter.write( cuboid );
      vtmWriter.write( rank );
      vtmWriter.write( discreteNormal );
      vtmWriter.write( discreteNormalType );

      vtmWriter.createMasterFile();
    }

    // Writes the vtm files and profile text file
    if ( iT%vtmIter==0 || hasConverged ) {
      sLattice.setProcessingContext(ProcessingContext::Evaluation);

      vtmWriter.write( iT );

      SuperEuklidNorm2D<T, DESCRIPTOR> normVel( velocity );
      BlockReduction2D2D<T> planeReduction( normVel, 600, BlockDataSyncMode::ReduceOnly );
      // write output as JPEG
      heatmap::write(planeReduction, iT);
    }
  }

  // Output on the console
  const int statIter = converter.getLatticeTime( maxPhysT/20. );
  if ( iT%statIter==0 || hasConverged || (iT == converter.getLatticeTime( maxPhysT ) - 1) ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );

    // Error norms
    if ((boundaryType != freeSlip) && (boundaryType != partialSlip)) {
      if ( (!eoc) || hasConverged || (iT == converter.getLatticeTime( maxPhysT ) - 1) ) {
        error( superGeometry, sLattice, converter );
      }
    }
  }

  // Gnuplot output
  if ((boundaryType != freeSlip) && (boundaryType != partialSlip)) {
    if ( hasConverged || (iT == converter.getLatticeTime( maxPhysT ) - 1) ) {
      sLattice.setProcessingContext(ProcessingContext::Evaluation);
      if ( eoc ) {
        if (flowType == nonForced){
          gplot.setData (
            T(converter.getResolution()),
            { velocityL1AbsError, velocityL2AbsError, velocityLinfAbsError,
              strainRateL1AbsError, strainRateL2AbsError, strainRateLinfAbsError,
              pressureL1AbsError, pressureL2AbsError, pressureLinfAbsError },
            { "velocity L1 abs Error","velocity L2 abs Error",
              "velocity Linf abs error","strain rate L1 abs error",
              "strain rate L2 abs error", "strain rate Linf abs error",
              "pressure L1 abs error", "pressure L2 abs error",
              "pressure Linf abs error" },
            "top right",
            { 'p','p','p','p','p','p','p','p','p' } );
        } else {
          // same as above, but without pressure computation
          gplot.setData (
            T(converter.getResolution()),
            { velocityL1AbsError, velocityL2AbsError, velocityLinfAbsError,
              strainRateL1AbsError, strainRateL2AbsError, strainRateLinfAbsError },
            { "velocity L1 abs Error","velocity L2 abs Error",
              "velocity Linf abs error","strain rate L1 abs error",
              "strain rate L2 abs error", "strain rate Linf abs error"},
            "top right",
            { 'p','p','p','p','p', 'p' } );
        }
      }
      else {  // if !eoc
        // plot velocity magnitude over line through the center of the simulation domain
        const T maxVelocity = converter.getPhysVelocity( converter.getCharLatticeVelocity() );
        T dx = 1. / T(converter.getResolution());
        T Ly = ly / converter.getConversionFactorLength();
        const T radius = (boundaryType == bounceBack) ?
          T(0.5) * (ly - converter.getPhysDeltaX()) : T(0.5) * ly;
        std::vector<T> axisPoint{ lx/2., ly/2. };
        std::vector<T> axisDirection{ 1, 0 };
        Poiseuille2D<T> uSol( axisPoint, axisDirection, maxVelocity, radius );
        SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity( sLattice, converter );
        AnalyticalFfromSuperF2D<T> intpolateVelocity( velocity, true );
        T point[2] { };
        point[0] = lx/2.;
        T analytical[2] { };
        T numerical[2] { };
        for ( int iY=0; iY<=Ly; ++iY ) {
          point[1] = ( T )iY/Ly;
          uSol( analytical,point );
          intpolateVelocity( numerical,point );
          gplot.setData( iY*dx, {analytical[0],numerical[0]}, {"analytical","numerical"} );
        }
        gplot.writePNG();
      }
    }
  }
}

// simulation method as former main method, runs the simulation for the parameter N
// decided whether eoc anlysis or not
void simulatePoiseuille( int N, Gnuplot<T>& gplot, bool eoc )
{
    OstreamManager clout( std::cout,"simulatePoiseuille" );

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
    int {N},     // resolution: number of voxels per charPhysL
    (T)   0.8,   // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   1,     // charPhysLength: reference length of simulation geometry
    (T)   1,     // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   1./Re, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.0    // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("poiseuille2d");


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
  CuboidGeometry2D<T> cuboidGeometry(
    cuboid, converter.getConversionFactorLength(), noOfCuboids );


  if (flowType == forced) {
    // Periodic boundaries in x-direction
    cuboidGeometry.setPeriodicity( true, false );
  }

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  const int overlap = (flowType == forced) ? 2 : 3;
  SuperGeometry<T,2> superGeometry( cuboidGeometry, loadBalancer, overlap );

  prepareGeometry( converter, superGeometry );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice( superGeometry );

  // Prepare lattice and set boundary conditions
  prepareLattice(converter, sLattice, superGeometry);

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( converter.getLatticeTime( maxPhysT ),
    superGeometry.getStatistics().getNvoxel() );
  util::ValueTracer<T> converge( converter.getLatticeTime( physInterval ), residuum );
  timer.start();

  for ( std::size_t iT = 0; iT < converter.getLatticeTime( maxPhysT ); ++iT ) {
    if ( converge.hasConverged() ) {
      clout << "Simulation converged." << std::endl;
      getResults( sLattice, converter, iT, superGeometry, timer, converge.hasConverged() ,gplot, eoc );

      break;
    }

    // === 5th Step: Definition of Initial and Boundary Conditions ===
    // in this application no boundary conditions have to be adjusted

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, iT,
      superGeometry, timer, converge.hasConverged() ,gplot, eoc );
    converge.takeValue( sLattice.getStatistics().getMaxU(), false );
  }

  timer.stop();
  timer.printSummary();
}
