/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2006, 2007, 2012 Jonas Latt, Mathias J. Krause,
 *  Louis Kronberg, Christian Vorwerk, Bastian Schäffauer
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

/* bstep2d.cpp:
 * The implementation of a backward facing step. It is furthermore
 * shown how to use checkpointing to save the state of the
 * simulation regularly.
 * The geometry of the step is based on the experiment described in
 * [Armaly, B.F., Durst, F., Pereira, J. C. F. and Schönung, B. Experimental
 * and theoretical investigation of backward-facing step flow. 1983.
 * J. Fluid Mech., vol. 127, pp. 473-496, DOI: 10.1017/S0022112083002839]
 */

#include "olb2D.h"
#include "olb2D.hh"


using namespace olb;
using namespace olb::descriptors;
using T = double;
using DESCRIPTOR = D2Q9<>;

using BulkDynamics = BGKdynamics<T,DESCRIPTOR>;

// Parameters for the simulation setup
const T lengthStep        = 0.2;                           // length of step in meter
const T heightStep        = 0.0049;                        // height of step in meter
const T lengthChannel     = 0.7;                           // length of channel in meter
const T heightChannel     = 0.0101;                        // height of channel in meter
const T heightInlet       = heightChannel - heightStep;    // height of inlet channel in meter
const T charL             = 2 * heightInlet;               // characteristic length
const int N               = 60;                            // resolution of the model
const T maxPhysT          = 2.;                            // max. simulation time in s, SI unit
const T relaxationTime    = 0.518;


// Stores geometry information in form of material numbers
void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter,
                      SuperGeometry<T,2>& superGeometry,
                      std::shared_ptr<IndicatorF2D<T>> channel,
                      std::shared_ptr<IndicatorF2D<T>> step )
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  // material numbers from zero to 2 inside geometry defined by indicator
  superGeometry.rename(0,2, channel-step);
  superGeometry.rename(2,1,{1,1} );

  Vector<T,2> extendBC_out( 0 + 1.*converter.getPhysDeltaX(),heightChannel );
  Vector<T,2> extendBC_in( 0, heightInlet );
  Vector<T,2> originBC_out( lengthChannel - 1.*converter.getPhysDeltaX(),0 );
  Vector<T,2> originBC_in( 0, heightStep);

  IndicatorCuboid2D<T> inflow( extendBC_in, originBC_in );
  // Set material number for inflow
  superGeometry.rename( 2,3,1,inflow );

  IndicatorCuboid2D<T> outflow( extendBC_out, originBC_out );
  // Set material number for outflow
  superGeometry.rename( 2,4,1,outflow );

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.getStatistics().print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice( UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperLattice<T,DESCRIPTOR>& sLattice,
                     SuperGeometry<T,2>& superGeometry )
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  auto bulkIndicator = superGeometry.getMaterialIndicator({1, 3, 4});

  // Material=0 -->do nothing
  sLattice.defineDynamics<NoDynamics>(superGeometry, 0);
  // Material=1 -->bulk dynamics
  // Material=3 -->bulk dynamics (inflow)
  // Material=4 -->bulk dynamics (outflow)
  sLattice.defineDynamics<BulkDynamics>(bulkIndicator);
  // Material=2 -->bounce back
  sLattice.defineDynamics<BounceBack>(superGeometry, 2);

  //if boundary conditions are chosen to be local
  setLocalVelocityBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 3);
  setLocalPressureBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 4);

  //if boundary conditions are chosen to be interpolated
  //setInterpolatedVelocityBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 3);
  //setInterpolatedPressureBoundary<T,DESCRIPTOR>(sLattice, omega, superGeometry, 4);

  // Initial conditions
  AnalyticalConst2D<T,T> ux( 0. );
  AnalyticalConst2D<T,T> uy( 0. );
  AnalyticalConst2D<T,T> rho( 1. );
  AnalyticalComposed2D<T,T> u( ux,uy );

  //Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU( bulkIndicator, rho, u );
  sLattice.iniEquilibrium( bulkIndicator, rho, u );

  sLattice.setParameter<descriptors::OMEGA>(omega);

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setBoundaryValues( UnitConverter<T,DESCRIPTOR> const& converter,
                        SuperLattice<T,DESCRIPTOR>& sLattice, int iT,
                        SuperGeometry<T,2>& superGeometry )
{
  OstreamManager clout( std::cout,"setBoundaryValues" );

  // time for smooth start-up
  int iTmaxStart = converter.getLatticeTime( maxPhysT*0.2 );
  int iTupdate = 100;

  if ( iT%iTupdate == 0 && iT<= iTmaxStart ) {
    // Smooth start curve, sinus
    // SinusStartScale<T,int> StartScale(iTmaxStart, T(1));
    // Smooth start curve, polynomial
    PolynomialStartScale<T,int> StartScale( iTmaxStart, T( 1 ) );
    // Creates and sets the Poiseuille inflow profile using functors
    int iTvec[1] = {iT};
    T frac[1] = {};
    StartScale( frac,iTvec );
    T maxVelocity = converter.getCharLatticeVelocity()*3./2.*frac[0];
    T distance2Wall = converter.getConversionFactorLength()/2.;
    Poiseuille2D<T> poiseuilleU( superGeometry, 3, maxVelocity, distance2Wall );
    // define lattice speed on inflow
    sLattice.defineU( superGeometry, 3, poiseuilleU );

    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
}

// write data to termimal and file system
void getResults( SuperLattice<T,DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                 SuperGeometry<T,2>& superGeometry, util::Timer<T>& timer,
                 SuperPlaneIntegralFluxVelocity2D<T>& velocityFlux,
                 SuperPlaneIntegralFluxPressure2D<T>& pressureFlux )
{
  OstreamManager clout( std::cout,"getResults" );
  SuperVTMwriter2D<T> vtmWriter( "bstep2d" );

  if ( iT==0 ) {
    // Writes geometry, cuboid no. and rank no. to file system
    SuperLatticeGeometry2D<T,DESCRIPTOR> geometry( sLattice, superGeometry );
    SuperLatticeCuboid2D<T,DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank2D<T,DESCRIPTOR> rank( sLattice );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }

  // Writes every 0.2  seconds
  // Writes every 0.1 simulated
  if ( iT%converter.getLatticeTime( 0.1 )==0 ) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    velocityFlux.print();
    pressureFlux.print();

    // write to terminal
    timer.update( iT );
    timer.printStep();
    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );
  }

  if ( iT%converter.getLatticeTime( 0.2 )==0 ) {
    SuperLatticePhysVelocity2D<T,DESCRIPTOR> velocity( sLattice, converter );
    SuperLatticePhysPressure2D<T,DESCRIPTOR> pressure( sLattice, converter );
    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( pressure );
    // write vtk to file system
    vtmWriter.write( iT );
    SuperEuklidNorm2D<T,DESCRIPTOR> normVel( velocity );
    BlockReduction2D2D<T> planeReduction( normVel, 1200, BlockDataSyncMode::ReduceOnly );
    // write output as JPEG
    heatmap::plotParam<T> jpeg_Param;
    jpeg_Param.maxValue = converter.getCharPhysVelocity() * 3./2.;
    jpeg_Param.minValue = 0.0;
    jpeg_Param.fullScreenPlot = true;
    heatmap::write(planeReduction, iT, jpeg_Param);
  }

  // Saves lattice data
  if ( iT%converter.getLatticeTime( 1 )==0 && iT>0 ) {
    //clout << "Checkpointing the system at t=" << iT << std::endl;
    //sLattice.save( "bstep2d.checkpoint" );
    // The data can be reloaded using
    //     sLattice.load("bstep2d.checkpoint");
  }
}


int main( int argc, char* argv[] )
{
  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );  // set output directory
  OstreamManager clout( std::cout, "main" );

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> converter(
    (T)   N,                 // resolution
    (T)   relaxationTime,    // relaxation time
    (T)   charL,             // charPhysLength: reference length of simulation geometry
    (T)   1.,                // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   1./19230.76923,    // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.                 // physDensity: physical density in __kg / m^3__
  );

  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("bstep2d");

  // === 2nd Step: Prepare Geometry ===

  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 6;
#endif

  // setup channel
  Vector<T,2> extendChannel( lengthChannel, heightChannel );
  Vector<T,2> originChannel( 0, 0 );
  std::shared_ptr<IndicatorF2D<T>> channel = std::make_shared<IndicatorCuboid2D<T>>( extendChannel, originChannel );

  // setup step
  Vector<T,2> extendStep( lengthStep, heightStep );
  Vector<T,2> originStep( 0, 0);
  std::shared_ptr<IndicatorF2D<T>> step = std::make_shared<IndicatorCuboid2D<T>>( extendStep, originStep );

  CuboidGeometry2D<T> cuboidGeometry( *(channel-step), converter.getConversionFactorLength(), noOfCuboids );

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  SuperGeometry<T,2> superGeometry( cuboidGeometry, loadBalancer );

  prepareGeometry( converter, superGeometry, channel, step );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T,DESCRIPTOR> sLattice( superGeometry );

  //prepare Lattice and set boundaryConditions
  prepareLattice( converter, sLattice, superGeometry );

  // instantiate reusable functors
  SuperPlaneIntegralFluxVelocity2D<T> velocityFlux( sLattice,
      converter,
      superGeometry,
  {lengthStep/2.,  heightInlet / 2.},
  {0.,  1.} );

  SuperPlaneIntegralFluxPressure2D<T> pressureFlux( sLattice,
      converter,
      superGeometry,
  {lengthStep/2.,  heightInlet / 2. },
  {0.,  1.} );

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for ( std::size_t iT = 0; iT < converter.getLatticeTime( maxPhysT ); ++iT ) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( converter, sLattice, iT, superGeometry );
    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();
    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, iT, superGeometry, timer, velocityFlux, pressureFlux );
  }

  timer.stop();
  timer.printSummary();
}
