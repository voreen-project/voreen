/*  Lattice Boltzmann sample, written in C++, using the OpenLBlibrary
 *
 *  Copyright (C) 2006-2019 Jonas Latt, Mathias J. Krause,
 *  Vojtech Cvrcek, Peter Weisbrod, Adrian Kummerländer
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

/* cylinder2d.cpp:
 * This example examines a steady flow past a cylinder placed in a channel.
 * The cylinder is offset somewhat from the center of the flow to make the
 * steady-state symmetrical flow unstable. At the inlet, a Poiseuille profile is
 * imposed on the velocity, whereas the outlet implements a Dirichlet pressure
 * condition set by p = 0.
 * Inspired by "Benchmark Computations of Laminar Flow Around
 * a Cylinder" by M.Schäfer and S.Turek. For high resolution, low
 * latticeU, and enough time to converge, the results for pressure drop, drag
 * and lift lie within the estimated intervals for the exact results.
 * An unsteady flow with Karman vortex street can be created by changing the
 * Reynolds number to Re=100.
 */

#include "olb2D.h"
#include "olb2D.hh"


using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = double;
using DESCRIPTOR = D2Q9<>;

#ifndef PLATFORM_GPU_CUDA
#define BOUZIDI
#endif

// Parameters for the simulation setup
const int N = 10;       // resolution of the model
const T Re = 20.;       // Reynolds number
const T maxPhysT = 16;  // max. simulation time in s, SI unit
const T L = 0.1/N;      // latticeL
const T lengthX = 2.2;
const T lengthY = .41+L;
const T centerCylinderX = 0.2;
const T centerCylinderY = 0.2+L/2.;
const T radiusCylinder = 0.05;


// Stores geometry information in form of material numbers
void prepareGeometry( UnitConverter<T, DESCRIPTOR> const& converter,
                      SuperGeometry<T,2>& superGeometry,
                      std::shared_ptr<IndicatorF2D<T>> circle)
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  Vector<T,2> extend( lengthX,lengthY );
  Vector<T,2> origin;

  superGeometry.rename( 0,2 );

  superGeometry.rename( 2,1,{1,1} );

  // Set material number for inflow
  extend[0] = 2.*L;
  origin[0] = -L;
  IndicatorCuboid2D<T> inflow( extend, origin );
  superGeometry.rename( 2,3,1,inflow );
  // Set material number for outflow
  origin[0] = lengthX-L;
  IndicatorCuboid2D<T> outflow( extend, origin );
  superGeometry.rename( 2,4,1,outflow );
  // Set material number for cylinder
  superGeometry.rename( 1,5, circle );

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice( SuperLattice<T,DESCRIPTOR>& sLattice,
                     UnitConverter<T, DESCRIPTOR> const& converter,
                     SuperGeometry<T,2>& superGeometry,
                     std::shared_ptr<IndicatorF2D<T>> circle)
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=0 -->do nothing
  sLattice.defineDynamics<NoDynamics>(superGeometry, 0);

  // Material=1 -->bulk dynamics
  auto bulkIndicator = superGeometry.getMaterialIndicator({1});
  sLattice.defineDynamics<BGKdynamics>(bulkIndicator);

  // Material=2 -->bounce back
  sLattice.defineDynamics<BounceBack>(superGeometry, 2);

  // Setting of the boundary conditions

  //if boundary conditions are chosen to be local
  //setLocalVelocityBoundary(sLattice, omega, superGeometry, 3);
  //setLocalPressureBoundary(sLattice, omega, superGeometry, 4);

  //if boundary conditions are chosen to be interpolated
  setInterpolatedVelocityBoundary(sLattice, omega, superGeometry, 3);
  setInterpolatedPressureBoundary(sLattice, omega, superGeometry, 4);

  // Material=5 -->bouzidi / bounce back

  #ifdef BOUZIDI
  sLattice.defineDynamics<NoDynamics>(superGeometry, 5);
  setBouzidiZeroVelocityBoundary(sLattice, superGeometry, 5, *circle);
  #else
  sLattice.defineDynamics<BounceBack>(superGeometry, 5);
  #endif

  // Initial conditions
  AnalyticalConst2D<T,T> rhoF( 1 );
  std::vector<T> velocity( 2,T( 0 ) );
  AnalyticalConst2D<T,T> uF( velocity );

  // Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU( bulkIndicator, rhoF, uF );
  sLattice.iniEquilibrium( bulkIndicator, rhoF, uF );

  sLattice.setParameter<descriptors::OMEGA>(omega);

  // Make the lattice ready for simulation
  sLattice.initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}

// Generates a slowly increasing inflow for the first iTMaxStart timesteps
void setBoundaryValues( SuperLattice<T, DESCRIPTOR>& sLattice,
                        UnitConverter<T, DESCRIPTOR> const& converter, int iT,
                        SuperGeometry<T,2>& superGeometry )
{

  OstreamManager clout( std::cout,"setBoundaryValues" );

  // No of time steps for smooth start-up
  int iTmaxStart = converter.getLatticeTime( maxPhysT*0.4 );
  int iTupdate = 5;

  if ( iT%iTupdate==0 && iT<= iTmaxStart ) {
    // Smooth start curve, sinus
    // SinusStartScale<T,int> StartScale(iTmaxStart, T(1));

    // Smooth start curve, polynomial
    PolynomialStartScale<T,T> StartScale( iTmaxStart, T( 1 ) );

    // Creates and sets the Poiseuille inflow profile using functors
    T iTvec[1] = {T( iT )};
    T frac[1] = {};
    StartScale( frac,iTvec );
    T maxVelocity = converter.getCharLatticeVelocity()*3./2.*frac[0];
    T distance2Wall = L/2.;
    Poiseuille2D<T> poiseuilleU( superGeometry, 3, maxVelocity, distance2Wall );

    sLattice.defineU( superGeometry, 3, poiseuilleU );

    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
}

// Computes the pressure drop between the voxels before and after the cylinder
void getResults( SuperLattice<T, DESCRIPTOR>& sLattice,
                 UnitConverter<T, DESCRIPTOR> const& converter, std::size_t iT,
                 SuperGeometry<T,2>& superGeometry, util::Timer<T>& timer )
{
  OstreamManager clout( std::cout,"getResults" );

  // Gnuplot constructor (must be static!)
  // for real-time plotting: gplot("name", true) // experimental!
  static Gnuplot<T> gplot( "drag" );

  SuperVTMwriter2D<T> vtmWriter( "cylinder2d" );
  SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure2D<T, DESCRIPTOR> pressure( sLattice, converter );
  SuperLatticeGeometry2D<T, DESCRIPTOR> materials( sLattice, superGeometry );
  SuperLatticeRefinementMetricKnudsen2D<T, DESCRIPTOR> quality( sLattice, converter);
  SuperRoundingF2D<T> roundedQuality ( quality, RoundingMode::NearestInteger);
  SuperDiscretizationF2D<T> discretization ( roundedQuality, 0., 2. );

  vtmWriter.addFunctor( materials );
  vtmWriter.addFunctor( quality );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );

  const int vtkIter  = converter.getLatticeTime( .3 );
  const int statIter = converter.getLatticeTime( .1 );

  T point[2] = {};
  point[0] = centerCylinderX + 3*radiusCylinder;
  point[1] = centerCylinderY;
  AnalyticalFfromSuperF2D<T> intpolateP( pressure, true );
  T p;
  intpolateP( &p,point );

  if ( iT == 0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry2D<T, DESCRIPTOR> geometry( sLattice, superGeometry );
    SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid( sLattice );
    SuperLatticeRank2D<T, DESCRIPTOR> rank( sLattice );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );

    vtmWriter.createMasterFile();
  }

  if ( iT%statIter == 0 ) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );

    // Drag, lift, pressure drop
    AnalyticalFfromSuperF2D<T> intpolatePressure( pressure, true );
    SuperLatticePhysDrag2D<T,DESCRIPTOR> drag( sLattice, superGeometry, 5, converter );

    T point1[2] = {};
    T point2[2] = {};

    point1[0] = centerCylinderX - radiusCylinder;
    point1[1] = centerCylinderY;

    point2[0] = centerCylinderX + radiusCylinder;
    point2[1] = centerCylinderY;

    T p1, p2;
    intpolatePressure( &p1,point1 );
    intpolatePressure( &p2,point2 );

    clout << "pressure1=" << p1;
    clout << "; pressure2=" << p2;

    T pressureDrop = p1-p2;
    clout << "; pressureDrop=" << pressureDrop;

    int input[3] = {};
    T _drag[drag.getTargetDim()];
    drag( _drag,input );
    clout << "; drag=" << _drag[0] << "; lift=" << _drag[1] << std::endl;

    // set data for gnuplot: input={xValue, yValue(s), names (optional), position of key (optional)}
    gplot.setData( converter.getPhysTime( iT ), {_drag[0], 5.58}, {"drag(openLB)", "drag(schaeferTurek)"}, "bottom right", {'l','l'} );

    // every (iT%vtkIter) write an png of the plot
    if ( iT%( vtkIter ) == 0 ) {
      // writes pngs: input={name of the files (optional), x range for the plot (optional)}
      gplot.writePNG( iT, maxPhysT );
    }
  }

  // Writes the vtk files
  if ( iT%vtkIter == 0 && iT > 0 ) {
    vtmWriter.write( iT );

    {
      SuperEuklidNorm2D<T, DESCRIPTOR> normVel( velocity );
      BlockReduction2D2D<T> planeReduction( normVel, 600, BlockDataSyncMode::ReduceOnly );
      // write output as JPEG
      heatmap::write(planeReduction, iT);
    }
    {
      BlockReduction2D2D<T> planeReduction( discretization, 600, BlockDataSyncMode::ReduceOnly );
      heatmap::plotParam<T> jpeg_scale;
      jpeg_scale.name = "quality";
      jpeg_scale.colour = "blackbody";
      heatmap::write( planeReduction, iT, jpeg_scale );
    }
  }

  // write pdf at last time step
  if ( iT == converter.getLatticeTime( maxPhysT )-1 ) {
    // writes pdf
    gplot.writePDF();
  }
}

int main( int argc, char* argv[] )
{
  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
    int {N},                        // resolution: number of voxels per charPhysL
    (T)   0.56,                     // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   2.0*radiusCylinder,       // charPhysLength: reference length of simulation geometry
    (T)   0.2,                      // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   0.2*2.*radiusCylinder/Re, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.0                       // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("cylinder2d");

  // === 2rd Step: Prepare Geometry ===
  Vector<T,2> extend( lengthX,lengthY );
  Vector<T,2> origin;
  IndicatorCuboid2D<T> cuboid( extend, origin );

  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 7;
#endif
  CuboidGeometry2D<T> cuboidGeometry( cuboid, L, noOfCuboids );

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  SuperGeometry<T,2> superGeometry( cuboidGeometry, loadBalancer );

  Vector<T,2> center( centerCylinderX,centerCylinderY );
  std::shared_ptr<IndicatorF2D<T>> circle = std::make_shared<IndicatorCircle2D<T>>( center, radiusCylinder );

  prepareGeometry( converter, superGeometry, circle );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice( superGeometry );

  //prepareLattice and set boundaryConditions
  prepareLattice( sLattice, converter, superGeometry, circle );

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for ( std::size_t iT = 0; iT < converter.getLatticeTime( maxPhysT ); ++iT ) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( sLattice, converter, iT, superGeometry );

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, iT, superGeometry, timer );
  }

  timer.stop();
  timer.printSummary();
}
