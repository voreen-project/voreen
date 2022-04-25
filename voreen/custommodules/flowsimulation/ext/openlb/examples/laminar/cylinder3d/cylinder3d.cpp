/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2011-2013 Mathias J. Krause, Thomas Henn, Tim Dornieden
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

/* cylinder3d.cpp:
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
 * It also shows the usage of the STL-reader and explains how
 * to set boundary conditions automatically.
 */


#include "olb3D.h"
#include "olb3D.hh"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = double;
using DESCRIPTOR = D3Q19<>;

#ifndef PLATFORM_GPU_CUDA
#define BOUZIDI
#endif

// Parameters for the simulation setup
const int N = 10;        // resolution of the model
const T Re = 20.;       // Reynolds number
const T maxPhysT = 16.; // max. simulation time in s, SI unit

// Stores data from stl file in geometry in form of material numbers
void prepareGeometry( UnitConverter<T,DESCRIPTOR> const& converter, IndicatorF3D<T>& indicator,
                      STLreader<T>& stlReader, SuperGeometry<T,3>& superGeometry )
{

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( 0,2,indicator );
  superGeometry.rename( 2,1,stlReader );
  superGeometry.clean();

  Vector<T,3> origin = superGeometry.getStatistics().getMinPhysR( 2 );
  origin[1] += converter.getConversionFactorLength()/2.;
  origin[2] += converter.getConversionFactorLength()/2.;

  Vector<T,3> extend = superGeometry.getStatistics().getMaxPhysR( 2 );
  extend[1] = extend[1]-origin[1]-converter.getConversionFactorLength()/2.;
  extend[2] = extend[2]-origin[2]-converter.getConversionFactorLength()/2.;

  // Set material number for inflow
  origin[0] = superGeometry.getStatistics().getMinPhysR( 2 )[0]-converter.getConversionFactorLength();
  extend[0] = 2*converter.getConversionFactorLength();
  IndicatorCuboid3D<T> inflow( extend,origin );
  superGeometry.rename( 2,3,inflow );

  // Set material number for outflow
  origin[0] = superGeometry.getStatistics().getMaxPhysR( 2 )[0]-converter.getConversionFactorLength();
  extend[0] = 2*converter.getConversionFactorLength();
  IndicatorCuboid3D<T> outflow( extend,origin );
  superGeometry.rename( 2,4,outflow );

  // Set material number for cylinder
  origin[0] = superGeometry.getStatistics().getMinPhysR( 2 )[0]+converter.getConversionFactorLength();
  extend[0] = ( superGeometry.getStatistics().getMaxPhysR( 2 )[0]-superGeometry.getStatistics().getMinPhysR( 2 )[0] )/2.;
  std::shared_ptr<IndicatorF3D<T>> cylinder = std::make_shared<IndicatorCuboid3D<T>>( extend, origin );
  superGeometry.rename( 2,5, cylinder );

  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

// Set up the geometry of the simulation
void prepareLattice( SuperLattice<T,DESCRIPTOR>& sLattice,
                     UnitConverter<T,DESCRIPTOR> const& converter,
                     STLreader<T>& stlReader,
                     SuperGeometry<T,3>& superGeometry )
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

  //if local boundary conditions are chosen
  //setLocalVelocityBoundary(sLattice, omega, superGeometry, 3);
  //setLocalPressureBoundary(sLattice, omega, superGeometry, 4);

  //if interpolated boundary conditions are chosen
  setInterpolatedVelocityBoundary(sLattice, omega, superGeometry, 3);
  setInterpolatedPressureBoundary(sLattice, omega, superGeometry, 4);

  // Material=5 -->bouzidi / bounce back
  #ifdef BOUZIDI
  sLattice.defineDynamics<NoDynamics>(superGeometry, 5);
  setBouzidiZeroVelocityBoundary<T,DESCRIPTOR>(sLattice, superGeometry, 5, stlReader);
  #else
  sLattice.defineDynamics<BounceBack>(superGeometry, 5);
  #endif

  // Initial conditions
  AnalyticalConst3D<T,T> rhoF( 1 );
  Vector<T,3> velocityV;
  AnalyticalConst3D<T,T> uF(velocityV);

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
                        UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                        SuperGeometry<T,3>& superGeometry )
{
  OstreamManager clout( std::cout,"setBoundaryValues" );

  // No of time steps for smooth start-up
  int iTmaxStart = converter.getLatticeTime( maxPhysT*0.4 );
  int iTupdate = 30;

  if ( iT%iTupdate == 0 && iT <= iTmaxStart ) {
    // Smooth start curve, sinus
    // SinusStartScale<T,int> StartScale(iTmaxStart, T(1));

    // Smooth start curve, polynomial
    PolynomialStartScale<T,int> StartScale( iTmaxStart, T( 1 ) );

    // Creates and sets the Poiseuille inflow profile using functors
    int iTvec[1] = {iT};
    T frac[1] = {};
    StartScale( frac,iTvec );
    std::vector<T> maxVelocity( 3,0 );
    maxVelocity[0] = 2.25*frac[0]*converter.getCharLatticeVelocity();

    T distance2Wall = converter.getConversionFactorLength()/2.;
    RectanglePoiseuille3D<T> poiseuilleU( superGeometry, 3, maxVelocity, distance2Wall, distance2Wall, distance2Wall );
    sLattice.defineU( superGeometry, 3, poiseuilleU );

    clout << "step=" << iT << "; maxVel=" << maxVelocity[0] << std::endl;

    sLattice.setProcessingContext<Array<momenta::FixedVelocityMomentumGeneric::VELOCITY>>(
      ProcessingContext::Simulation);
  }
}

// Computes the pressure drop between the voxels before and after the cylinder
void getResults( SuperLattice<T, DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                 SuperGeometry<T,3>& superGeometry, util::Timer<T>& timer,
                 STLreader<T>& stlReader )
{

  OstreamManager clout( std::cout,"getResults" );

  SuperVTMwriter3D<T> vtmWriter( "cylinder3d" );
  SuperLatticePhysVelocity3D<T, DESCRIPTOR> velocity( sLattice, converter );
  SuperLatticePhysPressure3D<T, DESCRIPTOR> pressure( sLattice, converter );
  SuperLatticeYplus3D<T, DESCRIPTOR> yPlus( sLattice, converter, superGeometry, stlReader, 5 );
  SuperLatticeRefinementMetricKnudsen3D<T, DESCRIPTOR> quality( sLattice, converter );
  SuperRoundingF3D<T, T> roundedQuality ( quality, RoundingMode::NearestInteger );
  SuperDiscretizationF3D<T> discretization ( roundedQuality, 0., 2. );

  vtmWriter.addFunctor( quality );
  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );
  vtmWriter.addFunctor( yPlus );

  const int vtkIter  = converter.getLatticeTime( .3 );
  const int statIter = converter.getLatticeTime( .1 );

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

  // Writes output on the console
  if ( iT%statIter == 0 ) {
    sLattice.setProcessingContext(ProcessingContext::Evaluation);

    // Timer console output
    timer.update( iT );
    timer.printStep();

    // Lattice statistics console output
    sLattice.getStatistics().print( iT,converter.getPhysTime( iT ) );

    // Drag, lift, pressure drop
    AnalyticalFfromSuperF3D<T> intpolatePressure( pressure, true );
    SuperLatticePhysDrag3D<T,DESCRIPTOR> drag( sLattice, superGeometry, 5, converter );

    olb::Vector<T, 3> point1V = superGeometry.getStatistics().getCenterPhysR( 5 );
    olb::Vector<T, 3> point2V = superGeometry.getStatistics().getCenterPhysR( 5 );
    T point1[3] = {};
    T point2[3] = {};
    for ( int i = 0; i<3; i++ ) {
      point1[i] = point1V[i];
      point2[i] = point2V[i];
    }
    point1[0] = superGeometry.getStatistics().getMinPhysR( 5 )[0] - converter.getConversionFactorLength();
    point2[0] = superGeometry.getStatistics().getMaxPhysR( 5 )[0] + converter.getConversionFactorLength();

    T p1, p2;
    intpolatePressure( &p1,point1 );
    intpolatePressure( &p2,point2 );

    clout << "pressure1=" << p1;
    clout << "; pressure2=" << p2;

    T pressureDrop = p1-p2;
    clout << "; pressureDrop=" << pressureDrop;

    T dragA[3];
    int input1[0];
    drag( dragA, input1 );
    clout << "; drag=" << dragA[0] << "; lift=" << dragA[1] << std::endl;

    int input[4] = {};
    SuperMax3D<T> yPlusMaxF( yPlus, superGeometry, 1 );
    T yPlusMax[1];
    yPlusMaxF( yPlusMax,input );
    clout << "yPlusMax=" << yPlusMax[0] << std::endl;
  }

  // Writes the vtk files
  if ( iT%vtkIter == 0 ) {
    vtmWriter.write( iT );

    {
      SuperEuklidNorm3D<T, DESCRIPTOR> normVel( velocity );
      BlockReduction3D2D<T> planeReduction( normVel, Vector<T,3>({0, 0, 1}) );
      // write output as JPEG
      heatmap::write(planeReduction, iT);
    }

    {
      BlockReduction3D2D<T> planeReduction( discretization, Vector<T,3>({0, 0, 1}) );
      heatmap::plotParam<T> jpeg_scale;
      jpeg_scale.colour = "blackbody";
      jpeg_scale.name = "quality";
      heatmap::write( planeReduction, iT, jpeg_scale );
    }
  }
}

int main( int argc, char* argv[] )
{

  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );
  // display messages from every single mpi process
  //clout.setMultiOutput(true);

  UnitConverterFromResolutionAndRelaxationTime<T, DESCRIPTOR> const converter(
    int {N},              // resolution: number of voxels per charPhysL
    (T)   0.53,           // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   0.1,            // charPhysLength: reference length of simulation geometry
    (T)   0.2,            // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   0.2*2.*0.05/Re, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.0             // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("cylinder3d");


  // === 2nd Step: Prepare Geometry ===

  // Instantiation of the STLreader class
  // file name, voxel size in meter, stl unit in meter, outer voxel no., inner voxel no.
  STLreader<T> stlReader( "cylinder3d.stl", converter.getConversionFactorLength(), 0.001 );
  IndicatorLayer3D<T> extendedDomain( stlReader, converter.getConversionFactorLength() );

  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 7;
#endif
  CuboidGeometry3D<T> cuboidGeometry( extendedDomain, converter.getConversionFactorLength(), noOfCuboids );

  // Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );

  // Instantiation of a superGeometry
  SuperGeometry<T,3> superGeometry( cuboidGeometry, loadBalancer );

  prepareGeometry( converter, extendedDomain, stlReader, superGeometry );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice( superGeometry );

  //prepareLattice and set boundaryCondition
  prepareLattice( sLattice, converter, stlReader, superGeometry );

  // === 4th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( converter.getLatticeTime( maxPhysT ), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for (std::size_t iT = 0; iT < converter.getLatticeTime( maxPhysT ); ++iT) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( sLattice, converter, iT, superGeometry );

    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, iT, superGeometry, timer, stlReader );
  }

  timer.stop();
  timer.printSummary();
}

