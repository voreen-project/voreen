/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2011-2016 Robin Trunk, Mathias J. Krause
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

/* bifurcation3d.cpp:
 * This example examines a steady particulate flow past a bifurcation. At the inlet,
 * an inflow condition with grad_n u = 0 and rho = 1 is implemented.
 * At both outlets, a Poiseuille profile is imposed on the velocity.
 * After a start time, particles are put into the bifurcation by imposing
 * a inflow condition with rho = 1 on the second euler phase at the inlet.
 * The particles are simulated as a continuum with a advection-diffusion equation
 * and experience a stokes drag force.
 *
 * A publication using the same geometry can be found here:
 * http://link.springer.com/chapter/10.1007/978-3-642-36961-2_5
 *  *
 */

#include "olb3D.h"
#include "olb3D.hh"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::util;

typedef double T;
typedef D3Q19<> NSDESCRIPTOR;
typedef D3Q7<VELOCITY,VELOCITY2> ADDESCRIPTOR;

const T Re = 50;               // Reynolds number
const int N = 19;              // resolution of the model
const int iTperiod = 100;      // amount of timesteps when new boundary conditions are reset and results are visualized
const T diffusion = 1.e-6;     // diffusion coefficient for advection-diffusion equation
const T radius = 1.5e-04;      // particles radius
const T partRho = 998.2;       // particles density
const T maxPhysT = 10.;        // max. simulation time in s, SI unit

// center of inflow and outflow regions [m]
Vector<T,3> inletCenter( T(), T(), 0.0786395 );
Vector<T,3> outletCenter0( -0.0235929682287551, -0.000052820468762797, -0.021445708949909 );
Vector<T,3> outletCenter1( 0.0233643529416147, 0.00000212439067050152, -0.0211994104877918 );

// radii of inflow and outflow regions [m]
T inletRadius = 0.00999839;
T outletRadius0 = 0.007927;
T outletRadius1 = 0.00787134;

// normals of inflow and outflow regions
Vector<T,3> inletNormal( T(), T(), T( -1 ) );
Vector<T,3> outletNormal0( 0.505126, -0.04177, 0.862034 );
Vector<T,3> outletNormal1( -0.483331, -0.0102764, 0.875377 );

void prepareGeometry( UnitConverter<T,NSDESCRIPTOR> const& converter,
                      IndicatorF3D<T>& indicator, STLreader<T>& stlReader,
                      SuperGeometry<T,3>& superGeometry )
{

  OstreamManager clout( std::cout, "prepareGeometry" );

  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( 0, 2, indicator );
  superGeometry.rename( 2, 1, stlReader );

  superGeometry.clean();

  // rename the material at the inlet
  IndicatorCircle3D<T> inletCircle( inletCenter, inletNormal, inletRadius );
  IndicatorCylinder3D<T> inlet( inletCircle, 2 * converter.getConversionFactorLength() );
  superGeometry.rename( 2, 3, 1, inlet );

  // rename the material at the outlet0
  IndicatorCircle3D<T> outletCircle0( outletCenter0, outletNormal0, 0.95*outletRadius0 );
  IndicatorCylinder3D<T> outlet0( outletCircle0, 4 * converter.getConversionFactorLength() );
  superGeometry.rename( 2, 4, outlet0 );

  // rename the material at the outlet1
  IndicatorCircle3D<T> outletCircle1( outletCenter1, outletNormal1, 0.95*outletRadius1 );
  IndicatorCylinder3D<T> outlet1( outletCircle1, 4 * converter.getConversionFactorLength() );
  superGeometry.rename( 2, 5, outlet1 );

  IndicatorCircle3D<T> inletCircleExtended( inletCenter, inletNormal, inletRadius + 2 * converter.getConversionFactorLength() );
  IndicatorCylinder3D<T> inletExtended( inletCircleExtended, 2 * converter.getConversionFactorLength() );
  superGeometry.rename(2, 6, inletExtended);

  superGeometry.clean();
  superGeometry.innerClean( true );
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
  return;
}

void prepareLattice( SuperLattice<T, NSDESCRIPTOR>& sLatticeNS,
                     SuperLattice<T, ADDESCRIPTOR>& sLatticeAD,
                     UnitConverter<T,NSDESCRIPTOR> const& converter,
                     SuperGeometry<T,3>& superGeometry,
                     T omegaAD )
{

  OstreamManager clout( std::cout, "prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeRelaxationFrequency();

  // Material=0 --> do nothing
  sLatticeNS.defineDynamics<NoDynamics>(superGeometry, 0);
  sLatticeAD.defineDynamics<NoDynamics>(superGeometry, 0);

  // Material=1 --> bulk dynamics
  // Material=3 --> bulk dynamics (inflow)
  auto inflowIndicator = superGeometry.getMaterialIndicator({1, 3});
  sLatticeNS.defineDynamics<BGKdynamics>(inflowIndicator);
  sLatticeAD.defineDynamics<ParticleAdvectionDiffusionBGKdynamics>(inflowIndicator);

  // Material=2 --> bounce-back / zero distribution dynamics
  sLatticeNS.defineDynamics<BounceBack>(superGeometry, 2);
  sLatticeAD.defineDynamics<ZeroDistributionDynamics>(superGeometry, 2);

  // Material=4,5 -->bulk dynamics / do-nothing (outflow)
  auto outflowIndicator = superGeometry.getMaterialIndicator({4, 5});
  sLatticeNS.defineDynamics<BGKdynamics>(outflowIndicator);
  sLatticeAD.defineDynamics<NoDynamics>(outflowIndicator);

  // Material=6 --> bounce-back / bounce-back
  sLatticeNS.defineDynamics<BounceBack>(superGeometry, 6);
  sLatticeAD.defineDynamics<BounceBack>(superGeometry, 6);

  // Setting of the boundary conditions
  setInterpolatedPressureBoundary<T,NSDESCRIPTOR>(sLatticeNS, omega, superGeometry, 3);
  setInterpolatedVelocityBoundary<T,NSDESCRIPTOR>(sLatticeNS, omega, outflowIndicator);
  setZeroDistributionBoundary<T,ADDESCRIPTOR>(sLatticeAD, superGeometry, 2);
  setAdvectionDiffusionTemperatureBoundary<T,ADDESCRIPTOR>(sLatticeAD, omegaAD, superGeometry, 3);
  setAdvectionDiffusionConvectionBoundary<T,ADDESCRIPTOR>(sLatticeAD, outflowIndicator);
  setExtFieldBoundary<T,ADDESCRIPTOR,descriptors::VELOCITY,descriptors::VELOCITY2>(sLatticeAD, superGeometry.getMaterialIndicator({2, 3, 4, 5, 6}));

  // Initial conditions
  AnalyticalConst3D<T,T> rho1( 1. );
  AnalyticalConst3D<T,T> rho0( 1.e-8 );
  std::vector<T> velocity( 3,T() );
  AnalyticalConst3D<T,T> u0( velocity );

  // Initialize all values of distribution functions to their local equilibrium
  sLatticeNS.defineRhoU( superGeometry.getMaterialIndicator({1, 2, 3, 4, 5, 6}), rho1, u0 );
  sLatticeNS.iniEquilibrium( superGeometry.getMaterialIndicator({1, 2, 3, 4, 5, 6}), rho1, u0 );
  sLatticeAD.defineRho( superGeometry, 3, rho1 );
  sLatticeAD.iniEquilibrium( superGeometry.getMaterialIndicator({1, 2, 4, 5, 6}), rho0, u0 );

  sLatticeNS.setParameter<descriptors::OMEGA>(omega);
  sLatticeAD.setParameter<descriptors::OMEGA>(omegaAD);

  // Lattice initialize
  sLatticeNS.initialize();
  sLatticeAD.initialize();

  {
    auto& communicator = sLatticeAD.getCommunicator(PostCoupling());
    communicator.requestField<descriptors::VELOCITY>();
    communicator.requestOverlap(sLatticeAD.getOverlap());
    communicator.exchangeRequests();
  }

  clout << "Prepare Lattice ... OK" << std::endl;
  return;
}

void setBoundaryValues( SuperLattice<T, NSDESCRIPTOR>& sLatticeNS,
                        UnitConverter<T,NSDESCRIPTOR> const& converter, int iT,
                        SuperGeometry<T,3>& superGeometry )
{

  OstreamManager clout( std::cout, "setBoundaryValues" );

  // No of time steps for smooth start-up
  int iTmaxStart = converter.getLatticeTime( 0.8*maxPhysT );
  // Set inflow velocity
  T maxVelocity = converter.getCharLatticeVelocity() * 3. / 4. * util::pow(
                    inletRadius, 2 ) / util::pow( outletRadius0, 2 );
  if ( iT % iTperiod == 0 ) {
    if ( iT <= iTmaxStart ) {
      SinusStartScale<T, int> startScale( iTmaxStart, T( 1 ) );
      int iTvec[1] = { iT };
      T frac[1] = { T( 0 ) };
      startScale( frac, iTvec );
      maxVelocity *= frac[0];
    }

    CirclePoiseuille3D<T> poiseuilleU4( outletCenter0[0], outletCenter0[1],
                                        outletCenter0[2], outletNormal0[0],
                                        outletNormal0[1], outletNormal0[2],
                                        outletRadius0 * 0.95, -maxVelocity );

    CirclePoiseuille3D<T> poiseuilleU5( outletCenter1[0], outletCenter1[1],
                                        outletCenter1[2], outletNormal1[0],
                                        outletNormal1[1], outletNormal1[2],
                                        outletRadius1 * 0.95, -maxVelocity );

    sLatticeNS.defineU( superGeometry, 4, poiseuilleU4 );
    sLatticeNS.defineU( superGeometry, 5, poiseuilleU5 );
  }
}

void getResults( SuperLattice<T, NSDESCRIPTOR>& sLatticeNS,
                 SuperLattice<T, ADDESCRIPTOR>& sLatticeAD,
                 UnitConverter<T,NSDESCRIPTOR> const& converter, int iT,
                 SuperGeometry<T,3>& superGeometry,
                 Timer<double>& timer )
{

  OstreamManager clout( std::cout, "getResults" );
  SuperVTMwriter3D<T> vtmWriter( "bifurcation3d_fluid" );
  SuperVTMwriter3D<T> vtmWriterAD( "bifurcation3d_particle" );
  SuperLatticePhysVelocity3D<T, NSDESCRIPTOR> velocity( sLatticeNS, converter );
  SuperLatticeVelocity3D<T, NSDESCRIPTOR> latticeVelocity( sLatticeNS);

  SuperLatticePhysPressure3D<T, NSDESCRIPTOR> pressure( sLatticeNS, converter );
  SuperLatticeDensity3D<T, ADDESCRIPTOR> particles( sLatticeAD );
  SuperLatticePhysField3D<T, ADDESCRIPTOR, descriptors::VELOCITY> extField(
    sLatticeAD, converter.getConversionFactorVelocity());

  vtmWriter.addFunctor( velocity );
  vtmWriter.addFunctor( pressure );
  vtmWriterAD.addFunctor( particles );
  vtmWriterAD.addFunctor( extField );

  if ( iT == 0 ) {
    SuperLatticeGeometry3D<T, NSDESCRIPTOR> geometry( sLatticeNS, superGeometry );
    SuperLatticeCuboid3D<T, NSDESCRIPTOR> cuboid( sLatticeNS );
    SuperLatticeRank3D<T, NSDESCRIPTOR> rank( sLatticeNS );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
    vtmWriterAD.createMasterFile();

    // Print some output of the chosen simulation setup
    clout << "N=" << N << "; maxTimeSteps=" << converter.getLatticeTime( maxPhysT )
          << "; noOfCuboid=" << superGeometry.getCuboidGeometry().getNc() << "; Re=" << Re
          << "; St=" << ( 2.*partRho*radius*radius*converter.getCharPhysVelocity() ) / ( 9.*converter.getPhysViscosity()*converter.getPhysDensity()*converter.getCharPhysLength() )
          << std::endl;
  }

  if ( iT % iTperiod == 0 ) {
    // Writes the vtk files
    vtmWriter.write( iT );
    vtmWriterAD.write( iT );

    // GIF Writer
    SuperEuklidNorm3D<T, NSDESCRIPTOR> normVel( velocity );
    HyperplaneLattice3D<T> gifLattice(
      superGeometry.getCuboidGeometry(),
      Hyperplane3D<T>()
      .centeredIn(superGeometry.getCuboidGeometry().getMotherCuboid())
      .normalTo({0, -1, 0}),
      600);
    BlockReduction3D2D<T> planeReductionVelocity( normVel, gifLattice, BlockDataSyncMode::ReduceOnly );
    BlockReduction3D2D<T> planeReductionParticles( particles, gifLattice, BlockDataSyncMode::ReduceOnly );
    // write output as JPEG
    heatmap::write(planeReductionVelocity, iT);
    heatmap::write(planeReductionParticles, iT);

    // Writes output on the console
    timer.update( iT );
    timer.printStep();
    sLatticeNS.getStatistics().print( iT, iT * converter.getCharLatticeVelocity() / T(converter.getResolution()) );

    // preparation for flux computations
    const std::vector<int> materials = { 1, 3, 4, 5 };
    IndicatorCircle3D<T> inlet( inletCenter + 2. * converter.getConversionFactorLength() * inletNormal,
                                inletNormal,
                                inletRadius + 2. * converter.getConversionFactorLength() );
    IndicatorCircle3D<T> outlet0( outletCenter0 + 2. * converter.getConversionFactorLength() * outletNormal0,
                                  outletNormal0,
                                  outletRadius0 + 2. * converter.getConversionFactorLength() );
    IndicatorCircle3D<T> outlet1( outletCenter1 + 2. * converter.getConversionFactorLength() * outletNormal1,
                                  outletNormal1,
                                  outletRadius1 + 2. * converter.getConversionFactorLength() );

    // Flux of the fluid at the inlet and outlet regions
    SuperPlaneIntegralFluxVelocity3D<T> vFluxInflow( sLatticeNS, converter, superGeometry, inlet, materials );
    vFluxInflow.print( "inflow", "ml/s" );
    SuperPlaneIntegralFluxPressure3D<T> pFluxInflow( sLatticeNS, converter, superGeometry, inlet, materials );
    pFluxInflow.print( "inflow", "N", "Pa" );
    SuperPlaneIntegralFluxVelocity3D<T> vFluxOutflow0( sLatticeNS, converter, superGeometry, outlet0, materials );
    vFluxOutflow0.print( "outflow0", "ml/s" );
    SuperPlaneIntegralFluxPressure3D<T> pFluxOutflow0( sLatticeNS, converter, superGeometry, outlet0, materials );
    pFluxOutflow0.print( "outflow0", "N", "Pa" );
    SuperPlaneIntegralFluxVelocity3D<T> vFluxOutflow1( sLatticeNS, converter, superGeometry, outlet1, materials );
    vFluxOutflow1.print( "outflow1", "ml/s" );
    SuperPlaneIntegralFluxPressure3D<T> pFluxOutflow1( sLatticeNS, converter, superGeometry, outlet1, materials );
    pFluxOutflow1.print( "outflow1", "N", "Pa" );

    int input[4] = {0};
    T mFlux[5] = {0.}, mFlux0[5] = {0.}, mFlux1[5] = {0.};
    // Flux of particles at the inlet and outlet regions: Inflow, Outflow0 and Outlfow1
    SuperPlaneIntegralFluxMass3D<T> mFluxInflow(latticeVelocity,particles,
        superGeometry, converter.getConversionFactorMass(),
        converter.getConversionFactorTime(), inlet, materials);
    SuperPlaneIntegralFluxMass3D<T> mFluxOutflow0(latticeVelocity, particles,
        superGeometry, converter.getConversionFactorMass(),
        converter.getConversionFactorTime(),outlet0, materials);
    SuperPlaneIntegralFluxMass3D<T> mFluxOutflow1(latticeVelocity, particles,
        superGeometry, converter.getConversionFactorMass(),
        converter.getConversionFactorTime(), outlet1, materials);

    mFluxInflow( mFlux,input );
    mFluxOutflow0( mFlux0,input );
    mFluxOutflow1( mFlux1,input );

    // Since more diffusion is added to ensure stability the computed escaperate falls short of the real value,
    // therefore it is scaled by the factor 1.4 computed by a simulation without drag force. This value is computed
    // for this specific setup. For further information see R.Trunk, T.Henn, W.Dörfler, H.Nirschl, M.J.Krause,
    // "Inertial Dilute Particulate Fluid Flow Simulations with an Euler-Euler Lattice Boltzmann Method"
    T escr = -( mFlux0[0]+mFlux1[0] )/mFlux[0]*1.4;
    clout << "escapeRate=" << escr << "; captureRate="<< 1-escr << std::endl;
  }
}

int main( int argc, char* argv[] )
{

  // === 1st Step: Initialization ===

  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout, "main" );

  UnitConverterFromResolutionAndRelaxationTime<T,NSDESCRIPTOR> const converter(
    int {N},                           // resolution: number of voxels per charPhysL
    (T)   0.557646,                    // latticeRelaxationTime: relaxation time, have to be greater than 0.5!
    (T)   inletRadius*2.,              // charPhysLength: reference length of simulation geometry
    (T)   Re*1.5e-5/( inletRadius*2 ), // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   1.5e-5,                      // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1.225                        // physDensity: physical density in __kg / m^3__
  );
  // Prints the converter log as console output
  converter.print();
  // Writes the converter log in a file
  converter.write("bifurcation3d");

  // compute relaxation parameter to solve the advection-diffusion equation in the lattice Boltzmann context
  T omegaAD = converter.getLatticeRelaxationFrequencyFromDiffusivity<ADDESCRIPTOR>( diffusion );

  // === 2nd Step: Prepare Geometry ===

  STLreader<T> stlReader( "../bifurcation3d.stl",converter.getConversionFactorLength() );
  IndicatorLayer3D<T> extendedDomain( stlReader,converter.getConversionFactorLength() );

  // Instantiation of an empty cuboidGeometry
  int noOfCuboids = util::max( 20, singleton::mpi().getSize() );
  CuboidGeometry3D<T> cuboidGeometry( extendedDomain, converter.getConversionFactorLength(),
                                      noOfCuboids, "weight" );
  cuboidGeometry.printExtended();
  clout << "min / max ratio (volume) = " << (T) cuboidGeometry.getMinLatticeVolume()
        / cuboidGeometry.getMaxLatticeVolume() << std::endl;
  clout << "min / max ratio (weight) = " << (T) cuboidGeometry.getMinLatticeWeight()
        / cuboidGeometry.getMaxLatticeWeight() << std::endl;

  // Instantiation of an empty loadBalancer
  HeuristicLoadBalancer<T> loadBalancer( cuboidGeometry );
  // Default instantiation of superGeometry
  SuperGeometry<T,3> superGeometry( cuboidGeometry, loadBalancer, 3 );

  prepareGeometry( converter, extendedDomain, stlReader, superGeometry );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, NSDESCRIPTOR> sLatticeNS( superGeometry );
  SuperLattice<T, ADDESCRIPTOR> sLatticeAD( superGeometry );


  //prepareLattice and setBoundaryConditions
  prepareLattice( sLatticeNS, sLatticeAD, converter, superGeometry, omegaAD );

  AdvectionDiffusionParticleCouplingGenerator3D<T,NSDESCRIPTOR,ADDESCRIPTOR,VELOCITY,VELOCITY2> coupling;

  AdvDiffDragForce3D<T, NSDESCRIPTOR, ADDESCRIPTOR> dragForce( converter,radius,partRho );
  coupling.addForce( dragForce );
  sLatticeNS.addLatticeCoupling( superGeometry, 1, coupling, sLatticeAD );

  // === 4th Step: Main Loop with Timer ===
  Timer<double> timer( converter.getLatticeTime( maxPhysT ),
                       superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for (std::size_t iT = 0; iT <= converter.getLatticeTime(maxPhysT); ++iT) {
    sLatticeAD.setParameter<descriptors::LATTICE_TIME>(iT);
    getResults( sLatticeNS, sLatticeAD, converter, iT, superGeometry, timer );
    setBoundaryValues( sLatticeNS, converter, iT, superGeometry );
    sLatticeNS.executeCoupling();
    sLatticeAD.executeCoupling();
    sLatticeNS.collideAndStream();
    sLatticeAD.collideAndStream();
  }
  timer.stop();
  timer.printSummary();

}
