/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2022 Davide Dapelo, Stephan Simonis
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

/*
 * This example simulates the a concentration transport through a pipe with a
 * square section. Whereas the fluid flow is simulated via approximating
 * Navier-Stokes with LBM (LatticeNS), the solute transport is modeled via a
 * one-way coupled advection-diffusion equation (LatticeAD) which is
 * approximated with finite differences (FD). Several FD schemes are
 * implemented below and tested in Dapelo et al., Journal of Computational
 * Science (2021) 51:101363, DOI: https://doi.org/10.1016/j.jocs.2021.101363.
 */

#define QUIET // disable BC's verbose warning
//#define NOVTM // disable vtm output
//#define DEBUG_LINEAR // gdb linear
//#define DEBUG_PARALLEL // gdb parallel
//#define TESTING // printed output as in gdb macros, but without gdb

#include "olb3D.h"
#include "olb3D.hh"

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

typedef double T;
typedef D3Q19<AD_FIELD,EUL2LAGR> DESCRIPTOR_NS;
typedef D3Q7<VELOCITY> DESCRIPTOR_AD;
typedef fd::tag::UPWIND  SCHEME_ADV_UPWIND1;
typedef fd::tag::UPWIND_2_ORDER  SCHEME_ADV_UPWIND2;
typedef fd::tag::CENTRAL  SCHEME_ADV_CENTRAL;
typedef fd::tag::CENTRAL SCHEME_DIFF;


///////////////////////////////////////////////////////////////////////////////////////////////////

size_t iT = 0; // global timestep

// Flags for advection-diffusion simulations
int  advectionScheme = 1;     // 1 = central; 2 = upwind 1st order; 3 = upwind 1st order with anti-diffusion; 4 = upwind 2nd order
bool fromCentre      = true; // Seeding point at the centre?

// (dimensionless) Parameters for the simulation setup
T    antiDiffusionTuning = 1.; // Anti-diffusion coefficient (0=0%; 1=100%)
int  nx         = 500;  // Number of internal lattice points along the x direction
int  ny         = 7;   // Number of internal lattice points along the y direction
size_t  tMax       = 5001; // Total number of lattice updates
size_t  tVtm       = 500;  // Number of timesteps before producing output

// Dimensional parameters
T charU       = 0.05;   // Flow velocity (m/s)
T deltaX      = 0.1;  // Lattice spacing (m)
T deltaT      = 0.01; // Timestep (s)
T diffusivity = 0.005; // Diffusivity (m^2/s)

// Parameters computed from the input
int nz = 11; // Number of internal lattice points along the z direction

// Other parameters
bool noPlots = false; // Prevent from plots (for scaling testing)?


///////////////////////////////////////////////////////////////////////////////////////////////////
// Stores geometry information in form of material numbers
void prepareGeometry(SuperGeometry<T,3>& superGeometry)
{
  /* MAT NUM | GEOMETRY     | NAVIER-STOKES | FINITE-DIFFERENCE | lATTICE-BOLTZMANN
   * 1       | Near-Border  | Bulk          | Bulk              | Bulk
   * 2       | Wall         | Free-slip     | No-penetration    | Bounce-back
   * 3       | Inlet-Outlet | Bulk          | No-penetration    | Bounce-back
   * 4       | Seed         | Bulk          | Bulk              | Bulk
   * 5       | Bulk         | Bulk          | Bulk              | Bulk
   */
  OstreamManager clout(std::cout,"prepareGeometry");
  if (! noPlots) clout << "Prepare Geometry ..." << std::endl;

  std::vector<T> origin { -(0.5*nx + 1)*deltaX, -(0.5*ny + 1)*deltaX, -(0.5*nz + 1)*deltaX };
  std::vector<T> extend { deltaX, (ny + 2.)*deltaX, (nz + 2.)*deltaX };
  IndicatorCuboid3D<T> inlet(extend, origin);

  origin[0] = 0.5*nx*deltaX;
  IndicatorCuboid3D<T> outlet(extend, origin);

  IndicatorCuboid3D<T> seed( {deltaX, (ny + 2.)*deltaX, (nz + 2.)*deltaX},
                             {(-0.25*(fromCentre ? 1. : nx) - 0.5)*deltaX, -(0.5*ny + 1)*deltaX, -(0.5*nz + 1)*deltaX} );

  superGeometry.rename(0,2);
  superGeometry.rename( 2,1,{1,1,1} );
  superGeometry.rename(2,3,inlet);
  superGeometry.rename(2,3,outlet);
  superGeometry.rename( 1,5,{1,1,1} );
  superGeometry.rename(1,4,seed);
  superGeometry.rename(5,4,seed);

  superGeometry.checkForErrors();
  superGeometry.print();

  if (! noPlots) clout << "Prepare Geometry ... OK" << std::endl;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Set up the geometry of the simulation
void prepareLattice( SuperLattice<T, DESCRIPTOR_NS>& sLatticeNS,
                     SuperLattice<T, DESCRIPTOR_AD>& sLatticeAD,
                     UnitConverter<T,DESCRIPTOR_NS>& converter,
                     SuperGeometry<T,3>& superGeometry )
{
  /* MAT NUM | GEOMETRY     | NAVIER-STOKES | FINITE-DIFFERENCE | lATTICE-BOLTZMANN
   * 1       | Near-Border  | Bulk          | Bulk              | Bulk
   * 2       | Wall         | Free-slip     | No-penetration    | Bounce-back
   * 3       | Inlet-Outlet | Bulk          | No-penetration    | Bounce-back
   * 4       | Seed         | Bulk          | Bulk              | Bulk
   * 5       | Bulk         | Bulk          | Bulk              | Bulk
   */
  OstreamManager clout( std::cout,"prepareLattice" );

  if (! noPlots) clout << "Prepare Lattice ..." << std::endl;
  sLatticeNS.defineDynamics<BGKdynamics>(superGeometry.getMaterialIndicator({1,2,3,4,5}));
  sLatticeNS.setParameter<descriptors::OMEGA>(0.8);

  std::shared_ptr<FdModel<T,DESCRIPTOR_NS>> adModelBulk;
  std::shared_ptr<FdModel<T,DESCRIPTOR_NS>> adModelNearBorder;
  if (advectionScheme == 1) {
    adModelBulk = std::make_shared<FdAdvectionDiffusionModel<T,DESCRIPTOR_NS,SCHEME_ADV_CENTRAL,SCHEME_DIFF>> (
                  diffusivity * converter.getPhysDeltaT() / ( converter.getPhysDeltaX() * converter.getPhysDeltaX() ) );
    adModelNearBorder = std::make_shared<FdAdvectionDiffusionModel<T,DESCRIPTOR_NS,SCHEME_ADV_CENTRAL,SCHEME_DIFF>> (
                        diffusivity * converter.getPhysDeltaT() / ( converter.getPhysDeltaX() * converter.getPhysDeltaX() ) );
  }
  else if (advectionScheme == 2) {
    adModelBulk       = std::make_shared<FdAdvectionDiffusionModel<T,DESCRIPTOR_NS,SCHEME_ADV_UPWIND1,SCHEME_DIFF>> (
                        diffusivity * converter.getPhysDeltaT() / ( converter.getPhysDeltaX() * converter.getPhysDeltaX() ) );
    adModelNearBorder = std::make_shared<FdAdvectionDiffusionModel<T,DESCRIPTOR_NS,SCHEME_ADV_UPWIND1,SCHEME_DIFF>> (
                        diffusivity * converter.getPhysDeltaT() / ( converter.getPhysDeltaX() * converter.getPhysDeltaX() ) );
  }
  else if (advectionScheme == 3) {
    adModelBulk       = std::make_shared<FdAdvectionDiffusionModelWithAntiDiffusion<T,DESCRIPTOR_NS,SCHEME_ADV_UPWIND1,SCHEME_DIFF>> (
                        diffusivity * converter.getPhysDeltaT() / ( converter.getPhysDeltaX() * converter.getPhysDeltaX() ), antiDiffusionTuning );
    adModelNearBorder = std::make_shared<FdAdvectionDiffusionModelWithAntiDiffusion<T,DESCRIPTOR_NS,SCHEME_ADV_UPWIND1,SCHEME_DIFF>> (
                        diffusivity * converter.getPhysDeltaT() / ( converter.getPhysDeltaX() * converter.getPhysDeltaX() ), antiDiffusionTuning );
  }
  else if (advectionScheme == 4) {
    adModelBulk       = std::make_shared<FdAdvectionDiffusionModel<T,DESCRIPTOR_NS,SCHEME_ADV_UPWIND2,SCHEME_DIFF>> (
                        diffusivity * converter.getPhysDeltaT() / ( converter.getPhysDeltaX() * converter.getPhysDeltaX() ) );
    adModelNearBorder = std::make_shared<FdAdvectionDiffusionModel<T,DESCRIPTOR_NS,SCHEME_ADV_UPWIND1,SCHEME_DIFF>> (
                        diffusivity * converter.getPhysDeltaT() / ( converter.getPhysDeltaX() * converter.getPhysDeltaX() ) );
  }
  else {
    throw std::out_of_range("The order of the finite-difference scheme must only be 1 or 2.");
  }

  FdPostProcessorGenerator3D<T,DESCRIPTOR_NS,descriptors::AD_FIELD> adPostPbulk(iT, adModelBulk);
  FdPostProcessorGenerator3D<T,DESCRIPTOR_NS,descriptors::AD_FIELD> adPostPnearBorder(iT, adModelNearBorder);
  sLatticeNS.addPostProcessor(superGeometry, 1, adPostPnearBorder);
  sLatticeNS.addPostProcessor(superGeometry, 3, adPostPnearBorder);
  sLatticeNS.addPostProcessor(superGeometry, 4, adPostPbulk);
  sLatticeNS.addPostProcessor(superGeometry, 5, adPostPbulk);

  if (advectionScheme == 1) {
    setFdNeumannZeroBoundary<T,DESCRIPTOR_NS,SCHEME_ADV_CENTRAL,descriptors::AD_FIELD>(sLatticeNS, iT, adModelBulk, superGeometry, 2);
  }
  else if (advectionScheme == 2 || advectionScheme == 3 || advectionScheme == 4) {
    setFdNeumannZeroBoundary<T,DESCRIPTOR_NS,SCHEME_ADV_UPWIND1,descriptors::AD_FIELD>(sLatticeNS, iT, adModelBulk, superGeometry, 2);
  }
  else {
    throw std::out_of_range (
          "prepareLattice(...): The order of the finite-difference scheme must only be 1, 2 or 3, but instead is" + std::to_string(advectionScheme) );
  }

  if (! noPlots) clout << "Prepare Lattice ... OK" << std::endl;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Set up the inter-lattice coupling
void prepareCoupling ( SuperLattice<T, DESCRIPTOR_NS>& sLatticeNS,
                       SuperLattice<T, DESCRIPTOR_AD>& sLatticeAD,
                       SuperGeometry<T,3>& superGeometry )
{
  /* MAT NUM | GEOMETRY     | NAVIER-STOKES | FINITE-DIFFERENCE | lATTICE-BOLTZMANN
   * 1       | Near-Border  | Bulk          | Bulk              | Bulk
   * 2       | Wall         | Free-slip     | No-penetration    | Bounce-back
   * 3       | Inlet-Outlet | Bulk          | No-penetration    | Bounce-back
   * 4       | Seed         | Bulk          | Bulk              | Bulk
   * 5       | Bulk         | Bulk          | Bulk              | Bulk
   */
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Coupling ..." << std::endl;

  auto& commFields = sLatticeNS.getCommunicator(PostPostProcess());
  commFields.requestField<AD_FIELD>();
  commFields.requestOverlap(sLatticeNS.getOverlap());
  commFields.exchangeRequests();

  clout << "Prepare Coupling ... OK" << std::endl;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Sets values at the boundary and external velocity field following linear Couette at iT=0
void setBoundaryValues( SuperLattice<T, DESCRIPTOR_NS>& sLatticeNS,
                        SuperLattice<T, DESCRIPTOR_AD>& sLatticeAD,
                        UnitConverter<T,DESCRIPTOR_NS>& converter,
                        SuperGeometry<T,3>& superGeometry )
{
  /* MAT NUM | GEOMETRY     | NAVIER-STOKES | FINITE-DIFFERENCE | lATTICE-BOLTZMANN
   * 1       | Near-Border  | Bulk          | Bulk              | Bulk
   * 2       | Wall         | Free-slip     | No-penetration    | Bounce-back
   * 3       | Inlet-Outlet | Bulk          | No-penetration    | Bounce-back
   * 4       | Seed         | Bulk          | Bulk              | Bulk
   * 5       | Bulk         | Bulk          | Bulk              | Bulk
   */
  OstreamManager clout( std::cout,"setBoundaryValues" );

  AnalyticalConst<3,T,T> zero   {0.};
  AnalyticalConst<3,T,T> zeroLB {0.};
  AnalyticalConst<3,T,T> u0 {converter.getCharLatticeVelocity(), T(), T()};

  sLatticeNS.defineU( superGeometry, 2, u0 );
  sLatticeNS.defineU( superGeometry, 3, u0 );

  if (iT==0) {
    AnalyticalConst<3,T,T> one   {1.};
    AnalyticalConst<3,T,T> oneLB {1.};
    sLatticeNS.defineRhoU( superGeometry, 1, one, u0 );
    sLatticeNS.defineRhoU( superGeometry, 4, one, u0 );
    sLatticeNS.defineRhoU( superGeometry, 5, one, u0 );
    sLatticeNS.template defineField<AD_FIELD>( superGeometry, 1, zero );
    sLatticeNS.template defineField<AD_FIELD>( superGeometry, 2, zero );
    sLatticeNS.template defineField<AD_FIELD>( superGeometry, 3, zero );
    sLatticeNS.template defineField<AD_FIELD>( superGeometry, 4, one  );
    sLatticeNS.template defineField<AD_FIELD>( superGeometry, 5, zero );
    sLatticeNS.initialize();
  }
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Plots the results
void getResults( SuperLattice<T, DESCRIPTOR_NS>& sLatticeNS,
                 SuperLattice<T, DESCRIPTOR_AD>& sLatticeAD,
                 UnitConverter<T,DESCRIPTOR_NS>& converter,
                 SuperGeometry<T,3>& superGeometry, util::Timer<T>& timer )
{
  OstreamManager clout( std::cout,"getResults" );

  SuperVTMwriter3D<T> vtmWriterNS( "advectionDiffusionPipe3d" );
  SuperLatticeDensity3D<T, DESCRIPTOR_NS>            density( sLatticeNS );
  SuperLatticePhysVelocity3D<T, DESCRIPTOR_NS>       velocity( sLatticeNS, converter );
  SuperLatticeEul2LagrDensity3D<T, DESCRIPTOR_NS>    particleDensity( sLatticeNS );
  SuperLatticeExternal3D<T, DESCRIPTOR_NS, AD_FIELD> fdField( sLatticeNS, iT );
  SuperLatticeGeometry3D<T, DESCRIPTOR_NS>           geometry( sLatticeNS, superGeometry );

  SuperVTMwriter3D<T> vtmWriterAD( "writerAD" );
  SuperLatticeDensity3D<T, DESCRIPTOR_AD>            lbField( sLatticeAD );

#ifndef NOVTM
  vtmWriterNS.addFunctor( density );
  vtmWriterNS.addFunctor( velocity );
  vtmWriterNS.addFunctor( fdField );
#endif

  if ( iT == 0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeCuboid3D<T, DESCRIPTOR_NS> cuboid( sLatticeNS );
    SuperLatticeRank3D<T, DESCRIPTOR_NS> rank( sLatticeNS );
    vtmWriterNS.write( geometry );
    vtmWriterNS.write( cuboid );
    vtmWriterNS.write( rank );
    vtmWriterNS.createMasterFile();
  }

#ifndef NOVTM
  /// Vtm:
  if (iT % tVtm == 0) {
    clout << "Writing vtm files at iT=" << iT << std::endl;
    vtmWriterNS.write(iT);
    vtmWriterAD.write(iT);
  }
#endif

  /// ImageData:
  if (iT % tVtm == 0) {
    Vector<T,3> centre = {0., 0., 0.};
    Vector<T,3> u = {1.,0.,  0.};
    Vector<T,3> v = {0., 0., 1.};

    BlockReduction3D2D<T> planeReduction_fdField( fdField, centre, u, v, 200, BlockDataSyncMode::ReduceOnly );
    heatmap::plotParam<T> plotParam_fdField;
    plotParam_fdField.name = "fdField_";
    heatmap::write(planeReduction_fdField, iT);
    SuperEuklidNorm3D<T, DESCRIPTOR_NS> normVel(velocity);
    BlockReduction3D2D<T> planeReduction_normVel( normVel, centre, u, v, 200, BlockDataSyncMode::ReduceOnly );
    heatmap::plotParam<T> plotParam_normVel;
    plotParam_normVel.name = "normVel_";
    heatmap::write(planeReduction_normVel, iT);
  }

  /// Gnuplot:
  if (iT % tVtm == 0) {
    T sigma = util::sqrt(2*iT*diffusivity/(deltaX*deltaX/deltaT)); // dimensionless
    Gnuplot<T> gplot( "distribution_" + std::to_string(iT) );
    T analyticalMeanPointZero[] { -0.25*(fromCentre ? 0. : nx) - 0.5*(1 - nx%2), 0., 0. };
    T analyticalMeanPoint[] { analyticalMeanPointZero[0] + iT*converter.getCharLatticeVelocity(), 0., 0. };
    T numericalMeanPoint_fd {analyticalMeanPointZero[0]};
    T numericalMeanValue_fd {T()};
    int point[] {0, 0, 0};
    T pointP[] {0., 0., 0.};
    T analytical[] {0., 0., 0.};
    T numericalFd[] {0., 0., 0.};
    AnalyticalNormal<3,T,int> sol({analyticalMeanPoint[0], analyticalMeanPoint[1], analyticalMeanPoint[2]}, sigma);
    AnalyticalFfromSuperF3D<T> numFd( fdField, true, 1 );
    AnalyticalFfromSuperF3D<T> numLB( lbField, true, 1 );
    AnalyticalFfromSuperF3D<T> numPa( particleDensity, true, 1 );
    T normFd2 = T();
    T normalization = T();
    for (int iX=-(nx+1)/2+1; iX<=(nx+1)/2-1; ++iX) {
      point[0] = iX;
      pointP[0] = iX*deltaX;
      sol(analytical,point);
      if (diffusivity>0.0) normalization += analytical[0] * analytical[0];
      std::vector<T> gplot_values {analytical[0]};
      std::vector<std::string> gplot_names {"analytical"};
      numFd(numericalFd,pointP);
      numericalFd[0] *= (ny + 1.)*(nz + 1.) / (ny*nz);
      if (numericalFd[0] > numericalMeanValue_fd) {
        numericalMeanValue_fd = numericalFd[0];
        numericalMeanPoint_fd = pointP[0];
      }
      normFd2 += util::pow(numericalFd[0] - (diffusivity>0.0 ? analytical[0] : 0.0), 2 );
      gplot_values.push_back(numericalFd[0]);
      gplot_names.push_back("finite-difference");
      gplot.setData( iX*deltaX, gplot_values, gplot_names );
    }
    gplot.writePNG();
    if (diffusivity > 0.0) {
      clout << "norm2=" << util::sqrt(normFd2 / normalization)
            << "   group_velocity_error=" << (analyticalMeanPoint[0] - numericalMeanPoint_fd/converter.getPhysDeltaX())
                                           / (analyticalMeanPoint[0] - analyticalMeanPointZero[0])
            << std::endl
            << "mean_value_analytical=" << analyticalMeanPoint[0]*converter.getPhysDeltaX()
            << "   mean_value_numerical=" << numericalMeanPoint_fd // <---
            << std::endl;
    }
    // || f(x) - delta(x-x0) ||_L2 = util::sqrt [ || f ||_L2^2 + 1 - 2 * f(x0) ]
    else {
      numFd(numericalFd,analyticalMeanPoint);
      clout << "norm2=" << util::sqrt(normFd2 + 1.0 - 2.0*numericalFd[0])
            << "   group_velocity_error=" << (analyticalMeanPoint[0] - numericalMeanPoint_fd/converter.getPhysDeltaX())
                                           / (analyticalMeanPoint[0] - analyticalMeanPointZero[0])
            << std::endl;
    }
    timer.update(iT);
  }
}


///////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir("./tmp/");

  std::string xmlName( "input.xml" );
  if (argc>1) std::string xmlName = argv[1];

  XMLreader config(xmlName);

  // Boolean flags for advection-diffusion simulations
  config["advectionScheme"].read(advectionScheme);
  config["fromCentre"].read(fromCentre);

  // (dimensionless) Parameters for the simulation setup
  config["antiDiffusionTuning"].read(antiDiffusionTuning);
  config["nx"].read(nx);
  config["ny"].read(ny);
  config["tMax"].read(tMax);
  config["tVtm"].read(tVtm);

  // Dimensional parameters
  config["charU"].read(charU);
  config["deltaX"].read(deltaX);
  config["deltaT"].read(deltaT);
  config["diffusivity"].read(diffusivity);

  // Other parameters
  config["noPlots"].read(noPlots);

  nz = ny;
//  singleton::directories().setOutputDir( "./" + xmlName.substr(0, xmlName.find(".xml")) + "/" );
  singleton::directories().setOutputDir("./tmp/");
  OstreamManager clout( std::cout,"main" );

  UnitConverter<T,DESCRIPTOR_NS> converter (
    ( T )   deltaX,      // physDeltaX
    ( T )   deltaT,      // physDeltaT
    ( T )   ny*deltaX,   // charPhysLength
    ( T )   charU,       // charPhysVelocity
    ( T )   0.1 * (deltaX * deltaX) / deltaT,         // physViscosity
	  ( T )   1.           // physDensity
    );
  if (! noPlots) clout << "---------- Input data: ------------" << std::endl
        << "advectionScheme    = " << advectionScheme << std::endl
        << "nx         = " << nx << std::endl
        << "ny         = " << ny << std::endl
        << "tMax       = " << tMax << std::endl
        << "tVtm       = " << tVtm << std::endl
        << "charU      = " << charU << std::endl
        << "deltaX     = " << deltaX << std::endl
        << "deltaT     = " << deltaT << std::endl
        << "diffusivity= " << diffusivity << std::endl
        << "-----------------------------------" << std::endl;
  if (! noPlots) clout << "LB advection-diffusion relaxation time="
        << 1. / converter.template getLatticeRelaxationFrequencyFromDiffusivity<DESCRIPTOR_AD> (
           diffusivity ) << std::endl;
           //diffusivity * (deltaX*deltaX/deltaT) ) << std::endl;
  if (! noPlots) clout << "---------- Stability criteria: ------------" << std::endl
        << "Peclet number:          Pe = " << ny * deltaX * charU / diffusivity << "  (it should be not too big for Lattice-Boltzmann)" << std::endl
        << "Diffusivity-spacing-timestep:   deltaX * deltaX / (2. * diffusivity * deltaT) = " << deltaX * deltaX / (2. * diffusivity * deltaT) << "  (it should be > 1 for finite-difference)" << std::endl
        << "Velocity-diffusivity-spacing:   2. * diffusivity / (deltaX * charU) = " << 2. * diffusivity / (deltaX * charU) << "  (it should be > 1 for 2nd-order finite-difference)" << std::endl
        << "3sigma: deltaX / (3. * sqrt(2. * diffusivity * deltaT)) = " << deltaX / (3. * util::sqrt(2. * diffusivity * deltaT)) << std::endl
        << "-------------------       -----------------" << std::endl;
  if (! noPlots) converter.print();

  /// === 2rd Step: Prepare Geometry ===
  std::vector<T> origin { -0.5*(nx + 1.)*deltaX, -0.5*(ny + 1.)*deltaX, -0.5*(nz + 1.)*deltaX };
  std::vector<T> extend { (nx + 1)*deltaX, (ny + 1)*deltaX, (nz + 1)*deltaX };
  IndicatorCuboid3D<T> cuboid(extend, origin);

  /// Instantiation of a cuboidGeometry with weights
#ifdef DEBUG_LINEAR
  const int noOfCuboids = 8;
#else
#ifdef PARALLEL_MODE_MPI
  //const int noOfCuboids = 2 * singleton::mpi().getSize();
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif
#endif
  CuboidGeometry3D<T> cuboidGeometry(cuboid, deltaX, noOfCuboids);
  cuboidGeometry.setPeriodicity(false, false, false);
  cuboidGeometry.print();

  /// Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

  /// Instantiation of a superGeometry
  SuperGeometry<T,3> superGeometry(cuboidGeometry, loadBalancer, 2);

  prepareGeometry(superGeometry);

  /// === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR_NS> sLatticeNS( superGeometry );
  SuperLattice<T, DESCRIPTOR_AD> sLatticeAD( superGeometry );
  prepareLattice( sLatticeNS, sLatticeAD, converter, superGeometry );

  // Communicator
  prepareCoupling( sLatticeNS, sLatticeAD, superGeometry);

  // === 4th Step: Main Loop with Timer ===
  if (! noPlots) clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( tMax, superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for ( iT = 0; iT <= tMax; ++iT ) {
    // === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues( sLatticeNS, sLatticeAD, converter, superGeometry );

    // === 6th Step: Collide and Stream Execution ===
    sLatticeNS.collideAndStream();

    // === 7th Step: Computation and Output of the Results ===
    if (! noPlots ) getResults( sLatticeNS, sLatticeAD, converter, superGeometry, timer );
  }

  if (noPlots) timer.update(iT);
  timer.stop();
  timer.printSummary();
}
