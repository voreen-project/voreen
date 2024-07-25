/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2008 Orestis Malaspinas, Andrea Parmigiani
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

/* contactAngle3d.cpp
 * In this example a semi-spherical droplet of fluid is initialised
 * within a different fluid at a solid boundary. The contact angle
 * is measured as the droplet comes to equilibrium. This is compared
 * with the analytical angle (100 degrees) predicted by the
 * parameters set for the boundary.
 *
 * This example demonstrates how to use the wetting solid boundaries
 * for the free-energy model with two fluid components.
 */

#include "olb3D.h"
#include "olb3D.hh"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

typedef double T;
typedef D3Q19<CHEM_POTENTIAL,FORCE> DESCRIPTOR;

// Parameters for the simulation setup
const int N  = 75;
const T nxy  = 75.;
const T nz   = 50.;
const T radius = 0.25 * nxy;

const T alpha = 1.;      // Interfacial width         [lattice units]
const T kappa1 = 0.005;  // For surface tensions      [lattice units]
const T kappa2 = 0.005;  // For surface tensions      [lattice units]
const T gama = 10.;      // For mobility of interface [lattice units]
const T h1 =  0.0001448; // Contact angle  80 degrees [lattice units]
const T h2 = -0.0001448; // Contact angle 100 degrees [lattice units]

const int maxIter = 10000;
const int vtkIter  = 200;
const int statIter = 200;
const bool calcAngle = true;

T angle_prev = 90.;


void prepareGeometry( SuperGeometry<T,3>& superGeometry,
                      UnitConverter<T, DESCRIPTOR>& converter)
{

  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename( 0,2 );

  Vector<T,3> extend(nxy+2., nxy+2., nz-1.*converter.getPhysDeltaX() );
  Vector<T,3> origin( -1., -1., 0.5*converter.getPhysDeltaX() );
  IndicatorCuboid3D<T> inner ( extend, origin );
  superGeometry.rename( 2,1,inner );

  superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}


void prepareLattice( SuperLattice<T, DESCRIPTOR>& sLattice1,
                     SuperLattice<T, DESCRIPTOR>& sLattice2,
                     UnitConverter<T, DESCRIPTOR>& converter,
                     SuperGeometry<T,3>& superGeometry)
{

  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  // Define lattice Dynamics
  sLattice1.defineDynamics<NoDynamics>(superGeometry, 0);
  sLattice2.defineDynamics<NoDynamics>(superGeometry, 0);
  sLattice1.defineDynamics<ForcedBGKdynamics>(superGeometry, 1);
  sLattice2.defineDynamics<FreeEnergyBGKdynamics>( superGeometry, 1);
  sLattice1.defineDynamics<NoDynamics>(superGeometry, 2);
  sLattice2.defineDynamics<NoDynamics>(superGeometry, 2);

  // Add wall boundary

  setFreeEnergyWallBoundary<T,DESCRIPTOR>(sLattice1, superGeometry, 2, alpha, kappa1, kappa2, h1, h2, 1);
  setFreeEnergyWallBoundary<T,DESCRIPTOR>(sLattice2, superGeometry, 2, alpha, kappa1, kappa2, h1, h2, 2);


  // Bulk initial conditions
  // Define spherical domain for fluid 2
  std::vector<T> v( 3,T() );
  AnalyticalConst3D<T,T> zeroVelocity( v );

  AnalyticalConst3D<T,T> one( 1.0 );
  IndicatorSphere3D<T> sphere( {nxy/2., nxy/2., 0.}, radius );
  SmoothIndicatorSphere3D<T,T> smoothSphere( sphere, 10.*alpha );

  AnalyticalIdentity3D<T,T> rho( one );
  AnalyticalIdentity3D<T,T> phi( one - smoothSphere - smoothSphere );

  sLattice1.iniEquilibrium( superGeometry, 1, rho, zeroVelocity );
  sLattice2.iniEquilibrium( superGeometry, 1, phi, zeroVelocity );

  sLattice1.iniEquilibrium( superGeometry, 2, rho, zeroVelocity );
  sLattice2.iniEquilibrium( superGeometry, 2, phi, zeroVelocity );

  sLattice1.setParameter<descriptors::OMEGA>( converter.getLatticeRelaxationFrequency() );
  sLattice2.setParameter<descriptors::OMEGA>( converter.getLatticeRelaxationFrequency() );

  sLattice1.initialize();
  sLattice2.initialize();

  sLattice1.communicate();
  sLattice2.communicate();

  {
    auto& communicator = sLattice1.getCommunicator(PostPostProcess());
    communicator.requestField<POPULATION>();
    communicator.requestOverlap(sLattice1.getOverlap());
    communicator.exchangeRequests();
  }
  {
    auto& communicator = sLattice2.getCommunicator(PostPostProcess());
    communicator.requestField<POPULATION>();
    communicator.requestOverlap(sLattice2.getOverlap());
    communicator.exchangeRequests();
  }

  clout << "Prepare Lattice ... OK" << std::endl;
}


void prepareCoupling( SuperLattice<T, DESCRIPTOR>& sLattice1,
                      SuperLattice<T, DESCRIPTOR>& sLattice2,
                      SuperGeometry<T,3>& superGeometry )
{

  OstreamManager clout( std::cout,"prepareCoupling" );
  clout << "Add lattice coupling" << std::endl;

  // Add the lattice couplings (not to the solid nodes)
  // The chemical potential coupling must come before the force coupling
  FreeEnergyChemicalPotentialGenerator3D<T, DESCRIPTOR> coupling1(
    alpha, kappa1, kappa2);
  FreeEnergyForceGenerator3D<T, DESCRIPTOR> coupling2;

  // Suppress compiler warnings
  coupling1.shift(0, 0, 0);
  coupling2.shift(0, 0, 0);

  sLattice1.addLatticeCoupling( superGeometry, 1, coupling1, sLattice2 );
  sLattice2.addLatticeCoupling( superGeometry, 1, coupling2, sLattice1 );

  {
    auto& communicator = sLattice1.getCommunicator(PostCoupling());
    communicator.requestField<CHEM_POTENTIAL>();
    communicator.requestOverlap(sLattice1.getOverlap());
    communicator.exchangeRequests();
  }
  {
    auto& communicator = sLattice2.getCommunicator(PreCoupling());
    communicator.requestField<CHEM_POTENTIAL>();
    communicator.requestOverlap(sLattice2.getOverlap());
    communicator.exchangeRequests();
  }

  clout << "Add lattice coupling ... OK!" << std::endl;
}


void getResults( SuperLattice<T, DESCRIPTOR>& sLattice1,
                 SuperLattice<T, DESCRIPTOR>& sLattice2, int iT,
                 SuperGeometry<T,3>& superGeometry, util::Timer<T>& timer,
                 UnitConverter<T, DESCRIPTOR> converter )
{

  OstreamManager clout( std::cout,"getResults" );
  SuperVTMwriter3D<T> vtmWriter( "contactAngle3d" );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry3D<T, DESCRIPTOR> geometry( sLattice1, superGeometry );
    SuperLatticeCuboid3D<T, DESCRIPTOR> cuboid( sLattice1 );
    SuperLatticeRank3D<T, DESCRIPTOR> rank( sLattice1 );
    vtmWriter.write( geometry );
    vtmWriter.write( cuboid );
    vtmWriter.write( rank );
    vtmWriter.createMasterFile();
  }

  // Get statistics
  if ( iT%statIter==0 ) {
    // Timer console output
    timer.update( iT );
    timer.printStep();
    sLattice1.getStatistics().print( iT, converter.getPhysTime(iT) );
    sLattice2.getStatistics().print( iT, converter.getPhysTime(iT) );
  }

  // Writes the VTK files
  if ( iT%vtkIter==0 ) {
    AnalyticalConst3D<T,T> half_( 0.5 );
    SuperLatticeFfromAnalyticalF3D<T, DESCRIPTOR> half(half_, sLattice1);

    SuperLatticeVelocity3D<T, DESCRIPTOR> velocity( sLattice1 );
    SuperLatticeDensity3D<T, DESCRIPTOR> rho( sLattice1 );
    rho.getName() = "rho";
    SuperLatticeDensity3D<T, DESCRIPTOR> phi( sLattice2 );
    phi.getName() = "phi";

    SuperIdentity3D<T,T> c1 (half*(rho+phi));
    c1.getName() = "density-fluid-1";
    SuperIdentity3D<T,T> c2 (half*(rho-phi));
    c2.getName() = "density-fluid-2";

    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( rho );
    vtmWriter.addFunctor( phi );
    vtmWriter.addFunctor( c1 );
    vtmWriter.addFunctor( c2 );
    vtmWriter.write( iT );

    // Evaluation of contact angle
    if (calcAngle) {
      int Nz = (int)( N * nz / nxy );
      T dx = converter.getPhysDeltaX();
      AnalyticalFfromSuperF3D<T,T> interpolPhi( phi, true, 1 );

      T base1 = 0.;
      T base2 = 0.;
      T height1 = 0.;
      T height2 = 0.;

      double pos[3] = {0., nxy/2., dx};
      for (int ix=0; ix<N; ix++) {
        T phi1, phi2;
        pos[0] = ix * dx;
        interpolPhi( &phi1, pos );
        if (phi1 < 0.) {
          pos[0] = (ix-1) * dx;
          interpolPhi( &phi2, pos );
          base1 = 2. * ( 0.5*N - ix + phi1/(phi1-phi2) );
          break;
        }
      }

      pos[2] = 3.*dx;
      for (int ix=0; ix<N; ix++) {
        T phi1, phi2;
        pos[0] = ix * dx;
        interpolPhi( &phi1, pos );
        if (phi1 < 0.) {
          pos[0] = (ix-1) * dx;
          interpolPhi( &phi2, pos );
          base2 = 2. * ( 0.5*N - ix + phi1/(phi1-phi2) );
          break;
        }
      }

      pos[0] = nxy / 2.;
      for (int iz=2; iz<Nz; iz++) {
        T phi1, phi2;
        pos[2] = iz * dx;
        interpolPhi( &phi1, pos );
        if (phi1 > 0.) {
          pos[2] = (iz-1) * dx;
          interpolPhi( &phi2, pos );
          height1 = iz - 1. - phi1/(phi1-phi2);
          height2 = iz - 3. - phi1/(phi1-phi2);
          break;
        }
      }

      // Calculate simulated contact angle
      T pi = 3.14159265;
      T height = height1 + 1.;
      T base = base1 + 2 * (radius - height1) / base1;
      T radius = (4.*height2*height2 + base2*base2) / ( 8.*height2 );
      T angle_rad = pi + util::atan( 0.5*base / (radius - height) );
      T angle = angle_rad * 180. / pi;
      if ( angle > 180. ) {
        angle -= 180.;
      }

      // Calculate theoretical contact angle
      T ak1 = alpha * kappa1;
      T ak2 = alpha * kappa2;
      T k12 = kappa1 + kappa2;
      T num1 = util::pow(ak1 + 4 * h1, 1.5) - util::pow(ak1 - 4 * h1, 1.5);
      T num2 = util::pow(ak2 + 4 * h2, 1.5) - util::pow(ak2 - 4 * h2, 1.5);
      T angle_an = 180 / pi * util::acos(num2 / (2 * k12 * util::sqrt(ak2)) - \
                                         num1 / (2 * k12 * util::sqrt(ak1)));

      clout << "----->>>>> Contact angle: " << angle << " ; ";
      clout << "Analytical contact angle: " << angle_an <<  std::endl;
      clout << "----->>>>> Difference to previous: " << angle-angle_prev << std::endl;
      angle_prev = angle;
    }
  }
}


int main( int argc, char *argv[] )
{

  // === 1st Step: Initialization ===

  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );

  UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR> converter(
    (T)   N, // resolution
    (T)   1., // lattice relaxation time (tau)
    (T)   nxy, // charPhysLength: reference length of simulation geometry
    (T)   0.0001, // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   1.002e-8, // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   1. // physDensity: physical density in __kg / m^3__
  );

  // Prints the converter log as console output
  converter.print();

  // === 2nd Step: Prepare Geometry ===
  std::vector<T> extend = { nxy, nxy, nz };
  std::vector<T> origin = { 0., 0., 0. };
  IndicatorCuboid3D<T> cuboid(extend,origin);
#ifdef PARALLEL_MODE_MPI
  CuboidGeometry3D<T> cGeometry( cuboid, converter.getPhysDeltaX(), singleton::mpi().getSize() );
#else
  CuboidGeometry3D<T> cGeometry( cuboid, converter.getPhysDeltaX() );
#endif

  // Set periodic boundaries to the domain
  cGeometry.setPeriodicity( true, true, false );

  // Instantiation of loadbalancer
  HeuristicLoadBalancer<T> loadBalancer( cGeometry );
  loadBalancer.print();

  // Instantiation of superGeometry
  SuperGeometry<T,3> superGeometry( cGeometry,loadBalancer );

  prepareGeometry( superGeometry, converter );

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice1( superGeometry );
  SuperLattice<T, DESCRIPTOR> sLattice2( superGeometry );

  //prepareLattice and set boundaryConditions
  prepareLattice( sLattice1, sLattice2, converter, superGeometry);

  prepareCoupling( sLattice1, sLattice2, superGeometry);

  // === 4th Step: Main Loop with Timer ===
  int iT = 0;
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( maxIter, superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for ( iT=0; iT<=maxIter; ++iT ) {
    // Computation and output of the results
    getResults( sLattice1, sLattice2, iT, superGeometry, timer, converter );

    // Collide and stream execution
    sLattice1.collideAndStream();
    sLattice2.collideAndStream();

    // Execute coupling between the two lattices
    sLattice1.executeCoupling();
    sLattice2.executeCoupling();
  }

  timer.stop();
  timer.printSummary();

}

