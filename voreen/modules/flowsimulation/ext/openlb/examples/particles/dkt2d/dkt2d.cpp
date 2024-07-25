/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2006-2021 Nicolas Hafen, Robin Trunk, Mathias J. Krause
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



/* dkt2d.cpp:
 * The case examines the settling of two circles under gravity
 * in a surrounding fluid. The rectangular domain is limited
 * by no-slip boundary conditions.
 * For the calculation of forces a DNS approach is chosen
 * which also leads to a back-coupling of the particle on the fluid,
 * inducing a flow.
 * The simulation is based on the homogenised lattice Boltzmann approach
 * (HLBM) introduced by Krause et al. in "Particle flow simulations
 * with homogenised lattice Boltzmann methods".
 * The drafting-kissing-tumbling benchmark case is e.g. described
 * in "Drafting, kissing and tumbling process of two particles
 * with different sizes" by Wang et al.
 * or "The immersed boundary-lattice Boltzmann method
 * for solving fluid-particles interaction problems" by Feng and Michaelides.
 * The example demonstrates the usage of HLBM in the OpenLB framework
 * as well as the utilisation of the Gnuplot-writer
 * to print simulation results.
 */

#define NEW_FRAMEWORK


#include "olb2D.h"
#include "olb2D.hh"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace olb::particles;
using namespace olb::particles::dynamics;
using namespace olb::util;

typedef double T;

//Define lattice type
typedef PorousParticleD2Q9Descriptor DESCRIPTOR;

//Define particle type
typedef ResolvedCircle2D PARTICLETYPE;

#define WriteVTK
#define WriteGnuPlot

std::string gnuplotFilename = "gnuplot.dat";

// Parameters for the simulation setup
int N = 1;
int M = N;

T eps = 0.5;      // eps*latticeL: width of transition area

T maxPhysT = 6.;  // max. simulation time in s, SI unit
T iTwrite = 0.125;  //converter.getLatticeTime(.3);

T lengthX = 0.02;
T lengthY = 0.08;

T centerX1 = 0.01;
T centerY1 = 0.068;
Vector<T,2> center1 = {centerX1,centerY1};
T centerX2 = 0.00999;
T centerY2 = 0.072;
Vector<T,2> center2 = {centerX2,centerY2};

T rhoP = 1010.;
T radiusP = 0.001;
Vector<T,2> accExt = {.0, -9.81 * (1. - 1000. / rhoP)};

void prepareGeometry(UnitConverter<T,DESCRIPTOR> const& converter,
                     SuperGeometry<T,2>& superGeometry)
{
  OstreamManager clout(std::cout, "prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0, 2);
  superGeometry.rename(2, 1, {1, 1});

  superGeometry.clean();
  superGeometry.innerClean();

  superGeometry.checkForErrors();
  superGeometry.getStatistics().print();
  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice(
  SuperLattice<T, DESCRIPTOR>& sLattice, UnitConverter<T,DESCRIPTOR> const& converter,
  SuperGeometry<T,2>& superGeometry)
{
  OstreamManager clout(std::cout, "prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  /// Material=0 -->do nothing
  sLattice.defineDynamics<NoDynamics>(superGeometry, 0);
  sLattice.defineDynamics<PorousParticleBGKdynamics>(superGeometry, 1);
  sLattice.defineDynamics<BounceBack>(superGeometry, 2);

  sLattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());

  clout << "Prepare Lattice ... OK" << std::endl;
}

void setBoundaryValues(SuperLattice<T, DESCRIPTOR>& sLattice,
                       UnitConverter<T,DESCRIPTOR> const& converter,
                       SuperGeometry<T,2>& superGeometry)
{
  OstreamManager clout(std::cout, "setBoundaryValues");

  AnalyticalConst2D<T, T> one(1.);
  sLattice.defineField<POROSITY>(superGeometry.getMaterialIndicator({1,2}), one);

  // Set initial condition
  AnalyticalConst2D<T, T> ux(0.);
  AnalyticalConst2D<T, T> uy(0.);
  AnalyticalConst2D<T, T> rho(1.);
  AnalyticalComposed2D<T, T> u(ux, uy);

  //Initialize all values of distribution functions to their local equilibrium
  sLattice.defineRhoU(superGeometry, 1, rho, u);
  sLattice.iniEquilibrium(superGeometry, 1, rho, u);

  // Make the lattice ready for simulation
  sLattice.initialize();
}

void getResults(SuperLattice<T, DESCRIPTOR>& sLattice,
                UnitConverter<T,DESCRIPTOR> const& converter, int iT,
                SuperGeometry<T,2>& superGeometry, Timer<double>& timer,
                ParticleSystem<T,PARTICLETYPE>& particleSystem )
{
  OstreamManager clout(std::cout, "getResults");

#ifdef WriteVTK
  SuperVTMwriter2D<T> vtkWriter("sedimentation");
  SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity(sLattice, converter);
  SuperLatticePhysPressure2D<T, DESCRIPTOR> pressure(sLattice, converter);
  SuperLatticePhysExternalPorosity2D<T, DESCRIPTOR> externalPor(sLattice, converter);
  SuperLatticeMomentumExchangeForceLocal<T, DESCRIPTOR, PARTICLETYPE> momentumExchange(
    sLattice, converter, superGeometry, particleSystem);

  vtkWriter.addFunctor(velocity);
  vtkWriter.addFunctor(pressure);
  vtkWriter.addFunctor(externalPor);
  vtkWriter.addFunctor(momentumExchange);

  if (iT == 0) {
    converter.write("dkt");
    SuperLatticeGeometry2D<T, DESCRIPTOR> geometry(sLattice, superGeometry);
    SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid(sLattice);
    SuperLatticeRank2D<T, DESCRIPTOR> rank(sLattice);
    vtkWriter.write(geometry);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);
    vtkWriter.createMasterFile();
  }

  if (iT % converter.getLatticeTime(iTwrite) == 0) {
    vtkWriter.write(iT);
  }
#endif


  auto particleA = particleSystem.get( 0 );
  auto particleB = particleSystem.get( 1 );

#ifdef WriteGnuPlot
  if (iT % converter.getLatticeTime(iTwrite) == 0) {
    if (singleton::mpi().getRank() == 0) {

      std::ofstream myfile;
      myfile.open (gnuplotFilename.c_str(), std::ios::app);
      T p2PosY = particleB.getField<GENERAL,POSITION>()[1];
      T p1PosY = particleA.getField<GENERAL,POSITION>()[1];
      T p2PosX = particleB.getField<GENERAL,POSITION>()[0];
      T p1PosX = particleA.getField<GENERAL,POSITION>()[0];
      myfile
          << converter.getPhysTime(iT) << " "
          << std::setprecision(9)
          << p2PosY << " "
          << p1PosY << " "
          << p2PosX << " "
          << p1PosX << std::endl;
      myfile.close();
    }
  }
#endif

  /// Writes output on the console
  if (iT % converter.getLatticeTime(iTwrite) == 0) {
    timer.update(iT);
    timer.printStep();
    sLattice.getStatistics().print(iT, converter.getPhysTime(iT));
    for (std::size_t iP=0; iP<particleSystem.size(); ++iP) {
      auto particle = particleSystem.get(iP);
      io::printResolvedParticleInfo(particle);
    }
  }

  return;
}

int main(int argc, char* argv[])
{
  /// === 1st Step: Initialization ===
  olbInit(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");
  OstreamManager clout(std::cout, "main");

  UnitConverter<T,DESCRIPTOR> converter(
    ( T )   0.0001/ N, //physDeltaX
    ( T )   5.e-4/(N*M), //physDeltaT,
    ( T )   .002, //charPhysLength
    ( T )   0.2, //charPhysVelocity
    ( T )   1E-6, //physViscosity
    ( T )   1000. //physDensity
  );
  converter.print();

  /// === 2nd Step: Prepare Geometry ===
  std::vector<T> extend(2, T());
  extend[0] = lengthX;
  extend[1] = lengthY;
  std::vector<T> origin(2, T());
  IndicatorCuboid2D<T> cuboid(extend, origin);

#ifdef PARALLEL_MODE_MPI
  CuboidGeometry2D<T> cuboidGeometry(cuboid, converter.getConversionFactorLength(), singleton::mpi().getSize());
#else
  CuboidGeometry2D<T> cuboidGeometry(cuboid, converter.getConversionFactorLength(), 1);
#endif

  HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);
  SuperGeometry<T,2> superGeometry(cuboidGeometry, loadBalancer, 2);
  prepareGeometry(converter, superGeometry);

  /// === 3rd Step: Prepare Lattice ===
  SuperLattice<T, DESCRIPTOR> sLattice(superGeometry);

  prepareLattice(sLattice, converter, superGeometry);

  /// === 4th Step: Main Loop with Timer ===
  Timer<double> timer(converter.getLatticeTime(maxPhysT), superGeometry.getStatistics().getNvoxel());
  timer.start();

  // Create ParticleSystem
  ParticleSystem<T,PARTICLETYPE> particleSystem;

  //Create particle manager handling coupling, gravity and particle dynamics
  ParticleManager<T,DESCRIPTOR,PARTICLETYPE> particleManager(
    particleSystem, superGeometry, sLattice, converter, accExt);

  T epsilon = eps*converter.getConversionFactorLength();
  T radius = radiusP;

  // Create Particle 1
  creators::addResolvedCircle2D( particleSystem, center1,
                                 radius, epsilon, rhoP );

  // Create Particle 2
  creators::addResolvedCircle2D( particleSystem, center2,
                                 radius, epsilon, rhoP );

  // Create and assign resolved particle dynamics
  VerletParticleDynamics<T,PARTICLETYPE> particleDynamics;
  for (std::size_t iP=0; iP<particleSystem.size(); ++iP) {
    particleSystem.defineDynamics( iP, &particleDynamics );
  }

  /// === 5th Step: Definition of Initial and Boundary Conditions ===
  setBoundaryValues(sLattice, converter, superGeometry);

  {
    auto& communicator = sLattice.getCommunicator(PostPostProcess());
    communicator.requestOverlap(sLattice.getOverlap());
    communicator.requestFields<POROSITY,VELOCITY_NUMERATOR,VELOCITY_DENOMINATOR>();
    communicator.exchangeRequests();
  }

  clout << "MaxIT: " << converter.getLatticeTime(maxPhysT) << std::endl;
  for (std::size_t iT = 0; iT < converter.getLatticeTime(maxPhysT)+10; ++iT) {

    // Execute particle manager
    particleManager.execute<
      couple_lattice_to_particles<T,DESCRIPTOR,PARTICLETYPE>,
      apply_gravity<T,PARTICLETYPE>,
      process_dynamics<T,PARTICLETYPE>,
      couple_particles_to_lattice<T,DESCRIPTOR,PARTICLETYPE>
    >();

    // Get Results
    getResults(sLattice, converter, iT, superGeometry, timer, particleSystem);

    // Collide and stream
    sLattice.collideAndStream();
  }

  // Run Gnuplot
  if (singleton::mpi().getRank() == 0) {
    if (!system(NULL)) {
      exit (EXIT_FAILURE);
    }
    int ret = system("gnuplot dkt.p");
    if (ret == -1) {
      clout << "Writing Gnuplot failed!" << std::endl;
    }
  }

  timer.stop();
  timer.printSummary();
}
