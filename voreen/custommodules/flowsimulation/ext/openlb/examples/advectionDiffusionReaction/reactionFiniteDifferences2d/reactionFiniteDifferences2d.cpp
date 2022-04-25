/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2006-2021 Davide Dapelo
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

/* advectionDiffusionReaction2d:
 * simulating a simple domain with no fluid motion and homogeneous
 * species concentration. Reproducing the chemical reaction |a|A -> |b|B:
 * reaction rate nu = [A]/t0;
 * initial conditions [A](t=0)=1; [B](t=0)=0.
 * Analytical solution:
 * [A](t) = util::exp( -|a|*t/t0 );
 * [B](t) = |b/a|*( 1 - util::exp( -|a|*t/t0 ) ).
 * t0 is a time conversion factor.
 */

#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>

#include "olb2D.h"
#include "olb2D.hh"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

namespace olb {
namespace descriptors {
struct A_FIELD  : public FIELD_BASE<2,  0, 0> { };
struct A_SOURCE : public FIELD_BASE<1,  0, 0> { };
struct B_FIELD  : public FIELD_BASE<2,  0, 0> { };
struct B_SOURCE : public FIELD_BASE<1,  0, 0> { };
}
}

typedef double T;
typedef D2Q9<A_FIELD,A_SOURCE,B_FIELD,B_SOURCE> DESCRIPTOR;
typedef fd::tag::UPWIND  SCHEME_ADV;
typedef fd::tag::CENTRAL SCHEME_DIFF;


///////////////////////////////////////////////////////////////////////////////////////////////////
// Functional for the analytical solution(s)
template <unsigned D, typename T, typename S>
class AnalyticalSol final: public AnalyticalF<D,T,S> {
private:
  T _t, _t0, _a, _b;
public:
  AnalyticalSol(T t, T t0, T a, T b) : AnalyticalF<D,T,S>(2), _t(t), _t0(t0), _a(a), _b(b)
  {
    this->getName() = "sol";
  }
  bool operator() (T output[], const S x[]) override
  {
    output[0] = util::exp(_a*_t/_t0);
    output[1] = -_b/_a*(1. - util::exp(_a*_t/_t0));
    return true;
  }
};


///////////////////////////////////////////////////////////////////////////////////////////////////
std::size_t iT   = 0; // global timestep

// (dimensionless) Parameters for the simulation setup
int nx   = 50;  // Number of internal lattice points along the x direction
int ny   = 50;   // Number of internal lattice points along the y direction
std::size_t tMax = 2000; // Total number of lattice updates
std::size_t tVtm = 100;  // Number of timesteps before producing output
std::size_t tGnu = 50; // Number of timesteps before producing Gnuplot output
T   a    = -3.; // A specie's stoichiometric coefficient. Negative because it is the reagent.
T   b    =  2.; // B specie's stoichiometric coefficient. Positive because it is the product.

// Dimensional parameters
T deltaX = 0.1;  // Lattice spacing (m)
T deltaT = 0.01; // Timestep (s)
T t0     = 10;  // time conversion factor for reaction rate


///////////////////////////////////////////////////////////////////////////////////////////////////
// Stores geometry information in form of material numbers
void prepareGeometry(SuperGeometry<T,2>& superGeometry)
{
  /* MAT NUM | GEOMETRY     | FINITE-DIFFERENCE | lATTICE-BOLTZMANN
   * 1       | Bulk         | Bulk              | Bulk
   * 2       | Wall         | No-penetration    | Bulk
   */
  OstreamManager clout(std::cout,"prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0,2);
  superGeometry.rename(2,1,{1,1});

  superGeometry.checkForErrors();
  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Set up the geometry of the simulation
void prepareLattice( SuperLattice<T,DESCRIPTOR>& sLattice,
                     UnitConverter<T,DESCRIPTOR>& converter,
                     SuperGeometry<T,2>& superGeometry )
{
  /* MAT NUM | GEOMETRY     | FINITE-DIFFERENCE | lATTICE-BOLTZMANN
   * 1       | Bulk         | Bulk              | Bulk
   * 2       | Wall         | No-penetration    | Bulk
   */
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  auto adModel = std::make_shared<FdAdvectionDiffusionModel<T,DESCRIPTOR,SCHEME_ADV,SCHEME_DIFF>>(0.);

  ReactionGenerator2D<T,DESCRIPTOR> reaction (
  {  new ExpOn1stSpecieRate<T>(converter.getLatticeTime(t0)) }, {
    {
      new FiniteDifferenceReactingSpecies2D<T,DESCRIPTOR,A_FIELD,A_SOURCE> (superGeometry, sLattice, a, iT),
      new FiniteDifferenceReactingSpecies2D<T,DESCRIPTOR,B_FIELD,B_SOURCE> (superGeometry, sLattice, b, iT)
    }
  } );
  FdPostProcessorGenerator2D<T,DESCRIPTOR,A_FIELD,A_SOURCE> adPostP_A(iT, adModel);
  FdPostProcessorGenerator2D<T,DESCRIPTOR,B_FIELD,B_SOURCE> adPostP_B(iT, adModel);

  sLattice.addPostProcessor(reaction);
  sLattice.addPostProcessor(superGeometry, 1, adPostP_A);
  sLattice.addPostProcessor(superGeometry, 1, adPostP_B);

  setFdNeumannZeroBoundary<T,DESCRIPTOR,SCHEME_ADV,A_FIELD,A_SOURCE>(sLattice, iT, adModel, superGeometry, 2);
  setFdNeumannZeroBoundary<T,DESCRIPTOR,SCHEME_ADV,B_FIELD,B_SOURCE>(sLattice, iT, adModel, superGeometry, 2);

  sLattice.defineDynamics<BGKdynamics>(superGeometry.getMaterialIndicator({1,2}));
  sLattice.setParameter<descriptors::OMEGA>(converter.getLatticeRelaxationFrequency());

  auto& commFields = sLattice.getCommunicator(PostPostProcess());
  commFields.requestField<A_FIELD>();
  commFields.requestField<A_SOURCE>();
  commFields.requestField<B_FIELD>();
  commFields.requestField<B_SOURCE>();
  commFields.requestOverlap(sLattice.getOverlap());
  commFields.exchangeRequests();
  clout << "Prepare Lattice ... OK" << std::endl;
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Sets values at the boundary and external velocity field following linear Couette at iT=0
void setBoundaryValues( SuperLattice<T, DESCRIPTOR>& sLattice,
                        UnitConverter<T,DESCRIPTOR>& converter,
                        SuperGeometry<T,2>& superGeometry )
{
  /* MAT NUM | GEOMETRY     | FINITE-DIFFERENCE | lATTICE-BOLTZMANN
   * 1       | Bulk         | Bulk              | Bulk
   * 2       | Wall         | No-penetration    | Bulk
   */
  OstreamManager clout( std::cout,"setBoundaryValues" );

  AnalyticalConst<2,T,T> one1  {1.};
  AnalyticalConst<2,T,T> zero2 {0., 0.};
  AnalyticalConst<2,T,T> one2  {1., 1.};
  AnalyticalConst<2,T,T> u0 {T(), T()};

  sLattice.defineRhoU( superGeometry, 0, one1, u0 );
  sLattice.template defineField<A_FIELD> ( superGeometry, 1, one2  );
  sLattice.template defineField<A_FIELD> ( superGeometry, 2, one2  );
  sLattice.template defineField<A_SOURCE>( superGeometry, 1, zero2 );
  sLattice.template defineField<A_SOURCE>( superGeometry, 2, zero2 );
  sLattice.template defineField<B_FIELD> ( superGeometry, 1, zero2 );
  sLattice.template defineField<B_FIELD> ( superGeometry, 2, zero2 );
  sLattice.template defineField<B_SOURCE>( superGeometry, 1, zero2 );
  sLattice.template defineField<B_SOURCE>( superGeometry, 2, zero2 );

  sLattice.initialize();
}


///////////////////////////////////////////////////////////////////////////////////////////////////
// Plots the results
void getResults( SuperLattice<T, DESCRIPTOR>& sLattice,
                 UnitConverter<T,DESCRIPTOR>& converter,
                 SuperGeometry<T,2>& superGeometry, util::Timer<T>& timer )
{
  OstreamManager clout( std::cout,"getResults" );

  static Gnuplot<T> gplot("concentrations");

  AnalyticalSol<2,T,int> sol (converter.getPhysTime(iT), t0, a, b);
  SuperLatticeExternal2D<T,DESCRIPTOR,A_FIELD> slA ( sLattice, iT );
  SuperLatticeExternal2D<T,DESCRIPTOR,B_FIELD> slB ( sLattice, iT );
  AnalyticalFfromSuperF2D<T> A( slA, true, 1 );
  AnalyticalFfromSuperF2D<T> B( slB, true, 1 );

  int point[] { (int)(nx/2)+1, (int)(ny/2)+1 };
  T pointP[] { point[0]*deltaX, point[1]*deltaX };
  T analytical[] { T(), T() };
  T numericalA[] { T() };
  T numericalB[] { T() };

  sol(analytical,point);
  A(numericalA,pointP);
  B(numericalB,pointP);

  if (iT % tGnu ==0) {
    gplot.setData(converter.getPhysTime(iT),
              {analytical[0], numericalA[0], analytical[1], numericalB[0]},
              {"A analytical", "A numerical", "B analytical", "B numerical"},
              "bottom right", {'l', 'p', 'l', 'p'});
    gplot.writePNG();
  }

  if (iT % tVtm == 0) {
    timer.update( iT );
    timer.printStep();
    sLattice.getStatistics().print(iT, converter.getPhysTime( iT ));
    clout << "A:  analytical=" << analytical[0] << "  numerical=" << numericalA[0]<< "  error=" << util::abs(analytical[0]-numericalA[0])/numericalA[0] << std::endl
          << "B:  analytical=" << analytical[1] << "  numerical=" << numericalB[0]<< "  error=" << util::abs(analytical[1]-numericalB[0])/numericalB[0] << std::endl
          << std::endl;
  }

  if (iT == tMax) {
    gplot.writePDF();
  }

}


///////////////////////////////////////////////////////////////////////////////////////////////////
int main(int argc, char* argv[])
{
  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );

  OstreamManager clout( std::cout,"main" );

  UnitConverter<T,DESCRIPTOR> converter (
    (T) deltaX,    // physDeltaX
    (T) deltaT,    // physDeltaT
    (T) nx*deltaX, // charPhysLength
    (T) 0.,        // charPhysVelocity
    (T) 0.1,       // physViscosity
    (T) 1.         // physDensity
  );
  clout << "---------- Input data: ------------" << std::endl
        << "nx         = " << nx << std::endl
        << "ny         = " << ny << std::endl
        << "tMax       = " << tMax << std::endl
        << "tVtm       = " << tVtm << std::endl
        << "deltaX     = " << deltaX << std::endl
        << "deltaT     = " << deltaT << std::endl
        << "a          = " << a << std::endl
        << "b          = " << b << std::endl
        << "t0         = " << t0 << std::endl
        << "------------------------------------" << std::endl;
  converter.print();

  /// === 2rd Step: Prepare Geometry ===
  std::vector<T> origin { 0., 0. };
  std::vector<T> extend { (nx-1)*deltaX, (ny-1)*deltaX };
  IndicatorCuboid2D<T> cuboid(extend, origin);

  /// Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif
  CuboidGeometry2D<T> cuboidGeometry(cuboid, deltaX, noOfCuboids);
  cuboidGeometry.setPeriodicity(false, false);
  cuboidGeometry.print();

  /// Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

  /// Instantiation of a superGeometry
  SuperGeometry<T,2> superGeometry(cuboidGeometry, loadBalancer, 2);

  prepareGeometry(superGeometry);

  /// === 3rd Step: Prepare Lattice ===
  SuperLattice<T,DESCRIPTOR> sLattice( superGeometry );
  prepareLattice( sLattice, converter, superGeometry );

  // === 4th Step: Definition of Initial and Boundary Conditions ===
  setBoundaryValues( sLattice, converter, superGeometry );

  // === 5th Step: Main Loop with Timer ===
  clout << "starting simulation..." << std::endl;
  util::Timer<T> timer( tMax, superGeometry.getStatistics().getNvoxel() );
  timer.start();

  for ( iT = 0; iT <= tMax; ++iT ) {
    // === 6th Step: Collide and Stream Execution ===
    sLattice.collideAndStream();

    // === 8th Step: Computation and Output of the Results ===
    getResults( sLattice, converter, superGeometry, timer );
  }

  timer.stop();
  timer.printSummary();
}
