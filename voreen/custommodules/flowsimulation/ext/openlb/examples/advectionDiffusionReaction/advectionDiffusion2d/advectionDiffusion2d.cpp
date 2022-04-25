
/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2020 Louis Kronberg, Stephan Simonis
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

/*  adePeriodic2d:
 *  The solution to a linear, scalar, two-dimensional advection-diffusion
 *  equation is approximated.
 *  The numerical setup and the analytical solution are taken from
 *  [Simonis, S., Frank, M., and Krause, M. J. 2020. Phil. Trans. R. Soc. A378:
 *  20190400. DOI: 10.1098/rsta.2019.0400]
 *  Error norms are calculated for three subsequent resolutions and stored
 *  in the respective /tmp folders. A python script is provided to calculate
 *  the experimental order of convergence towards the analytical solution.
 */


#include "olb2D.h"
#include "olb2D.hh"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

typedef double T;

typedef D2Q5<VELOCITY> TDESCRIPTOR;

/// uncomment for console output of other norms
// #define allNorms

const int runs = 3;                // # simulations with increasing resolution

const int N0 = 50;                 // initial # discrete points per dimension
const int statIter0 = N0;          // initial # lattice output timesteps
const T mue0 = 0.05;                // physical diffusivity
const T peclet0 = 100.;            // Peclet number (Pe = u*L/mue)
const T physLength0 = 2.;          // physical domain length in each dimension
const int maxN = util::pow(2,runs-1)*N0; // maximum resolution


template <typename T>
class AdePhysTemp2D : public AnalyticalF2D<T,T> {

protected:
  T t;
  T mue;
  T uMag;
  T res;
public:
  AdePhysTemp2D(T time, AdeUnitConverter<T, TDESCRIPTOR> converter) : AnalyticalF2D<T,T>(1),
    t(time),
    mue(converter.getPhysDiffusivity()),
    uMag(converter.getCharPhysVelocity()),
    res(converter.getResolution())

  {};

  bool operator()(T output[], const T input[]) override
  {
    T x = input[0];
    T y = input[1];
    T gf = res/(res+1.);
    T uMagX = uMag;
    T uMagY = uMag;

    // initial condition (2D)
    output[0] = util::sin(M_PI*(x-uMagX*t)*gf);
    output[0] *= util::sin(M_PI*(y-uMagY*t)*gf);
    output[0] *= util::exp(-2*mue*util::pow(M_PI,2)*t*util::pow(gf,2));
    output[0] += 1.;

    return true;
  };
};


void prepareGeometry(SuperGeometry<T,2>& superGeometry,
                     IndicatorF2D<T>& indicator )
{
  OstreamManager clout(std::cout,"prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0, 1); // , indicator);
  superGeometry.communicate();

  superGeometry.clean();
  /// Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}


void prepareLattice(  SuperLattice<T, TDESCRIPTOR>& ADlattice,
                      SuperGeometry<T,2>& superGeometry,
                      AdeUnitConverter<T, TDESCRIPTOR> converter )
{
  OstreamManager clout(std::cout,"prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  // buffer layer
  ADlattice.defineDynamics<NoDynamics>(superGeometry.getMaterialIndicator({0}));
  // bulk
  ADlattice.defineDynamics<AdvectionDiffusionBGKdynamics>(superGeometry.getMaterialIndicator({1}));

  AnalyticalConst2D<T,T> u0( converter.getCharLatticeVelocity(), converter.getCharLatticeVelocity() );

  AdePhysTemp2D<T> Tinit( 0.0, converter );

  auto bulkIndicator = superGeometry.getMaterialIndicator({0,1});

  ADlattice.defineField<descriptors::VELOCITY>( bulkIndicator, u0 );
  ADlattice.defineRho( bulkIndicator, Tinit );
  ADlattice.iniEquilibrium( bulkIndicator, Tinit, u0 );

  ADlattice.setParameter<descriptors::OMEGA>( converter.getLatticeAdeRelaxationFrequency() );

  /// Make the lattice ready for simulation
  ADlattice.initialize();
  clout << "Prepare Lattice ... OK" << std::endl;
}

// compute absolute and relative errors with different norms
T error( SuperLattice<T, TDESCRIPTOR>& ADlattice,
         SuperGeometry<T,2>& superGeometry,
         int iT,
         AdeUnitConverter<T, TDESCRIPTOR> converter )
{
  OstreamManager clout(std::cout,"error");
  Gnuplot<T> plt("UNUSED_FILE");
  CSV<T> csvWriter("UNUSED_CSV_FILE");

  T result[2] = {T(), T()};
  int tmp[] = {int()};

  SuperLatticeDensity2D<T, TDESCRIPTOR> temperature(ADlattice);
  AdePhysTemp2D<T> temperatureSol( converter.getPhysTime(iT), converter);

  auto indicatorF = superGeometry.getMaterialIndicator({1});

  T l2rel_error;

  SuperRelativeErrorL2Norm2D<T> relTemperatureErrorL2Norm(temperature, temperatureSol, indicatorF);
  SuperAbsoluteErrorL2Norm2D<T> absTemperatureErrorL2Norm(temperature, temperatureSol, indicatorF);

  /// Console output of other error norms
#ifdef allNorms
  SuperAbsoluteErrorL1Norm2D<T> absTemperatureErrorL1Norm(temperature, temperatureSol, indicatorF);
  absTemperatureErrorL1Norm(result, tmp);
  clout << "temperature-L1-absolute-error=" << result[0] << std::endl;
  csvWriter.writeDataFile(iT, result[0], "l1-abs-error");

  SuperRelativeErrorL1Norm2D<T> relTemperatureErrorL1Norm(temperature, temperatureSol, indicatorF);
  relTemperatureErrorL1Norm(result, tmp);
  clout << "temperature-L1-relative-error=" << result[0] << std::endl;
  csvWriter.writeDataFile(iT, result[0],  "l1-rel-error");

  SuperAbsoluteErrorLinfNorm2D<T> absTemperatureErrorLinfNorm(temperature, temperatureSol, indicatorF);
  absTemperatureErrorLinfNorm(result, tmp);
  clout << "temperature-Linf-absolute-error=" << result[0] << std::endl;
  csvWriter.writeDataFile(iT, result[0],  "linf-abs-error");

  SuperRelativeErrorLinfNorm2D<T> relTemperatureErrorLinfNorm(temperature, temperatureSol, indicatorF);
  relTemperatureErrorLinfNorm(result, tmp);
  clout << "temperature-Linf-relative-error=" << result[0] << std::endl;
  csvWriter.writeDataFile(iT, result[0], "linf-rel-error");
#endif

  absTemperatureErrorL2Norm(result, tmp);
  clout << "temperature-L2-absolute-error=" << result[0] << std::endl;
  csvWriter.writeDataFile(iT, result[0],  "l2-abs-error");

  relTemperatureErrorL2Norm(result, tmp);
  clout << "temperature-L2-relative-error=" << result[0] << std::endl;
  csvWriter.writeDataFile(iT, result[0], "l2-rel-error");

  l2rel_error = result[0];

  return l2rel_error;
}

// will always return ZERO for iT == 0
T getResults( SuperLattice<T, TDESCRIPTOR>& ADlattice,
              int iT,
              int statIter,
              AdeUnitConverter<T, TDESCRIPTOR> converter,
              SuperGeometry<T,2>& superGeometry,
              util::Timer<T>& timer)
{
  OstreamManager clout(std::cout,"getResults");
  SuperVTMwriter2D<T> vtkWriter("advectionDiffusion2d");

  AdePhysTemp2D<T> temperatureSol( converter.getPhysTime(iT), converter );

  T avg = 0;

  if (iT == 0) {
    /// Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry2D<T, TDESCRIPTOR> geometry(ADlattice, superGeometry);
    SuperLatticeCuboid2D<T, TDESCRIPTOR> cuboid(ADlattice);
    SuperLatticeRank2D<T, TDESCRIPTOR> rank(ADlattice);
    vtkWriter.write(geometry);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);

    vtkWriter.createMasterFile();

    SuperLatticeDensity2D<T, TDESCRIPTOR> temperature(ADlattice);
    SuperLatticeFfromAnalyticalF2D<T, TDESCRIPTOR> solution(temperatureSol, ADlattice);

    vtkWriter.addFunctor( temperature );
    vtkWriter.addFunctor( solution );
    vtkWriter.write( iT );

    /// Creates contour plots as jpeg files for highest resolution run
    if (converter.getResolution() == maxN) {
      heatmap::plotParam<T> plotParam_temp0;
      plotParam_temp0.minValue = 0;
      plotParam_temp0.maxValue = 2;
      plotParam_temp0.contourlevel = 10;
      BlockReduction2D2D<T> planeReduction0( temperature, converter.getResolution(), BlockDataSyncMode::ReduceOnly );
      heatmap::write( planeReduction0, iT, plotParam_temp0 );
    }

    timer.update(iT);
    timer.printStep();
    ADlattice.getStatistics().print( iT, converter.getPhysTime(iT) );
  }
  else if (iT % statIter == 0 ) {
    /// Writes the VTK files
    SuperLatticeDensity2D<T, TDESCRIPTOR> temperature(ADlattice);
    SuperLatticeFfromAnalyticalF2D<T, TDESCRIPTOR> solution(temperatureSol, ADlattice);

    vtkWriter.addFunctor( temperature );
    vtkWriter.addFunctor( solution );
    vtkWriter.write( iT );

    /// Creates contour plots as jpeg files for highest resolution run
    if (converter.getResolution() == maxN) {
      heatmap::plotParam<T> plotParam_temp;
      plotParam_temp.minValue = 0;
      plotParam_temp.maxValue = 2;
      plotParam_temp.contourlevel = 10;
      BlockReduction2D2D<T> planeReduction( temperature, converter.getResolution(), BlockDataSyncMode::ReduceOnly );
      heatmap::write( planeReduction, iT, plotParam_temp );
    }

    /// ADLattice statistics console output
    timer.update(iT);
    timer.printStep();
    ADlattice.getStatistics().print(iT, converter.getPhysTime( iT ));

    // compute relative error
    avg = error(ADlattice, superGeometry, iT, converter);
    clout << "Relative L2-error norm: "  << avg << std::endl;
  }

  return avg;
}


void simulate(int N, int statIter, T mue, T peclet, T physLength)
{

  OstreamManager clout(std::cout,"simulate");
  clout << "Executing the simulation with N=" << std::to_string(N) << std::endl;

  AdeUnitConverter<T,TDESCRIPTOR> converter(
    physLength/N,            // physDeltaX
    util::pow(physLength/N, 2),    // physDeltaT (diffusive scaling)
    physLength,              // charPhysLength
    peclet*mue/physLength,   // charPhysVelocity from Peclet
    mue,                     // physDiffusivity
    1                        // physDensity,
  );

  converter.print();

  // switch outdirectory only if there are multiple simulation runs
  if (runs > 1) {
    singleton::directories().setOutputDir("./tmp/N" + std::to_string(N) + "/");
  }

  /// === 2nd Step: Prepare Geometry ===
  std::vector<T> extend(2,T());
  std::vector<T> origin(2,T());

  extend[0] = converter.getCharPhysLength();
  extend[1] = converter.getCharPhysLength();

  origin[0] = - converter.getCharPhysLength()/2;
  origin[1] = - converter.getCharPhysLength()/2;

  IndicatorCuboid2D<T> cuboid(extend, origin);

  /// Instantiation of an empty cuboidGeometry
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif
  CuboidGeometry2D<T> cuboidGeometry(cuboid, converter.getPhysDeltaX(), noOfCuboids);
  cuboidGeometry.setPeriodicity(true, true);

  /// Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

  /// Instantiation of a superGeometry
  SuperGeometry<T,2> superGeometry(cuboidGeometry, loadBalancer, 2);
  prepareGeometry(superGeometry, cuboid);

  /// === 3rd Step: Prepare Lattice ===
  SuperLattice<T, TDESCRIPTOR> ADlattice(superGeometry);

  prepareLattice(ADlattice, superGeometry, converter);

  /// === 4th Step: Main Loop with Timer ===
  util::Timer<T> timer( superGeometry.getStatistics().getNvoxel() );
  timer.start();

  T pulseDiffBound = 1e-1;
  int timeCount = util::ceil( -1./converter.getPhysDeltaT() * util::log(pulseDiffBound)/(2.*converter.getPhysDiffusivity()*util::pow(M_PI,2)) /statIter );
  int iTmax = timeCount * statIter;

  int iT;
  T simulationAverage = .0;
  for (iT = 0; iT < iTmax; ++iT) {
    simulationAverage += getResults(ADlattice, iT, statIter, converter, superGeometry, timer);

    /// === 6th Step: Collide and Stream Execution ===
    ADlattice.collideAndStream();
  }

  simulationAverage /= timeCount;

  // this outputs into ./tmp/gnuplotData/data/averageSimL2RelErr
  singleton::directories().setOutputDir("./tmp/");
  Gnuplot<T> plt("UNUSED");
  CSV<T> csvWriter("UNUSED_CSV");
  csvWriter.writeDataFile(N, simulationAverage, "averageSimL2RelErr");

  clout << "Simulation Average Relative L2 Error: " << simulationAverage << std::endl;

  timer.stop();
  timer.printSummary();
}


int main(int argc, char *argv[])
{
  OstreamManager clout(std::cout,"main");
  olbInit(&argc, &argv);

  singleton::directories().setOutputDir("./tmp/");

  for (int i = 0; i < runs; ++i) {
    simulate( util::pow(2,i) * N0,
              util::pow(4,i) * statIter0,
              mue0,
              peclet0,
              physLength0 );
  }
}
