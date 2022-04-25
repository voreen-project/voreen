
/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2021 Stephan Simonis, Davide Dapelo
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

/*  adePeriodic3d:
 *  The solution to a linear, scalar, three-dimensional advection-diffusion
 *  equation is approximated.
 *  The numerical setup and the analytical solution for the smooth initial value
 *  problem (IVP) are extruded to three dimensions from the one in
 *  Simonis et al. Phil. Trans. R. Soc. A (2020) 378:20190400,
 *  DOI: https://doi.org/10.1098/rsta.2019.0400.
 *  The numerical setup and the analytical solution for the unsmooth IVP are
 *  proposed as an initial Dirac comb profile in
 *  Dapelo et al., Journal of Computational Science (2021) 51:101363,
 *  DOI: https://doi.org/10.1016/j.jocs.2021.101363.
 *  Error norms are calculated for three subsequent resolutions and stored
 *  in the respective /tmp folders. A python script is provided to calculate
 *  the experimental order of convergence towards the analytical solution.
 */


#include "olb3D.h"
#include "olb3D.hh"
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>

#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

using T = double;
using TDESCRIPTOR = D3Q7<VELOCITY>;

/// uncomment for console output of other norms
#define allNorms

const int nonsmooth = 1;           //switch for initial condition
const int runs = 3;                // # simulations with increasing resolution
const int N0 = 20;                 // initial # discrete points per dimension
const int statIter0 = 20;          // initial # lattice output timesteps

// Note: the peclet number can also be passed as an argument
T peclet0 = 100.;
const T physVel0 = 2.5;
const T physLength0 = 2.;          // physical domain length in each dimension
const int maxN = pow(2,runs-1)*N0; // maximum resolution

template <typename T>
class AdePhysTemp3D : public AnalyticalF3D<T,T> {

protected:
    T t;
    T mue;
    T uMag;
    T res;
    T deltaT;

public:
    AdePhysTemp3D(T time, AdeUnitConverter <T, TDESCRIPTOR> converter) :
      AnalyticalF3D<T, T>(1),
      t(time),
      mue(converter.getPhysDiffusivity()),
      uMag(converter.getCharPhysVelocity()),
      res(converter.getResolution()),
      deltaT(converter.getPhysDeltaT()) {};

    bool operator()(T output[], const T input[]) override {
        T x = input[0];
        T y = input[1];
        T z = input[2];
        T gf = res/(res+1.); // domain correction
        if (nonsmooth == 1){
          gf = 1.;
        }
        T uMagX = uMag;
        T uMagY = uMag;
        T uMagZ = uMag;

        if (nonsmooth == 1) {
          /// dirac comb (see https://doi.org/10.1016/j.jocs.2021.101363)
          T startingPt = 0.;
          T tShift = deltaT;
          int numSpikes = 10;
          output[0] = 0.;
          for (int i=-numSpikes; i<=numSpikes; ++i){
            output[0] += exp( -pow(gf*(x-uMagX*(t+tShift)-startingPt-i*(res+1)*sqrt(deltaT)),2) / (4*(t+tShift)*mue*pow(gf,2)) );
          }
          T initialConcentr = 1./sqrt(4*M_PI*mue);
          output[0] *= initialConcentr * ( 1 / sqrt((t+tShift)*pow(gf,2)) );
          output[0] += 1.;
        }
        else {
          /// smooth modes (3D extrude of https://doi.org/10.1098/rsta.2019.0400)
          output[0] = sin(M_PI*(x-uMagX*t)*gf);
          output[0] *= sin(M_PI*(y-uMagY*t)*gf);
          output[0] *= sin(M_PI*(z-uMagZ*t)*gf);
          output[0] *= exp(-3*mue*pow(M_PI,2)*t*pow(gf,2));
          output[0] += 1.;
        }
        return true;
    };
};

template<typename T>
class AdePhysTemp3Dinit : public AnalyticalF3D<T, T> {

protected:
    T t;
    T mue;
    T uMag;
    T res;
    T deltaT;

public:
    AdePhysTemp3Dinit(T time, AdeUnitConverter <T, TDESCRIPTOR> converter) :
      AnalyticalF3D<T, T>(1),
      t(time),
      mue(converter.getPhysDiffusivity()),
      uMag(converter.getCharPhysVelocity()),
      res(converter.getResolution()),
      deltaT(converter.getPhysDeltaT()) {};

    bool operator()(T output[], const T input[]) override {

        if (nonsmooth==1) {
          T x = input[0];
          T gf = 1.;
          /// naive dirac comb (single spike)
          if (std::fabs(x) <= 1. / res) {
              output[0] = 1 / sqrt(deltaT * pow(gf, 2));
              output[0] += 1.;
          } else {
              output[0] = 1.;
          }
        }
        else {
          T x = input[0];
          T y = input[1];
          T z = input[2];
          T gf = res/(res+1.); // domain correction
          /// smooth initial condition (3D)
          T uMagX = uMag;
          T uMagY = uMag;
          T uMagZ = uMag;
          output[0] = sin(M_PI*(x-uMagX*t)*gf);
          output[0] *= sin(M_PI*(y-uMagY*t)*gf);
          output[0] *= sin(M_PI*(z-uMagZ*t)*gf);
          output[0] *= exp(-3*mue*pow(M_PI,2)*t*pow(gf,2));
          output[0] += 1.;
        }

        return true;
    };
};

// geometry preparation for cube
void prepareGeometry(SuperGeometry <T,3> &superGeometry,
                     IndicatorF3D <T> &indicator) {
    OstreamManager clout(std::cout, "prepareGeometry");
    clout << "Prepare Geometry ..." << std::endl;

    superGeometry.rename(0, 1);
    superGeometry.communicate();
    superGeometry.print();

    clout << "Prepare Geometry ... OK" << std::endl;
}


void prepareLattice(SuperLattice <T, TDESCRIPTOR> &ADlattice,
                    SuperGeometry <T,3> &superGeometry,
                    AdeUnitConverter <T, TDESCRIPTOR> converter) {
    OstreamManager clout(std::cout, "prepareLattice");
    clout << "Prepare Lattice ..." << std::endl;

    // buffer layer
    ADlattice.defineDynamics<NoDynamics>(superGeometry.getMaterialIndicator({0}));
    // bulk
    ADlattice.defineDynamics<AdvectionDiffusionBGKdynamics>(superGeometry.getMaterialIndicator({1}));
    // Initial advection velocity
    AnalyticalConst3D<T,T> u0( converter.getCharLatticeVelocity(),
                               converter.getCharLatticeVelocity(),
                               converter.getCharLatticeVelocity() );
    if(nonsmooth==1){
        AnalyticalConst3D <T, T> u0(converter.getCharLatticeVelocity(), 0.0, 0.0);
    }

    AdePhysTemp3Dinit<T> Tinit( 0.0, converter );

    auto bulkIndicator = superGeometry.getMaterialIndicator({0, 1});

    ADlattice.defineField<descriptors::VELOCITY>(bulkIndicator, u0);
    ADlattice.defineRho(bulkIndicator, Tinit);
    ADlattice.iniEquilibrium(bulkIndicator, Tinit, u0);

    ADlattice.setParameter<descriptors::OMEGA>( converter.getLatticeAdeRelaxationFrequency() );
    ADlattice.initialize();
    clout << "Prepare Lattice ... OK" << std::endl;
}

// compute absolute and relative errors with different norms
T error(SuperLattice<T, TDESCRIPTOR>& ADlattice,
        SuperGeometry<T,3>& superGeometry,
        int iT,
        AdeUnitConverter <T, TDESCRIPTOR> converter) {
    OstreamManager clout(std::cout, "error");
    Gnuplot <T> plt("UNUSED_FILE");
    CSV<T> csvWriter("UNUSED_CSV_FILE");

    T result[3] = {T(), T(), T()};
    int tmp[] = {int()};

    SuperLatticeDensity3D <T, TDESCRIPTOR> temperature(ADlattice);
    AdePhysTemp3D<T> temperatureSol(converter.getPhysTime(iT), converter);

    auto indicatorF = superGeometry.getMaterialIndicator({1});

    T l2rel_error;

    SuperRelativeErrorL2Norm3D <T> relTemperatureErrorL2Norm(temperature, temperatureSol, indicatorF);
    SuperAbsoluteErrorL2Norm3D <T> absTemperatureErrorL2Norm(temperature, temperatureSol, indicatorF);

    /// Console output of other error norms
#ifdef allNorms
    SuperAbsoluteErrorL1Norm3D<T> absTemperatureErrorL1Norm(temperature, temperatureSol, indicatorF);
    absTemperatureErrorL1Norm(result, tmp);
    clout << "temperature-L1-absolute-error=" << result[0] << std::endl;
    csvWriter.writeDataFile(iT, result[0], "l1-abs-error");

    SuperRelativeErrorL1Norm3D<T> relTemperatureErrorL1Norm(temperature, temperatureSol, indicatorF);
    relTemperatureErrorL1Norm(result, tmp);
    clout << "temperature-L1-relative-error=" << result[0] << std::endl;
    csvWriter.writeDataFile(iT, result[0],  "l1-rel-error");

    SuperAbsoluteErrorLinfNorm3D<T> absTemperatureErrorLinfNorm(temperature, temperatureSol, indicatorF);
    absTemperatureErrorLinfNorm(result, tmp);
    clout << "temperature-Linf-absolute-error=" << result[0] << std::endl;
    csvWriter.writeDataFile(iT, result[0],  "linf-abs-error");

    SuperRelativeErrorLinfNorm3D<T> relTemperatureErrorLinfNorm(temperature, temperatureSol, indicatorF);
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

T getResults(SuperLattice<T, TDESCRIPTOR> &ADlattice,
             int iT,
             int statIter,
             AdeUnitConverter <T, TDESCRIPTOR> converter,
             SuperGeometry<T,3> &superGeometry,
             util::Timer <T> &timer) {
    OstreamManager clout(std::cout, "getResults");
    SuperVTMwriter3D <T> vtkWriter("advectionDiffusion3d");

    AdePhysTemp3D<T> temperatureSol(converter.getPhysTime(iT), converter);

    T avg = 0;

    if (iT == 0) {
      /// Writes the geometry, cuboid no. and rank no. as vti file for visualization
      SuperLatticeGeometry3D <T, TDESCRIPTOR> geometry(ADlattice, superGeometry);
      SuperLatticeCuboid3D <T, TDESCRIPTOR> cuboid(ADlattice);
      SuperLatticeRank3D <T, TDESCRIPTOR> rank(ADlattice);
      vtkWriter.write(geometry);
      vtkWriter.write(cuboid);
      vtkWriter.write(rank);

      vtkWriter.createMasterFile();

      SuperLatticeDensity3D <T, TDESCRIPTOR> temperature(ADlattice);
      SuperLatticeFfromAnalyticalF3D <T, TDESCRIPTOR> solution(temperatureSol, ADlattice);

      vtkWriter.addFunctor(temperature);
      vtkWriter.addFunctor(solution);
      vtkWriter.write(iT);

      timer.update(iT);
      timer.printStep();
      ADlattice.getStatistics().print(iT, converter.getPhysTime(iT));
    } else if (iT % statIter == 0 && iT != 0) {
      ADlattice.setProcessingContext(ProcessingContext::Evaluation);
      /// Writes the VTK files
      SuperLatticeDensity3D <T, TDESCRIPTOR> temperature(ADlattice);
      SuperLatticeFfromAnalyticalF3D <T, TDESCRIPTOR> solution(temperatureSol, ADlattice);

      vtkWriter.addFunctor(temperature);
      vtkWriter.addFunctor(solution);
      vtkWriter.write(iT);

      /// ADLattice statistics console output
      timer.update(iT);
      timer.printStep();
      ADlattice.getStatistics().print(iT, converter.getPhysTime(iT));

      // compute relative error
      avg = error(ADlattice, superGeometry, iT, converter);
      clout << "Relative L2-error norm: " << avg << std::endl;
    }

    return avg;
}


void simulate(int N, int statIter, T physVel, T peclet, T physLength) {

    OstreamManager clout(std::cout, "simulate");
    clout << "Executing the simulation with N=" << std::to_string(N) << std::endl;

    AdeUnitConverter<T,TDESCRIPTOR> converter(
      physLength/N,                   // physDeltaX
      util::pow(physLength/N, 2),     // physDeltaT (diffusive scaling)
      physLength,                     // charPhysLength
      physVel,                        // charPhysVelocity
      (physLength*physVel)/peclet,    // physDiffusivity from peclet
      1                               // physDensity,
    );

    converter.print();

    // switch outdirectory if there multiple simulation runs
    // store results in subfolder corresponding to peclet number "p_XX"
    if (runs > 1) {
        singleton::directories().setOutputDir("./tmp/p_" + std::to_string((int)peclet) + "/N" + std::to_string(N) + "/");
    }

    /// === 2nd Step: Prepare Geometry ===
    std::vector <T> extend(3, T());
    std::vector <T> origin(3, T());

    T eps = std::numeric_limits<T>::epsilon();

    extend[0] = converter.getCharPhysLength() + 2*eps;
    extend[1] = converter.getCharPhysLength() + 2*eps;
    extend[2] = converter.getCharPhysLength() + 2*eps;

    origin[0] = -converter.getCharPhysLength() / 2.0 - eps;
    origin[1] = -converter.getCharPhysLength() / 2.0 - eps;
    origin[2] = -converter.getCharPhysLength() / 2.0 - eps;

    IndicatorCuboid3D <T> cuboid(extend, origin);

    /// Instantiation of an empty cuboidGeometry
#ifdef PARALLEL_MODE_MPI
    const int noOfCuboids = singleton::mpi().getSize();
#else
    const int noOfCuboids = 1;
#endif
    CuboidGeometry3D <T> cuboidGeometry(cuboid, converter.getPhysDeltaX(), noOfCuboids);
    cuboidGeometry.setPeriodicity(true, true, true);

    /// Instantiation of a loadBalancer
    HeuristicLoadBalancer <T> loadBalancer(cuboidGeometry);

    /// Instantiation of a superGeometry
    SuperGeometry<T,3> superGeometry(cuboidGeometry, loadBalancer, 2);
    prepareGeometry(superGeometry, cuboid);

    /// === 3rd Step: Prepare Lattice ===
    SuperLattice<T, TDESCRIPTOR> ADlattice(superGeometry);

    prepareLattice(ADlattice, superGeometry, converter);

    T physTmax = 1.52;
    T pulseDiffBound = 1e-1;
    int iTmax;
    int timeCount;

    if (nonsmooth==1){
      iTmax = converter.getLatticeTime(physTmax);
      timeCount = ceil(iTmax / statIter);
    } else {
       timeCount = ceil( -1./converter.getPhysDeltaT() * log(pulseDiffBound)/(3.*converter.getPhysDiffusivity()*util::pow(M_PI,2)) /statIter );
       iTmax = timeCount * statIter;
    }

    /// === 4th Step: Main Loop with Timer ===
    util::Timer<T> timer(iTmax, superGeometry.getStatistics().getNvoxel());
    timer.start();

    int iT;
    T simulationAverage = .0;
    for (iT = 0; iT < iTmax; ++iT) {
        simulationAverage += getResults(ADlattice, iT, statIter, converter, superGeometry, timer);

        /// === 6th Step: Collide and Stream Execution ===
        ADlattice.collideAndStream();
    }

    simulationAverage /= timeCount;

    // this outputs into ./tmp/gnuplotData/data/averageSimL2RelErr
    singleton::directories().setOutputDir("./tmp/p_" + std::to_string((int)peclet) + "/");
    Gnuplot <T> plt("UNUSED");
    CSV<T> csvWriter("UNUSED_CSV");
    csvWriter.writeDataFile(N, simulationAverage, "averageSimL2RelErr");

    clout << "Simulation Average Relative L2 Error: " << simulationAverage << std::endl;

    timer.stop();
    timer.printSummary();
}


int main(int argc, char *argv[]) {
    OstreamManager clout(std::cout, "main");
    olbInit(&argc, &argv);

    // Get peclet number passed as argument
    if (argc > 1) peclet0 = atof(argv[1]);
    singleton::directories().setOutputDir("./tmp/p_" + std::to_string((int)peclet0) + "/");

    for (int i = 0; i < runs; ++i) {
        clout << "<------- Simulate with -------> " << std::endl;
        clout << "Resolution: " << util::pow(2,i)*N0 << std::endl;
        clout << "Grid Peclet number: " << peclet0/(util::pow(2,i)*N0) << std::endl;
        clout << "Courant number: " << physLength0*physVel0/(util::pow(2,i)*N0) << std::endl;
        simulate(pow(2, i) * N0,
                 pow(4, i) * statIter0,
                 physVel0,
                 peclet0,
                 physLength0);
    }

}
