/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2008 Orestis Malaspinas
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

/* rayleighBenard2d.cpp:
 * Rayleigh-Benard convection rolls in 2D, simulated with
 * the thermal LB model by Z. Guo e.a., between a hot plate at
 * the bottom and a cold plate at the top.
 */


#include "olb2D.h"
#include "olb2D.hh"   // use only generic version!
#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>


using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;
using namespace std;

typedef double T;

#define NSDESCRIPTOR ForcedD2Q9Descriptor
#define TDESCRIPTOR AdvectionDiffusionD2Q5Descriptor

// const int maxIter  = 1000000;
// const int saveIter = 5000;

// Parameters for the simulation setup
const T lx  = 1.;      // length of the channel
const T ly  = 1.;      // height of the channel
int     N = 20;        // resolution of the model
T       tau = 1.;      // relaxation time
const T Re = 5.;       // Reynolds number
const T Ra = 100.;     // Rayleigh number
const T Pr = 0.71;     // Prandtl number
const T maxPhysT = 1e4; // max. simulation time in s, SI unit
const T epsilon = 1.e-5; // precision of the convergence (residuum)

const T Tcold = 273.15;
const T Thot = 274.15;

// analytical solution from point light source in infinte domain
// appliacation from R3 to R1.
// effective for x in R3, only the distance to (0,0) is needed.
// documentation e.g. Biomedical Optics, Lihong V. Wang Hsin-I Wu
template <typename T, typename S>
class AnalyticalVelocityPorousPlate2D : public AnalyticalF2D<T, S> {
private:
  T _Re;
  T _u0;
  T _v0;
  T _ly;
public:
  AnalyticalVelocityPorousPlate2D(T Re, T u0, T v0, T ly) : AnalyticalF2D<T, S>(2),
   _Re(Re), _u0(u0), _v0(v0), _ly(ly)
  {
    this->getName() = "AnalyticalVelocityPorousPlate2D";
  };

  bool operator()(T output[2], const S x[2])
  {
    output[0] = _u0*((exp(_Re* x[1] / _ly) - 1) / (exp(_Re) - 1));
    output[1] = _v0;
    return true;
  };
};

template <typename T, typename S>
class AnalyticalTemperaturePorousPlate2D : public AnalyticalF2D<T, S> {
private:
  T _Re;
  T _Pr;
  T _ly;
  T _T0;
  T _deltaT;
public:
  AnalyticalTemperaturePorousPlate2D(T Re, T Pr, T ly, T T0, T deltaT) : AnalyticalF2D<T, S>(1),
    _Re(Re), _Pr(Pr), _ly(ly), _T0(T0), _deltaT(deltaT)
  {
    this->getName() = "AnalyticalTemperaturePorousPlate2D";
  };

  bool operator()(T output[1], const S x[2])
  {
    output[0] = _T0 + _deltaT*((exp(_Pr*_Re*x[1] / _ly) - 1) / (exp(_Pr*_Re) - 1));
    return true;
  };
};

void error( SuperGeometry2D<T>& superGeometry,
            SuperLattice2D<T, NSDESCRIPTOR>& NSlattice,
            SuperLattice2D<T, TDESCRIPTOR>& ADlattice,
            ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter,
            T Re )
{
  OstreamManager clout( std::cout, "error" );

  int input[1] = {0};
  T normAnaSol[1], absErr[1], relErr[1];

  T u_Re = Re * converter.getPhysViscosity() / converter.getCharPhysLength();
  AnalyticalVelocityPorousPlate2D<T,T> uSol(Re, converter.getCharPhysVelocity(), u_Re, converter.getCharPhysLength());
  SuperLatticePhysVelocity2D<T,NSDESCRIPTOR> u1(NSlattice,converter);
  SuperLatticeFfromAnalyticalF2D<T,NSDESCRIPTOR> uSolLattice(uSol,NSlattice);

  SuperL2Norm2D<T> uL2Norm( uSolLattice - u1, superGeometry, 1 );
  SuperL2Norm2D<T> uSolL2Norm( uSolLattice, superGeometry, 1 );
  uL2Norm( absErr, input );
  uSolL2Norm( normAnaSol, input );
  relErr[0] = absErr[0] / normAnaSol[0];
  clout << "velocity-L2-error(abs)=" << absErr[0] << "; velocity-L2-error(rel)=" << relErr[0] << std::endl;

  int inputT[1] = {0};
  T normAnaSolT[1], absErrT[1], relErrT[1];

  AnalyticalTemperaturePorousPlate2D<T,T> TSol(Re, Pr, converter.getCharPhysLength(), converter.getCharPhysLowTemperature(), converter.getCharPhysTemperatureDifference());
  SuperLatticePhysTemperature2D<T,NSDESCRIPTOR,TDESCRIPTOR> T1(ADlattice,converter);
  SuperLatticeFfromAnalyticalF2D<T,TDESCRIPTOR> TSolLattice(TSol,ADlattice);

  SuperL2Norm2D<T> TL2Norm( TSolLattice - T1, superGeometry, 1 );
  SuperL2Norm2D<T> TSolL2Norm( TSolLattice, superGeometry, 1 );
  TL2Norm( absErrT, inputT );
  TSolL2Norm( normAnaSolT, inputT );
  relErrT[0] = absErrT[0] / normAnaSolT[0];
  clout << "temperature-L2-error(abs)=" << absErrT[0] << "; temperature-L2-error(rel)=" << relErrT[0] << std::endl;

}

/// Stores geometry information in form of material numbers
void prepareGeometry(SuperGeometry2D<T>& superGeometry,
                     ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter)
{

  OstreamManager clout(std::cout,"prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  superGeometry.rename(0,2);
  superGeometry.rename(2,1,0,1);

  std::vector<T> extend( 2, T(0) );
  extend[0] = lx+2*converter.getPhysLength(1);
  extend[1] = converter.getPhysLength(1);
  std::vector<T> origin( 2, T(0) );
  origin[0] = -converter.getPhysLength(1);
  IndicatorCuboid2D<T> bottom(extend, origin);
  /// Set material number for bottom
  superGeometry.rename(2,3,1,bottom);

  /// Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  /// Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

void prepareLattice( ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter,
                     SuperLattice2D<T, NSDESCRIPTOR>& NSlattice,
                     SuperLattice2D<T, TDESCRIPTOR>& ADlattice,
                     Dynamics<T, NSDESCRIPTOR> &bulkDynamics,
                     Dynamics<T, TDESCRIPTOR>& advectionDiffusionBulkDynamics,
                     sOnLatticeBoundaryCondition2D<T,NSDESCRIPTOR>& NSboundaryCondition,
                     sOnLatticeBoundaryCondition2D<T,TDESCRIPTOR>& TboundaryCondition,
                     SuperGeometry2D<T>& superGeometry )
{

  OstreamManager clout(std::cout,"prepareLattice");

  T Tomega  = converter.getLatticeThermalRelaxationFrequency();
  T NSomega = converter.getLatticeRelaxationFrequency();

  /// define lattice Dynamics
  clout << "defining dynamics" << endl;

  ADlattice.defineDynamics(superGeometry, 0, &instances::getNoDynamics<T, TDESCRIPTOR>());
  NSlattice.defineDynamics(superGeometry, 0, &instances::getNoDynamics<T, NSDESCRIPTOR>());

  ADlattice.defineDynamics(superGeometry, 1, &advectionDiffusionBulkDynamics);
  ADlattice.defineDynamics(superGeometry, 2, &advectionDiffusionBulkDynamics);
  ADlattice.defineDynamics(superGeometry, 3, &advectionDiffusionBulkDynamics);
  NSlattice.defineDynamics(superGeometry, 1, &bulkDynamics);
  NSlattice.defineDynamics(superGeometry, 2, &bulkDynamics);
  NSlattice.defineDynamics(superGeometry, 3, &bulkDynamics);


  /// sets boundary
  NSboundaryCondition.addVelocityBoundary(superGeometry, 2, NSomega);
  NSboundaryCondition.addVelocityBoundary(superGeometry, 3, NSomega);
  TboundaryCondition.addTemperatureBoundary(superGeometry, 2, Tomega);
  TboundaryCondition.addTemperatureBoundary(superGeometry, 3, Tomega);
}

void setBoundaryValues(ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter,
                       SuperLattice2D<T, NSDESCRIPTOR>& NSlattice,
                       SuperLattice2D<T, TDESCRIPTOR>& ADlattice,
                       int iT, SuperGeometry2D<T>& superGeometry)
{

  if (iT == 0) {

    // typedef advectionDiffusionLbHelpers<T,TDESCRIPTOR> TlbH;

    /// for each material set the defineRhoU and the Equilibrium
    std::vector<T> zero(2,T());
    AnalyticalConst2D<T,T> u(zero);
    AnalyticalConst2D<T,T> rho(1.);
    AnalyticalConst2D<T,T> force(zero);

    T u_Re = converter.getLatticeVelocity( Re * converter.getPhysViscosity() / converter.getCharPhysLength() );

    AnalyticalConst2D<T,T> u_top(converter.getCharLatticeVelocity(), u_Re);
    AnalyticalConst2D<T,T> u_bot(0.0, u_Re);

    NSlattice.defineRhoU(superGeometry, 1, rho, u);
    NSlattice.iniEquilibrium(superGeometry, 1, rho, u);
    NSlattice.defineExternalField(superGeometry, 1,
                                  NSDESCRIPTOR<T>::ExternalField::forceBeginsAt,
                                  NSDESCRIPTOR<T>::ExternalField::sizeOfForce, force );
    NSlattice.defineRhoU(superGeometry, 2, rho, u_top);
    NSlattice.iniEquilibrium(superGeometry, 2, rho, u_top);
    NSlattice.defineExternalField(superGeometry, 2,
                                  NSDESCRIPTOR<T>::ExternalField::forceBeginsAt,
                                  NSDESCRIPTOR<T>::ExternalField::sizeOfForce, force );
    NSlattice.defineRhoU(superGeometry, 3, rho, u_bot);
    NSlattice.iniEquilibrium(superGeometry, 3, rho, u_bot);
    NSlattice.defineExternalField(superGeometry, 3,
                                  NSDESCRIPTOR<T>::ExternalField::forceBeginsAt,
                                  NSDESCRIPTOR<T>::ExternalField::sizeOfForce, force );

    AnalyticalConst2D<T,T> Cold(converter.getLatticeTemperature(Tcold));
    AnalyticalConst2D<T,T> Hot(converter.getLatticeTemperature(Thot));

    ADlattice.defineRho(superGeometry, 1, Cold);
    ADlattice.iniEquilibrium(superGeometry, 1, Cold, u);
    ADlattice.defineExternalField(superGeometry, 1,
                                  TDESCRIPTOR<T>::ExternalField::velocityBeginsAt,
                                  TDESCRIPTOR<T>::ExternalField::sizeOfVelocity, u );
    ADlattice.defineRho(superGeometry, 3, Cold);
    ADlattice.iniEquilibrium(superGeometry, 3, Cold, u);
    ADlattice.defineExternalField(superGeometry, 3,
                                  TDESCRIPTOR<T>::ExternalField::velocityBeginsAt,
                                  TDESCRIPTOR<T>::ExternalField::sizeOfVelocity, u );
    ADlattice.defineRho(superGeometry, 2, Hot);
    ADlattice.iniEquilibrium(superGeometry, 2, Hot, u);
    ADlattice.defineExternalField(superGeometry, 2,
                                  TDESCRIPTOR<T>::ExternalField::velocityBeginsAt,
                                  TDESCRIPTOR<T>::ExternalField::sizeOfVelocity, u );

    /// Make the lattice ready for simulation
    NSlattice.initialize();
    ADlattice.initialize();
  }
}

void getResults(ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const& converter,
                SuperLattice2D<T, NSDESCRIPTOR>& NSlattice,
                SuperLattice2D<T, TDESCRIPTOR>& ADlattice, int iT,
                SuperGeometry2D<T>& superGeometry,
                Timer<T>& timer,
                bool converged)
{

  OstreamManager clout(std::cout,"getResults");

  SuperVTMwriter2D<T> vtkWriter("thermalPorousPlate2d");
  SuperLatticePhysVelocity2D<T, NSDESCRIPTOR> velocity(NSlattice, converter);
  SuperLatticePhysPressure2D<T, NSDESCRIPTOR> pressure(NSlattice, converter);
  SuperLatticePhysTemperature2D<T, NSDESCRIPTOR, TDESCRIPTOR> temperature(ADlattice, converter);
  vtkWriter.addFunctor( pressure );
  vtkWriter.addFunctor( velocity );
  vtkWriter.addFunctor( temperature );

  const int vtkIter = converter.getLatticeTime(100.);

  if (iT == 0) {
    /// Writes the converter log file
    // writeLogFile(converter,"thermalPorousPlate2d");

    /// Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry2D<T, NSDESCRIPTOR> geometry(NSlattice, superGeometry);
    SuperLatticeCuboid2D<T, NSDESCRIPTOR> cuboid(NSlattice);
    SuperLatticeRank2D<T, NSDESCRIPTOR> rank(NSlattice);
    vtkWriter.write(geometry);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);

    vtkWriter.createMasterFile();
  }

  /// Writes the VTK files
  if (iT%vtkIter == 0 || converged) {
    NSlattice.getStatistics().print(iT,converter.getPhysTime(iT));
    timer.print(iT);
    error(superGeometry, NSlattice, ADlattice, converter, Re);

    vtkWriter.write(iT);

    ///writes Jpeg
    //SuperEuklidNorm2D<T, DESCRIPTOR> normVel(velocity);
    BlockReduction2D2D<T> planeReduction( temperature, 600, BlockDataSyncMode::ReduceOnly );
    // write output of velocity as JPEG
    heatmap::plotParam<T> jpeg_Param;
    jpeg_Param.maxValue = Thot;
    jpeg_Param.minValue = Tcold;
    heatmap::write(planeReduction, iT, jpeg_Param);
  }

}

int main(int argc, char *argv[])
{

  /// === 1st Step: Initialization ===
  OstreamManager clout(std::cout,"main");
  olbInit(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");

  if (argc >= 2) {
    N = atoi(argv[1]);
  }
  if (argc == 3) {
    tau = atof(argv[2]);
  }

  ThermalUnitConverter<T, NSDESCRIPTOR, TDESCRIPTOR> const converter(
    (T) 1.0 / N, // physDeltaX
    (T) 1.0 / N * 1.0 / 1e-3 * (tau - 0.5) / 3 / N, // physDeltaT
    (T) 1.0, // charPhysLength
    (T) sqrt( 9.81 * Ra * 1e-3 * 1e-3 / Pr / 9.81 / (Thot - Tcold) / pow(1.0, 3) * (Thot - Tcold) * 1.0 ), // charPhysVelocity
    (T) 1e-3, // physViscosity
    (T) 1.0, // physDensity
    (T) 0.03, // physThermalConductivity
    (T) Pr * 0.03 / 1e-3 / 1.0, // physSpecificHeatCapacity
    (T) Ra * 1e-3 * 1e-3 / Pr / 9.81 / (Thot - Tcold) / pow(1.0, 3), // physThermalExpansionCoefficient
    (T) Tcold, // charPhysLowTemperature
    (T) Thot // charPhysHighTemperature
  );
  converter.print();

  /// === 2nd Step: Prepare Geometry ===
  std::vector<T> extend(2,T());
  extend[0] = lx;
  extend[1] = ly;
  std::vector<T> origin(2,T());
  IndicatorCuboid2D<T> cuboid(extend, origin);

  /// Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif
  CuboidGeometry2D<T> cuboidGeometry(cuboid, converter.getPhysDeltaX(), noOfCuboids);
  cuboidGeometry.setPeriodicity(true, false);

  /// Instantiation of a loadBalancer
  HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

  /// Instantiation of a superGeometry
  SuperGeometry2D<T> superGeometry(cuboidGeometry, loadBalancer, 2);

  prepareGeometry(superGeometry, converter);

  /// === 3rd Step: Prepare Lattice ===

  SuperLattice2D<T, TDESCRIPTOR> ADlattice(superGeometry);
  SuperLattice2D<T, NSDESCRIPTOR> NSlattice(superGeometry);

  sOnLatticeBoundaryCondition2D<T, NSDESCRIPTOR> NSboundaryCondition(NSlattice);
  createLocalBoundaryCondition2D<T, NSDESCRIPTOR>(NSboundaryCondition);

  sOnLatticeBoundaryCondition2D<T, TDESCRIPTOR> TboundaryCondition(ADlattice);
  createAdvectionDiffusionBoundaryCondition2D<T, TDESCRIPTOR>(TboundaryCondition);

  ForcedBGKdynamics<T, NSDESCRIPTOR> NSbulkDynamics(
    converter.getLatticeRelaxationFrequency(),
    instances::getBulkMomenta<T,NSDESCRIPTOR>());

  AdvectionDiffusionBGKdynamics<T, TDESCRIPTOR> TbulkDynamics (
    converter.getLatticeThermalRelaxationFrequency(),
    instances::getAdvectionDiffusionBulkMomenta<T,TDESCRIPTOR>()
  );

  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//
  // This coupling must be necessarily be put on the Navier-Stokes lattice!!
  // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!//

  std::vector<T> dir{0.0, 1.0};

  T boussinesqForcePrefactor = 9.81 / converter.getConversionFactorVelocity() * converter.getConversionFactorTime() *
    converter.getCharPhysTemperatureDifference() * converter.getPhysThermalExpansionCoefficient();

  NavierStokesAdvectionDiffusionCouplingGenerator2D<T,NSDESCRIPTOR>
    coupling(0, converter.getLatticeLength(lx), 0, converter.getLatticeLength(ly),
        boussinesqForcePrefactor, converter.getLatticeTemperature(Tcold), 1., dir);

  NSlattice.addLatticeCoupling(superGeometry, 1, coupling, ADlattice);
  NSlattice.addLatticeCoupling(superGeometry, 2, coupling, ADlattice);
  NSlattice.addLatticeCoupling(superGeometry, 3, coupling, ADlattice);

  prepareLattice(converter,
                 NSlattice, ADlattice,
                 NSbulkDynamics, TbulkDynamics,
                 NSboundaryCondition, TboundaryCondition, superGeometry );

  /// === 4th Step: Main Loop with Timer ===
  Timer<T> timer(converter.getLatticeTime(maxPhysT), superGeometry.getStatistics().getNvoxel() );
  timer.start();

  util::ValueTracer<T> converge(converter.getLatticeTime(1.0),epsilon);
  for (int iT = 0; iT < converter.getLatticeTime(maxPhysT); ++iT) {

    if (converge.hasConverged()) {
      clout << "Simulation converged." << endl;
      getResults(converter, NSlattice, ADlattice, iT, superGeometry, timer, converge.hasConverged());
      break;
    }

    /// === 5th Step: Definition of Initial and Boundary Conditions ===
    setBoundaryValues(converter, NSlattice, ADlattice, iT, superGeometry);

    /// === 6th Step: Collide and Stream Execution ===
    ADlattice.collideAndStream();
    NSlattice.collideAndStream();

    NSlattice.executeCoupling();

    /// === 7th Step: Computation and Output of the Results ===
    getResults(converter, NSlattice, ADlattice, iT, superGeometry, timer, converge.hasConverged());
    converge.takeValue(ADlattice.getStatistics().getAverageEnergy());
  }

  timer.stop();
  timer.printSummary();
}
