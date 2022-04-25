/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2021 Johanna Moedl, Julius Jessberger
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

/* advectionDiffusionReaction2d.cpp:
 * This example illustrates a steady-state chemical reaction in a plug flow reactor.
 * One can choose two types of reaction, A->C and A<->C.
 * The concentration and analytical solution along the centerline of the rectangle domain
 * is given in tmp/N<resolution>/gnuplotData as well as the error plot for the concentration
 * along the centerline.
 * There are runs=3 simulation runs executed and the average L2 error over the centerline
 * is computed for each resolution. An EOC plot can be found in tmp/gnuplotData.
 */

#include "olb2D.h"
#include "olb2D.hh"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

typedef double T;
typedef D2Q5<VELOCITY,SOURCE> TDESCRIPTOR;

const int runs = 3;              // # simulations with increasing resolution

// Parameters for the simulation setup
const T lx  = 10.;               // length of the channel
const int N0 = 200;              // resolution of the model
const int statIter0 = N0*N0;     // initial # lattice output timesteps
const T mue = 0.1;               // physical diffusivity
const T peclet = 100.;           // Peclet number (Pe = u*L/mue)
const T physLength = 2.;         // physical domain length in each dimension
const T maxPhysT = 50.;          // max. simulation time in s, SI unit
const T flow_rate = 0.5;         // flow rate inside ractor
const T epsilon = 10e-8;         // stationary check
T C_0_0_phys = 0.;
T ly  = lx/N0;                   // height of the channel

// Verification type
bool AnalySol = true;                  // comparison with analytical solution

/// Initialize gnuplot
static Gnuplot<T> gploteoc( "concentration_eoc", false,"set terminal png size 720, 720 font 'Arial,10'",Gnuplot<T>::LOGLOG, Gnuplot<T>::LINREG);

// Reaction parameters
enum ReactionType {a2c, a2cAndBack};
constexpr ReactionType reactionType = a2c;

template<typename T, ReactionType type>
struct ReactionData
{
  int numReactions;                          // number of reactions
  std::vector<T> physReactionCoeff;          // reaction rate coefficients (Equilibrium reactions are defined as two separate reactions)
  int numComponents;                         // number of components
  std::vector<std::string> names;            // name of each component
  std::vector<T> stochCoeff;                 // stoichiometric coefficient according to {Reaction1: comp1, comp2, comp3,..., Reaction2: comp1, comp2, comp3,...}
                                             // component for which the kinetic was defined is required to have the coefficient of 1
  std::vector<T> reactionOrders;             // reaction orders per component and reaction for reaction rate power law
  std::vector<T> physInitialConcentrations;  // concentration per component at inlet

  ReactionData(){
    if constexpr (type == a2c) {
      numReactions = 1;
      physReactionCoeff = std::vector<T>({0.25});
      numComponents = 2;
      names = std::vector<std::string>({"Hydrogen","Hydrogene Iodine"});
      stochCoeff = std::vector<T>({-1,1.});
      reactionOrders = std::vector<T>({1., 0.});
      physInitialConcentrations = std::vector<T>({50.0, 0.0});

    }
    else if constexpr (type == a2cAndBack) {
      numReactions = 2;
      physReactionCoeff = std::vector<T>({0.25, 0.05});
      numComponents = 2;
      names = std::vector<std::string>({"Hydrogen","Hydrogene Iodine"});
      stochCoeff = std::vector<T>({-1, 1., 1., -1.});
      reactionOrders = std::vector<T>({1., 0., 0., 1.});
      physInitialConcentrations = std::vector<T>({50.0, 0.0});
    }
  }
};

ReactionData<T,ReactionType::a2c> reactionData;   //chooses reactionType


template <typename T>
class AdePhysConc2D : public AnalyticalF2D<T,T> {       //analytical solution
protected:
  T mue;
  T CA_0_phys = reactionData.physInitialConcentrations[0];
public:
  AdePhysConc2D(AdeUnitConverter<T, TDESCRIPTOR> converter) : AnalyticalF2D<T,T>(1),
    mue(converter.getPhysDiffusivity())
  { }

  bool operator()(T output[], const T input[]) override
  {
    T x = input[0];
    if constexpr (reactionType == a2c){
      T uD = flow_rate / mue;
      T lamda2 = 0.5 * (uD - util::sqrt(uD * uD + 4 * reactionData.physReactionCoeff[0]/mue));
      output[0] = CA_0_phys * util::exp(lamda2 * x);
    }
    else if constexpr (reactionType == a2cAndBack) {
      T Constant = reactionData.physReactionCoeff[1] * CA_0_phys / (reactionData.physReactionCoeff[0] + reactionData.physReactionCoeff[1]);
      T uD = flow_rate / mue;
      T lamda2 = 0.5 * (uD - util::sqrt(uD*uD + 4 * (reactionData.physReactionCoeff[0] + reactionData.physReactionCoeff[1]) / mue));
      output[0] = (CA_0_phys - Constant) * util::exp(lamda2*x) + Constant;
    }
    return true;
  }
};



// Stores geometry information in form of material numbers
void prepareGeometry(SuperGeometry<T,2>& superGeometry,
                     AdeUnitConverter<T, TDESCRIPTOR> &converter)
{
  OstreamManager clout(std::cout,"prepareGeometry");
  clout << "Prepare Geometry ..." << std::endl;

  Vector<T,2> extend;
  Vector<T,2> origin;

  superGeometry.rename( 0,2 );

  superGeometry.rename( 2,1,{1,0} );

  // defining size of inflow and outflow
  extend[1] = ly;
  extend[0]= ly/2;

  // Set material number for inflow
  origin[0] -= ly/4;
  IndicatorCuboid2D<T> inflow( extend, origin );
  superGeometry.rename( 2,3,1,inflow );

  // Set material number for outflow
  origin[0] = lx- ly/4;
  IndicatorCuboid2D<T> outflow( extend, origin );
  superGeometry.rename( 2,4,1,outflow );


  // Removes all not needed boundary voxels outside the surface
  superGeometry.clean();
  // Removes all not needed boundary voxels inside the surface
  superGeometry.innerClean();
  superGeometry.checkForErrors();

  superGeometry.print();

  clout << "Prepare Geometry ... OK" << std::endl;
}

// initiates lattice, initial conditions and boundary conditions
void prepareLattice(  SuperLattice<T, TDESCRIPTOR>*& ADlattice,
                      SuperGeometry<T,2>& superGeometry, T rho, T rhoEnd,
                      AdeUnitConverter<T, TDESCRIPTOR> converter )
{
  OstreamManager clout(std::cout,"prepareLattice");
  clout << "Prepare Lattice ..." << std::endl;

  const T omega = converter.getLatticeAdeRelaxationFrequency();

  // buffer layer
  ADlattice->defineDynamics<NoDynamics>(superGeometry, 0);

  // bulk
  ADlattice->defineDynamics<SourcedAdvectionDiffusionBGKdynamics>(superGeometry, 1);

  AnalyticalConst2D<T,T> u(converter.getLatticeVelocity(flow_rate),converter.getLatticeVelocity(0.0));
  AnalyticalConst2D<T,T> rho0_(converter.getLatticeDensity(C_0_0_phys));
  AnalyticalConst2D<T,T> rho_(converter.getLatticeDensity(rho));
  AnalyticalConst2D<T,T> rhoEnd_(converter.getLatticeDensity(rhoEnd));

  // Setting of the boundary conditions, inflow and outflow with Dirichlet condition according to analytical solution
  setAdvectionDiffusionTemperatureBoundary<T, TDESCRIPTOR>(*ADlattice, omega, superGeometry.getMaterialIndicator({3}));
  setAdvectionDiffusionTemperatureBoundary<T, TDESCRIPTOR>(*ADlattice, omega, superGeometry.getMaterialIndicator({4}));

  ADlattice->setParameter<descriptors::OMEGA>(converter.getLatticeAdeRelaxationFrequency());

  auto bulkIndicator = superGeometry.getMaterialIndicator({1});
  auto everywhere = superGeometry.getMaterialIndicator({1,2,3,4});

  // setting flowrate inside reactor, as well as initial concentrations
  ADlattice->defineField<descriptors::VELOCITY>(everywhere, u );
  ADlattice->defineField<descriptors::SOURCE>(bulkIndicator, rho0_);
  ADlattice->defineRho( bulkIndicator, rho0_);
  ADlattice->iniEquilibrium( superGeometry,1, rho0_, u );

  // setting actual values for the boundary
  ADlattice->defineRho( superGeometry, 3, rho_);
  ADlattice->iniEquilibrium( superGeometry, 3, rho_, u );
  ADlattice->defineRho( superGeometry, 4, rhoEnd_);
  ADlattice->iniEquilibrium( superGeometry, 4, rhoEnd_, u );

  // Make the lattice ready for simulation
  ADlattice->initialize();

  clout << "Prepare Lattice ... OK" << std::endl;
}


// calcualtes the error between numerical and analytical solution along the centerline
// return value is the average L2 error along the centerline
T errorOverLine( std::vector<SuperLatticeDensity2D<T, TDESCRIPTOR>*>& densities,
                 SuperGeometry<T,2>& superGeometry,
                 int iT,
                 AdeUnitConverter<T, TDESCRIPTOR> converter)
{
  OstreamManager clout(std::cout,"error calculation");
  Gnuplot<T> gplot( "centerConcentration" );
  gplot.setLabel("L / m", "concentration / mol/m^2");
  Gnuplot<T> gplot2( "Error Concentration" );
  gplot2.setLabel("L / m", "concentration / mol/m^2");

  T dist = converter.getPhysDeltaX();

  AnalyticalFfromSuperF2D<T> interpolation_Concentration_A( *densities[0], true, 1 ); //numerical density
  AnalyticalFfromSuperF2D<T> interpolation_Concentration_C( *densities[1], true, 1 );
  AdePhysConc2D<T> concentration_analytical( converter);

  T resultA[1] = {0.};    // will contain result of analytical solution
  T resultC[1] = {0.};
  T resultSimA[1] = {0.}; // will contain result of lattice evaluation
  T resultSimC[1] = {0.};

  T ndatapoints = converter.getResolution(); // number of data points on line
  std::vector<T> Error(2,T());
  T tempMSEA = .0;                              // mean squared error for A
  T tempSquaredSumA = .0;                       // sum of analytical concentration for A

  T tempMSEC = .0;                              // mean squared error for C
  T tempSquaredSumC = .0;                       // sum of analytical concentration for C
  for (int i = 0; i <= ndatapoints; i++) {
    T input[2] =  {i*dist, ly/2};
    interpolation_Concentration_A(resultSimA, input);
    interpolation_Concentration_C(resultSimC, input);
    concentration_analytical(resultA, input);
    resultC[0]= reactionData.physInitialConcentrations[0]-resultA[0];
    tempMSEA += (resultA[0]  - resultSimA[0]) * (resultA[0]  - resultSimA[0]);
    tempSquaredSumA += (resultA[0]) * (resultA[0]);
    tempMSEC += (resultC[0]  - resultSimC[0]) * (resultC[0]  - resultSimC[0]);
    tempSquaredSumC += (resultC[0]) * (resultC[0]);

    gplot.setData( input[0], {resultSimA[0], resultA[0], resultSimC[0], resultC[0]},
          {reactionData.names[0] + " simulated",reactionData.names[0] + " analytical", reactionData.names[1]+" simulated" ,reactionData.names[1] + " analytical"});
    gplot2.setData( input[0],{ (util::abs(resultSimA[0]-resultA[0])), (util::abs(resultSimC[0]-resultC[0]))}, {"absolute Error Concentration" + reactionData.names[0], "absolute Error Concentration" + reactionData.names[1]});
  }
  Error[0]= util::sqrt(tempMSEA/tempSquaredSumA);
  Error[1]= util::sqrt(tempMSEC/tempSquaredSumC);

  gplot2.writePDF({"Plot2"});
  gplot.writePDF({"Plot1"}); // plot is generated
  singleton::directories().setOutputDir("./tmp/");
  CSV<T> csvWriterErr;
  csvWriterErr.writeDataFile( converter.getResolution(), Error, "averageL2RelError", 16);

  clout << "Average relative L2 error for species A " << std::to_string(util::sqrt(tempMSEA/tempSquaredSumA)) << std::endl;
  clout << "Average relative L2 error for species C " << std::to_string(util::sqrt(tempMSEC/tempSquaredSumC)) << std::endl;

  return util::sqrt(tempMSEA/tempSquaredSumA);
}

// creating vtk output as well as verification plots
void getResults(int statIter, AdeUnitConverter<T, TDESCRIPTOR> converter,
                std::vector<SuperLattice<T, TDESCRIPTOR>*>& adlattices,
                int iT,
                SuperGeometry<T,2>& superGeometry,
                util::Timer<T>& timer , int iTmax, bool Conc_stationary)
{
  if (iT%statIter == 0 || Conc_stationary == true || iT == iTmax) {
    OstreamManager clout(std::cout,"getResults");

    // initiate vtk Data
    SuperVTMwriter2D<T> vtkWriter("advectionDiffusionReaction2d");

    // insert concentrations and velocoties for each component into vector and vtk
    std::vector<SuperLatticePhysVelocity2D<T, TDESCRIPTOR>*> velocities;
    std::vector<SuperLatticeDensity2D<T, TDESCRIPTOR>*> densities;

    for (int i = 0; i<reactionData.numComponents;i++){
      adlattices[i]->communicate();
      velocities.emplace_back( new SuperLatticePhysVelocity2D<T, TDESCRIPTOR> (*adlattices[i], converter));
      densities.emplace_back( new SuperLatticeDensity2D<T, TDESCRIPTOR> (*adlattices[i]));
      velocities[i]->getName() = "velocity " + reactionData.names[i];
      densities[i]->getName() = "Concentration " + reactionData.names[i];
      vtkWriter.addFunctor(*velocities[i]);
      vtkWriter.addFunctor(*densities[i]);
    }

  if (iT == 0) {
    /// Writes the converter log file
    // writeLogFile(converter,"rayleighBenard2d");

    /// Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry2D<T, TDESCRIPTOR> geometry(*adlattices[0], superGeometry);
    SuperLatticeCuboid2D<T, TDESCRIPTOR> cuboid(*adlattices[0]);
    SuperLatticeRank2D<T, TDESCRIPTOR> rank(*adlattices[0]);

    vtkWriter.write(geometry);
    vtkWriter.write(cuboid);
    vtkWriter.write(rank);

    vtkWriter.createMasterFile();
  }

    /// Writes the VTK files and prints statistics
    if (iT%statIter == 0 || Conc_stationary == true) {
      /// Timer console output
      timer.update(iT);
      timer.printStep();
      /// Lattice statistics console output
      adlattices[0]->getStatistics().print(iT,converter.getPhysTime(iT));

      if (AnalySol == true){
        AdePhysConc2D<T> concentration_analytical( converter);
        SuperLatticeFfromAnalyticalF2D<T, TDESCRIPTOR> analyticalA_lattice(concentration_analytical, *adlattices[0]);
        AnalyticalConst2D<T,T> rho0_(converter.getLatticeDensity(reactionData.physInitialConcentrations[0]));
        SuperLatticeFfromAnalyticalF2D<T, TDESCRIPTOR> rho_lattice(rho0_, *adlattices[0]);
        SuperLatticeField2D<T, TDESCRIPTOR, SOURCE> source (*adlattices[0]);
        vtkWriter.addFunctor(analyticalA_lattice, "Analytical Concentration " + reactionData.names[0]);
        vtkWriter.addFunctor(rho_lattice - analyticalA_lattice,"Analytical Conentration " + reactionData.names[1]);
        vtkWriter.addFunctor(analyticalA_lattice - *densities[0], "Error " + reactionData.names[0]);
        vtkWriter.addFunctor(rho_lattice - analyticalA_lattice - *densities[1], "Error " + reactionData.names[1]);
        vtkWriter.addFunctor(source);

        vtkWriter.write(iT);
      }
      else{
        vtkWriter.write(iT);
     }
    }

    /// create Gnuplot
    if ( iT == iTmax || Conc_stationary == true){
      OstreamManager clout(std::cout,"Gnuplot");
      clout << "Gnuplot ..." << std::endl;

      if (AnalySol == true ){
        T simulationAverage = .0;
        simulationAverage = errorOverLine(densities, superGeometry, iT, converter); //computes error along the middle of the PRF and saves also the concetrations along it
       if (runs >1){  // EOC Plot
           singleton::directories().setOutputDir("./tmp/");
           gploteoc.setData( T(converter.getResolution()) , {simulationAverage}, {"average L2  Error"}, "top right", {'p'});
           gploteoc.writePNG();
       }
      }

      else {         // AnalySol == false means no comparison with the analytical solution (no error plots)
        Gnuplot<T> gplot( "centerConcentration" );
        gplot.setLabel("L / m", "concentration / mol/m^2");
        // Vectors for simulated solution
        T Concentration_i[1] = {T()};
        std::vector<T>yValues;
        std::vector<std::string>names_Gnuplot;
        for (int i=0; i< reactionData.numComponents; i++){
          yValues.emplace_back(0.);
          names_Gnuplot.emplace_back(reactionData.names[i]);
        }
        T dist = converter.getPhysDeltaX();
        T ndatapoints = converter.getResolution(); // number of data points on line
        CSV<T> csvWriterConcentration;
         //save concentration along the middle of the PRF
        for (int n = 0; n <= ndatapoints; n++) {
          T input[2] =  {n*dist, ly/2};
          for (int i=0; i< reactionData.numComponents; i++){
            AnalyticalFfromSuperF2D<T> interpolation_Concentration( *densities[i], true, 1 );
            interpolation_Concentration(Concentration_i, input);
            yValues[i] = Concentration_i[0];
            csvWriterConcentration.writeDataFile( input[0], Concentration_i[0], "simulation" + reactionData.names[i] , 16);
          }

          gplot.setData( input[0], yValues,
          names_Gnuplot);
        }
        gplot.writePDF({"Plot1"}); // plot is generated
      }

    clout << "Gnuplot done" << std::endl;
    }
  }
}

void simulate(int N, int statIter, T ly, T mue, T peclet, Gnuplot<T>& gploteoc)
{
OstreamManager clout(std::cout,"simulate");
  clout << "Executing the simulation with N=" << std::to_string(N) << std::endl;
   if (runs > 1) {
    singleton::directories().setOutputDir("./tmp/N" + std::to_string(N)  + "/");
  }

  AdeUnitConverter<T,TDESCRIPTOR> converter(
    lx/N,            // physDeltaX
    util::pow(lx/N, 2),    // physDeltaT (diffusive scaling)
    lx,              // charPhysLength
    peclet*mue/lx,   // charPhysVelocity from Peclet */

    mue,             // physDiffusivity
    1                // physDensity,
  );
  converter.print();

  std::vector<T> latticeReactionCoeff(reactionData.numReactions);
  std::transform(reactionData.physReactionCoeff.begin(), reactionData.physReactionCoeff.end(), latticeReactionCoeff.begin(), [&](auto& c){
    return c * converter.getConversionFactorTime();
  });

  /// === 2nd Step: Prepare Geometry ===
  std::vector<T> extend(2,T());
  extend[0] = lx;
  extend[1] = ly;
  std::vector<T> origin(2,T());
  IndicatorCuboid2D<T> cuboid(extend, origin);

  // Instantiation of a cuboidGeometry with weights
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif
  CuboidGeometry2D<T> cuboidGeometry(cuboid, converter.getPhysDeltaX(), noOfCuboids);

  // periodic boundary in y-direction -> pipe of uncertain diameter
  cuboidGeometry.setPeriodicity(false, true);

  HeuristicLoadBalancer<T> loadBalancer(cuboidGeometry);

  SuperGeometry<T,2> superGeometry(cuboidGeometry, loadBalancer, 2);

  prepareGeometry(superGeometry, converter);

  // === 3rd Step: Prepare Lattice ===
  SuperLattice<T, TDESCRIPTOR> ADlattice(superGeometry);
  std::vector<T> boundaryEnd(2, T());
  AdePhysConc2D<T> concentration_analytical( converter);
  T input[2]= {lx,ly};
  T tmpA[1] = {0.};
  concentration_analytical(tmpA, input);
  boundaryEnd[0]= tmpA[0];
  boundaryEnd[1] = reactionData.physInitialConcentrations[0]-boundaryEnd[0];

  // initiate reference lattice for coupling, this lattice will be the first component
  std::vector<SuperLattice<T, TDESCRIPTOR>*> adlattices;
  adlattices.push_back(&ADlattice);

  // initiate further lattices for the remaining components
  for (int i =1; i<reactionData.numComponents;i++){
    adlattices.emplace_back( new SuperLattice<T, TDESCRIPTOR> (superGeometry));
  }
  // partners are used for coupling; coupled with adlattice[0]
  std::vector<SuperLattice<T, TDESCRIPTOR>*> partners;
  for (int i =1; i<reactionData.numComponents; i++){
    partners.emplace_back(adlattices[i]);
  }

  // prepare all lattices with the respective concentration
  for(int i = 0; i<reactionData.numComponents; i++){
    prepareLattice(adlattices[i],
      superGeometry,reactionData.physInitialConcentrations[i], boundaryEnd[i],
      converter);
  }

  /// === 3.1 Step: Prepare Coupling ===
  ConcentrationAdvectionDiffusionCouplingGenerator2D<T,TDESCRIPTOR> coupling(
    0, converter.getLatticeLength(lx), 0, converter.getLatticeLength(ly),
    reactionData.stochCoeff, latticeReactionCoeff, reactionData.reactionOrders);

  //coupling of lattice[0] with the other lattices (=partners)
  ADlattice.addLatticeCoupling(coupling, partners);
  std::size_t iTmax = converter.getLatticeTime(maxPhysT);

  util::Timer<T> timer( iTmax, superGeometry.getStatistics().getNvoxel() );
  timer.start();

  std::vector<util::ValueTracer<T>*>converge;


  for (int i=0; i< reactionData.numComponents; i++){
    converge.emplace_back(new util::ValueTracer<T> (converter.getLatticeTime(1.0),epsilon)) ;
    }

  for (std::size_t iT = 0; iT <= iTmax; ++iT) {
    /// === 6th Step: Collide and Stream Execution ===
    ADlattice.executeCoupling();
    for (int i = 0; i<reactionData.numComponents; i++){
      adlattices[i]->collideAndStream();
    }

    bool hasConverged = true;
    //hasConverged() has to be true for ALL lattices when the simulation ends before iTmax
    for (int i=0; i<reactionData.numComponents && hasConverged == true; i++) {
      hasConverged = converge[i]->hasConverged();
    }
    //stationary Check
    if (hasConverged) {
      getResults(statIter, converter, adlattices, iT, superGeometry, timer, iTmax, hasConverged);
      OstreamManager clout(std::cout,"stationary check");
      clout << "stationary condition reached -> ending simulation" << std::endl;
      break;
    }

    getResults(statIter, converter, adlattices, iT, superGeometry, timer, iTmax, hasConverged);

    for (int i=0; i<reactionData.numComponents; i++){
      converge[i]->takeValue(adlattices[i]->getStatistics().getAverageRho()); //stationary Check
    }
  }
}

int main(int argc, char *argv[])
{
  /// === 1st Step: Initialization ===
  OstreamManager clout(std::cout,"main");
  olbInit(&argc, &argv);
  singleton::directories().setOutputDir("./tmp/");

  gploteoc.setLabel("Resolution test", "average Error");

  for (int i = 0; i < runs; ++i) {
    ly = (lx/(N0+50*i));
    //increasing resolution +50 for each run
    simulate(N0+50*i,
             (N0+50*i)*5,
             ly,
             mue,
             peclet, gploteoc);
  }
}
