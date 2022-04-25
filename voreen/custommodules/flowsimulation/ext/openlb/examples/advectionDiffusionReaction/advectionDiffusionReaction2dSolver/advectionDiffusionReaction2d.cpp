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
using namespace olb::names;

using TDESCRIPTOR = descriptors::D2Q5<VELOCITY, SOURCE>;

constexpr unsigned numComponents = 2;

using LATTICES = meta::map<
  Concentration<0>, TDESCRIPTOR,
  Concentration<1>, TDESCRIPTOR
>;

// Reaction parameters
enum ReactionType {a2c, a2cAndBack};              // decide between forward and reverse reaction
constexpr ReactionType reactionType = ReactionType::a2cAndBack;

template<typename T>
struct SimulationErrors : public parameters::ResultsBase
{
  T averageRelativeL2Error   {0};
};


template<typename T>
struct ReactionData
{
  unsigned numReactions;                     // number of reactions
  std::vector<T> physReactionCoeff;          // reaction rate coefficients (Equilibrium reactions are defined as two separate reactions)
  unsigned numComponents;                    // number of components
  std::vector<std::string> names;            // name of each component
  std::vector<T> stoichiometricCoeff;        // stoichiometric coefficient according to {Reaktion1: comp1, comp2, comp3,..., Reaktion2: comp1, comp2, comp3,...}
                                             // component for which the kinetic was defined is required to have the coefficient of 1
  std::vector<T> reactionOrders;             // reaction orders per component and reaction for reaction rate power law
  std::vector<T> physInitialConcentrations;  // concentration per component at inlet and inside reactor at the beginning


  ReactionData(ReactionType type) {
    if (type == a2c) {
      numReactions = 1;
      physReactionCoeff = std::vector<T>({0.25});
      numComponents = 2;
      names = std::vector<std::string>({"Hydrogen","Hydrogene Iodine"});
      stoichiometricCoeff = std::vector<T>({-1,1.});
      reactionOrders = std::vector<T>({1., 0.});
      physInitialConcentrations = std::vector<T>({50.0, 0.0});
      }

    else if (type == a2cAndBack) {
      numReactions = 2;
      physReactionCoeff = std::vector<T>({0.25, 0.05});
      numComponents = 2;
      names = std::vector<std::string>({"Hydrogen","Hydrogene Iodine"});
      stoichiometricCoeff = std::vector<T>({-1, 1., 1., -1.});
      reactionOrders = std::vector<T>({1., 0., 0., 1.});
      physInitialConcentrations = std::vector<T>({50.0, 0.0});
      }

  }

};

template<typename T>
struct ReactionSimulationParameters
 : public parameters::SimulationBase<T>
{
  const T lx  = 10.;                  // length of the channel
  int N;                              // resolution of the model
  //const int statIter0;           // initial # lattice output timesteps
  const T mue = 0.1;//2.3e-4;                 // physical diffusivity
  const T peclet = 100.;              // Peclet number (Pe = u*L/mue)
  const T physLength = 2.;            // physical domain length in each dimension
  const T flow_rate = 0.5;              // flow rate inside ractor
  T ly;                              // height of the channel
  const int runs = 1;        // # simulations with increasing resolution


  std::shared_ptr<ReactionData<T> const> reactionData;

  std::shared_ptr<AdeUnitConverter<T, TDESCRIPTOR> const>
    converter;

  ReactionSimulationParameters( ReactionType type, int N_=200)
   : N(N_), ly(lx/N_)
  {

   reactionData = std::make_shared<ReactionData<T> const> (type);

   this->converter = std::make_shared<AdeUnitConverter<T,TDESCRIPTOR> const> (
      lx/N,            // physDeltaX
      util::pow(lx/N, 2),    // physDeltaT (diffusive scaling)
      lx,              // charPhysLength
      peclet*mue/lx,   // charPhysVelocity from Peclet
      mue,             // physDiffusivity
      1                // physDensity,
    );

  this->maxTime = 50.0;
  }
};

template<typename T>
using Parameters = meta::map<
  Simulation, ReactionSimulationParameters<T>,
  Stationarity, parameters::Stationarity<T,Concentration<0>, Concentration<1>>,
  Output, parameters::OutputBase<T,Concentration<0>>,
  Errors, SimulationErrors<T>
>;

// Analytical solution
template <typename T>
class AdePhysConc2D : public AnalyticalF2D<T,T> {
protected:
  ReactionType _type;
  const ReactionData<T>& _data;
  T mue;
  T _flow_rate;
  T CA_0_phys = _data.physInitialConcentrations[0];
public:
  AdePhysConc2D(
    const AdeUnitConverter<T,TDESCRIPTOR>& converter, ReactionType type, const ReactionData<T>& data, T flow_rate)
   : AnalyticalF2D<T,T>(1), _type(type), _data(data), mue(converter.getPhysDiffusivity()), _flow_rate(flow_rate),
    CA_0_phys(data.physInitialConcentrations[0])
  { }

  bool operator()(T output[], const T input[]) override
  {
    if (_type == a2c){
      T uD =_flow_rate / mue;
      T lamda2 = 0.5 * (uD - util::sqrt(uD * uD + 4 * _data.physReactionCoeff[0]/mue));
      output[0] = CA_0_phys * util::exp(lamda2 * input[0]);
    }
    else if (_type == a2cAndBack) {
      T Constant = _data.physReactionCoeff[1] * CA_0_phys / (_data.physReactionCoeff[0] + _data.physReactionCoeff[1]);
      T uD = _flow_rate/ mue;
      T lamda2 = 0.5 * (uD - util::sqrt(uD*uD + 4 * (_data.physReactionCoeff[0] + _data.physReactionCoeff[1]) / mue));
      output[0] = (CA_0_phys - Constant) * util::exp(lamda2*input[0]) + Constant;
    }
    return true;
  }
};


template<typename T>
class Reaction2dSolver : public LBSolver<T, Parameters<T>, LATTICES> {
 private:
  mutable OstreamManager            clout {std::cout, "Reaction2dSolver"};

public:
  Reaction2dSolver(utilities::TypeIndexedTuple<Parameters<T>> params)
   : Reaction2dSolver::LBSolver(params)
   { }


protected:
  void prepareGeometry() override
  {
    std::vector<T> extend(2,T());
    extend[0] = this->parameters(Simulation()).lx;
    extend[1] = this->parameters(Simulation()).ly;
    const std::vector<T> origin(2,T());
    IndicatorCuboid2D<T> cuboid(extend, origin);

    // Instantiation of a cuboidGeometry with weights
  #ifdef PARALLEL_MODE_MPI
    const int noOfCuboids = this->parameters(Simulation()).noC * singleton::mpi().getSize();
  #else
    const int noOfCuboids = 1;
  #endif
    this->_cGeometry = std::make_shared<CuboidGeometry2D<T>>(
      cuboid,
      this->converter().getPhysDeltaX(),
      noOfCuboids);

    // periodic boundary in y-direction -> pipe of uncertain diameter
    this->_cGeometry->setPeriodicity(false, true);

    this->_loadBalancer
     = std::make_shared<HeuristicLoadBalancer<T>> (*this->_cGeometry);

    this->_sGeometry = std::make_shared<SuperGeometry<T,2>> (
      *this->_cGeometry,
      *this->_loadBalancer,
      this->parameters(Simulation()).overlap);

    this->geometry().rename( 0,2 );
    this->geometry().rename( 2,1,{1,0} );

    // defining size of inflow and outflow
    Vector<T,2> extend_bdry;
    Vector<T,2> origin_bdry;
    extend_bdry[1] = this->parameters(Simulation()).ly;
    extend_bdry[0]= this->parameters(Simulation()).ly / 2;

    // Set material number for inflow
    origin_bdry[0] -= this->parameters(Simulation()).ly / 4;
    IndicatorCuboid2D<T> inflow( extend_bdry, origin_bdry );
    this->geometry().rename( 2,3,1,inflow );

    // Set material number for outflow
    origin_bdry[0]
     = this->parameters(Simulation()).lx- this->parameters(Simulation()).ly / 4;
    IndicatorCuboid2D<T> outflow( extend_bdry, origin_bdry );
    this->geometry().rename( 2,4,1,outflow );
  }

  void prepareLattices() override
  {
    const auto& params = this->parameters(Simulation());
    const auto& reactionParams = *(params.reactionData);


    std::vector<T> latticeReactionCoeff(
      reactionParams.numReactions);
    std::transform(
      reactionParams.physReactionCoeff.begin(),
      reactionParams.physReactionCoeff.end(),
      latticeReactionCoeff.begin(), [&](auto& c){
      return c * this->converter().getConversionFactorTime();
    });

    //Coupling of Lattices for reaction rate calculation
    ConcentrationAdvectionDiffusionCouplingGenerator2D<T,TDESCRIPTOR> coupling(
      0, this->converter().getLatticeLength(params.lx),
      0, this->converter().getLatticeLength(params.ly),
      reactionParams.stoichiometricCoeff, latticeReactionCoeff, reactionParams.reactionOrders);
    std::vector <SuperLattice<T,TDESCRIPTOR>*> lattices;
    meta::tuple_for_each(this->_sLattices, [&](auto& lattice, unsigned iLattice) {
      if (iLattice != 0) {
        lattices.emplace_back(lattice.get());
        }
    });
    this->lattice(Concentration<0>()).addLatticeCoupling(coupling, lattices);

    const T omega = this->converter().getLatticeRelaxationFrequency();

    meta::tuple_for_each(this->_sLattices, [&](auto& lattice, unsigned iLattice) {
      // buffer layer
      lattice->template defineDynamics<NoDynamics>(this->geometry(), 0);
      // Dynamics for the use of a source term
      // Material=1 -->bulk dynamics
      lattice->template defineDynamics<SourcedAdvectionDiffusionBGKdynamics>(
      this->geometry().getMaterialIndicator({1}));

      // Setting of the boundary conditions, inflow with Dirichlet condition and outflow calculated from analytical solution
      setAdvectionDiffusionTemperatureBoundary<T,TDESCRIPTOR>(
        *lattice, omega, this->geometry().getMaterialIndicator({3}));
      setAdvectionDiffusionTemperatureBoundary<T,TDESCRIPTOR>(
        *lattice, omega, this->geometry().getMaterialIndicator({4}));

      lattice->template setParameter<descriptors::OMEGA>(this->converter().getLatticeAdeRelaxationFrequency());
    });
}

  void setInitialValues() override
  {
    const auto& converter = this->converter();
    const auto& params = this->parameters(Simulation());
    const auto& reactionParams = *(params.reactionData);
    std::vector<T> boundaryEnd(2, T());
      AdePhysConc2D<T> concentration_analytical(
        converter, reactionType, reactionParams, params.flow_rate);
      T input[2]= {params.lx,params.ly};
      T tmpA[1] = {0.};
      concentration_analytical(tmpA, input);
      boundaryEnd[0]= tmpA[0];
      boundaryEnd[1] = reactionParams.physInitialConcentrations[0]-boundaryEnd[0];


 meta::tuple_for_each(this->_sLattices, [&](auto& lattice, unsigned iLattice) {
      AnalyticalConst2D<T,T> u(converter.getLatticeVelocity(params.flow_rate),
        converter.getLatticeVelocity(0.0));
      AnalyticalConst2D<T,T> rho0_(converter.getLatticeDensity(0.0));
      AnalyticalConst2D<T,T> rho_(converter.getLatticeDensity(
        reactionParams.physInitialConcentrations[iLattice]));
      AnalyticalConst2D<T,T> rhoEnd_(converter.getLatticeDensity(boundaryEnd[iLattice]));
      AnalyticalConst2D<T,T> id_lattice (iLattice);

      auto bulkIndicator = this->geometry().getMaterialIndicator({1});
      auto everywhere = this->geometry().getMaterialIndicator({1,2,3,4});

      // setting flowrate inside reactor, as well as initial concentrations
      lattice->template defineField<descriptors::VELOCITY>( everywhere, u );
      lattice->template defineField<descriptors::SOURCE>( everywhere, rho0_ );
      lattice->defineRho( bulkIndicator, rho0_);
      lattice->iniEquilibrium( this->geometry(),1, rho0_, u );

      lattice->defineRho( this->geometry(), 3, rho_);
      lattice->iniEquilibrium( this->geometry(), 3, rho_, u );

      lattice->defineRho( this->geometry(), 4, rhoEnd_);
      lattice->iniEquilibrium( this->geometry(), 4, rhoEnd_, u );

  });

 }

  void setBoundaryValues(std::size_t iT) override { }

  void writeImages (std::size_t iT) const override { }

   void writeVTK(std:: size_t iT) const override
   {
    const auto& params = this->parameters(Simulation());
    SuperVTMwriter2D<T> vtkWriter(this->parameters(Output()).vtkFilename);

    std::vector<SuperLatticeDensity2D<T,TDESCRIPTOR>*> densities;

    unsigned counter (0);
    meta::tuple_for_each(this->_sLattices, [&](auto& lattice, unsigned iLattice) {
      lattice->communicate();
      densities.emplace_back( new SuperLatticeDensity2D (*lattice));
      densities[counter]->getName() = "Concentration " + params.reactionData->names[iLattice];

      vtkWriter.addFunctor(*densities[counter]);
      ++counter;
    });

      auto& sLattice = this->lattice(Concentration<0>());

      AdePhysConc2D<T> concentration_analytical (
        this->converter(), reactionType,
        *(params.reactionData),
        params.flow_rate);
      SuperLatticeFfromAnalyticalF2D<T, TDESCRIPTOR> analyticalA_lattice(concentration_analytical, sLattice);
      AnalyticalConst2D<T,T> rho0_(this->converter().getLatticeDensity(params.reactionData->physInitialConcentrations[0]));
      SuperLatticeFfromAnalyticalF2D<T, TDESCRIPTOR> rho_lattice(rho0_, sLattice);
      SuperLatticeField2D<T, TDESCRIPTOR, SOURCE> source (this->lattice(Concentration<1>()));

      vtkWriter.addFunctor(analyticalA_lattice, "Analytical Concentration " + params.reactionData->names[0]);
      vtkWriter.addFunctor(rho_lattice - analyticalA_lattice,"Analytical Conentration " + params.reactionData->names[1]);
      vtkWriter.addFunctor(analyticalA_lattice - *densities[0], "Error " + params.reactionData->names[0]);
      vtkWriter.addFunctor(rho_lattice - analyticalA_lattice - *densities[1], "Error " + params.reactionData->names[1]);
      vtkWriter.addFunctor(source);

      vtkWriter.write( iT );

    }

  void writeGnuplot() const override {
    //set up to use Gnuplot
    const auto& params = this->parameters(Simulation());
    Gnuplot<T> gplot( "centerConcentration" );
    gplot.setLabel("L / m", "concentration / mol/m^2");

    T dist = this->converter().getPhysDeltaX();
    T ndatapoints = this->converter().getResolution(); // number of data points on line

    std::vector<SuperLatticeDensity2D<T,TDESCRIPTOR>*> densities;
    unsigned counter (0);
    meta::tuple_for_each(this->_sLattices, [&](auto& lattice, unsigned iLattice) {
      densities.emplace_back( new SuperLatticeDensity2D (*lattice));
      densities[counter]->getName() = "Concentration " + params.reactionData->names[iLattice];
      ++counter;
    });


        Gnuplot<T> plt("UNUSED_FILE");
        Gnuplot<T> gplot2( "Error Concentration" );
        gplot2.setLabel("L / m", "concentration / mol/m^2");

        AnalyticalFfromSuperF2D<T> interpolation_Concentration_A( *densities[1], true, 1 );
        AnalyticalFfromSuperF2D<T> interpolation_Concentration_C( *densities[0], true, 1 );
        AdePhysConc2D<T> concentration_analytical (
          this->converter(),
          reactionType,
          *(this->parameters(Simulation()).reactionData),
           this->parameters(Simulation()).flow_rate);

        T resultA[1] = {0.};    // will contain result of analytical solution
        T resultC[1] = {0.};    //
        T resultSimA[1] = {0.}; // will contain result of lattice evaluation
        T resultSimC[1] = {0.};
        T tempMSEA = .0;                              // mean squared error for species A
        T tempSquaredSumA = .0;                       // sum of analytical concentration for species A
        T tempMSEC = .0;                              // mean squared error for species C
        T tempSquaredSumC = .0;                       // sum of analytical concentration for species C

       for (int i = 0; i <= ndatapoints; ++i) {
          T input[2] =  {i*dist, params.ly/2};

          interpolation_Concentration_A(resultSimA, input);
          interpolation_Concentration_C(resultSimC, input);
          concentration_analytical(resultA, input);
          resultC[0]=  params.reactionData->physInitialConcentrations[0]-resultA[0];
          tempMSEA += (resultA[0]  - resultSimA[0]) * (resultA[0]  - resultSimA[0]);
          tempSquaredSumA += (resultA[0]) * (resultA[0]);
          tempMSEC += (resultC[0]  - resultSimC[0]) * (resultC[0]  - resultSimC[0]);
          tempSquaredSumC += (resultC[0]) * (resultC[0]);

          gplot.setData( input[0], {resultSimA[0], resultA[0], resultSimC[0], resultC[0]},
          {params.reactionData->names[0] + " simulated", params.reactionData->names[0] + " analytical", params.reactionData->names[1]+" simulated" , params.reactionData->names[1] + " analytical"});
          gplot2.setData( input[0],{ (util::abs(resultSimA[0]-resultA[0])), (util::abs(resultSimC[0]-resultC[0]))}, {"absolute Error Concentration" + params.reactionData->names[0], "absolute Error Concentration" + params.reactionData->names[1]});

        }

        CSV<T> csvWriterErr;
        csvWriterErr.writeDataFile(this->converter().getResolution(),
                                   {util::sqrt(tempMSEA/tempSquaredSumA), util::sqrt(tempMSEC/tempSquaredSumC)}, "averageL2RelError", 16);
        gplot2.writePDF({"Plot2"});
        gplot.writePDF({"Plot1"}); // plot is generated

        clout << "Average relative L2 error for species A " << std::to_string(util::sqrt(tempMSEA/tempSquaredSumA)) << std::endl;
        clout << "Average relative L2 error for species C " << std::to_string(util::sqrt(tempMSEC/tempSquaredSumC)) << std::endl;

        this->parameters(Errors()).averageRelativeL2Error = util::sqrt(tempMSEA/tempSquaredSumA);

  }

};

int main(int argc, char *argv[])
{
  using T = double;

  olbInit (&argc, &argv );

  OstreamManager clout (std::cout, "Main Loop EOC");

  Gnuplot<T> gploteoc(
    "concentration_eoc", false,
    "set terminal png size 720, 720 font 'Arial,10'",
     Gnuplot<T>::LOGLOG, Gnuplot<T>::LINREG);
  gploteoc.setLabel("Resolution test", "average Error");

  T runs = 3;         // for EOC plot, number of runs with increasing resolution
  int N;

  for (int i = 0; i < runs; ++i) {
    N = 200+50*i;               // increasing spatial resolution
    clout << "Executing the simulation with N=" << std::to_string(N) << std::endl;

    std::string output("./tmp/N" + std::to_string(N)  + "/");
    singleton::directories().setOutputDir(output);

    utilities::TypeIndexedTuple<Parameters<T>> params;
    params.template get<Simulation>()= std::make_shared<ReactionSimulationParameters<T>>(ReactionType::a2cAndBack, N);
    params.template get<Stationarity>() = std::make_shared<parameters::Stationarity<T,Concentration<0>,Concentration<1>>>(
      parameters::Stationarity<T,Concentration<0>, Concentration<1>>:: AverageRho, 1.0, 1e-8);
    params.template get<Errors>() = std::make_shared<SimulationErrors<T>>();
    params.template get<Output>() = std::make_shared<parameters::OutputBase<T, names::Concentration<0>>>(
      "advectionDiffusionReaction2d", "../../../", output,
      true, true, 1., 0,
      true, "advectionDiffusionReaction2d", 1.,
      false, "", 0.,
      true, "advectionDiffusionReaction2d", 1.
    );

    Reaction2dSolver<T> reaction2d(params);

    reaction2d.solve();

    const auto& error = reaction2d.parameters(Errors());
    singleton::directories().setOutputDir("./tmp/");

    gploteoc.setData( T(N) , {error.averageRelativeL2Error},
           {"average L2  Error"}, "top right", {'p'});

    gploteoc.writePNG();


  }

 }

