/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2020 Johannes Nguyen, Stephan Simonis, Sam Avis
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
 
/* fourRollMill2d.cpp
 * A spherical domain of fluid phase II is immersed in a cuboid
 * filled with fluid phase I. The bottom left and top right
 * cylinders begin to spin in counterclock-wise direction.
 * The top left and bottom right cylinders spin in clock-wise direction-
 * The velocity field of extensional type deforms the droplet accordingly.
 *
 * The droplet breakup case should be run on parallel mode
 * because of increased computation time.
 */

#include "olb2D.h"
#include "olb2D.hh"

using namespace olb;
using namespace olb::descriptors;
using namespace olb::graphics;

typedef double T;
#define DESCRIPTOR D2Q9<CHEM_POTENTIAL,FORCE>

// Guideline parameters:
// Ca = 0.1  for steady state
// Ca = 0.42 for drop breakup
// Recommended xi = 1.14, min 1.0

//Set physical input parameters
const T radius = 1.;                                // [pu] physical droplet radius
const T nx = 20.;                                   // [pu] domain length H = 20
const T cylinderOffset = nx / 2. / 3.;              // [pu] cylinder offset from x- and y- axis, symmetrical!
const T radiusCylinder = cylinderOffset * 0.625;    // [pu] radius of the cylinder

const T contVisc = 1.;                              // [pu] physical(!) kinematic viscosity of continuous phase
const T physDensity = 1.;                           // [pu] physical density

// Set numerical parameters
const T N = 10;                                     // droplet radius [lu] = N/radius
const T latticeRelaxationTime = 1.;                 // [lu] tau for first lattice (f)
const T latticeRelaxationTime_g = 1.;               // [lu] caution: not used..., i.e. tau_f = tau_g

// Set (model) dimensionless numbers
const T dropRe = 0.0625;                            // definition: extensionRate * radius^2 / continuousViscosity
//T capillaryNr = 0.1;                              // definition: radius * extensionRate * viscosity(non kinematic) / surfTen
const T capillaryNr = 0.42;                         // capillary number for breakup

// Set (numerical) dimensionless numbers
const T cahnNr = 0.114;                             // Cahn Number = xi / radius;
const T pecletNr = 0.43;                            // Pe = extensionRate * radius * xi / ( M A )

// Set time interval to simulate and normalized output
const T maxPhysTime = 160.;                         // normalized time normT = T * extensionRate [either in lu or pu]
const T physVtkIter = 1.;                           // normalized vtk output (i.e. per normalized time step, one output!)
const T physStatIter = 1.;                          // normalized vtk output (i.e. per normalized time step, one output!)
const T physInterval = 1.;                          // interval for the convergence check in s
const T residuum = 7e-5;                            // residuum for the convergence check
const T defResiduum = 0.001;						// deformation residuum criterion

// Set contact angles
// no boundary touching here...
// h_i= Parameter related to resulting contact angle of the boundary. [lattice units]
const T h1 = 0.;                                    // Contact angle 90 degrees   [lattice units]
const T h2 = 0.;                                    // Contact angle 90 degrees   [lattice units]

// Center of rotating cylinders
const T cylinderOffsetMinus = (nx/2. - cylinderOffset);
const T cylinderOffsetPlus = (nx/2. + cylinderOffset);
std::vector<T> center1 = {cylinderOffsetMinus, cylinderOffsetMinus};
std::vector<T> center2 = {cylinderOffsetPlus, cylinderOffsetMinus};
std::vector<T> center3 = {cylinderOffsetPlus, cylinderOffsetPlus};
std::vector<T> center4 = {cylinderOffsetMinus, cylinderOffsetPlus};

// Initialize remaining parameters (conversion, computation, free energy)
T strainRate;       // [pu] physical shear rate [1/s]
T rotVelo;          // [pu] physical shear speed corr. to dropRe [m/s]

T surfTen;          // [pu] physical surface tension
T surfTenLatt;      // [lu]

T kappa;            // [lu] this is kappa1+kappa2, to simplify, and meet komrakova model. this is NOT kappa in Komrakova2013!!
T kappa1;           // [lu]
T kappa2;           // [lu]
T alpha;            // [lu] free energy alpha in surface tension
T gama;             // [lu] gama "Diffusivity of the interface" same as in Komrakova2013

T xiThickness;      // [lu] characteristical thickness of interface (xi from Komrakova2013)
T kappaKomra;       // [lu] this is kappa from Komrakova2013
T aKomra;           // [lu] this is A from Komrakova2013
T transportCoeff;   // [lu] this is M from Komrakova2013

int maxIter;        // [lu] maximum lattice time
int vtkIter;        // [lu] lattice time step, when output happens
int statIter;       // [lu] lattice timestep when statistics apper in terminal


// Functors for rotating cylinders;
// cylinders rotate with tangential speed rotVel[lu]
template <typename T>
class RotatingMill2D : public AnalyticalF2D<T,T> {
  public:
    RotatingMill2D(std::vector<T> center, T radius, T rotVel, bool counterClock);
    bool operator()(T output[], const T input[]) override;
  private:
    T _center[2];
    T _radius;
    T _rotVel;
    bool _counterClock;
};

template <typename T>
RotatingMill2D<T>::RotatingMill2D(std::vector<T> center, T radius, T rotVel, bool counterClock) : AnalyticalF2D<T,T>(2){
  this->_center[0] = center[0];
  this->_center[1] = center[1];
  this->_radius = radius;
  this->_rotVel = rotVel;
  this->_counterClock = counterClock;
}

template <typename T>
bool RotatingMill2D<T>::operator()(T output[], const T input[]){
  output[0] = -(input[1]-_center[1])*_rotVel/_radius;
  output[1] =  (input[0]-_center[0])*_rotVel/_radius;
  if(!_counterClock){
    output[0] =  (input[1]-_center[1])*_rotVel/_radius;
    output[1] = -(input[0]-_center[0])*_rotVel/_radius;
  }
  return true;
}

// Measurement of droplet: interpolation in x-direction
// Simple interpolation on y = ny/2 with assumption of droplet center in domain center
T horizontalMeasure(SuperLattice<T, DESCRIPTOR>& sLattice2, int points){
  SuperLattice<T, DESCRIPTOR> *sLatticeTest = &sLattice2;;
  SuperLatticeDensity2D<T, DESCRIPTOR> density2( *sLatticeTest );

  AnalyticalFfromSuperF2D<T> aDensityH(density2, true, 1);

  T dropletLength = 0.;
  T coordinates[2] = {T(0),nx/2.};	//interpolation coordinates, y fixed

  for(int i=0; i<=points-1; i++){
    coordinates[0] = nx/2. + ((nx/2.) * (T) i/ ((T) (points-1.)));
    T output[1];
    aDensityH( output, coordinates );

    if(output[0] >= 0.0){			//fluid equals continuous phase
      dropletLength = coordinates[0] - nx/2.;
      break;
    }
  }
  return dropletLength;
}

// Measurement of droplet: interpolation in y-direction
// Simple interpolation on x = nx/2 with assumption of droplet center in domain center
T verticalMeasure(SuperLattice<T, DESCRIPTOR>& sLattice2, int points){
  SuperLattice<T, DESCRIPTOR> *sLatticeTest = &sLattice2;;
  SuperLatticeDensity2D<T, DESCRIPTOR> density2( *sLatticeTest );

  AnalyticalFfromSuperF2D<T> aDensityV(density2, true, 1);

  T dropletLength = 0.;
  T coordinates[2] = {nx/2.,T(0)};	//interpolation coordinates, x fixed

  for(int i=0; i<=points-1; i++){
    coordinates[1] = nx/2. + ((nx/2.) * (T) i/ ((T) (points-1.)));
    T output[1];
    aDensityV( output, coordinates );

    if(output[0] >= 0.0){			//fluid equals continuous phase
      dropletLength = coordinates[1] - nx/2.;
      break;
    }
  }
  return dropletLength;
}

SuperGeometry<T,2> prepareGeometry(UnitConverter<T,DESCRIPTOR> const& converter)
{
  OstreamManager clout( std::cout,"prepareGeometry" );
  clout << "Prepare Geometry ..." << std::endl;

  std::vector<T> extend = { nx, nx };
  std::vector<T> origin = { 0, 0 };
  std::shared_ptr<IndicatorF2D<T>> cuboid = std::make_shared<IndicatorCuboid2D<T>>( extend, origin );
#ifdef PARALLEL_MODE_MPI
  const int noOfCuboids = singleton::mpi().getSize();
#else
  const int noOfCuboids = 1;
#endif
  std::shared_ptr<IndicatorF2D<T>> cylind1 = std::make_shared<IndicatorCircle2D<T>>( center1, radiusCylinder );
  std::shared_ptr<IndicatorF2D<T>> cylind2 = std::make_shared<IndicatorCircle2D<T>>( center2, radiusCylinder );
  std::shared_ptr<IndicatorF2D<T>> cylind3 = std::make_shared<IndicatorCircle2D<T>>( center3, radiusCylinder );
  std::shared_ptr<IndicatorF2D<T>> cylind4 = std::make_shared<IndicatorCircle2D<T>>( center4, radiusCylinder );
  IndicatorCircle2D<T> cylinder1( center1, radiusCylinder );
  IndicatorCircle2D<T> cylinder2( center2, radiusCylinder );
  IndicatorCircle2D<T> cylinder3( center3, radiusCylinder );
  IndicatorCircle2D<T> cylinder4( center4, radiusCylinder );
  CuboidGeometry2D<T>* cGeometry = new CuboidGeometry2D<T>( *(cuboid-(cylind1+cylind2+cylind3+cylind4)),
                                converter.getPhysDeltaX(), noOfCuboids );
  cGeometry->setPeriodicity( true, true );
  HeuristicLoadBalancer<T>* loadBalancer = new HeuristicLoadBalancer<T>( *cGeometry );
  SuperGeometry<T,2> superGeometry( *cGeometry,*loadBalancer );

  superGeometry.rename( 0,2 );
  superGeometry.rename( 2,1,{0,0} );
  superGeometry.rename( 1,3,cylinder1 );	//Material number 3: cylinder bottom left
  superGeometry.rename( 1,4,cylinder2 );	//Material number 4: cylinder bottom right
  superGeometry.rename( 1,5,cylinder3 );	//Material number 5: cylinder top right
  superGeometry.rename( 1,6,cylinder4 );	//Material number 6: cylinder top left

  // clean up
  superGeometry.clean();
  superGeometry.innerClean();
  superGeometry.checkForErrors();
  superGeometry.getStatistics().print();

  clout << "Prepare Geometry ... OK" << std::endl;
  return superGeometry;
}

void prepareLattice( SuperLattice<T, DESCRIPTOR>& sLattice1,
                     SuperLattice<T, DESCRIPTOR>& sLattice2,
                     UnitConverter<T, DESCRIPTOR>& converter,
                     SuperGeometry<T,2>& superGeometry)
{
  OstreamManager clout( std::cout,"prepareLattice" );
  clout << "Prepare Lattice ..." << std::endl;

  T omega = converter.getLatticeRelaxationFrequency();

  // define lattice Dynamics
  sLattice1.defineDynamics<NoDynamics>( superGeometry, 0 );
  sLattice2.defineDynamics<NoDynamics>( superGeometry, 0 );

  sLattice1.defineDynamics<ForcedBGKdynamics>( superGeometry, 1 );
  sLattice2.defineDynamics<FreeEnergyBGKdynamics>( superGeometry, 1 );

  sLattice1.defineDynamics<BounceBack>( superGeometry, 2 );
  sLattice2.defineDynamics<BounceBack>( superGeometry, 2 );

  sLattice1.defineDynamics<NoDynamics>( superGeometry, 3 );
  sLattice2.defineDynamics<NoDynamics>( superGeometry, 3 );

  auto cylinderIndicator1 = superGeometry.getMaterialIndicator(3);
  setFreeEnergyInletBoundary<T, DESCRIPTOR>(sLattice1, omega, cylinderIndicator1, "velocity", 1);
  setFreeEnergyInletBoundary<T, DESCRIPTOR>(sLattice2, omega, cylinderIndicator1, "velocity", 2);
  auto cylinderIndicator2 = superGeometry.getMaterialIndicator(4);
  setFreeEnergyInletBoundary<T, DESCRIPTOR>(sLattice1, omega, cylinderIndicator2, "velocity", 1);
  setFreeEnergyInletBoundary<T, DESCRIPTOR>(sLattice2, omega, cylinderIndicator2, "velocity", 2);
  auto cylinderIndicator3 = superGeometry.getMaterialIndicator(5);
  setFreeEnergyInletBoundary<T, DESCRIPTOR>(sLattice1, omega, cylinderIndicator3, "velocity", 1);
  setFreeEnergyInletBoundary<T, DESCRIPTOR>(sLattice2, omega, cylinderIndicator3, "velocity", 2);
  auto cylinderIndicator4 = superGeometry.getMaterialIndicator(6);
  setFreeEnergyInletBoundary<T, DESCRIPTOR>(sLattice1, omega, cylinderIndicator4, "velocity", 1);
  setFreeEnergyInletBoundary<T, DESCRIPTOR>(sLattice2, omega, cylinderIndicator4, "velocity", 2);

  sLattice1.setParameter<OMEGA>(omega);
  sLattice2.setParameter<OMEGA>(omega);
  sLattice2.setParameter<collision::FreeEnergy::GAMMA>(gama);

  // bulk initial conditions
  // define circular domain for fluid 2
  std::vector<T> v( 2,T() );
  AnalyticalConst2D<T,T> zeroVelocity( v );

  AnalyticalConst2D<T,T> zero( 0. );
  AnalyticalConst2D<T,T> one ( 1. );

  SmoothIndicatorCircle2D<T,T> circle( {nx/2., nx/2.}, radius, converter.getPhysLength(alpha) );
  AnalyticalIdentity2D<T,T> rho( one );
  AnalyticalIdentity2D<T,T> phi( one - circle - circle );

  RotatingMill2D<T> rotator1 ( center1, radiusCylinder, converter.getCharLatticeVelocity(), true);
  RotatingMill2D<T> rotator2 ( center2, radiusCylinder, converter.getCharLatticeVelocity(), false);
  RotatingMill2D<T> rotator3 ( center3, radiusCylinder, converter.getCharLatticeVelocity(), true);
  RotatingMill2D<T> rotator4 ( center4, radiusCylinder, converter.getCharLatticeVelocity(), false);

  sLattice1.defineRhoU( cylinderIndicator1, rho, rotator1 );
  sLattice2.defineRho ( cylinderIndicator1, phi );
  sLattice1.defineRhoU( cylinderIndicator2, rho, rotator2 );
  sLattice2.defineRho ( cylinderIndicator2, phi );
  sLattice1.defineRhoU( cylinderIndicator3, rho, rotator3 );
  sLattice2.defineRho ( cylinderIndicator3, phi );
  sLattice1.defineRhoU( cylinderIndicator4, rho, rotator4 );
  sLattice2.defineRho ( cylinderIndicator4, phi );

  // equilibrium population initialization
  sLattice1.iniEquilibrium( superGeometry, 1, rho, zeroVelocity );
  sLattice2.iniEquilibrium( superGeometry, 1, phi, zeroVelocity );
  sLattice1.iniEquilibrium( superGeometry, 3, rho, rotator1 );
  sLattice2.iniEquilibrium( superGeometry, 3, phi, rotator1 );
  sLattice1.iniEquilibrium( superGeometry, 4, rho, rotator2 );
  sLattice2.iniEquilibrium( superGeometry, 4, phi, rotator2 );
  sLattice1.iniEquilibrium( superGeometry, 5, rho, rotator3 );
  sLattice2.iniEquilibrium( superGeometry, 5, phi, rotator3 );
  sLattice1.iniEquilibrium( superGeometry, 6, rho, rotator4 );
  sLattice2.iniEquilibrium( superGeometry, 6, phi, rotator4 );


  sLattice1.initialize();
  sLattice2.initialize();

  sLattice1.communicate();
  sLattice2.communicate();

  clout << "Prepare Lattice ... OK" << std::endl;
}

void prepareCoupling(SuperLattice<T, DESCRIPTOR>& sLattice1,
                     SuperLattice<T, DESCRIPTOR>& sLattice2) {

  OstreamManager clout( std::cout,"prepareCoupling" );
  clout << "Add lattice coupling" << std::endl;

  // Add the lattice couplings
  // The chemical potential coupling must come before the force coupling
  FreeEnergyChemicalPotentialGenerator2D<T, DESCRIPTOR> coupling1(
    alpha, kappa1, kappa2);
  FreeEnergyForceGenerator2D<T, DESCRIPTOR> coupling2;

  sLattice1.addLatticeCoupling( coupling1, sLattice2 );
  sLattice2.addLatticeCoupling( coupling2, sLattice1 );

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

void getResults( SuperLattice<T, DESCRIPTOR>& sLattice2,
                 SuperLattice<T, DESCRIPTOR>& sLattice1, int iT,
                 SuperGeometry<T,2>& superGeometry, util::Timer<T>& timer,
                 UnitConverter<T, DESCRIPTOR> converter) {

  OstreamManager clout( std::cout,"getResults" );
  SuperVTMwriter2D<T> vtmWriter( "fourRollMill2D" );

  if ( iT==0 ) {
    // Writes the geometry, cuboid no. and rank no. as vti file for visualization
    SuperLatticeGeometry2D<T, DESCRIPTOR> geometry( sLattice1, superGeometry );
    SuperLatticeCuboid2D<T, DESCRIPTOR> cuboid( sLattice1 );
    SuperLatticeRank2D<T, DESCRIPTOR> rank( sLattice1 );
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
   if ( iT%statIter==0 ) {
    AnalyticalConst2D<T,T> half_( 0.5 );
    SuperLatticeFfromAnalyticalF2D<T, DESCRIPTOR> half(half_, sLattice1);

    SuperLatticeDensity2D<T, DESCRIPTOR> density1( sLattice1 );
    density1.getName() = "rho";
    SuperLatticeDensity2D<T, DESCRIPTOR> density2( sLattice2 );
    density2.getName() = "phi";

    SuperIdentity2D<T,T> c1 (half*(density1+density2));
    c1.getName() = "density-fluid-1";
    SuperIdentity2D<T,T> c2 (half*(density1-density2));
    c2.getName() = "density-fluid-2";

    // velocity (seperately and combined)
    SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity1( sLattice1, converter );
    SuperLatticePhysVelocity2D<T, DESCRIPTOR> velocity2( sLattice2, converter );
    SuperIdentity2D<T,T> velocityMittel (half*(velocity1+velocity2));
    velocityMittel.getName() = "velocity-Mittel";

    vtmWriter.addFunctor( density1 );
    vtmWriter.addFunctor( density2 );
    vtmWriter.addFunctor( c1 );
    vtmWriter.addFunctor( c2 );
    vtmWriter.addFunctor( velocityMittel );
    vtmWriter.write( iT );
  }
}

int main( int argc, char *argv[] )
{
  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );

  // use remaining equations from Komrakova2013 to calculate missing parameters:
  // 0. compute velocity from Re (first strainRate, then velocity)
  // formula: https://doi.org/10.1016/S0377-0257(01)00123-9
  // (do that before converter, since it is required as argument!)
  strainRate = contVisc * dropRe * 0.35/ util::pow(radius, 2);
  // apparatus specific constant (~0.35) is only for this simulations setup with this cylinderRadius and offset
  // has to be calculated empirically: https://doi.org/10.1098/rspa.1934.0169
  rotVelo = strainRate / 0.35;

  UnitConverterFromResolutionAndRelaxationTime<T,DESCRIPTOR> converter(
    (T)   N,                     // resolution
    (T)   latticeRelaxationTime, // lattice relaxation time
    (T)   util::pow(radius, 2.), // charPhysLength: reference length of simulation geometry (dropRe!!)
    (T)   rotVelo,               // charPhysVelocity: maximal/highest expected velocity during simulation in __m / s__
    (T)   contVisc,              // physViscosity: physical kinematic viscosity in __m^2 / s__
    (T)   physDensity            // physDensity: physical density in __kg / m^3__
  );

  T deborahNumber = strainRate * converter.getLatticeRelaxationTime();

  // 1. compute surface tension from Ca
  surfTen = radius * strainRate * physDensity * contVisc / capillaryNr;

  surfTenLatt = surfTen / (converter.getConversionFactorPressure() * converter.getConversionFactorLength());

  // 2. compute interface thickness from Ch
  xiThickness = cahnNr * converter.getLatticeLength(radius);

  // 3. compute kappa from xi (Komrakova2013 equation (12))
  kappaKomra = (3./4.) * surfTenLatt * xiThickness;

  // 4. compute A from kappa (Komrakova2013 equation (11), neglect minus sign, since solely multiplications here)
  aKomra = 2. * kappaKomra / util::pow(xiThickness, 2);

  // 5. compute M from Pe (Komrakova2013 equation (17), there's a typo in the paper, should be a minus on the right hand side... we dropped that anyway.)
  transportCoeff = strainRate * converter.getConversionFactorTime() * converter.getLatticeLength(radius) * xiThickness / (pecletNr * aKomra);

  // 6. compute gamma from M (Komrakova2013 equation (8), caution: \Delta t is unity in lattice units!)
  gama = transportCoeff / (converter.getLatticeRelaxationTime() - .5);

  // Parameter fit from Komrakova to OpenLB (Semprebon2016):
  kappa1 = 4. * aKomra;
  kappa2 = kappa1;
  alpha = xiThickness / 2.;

  // transform normalized time to physical and then to lattice:
  maxIter = converter.getLatticeTime( maxPhysTime );
  vtkIter = converter.getLatticeTime( physVtkIter );
  statIter = converter.getLatticeTime( physStatIter );;

  // Prints the converter log as console output
  converter.print();

  // now print everything additionally to the unit converter (helps for double-checking)
  clout << "----------------------------------------------------------------------" << std::endl;
  clout << "----------------------------------------------------------------------" << std::endl;
  clout << "strain rate           " << strainRate     << "[pu], " << converter.getLatticeVelocity(strainRate) << "[lu]" << std::endl;
  clout << "Re_drop               " << dropRe         << std::endl;
  clout << "Deborah               " << deborahNumber  << std::endl;
  clout << "Ca                    " << capillaryNr    << std::endl;
  clout << "radius                " << radius         << "[pu], " << converter.getLatticeLength(radius) << "[lu]" << std::endl;
  clout << "surface tension       " << surfTen        << "[pu], " << surfTenLatt << "[lu]" << std::endl;
  clout << "Cahn Nr.              " << cahnNr         << std::endl;
  clout << "xi                    " << xiThickness    << "[lu]" << std::endl;
  clout << "kappa1                " << kappa1         << "[lu]" << std::endl;
  clout << "kappa2                " << kappa2         << "[lu]" << std::endl;
  clout << "Peclet Nr.            " << pecletNr       << std::endl;
  clout << "A                     " << aKomra         << "[lu]" << std::endl;
  clout << "Mobility M            " << transportCoeff << "[lu]" << std::endl;
  clout << "Gamma                 " << gama           << "[lu]" << std::endl;
  clout << "cylinderRadius        " << radiusCylinder << "[pu], " << converter.getLatticeLength(radiusCylinder) << "[lu]" << std::endl;
  clout << "vtkIter               " << vtkIter        << "[lu]" << std::endl;
  clout << "statIter              " << statIter       << "[lu]" << std::endl;
  clout << "maxIter               " << maxIter        << "[lu]" << std::endl;
  clout << "----------------------------------------------------------------------" << std::endl;
  clout << "----------------------------------------------------------------------" << std::endl;

  // === 2nd Step: Prepare Geometry ===
  // Instantiation of superGeometry
  SuperGeometry<T,2> superGeometry( prepareGeometry( converter ) );

  // === 3rd Step: Prepare Lattice ===ca
  SuperLattice<T, DESCRIPTOR> sLattice1( superGeometry );
  SuperLattice<T, DESCRIPTOR> sLattice2( superGeometry );

  prepareLattice( sLattice1, sLattice2, converter, superGeometry );

  prepareCoupling( sLattice1, sLattice2 );

  // === 4th Step: Main Loop with Timer ===
  int iT = 0;
  clout << "starting simulation..." << std::endl;
  util::ValueTracer<T> converge( converter.getLatticeTime( physInterval ), residuum );
  util::Timer<T> timer( maxIter, superGeometry.getStatistics().getNvoxel() );
  timer.start();

  T old_deformation = 0;
  T deformation = 0;

  for ( iT=0; iT<=maxIter; ++iT ) {
    // doplet measurement
  if ( iT%statIter==0 ) {
      clout<<"--------------- t: "<< converter.getPhysTime(iT) <<" ----------------"<<std::endl;
    old_deformation = deformation;
      T lengthHorizontal = horizontalMeasure( sLattice2, N * 1000);
      T lengthVertical   = verticalMeasure(   sLattice2, N * 1000);
    deformation = (lengthHorizontal - lengthVertical) / (lengthHorizontal + lengthVertical);
      clout << "Length(horizontal): " << lengthHorizontal << std::endl;
      clout << "Length(vertical):   " << lengthVertical   << std::endl;
      clout << "Deformation:        " << deformation << std::endl;
      clout << "Length/radius       " << lengthHorizontal / radius << std::endl;
  }
    // Convergence check
    if ( converge.hasConverged() && util::abs(deformation - old_deformation) < defResiduum ) {
      clout << "Simulation converged." << std::endl;
      getResults( sLattice2, sLattice1, iT, superGeometry, timer, converter );
      break;
    }
    // Computation and output of the results
    getResults( sLattice2, sLattice1, iT, superGeometry, timer, converter );

    // Collide and stream execution
    sLattice1.collideAndStream();
    sLattice2.collideAndStream();

    // MPI communication for lattice data
    sLattice1.communicate();
    sLattice2.communicate();

    // Execute coupling between the two lattices
    sLattice1.executeCoupling();
    sLattice2.executeCoupling();

    converge.takeValue( sLattice2.getStatistics().getAverageRho(), true );
  }

  timer.stop();
  timer.printSummary();
}
