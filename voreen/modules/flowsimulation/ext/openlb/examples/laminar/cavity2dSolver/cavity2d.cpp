/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2006 - 2012 Mathias J. Krause, Jonas Fietz,
 *  Jonas Latt, Jonas Kratzke
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

/* cavity2d.cpp:
 * This example illustrates a flow in a cuboid, lid-driven cavity.
 * It is similar to the other cavity2d example and illustrates the
 * application a solver class which encapsulates large parts of the simulation.
 */

#include "olb2D.h"
#include "olb2D.hh"

using namespace olb;
using namespace olb::names;

using Descriptor = descriptors::D2Q9<>;
using Lattices = meta::map<
  NavierStokes, Descriptor
>;
template<typename T>
using Parameters = meta::map<
  Simulation, parameters::XmlSimulation<T,Lattices>,
  Stationarity, parameters::Stationarity<T>,
  Output, parameters::OutputBase<T>
>;

template<typename T>
class Cavity2dSolver : public LBSolver<T,Parameters<T>,Lattices> {
private:
  mutable OstreamManager            clout {std::cout, "Cavity2dSolver"};

public:
  Cavity2dSolver(utilities::TypeIndexedTuple<Parameters<T>> params)
   : Cavity2dSolver::LBSolver(params)
  { }

protected:
  void prepareGeometry() override
  {
    Vector<T,2> extend( 1,1 );
    Vector<T,2> origin( 0,0 );
    IndicatorCuboid2D<T> cuboid( extend, origin );

    this->_cGeometry = std::make_shared<CuboidGeometry2D<T>> (
      cuboid, 
      this->converter().getConversionFactorLength(),
      this->parameters(Simulation()).noC * singleton::mpi().getSize() );

    this->_cGeometry->print();

    this->_loadBalancer = std::make_shared<HeuristicLoadBalancer<T>> (
      *this->_cGeometry);
    this->_sGeometry = std::make_shared<SuperGeometry<T,2>> (
      *this->_cGeometry,
      *this->_loadBalancer, 
      this->parameters(Simulation()).overlap);

    this->geometry().rename( 0,2 );
    this->geometry().rename( 2,1,{1,1} );

    T eps = this->converter().getConversionFactorLength();
    Vector<T,2> extendShifted( T( 1 ) + 2*eps, 2*eps );
    Vector<T,2> originShifted( T() - eps, T( 1 ) - eps );
    IndicatorCuboid2D<T> lid( extendShifted, originShifted );
    this->geometry().rename( 2,3,1,lid );
  }

  void prepareLattices() override
  {
    // link lattice with dynamics for collision step
    // Material=0 -->do nothing
    this->lattice().template defineDynamics<NoDynamics<T,Descriptor>>(this->geometry(), 0);

    // Material=1 -->bulk dynamics
    // Material=2,3 -->bulk dynamics, velocity boundary
    using BulkDynamics = ConstRhoBGKdynamics<T,Descriptor>;
    this->lattice().template defineDynamics<BulkDynamics>(
      this->geometry().getMaterialIndicator({1,2,3}));

    const T omega = this->converter().getLatticeRelaxationFrequency();
    setInterpolatedVelocityBoundary<T,Descriptor,BulkDynamics>(
      this->lattice(), omega, this->geometry(), 2);
    setInterpolatedVelocityBoundary<T,Descriptor,BulkDynamics>(
      this->lattice(), omega, this->geometry(), 3);

    this->lattice().template setParameter<descriptors::OMEGA>(omega);
  }

  void setInitialValues() override
  {
    // set initial values: v = [0,0]
    AnalyticalConst2D<T,T> rhoF( 1 );
    std::vector<T> velocity( 2,T() );
    AnalyticalConst2D<T,T> uF( velocity );

    const auto& bulkIndicator = this->geometry().getMaterialIndicator({1, 2, 3});
    this->lattice().iniEquilibrium( bulkIndicator, rhoF, uF );
    this->lattice().defineRhoU( bulkIndicator, rhoF, uF );

    // set non-zero velocity for upper boundary cells
    velocity[0] = this->converter().getCharLatticeVelocity();
    AnalyticalConst2D<T,T> u( velocity );
    this->lattice().defineU( this->geometry(), 3, u );
  }

  void setBoundaryValues(std::size_t iT) override { }

  void prepareVTK() const override {
    SuperVTMwriter2D<T> vtmWriter(this->parameters(Output()).vtkFilename);
    auto sLattice = &this->lattice();

    SuperLatticeGeometry<T,Descriptor> geometry(*sLattice, this->geometry());
    SuperLatticeCuboid<T,Descriptor> cuboid(*sLattice);
    SuperLatticeRank<T,Descriptor> rank(*sLattice);
    SuperLatticeDiscreteNormal2D discreteNormal(
      *sLattice, this->geometry(), this->geometry().getMaterialIndicator({2, 3}) );
    SuperLatticeDiscreteNormalType2D discreteNormalType(
      *sLattice, this->geometry(), this->geometry().getMaterialIndicator({2, 3}) );

    vtmWriter.write( cuboid );
    vtmWriter.write( geometry );
    vtmWriter.write( rank );
    vtmWriter.write( discreteNormal );
    vtmWriter.write( discreteNormalType );

    vtmWriter.createMasterFile();
  }

  void writeImages(std::size_t iT) const override {
    SuperLatticePhysVelocity2D velocity( this->lattice(), this->converter() );
    SuperEuklidNorm2D normVel( velocity );
    BlockReduction2D2D<T> planeReduction( normVel, 600, BlockDataSyncMode::ReduceOnly );
    // write output of velocity as JPEG
    heatmap::write(planeReduction, iT);
  }

  void writeVTK(std::size_t iT) const override {
    SuperVTMwriter2D<T> vtmWriter( this->parameters(Output()).vtkFilename );

    SuperLatticePhysVelocity2D velocity( this->lattice(), this->converter() );
    SuperLatticePhysPressure2D pressure( this->lattice(), this->converter() );
    
    vtmWriter.addFunctor( velocity );
    vtmWriter.addFunctor( pressure );
    vtmWriter.write( iT );
  }

  void writeGnuplot() const override {
    // Gives access to velocity information on lattice
    SuperLatticePhysVelocity2D velocityField( this->lattice(), this->converter() );
    // Interpolation functor with velocityField information
    AnalyticalFfromSuperF2D<T> interpolation( velocityField, true, 1 );

    Vector<int,17> y_coord( {128, 125, 124, 123, 122, 109, 94, 79, 64, 58, 36, 22, 13, 9, 8, 7, 0} );
    // Ghia, Ghia and Shin, 1982: "High-Re Solutions for Incompressible Flow 
    // Using the Navier-Stokes Equations and a Multigrid Method";  Table 1
    Vector<T,17> vel_ghia_RE1000( { 1.0,     0.65928, 0.57492, 0.51117, 0.46604,
                                    0.33304, 0.18719, 0.05702,-0.06080,-0.10648,
                                    -0.27805,-0.38289,-0.29730,-0.22220,-0.20196,
                                    -0.18109, 0.0
                                  } );
    Vector<T,17> vel_ghia_RE100( {1.0,     0.84123, 0.78871, 0.73722, 0.68717,
                                  0.23151, 0.00332,-0.13641,-0.20581,-0.21090,
                                  -0.15662,-0.10150,-0.06434,-0.04775,-0.04192,
                                  -0.03717, 0.0
                                } );
    Vector<T,17> vel_simulation;

    // Gnuplot interface to create plots
    static Gnuplot<T> gplot( "centerVelocityX" );
    // Define comparison values
    Vector<T,17> comparison = vel_ghia_RE1000;

    for ( int nY = 0; nY < 17; ++nY ) {
      // 17 data points evenly distributed between 0 and 1 (height)
      T position[2] = {0.5, y_coord[nY]/128.0};
      T velocity[2] = {T(), T()};
      // Interpolate velocityField at "position" and save it in "velocity"
      interpolation( velocity, position );
      // Save value of velocity (in x-direction) in "vel_simulation" for every position "nY"
      vel_simulation[nY] = velocity[0];
      // Set data for plot output
      gplot.setData( position[1], {vel_simulation[nY],comparison[nY]}, {"simulated","Ghia"} );
    }
    // Create PNG file
    gplot.writePNG();
    // Console output with results
    clout
      << "absoluteErrorL2(line)=" << norm(vel_simulation - comparison) / 17. 
      << "; relativeErrorL2(line)=" << norm(vel_simulation - comparison) / norm(comparison)
      << std::endl;
  }
};


int main( int argc, char* argv[] )
{
  using T = double;

  olbInit( &argc, &argv );

  XMLreader config( "cavity2d.xml" );
  auto cavity2d = createLBSolver<Cavity2dSolver<T>> (config);

  cavity2d->solve();
}
