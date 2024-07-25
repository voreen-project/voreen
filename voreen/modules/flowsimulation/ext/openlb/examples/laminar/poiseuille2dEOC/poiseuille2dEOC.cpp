/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2007, 2012 Jonas Latt, Mathias J. Krause
 *  Vojtech Cvrcek, Peter Weisbrod
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

/* poiseuille2dEOC.cpp:
 * This example examines a 2D Poseuille flow
 * It illustrates the error analysis.
 */


#include "../poiseuille2d/poiseuille2d.h"

/// Initialize gnuplot
static Gnuplot<T> gplot(
  "Velocity_and_StrainRate_eoc",
  false,
  "set terminal png size 720, 720 font 'Arial,10'",
  Gnuplot<T>::LOGLOG,
  Gnuplot<T>::LINREG);

int main( int argc, char* argv[] )
{
  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );

  /// Simulation Parameter
  int startN = 50;

  if (argc > 1) {
    if (argv[1][0]=='-'&&argv[1][1]=='h') {
    OstreamManager clout( std::cout,"help" );
    clout<<"Usage: program [Resolution] [FlowType] [BoundaryType]"<<std::endl;
    clout<<"FlowType: 0=forced, 1=nonForced"<<std::endl;
    clout<<"BoundaryType: 0=bounceBack, 1=local, 2=interpolated"<<std::endl;
    clout<<"Default: Resolution=" << startN << ", FlowType=forced, BoundaryType=interpolated"<<std::endl;
    return 0;
    }
  }

  if (argc > 1) {
    startN = atoi(argv[1]);
    if (startN < 1) {
    std::cerr << "Fluid domain is too small" << std::endl;
    return 1;
    }
  }

  if (argc > 2) {
    int flowTypeNumber = atoi(argv[2]);
    if (flowTypeNumber < 0 || flowTypeNumber > (int)nonForced) {
    std::cerr << "Unknown fluid flow type" << std::endl;
    return 2;
    }
    flowType = (FlowType) flowTypeNumber;
  }

  if (argc > 3) {
    int boundaryTypeNumber = atoi(argv[3]);
    if (boundaryTypeNumber < 0 || boundaryTypeNumber > (int) partialSlip) {
    std::cerr << "Unknown boundary type" << std::endl;
    return 3;
    }
    boundaryType = (BoundaryType) boundaryTypeNumber;
  }

  if ((boundaryType == freeSlip) || (boundaryType == partialSlip)) {
    throw std::invalid_argument(
      "eoc computation is currently not supported for slip boundary conditions");
  }

  int maxN = startN + 41;
  bool eoc = true;

  // set the labels for the plot
  gplot.setLabel("Resolution test", "average Error");
    

  // loop over the different simulations
  for(int N = startN; N < maxN; N += 10){

    /// Run the simulations
    clout << "Starting next simulation with N = " << N << std::endl;
    simulatePoiseuille(N, gplot, eoc);
  }

  gplot.writePNG();

  return 0;
}