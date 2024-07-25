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

/* poiseuille2d.cpp:
 * This example examines a 2D Poseuille flow
 * It illustrates the computation of error norms.
 */

// the main code of the simulation is in poiseuille2d.h as it is also used by the
// example poiseuille2dEOC
#include "poiseuille2d.h"
//Initialize Gnuplot
static Gnuplot<T> gplot("centerVelocity");

int main( int argc, char* argv[] )
{
  // === 1st Step: Initialization ===
  olbInit( &argc, &argv );
  singleton::directories().setOutputDir( "./tmp/" );
  OstreamManager clout( std::cout,"main" );

  // Simulation parameter
  int N = 50;
  bool eoc = false;

  if (argc > 1) {
    if (argv[1][0]=='-'&&argv[1][1]=='h') {
      OstreamManager clout( std::cout,"help" );
      clout<<"Usage: program [Resolution] [FlowType] [BoundaryType]"<<std::endl;
      clout<<"FlowType: 0=forced, 1=nonForced"<<std::endl;
      clout<<"BoundaryType: 0=bounceBack, 1=local, " 
           <<"2=interpolated, 3=freeSlip, 4=partialSlip"<<std::endl;
      clout<<"Default: Resolution=50, FlowType=forced, "
           <<"BoundaryType=interpolated"<<std::endl;
      return 0;
    }
  }

  if (argc > 1) {
    N = atoi(argv[1]);
    if (N < 1) {
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
  
  simulatePoiseuille(N, gplot, eoc);
}
