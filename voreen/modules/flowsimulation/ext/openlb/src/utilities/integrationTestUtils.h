/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Nicolas Hafen
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

#ifndef INTEGRATION_TEST_UTILS_H
#define INTEGRATION_TEST_UTILS_H

namespace olb {

namespace util {

//Test function for final values
//TODO: Adapt to parallel runs
template<typename T>
int evaluateIntegration( std::vector<T>& testValues, bool print=false ){
  OstreamManager clout( std::cout,"Integrationtest" );
  bool allPassed = true;
  for (int i=0; i<testValues.size(); i+=3){
    T val = testValues[i];
    T ref = testValues[i+1];
    T maxErr = testValues[i+2];
    T err = std::abs(ref-val)/ref;
    bool passed = err <= maxErr;
    if (print || !passed ){ 
      clout << std::setprecision(16);
      clout << "i:|"     << "Value:                 |" 
                         << "Reference:             |" 
                         << "Error:                 |" 
                         << "MaxError:              |"
                         << "Passed:" <<  std::endl;
      clout << i << " |" << std::setw(22) << val
                 << " |" << std::setw(22) << ref
                 << " |" << std::setw(22) << err 
                 << " |" << std::setw(22) << maxErr 
                 << " |" << passed << std::endl;
      clout << std::setprecision(6);
    }
    if (!passed){ allPassed=false; };
  }
  if (allPassed){
    return 0;
  } else {
    return 1;
  }
}

}//namespace util

}//namespace olb

#endif
