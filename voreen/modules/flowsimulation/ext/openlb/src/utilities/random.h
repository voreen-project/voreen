/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2021 Nicolas Hafen, Mathias J. Krause
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

#ifndef RANDOM_H
#define RANDOM_H


namespace olb {

namespace util {

/* The present class is supposed to represent a random scalar generator that can
 * be fed from stored values also. This is supposed to enable a semi-random approach,
 * which will enable reproducible results with random-like values. This is primarily
 * intendet to be used in the random particle seeding.
 * When passing useStored=true, 'generateScalar' will use stored values provided in
 * the constructor or by a specified file.
 */
template<typename T,bool useStored=false>
class Randomizer{
private:
  std::vector<T> _storedSequence;
  std::size_t _count = 0;

public:
  /// Constructor for (useStored=false)
  Randomizer();

  /// Constructor with passed sequence 
  Randomizer(std::vector<T> sequence);

  /// Constroctor with filePath to stored sequence
  Randomizer(std::string filePathSequence);

  /// Generate scalar or vector filled with scalars
  template<typename O=T>
  O generate();

  /// Generate scalar leading to normal distribution based on the Box-Muller approach
  T generateScalarNormal(T avg, T stdDev );
  T generateScalarNormal(T avg, T stdDev, T cutoff );

  /// Write sequence to file for later retrieval
  void writeSequence( std::size_t numOfValues, 
                      std::string filePathSequence = "./randomSequence.dat",
                      int precision = 5 );

private:

  /// Generate scalar based on rand
  T generateScalarRandom();

  /// Generate scalar based on stored sequence
  T generateScalarStored();

  /// Generate scalar either by rand or stored sequence defined by useStored
  T generateScalar();
};


} //namespace util

} //namespace olb

#endif
