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

#ifndef RANDOM_HH
#define RANDOM_HH


#include <type_traits>

namespace olb {

namespace util {

//TODO: disable constructor for <T,true> ?
template<typename T,bool useStored>
Randomizer<T,useStored>::Randomizer()
{
  static_assert(!useStored, "ERROR: no sequence provided");
}

template<typename T,bool useStored>
Randomizer<T,useStored>::Randomizer(std::vector<T> sequence)
  : _storedSequence(sequence)
{}


template<typename T,bool useStored>
Randomizer<T,useStored>::Randomizer(std::string filePathSequence)
{
  std::size_t count = 0;
  int rank = singleton::mpi().getRank();
  if (rank == 0){
    std::string strVal;
    std::ifstream file(filePathSequence);
    while (std::getline(file, strVal))
    {
      _storedSequence.push_back(std::stod(strVal));
      ++count;
    }
  }

#ifdef PARALLEL_MODE_MPI
  //Align count on all cores
  singleton::mpi().bCast<std::size_t>( &count, 1 );

  //Align vector size on all cores
  if (rank != 0){
    _storedSequence.resize(count);
  }

  //Align vector content on all cores
  singleton::mpi().bCast( _storedSequence.data(), count );
#endif
}



template<typename T,bool useStored>
template<typename O>
O Randomizer<T,useStored>::generate()
{
  if constexpr (std::is_arithmetic<O>::value) {
    O output = generateScalar();
    return output;
  } else {
    O output;
    for(unsigned iDim=0; iDim<O::d; ++iDim){
      output[iDim] = generateScalar();
    }
    return output;
  }
}

template<typename T,bool useStored>
T Randomizer<T,useStored>::generateScalarNormal(T avg, T stdDev )
{
  T A = generateScalar();
  T B = generateScalar();
  T x = util::cos(2 * M_PI * A) * util::sqrt(-2 * util::log(B));
  T output = avg + x * stdDev;
  return output;
}

template<typename T,bool useStored>
T Randomizer<T,useStored>::generateScalarNormal(T avg, T stdDev, T cutoff )
{
  T scalar = generateScalarNormal( avg, stdDev );
  if (scalar<(avg-cutoff) || scalar>(avg+cutoff)){
    scalar = generateScalarNormal(avg, stdDev, cutoff );
  }
  return scalar;
}

template<typename T,bool useStored>
void Randomizer<T,useStored>::writeSequence( std::size_t numOfValues,
                                             std::string filePathSequence,
                                             int precision )
{
  int rank = singleton::mpi().getRank();
  if (rank == 0){
    std::ofstream writer;
    writer.open( filePathSequence );
    writer << std::setprecision(precision);
    for (std::size_t i=0; i<numOfValues; ++i){
      writer << generateScalarRandom() << " " << std::endl;
    }
    writer.close();
  }
}


template<typename T,bool useStored>
T Randomizer<T,useStored>::generateScalarRandom()
{
  srand(clock());
  return (T) (rand() % 100000) / 100000. ; //precision 5
}

template<typename T,bool useStored>
T Randomizer<T,useStored>::generateScalarStored()
{
  static_assert(useStored, "ERROR: no sequence provided");
  auto idx = _count%_storedSequence.size();
  T scalar = _storedSequence[idx];
  ++_count;
  return scalar;
}

template<typename T,bool useStored>
T Randomizer<T,useStored>::generateScalar()
{
  if constexpr (useStored){
    return generateScalarStored();
  } else {
    T scalar = generateScalarRandom();
#ifdef PARALLEL_MODE_MPI
  singleton::mpi().bCast<T>( &scalar, 1 );
#endif
    return scalar;
  }
}

} //namespace util

} //namespace olb

#endif
