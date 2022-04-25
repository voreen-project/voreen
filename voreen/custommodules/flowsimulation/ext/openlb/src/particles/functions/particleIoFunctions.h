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


/* This file contains functions for console and vtk output.
 * Those include spefic and generic version.
*/

#ifndef PARTICLE_IO_FUNCTIONS_H
#define PARTICLE_IO_FUNCTIONS_H


namespace olb {

namespace particles {

namespace io {

/// Generic vector print method
template<typename VECTOR>
std::string strFormatVector(VECTOR vector, int strWidth=13)
{
  const int dim=vector.getDim();
  std::stringstream stream;
  stream << "(";
  for( int iDim=0; iDim<dim-1; ++iDim ){  
    stream << std::setw(strWidth) << vector[iDim] << ", ";
  }
  stream << std::setw(strWidth) << vector[dim-1] << ")";
  return stream.str();
}


/// Generic printing method for particles, which automatically traverses through all provided fields
template<typename T, typename PARTICLETYPE, bool multiOutput=false>
void printGenericParticleInfo( DynamicFieldGroupsD<T, typename PARTICLETYPE::fieldList>& dynamicFieldGroups, std::size_t iP )
{
  OstreamManager clout( std::cout, "Particle "+std::to_string(iP) );
  clout.setMultiOutput(multiOutput);
  clout << "================================================" << std::endl;
  //Define output lambda expression
  typedef std::function<bool(const std::type_info&,int,std::string)> FunctionType;
  FunctionType printFunction = [](const std::type_info& typeInfo, int fieldSize, std::string fieldContentStr) {
    OstreamManager clout( std::cout,"Field" );
    clout << std::setw(45) << std::string(typeInfo.name()) << " (" << fieldSize << ") : " << fieldContentStr << std::endl;
    return false; //resetField=false
  };
  //Call recursive field traversal function with lambda expression
  descriptors::access_field_content<FunctionType,T,PARTICLETYPE,typename PARTICLETYPE::fieldList>::fieldsL2(
    printFunction, dynamicFieldGroups, iP );
  clout << "================================================" << std::endl;
  clout.setMultiOutput(false);
}


/// Simple printing method for resolved particles
template<typename T, typename PARTICLETYPE>
void printResolvedParticleInfoSimple( Particle<T, PARTICLETYPE>& particle )
{
  using namespace descriptors;
  constexpr unsigned D = PARTICLETYPE::d;
  const unsigned Drot = utilities::dimensions::convert<D>::rotation;
  //Retrieve Data
  Vector<T,D> position( particle.template getField<GENERAL,POSITION>() );
  Vector<T,D> velocity( particle.template getField<MOBILITY,VELOCITY>() );
  Vector<T,Drot> angularVelocity( particle.template getField<MOBILITY,ANG_VELOCITY>() );
  Vector<T,Drot> angle( particle.template getField<SURFACE,ANGLE>() );
  auto indicatorPtr = particle.template getField<SURFACE,SINDICATOR>();
  //Print data
  OstreamManager clout( std::cout,"Particle" );
  clout << "================================================" << std::endl;
  clout << "Position(m)=           " << std::setw(13) << position << std::endl;
  clout << "Velocity(m/s)=         " << std::setw(13) << velocity << std::endl;
  clout << "Ang. Velocity(m/s)=    " << std::setw(13) << angularVelocity << std::endl;
  clout << "Angle(rad)=            " << std::setw(13) << angle << std::endl;
  clout << "Indi Pos=              " << std::setw(13) << indicatorPtr->getPos() << std::endl;
  clout << "Indi Min=              " << std::setw(13) << indicatorPtr->getMin() << std::endl;
  clout << "Indi Max=              " << std::setw(13) << indicatorPtr->getMax() << std::endl;
  clout << "================================================" << std::endl;
}



/// Printing method adapted to the original output
template<typename T, typename PARTICLETYPE>
void printResolvedParticleInfo( Particle<T, PARTICLETYPE>& particle, const std::string& streamName="ParticleInfo" )
{

  //Set up environment
  using namespace descriptors;
  constexpr unsigned D = PARTICLETYPE::d;
  const unsigned Drot = utilities::dimensions::convert<D>::rotation;

  //Retrieve data
  auto indicatorPtr = particle.template getField<SURFACE,SINDICATOR>();
  T sIndiCircumRadius = indicatorPtr->getCircumRadius();
  T mass = particle.template getField<PHYSPROPERTIES,MASS>();
  Vector<T,D> position( particle.template getField<GENERAL,POSITION>() );
  Vector<T,D> velocity( particle.template getField<MOBILITY,VELOCITY>() );
  Vector<T,Drot> angularVelocity( particle.template getField<MOBILITY,ANG_VELOCITY>() );
  Vector<T,Drot> angle( particle.template getField<SURFACE,ANGLE>() );
  Vector<T,D> force( particle.template getField<FORCING,FORCE>() );
  Vector<T,D> acceleration( particle.template getField<MOBILITY,ACCELERATION_STRD>() );
  Vector<T,Drot> angularAcceleration( particle.template getField<MOBILITY,ANG_ACC_STRD>() );

  //Translate angle from radian to degree
  angle               *= (180/M_PI);
  angularVelocity     *= (180/M_PI);
  angularAcceleration *= (180/M_PI);

  //Print data
  OstreamManager clout( std::cout,streamName );
  clout << "Particle " << particle.getId() << std::endl;
  clout << " |Circum radius(m)=     " << std::setw(13) << sIndiCircumRadius << std::endl;
  clout << " |Mass(kg)=             " << std::setw(13) << mass << std::endl;
  clout << " |Position(m)=         " << strFormatVector<Vector<T,D>>(position) << std::endl;
  clout << " |Angle(°)=            " << strFormatVector<Vector<T,Drot>>(angle) << std::endl;
  clout << " |Velocity(m/s)=       " << strFormatVector<Vector<T,D>>(velocity) << std::endl;
  clout << " |Ang. Velocity(°/s)=  " << strFormatVector<Vector<T,Drot>>(angularVelocity) << std::endl;
  clout << " |Force(N)=            " << strFormatVector<Vector<T,D>>(force) << std::endl;
  clout << " |Acceleration(m/s^2)= " << strFormatVector<Vector<T,D>>(acceleration) << std::endl;
  clout << " |Ang. acc.(°/s^2)=    " << strFormatVector<Vector<T,Drot>>(angularAcceleration) << std::endl;
}


} //namespace io

} //namespace particles

} //namespace olb

#endif
