/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt
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

/** \file
 * Descriptor for all types of 2D and 3D lattices. In principle, thanks
 * to the fact that the OpenLB code is generic, it is sufficient to
 * write a new descriptor when a new type of lattice is to be used.
 *  -- header file
 */
#ifndef LATTICE_DESCRIPTORS_H
#define LATTICE_DESCRIPTORS_H

#include <vector>
#include "core/olbDebug.h"

namespace olb {

/// Descriptors for the 2D and 3D lattices.
/** \warning Attention: The lattice directions must always be ordered in
 * such a way that c[i] = -c[i+(q-1)/2] for i=1..(q-1)/2, and c[0] = 0 must
 * be the rest velocity. Furthermore, the velocities c[i] for i=1..(q-1)/2
 * must verify
 *  - in 2D: (c[i][0]<0) || (c[i][0]==0 && c[i][1]<0)
 *  - in 3D: (c[i][0]<0) || (c[i][0]==0 && c[i][1]<0)
 *                       || (c[i][0]==0 && c[i][1]==0 && c[i][2]<0)
 * Otherwise some of the code will work erroneously, because the
 * aformentioned relations are taken as given to enable a few
 * optimizations.
*/
namespace descriptors {

/// descriptor base
/*template <unsigned Size> class DescriptorBase
{
  enum { size = Size };       ///< number of dimensions
};

template <typename T, unsigned Size> class Vector
    : public DescriptorBase<Size>
{
    T data[Size];  ///< number of dimensions
};

template <typename T> class Scalar
    : public Vector<T,1>
{
};

template <typename T> class Vector2D
    : public Vector<T,2>
{
};

template <typename T> class Vector3D
    : public Vector<T,1>
{
};*/


struct NoExternalField {
  static const int numScalars = 0;
  static const int numSpecies = 0;
  static const int forceBeginsAt = 0;
  static const int sizeOfForce   = 0;
  static const int velocityBeginsAt = 0;
  static const int sizeOfVelocity   = 0;
};

struct NoExternalFieldBase {
  typedef NoExternalField ExternalField;
};

struct Force2dDescriptor {
  static const int numScalars = 2;
  static const int numSpecies = 1;
  static const int forceBeginsAt = 0;
  static const int sizeOfForce   = 2;
};

struct Force2dDescriptorBase {
  typedef Force2dDescriptor ExternalField;
};

struct Force3dDescriptor {
  static const int numScalars = 3;
  static const int numSpecies = 1;
  static const int forceBeginsAt = 0;
  static const int sizeOfForce   = 3;
};

struct Force3dDescriptorBase {
  typedef Force3dDescriptor ExternalField;
};

struct V6Force2dDescriptor {
  static const int numScalars = 8;
  static const int numSpecies = 2;
  static const int forceBeginsAt = 0;
  static const int sizeOfForce   = 2;
  static const int vBeginsAt = 2;
  static const int sizeOfV   = 6;
};

struct V6Force2dDescriptorBase {
  typedef V6Force2dDescriptor ExternalField;
};

struct V12Force3dDescriptor {
  static const int numScalars = 15;
  static const int numSpecies = 2;
  static const int forceBeginsAt = 0;
  static const int sizeOfForce   = 3;
  static const int vBeginsAt = 3;
  static const int sizeOfV   = 12;
};

struct V12Force3dDescriptorBase {
  typedef V12Force3dDescriptor ExternalField;
};

template<typename T, typename ExternalField>
class ExternalFieldArray {
public:
  T* get(int index)
  {
    OLB_PRECONDITION( index < ExternalField::numScalars );
    return data+index;
  }
  T const* get(int index) const
  {
    OLB_PRECONDITION( index < ExternalField::numScalars );
    return data+index;
  }
private:
  T data[ExternalField::numScalars];
};

template<typename T>
class ExternalFieldArray<T,descriptors::NoExternalField> {
public:
  T* get(unsigned index)
  {
    OLB_PRECONDITION( false );
    static T data = T();
    return &data;
  }
  T const* get(unsigned index) const
  {
    OLB_PRECONDITION( false );
    static T data = T();
    return &data;
  }
};

struct Velocity2dDescriptor {
  static const int numScalars = 2;
  static const int numSpecies = 1;
  static const int velocityBeginsAt = 0;
  static const int sizeOfVelocity   = 2;
};

struct Velocity2dBase {
  typedef Velocity2dDescriptor ExternalField;
};

struct Velocity3dDescriptor {
  static const int numScalars = 3;
  static const int numSpecies = 1;
  static const int velocityBeginsAt = 0;
  static const int sizeOfVelocity   = 3;
};

struct Velocity3dBase {
  typedef Velocity3dDescriptor ExternalField;
};

struct ParticleAdvectionDiffusion3dDescriptor {
  static const int numScalars = 6;
  static const int numSpecies = 2;
  static const int velocityBeginsAt = 0;
  static const int sizeOfVelocity   = 3;
  static const int velocity2BeginsAt = 3;
  static const int sizeOfVelocity2   = 3;
};

struct ParticleAdvectionDiffusion3dDescriptorBase {
  typedef ParticleAdvectionDiffusion3dDescriptor ExternalField;
};

/// D2Q5 lattice
template <typename T> struct D2Q5DescriptorBase {
  typedef D2Q5DescriptorBase<T> BaseDescriptor;
  enum { d = 2, q = 5 };      ///< number of dimensions/distr. functions
  static const int vicinity;  ///< size of neighborhood
  static const int c[q][d];   ///< lattice directions
  static const int opposite[q]; ///< opposite entry
  static const T t[q];        ///< lattice weights
  static const T invCs2;      ///< inverse square of speed of sound
};

/// D2Q5 lattice lattice for advection-diffusion problems (MRT)
template <typename T>
struct AdvectionDiffusionMRTD2Q5DescriptorBase {
  typedef AdvectionDiffusionMRTD2Q5DescriptorBase<T> BaseDescriptor;
  enum { d = 2, q = 5 };          // number of dimensions/distr. functions
  static const int vicinity;      // size of neighborhood
  static const int c[q][d];       // lattice directions
  static const T t[q];            // lattice weights
  static const int opposite[q];   // opposite entry
  static const T M[q][q];         // Matrix of base change between f and moments : moments=M.f
  static const T invM[q][q];      // inverse of base change matrix : f=invM.moments
  static const T S[q];            // relaxation times
  enum { shearIndexes = 2 };
  static const int shearViscIndexes[shearIndexes]; // relevant indexes of r. t. for shear viscosity
  static const int bulkViscIndex = 2; // relevant index of r. t. for bulk viscosity
  static const T invCs2;              // inverse square of speed of sound
};

template <typename T> struct AdvectionDiffusionMRTD2Q5Descriptor
    : public AdvectionDiffusionMRTD2Q5DescriptorBase<T>, public Velocity2dBase {
};


/// D2Q9 lattice
template <typename T> struct D2Q9DescriptorBase {
  typedef D2Q9DescriptorBase<T> BaseDescriptor;
  enum { d = 2, q = 9 };      ///< number of dimensions/distr. functions
  static const int vicinity;  ///< size of neighborhood
  static const int c[q][d];   ///< lattice directions
  static const int opposite[q]; ///< opposite entry
  static const T t[q];        ///< lattice weights
  static const T invCs2;      ///< inverse square of speed of sound
};

/// D3Q7 lattice
template <typename T> struct D3Q7DescriptorBase {
  typedef D3Q7DescriptorBase<T> BaseDescriptor;
  enum { d = 3, q = 7 };     ///< number of dimensions/distr. functions
  static const int vicinity;  ///< size of neighborhood
  static const int c[q][d];   ///< lattice directions
  static const int opposite[q]; ///< opposite entry
  static const T t[q];        ///< lattice weights
  static const T invCs2;      ///< inverse square of speed of sound
};

/// D3Q7 lattice for advection-diffusion problems (MRT)
template <typename T>
struct AdvectionDiffusionMRTD3Q7DescriptorBase {
  typedef AdvectionDiffusionMRTD3Q7DescriptorBase<T> BaseDescriptor;
  enum { d = 3, q = 7 };          // number of dimensions/distr. functions
  static const int vicinity;      // size of neighborhood
  static const int c[q][d];       // lattice directions
  static const T t[q];            // lattice weights
  static const int opposite[q];   // opposite entry
  static const T M[q][q];         // Matrix of base change between f and moments : moments=M.f
  static const T invM[q][q];      // inverse of base change matrix : f=invM.moments
  static const T S[q];            // relaxation times
  static const T S_2[q];          // relaxation times
  enum { shearIndexes = 3 };
  static const int shearViscIndexes[shearIndexes]; // relevant indexes of r. t. for shear viscosity
  static const int bulkViscIndex = 1; // relevant index of r. t. for bulk viscosity
  static const T invCs2;              // inverse square of speed of sound
};

template <typename T> struct AdvectionDiffusionMRTD3Q7Descriptor
    : public AdvectionDiffusionMRTD3Q7DescriptorBase<T>, public Velocity3dBase {
};

template <typename T> struct ParticleAdvectionDiffusionMRTD3Q7Descriptor
    : public D3Q7DescriptorBase<T>, public ParticleAdvectionDiffusion3dDescriptorBase {
};

/// D3Q13 lattice
template <typename T> struct D3Q13DescriptorBase {
  typedef D3Q13DescriptorBase<T> BaseDescriptor;
  enum { d = 3, q = 13 };     ///< number of dimensions/distr. functions
  static const int vicinity;  ///< size of neighborhood
  static const int c[q][d];   ///< lattice directions
  static const int opposite[q]; ///< opposite entry
  static const T t[q];        ///< lattice weights
  static const T invCs2;      ///< inverse square of speed of sound
  static const T lambda_e;    ///< relaxation parameter for the bulk stress
  static const T lambda_h;    ///< additional relaxation parameter
};

/// D3Q15 lattice
template <typename T> struct D3Q15DescriptorBase {
  typedef D3Q15DescriptorBase<T> BaseDescriptor;
  enum { d = 3, q = 15 };     ///< number of dimensions/distr. functions
  static const int vicinity;  ///< size of neighborhood
  static const int c[q][d];   ///< lattice directions
  static const int opposite[q]; ///< opposite entry
  static const T t[q];        ///< lattice weights
  static const T invCs2;      ///< inverse square of speed of sound
};

/// D3Q19 lattice
template <typename T> struct D3Q19DescriptorBase {
  typedef D3Q19DescriptorBase<T> BaseDescriptor;
  enum { d = 3, q = 19 };     ///< number of dimensions/distr. functions
  static const int vicinity;  ///< size of neighborhood
  static const int c[q][d];   ///< lattice directions
  static const int opposite[q]; ///< opposite entry
  static const T t[q];        ///< lattice weights
  static const T invCs2;      ///< inverse square of speed of sound
};

/// D3Q27 lattice
template <typename T> struct D3Q27DescriptorBase {
  typedef D3Q27DescriptorBase<T> BaseDescriptor;
  enum { d = 3, q = 27 };     ///< number of dimensions/distr. functions
  static const int vicinity;  ///< size of neighborhood
  static const int c[q][d];   ///< lattice directions
  static const int opposite[q]; ///< opposite entry
  static const T t[q];        ///< lattice weights
  static const T invCs2;      ///< inverse square of speed of sound
};

template <typename T>struct AdvectionDiffusionD2Q5Descriptor
  : public D2Q5DescriptorBase<T>, public Velocity2dBase {
};

template <typename T> struct D2Q9Descriptor
  : public D2Q9DescriptorBase<T>, public NoExternalFieldBase {
};

template <typename T> struct ForcedD2Q9Descriptor
  : public D2Q9DescriptorBase<T>, public Force2dDescriptorBase {
};

template <typename T> struct V6ForcedD2Q9Descriptor
  : public D2Q9DescriptorBase<T>, public V6Force2dDescriptorBase {
};

template <typename T> struct AdvectionDiffusionD3Q7Descriptor
  : public D3Q7DescriptorBase<T>, public Velocity3dBase {
};

template <typename T> struct ParticleAdvectionDiffusionD3Q7Descriptor
  : public D3Q7DescriptorBase<T>, public ParticleAdvectionDiffusion3dDescriptorBase {
};

template <typename T> struct D3Q13Descriptor
  : public D3Q13DescriptorBase<T>, public NoExternalFieldBase {
};

template <typename T> struct ForcedD3Q13Descriptor
  : public D3Q13DescriptorBase<T>, public Force3dDescriptorBase {
};

template <typename T> struct D3Q15Descriptor
  : public D3Q15DescriptorBase<T>, public NoExternalFieldBase {
};

template <typename T> struct ForcedD3Q15Descriptor
  : public D3Q15DescriptorBase<T>, public Force3dDescriptorBase {
};

template <typename T> struct D3Q19Descriptor
  : public D3Q19DescriptorBase<T>, public NoExternalFieldBase {
};

template <typename T> struct ForcedD3Q19Descriptor
  : public D3Q19DescriptorBase<T>, public Force3dDescriptorBase {
};

template <typename T> struct V12ForcedD3Q19Descriptor
  : public D3Q19DescriptorBase<T>, public V12Force3dDescriptorBase {
};

template <typename T> struct ParticleAdvectionDiffusionD3Q19Descriptor
  : public D3Q19DescriptorBase<T>, public ParticleAdvectionDiffusion3dDescriptorBase {
};


template <typename T> struct D3Q27Descriptor
  : public D3Q27DescriptorBase<T>, public NoExternalFieldBase {
};

template <typename T> struct ForcedD3Q27Descriptor
  : public D3Q27DescriptorBase<T>, public Force3dDescriptorBase {
};

}  // namespace descriptors

}  // namespace olb

#endif
