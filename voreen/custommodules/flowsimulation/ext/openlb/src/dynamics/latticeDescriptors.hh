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
 *  -- generic code
 */
#ifndef LATTICE_DESCRIPTORS_HH
#define LATTICE_DESCRIPTORS_HH

#include "latticeDescriptors.h"

namespace olb {

namespace descriptors {

// AdvectionDiffusion D2Q5 //////////////////////////////////////////////

template<typename T>
const int D2Q5DescriptorBase<T>::c[D2Q5DescriptorBase<T>::q][D2Q5DescriptorBase<T>::d] = {
  { 0, 0},
  {-1, 0}, {0, -1}, {1,0}, { 0,1}
};

template<typename T>
const int D2Q5DescriptorBase<T>::opposite[D2Q5DescriptorBase<T>::q] = {
  0, 3, 4, 1, 2
};

template<typename T>
const int D2Q5DescriptorBase<T>::vicinity = 1;

template<typename T>
const T D2Q5DescriptorBase<T>::invCs2 = (T)3;

template<typename T>
const T D2Q5DescriptorBase<T>::t[D2Q5DescriptorBase<T>::q] = {
  (T)1-(T)2/invCs2,
  (T)1/(invCs2*(T)2), (T)1/(invCs2*(T)2),
  (T)1/(invCs2*(T)2), (T)1/(invCs2*(T)2)
};

// D2Q5 AdvectionDiffusionMRT ////////////////////////////////////////////////////////////
/*
 * Based on: Liu, Q., & He, Y. L. (2015). Double multiple-relaxation-time lattice Boltzmann model
 *           for solid–liquid phase change with natural convection in porous media.
 *           Physica A: Statistical Mechanics and its Applications, 438, 94-106.
 */

template<typename T>
const int AdvectionDiffusionMRTD2Q5DescriptorBase<T>::c
[AdvectionDiffusionMRTD2Q5DescriptorBase<T>::q][AdvectionDiffusionMRTD2Q5DescriptorBase<T>::d] = {
  { 0, 0},
  {-1, 0},
  { 0,-1},
  { 1, 0},
  { 0, 1}
};

template<typename T>
const int AdvectionDiffusionMRTD2Q5DescriptorBase<T>::opposite[AdvectionDiffusionMRTD2Q5DescriptorBase<T>::q] = {
  0, 3, 4, 1, 2
};

template<typename T>
const int AdvectionDiffusionMRTD2Q5DescriptorBase<T>::vicinity = 1;

template<typename T>
const T AdvectionDiffusionMRTD2Q5DescriptorBase<T>::M[AdvectionDiffusionMRTD2Q5DescriptorBase<T>::q][AdvectionDiffusionMRTD2Q5DescriptorBase<T>::q] = {
  {(T)1 , (T)1, (T)1, (T)1, (T)1},
  {T()  ,-(T)1, T() , (T)1, T() },
  {T()  , T() ,-(T)1, T() , (T)1},
  {-(T)4, (T)1, (T)1, (T)1, (T)1},
  {T()  , (T)1,-(T)1, (T)1,-(T)1}
};

template<typename T>
const T AdvectionDiffusionMRTD2Q5DescriptorBase<T>::invM[AdvectionDiffusionMRTD2Q5DescriptorBase<T>::q][AdvectionDiffusionMRTD2Q5DescriptorBase<T>::q] = {
  {(T)1/(T)5,   T(),         T(),        -(T)1/(T)5,    T()},
  {(T)1/(T)5,  -(T)1/(T)2,   T(),         (T)1/(T)20,   (T)1/(T)4},
  {(T)1/(T)5,   T(),        -(T)1/(T)2,   (T)1/(T)20,  -(T)1/(T)4},
  {(T)1/(T)5,   (T)1/(T)2,   T(),         (T)1/(T)20,   (T)1/(T)4},
  {(T)1/(T)5,   T(),         (T)1/(T)2,   (T)1/(T)20,  -(T)1/(T)4}
};

template<typename T>
const T AdvectionDiffusionMRTD2Q5DescriptorBase<T>::S[AdvectionDiffusionMRTD2Q5DescriptorBase<T>::q] =
{T(), T(), T(), (T)1.5, (T)1.5};

template<typename T>
const int AdvectionDiffusionMRTD2Q5DescriptorBase<T>::shearViscIndexes[AdvectionDiffusionMRTD2Q5DescriptorBase<T>::shearIndexes] = {1, 2};

template<typename T>
const T AdvectionDiffusionMRTD2Q5DescriptorBase<T>::invCs2 = (T)3;

template<typename T>
const T AdvectionDiffusionMRTD2Q5DescriptorBase<T>::t[AdvectionDiffusionMRTD2Q5DescriptorBase<T>::q] = {
  (T)1-(T)2/invCs2,
  (T)1/(invCs2*(T)2), (T)1/(invCs2*(T)2),
  (T)1/(invCs2*(T)2), (T)1/(invCs2*(T)2)
};

// D2Q9 ////////////////////////////////////////////////////////////

template<typename T>
const int D2Q9DescriptorBase<T>::vicinity = 1;

template<typename T>
const int D2Q9DescriptorBase<T>::c
[D2Q9DescriptorBase<T>::q][D2Q9DescriptorBase<T>::d] = {
  { 0, 0},
  {-1, 1}, {-1, 0}, {-1,-1}, { 0,-1},
  { 1,-1}, { 1, 0}, { 1, 1}, { 0, 1}
};

template<typename T>
const int D2Q9DescriptorBase<T>::opposite[D2Q9DescriptorBase<T>::q] = {
  0, 5, 6, 7, 8, 1, 2, 3, 4
};


template<typename T>
const T D2Q9DescriptorBase<T>::t[D2Q9DescriptorBase<T>::q] = {
  (T)4/(T)9, (T)1/(T)36, (T)1/(T)9, (T)1/(T)36, (T)1/(T)9,
  (T)1/(T)36, (T)1/(T)9, (T)1/(T)36, (T)1/(T)9
};

template<typename T>
const T D2Q9DescriptorBase<T>::invCs2 = (T)3;

// AdvectionDiffusion D3Q7 ////////////////////////////////////////////////////

template<typename T>
const int D3Q7DescriptorBase<T>::c[D3Q7DescriptorBase<T>::q][D3Q7DescriptorBase<T>::d] = {
  { 0, 0, 0},

  {-1, 0, 0}, {0,-1, 0},
  { 0, 0,-1}, {1, 0, 0},
  { 0, 1, 0}, {0, 0, 1},
};

template<typename T>
const int D3Q7DescriptorBase<T>::opposite[D3Q7DescriptorBase<T>::q] = {
  0, 4, 5, 6, 1, 2, 3
};

template<typename T>
const int D3Q7DescriptorBase<T>::vicinity = 1;

template<typename T>
const T D3Q7DescriptorBase<T>::invCs2 = (T)4;

template<typename T>
const T D3Q7DescriptorBase<T>::t[D3Q7DescriptorBase<T>::q] = {
  (T)1 -(T)3 / invCs2,

  (T)1/(invCs2*(T)2), (T)1/(invCs2*(T)2), (T)1/(invCs2*(T)2),
  (T)1/(invCs2*(T)2), (T)1/(invCs2*(T)2), (T)1/(invCs2*(T)2)
};



// D3Q7 AdvectionDiffusionMRT ////////////////////////////////////////////////////////////
/*
 * Based on: Wu, H., Wang, J., & Tao, Z. (2011). Passive heat transfer in a turbulent
 *           channel flow simulation using large eddy simulation based on the lattice
 *           Boltzmann method framework.
 *           International Journal of Heat and Fluid Flow, 32(6), 1111-1119.
 *
 * There are some differences in respect to the order of the columns based on the lattice directions
 *
 * TODO @AP: Check the D3Q7 M and invM matrices, if they are consistent to the OpenLB lattice directions
 */
template<typename T>
const int AdvectionDiffusionMRTD3Q7DescriptorBase<T>::c
[AdvectionDiffusionMRTD3Q7DescriptorBase<T>::q][AdvectionDiffusionMRTD3Q7DescriptorBase<T>::d] = {
      { 0, 0, 0},
      {-1, 0, 0},
      { 0,-1, 0},
      { 0, 0,-1},
      { 1, 0, 0},
      { 0, 1, 0},
      { 0, 0, 1},
};

template<typename T>
const int AdvectionDiffusionMRTD3Q7DescriptorBase<T>::opposite[AdvectionDiffusionMRTD3Q7DescriptorBase<T>::q] = {
  0, 4, 5, 6, 1, 2, 3
};

template<typename T>
const int AdvectionDiffusionMRTD3Q7DescriptorBase<T>::vicinity = 1;

template<typename T>
const T AdvectionDiffusionMRTD3Q7DescriptorBase<T>::M[AdvectionDiffusionMRTD3Q7DescriptorBase<T>::q][AdvectionDiffusionMRTD3Q7DescriptorBase<T>::q] = {
//  Maria Lloret OpenLB Guide: In my oppinion the matrix is wrong taking the OpenLB definition for lattice directions
//  {(T)1,  (T)1,   (T)1,   (T)1,   (T)1,   (T)1,   (T)1},
//  {T(),   -(T)1,  T(),    T(),    (T)1,   T(),    T() },
//  {T(),   T(),    -(T)1,  T(),    T(),    (T)1,   T() },
//  {T(),   T(),    T(),    -(T)1,  T(),    T(),    (T)1},
//  {-(T)6, (T)1,   (T)1,   (T)1,   (T)1,   (T)1,   (T)1},
//  {T(),   (T)1,   -(T)1,  T(),    (T)1,   -(T)1,   T()},
//  {T(),   (T)1,   (T)1,   -(T)2,  (T)1,   (T)1,   -(T)2}

//  Wu et.al 2011: The directions are modified for the OpenLB definition // might not be correct
//    {(T)1,  (T)1,   (T)1,   (T)1,   (T)1,   (T)1,   (T)1 },
//    {T(),   (T)1,   T(),   -(T)1,   T(),    T(),    T()  },
//    {T(),   T(),   -(T)1,   T(),    T(),    (T)1,   T()  },
//    {T(),   T(),    T(),    T(),    (T)1,   T(),    -(T)1},
//    {-(T)6, (T)1,   (T)1,   (T)1,   (T)1,   (T)1,   (T)1 },
//    {T(),   (T)1,  -(T)1,   (T)1,   T(),    -(T)1,  T()  },
//    {T(),   (T)1,   (T)1,   (T)1,   (T)2,   (T)1,   -(T)2}

//  Li, Yang et al 2016: The directions are modified for the OpenLB definition
    {(T)1,  (T)1,   (T)1,   (T)1,   (T)1,   (T)1,   (T)1 },
    {T(),   (T)1,   T(),    -(T)1,  T(),    T(),    T()  },
    {T(),   T(),    -(T)1,  T(),    T(),    (T)1,   T()  },
    {T(),   T(),    T(),    T(),    (T)1,   T(),    -(T)1},
    {(T)6,  -(T)1,  -(T)1,  -(T)1,  -(T)1,  -(T)1,  -(T)1},
    {T(),   (T)2,   -(T)1,  (T)2,   -(T)1,  -(T)1,  -(T)1},
    {T(),   T(),    (T)1,   T(),    -(T)1,  (T)1,   -(T)1}

};

template<typename T>
const T AdvectionDiffusionMRTD3Q7DescriptorBase<T>::invM[AdvectionDiffusionMRTD3Q7DescriptorBase<T>::q][AdvectionDiffusionMRTD3Q7DescriptorBase<T>::q] = {
//  Maria Lloret OpenLB Guide: In my oppinion the matrix is wrong taking the OpenLB definition for lattice directions
//  {(T)1/(T)7,   T(),          T(),          T(),        -(T)1/(T)7,     T(),          T()},
//  {(T)1/(T)7,   -(T)1/(T)2,   T(),          T(),        (T)1/(T)42,     (T)1/(T)4,    (T)1/(T)12},
//  {(T)1/(T)7,   T(),          -(T)1/(T)2,   T(),        (T)1/(T)42,     -(T)1/(T)4,   (T)1/(T)12},
//  {(T)1/(T)7,   T(),          T(),          -(T)1/(T)2, (T)1/(T)42,     T(),          -(T)1/(T)6},
//  {(T)1/(T)7,  (T)1/(T)2,     T(),          T(),        (T)1/(T)42,     (T)1/(T)4,    (T)1/(T)12},
//  {(T)1/(T)7,   T(),          (T)1/(T)2,    T(),        (T)1/(T)42,     -(T)1/(T)4,   (T)1/(T)12},
//  {(T)1/(T)7,   T(),          T(),          (T)1/(T)2,  (T)1/(T)42,     T(),          -(T)1/(T)6}

//  Wu et.al 2011: The directions are modified for the OpenLB definition // might not be correct
//    {(T)1/(T)7,   T(),          T(),          T(),        -(T)1/(T)7,     T(),              T()},
//    {(T)1/(T)7,   (T)1/(T)2,    T(),          T(),        (T)1/(T)42,     (T)1/(T)4,        (T)1/(T)12},
//    {(T)1/(T)7,   T(),          -(T)1/(T)2,   T(),        (T)1/(T)42,     -(T)1/(T)4,       (T)1/(T)12},
//    {(T)1/(T)7,   -(T)1/(T)2,   T(),          T(),        (T)1/(T)42,     (T)1/(T)4,        (T)1/(T)12},
//    {(T)1/(T)7,   T(),          T(),          (T)1/(T)2,  (T)1/(T)42,     T(),              -(T)1/(T)6},
//    {(T)1/(T)7,   T(),          (T)1/(T)2,    T(),        (T)1/(T)42,     -(T)1/(T)4,       (T)1/(T)12},
//    {(T)1/(T)7,   T(),          T(),          -(T)1/(T)2, (T)1/(T)42,     T(),              -(T)1/(T)6}

//  Li, Yang et al 2016: The directions are modified for the OpenLB definition
    {(T)1/(T)7,   T(),          T(),          T(),        (T)1/(T)7,       T(),              T()},
    {(T)1/(T)7,   (T)1/(T)2,    T(),          T(),        -(T)1/(T)42,     (T)1/(T)6,        T()},
    {(T)1/(T)7,   T(),          -(T)1/(T)2,   T(),        -(T)1/(T)42,     -(T)1/(T)12,      (T)1/(T)4},
    {(T)1/(T)7,   -(T)1/(T)2,   T(),          T(),        -(T)1/(T)42,     (T)1/(T)6,        T()},
    {(T)1/(T)7,   T(),          T(),          (T)1/(T)2,  -(T)1/(T)42,     -(T)1/(T)12,      -(T)1/(T)4},
    {(T)1/(T)7,   T(),          (T)1/(T)2,    T(),        -(T)1/(T)42,     -(T)1/(T)12,      (T)1/(T)4},
    {(T)1/(T)7,   T(),          T(),          -(T)1/(T)2, -(T)1/(T)42,     -(T)1/(T)12,      -(T)1/(T)4}

};

template<typename T>
const T AdvectionDiffusionMRTD3Q7DescriptorBase<T>::S[AdvectionDiffusionMRTD3Q7DescriptorBase<T>::q] = {
  // Original MRT Relaxation times
  /*s0*/  T(),  // rho (conserved)
  /*s1*/  T(),  // Function of the thermal diffusivity: S_a = 1/t_a = 1/(4*a + 1/2)
  /*s2*/  T(),  // Function of the thermal diffusivity: S_a = 1/t_a = 1/(4*a + 1/2)
  /*s3*/  T(),  // Function of the thermal diffusivity: S_a = 1/t_a = 1/(4*a + 1/2)
  /*s4*/  (T)1.9,
  /*s5*/  (T)1.9,
  /*s6*/  (T)1.9
};

template<typename T>
const int AdvectionDiffusionMRTD3Q7DescriptorBase<T>::shearViscIndexes[AdvectionDiffusionMRTD3Q7DescriptorBase<T>::shearIndexes] = {1, 2, 3};

template<typename T>
const T AdvectionDiffusionMRTD3Q7DescriptorBase<T>::invCs2 = (T)4;

template<typename T>
const T AdvectionDiffusionMRTD3Q7DescriptorBase<T>::t[AdvectionDiffusionMRTD3Q7DescriptorBase<T>::q] = {
  (T)1-(T)3/invCs2,
  (T)1/(invCs2*(T)2), (T)1/(invCs2*(T)2), (T)1/(invCs2*(T)2),
  (T)1/(invCs2*(T)2), (T)1/(invCs2*(T)2), (T)1/(invCs2*(T)2)
};

// D3Q13 ///////////////////////////////////////////////////////////

template<typename T>
const int D3Q13DescriptorBase<T>::vicinity = 1;

template<typename T>
const int D3Q13DescriptorBase<T>::c
[D3Q13DescriptorBase<T>::q][D3Q13DescriptorBase<T>::d] = {
  { 0, 0, 0},

  {-1,-1, 0}, {-1, 1, 0}, {-1, 0,-1},
  {-1, 0, 1}, { 0,-1,-1}, { 0,-1, 1},

  { 1, 1, 0}, { 1,-1, 0}, { 1, 0, 1},
  { 1, 0,-1}, { 0, 1, 1}, { 0, 1,-1}
};

template<typename T>
const int D3Q13DescriptorBase<T>::opposite[D3Q13DescriptorBase<T>::q] = {
  0, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6
};

template<typename T>
const T D3Q13DescriptorBase<T>::t[D3Q13DescriptorBase<T>::q] = {
  (T)1/(T)2,

  (T)1/(T)24, (T)1/(T)24, (T)1/(T)24,
  (T)1/(T)24, (T)1/(T)24, (T)1/(T)24,

  (T)1/(T)24, (T)1/(T)24, (T)1/(T)24,
  (T)1/(T)24, (T)1/(T)24, (T)1/(T)24
};


/** This parameter is chosen to enhance numerical stability */
template<typename T>
const T D3Q13DescriptorBase<T>::invCs2 = (T)3;

/** This parameter is chosen to enhance numerical stability */
template<typename T>
const T D3Q13DescriptorBase<T>::lambda_e = (T)1.5;

/** This parameter is chosen to enhance numerical stability */
template<typename T>
const T D3Q13DescriptorBase<T>::lambda_h = (T)1.8;


// D3Q15 ///////////////////////////////////////////////////////////

template<typename T>
const int D3Q15DescriptorBase<T>::vicinity = 1;

template<typename T>
const int D3Q15DescriptorBase<T>::c
[D3Q15DescriptorBase<T>::q][D3Q15DescriptorBase<T>::d] = {
  { 0, 0, 0},

  {-1, 0, 0}, { 0,-1, 0}, { 0, 0,-1},
  {-1,-1,-1}, {-1,-1, 1}, {-1, 1,-1}, {-1, 1, 1},

  { 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1},
  { 1, 1, 1}, { 1, 1,-1}, { 1,-1, 1}, { 1,-1,-1},

};

template<typename T>
const int D3Q15DescriptorBase<T>::opposite[D3Q15DescriptorBase<T>::q] = {
  0, 8, 9, 10, 11, 12, 13, 14, 1, 2, 3, 4, 5, 6, 7
};

template<typename T>
const T D3Q15DescriptorBase<T>::t[D3Q15DescriptorBase<T>::q] = {
  (T)2/(T)9,

  (T)1/(T)9, (T)1/(T)9, (T)1/(T)9,
  (T)1/(T)72, (T)1/(T)72, (T)1/(T)72, (T)1/(T)72,

  (T)1/(T)9, (T)1/(T)9, (T)1/(T)9,
  (T)1/(T)72, (T)1/(T)72, (T)1/(T)72, (T)1/(T)72
};

template<typename T>
const T D3Q15DescriptorBase<T>::invCs2 = (T)3;


// D3Q19 ///////////////////////////////////////////////////////////

template<typename T>
const int D3Q19DescriptorBase<T>::vicinity = 1;

template<typename T>
const int D3Q19DescriptorBase<T>::c
[D3Q19DescriptorBase<T>::q][D3Q19DescriptorBase<T>::d] = {
  { 0, 0, 0},

  {-1, 0, 0}, { 0,-1, 0}, { 0, 0,-1},
  {-1,-1, 0}, {-1, 1, 0}, {-1, 0,-1},
  {-1, 0, 1}, { 0,-1,-1}, { 0,-1, 1},

  { 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1},
  { 1, 1, 0}, { 1,-1, 0}, { 1, 0, 1},
  { 1, 0,-1}, { 0, 1, 1}, { 0, 1,-1}
};

template<typename T>
const int D3Q19DescriptorBase<T>::opposite[D3Q19DescriptorBase<T>::q] = {
  0, 10, 11, 12, 13, 14, 15, 16, 17, 18, 1, 2, 3, 4, 5, 6, 7, 8, 9
};

template<typename T>
const T D3Q19DescriptorBase<T>::t[D3Q19DescriptorBase<T>::q] = {
  (T)1/(T)3,

  (T)1/(T)18, (T)1/(T)18, (T)1/(T)18,
  (T)1/(T)36, (T)1/(T)36, (T)1/(T)36,
  (T)1/(T)36, (T)1/(T)36, (T)1/(T)36,

  (T)1/(T)18, (T)1/(T)18, (T)1/(T)18,
  (T)1/(T)36, (T)1/(T)36, (T)1/(T)36,
  (T)1/(T)36, (T)1/(T)36, (T)1/(T)36
};

template<typename T>
const T D3Q19DescriptorBase<T>::invCs2 = (T)3;


// D3Q27 ///////////////////////////////////////////////////////////

template<typename T>
const int D3Q27DescriptorBase<T>::vicinity = 1;

template<typename T>
const int D3Q27DescriptorBase<T>::c
[D3Q27DescriptorBase<T>::q][D3Q27DescriptorBase<T>::d] = {
  { 0, 0, 0},

  {-1, 0, 0}, { 0,-1, 0}, { 0, 0,-1},
  {-1,-1, 0}, {-1, 1, 0}, {-1, 0,-1},
  {-1, 0, 1}, { 0,-1,-1}, { 0,-1, 1},
  {-1,-1,-1}, {-1,-1, 1}, {-1, 1,-1}, {-1, 1, 1},

  { 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1},
  { 1, 1, 0}, { 1,-1, 0}, { 1, 0, 1},
  { 1, 0,-1}, { 0, 1, 1}, { 0, 1,-1},
  { 1, 1, 1}, { 1, 1,-1}, { 1,-1, 1}, { 1,-1,-1}
};

template<typename T>
const int D3Q27DescriptorBase<T>::opposite[D3Q27DescriptorBase<T>::q] = {
  0, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26,
  1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13
};

template<typename T>
const T D3Q27DescriptorBase<T>::t[D3Q27DescriptorBase<T>::q] = {
  (T)8/(T)27,

  (T)2/(T)27, (T)2/(T)27, (T)2/(T)27,
  (T)1/(T)54, (T)1/(T)54, (T)1/(T)54,
  (T)1/(T)54, (T)1/(T)54, (T)1/(T)54,
  (T)1/(T)216, (T)1/(T)216, (T)1/(T)216, (T)1/(T)216,

  (T)2/(T)27, (T)2/(T)27, (T)2/(T)27,
  (T)1/(T)54, (T)1/(T)54, (T)1/(T)54,
  (T)1/(T)54, (T)1/(T)54, (T)1/(T)54,
  (T)1/(T)216, (T)1/(T)216, (T)1/(T)216, (T)1/(T)216
};

template<typename T>
const T D3Q27DescriptorBase<T>::invCs2 = (T)3;


}  // namespace descriptors

}  // namespace olb

#endif
