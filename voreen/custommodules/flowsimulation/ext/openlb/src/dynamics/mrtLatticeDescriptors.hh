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
#ifndef MRT_LATTICE_DESCRIPTORS_HH
#define MRT_LATTICE_DESCRIPTORS_HH

namespace olb {
namespace descriptors {

// MRT D2Q9 ////////////////////////////////////////////////////////////

template<typename T>
const T MRTD2Q9DescriptorBase<T>::M[MRTD2Q9DescriptorBase<T>::q_][MRTD2Q9DescriptorBase<T>::q_] = {
  {(T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1},
  {-(T)4,(T)2,-(T)1, (T)2,-(T)1, (T)2,-(T)1, (T)2,-(T)1},
  {(T)4, (T)1,-(T)2, (T)1,-(T)2, (T)1,-(T)2, (T)1,-(T)2},
  {T(),-(T)1,-(T)1,-(T)1, T(), (T)1, (T)1, (T)1, T()},
  {T(),-(T)1, (T)2,-(T)1, T(), (T)1,-(T)2, (T)1, T()},
  {T(), (T)1, T(),-(T)1,-(T)1,-(T)1, T(), (T)1, (T)1},
  {T(), (T)1, T(),-(T)1, (T)2,-(T)1, T(), (T)1,-(T)2},
  {T(), T(), (T)1, T(),-(T)1, T(), (T)1, T(),-(T)1},
  {T(),-(T)1, T(), (T)1, T(),-(T)1, T(), (T)1, T()}
};

template<typename T>
const T MRTD2Q9DescriptorBase<T>::invM[MRTD2Q9DescriptorBase<T>::q_][MRTD2Q9DescriptorBase<T>::q_] = {
  {(T)1/(T)9, -(T)1/(T)9, (T)1/(T)9, T(), T(), T(), T(), T(), T()},
  {(T)1/(T)9, (T)1/(T)18, (T)1/(T)36, -(T)1/(T)6, -(T)1/(T)12, (T)1/(T)6, (T)1/(T)12, T(), -(T)1/(T)4},
  {(T)1/(T)9, -(T)1/(T)36, -(T)1/(T)18, -(T)1/(T)6, (T)1/(T)6, T(), T(), (T)1/(T)4, T()},
  {(T)1/(T)9, (T)1/(T)18, (T)1/(T)36, -(T)1/(T)6, -(T)1/(T)12, -(T)1/(T)6, -(T)1/(T)12, T(), (T)1/(T)4},
  {(T)1/(T)9, -(T)1/(T)36, -(T)1/(T)18, T(), T(), -(T)1/(T)6, (T)1/(T)6, -(T)1/(T)4, T()},
  {(T)1/(T)9, (T)1/(T)18, (T)1/(T)36, (T)1/(T)6, (T)1/(T)12, -(T)1/(T)6, -(T)1/(T)12, T(), -(T)1/(T)4},
  {(T)1/(T)9, -(T)1/(T)36, -(T)1/(T)18, (T)1/(T)6, -(T)1/(T)6, T(), T(), (T)1/(T)4, T()},
  {(T)1/(T)9, (T)1/(T)18, (T)1/(T)36, (T)1/(T)6, (T)1/(T)12, (T)1/(T)6, (T)1/(T)12, T(), (T)1/(T)4},
  {(T)1/(T)9, -(T)1/(T)36, -(T)1/(T)18, T(), T(), (T)1/(T)6, -(T)1/(T)6, -(T)1/(T)4, T()}
};

template<typename T>
const T MRTD2Q9DescriptorBase<T>::S[MRTD2Q9DescriptorBase<T>::q_] =
{ T(), (T)1.1, (T)1.1, T(), (T)1.1, T(), (T)1.1, T(), T() };
// s7=s8 to have a shear viscosity nu
// and the bulk viscosity depends on s2.


template<typename T>
const int MRTD2Q9DescriptorBase<T>::shearViscIndexes[MRTD2Q9DescriptorBase<T>::shearIndexes] = {7, 8};


// MRT D3Q19 ////////////////////////////////////////////////////////////

template<typename T>
const T MRTD3Q19DescriptorBase<T>::M[MRTD3Q19DescriptorBase<T>::q_][MRTD3Q19DescriptorBase<T>::q_] = {
  /*0*/      {
    (T)1,
    (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1,
    (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1
  },
  /*1*/       {
    (T)-30,
    (T)-11, (T)-11, (T)-11, (T)8, (T)8, (T)8, (T)8, (T)8, (T)8,
    (T)-11, (T)-11, (T)-11, (T)8, (T)8, (T)8, (T)8, (T)8, (T)8
  },
  /*2*/       {
    (T)12,
    (T)-4, (T)-4, (T)-4, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1,
    (T)-4, (T)-4, (T)-4, (T)1, (T)1, (T)1, (T)1, (T)1, (T)1
  },
  /*3*/       {
    T(),
    (T)-1, T(), T(), (T)-1, (T)-1, (T)-1, (T)-1, T(), T(),
    (T)1, T(), T(), (T)1, (T)1, (T)1, (T)1, T(), T()
  },
  /*4*/       {
    T(),
    (T)4, T(), T(), (T)-1, (T)-1, (T)-1, (T)-1, T(), T(),
    (T)-4, T(), T(), (T)1, (T)1, (T)1, (T)1, T(), T()
  },
  /*5*/       {
    T(),
    T(), (T)-1, T(), (T)-1, (T)1, T(), T(), (T)-1, (T)-1,
    T(), (T)1, T(), (T)1, (T)-1, T(), T(), (T)1, (T)1
  },
  /*6*/       {
    T(),
    T(), (T)4, T(), (T)-1, (T)1, T(), T(), (T)-1, (T)-1,
    T(), (T)-4, T(), (T)1, (T)-1, T(), T(), (T)1, (T)1
  },
  /*7*/       {
    T(),
    T(), T(), (T)-1, T(), T(), (T)-1, (T)1, (T)-1, (T)1,
    T(), T(), (T)1, T(), T(), (T)1, (T)-1, (T)1, (T)-1
  },
  /*8*/       {
    T(),
    T(), T(), (T)4, T(), T(), (T)-1, (T)1, (T)-1, (T)1,
    T(), T(), (T)-4, T(), T(), (T)1, (T)-1, (T)1, (T)-1
  },
  /*9*/       {
    T(),
    (T)2, (T)-1, (T)-1, (T)1, (T)1, (T)1, (T)1, (T)-2, (T)-2,
    (T)2, (T)-1, (T)-1, (T)1, (T)1, (T)1, (T)1, (T)-2, (T)-2
  },
  /*10*/      {
    T(),
    (T)-4, (T)2, (T)2, (T)1, (T)1, (T)1, (T)1, (T)-2, (T)-2,
    (T)-4, (T)2, (T)2, (T)1, (T)1, (T)1, (T)1, (T)-2, (T)-2
  },
  /*11*/      {
    T(),
    T(), (T)1, (T)-1, (T)1, (T)1, (T)-1, (T)-1, T(), T(),
    T(), (T)1, (T)-1, (T)1, (T)1, (T)-1, (T)-1, T(), T()
  },
  /*12*/      {
    T(),
    T(), (T)-2, (T)2, (T)1, (T)1, (T)-1, (T)-1, T(), T(),
    T(), (T)-2, (T)2, (T)1, (T)1, (T)-1, (T)-1, T(), T()
  },
  /*13*/      {
    T(),
    T(), T(), T(), (T)1, (T)-1, T(), T(), T(), T(),
    T(), T(), T(), (T)1, (T)-1, T(), T(), T(), T()
  },
  /*14*/      {
    T(),
    T(), T(), T(), T(), T(), T(), T(), (T)1, (T)-1,
    T(), T(), T(), T(), T(), T(), T(), (T)1, (T)-1
  },
  /*15*/      {
    T(),
    T(), T(), T(), T(), T(), (T)1, (T)-1, T(), T(),
    T(), T(), T(), T(), T(), (T)1, (T)-1, T(), T()
  },
  /*16*/      {
    T(),
    T(), T(), T(), (T)-1, (T)-1, (T)1, (T)1, T(), T(),
    T(), T(), T(), (T)1, (T)1, (T)-1, (T)-1, T(), T()
  },
  /*17*/      {
    T(),
    T(), T(), T(), (T)1, (T)-1, T(), T(), (T)-1, (T)-1,
    T(), T(), T(), (T)-1, (T)1, T(), T(), (T)1, (T)1
  },
  /*18*/      {
    T(),
    T(), T(), T(), T(), T(), (T)-1, (T)1, (T)1, (T)-1,
    T(), T(), T(), T(), T(), (T)1, (T)-1, (T)-1, (T)1
  }
};

template<typename T>
const T MRTD3Q19DescriptorBase<T>::invM[MRTD3Q19DescriptorBase<T>::q_][MRTD3Q19DescriptorBase<T>::q_] = {
  /*0*/       {
    (T)1/(T)19,
    -(T)5/(T)399,(T)1/(T)21,T(),T(),T(),T(),T(),T(),
    T(),T(),T(),T(),T(),T(),T(),T(),T(),T()
  },

  /*1*/       {
    (T)1/(T)19,
    -(T)11/(T)2394,-(T)1/(T)63,-(T)1/(T)10,(T)1/(T)10,T(),T(),T(),T(),(T)1/(T)18,
    -(T)1/(T)18,T(),T(),T(),T(),T(),T(),T(),T()
  },

  /*2*/       {
    (T)1/(T)19,
    -(T)11/(T)2394,-(T)1/(T)63,T(),T(),-(T)1/(T)10,(T)1/(T)10,T(),T(),-(T)1/(T)36,
    (T)1/(T)36,(T)1/(T)12,-(T)1/(T)12,T(),T(),T(),T(),T(),T()
  },

  /*3*/       {
    (T)1/(T)19,
    -(T)11/(T)2394,-(T)1/(T)63,T(),T(),T(),T(),-(T)1/(T)10,(T)1/(T)10,-(T)1/(T)36,
    (T)1/(T)36,-(T)1/(T)12,(T)1/(T)12,T(),T(),T(),T(),T(),T()
  },

  /*4*/       {
    (T)1/(T)19,
    (T)4/(T)1197,(T)1/(T)252,-(T)1/(T)10,-(T)1/(T)40,-(T)1/(T)10,-(T)1/(T)40,T(),T(),(T)1/(T)36,
    (T)1/(T)72,(T)1/(T)12,(T)1/(T)24,(T)1/(T)4,T(),T(),-(T)1/(T)8,(T)1/(T)8,T()
  },

  /*5*/       {
    (T)1/(T)19,
    (T)4/(T)1197,(T)1/(T)252,-(T)1/(T)10,-(T)1/(T)40,(T)1/(T)10,(T)1/(T)40,T(),T(),(T)1/(T)36,
    (T)1/(T)72,(T)1/(T)12,(T)1/(T)24,-(T)1/(T)4,T(),T(),-(T)1/(T)8,-(T)1/(T)8,T()
  },

  /*6*/       {
    (T)1/(T)19,
    (T)4/(T)1197,(T)1/(T)252,-(T)1/(T)10,-(T)1/(T)40,T(),T(),-(T)1/(T)10,-(T)1/(T)40,(T)1/(T)36,
    (T)1/(T)72,-(T)1/(T)12,-(T)1/(T)24,T(),T(),(T)1/(T)4,(T)1/(T)8,T(),-(T)1/(T)8
  },

  /*7*/       {
    (T)1/(T)19,
    (T)4/(T)1197,(T)1/(T)252,-(T)1/(T)10,-(T)1/(T)40,T(),T(),(T)1/(T)10,(T)1/(T)40,(T)1/(T)36,
    (T)1/(T)72,-(T)1/(T)12,-(T)1/(T)24,T(),T(),-(T)1/(T)4,(T)1/(T)8,T(),(T)1/(T)8
  },

  /*8*/       {
    (T)1/(T)19,
    (T)4/(T)1197,(T)1/(T)252,T(),T(),-(T)1/(T)10,-(T)1/(T)40,-(T)1/(T)10,-(T)1/(T)40,-(T)1/(T)18,
    -(T)1/(T)36,T(),T(),T(),(T)1/(T)4,T(),T(),-(T)1/(T)8,(T)1/(T)8
  },

  /*9*/       {
    (T)1/(T)19,
    (T)4/(T)1197,(T)1/(T)252,T(),T(),-(T)1/(T)10,-(T)1/(T)40,(T)1/(T)10,(T)1/(T)40,-(T)1/(T)18,
    -(T)1/(T)36,T(),T(),T(),-(T)1/(T)4,T(),T(),-(T)1/(T)8,-(T)1/(T)8
  },

  /*10*/      {
    (T)1/(T)19,
    -(T)11/(T)2394,-(T)1/(T)63,(T)1/(T)10,-(T)1/(T)10,T(),T(),T(),T(),(T)1/(T)18,
    -(T)1/(T)18,T(),T(),T(),T(),T(),T(),T(),T()
  },

  /*11*/      {
    (T)1/(T)19,
    -(T)11/(T)2394,-(T)1/(T)63,T(),T(),(T)1/(T)10,-(T)1/(T)10,T(),T(),-(T)1/(T)36,(T)1/(T)36,
    (T)1/(T)12,-(T)1/(T)12,T(),T(),T(),T(),T(),T()
  },

  /*12*/      {
    (T)1/(T)19,
    -(T)11/(T)2394,-(T)1/(T)63,T(),T(),T(),T(),(T)1/(T)10,-(T)1/(T)10,-(T)1/(T)36,
    (T)1/(T)36,-(T)1/(T)12,(T)1/(T)12,T(),T(),T(),T(),T(),T()
  },

  /*13*/      {
    (T)1/(T)19,
    (T)4/(T)1197,(T)1/(T)252,(T)1/(T)10,(T)1/(T)40,(T)1/(T)10,(T)1/(T)40,T(),T(),(T)1/(T)36,
    (T)1/(T)72,(T)1/(T)12,(T)1/(T)24,(T)1/(T)4,T(),T(),(T)1/(T)8,-(T)1/(T)8,T()
  },

  /*14*/      {
    (T)1/(T)19,
    (T)4/(T)1197,(T)1/(T)252,(T)1/(T)10,(T)1/(T)40,-(T)1/(T)10,-(T)1/(T)40,T(),T(),(T)1/(T)36,
    (T)1/(T)72,(T)1/(T)12,(T)1/(T)24,-(T)1/(T)4,T(),T(),(T)1/(T)8,(T)1/(T)8,T()
  },

  /*15*/      {
    (T)1/(T)19,
    (T)4/(T)1197,(T)1/(T)252,(T)1/(T)10,(T)1/(T)40,T(),T(),(T)1/(T)10,(T)1/(T)40,(T)1/(T)36,
    (T)1/(T)72,-(T)1/(T)12,-(T)1/(T)24,T(),T(),(T)1/(T)4,-(T)1/(T)8,T(),(T)1/(T)8
  },

  /*16*/      {
    (T)1/(T)19,
    (T)4/(T)1197,(T)1/(T)252,(T)1/(T)10,(T)1/(T)40,T(),T(),-(T)1/(T)10,-(T)1/(T)40,(T)1/(T)36,
    (T)1/(T)72,-(T)1/(T)12,-(T)1/(T)24,T(),T(),-(T)1/(T)4,-(T)1/(T)8,T(),-(T)1/(T)8
  },

  /*17*/      {
    (T)1/(T)19,
    (T)4/(T)1197,(T)1/(T)252,T(),T(),(T)1/(T)10,(T)1/(T)40,(T)1/(T)10,(T)1/(T)40,-(T)1/(T)18,
    -(T)1/(T)36,T(),T(),T(),(T)1/(T)4,T(),T(),(T)1/(T)8,-(T)1/(T)8
  },

  /*18*/      {
    (T)1/(T)19,
    (T)4/(T)1197,(T)1/(T)252,T(),T(),(T)1/(T)10,(T)1/(T)40,-(T)1/(T)10,-(T)1/(T)40,-(T)1/(T)18,
    -(T)1/(T)36,T(),T(),T(),-(T)1/(T)4,T(),T(),(T)1/(T)8,(T)1/(T)8
  }
};

template<typename T>
const T MRTD3Q19DescriptorBase<T>::S[MRTD3Q19DescriptorBase<T>::q_] = {


  // Original MRT Relaxation times
  /*s0*/  T(), // rho (conserved)
  /*s1*/  (T)1.19,
  /*s2*/  (T)1.4,
  /*s3*/  T(), // rho*ux (conserved)
  /*s4*/  (T)1.2,
  /*s5*/  T(), // rho*uy (conserved)
  /*s6*/  (T)1.2, // = s4
  /*s7*/  T(), // rho*uz (conserved)
  /*s8*/  (T)1.2, // = s4
  /*s9*/  T(), //should be equal to s13, used to define nu
  /*s10*/ (T)1.4,
  /*s11*/ T(), // = s9,
  /*s12*/ (T)1.4,
  /*s13*/ T(), //should be equal to s9, used to define nu
  /*s14*/ T(), // = s13,
  /*s15*/ T(), // = s13,
  /*s16*/ (T)1.98,
  /*s17*/ (T)1.98, // = s16,
  /*s18*/ (T)1.98  // = s16,


};

template<typename T>
const T MRTD3Q19DescriptorBase<T>::S_2[MRTD3Q19DescriptorBase<T>::q_] = {

  //  Use these relaxation time for higher stability
  /*s0*/  T(), // rho (conserved)
  /*s1*/  (T)1.0,
  /*s2*/  (T)1.0,
  /*s3*/  T(), // rho*ux (conserved)
  /*s4*/  (T)1.0,
  /*s5*/  T(), // rho*uy (conserved)
  /*s6*/  (T)1.0, // = s4
  /*s7*/  T(), // rho*uz (conserved)
  /*s8*/  (T)1.0, // = s4
  /*s9*/  T(), //should be equal to s13, used to define nu
  /*s10*/ (T)1.0,
  /*s11*/ T(), // = s9,
  /*s12*/ (T)1.0,
  /*s13*/ T(), //should be equal to s9, used to define nu
  /*s14*/ T(), // = s13,
  /*s15*/ T(), // = s13,
  /*s16*/ (T)1.0,
  /*s17*/ (T)1.0, // = s16,
  /*s18*/ (T)1.0  // = s16,

};

template<typename T>
const int MRTD3Q19DescriptorBase<T>::shearViscIndexes[MRTD3Q19DescriptorBase<T>::shearIndexes] =
{9, 11, 13, 14, 15};

}  // namespace descriptors

}  // namespace olb

#endif
