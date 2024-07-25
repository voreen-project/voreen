/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Davide Dapelo
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
 * 3D postprocessor to locally update the finite-diffusion advection-diffusion external field.
 *  -- generic implementation
 */
#ifndef AD_POST_PROCESSOR_3D_HH
#define AD_POST_PROCESSOR_3D_HH

#include <string>
namespace olb {


////////  FdBasePostProcessor3D ///////////////////////////////////

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
FdBasePostProcessor3D<T,DESCRIPTOR,FIELD,SOURCE>::FdBasePostProcessor3D(size_t& iT)
  : _iT(iT)
{
  static_assert(DESCRIPTOR::template size<FIELD>()  == 2, "FIELD must have size 2." );
  this->getName() = "FdBasePostProcessor3D";
}

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
int FdBasePostProcessor3D<T,DESCRIPTOR,FIELD,SOURCE>::getIT() const
{
  return _iT;
}

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
std::vector<std::vector<T>> FdBasePostProcessor3D <T,DESCRIPTOR,FIELD,SOURCE>::
                         initMatrix(int size0)
{
  return std::vector<std::vector<T>> ( size0, std::vector<T>(DESCRIPTOR::d, T()) );
}

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
template<typename SOURCE_CHK>
std::enable_if_t<  std::is_void<SOURCE_CHK>::value, void> FdBasePostProcessor3D <T,DESCRIPTOR,FIELD,SOURCE>::
applySourceTerm(T* fNew, Cell<T,DESCRIPTOR>& cell)
{ }

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
template<typename SOURCE_CHK>
std::enable_if_t<! std::is_void<SOURCE_CHK>::value, void> FdBasePostProcessor3D <T,DESCRIPTOR,FIELD,SOURCE>::
applySourceTerm(T* fNew, Cell<T,DESCRIPTOR>& cell)
{
  //T before = *fNew; // <---
  *fNew += cell.template getFieldPointer<SOURCE_CHK>()[0];
  //T after = *fNew; // <---
  //std::cout << "BEFORE=" << before << "     AFTER=" << after << std::endl; // <---
}


//////// FdBasePostProcessorGenerator3D ///////////////////////////////////

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
FdBasePostProcessorGenerator3D<T,DESCRIPTOR,FIELD,SOURCE>::FdBasePostProcessorGenerator3D (
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, std::size_t& iT)
  : PostProcessorGenerator3D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_, z0_, z1_),
    _iT(iT)
{ }


////////  FdPostProcessor3D ///////////////////////////////////

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
FdPostProcessor3D <T,DESCRIPTOR,FIELD,SOURCE>::
FdPostProcessor3D ( int x0, int x1, int y0, int y1, int z0, int z1, std::size_t& iT, std::shared_ptr<FdModel<T,DESCRIPTOR>> model)
  :  FdBasePostProcessor3D<T,DESCRIPTOR,FIELD,SOURCE>(iT),
     _x0(x0), _x1(x1), _y0(y0), _y1(y1), _z0(z0), _z1(z1),
     _model(model),
     f{this->initMatrix(_model->extent())},
     F{this->initMatrix(_model->extent())}
{ }

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
FdPostProcessor3D <T,DESCRIPTOR,FIELD,SOURCE>::
FdPostProcessor3D ( std::size_t& iT, std::shared_ptr<FdModel<T,DESCRIPTOR>> model)
  :  FdPostProcessor3D<T,DESCRIPTOR,FIELD,SOURCE>(0,0,0,0,0,0,iT,model)
{ }

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
void FdPostProcessor3D<T,DESCRIPTOR,FIELD,SOURCE>::
processSubDomain( BlockLattice<T,DESCRIPTOR>& blockLattice,
                  int x0, int x1, int y0, int y1, int z0, int z1 )
{
  int newX0, newX1, newY0, newY1, newZ0, newZ1;
  if ( util::intersect (
         x0, x1, y0, y1, z0, z1,
         _x0, _x1, _y0, _y1, _z0, _z1,
         newX0, newX1, newY0, newY1, newZ0, newZ1 ) ) {

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        for (int iZ=newZ0; iZ<=newZ1; ++iZ) {
          T* fNew = fd::accessNew<T,DESCRIPTOR,FIELD>( blockLattice.get ( iX, iY, iZ    ), this->getIT() );
          T* f0   = fd::accessOld<T,DESCRIPTOR,FIELD>( blockLattice.get ( iX, iY, iZ    ), this->getIT() );
          for (int iN=1; iN<=_model->extent(); ++iN) {
            T* fx = fd::accessOld<T,DESCRIPTOR,FIELD>( blockLattice.get ( iX-iN, iY, iZ    ), this->getIT() );
            T* fX = fd::accessOld<T,DESCRIPTOR,FIELD>( blockLattice.get ( iX+iN, iY, iZ    ), this->getIT() );
            T* fy = fd::accessOld<T,DESCRIPTOR,FIELD>( blockLattice.get ( iX, iY-iN, iZ    ), this->getIT() );
            T* fY = fd::accessOld<T,DESCRIPTOR,FIELD>( blockLattice.get ( iX, iY+iN, iZ    ), this->getIT() );
            T* fz = fd::accessOld<T,DESCRIPTOR,FIELD>( blockLattice.get ( iX, iY, iZ-iN ), this->getIT() );
            T* fZ = fd::accessOld<T,DESCRIPTOR,FIELD>( blockLattice.get ( iX, iY, iZ+iN ), this->getIT() );
            f[iN-1][0] = *fx;
            F[iN-1][0] = *fX;
            f[iN-1][1] = *fy;
            F[iN-1][1] = *fY;
            f[iN-1][2] = *fz;
            F[iN-1][2] = *fZ;
          }
          Cell<T,DESCRIPTOR> cell = blockLattice.get(iX,iY,iZ);
          _model->operator()(fNew, f0, f, F, cell);
          this->applySourceTerm(fNew, cell);
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
void FdPostProcessor3D<T,DESCRIPTOR,FIELD,SOURCE>::
process(BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, _x0, _x1, _y0, _y1, _z0, _z1);
}


//////// FdPostProcessorGenerator3D ///////////////////////////////////

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
FdPostProcessorGenerator3D<T,DESCRIPTOR,FIELD,SOURCE>::FdPostProcessorGenerator3D (
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, std::size_t& iT, std::shared_ptr<FdModel<T,DESCRIPTOR>> model )
  : FdBasePostProcessorGenerator3D<T,DESCRIPTOR,FIELD,SOURCE>(x0_, x1_, y0_, y1_, z0_, z1_, iT),
    _model(model)
{ }

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
FdPostProcessorGenerator3D<T,DESCRIPTOR,FIELD,SOURCE>::FdPostProcessorGenerator3D(size_t& iT, std::shared_ptr<FdModel<T,DESCRIPTOR>> model)
  : FdPostProcessorGenerator3D<T,DESCRIPTOR,FIELD,SOURCE>(0, 0, 0, 0, 0, 0, iT, model)
{ }

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
PostProcessor3D<T,DESCRIPTOR>* FdPostProcessorGenerator3D<T,DESCRIPTOR,FIELD,SOURCE>::generate() const
{
  return new FdPostProcessor3D<T,DESCRIPTOR,FIELD,SOURCE>(this->x0,this->x1,this->y0,this->y1,this->z0,this->z1,
         this->_iT,_model);
}

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
PostProcessorGenerator3D<T,DESCRIPTOR>* FdPostProcessorGenerator3D<T,DESCRIPTOR,FIELD,SOURCE>::clone() const
{
  return new FdPostProcessorGenerator3D<T,DESCRIPTOR,FIELD,SOURCE>(*this);
}


}  // namespace olb

#endif
