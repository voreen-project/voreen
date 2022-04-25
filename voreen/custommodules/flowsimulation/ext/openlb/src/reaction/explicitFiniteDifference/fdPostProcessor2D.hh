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
 * 2D postprocessor to locally update the finite-diffusion advection-diffusion external field.
 *  -- generic implementation
 */
#ifndef AD_POST_PROCESSOR_2D_HH
#define AD_POST_PROCESSOR_2D_HH

namespace olb {


////////  FdBasePostProcessor2D ///////////////////////////////////

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
FdBasePostProcessor2D<T,DESCRIPTOR,FIELD,SOURCE>::FdBasePostProcessor2D(size_t& iT)
  : _iT(iT)
{
  static_assert(DESCRIPTOR::template size<FIELD>()  == 2, "FIELD must have size 2." );
  this->getName() = "FdBasePostProcessor2D";
}

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
int FdBasePostProcessor2D<T,DESCRIPTOR,FIELD,SOURCE>::getIT() const
{
  return _iT;
}

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
std::vector<std::vector<T>> FdBasePostProcessor2D <T,DESCRIPTOR,FIELD,SOURCE>::
                         initMatrix(int size0)
{
  return std::vector<std::vector<T>> ( size0, std::vector<T>(DESCRIPTOR::d, T()) );
}

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
template<typename SOURCE_CHK>
std::enable_if_t<  std::is_void<SOURCE_CHK>::value, void> FdBasePostProcessor2D <T,DESCRIPTOR,FIELD,SOURCE>::
applySourceTerm(T* fNew, Cell<T,DESCRIPTOR>& cell)
{ }

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
template<typename SOURCE_CHK>
std::enable_if_t<! std::is_void<SOURCE_CHK>::value, void> FdBasePostProcessor2D <T,DESCRIPTOR,FIELD,SOURCE>::
applySourceTerm(T* fNew, Cell<T,DESCRIPTOR>& cell)
{
  *fNew += cell.template getFieldPointer<SOURCE_CHK>()[0];
}


//////// FdBasePostProcessorGenerator2D ///////////////////////////////////

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
FdBasePostProcessorGenerator2D<T,DESCRIPTOR,FIELD,SOURCE>::FdBasePostProcessorGenerator2D (
  int x0_, int x1_, int y0_, int y1_, std::size_t& iT )
  : PostProcessorGenerator2D<T,DESCRIPTOR>(x0_, x1_, y0_, y1_),
    _iT(iT)
{ }


////////  FdPostProcessor2D ///////////////////////////////////

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
FdPostProcessor2D <T,DESCRIPTOR,FIELD,SOURCE>::
FdPostProcessor2D ( int x0, int x1, int y0, int y1, std::size_t& iT, std::shared_ptr<FdModel<T,DESCRIPTOR>> model)
  :  FdBasePostProcessor2D<T,DESCRIPTOR,FIELD,SOURCE>(iT),
     _x0(x0), _x1(x1), _y0(y0), _y1(y1),
     _model(model),
     f{this->initMatrix(_model->extent())},
     F{this->initMatrix(_model->extent())}
{ }

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
FdPostProcessor2D <T,DESCRIPTOR,FIELD,SOURCE>::
FdPostProcessor2D ( std::size_t& iT, std::shared_ptr<FdModel<T,DESCRIPTOR>> model)
  :  FdPostProcessor2D<T,DESCRIPTOR,FIELD,SOURCE>(0,0,0,0,iT,model)
{ }

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
void FdPostProcessor2D<T,DESCRIPTOR,FIELD,SOURCE>::
processSubDomain( BlockLattice<T,DESCRIPTOR>& blockLattice,
                  int x0, int x1, int y0, int y1 )
{
  int newX0, newX1, newY0, newY1;
  if ( util::intersect (
         x0, x1, y0, y1,
         _x0, _x1, _y0, _y1,
         newX0, newX1, newY0, newY1 ) ) {

    for (int iX=newX0; iX<=newX1; ++iX) {
      for (int iY=newY0; iY<=newY1; ++iY) {
        T* fNew = fd::accessNew<T,DESCRIPTOR,FIELD>( blockLattice.get ( iX, iY    ), this->getIT() );
        T* f0   = fd::accessOld<T,DESCRIPTOR,FIELD>( blockLattice.get ( iX, iY    ), this->getIT() );
        for (int iN=1; iN<=_model->extent(); ++iN) {
          T* fx = fd::accessOld<T,DESCRIPTOR,FIELD>( blockLattice.get ( iX-iN, iY    ), this->getIT() );
          T* fX = fd::accessOld<T,DESCRIPTOR,FIELD>( blockLattice.get ( iX+iN, iY    ), this->getIT() );
          T* fy = fd::accessOld<T,DESCRIPTOR,FIELD>( blockLattice.get ( iX, iY-iN ), this->getIT() );
          T* fY = fd::accessOld<T,DESCRIPTOR,FIELD>( blockLattice.get ( iX, iY+iN ), this->getIT() );
          f[iN-1][0] = *fx;
          F[iN-1][0] = *fX;
          f[iN-1][1] = *fy;
          F[iN-1][1] = *fY;
        }
        Cell<T,DESCRIPTOR> cell = blockLattice.get(iX,iY);
        _model->operator()(fNew, f0, f, F, cell);
        this->applySourceTerm(fNew, cell);
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
void FdPostProcessor2D<T,DESCRIPTOR,FIELD,SOURCE>::
process(BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, _x0, _x1, _y0, _y1);
}


//////// FdPostProcessorGenerator2D ///////////////////////////////////

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
FdPostProcessorGenerator2D<T,DESCRIPTOR,FIELD,SOURCE>::FdPostProcessorGenerator2D (
  int x0_, int x1_, int y0_, int y1_, std::size_t& iT, std::shared_ptr<FdModel<T,DESCRIPTOR>> model )
  : FdBasePostProcessorGenerator2D<T,DESCRIPTOR,FIELD,SOURCE>(x0_, x1_, y0_, y1_, iT),
    _model(model)
{ }

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
FdPostProcessorGenerator2D<T,DESCRIPTOR,FIELD,SOURCE>::FdPostProcessorGenerator2D(size_t& iT, std::shared_ptr<FdModel<T,DESCRIPTOR>> model)
  : FdPostProcessorGenerator2D<T,DESCRIPTOR,FIELD,SOURCE>(0, 0, 0, 0, iT, model)
{ }

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
PostProcessor2D<T,DESCRIPTOR>* FdPostProcessorGenerator2D<T,DESCRIPTOR,FIELD,SOURCE>::generate() const
{
  return new FdPostProcessor2D<T,DESCRIPTOR,FIELD,SOURCE>(this->x0,this->x1,this->y0,this->y1,
         this->_iT,_model);
}

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
PostProcessorGenerator2D<T,DESCRIPTOR>* FdPostProcessorGenerator2D<T,DESCRIPTOR,FIELD,SOURCE>::clone() const
{
  return new FdPostProcessorGenerator2D<T,DESCRIPTOR,FIELD,SOURCE>(*this);
}


}  // namespace olb

#endif
