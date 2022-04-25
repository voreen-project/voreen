/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2006, 2007 Jonas Latt, Davide Dapelo
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

#ifndef FD_NO_PENETRATION_BOUNDARIES_2D_HH
#define FD_NO_PENETRATION_BOUNDARIES_2D_HH

namespace olb {


////////  FdBoundaryPostProcessor2D ///////////////////////////////////

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
FdBoundaryPostProcessor2D <T,DESCRIPTOR,FIELD,SOURCE>::
FdBoundaryPostProcessor2D ( int x0, int x1, int y0, int y1, std::size_t& iT, int normalX, int normalY,
    std::shared_ptr<FdModel<T,DESCRIPTOR>> model, std::shared_ptr<fd::AdBoundarySchemeBase<2,T>> boundaryScheme )
  :  FdBasePostProcessor2D<T,DESCRIPTOR,FIELD,SOURCE>(iT),
     _x0(x0), _x1(x1), _y0(y0), _y1(y1),
     _normalX(normalX), _normalY(normalY),
     _model(model),
     _boundaryScheme(boundaryScheme),
     f{this->initMatrix(_model->extent())},
     F{this->initMatrix(_model->extent())},
     fNormal{this->initMatrix(_model->extent()+_boundaryScheme->getExtraExtent())},
     fGhost{this->initMatrix(_model->extent())}
{ }

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
FdBoundaryPostProcessor2D <T,DESCRIPTOR,FIELD,SOURCE>::
FdBoundaryPostProcessor2D ( std::size_t& iT, int normalX, int normalY,
    std::shared_ptr<FdModel<T,DESCRIPTOR>> model, std::shared_ptr<fd::AdBoundarySchemeBase<2,T>> boundaryScheme )
  :  FdBoundaryPostProcessor2D<T,DESCRIPTOR,FIELD,SOURCE>(0,0,0,0,iT,normalX,normalY,model,boundaryScheme)
{ }

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
void FdBoundaryPostProcessor2D<T,DESCRIPTOR,FIELD,SOURCE>::
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
        std::vector<int> normal {_normalX, _normalY};
        T* fNew = fd::accessNew<T,DESCRIPTOR,FIELD>( blockLattice.get ( iX, iY             ), this->getIT() );
        T* f0   = fd::accessOld<T,DESCRIPTOR,FIELD>( blockLattice.get ( iX, iY             ), this->getIT() );
        for (int iN=1; iN<=_model->extent()+_boundaryScheme->getExtraExtent(); ++iN) {
          T* fnX= fd::accessOld<T,DESCRIPTOR,FIELD>( blockLattice.get ( iX+_normalX*iN, iY             ), this->getIT() );
          T* fnY= fd::accessOld<T,DESCRIPTOR,FIELD>( blockLattice.get ( iX, iY+_normalY*iN ), this->getIT() );
          fNormal[iN-1][0] = *fnX;
          fNormal[iN-1][1] = *fnY;
        }
        T uArr[2];
        Cell<T,DESCRIPTOR> cell = blockLattice.get(iX,iY);
        cell.computeU(uArr);
        std::vector<T> u {uArr[0], uArr[1]};
        _boundaryScheme->operator()(fGhost, *f0, fNormal, normal, u);
        for (int iN=1; iN<=_model->extent(); ++iN) {
          T* fx = fd::accessOld<T,DESCRIPTOR,FIELD>( blockLattice.get ( iX-iN, iY             ), this->getIT() );
          T* fX = fd::accessOld<T,DESCRIPTOR,FIELD>( blockLattice.get ( iX+iN, iY             ), this->getIT() );
          T* fy = fd::accessOld<T,DESCRIPTOR,FIELD>( blockLattice.get ( iX, iY-iN          ), this->getIT() );
          T* fY = fd::accessOld<T,DESCRIPTOR,FIELD>( blockLattice.get ( iX, iY+iN          ), this->getIT() );
          f[iN-1][0] = _normalX<=0 ? *fx : fGhost[iN-1][0];
          F[iN-1][0] = _normalX>=0 ? *fX : fGhost[iN-1][0];
          f[iN-1][1] = _normalY<=0 ? *fy : fGhost[iN-1][1];
          F[iN-1][1] = _normalY>=0 ? *fY : fGhost[iN-1][1];
        }
        _model->operator()(fNew, f0, f, F, cell);
        this->applySourceTerm(fNew, cell);
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
void FdBoundaryPostProcessor2D<T,DESCRIPTOR,FIELD,SOURCE>::
process(BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, _x0, _x1, _y0, _y1);
}


//////// FdBaseBoundaryPostProcessorGenerator2D ///////////////////////////////////

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
FdBaseBoundaryPostProcessorGenerator2D<T,DESCRIPTOR,FIELD,SOURCE>::FdBaseBoundaryPostProcessorGenerator2D (
  int x0_, int x1_, int y0_, int y1_, std::size_t& iT, int normalX, int normalY,
  std::shared_ptr<FdModel<T,DESCRIPTOR>> model, std::shared_ptr<fd::AdBoundarySchemeBase<2,T>> boundaryScheme )
  : FdBasePostProcessorGenerator2D<T,DESCRIPTOR,FIELD,SOURCE>(x0_, x1_, y0_, y1_, iT),
    _normalX(normalX), _normalY(normalY),
    _model(model), _boundaryScheme(boundaryScheme)
{ }

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
FdBaseBoundaryPostProcessorGenerator2D<T,DESCRIPTOR,FIELD,SOURCE>::
FdBaseBoundaryPostProcessorGenerator2D(size_t& iT, int normalX, int normalY,
      std::shared_ptr<FdModel<T,DESCRIPTOR>> model, std::shared_ptr<fd::AdBoundarySchemeBase<2,T>> boundaryScheme )
  : FdBaseBoundaryPostProcessorGenerator2D<T,DESCRIPTOR,FIELD,SOURCE>(0, 0, 0, 0, iT, normalX, normalY, model, boundaryScheme)
{ }

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
PostProcessor2D<T,DESCRIPTOR>* FdBaseBoundaryPostProcessorGenerator2D<T,DESCRIPTOR,FIELD,SOURCE>::generate() const
{
  return new FdBoundaryPostProcessor2D<T,DESCRIPTOR,FIELD,SOURCE>(this->x0,this->x1,this->y0,this->y1,
         this->_iT,this->_normalX,this->_normalY,this->_model,this->_boundaryScheme);
}


//////// FdNeumannZeroBoundaryPostProcessorGenerator2D ///////////////////////////////////

template<typename T, typename DESCRIPTOR, typename SCHEME_ADV, typename FIELD, typename SOURCE>
FdNeumannZeroBoundaryPostProcessorGenerator2D<T,DESCRIPTOR,SCHEME_ADV,FIELD,SOURCE>::FdNeumannZeroBoundaryPostProcessorGenerator2D (
  int x0_, int x1_, int y0_, int y1_, std::size_t& iT, int normalX, int normalY, std::shared_ptr<FdModel<T,DESCRIPTOR>> model)
  : FdBaseBoundaryPostProcessorGenerator2D<T,DESCRIPTOR,FIELD,SOURCE>(x0_, x1_, y0_, y1_, iT, normalX, normalY, model,
    std::make_shared<fd::AdNeumannZeroBoundaryScheme<2,T,SCHEME_ADV>>())
{ }

template<typename T, typename DESCRIPTOR, typename SCHEME_ADV, typename FIELD, typename SOURCE>
FdNeumannZeroBoundaryPostProcessorGenerator2D<T,DESCRIPTOR,SCHEME_ADV,FIELD,SOURCE>::
FdNeumannZeroBoundaryPostProcessorGenerator2D(size_t& iT, int normalX, int normalY, std::shared_ptr<FdModel<T,DESCRIPTOR>> model)
  : FdNeumannZeroBoundaryPostProcessorGenerator2D<T,DESCRIPTOR,SCHEME_ADV,FIELD,SOURCE>(0, 0, 0, 0, iT, normalX, normalY, model)
{ }

template<typename T, typename DESCRIPTOR, typename SCHEME_ADV, typename FIELD, typename SOURCE>
PostProcessorGenerator2D<T,DESCRIPTOR>* FdNeumannZeroBoundaryPostProcessorGenerator2D<T,DESCRIPTOR,SCHEME_ADV,FIELD,SOURCE>::clone() const
{
  return new FdNeumannZeroBoundaryPostProcessorGenerator2D<T,DESCRIPTOR,SCHEME_ADV,FIELD,SOURCE>(*this);
}


}  // namespace olb

#endif
