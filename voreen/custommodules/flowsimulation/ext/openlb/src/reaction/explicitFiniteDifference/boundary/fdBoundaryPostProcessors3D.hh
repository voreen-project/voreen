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

#ifndef FD_NO_PENETRATION_BOUNDARY_POST_PROCESSORS_3D_HH
#define FD_NO_PENETRATION_BOUNDARY_POST_PROCESSORS_3D_HH

namespace olb {


////////  FdBoundaryPostProcessor3D ///////////////////////////////////

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
FdBoundaryPostProcessor3D <T,DESCRIPTOR,FIELD,SOURCE>::
FdBoundaryPostProcessor3D ( int x0, int x1, int y0, int y1, int z0, int z1, std::size_t& iT, int normalX, int normalY, int normalZ,
    std::shared_ptr<FdModel<T,DESCRIPTOR>> model, std::shared_ptr<fd::AdBoundarySchemeBase<3,T>> boundaryScheme )
  :  FdBasePostProcessor3D<T,DESCRIPTOR,FIELD,SOURCE>(iT),
     _x0(x0), _x1(x1), _y0(y0), _y1(y1), _z0(z0), _z1(z1),
     _normalX(normalX), _normalY(normalY), _normalZ(normalZ),
     _model(model),
     _boundaryScheme(boundaryScheme),
     f{this->initMatrix(_model->extent())},
     F{this->initMatrix(_model->extent())},
     fNormal{this->initMatrix(_model->extent()+_boundaryScheme->getExtraExtent())},
     fGhost{this->initMatrix(_model->extent())}
{ }

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
FdBoundaryPostProcessor3D <T,DESCRIPTOR,FIELD,SOURCE>::
FdBoundaryPostProcessor3D ( std::size_t& iT, int normalX, int normalY, int normalZ,
    std::shared_ptr<FdModel<T,DESCRIPTOR>> model, std::shared_ptr<fd::AdBoundarySchemeBase<3,T>> boundaryScheme )
  :  FdBoundaryPostProcessor3D<T,DESCRIPTOR,FIELD,SOURCE>(0,0,0,0,0,0,iT,normalX,normalY,normalZ,model,boundaryScheme)
{ }

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
void FdBoundaryPostProcessor3D<T,DESCRIPTOR,FIELD,SOURCE>::
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
          std::vector<int> normal {_normalX, _normalY, _normalZ};
          //int normal[] {_normalX, _normalY, _normalX}; // <--- wrong??
          T* fNew = fd::accessNew<T,DESCRIPTOR,FIELD>( blockLattice.get(iX, iY, iZ            ), this->getIT() );
          T* f0   = fd::accessOld<T,DESCRIPTOR,FIELD>( blockLattice.get(iX, iY, iZ            ), this->getIT() );
          for (int iN=1; iN<=_model->extent()+_boundaryScheme->getExtraExtent(); ++iN) {
            T* fnX= fd::accessOld<T,DESCRIPTOR,FIELD>( blockLattice.get(iX+_normalX*iN, iY, iZ            ), this->getIT() );
            T* fnY= fd::accessOld<T,DESCRIPTOR,FIELD>( blockLattice.get(iX, iY+_normalY*iN, iZ            ), this->getIT() );
            T* fnZ= fd::accessOld<T,DESCRIPTOR,FIELD>( blockLattice.get(iX, iY, iZ+_normalZ*iN), this->getIT() );
            fNormal[iN-1][0] = *fnX;
            fNormal[iN-1][1] = *fnY;
            fNormal[iN-1][2] = *fnZ;
          }
          T uArr[3];
          Cell<T,DESCRIPTOR> cell = blockLattice.get(iX,iY,iZ);
          cell.computeU(uArr);
          std::vector<T> u {uArr[0], uArr[1], uArr[2]};
          _boundaryScheme->operator()(fGhost, *f0, fNormal, normal, u);
          for (int iN=1; iN<=_model->extent(); ++iN) {
            T* fx = fd::accessOld<T,DESCRIPTOR,FIELD>( blockLattice.get(iX-iN, iY, iZ            ), this->getIT() );
            T* fX = fd::accessOld<T,DESCRIPTOR,FIELD>( blockLattice.get(iX+iN, iY, iZ            ), this->getIT() );
            T* fy = fd::accessOld<T,DESCRIPTOR,FIELD>( blockLattice.get(iX, iY-iN, iZ            ), this->getIT() );
            T* fY = fd::accessOld<T,DESCRIPTOR,FIELD>( blockLattice.get(iX, iY+iN, iZ            ), this->getIT() );
            T* fz = fd::accessOld<T,DESCRIPTOR,FIELD>( blockLattice.get(iX, iY, iZ-iN         ), this->getIT() );
            T* fZ = fd::accessOld<T,DESCRIPTOR,FIELD>( blockLattice.get(iX, iY, iZ+iN         ), this->getIT() );
            f[iN-1][0] = _normalX<=0 ? *fx : fGhost[iN-1][0];
            F[iN-1][0] = _normalX>=0 ? *fX : fGhost[iN-1][0];
            f[iN-1][1] = _normalY<=0 ? *fy : fGhost[iN-1][1];
            F[iN-1][1] = _normalY>=0 ? *fY : fGhost[iN-1][1];
            f[iN-1][2] = _normalZ<=0 ? *fz : fGhost[iN-1][2];
            F[iN-1][2] = _normalZ>=0 ? *fZ : fGhost[iN-1][2];
          }
          _model->operator()(fNew, f0, f, F, cell);
          this->applySourceTerm(fNew, cell);
        }
      }
    }
  }
}

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
void FdBoundaryPostProcessor3D<T,DESCRIPTOR,FIELD,SOURCE>::
process(BlockLattice<T,DESCRIPTOR>& blockLattice)
{
  processSubDomain(blockLattice, _x0, _x1, _y0, _y1, _z0, _z1);
}


//////// FdBaseBoundaryPostProcessorGenerator3D ///////////////////////////////////

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
FdBaseBoundaryPostProcessorGenerator3D<T,DESCRIPTOR,FIELD,SOURCE>::FdBaseBoundaryPostProcessorGenerator3D (
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, std::size_t& iT, int normalX, int normalY, int normalZ,
  std::shared_ptr<FdModel<T,DESCRIPTOR>> model, std::shared_ptr<fd::AdBoundarySchemeBase<3,T>> boundaryScheme )
  : FdBasePostProcessorGenerator3D<T,DESCRIPTOR,FIELD,SOURCE>(x0_, x1_, y0_, y1_, z0_, z1_, iT),
    _normalX(normalX), _normalY(normalY), _normalZ(normalZ),
    _model(model), _boundaryScheme(boundaryScheme)
{ }

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
FdBaseBoundaryPostProcessorGenerator3D<T,DESCRIPTOR,FIELD,SOURCE>::
FdBaseBoundaryPostProcessorGenerator3D(size_t& iT, int normalX, int normalY, int normalZ,
  std::shared_ptr<FdModel<T,DESCRIPTOR>> model, std::shared_ptr<fd::AdBoundarySchemeBase<3,T>> boundaryScheme )
  : FdBaseBoundaryPostProcessorGenerator3D<T,DESCRIPTOR,FIELD,SOURCE>(0, 0, 0, 0, 0, 0, iT, normalX, normalY, normalZ, model, boundaryScheme)
{ }

template<typename T, typename DESCRIPTOR, typename FIELD, typename SOURCE>
PostProcessor3D<T,DESCRIPTOR>* FdBaseBoundaryPostProcessorGenerator3D<T,DESCRIPTOR,FIELD,SOURCE>::generate() const
{
  return new FdBoundaryPostProcessor3D<T,DESCRIPTOR,FIELD,SOURCE>(this->x0,this->x1,this->y0,this->y1,this->z0,this->z1,
         this->_iT,this->_normalX,this->_normalY,this->_normalZ,this->_model,this->_boundaryScheme);
}


//////// FdNeumannZeroBoundaryPostProcessorGenerator3D ///////////////////////////////////

template<typename T, typename DESCRIPTOR, typename SCHEME_ADV, typename FIELD, typename SOURCE>
FdNeumannZeroBoundaryPostProcessorGenerator3D<T,DESCRIPTOR,SCHEME_ADV,FIELD,SOURCE>::FdNeumannZeroBoundaryPostProcessorGenerator3D (
  int x0_, int x1_, int y0_, int y1_, int z0_, int z1_, std::size_t& iT, int normalX, int normalY, int normalZ, std::shared_ptr<FdModel<T,DESCRIPTOR>> model)
  : FdBaseBoundaryPostProcessorGenerator3D<T,DESCRIPTOR,FIELD,SOURCE>(x0_, x1_, y0_, y1_, z0_, z1_, iT, normalX, normalY, normalZ, model,
    std::make_shared<fd::AdNeumannZeroBoundaryScheme<3,T,SCHEME_ADV>>())
{ }

template<typename T, typename DESCRIPTOR, typename SCHEME_ADV, typename FIELD, typename SOURCE>
FdNeumannZeroBoundaryPostProcessorGenerator3D<T,DESCRIPTOR,SCHEME_ADV,FIELD,SOURCE>::
FdNeumannZeroBoundaryPostProcessorGenerator3D(size_t& iT, int normalX, int normalY, int normalZ, std::shared_ptr<FdModel<T,DESCRIPTOR>> model)
  : FdNeumannZeroBoundaryPostProcessorGenerator3D<T,DESCRIPTOR,SCHEME_ADV,FIELD,SOURCE>(0, 0, 0, 0, 0, 0, iT, normalX, normalY, normalZ, model)
{ }

template<typename T, typename DESCRIPTOR, typename SCHEME_ADV, typename FIELD, typename SOURCE>
PostProcessorGenerator3D<T,DESCRIPTOR>* FdNeumannZeroBoundaryPostProcessorGenerator3D<T,DESCRIPTOR,SCHEME_ADV,FIELD,SOURCE>::clone() const
{
  return new FdNeumannZeroBoundaryPostProcessorGenerator3D<T,DESCRIPTOR,SCHEME_ADV,FIELD,SOURCE>(*this);
}


}  // namespace olb

#endif
