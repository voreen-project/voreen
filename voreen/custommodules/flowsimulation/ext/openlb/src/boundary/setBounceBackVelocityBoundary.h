#ifndef SET_BOUNCE_BACK_VELOCITY_BOUNDARY_H
#define SET_BOUNCE_BACK_VELOCITY_BOUNDARY_H

#include <vector>
#include "utilities/functorPtr.h"
#include "extendedFiniteDifferenceBoundary3D.h"
#include "geometry/superGeometry3D.h"
#include "extendedFiniteDifferenceBoundary3D.h"
#include "core/superLattice3D.h"
#include "functors/lattice/indicator/superIndicatorBaseF3D.h"
#include "dynamics/dynamics.h"
#include "functors/lattice/indicator/blockIndicatorF3D.h"
#include "momentaOnBoundaries3D.h"
#include "io/ostreamManager.h"
#include "functors/lattice/indicator/blockIndicatorF3D.h"
#include "dynamics/freeEnergyDynamics.h"

namespace olb {
//1st setVeloctiyBoundary member-function out of sOnLatticeBoundaryCondition
//eventually MixinDynamics=RH? Dynamics question
template<typename T, typename DESCRIPTOR, typename MixinDynamics=BGKdynamics<T,DESCRIPTOR>>
void setBounceBackVelocityBoundary(SuperGeometry3D<T>& superGeometry, int
                                   material, T omega, SuperLattice3D<T, DESCRIPTOR>& sLattice);

//2nd setLocalVeloctiyBoundary Member function out of sOnLatticeBoundaryCondition
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setBounceBackVelocityBoundary(FunctorPtr<SuperIndicatorF3D<T>>&& indicator, T omega,SuperLattice3D<T, DESCRIPTOR>& sLattice);


//1st setLocalVelocityBoundary member-function out of BoundaryConditionInstantiator
template<typename T, typename DESCRIPTOR, typename MixinDynamics>
void setBounceBackVelocityBoundary(BlockIndicatorF3D<T>& indicator, T omega, bool includeOuterCells,BlockLatticeStructure3D<T,DESCRIPTOR>& _block);

//call of setLocalVelocityBoundary out of BoundaryInstantiator
/*template<int direction, int orientation,typename T, typename DESCRIPTOR, typename MixinDynamics>
void setBounceBackVelocityBoundary(int iX,int iY,int iZ, T omega, BlockLatticeStructure3D<T,DESCRIPTOR>& _block,
    Momenta<T,DESCRIPTOR>* momenta,Dynamics<T,DESCRIPTOR>* dynamics,PostProcessorGenerator3D<T,DESCRIPTOR>* postProcessor);*/
}

#endif
