#ifndef SET_BOUNCE_BACK_VELOCITY_BOUNDARY_HH
#define SET_BOUNCE_BACK_VELOCITY_BOUNDARY_HH

#include "setLocalVelocityBoundary3D.h"

namespace olb {
template<typename T,typename DESCRIPTOR, class MixinDynamics>
void setBounceBackVelocityBoundary(SuperGeometry3D<T>& superGeometry, int material, T omega, SuperLattice3D<T, DESCRIPTOR>& sLattice)
{
  //Arguments with sLattice adjusted
  setBounceBackVelocityBoundary<T,DESCRIPTOR,MixinDynamics>(superGeometry.getMaterialIndicator(material),omega,sLattice);
}
//2nd addLocalVeloctiyBoundary Member function out of sOnLatticeBoundaryCondition
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setBounceBackVelocityBoundary(FunctorPtr<SuperIndicatorF3D<T>>&& indicator, T omega,SuperLattice3D<T, DESCRIPTOR>& sLattice)
{
  OstreamManager clout(std::cout, "BounceBackVelocityBoundary");
  int _overlap = indicator->getSuperGeometry().getOverlap();
  bool includeOuterCells = false;
  if (indicator->getSuperGeometry().getOverlap() == 1) {
    includeOuterCells = true;
    clout << "WARNING: overlap == 1, boundary conditions set on overlap despite unknown neighbor materials" << std::endl;
  }
  // addPoints2CommBC(sLattice,std::forward<decltype(indicator)>(indicator), _overlap);
  clout << sLattice.getLoadBalancer().size() <<"sLattice.getLoadBalancer.size()" << std::endl;
  for (int iC = 0; iC < sLattice.getLoadBalancer().size(); ++iC) {
    setBounceBackVelocityBoundary<T,DESCRIPTOR,MixinDynamics>(indicator->getExtendedBlockIndicatorF(iC),omega,includeOuterCells,sLattice.getExtendedBlockLattice(iC));
  }
}


//Maybe this is the only function needed to distinguish between local and interpolated because here the interpolated or local boundary manager would've been called

//1st addLocalVelociyBoundary member-function out of OnLatticeBoundaryCondition3d
template<typename T, typename DESCRIPTOR, class MixinDynamics>
void setBounceBackVelocityBoundary(BlockIndicatorF3D<T>& indicator, T omega, bool includeOuterCells,BlockLatticeStructure3D<T,DESCRIPTOR>& _block)
{
  auto& blockGeometryStructure = indicator.getBlockGeometryStructure();
  const int margin = includeOuterCells ? 0 : 1;
  int x0 = margin;
  int y0 = margin;
  int z0 = margin;
  int x1 = blockGeometryStructure.getNx()-1 -margin;
  int y1 = blockGeometryStructure.getNy()-1 -margin;
  int z1 = blockGeometryStructure.getNz()-1 -margin;
  std::vector<int> discreteNormal(4,0);
  T default_rho = 1.0;
  T default_u[] = {0.0, 0.0, 0.0};
  for (int iX = x0; iX <= x1; iX++) {
    for (int iY = y0; iY <= y1; iY++) {
      for (int iZ = z0; iZ <= z1; iZ++) {
        Momenta<T,DESCRIPTOR>* momenta = nullptr;
        Dynamics<T,DESCRIPTOR>* dynamics = nullptr;
        PostProcessorGenerator3D<T,DESCRIPTOR>* postProcessor = nullptr;
        if (indicator(iX, iY, iZ)) {
          momenta = nullptr;
          dynamics = new BounceBackVelocity<T,DESCRIPTOR>(default_rho, default_u);
          postProcessor = nullptr;
          //Defined in local velocity boundary
          setBoundary<T, DESCRIPTOR, MixinDynamics>(_block, iX,iY,iZ, omega, momenta, dynamics, postProcessor);

        }
      }
    }
  }

}

}
#endif
