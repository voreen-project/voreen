/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2020 Alexander Schulz
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

#include "setLocalVelocityBoundary3D.h"
#include "setLocalVelocityBoundary3D.hh"
#include "setInterpolatedVelocityBoundary3D.h"
#include "setInterpolatedVelocityBoundary3D.hh"
#include "setLocalPressureBoundary3D.h"
#include "setLocalPressureBoundary3D.hh"
#include "setInterpolatedPressureBoundary3D.h"
#include "setInterpolatedPressureBoundary3D.hh"
#include "setSlipBoundary3D.h"
#include "setSlipBoundary3D.hh"
#include "setPartialSlipBoundary3D.h"
#include "setPartialSlipBoundary3D.hh"
#include "setFreeEnergyWallBoundary3D.h"
#include "setFreeEnergyWallBoundary3D.hh"
#include "setFreeEnergyInletBoundary3D.h"
#include "setFreeEnergyInletBoundary3D.hh"
#include "setFreeEnergyOutletBoundary3D.h"
#include "setFreeEnergyOutletBoundary3D.hh"
#include "setWallFunctionBoundary3D.h"
#include "setWallFunctionBoundary3D.hh"
#include "setLocalConvectionBoundary3D.h"
#include "setLocalConvectionBoundary3D.hh"
#include "setInterpolatedConvectionBoundary3D.h"
#include "setInterpolatedConvectionBoundary3D.hh"
#include "setAdvectionDiffusionConvectionBoundary3D.h"
#include "setAdvectionDiffusionConvectionBoundary3D.hh"
#include "setAdvectionDiffusionTemperatureBoundary3D.h"
#include "setAdvectionDiffusionTemperatureBoundary3D.hh"
//#include "setExtFieldBoundary3D.h"
//#include "setExtFieldBoundary3D.hh"
#include "setZeroDistributionBoundary3D.h"
#include "setZeroDistributionBoundary3D.hh"
#include "setBouzidiVelocityBoundary3D.h"
#include "setBouzidiVelocityBoundary3D.hh"
#include "setBouzidiZeroVelocityBoundary3D.h"
#include "setBouzidiZeroVelocityBoundary3D.hh"
#include "setZouHeVelocityBoundary3D.h"
#include "setZouHeVelocityBoundary3D.hh"
#include "setZouHePressureBoundary3D.h"
#include "setZouHePressureBoundary3D.hh"
#include "setRtlbmDiffuseTemperatureBoundary3D.h"
#include "setRtlbmDiffuseTemperatureBoundary3D.hh"
#include "setRtlbmDiffuseConstTemperatureBoundary3D.h"
#include "setRtlbmDiffuseConstTemperatureBoundary3D.hh"
#include "setRtlbmDirectedTemperatureBoundary3D.h"
#include "setRtlbmDirectedTemperatureBoundary3D.hh"
#include "defineU3D.h"
#include "defineU3D.hh"

#include "boundary3D.h"

namespace olb {

template void setLocalVelocityBoundary<double, descriptors::D3Q19<>, RLBdynamics<double, descriptors::D3Q19<>>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, double omega, SuperGeometry3D<double>& superGeometry, int material);

template void setLocalVelocityBoundary<double, descriptors::D3Q19<>, RLBdynamics<double, descriptors::D3Q19<>>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, double omega, FunctorPtr<SuperIndicatorF3D<double>>&& indicator);

template void setLocalVelocityBoundary<double, descriptors::D3Q19<>, RLBdynamics<double, descriptors::D3Q19<>>>
(BlockLatticeStructure3D<double, descriptors::D3Q19<>>& _block, double omega, BlockIndicatorF3D<double>& indicator, bool includeOuterCells);

template void addPoints2CommBC<double, descriptors::D3Q19<>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, FunctorPtr<SuperIndicatorF3D<double>>&& indicator, int _overlap);

template void setBoundary<double, descriptors::D3Q19<>, RLBdynamics<double, descriptors::D3Q19<>>>
(BlockLatticeStructure3D<double, descriptors::D3Q19<>>& _block, double omega, int iX, int iY, int iZ,
		Momenta<double, descriptors::D3Q19<>>* momenta, Dynamics<double,  descriptors::D3Q19<>>* dynamics,
		PostProcessorGenerator3D<double,  descriptors::D3Q19<>>* postProcessor);



template void setInterpolatedVelocityBoundary<double, descriptors::D3Q19<>, BGKdynamics<double, descriptors::D3Q19<>>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, double omega, SuperGeometry3D<double>& superGeometry, int material);

template void setInterpolatedVelocityBoundary<double, descriptors::D3Q19<>, ConstRhoBGKdynamics<double, descriptors::D3Q19<>>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, double omega, SuperGeometry3D<double>& superGeometry, int material);

template void setInterpolatedVelocityBoundary<double, descriptors::D3Q19<>, BGKdynamics<double, descriptors::D3Q19<>>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, double omega, FunctorPtr<SuperIndicatorF3D<double>>&& indicator);

template void setInterpolatedVelocityBoundary<double, descriptors::D3Q19<>, BGKdynamics<double, descriptors::D3Q19<>>>
(BlockLatticeStructure3D<double, descriptors::D3Q19<>>& _block, double omega, BlockIndicatorF3D<double>& indicator, bool includeOuterCells);



template void setLocalPressureBoundary<double, descriptors::D3Q19<>, RLBdynamics<double, descriptors::D3Q19<>>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, double omega, SuperGeometry3D<double>& superGeometry, int material);

template void setLocalPressureBoundary<double, descriptors::D3Q19<>, RLBdynamics<double, descriptors::D3Q19<>>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, double omega, FunctorPtr<SuperIndicatorF3D<double>>&& indicator);

template void setLocalPressureBoundary<double, descriptors::D3Q19<>, RLBdynamics<double, descriptors::D3Q19<>>>
(BlockLatticeStructure3D<double, descriptors::D3Q19<>>& _block, double omega, BlockIndicatorF3D<double>& indicator, bool includeOuterCells);



template void setInterpolatedPressureBoundary<double, descriptors::D3Q19<>, BGKdynamics<double, descriptors::D3Q19<>>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, double omega, SuperGeometry3D<double>& superGeometry, int material);

template void setInterpolatedPressureBoundary<double, descriptors::D3Q19<>, BGKdynamics<double, descriptors::D3Q19<>>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, double omega, FunctorPtr<SuperIndicatorF3D<double>>&& indicator);

template void setInterpolatedPressureBoundary<double, descriptors::D3Q19<>, BGKdynamics<double, descriptors::D3Q19<>>>
(BlockLatticeStructure3D<double, descriptors::D3Q19<>>& _block, double omega, BlockIndicatorF3D<double>& indicator, bool includeOuterCells);



template void setSlipBoundary<double, descriptors::D3Q19<>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, SuperGeometry3D<double>& superGeometry, int material);

template void setSlipBoundary<double, descriptors::D3Q19<>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, FunctorPtr<SuperIndicatorF3D<double>>&& indicator);

template void setSlipBoundary<double, descriptors::D3Q19<>>
(BlockLatticeStructure3D<double, descriptors::D3Q19<>>& block, BlockIndicatorF3D<double>& indicator, bool includeOuterCells);



template void setPartialSlipBoundary<double, descriptors::D3Q19<>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, double tuner, SuperGeometry3D<double>& superGeometry, int material);

template void setPartialSlipBoundary<double, descriptors::D3Q19<>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, double tuner, FunctorPtr<SuperIndicatorF3D<double>>&& indicator);

template void setPartialSlipBoundary<double, descriptors::D3Q19<>>
(BlockLatticeStructure3D<double, descriptors::D3Q19<>>& block, double tuner, BlockIndicatorF3D<double>& indicator, bool includeOuterCells);



template void setFreeEnergyWallBoundary<double, descriptors::D3Q19<>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, SuperGeometry3D<double>& superGeometry, int material, double alpha, double kappa1,
		double kappa2, double h1, double h2, int latticeNumber);

template void setFreeEnergyWallBoundary<double, descriptors::D3Q19<>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, FunctorPtr<SuperIndicatorF3D<double>>&& indicator, double alpha, double kappa1,
		double kappa2, double h1, double h2, int latticeNumber);

template void setFreeEnergyWallBoundary<double, descriptors::D3Q19<>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, SuperGeometry3D<double>& superGeometry, int material, double alpha, double kappa1,
		double kappa2, double kappa3, double h1, double h2, double h3, int latticeNumber);

template void setFreeEnergyWallBoundary<double, descriptors::D3Q19<>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, FunctorPtr<SuperIndicatorF3D<double>>&& indicator, double alpha, double kappa1,
		double kappa2, double kappa3, double h1, double h2, double h3, int latticeNumber);



template void setFreeEnergyInletBoundary<double, descriptors::D3Q19<>, RLBdynamics<double, descriptors::D3Q19<>>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, double omega, SuperGeometry3D<double>& superGeometry, int material,
		std::string type, int latticeNumber);

template void setFreeEnergyInletBoundary<double, descriptors::D3Q19<>, RLBdynamics<double, descriptors::D3Q19<>>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, double omega, FunctorPtr<SuperIndicatorF3D<double>>&& indicator, std::string type,
		int latticeNumber);

template void setFreeEnergyInletBoundary<double, descriptors::D3Q19<>, RLBdynamics<double, descriptors::D3Q19<>>>
(BlockLatticeStructure3D<double, descriptors::D3Q19<>>& _block, double omega, BlockIndicatorF3D<double>& indicator, std::string type,
		int latticeNumber, bool includeOuterCells);



template void setFreeEnergyOutletBoundary<double, descriptors::D3Q19<>, RLBdynamics<double, descriptors::D3Q19<>>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, double omega, SuperGeometry3D<double>& superGeometry, int material,
	  std::string type, int latticeNumber);

template void setFreeEnergyOutletBoundary<double, descriptors::D3Q19<>, RLBdynamics<double, descriptors::D3Q19<>>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, double omega, FunctorPtr<SuperIndicatorF3D<double>>&& indicator,
    std::string type, int latticeNumber);

template void setFreeEnergyOutletBoundary<double, descriptors::D3Q19<>, RLBdynamics<double, descriptors::D3Q19<>>>
(BlockLatticeStructure3D<double, descriptors::D3Q19<>>& _block, BlockIndicatorF3D<double>& indicator, double omega, std::string type,
		int latticeNumber, bool includeOuterCells);



//template void setWallFunctionBoundary<double, descriptors::D3Q19<>, RLBdynamics<double, descriptors::D3Q19<>>>
//(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, SuperGeometry3D<double>& superGeometry, int material,
//		UnitConverter<double, descriptors::D3Q19<>> const& converter, wallFunctionParam<double> const& wallFunctionParam,
//		IndicatorF3D<double>* geoIndicator);
//
//template void setWallFunctionBoundary<double, descriptors::D3Q19<>, RLBdynamics<double, descriptors::D3Q19<>>>
//(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, FunctorPtr<SuperIndicatorF3D<double>>&& indicator,
//		UnitConverter<double, descriptors::D3Q19<>> const& converter, wallFunctionParam<double> const& wallFunctionParam,
//		IndicatorF3D<double>* geoIndicator);
//
//template void setWallFunctionBoundary<double, descriptors::D3Q19<>, RLBdynamics<double, descriptors::D3Q19<>>>
//(BlockLatticeStructure3D<double, descriptors::D3Q19<>>& _block, BlockIndicatorF3D<double>& indicator,
//		UnitConverter<double, descriptors::D3Q19<>> const& converter,  wallFunctionParam<double> const& wallFunctionParam,
//		IndicatorF3D<double>* geoIndicator, bool includeOuterCells);



template void setLocalConvectionBoundary<double, descriptors::D3Q19<>, RLBdynamics<double, descriptors::D3Q19<>>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, double omega, SuperGeometry3D<double>& superGeometry, int material,
		double* uAv);

template void setLocalConvectionBoundary<double, descriptors::D3Q19<>, RLBdynamics<double, descriptors::D3Q19<>>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, double omega, FunctorPtr<SuperIndicatorF3D<double>>&& indicator,
		double* uAv);

template void setLocalConvectionBoundary<double, descriptors::D3Q19<>, RLBdynamics<double, descriptors::D3Q19<>>>
(BlockLatticeStructure3D<double, descriptors::D3Q19<>>& _block, double omega, BlockIndicatorF3D<double>& indicator,
		double* uAv, bool includeOuterCells);



template void setInterpolatedConvectionBoundary<double, descriptors::D3Q19<>, BGKdynamics<double, descriptors::D3Q19<>>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, double omega, SuperGeometry3D<double>& superGeometry, int material,
		double* uAv);

template void setInterpolatedConvectionBoundary<double, descriptors::D3Q19<>, BGKdynamics<double, descriptors::D3Q19<>>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, double omega, FunctorPtr<SuperIndicatorF3D<double>>&& indicator,
		double* uAv);

template void setInterpolatedConvectionBoundary<double, descriptors::D3Q19<>, BGKdynamics<double, descriptors::D3Q19<>>>
(BlockLatticeStructure3D<double, descriptors::D3Q19<>>& _block, double omega, BlockIndicatorF3D<double>& indicator,
		double* uAv, bool includeOuterCells);



//template void setAdvectionDiffusionConvectionBoundary<double, descriptors::D3Q19<>, AdvectionDiffusionRLBdynamics<double, descriptors::D3Q19<>>>
//(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, SuperGeometry3D<double>& superGeometry, int material);
//
//template void setAdvectionDiffusionConvectionBoundary<double, descriptors::D3Q19<>, AdvectionDiffusionRLBdynamics<double, descriptors::D3Q19<>>>
//(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, FunctorPtr<SuperIndicatorF3D<double>>&& indicator);
//
//template void setAdvectionDiffusionConvectionBoundary<double, descriptors::D3Q19<>, AdvectionDiffusionRLBdynamics<double, descriptors::D3Q19<>>>
//(BlockLatticeStructure3D<double, descriptors::D3Q19<>>& _block, BlockIndicatorF3D<double>& indicator, bool includeOuterCells);
//
//
//
//template void setAdvectionDiffusionTemperatureBoundary<double, descriptors::D3Q19<>, AdvectionDiffusionRLBdynamics<double, descriptors::D3Q19<>>>
//(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, double omega, SuperGeometry3D<double>& superGeometry, int material);
//
//template void setAdvectionDiffusionTemperatureBoundary<double, descriptors::D3Q19<>, AdvectionDiffusionRLBdynamics<double, descriptors::D3Q19<>>>
//(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, double omega, FunctorPtr<SuperIndicatorF3D<double>>&& indicator);
//
//template void setAdvectionDiffusionTemperatureBoundary<double, descriptors::D3Q19<>, AdvectionDiffusionRLBdynamics<double, descriptors::D3Q19<>>>
//(BlockLatticeStructure3D<double, descriptors::D3Q19<>>& _block, BlockIndicatorF3D<double>& indicator, double omega, bool includeOuterCells);



//template void setExtFieldBoundary<double, descriptors::D3Q19<>, AdvectionDiffusionRLBdynamics<double, descriptors::D3Q19<>>>
//(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, SuperGeometry3D<double>& superGeometry, int material);
//
//template void setExtFieldBoundary<double, descriptors::D3Q19<>, AdvectionDiffusionRLBdynamics<double, descriptors::D3Q19<>>>
//(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, FunctorPtr<SuperIndicatorF3D<double>>&& indicator);
//
//template void setExtFieldBoundary<double, descriptors::D3Q19<>,  AdvectionDiffusionRLBdynamics<double, descriptors::D3Q19<>>>
//(BlockLatticeStructure3D<double, descriptors::D3Q19<>>& _block, BlockIndicatorF3D<double>& indicator, int offset, bool includeOuterCells);



//template void setZeroDistributionBoundary<double, descriptors::D3Q19<>, AdvectionDiffusionRLBdynamics<double, descriptors::D3Q19<>>>
//(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, SuperGeometry3D<double>& superGeometry, int material);
//
//template void setZeroDistributionBoundary<double, descriptors::D3Q19<>, AdvectionDiffusionRLBdynamics<double, descriptors::D3Q19<>>>
//(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, FunctorPtr<SuperIndicatorF3D<double>>&& indicator);
//
//template void setZeroDistributionBoundary<double, descriptors::D3Q19<>, AdvectionDiffusionRLBdynamics<double, descriptors::D3Q19<>>>
//(BlockLatticeStructure3D<double, descriptors::D3Q19<>>& _block, BlockIndicatorF3D<double>& indicator, bool includeOuterCells);



template void setBouzidiVelocityBoundary<double, descriptors::D3Q19<>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, SuperGeometry3D<double>& superGeometry, int material, IndicatorF3D<double>& indicator,
		std::vector<int> bulkMaterials);

template void setBouzidiVelocityBoundary<double, descriptors::D3Q19<>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, FunctorPtr<SuperIndicatorF3D<double>>&& boundaryIndicator,
		FunctorPtr<SuperIndicatorF3D<double>>&& bulkIndicator, IndicatorF3D<double>& geometryIndicator);

template void setBouzidiVelocityBoundary<double, descriptors::D3Q19<>>
(BlockLatticeStructure3D<double, descriptors::D3Q19<>>& block, BlockIndicatorF3D<double>& boundaryIndicator, BlockIndicatorF3D<double>& bulkIndicator,
		IndicatorF3D<double>& geometryIndicator, double _epsFraction);

template void setBouzidiVelocityBoundary<double, descriptors::D3Q19<>>
(BlockLatticeStructure3D<double, descriptors::D3Q19<>>& block, BlockGeometryStructure3D<double>& blockGeometryStructure, int iX, int iY, int iZ,
		IndicatorF3D<double>& geometryIndicator, BlockIndicatorF3D<double>& bulkIndicator, double _epsFraction);

template void setBouzidiVelocityBoundary<double, descriptors::D3Q19<>>
(BlockLatticeStructure3D<double, descriptors::D3Q19<>>& block, BlockGeometryStructure3D<double>& blockGeometryStructure, int x, int y, int z,
		double distances[descriptors::D3Q19<>::q]);



template void setBouzidiZeroVelocityBoundary<double, descriptors::D3Q19<>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, SuperGeometry3D<double>& superGeometry, int material, IndicatorF3D<double>& geometryIndicator,
		std::vector<int> bulkMaterials);

template void setBouzidiZeroVelocityBoundary<double, descriptors::D3Q19<>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, FunctorPtr<SuperIndicatorF3D<double>>&& boundaryIndicator,
		IndicatorF3D<double>& geometryIndicator, std::vector<int> bulkMaterials);

template void setBouzidiZeroVelocityBoundary<double, descriptors::D3Q19<>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, FunctorPtr<SuperIndicatorF3D<double>>&& boundaryIndicator,
		FunctorPtr<SuperIndicatorF3D<double>>&& bulkIndicator, IndicatorF3D<double>& geometryIndicator);

template void setBouzidiZeroVelocityBoundary<double, descriptors::D3Q19<>>
(BlockLatticeStructure3D<double, descriptors::D3Q19<>>& block, BlockIndicatorF3D<double>& boundaryIndicator, BlockIndicatorF3D<double>& bulkIndicator,
		IndicatorF3D<double>& geometryIndicator, double _epsFraction);

template void setBouzidiZeroVelocityBoundary1<double, descriptors::D3Q19<>>
(BlockLatticeStructure3D<double, descriptors::D3Q19<>>& block, BlockGeometryStructure3D<double>& blockGeometryStructure, int iX, int iY, int iZ,
		IndicatorF3D<double>& geometryIndicator, BlockIndicatorF3D<double>& bulkIndicator, double _epsFraction);

template void setBouzidiZeroVelocityBoundary<double, descriptors::D3Q19<>>
(BlockLatticeStructure3D<double, descriptors::D3Q19<>>& block, BlockGeometryStructure3D<double>& blockGeometryStructure, int x, int y, int z,
		double distances[descriptors::D3Q19<>::q]);

template void setBouzidiZeroVelocityBoundary<double, descriptors::D3Q19<>>
(BlockLatticeStructure3D<double, descriptors::D3Q19<>>& block, BlockGeometryStructure3D<double>& blockGeometryStructure, int x, int y, int z, int iPop,
		double dist);


template void setZouHeVelocityBoundary<double, descriptors::D3Q19<>, BGKdynamics<double, descriptors::D3Q19<>>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, double omega, SuperGeometry3D<double>& superGeometry, int material);

template void setZouHeVelocityBoundary<double, descriptors::D3Q19<>, BGKdynamics<double, descriptors::D3Q19<>>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, double omega, FunctorPtr<SuperIndicatorF3D<double>>&& indicator);

template void setZouHeVelocityBoundary<double, descriptors::D3Q19<>, BGKdynamics<double, descriptors::D3Q19<>>>
(BlockLatticeStructure3D<double, descriptors::D3Q19<>>& _block, double omega, BlockIndicatorF3D<double>& indicator, bool includeOuterCells);




template void setZouHePressureBoundary<double, descriptors::D3Q19<>, BGKdynamics<double, descriptors::D3Q19<>>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, double omega, SuperGeometry3D<double>& superGeometry, int material);

template void setZouHePressureBoundary<double, descriptors::D3Q19<>, BGKdynamics<double, descriptors::D3Q19<>>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, double omega, FunctorPtr<SuperIndicatorF3D<double>>&& indicator);

template void setZouHePressureBoundary<double, descriptors::D3Q19<>, BGKdynamics<double, descriptors::D3Q19<>>>
(BlockLatticeStructure3D<double, descriptors::D3Q19<>>& _block, double omega, BlockIndicatorF3D<double>& indicator, bool includeOuterCells);



template void setRtlbmDiffuseTemperatureBoundary<double, descriptors::D3Q19<>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, double omega, SuperGeometry3D<double>& superGeometry, int material);

template void setRtlbmDiffuseTemperatureBoundary<double, descriptors::D3Q19<>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, double omega, FunctorPtr<SuperIndicatorF3D<double>>&& indicator);

template void setRtlbmDiffuseTemperatureBoundary<double, descriptors::D3Q19<>>
(BlockLatticeStructure3D<double, descriptors::D3Q19<>>& _block, BlockIndicatorF3D<double>& indicator, double omega, bool includeOuterCells);



template void setRtlbmDiffuseConstTemperatureBoundary<double, descriptors::D3Q19<>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, double omega, SuperGeometry3D<double>& superGeometry, int material);

template void setRtlbmDiffuseConstTemperatureBoundary<double, descriptors::D3Q19<>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, double omega, FunctorPtr<SuperIndicatorF3D<double>>&& indicator);

template void setRtlbmDiffuseConstTemperatureBoundary<double, descriptors::D3Q19<>>
(BlockLatticeStructure3D<double, descriptors::D3Q19<>>& _block, BlockIndicatorF3D<double>& indicator, double omega, bool includeOuterCells);



template void setRtlbmDirectedTemperatureBoundary<double, descriptors::D3Q19<>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, double omega, SuperGeometry3D<double>& superGeometry, int material);

template void setRtlbmDirectedTemperatureBoundary<double, descriptors::D3Q19<>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, double omega, FunctorPtr<SuperIndicatorF3D<double>>&& indicator);

template void setRtlbmDirectedTemperatureBoundary<double, descriptors::D3Q19<>>
(BlockLatticeStructure3D<double, descriptors::D3Q19<>>& _block, BlockIndicatorF3D<double>& indicator, double omega, bool includeOuterCells);



template void defineUBouzidi<double, descriptors::D3Q19<>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, SuperGeometry3D<double>& superGeometry, int material, AnalyticalF3D<double, double>& u, std::vector<int> bulkMaterials);

template void defineUBouzidi<double, descriptors::D3Q19<>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, FunctorPtr<SuperIndicatorF3D<double>>&& boundaryIndicator, AnalyticalF3D<double, double>& u, std::vector<int> bulkMaterials);

template void defineUBouzidi<double, descriptors::D3Q19<>>
(SuperLattice3D<double, descriptors::D3Q19<>>& sLattice, FunctorPtr<SuperIndicatorF3D<double>>&& boundaryIndicator, FunctorPtr<SuperIndicatorF3D<double>>&& bulkIndicator, AnalyticalF3D<double, double>& u);

template void defineUBouzidi<double, descriptors::D3Q19<>>
(BlockLatticeStructure3D<double, descriptors::D3Q19<>>& _block, BlockIndicatorF3D<double>& indicator, BlockIndicatorF3D<double>& bulkIndicator, AnalyticalF3D<double, double>& u);

template void defineUBouzidi<double, descriptors::D3Q19<>>
(BlockLatticeStructure3D<double, descriptors::D3Q19<>>& _block, int iX, int iY, int iZ, int iPop, const double u[descriptors::D3Q19<>::d]);

}  // namespace olb

