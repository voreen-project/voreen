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

#include "setLocalVelocityBoundary2D.h"
#include "setLocalVelocityBoundary2D.hh"
#include "setInterpolatedVelocityBoundary2D.h"
#include "setInterpolatedVelocityBoundary2D.hh"
#include "setLocalPressureBoundary2D.h"
#include "setLocalPressureBoundary2D.hh"
#include "setInterpolatedPressureBoundary2D.h"
#include "setInterpolatedPressureBoundary2D.hh"
#include "setSlipBoundary2D.h"
#include "setSlipBoundary2D.hh"
#include "setPartialSlipBoundary2D.h"
#include "setPartialSlipBoundary2D.hh"
#include "setFreeEnergyWallBoundary2D.h"
#include "setFreeEnergyWallBoundary2D.hh"
#include "setFreeEnergyInletBoundary2D.h"
#include "setFreeEnergyInletBoundary2D.hh"
#include "setFreeEnergyOutletBoundary2D.h"
#include "setFreeEnergyOutletBoundary2D.hh"
#include "setRegularizedTemperatureBoundary2D.h"
#include "setRegularizedTemperatureBoundary2D.hh"
#include "setAdvectionDiffusionTemperatureBoundary2D.h"
#include "setAdvectionDiffusionTemperatureBoundary2D.hh"
#include "setRegularizedHeatFluxBoundary2D.h"
#include "setRegularizedHeatFluxBoundary2D.hh"
#include "setLocalConvectionBoundary2D.h"
#include "setLocalConvectionBoundary2D.hh"
#include "setInterpolatedConvectionBoundary2D.h"
#include "setInterpolatedConvectionBoundary2D.hh"
#include "setBouzidiVelocityBoundary2D.h"
#include "setBouzidiVelocityBoundary2D.hh"
#include "setBounceBackVelocityBoundary2D.h"
#include "setBounceBackVelocityBoundary2D.hh"
#include "setBouzidiZeroVelocityBoundary2D.h"
#include "setBouzidiZeroVelocityBoundary2D.hh"
#include "setBounceBackZeroVelocityBoundary2D.h"
#include "setBounceBackZeroVelocityBoundary2D.hh"
#include "setZouHeVelocityBoundary2D.h"
#include "setZouHeVelocityBoundary2D.hh"
#include "setZouHePressureBoundary2D.h"
#include "setZouHePressureBoundary2D.hh"


#include "boundary2D.h"

namespace olb {

template void setLocalVelocityBoundary<double,descriptors::D2Q9<>, RLBdynamics<double, descriptors::D2Q9<>>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice,double omega, SuperGeometry2D<double>& superGeometry, int material);

template void setLocalVelocityBoundary<double, descriptors::D2Q9<>, RLBdynamics<double, descriptors::D2Q9<>>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, double omega, FunctorPtr<SuperIndicatorF2D<double>>&& indicator);

template void setLocalVelocityBoundary<double, descriptors::D2Q9<>,  RLBdynamics<double, descriptors::D2Q9<>>>
(BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block, double omega, BlockIndicatorF2D<double>& indicator, bool includeOuterCells);

template void addPoints2CommBC<double,descriptors::D2Q9<>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, FunctorPtr<SuperIndicatorF2D<double>>&& indicator, int _overlap);

template void setBoundary<double,descriptors::D2Q9<>, RLBdynamics<double, descriptors::D2Q9<>>>
(BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block, double omega, int iX, int iY, Momenta<double, descriptors::D2Q9<>>* momenta,
		Dynamics<double, descriptors::D2Q9<>>* dynamics, PostProcessorGenerator2D<double, descriptors::D2Q9<>>* postProcessor);




template void setInterpolatedVelocityBoundary<double, descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice,double omega, SuperGeometry2D<double>& superGeometry, int material);

template void setInterpolatedVelocityBoundary<double, descriptors::D2Q9<>, ConstRhoBGKdynamics<double, descriptors::D2Q9<>>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice,double omega, SuperGeometry2D<double>& superGeometry, int material);

template void setInterpolatedVelocityBoundary<double, descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, double omega, FunctorPtr<SuperIndicatorF2D<double>>&& indicator);

template void setInterpolatedVelocityBoundary<double, descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block, double omega, BlockIndicatorF2D<double>& indicator, bool includeOuterCells);




template void setLocalPressureBoundary<double,descriptors::D2Q9<>, RLBdynamics<double, descriptors::D2Q9<>>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice,double omega, SuperGeometry2D<double>& superGeometry, int material);

template void setLocalPressureBoundary<double,descriptors::D2Q9<>, RLBdynamics<double, descriptors::D2Q9<>>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, double omega, FunctorPtr<SuperIndicatorF2D<double>>&& indicator);

template void setLocalPressureBoundary<double,descriptors::D2Q9<>, RLBdynamics<double, descriptors::D2Q9<>>>
(BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block, double omega, BlockIndicatorF2D<double>& indicator, bool includeOuterCells);



template void setInterpolatedPressureBoundary<double, descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice,double omega, SuperGeometry2D<double>& superGeometry, int material);

template void setInterpolatedPressureBoundary<double, descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, double omega, FunctorPtr<SuperIndicatorF2D<double>>&& indicator);

template void setInterpolatedPressureBoundary<double, descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block, double omega, BlockIndicatorF2D<double>& indicator, bool includeOuterCells);



template void setSlipBoundary<double,descriptors::D2Q9<>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, SuperGeometry2D<double>& superGeometry, int material);

template void setSlipBoundary<double,descriptors::D2Q9<>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, FunctorPtr<SuperIndicatorF2D<double>>&& indicator);

template void setSlipBoundary<double,descriptors::D2Q9<>>
(BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block, BlockIndicatorF2D<double>& indicator, bool includeOuterCells);



template void setPartialSlipBoundary<double,descriptors::D2Q9<>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice,double tuner, SuperGeometry2D<double>& superGeometry, int material);

template void setPartialSlipBoundary<double,descriptors::D2Q9<>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, double tuner, FunctorPtr<SuperIndicatorF2D<double>>&& indicator);

template void setPartialSlipBoundary<double,descriptors::D2Q9<>>
(BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block, double tuner, BlockIndicatorF2D<double>& indicator, bool includeOuterCells);



template void setFreeEnergyWallBoundary<double,descriptors::D2Q9<>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, SuperGeometry2D<double>& superGeometry, int material, double alpha,
		double kappa1, double kappa2, double h1, double h2, int latticeNumber);

template void setFreeEnergyWallBoundary<double,descriptors::D2Q9<>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, FunctorPtr<SuperIndicatorF2D<double>>&& indicator, double alpha,
		double kappa1, double kappa2, double h1, double h2, int latticeNumber);

template void setFreeEnergyWallBoundary<double,descriptors::D2Q9<>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, SuperGeometry2D<double>& superGeometry, int material, double alpha,
		double kappa1, double kappa2, double kappa3, double h1, double h2, double h3, int latticeNumber);

template void setFreeEnergyWallBoundary<double,descriptors::D2Q9<>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, FunctorPtr<SuperIndicatorF2D<double>>&& indicator, double alpha,
		double kappa1, double kappa2, double kappa3, double h1, double h2, double h3, int latticeNumber);

template void setFreeEnergyWallBoundary<double,descriptors::D2Q9<>>
(BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block, BlockIndicatorF2D<double>& indicator, double addend, int latticeNumber,
		bool includeOuterCells);



template void setFreeEnergyInletBoundary<double,descriptors::D2Q9<>, RLBdynamics<double, descriptors::D2Q9<>>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, double omega, SuperGeometry2D<double>& superGeometry, int material,
		std::string type, int latticeNumber);

template void setFreeEnergyInletBoundary<double,descriptors::D2Q9<>, RLBdynamics<double, descriptors::D2Q9<>>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, double omega, FunctorPtr<SuperIndicatorF2D<double>>&& indicator,
		std::string type, int latticeNumber);

template void setFreeEnergyInletBoundary<double,descriptors::D2Q9<>, RLBdynamics<double, descriptors::D2Q9<>>>
(BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block, double omega, BlockIndicatorF2D<double>& indicator,
		std::string type, int latticeNumber, bool includeOuterCells);



template void setFreeEnergyOutletBoundary<double,descriptors::D2Q9<>, RLBdynamics<double, descriptors::D2Q9<>>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, double omega, SuperGeometry2D<double>& superGeometry, int material,
		std::string type, int latticeNumber);

template void setFreeEnergyOutletBoundary<double,descriptors::D2Q9<>, RLBdynamics<double, descriptors::D2Q9<>>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, double omega, FunctorPtr<SuperIndicatorF2D<double>>&& indicator,
		std::string type, int latticeNumber);

template void setFreeEnergyOutletBoundary<double,descriptors::D2Q9<>, RLBdynamics<double, descriptors::D2Q9<>>>
(BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block, double omega, BlockIndicatorF2D<double>& indicator,
		std::string type, int latticeNumber, bool includeOuterCells);



//template void setRegularizedTemperatureBoundary<double,descriptors::D2Q9<>, AdvectionDiffusionRLBdynamics<double, descriptors::D2Q9<>>>
//(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, double omega, SuperGeometry2D<double>& superGeometry, int material);
//
//template void setRegularizedTemperatureBoundary<double,descriptors::D2Q9<>, AdvectionDiffusionRLBdynamics<double, descriptors::D2Q9<>>>
//(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, double omega, FunctorPtr<SuperIndicatorF2D<double>>&& indicator);
//
//template void setRegularizedTemperatureBoundary<double,descriptors::D2Q9<>, AdvectionDiffusionRLBdynamics<double, descriptors::D2Q9<>>>
//(BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block, double omega, BlockIndicatorF2D<double>& indicator);
//
//
//
//template void setTemperatureBoundary<double,descriptors::D2Q9<>, AdvectionDiffusionRLBdynamics<double, descriptors::D2Q9<>>>
//(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, double omega, SuperGeometry2D<double>& superGeometry, int material);
//
//template void setTemperatureBoundary<double,descriptors::D2Q9<>, AdvectionDiffusionRLBdynamics<double, descriptors::D2Q9<>>>
//(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, double omega, FunctorPtr<SuperIndicatorF2D<double>>&& indicator);
//
//template void setTemperatureBoundary<double,descriptors::D2Q9<>, AdvectionDiffusionRLBdynamics<double, descriptors::D2Q9<>>>
//(BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block, double omega, BlockIndicatorF2D<double>& indicator, bool includeOuterCells);
//
//
//
//template void setRegularizedHeatFluxBoundary<double,descriptors::D2Q9<>, AdvectionDiffusionRLBdynamics<double, descriptors::D2Q9<>>>
//(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, double omega, SuperGeometry2D<double>& superGeometry, int material, double *heatFlux);
//
//template void setRegularizedHeatFluxBoundary<double,descriptors::D2Q9<>, AdvectionDiffusionRLBdynamics<double, descriptors::D2Q9<>>>
//(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, double omega, FunctorPtr<SuperIndicatorF2D<double>>&& indicator, double *heatFlux);
//
//template void setRegularizedHeatFluxBoundary<double,descriptors::D2Q9<>, AdvectionDiffusionRLBdynamics<double, descriptors::D2Q9<>>>
//(BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block, double omega, BlockIndicatorF2D<double>& indicator, double *heatFlux);



template void setLocalConvectionBoundary<double,descriptors::D2Q9<>, RLBdynamics<double, descriptors::D2Q9<>>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, double omega, SuperGeometry2D<double>& superGeometry, int material, double* uAv);

template void setLocalConvectionBoundary<double,descriptors::D2Q9<>, RLBdynamics<double, descriptors::D2Q9<>>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, double omega, FunctorPtr<SuperIndicatorF2D<double>>&& indicator, double* uAv);

template void setLocalConvectionBoundary<double,descriptors::D2Q9<>, RLBdynamics<double, descriptors::D2Q9<>>>
(BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block, double omega, BlockIndicatorF2D<double>& indicator, double* uAv,
		bool includeOuterCells);



template void setInterpolatedConvectionBoundary<double,descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, double omega, SuperGeometry2D<double>& superGeometry, int material, double* uAv);

template void setInterpolatedConvectionBoundary<double,descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, double omega, FunctorPtr<SuperIndicatorF2D<double>>&& indicator, double* uAv);

template void setInterpolatedConvectionBoundary<double,descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block, double omega, BlockIndicatorF2D<double>& indicator, double* uAv,
		bool includeOuterCells);



template void setBouzidiVelocityBoundary<double,descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, SuperGeometry2D<double>& superGeometry, int material,
		IndicatorF2D<double>& geometryIndicator, std::vector<int> bulkMaterials);

template void setBouzidiVelocityBoundary<double,descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, FunctorPtr<SuperIndicatorF2D<double>>&& boundaryIndicator,
		FunctorPtr<SuperIndicatorF2D<double>>&& bulkIndicator, IndicatorF2D<double>& geometryIndicator);

template void setBouzidiVelocityBoundary<double,descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block, BlockIndicatorF2D<double>& boundaryIndicator,
		BlockIndicatorF2D<double>& bulkIndicator, IndicatorF2D<double>& geometryIndicator);

template void setBouzidiVelocityBoundary<double,descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block, BlockGeometryStructure2D<double>& blockGeometryStructure,
		int iX, int iY, BlockIndicatorF2D<double>& bulkIndicator, IndicatorF2D<double>& geometryIndicator);

template void setBouzidiVelocityBoundary<double,descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block, BlockGeometryStructure2D<double>& blockGeometryStructure,
		int x, int y, double distances[descriptors::D2Q9<>::q]);

template void setBouzidiVelocityBoundary<double,descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block, BlockGeometryStructure2D<double>& blockGeometryStructure,
		int x, int y, int iPop, double dist);

template void setOffDynamics<double,descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block, int x, int y, double location[descriptors::D2Q9<>::d],
		double distances[descriptors::D2Q9<>::q]);



template void setBounceBackVelocityBoundary<double,descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, SuperGeometry2D<double>& superGeometry, int material,
		IndicatorF2D<double>& geometryIndicator, std::vector<int> bulkMaterials);

template void setBounceBackVelocityBoundary<double,descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, FunctorPtr<SuperIndicatorF2D<double>>&& boundaryIndicator,
		FunctorPtr<SuperIndicatorF2D<double>>&& bulkIndicator, IndicatorF2D<double>& geometryIndicator);

template void setBounceBackVelocityBoundary<double,descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block, BlockIndicatorF2D<double>& boundaryIndicator,
		BlockIndicatorF2D<double>& bulkIndicator, IndicatorF2D<double>& geometryIndicator);

template void setBounceBackVelocityBoundary<double,descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block, BlockGeometryStructure2D<double>& blockGeometryStructure,
		int iX, int iY, BlockIndicatorF2D<double>& bulkIndicator, IndicatorF2D<double>& geometryIndicator);

template void setBounceBackVelocityBoundary<double,descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block, BlockGeometryStructure2D<double>& blockGeometryStructure,
		int x, int y, double distances[descriptors::D2Q9<>::q]);

template void setBounceBackVelocityBoundary<double,descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block, BlockGeometryStructure2D<double>& blockGeometryStructure,
		int x, int y, int iPop, double dist);



template void setBouzidiZeroVelocityBoundary<double,descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, SuperGeometry2D<double>& superGeometry, int material,
		IndicatorF2D<double>& geometryIndicator, std::vector<int> bulkMaterials);

template void setBouzidiZeroVelocityBoundary<double,descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, FunctorPtr<SuperIndicatorF2D<double>>&& boundaryIndicator,
		FunctorPtr<SuperIndicatorF2D<double>>&& bulkIndicator, IndicatorF2D<double>& geometryIndicator);

template void setBouzidiZeroVelocityBoundary<double,descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block, BlockIndicatorF2D<double>& boundaryIndicator,
		BlockIndicatorF2D<double>& bulkIndicator, IndicatorF2D<double>& geometryIndicator);

template void setBouzidiZeroVelocityBoundary<double,descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block, BlockGeometryStructure2D<double>& blockGeometryStructure,
		int iX, int iY, BlockIndicatorF2D<double>& bulkIndicator, IndicatorF2D<double>& geometryIndicator);

template void setBouzidiZeroVelocityBoundary<double,descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block, BlockGeometryStructure2D<double>& blockGeometryStructure,
		int x, int y, double distances[descriptors::D2Q9<>::q]);

template void setBouzidiZeroVelocityBoundary<double,descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block, BlockGeometryStructure2D<double>& blockGeometryStructure,
		int x, int y, int iPop, double dist);



template void setBounceBackZeroVelocityBoundary<double,descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, SuperGeometry2D<double>& superGeometry, int material,
		IndicatorF2D<double>& geometryIndicator, std::vector<int> bulkMaterials);

template void setBounceBackZeroVelocityBoundary<double,descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, FunctorPtr<SuperIndicatorF2D<double>>&& boundaryIndicator,
		FunctorPtr<SuperIndicatorF2D<double>>&& bulkIndicator, IndicatorF2D<double>& geometryIndicator);

template void setBounceBackZeroVelocityBoundary<double,descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block, BlockIndicatorF2D<double>& boundaryIndicator,
		BlockIndicatorF2D<double>& bulkIndicator, IndicatorF2D<double>& geometryIndicator);

template void setBounceBackZeroVelocityBoundary<double,descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block, BlockGeometryStructure2D<double>& blockGeometryStructure,
		int iX, int iY, BlockIndicatorF2D<double>& bulkIndicator, IndicatorF2D<double>& geometryIndicator);

template void setBounceBackZeroVelocityBoundary<double,descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block, BlockGeometryStructure2D<double>& blockGeometryStructure,
		int x, int y, double distances[descriptors::D2Q9<>::q]);

template void setBounceBackZeroVelocityBoundary<double,descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(BlockLatticeStructure2D<double,descriptors::D2Q9<>>& block, BlockGeometryStructure2D<double>& blockGeometryStructure,
		int x, int y, int iPop, double dist);



template void setZouHePressureBoundary<double, descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, double omega, SuperGeometry2D<double>& superGeometry, int material);

template void setZouHePressureBoundary<double, descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, double omega, FunctorPtr<SuperIndicatorF2D<double>>&& indicator);

template void setZouHePressureBoundary<double, descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(BlockLatticeStructure2D<double, descriptors::D2Q9<>>& block, double omega, BlockIndicatorF2D<double>& indicator, bool includeOuterCells);



template void setZouHeVelocityBoundary<double, descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, double omega, SuperGeometry2D<double>& superGeometry, int material);

template void setZouHeVelocityBoundary<double, descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(SuperLattice2D<double, descriptors::D2Q9<>>& sLattice, double omega, FunctorPtr<SuperIndicatorF2D<double>>&& indicator);

template void setZouHeVelocityBoundary<double, descriptors::D2Q9<>, BGKdynamics<double, descriptors::D2Q9<>>>
(BlockLatticeStructure2D<double, descriptors::D2Q9<>>& block, double omega, BlockIndicatorF2D<double>& indicator, bool includeOuterCells);




}  // namespace olb

