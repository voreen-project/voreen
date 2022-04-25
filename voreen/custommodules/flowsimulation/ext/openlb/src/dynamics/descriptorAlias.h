/*  This file is part of the OpenLB library
 *
 *  Copyright (C) 2019 Mathias J. Krause, Tom Braun, Jonas Latt,
 *  Adrian Kummerländer, Patrick Nathen, Andrea Parmigiani,
 *  Orestis Malaspinas, Albert Mink, Marc Haussmann, Davide Dapelo
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

#ifndef DESCRIPTOR_ALIAS_H
#define DESCRIPTOR_ALIAS_H


#include "latticeDescriptors.h"
#include "mrtLatticeDescriptors.h"
#include "rtlbmDescriptors.h"

#include <vector>

namespace olb {

namespace descriptors {


/// D2Q5 descriptor
using D2Q5Descriptor                                = D2Q5<>;
using AdvectionDiffusionD2Q5Descriptor              = D2Q5<VELOCITY>;
using SourcedAdvectionDiffusionD2Q5Descriptor       = D2Q5<SOURCE,VELOCITY>;
using SmagorinskyAdvectionDiffusionD2Q5Descriptor   = D2Q5<VELOCITY,TAU_EFF>;
using MixedScaleAdvectionDiffusionD2Q5Descriptor    = D2Q5<VELOCITY,TAU_EFF,CUTOFF_HEAT_FLUX>;


/// MRTD2Q5 descriptor
using MRTD2Q5Descriptor                             = D2Q5<tag::MRT>;
using AdvectionDiffusionMRTD2Q5Descriptor           = D2Q5<tag::MRT,VELOCITY>;


/// D2Q9 descriptor
using D2Q9Descriptor                                = D2Q9<>;
using AdvectionDiffusionD2Q9Descriptor              = D2Q9<VELOCITY>;
using ForcedD2Q9Descriptor                          = D2Q9<FORCE>;
using SmagorinskyForcedD2Q9Descriptor               = D2Q9<FORCE,TAU_EFF>;
using MixedScaleForcedD2Q9Descriptor                = D2Q9<FORCE,TAU_EFF,CUTOFF_KIN_ENERGY>;
using FreeEnergyD2Q9Descriptor                      = D2Q9<CHEM_POTENTIAL,FORCE>;
using V6ForcedD2Q9Descriptor                        = D2Q9<FORCE,V6>;
using DynOmegaD2Q9Descriptor                        = D2Q9<OMEGA>;
using ForcedDynOmegaD2Q9Descriptor                  = D2Q9<FORCE,OMEGA>;
using DynOmegaPorousParticleD2Q9Descriptor          = D2Q9<POROSITY,VELOCITY_NUMERATOR,VELOCITY_DENOMINATOR,OMEGA>;
using ShearSmagorinskyD2Q9Descriptor                = D2Q9<AV_SHEAR>;
using ShanChenForcedD2Q9Descriptor                  = D2Q9<VELOCITY,FORCE>;
using ShanChenDynOmegaD2Q9Descriptor                = D2Q9<VELOCITY,FORCE,OMEGA>;
using ShanChenDynOmegaForcedD2Q9Descriptor          = D2Q9<VELOCITY,FORCE,EXTERNAL_FORCE,OMEGA>;
using ShanChenDynGD2Q9Descriptor                    = D2Q9<VELOCITY,FORCE,G>;
using ShanChenDynGForcedD2Q9Descriptor              = D2Q9<VELOCITY,FORCE,EXTERNAL_FORCE,G>;
using DynSmagorinskyD2Q9Descriptor                  = D2Q9<SMAGO_CONST>;
using PorousD2Q9Descriptor                          = D2Q9<POROSITY>;
using ExtendedPorousD2Q9Descriptor                  = D2Q9<POROSITY,LOCAL_DRAG>;
using PorousParticleD2Q9Descriptor                  = D2Q9<POROSITY,VELOCITY_NUMERATOR,VELOCITY_DENOMINATOR>;
using PorousParticleWithCollisionD2Q9Descriptor     = D2Q9<POROSITY,VELOCITY_NUMERATOR,VELOCITY_DENOMINATOR,COLLISION_DETECTION>;
using PSMD2Q9Descriptor                             = D2Q9<POROSITY,VELOCITY_SOLID>;
using GuoZhaoD2Q9Descriptor                         = D2Q9<FORCE,EPSILON,K,NU,BODY_FORCE>;
using ADMD2Q9Descriptor                             = D2Q9<FIL_RHO,LOCAL_FIL_VEL_X,LOCAL_FIL_VEL_Y>;
using ForcedADMD2Q9Descriptor                       = D2Q9<FORCE,FIL_RHO,LOCAL_FIL_VEL_X,LOCAL_FIL_VEL_Y>;


/// MRTD2Q9 descriptor
using MRTD2Q9Descriptor                             = D2Q9<tag::MRT>;
using ForcedMRTD2Q9Descriptor                       = D2Q9<tag::MRT,FORCE>;


/// D3Q7 descriptor
using D3Q7Descriptor                                = D3Q7<>;
using AdvectionDiffusionD3Q7Descriptor              = D3Q7<VELOCITY>;
using SourcedAdvectionDiffusionD3Q7Descriptor       = D3Q7<SOURCE,VELOCITY>;
using SmagorinskyAdvectionDiffusionD3Q7Descriptor   = D3Q7<VELOCITY,TAU_EFF>;
using ParticleAdvectionDiffusionD3Q7Descriptor      = D3Q7<VELOCITY,VELOCITY2>;
using ParticleAdvectionDiffusionMRTD3Q7Descriptor   = D3Q7<VELOCITY,VELOCITY2>;
using PorousAdvectionDiffusionD3Q7Descriptor        = D3Q7<POROSITY,VELOCITY>;
using D3Q7DescriptorLebedev                         = D3Q7<tag::RTLBM>;


/// MRTD3Q7 descriptor
using MRTD3Q7Descriptor                             = D3Q7<tag::MRT>;
using AdvectionDiffusionMRTD3Q7Descriptor           = D3Q7<tag::MRT,VELOCITY>;


/// D3Q13 descriptor
using ForcedD3Q13Descriptor                         = D3Q13<FORCE>;


/// D3Q15 descriptor
using D3Q15Descriptor                               = D3Q15<>;
using ForcedD3Q15Descriptor                         = D3Q15<FORCE>;
using D3Q15DescriptorLebedev                        = D3Q15<tag::RTLBM>;


/// D3Q19 descriptor
using D3Q19Descriptor                               = D3Q19<>;
using ForcedD3Q19Descriptor                         = D3Q19<FORCE>;
using SmagorinskyForcedD3Q19Descriptor              = D3Q19<FORCE,TAU_EFF>;
using V12ForcedD3Q19Descriptor                      = D3Q19<FORCE,V12>;
using FreeEnergyD3Q19Descriptor                     = D3Q19<CHEM_POTENTIAL,FORCE>;
using ParticleAdvectionDiffusionD3Q19Descriptor     = D3Q19<VELOCITY,VELOCITY2>;
using DynOmegaD3Q19Descriptor                       = D3Q19<OMEGA>;
using ForcedDynOmegaD3Q19Descriptor                 = D3Q19<OMEGA,FORCE>;
using DynOmegaPorousParticleD3Q19Descriptor         = D3Q19<POROSITY,VELOCITY_NUMERATOR,VELOCITY_DENOMINATOR,OMEGA>;
using ShearSmagorinskyD3Q19Descriptor               = D3Q19<AV_SHEAR>;
using ShearSmagorinskyForcedD3Q19Descriptor         = D3Q19<AV_SHEAR,FORCE>;
using ForcedShearWallSmagorinskyD3Q19Descriptor     = D3Q19<AV_SHEAR,FORCE,TAU_W>;
using FDKalmanShearSmagorinskyD3Q19Descriptor       = D3Q19<ERROR_COVARIANCE,VARIANCE,VELOCITY,FILTERED_VEL_GRAD,VELO_GRAD>;
using FDKalmanShearSmagorinskyForcedD3Q19Descriptor = D3Q19<ERROR_COVARIANCE,VARIANCE,VELOCITY,FILTERED_VEL_GRAD,VELO_GRAD,FORCE>;

// Kalman filter : Adaptive exponential smoothing //
// 3D Descriptors for flow with Shear-Improved Smagorinsky - Kalman Filter
// Boudet et al. (2016) A Kalman filter adapted of the estimation of mean gradients
//   in the a large-eddy simulation of unsteady turbulent flows.
using KalmanShearSmagorinskyD3Q19Descriptor         = D3Q19<ERROR_COVARIANCE,VARIANCE,TAU_SGS,FILTERED_POPULATION>;

using WALED3Q19Descriptor                           = D3Q19<EFFECTIVE_OMEGA,VELO_GRAD>;
using WALEForcedD3Q19Descriptor                     = D3Q19<EFFECTIVE_OMEGA,VELO_GRAD,FORCE>;
using ShanChenForcedD3Q19Descriptor                 = D3Q19<VELOCITY,FORCE>;
using ShanChenDynOmegaD3Q19Descriptor               = D3Q19<VELOCITY,FORCE,OMEGA>;
using ShanChenDynOmegaForcedD3Q19Descriptor         = D3Q19<VELOCITY,FORCE,EXTERNAL_FORCE,OMEGA>;
using WallFunctionD3Q19Descriptor                   = D3Q19<TAU_W,TAU_EFF>;
using WallFunctionForcedD3Q19Descriptor             = D3Q19<TAU_W,TAU_EFF,FORCE>;
using DynSmagorinskyD3Q19Descriptor                 = D3Q19<SMAGO_CONST>;
using PorousD3Q19Descriptor                         = D3Q19<POROSITY>;
using PorousForcedD3Q19Descriptor                   = D3Q19<POROSITY,FORCE>;
using ExtendedPorousD3Q19Descriptor                 = D3Q19<POROSITY,LOCAL_DRAG>;
using PorousParticleD3Q19Descriptor                 = D3Q19<POROSITY,VELOCITY_NUMERATOR,VELOCITY_DENOMINATOR>;
using PorousParticleWithCollisionD3Q19Descriptor    = D3Q19<POROSITY,VELOCITY_NUMERATOR,VELOCITY_DENOMINATOR,COLLISION_DETECTION>;
using PSMD3Q19Descriptor                            = D3Q19<POROSITY,VELOCITY_SOLID>;
using ForcedPSMD3Q19Descriptor                      = D3Q19<POROSITY,VELOCITY_SOLID,FORCE>;
using GuoZhaoD3Q19Descriptor                        = D3Q19<FORCE,EPSILON,K,NU,BODY_FORCE>;
using ADMD3Q19Descriptor                            = D3Q19<FIL_RHO,LOCAL_FIL_VEL_X,LOCAL_FIL_VEL_Y,LOCAL_FIL_VEL_Z>;
using ForcedADMD3Q19Descriptor                      = D3Q19<FORCE,FIL_RHO,LOCAL_FIL_VEL_X,LOCAL_FIL_VEL_Y,LOCAL_FIL_VEL_Z>;
using ForcedAdaptiveADMD3Q19Descriptor              = D3Q19<FORCE,FIL_RHO,LOCAL_FIL_VEL_X,LOCAL_FIL_VEL_Y,LOCAL_FIL_VEL_Z,LOCAL_AV_DISS,LOCAL_AV_TKE,LOCAL_SIGMA_ADM,LOCAL_NU_EDDY,TAU_W>;


/// MRTD3Q19 descriptor
using MRTD3Q19Descriptor                            = D3Q19<tag::MRT>;
using ForcedMRTD3Q19Descriptor                      = D3Q19<tag::MRT,FORCE>;


/// D3Q27 descriptor
using D3Q27Descriptor                               = D3Q27<>;
using ForcedD3Q27Descriptor                         = D3Q27<FORCE>;
using DynOmegaD3Q27Descriptor                       = D3Q27<OMEGA>;
using ForcedDynOmegaD3Q27Descriptor                 = D3Q27<OMEGA,FORCE>;
using DynOmegaPorousParticleD3Q27Descriptor         = D3Q27<POROSITY,VELOCITY_NUMERATOR,VELOCITY_DENOMINATOR,OMEGA>;
using WALED3Q27Descriptor                           = D3Q27<EFFECTIVE_OMEGA,VELO_GRAD>;
using WALEForcedD3Q27Descriptor                     = D3Q27<EFFECTIVE_OMEGA,VELO_GRAD,FORCE>;
using WallFunctionForcedD3Q27Descriptor             = D3Q27<TAU_W,TAU_EFF,FORCE>;
using D3Q27DescriptorLebedev                        = D3Q27<tag::RTLBM>;

}

}

#endif