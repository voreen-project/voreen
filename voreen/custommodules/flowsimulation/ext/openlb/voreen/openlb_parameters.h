/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
 * Department of Computer Science.                                                 *
 * For a list of authors please refer to the file "CREDITS.txt".                   *
 *                                                                                 *
 * This file is part of the Voreen software package. Voreen is free software:      *
 * you can redistribute it and/or modify it under the terms of the GNU General     *
 * Public License version 2 as published by the Free Software Foundation.          *
 *                                                                                 *
 * Voreen is distributed in the hope that it will be useful, but WITHOUT ANY       *
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR   *
 * A PARTICULAR PURPOSE. See the GNU General Public License for more details.      *
 *                                                                                 *
 * You should have received a copy of the GNU General Public License in the file   *
 * "LICENSE.txt" along with this file. If not, see <http://www.gnu.org/licenses/>. *
 *                                                                                 *
 * For non-commercial academic use see the license exception specified in the file *
 * "LICENSE-academic.txt". To get information about commercial licensing please    *
 * contact the authors.                                                            *
 *                                                                                 *
 ***********************************************************************************/

#ifndef VRN_OPENLB_PARAMETERS_H
#define VRN_OPENLB_PARAMETERS_H

// Note: no OpenLB dependency shall be added here!

#include <map>
#include <string>
#include <vector>

#ifdef VRN_MODULE_FLOWSIMULATION
#include "voreen/core/voreencoreapi.h"
#include "tgt/vector.h"
using vec3 = tgt::vec3;
using vec4 = tgt::vec4;
#else
#define VRN_CORE_API
typedef struct { float x, y, z; } vec3;
typedef struct { float x, y, z, w; } vec4;
#endif

namespace voreen {

static const float VOREEN_LENGTH_TO_SI = 0.001f;
static const float VOREEN_TIME_TO_SI = 0.001f;
const std::string META_DATA_NAME_OFFSET = "Offset";
const std::string META_DATA_NAME_SPACING = "Spacing";
const std::string META_DATA_NAME_TIMESTEP = "Timestep";
const std::string META_DATA_NAME_REAL_WORLD_MAPPING = "RealWorldMapping";
const std::string META_DATA_NAME_MODALITY = "Modality";


// List of material ids that are basically predefined by OpenLB.
enum Material {
    MAT_EMPTY = 0,
    MAT_FLUID = 1,
    MAT_WALL = 2,
    MAT_COUNT,
};

// List of flow features which can be extracted during simulation.
// Values have to be power of two (bitfield).
enum FlowFeatures {
    FF_NONE = 0, ///< No flow feature
    FF_VELOCITY = 1, ///< Velocity vector field
    FF_MAGNITUDE = 2, ///< Magnitude scalar field (from velocity vector field)
    FF_PRESSURE = 4, ///< Pressure scalar field
    FF_WALLSHEARSTRESS = 8, ///< Wall shear stress scalar field
};

enum FlowIndicatorType {
    FIT_INVALID = -1, ///< Denotes an invalid indicator.
    FIT_CANDIDATE = 0, ///< This indicator is just a candidate and has no function yet.
    FIT_VELOCITY = 1, ///< This indicator is a velocity boundary condition.
    FIT_PRESSURE = 2, ///< This indicator is a pressure boundary condition.
    FIT_MEASURE = 3, ///< This indicator serves as a flux measure.
};

enum FlowProfile {
    FP_NONE = 0, ///< No flow profile
    FP_POISEUILLE = 1, ///< Poiseuille flow profile
    FP_POWERLAW = 2, ///< Power law flow profile
    FP_CONSTANT = 3, ///< constant flow profile
    FP_VOLUME = 4, ///< samples from volume
};

enum FlowBoundaryCondition {
    FBC_NONE = 0,
    FBC_BOUNCE_BACK = 1,
    FBC_BOUZIDI = 2,
};

enum FlowTurbulenceModel {
    FTM_NONE = 0,
    FTM_SMAGORINSKY = 1,
    FTM_SMAGORINSKY_SHEAR_IMPROVED = 2,
    FTM_SMAGORINSKY_CONSISTENT = 3,
    FTM_SMAGORINSKY_CONSISTENT_STRAIN = 4,
    FTM_BGK = 5,
};

class VRN_CORE_API VelocityCurve {
public:

    VelocityCurve();

    float operator()(float t) const;

    float& operator[](float t);

    void setPeriodic(bool enabled);
    bool isPeriodic() const;

    void setScale(float scale);
    float getScale() const;

    float getMinVelocity() const;
    float getMaxVelocity() const;

    float getStartTime() const;
    float getEndTime() const;
    float getDuration() const;

    static VelocityCurve createConstantCurve(float value);
    static VelocityCurve createLinearCurve(float duration, float maxValue);
    static VelocityCurve createSinusoidalCurve(float duration, float maxValue, int steps = 30);
    static VelocityCurve createHumanHeartBeat();
    static VelocityCurve createFromCSV(const std::string& file);
    static VelocityCurve createFromMap(const std::map<float, float>& map);

protected:
    std::map<float, float> peakVelocities_;
    bool periodic_;
    float scale_;
};

// Indicates flux through an arbitrary, circle-shaped area.
struct FlowIndicator {

    FlowIndicatorType type_{FIT_INVALID};       ///< FlowIndicator type, @see FlowIndicatorType.
    int id_{-1};                                ///< Unique identifier. Also used by OpenLB to indicate material.
    std::string name_{};                        ///< Optional name.
    vec4 color_{1, 1, 1, 1};                    ///< Optional color.

    vec3 center_{0, 0, 0};                      ///< Center position of the circle shaped area in world space.
    vec3 normal_{0, 0, 0};                      ///< (Normalized) Normal vector defining the orientation.
    float radius_{0};                           ///< Radius of the disk.
    float length_{0};                           ///< Length of the disk/cylinder.
    bool roleSwapped_{false};                   ///< Determines if this inlet serves as outlet and vice-versa.

    // Used by generating flow indicators:
    FlowProfile flowProfile_{FP_NONE};          ///< Flow profile, @see FlowProfile.
    VelocityCurve velocityCurve_{};             ///< Velocity curve mapping time points to velocities.
};

/**
 * Datastructure used to represent flow parameters for setting up a flow simulation. It is used in the flowsimulation module.
 */
class VRN_CORE_API Parameters {
public:

    /** Constructor */
    Parameters();

    virtual ~Parameters() = default;

    /**
     * This function generates a unique and distinguishable name for each parametrization.
     */
    std::string name_;

    /**
     * Returns the spatial resolution in voxels, the largest vessel diameter should be divided into.
     * A high resolution is important for simulation accuracy.
     */
    int spatialResolution_;

    /**
     * Returns the relaxation time parameter.
     */
    float relaxationTime_;

    /**
     * Returns the max expected length in m within the simulation geometry.
     * E.g., the largest diameter of all contained vessels.
     */
    float characteristicLength_;

    /**
     * Returns the highest expected velocity in m/s.
     */
    float characteristicVelocity_;

    /**
     * Returns the kinematic viscosity in m^2/s.
     */
    float viscosity_;

    /**
     * Returns the fluid mass density in kg/m^3
     */
    float density_;

    /**
     * Returns the turbulence model.
     */
    FlowTurbulenceModel turbulenceModel_;

    /**
     * Returns the constant used by the turbulence model.
     */
    float smagorinskyConstant_;

    /**
     * Returns the wall boundary condition.
     */
    FlowBoundaryCondition wallBoundaryCondition_;

    /**
     * Returns if the lattice shall be perturbed.
     * This is able to enforce the turbulent regime in an otherwise laminar simulation.
     */
    bool latticePerturbation_;

    /**
     * Returns a factor that is multiplied with the inlet velocity.
     */
    float inletVelocityMultiplier_;

    /**
     * Returns Reynolds number.
     */
    float getReynoldsNumber() const;

    /**
     * Determines, if this parameter set operates in the stable limit.
     */
    bool isValid() const;
};

}

#endif
