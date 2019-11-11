/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#ifndef VRN_FLOWPARAMETERS_H
#define VRN_FLOWPARAMETERS_H

#include "voreen/core/voreencoreapi.h"

#include "voreen/core/datastructures/datainvalidationobserver.h"
#include "voreen/core/io/serialization/serializable.h"

#include "tgt/vector.h"

namespace voreen {

// List of flow features which can be extracted during simulation.
// Values have to be power of two (bitfield).
enum FlowFeatures {
    FF_NONE             = 0,
    FF_VELOCITY         = 1,
    FF_MAGNITUDE        = 2,
    FF_PRESSURE         = 4,
    FF_WALLSHEARSTRESS  = 8,
};

enum FlowIndicatorType {
    FIT_CANDIDATE = 0, ///< This indicator is just a candidate and has no function yet.
    FIT_GENERATOR = 1, ///< This indicator is flow generating.
    FIT_PRESSURE  = 2, ///< This indicator is a pressure boundary.
    FIT_MEASURE   = 3, ///< This indicator serves as a flux measure.
};

enum FlowProfile {
    FP_NONE       = 0,
    FP_POISEUILLE = 1,
    FP_POWERLAW   = 2,
    FP_CONSTANT   = 3,
};

enum FlowStartPhase {
    FSP_NONE     = 0,
    FSP_CONSTANT = 1,
    FSP_SINUS    = 2,
};

// Indicates flux through an arbitrary, circle-shaped area.
struct VRN_CORE_API FlowIndicator : public Serializable {

    FlowIndicatorType   type_;
    tgt::vec3           center_;
    tgt::vec3           normal_;
    float               radius_;

    // Used by generating flow indicators:
    FlowProfile         flowProfile_;

    FlowStartPhase      startPhaseFunction_;
    float               startPhaseDuration_;

    float               targetVelocity_; // TODO: remove.

    FlowIndicator();

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);
};

/**
 * Datastructure used to represent flow parameters for setting up a flow simulation. It is used in the flowsimulation module.
 */
class VRN_CORE_API FlowParameterSet : public Serializable {
public:

    /** Constructor */
    FlowParameterSet(); // For deserialization only.
    explicit FlowParameterSet(const std::string& name);

    /**
     * This function generates a unique and distinguishable name for each parametrization.
     */
    const std::string& getName() const;

    /**
     * Returns the spatial resolution in voxels, the largest vessel diameter should be divided into.
     * A high resolution is important for simulation accuracy.
     */
    int getSpatialResolution() const;
    void setSpatialResolution(int spatialResolution);

    /**
     * Returns the temporal granularity, i.e., the duration of each time step in s.
     */
    float getTemporalResolution() const;
    void setTemporalResolution(float temporalResolution);

    /**
     * Returns the max expected length in mm within the simulation geometry.
     * E.g., the largest diameter of all contained vessels.
     */
    float getCharacteristicLength() const;
    void setCharacteristicLength(float characteristicLength);

    /**
     * Returns the highest expected velocity in mm/s.
     */
    float getCharacteristicVelocity() const;
    void setCharacteristicVelocity(float characteristicVelocity);

    /**
     * Returns the kinematic viscosity in 10^-3 m^2/s.
     * Note: in order to achieve the correct physical value, multiply by 0.001
     */
    float getViscosity() const;
    void setViscosity(float viscosity);

    /**
     * Returns the mass density in kg/m^3
     */
    float getDensity() const;
    void setDensity(float density);

    /**
     * Returns the constant for the Smagorinsky turbulence model.
     */
    float getSmagorinskyConstant() const;
    void setSmagorinskyConstant(float smagorinskyConstant);

    /**
     * Returns whether Bouzidi boundary condition should be used for unaligned simulation geometries.
     */
    bool getBouzidi() const;
    void setBouzidi(bool bouzidi);

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

private:

    // Identifier of the parameter set.
    std::string name_;

    // All other relevant parameters.
    int spatialResolution_;        ///< spatial resolution in voxels (per dimension and characteristic length)
    float temporalResolution_;     ///< temporal resolution in seconds
    float characteristicLength_;    ///< characteristic length in mm
    float characteristicVelocity_;  ///< characteristic velocity in mm/s
    float viscosity_;               ///< viscosity in 10^-3 m^2/s
    float density_;                 ///< density in kg/m^3
    float smagorinskyConstant_;     ///< constant for Smagorinsky turbulence model
    bool bouzidi_;                  ///< bouzidi boundary condition
};

/**
 * Parametrization List, implementing thread safety for AsyncComputeProcessor.
 */
class VRN_CORE_API FlowParameterSetEnsemble : public DataInvalidationObservable, public Serializable {

    static const int VERSION;

public:

    static const size_t ALL_PARAMETER_SETS;

    explicit FlowParameterSetEnsemble(const std::string& name);
    FlowParameterSetEnsemble(const FlowParameterSetEnsemble& origin);

    const std::string& getName() const;

    /**
     * Returns the time in s which should be simulated.
     */
    float getSimulationTime() const;
    void setSimulationTime(float simulationTime);

    /**
     * Returns the number of time steps (intermediate results), the simulation should store.
     * Note: The initial time step will be stored additionally and in any case.
     */
    int getNumTimeSteps() const;
    void setNumTimeSteps(int numTimeSteps);

    /**
     * Returns the output resolution of the intermediate time steps for each volume and their dimension.
     * This enforces basically a resampling of the simulation domain.
     * Note: this currently acts as the max(!) resolution. If all the features can be captured by a lower resolution
     *       the latter will be taken. TODO: rename, if tested properly!
     */
    int getOutputResolution() const;
    void setOutputResolution(int outputResolution);

    /**
     * Returns the flow features as bitmask.
     * To test for a single feature, binarily test for the available features.
     */
    int getFlowFeatures() const;
    void setFlowFeatures(int flowFeatures);

    void addFlowIndicator(const FlowIndicator& flowIndicator);
    const std::vector<FlowIndicator>& getFlowIndicators() const;

    void addFlowParameterSet(const FlowParameterSet& parameters);
    const std::vector<FlowParameterSet>& getFlowParameterSets() const;

    // Shortcuts
    bool empty() const;
    size_t size() const;
    const FlowParameterSet& at(size_t index) const;

    /** Used to export parametrization file. */
    std::string toJSONString(size_t param = ALL_PARAMETER_SETS) const;
    std::string toXMLString(size_t param = ALL_PARAMETER_SETS) const;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

private:

    void serializeInternal(Serializer& s, size_t param) const;

    // Ensemble name.
    std::string name_;

    // Configuration.
    float simulationTime_;         ///< simulation time in seconds
    int numTimeSteps_;             ///< number of time steps of output
    int outputResolution_;         ///< spatial resolution of output in voxels (per dimension)
    int flowFeatures_;             ///< bitmask storing flow features

    // Flow indication (in-/out flow).
    std::vector<FlowIndicator> flowIndicators_;

    // Actual parameters.
    std::vector<FlowParameterSet> flowParameterSets;
};

}   // namespace

#endif
