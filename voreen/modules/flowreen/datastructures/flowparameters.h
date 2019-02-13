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

enum FlowDirection {
    NONE = -1,
    IN   =  0,
    OUT  =  1,
};

// Indicates flux through an arbitrary, circle-shaped area.
struct VRN_CORE_API FlowIndicator {
    FlowDirection   direction_      { NONE };
    tgt::vec3       center_         { tgt::vec3::zero };
    tgt::vec3       normal_         { tgt::vec3::zero };
    float           radius_         { 0.0f };
};

/**
 * Datastructure used to represent flow parameters for setting up a flow simulation. It is used in the flowreen module.
 */
class VRN_CORE_API FlowParameters : public Serializable {
public:

    /** Constructor */
    explicit FlowParameters(const std::string& name);

    /**
     * This function generates a unique and distinguishable name for each parametrization.
     */
    const std::string& getName() const;

    float getCharacteristicLength() const;
    void setCharacteristicLength(float characteristicLength);

    float getViscosity() const;
    void setViscosity(float viscosity);

    float getDensity() const;
    void setDensity(float density);

    bool getBouzidi() const;
    void setBouzidi(bool bouzidi);

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

private:

    // Identifier of the parameterization.
    std::string name_;

    // All other relevant parameters.
    float characteristicLength_;
    float viscosity_;
    float density_;
    bool bouzidi_;
};

/**
 * Parametrization List, implementing thread safety for AsyncComputeProcessor.
 */
class VRN_CORE_API FlowParametrizationList : public DataInvalidationObservable, public Serializable {

    static const int VERSION;

public:

    explicit FlowParametrizationList(const std::string& name);
    FlowParametrizationList(const FlowParametrizationList& origin);

    const std::string& getName() const;

    float getSimulationTime() const;
    void setSimulationTime(float simulationTime);

    float getTemporalResolution() const;
    void setTemporalResolution(float temporalResolution);

    void addFlowIndicator(const FlowIndicator& flowIndicator);
    const std::vector<FlowIndicator>& getFlowIndicators() const;

    void addFlowParameters(const FlowParameters& parameters);
    const std::vector<FlowParameters>& getFlowParametrizations() const;

    // Shortcuts
    bool empty() const;
    size_t size() const;
    const FlowParameters& at(size_t index) const;

    /** Used to save as CSV file. */
    std::string toCSVString() const;
    std::string toJSONString() const;
    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

private:

    // Ensemble name.
    std::string name_;

    // Time specification.
    float simulationTime_;
    float temporalResolution_;

    // Flow indication (in-/out flow).
    std::vector<FlowIndicator> flowIndicators_;

    // Actual parameters.
    std::vector<FlowParameters> flowParametrizations_;
};

}   // namespace

#endif
