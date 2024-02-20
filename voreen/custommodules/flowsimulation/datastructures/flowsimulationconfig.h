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

#ifndef VRN_FLOWPARAMETERS_H
#define VRN_FLOWPARAMETERS_H

#include "../ext/openlb/voreen/openlb_parameters.h"

#include "voreen/core/voreencoreapi.h"

#include "voreen/core/datastructures/datainvalidationobserver.h"
#include "voreen/core/io/serialization/serializable.h"

#include "tgt/vector.h"
#include "tgt/matrix.h"

#include <map>

namespace voreen {

/**
 * Parametrization List, implementing thread safety for AsyncComputeProcessor.
 */
class VRN_CORE_API FlowSimulationConfig : public DataInvalidationObservable, public Serializable {

    static const int VERSION;

public:

    static const size_t ALL_PARAMETER_SETS;

    /**
     * Returns the offset used to generate flow indicator ids.
     * Note: the offset depends on the simulation framework.
     * E.g. OpenLB requires an offset of 3.
     */
    static int getFlowIndicatorIdOffset() { return MAT_COUNT; }


    explicit FlowSimulationConfig(const std::string& name);
    FlowSimulationConfig(const FlowSimulationConfig& origin);

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
     * Note: this currently acts as the max(!) resolution. If all the features can be captured at a lower resolution
     *       the latter will be taken.
     */
    int getOutputResolution() const;
    void setOutputResolution(int outputResolution);

    /**
     * Returns the output file format.
     */
    const std::string& getOutputFileFormat() const;
    void setOutputFileFormat(const std::string& format);

    /**
     * Returns the flow features as bitmask.
     * To test for a single feature, perform a bit test for the available features.
     */
    int getFlowFeatures() const;
    void setFlowFeatures(int flowFeatures);

    /**
     * Returns the Transformation matrix for the entire domain.
     * This includes flow indicators.
     */
    tgt::mat4 getTransformationMatrix() const;
    void setTransformationMatrix(tgt::mat4 transformation);
    tgt::mat4 getInvertedTransformationMatrix() const;

    /**
     * Add a flow indicator to the internal list.
     * Note: This will set the unique id within the parameter set ensemble.
     */
    void addFlowIndicator(const FlowIndicator& flowIndicator);
    std::vector<FlowIndicator> getFlowIndicators(bool transformed = false) const;

    void addFlowParameterSet(const Parameters& parameters);
    const std::vector<Parameters>& getFlowParameterSets() const;

    // Shortcuts
    bool empty() const;
    size_t size() const;
    const Parameters& at(size_t index) const;

    /** Used to export parametrization file. */
    std::string toJSONString(size_t param = ALL_PARAMETER_SETS) const;
    std::string toXMLString(size_t param = ALL_PARAMETER_SETS) const;

    void serialize(Serializer& s) const;
    void deserialize(Deserializer& s);

private:

    void serializeInternal(Serializer& s, size_t param) const;
    int generateIndicatorId() const;

    // Ensemble name.
    std::string name_;

    // Configuration.
    float simulationTime_;         ///< simulation time in seconds
    int numTimeSteps_;             ///< number of time steps of output
    int outputResolution_;         ///< spatial resolution of output in voxels (per dimension)
    std::string outputFileFormat_; ///< output file format
    int flowFeatures_;             ///< bitmask storing flow features
    tgt::mat4 transformation_;     ///< transformation matrix for the domain (geometry, indicators, ...)

    // Flow indication (in-/out flow).
    std::vector<FlowIndicator> flowIndicators_;

    // Actual parameters.
    std::vector<Parameters> flowParameters_;
};

}   // namespace

#endif
