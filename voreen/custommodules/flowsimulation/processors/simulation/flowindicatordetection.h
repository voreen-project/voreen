/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#ifndef VRN_FLOWINDICATORDETECTION_H
#define VRN_FLOWINDICATORDETECTION_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/properties/string/stringtableproperty.h"

#include "../../datastructures/flowparameters.h"
#include "../../ports/flowparametrizationport.h"

#include "modules/vesselnetworkanalysis/ports/vesselgraphport.h"

namespace voreen {

/**
 * This processor is being used to select in and out flow.
 */
class VRN_CORE_API FlowIndicatorDetection : public Processor {
public:
    FlowIndicatorDetection();
    virtual Processor* create() const         { return new FlowIndicatorDetection();    }

    virtual std::string getClassName() const  { return "FlowIndicatorDetection";        }
    virtual std::string getCategory() const   { return "Simulation";                    }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL;         }

    virtual bool isReady() const;
    virtual void process();

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

protected:

    virtual void adjustPropertiesToInput();

    virtual void setDescriptions() {
        setDescription("This processor is being used to select in and out flow.");
    }

private:

    /**
     * This helper struct stores the settings of each indicator
     * that they don't store directly.
     */
    struct FlowIndicatorSettings : public Serializable {
        // With respect to current vessel graph.
        VGNodeID nodeId_;
        VGEdgeID edgeId_;
        int centerlinePosition_;
        float relativeRadiusCorrection_;
        bool invertDirection_;
        bool forceAxisAlignment_;

        // Other settings.
        std::string velocityCurveType_;
        float velocityCurveDuration_;
        float targetVelocity_;
        std::string velocityCurveFile_;
        bool velocityCurvePeriodic_;

        FlowIndicatorSettings();
        FlowIndicatorSettings(VGNodeID nodeId, VGEdgeID edgeId);

        virtual void serialize(Serializer& s) const;
        virtual void deserialize(Deserializer& s);
    };

    /**
     * Callback that is triggered as soon as the vessel graph changes.
     * This will add new flow indicator candidates (and estimate their type).
     */
    void detectFlowIndicators(bool forced);

    /**
     * Callback that is triggered as soon as the selection was changed.
     * This will set the related property values accordingly.
     */
    void updateIndicatorUI();

    /**
     * Callback that is triggered as soon as properties of the selected indicator
     * were changed. This will update the indicator accordingly.
     */
    void onIndicatorConfigChange(bool needReinitialization);
    inline void onIndicatorConfigChange() { onIndicatorConfigChange(false); }

    /**
     * Callback for cloning the currently selected flow indicator (if any).
     */
    void onCloneFlowIndicator();

    /**
     * Callback for removing the currently selected flow indicator (if any).
     */
    void onRemoveFlowIndicator();

    /**
     * Estimates the indicator type.
     */
    FlowIndicatorType estimateType(const FlowIndicator& indicator, const tgt::vec3& velocity) const;

    /**
     * Initializes a flow indicator according to the specified settings.
     * This includes position and normal and radius.
     */
    FlowIndicator initializeIndicator(FlowIndicatorSettings& settings);

    /**
     * Creates a velocity curve according to the specified settings.
     */
    VelocityCurve createCurveFromSettings(FlowIndicatorSettings& settings);

    /**
     * Creates the overview table from the current flow indicator list.
     * Note: This might remove the widget's focus!
     */
    void buildTable();

    FlowParametrizationPort parameterInport_;
    VesselGraphPort vesselGraphPort_;
    VolumePort volumePort_;
    FlowParametrizationPort parameterOutport_;

    StringTableProperty flowIndicatorTable_;
    ButtonProperty cloneFlowIndicator_;
    ButtonProperty removeFlowIndicator_;
    ButtonProperty resetFlowIndicators_;
    IntProperty angleThreshold_;

    StringProperty indicatorName_;
    IntProperty centerlinePosition_;
    FloatProperty radius_;
    FloatProperty relativeRadiusCorrection_;
    BoolProperty invertDirection_;
    BoolProperty forceAxisAlignment_;
    OptionProperty<FlowIndicatorType> indicatorType_;

    OptionProperty<FlowProfile> flowProfile_;
    StringOptionProperty velocityCurveType_;
    FloatProperty velocityCurveDuration_;
    FloatProperty targetVelocity_;
    FileDialogProperty velocityCurveFile_;
    BoolProperty velocityCurvePeriodicity_;

    std::vector<FlowIndicator> flowIndicators_;
    std::vector<FlowIndicatorSettings> flowIndicatorSettings_;

    bool triggertBySelection_;

    std::string vesselGraphHash_;

    static const std::string loggerCat_;
};

}   //namespace

#endif
