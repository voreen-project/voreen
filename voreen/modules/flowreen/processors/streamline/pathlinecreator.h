/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2017 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#ifndef VRN_PATHLINECREATOR_H
#define VRN_PATHLINECREATOR_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/utils/backgroundthread.h"

#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/progressproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"

#include "voreen/core/ports/volumeport.h"
#include "../../ports/streamlinelistport.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"

#include <random>

namespace voreen {

    class PathlineCreatorBackgroundThread;
    class StreamlineBundleDetectorBackgroundThread;
    class StreamlineList;

/**
 * This processor is used to create streamlines from a vec3 volume.
 * It can be used with the PathlineRenderer3D. At the moment only RAM volumes are supported.
 *
 * @Note: It uses a background thread to handle changed parameters during calculation.
 */
class VRN_CORE_API PathlineCreator : public Processor, public PortObserver {
    friend class PathlineCreatorBackgroundThread;
public:
    PathlineCreator();
    virtual ~PathlineCreator();

    virtual Processor* create() const { return new PathlineCreator(); }

    virtual std::string getCategory() const { return "Pathline Processing"; }
    virtual std::string getClassName() const { return "PathlineCreator"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_EXPERIMENTAL; }
protected:
    virtual void setDescriptions() {
    }

    virtual bool usesExpensiveComputation() const { return true; }

    /** Triggers the backgound thread. */
    virtual void process();
    /** Invalidates the processor, if wait for thread is active. */
    virtual void afterProcess();

    //------------------
    //  Observer
    //------------------
    virtual void afterConnectionAdded(const Port* source, const Port* connectedPort);
    virtual void beforeConnectionRemoved(const Port* source, const Port*);

    //------------------
    //  Thread handling
    //------------------
    void stopBackgroundThread();

    //------------------
    //  Callbacks
    //------------------
    /** Adjusts the threshold property values on inport changes. */
    //void volumePortHasChanged();
    /** Adjusts the relative threshold according to the absolute one. */
    void adjustRelativeThreshold();
    /** Generates new default seeds on button press. */
    void resetSeedsOnChange();
    /** Starts a new calculation on button press. */
    void calculatePathlinesOnChange();
    /** PathlineSettings have been changed. */
    void streamlineSettingsHaveBeenChanged();
    /** PathlineBundleSettings have been changed. */
    void streamlineBundleSettingsHaveBeenChanged();
    /** The toggle for bundle detection has actually been toggled. */
    void generateBundlesHasBeenChanged();

    //------------------
    //  Members
    //------------------
private:
    enum FilterMode {
        NEAREST,
        LINEAR
    };

    enum DirtyFlag {
        NONE,
        PATHLINES,
        PATHLINEBUNDLES
    };

    //general configuration
    ButtonProperty calculatePathlinesProp_;
    BoolProperty autoGenerateProp_;
    BoolProperty waitForThreadFinishedProp_;
    ProgressProperty progressProp_;                         ///< used for feedback in application mode

    //streamline settings
    IntProperty maxNumPathlinesProp_;                     ///< maximal number of streamlines
    IntIntervalProperty streamlineLengthThresholdProp_;         ///< streamline length must be in this interval
    FloatIntervalProperty absoluteMagnitudeThresholdProp_;      ///< only magnitudes in this intervall are used
    FloatIntervalProperty relativeMagnitudeThresholdProp_;      ///< debug output
    OptionProperty<FilterMode> filterModeProp_;             ///< filtering inside the dataset

    //streamline bundle settings
    BoolProperty detectPathlineBundlesProp_;              ///< determines if bundles shall be detected
    ProgressProperty bundleDetectionProgressProp_;          ///< progress bar for the bundle-detection thread
    FloatProperty maxAverageDistanceThresholdProp_;         ///< distance threshold for the bundle algorithm
    FloatProperty minNumPathlinesPerBundleProp_;          ///< bundle must contain more than this percentage of streamlines
    IntProperty resampleSizeProp_;                          ///< streamlines are resampled to this value

    //seed points
    IntProperty seedTimeProp_;          ///< time used to init the rand calcualtion
    ButtonProperty resetSeedProp_;      ///< resets the seedTimeProp_;

    //ports
    VolumeListPort volInport_;
    StreamlineListPort streamlineOutport_;

    //background thread
    ProcessorBackgroundThread<PathlineCreator>* backgroundThread_;///< thread used to calcualte the streamlines
    StreamlineList* streamlineListThreadOutput_;                    ///< thread stores calcualtion here
    VolumeRAM_3xFloat* referenceVolume_;
    DirtyFlag dirtyFlag_;                                           ///< flag holding which data should be recalculated
    bool backgroundThreadIsPathlineCreator_;                      ///< determines whether the active background thread is the streamlinecreator
};

}   // namespace

#endif  // VRN_PATHLINECREATOR_H
