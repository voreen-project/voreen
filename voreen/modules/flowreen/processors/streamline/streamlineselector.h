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

#ifndef VRN_STREAMLINESELECTOR_H
#define VRN_STREAMLINESELECTOR_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/utils/backgroundthread.h"

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/boundingboxproperty.h"
#include "voreen/core/properties/color/colorswitchproperty.h"
#include "voreen/core/properties/progressproperty.h"

#include "../../ports/streamlinelistport.h"
#include "voreen/core/ports/geometryport.h"

namespace voreen {

        class StreamlineSelectorBackgroundThread;
        class StreamlineList;

/**
 * Used to select(clip) streamlines according to a region of interest.
 * A Box can be defined, in which streamlines must/must not begin/end/intersect.
 *
 * @Note: It uses a background thread to handle changed parameters during calculation.
 */
class StreamlineSelector : public Processor, public PortObserver {

    friend class StreamlineSelectorBackgroundThread;

    /** Enum used to define intersection mode. */
    enum StreamlineSelectionMode {
        STREAM_BEGIN,
        STREAM_INTERSECT,
        STREAM_END
    };

public:
    StreamlineSelector();
    virtual ~StreamlineSelector();

    virtual Processor* create() const { return new StreamlineSelector(); }

    virtual std::string getCategory() const { return "Streamline Processing"; }
    virtual std::string getClassName() const { return "StreamlineSelector"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_EXPERIMENTAL; }
protected:
    virtual void setDescriptions() {
        setDescription("Used to select(clip) streamlines according to a region of interest.  Box can be defined, " \
                       "in which streamlines must (not) begin/end/intersect.");
        //ports
        streamlineInport_.setDescription("Streamlines, which should be filtered.");
        streamlineOutport_.setDescription("Filtered streamlines.");
        geometryOutport_.setDescription("LLF and URB of the ROI. Can be used with the BoundingBoxRenderer to visualize the ROI.");
        //properties
        enableProp_.setDescription("If disabled, the input is not modified.");
        selectStreamlinesProp_.setDescription("This button triggers a new streamline selection. This process can take a few minutes. " \
                                              "A progress bar will indicate the current progress. Any input or property change will abort the selection.");
        autoGenerateProp_.setDescription("If checked, an input or property change will trigger a recalculation. It can be used to calculate streamlines automaticaly " \
                                         "if the workspace is beeing opened. During parameter adjustment this option should not be checked.");
        clearSelectionProp_.setDescription("Clears the selection without changing the ROI.");
        insideProp_.setDescription("Defines, if streamlines must (not) be in the region of interest.");
        selectionModeProp_.setDescription("Defines, if streamlines must begin/end in the region of interest or must intersect it.");
        roiProp_.setDescription("The region of interest where streamlines must (not) begin/end or intersect.");
        showGeometryProp_.setDescription("Enables the geometry outport to visualize the current region of interest.");
        colorProp_.setDescription("The color of the ROI.");
        resetToLastUsedGeometryProp_.setDescription("Reset ROI to last used configuration.");
    }

    virtual void process();
    virtual bool usesExpensiveComputation() const {return true;}

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
    void inportHasChanged();
    /** Toggles the calculation button */
    void enableOnChange();
    /** Starts a new selection on button press. */
    void selectStreamlinesOnChange();
    /** Anything has been changed. */
    void anythingHasBeenChanged();
    /** Triggered by the button with the same name. */
    void clearSelectionOnChange();
    /** Triggered by the button with the same name. */
    void resetToLastUsedOnChange();

    //------------------
    //  Members
    //------------------
private:
    // ports
    StreamlineListPort streamlineInport_;
    StreamlineListPort streamlineOutport_;
    GeometryPort geometryOutport_;
    // properties
        // enable
    BoolProperty enableProp_;                                   ///< toggles the processor on and of
        // start configuration
    ButtonProperty selectStreamlinesProp_;                      ///< starts a new selection
    BoolProperty autoGenerateProp_;                             ///< starts a selection on any change
    ButtonProperty clearSelectionProp_;                         ///< clears the current selection
    ProgressProperty progressProp_;                             ///< used to show the progress in application mode
        // config
    OptionProperty<bool> insideProp_;                           ///< streamlines must be (not) inside the ROI
    OptionProperty<StreamlineSelectionMode> selectionModeProp_; ///< streamlines must start/end/intersect ROI
    FloatBoundingBoxProperty roiProp_;                          ///< ROI for selecting/clipping
        //roi representation (advanced)
    BoolProperty showGeometryProp_;                             ///< enables the geometry outport
    ColorSwitchProperty colorProp_;                             ///< color for handling de bounding box color
    ButtonProperty resetToLastUsedGeometryProp_;                    ///< reset to last used geometry

    //background thread
    StreamlineSelectorBackgroundThread* backgroundThread_;  ///< thread used to calcualte the streamlines
    StreamlineList* streamlineListThreadOutput_;            ///< thread stores calcualtion here
    StreamlineList* lastSelectedList_;                      ///< last selected list
    tgt::Bounds*    lastUsedGeometry_;                      ///< last used geometry (in voxel space)
    bool reselectStreamlines_;                              ///< process is only executed, if this is true
};


/**
 * Background thread used to select the stream lines.
 */
class VRN_CORE_API StreamlineSelectorBackgroundThread : public ProcessorBackgroundThread<StreamlineSelector> {
    friend class StreamlineSelector;
public:
    /** Constructor */
    StreamlineSelectorBackgroundThread(StreamlineSelector* processor, StreamlineList* output,
                                      bool inside, StreamlineSelector::StreamlineSelectionMode selectionMode, tgt::Bounds roi);
    /** Destructor */
    virtual ~StreamlineSelectorBackgroundThread();
protected:
    /** Main-Function used to calculate the streamlines */
    virtual void threadMain();
    /** Used to clean up */
    virtual void handleInterruption();
    //-----------------
    //  Helpers
    //-----------------
    /** Used for "begin" selection. */
    void selectBegin();
    /** Used for "end" selection. */
    void selectEnd();
    /** Used for "intersect" selection. */
    void selectIntersect();
private:
    //-----------
    //  Members
    //-----------
    //derived from properties
    StreamlineList* output_;     ///< output, which will be used by the processor
    bool inside_;                ///<  (not) in the roi?
    StreamlineSelector::StreamlineSelectionMode selectionMode_; ///< current selection modestreamline length must be in this interval
    tgt::Bounds roi_; ///< the region of selection interest
};

}   // namespace

#endif  // VRN_STREAMLINESELCETOR_H
