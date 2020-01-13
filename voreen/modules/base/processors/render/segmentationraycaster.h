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

#ifndef VRN_SEGMENTATIONRAYCASTER_H
#define VRN_SEGMENTATIONRAYCASTER_H

#include "voreen/core/processors/volumeraycaster.h"

#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/buttonproperty.h"

#include "voreen/core/ports/volumeport.h"

namespace voreen {

/**
 * This Raycaster can be used to render a segmented data set where every segment can have its own transfer function.
 * Currently only uint8 segment data sets are supported.
 */
class VRN_CORE_API SegmentationRaycaster : public VolumeRaycaster {
public:
    const static int TRANSFUNC_TEXTURE_SIZE = 1024;

    SegmentationRaycaster();
    virtual ~SegmentationRaycaster();

    virtual std::string getClassName() const  { return "SegmentationRaycaster"; }
    virtual std::string getCategory() const   { return "Raycasting"; }
    virtual CodeState getCodeState() const    { return CODE_STATE_STABLE; }
    virtual Processor* create() const         { return new SegmentationRaycaster(); }

    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        //properties (volume raycaster)
        VolumeRaycaster::setDescriptions();
        //general
        setDescription("This Raycaster can be used to render a segmented data set where every segment can have its own transfer function. \
                        Currently only uint8 segmentation data sets are supported. If no segmentation data set is connected, it can be \
                        used as a normal raycaster by using the \"Default Transfer Function\" if \"Apply Segmentation\" is disabled.");
        //ports
        volumeInport_.setDescription("Connect the original volume the segmentation has been done of. If no segemntation volume is connected, the \
                                      volume can be rendered normally by disabling \"Apply Segmentation\" and using the \"Default Transfer Function\"");
        segmentationInport_.setDescription("Connect the segmentation volume. Supported formats are uint8 at the moment.");
        entryPort_.setDescription("Connect the entry points, e.g., by the MeshEntryExitPoints processor.");
        exitPort_.setDescription("Connect the exit points, e.g., by the MeshEntryExitPoints processor.");
        outport1_.setDescription("First rendering result. This processor can render three images at the same time. For each outport \
                                  different compositing modes can be choosen.");
        outport2_.setDescription("Second rendering result.");
        outport3_.setDescription("Third rendering result.");
        //properties
        defaultTransFuncProp_.setDescription("Transfer function used to render the input volume. The function can be replaced for each \
                                              segment in the segmentation volume.");
        compositingMode2_.setDescription("Like Compositing, but for the second outport.");
        compositingMode3_.setDescription("Like Compositing, but for the third outport.");
        applySegmentationProp_.setDescription("If true, the segmentation color settings are applied. Otherwise, the input volume will be \
                                               rendered using the default transfer function.");
        segmentIndexProp_.setDescription("Current segment, which can be configurated.");
        isDefaultTFProp_.setDescription("If true, the default transfer function is used for this segment. Otherwise, the \
                                        function specified below will be used.");
        segmentTransFuncProp_.setDescription("Transfer function used for the currently selected segment, if the default transfer function \
                                              is not been used.");
        resetCurrentSegmentProp_.setDescription("Resets the current segment TF to the default one.");
        resetAllSegmentsProp_.setDescription("Resets all segment TFs to the default one.");
    }

    /** @see VolumeRaycaster */
    virtual void initialize() override;
    /** @see VolumeRaycaster */
    virtual void deinitialize() override;
    /** @see VolumeRaycaster */
    virtual void beforeProcess() override;
    /** @see VolumeRaycaster */
    virtual void process() override;
    /** @see VolumeRaycaster */
    virtual std::string generateHeader(const tgt::GpuCapabilities::GlVersion* version = 0) override;
    /** Rebuilds the shader. */
    virtual void compile();

private:
    //--------------------------------------//
    //  Handle Segmentation                 //
    //--------------------------------------//
    /**
     * This Vector contains the currently used transfer functions for each segment.
     * If the vector contains a null pointer, it is interpreted as using the default transfer function.
     */
    std::vector<TransFunc1DKeys*> segmentTransFuncVector_;  ///< vector since arrays cant be serialized
    /**
     * Used to color each segment during raycsting. Each segment TF is stored in Y-order in 3 rows thick slices.
     */
    tgt::Texture* segmentationTransFuncTex_;  ///< 2D texture that contains the segments' 1D transfer functions in row-order
    /**
     * Determines whether the segmentationTransFuncTex_ has to be uploaded before next use i.e. in process().
     * For performance reasons it is done lazy and not on TF change.
     */
    bool segmentationTransFuncTexValid_;
    /**
     * Domains of each transfer segment transfer function.
     */
    tgt::vec2 segmentationTransferFuncDomains_[256]; ///< must be the MAX_SEGMENTATIONS value (see .cpp)

    /**
     * Updates the entire segmentationTransFuncTex_ for each segment. It sets the segmentationTransFuncTexValid_ flag to false to enable the lazy upload.
     */
    void updateEntireSegmentationTransFuncTexture();
    /**
     * Updates only one segment of the segmentationTransFuncTex_ (3 rows). It sets the segmentationTransFuncTexValid_ flag to false to enable the lazy upload.
     */
    void updateSegmentationTransFuncTexture(int segment);

    /**
     * Resample a 1D Texture if insize and outsize are a multiple of each other
     */
    void resampleTexture1D(uint32_t *out, int outsize, const uint32_t *in, int insize);


    /** Serialize segmentTransFuncVector_. */
    virtual void serialize(Serializer& s) const override;
    /** Serialize segmentTransFuncVector_. */
    virtual void deserialize(Deserializer& d) override;

    //--------------------------------------//
    //  Callbacks                           //
    //--------------------------------------//
    /** Used to set up the default raycasting properties. */
    void adjustPropertyVisibilities();
    /** Called on port change. */
    void volumeInportOnChange();
    /** Called on property change. */
    void defaultTransFuncOnChange();
    /** Called on property change. */
    void applySegmentationOnChange();
    /** Called from reset button. */
    void resetCurrentSegment();
    /** Called from reset button. */
    void resetAllSegments();
    /** Called on property change. */
    void segmentIndexOnChange();
    /** Called on property change. */
    void segmentTransFuncOnChange();

    //--------------------------------------//
    //  Members                             //
    //--------------------------------------//
        //ports
    VolumePort segmentationInport_;     ///< port for the segmentation volume (supported are uint8/12/16)

        //properties
    //raycasting
        //see VolumeRaycaster
    //rendering
    TransFunc1DKeysProperty defaultTransFuncProp_;  ///< The transfer function property used as default tf

    //segmentation
    BoolProperty applySegmentationProp_;            ///< segmentation is only applied, it this is true
    OptionProperty<int> segmentIndexProp_;          ///< current segment, which will be modified
    BoolProperty isDefaultTFProp_;                  ///< flag, if the default function is used for this segemnt
    ButtonProperty resetCurrentSegmentProp_;        ///< resets the transfer function to the default one
    TransFunc1DKeysProperty segmentTransFuncProp_;  ///< TF used for the selected segment
    ButtonProperty resetAllSegmentsProp_;           ///< restes all segment functions to the default one
    //light
        //see VolumeRenderer
    //debug
    ShaderProperty shaderProp_;     ///< The shader property used by this raycaster
    CameraProperty cameraProp_;     ///< The used camera
};

} // namespace

#endif // VRN_SEGMENTATIONRAYCASTER_H
