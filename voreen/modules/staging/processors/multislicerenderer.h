/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
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

#ifndef VRN_MULTISLICERENDERER_H
#define VRN_MULTISLICERENDERER_H

#include "voreen/core/processors/volumerenderer.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/interaction/slicecamerainteractionhandler.h"

namespace voreen {

class MultiSliceRenderer : public VolumeRenderer {
public:
    MultiSliceRenderer();

    virtual ~MultiSliceRenderer();

    virtual std::string getCategory() const { return "Slice Rendering"; }
    virtual std::string getClassName() const { return "MultiSliceRenderer"; }
    virtual CodeState getCodeState() const { return CODE_STATE_TESTING; }
    virtual Processor* create() const { return new MultiSliceRenderer(); }

    virtual void initialize();
    virtual void deinitialize();

    virtual void process();

    virtual bool isReady() const;
protected:
    virtual void setDescriptions() {
        setDescription("Multi-modal slice-renderer, based on a texture with (world-space) sampling positions (think: entry-exit points). Use in combination with AlignedSliceProxyGeometry or SliceProxyGeometry.");
    }

    virtual std::string generateHeader(const tgt::GpuCapabilities::GlVersion* version = 0);
    virtual void compile();

private:
    tgt::Shader* raycastPrg_;

    CameraProperty camera_; //temporary?
    SliceCameraInteractionHandler* interactionHandler_;

    TransFunc1DKeysProperty transferFunc1_;
    TransFunc1DKeysProperty transferFunc2_;
    TransFunc1DKeysProperty transferFunc3_;
    TransFunc1DKeysProperty transferFunc4_;

    StringOptionProperty blendingMode1_;
    StringOptionProperty blendingMode2_;
    StringOptionProperty blendingMode3_;
    StringOptionProperty blendingMode4_;

    FloatProperty blendingFactor1_;
    FloatProperty blendingFactor2_;
    FloatProperty blendingFactor3_;
    FloatProperty blendingFactor4_;

    IntOptionProperty texFilterMode1_;      ///< filter mode for volume 1
    IntOptionProperty texFilterMode2_;      ///< filter mode for volume 2
    IntOptionProperty texFilterMode3_;      ///< filter mode for volume 3
    IntOptionProperty texFilterMode4_;      ///< filter mode for volume 4

    GLEnumOptionProperty texClampMode1_;    ///< clamp mode for volume 1
    GLEnumOptionProperty texClampMode2_;    ///< clamp mode for volume 2
    GLEnumOptionProperty texClampMode3_;    ///< clamp mode for volume 3
    GLEnumOptionProperty texClampMode4_;    ///< clamp mode for volume 4
    FloatProperty texBorderIntensity_;      ///< border intensity for all volumes

    VolumePort volumeInport1_;
    VolumePort volumeInport2_;
    VolumePort volumeInport3_;
    VolumePort volumeInport4_;

    RenderPort entryPort_;
    RenderPort outport_;
};

} // namespace voreen

#endif // VRN_MULTISLICERENDERER_H
