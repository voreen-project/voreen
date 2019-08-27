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

#ifndef VRN_INTERACTIVEPROJECTIONLABELING_H
#define VRN_INTERACTIVEPROJECTIONLABELING_H

#include "voreen/core/processors/renderprocessor.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/renderport.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/boundingboxproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "modules/bigdataimageprocessing/datastructures/lz4slicevolume.h"

#include "tgt/vector.h"

namespace voreen {

struct LabelProjection;

struct LabelGuard {
public:
    LabelGuard(LabelProjection& labelProjection);
    ~LabelGuard();
    uint8_t& at(tgt::svec3);
private:
    LabelProjection& labelProjection_;
};

struct LabelProjection {

    enum LabelVoxelState {
        UNLABELED = 0,
        FOREGROUND = 1,
        BACKGROUND = 2,
        SUGGESTED_FOREGROUND = 3,
    };

    void withLabels(std::function<void(VolumeAtomic<uint8_t>&)>);

    LabelProjection();
    LabelProjection(tgt::svec3 dimensions, tgt::svec3 dimensionMask, tgt::mat3 realToProjectedMat_);

    const VolumeAtomic<uint8_t>& labels() const {
        return labels_;
    }
    LabelGuard labels_mut() {
        return LabelGuard { *this };
    }
    void bindProjectionTexture();
    void bindLabelTexture();

    tgt::mat3 realToProjected() const {
        return realToProjectedMat_;
    }
    tgt::mat3 projectedToReal() const {
        return tgt::transpose(realToProjectedMat_);
    }


private:
    friend struct LabelGuard;
    void ensureTexturesPresent();
    size_t getVoxelIndex(tgt::svec3 pos3d);

    VolumeAtomic<uint8_t> labels_;
    boost::optional<tgt::Texture> projectionTexture_;
    boost::optional<tgt::Texture> labelTexture_;
    tgt::svec3 dimensionMask_;
    tgt::mat3 realToProjectedMat_;
};

class InteractiveProjectionLabeling : public RenderProcessor {
public:
    InteractiveProjectionLabeling();
    virtual ~InteractiveProjectionLabeling();

    virtual std::string getClassName() const         { return "InteractiveProjectionLabeling"; }
    virtual std::string getCategory() const       { return "Volume Processing"; }
    virtual std::string setDescriptions() const       { return "Volume Processing"; }
    virtual VoreenSerializableObject* create() const;
    virtual void setDescriptions() { setDescription( "TODO"); }
    virtual CodeState getCodeState() const        { return CODE_STATE_EXPERIMENTAL;   }
    virtual void process();
    virtual void adjustPropertiesToInput();
    virtual void onPortEvent(tgt::Event* e, Port* port);


private:
    void drawEvent(tgt::MouseEvent* e, LabelProjection& p);
    void renderToPort(RenderPort& port, LabelProjection& p);
    void syncLabels();
    void withOutputVolume(std::function<void(LZ4SliceVolume<uint8_t>&)>);

    VolumePort inport_;
    VolumePort labelVolume_;
    RenderPort xyProjectionOutport_;
    RenderPort xzProjectionOutport_;
    RenderPort yzProjectionOutport_;

    LabelProjection xy_;
    LabelProjection xz_;
    LabelProjection yz_;

    IntBoundingBoxProperty projectionRegion_;

    boost::optional<LZ4SliceVolume<uint8_t>> outputVolume_;

    EventProperty<InteractiveProjectionLabeling> mouseDrawEvent_;

    ShaderProperty shader_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_INTERACTIVEPROJECTIONLABELING_H
