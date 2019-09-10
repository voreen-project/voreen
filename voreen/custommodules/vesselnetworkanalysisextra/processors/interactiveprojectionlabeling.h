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
    void set(size_t x, size_t y, tgt::svec2 range, uint8_t val);
private:
    LabelProjection& labelProjection_;
};

struct LabelProjection {
    LabelProjection();

    const VolumeAtomic<uint8_t>& projection() const {
        return projection_;
    }
    LabelGuard projection_mut() {
        return LabelGuard { *this };
    }
    void bindTexture();
private:
    friend struct LabelGuard;
    void ensureTexturesPresent();

    VolumeAtomic<uint8_t> projection_;
    boost::optional<tgt::Texture> projectionTexture_;
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
    virtual void initialize();
    virtual void deinitialize();


private:
    void updateSizes();
    void renderOverlay();
    void renderProjection();
    void withOutputVolume(std::function<void(LZ4SliceVolume<uint8_t>&)>);

    void projectionEvent(tgt::MouseEvent* e);
    void overlayEvent(tgt::MouseEvent* e);

    VolumePort inport_;
    VolumePort labelVolume_;
    RenderPort overlayInput_;
    RenderPort overlayOutput_;
    RenderPort projectionOutput_;

    RenderPort fhp_;
    RenderPort lhp_;

    boost::optional<LZ4SliceVolume<uint8_t>> outputVolume_;

    FloatProperty interpolationValue_;

    static const std::string loggerCat_;

    tgt::Shader* copyShader_;
    ShaderProperty projectionShader_;

    std::vector<tgt::vec2> displayLine_;
};

} // namespace voreen

#endif // VRN_INTERACTIVEPROJECTIONLABELING_H
