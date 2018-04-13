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

#ifndef VRN_DISTANCEMEASURE_H
#define VRN_DISTANCEMEASURE_H

#include "voreen/core/processors/imageprocessor.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/utils/stringutils.h"
#include "voreen/core/datastructures/geometry/glmeshgeometry.h"

#include "voreen/core/ports/volumeport.h"

#include "tgt/font.h"
#include "tgt/glmath.h"
#include "tgt/immediatemode/immediatemode.h"

namespace voreen {

class VRN_CORE_API DistanceMeasure : public ImageProcessor {
public:
    DistanceMeasure();
    ~DistanceMeasure();
    virtual Processor* create() const;

    virtual std::string getCategory() const  { return "Utility";         }
    virtual std::string getClassName() const { return "DistanceMeasure"; }
    virtual CodeState getCodeState() const   { return CODE_STATE_STABLE; }
    virtual bool isUtility() const           { return true; }

    virtual bool isReady() const;

    void measure(tgt::MouseEvent* e);

protected:
    virtual void setDescriptions() {
        setDescription(
                "Allows to interactively measure distances on the surface rendered volumes. "
                "This processor expects the rendered volume, a rendering of the first hitpoints "
                "(use FHP compositing in a SingleVolumeRaycaster), and the volume that is currently "
                "being rendered (for scaling information) as input."
                );
    }

    void process();

private:
    tgt::ivec2 clampToViewport(tgt::ivec2 mousePos);

    RenderPort imgInport_;
    RenderPort fhpInport_;
    VolumePort refInport_;
    RenderPort outport_;

    EventProperty<DistanceMeasure> mouseEventProp_;
    CameraProperty camera_;
    BoolProperty renderSpheres_;

    tgt::ivec2 mouseCurPos2D_;
    tgt::vec4 mouseCurPos3D_;
    tgt::ivec2 mouseStartPos2D_;
    tgt::vec4 mouseStartPos3D_;
    bool mouseDown_;

    float distance_;

    tgt::Font font_;
    GlMeshGeometryUInt16Normal mesh_;
    tgt::ImmediateMode::LightSource lightSource_;
    tgt::ImmediateMode::Material material_;
};

} // namespace

#endif // VRN_DISTANCEMEASURE_H
