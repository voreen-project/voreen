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

#ifndef VRN_VOLUMEPICKING_H
#define VRN_VOLUMEPICKING_H

#include "voreen/core/processors/renderprocessor.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/datastructures/geometry/pointlistgeometry.h"

#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/vectorproperty.h"

#include "voreen/core/ports/volumeport.h"

namespace voreen {

class VRN_CORE_API VolumePicking : public RenderProcessor {
public:
    VolumePicking();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "VolumePicking";    }
    virtual std::string getCategory() const  { return "Utility";          }
    virtual CodeState getCodeState() const   { return CODE_STATE_TESTING; }
    virtual bool isUtility() const           { return true; }

    virtual void invalidate(int inv = INVALID_RESULT);

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& d);

protected:
    virtual void setDescriptions() {
        setDescription("Converts 2D mouse coordinates into a selectable volume coordinate system. The picked 3D positions are put out as PointListGeometry. Additionally, the last picked position is provided through a property. The RenderPort needs to be connected to a first-hit rendering of the picked volume (currently expecting normalized texture coordinates).");
    }

    virtual void process();

private:
    tgt::vec3 getPickedPosition(tgt::ivec2 mousePosition);
    PointListGeometryVec3* getTransformedPositions() const;

    void mouseEvent(tgt::MouseEvent* e);
    void clearLast();
    void clearAll();
    void updateProperties();

    RenderPort inportFirstHit_;
    VolumePort inportVolume_;
    GeometryPort outportPicked_;

    StringOptionProperty outCoordinateSystem_;
    ButtonProperty clearLast_;
    ButtonProperty clearAll_;
    IntProperty numPicked_;
    FloatVec3Property lastPicked_;

    EventProperty<VolumePicking> pickEvent_;

    std::vector<tgt::vec3> pickedPositions_;
    tgt::ivec2 lastMousePosition_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif //VRN_VOLUMEPICKING_H
