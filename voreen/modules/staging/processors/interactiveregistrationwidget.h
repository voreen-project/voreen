/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#ifndef VRN_INTERACTIVEREGISTATIONWIDGET_H
#define VRN_INTERACTIVEREGISTATIONWIDGET_H

#include <string>
#include "voreen/core/processors/renderprocessor.h"

#include "voreen/core/ports/renderport.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/textport.h"

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/properties/matrixproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/colorproperty.h"

namespace voreen {

class Volume;

class VRN_CORE_API InteractiveRegistrationWidget : public RenderProcessor {
public:
    InteractiveRegistrationWidget();
    virtual ~InteractiveRegistrationWidget();
    virtual Processor* create() const;

    virtual std::string getCategory() const   { return "Volume Processing"; }
    virtual std::string getClassName() const  { return "InteractiveRegistrationWidget"; }
    virtual CodeState getCodeState() const    { return CODE_STATE_TESTING; }

    virtual void onEvent(tgt::Event* e);
    virtual bool isReady() const;
protected:
    virtual void setDescriptions() {
        setDescription("Renders a widget to rotate and translate volume for co-registration purposes.");
    }

    virtual void process();

    virtual void initialize();
    virtual void deinitialize();

    void translate(tgt::vec3 v);
    void rotate(tgt::vec3 center, tgt::vec3 v, float angle);

    void planeChanged();
    void centerPoint();

private:
    RenderPort inport_;
    RenderPort outport_;
    RenderPort pickingPort_;
    TextPort textPort_;

    FloatMat4Property transformMatrix_;
    FloatVec3Property point_;
    FloatVec3Property plane_;
    FloatProperty planeDist_;
    BoolProperty render_;
    CameraProperty camera_;
    FloatProperty sphereRadius_;
    FloatProperty ringRadius_;
    ColorProperty ringColor_;
    ColorProperty sphereColor_;
    ButtonProperty centerPoint_;

    tgt::Shader* copyShader_;
    tgt::Shader* sphereShader_;
    tgt::ivec2 lastCoord_;
    tgt::vec3 startDragCoord_;
    tgt::vec3 curDragCoord_;
    float rotAngle_;
    int mouseDown_;

    static const std::string loggerCat_; ///< category used in logging
};

}   //namespace

#endif
