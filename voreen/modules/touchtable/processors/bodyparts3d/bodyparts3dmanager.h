/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#ifndef VRN_BODYPARTS3DMANAGER_H
#define VRN_BODYPARTS3DMANAGER_H

#include "voreen/core/ports/textport.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "modules/touchtable/datastructures/geometry/trianglemeshgeometrysinglecolor.h"
#include "modules/touchtable/ports/bodyparts3d/bodypartstextureport.h"
#include "voreen/core/processors/renderprocessor.h"
#include "voreen/core/interaction/idmanager.h"
#include "voreen/core/interaction/camerainteractionhandler.h"
#include "tgt/shadermanager.h"

namespace voreen{

class VRN_CORE_API BodyParts3DManager : public RenderProcessor {
public:
    BodyParts3DManager();
    virtual ~BodyParts3DManager();

    virtual Processor* create() const {return new BodyParts3DManager();}
    virtual std::string getClassName() const { return "BodyParts3DManager";}
    virtual std::string getCategory() const  { return "TBD";}
    virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL;}

    virtual bool isReady() const;

    void updatePropertyVisibilities();
    void adjustPickingRenderSize();
protected:

    virtual void setDescriptions(){
        setDescription("This is the Processor responsible for managing BodyParts3D");
    };

    virtual void initialize();
    virtual void beforeProcess();
    virtual void process();
    virtual void deinitialize();

    tgt::Texture* renderGeometryTexture(TriangleMeshGeometryBodyParts3D* geom);

    void onEvent(tgt::Event* e);
    void highlightGeometry(tgt::MouseEvent* e);
    void removeGeometry(tgt::MouseEvent* e);
    void addGeometry();
    void restore();
    void resetHighlight();

private:
    GeometryPort inPort_;
    BodyPartsTexturePort textureOutPort_;
    GeometryPort outPort_;
    TextPort textOutPort_;
    RenderPort pickingPort_;
    RenderPort privateRenderPort_;

    tgt::Shader* pickingShader_;
    tgt::Shader* textureShader_;

    CameraProperty camera_;
    CameraInteractionHandler* cameraHandler_;

    IDManager idManager_;

    IntVec2Property canvasSize_;
    BoolProperty action_;
    ButtonProperty revert_;

    // clipping
    BoolProperty enableClipping_;
    FloatVec3Property planeNormal_;
    FloatProperty planeDistance_;
    BoolProperty invertPlane_;

    // lighting
    ColorProperty lightAmbient_;
    ColorProperty lightDiffuse_;
    ColorProperty lightSpecular_;
    FloatProperty shininess_;

    int highlight_;

    TriangleMeshGeometryCollectionBodyParts3D* portData_;
    TriangleMeshGeometryCollectionBodyParts3D* removedMeshes_;

    std::vector<std::pair<tgt::Texture*,std::string> > textures_;

    static const std::string loggerCat_;
};

}

#endif //VRN_BODYPARTS3DMANAGER_H
