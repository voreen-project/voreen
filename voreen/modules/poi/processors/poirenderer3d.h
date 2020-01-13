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

#ifndef VRN_POIRENDERER3D_H
#define VRN_POIRENDERER3D_H

#include "voreen/core/processors/renderprocessor.h"
#include "poistorage.h"
#include "voreen/core/properties/stringexpressionproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/ports/genericcoprocessorport.h"
#include "voreen/core/ports/textport.h"
#include "voreen/core/interaction/camerainteractionhandler.h"
#include "voreen/core/ports/renderport.h"
#include "voreen/core/ports/volumeport.h"

namespace voreen {

/**
 * This renders POIs as 3d Spheres by creating every sphere from a quad in a shader
 */
class VRN_CORE_API POIRenderer3d : public RenderProcessor {
public:
    POIRenderer3d();

    virtual Processor* create() const;
    virtual std::string getClassName() const;
    virtual std::string getCategory() const;
    virtual bool isReady() const override;
protected:
    virtual void setDescriptions();
    virtual void process();

    void renderPOIs();

    virtual void initialize() override;
    
    // mouse events
    void mouseMoveEvent(tgt::MouseEvent* ev);
    void mouseClickEvent(tgt::MouseEvent* ev);
    void endSelecting(tgt::MouseEvent* ev);
    void startSelecting(tgt::MouseEvent* ev);

    void autodjustRadius();
    void renderSelection(tgt::vec2 begin, tgt::vec2 end);
private:
    GenericCoProcessorPort<POIStorage> cpPort_;
    VolumePort volumeInport_;
    RenderPort renderPort_; 

    // shaders
    ShaderProperty shader_; ///< shader for normal rendering of pois

    CameraProperty camera_;
    BoolProperty autoadjustRadius_;
    FloatProperty radius_; ///< radius of the poi spheres in world space
    FloatProperty radiusDivieder_;
    ColorProperty distanceMeasureColor_; ///< color of distance mesure line

    OptionProperty<POIStorage::SelectionMode> selectionMode_;
    std::unique_ptr<EventProperty<POIRenderer3d> > mouseMoveEventProp_;

    std::unique_ptr<CameraInteractionHandler> cameraHandler_;
    std::unique_ptr<EventProperty<POIRenderer3d> > mouseClickEventProp_;
    std::unique_ptr<EventProperty<POIRenderer3d> > mousePressEventProp_;
    std::unique_ptr<EventProperty<POIRenderer3d> > mouseReleasedEventProp_;

    /**
     * Gets the id of the POI at position in screen space.
     * @parm pos Position in screen space
     */
    POIPointID idAtPosition(tgt::ivec2 pos);    

    // Vertexbuffer for current points of interest
    void buildBuffers();
    GLuint vao_;
    GLuint vbo_;
    int primitiveCount_;
    bool shouldRebuildBuffers_;
    bool selecting_;
    tgt::vec2 selectionBegin_;
    tgt::vec2 selectionEnd_;
    
};
} // namespace

#endif
