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

#ifndef VRN_CMHALODISTANCEOVERLAY_H
#define VRN_CMHALODISTANCEOVERLAY_H

#include "voreen/core/processors/imageprocessor.h"

#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/fontproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "tgt/texturemanager.h"
#include "../ports/cmhaloport.h"

namespace voreen {

/**
 * Overlays a simple ruler like structure between (the position of) the currently selected halo
 * and (the position of) the halo hovered over with the mouse. The distance between those halos
 * will be displayed as well.
 */
class CMHaloDistanceOverlay : public ImageProcessor {
public:
    CMHaloDistanceOverlay();
    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "CMHaloDistanceOverlay"; }
    virtual std::string getCategory() const   { return "Image Processing"; }
    virtual CodeState getCodeState() const    { return CODE_STATE_TESTING;  }


protected:
    virtual void setDescriptions() {
        setDescription("Provides an overlay that visualises the distance between two halos.");
    }

    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

private:

    /**
     * Creates and sets up opengl buffers for use.
     */
    void setupBuffers();

    /**
     * Used to update interal datastructures if the hovered over or the selected halo changes.
     */
    void halosChanged();

    /// Inport for the image that this overlay will be put upon.
    RenderPort imageInport_;
    /// Inport for halo data
    CMHaloPort inport_;
    /// Outport for the inport image overlayed with halo distances
    RenderPort outport_;

    /// The camera (should be the same camera used to generate the inport image)
    CameraProperty camera_;
    /// Font used to render distance information
    FontProperty font_;
    /// ID of the selected halo
    IntProperty selectedHaloIDProp_;
    /// ID of the hovered over halo
    IntProperty mouseOverHaloIDProp_;
    /// Shader used for rendering the overlay
    ShaderProperty shaderProp_;
    /// Holds the 3d positions of the selected and hovered over halo
    GLuint vbo_;
    /// vao used for rendering
    GLuint vao_;
};

} // namespace

#endif
