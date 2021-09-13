/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2014 University of Muenster, Germany.                        *
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

#ifndef VRN_CMHALOSUBSTRUCTURERENDERER_H
#define VRN_CMHALOSUBSTRUCTURERENDERER_H


#include "voreen/core/processors/renderprocessor.h"

#include "voreen/core/ports/renderport.h"

#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"


#include "voreen/core/interaction/camerainteractionhandler.h"
#include "../interaction/cmpropertyanimator.h"
#include "../interaction/cmselectionmanager.h"

#include "../ports/cmhaloport.h"

//use namespace voreen
namespace voreen {

/**
 * Renderprocessor that renders the currently selected halo and its direct satellite halos with correct radii.
 */
class VRN_CORE_API CMHaloSubstructureRenderer : public RenderProcessor {

public:
    /**
     * Constructor
     */
    CMHaloSubstructureRenderer();
    ~CMHaloSubstructureRenderer();

    //------------------------------------------
    //  Pure virtual functions of base classes
    //------------------------------------------
    virtual Processor* create() const { return new CMHaloSubstructureRenderer();     }
    virtual std::string getClassName() const { return "CMHaloSubstructureRenderer";  }
    virtual std::string getCategory() const  { return "Viscontest2015";         }

protected:

    virtual void setDescriptions() { setDescription("Renderprocessor that renders the currently selected halo and its direct satellite halos with correct radii."); }
    virtual void process();


    /**
     * Overwrites the base implementation of this function.
     * It is used to load the needed shader.
     * @see Processor
     */
    virtual void initialize();

    /**
     * Overwrites the base implementation of this function.
     * It is used to free the used shader.
     * @see Processor
     */
    virtual void deinitialize();

private:

    /**
     * Struct used represent halos when passed to openGL
     */
    struct HaloVertexLayout{
        tgt::vec3 pos;
        //tgt::vec3 color;
        float bodyValue;
        tgt::vec3 vel;
        tgt::vec3 angularMomenta;
        float radius;
        float spinParameter;
        GLuint haloID;
    };

    /**
     * Different modes used to render the halo bodies.
     */
    enum RenderMode{
        MASS_TRANSFER_FUNC = 0,
        DIRECT_VELOCITY_COLOR = 1,
        BLACK_HOLE_MASS_TRANSFER_FUNC = 2,
        BLACK_HOLE_SPIN_TRANSFER_FUNC = 3,
        SPHEROID_RADIUS_TRANSFER_FUNC = 4,
        SPHEROID_MASS_GAS_TRANSFER_FUNC = 5,
        SPHEROID_VELOCITY_TRANSFER_FUNC = 6,
        DISK_RADIUS_TRANSFER_FUNC = 7,
        DISK_MASS_GAS_TRANSFER_FUNC = 8,
        DISK_VELOCITY_TRANSFER_FUNC = 9
    };

    /**
     * Render the current halo and its direct satellite halos to screen.
     */
    void processHalo(const tgt::mat4& projectionMatrix, const tgt::mat4& viewMatrix);

    /**
     * Render the current halo and its direct satellite halos to the framebuffer
     * provided by selectionManager_.
     */
    void processHaloSelection(const tgt::mat4& projectionMatrix, const tgt::mat4& viewMatrix);

    /**
     * Creates and sets up opengl buffers for use.
     */
    void setupBuffers();

    /**
     * Will be called when the current render mode changes.
     */
    void changedRenderMode();

    /**
     * Will be called when the alpha factor property changes.
     */
    void changeUseAlpha();

    /**
     * Will be called when the currently selected halo changes.
     */
    void selectedHaloChanged();

    /**
     * Will be called when the time step changes.
     */
    void timeStepChanged();

    /**
     * Will be called when the halo data to render changes.
     */
    void inputDataChanged();

    /**
     * Get the value of the halo that is currently selected to decide the color of the halo body
     * according to renderMode_.
     */
    float getSelectedBodyValue(const CMHalo* h) const;

    /**
     * Prepare hostHalo and all its direct satellites for rendering, i.e. append them
     * to parameter verts.
     */
    void collectHalos(const CMHalo* hostHalo, std::vector<HaloVertexLayout>& verts) const;

    //-------------
    //  members
    //-------------
    /// Output of the rendered image
    RenderPort outport_;
    /// Halo data input
    CMHaloPort inport_;
    /// Halo radii will be divided by this factor
    FloatProperty radiusDivProp_;
    /// Thickness of the ring that represents the currently selected halo in rendering
    FloatProperty hostHaloRingThicknessProp_;
    /// The camera
    CameraProperty camera_;
    /// Whether or not to handle camera movement if the currently selected halo changes
    BoolProperty handleCameraMovement_;
    /// Initial distance between camera and currently selected halo on changed targets
    FloatProperty initialZoomDistanceProp_;
    /// Duration of animations (in ms) on changed selectedHaloIDProp_
    IntProperty animationDurationProp_;
    /// Shader used for rendering bodies when visualizing scalar halo attributes
    ShaderProperty shaderProp_;
    /// Shader used for rendering bodies when visualizing vector halo attributes
    ShaderProperty shaderPropDVC_;
    /// Shader used to render to selection frame buffer
    ShaderProperty shaderPropSelection_;
    /// TF used for visualizing scalar halo attributes
    TransFunc1DKeysProperty transFunc_;
    /// Render mode used for halo bodies
    OptionProperty<RenderMode> renderMode_;
    /// Whether or not to use alpha
    BoolProperty useAlpha_;
    /// Alpha factor
    FloatProperty alphaFactor_;
    /// ID of the currently selected halo
    IntProperty selectedHaloIDProp_;
    /// ID of a hovered over halo
    IntProperty mouseOverHaloIDProp_;
    /// Time step of the currently selected halo
    FloatProperty timeStep_;
    /// Currently renderd time step
    FloatProperty displayTimeStep_;

    /// Handles camera interaction
    CameraInteractionHandler* cameraHandler_;
    /// Used to handle selection via picking
    CMSelectionManager* selectionManager_;
    /// Used to handle time animation
    CMPropertyAnimator timeAnimator_;
    /// Used to handle camera movement animation
    CMPropertyAnimator cameraAnimator_;

    /// VBO that holds halo data
    GLuint vbo_;
    /// VAO used for rendering
    GLuint vao_;
    /// VAO used to render halo ids for selection
    GLuint selectionVao_;

    /// Number of halos to render
    int particleCount_;

    /// Histogram of the currently visualized halo properties
    const Volume* bodyValues_;
};

} // namespace

#endif // VRN_CMHALOSUBSTRUCTURERENDERER_H
