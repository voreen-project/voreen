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

#ifndef VRN_COSMOLOGYHALORENDERERER_H
#define VRN_COSMOLOGYHALORENDERERER_H


#include "voreen/core/processors/renderprocessor.h"

#include "voreen/core/ports/renderport.h"

#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/datastructures/volume/volume.h"



#include "voreen/core/interaction/camerainteractionhandler.h"
#include "../interaction/cmselectionmanager.h"
#include "../interaction/cmpropertyanimator.h"

#include "../ports/cmhaloport.h"

//use namespace voreen
namespace voreen {

/**
 * Sample render processor, which gray-scales an input image based on a user-defined parameter.
 * VRN_CORE_API is a macro needed for shared libs on windows (see voreencoreapi.h)
 */
class VRN_CORE_API CosmologyHaloRenderer : public RenderProcessor {

public:
    /**
     * Constructor
     */
    CosmologyHaloRenderer();
    ~CosmologyHaloRenderer();

    //------------------------------------------
    //  Pure virtual functions of base classes
    //------------------------------------------
    virtual Processor* create() const { return new CosmologyHaloRenderer();     }
    virtual std::string getClassName() const { return "CosmologyHaloRenderer";  }
    virtual std::string getCategory() const  { return "Viscontest2015";         }

protected:

    virtual void setDescriptions() { setDescription("Sample render processor that transforms its input image to gray scale."); }
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

    enum FocusMode {
        HALO = 0,
        OVERVIEW = 1,
    };

    void processHalo(const tgt::mat4& projectionMatrix, const tgt::mat4& viewMatrix, float alphaFactor, GLint first, GLsizei count);
    void processHaloSelection(const tgt::mat4& projectionMatrix, const tgt::mat4& viewMatrix);
    void processRotationAxis(const tgt::mat4& projectionMatrix, const tgt::mat4& viewMatrix, float alphaFactor, GLint first, GLsizei count);
    void processRotationOrbit(const tgt::mat4& projectionMatrix, const tgt::mat4& viewMatrix, float alphaFactor, GLint first, GLsizei count);
    void processSatelliteLinks(const tgt::mat4& projectionMatrix, const tgt::mat4& viewMatrix, float alphaFactor);
    void setupBuffers();
    void changedRenderMode();
    void changeUseAlpha();
    void selectedHaloChanged();
    void mouseOverHaloChanged();
    void timeStepChanged();
    void focusModeChanged();
    void inputDataChanged();
    float getSelectedBodyValue(const CMHalo* h) const;

    void collectHalosAndLinks(const CMHalo* hostHalo, int hostVertID, std::vector<HaloVertexLayout>& verts, std::vector<GLushort>& satelliteLinks) const;

    //-------------
    //  members
    //-------------
    RenderPort outport_;            ///< output of the modified image
    CMHaloPort inport_;
	BoolProperty enable_;
    FloatProperty radiusDivProp_;
    CameraProperty camera_;
    BoolProperty handleCameraMovement_;
    OptionProperty<FocusMode> focusMode_;
    FloatProperty initialZoomDistanceProp_;
    IntProperty animationDurationProp_;
    ShaderProperty shaderProp_;
    ShaderProperty shaderPropSelection_;
    ShaderProperty shaderPropRotation_;
    ShaderProperty shaderPropDVC_;
    ShaderProperty shaderPropSatelliteLinks_;
    ShaderProperty rotationAxisProp_;
    ShaderProperty rotationOrbitProp_;
    TransFunc1DKeysProperty transFunc_;
    TransFunc1DKeysProperty axisTransFunc_;
    TransFunc1DKeysProperty orbitTransFunc_;
    OptionProperty<RenderMode> renderMode_;
    BoolProperty useAlpha_;
    FloatProperty alphaFactor_;
    FloatProperty previewHaloAlphaFactor_;
    IntProperty selectedHaloIDProp_;
    IntProperty mouseOverHaloIDProp_;
    FloatProperty timeStep_;
    FloatVec3Property overviewFocusProp_;
    FloatProperty overviewZoomProp_;
    FloatProperty satelliteLinkWidthProp_;
    FloatProperty satelliteLinkDistanceProp_;
    BoolProperty  transferFunctionVisibleProp_;
    StringProperty bodyUnitProp_;

    CameraInteractionHandler* cameraHandler_;
    CMSelectionManager* selectionManager_;
    CMPropertyAnimator cameraAnimator_;
    CMPropertyAnimator timeAnimator_;

    GLuint vbo_;
    GLuint vao_;
    GLuint selectionVao_;
    GLuint linkIbo_;

    int particleCount_;
    int linkCount_;

    const Volume* velocityValues_;
    const Volume* bodyValues_;
    const Volume* spinValues_;
};

} // namespace

#endif // VRN_COSMOLOGYHALORENDERERER_H
