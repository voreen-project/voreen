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

#ifndef VRN_COSMOLOGYPARTICLERENDERER_H
#define VRN_COSMOLOGYPARTICLERENDERER_H


#include "voreen/core/processors/renderprocessor.h"

#include "voreen/core/ports/renderport.h"

#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/boundingboxproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/matrixproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"



#include "voreen/core/interaction/camerainteractionhandler.h"

#include "../ports/cmparticleport.h"

namespace voreen {

class VRN_CORE_API CosmologyParticleRenderer : public RenderProcessor {

public:
    CosmologyParticleRenderer();
    ~CosmologyParticleRenderer();

    virtual Processor*  create() const       { return new CosmologyParticleRenderer(); }
    virtual std::string getClassName() const { return "CosmologyParticleRenderer";     }
    virtual std::string getCategory() const  { return "Viscontest201x";                }

protected:
    virtual void setDescriptions() { setDescription("Processor that draw particles as cirular sprites."); }
    virtual void process();

    virtual void initialize(); 
    virtual void deinitialize();

private:
    enum RenderMode{
        VELOCITY_TRANSFER_FUNC = 0,
        DIRECT_VELOCITY_COLOR  = 1,
        PHI_TRANSFER_FUNC      = 2,
		INT_EGY_TRANSFER_FUNC  = 3,
		MOL_WGT_TRANSFER_FUNC  = 4,
        MASS                   = 5,
		PARTICLE_TYPE		   = 6,
		SPH_DEN_TRANSFER_FUNC  = 7,
		SPH_LEN_TRANSFER_FUNC  = 8,
        MAX_RENDER_MODE
    };

    void setupBuffers();
    void changedRenderMode();
    void changeUseAlpha();
    void changedAlphaFactor();
    void changeUseRelativeRadius();
    void invalidateBuffers();

    bool isScalarRenderMode(RenderMode rendermode){
        return (rendermode == VELOCITY_TRANSFER_FUNC) || (rendermode == PHI_TRANSFER_FUNC) || (rendermode == INT_EGY_TRANSFER_FUNC) || (rendermode == MOL_WGT_TRANSFER_FUNC) || (rendermode == SPH_DEN_TRANSFER_FUNC) || (rendermode == SPH_DEN_TRANSFER_FUNC);
    }

    //-------------
    //  members
    //-------------
    RenderPort outport_;            ///< output of the modified image
    CMParticlePort inport_;
    OptionProperty<bool> useRelativeRadius_;
    FloatProperty radiusProp_;
    FloatProperty relativRadiusProp_;
    CameraProperty camera_;
    ShaderProperty shaderProp_;
    ShaderProperty shaderPropDVC_;
	ShaderProperty shaderPropPTP_;
    ShaderProperty shaderPropMASS_;
    TransFunc1DKeysProperty transFunc_;
    TransFunc1DKeysProperty massTransFunc_;
    OptionProperty<RenderMode> renderMode_;
    BoolProperty useAlpha_;
    FloatProperty alphaFactor_;
    FloatProperty timeStep_;
    ColorProperty massColor_;

    BoolProperty  parallelizeProp_;
    BoolProperty adaptCamera_;
    FloatBoundingBoxProperty universeDimensions_;
    FloatMat4Property universeMatrix_;
	StringProperty unitDisplayed_;


    BoolProperty visualizesDirection_;
	BoolProperty visualizesPType_;
    BoolProperty usesTf_;
    BoolProperty visualizesMass_;

    CameraInteractionHandler* cameraHandler_;

    GLuint vbo_[MAX_RENDER_MODE];
    GLuint posvbo_;
    GLuint vao_;

    int particleCount_;
    tgt::vec2 minmaxVel_;
    tgt::vec2 minmaxPhi_;
	//tgt::vec2 minmaxRho_;
	//tgt::vec2 minmaxHh_;
	tgt::vec2 minmaxMu_;
	tgt::vec2 minmaxUu_;

    tgt::Bounds bounds_;
    bool buffersInvalid_;

    const Volume* phiValues_;
    const Volume* velValues_;
	//const Volume* rhoValues_;
	//const Volume* hhValues_;
	const Volume* uuValues_;
	const Volume* muValues_;
	const Volume* ptpValues_;

    static const std::string loggerCat_; ///< category used in logging
};

} // namespace

#endif // VRN_COSMOLOGYPARTICLERENDERER_H
