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

#ifndef VRN_COSMOLOGYTIMESTEPOVERLAY_H
#define VRN_COSMOLOGYTIMESTEPOVERLAY_H

#include "voreen/core/processors/imageprocessor.h"

#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/fontproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "tgt/texturemanager.h"
#include "../ports/cmparticleport.h"
namespace voreen {

/**
 * Overlays a transfer function on top of the rendering.
 */
class CosmologyTimeStepOverlay : public ImageProcessor {
public:
    CosmologyTimeStepOverlay();
    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "CosmologyTimeStepOverlay"; }
    virtual std::string getCategory() const   { return "Image Processing"; }
    virtual CodeState getCodeState() const    { return CODE_STATE_TESTING;  }


protected:
    virtual void setDescriptions() {
        setDescription("Provides an overlay that renders a transfer function.");
    }

    virtual void process();
    virtual void initialize();
    virtual void deinitialize();
    
    virtual void onEvent(tgt::Event *ev) override;

private:
    void onChangeUsePixelCoordinates();

    RenderPort imageInport_;
    CMParticlePort inport_;
    RenderPort outport_;

	BoolProperty enabled_;
	FloatProperty timeStep_;
    FontProperty font_;
    IntVec2Property offset_;
    BoolProperty readOnly_;

	IntProperty currentSelectedView_;
	BoolProperty multiView_;
	BoolProperty overView_;
	BoolProperty view1_;
	BoolProperty view2_;
	BoolProperty view3_;
	BoolProperty view4_;

	

    bool isDragging_;


};

} // namespace

#endif
