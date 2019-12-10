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

#ifndef VRN_PBLINKCONTROL_H
#define VRN_PBLINKCONTROL_H

#include "voreen/core/processors/processor.h"

#include "modules/flowreen/processors/streamlinerenderer3d.h"

#include "voreen/core/ports/volumeport.h"

#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/boolproperty.h"

namespace voreen {

    /**
     * Utility processor for linking the string option of the 3d streamline rendering to display the different overlays 
     * (magnitude transfunc overlay or streamline direction color coding overlay).
     */
class PBLinkControl : public Processor {
public:
    PBLinkControl();
    virtual Processor* create() const;
    virtual std::string getClassName() const { return "PBLinkControl";  }
    virtual std::string getCategory() const  { return "Utility";     }
    virtual CodeState getCodeState() const   { return CODE_STATE_TESTING; }

protected:
    virtual void setDescriptions() {
        setDescription("Use to configure (do magic) in the streamline workspaces. <br> \
                        The overlay bool properties should be linked to all <i>FlowDirectionOverlay</i>s.");
    }

    virtual void process();

    /**
     * Callback for overlay settings.
     */
    virtual void updateOverlayCallback();


    //--------------//
    //  Properties  //
    //--------------//
    //overlay
    OptionProperty<StreamlineRenderer3D::StreamlineColorCoding> colorProp_; ///< toogle overlay options
        //hidden properties
    BoolProperty renderMagnitudeOverlay_;
    BoolProperty renderColorCodingOverlay_;

    static const std::string loggerCat_; //< static member for output LERROR etc...
};

} //namespace voreen

#endif
