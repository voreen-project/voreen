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

#ifndef VRN_RADARGLYPHRENDERER3D_H
#define VRN_RADARGLYPHRENDERER3D_H

//super class
#include "radarglyphrendererbase.h"
//ports
#include "voreen/core/ports/volumeport.h"
//properties
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/interaction/camerainteractionhandler.h"

namespace voreen {

class RadarGlyphRenderer3D : public RadarGlyphRendererBase
{
public:
    RadarGlyphRenderer3D();
    virtual ~RadarGlyphRenderer3D();

    virtual Processor* create() const { return new RadarGlyphRenderer3D(); }

    virtual std::string getClassName() const { return "RadarGlyphRenderer3D"; }
    virtual void process();
private:

    enum RenderMode3D {
        NORMAL_RADAR_GLYPHS_3D = 0
    };
    friend class OptionProperty<int>;

    //3D Control
    CameraProperty cameraProp_;
    CameraInteractionHandler* cameraHandler_;

    //onChange
    void renderModeOnChange();
    void inportOnChange();

    //calculate glyphs
    void fillGlyphVBO();
    //void fillIndexVBO();
};

}   // namespace

#endif  // VRN_RADARGLYPHRENDERER3D_H
