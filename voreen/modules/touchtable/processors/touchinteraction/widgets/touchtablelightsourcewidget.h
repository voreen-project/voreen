/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#ifndef VRN_TOUCHTABLELIGHTSOURCEWIDGET_H
#define VRN_TOUCHTABLELIGHTSOURCEWIDGET_H

#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/lightsourceproperty.h"

#include "../touchtablecontrolelement.h"
#include "../touchtableoverlay.h"
#include "../touchtableslider.h"

#include "touchtablemenuwidget.h"

namespace voreen {

    class TouchTableLightSourceWidget : public TouchTableMenuWidget {
    public:
        TouchTableLightSourceWidget();

        Processor* create() const{
            return new TouchTableLightSourceWidget();
        }

        virtual bool isReady() const{
            return true;
        }

        virtual std::string getClassName() const{
            return "TouchTableLightSourceWidget";
        }

        virtual void handleTouchPoints(const std::deque<tgt::TouchPoint>& tp);

    protected:

        virtual void initialize();

        virtual void deinitialize();

        virtual void renderComponents();

        virtual void process();

        virtual void updateComponents();

        void updateProperties();

        bool hitsQuad(const tgt::ivec2& pos);

        void adjustLightToCamera();

        void updateLightSourcePosition(const tgt::ivec2& pos);


        tgt::vec3 oldCameraNormal_;

        int sliderID_;
        int quadID_;

        tgt::ivec2 quadUR_;
        tgt::ivec2 quadLL_;

        TouchTableControlElement followCamera_;
        TouchTableControlElement invertLight_;
        TouchTableControlElement lightElem_;                ///< control element for showing the light source

        BoolProperty renderLightSource_;                    ///< determines if the light source is rendered and may be moved
        FloatProperty shininess_;
        ColorProperty diffuseLight_;
        ColorProperty ambientLight_;
        ColorProperty specularLight_;
        CameraProperty camera_;
        LightSourceProperty lightPosition_;
        TouchTableSlider lightDistance_;
        tgt::Texture* followCameraTex_;
        tgt::Texture* invertLightTex_;

        tgt::Texture* lightTex_;

        static const std::string loggerCat_;
    };

} // namespace

#endif // VRN_TOUCHTABLELIGHTSOURCEWIDGET_H

