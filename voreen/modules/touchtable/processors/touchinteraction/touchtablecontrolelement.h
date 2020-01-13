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

#ifndef VRN_TOUCHTABLECONTROLELEMENT_H
#define VRN_TOUCHTABLECONTROLELEMENT_H

#include "tgt/vector.h"
#include "tgt/texture.h"

namespace voreen {


    /**
    * Class for control elements in a menu of a renderer    */
    class TouchTableControlElement {

    public:

        TouchTableControlElement();
        TouchTableControlElement(int radius, tgt::ivec2 pos, tgt::Texture* symbolTexture = 0, float opacity = 1.f, tgt::vec4 colorMod = tgt::vec4(0.f), tgt::vec4 propColor = tgt::vec4(0.f), int id = -1, bool usePropColor = false, bool isActivatable = true);
        virtual bool isActive() const;
        virtual void setActive(bool active);
        virtual tgt::Texture* getSymbolTexture() const;
        virtual void setSymbolTexture(tgt::Texture* symbolTexture);
        virtual float getOpacity() const;
        virtual void setOpacity(float opacity);
        virtual tgt::ivec2 getPos() const;
        virtual int getRadius() const;
        virtual tgt::vec4 getColorMod() const;
        virtual void setColorMod(tgt::vec4 colorMod);
        virtual int getId() const;

        virtual tgt::ivec2 getPosition() const;
        virtual void setPosition(tgt::ivec2 pos);
        virtual void setRadius(int radius);

        virtual bool isActivatable() const;

    protected:

        tgt::ivec2 pos_;
        int radius_;
        tgt::Texture* symbolTexture_;
        float opacity_;
        tgt::vec4 colorMod_;
        tgt::vec4 propColor_;
        bool usePropColor_;
        bool active_;
        int id_;
        bool isActivatable_;

    };

} // namespace

#endif // VRN_TOUCHTABLECONTROLELEMENT_H
