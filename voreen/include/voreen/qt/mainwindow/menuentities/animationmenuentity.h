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

#ifndef VRN_ANIMATIONMENUENTITY_H
#define VRN_ANIMATIONMENUENTITY_H

#include "voreenqtmenuentity.h"

namespace voreen {

class RenderTargetViewer;
class AnimationEditor;

/**
 * Menu entity for the animation editor
 */
class AnimationMenuEntity : public VoreenQtMenuEntity {

    friend class VoreenVEMainWindow;

public:
    AnimationMenuEntity();
    virtual ~AnimationMenuEntity();
    virtual QIcon getIcon() const {return QIcon(":/qt/icons/video.png");}
    virtual std::string getName() const {return "Animation";}
    virtual MenuCategory getMenuCategory() const {return MC_BASIC_TOOL;}
    virtual Qt::DockWidgetAreas getAllowedDockWidgetAreas() const {return Qt::NoDockWidgetArea;}
    virtual Qt::DockWidgetArea getInitialDockWidgetArea() const {return Qt::NoDockWidgetArea;}
    virtual bool getDefaultVisibilityInApplicationMode() const {return true;}
protected:
    virtual QWidget* createWidget() const;
    virtual void deinitialize();
private:
    mutable AnimationEditor* editor_;
};

} // namespace

#endif // VRN_ANIMATIONMENUENTITY_H
