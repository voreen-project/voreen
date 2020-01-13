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

#ifndef VRN_RENDERTARGETVIEWERMENUENTITY_H
#define VRN_RENDERTARGETVIEWERMENUENTITY_H

#include "../voreenqtmenuentity.h"

namespace voreen {

class RenderTargetViewer;

/**
 * Menu entity for the render target viewer
 */
class RenderTargetViewerMenuEntity : public VoreenQtMenuEntity {

    friend class VoreenQtMainWindow;

public:
    RenderTargetViewerMenuEntity();
    virtual ~RenderTargetViewerMenuEntity();
    virtual QIcon getIcon() const {return QIcon(":/qt/icons/rendertargets.png");}
    virtual std::string getName() const {return "Render Target Viewer";}
    virtual std::string getShortCut() const {return "Ctrl+T";}
    virtual MenuCategory getMenuCategory() const {return MC_BASIC_TOOL;}
    virtual bool getDefaultVisibilityInApplicationMode() const {return false;}
protected:
    virtual QWidget* createWidget() const;
    virtual void initialize();
    virtual void deinitialize();
private:
    mutable RenderTargetViewer* renderTargetViewer_;
};

} // namespace

#endif // VRN_RENDERTARGETVIEWERENTITY_H
