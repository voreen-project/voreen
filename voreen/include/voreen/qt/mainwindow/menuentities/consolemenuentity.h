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

#ifndef VRN_CONSOLEMENUENTITY_H
#define VRN_CONSOLEMENUENTITY_H

#include "voreenqtmenuentity.h"

namespace voreen {

class ConsolePlugin;

/**
 * Menu entity for the debug console
 */
class ConsoleMenuEntity : public VoreenQtMenuEntity {

    friend class VoreenVEMainWindow;

public:
    ConsoleMenuEntity(ConsolePlugin* console);
    virtual ~ConsoleMenuEntity();
    virtual QIcon getIcon() const {return QIcon(":/qt/icons/console.png");}
    virtual std::string getName() const {return "Debug Console";}
    virtual std::string getShortCut() const {return "";}
    virtual MenuCategory getMenuCategory() const {return MC_BASIC_SIDEBAR;}
    virtual Qt::DockWidgetAreas getAllowedDockWidgetAreas() const {return Qt::BottomDockWidgetArea;}
    virtual Qt::DockWidgetArea getInitialDockWidgetArea() const {return Qt::BottomDockWidgetArea;}
    virtual bool getDefaultVisibilityInApplicationMode() const {return false;}
protected:
    virtual QWidget* createWidget() const;
private:
    mutable ConsolePlugin* consolePlugin_;
};

} // namespace

#endif // VRN_RENDERTARGETVIEWERENTITY_H
