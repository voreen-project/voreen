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

#ifndef VRN_APPLICATIONMODECONFIGMENUENTITY_H
#define VRN_APPLICATIONMODECONFIGMENUENTITY_H

#include "voreenqtmenuentity.h"

namespace voreen {

class ApplicationModeConfigDialog;

/**
 * Menu entity for the application mode settings
 */
class ApplicationModeConfigMenuEntity : public VoreenQtMenuEntity {

    friend class VoreenVEMainWindow;

public:
    ApplicationModeConfigMenuEntity();
    virtual ~ApplicationModeConfigMenuEntity();
    virtual QIcon getIcon() const {return QIcon(":/qt/icons/appmodeconfig.png");}
    virtual std::string getName() const {return "Application Mode Config";}
    virtual std::string getShortCut() const {return "";}
    virtual MenuCategory getMenuCategory() const {return MC_BASIC_TOOL;}
    virtual bool getDefaultVisibilityInApplicationMode() const {return false;}
protected:
    virtual QWidget* createWidget() const;
    void deinitialize();
private:
    mutable ApplicationModeConfigDialog* dialog_;
};

} // namespace

#endif // VRN_APPLICATIONMODECONFIGMENUENTITY_H
