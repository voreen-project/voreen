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

#ifndef VRN_NETWORKMODESETTINGSMENUENTITY_H
#define VRN_NETWORKMODESETTINGSMENUENTITY_H

#include "voreenqtmenuentity.h"

namespace voreen {

class NetworkModeConfigDialog;

/**
 * Menu entity for the network mode settings
 */
class NetworkModeSettingsMenuEntity : public VoreenQtMenuEntity {

    friend class VoreenQtMainWindow;

public:
    NetworkModeSettingsMenuEntity();
    virtual ~NetworkModeSettingsMenuEntity();
    virtual QIcon getIcon() const {return QIcon(":/qt/icons/network_settings.png");}
    virtual std::string getName() const {return "Network Editor Settings...";}
    virtual std::string getShortCut() const {return "";}
    virtual MenuCategory getMenuCategory() const {return MC_SETTINGS;}
    virtual bool getDefaultVisibilityInApplicationMode() const {return false;}
    //networksettings should not be displayed in toolbar
    virtual QAction* createToolbarAction() {return 0;}
protected:
    virtual QWidget* createWidget() const;
    virtual void initialize();
private:
    mutable NetworkModeConfigDialog* dialog_;
};

} // namespace

#endif // VRN_NETWORKMODESETTINGSMENUENTITY_H
