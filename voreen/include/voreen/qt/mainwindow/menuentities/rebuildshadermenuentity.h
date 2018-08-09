/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#ifndef VRN_REBUILDSHADERMENUENTITY_H
#define VRN_REBUILDSHADERMENUENTITY_H

#include "voreenqtmenuentity.h"

#include <QObject>

namespace voreen {

/**
 * Menu entity for the rebuilding the shaders
 */
class RebuildShaderMenuEntity : public QObject, public VoreenQtMenuEntity {
    Q_OBJECT
    friend class VoreenQtMainWindow;

public:
    RebuildShaderMenuEntity();
    virtual ~RebuildShaderMenuEntity();
    virtual QIcon getIcon() const {return QIcon(":/qt/icons/rebuild.png");}
    virtual std::string getName() const {return "Rebuild Shaders";}
    virtual std::string getShortCut() const {return "F5";}
    virtual MenuCategory getMenuCategory() const {return MC_BASIC_TOOL;}
    virtual bool getDefaultVisibilityInApplicationMode() const {return false;}
protected:
    virtual QAction* createMenuAction();
protected slots:
    void rebuildShaders();
};

} // namespace

#endif // VRN_REBUILDSHADERMENUENTITY_H
