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

#ifndef VRN_VOREENMODULEQT_H
#define VRN_VOREENMODULEQT_H

#include "voreen/core/voreenmodule.h"

#include <string>
#include <vector>
#include "voreen/core/utils/exception.h"
#include "voreen/qt/voreenqtapi.h"

namespace voreen {

    class VoreenQtMenuEntity;

/**
 * Marker interface for Voreen Qt module classes.
 *
 * A qt module class is responsible for registering the
 * module's widgets at runtime:
 *  - ProcessorWidgets
 *  - PropertyWidgets
 *
 * The registration is to be done in the constructor.
 *
 * @see http://www.voreen.org/402-Adding-a-Module.html
 *
 */
class VRN_QT_API VoreenModuleQt : public VoreenModule {

    friend class VoreenApplicationQt;

public:
    VoreenModuleQt(const std::string& modulePath);
    ~VoreenModuleQt();

    virtual std::string getDescription() const { return "Provides Qt-based widgets."; }

    /**
     * Returns the module's custom menu entities.
     *
     * @see registerCustomMenuEntity
     */
    const std::vector<VoreenQtMenuEntity*>& getCustomMenuEntities() const;
protected:
    /**
     * Used to register custom menu entities.
     * @param entity entity to be registered
     */
    void registerCustomMenuEntity(VoreenQtMenuEntity* entity);

    std::vector<VoreenQtMenuEntity*> customMenuEntities_; ///< contains all registered entities

    static const std::string loggerCat_;
};

} // namespace

#endif // VRN_VOREENMODULEQT_H
