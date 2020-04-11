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

#ifndef VRN_ENSEMBLEANALYSISMODULE_H
#define VRN_ENSEMBLEANALYSISMODULE_H

#include "voreen/core/voreenmodule.h"
#include "voreen/core/properties/boolproperty.h"

namespace voreen {

class EnsembleAnalysisModule: public VoreenModule {

public:
    EnsembleAnalysisModule(const std::string& modulePath);

    virtual std::string getDescription() const {
        return "Module for ensemble analysis, initiated for the SciVis Contest 2018";
    }

    /**
     * Enables or disables forcing a disk representation of loaded volumes, even
     * if they have no native one. If enabled, a hdf5 disk representation will be used.
     * Default: enabled
     */
    void setForceDiskRepresentation(bool enabled);

    /**
     * Returns whether disk representations are forced.
     */
    bool getForceDiskRepresentation() const;

    /**
     * Returns the global instance of this class.
     *
     * @note Does not create the instance. If the module class has not been
     *       instantiated yet, the null pointer is returned.
     */
    static EnsembleAnalysisModule* getInstance();

private:

    static EnsembleAnalysisModule* instance_;

    BoolProperty forceDiskRepresentation_;
};

} // namespace

#endif
