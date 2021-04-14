/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#ifndef VRN_ENSEMBLETIMESTEPLINKER_H
#define VRN_ENSEMBLETIMESTEPLINKER_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"

#include "../ports/ensembledatasetport.h"

namespace voreen {

/**
 * Allows to link discrete time steps to the common and global time range of the ensemble.
 */
class VRN_CORE_API EnsembleTimeStepLinker : public Processor {
public:
    EnsembleTimeStepLinker();
    virtual Processor* create() const;

    virtual std::string getClassName() const    { return "EnsembleTimeStepLinker"; }
    virtual std::string getCategory() const     { return "Utility"; };
    virtual CodeState getCodeState() const      { return CODE_STATE_EXPERIMENTAL; }
    virtual bool isUtility() const              { return true; }

protected:
    virtual void setDescriptions() {
        setDescription("Allows to link discrete time steps to the common and global time range of the ensemble.\n"
                       "This processors main intention is to be used to link a selected time step in the parallel coordinates plot "
                       "to one or multiple EnsembleFilter processor time ranges (using common range).");
    }

    virtual void process() {};

private:

    void onTimeStepChange();

    EnsembleDatasetPort ensemblePort_;

    IntProperty inTimeStep_;
    FloatIntervalProperty outFromCommonTimeRange_;
    FloatIntervalProperty outFromGlobalTimeRange_;
};

} // namespace

#endif // VRN_ENSEMBLETIMESTEPLINKER_H
