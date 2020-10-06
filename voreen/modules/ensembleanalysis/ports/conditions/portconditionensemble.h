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

#ifndef VRN_PORTCONDITION_ENSEMBLE_H
#define VRN_PORTCONDITION_ENSEMBLE_H

#include "voreen/core/ports/conditions/portcondition.h"
#include "../ensembledatasetport.h"

namespace voreen {

class VRN_CORE_API PortConditionEnsemble : public PortCondition {
public:
    PortConditionEnsemble(const std::string& message);

    virtual ~PortConditionEnsemble();

protected:
    virtual void setCheckedPort(const Port* checkedPort);

    const EnsembleDatasetPort* ensemblePort_;
};


class VRN_CORE_API PortConditionEnsembleFieldName : public PortConditionEnsemble {
public:
    PortConditionEnsembleFieldName(const std::string& fieldName);
    PortConditionEnsembleFieldName(const std::vector<std::string>& fieldNames);
    virtual ~PortConditionEnsembleFieldName();

    virtual bool acceptsPortData() const;

protected:
    const std::vector<std::string> fieldNames_;
};

class VRN_CORE_API PortConditionEnsembleFieldCount : public PortConditionEnsemble {
public:
    PortConditionEnsembleFieldCount(size_t fieldCount);
    virtual ~PortConditionEnsembleFieldCount();

    virtual bool acceptsPortData() const;

protected:
    const size_t fieldCount_;
};

class VRN_CORE_API PortConditionEnsembleSingleTimeStep : public PortConditionEnsemble {
public:
    PortConditionEnsembleSingleTimeStep();
    virtual ~PortConditionEnsembleSingleTimeStep();

    virtual bool acceptsPortData() const;
};


} // namespace

#endif // VRN_PORTCONDITION_ENSEMBLE_H
