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

#ifndef VRN_PORTCONDITION_VOLUMELIST_H
#define VRN_PORTCONDITION_VOLUMELIST_H

#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include "voreen/core/ports/genericport.h"

namespace voreen {

/**
 * Base class for port conditions that check for a specific volume list.
 */
class VRN_CORE_API PortConditionVolumeList : public PortCondition {
public:
    PortConditionVolumeList(const std::string& description);
    virtual ~PortConditionVolumeList();

protected:
    virtual void setCheckedPort(const Port* checkedPort);
    const VolumeListPort* volumeListPort_;
};

// ------------------------------------------------------------------------------------------------

/**
 * Port condition that allows volume lists with a specific number of volumes.
 */
class VRN_CORE_API PortConditionVolumeListElementCount : public PortConditionVolumeList {
public:
    PortConditionVolumeListElementCount(size_t numVolumes);

    virtual bool acceptsPortData() const;

protected:
    size_t numVolumes_;
};

// ------------------------------------------------------------------------------------------------

/**
 * Port condition that allows volume lists containing volumes with identical meta data alltogether.
 */
class VRN_CORE_API PortConditionVolumeListEnsemble : public PortConditionVolumeList {
public:
    PortConditionVolumeListEnsemble();

    virtual bool acceptsPortData() const;
};

// ------------------------------------------------------------------------------------------------

/**
 * Port condition that allows volume lists that fulfill the passed condition for each contained volume.
 * Node: The condition takes ownership of the specified condition.
 */
class VRN_CORE_API PortConditionVolumeListAdapter : public PortConditionVolumeList {
public:
    PortConditionVolumeListAdapter(PortCondition* condition);

    virtual bool acceptsPortData() const;

protected:
    std::unique_ptr<PortCondition> condition_;
    std::unique_ptr<VolumePort> tmpPort_;
};


} // namespace

#endif // VRN_PORTCONDITION_VOLUMELIST_H
