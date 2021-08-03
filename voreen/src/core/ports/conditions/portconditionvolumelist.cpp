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

#include "voreen/core/ports/conditions/portconditionvolumelist.h"

#include "voreen/core/utils/stringutils.h"

namespace voreen {

PortConditionVolumeList::PortConditionVolumeList(const std::string& description)
    : PortCondition(description)
{}

PortConditionVolumeList::~PortConditionVolumeList()
{}

void PortConditionVolumeList::setCheckedPort(const Port* checkedPort) {
    if (!dynamic_cast<const VolumeListPort *>(checkedPort)) {
        LERRORC("voreen.PortConditionVolumeList", "Assigned port is not a volume list port");
    } else {
        volumeListPort_ = static_cast<const VolumeListPort*>(checkedPort);
    }
}

// ----------------------------------------------------------------------------

PortConditionVolumeListElementCount::PortConditionVolumeListElementCount(size_t numVolumes)
    : PortConditionVolumeList("Expected list of size: " + std::to_string(numVolumes))
    , numVolumes_(numVolumes)
{
}

bool PortConditionVolumeListElementCount::acceptsPortData() const {
    if (!volumeListPort_)
        return false;

    const VolumeList* volumeList = volumeListPort_->getData();
    return (volumeList && volumeList->size() == numVolumes_);
}

// ----------------------------------------------------------------------------

PortConditionVolumeListEnsemble::PortConditionVolumeListEnsemble()
    : PortConditionVolumeList("Ensemble expected (identical meta data)")
{
}

bool PortConditionVolumeListEnsemble::acceptsPortData() const {
    if (!volumeListPort_)
        return false;

    const VolumeList* volumeList = volumeListPort_->getData();
    if(!volumeList)
        return false;

    // We define an empty list to be fulfilling the condition.
    if(volumeList->empty())
        return true;

    // Retrieve reference information.
    std::string type            = volumeList->first()->getFormat();
    tgt::svec3 dimensions       = volumeList->first()->getDimensions();
    tgt::vec3 spacing           = volumeList->first()->getSpacing();
    tgt::vec3 offset            = volumeList->first()->getOffset();
    tgt::mat4 transformation    = volumeList->first()->getPhysicalToWorldMatrix();

    for(size_t i = 1; i < volumeList->size(); i++) {
        if (type != volumeList->at(i)->getFormat()) {
            return false;
        }

        if (dimensions != volumeList->at(i)->getDimensions()) {
            return false;
        }

        if (spacing != volumeList->at(i)->getSpacing()) {
            return false;
        }

        if (offset != volumeList->at(i)->getOffset()) {
            return false;
        }

        if (transformation != volumeList->at(i)->getPhysicalToWorldMatrix()) {
            return false;
        }
    }
    return true;
}

// ----------------------------------------------------------------------------

PortConditionVolumeListAdapter::PortConditionVolumeListAdapter(PortCondition* condition)
    : PortConditionVolumeList(condition->getDescription())
    , condition_(condition)
    , tmpPort_(new VolumePort(Port::OUTPORT, "tmp"))
{
    condition->setCheckedPort(tmpPort_.get());
}

bool PortConditionVolumeListAdapter::acceptsPortData() const {
    if (!volumeListPort_)
        return false;

    const VolumeList* volumeList = volumeListPort_->getData();
    if(!volumeList)
        return false;

    // We define an empty list to be fulfilling the condition.
    if(volumeList->empty())
        return true;

    for(size_t i = 0; i < volumeList->size(); i++) {
        tmpPort_->setData(volumeList->at(i), false);
        if(!condition_->acceptsPortData()) {
            return false;
        }
    }
    return true;
}

} // namespace
