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

#include "ensemblesamplepointmapping.h"

#include "voreen/core/io/serialization/serializer.h"
#include "voreen/core/io/serialization/deserializer.h"

namespace voreen {

/*************************************************
 *          EnsembleSamplePointMapping           *
 *************************************************/

EnsembleSamplePointMapping::EnsembleSamplePointMapping()
{
}

EnsembleSamplePointMapping::EnsembleSamplePointMapping(const EnsembleDataset& ensemble)
{
    // Add a POIList for every ensemble-member
    for (EnsembleMember member : ensemble.getMembers()) {
        const std::string& name = member.getName();
        
        if (samplePointMapping_.find(name) != samplePointMapping_.end()) {
            // Member-names not unique
            throw VoreenException("Member names in ensemble are not unique!");
            continue;
        }
        
        samplePointMapping_.insert({ name, POIList() });
    }
}

const POIList* EnsembleSamplePointMapping::getPOIList(const std::string& memberName) const
{
    auto res = samplePointMapping_.find(memberName);
    if (res != samplePointMapping_.end()) {
        return &(res->second);
    }

    return nullptr;
}

void EnsembleSamplePointMapping::setPoiList(const std::string& memberName, const POIList* poiList)
{
    auto res = samplePointMapping_.find(memberName);
    if (res == samplePointMapping_.end()) {
        // Unknown membername
        return;
    }

    samplePointMapping_.erase(memberName);
    samplePointMapping_.insert({ memberName, *poiList });
}


void voreen::EnsembleSamplePointMapping::addDefaultPOIPoint(const SamplePointConfig& config) {

    // Add Samplepoint for each member
    for (auto iter = samplePointMapping_.begin(); iter != samplePointMapping_.end(); ++iter) {
        const std::string memberName = iter->first;
        POIList& poiList = iter->second;
        POIGroupID gId = poiList.addGroup(config.name_);
        poiList.addPoint(config.defaultPos_, gId);
    }
}

void voreen::EnsembleSamplePointMapping::addDefaultPOIPoint(std::vector<SamplePointConfig>& config){
    for (auto it = config.begin(); it != config.end(); ++it) {
        addDefaultPOIPoint(*it);
    }
}

void EnsembleSamplePointMapping::serialize(Serializer& s) const {
    s.serialize("spMapping", samplePointMapping_);
}
void EnsembleSamplePointMapping::deserialize(Deserializer& s) {
    s.deserialize("spMapping", samplePointMapping_);
}

const std::vector<std::string> voreen::EnsembleSamplePointMapping::getMappedMemberNames() const{
    std::vector<std::string> names;
    for (auto const& entry : samplePointMapping_) {
        names.push_back(entry.first);
    }

    return names;
}

}