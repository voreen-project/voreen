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

#ifndef VRN_ENSEMBLESAMPLEPOINTMAPPING_H
#define VRN_ENSEMBLESAMPLEPOINTMAPPING_H

#include "voreen/core/io/serialization/serializable.h"

#include "../modules/ensembleanalysis/datastructures/ensembledataset.h"
#include "../modules/poi/datastructures/poilist.h"
#include "../utils/samplepointconfigloader.h"

#include "tgt/vector.h"

#include <vector>
#include <map>

namespace voreen {

class EnsembleSamplePointMapping : public Serializable{

public:
    EnsembleSamplePointMapping();

    /** Creates a new EnsembleSamplePointMapping, which contains an empty entry 
     *  for every member in the given Ensemble.
     */
    EnsembleSamplePointMapping(const EnsembleDataset& ensemble);

    /**
     *  Returns the POIList which is mapped to the member
     *  with the given member-name. If the member-name is not part of the mapping
     *  a nullptr is returned.
     */
    const POIList* getPOIList(const std::string& memberName) const;

    /**
     *  Assigns a copy of the given POIList to the member
     *  with the given member-name. If the member-name is not part of the mapping
     *  nothing happens
     */
    void setPoiList(const std::string& memberName, const POIList *poiList);

    /**
     * Adds a samplepoint to every member of the mapping.
     * The Points have the label and position specified by the passed SamplePointConfig-Object.
     * If a samplepoint with the same label as specified in 'config' is already mapped to a member, 
     * the already mapped samplepoint is not modified.
     */
    void addDefaultPOIPoint(const SamplePointConfig& config);

    /**
      * Calls addDefaultSamplePoint for every element in the vector.
      */
    void addDefaultPOIPoint(std::vector<SamplePointConfig>& config);

    /**
     *  Returns a list of all mapped member-names
     **/
    const std::vector<std::string> getMappedMemberNames() const;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);
private:
    // Maps the name of ensemble members to a vector of sample points
    std::map<std::string, POIList> samplePointMapping_;
};

}

#endif  