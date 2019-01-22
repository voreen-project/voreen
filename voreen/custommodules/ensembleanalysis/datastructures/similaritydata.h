/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2017 University of Muenster, Germany.                        *
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

#ifndef VRN_SIMILARITYDATA_H
#define VRN_SIMILARITYDATA_H

#include "voreen/core/voreencoreapi.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumelist.h"
#include "voreen/core/io/serialization/serializable.h"
#include "voreen/core/io/serialization/xmlserializer.h"
#include "voreen/core/io/serialization/xmldeserializer.h"
#include "voreen/core/ports/port.h"

#include "tgt/vector.h"

#include <map>

namespace voreen {

/**
 * Datastructure used to represent the actual field data in a
 */
class VRN_CORE_API SimilarityData : public DataInvalidationObservable {
public:

    /** Constructor */
    explicit SimilarityData();
    /** Destructor */
    ~SimilarityData();

public:

    //----------------
    //  Access
    //----------------
    const std::vector<std::vector<float>>& getData() const;
    const tgt::ivec3 getDimensions() const;
    const std::vector<std::string> getRuns() const;

private:

    //----------------
    //  Members
    //----------------

    std::vector<std::vector<float>> data_;

};

}   // namespace

#endif //VRN_SIMILARITYDATA_H
