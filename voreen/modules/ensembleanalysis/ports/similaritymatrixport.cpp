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

#include "similaritymatrixport.h"

namespace voreen {

SimilarityMatrixPort::SimilarityMatrixPort(PortDirection direction, const std::string& id, const std::string& guiName, bool allowMultipleConnections, Processor::InvalidationLevel invalidationLevel)
        : GenericPort<SimilarityMatrixList>(direction, id, guiName, allowMultipleConnections, invalidationLevel) {
}

std::string SimilarityMatrixPort::getClassName() const {
    return "SimilarityMatrixPort";
}

Port* SimilarityMatrixPort::create(PortDirection direction, const std::string& id, const std::string& guiName) const {
    return new SimilarityMatrixPort(direction,id,guiName);
}

tgt::col3 SimilarityMatrixPort::getColorHint() const {
    return tgt::col3(150, 72, 190);
}

std::string SimilarityMatrixPort::getContentDescription() const {
    std::stringstream strstr;
    strstr << Port::getContentDescription();
    if(hasData()) {
        strstr << std::endl << "Size: " << getData()->getSize();
#ifdef VRN_DEBUG
        strstr << std::endl << "Hash: " << getData()->getHash();
#endif
    }
    return strstr.str();
}

std::string SimilarityMatrixPort::getContentDescriptionHTML() const {
    std::stringstream strstr;
    strstr << Port::getContentDescriptionHTML();
    if(hasData()) {
        strstr << "<br>" << "Size: " << getData()->getSize();
#ifdef VRN_DEBUG
        strstr << "<br>" << "Hash: " << getData()->getHash();
#endif
    }
    return strstr.str();
}

} // namespace
