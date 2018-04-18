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

#include "similaritydataport.h"

namespace voreen {

SimilarityDataPort::SimilarityDataPort(PortDirection direction, const std::string& id, const std::string& guiName, bool allowMultipleConnections, Processor::InvalidationLevel invalidationLevel)
    : GenericPort<SimilarityData>(direction, id, guiName, allowMultipleConnections, invalidationLevel) {
}

std::string SimilarityDataPort::getClassName() const {
    return "SimilarityDataPort";
}

Port* SimilarityDataPort::create(PortDirection direction, const std::string& id, const std::string& guiName) const {
    return new SimilarityDataPort(direction,id,guiName);
}

tgt::col3 SimilarityDataPort::getColorHint() const {
    return tgt::col3(25, 172, 204);
}

std::string SimilarityDataPort::getContentDescription() const {
    std::stringstream strstr;
    strstr << Port::getContentDescription();
    if(hasData()) {
        strstr << std::endl << "TODO";
    }
    return strstr.str();
}

std::string SimilarityDataPort::getContentDescriptionHTML() const {
    std::stringstream strstr;
    strstr << Port::getContentDescriptionHTML();
    if(hasData()) {
        strstr << "<br>" << "TODO";
    }
    return strstr.str();
}

} // namespace
