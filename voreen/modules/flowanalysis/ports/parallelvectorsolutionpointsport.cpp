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

#include "parallelvectorsolutionpointsport.h"
#include <sstream>

namespace voreen {

ParallelVectorSolutionPointsPort::ParallelVectorSolutionPointsPort(PortDirection direction, const std::string& id, const std::string& guiName, bool allowMultipleConnections, Processor::InvalidationLevel invalidationLevel)
    : GenericPort<ParallelVectorSolutions>(direction, id, guiName, allowMultipleConnections, invalidationLevel)
{
}

Port *ParallelVectorSolutionPointsPort::create(PortDirection direction, const std::string& id, const std::string& guiName) const {
    return new ParallelVectorSolutionPointsPort(direction, id, guiName);
}
std::string ParallelVectorSolutionPointsPort::getClassName() const {
    return "ParallelVectorSolutionPointsPort";
}

tgt::col3 ParallelVectorSolutionPointsPort::getColorHint() const {
    return GenericPort::getColorHint();
}

std::string ParallelVectorSolutionPointsPort::getContentDescription() const {
    auto stream = std::stringstream();
    stream << Port::getContentDescription();

    if (hasData()) {
        stream << std::endl << "Size: " << getData()->solutions.size();
    }

    return stream.str();
}

std::string ParallelVectorSolutionPointsPort::getContentDescriptionHTML() const {
    auto stream = std::stringstream();
    stream << Port::getContentDescriptionHTML();

    if (hasData()) {
        stream << "<br/>" << "Size: " << getData()->solutions.size();
    }

    return stream.str();
}

} // namespace voreen
