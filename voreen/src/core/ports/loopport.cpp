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

#include "voreen/core/ports/loopport.h"

namespace voreen {


LoopPort::LoopPort(PortDirection direction, const std::string& name, const std::string& guiName,
                   bool allowMultipleConnections, Processor::InvalidationLevel invalidationLevel)
    : Port(direction, name, guiName, allowMultipleConnections, invalidationLevel)
{
    setLoopPort(true);
}

std::string LoopPort::getContentDescription() const {
    std::stringstream strstr;
    strstr << getGuiName() << std::endl
           << "Type: " << getClassName()
           << "Loop-Iterations: " << getNumLoopIterations();
    return strstr.str();
}

std::string LoopPort::getContentDescriptionHTML() const {
    std::stringstream strstr;
    strstr << "<center><font><b>" << getGuiName() << "</b></font></center>"
           << "Type: " << getClassName()
           << "Loop-Iterations: " << getNumLoopIterations();
    return strstr.str();
}

} // namespace
