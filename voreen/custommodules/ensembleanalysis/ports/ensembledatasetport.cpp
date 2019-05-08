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

#include "ensembledatasetport.h"

namespace voreen {

EnsembleDatasetPort::EnsembleDatasetPort(PortDirection direction, const std::string& id, const std::string& guiName, bool allowMultipleConnections, Processor::InvalidationLevel invalidationLevel)
    : GenericPort<EnsembleDataset>(direction, id, guiName, allowMultipleConnections, invalidationLevel) {
}

std::string EnsembleDatasetPort::getClassName() const {
    return "EnsembleDatasetPort";
}

Port* EnsembleDatasetPort::create(PortDirection direction, const std::string& id, const std::string& guiName) const {
    return new EnsembleDatasetPort(direction,id,guiName);
}

tgt::col3 EnsembleDatasetPort::getColorHint() const {
    return tgt::col3(24, 72, 124);
}

std::string EnsembleDatasetPort::getContentDescription() const {
    std::stringstream strstr;
    strstr << Port::getContentDescription();
    if(hasData()) {
        if (!getData()->getRuns().empty()) {
            strstr << std::endl << "Number of runs: " << getData()->getRuns().size();
            strstr << std::endl << "Number of unique channels: " << getData()->getUniqueChannels().size();
            strstr << std::endl << "Number of common Channels: " << getData()->getCommonChannels().size();
            strstr << std::endl << "Min Number of Time Steps: " << getData()->getMinNumTimeSteps();
            strstr << std::endl << "Max Number of Time Steps: " << getData()->getMaxNumTimeSteps();
            strstr << std::endl << "Min Time Step Duration: " << getData()->getMinTimeStepDuration();
            strstr << std::endl << "Max Time Step Duration: " << getData()->getMaxTimeStepDuration();
            strstr << std::endl << "Start Time: " << getData()->getStartTime();
            strstr << std::endl << "End Time: " << getData()->getEndTime();
            strstr << std::endl << "Max Total Duration: " << getData()->getMaxTotalDuration();

            tgt::vec3 llf = getData()->getRoi().getLLF();
            strstr << std::endl << "ROI LLF: (" << llf.x << ", " << llf.y << ", " << llf.z << ")";
            tgt::vec3 urb = getData()->getRoi().getURB();
            strstr << std::endl << "ROI URB: (" << urb.x << ", " << urb.y << ", " << urb.z << ")";
        }
        else
            strstr << std::endl << "Empty Ensemble Dataset";
    }
    return strstr.str();
}

std::string EnsembleDatasetPort::getContentDescriptionHTML() const {
    std::stringstream strstr;
    strstr << Port::getContentDescriptionHTML();
    if(hasData()) {
        if (!getData()->getRuns().empty()) {
            strstr << "<br>" << "Number of runs: " << getData()->getRuns().size();
            strstr << "<br>" << "Number of unique Channels: " << getData()->getUniqueChannels().size();
            strstr << "<br>" << "Number of common Channels: " << getData()->getCommonChannels().size();
            strstr << "<br>" << "Min Number of Time Steps: " << getData()->getMinNumTimeSteps();
            strstr << "<br>" << "Max Number of Time Steps: " << getData()->getMaxNumTimeSteps();
            strstr << "<br>" << "Min Time Step Duration: " << getData()->getMinTimeStepDuration();
            strstr << "<br>" << "Max Time Step Duration: " << getData()->getMaxTimeStepDuration();
            strstr << "<br>" << "Start Time: " << getData()->getStartTime();
            strstr << "<br>" << "End Time: " << getData()->getEndTime();
            strstr << "<br>" << "Max Total Duration: " << getData()->getMaxTotalDuration();

            tgt::vec3 llf = getData()->getRoi().getLLF();
            strstr << "<br>" << "ROI LLF: (" << llf.x << ", " << llf.y << ", " << llf.z << ")";
            tgt::vec3 urb = getData()->getRoi().getURB();
            strstr << "<br>" << "ROI URB: (" << urb.x << ", " << urb.y << ", " << urb.z << ")";
        }
        else
            strstr << "<br>" << "Empty Ensemble Dataset";
    }
    return strstr.str();
}

} // namespace
