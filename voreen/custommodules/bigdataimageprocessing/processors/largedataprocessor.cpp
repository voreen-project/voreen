/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
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

#include "largedataprocessor.h"

#include "voreen/core/datastructures/callback/lambdacallback.h"

namespace voreen {

const char* ProcessingInterruptedException::what() const throw() {
    return "Computation interrupted.";
}


InterruptableProgressReporter::InterruptableProgressReporter()
    : supervisingReporters_()
    , interrupted_(false)
{
}

void InterruptableProgressReporter::addSupervisor(ProgressReporter* supervisor) {
    supervisingReporters_.push_back(supervisor);
}

void InterruptableProgressReporter::signalInterruption() {
    interrupted_ = true;
}

void InterruptableProgressReporter::setProgress(float progress) {
    ProgressReporter::setProgress(progress);
    for(auto supervisor : supervisingReporters_) {
        supervisor->setProgress(progress);
    }
    if(interrupted_) {
        interrupted_ = false;
        throw ProcessingInterruptedException();
    }
}
void InterruptableProgressReporter::setProgressRange(const tgt::vec2& progressRange) {
    for(auto supervisor : supervisingReporters_) {
        supervisor->setProgressRange(progressRange);
    }
}
void InterruptableProgressReporter::setProgressMessage(const std::string& message) {
    for(auto supervisor : supervisingReporters_) {
        supervisor->setProgressMessage(message);
    }
}
std::string InterruptableProgressReporter::getProgressMessage() const {
    if(supervisingReporters_.empty()) {
        return "";
    } else {
        return supervisingReporters_[0]->getProgressMessage();
    }
}


LargeDataProcessor::LargeDataProcessor()
    : Processor()
    , continuousUpdate_("continuousUpdate_", "Continuous Update", false, Processor::INVALID_RESULT, Property::LOD_ADVANCED)
    , manualUpdateButton_("manualUpdateButton_", "Start", Processor::INVALID_RESULT, Property::LOD_DEFAULT)
    , stopUpdateButton_("stopUpdateButton", "Stop", Processor::VALID, Property::LOD_DEFAULT)
    , progressDisplay_("progressDisplay", "Progress")
    , progressReporter_()
    , updateForced_(false)
{
    addProperty(continuousUpdate_);
        continuousUpdate_.setGroupID("ldp_processing");
    addProperty(manualUpdateButton_);
        manualUpdateButton_.setGroupID("ldp_processing");
        ON_CHANGE_LAMBDA(manualUpdateButton_, [this] () {
                updateForced_ = true;
                //invalidate();
                });
    addProperty(stopUpdateButton_);
        stopUpdateButton_.setGroupID("ldp_processing");
        stopUpdateButton_.setVisibleFlag(false);
        ON_CHANGE_LAMBDA(stopUpdateButton_, [this] () {
                progressReporter_.signalInterruption();
                });
    addProperty(progressDisplay_);
        progressDisplay_.setGroupID("ldp_processing");
    setPropertyGroupGuiName("ldp_processing", "Processing");

    progressReporter_.addSupervisor(this);
    progressReporter_.addSupervisor(&progressDisplay_);
}

LargeDataProcessor::~LargeDataProcessor() {
}

void LargeDataProcessor::beforeProcess() {
    Processor::beforeProcess();

    for(Property* property : getProperties()) {
        if(property != &manualUpdateButton_
                && property != &stopUpdateButton_
                && property != &progressDisplay_) {
            propertyWritablilityMap_[property] = property->isReadOnlyFlagSet();
            property->setReadOnlyFlag(true);
        }
    }

    manualUpdateButton_.setVisibleFlag(false);
    stopUpdateButton_.setVisibleFlag(true);
}
void LargeDataProcessor::afterProcess() {
    for(Property* property : getProperties()) {
        if(property != &manualUpdateButton_
                && property != &stopUpdateButton_
                && property != &progressDisplay_) {
            tgtAssert(propertyWritablilityMap_.find(property) != propertyWritablilityMap_.end(), "Invalid property. Did you add properties during computation?");
            property->setReadOnlyFlag(propertyWritablilityMap_[property]);
        }
    }

    manualUpdateButton_.setVisibleFlag(true);
    stopUpdateButton_.setVisibleFlag(false);

    updateForced_ = false;

    Processor::afterProcess();
}

bool LargeDataProcessor::isValid() const {
    return !continuousUpdate_.get() && !updateForced_ || Processor::isValid();
}

ProgressReporter& LargeDataProcessor::getProgressReporter() {
    return progressReporter_;
}

} // namespace voreen
