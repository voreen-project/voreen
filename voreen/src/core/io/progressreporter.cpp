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

#include "voreen/core/io/progressreporter.h"

#include "tgt/assert.h"
#include "tgt/logmanager.h"

namespace voreen {

std::string ProgressReporter::loggerCat_ = "voreen.ProgressReporter";

ProgressReporter::ProgressReporter()
    : progress_(0)
    , progressRange_(0.f, 1.f)
    , printedErrorMessage_(false)
{ }

void ProgressReporter::setProgress(float progress) {
    tgtAssert(progressRange_.x <= progressRange_.y && progressRange_.x >= 0.f && progressRange_.y <= 1.f, "invalid progress range");

    if (progress < 0.f) {
        if (!printedErrorMessage_)
            LWARNING("progress value " << progress << " out of valid range [0,1]");
        printedErrorMessage_ = true;
        progress = 0.f;
    }
    else if (progress > 1.f) {
        if (!printedErrorMessage_)
            LWARNING("progress value " << progress << " out of valid range [0,1]");
        printedErrorMessage_ = true;
        progress = 1.f;
    }
    else {
        printedErrorMessage_ = false;
    }

    progress_ = progressRange_.x + progress*(progressRange_.y-progressRange_.x);

    update();
}

float ProgressReporter::getProgress() const {
    return progress_;
}

void ProgressReporter::setProgressRange(const tgt::vec2& progressRange) {
    tgtAssert(progressRange.x <= progressRange.y && progressRange.x >= 0.f && progressRange.y <= 1.f, "invalid progress range");
    progressRange_ = progressRange;
}

tgt::vec2 ProgressReporter::getProgressRange() const {
    return progressRange_;
}

SubtaskProgressReporter::SubtaskProgressReporter(ProgressReporter& parent, const tgt::vec2& subtaskProgressRange)
    : parent_(parent)
    , subtaskProgressRange_(subtaskProgressRange)
{
    tgtAssert(subtaskProgressRange_.x < subtaskProgressRange_.y, "Invalid subtask progress range");
}

void SubtaskProgressReporter::setProgress(float progress) {
    float progressAccordingToLocalRange = progressRange_.x + progress*(progressRange_.y-progressRange_.x);

    float parentProgress = subtaskProgressRange_.x + progressAccordingToLocalRange*(subtaskProgressRange_.y-subtaskProgressRange_.x);
    parent_.setProgress(parentProgress);
}
float SubtaskProgressReporter::getProgress() const {
    float parentProgress = parent_.getProgress();
    return (parentProgress - subtaskProgressRange_.x) / (subtaskProgressRange_.y-subtaskProgressRange_.x);
}

void SubtaskProgressReporter::setProgressMessage(const std::string& message) {
    parent_.setProgressMessage(message);
}

std::string SubtaskProgressReporter::getProgressMessage() const {
    return parent_.getProgressMessage();
}

ThreadedTaskProgressReporter::ThreadedTaskProgressReporter(ProgressReporter& parent, size_t totalNumberOfSteps)
    : parent_(parent)
    , mutex_()
    , totalNumberOfSteps_(totalNumberOfSteps)
    , step_(0)
{
    parent_.setProgress(0);
}

bool ThreadedTaskProgressReporter::reportStepDone() {
    boost::lock_guard<boost::mutex> lock(mutex_);

    ++step_;
    tgtAssert(step_ <= totalNumberOfSteps_, "Invalid step");

    try {
        parent_.setProgress(static_cast<float>(step_)/totalNumberOfSteps_);
    } catch(boost::thread_interrupted& e) {
        return true;
    }

    return false;
}

}
