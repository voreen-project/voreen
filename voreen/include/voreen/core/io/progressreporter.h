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

#ifndef VRN_PROGRESSREPORTER_H
#define VRN_PROGRESSREPORTER_H

#include "voreen/core/voreencoreapi.h"
#include "tgt/vector.h"
#include "tgt/assert.h"
#include <string>
#include <array>
#include <memory>
#include <boost/thread.hpp>
#include <numeric>

namespace voreen {

class SubtaskProgressReporter;

/**
 * Abstract base class for classes that report a task progress to the user.
 */
class VRN_CORE_API ProgressReporter {
public:

    ProgressReporter();
    virtual ~ProgressReporter() {}

    /**
     * Assigns the current progress state.
     *
     * @param progress The current progress. Must lie in [0, 1]
     */
    virtual void setProgress(float progress);

    /// Returns the current progress value.
    virtual float getProgress() const;

    /**
     * Assigns the range into which the progress value will be transformed,
     * i.e., actualProgress = progressRange.x + progress*(progressRange.y-progressRange.x)
     * The progress range must be a subrange of [0.f;1.f].
     * The default range is [0.f;1.f].
     */
    virtual void setProgressRange(const tgt::vec2& progressRange);

    /// Returns the current progress range.
    virtual tgt::vec2 getProgressRange() const;

    /**
     * Assigns a message that is to be displayed, e.g. in a progress dialog.
     */
    virtual void setProgressMessage(const std::string& message) = 0;

    /// Returns the current progress message;
    virtual std::string getProgressMessage() const = 0;

    /// update function (e.g., for GUI), is called on progress, e.g., by setProgress(). default implementation is empty
    virtual void update() { }

protected:

    float progress_;
    tgt::vec2 progressRange_;

    bool printedErrorMessage_;

    static std::string loggerCat_;
};

/**
 * A class that maps its progress to a subrange of the parent reporter.
 * This can be seen as a more object oriented approach to the setProgressRange method of ProgressReporter
 * (so you should probably use this instead).
 */
class VRN_CORE_API SubtaskProgressReporter : public ProgressReporter {
public:
    SubtaskProgressReporter(ProgressReporter& parent, const tgt::vec2& subtaskProgressRange);
    virtual ~SubtaskProgressReporter() {}

    /**
     * Assigns the current progress state.
     *
     * @param progress The current progress. Must lie in [0, 1]
     */
    virtual void setProgress(float progress);

    /// Returns the current progress value.
    virtual float getProgress() const;

    /**
     * Assigns a message that is to be displayed, e.g. in a progress dialog.
     * Passes the message through to the parent
     */
    virtual void setProgressMessage(const std::string& message);

    /**
     * Returns the current progress message (Taken from parent)
     */
    virtual std::string getProgressMessage() const;

protected:
    ProgressReporter& parent_;
    tgt::vec2 subtaskProgressRange_;
};

/**
 * A collection of progress reporters with a static size n.
 * The range of progress is devided equally amoung the n subreporters.
 */
template<uint8_t n>
struct SubtaskProgressReporterCollection {
    SubtaskProgressReporterCollection(ProgressReporter& parent)
        : reporters_()
    {
        for(uint8_t i=0; i < n; ++i) {
            reporters_[i].reset(new SubtaskProgressReporter(parent, tgt::vec2(static_cast<float>(i)/n, static_cast<float>(i+1)/n)));
        }
    }

    SubtaskProgressReporterCollection(ProgressReporter& parent, const std::array<float, n>& subtaskWorkWeight)
        : reporters_()
    {
        float weightSum = std::accumulate(subtaskWorkWeight.begin(), subtaskWorkWeight.end(), 0.0f);
        float progressBegin = 0.0f;
        for(uint8_t i=0; i < n; ++i) {
            float progressEnd = progressBegin + subtaskWorkWeight[i]/weightSum;
            reporters_[i].reset(new SubtaskProgressReporter(parent, tgt::vec2(progressBegin, progressEnd)));
            progressBegin = progressEnd;
        }
        tgtAssert(progressBegin == 1.0f, "Invalid progress end");
    }

    template<uint8_t i>
    SubtaskProgressReporter& get() {
        static_assert(i < n, "Invalid SubtaskProgressReporter index.");
        return *reporters_[i];
    }

    std::array<std::unique_ptr<SubtaskProgressReporter>, n> reporters_;
};

/**
 * An abstraction over ProgressReporter that offers less functionality (and thus is not
 * directly derived from it), but can safely be used from multiple threads at a time.
 *
 * Provide the total number of computation steps that need to be done in the constructor and
 * call reportStepDone from the worker threads once one of the step is finished.
 *
 * Do not call reportStepDone more often than totalNumberOfSteps!
 */
struct ThreadedTaskProgressReporter {
    ThreadedTaskProgressReporter(ProgressReporter& parent, size_t totalNumberOfSteps);

    // Returns true if the execution was interrupted, and further computations should be canceled.
    bool reportStepDone();

    ProgressReporter& parent_;
    boost::mutex mutex_;
    size_t totalNumberOfSteps_;
    size_t step_;
};


} // namespace voreen

#endif // VRN_PROGRESSREPORTER_H
