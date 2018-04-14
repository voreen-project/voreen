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

#ifndef VRN_LARGE_DATA_PROCESSOR_H
#define VRN_LARGE_DATA_PROCESSOR_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/progressproperty.h"

#include <exception>

namespace voreen {

class ProcessingInterruptedException : public std::exception {
public:
    virtual ~ProcessingInterruptedException() {}
    virtual const char* what() const throw();
};

/*
 * ProgressReporter that can be used
 *  1. as a proxy to multiple others and
 *  2. to interrupt processing
 */
class InterruptableProgressReporter : public ProgressReporter {
public:
    InterruptableProgressReporter();
    virtual ~InterruptableProgressReporter() {}

    /**
     * Add a ProgressReporter that will be proxied by this one
     */
    void addSupervisor(ProgressReporter* supervisor);

    /**
     * Signals the reporter that an interruption has occured.
     * The reporter will throw a ProcessingInterruptedException during the next call of setProgress
     */
    void signalInterruption();

    virtual void setProgress(float progress);
    virtual void setProgressRange(const tgt::vec2& progressRange);
    virtual void setProgressMessage(const std::string& message);
    virtual std::string getProgressMessage() const;

protected:
    bool interrupted_;
    std::vector<ProgressReporter*> supervisingReporters_;
};


/**
 * Processor that can be subclassed for Processors that do expensive
 * computation in their process() method.
 * Subclassing processors should report the progress of their computation
 * using the ProgressReporter returned by getProgressReporter().
 * The LargeDataProcessor will provide an interface that shows the progress
 * in the GUI and means to start/stop computation.
 *
 * Attention: The setProgress() method can throw a ProcessingInterruptedException
 * which has to be caught in process() at the latest.
 */
class LargeDataProcessor : public Processor {
public:
    LargeDataProcessor();
    virtual ~LargeDataProcessor();

    virtual void beforeProcess();
    virtual void afterProcess();
    virtual bool isValid() const;

    virtual bool usesExpensiveComputation() const {
        return true;
    }

protected:
    ProgressReporter& getProgressReporter();

private:
    BoolProperty continuousUpdate_;
    ButtonProperty manualUpdateButton_;
    ButtonProperty stopUpdateButton_;
    ProgressProperty progressDisplay_;

    /// Used to save/restore property writability before/after expensive computation
    std::map<Property*, bool> propertyWritablilityMap_;

    InterruptableProgressReporter progressReporter_;

    bool updateForced_;
};

} // namespace voreen

#endif // VRN_LARGE_DATA_PROCESSOR_H
