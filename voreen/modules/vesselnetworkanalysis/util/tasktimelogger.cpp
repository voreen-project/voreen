/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "tasktimelogger.h"
#include "voreen/core/utils/stringutils.h"

namespace voreen {

TaskTimeLogger::TaskTimeLogger(std::string taskName, tgt::LogLevel loglevel)
    : taskName_(taskName)
    , loglevel_(loglevel)
    , startTime_(TaskTimeLogger::Clock::now())
{
}

TaskTimeLogger::~TaskTimeLogger() {
    TimePoint now = Clock::now();
    auto totaltime = now - startTime_;
    auto totaltime_ms = std::chrono::duration_cast<std::chrono::milliseconds>(totaltime);

    std::ostringstream msg;
    msg << taskName_;
    msg << ": ";
    msg << formatTime(totaltime_ms.count());
    LogMgr.log("TaskTimeLogger", loglevel_, msg.str());
}

};
