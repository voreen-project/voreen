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

#include "timeserieslistsave.h"

#include <fstream>

namespace voreen {

    const std::string TimeSeriesListSave::loggerCat_("voreen.core.TimeSeriesListSave");

    TimeSeriesListSave::TimeSeriesListSave()
        : Processor()
        , inport_(Port::INPORT, "inport", "Time Series List Input")
        , fileProp_("file", "Time Series List File", "Select Time Series List file...", "./",
            "Time Series List files (*.tsl)", FileDialogProperty::SAVE_FILE, Processor::INVALID_PATH)
        , saveButton_("save", "Save")
        , continousSave_("continousSave", "Save continuously", false)
    {
        addPort(inport_);

        saveButton_.onChange(MemberFunctionCallback<TimeSeriesListSave>(this, &TimeSeriesListSave::saveFile));
        addProperty(fileProp_);
        addProperty(saveButton_);
        addProperty(continousSave_);
    }

    Processor* TimeSeriesListSave::create() const {
        return new TimeSeriesListSave();
    }

    void TimeSeriesListSave::invalidate(int inv) {
        Processor::invalidate(inv);
        //auto save on path change
        if (inv == Processor::INVALID_PATH && isInitialized())
            saveFile();
    }

    void TimeSeriesListSave::process() {
        if (inport_.hasChanged() && continousSave_.get())
            saveFile();
    }

    void TimeSeriesListSave::saveFile() {
        const TimeSeriesList* timeSeriesList = inport_.getData();
        if (!timeSeriesList)
            return;

        std::string filename = fileProp_.get();
        if (filename.empty()) {
            LWARNING("Could not save time series list: filename is empty");
            return;
        }

        LINFO("Saving Time Series List to file: " << filename);

        std::fstream stream(filename.c_str(), std::ios::out);
        if (stream.fail()) {
            LERROR("Failed to open file for writing: " << filename);
        }
        else {
            XmlSerializer s;
            try {
                s.serialize("TimeSeriesList", timeSeriesList);
            }
            catch (VoreenException & e) {
                LERRORC("voreen.core.TimeSeriesListSave", "Failed to serialize time series: " + std::string(e.what()));
                return;
            }
            s.write(stream);
            stream.close();
        }
    }

    void TimeSeriesListSave::adjustPropertiesToInput() {
        fileProp_.setFileFilter("Time Series List files (*.tsl)");
    }

} // namespace voreen
