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

#include "timeserieslistsource.h"

namespace voreen {

const std::string TimeSeriesListSource::loggerCat_("voreen.scivis21.TimeSeriesListSource");

TimeSeriesListSource::TimeSeriesListSource()
    : Processor()
    , outport_(Port::OUTPORT, "outport", "Timeseries List Output", false)
    , filenameProp_("filenameprop", "Load Time Series List File from", "Select file...", VoreenApplication::app()->getUserDataPath(), "Time Series List files (*.tsl)", FileDialogProperty::OPEN_FILE, Processor::INVALID_PATH)
    , loadButton_("loadButton", "Load", INVALID_PATH)
    , loadTimeSeriesList_(true)
{
    addPort(outport_);
    addProperty(filenameProp_);
    addProperty(loadButton_);
}

void TimeSeriesListSource::invalidate(int inv) {
    Processor::invalidate(inv);

    if (inv == Processor::INVALID_PATH && isInitialized()) {
        loadTimeSeriesList_ = true;
    }
}

void TimeSeriesListSource::process() {
    if (loadTimeSeriesList_){
        loadTimeSeriesList();
        loadTimeSeriesList_ = false;
    }
}

void TimeSeriesListSource::loadTimeSeriesList() {
    if (!isInitialized())
        return;

    outport_.setData(nullptr);

    if (filenameProp_.get().empty()) {
        LWARNING("no filename specified");
        return;
    }

    try {
        std::unique_ptr<TimeSeriesList> timeseries(new TimeSeriesList());

        std::ifstream stream(filenameProp_.get());
        XmlDeserializer json;
        json.read(stream);
        Deserializer s(json);
        s.deserialize("TimeSeriesList", *timeseries);
        outport_.setData(timeseries.release(), true);
        LINFO(filenameProp_.get() << " loaded sucessfully!");
    } catch(tgt::Exception& e) {
        LERROR(e.what());
        filenameProp_.set("");
    }
}

}