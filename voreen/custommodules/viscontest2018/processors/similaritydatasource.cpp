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

#include "similaritydatasource.h"

#include "voreen/core/datastructures/volume/volumelist.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/io/volumereader.h"
#include "voreen/core/io/volumeserializer.h"

#include "../datastructures/similaritydata.h"
#include "../utils/colorpool.h"

namespace voreen {

const std::string SimilarityDataSource::SCALAR_FIELD_NAME = "Scalar";
const std::string SimilarityDataSource::SIMULATED_TIME_NAME = "simulated_time";
const std::string SimilarityDataSource::loggerCat_("voreen.viscontest2018.SimilarityDataSource");

SimilarityDataSource::SimilarityDataSource()
    : Processor()
    , similarityPath_("ensemblepath", "Ensemble Path", "Select Ensemble root folder", "", "", FileDialogProperty::DIRECTORY)
    , loadSimilarityDataButton_("loadDataset", "Load Similarity Data")
    , outport_(Port::OUTPORT, "similaritydata", "SimilarityData Output", false)
    , loadedChannels_("loadedChannels", "Loaded Channels")
    , selectedTimeStep_("selectedTimeStep", "Selected time step", 0)
{
    addPort(outport_);
    addProperty(similarityPath_);
    addProperty(loadSimilarityDataButton_);
    addProperty(loadedChannels_);
    addProperty(selectedTimeStep_);

    loadSimilarityDataButton_.onChange(MemberFunctionCallback<SimilarityDataSource>(this, &SimilarityDataSource::buildSimilarityData));
}

SimilarityDataSource::~SimilarityDataSource() {
    clearSimilarityData();
}

Processor* SimilarityDataSource::create() const {
    return new SimilarityDataSource();
}

void SimilarityDataSource::process() {
    // Nothing to do.
}

void SimilarityDataSource::initialize() {
    Processor::initialize();
}

void SimilarityDataSource::clearSimilarityData() {
    // TODO clear it
    loadedChannels_.reset();
}

void SimilarityDataSource::buildSimilarityData() {

    // Delete old data.
    clearSimilarityData();

    if(similarityPath_.get().empty())
        return;

    loadSimilarityDataButton_.setReadOnlyFlag(true);

    SimilarityData* similarityData = new SimilarityData();

    std::vector<std::string> channels = tgt::FileSystem::listSubDirectories(similarityPath_.get(), true);

    int numTimeSteps = 0;

    for(const std::string& channelName : channels) {
        std::string channelPath = similarityPath_.get() + "/" + channelName;
        std::vector<std::string> timeStepNames = tgt::FileSystem::readDirectory(channelPath, true, false);
        loadedChannels_.addRow(channelName);
        numTimeSteps = std::max(numTimeSteps, static_cast<int>(timeStepNames.size()));
    }

    selectedTimeStep_.setMinValue(0);
    selectedTimeStep_.setMaxValue(numTimeSteps);
    selectedTimeStep_.set(0);

    outport_.setData(similarityData, true);

    loadSimilarityDataButton_.setReadOnlyFlag(false);
}

} // namespace
