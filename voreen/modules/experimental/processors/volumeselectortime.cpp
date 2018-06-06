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

#include "volumeselectortime.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumelist.h"
#include "voreen/core/processors/processorwidgetfactory.h"

namespace voreen {

const std::string VolumeSelectorTime::loggerCat_("voreen.VolumeSelectorTime");

VolumeSelectorTime::VolumeSelectorTime()
    : Processor(),
      volumeID_("volumeID", "Selected volume", 0, -99999999, std::numeric_limits<int>::max() ),
      inport_(Port::INPORT, "volumecollection", "Volume List", false),
      outport_(Port::OUTPORT, "volumehandle.volumehandle", "volumehandle.volumehandle", false)
{
    addPort(inport_);
    addPort(outport_);

    addProperty(volumeID_);
}

Processor* VolumeSelectorTime::create() const {
    return new VolumeSelectorTime();
}

void VolumeSelectorTime::process() {
    // nothing
}

void VolumeSelectorTime::initialize() {
    Processor::initialize();

    adjustToVolumeList();
}

void VolumeSelectorTime::invalidate(int /*inv = INVALID_RESULT*/) {
    adjustToVolumeList();
}

void VolumeSelectorTime::updateOutport(int t) {
    const VolumeList* collection = inport_.getData();
    if (collection && !collection->empty()) {
        size_t next = 0;

        for(size_t i=1; i<collection->size(); i++) {

            if(abs(t - collection->at(i)->getMetaDataValue<IntMetaData>("FrameTime", 0)) < abs(t - collection->at(next)->getMetaDataValue<IntMetaData>("FrameTime", 0)))
                next = i;
        }

        VolumeBase* vh = collection->at(next);
        if(vh != outport_.getData())
            outport_.setData(vh, false);
    }
}

void VolumeSelectorTime::adjustToVolumeList() {
    if (!outport_.isInitialized())
        return;

    const VolumeList* collection = inport_.getData();

    if (collection && !collection->empty()) {
        int minTime = collection->at(0)->getMetaDataValue<IntMetaData>("FrameTime", 0);
        int maxTime = collection->at(0)->getMetaDataValue<IntMetaData>("FrameTime", 0);

        for(size_t i=1; i<collection->size(); i++) {
            int t = collection->at(i)->getMetaDataValue<IntMetaData>("FrameTime", 0);
            if(minTime > t)
                minTime = t;
            else if(maxTime < t)
                maxTime = t;
        }

        volumeID_.setMinValue(minTime);
        volumeID_.setMaxValue(maxTime);
    }
    else {
        volumeID_.setMinValue(-999999999);
        volumeID_.setMaxValue(std::numeric_limits<int>::max());
    }
    volumeID_.updateWidgets();

    updateOutport(volumeID_.get());
}

} // namespace
