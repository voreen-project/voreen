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

#include "metadataadder.h"
#include <fstream>
namespace voreen {

MetaDataAdder::MetaDataAdder()
    : Processor()
    , inport_(Port::INPORT, "inport", "Volume List Inport")
    , outport_(Port::OUTPORT, "outport", "Volume List Outport")
    , addTime_("addTime", "Add Time?")
    , timeInformationFile_("timeInformationFile", "Time Information: ", "Select file with timestep information", "")
    , modalityString_("nameString", "Modality")
{
    addPort(inport_);
    addPort(outport_);
    addProperty(addTime_);
    ON_CHANGE_LAMBDA(addTime_, [this] {
       timeInformationFile_.setVisibleFlag(addTime_.get());
    });
    addProperty(timeInformationFile_);
    timeInformationFile_.setVisibleFlag(false);
    addProperty(modalityString_);
}

Processor* MetaDataAdder::create() const {
    return new MetaDataAdder();
}

void MetaDataAdder::setDescriptions() {
    setDescription("This processor can add explicit time steps to the data as well as " \
                   "a name. The time information should be given in a file where each line " \
                   "contains one number which is used as a timestep.");
}

void MetaDataAdder::process() {
    clearOutput();

    VolumeList* list = new VolumeList();
    std::vector<float>  timesteps;

    // Read time data
    if(addTime_.get() && !timeInformationFile_.get().empty()) {
        std::string filepath = timeInformationFile_.get();
        std::ifstream infile(filepath);
        float time;
        while (infile >> time) {
            timesteps.push_back(time);
        }
        infile.close();
    }

    // Add data to volumes
    for(size_t i = 0; i<inport_.getData()->size(); i++) {
        VolumeBase* volumeDec = new VolumeDecoratorIdentity(inport_.getData()->at(i));
        decorators_.push_back(std::unique_ptr<VolumeBase>(volumeDec));
        if(timesteps.size() > i) {
            volumeDec = new VolumeDecoratorReplaceTimestep(volumeDec, timesteps[i]);
            decorators_.push_back(std::unique_ptr<VolumeBase>(volumeDec));
        }
        if(!modalityString_.get().empty()) {
            volumeDec = new VolumeDecoratorReplace(volumeDec, VolumeBase::META_DATA_NAME_MODALITY, new StringMetaData(modalityString_.get()), true);
            decorators_.push_back(std::unique_ptr<VolumeBase>(volumeDec));
        }
        list->add(volumeDec);
    }
    outport_.setData(list);
}

void MetaDataAdder::clearOutput() {
    outport_.clear();
    decorators_.clear();
}

} //namespace
