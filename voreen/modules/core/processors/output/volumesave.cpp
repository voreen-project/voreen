/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "volumesave.h"
#include "voreen/core/io/volumeserializer.h"
#include "voreen/core/io/volumeserializerpopulator.h"
#include "voreen/core/io/volumewriter.h"
#include "voreen/core/utils/stringutils.h"
#include "voreen/core/io/progressbar.h"

namespace voreen {

const std::string VolumeSave::loggerCat_("voreen.core.VolumeSave");

VolumeSave::VolumeSave()
    : VolumeProcessor()
    , inport_(Port::INPORT, "volumehandle.input", "Volume Input", false)
    , filename_("outputFilename", "File", "Select file...", "",
            "Voreen Volume Data (*.vvd)", FileDialogProperty::SAVE_FILE, Processor::INVALID_PATH)
    , saveButton_("save", "Save")
    , saveOnPathChange_("saveOnPathChange","Save on path change",true)
    , continousSave_("continousSave", "Save on inport change", false)
    , volumeInfo_("volumeInfo","info")
    , saveVolume_(false)
    , volSerializerPopulator_(0)
{
    addPort(inport_);

    saveButton_.onChange(MemberFunctionCallback<VolumeSave>(this, &VolumeSave::saveVolume));

    addProperty(filename_);
    addProperty(saveButton_);
    addProperty(saveOnPathChange_);
    addProperty(continousSave_);
    addProperty(volumeInfo_);
}

VolumeSave::~VolumeSave() {
}

Processor* VolumeSave::create() const {
    return new VolumeSave();
}

void VolumeSave::initialize() {
    VolumeProcessor::initialize();

    tgtAssert(!volSerializerPopulator_, "VolumeSerializerPopulator already created");
    volSerializerPopulator_ = new VolumeSerializerPopulator();
    std::vector<std::string> filters = constructFilters();
    filename_.setFileFilter(strJoin(filters, ";;"));
}

void VolumeSave::deinitialize() {
    delete volSerializerPopulator_;
    volSerializerPopulator_ = 0;
    volumeInfo_.setVolume(0);
    VolumeProcessor::deinitialize();
}

void VolumeSave::invalidate(int inv) {
    Processor::invalidate(inv);

    if (inport_.hasChanged()) {
        volumeInfo_.setVolume(inport_.getData());
        if (continousSave_.get())
            saveVolume_ = true;
    }

    if(inv == Processor::INVALID_PATH && saveOnPathChange_.get() && isInitialized()) {
        saveVolume_ = true;
    }
}

void VolumeSave::process() {
    if (saveVolume_){
        saveVolume();
        saveVolume_ = false;
    }
}

void VolumeSave::saveVolume() {
    if (!isInitialized())
        return;
    if (!inport_.getData()) {
        LWARNING("no input volume");
        return;
    }
    if (filename_.get() == "") {
        LWARNING("no filename specified");
        return;
    }

    tgtAssert(volSerializerPopulator_, "VolumeSerializerPopulator not instantiated");
    tgtAssert(volSerializerPopulator_->getVolumeSerializer(), "no VolumeSerializer");
    VolumeSerializer* s = volSerializerPopulator_->getVolumeSerializer();
    //set the progress bar, if possible
    if (progressBars_.size() >= 1)
        s->setProgressBar(dynamic_cast<ProgressBar*>(progressBars_.at(0)));
    else
        s->setProgressBar(0);
    try {
        s->write(filename_.get(), inport_.getData());
    }
    catch(tgt::FileException& e) {
        LERROR(e.what());
        filename_.set("");
        s->setProgressBar(0);
    }
}

std::vector<std::string> VolumeSave::constructFilters() const {
    tgtAssert(volSerializerPopulator_, "VolumeSerializerPopualator not instantiated");
    std::vector<std::string> filters;

    const std::vector<VolumeWriter*> volumeWriters = volSerializerPopulator_->getVolumeSerializer()->getWriters();
    std::string vvdWriterFilter;
    for (size_t i=0; i<volumeWriters.size(); i++) {
        VolumeWriter* curWriter = volumeWriters.at(i);

        // extensions
        std::vector<std::string> extensionVec = curWriter->getSupportedExtensions();
        for (size_t j=0; j<extensionVec.size(); j++) {
            std::string extension = extensionVec.at(j);
            std::string filterStr = curWriter->getFormatDescription() + " (*." + extension + ")";
            if (extension == "vvd")
                vvdWriterFilter = filterStr;
            else
                filters.push_back(filterStr);
        }

        // filename
        std::vector<std::string> filenamesVec = curWriter->getSupportedFilenames();
        for (size_t j=0; j<filenamesVec.size(); j++) {
            std::string filename = filenamesVec.at(j);
            std::string filterStr = curWriter->getFormatDescription() + " (" + filename + ")";
            filters.push_back(filterStr);
        }

    }

    if (vvdWriterFilter != "")
        filters.insert(filters.begin(), vvdWriterFilter);

    return filters;
}

}   // namespace
