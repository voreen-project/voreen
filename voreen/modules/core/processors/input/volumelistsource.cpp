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

#include "volumelistsource.h"

#include "voreen/core/datastructures/volume/volumelist.h"
#include "voreen/core/io/volumeserializerpopulator.h"
#include "voreen/core/io/volumeserializer.h"

namespace voreen {

const std::string VolumeListSource::loggerCat_("voreen.core.VolumeListSource");

VolumeListSource::VolumeListSource()
    : Processor(),
      outport_(Port::OUTPORT, "volumecollection", "VolumeList Output", false),
      volumeURLList_("volumeURLList", "Volume URL List", std::vector<std::string>())
{
    addPort(outport_);
    addProperty(volumeURLList_);
}

VolumeListSource::~VolumeListSource() {
}

Processor* VolumeListSource::create() const {
    return new VolumeListSource();
}

void VolumeListSource::process() {
    // nothing
}

void VolumeListSource::initialize() {
    Processor::initialize();

    volumeURLList_.loadVolumes(false, true);

    outport_.setData(volumeURLList_.getVolumes(true), true);
}

void VolumeListSource::invalidate(int inv) {
    outport_.setData(volumeURLList_.getVolumes(true), true);

    Processor::invalidate(inv);
}

void VolumeListSource::setVolumeList(VolumeList* collection, bool owner, bool selected) {
    volumeURLList_.clear();
    if (!collection)
        return;
    for (size_t i=0; i<collection->size(); i++)
        volumeURLList_.addVolume(collection->at(i), owner, selected);

    invalidate();
}

VolumeList* VolumeListSource::getVolumeList() const {
    return volumeURLList_.getVolumes(false);
}

void VolumeListSource::loadVolumes(const std::string& url, bool selected /*= true*/, bool clear /*= false*/) {
    tgtAssert(!url.empty(), "passed url is empty");

    // load volumes
    VolumeSerializerPopulator serializerPopulator;
    VolumeList* volumeList = 0;
    try {
        volumeList = serializerPopulator.getVolumeSerializer()->read(url);
    }
    catch (tgt::FileException& e) {
        LWARNING(e.what());
        throw e;
    }
    catch (std::bad_alloc& e) {
        LWARNING("Bad allocation while loading volumes from: " << url);
        throw e;
    }
    catch (tgt::Exception& e) {
        LWARNING("Unknown exception while loading volumes from '" << url << "':" << e.what());
        throw tgt::FileException(e.what(), url);
    }
    tgtAssert(volumeList && !volumeList->empty(), "no volumes loaded");

    // Clear old volumes on demand.
    if(clear)
        volumeURLList_.clear();

    // add loaded volumes
    for (size_t i=0; i<volumeList->size(); i++) {
        tgtAssert(volumeList->at(i), "null pointer as volume");
        volumeURLList_.addVolume(volumeList->at(i), true, selected);
    }
    delete volumeList;

    invalidate();
}

} // namespace
