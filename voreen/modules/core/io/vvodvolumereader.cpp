/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "vvodvolumereader.h"
#include "vvodformat.h"

#include <fstream>
#include <iostream>

#include "tgt/exception.h"
#include "tgt/vector.h"
#include "tgt/filesystem.h"

#include "voreen/core/utils/stringutils.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/octree/volumeoctree.h"
#include "voreen/core/datastructures/volume/volumepreview.h"
#include "voreen/core/datastructures/volume/volumehash.h"
#include "voreen/core/datastructures/octree/octreebrickpoolmanagerdisk.h"
#include "voreen/core/voreenapplication.h"

#include "voreen/core/utils/hashing.h"

using tgt::vec3;
using tgt::ivec3;
using tgt::hor;
using tgt::lessThanEqual;

namespace voreen {

const std::string VvodVolumeReader::loggerCat_ = "voreen.io.VolumeReader.vvod";

VvodVolumeReader::VvodVolumeReader(ProgressBar* progress)
    : VolumeReader(progress)
{
    extensions_.push_back("vvod");
}

VolumeList* VvodVolumeReader::read(const std::string &url) {
    // create volume url
    VolumeURL origin(url);
    std::string fileName = origin.getPath();

    // try to deserialize octree
    if (!tgt::FileSystem::fileExists(fileName))
        throw tgt::FileException("File does not exist: " + fileName);

    // open file for reading
    std::fstream fileStream(fileName.c_str(), std::ios_base::in);
    if (fileStream.fail())
        throw tgt::FileException("Failed to open cached octree file '" + fileName + "' for reading.");

    // read data stream into deserializer
    XmlDeserializer d(fileName);
    try {
        d.read(fileStream);
    }
    catch (std::exception& e) {
        fileStream.close();
        throw tgt::FileException("Failed to read serialization data stream from cached octree file '" + fileName + "': " + e.what());
    }

    fileStream.close();

    // deserialize octree and meta data from data stream
    VvodStorageObject storageObject;
    try {
        d.deserialize("OctreeVolume", storageObject);
    }
    catch (std::exception& e) {
        throw tgt::FileException("Failed to restore cached octree from file '" + fileName + "': " + e.what());
    }

    if (!storageObject.volumeOctree_)
        throw tgt::FileException("Could not restore octree from file: " + fileName);

    if (!storageObject.volumeOctree_->getRootNode())
        throw tgt::FileException("Deserialized octree has no root node: " + fileName);

    VolumeOctree* octree = storageObject.volumeOctree_;

    // octree has been deserialized, now do some additional stuff and set octree representation to volume

    // assign RAM limit
    size_t ramLimit = VoreenApplication::app()->getCpuRamLimit();
    if (const OctreeBrickPoolManagerDisk* brickPoolManager =
        dynamic_cast<const OctreeBrickPoolManagerDisk*>(static_cast<VolumeOctree*>(octree)->getBrickPoolManager())) {
            const_cast<OctreeBrickPoolManagerDisk*>(brickPoolManager)->setRAMLimit(ramLimit);
    }

    octree->logDescription();

    // min/max values
    std::vector<float> minValues, maxValues, minNormValues, maxNormValues;
    for (size_t i=0; i<octree->getNumChannels(); i++) {
        float minNorm = octree->getRootNode()->getMinValue(i) / 65535.f;
        float maxNorm = octree->getRootNode()->getMaxValue(i) / 65535.f;
        tgtAssert(minNorm <= maxNorm, "invalid min/max values");
        // get realworld mapping
        float min = minNorm;
        float max = maxNorm;
        if (storageObject.metaData_.hasMetaData("RealWorldMapping")) {
            RealWorldMapping rwm = (dynamic_cast<RealWorldMappingMetaData*>(storageObject.metaData_.getMetaData("RealWorldMapping")))->getValue();
            min = rwm.normalizedToRealWorld(minNorm);
            max = rwm.normalizedToRealWorld(maxNorm);
        }
        else
            LWARNING("Found no RealWorldMapping meta data, using normalized values");

        minValues.push_back(min);
        maxValues.push_back(max);
        minNormValues.push_back(minNorm);
        maxNormValues.push_back(maxNorm);
    }
    VolumeMinMax* volumeMinMax = new VolumeMinMax(minValues, maxValues, minNormValues, maxNormValues);

    // histograms
    std::vector<Histogram1D> histograms;
    for (size_t i=0; i<octree->getNumChannels(); i++) {
        histograms.push_back(Histogram1D(*(octree->getHistogram(i))));
    }
    VolumeHistogramIntensity* histogramData = new VolumeHistogramIntensity(histograms);

    // create volume
    VolumeBase* outputVolume = new Volume(octree, &(storageObject.metaData_));
    outputVolume->setOrigin(origin);

    outputVolume->addDerivedData(volumeMinMax);
    outputVolume->addDerivedData(histogramData);
    outputVolume->getDerivedData<VolumePreview>(); //< prevent creation in background thread (brick pool access)

    VolumeHash* volumeHash = new VolumeHash();
    volumeHash->setHash(VoreenHash::getHash(octree->getOctreeConfigurationHash()));
    outputVolume->addDerivedData(volumeHash);

    // create VolumeList with one volume
    VolumeList* vl = new VolumeList();
    vl->add(outputVolume);

    return vl;
}

VolumeReader* VvodVolumeReader::create(ProgressBar* progress) const {
    return new VvodVolumeReader(progress);
}

} // namespace voreen
