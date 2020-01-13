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

#include "octreequantificationresults.h"

#include "tgt/vector.h"
#include <algorithm>

namespace voreen {

    OctreeQuantificationResults::OctreeQuantificationResults(size_t numChannels, size_t numROIs)
        : numChannels_(numChannels)
        , numROIs_(numROIs)
    {
        // initialize data structures
        initializeData();
    }

    void OctreeQuantificationResults::initializeData() {
        // clear old data and initialize num voxels
        numVoxelsInROIs_.clear();
        intensityHistogramsInROIs_.clear();

        numVoxelsInVolume_ = 0;

        // allocate data structures
        numVoxelsInROIs_.resize(numROIs_, 0);
        
        // first, we have to have the right number of channels
        intensityHistogramsInROIs_.resize(numChannels_);
        for (auto& channel : intensityHistogramsInROIs_) {
            // for each channel, we have to have the right number of ROIs
            channel.resize(numROIs_);
            for (auto& roi : channel) {
                // for each (channel, roi) the intensity histogram has entries 0 .. uint16_max (all initialized with zero)
                roi.resize(static_cast<size_t>(std::numeric_limits<uint16_t>::max()) + 1, 0);
            }
        }
    }

    void OctreeQuantificationResults::clear() {
        initializeData();
    }

    void OctreeQuantificationResults::addResults(OctreeQuantificationResults& results) {
        tgtAssert(numChannels_ == results.numChannels_, "Wrong number of channels");
        tgtAssert(numROIs_ == results.numROIs_, "Wrong number of ROIs");
        tgtAssert(results.numVoxelsInROIs_.size() == numROIs_, "Wrong number of entries in ROI voxel results");
    
        // add the actual results
        numVoxelsInVolume_ += results.numVoxelsInVolume_;
        
        for (size_t r = 0; r < numROIs_; ++r) 
            numVoxelsInROIs_.at(r) += results.numVoxelsInROIs_.at(r);

        for (size_t channel = 0; channel < numChannels_; ++channel) {
            for (size_t roi = 0; roi < numROIs_; ++roi) {
                // use transform to add up all elements
                std::vector<size_t>& myIntensities = intensityHistogramsInROIs_.at(channel).at(roi);
                std::vector<size_t>& theirIntensities = results.intensityHistogramsInROIs_.at(channel).at(roi);
                std::transform(myIntensities.begin(), myIntensities.end(), theirIntensities.begin(), myIntensities.begin(), std::plus<size_t>());

                /*for (size_t intensity = 0; intensity <= static_cast<size_t>(std::numeric_limits<uint16_t>::max()); ++intensity) {       
                    intensityHistogramsInROIs_.at(channel).at(roi).at(intensity) += results.intensityHistogramsInROIs_.at(channel).at(roi).at(intensity);
                }*/
            }
        }
    }
        

    void OctreeQuantificationResults::addVoxel(std::vector<uint16_t> intensities, std::vector<bool> inROIs) {
        tgtAssert(intensities.size() == numChannels_, "Wrong number of channels");
        tgtAssert(inROIs.size() == numROIs_, "Wrong number of ROIs");

        // count voxels
        numVoxelsInVolume_++;
        // count voxels for each ROI
        for (size_t r = 0; r < numROIs_; ++r) {
            if (inROIs.at(r))
                numVoxelsInROIs_.at(r)++;
        }
        // count intensities for each roi and channel separately
        for (size_t channel = 0; channel < numChannels_; ++channel) {
            for (size_t roi = 0; roi < numROIs_; ++roi) {
                if (inROIs.at(roi)) {
                    intensityHistogramsInROIs_.at(channel).at(roi).at(intensities.at(channel))++;
                }
            }
        }
    }

    size_t OctreeQuantificationResults::getNumVoxelsInVolume() {
        return numVoxelsInVolume_;
    }
    
    std::vector<size_t>& OctreeQuantificationResults::getNumVoxelsInROIs() {
        return numVoxelsInROIs_;
    }

    std::vector<std::vector<std::vector<size_t> > >& OctreeQuantificationResults::getIntensityHistogramsInROIs() {
        return intensityHistogramsInROIs_;
    }

    size_t OctreeQuantificationResults::getNumChannels() const {
        return numChannels_;
    }

    size_t OctreeQuantificationResults::getNumROIs() const {
        return numROIs_;
    }

}   // namespace
