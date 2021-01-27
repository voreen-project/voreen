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

#ifndef VRN_OCTREEQUANTIFICATIONRESULTS_H
#define VRN_OCTREEQUANTIFICATIONRESULTS_H

#include <vector>
#include "voreen/core/datastructures/octree/octreeutils.h" 

namespace voreen {

class OctreeQuantificationResults {

public:

    OctreeQuantificationResults(size_t numChannels, size_t numROIs);

    // set everything to zero
    void clear();

    // add results (e.g., from a local thread)
    void addResults(OctreeQuantificationResults& results);

    /**
     * Adds a voxel to the result. 
     * @param intensities the raw intensity values of all channels (must contain numChannels entries)
     * @param inROIs boolean value for each of the ROIS if the voxel lies within that ROI (must contain numROIs entries)
     */
    void addVoxel(std::vector<uint16_t> intensities, std::vector<bool> inROIs);

    // getter for results
    size_t getNumVoxelsInVolume();
    std::vector<size_t>& getNumVoxelsInROIs();
    std::vector<std::vector<std::vector<size_t> > >& getIntensityHistogramsInROIs();

    size_t getNumChannels() const;
    size_t getNumROIs() const;

private:

    OctreeQuantificationResults();

    void initializeData();

    size_t numVoxelsInVolume_;    
    std::vector<size_t> numVoxelsInROIs_;   // one entry for each ROI (if a ROI in between is not present: zero)
    std::vector<std::vector<std::vector<size_t> > > intensityHistogramsInROIs_;   // inner vector contains an intensity histogram (one entry for each raw uint16 value), middle vector one for each ROI, outer vector one for each channel

    size_t numChannels_;
    size_t numROIs_;
};

} // namespace voreen

#endif
