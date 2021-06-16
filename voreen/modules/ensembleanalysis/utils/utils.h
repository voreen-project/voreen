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

#ifndef VRN_ENSEMBLE_UTILS_H
#define VRN_ENSEMBLE_UTILS_H

#include "voreen/core/voreencoreapi.h"
#include "voreen/core/io/volumereader.h"
#include "voreen/core/io/volumeserializerpopulator.h"

#include "tgt/vector.h"

namespace voreen {

class VolumeBase;
class VolumeURL;


/**
 * Simple helper class for loading volumes as part of an ensemble.
 * This class can be used to incorporate custom readers, if desired.
 * E.g. multi-channel volumes stored in HDF5 files will not be split into multiple volumes using this function.
 * This class can be used as other volume readers, however, this is currently not the intended way.
 * Hence, the reader is not registered by the module.
 */
class EnsembleVolumeReader : public VolumeReader {
public:

    /**
     * Constructor.
     *
     * @param progressBar Optional progress bar to assign to the volume readers and writers.
     *          Is <emph>not</emph> deleted by the populator's destructor.
     */
    EnsembleVolumeReader(ProgressBar* progressBar = nullptr);
    VolumeReader* create(ProgressBar* progressBar = nullptr) const;

    virtual std::string getClassName() const { return "EnsembleVolumeReader"; }
    virtual std::string getFormatDescription() const { return "Voreen ensemble volume reader"; }

    /**
     * Returns true if a suitable reader is available for the given path, false otherwise.
     */
    bool canRead(const std::string& path);

    //! @see VolumeReader
    VolumeList* read(const std::string& url);
    std::vector<VolumeURL> listVolumes(const std::string& url) const;

    /**
     * Reads the volume at the given URL.
     * In contrast to the VolumeReader specification, this implementation is more forgiving and
     * will return a nullptr instead of throwing an exception if the volume could not be loaded.
     */
    VolumeBase* read(const VolumeURL& origin);

    /**
     * Tries to find an appropriate reader for the given path.
     * @param path path to the volume file
     * @return a volume reader, if one was found or nullptr otherwise
     */
    VolumeReader* getVolumeReader(const std::string& path) const;

private:

    VolumeSerializerPopulator volumeSerializerPopulator_;
};

/**
 * Utility function mapping a value within range A to the equivalent value in range B.
 */
template<typename T, typename S>
S mapRange(const T& valA, const T& minA, const T& maxA, const S& minB, const S& maxB) {
    //tgtAssert(valA >= minA && valA <= maxA, "value out of range"); // value may lay outside intentionally!
    // Cast into receiver type.
    return S(minB + (maxB - minB) * (valA - minA) / (maxA - minA));
}

template<typename T, typename S>
S mapRange(const T& valA, const tgt::Vector2<T>& rangeA, const tgt::Vector2<S>& rangeB) {
    return S(rangeB.x + (rangeB.y - rangeB.x) * (valA - rangeA.x) / (rangeA.y - rangeA.x));
}

}

#endif
