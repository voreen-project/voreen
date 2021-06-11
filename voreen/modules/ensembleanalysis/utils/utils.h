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
#include "voreen/core/io/volumeserializerpopulator.h"

#include "tgt/vector.h"

namespace voreen {

class VolumeBase;

/**
 * Simple helper class for loading volumes as part of an ensemble.
 * Basically wraps a VolumeSerializerPopulator and allows to incorporate custom readers.
 */
class EnsembleVolumeReaderPopulator {
public:

    /**
     * Constructor.
     *
     * @param progressBar Optional progress bar to assign to the volume readers and writers.
     *          Is <emph>not</emph> deleted by the populator's destructor.
     */
    EnsembleVolumeReaderPopulator(ProgressBar* progressBar = nullptr);

    /**
     * Returns a volume reader for the given path or nullptr, if no suitable reader was found.
     * This function can be used to incorporate custom readers, if desired.
     * E.g. multi-channel volumes stored in HDF5 files will not be split into multiple volumes using this function.
     * @note EnsembleVolumeReaderPopulator own the returned reader!
     * @see TimeStep
     * @see EnsembleDataSource
     */
    VolumeReader* getVolumeReader(const std::string& path) const;

private:
    VolumeSerializerPopulator volumeSerializerPopulator_;
};

/**
 * If a volume reader (or file format) does not support a disk representation, a Swap disk can be added.
 * If a RAM representation is requested on a volume with a swap representation, it will be loaded from disk
 * on demand, which typically takes longer than just loading a disk representation, but serves as
 * (experimental) workaround for missing disk representations.
 */
class VolumeRAMSwap {
public:
    static bool tryAddVolumeSwap(VolumeBase* volumeBase);
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
