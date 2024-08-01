/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#ifndef VRN_QUADHIDACVOLUMEREADER_H
#define VRN_QUADHIDACVOLUMEREADER_H

#include "voreen/core/io/volumereader.h"

namespace voreen {

class Volume;

struct VRN_CORE_API QuadHidacHeader {
    std::string title;
    int subsets;
    int ntf;
    int nbg;
    int nbrg;
    float dxyr;
    float rho;
    std::vector<float> dt;
    tgt::svec3 dims;
    unsigned int frameTime;
    tgt::mat4 transformationMatrix;

    size_t dataOffset; //Start of actual data

    QuadHidacHeader() : title("")
                      , subsets(0)
                      , ntf(1)
                      , nbg(1)
                      , nbrg(1)
                      , dxyr(0.0f)
                      , rho(0.0f)
                      , dims(tgt::svec3(static_cast<size_t>(0)))
                      , frameTime(0)
                      , transformationMatrix(tgt::mat4::identity)
                      , dataOffset(0)
    {
    }
};

/**
 * Reads a volume dataset from a file produced by the quadHIDAC small-animal PET.
 */
class VRN_CORE_API QuadHidacVolumeReader : public VolumeReader {
public:
    QuadHidacVolumeReader(ProgressBar* progress = 0);
    virtual VolumeReader* create(ProgressBar* progress = 0) const;

    virtual std::string getClassName() const   { return "QuadHidacVolumeReader"; }
    virtual std::string getFormatDescription() const { return "quadHIDAC small-animal PET"; }

    QuadHidacHeader readHeader(const std::string& filename) const;

    virtual VolumeList* read(const std::string& url);

    ///Passes the origin's URL to read(url).
    virtual VolumeBase* read(const VolumeURL& origin);

    std::vector<VolumeURL> listVolumes(const std::string& url) const;

    Volume* loadVolume(const std::string& filename, int i, const QuadHidacHeader& header) const;

private:
    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_QUADHIDACVOLUMEREADER_H
