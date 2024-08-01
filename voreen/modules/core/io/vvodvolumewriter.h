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

#ifndef VRN_VVODVOLUMEWRITER_H
#define VRN_VVODVOLUMEWRITER_H

#include "voreen/core/io/volumewriter.h"
#include "voreen/core/datastructures/octree/volumeoctree.h"
#include <sstream>
#include "voreen/core/io/progressbar.h"

namespace voreen {

/**
 * Writes the volume into a .vvd and a .raw file (Voreen Volume Data, new Voreen format).
 */
class VRN_CORE_API VvodVolumeWriter : public VolumeWriter {
public:
    VvodVolumeWriter(ProgressReporter* p = 0);
    virtual VolumeWriter* create(ProgressBar* progress = 0) const;

    virtual std::string getClassName() const   { return "VvodVolumeWriter"; }
    virtual std::string getFormatDescription() const { return "Voreen Volume Octree Data"; }

    /**
     * Writes the data of a volume octree representation into a vvod-file and several brick buffers within a subfolder.
     *
     * @param fileName File name to be written
     * @param volume Volume dataset
     */
    virtual void write(const std::string& filename, const VolumeBase* volumeHandle);

protected:
    std::string getOctreeStoragePath(const VolumeOctree* octree) const;

private:

    ProgressReporter* progressBar_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif
