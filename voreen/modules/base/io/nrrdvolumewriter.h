/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
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

#ifndef VRN_NRRDVOLUMEWRITER_H
#define VRN_NRRDVOLUMEWRITER_H

#include "voreen/core/io/volumewriter.h"

namespace voreen {

/**
 * Writer for <tt>.nrrd</tt> volume files (nearly raw raster data).
 * Writes the volume into a .nhdr and a .raw file.
 *
 * See http://teem.sourceforge.net/nrrd/ for details about the file format.
 */
class VRN_CORE_API NrrdVolumeWriter : public VolumeWriter {
public:
    NrrdVolumeWriter();
    virtual VolumeWriter* create(ProgressBar* progress = 0) const;

    virtual std::string getClassName() const   { return "NrrdVolumeWriter"; }
    virtual std::string getFormatDescription() const { return "Nearly Raw Raster Data"; }

    /**
     * Writes the data of a volume into a nrrd- and a raw-file.
     *
     * @param fileName File name to be written
     * @param volume Volume dataset
     */
    virtual void write(const std::string& filename, const VolumeBase* volume);

private:
    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_NRRDVOLUMEWRITER_H
