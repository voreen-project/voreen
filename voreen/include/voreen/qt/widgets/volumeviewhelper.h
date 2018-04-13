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

#ifndef VRN_VOLUMEVIEW_H
#define VRN_VOLUMEVIEW_H

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/qt/voreenqtapi.h"

#include <QPixmap>

namespace voreen {

class VRN_QT_API VolumeViewHelper {
public:
    /**
     * Returns a Preview with the height of the given int
     * @param the height of the preview
     */
    static QPixmap generatePreview(const VolumeBase*, int);

    /**
     * Returns a quadratic Preview with the height of the given int and white border
     * @param the height of the preview
     * @param the size of the border
     */
    static QPixmap generateBorderedPreview(const VolumeBase*, int /*height*/, int /* border*/);

    /**
     * Returns a formatted std::string wich contains Volume Information
     */
    static std::string volumeInfoString(const VolumeBase*);

    /**
     * Returns the Volumetype
     */
    static std::string getVolumeType(const VolumeBase*);

    /**
     * Returns the Volumename (incl. path)
     */
    static std::string getVolumeName(const VolumeBase*);

    /**
     * Returns the Volumename (excl. path)
     */
    static std::string getStrippedVolumeName(const VolumeBase*);

    /**
     * Returns the Volumepath
     */
    static std::string getVolumePath(const VolumeBase*);

    /**
     * Returns the dimension
     */
    static std::string getVolumeDimension(const VolumeBase*);

    /**
     * Returns the spacing
     */
    static std::string getVolumeSpacing(const VolumeBase*);

    /**
     * Returns a string desribing the amount of memory
     * that is actually occupied by the volume.
     */
    static std::string getVolumeMemorySize(const VolumeBase*);

    static uint64_t getVolumeMemorySizeByte(const VolumeBase*);

};

} // namespace

#endif
