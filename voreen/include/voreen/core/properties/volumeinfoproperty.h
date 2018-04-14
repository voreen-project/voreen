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

#ifndef VRN_VOLUMEINFOPROPERTY_H
#define VRN_VOLUMEINFOPROPERTY_H

#include "property.h"

namespace voreen {

class VolumeBase;

/**
 * Property used to visualize attributes of a given volume.
 * 'This' is used for instance in VolumeSource or VolumeSave.
 */
class VRN_CORE_API VolumeInfoProperty : public Property {
public:
    VolumeInfoProperty(const std::string& id, const std::string& guiText,
                       int invalidationLevel = Processor::INVALID_PARAMETERS, Property::LevelOfDetail lod = Property::LOD_DEFAULT);
    VolumeInfoProperty();

    virtual Property* create() const;

    virtual std::string getClassName() const       { return "VolumeInfoProperty"; }
    virtual std::string getTypeDescription() const { return "VolumeInfo"; }
    void reset();

    /**
     * Sets the volume which attributes will be visualized. It can be null.
     */
    void setVolume(const VolumeBase* handle);

    /**
     * Returns the referenced volume, may be null.
     */
    const VolumeBase* getVolume() const;

private:
    const VolumeBase* volume_; ///< the volume, which infos are displayed

    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_VOLUMEINFOPROPERTY_H
