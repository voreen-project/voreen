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

#ifndef VRN_VOLUMEREPRESENTATION_H
#define VRN_VOLUMEREPRESENTATION_H

#include "voreen/core/voreencoreapi.h"
#include "tgt/vector.h"
#include <vector>

namespace voreen {

/*
 * A VolumeRepresentation is a class storing the volume data (the voxels) in some form.
 * Code should generally work with Volumes instead of VolumeRepresentations.
 */
class VRN_CORE_API VolumeRepresentation {
public:
    VolumeRepresentation(const tgt::svec3& dimensions);
    virtual ~VolumeRepresentation() {}

    /// Returns the VolumeRepresentation's dimensions in voxel coordinates.
    virtual tgt::svec3 getDimensions() const;

    /// Returns the number of voxels contained by the VolumeRepresentation.
    virtual size_t getNumVoxels() const;

    virtual size_t getNumChannels() const = 0;
    virtual size_t getBytesPerVoxel() const = 0;

    /// Returns the format of the volume as string (e.g., "uint8" or "Vector3(float)", @see VolumeFactory).
    virtual std::string getFormat() const = 0;

    /// Returns the base type (e.g., "float" for a representation of format "Vector3(float)").
    virtual std::string getBaseType() const = 0;

    /**
     * Should be called when the volume was modified.
     */
    virtual void invalidate() const {};

protected:
    // protected default constructor
    VolumeRepresentation() {}

    tgt::svec3  dimensions_;
    size_t      numVoxels_;

    static const std::string loggerCat_;
};

//---------------------------------------------------------
//TODO: move to separate file...

class VRN_CORE_API RepresentationConverterBase {
public:
    virtual bool canConvert(const VolumeRepresentation* /*source*/) const = 0;
    virtual VolumeRepresentation* convert(const VolumeRepresentation* /*source*/) const = 0;
    virtual ~RepresentationConverterBase() {}
};

template<class T>
class RepresentationConverter : public RepresentationConverterBase {
};

class VRN_CORE_API ConverterFactory {
public:
    ConverterFactory();
    ~ConverterFactory();

    void addConverter(RepresentationConverterBase* conv) {
        converters_.push_back(conv);
        //TODO: check for duplicates using RTI
    }

    template<class T>
    RepresentationConverter<T>* findConverter(const VolumeRepresentation* source) {
        for(size_t i=0; i<converters_.size(); i++) {
            RepresentationConverter<T>* test = dynamic_cast<RepresentationConverter<T>*> (converters_[i]);

            if(test) {
                //TODO priorities
                if(test->canConvert(source))
                    return test;
            }
        }
        return 0;
    }

protected:
    std::vector<RepresentationConverterBase*> converters_;
};

} // namespace voreen

#endif
