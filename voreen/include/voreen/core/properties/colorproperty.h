/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#ifndef VRN_COLORPROPERTY_H
#define VRN_COLORPROPERTY_H

#include "voreen/core/properties/vectorproperty.h"

namespace voreen {

    /**
     * Used to represent colors in Voreen. A color is represented in the normalized RGBA color space between 0.f - 1.f.
     * The constructor contains a parameter to enable/disable the alpha channel. A disabled alpha channel will always return
     * 1.f as alpha value.
     */
class VRN_CORE_API ColorProperty : public FloatVec4Property {
    friend class ColorPropertyWidget;
public:
    /** Constructor */
    ColorProperty(const std::string& id, const std::string& guiText,    ///< id and gui name
                tgt::vec4 value = tgt::vec4(1.f),                       ///< normalized initial value
                int invalidationLevel=Processor::INVALID_RESULT,        ///< invalidation level
                Property::LevelOfDetail lod = Property::LOD_DEFAULT,    ///< property visibility
                bool useAlphaChannel = true);                           ///< use alpha? If false, alpha will be always 1.f
    /** Default Constructor */
    ColorProperty();

    /** @override Property::create() */
    virtual Property* create() const override { return new ColorProperty(); }
    /** @override Property::getClassName() */
    virtual std::string getClassName() const override       { return "ColorProperty"; }
    /** @override Property::getTypeDescription() */
    virtual std::string getTypeDescription() const override { return "FloatVector4"; }

    /**
     * @override FloatVec4Property::set()
     * Adjust alpha to 1.f, if alpha channel is not used
     */
    virtual void set(const tgt::vec4& value) override;

    /**
     * Since get() returns a normalized color vector, this function returns a uint_8 vector.
     * @return retuns a uint_8 color vector. Use get() for a normalized vec4() color vector.
     */
    tgt::col4 getRGBA() const;

    //------------------
    //  Members
    //------------------
protected:
    bool useAlphaChannel_; ///< stores, if the alpha channel should be used or not
};

}   // namespace

#endif
