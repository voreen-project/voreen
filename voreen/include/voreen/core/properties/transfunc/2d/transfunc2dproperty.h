/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#ifndef VRN_TRANSFUNC2DPROPERTY_H
#define VRN_TRANSFUNC2DPROPERTY_H

#include "voreen/core/properties/transfunc/transfuncpropertybase.h"
#include "voreen/core/datastructures/transfunc/2d/transfunc2d.h"

namespace voreen {

/**
 * Generic property for all transfer functions.
 * Used to get rid of dynamic casts.
 */
class VRN_CORE_API TransFunc2DProperty : public TransFuncPropertyBase {
public:

    /**
     * Constructor
     *
     * @param ident identifier that is used in serialization
     * @param guiText text that is shown in the gui
     * @param invalidationLevel The owner is invalidated with this InvalidationLevel upon change.
     * @param lod level of detail in the gui representation
     */
    TransFunc2DProperty(const std::string& ident, const std::string& guiText, int invalidationLevel = Processor::INVALID_RESULT,
                          Property::LevelOfDetail lod = Property::LOD_DEFAULT);
    TransFunc2DProperty();
    virtual ~TransFunc2DProperty();

    //-------------------------------------------------
    //  Property Functions
    //-------------------------------------------------
    /** @see Property::deserialize */
    virtual void deinitialize();
    /**
     * Sets the stored transfer function to the given one.
     *
     * @note The TransFuncProperty takes ownership of the passed
     *  object. Therefore, the caller must not delete it.
     *
     * @param tf transfer function the property is set to
     */
    void set2D(TransFunc2D* tf);

    /**
     * Returns the current transfer function.
     */
    virtual TransFunc2D* get() const;

    //-------------------------------------------------
    //  Member
    //-------------------------------------------------
protected:
    TransFunc2D* transFunc2D_;  ///< pointer to the current function.

    //-------------------------------------------------
    //  Hide base pointer of super class
    //-------------------------------------------------
private:
    using TransFuncPropertyBase::baseFunction_;
    virtual void set(TransFuncBase* tf) {TransFuncPropertyBase::set(tf);}
};

} // namespace voreen

#endif // VRN_TRANSFUNC2DPROPERTY_H
