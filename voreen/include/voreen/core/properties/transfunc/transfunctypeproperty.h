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

#ifndef VRN_TRANSFUNCTYPEPROPERTY_H
#define VRN_TRANSFUNCTYPEPROPERTY_H

#include "voreen/core/properties/optionproperty.h"

#include <set>

namespace voreen {

    class TransFuncPropertyBase;

    /**
     * Enum of all supported tf types.
     * New tf types have to be added here.
     */
    enum TransFuncType {
        TFT_UNSUPPORTED,    ///< default value. should not be used
        TFT_1DKEYS,         ///< support of TransFunc1DKeys
        TFT_1DGAUSSIAN,     ///< support of TransFunc1DGaussian
        TFT_2DPRIMITIVES    ///< support of TransFunc2DPrimitives
    };

/**
 * This property is used to switch between tf property visibilities sorted by class type.
 * If a processor supports more kinds of tf types, it should add this property and add all its tf properties
 * by calling addTransFuncProperty(..). This automatically adds a new option and switches the tf visibilities.
 */
class VRN_CORE_API TransFuncTypeProperty : public OptionProperty<TransFuncType> {
public:
    TransFuncTypeProperty(const std::string& id, const std::string& guiText, int invalidationLevel = Processor::INVALID_PROGRAM, Property::LevelOfDetail lod = Property::LOD_DEFAULT);
    TransFuncTypeProperty();
    virtual ~TransFuncTypeProperty();

    //------------------------------------
    //  important functions
    //------------------------------------
public:
    /**
     * Adds a new tf property. If the tf class is not supported, an assert will be thrown.
     * This function adds a new option to the property according the the passed tf type.
     */
    void addTransFuncProperty(TransFuncPropertyBase* transFuncProperty);
    /**
     * This function do not have to be called, since it will be done in addTransFuncProperty().
     * But sometimes it can be useful to have the option in this property without having any tf properties.
     */
    void addTransFuncTypeOption(TransFuncType option);

    //------------------------------------
    //  override functions
    //------------------------------------
    virtual Property* create() const { return new TransFuncTypeProperty(); }
    virtual std::string getClassName() const       { return "TransFuncTypeProperty"; }
    virtual std::string getTypeDescription() const { return "TransFuncType"; }
    virtual void select(const std::string& key);
    virtual void selectByKey(const std::string& key);
    virtual void invalidate();

    //------------------------------------
    //  helper and member
    //------------------------------------
private:
    /**
     * Helper to update all visible flags.
     */
    void applyTypeSettings();

    std::set<TransFuncPropertyBase*> properties_;    ///< list of added tf properties
};


} // namespace voreen

#endif // VRN_TRANSFUNCTYPEPROPERTY_H
