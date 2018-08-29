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

#ifndef VRN_METADATAEXTRACTOR_H
#define VRN_METADATAEXTRACTOR_H

#include <string>
#include <vector>
#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/textport.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/properties/stringexpressionproperty.h"
#include "voreen/core/properties/buttonproperty.h"

namespace voreen {

/**
 * Provides a textport output with content defined by the Meta Data of the input Volume.
 */
class VRN_CORE_API MetaDataExtractor : public Processor {
public:
    MetaDataExtractor();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "MetaDataExtractor"; }
    virtual std::string getCategory() const  { return "Utility";           }
    virtual CodeState getCodeState() const   { return CODE_STATE_STABLE;   }
    virtual bool isUtility() const           { return true;                }

    virtual bool isReady() const;

    virtual void dumpToConsole();

    /**
     * Set the list of MetaData from the Volume to the StringExpressionProperty
     */
    void updateMetaDataList();

protected:
    virtual void setDescriptions() {
        setDescription("Provides functionality to get the Meta Information of a Volume.");
        outport_.setDescription("An Output String constructed of the added Meta Data.");
        inport_.setDescription("The Volume that provides the Meta Data");
        expressionProperty_.setDescription("Contains the editable text field, a list of placeholders and a button to clear the text field. New placeholders may be added by selecting them from the list.<br> The text field contains normal text as well as placeholders for the values of the Meta Information that have been added. Placeholders support use of the dot operator to get single components of compound meta data, e.g. {Offset.x}. <ul> <li> For vectors, the components \"x\", \"y\", \"z\", \"w\" are allowed. </li> <li> RealWorldMapping meta data supports \"scale\", \"offset\" and \"unit\". </li> <li> For date and time, the following components are allowed: <ul> <li> \"year\", \"month\", \"day\", \"hour\", \"minute\", \"second\", \"millisecond\" (for unformatted output) </li> <li> \"YYYY\", \"MM\", \"DD\", \"hh\", \"mm\", \"ss\", \"fff\" (formatted according to the input).\
 </li> </ul></li> <li> For Memory Size, the size can be converted by using \"B\" for bytes, \"kB\" for kilobytes and \"MB\" for megabytes. </li> </ul>");
    }

    virtual void process();

private:
    VolumePort inport_; ///< input Volume
    TextPort outport_;  ///< put out the Meta Information of the input Volume

    StringExpressionProperty expressionProperty_;
    ButtonProperty dumpButtonProperty_;
    ButtonProperty updateButtonProperty_;

    /**
     * Gets the placeholders used in the textfield of the StringExpressionProperty,
     * sets the corresponding replacements and returns the text with the replaced placeholders.
     */
    std::string replaceMetaDataAndGetString() const;

    static const std::string loggerCat_; ///< category used in logging
};

}   //namespace

#endif
