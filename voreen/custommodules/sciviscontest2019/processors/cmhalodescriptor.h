/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2015 University of Muenster, Germany.                        *
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

#ifndef VRN_CMHALODESCRIPTOR_H
#define VRN_CMHALODESCRIPTOR_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/ports/textport.h"
#include "voreen/core/properties/stringexpressionproperty.h"
#include "voreen/core/properties/buttonproperty.h"

#include "../ports/cmhaloport.h"

//use namespace voreen
namespace voreen {

/**
 * Processor that outputs a textual description of the currently selected halo.
 * The content of the description can be changed dynamically much like the output of
 * MetaDataExtractor for Volumes.
 */
class VRN_CORE_API CMHaloDescriptor : public Processor {

public:
    /**
     * Constructor
     */
    CMHaloDescriptor();
    ~CMHaloDescriptor();

    //------------------------------------------
    //  Pure virtual functions of base classes
    //------------------------------------------
    virtual Processor* create() const { return new CMHaloDescriptor();     }
    virtual std::string getClassName() const { return "CMHaloDescriptor";  }
    virtual std::string getCategory() const  { return "Viscontest2015";        }

    virtual void dumpToConsole();
    /**
     * Set the list of MetaData from the Volume to the StringExpressionProperty
     */
    void updateMetaDataList();

protected:

    virtual void setDescriptions() { setDescription("Processor that outputs Text which describes a selected Halo"); }
    virtual void process();

    /**
     * Overwrites the base implementation of this function.
     * It is used to load the needed shader.
     * @see Processor
     */
    virtual void initialize();

    /**
     * Overwrites the base implementation of this function.
     * It is used to free the used shader.
     * @see Processor
     */
    virtual void deinitialize();

private:
    const CMHalo* getSelectedHalo() const;

    //-------------
    //  members
    //-------------

    /// Inport for halo data
    CMHaloPort inport_;
    /// Outport for description of the focused halo
    TextPort outport_;
    /// ID of the currently selected halo
    IntProperty selectedHaloIDProp_;

    /// Holds the expressions that can be exchanged in the text
    StringExpressionProperty expressionProperty_;
    /// Used to dump the outport text to console for debugging purposes
    ButtonProperty dumpButtonProperty_;
    /// Used to update the meta data list
    ButtonProperty updateButtonProperty_;

    /**
     * Gets the placeholders used in the textfield of the StringExpressionProperty,
     * sets the corresponding replacements and returns the text with the replaced placeholders.
     */
    std::string replaceMetaDataAndGetString() const;
};

} // namespace

#endif // VRN_VRN_CMHALODESCRIPTOR_H
