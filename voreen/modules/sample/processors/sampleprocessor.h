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

#ifndef VRN_SAMPLEPROCESSOR_H
#define VRN_SAMPLEPROCESSOR_H

//add base class header
#include "voreen/core/processors/processor.h"

//add used port headers
#include "voreen/core/ports/textport.h"

//add used property headers
#include "voreen/core/properties/stringproperty.h"

//use namespace voreen
namespace voreen {

/**
 * Sample processor, which adds a user-defined prefix to a given text.
 * VRN_CORE_API is a macro needed for shared libs on windows (see voreencoreapi.h)
 */
class VRN_CORE_API SampleProcessor : public Processor {

public:

    /**
     *  Constructor
     */
    SampleProcessor();

    //------------------------------------------
    //  Pure virtual functions of base classes
    //------------------------------------------

    /**
     * Virtual constructor
     * @see VoreenSerializableObject
     */
    virtual Processor* create() const;

    /**
     * Function to get the class name, which is not
     * directly supported by c++.
     * @see VoreenSerializableObject
     */
    virtual std::string getClassName() const;

    /**
     * Function to return the catagory of the processor.
     * It will be shown in the VoreenVE GUI.
     * @see Processor
     */
    virtual std::string getCategory() const;

protected:

    /**
     * Function to set a description of the processors functionality.
     * It will be shown in the VoreenVE GUI.
     * @see Processor
     */
    virtual void setDescriptions();

    /**
     * The main function of each processor.
     * The main functionality should be implemented here.
     * @see Processor
     */
    virtual void process();

private:

    //-------------
    //  members
    //-------------
    TextPort inport_;           ///< inport used to receive the text to be modified
    TextPort outport_;          ///< outport used to pass the modified text
    StringProperty prefixProp_; ///< property for the user-defined prefix
};

} // namespace

#endif // VRN_SAMPLEPROCESSOR_H
