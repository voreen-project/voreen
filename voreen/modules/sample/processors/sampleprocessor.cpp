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

//include header file
#include "sampleprocessor.h"

//using namespace voreen
namespace voreen {

SampleProcessor::SampleProcessor()
    //constructor of base class
    : Processor()
    //constructors of ports
    , inport_(Port::INPORT,             ///< port type, i.e inport
             "inport",                  ///< unique port ID (unique for each processor)
             "Unmodified Text Inport")  ///< port name used in the VoreenVE GUI
    , outport_(Port::OUTPORT,           ///< port type, i.e outport
              "outport",                ///< unique port ID (unique for each processor)
              "Modified Text Outport")  ///< port name used in the VoreenVE GUI
    //constructor of property
    , prefixProp_("prefixProp",         ///< unique property ID (unique for each processor)
                  "Prefix",             ///< property name used in the VoreenVE GUI
                  "Simon says: ")       ///< default value of the property
{
    //register ports
    addPort(inport_);
    addPort(outport_);
    //register properties
    addProperty(prefixProp_);
}

Processor* SampleProcessor::create() const {
    return new SampleProcessor();
}

std::string SampleProcessor::getClassName() const {
    return "SampleProcessor";
}

std::string SampleProcessor::getCategory() const {
    return "Text Processing";
}

void SampleProcessor::setDescriptions() {
    setDescription("Minimal sample processor that appends a user-defined prefix to a given text.");
}

void SampleProcessor::process() {
    //get inport data
    std::string inString = inport_.getData();
    //get prefix string
    std::string prefixString = prefixProp_.get();
    //combine both strings
    std::string outString = prefixString + inString;
    //set outport data
    outport_.setData(outString);
}

} // namespace
