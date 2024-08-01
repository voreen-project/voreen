/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include "voreen/core/io/serialization/serializable.h"
#include "voreen/core/io/serialization/serialization.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/utils/stringutils.h"

#include "baseclasses.h"
#include "template.h"
#include "itk_wrapper.h"

using namespace voreen;
using namespace std;

/***********************************************************************************
 * Class Argument
 ***********************************************************************************/
void Argument::serialize(Serializer& s) const {
    s.serialize("name", name_);
    s.serialize("argumenttype", type_);
    s.serialize("inputArgument", inputArgument_);

    if (type_ != "KernelType" && type_ != "AddSeedType"
        && type_ != "SetSeedType" && type_ != "VectorSeedType" && inputArgument_) {
        s.serialize("defaultValue", defaultValue_);
        if (type_ != "Bool") {
            s.serialize("minValue", minValue_);
            s.serialize("maxValue", maxValue_);
        }
    }
}

void Argument::deserialize(Deserializer& s) {
    try {
        s.deserialize("name", name_);
        s.deserialize("argumenttype", type_);
    }
    catch (SerializationNoSuchDataException e) {
    }
    s.optionalDeserialize<bool>("inputArgument", inputArgument_, true);

    if (type_ != "KernelType" && type_ != "AddSeedType"
        && type_ != "SetSeedType" && type_ != "VectorSeedType") {
        if (type_ != "PixelType") {
            try {
                s.deserialize("defaultValue", defaultValue_);
            }
            catch (SerializationNoSuchDataException e) {
            }
        }
        if (type_ != "Bool") {
            try {
                s.deserialize("minValue", minValue_);
                s.deserialize("maxValue", maxValue_);
            }
            catch (SerializationNoSuchDataException e) {
            }
        }
    }
}

std::string Argument::getPropertyname() {
    std::string propertyname = name_ + "_";
    propertyname[0] = tolower(propertyname[0]);
    return propertyname;
}

std::string Argument::getProperty() {
    std::string result;
    result = result + "    " + getPropertyname() + "(\"" + getPropertyID()
        + "\", \"" + name_;

    if (type_ == "Float" || type_ == "Int") {
        if (inputArgument_) {
            result = result + "\", " + defaultValue_ + ", " + minValue_ + ", " + maxValue_ + ")";
        }
        else {
            if (type_ == "Float") {
                result = result + "\", 0.f, -FLT_MAX, FLT_MAX)";
            }
            else {
                result = result + "\", 0, INT_MIN, INT_MAX)";
            }
        }
    }
    else if (type_ == "Bool") {
        if (inputArgument_) {
            result = result + "\", " + defaultValue_ + ")";
        }
        else {
            result = result + "\", false)";
        }
    }
    else if (type_ == "PixelType") {
        if (inputArgument_) {
            result = result + "\")";
        }
        else {
            result = result + "\")";
        }
    }
    else if (type_ == "SizeType" || type_ == "IndexType") {
        if (inputArgument_) {
            result = result + "\", tgt::ivec3" + defaultValue_ + ", tgt::ivec3" + minValue_
                + ", tgt::ivec3" + maxValue_ + ")";
        }
        else {
            result = result + "\", tgt::ivec3 (0), tgt::ivec3 (INT_MIN)"
                + ", tgt::ivec3 (INT_MAX))";
        }
    }
    else if (type_ == "ArrayType") {
        if (inputArgument_) {
            result = result + "\", tgt::vec3" + defaultValue_ + ", tgt::vec3" + minValue_
                + ", tgt::vec3" + maxValue_ + ")";
        }
        else {
            result = result + "\", tgt::vec3 (0.f), tgt::vec3 (-FLT_MAX)"
                + ", tgt::vec3 (FLT_MAX))";
        }
    }
    else {
        result = "";
    }

    return result;
}

std::string Argument::getFilterSetOrGet() {
    std::string result;

    if (type_ == "Float" || type_ == "Int" || type_ == "Bool") {
        if (inputArgument_) {
            result = result + "    filter->Set" + name_ + "(" + getPropertyname() + ".get()" + ");\n";
        }
        else {
            result = result + "        " + getPropertyname() + ".set(filter->Get" + name_ + "());\n";
        }
    }
    else if (type_ == "PixelType") {
        if (inputArgument_) {
            result = result + "    " + "filter->Set" + name_+ "(" + getPropertyname()
                + ".getValue<T>()" + ");\n";
        }
        else {
            result = result + "        " + getPropertyname() + ".setValue<T>(filter->Get"
                + name_ + "());\n";
        }
    }
    else if (type_ == "SizeType" || type_ == "IndexType") {
        if (inputArgument_) {
            result = result + "    typename InputImageType1::" + type_ + " " + getTempName() + ";\n"
                + "    " + getTempName() + "[0] = " + getPropertyname() + ".get()" + ".x;\n"
                + "    " + getTempName() + "[1] = " + getPropertyname() + ".get()" + ".y;\n"
                + "    " + getTempName() + "[2] = " + getPropertyname() + ".get()" + ".z;\n"
                + "    " + "filter->Set" + name_ + "(" + getTempName() + ");\n";
        }
        else {
            result = result + "    typename InputImageType1::" + type_ + " " + getTempName() + ";\n"
                + "    " + getTempName() + " = filter->Get" + name_ + "();\n"
                + "    tgt::ivec3 help;\n"
                + "    help.x = " + getTempName() + "[0];\n"
                + "    help.y = " + getTempName() + "[1];\n"
                + "    help.z = " + getTempName() + "[2];\n"
                + "    " + getPropertyname() + ".set(help);\n";
        }
    }
    else if (type_ == "ArrayType") {
        if (inputArgument_) {
            result = result + "   typename FilterType::" + type_ + " " + getTempName() + ";\n"
                + "    " + getTempName() + "[0] = " + getPropertyname() + ".get()" + ".x;\n"
                + "    " + getTempName() + "[1] = " + getPropertyname() + ".get()" + ".y;\n"
                + "    " + getTempName() + "[2] = " + getPropertyname() + ".get()" + ".z;\n"
                + "    " + "filter->Set" + name_ + "(" + getTempName() + ");\n";
        }
        else {
            result = result + "   typename FilterType::" + type_ + " " + getTempName() + ";\n"
                + "    " + getTempName() + " = filter->Get" + name_ + "();\n"
                + "    tgt::vec3 help;\n"
                + "    help.x = " + getTempName() + "[0];\n"
                + "    help.y = " + getTempName() + "[1];\n"
                + "    help.z = " + getTempName() + "[2];\n"
                + "    " + getPropertyname() + ".set(help);\n";
        }
    }
    else if (type_ == "AddSeedType") {
        std::string help = "";
        if (name_ == "Seed1")
            help = "1";
        if (name_ == "Seed2")
            help = "2";

        result = result + "\n"
            + "    if (seedPointPort" + help + "_.hasChanged()) {\n"
            + "        const PointListGeometry<tgt::vec3>* pointList" + help
            + " = dynamic_cast< const PointListGeometry<tgt::vec3>* >(seedPointPort" + help
            + "_.getData());\n"
            + "        if (pointList" + help + ") {\n"
            + "            seedPoints" + help + " = pointList" + help + "->getData();\n"
            + "        }\n"
            + "    }\n\n"
            + "    filter->ClearSeeds" + help + "();\n"
            + "    typename InputImageType1::IndexType seed" + help + ";\n"
            + "    for (size_t i = 0; i < seedPoints" + help + ".size(); i++) {\n"
            + "        seed" + help + "[0] = seedPoints" + help + "[i].x;\n"
            + "        seed" + help + "[1] = seedPoints" + help + "[i].y;\n"
            + "        seed" + help + "[2] = seedPoints" + help + "[i].z;\n"
            + "        filter->Add" + name_ + "(seed" + help + ");\n"
            + "    }\n";
    }
    else if (type_ == "SetSeedType") {
        std::string help = "";
        if (name_ == "Seed1")
            help = "1";
        if (name_ == "Seed2")
            help = "2";

        result = result + "\n"
            + "    if (seedPointPort" + help + "_.hasChanged()) {\n"
            + "        const PointListGeometry<tgt::vec3>* pointList" + help
            + " = dynamic_cast< const PointListGeometry<tgt::vec3>* >(seedPointPort" + help
            + "_.getData());\n"
            + "        if (pointList" + help + ") {\n"
            + "            seedPoints" + help + " = pointList" + help + "->getData();\n"
            + "            if (!seedPoints" + help + ".empty()) {\n"
            + "                numSeedPoint" + help + "_.setMinValue(1);\n"
            + "                numSeedPoint" + help + "_.setMaxValue(seedPoints" + help + ".size());\n"
            + "                numSeedPoint" + help + "_.set(seedPoints" + help + ".size());\n"
            + "            }\n"
            + "            else {\n"
            + "                numSeedPoint" + help + "_.setMinValue(0);\n"
            + "                numSeedPoint" + help + "_.setMaxValue(0);\n"
            + "                numSeedPoint" + help + "_.set(0);\n"
            + "            }\n"
            + "        }\n"
            + "    }\n\n"
            + "    if(!seedPoints" + help + ".empty()) {\n"
            + "        seedPoint" + help + "_.set(seedPoints" + help + "[numSeedPoint" + help
            + "_.get()-1]);\n"
            + "    }\n"
            + "    else {\n"
            + "        seedPoint" + help + "_.set(tgt::vec3 (1));\n"
            + "    }\n"
            + "    typename InputImageType1::IndexType seed" + help + ";\n"
            + "    seed" + help + "[0] = seedPoint" + help + "_.get().x;\n"
            + "    seed" + help + "[1] = seedPoint" + help + "_.get().y;\n"
            + "    seed" + help + "[2] = seedPoint" + help + "_.get().z;\n\n"
            + "    filter->Set" + name_ + "(seed" + help + ");\n";
    }
    else if (type_ == "VectorSeedType") {
        std::string help;
        if (name_ == "SeedPoints1")
            help = "1";
        if (name_ == "SeedPoints2")
            help = "2";

        result = result + "\n"
            + "    if (seedPointPort" + help + "_.hasChanged()) {\n"
            + "        const PointListGeometry<tgt::vec3>* pointList" + help
            + " = dynamic_cast< const PointListGeometry<tgt::vec3>* >(seedPointPort" + help
            + "_.getData());\n"
            + "        if (pointList" + help + ") {\n"
            + "            seedPoints" + help + " = pointList" + help + "->getData();\n"
            + "        }\n"
            + "    }\n\n"
            + "    filter->Set" + name_ + "(seedPoints" + help + ");\n";
    }
    return result;
}

std::string Argument::getName() {
    return name_;
}

std::string Argument::getType() {
    return type_;
}

bool Argument::isInputArgument() {
    return inputArgument_;
}

std::string Argument::getDefaultValue() {
    return defaultValue_;
}

std::string Argument::getMinValue() {
    return minValue_;
}

std::string Argument::getMaxValue() {
    return maxValue_;
}

std::string Argument::getPropertyID() {
    std::string propertyID = name_;
    propertyID[0] = tolower(propertyID[0]);
    return propertyID;
}

std::string Argument::getTempName() {
    std::string tempName = name_;
    tempName[0] = tolower(tempName[0]);
    return tempName;
}


/***********************************************************************************
 * Class Ports
 ***********************************************************************************/
void Ports::serialize(Serializer& s) const {
    s.serialize("name", name_);
    s.serialize("filterTemplate", filterTemplate_);
    s.serialize("nameIsSetter", nameIsSetter_);
    s.serialize("transferRMW", transferRWM_);
    s.serialize("transferTransformation", transferTransformation_);
    s.serialize("possibleTypes", possibleTypes_, "type");
}

void Ports::deserialize(Deserializer& s) {
    s.deserialize("name", name_);
    s.optionalDeserialize<bool>("filterTemplate", filterTemplate_, true);
    s.optionalDeserialize<bool>("nameIsSetter", nameIsSetter_, false);
    s.optionalDeserialize<bool>("transferRMW", transferRWM_, true);
    s.optionalDeserialize<bool>("transferTransformation", transferTransformation_, true);
    try {
        s.deserialize("possibleTypes", possibleTypes_, "type");
    }
    catch (SerializationNoSuchDataException e) {
    }
}

std::string Ports::getITKTypename(int index) {
    std::string typeName = possibleTypes_[index];

    if (typeName == "UInt8" || typeName == "4xUInt8" || typeName == "3xUInt8"
        || typeName == "2xUInt8") {
        return "uint8_t";
    }
    else if (typeName == "Int8" || typeName == "4xInt8" || typeName == "3xInt8"
        || typeName == "2xInt8") {
        return "int8_t";
    }
    else if (typeName == "UInt16" || typeName == "4xUInt16" || typeName == "3xUInt16"
        || typeName == "2xUInt16") {
        return "uint16_t";
    }
    else if (typeName == "Int16" || typeName == "4xInt16" || typeName == "3xInt16"
        || typeName == "2xInt16") {
        return "int16_t";
    }
    else if (typeName == "UInt32") {
        return "uint32_t";
    }
    else if (typeName == "Int32") {
        return "int32_t";
    }
    else if (typeName == "UInt64") {
        return "uint64_t";
    }
    else if (typeName == "Int64") {
        return "int64_t";
    }
    else if (typeName == "Float" || typeName == "4xFloat" || typeName == "3xFloat"
        || typeName == "2xFloat") {
        return "float";
    }
    else if (typeName == "Double" || typeName == "4xDouble" || typeName == "3xDouble"
        || typeName == "2xDouble") {
        return "double";
    }
    return "";
}

std::string Ports::getName() {
    return name_;
}

void Ports::setName(std::string name) {
    name_ = name;
}

bool Ports::isFilterTemplate() {
    return filterTemplate_;
}

void Ports::setFilterTemplate(bool filterTemplate) {
    filterTemplate_ = filterTemplate;
}

bool Ports::isNameIsSetter() {
    return nameIsSetter_;
}

std::vector<std::string> Ports::getPossibleTypes() {
    return possibleTypes_;
}

void Ports::setNameIsSetter(bool nameIsSetter) {
    nameIsSetter_ = nameIsSetter;
}

bool Ports::transferRWM() {
    return transferRWM_;
}

void Ports::setTransferRMW(bool tRWM) {
    transferRWM_ = tRWM;
}

bool Ports::transferTransformation() {
    return transferTransformation_;
}

void Ports::setTransferTransformation(bool tTrafo) {
    transferTransformation_ = tTrafo;
}

void Ports::setPossibleTypes(std::vector<std::string> types) {
    possibleTypes_ = types;
}

void Ports::addPossibleType(std::string type) {
    possibleTypes_.push_back(type);
}

void Ports::deleteMetaType(size_t position) {
    possibleTypes_.erase(possibleTypes_.begin()+position);
}


/***********************************************************************************
 * Class Filter
 ***********************************************************************************/
void Filter::serialize(Serializer& s) const {
    s.serialize("name", name_);
    s.serialize("enabled", enabled_);
    s.serialize("description", description_);
    s.serialize("autoGenerated", autoGenerated_);
    s.serialize("codeState", codeState_);
    s.serialize("inports", inports_, "port");
    s.serialize("outports", outports_, "port");
    s.serialize("arguments", arguments_, "argument");
}

void Filter::deserialize(Deserializer& s) {
    s.deserialize("name", name_);
    s.optionalDeserialize<bool>("enabled", enabled_, true);
    s.optionalDeserialize<std::string>("description", description_, "");
    s.optionalDeserialize<bool>("autoGenerated", autoGenerated_, true);
    s.optionalDeserialize<std::string>("codeState", codeState_, "EXPERIMENTAL");

    try {
        s.deserialize("inports", inports_, "port");
    }
    catch (SerializationNoSuchDataException e) {
    }
    try {
        s.deserialize("outports", outports_, "port");
    }
    catch (SerializationNoSuchDataException e) {
    }
    try {
        s.deserialize("arguments", arguments_, "argument");
    }
    catch (SerializationNoSuchDataException e) {
    }
}

void Filter::initializeFilterPorts(std::vector<std::string> integerScalarTypes
                                   , std::vector<std::string> realScalarTypes
                                   , std::vector<std::string> scalarTypes
                                   , std::vector<std::string> integerVectorTypes
                                   , std::vector<std::string> realVectorTypes
                                   , std::vector<std::string> vectorTypes) {

    scalarVolumeTypesPossible_ = false;
    vec2VolumeTypesPossible_ = false;
    vec3VolumeTypesPossible_ = false;
    vec4VolumeTypesPossible_ = false;

    if (inports_.empty() && outports_.empty()) {
        Ports inport;
        inport.setName("InputImage");
        inport.setFilterTemplate(true);
        inport.setNameIsSetter(false);
        inports_.push_back(inport);
        Ports outport;
        outport.setName("OutputImage");
        outport.setFilterTemplate(true);
        outport.setNameIsSetter(false);
        outports_.push_back(outport);
    }

    for (size_t i=0; i<inports_.size(); ++i){
        if (inports_[i].getPossibleTypes().empty()) {
            inports_[i].setPossibleTypes(scalarTypes);
            scalarVolumeTypesPossible_ = true;
        }
        else {
            std::string volumeType;
            for (size_t j=0; j<inports_[i].getPossibleTypes().size(); ++j) {
                volumeType = inports_[i].getPossibleTypes()[j];
                if (volumeType == "Scalar") {
                    scalarVolumeTypesPossible_ = true;
                    for (size_t k=0; k<scalarTypes.size(); ++k) {
                        inports_[i].addPossibleType(scalarTypes[k]);
                    }
                    inports_[i].deleteMetaType(j);
                    --j;
                }
                else if (volumeType == "IntegerScalar") {
                    scalarVolumeTypesPossible_ = true;
                    for (size_t k=0; k<integerScalarTypes.size(); ++k) {
                        inports_[i].addPossibleType(integerScalarTypes[k]);
                    }
                    inports_[i].deleteMetaType(j);
                    --j;
                }
                else if (volumeType == "RealScalar") {
                    scalarVolumeTypesPossible_ = true;
                    for (size_t k=0; k<realScalarTypes.size(); ++k) {
                        inports_[i].addPossibleType(realScalarTypes[k]);
                    }
                    inports_[i].deleteMetaType(j);
                    --j;
                }
                else if (volumeType == "Vector") {
                    vec2VolumeTypesPossible_ = true;
                    vec3VolumeTypesPossible_ = true;
                    vec4VolumeTypesPossible_ = true;
                    for (size_t k=0; k<vectorTypes.size(); ++k) {
                        inports_[i].addPossibleType(vectorTypes[k]);
                    }
                    inports_[i].deleteMetaType(j);
                    --j;
                }
                else if (volumeType == "IntegerVector") {
                    vec2VolumeTypesPossible_ = true;
                    vec3VolumeTypesPossible_ = true;
                    vec4VolumeTypesPossible_ = true;
                    for (size_t k=0; k<integerVectorTypes.size(); ++k) {
                        inports_[i].addPossibleType(integerVectorTypes[k]);
                    }
                    inports_[i].deleteMetaType(j);
                    --j;
                }
                else if (volumeType == "RealVector") {
                    vec2VolumeTypesPossible_ = true;
                    vec3VolumeTypesPossible_ = true;
                    vec4VolumeTypesPossible_ = true;
                    for (size_t k=0; k<realVectorTypes.size(); ++k) {
                        inports_[i].addPossibleType(realVectorTypes[k]);
                    }
                    inports_[i].deleteMetaType(j);
                    --j;
                }
                else if (volumeType == "UInt8" || volumeType == "Int8"
                    || volumeType == "UInt16" || volumeType == "Int16"
                    || volumeType == "UInt32" || volumeType == "Int32"
                    || volumeType == "UInt64" || volumeType == "Int64"
                    || volumeType == "Float" || volumeType == "Double") {
                    scalarVolumeTypesPossible_ = true;
                }
                else if (volumeType == "2xUInt8" || volumeType == "2xInt8"
                    || volumeType == "2xUInt16" || volumeType == "2xInt16"
                    || volumeType == "2xFloat" || volumeType == "2xDouble") {
                    vec2VolumeTypesPossible_ = true;
                }
                else if (volumeType == "3xUInt8" || volumeType == "3xInt8"
                    || volumeType == "3xUInt16" || volumeType == "3xInt16"
                    || volumeType == "3xFloat" || volumeType == "3xDouble") {
                    vec3VolumeTypesPossible_ = true;
                }
                else if (volumeType == "4xUInt8" || volumeType == "4xInt8"
                    || volumeType == "4xUInt16" || volumeType == "4xInt16"
                    || volumeType == "4xFloat" || volumeType == "4xDouble") {
                    vec4VolumeTypesPossible_ = true;
                }
            }
        }
    }
}

void Filter::createCPP(std::string moduleName) {
    // replacements depending on input and output of the filter
    std::string rep_ports;
    std::string rep_portadd;

    // replacements depending on the arguments of the filter
    std::string rep_properties;
    std::string rep_propertyadd;

    // replacement for filter-methods
    std::string rep_filtermethods;

    // replacement for the processing methods
    std::string rep_processingmethods;

    // helpers
    std::string nextline;
    std::string number;
    std::string outputname;
    std::string volumeType;

    // create input-part of the replacements depending on input and output
    for (size_t i=0; i<inports_.size(); ++i) {
        nextline = "";
        if (i<inports_.size()-1 || outports_.size()>0)
            nextline = ",\n";

        number = itos(i+1);
        rep_ports = rep_ports + "    inport" + number + "_(Port::INPORT, \""
            + inports_[i].getName() + "\")" + nextline;

        rep_portadd = rep_portadd + "    addPort(inport" + number + "_);\n"
            + "    PortConditionLogicalOr* orCondition" + number
            + " = new PortConditionLogicalOr();\n";

        for (size_t j=0; j<inports_[i].getPossibleTypes().size(); ++j) {
            volumeType = inports_[i].getPossibleTypes()[j];
            rep_portadd = rep_portadd + "    orCondition" + number
                + "->addLinkedCondition(new PortConditionVolumeType" + volumeType + "());\n";
        }
        rep_portadd = rep_portadd + "    inport" + number + "_.addCondition(orCondition"
            + number + ");\n";
    }

    // create output-part of the replacements depending on input and output
    for (size_t i=0; i<outports_.size(); ++i) {
        number = itos(i+1);
        nextline = "";
        if (i<outports_.size()-1)
            nextline = ",\n";

        rep_ports = rep_ports + "    outport" + itos(i+1) + "_(Port::OUTPORT, \""
            + outports_[i].getName() + "\")" + nextline;
        rep_portadd = rep_portadd + "    addPort(outport" + itos(i+1) + "_);\n";
    }

    // for standard-filters with only one inport and one outport an enable-property will be added
    if (inports_.size() == 1 && outports_.size() == 1) {
        rep_properties = rep_properties + ",\n"
            + "    enableProcessing_(\"enabled\", \"Enable\", false)";
        rep_propertyadd = rep_propertyadd + "    addProperty(enableProcessing_);\n";
    }

    Argument argument;
    bool kernel = false;

    // create replacements depending on the arguments of the filter
    for (size_t j=0; j<arguments_.size(); ++j) {
        argument = arguments_[j];

        if (argument.getType() != "KernelType" && argument.getType() != "AddSeedType"
            && argument.getType() != "SetSeedType" && argument.getType() != "VectorSeedType") {
            rep_properties = rep_properties + ",\n";
        }
        if (argument.getType() == "KernelType") {
            kernel = true;
        }
        else if (argument.getType() == "AddSeedType" || argument.getType() == "SetSeedType"
            || argument.getType() == "VectorSeedType") {
            std::string help = "";
            if (argument.getName() == "Seed1" || argument.getName() == "SeedPoints1")
                help = "1";
            if (argument.getName() == "Seed2" || argument.getName() == "SeedPoints2")
                help = "2";

            rep_ports = rep_ports + ",\n    seedPointPort" + help
                + "_(Port::INPORT, \"seedPointInput" + help + "\")";
            rep_portadd = rep_portadd + "    addPort(seedPointPort" + help + "_);\n";

            if (argument.getType() == "SetSeedType") {
                rep_properties = rep_properties + ",\n    numSeedPoint" + help
                    + "_(\"numSeedPoint" + help + "\", \"Take Seed-Point\", 0, 0, 0, Processor::VALID),\n"
                    + "    seedPoint" + help + "_(\"SeedPoint" + help + "\", \"Seed-Point" + help
                    + "\", tgt::vec3(1),tgt::vec3(1),tgt::vec3(500))";
                rep_propertyadd = rep_propertyadd + "    addProperty(numSeedPoint" + help + "_);\n"
                    + "    addProperty(seedPoint" + help + "_);\n"
                    + "    numSeedPoint" + help + "_.onChange(MemberFunctionCallback<" + name_
                    + "ITK>(this, &" + name_ + "ITK::process));\n"
                    + "    seedPoint" + help + "_.setReadOnlyFlag(true);\n";
            }
        }
        else {
            rep_properties = rep_properties + argument.getProperty();
            rep_propertyadd = rep_propertyadd + "    addProperty("
                + argument.getPropertyname() + ");\n";
            if (!argument.isInputArgument()) {
                rep_propertyadd = rep_propertyadd + "    " + argument.getPropertyname()
                    + ".setReadOnlyFlag(true);\n";
            }
        }
    }

    // create replacement for filter-methods
    if (scalarVolumeTypesPossible_) {
        rep_filtermethods = rep_filtermethods + getFilterMethod("", kernel);
    }
    if (vec2VolumeTypesPossible_) {
        rep_filtermethods = rep_filtermethods + getFilterMethod("Vec2", kernel);
    }
    if (vec3VolumeTypesPossible_) {
        rep_filtermethods = rep_filtermethods + getFilterMethod("Vec3", kernel);
    }
    if (vec4VolumeTypesPossible_) {
        rep_filtermethods = rep_filtermethods + getFilterMethod("Vec4", kernel);
    }

    // create replacement for the processing methods
    rep_processingmethods = getProcessingMethods();

    // create template and replace the keywords by the generated replacements
    std::string templateFilename;
    if (kernel) {
        rep_properties = rep_properties + ",";
        templateFilename = getTemplatePath("itk_filter_cpp_KERNEL.txt");
    }
    else {
        templateFilename = getTemplatePath("itk_filter_cpp.txt");
    }

    Template t(templateFilename);
    t.replace("§LOWERFILTERNAME§", getFilename());
    t.replace("§FILTERNAME§", name_);

    t.replace("§PORTS§", rep_ports);
    t.replace("§PORTADD§", rep_portadd);

    t.replace("§PROPERTIES§", rep_properties);
    t.replace("§PROPERTYADD§", rep_propertyadd);

    t.replace("§FILTERMETHODS§", rep_filtermethods);
    t.replace("§PROCESSINGMETHODS§", rep_processingmethods);

    std::string filterFilename;
    filterFilename = filterFilename + getFilename() + ".cpp";
    t.writeToFile(getProcessorOutputPath(moduleName, filterFilename));
}

void Filter::createHeaderfile(std::string moduleName, std::string group) {
    Argument argument;

    // replacements depending on input and output of the filter
    std::string rep_ports;
    std::string rep_volumetypeswitch;

    // replacement for filter-methods
    std::string rep_filtermethods;

    // replacements depending on the arguments of the filter
    std::string rep_propertyinclude;
    std::string rep_propertydef;

    // helper
    std::string templates;


    // create replacements depending on input and output of the filter
    for (size_t i=0; i<inports_.size(); ++i) {
        rep_ports = rep_ports + "    VolumePort inport" + itos(i+1) + "_;\n";
    }
    for (size_t i=0; i<outports_.size(); ++i) {
        rep_ports = rep_ports + "    VolumePort outport" + itos(i+1) + "_;\n";
    }

    if (inports_.size() >= 2) {
        templates = ", class S";
        rep_volumetypeswitch = rep_volumetypeswitch + "    template<class T>\n"
            + "    void volumeTypeSwitch1();\n";

        if (inports_.size() >= 3) {
            templates = templates + ", class R";
            rep_volumetypeswitch = rep_volumetypeswitch + "\n    template<class T, class S>\n"
                + "    void volumeTypeSwitch2();\n";
        }
    }

    // create replacement for filter-methods
    if (scalarVolumeTypesPossible_) {
        rep_filtermethods = rep_filtermethods + "    template<class T" + templates + ">\n"
            + "    void " + getProcessFunctionName() + "ITK();\n";
    }
    if (vec2VolumeTypesPossible_) {
        rep_filtermethods = rep_filtermethods + "    template<class T" + templates + ">\n"
            + "    void " + getProcessFunctionName() + "Vec2ITK();\n";
    }
    if (vec3VolumeTypesPossible_) {
        rep_filtermethods = rep_filtermethods + "    template<class T" + templates + ">\n"
            + "    void " + getProcessFunctionName() + "Vec3ITK();\n";
    }
    if (vec4VolumeTypesPossible_) {
        rep_filtermethods = rep_filtermethods + "    template<class T" + templates + ">\n"
            + "    void " + getProcessFunctionName() + "Vec4ITK();\n";
    }

    bool boolInclude = false;
    bool intInclude = false;
    bool floatInclude = false;
    bool voxelInclude = false;
    bool vectorInclude = false;
    bool pointlistGeometryInclude = false;

    // create replacements depending on the arguments of the filter
    if (inports_.size() == 1 && outports_.size() == 1) {
        if (!boolInclude) {
            rep_propertyinclude = rep_propertyinclude
                + "#include \"voreen/core/properties/boolproperty.h\"\n";
            boolInclude = true;
        }
        rep_propertydef = rep_propertydef + "    BoolProperty enableProcessing_;\n";
    }

    for (size_t j=0 ; j<arguments_.size(); ++j) {
        argument = arguments_[j];

        if (argument.getType() == "KernelType") {
            rep_propertyinclude = rep_propertyinclude
                + "#include \"voreen/core/properties/optionproperty.h\"\n"
                + "#include \"voreen/core/properties/vectorproperty.h\"\n";
            rep_propertydef = rep_propertydef+"    StringOptionProperty structuringElement_;\n"
                + "    StringOptionProperty shape_;\n"
                + "    IntVec3Property radius_;\n";
        }
        else if (argument.getType() == "Bool") {
            if (!boolInclude) {
                rep_propertyinclude = rep_propertyinclude
                    + "#include \"voreen/core/properties/boolproperty.h\"\n";
                boolInclude = true;
            }
            rep_propertydef = rep_propertydef + "    BoolProperty "
                + argument.getPropertyname() + ";\n";
        }
        else if (argument.getType() == "Int") {
            if (!intInclude) {
                rep_propertyinclude = rep_propertyinclude
                    + "#include \"voreen/core/properties/intproperty.h\"\n";
                intInclude = true;
            }
            rep_propertydef = rep_propertydef + "    IntProperty "
                + argument.getPropertyname() + ";\n";
            }
        else if (argument.getType() == "Float") {
            if (!floatInclude) {
                rep_propertyinclude = rep_propertyinclude
                    + "#include \"voreen/core/properties/floatproperty.h\"\n";
                floatInclude = true;
            }
            rep_propertydef = rep_propertydef + "    FloatProperty "
                + argument.getPropertyname() + ";\n";
        }
        else if (argument.getType() == "PixelType") {
            if (!voxelInclude) {
                rep_propertyinclude = rep_propertyinclude
                    + "#include \"voreen/core/properties/voxeltypeproperty.h\"\n";
                voxelInclude = true;
            }
            rep_propertydef = rep_propertydef + "    VoxelTypeProperty "
                + argument.getPropertyname() + ";\n";
        }
        else if (argument.getType() == "SizeType" || argument.getType() == "IndexType") {
            if (!vectorInclude) {
                rep_propertyinclude = rep_propertyinclude
                    + "#include \"voreen/core/properties/vectorproperty.h\"\n";
                vectorInclude = true;
            }
            rep_propertydef = rep_propertydef + "    IntVec3Property "
                + argument.getPropertyname() + ";\n";
        }
        else if (argument.getType() == "ArrayType") {
            if (!vectorInclude) {
                rep_propertyinclude = rep_propertyinclude
                    + "#include \"voreen/core/properties/vectorproperty.h\"\n";
                vectorInclude = true;
            }
            rep_propertydef = rep_propertydef + "    FloatVec3Property "
                + argument.getPropertyname() + ";\n";
        }
        else if (argument.getType() == "AddSeedType" || argument.getType() == "SetSeedType"
            || argument.getType() == "VectorSeedType") {
            std::string help = "";
            if (argument.getName() == "Seed1" || argument.getName() == "SeedPoints1")
                help = "1";
            if (argument.getName() == "Seed2" || argument.getName() == "SeedPoints2")
                help = "2";

            rep_ports = rep_ports + "    GeometryPort seedPointPort" + help + "_;\n";
            if (!pointlistGeometryInclude) {
                rep_propertyinclude = rep_propertyinclude
                    + "#include \"voreen/core/datastructures/geometry/pointlistgeometry.h\"\n";
                pointlistGeometryInclude = true;
            }
            if (argument.getType() == "SetSeedType") {
                if (!intInclude) {
                    rep_propertyinclude = rep_propertyinclude
                        + "#include \"voreen/core/properties/intproperty.h\"\n";
                    intInclude = true;
                }
                if (!vectorInclude) {
                    rep_propertyinclude = rep_propertyinclude
                        + "#include \"voreen/core/properties/vectorproperty.h\"\n";
                    vectorInclude = true;
                }
                rep_propertydef = rep_propertydef + "    IntProperty numSeedPoint" + help + "_;\n"
                    + "    FloatVec3Property seedPoint" + help + "_;\n\n";
            }
            rep_propertydef = rep_propertydef + "    std::vector<tgt::vec3> seedPoints" + help + ";\n\n";
        }
    }

    // create template and replace the keywords by the generated replacements
    Template t(getTemplatePath("itk_filter_h.txt"));
    t.replace("§UPPERFILTERNAME§", getHeaderGuard());
    t.replace("§FILTERNAME§", name_);
    t.replace("§GROUP§", group);
    t.replace("§MODULENAME§", moduleName.substr(4));
    t.replace("§CODESTATE§", codeState_);
    t.replace("§PROPERTYINCLUDE§", rep_propertyinclude);
    t.replace("§FILTERMETHODS§", rep_filtermethods);
    t.replace("§VOLUMETYPESWITCH§", rep_volumetypeswitch);
    t.replace("§PORTS§", rep_ports);
    t.replace("§PROPERTYDEFINITION§", rep_propertydef);

    std::string filterFilename;
    filterFilename = filterFilename + getFilename() + ".h";
    t.writeToFile(getProcessorOutputPath(moduleName, filterFilename));
}

std::string Filter::getFilename() {
    return toLower(name_);
}

std::string Filter::getHeaderGuard() {
    return toUpper(name_);
}

std::string Filter::getTemplatePath(std::string templateFilename) {
    std::string path;
    path = path + VoreenApplication::app()->getBasePath()
        + "/apps/itk_wrapper/templates/" + templateFilename;
    return path;
}

std::string Filter::getProcessorOutputPath(std::string moduleName, std::string filterFilename) {
    std::string path;
    path = path + VoreenApplication::app()->getBasePath("modules/itk_generated/processors/" + moduleName + "/" + filterFilename);
    return path;
}

std::string Filter::getName() {
    return name_;
}

bool Filter::isEnabled() {
    return enabled_;
}

std::string Filter::getDescription() {
    return description_;
}

bool Filter::isAutoGenerated() {
    return autoGenerated_;
}

std::string Filter::getCodeState() {
    return codeState_;
}

std::vector<Ports> Filter::getInports() {
    return inports_;
}

std::vector<Ports> Filter::getOutports() {
    return outports_;
}

std::vector<Argument> Filter::getArguments() {
    return arguments_;
}

std::string Filter::getFilterMethod(std::string methodType, bool kernel) {
    // replacements depending on input and output of the filter
    std::string rep_templates;
    std::string rep_enableprocessor;
    std::string rep_itkimagetypedef;
    std::string rep_itkpointer;
    std::string rep_itkimagetypes;
    std::string rep_filterinputset;
    std::string rep_outputset;

    // replacements depending on the arguments of the filter
    std::string rep_voxelproperty;
    std::string rep_filterattributeset;
    std::string rep_filterattributeget;

    // helpers
    std::string templateLetters[] = {"T", "S", "R"};
    std::string komma;
    std::string number;
    std::string outputname;

    // create input-part of the replacements depending on input and output
    for (size_t i=0; i<inports_.size(); ++i) {
        komma = "";
        if (i>0)
            komma = ", ";
        number = itos(i+1);

        rep_templates = rep_templates + komma + "class " + templateLetters[i];

        if (inports_[i].isFilterTemplate()) {
            if (methodType == "Vec2") {
                rep_itkimagetypedef = rep_itkimagetypedef
                    + "    typedef itk::Image<itk::CovariantVector<" + templateLetters[i]
                    + ",2>, 3> InputImageType" + number+ ";\n";
            }
            else if (methodType == "Vec3") {
                rep_itkimagetypedef = rep_itkimagetypedef
                    + "    typedef itk::Image<itk::CovariantVector<" + templateLetters[i]
                    + ",3>, 3> InputImageType" + number+ ";\n";
            }
            else if (methodType == "Vec4") {
                rep_itkimagetypedef = rep_itkimagetypedef
                    + "    typedef itk::Image<itk::CovariantVector<" + templateLetters[i]
                    + ",4>, 3> InputImageType" + number+ ";\n";
            }
            else {
                rep_itkimagetypedef = rep_itkimagetypedef + "    typedef itk::Image<"
                    + templateLetters[i] + ", 3> InputImageType" + number+ ";\n";
            }
            rep_itkpointer = rep_itkpointer + "    typename InputImageType" + number
                + "::Pointer p" + number + " = voreen" + methodType + "ToITK" + methodType
                + "<" + templateLetters[i]+ ">(inport"+ number + "_.getData());\n";
            rep_itkimagetypes = rep_itkimagetypes + komma + "InputImageType" + number;
        }
        else {
            rep_itkpointer = rep_itkpointer + "    typename InputImageType1::Pointer p"
                + number + " = voreen" + methodType + "ToITK" + methodType + "<T>(inport"
                + number + "_.getData());\n";
        }

        if (inports_.size() > 1) {
            if (!inports_[i].isNameIsSetter()) {
                if (inports_[1].isNameIsSetter()) {
                    rep_filterinputset = rep_filterinputset
                    + "    filter->SetInput(p" + number + ");\n";
                }
                else {
                    rep_filterinputset = rep_filterinputset
                    + "    filter->SetInput" + number + "(p" + number + ");\n";
                }
            }
            else {
                rep_filterinputset = rep_filterinputset + "    filter->Set"
                    + inports_[i].getName() + "(p" + number + ");\n";
            }
        }
        else {
            if (!inports_[i].isNameIsSetter()) {
                rep_filterinputset = rep_filterinputset + "    filter->SetInput(p"
                    + number + ");\n";
            }
            else {
                rep_filterinputset = rep_filterinputset + "    filter->Set" + inports_[i].getName()
                    + "(p" + number + ");\n";
            }
        }
    }

    // create output-part of the replacements depending on input and output
    for (size_t i=0; i<outports_.size(); ++i) {
        number = itos(i+1);
        if (outports_.size() == 1) {
            outputname = "Output";
        }
        else {
            outputname = outports_[i].getName();
        }

        rep_outputset = rep_outputset + "    Volume* outputVolume" + number + " = 0;\n";
        if (!outports_[i].getPossibleTypes().empty()) {
            std::string volumeType = outports_[i].getPossibleTypes()[0];
            if (volumeType == "2xUInt8" || volumeType == "2xInt8"
                || volumeType == "2xUInt16" || volumeType == "2xInt16"
                || volumeType == "2xFloat" || volumeType == "2xDouble") {

                rep_itkimagetypedef = rep_itkimagetypedef
                    + "    typedef itk::Image<itk::CovariantVector<" + outports_[i].getITKTypename(0)
                    + ",2>, 3> OutputImageType" + number + ";\n";
                rep_outputset = rep_outputset + "    outputVolume" + number
                    + " = ITKVec2ToVoreenVec2Copy<" + outports_[i].getITKTypename(0)
                    + ">(filter->Get" + outputname + "());\n";
            }
            else if (volumeType == "3xUInt8" || volumeType == "3xInt8"
                || volumeType == "3xUInt16" || volumeType == "3xInt16"
                || volumeType == "3xFloat" || volumeType == "3xDouble") {
                rep_itkimagetypedef = rep_itkimagetypedef
                    + "    typedef itk::Image<itk::CovariantVector<" + outports_[i].getITKTypename(0)
                    + ",3>, 3> OutputImageType" + number + ";\n";
                rep_outputset = rep_outputset + "    outputVolume" + number
                    + " = ITKVec3ToVoreenVec3Copy<" + outports_[i].getITKTypename(0)
                    + ">(filter->Get" + outputname + "());\n";
            }
            else if (volumeType == "4xUInt8" || volumeType == "4xInt8"
                || volumeType == "4xUInt16" || volumeType == "4xInt16"
                || volumeType == "4xFloat" || volumeType == "4xDouble") {
                rep_itkimagetypedef = rep_itkimagetypedef
                    + "    typedef itk::Image<itk::CovariantVector<"
                    + outports_[i].getITKTypename(0) + ",4>, 3> OutputImageType" + number + ";\n";
                rep_outputset = rep_outputset + "    outputVolume" + number
                    + " = ITKVec4ToVoreenVec4Copy<" + outports_[i].getITKTypename(0)
                    + ">(filter->Get" + outputname + "());\n";
            }
            else {
                rep_itkimagetypedef = rep_itkimagetypedef + "    typedef itk::Image<"
                    + outports_[i].getITKTypename(0) + ", 3> OutputImageType" + number + ";\n";
                rep_outputset = rep_outputset + "    outputVolume" + number
                    + " = ITKToVoreenCopy<" + outports_[i].getITKTypename(0) +">(filter->Get"
                    + outputname + "());\n";
            }
        }
        else {
            if (methodType == "Vec2") {
                rep_itkimagetypedef = rep_itkimagetypedef
                    + "    typedef itk::Image<itk::CovariantVector<T,2>, 3> OutputImageType"
                    + number+ ";\n";
            }
            else if (methodType == "Vec3") {
                rep_itkimagetypedef = rep_itkimagetypedef
                    + "    typedef itk::Image<itk::CovariantVector<T,3>, 3> OutputImageType"
                    + number+ ";\n";
            }
            else if (methodType == "Vec4") {
                rep_itkimagetypedef = rep_itkimagetypedef
                    + "    typedef itk::Image<itk::CovariantVector<T,4>, 3> OutputImageType"
                    + number+ ";\n";
            }
            else {
                rep_itkimagetypedef = rep_itkimagetypedef
                    + "    typedef itk::Image<T, 3> OutputImageType" + number +";\n";
            }
            rep_outputset = rep_outputset + "    outputVolume" + number + " = ITK"
                + methodType + "ToVoreen" + methodType + "Copy<T>(filter->Get" + outputname + "());\n";
        }

        rep_outputset = rep_outputset + "\n    if (outputVolume" + number + ") {\n";

        if(outports_[i].transferRWM())
            rep_outputset += "        transferRWM(inport1_.getData(), outputVolume" + number + ");\n";

        if(outports_[i].transferTransformation())
            rep_outputset += "        transferTransformation(inport1_.getData(), outputVolume" + number + ");\n";

        rep_outputset += "        outport" + number + "_.setData(outputVolume" + number + ");\n";
        rep_outputset += "    } else\n";
        rep_outputset += "        outport" + number + "_.setData(0);\n\n";

        if (outports_[i].isFilterTemplate()) {
            rep_itkimagetypes = rep_itkimagetypes + ", OutputImageType" + number;
        }
    }

    Argument argument;

    // create replacements depending on the arguments of the filter
    for(size_t j=0; j<arguments_.size(); ++j) {
        argument = arguments_[j];

        if (argument.isInputArgument()) {
            rep_filterattributeset = rep_filterattributeset + argument.getFilterSetOrGet();
        }
        else {
            rep_filterattributeget = rep_filterattributeget + argument.getFilterSetOrGet();
        }

        if (argument.getType() == "PixelType") {
            rep_voxelproperty = rep_voxelproperty + "    " + argument.getPropertyname()
                + ".setVolume(inport1_.getData());\n";
            if (!argument.getMinValue().empty()) {
                rep_voxelproperty = rep_voxelproperty + "    " + argument.getPropertyname()
                    + ".setMinValue<T>(" + argument.getMinValue() + ");\n";
            }
            if (!argument.getMaxValue().empty()) {
                rep_voxelproperty = rep_voxelproperty + "    " + argument.getPropertyname()
                    + ".setMaxValue<T>(" + argument.getMaxValue() + ");\n";
            }
        }
    }

    if (inports_.size() == 1 && outports_.size() == 1) {
        rep_enableprocessor = rep_enableprocessor + "    if (!enableProcessing_.get()) {\n"
            + "        outport1_.setData(inport1_.getData(), false);\n"
            + "        return;\n"
            + "    }\n";
    }

    // create template and replace the keywords by the generated replacements
    std::string templateFilename;
    if (kernel) {
        templateFilename = getTemplatePath("filterMethodKERNEL.txt");
    }
    else {
        templateFilename = getTemplatePath("filterMethod.txt");
    }

    Template t(templateFilename);
    t.replace("§FILTERNAME§", name_);
    t.replace("§FIRSTLOWFILTERNAME§", getProcessFunctionName());
    t.replace("§METHODTYPE§", methodType);

    t.replace("§TEMPLATES§", rep_templates);
    t.replace("§ENABLEPROCESSOR§", rep_enableprocessor);
    t.replace("§ITKIMAGETYPEDEF§", rep_itkimagetypedef);
    t.replace("§ITKPOINTER§", rep_itkpointer);
    t.replace("§ITKIMAGETYPES§", rep_itkimagetypes);
    t.replace("§FILTERINPUTSET§", rep_filterinputset);
    t.replace("§OUTPUTSET§", rep_outputset);

    t.replace("§VOXELPROPERTY§", rep_voxelproperty);
    t.replace("§FILTERATTRIBUTESET§", rep_filterattributeset);
    t.replace("§FILTERATTRIBUTEGET§", rep_filterattributeget);

    std::string filterMethod = t.getTemplateString();
    filterMethod = filterMethod + "\n\n";

    return filterMethod;
}

std::string Filter::getProcessingMethods() {
    std::string processingMethods;

    // replacements
    std::string rep_templated;
    std::string rep_methodname;
    std::string rep_number;
    std::string rep_volumetypecases;

    // helpers
    std::string templates;

    // create replacements
    for (size_t i=0; i<inports_.size(); ++i) {
        rep_number = itos(i+1);
        rep_templated = "";
        rep_methodname = "";
        templates = "";

        if (i==0) {
            rep_methodname = "process";
        }
        else {
            if (i == 1) {
                rep_templated = "template <class T>";
                rep_methodname = "volumeTypeSwitch1";
                templates = "T, ";
            }
            else {
                rep_templated = "template <class T, class S>";
                rep_methodname = "volumeTypeSwitch2";
                templates = "T, S, ";
            }
        }

        rep_volumetypecases = "    ";
        if (!inports_[i].isFilterTemplate()) {
            if (i == (inports_.size()-1)) {
                rep_volumetypecases = rep_volumetypecases
                    + "if (!dynamic_cast<const VolumeAtomic<T>*>(inputVolume" + rep_number
                    + "))  {\n"
                    + "        LERROR(\"Volume " + rep_number
                    + " must be from the same type as Volume 1!\");\n"
                    + "    }\n";

                std::string volumeType;
                std::string methodType;
                for (size_t j=0; j<inports_[0].getPossibleTypes().size(); ++j) {
                    volumeType = inports_[0].getPossibleTypes()[j];
                    methodType = "";
                    if (volumeType == "2xUInt8" || volumeType == "2xInt8"
                        || volumeType == "2xUInt16" || volumeType == "2xInt16"
                        || volumeType == "2xFloat" || volumeType == "2xDouble") {
                        methodType = "Vec2";
                    }
                    if (volumeType == "3xUInt8" || volumeType == "3xInt8"
                        || volumeType == "3xUInt16" || volumeType == "3xInt16"
                        || volumeType == "3xFloat" || volumeType == "3xDouble") {
                        methodType = "Vec3";
                    }
                    if (volumeType == "4xUInt8" || volumeType == "4xInt8"
                        || volumeType == "4xUInt16" || volumeType == "4xInt16"
                        || volumeType == "4xFloat" || volumeType == "4xDouble") {
                        methodType = "Vec4";
                    }

                    rep_volumetypecases = rep_volumetypecases + "else if (dynamic_cast<const VolumeRAM_"
                        + inports_[i].getPossibleTypes()[j] + "*>(inputVolume" + rep_number + "))  {\n"
                        + "        " + getProcessFunctionName() + methodType + "ITK<"
                        + templates + inports_[i].getITKTypename(static_cast<int>(j)) + ">();\n    }\n";
                }
            }
            else {
                rep_volumetypecases = rep_volumetypecases
                    + "if (!dynamic_cast<const VolumeAtomic<T>*>(inputVolume" + rep_number
                    + "))  {\n"
                    + "        LERROR(\"Volume " + rep_number
                    + " must be from the same type as Volume 1!\");\n"
                    + "    }\n    else {\n"
                    + "        volumeTypeSwitch" + rep_number + "<" + templates + "T>();\n"
                    + "    }\n";
            }
        }
        else {
            if (i == (inports_.size()-1)) {
                std::string volumeType;
                std::string methodType;
                for (size_t j=0; j<inports_[i].getPossibleTypes().size(); ++j) {
                    volumeType = inports_[i].getPossibleTypes()[j];
                    methodType = "";
                    if (volumeType == "2xUInt8" || volumeType == "2xInt8"
                        || volumeType == "2xUInt16" || volumeType == "2xInt16"
                        || volumeType == "2xFloat" || volumeType == "2xDouble") {
                        methodType = "Vec2";
                    }
                    if (volumeType == "3xUInt8" || volumeType == "3xInt8"
                        || volumeType == "3xUInt16" || volumeType == "3xInt16"
                        || volumeType == "3xFloat" || volumeType == "3xDouble") {
                        methodType = "Vec3";
                    }
                    if (volumeType == "4xUInt8" || volumeType == "4xInt8"
                        || volumeType == "4xUInt16" || volumeType == "4xInt16"
                        || volumeType == "4xFloat" || volumeType == "4xDouble") {
                        methodType = "Vec4";
                    }

                    rep_volumetypecases = rep_volumetypecases + "if (dynamic_cast<const VolumeRAM_"
                        + inports_[i].getPossibleTypes()[j] + "*>(inputVolume" + rep_number + "))  {\n"
                        + "        " + getProcessFunctionName() + methodType + "ITK<"
                        + templates + inports_[i].getITKTypename(static_cast<int>(j)) + ">();\n"
                        + "    }\n    else ";
                }

                rep_volumetypecases = rep_volumetypecases + "{\n"
                    + "        LERROR(\"Inputformat of Volume " + rep_number
                    + " is not supported!\");\n"
                    + "    }";
            }
            else {
                for (size_t j=0; j<inports_[i].getPossibleTypes().size(); ++j){
                    rep_volumetypecases = rep_volumetypecases + "if (dynamic_cast<const VolumeRAM_"
                        + inports_[i].getPossibleTypes()[j] + "*>(inputVolume" + rep_number + "))  {\n"
                        + "        volumeTypeSwitch" + rep_number + "<" + templates
                        + inports_[i].getITKTypename(static_cast<int>(j)) + ">();\n"
                        + "    }\n    else ";
                }

                rep_volumetypecases = rep_volumetypecases + "{\n"
                    + "        LERROR(\"Inputformat of Volume " + rep_number
                    + " is not supported!\");\n"
                    + "    }";
            }
        }

        // create template and replace the keywords by the generated replacements
        Template processMethod(getTemplatePath("processBase.txt"));
        processMethod.replace("§TEMPLATED§", rep_templated);
        processMethod.replace("§FILTERNAME§", name_);
        processMethod.replace("§METHODNAME§", rep_methodname);
        processMethod.replace("§NUMBER§", rep_number);
        processMethod.replace("§VOLUMETYPECASES§", rep_volumetypecases);

        processingMethods = processingMethods + processMethod.getTemplateString() + "\n\n" ;

    }
    return processingMethods;
}

std::string Filter::getProcessFunctionName() {
    std::string processFunctionName = name_;
    processFunctionName[0] = tolower(processFunctionName[0]);
    return processFunctionName;
}

/***********************************************************************************
 * Class ITK_Module
 ***********************************************************************************/
void ITK_Module::serialize(Serializer& s) const {
    s.serialize("name", name_);
    s.serialize("group", group_);
    s.serialize("filterlist", filters_, "filter");
}

void ITK_Module::deserialize(Deserializer& s) {
    s.deserialize("name", name_);
    s.deserialize("group", group_);
    s.deserialize("filterlist", filters_, "filter");
}

std::string ITK_Module::getName() {
    return name_;
}

std::string ITK_Module::getGroup() {
    return group_;
}

std::vector<Filter> ITK_Module::getFilters() {
    return filters_;
}

