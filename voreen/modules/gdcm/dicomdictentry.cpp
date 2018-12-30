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

#include "dicomdictentry.h"

namespace voreen {

using std::string;

DicomDictEntry::DicomDictEntry()
{
    groupNumber_ = "0";
    elementNumber_ = "0";
    keyword_ = "none";
    name_ = "none";
    vm_ = "0";
    vr_ = "none";
    metaData_ = false;
}

DicomDictEntry::DicomDictEntry(string keyword, string name, string vm, string vr, bool meta){
    keyword_ = keyword;
    name_ = name;
    vm_ = vm;
    vr_ = vr;
    metaData_ = meta;
}

void DicomDictEntry::serialize(Serializer &s) const {
    s.serialize("group", groupNumber_);
    s.serialize("element", elementNumber_);
    s.serialize("keyword", keyword_);
    s.serialize("name", name_);
    s.serialize("vr", vr_);
    s.serialize("vm", vm_);
    s.serialize("metaData", metaData_);
}

void DicomDictEntry::deserialize(Deserializer &s) {
    s.deserialize("group", groupNumber_);
    s.deserialize("element", elementNumber_);
    s.deserialize("keyword", keyword_);
    s.deserialize("name", name_);
    s.deserialize("vr", vr_);
    s.deserialize("vm", vm_);
    s.deserialize("metaData", metaData_);
}


string DicomDictEntry::getKeyword() const{
    return keyword_;
}

string DicomDictEntry::getValueRepresentation() const{
    return vr_;
}

string DicomDictEntry::getValueMultiplicity()const{
    return vm_;
}

string DicomDictEntry::getName() const{
    return name_;
}

string DicomDictEntry::getGroupNumber() const{
    return groupNumber_;
}

string DicomDictEntry::getElementNumber() const{
    return elementNumber_;
}

uint16_t DicomDictEntry::getElementTagUint16() const{
    std::stringstream n;
    n << std::hex << elementNumber_;
    uint16_t group;
    n >> group;
    return group;
}

uint16_t DicomDictEntry::getGroupTagUint16() const{
    std::stringstream n;
    n << std::hex << groupNumber_;
    uint16_t elem;
    n >> elem;
    return elem;
}

bool DicomDictEntry::isMetaData() const {
    return metaData_;
}

void DicomDictEntry::setGroupNumber(string group) {
        //check if length is 4
        if (!(group.length() == 4)) {
            groupNumber_ = "0000";
            return;
        }

        //and if the string represents a hexadecimal value and may be converted into an unsigned int
        std::stringstream n;
        n << std::hex << group;
        uint16_t gnum;
        if (n >> gnum)
            groupNumber_ = group;
        else
            groupNumber_ = "0000";

}

void DicomDictEntry::setElementNumber(string element) {
        //check if length is 4
        if (!(element.length() == 4)) {
            elementNumber_ = "0000";
            return;
        }

        //and if the string represents a hexadecimal value and may be converted into an unsigned int
        std::stringstream n;
        n << std::hex << element;
        uint16_t elem;
        if (n >> elem)
            elementNumber_ = element;
        else
            elementNumber_ = "0000";

}

string DicomDictEntry::toString() const {
    string representation = "(" + groupNumber_ + "|" + elementNumber_ +")";
    representation.append("; Name: "+name_);
    representation.append("; Keyword: " + keyword_);
    representation.append("; VR: " + vr_);
    representation.append("; VM: " + vm_);
    representation.append("; MetaData: " + std::to_string(metaData_));

    return representation;
}

} //end namespace
