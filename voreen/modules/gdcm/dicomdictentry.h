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

#ifndef VRN_DICOMDICTENTRY_H
#define VRN_DICOMDICTENTRY_H

// include this before any GDCM header in order to make sure that C99 types are defined
#include "tgt/types.h"

#include "voreen/core/io/serialization/serializable.h"
#include "voreen/core/io/serialization/xmlserializer.h"
#include "voreen/core/io/serialization/xmldeserializer.h"
#include <string>

namespace voreen {

/**
 * Class to represent an entry in a Dicom Dictionary
 *
 * Used by GdcmVolumeReader to be loaded from XML files.
 * Includes information about Group Number, Element Number, Keyword, Name, Value Representation and Value Multiplicity of the entry
 * as well as an attribute to represent if this Entry is essentially required for the Dictionary
 */
class VRN_CORE_API DicomDictEntry : public Serializable
{
public:
    DicomDictEntry();

    DicomDictEntry(std::string keyword, std::string name, std::string vm, std::string vr, bool isMetaData);

    /**
     * Serializes all member variables using the given XmlSerializer
     *
     * @param s the XmlSerializer to be used
     */
    virtual void serialize(Serializer& s) const;

    /**
     * Deserializes all member variables using the given XmlDeserializer
     *
     * @param s the XmlDeserializer to be used
     */
    virtual void deserialize(Deserializer& s);

    /**
     * Sets the group number.
     *
     * @param group must be of the following form: gggg (where g is a hexadecimal digit), otherwise it is set to 0
     */
    void setGroupNumber(std::string group);

    /**
     * Sets the element number.
     *
     * @param element must be of the following form: eeee (where e is a hexadecimal digit), otherwise it is set to 0
     */
    void setElementNumber(std::string element);

    /**
     * Returns a std::string representation of the DictEntry
     */
    std::string toString() const;

    /**
     * get the keyword as a string
     */
    std::string getKeyword() const;

    /**
     * get VR as string
     */
    std::string getValueRepresentation() const;

    /**
     * get VM as string
     */
    std::string getValueMultiplicity() const;

    /**
     * get name as string
     */
    std::string getName() const;

    /**
     * get group number as (hex-)string
     */
    std::string getGroupNumber() const;

    /**
     * get element number as (hex-)string
     */
    std::string getElementNumber() const;

    /**
     * sould the value of this Entry be added to the Meta Information of the Volume?
     */
    bool isMetaData() const;

    /**
     * get group number as uint16_t
     */
    uint16_t getGroupTagUint16() const;

    /**
     * get group number as uint16_t
     */
    uint16_t getElementTagUint16() const;

private:
    std::string keyword_; ///< The keyword to find the entry
    std::string vm_;      ///< The Value Multiplicity
    std::string vr_;      ///< The Value Representation
    std::string name_;    ///< The Name of the Entry
    std::string groupNumber_; ///< The Tag Group Number
    std::string elementNumber_; ///< The Tag Element Number
    bool metaData_; ///< implies if this entry (ie. it's value in the loaded Files) should be added to the Volume's MetaDataContainer
};

} //end namespace

#endif // DICOMDICTENTRY_H
