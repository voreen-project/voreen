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

#ifndef VRN_DICOMDICT_H
#define VRN_DICOMDICT_H

#include "voreen/core/io/serialization/serializable.h"
#include "voreen/core/io/serialization/xmlserializer.h"
#include "voreen/core/io/serialization/xmldeserializer.h"
#include "tgt/filesystem.h"
#include "tgt/exception.h"
#include "./dicomdictentry.h"
#include <string>
#include <map>

namespace voreen {

/**
 * Class to represent a Dicom Dictionary, containing DicomDictEntry-Objects
 */
class VRN_CORE_API DicomDict : public Serializable
{
public:
    DicomDict();

    /**
     * Serializes all member variables using the given XmlSerializer
     *
     * @param s the XmlSerializer to be used
     */
    virtual void serialize(Serializer& s) const;

    /**
     * Deserializes all member variables using the given XmlDeserializer.
     * If a keyword occurs more than once, a tgt::Exception is thrown
     *
     * @param s the XmlDeserializer to be used
     */
    virtual void deserialize(Deserializer& s);

    /**
     * Returns the DicomDictEntry associated with the keyword.
     * If no entry with this keyword can be found, a tgt::FileException will be thrown
     *
     * @param keyword the keyword by which the entry is identified
     */
    const DicomDictEntry getDictEntryByKeyword(const std::string &keyword) const;

    /**
     * Loads a DicomDict from an XML file and returns a pointer to it.
     * If the file is corrupt, can not be found, or a keyword occurs more than once, a tgt::FileException is thrown
     *
     * @param fileName the file to be loaded
     */
    static DicomDict* loadFromFile(const std::string &fileName);

    /**
     * returns true if the Dict contains an entry with the given keyword
     *
     * @param key the keyword to be lokked for
     */
    bool containsKeyword(const std::string &key) const;

    /**
     * returns a std::vector containing all keywords in the DicomDict
     */
    const std::vector<std::string> getKeywordVector() const;

private:
    std::map<std::string, DicomDictEntry> entries_; ///< the actual map containing the entries with their corresponding keywords as keys
};

} //namespace

#endif // DICOMDICT_H
