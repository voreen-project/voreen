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

#include "dicomdict.h"

namespace voreen {

DicomDict::DicomDict()
{
}

void DicomDict::serialize(Serializer &s) const {
    //put DictEntries into a Vector
    std::vector<DicomDictEntry> entryVector;

    std::map<std::string,DicomDictEntry>::const_iterator it;
    for (it = entries_.begin(); it != entries_.end(); ++it) {
        entryVector.push_back(it->second);
    }

    //serialize this vector
    s.serialize("Dict", entryVector, "entry");
}

void DicomDict::deserialize(Deserializer &d) {
    //get Vector of DictEntries
    std::vector<DicomDictEntry> entryVector;

    //deserialize Vector
    //d.setUseAttributes(false);
    d.deserialize("Dict", entryVector, "entry");

    //iterate over Vector and put in map
    std::vector<DicomDictEntry>::iterator it;
    for (it = entryVector.begin(); it != entryVector.end(); ++it) {
        //check if keyword is already in entries_
        if (entries_.find(it->getKeyword()) != entries_.end()) {
            throw tgt::Exception("Keyword in Dictionary must be unique! Found keyword more than once: " + it->getKeyword());
        }

        entries_[it->getKeyword()] = *it;
    }
}

const DicomDictEntry DicomDict::getDictEntryByKeyword(const std::string &keyword) const {
    if (entries_.find(keyword) == entries_.end())
        throw tgt::FileException("Keyword could not be found in Dictionary: "+ keyword);
    else
        return entries_.find(keyword)->second;
}

DicomDict* DicomDict::loadFromFile(const std::string &fileName) {
    //Check the given file
    if (!(tgt::FileSystem::fileExists(fileName)))
            throw tgt::FileNotFoundException("Dictionary not found", fileName);

    if (!(tgt::FileSystem::fileExtension(fileName) == "xml"))
            throw tgt::FileAccessException("Dictionary is not an XML file", fileName);

    //get Data from the file
    std::string xmlData;
    tgt::FileSystem sys;
    tgt::File* file = sys.open(fileName);

    if (!file->good())
        throw tgt::FileAccessException("Dictionary file cannot be accessed", fileName);

    xmlData = file->getAsString();

    if (file->isOpen())
        file->close();

    delete file;
    file = 0;

    DicomDict* dict = new DicomDict();

    //Put Data into a Stream to be read by XmlDeserializer and try to deserialize it
    std::stringstream stream;
    stream << xmlData;

    XmlDeserializer d;
    d.setUseAttributes(true);
    d.read(stream);

    try {
        Deserializer deserializer(d);
        dict->deserialize(deserializer);
    }
    catch (tgt::Exception e) {
        throw tgt::FileAccessException("Could not load Dictionary: " + std::string(e.what()),fileName);
    }

    return dict;
}

bool DicomDict::containsKeyword(const std::string &key) const {
    return (entries_.count(key) == 1);
}

const std::vector<std::string> DicomDict::getKeywordVector() const {
    std::vector<std::string> v;

    std::map<std::string,DicomDictEntry>::const_iterator it;
    for (it = entries_.begin(); it != entries_.end(); ++it) {
        v.push_back(it->first);
    }

    return v;
}

} //namespace


