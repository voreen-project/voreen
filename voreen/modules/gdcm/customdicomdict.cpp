/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "customdicomdict.h"

namespace voreen {

CustomDicomDict::CustomDicomDict()
{
}

void CustomDicomDict::serialize(Serializer &s) const {
    s.serialize("conditions", conditions_);
    s.serialize("divisionkeys", subdivisionKeywords_);
    s.serialize("dict",dict_);
}

void CustomDicomDict::deserialize(Deserializer& s) {
    s.deserialize("conditions", conditions_);
    s.deserialize("divisionkeys", subdivisionKeywords_);
    s.deserialize("dict", dict_);
}

void CustomDicomDict::loadFromFile(const std::string &fileName) {
    //Check the given file
    if (!(tgt::FileSystem::fileExists(fileName)))
            throw tgt::FileNotFoundException("SpecialDictionary not found", fileName);

    if (!(tgt::FileSystem::fileExtension(fileName) == "xml"))
            throw tgt::FileAccessException("SpecialDictionary is not an XML file", fileName);

    //get Data from the file
    std::string xmlData;
    tgt::FileSystem sys;
    tgt::File* file = sys.open(fileName);

    if (!file->good())
        throw tgt::FileAccessException("SpecialDictionary file cannot be accessed", fileName);

    xmlData = file->getAsString();

    if (file->isOpen())
        file->close();

    delete file;
    file = 0;

    //Put Data into a Stream to be read by XmlDeserializer and try to deserialize it
    std::stringstream stream;
    stream << xmlData;

    XmlDeserializer d;
    d.setUseAttributes(true);
    d.read(stream);

    try {
        Deserializer deserializer(d);
        deserialize(deserializer);
    }
    catch (tgt::Exception e) {
        throw tgt::FileAccessException("Could not load Dictionary: " + std::string(e.what()),fileName);
    }
}

const DicomDict* CustomDicomDict::getDict() const {
    return &dict_;
}

const std::vector<std::string>* CustomDicomDict::getSubdivisionKeywords() const{
    return &subdivisionKeywords_;
}

const std::vector<std::pair<std::string, std::vector<std::string> > >* CustomDicomDict::getConditions() const{
    return &conditions_;
}

} //namespace

