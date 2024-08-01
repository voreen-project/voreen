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

#ifndef VRN_CUSTOMDICOMDICT_H
#define VRN_CUSTOMDICOMDICT_H

#include "dicomdict.h"

namespace voreen {

/**
 * This class represents an extended Dictionary. In addition to the DicomDict itself it contains conditions, for which DICOM files it is suitable,
 * which are represented by a list of Keywords for a Tag and the possible values for these DICOM Tags.
 * Furthermore, a list of subdivisionKeywords is included. The DICOM files of one SeriesInstanceUID should additionally be subdivided into
 * several Volumes by these Tags (e.g. for DTI DICOM Data Sets when one Series has to be subdivided by Gradient Diffusion Orientation).
 * The DicomDict of a CustomDicomDict necessarily has to contain ALL of the Keywords used in this CustomDicomDict, even if they are part of the Standard Dictionary!
 * (Although they may contain more Keywords for additional information)
 */
class VRN_CORE_API CustomDicomDict : public Serializable
{
public:
    CustomDicomDict();

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
     * Loads a CustomDicomDict from an XML file.
     * If the file is corrupt, can not be found, or another error occurs, a tgt::FileException is thrown
     *
     * @param fileName the file to be loaded
     */
    void loadFromFile(const std::string &fileName);

    /**
     * returns a pointer to the DicomDict within this CustomDicomDict
     */
    const DicomDict* getDict() const;

    /**
     * returns a pointer to the list of keywords for additional subdivision of one SeriesInstanceUID
     */
    const std::vector<std::string>* getSubdivisionKeywords() const;

    /**
     * get a pointer to the keyword / value possibility combinations
     */
    const std::vector<std::pair<std::string, std::vector<std::string> > >* getConditions() const;

private:
    DicomDict dict_; ///< the Dictionary
    std::vector<std::string> subdivisionKeywords_; ///< std::vector of keywords by which slices belonging to one SeriesInstanceUID should additionally be subdivided
    std::vector<std::pair<std::string, std::vector<std::string> > > conditions_; ///< keywords and corresponding value possibilities as conditions Dicom files should fulfill to use this CustomDicomDict... if the value possibilities-vector is empty, the attribute is required to be present in the files, but no specific value is required
};


} //namespace

#endif // CUSTOMDICOMDICT_H
