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

#ifndef VRN_DICOMDIRPARSER_H
#define VRN_DICOMDIRPARSER_H

#include "voreen/core/voreencoreapi.h"

// include this before any GDCM header in order to make sure that C99 types are defined
#include "tgt/types.h"
#include "tgt/exception.h"

#include <gdcmReader.h>
#include <gdcmStringFilter.h>

namespace voreen {

/**
 * Provides functionality to parse the Items in a DirectoryRecordSequence, identify their type and extract values.
 */
class VRN_CORE_API DirectoryRecordSequenceIterator {

public:
    DirectoryRecordSequenceIterator();

    /**
     * Sets the DataElement. If the DataElement is not a DirectoryRecordSequence or a DirectoryRecordSequence that contains 0 items, a tgt::FileException is thrown.
     * If the DataElement is successfully set, the iterator points is set to the first item in the DirectoryRecordSequence.
     */
    virtual void setDataElement(gdcm::DataElement data);

    /**
     * Sets iterator to the first item in the DirectoryRecordSequence
     */
    virtual void setIteratorToFirstElement() const;

    /**
     * Sets the iterator to the next item in the DirectoryRecordSequence and returns true.
     * If there is no item left, the iterator will still be set (trying to read the item afterwards results in a tgt::FileException).
     */
    virtual void setIteratorToNextElement() const;

    /**
     * Check is at the end of the DirectoryRecordSequence so that no item is left
     */
    virtual bool iteratorIsAtEnd() const;

    /**
     * returns the number of items in the DirectoryRecordSequence
     */
    virtual int getNumberOfItems() const;

    /**
     * get the name of the item the iterator is presently set to.
     */
    virtual std::string getPresentItemName() const;

    /**
     * check if a tag is present in the current item
     */
    virtual bool findTagInPresentItem(gdcm::Tag t) const;

    /**
     * Get the value of the current item for the given Tag.
     * If there is no value present for this Tag, an empty string will be returned.
     */
    virtual std::string getElementValue(gdcm::Tag t) const;

private:
    gdcm::DataElement data_; ///< the actual DirectoryRecordSequence
    gdcm::SmartPointer<gdcm::SequenceOfItems> sqi_; ///< gdcm::SmartPointer for parsing the DirectoryRecordSequence
    mutable unsigned int item_; ///< the item currently pointed to
};



/**
 * DicomDirParser uses functionality of GDCM library to open a DICOMDIR file
 * and parse its content in an abstract and simple way.
 */
class DicomDirParser
{

public:
    DicomDirParser();

    /**
     * Constructor that also sets filename, but does not open the file.
     * => Method open() still has to be called explicitly.
     */
    DicomDirParser(std::string filename);

    void setFilename(std::string filename);

    std::string getFilename() const;

    /**
     * Opens the specified DICOMDIR.
     * If no filename is set, if the file is not a DICOMDIR or could not be found, a tgt::FileException is thrown
     */
    virtual void open();

    virtual size_t getNumberOfDirectoryRecordSequences() const;

    /**
     * Returns the DirectoryRecordSequenceIterator of the DirectoryRecordSequence the iterator is presently referencing.
     * If the iterator is at the end or there is no DICOMDIR loaded, a tgt::FileAccessException is thrown.
     */
    virtual DirectoryRecordSequenceIterator getPresentDirectoryRecordSequence() const;

    /**
     * Check if there is no DirectoryRecordSequence left
     */
    virtual bool iteratorIsAtEnd() const;

    /**
     * Used to parse the content of a DICOMDIR:
     * DICOMDIR -> DirectoryRecordSequence (-> Items).
     * If there is no next DirectoryRecordSequence left, the iterator will still be set (and the method iteratorIsAtEnd() will return true).
     * Trying to get the present item will then result in an Exception.
     * If there has not already been opened a DICOMDIR correctly or no DirectoryRecordSequendes are available, a tgt::FileAccessException is thrown.
     */
    virtual void setIteratorToNextDirectoryRecordSequence() const;

    /**
     * Sets the Iterator to the first DirectoryRecordSequence found in the DICOMDIR, so that the structure might newly be parsed.
     * If there is no DirectoryRecordSequence found or if there has not already been opened a DICOMDIR correctly, a tgt::FileAccessException is thrown.
     */
    virtual void setIteratorToFirstDirectoryRecordSequence() const;

private:
    std::string filename_; ///< the filename (including the path to it) of the file to be opened / already opened
    gdcm::DataSet dataSet_; ///< the DataSet extracted from the opened file
    std::vector<gdcm::DataElement> directoryRecordSequences_; ///< pointer to the DataElements in the opened DICOMDIR that are a DirectoryRecordSequence
    mutable std::vector<gdcm::DataElement>::const_iterator it_; ///< iterator for directoryRecordSequences_
};

} //namespace

#endif // DICOMDIRPARSER_H
