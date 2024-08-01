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

#include "dicomdirparser.h"

#include "tgt/filesystem.h"

namespace voreen {

/*--------------------------------
  DirectoryRecordSequence
---------------------------------*/

DirectoryRecordSequenceIterator::DirectoryRecordSequenceIterator() {
    item_ = 1;
}

void DirectoryRecordSequenceIterator::setDataElement(gdcm::DataElement data) {
    //check if DataSet is DirectoryRecordSequence
    if (!(data.GetTag() == gdcm::Tag(0x0004,0x1220)))
        throw tgt::FileException("DataSet is not a DirectoryRecordSequence!");

    //check for length
    gdcm::SmartPointer<gdcm::SequenceOfItems> sqi = data.GetValueAsSQ();

    if (sqi->GetNumberOfItems() == 0)
        throw tgt::FileException("DirectoryRecordSequence contains 0 items!");

    data_ = data;
    sqi_ = data_.GetValueAsSQ();
    item_ = 1;
}

void DirectoryRecordSequenceIterator::setIteratorToFirstElement() const {
    item_ = 1;
}

void DirectoryRecordSequenceIterator::setIteratorToNextElement() const {
    if (item_ < (sqi_->GetNumberOfItems()+1))
        item_++;
}

bool DirectoryRecordSequenceIterator::iteratorIsAtEnd() const {
    return (item_ > sqi_->GetNumberOfItems());
}

int DirectoryRecordSequenceIterator::getNumberOfItems() const {
    return static_cast<int>(sqi_->GetNumberOfItems());
}

std::string DirectoryRecordSequenceIterator::getPresentItemName() const {
    std::stringstream strm;

    if (sqi_->GetItem(item_).FindDataElement(gdcm::Tag(0x0004,0x1430)))
      sqi_->GetItem(item_).GetDataElement(gdcm::Tag(0x0004,0x1430)).GetValue().Print(strm);

    std::string itemName(strm.str());
    return itemName;
}

bool DirectoryRecordSequenceIterator::findTagInPresentItem(gdcm::Tag t) const {
    return sqi_->GetItem(item_).FindDataElement(t);
}

std::string DirectoryRecordSequenceIterator::getElementValue(gdcm::Tag t) const {
    std::stringstream strm;

    const gdcm::DataElement& elem = sqi_->GetItem(item_).GetDataElement(t);
    if (!elem.IsEmpty()) {
        elem.GetValue().Print(strm);
    }

    return strm.str();
}


/*--------------------------------
  DicomDirParser
---------------------------------*/

DicomDirParser::DicomDirParser(){}

DicomDirParser::DicomDirParser(std::string filename) {
    filename_ = filename;
}

void DicomDirParser::setFilename(std::string filename) {
    filename_ = filename;
}

std::string DicomDirParser::getFilename() const {
    return filename_;
}

void DicomDirParser::open() {

    directoryRecordSequences_.clear(); //clear old Data

    //check if filename_ is set
    if (filename_.empty()) {
        it_ = directoryRecordSequences_.begin();
        throw tgt::FileException("No filename specified!");
    }


    //check if file could be opened as DICOM file
    gdcm::Reader reader;
    reader.SetFileName(filename_.c_str());
    if (!reader.Read()) {
        filename_.clear();
        it_ = directoryRecordSequences_.begin();
        throw tgt::FileException("File does not exist or is not a DICOM file: " + filename_, filename_);
    }


    //check if file is DICOMDIR
    gdcm::MediaStorage ms;
    ms.SetFromFile(reader.GetFile());

    if(!(ms == gdcm::MediaStorage::MediaStorageDirectoryStorage)) {
        filename_.clear();
        it_ = directoryRecordSequences_.begin();
        throw tgt::FileException("File is not a DICOMDIR!", filename_);
    }


    gdcm::StringFilter sf;
    sf.SetFile(reader.GetFile());

    //check Media Storage Sop Class UID
    std::string sop = sf.ToString(gdcm::Tag(0x0002, 0x0002));
    if (sop != "1.2.840.10008.1.3.10") {
        filename_.clear();
        it_ = directoryRecordSequences_.begin();
        throw tgt::FileException("DICOMDIR has wrong Media Storage SOP Class: " + sop, filename_);
    }

    //now get the actual Data
    dataSet_ = reader.GetFile().GetDataSet();

    std::set<gdcm::DataElement>::iterator it;

    //get all DataElements with DirectoryRecordSequence and add pointers to directoryRecordSequences_
    for(it = dataSet_.GetDES().begin() ; it != dataSet_.GetDES().end(); ++it){

        //check if iterator points to a DirectoryRecordSequence
        if (it->GetTag()==gdcm::Tag(0x0004,0x1220)){
            gdcm::DataElement element = *it;
            directoryRecordSequences_.push_back(element);
        }
    }

    it_ = directoryRecordSequences_.begin();
}

size_t DicomDirParser::getNumberOfDirectoryRecordSequences() const {
    return directoryRecordSequences_.size();
}

bool DicomDirParser::iteratorIsAtEnd() const {
    return (it_ == directoryRecordSequences_.end());
}

void DicomDirParser::setIteratorToNextDirectoryRecordSequence() const {
    if (filename_.empty())
        throw tgt::FileAccessException("No file specified to access DirectoryRecordSequence!");

    if (directoryRecordSequences_.empty())
        throw tgt::FileAccessException("Found no DirectoryRecordSequences!");

    if (iteratorIsAtEnd())
        throw tgt::FileAccessException("No DirectoryRecordSequence left!");

    ++it_;
}


void DicomDirParser::setIteratorToFirstDirectoryRecordSequence() const {
    if (filename_.empty())
        throw tgt::FileAccessException("No file specified to access DirectoryRecordSequence!");

    if (directoryRecordSequences_.empty())
        throw tgt::FileAccessException("Found no DirectoryRecordSequences!");

    it_ = directoryRecordSequences_.begin();
}

DirectoryRecordSequenceIterator DicomDirParser::getPresentDirectoryRecordSequence() const {
    if (it_ == directoryRecordSequences_.end())
        throw tgt::FileAccessException("No DirectoryRecordSequence left!");

    DirectoryRecordSequenceIterator seq;
    seq.setDataElement(*it_);
    return seq;
}



} //namespace

