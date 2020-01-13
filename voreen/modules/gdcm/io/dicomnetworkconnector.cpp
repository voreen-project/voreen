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

#include "dicomnetworkconnector.h"

namespace voreen {

const std::string DicomNetworkConnector::loggerCat_ = "voreen.gdcm.DicomNetworkConnector";

DicomNetworkConnector::DicomNetworkConnector(const DicomDict& dict)
{
    dict_ = dict;
    aeTitle_ = "ANY-SCP";
}

DicomNetworkConnector::DicomNetworkConnector(const DicomDict& dict, const std::string& aeTitle, uint16_t incomingPort)
{
    dict_ = dict;
    aeTitle_ = aeTitle;
    incomingPort_ = incomingPort;
}

void DicomNetworkConnector::setIncomingPort(uint16_t port) {
    incomingPort_ = port;
}

void DicomNetworkConnector::setAETitle(const std::string& title) {
    aeTitle_ = title;
}

void DicomNetworkConnector::setProgressBar(ProgressBar* progress) {
    progressBar_ = progress;
}

std::string DicomNetworkConnector::getAETitle() const {
    return aeTitle_;
}

std::string DicomNetworkConnector::getDataElementValue(const gdcm::DataSet& set, const std::string& keyword, const DicomDict* dict /*= 0*/) const {
    if (!dict)
        dict = &dict_;
    std::stringstream stream;
    if (set.FindDataElement(getTagFromDictEntry(dict->getDictEntryByKeyword(keyword)))) {
        if (!set.GetDataElement(getTagFromDictEntry(dict->getDictEntryByKeyword(keyword))).IsEmpty())
            set.GetDataElement(getTagFromDictEntry(dict->getDictEntryByKeyword(keyword))).GetValue().Print(stream);
    }
    return stream.str();
}

bool DicomNetworkConnector::myCMove( const std::string& remote, uint16_t portno, const gdcm::BaseRootQuery* query, const std::string& callee, const std::string& outputdir) const {
    //handle some wrong input parameters
    if(remote.empty())
        throw VoreenException("Empty Remote Server!");

    if(callee.empty())
        throw VoreenException("AE Title of SCP may not be empty!");

    if (portno < 1)
        throw VoreenException("Port Number is not possible: " + itos(portno));

    if (!query)
        throw VoreenException("Trying to use C-FIND with empty query!");

    if (aeTitle_.empty())
        throw VoreenException("AE Title is empty: please set AE Title before contacting any SCP");

    if (incomingPort_ <= 0)
        throw VoreenException("Incoming Port Number is 0: please set Incoming Port Number before contacting any SCP");

    const char* aetitle = aeTitle_.c_str();
    const char* call = callee.c_str();
    const char* dir;

    if (outputdir.empty())
        dir = ".";
    else
        dir = outputdir.c_str();

    // Generate the PresentationContext array from the query UID:
    gdcm::PresentationContextGenerator generator;
    if(!generator.GenerateFromUID(query->GetAbstractSyntaxUID())) {
        throw VoreenException("GDCM Problem: Failed to generate presentation context.");
    }

    //establish the connection
    gdcm::network::ULConnectionManager theManager;
    if (!theManager.EstablishConnectionMove(aetitle, call, remote, 0, portno, 1000, incomingPort_, generator.GetPresentationContexts())) {
        //break the connection
        theManager.BreakConnectionNow();
        throw VoreenException("Failed to establish connection.");
    }

    gdcm::network::ULWritingCallback theCallback;
    theCallback.SetDirectory(dir);
    theManager.SendMove( query, &theCallback );

    theManager.BreakConnection(-1);//wait for a while for the connection to break, ie, infinite
    return true;
}

void DicomNetworkConnector::moveImage(const std::string& remote, uint16_t portno, const std::string& callee, const std::string& patientID,const std::string& studyInstanceUID,const std::string& seriesInstanceUID,const std::string& sopInstanceUID, const std::string& outputDirectory) const {
    //create queryDataSet
    std::vector<std::pair<gdcm::Tag, std::string> > moveQueryDataSet;
    moveQueryDataSet.push_back(std::make_pair(getTagFromDictEntry(dict_.getDictEntryByKeyword("PatientID")),patientID)); //PATIENT ID
    moveQueryDataSet.push_back(std::make_pair(getTagFromDictEntry(dict_.getDictEntryByKeyword("StudyInstanceUID")),studyInstanceUID)); //StudyInstanceUID
    moveQueryDataSet.push_back(std::make_pair(getTagFromDictEntry(dict_.getDictEntryByKeyword("SeriesInstanceUID")),seriesInstanceUID)); //SeriesInstanceUID
    moveQueryDataSet.push_back(std::make_pair(getTagFromDictEntry(dict_.getDictEntryByKeyword("SOPInstanceUID")),sopInstanceUID)); //SOP Instance UID
    //create query
#if GDCM_MAJOR_VERSION == 3 || GDCM_MAJOR_VERSION == 2 && ((GDCM_MINOR_VERSION == 2 && GDCM_BUILD_VERSION >= 2) || (GDCM_MINOR_VERSION > 2))
    gdcm::BaseRootQuery* moveQuery = gdcm::CompositeNetworkFunctions::ConstructQuery(gdcm::eStudyRootType, gdcm::eImage, moveQueryDataSet);
#else
    gdcm::BaseRootQuery* moveQuery = gdcm::CompositeNetworkFunctions::ConstructQuery(gdcm::eStudyRootType, gdcm::eImageOrFrame, moveQueryDataSet);
#endif
    //use C-MOVE to execute query
    bool test;
    try {
        test = myCMove(remote, portno, moveQuery, callee, outputDirectory);
    }
    catch (VoreenException v) {
        delete moveQuery;
        throw v;
    }

    delete moveQuery;

    if (!test) {
        std::string error = "C-MOVE of image failed:\n";
        error.append("Patient ID: "+patientID +"\n");
        error.append("StudyInstanceUID: "+studyInstanceUID+"\n");
        error.append("SeriesInstanceUID: "+seriesInstanceUID+"\n");
        error.append("SOPInstanceUID: "+sopInstanceUID);
        throw VoreenException(error);
    }
}

void DicomNetworkConnector::moveSeries(VolumeURL origin, const std::string& path) const {
    //first: execute C-FIND query to get all images for the series
    if (progressBar_) {
        progressBar_->setTitle("Moving Series...");
        progressBar_->setProgress(0.f);
        progressBar_->setProgressMessage("Executing C-FIND to find images.");
        progressBar_->update();
    }
    std::set<gdcm::Tag> emptySet;
    std::vector<gdcm::DataSet> imageDataSetVector;

    try {
        imageDataSetVector = findImages(origin.getPath(),static_cast<uint16_t>(stoi(origin.getSearchParameter("port"))),origin.getSearchParameter("StudyInstanceUID"), origin.getSearchParameter("SeriesInstanceUID") , emptySet, origin.getSearchParameter("SCP-AET"));
    }
    catch (VoreenException v) {
        progressBar_->hide();
        throw v;
    }

    if (imageDataSetVector.empty()) {
        if (progressBar_)
            progressBar_->hide();
        throw VoreenException("Found no images for this Series!");
    }

    //now move all images of this Series
    int progressCount = 0;
    std::vector<gdcm::DataSet>::iterator imageIterator;
    for (imageIterator = imageDataSetVector.begin(); imageIterator != imageDataSetVector.end(); ++imageIterator) {
        if (progressBar_) {
            progressBar_->setTitle("Moving Series...");
            progressBar_->setProgress(static_cast<float>(progressCount)/static_cast<float>(imageDataSetVector.size()));
            progressBar_->setProgressMessage("Moving images...");
            progressBar_->update();
            progressCount++;
        }

        try {
            moveImage(origin.getPath(), static_cast<uint16_t>(stoi(origin.getSearchParameter("port"))), origin.getSearchParameter("SCP-AET"), origin.getSearchParameter("PatientID"),origin.getSearchParameter("StudyInstanceUID"), origin.getSearchParameter("SeriesInstanceUID") , getDataElementValue(*imageIterator, "SOPInstanceUID"), path);
        }
        catch (VoreenException v) {
            LERROR(v.what());
        }
    }
    if (progressBar_)
        progressBar_->hide();
}

std::vector<gdcm::DataSet> DicomNetworkConnector::myCFind( const std::string& remote, uint16_t portno, const gdcm::BaseRootQuery* query, const std::string& callee) const {
    //handle some wrong input parameters
    if(remote.empty())
        throw VoreenException("Empty Remote Server!");

    if(callee.empty())
        throw VoreenException("AE Title of SCP may not be empty!");

    if (portno < 1)
        throw VoreenException("Port Number is not possible: " + itos(portno));

    if (!query)
        throw VoreenException("Trying to use C-FIND with empty query!");

    if (aeTitle_.empty())
        throw VoreenException("AE Title is empty: please set AE Title before contacting any SCP");

    const char* aetitle = aeTitle_.c_str();
    const char* call = callee.c_str();

    // Generate the PresentationContext array from the query UID:
    gdcm::PresentationContextGenerator generator;
    if(!generator.GenerateFromUID(query->GetAbstractSyntaxUID())) {
        throw VoreenException("GDCM Problem: Failed to generate presentation context.");
    }

    //establish the connection
    gdcm::network::ULConnectionManager theManager;
    if (!theManager.EstablishConnection(aetitle, call, remote, 0, portno, 1000, generator.GetPresentationContexts())) {
        //break the connection
        theManager.BreakConnectionNow();
        throw VoreenException("Failed to establish connection.");
    }

    //send the actual query
    std::vector<gdcm::DataSet> retDataSets = theManager.SendFind(query);

    //break the connection
    theManager.BreakConnection(-1);//wait for a while for the connection to break, ie, infinite

    //return the vector
    return retDataSets;
}

std::vector<gdcm::DataSet> DicomNetworkConnector::findPatients(const std::string& remote, uint16_t portno, const std::string& callee, const std::string& patientID) const {
    //create queryDataSet
    std::vector<std::pair<gdcm::Tag, std::string> > patientQueryDataSet;
    patientQueryDataSet.push_back(std::make_pair(getTagFromDictEntry(dict_.getDictEntryByKeyword("PatientName")),"")); //PATIENT NAME
    patientQueryDataSet.push_back(std::make_pair(getTagFromDictEntry(dict_.getDictEntryByKeyword("PatientID")),patientID)); //PATIENT ID
    //create query
    gdcm::BaseRootQuery* patientFindQuery = gdcm::CompositeNetworkFunctions::ConstructQuery(gdcm::ePatientRootType, gdcm::ePatient, patientQueryDataSet);
    //use C-FIND to execute query
    std::vector<gdcm::DataSet> patientDataSetVector;
    try {
        patientDataSetVector = myCFind(remote, portno, patientFindQuery, callee);
    }
    catch (VoreenException v) {
        delete patientFindQuery;
        throw v;
    }
    delete patientFindQuery;
    return patientDataSetVector;
}

std::vector<gdcm::DataSet> DicomNetworkConnector::findStudies(const std::string& remote, uint16_t portno, const std::string& patientID, const std::string& callee, const std::string& studyInstanceUID) const {
    //query data set
    std::vector<std::pair<gdcm::Tag, std::string> > studyQueryDataSet;
    //add PatientID
    studyQueryDataSet.push_back(std::make_pair(getTagFromDictEntry(dict_.getDictEntryByKeyword("PatientID")),patientID)); //PATIENT ID
    //add the relevant study tags
    studyQueryDataSet.push_back(std::make_pair(getTagFromDictEntry(dict_.getDictEntryByKeyword("StudyDate")),"")); //STUDY DATE
    studyQueryDataSet.push_back(std::make_pair(getTagFromDictEntry(dict_.getDictEntryByKeyword("StudyTime")),"")); //STUDY TIME
    studyQueryDataSet.push_back(std::make_pair(getTagFromDictEntry(dict_.getDictEntryByKeyword("StudyDescription")),"")); //STUDY DESCRIPTION
    studyQueryDataSet.push_back(std::make_pair(getTagFromDictEntry(dict_.getDictEntryByKeyword("StudyInstanceUID")),studyInstanceUID)); //STUDY INSTANCE UID

    //create query
    gdcm::BaseRootQuery* studyQuery = gdcm::CompositeNetworkFunctions::ConstructQuery(gdcm::eStudyRootType, gdcm::eStudy, studyQueryDataSet);
    //use C-FIND to execute query
    std::vector<gdcm::DataSet> studyDataSetVector;
    try {
        studyDataSetVector = myCFind(remote, portno, studyQuery, callee);
    }
    catch (VoreenException v) {
        delete studyQuery;
        throw v;
    }

    delete studyQuery;
    return studyDataSetVector;
}

std::vector<gdcm::DataSet> DicomNetworkConnector::findSeries(const std::string& remote, uint16_t portno, const std::string& patientID, const std::string& studyInstanceUID,const std::string& callee, const std::string& seriesInstanceUID) const {
    //query data set
    std::vector<std::pair<gdcm::Tag, std::string> > seriesQueryDataSet;
    //add PatientID
    seriesQueryDataSet.push_back(std::make_pair(getTagFromDictEntry(dict_.getDictEntryByKeyword("PatientID")),patientID)); //PATIENT ID
    //add StudyInstanceUID
    seriesQueryDataSet.push_back(std::make_pair(getTagFromDictEntry(dict_.getDictEntryByKeyword("StudyInstanceUID")),studyInstanceUID)); //StudyInstanceUID

    //add the relevant series tags
    seriesQueryDataSet.push_back(std::make_pair(getTagFromDictEntry(dict_.getDictEntryByKeyword("SeriesInstanceUID")),seriesInstanceUID)); //SERIES INSTANCE UID
    seriesQueryDataSet.push_back(std::make_pair(getTagFromDictEntry(dict_.getDictEntryByKeyword("SeriesDescription")),"")); //SERIES DESCRIPTION
    seriesQueryDataSet.push_back(std::make_pair(getTagFromDictEntry(dict_.getDictEntryByKeyword("SeriesDate")),""));
    seriesQueryDataSet.push_back(std::make_pair(getTagFromDictEntry(dict_.getDictEntryByKeyword("SeriesTime")),""));

    //create query
    gdcm::BaseRootQuery* seriesQuery = gdcm::CompositeNetworkFunctions::ConstructQuery(gdcm::eStudyRootType, gdcm::eSeries, seriesQueryDataSet);
    //use C-FIND to execute query
    std::vector<gdcm::DataSet> seriesDataSetVector;
    try {
        seriesDataSetVector = myCFind(remote, portno, seriesQuery, callee);
    }
    catch (VoreenException v) {
        delete seriesQuery;
        throw v;
    }

    delete seriesQuery;

    return seriesDataSetVector;
}

std::vector<gdcm::DataSet> DicomNetworkConnector::findImages(const std::string& remote, uint16_t portno, const std::string& studyInstanceUID, const std::string& seriesInstanceUID, std::set<gdcm::Tag> imageLevelTags, const std::string& callee) const {
    //query data set
    std::vector<std::pair<gdcm::Tag, std::string> > imageQueryDataSet;
    //add StudyInstanceUID
    imageQueryDataSet.push_back(std::make_pair(getTagFromDictEntry(dict_.getDictEntryByKeyword("StudyInstanceUID")),studyInstanceUID)); //StudyInstanceUID
    //add SeriesInstanceUID
    imageQueryDataSet.push_back(std::make_pair(getTagFromDictEntry(dict_.getDictEntryByKeyword("SeriesInstanceUID")),seriesInstanceUID)); //SERIES INSTANCE UID

    //relevant Image tags
    imageQueryDataSet.push_back(std::make_pair(getTagFromDictEntry(dict_.getDictEntryByKeyword("SOPInstanceUID")),"")); //SOP Instance UID

    //add all the tags that should be added
    std::set<gdcm::Tag>::iterator tagIterator;
    for (tagIterator = imageLevelTags.begin(); tagIterator != imageLevelTags.end(); ++tagIterator) {
        imageQueryDataSet.push_back(std::make_pair(*tagIterator, ""));
    }

    //create query
#if GDCM_MAJOR_VERSION == 3 || GDCM_MAJOR_VERSION == 2 && ((GDCM_MINOR_VERSION == 2 && GDCM_BUILD_VERSION >= 2) || (GDCM_MINOR_VERSION > 2))
    gdcm::BaseRootQuery* imageQuery = gdcm::CompositeNetworkFunctions::ConstructQuery(gdcm::eStudyRootType, gdcm::eImage, imageQueryDataSet);
#else
    gdcm::BaseRootQuery* imageQuery = gdcm::CompositeNetworkFunctions::ConstructQuery(gdcm::eStudyRootType, gdcm::eImageOrFrame, imageQueryDataSet);
#endif
    //use C-FIND to execute query
    std::vector<gdcm::DataSet> imageDataSetVector;
    try {
        imageDataSetVector = myCFind(remote, portno, imageQuery, callee);
    }
    catch (VoreenException v) {
        delete imageQuery;
        throw v;
    }

    delete imageQuery;

    return imageDataSetVector;
}

std::vector<PatientInfo> DicomNetworkConnector::findNetworkPatients(const std::string& remote, const std::string& call, uint16_t portno) const {
    std::vector<PatientInfo> result;

    std::vector<gdcm::DataSet> patientDataSetVector = findPatients(remote,portno,call,"");

    std::vector<gdcm::DataSet>::const_iterator patientIterator;
    for (patientIterator = patientDataSetVector.begin(); patientIterator != patientDataSetVector.end(); ++patientIterator) {

        //check for PatientID
        std::string presentPatientID = getDataElementValue(*patientIterator,"PatientID");
        if (presentPatientID.empty()) {
            LWARNING("Found Patient with empty PatientID...skipping Patient.");
            continue;
        }

        PatientInfo p;
        p.patientID_ = presentPatientID;
        p.patientName_ = getDataElementValue(*patientIterator,"PatientName");
        result.push_back(p);
    }

    return result;
}

std::vector<StudyInfo> DicomNetworkConnector::findNetworkStudies(const std::string& remote, const std::string& call, uint16_t portno, const std::string& patientID) const {
    std::vector<StudyInfo> result;

    std::vector<gdcm::DataSet> studyDataSetVector = findStudies(remote,portno,patientID,call,"");

    std::vector<gdcm::DataSet>::const_iterator studyIterator;
    for (studyIterator = studyDataSetVector.begin(); studyIterator != studyDataSetVector.end(); ++studyIterator) {

        //check for PatientID
        std::string presentStudyID = getDataElementValue(*studyIterator,"StudyInstanceUID");
        if (presentStudyID.empty()) {
            LWARNING("Found Study with empty StudyInstanceUID...skipping Study.");
            continue;
        }

        StudyInfo s;
        s.studyInstanceUID_ = presentStudyID;
        s.studyDate_ = getDataElementValue(*studyIterator,"StudyDate");
        s.studyDescription_ = getDataElementValue(*studyIterator,"StudyDescription");
        s.studyTime_ = getDataElementValue(*studyIterator,"StudyTime");
        result.push_back(s);
    }

    return result;
}

std::vector<SeriesInfo> DicomNetworkConnector::findNetworkSeries(const std::string& remote, const std::string& call, uint16_t portno, const std::string& patientID, const std::string& studyID) const {
    std::vector<SeriesInfo> result;

    std::vector<gdcm::DataSet> seriesDataSetVector = findSeries(remote, portno,patientID, studyID, call, "");

    std::vector<gdcm::DataSet>::const_iterator seriesIterator;
    for (seriesIterator = seriesDataSetVector.begin(); seriesIterator != seriesDataSetVector.end(); ++seriesIterator) {

        //check for PatientID
        std::string presentSeriesID = getDataElementValue(*seriesIterator,"SeriesInstanceUID");
        if (presentSeriesID.empty()) {
            LWARNING("Found Series with empty SeriesInstanceUID...skipping Series.");
            continue;
        }

        SeriesInfo s;
        s.seriesInstanceUID_ = presentSeriesID;
        s.seriesDate_ = getDataElementValue(*seriesIterator,"SeriesDate");
        s.seriesDescription_ = getDataElementValue(*seriesIterator,"SeriesDescription");
        s.seriesTime_ = getDataElementValue(*seriesIterator,"SeriesTime");

        //find number of images
        std::set<gdcm::Tag> emptySet;
        std::vector<gdcm::DataSet> images = findImages(remote,portno,studyID,presentSeriesID,emptySet,call);
        s.numberOfImages_ = static_cast<int>(images.size());

        result.push_back(s);
    }

    return result;
}

/*std::vector<VolumeURL> DicomNetworkConnector::listVolumes(const std::string& remote, const std::string& call, uint16_t portno, const std::string& patientID, const std::string& studyInstanceUID, const std::string& seriesInstanceUID) const {

    //get all the tags required in any of the CustomDicomDicts, later needed for image level queries
    std::set<gdcm::Tag> customTags = getRelevantTagsFromCustomDicts();

    //hierarchical C-FIND queries: PATIENT level - STUDY level - SERIES level
    std::vector<VolumeURL> result; //result vector
    //std::vector<gdcm::DataSet> seriesDataSetVector; //series vector

    // PATIENT level
    std::vector<gdcm::DataSet> patientDataSetVector = findPatients(remote,portno,call,patientID);

    if (patientDataSetVector.empty())
        throw VoreenException("No Patients found!");

    // STUDY LEVEL
    //iterate over patients and try to find studies
    std::vector<gdcm::DataSet>::const_iterator patientIterator;
    //int patientCounter = 0; //for testing purposes
    for (patientIterator = patientDataSetVector.begin(); patientIterator != patientDataSetVector.end(); ++patientIterator) {

        //check for PatientID
        std::string presentPatientID = getDataElementValue(*patientIterator,"PatientID");
        if (presentPatientID.empty()) {
            LWARNING("Found Patient with empty PatientID...skipping Patient.");
            continue;
        }

        //find studies for this PatientID
        std::vector<gdcm::DataSet> tmpStudyDataSetVector = findStudies(remote,portno,presentPatientID,call,studyInstanceUID);

        // SERIES LEVEL
        //iterate over studies and try to find series
        std::vector<gdcm::DataSet>::const_iterator studyIterator;
        for (studyIterator = tmpStudyDataSetVector.begin(); studyIterator != tmpStudyDataSetVector.end(); ++studyIterator) {

            //check for PATIENT ID
            std::string presentPatientID = getDataElementValue(*studyIterator,"PatientID");
            if (presentPatientID.empty()) {
                LWARNING("Found Study with empty PatientID...skipping Study.");
                continue;
            }

            //check for STUDY INSTANCE UID
            std::string presentStudyInstanceUID = getDataElementValue(*studyIterator,"StudyInstanceUID");
            if (presentStudyInstanceUID.empty()) {
                LWARNING("Found Study with empty StudyInstanceUID...skipping Study.");
                continue;
            }

            //find series
            std::vector<gdcm::DataSet> tmpSeriesDataSetVector = findSeries(remote,portno,presentPatientID,presentStudyInstanceUID,call,seriesInstanceUID);

            //iterate over series
            std::vector<gdcm::DataSet>::const_iterator seriesIterator;
            for (seriesIterator = tmpSeriesDataSetVector.begin(); seriesIterator != tmpSeriesDataSetVector.end(); ++seriesIterator) {
                //check for SERIES INSTANCE UID
                std::string presentSeriesInstanceUID = getDataElementValue(*seriesIterator,"SeriesInstanceUID");
                if (presentSeriesInstanceUID.empty()) {
                    LWARNING("Found Series with empty SeriesInstanceUID...skipping Series.");
                    continue;
                }

                //find images for this series
                std::vector<gdcm::DataSet> tmpImageDataSetVector = findImages(remote,portno,presentStudyInstanceUID, presentSeriesInstanceUID , customTags, call);

                if (tmpImageDataSetVector.empty()) {
                    LWARNING("Series " + presentSeriesInstanceUID + " has no image files... skipping.");
                    continue;
                }

                //check if the images of this series have to be subdivided using a CustomDicomDict
                std::vector<CustomDicomDict> customDicts = getCustomDicts();
                std::vector<CustomDicomDict>::const_iterator dictIt;
                for (dictIt = customDicts.begin(); dictIt != customDicts.end(); ++dictIt) {
                    //check if the current CustomDicomDict is right for the given files
                    const std::vector<std::pair<std::string, std::vector<std::string> > >* conditions = dictIt->getConditions();

                    bool fits = true; //change if conditions do not apply

                    //iterate over all conditions and check if they apply for the given files
                    std::vector<std::pair<std::string, std::vector<std::string> > >::const_iterator condIterator;
                    for (condIterator = conditions->begin(); condIterator != conditions->end(); condIterator++) {
                        std::string value = getDataElementValue(tmpImageDataSetVector.at(0), condIterator->first, dictIt->getDict());
                        value = trim(value, " ");
                        //iterate over possible values for the given tag
                        fits = false;
                        std::vector<std::string>::const_iterator valueIt;
                        //LERROR(value);
                        for (valueIt = condIterator->second.begin(); valueIt != condIterator->second.end(); valueIt++) {
                            //LERROR(*valueIt);
                            if (*valueIt == value) {
                                fits = true;
                                break;
                            }
                        }
                        //if no values are given, the attribute must be present in the file (=> no empty string)
                        if (condIterator->second.empty() && !value.empty())
                            fits = true;


                        if (!fits)
                            break;
                    }

                    if (fits) {
                        //CustomDicomDict fits for these image files
                        LINFO("Image files of Series " + presentSeriesInstanceUID + " are subdivided by a CustomDicomDict.");

                        std::vector<gdcm::DataSet> subdivisionGroups = subdivideSeriesImagesByCustomDict(tmpImageDataSetVector, *dictIt);

                        //now iterate over the groups and add the corresponding VolumeURLs
                        std::vector<gdcm::DataSet>::const_iterator groupIter;
                        for (groupIter = subdivisionGroups.begin(); groupIter != subdivisionGroups.end(); ++groupIter) {
                            //create VolumeURL and add Search Parameters
                            VolumeURL origin("dicom-scp", remote);
                            origin.addSearchParameter("port", itos(portno));
                            origin.addSearchParameter("SCP-AET",call);
                            origin.addSearchParameter("PatientID",presentPatientID);
                            origin.addSearchParameter("StudyInstanceUID",presentStudyInstanceUID);
                            origin.addSearchParameter("SeriesInstanceUID",presentSeriesInstanceUID);

                            //get information on patient-, study-, and series-level and add it to MetaData of VolumeURL
                            origin.getMetaDataContainer().addMetaData("PatientName", new StringMetaData(getDataElementValue(*patientIterator,"PatientName")));
                            origin.getMetaDataContainer().addMetaData("PatientID", new StringMetaData(getDataElementValue(*patientIterator,"PatientID")));
                            origin.getMetaDataContainer().addMetaData("StudyDescription", new StringMetaData(getDataElementValue(*studyIterator,"StudyDescription")));
                            origin.getMetaDataContainer().addMetaData("StudyTime", new StringMetaData(getDataElementValue(*studyIterator, "StudyTime")));
                            origin.getMetaDataContainer().addMetaData("StudyDate", new StringMetaData(getDataElementValue(*studyIterator,"StudyDate")));
                            origin.getMetaDataContainer().addMetaData("StudyInstanceUID", new StringMetaData(getDataElementValue(*studyIterator,"StudyInstanceUID")));
                            origin.getMetaDataContainer().addMetaData("SeriesDescription", new StringMetaData(getDataElementValue(*seriesIterator,"SeriesDescription")));
                            origin.getMetaDataContainer().addMetaData("SeriesInstanceUID", new StringMetaData(getDataElementValue(*seriesIterator,"SeriesInstanceUID")));
                            origin.getMetaDataContainer().addMetaData("SeriesTime", new StringMetaData(getDataElementValue(*studyIterator, "SeriesTime")));
                            origin.getMetaDataContainer().addMetaData("SeriesDate", new StringMetaData(getDataElementValue(*studyIterator,"SeriesDate")));
                            origin.getMetaDataContainer().addMetaData("Number of images in Series", new StringMetaData(itos(tmpImageDataSetVector.size())));

                            //add subdivision keywords of the custom dict as search parameters
                            const std::vector<std::string> *subdivisionKeywords = dictIt->getSubdivisionKeywords();
                            std::vector<std::string>::const_iterator keyIterator;
                            for (keyIterator = subdivisionKeywords->begin(); keyIterator != subdivisionKeywords->end(); ++keyIterator) {
                                std::string value = getDataElementValue(*groupIter,*keyIterator, dictIt->getDict());
                                value = trim(value, " ");
                                if (!value.empty())
                                    origin.addSearchParameter(*keyIterator, value);
                            }
                        }

                        break; //stop iterating over custom dicts
                    }
                }

                //if no CustomDicom could be found: add VolumeURL for the series
                if (dictIt == customDicts.end()) {
                    //create VolumeURL and add Search Parameters
                    VolumeURL origin("dicom-scp", remote);
                    origin.addSearchParameter("port", itos(portno));
                    origin.addSearchParameter("SCP-AET",call);
                    origin.addSearchParameter("PatientID",presentPatientID);
                    origin.addSearchParameter("StudyInstanceUID",presentStudyInstanceUID);
                    origin.addSearchParameter("SeriesInstanceUID",presentSeriesInstanceUID);

                    //get information on patient-, study-, and series-level and add it to MetaData of VolumeURL
                    origin.getMetaDataContainer().addMetaData("PatientName", new StringMetaData(getDataElementValue(*patientIterator,"PatientName")));
                    origin.getMetaDataContainer().addMetaData("PatientID", new StringMetaData(getDataElementValue(*patientIterator,"PatientID")));
                    origin.getMetaDataContainer().addMetaData("StudyDescription", new StringMetaData(getDataElementValue(*studyIterator,"StudyDescription")));
                    origin.getMetaDataContainer().addMetaData("StudyTime", new StringMetaData(getDataElementValue(*studyIterator, "StudyTime")));
                    origin.getMetaDataContainer().addMetaData("StudyDate", new StringMetaData(getDataElementValue(*studyIterator,"StudyDate")));
                    origin.getMetaDataContainer().addMetaData("StudyInstanceUID", new StringMetaData(getDataElementValue(*studyIterator,"StudyInstanceUID")));
                    origin.getMetaDataContainer().addMetaData("SeriesDescription", new StringMetaData(getDataElementValue(*seriesIterator,"SeriesDescription")));
                    origin.getMetaDataContainer().addMetaData("SeriesInstanceUID", new StringMetaData(getDataElementValue(*seriesIterator,"SeriesInstanceUID")));
                    origin.getMetaDataContainer().addMetaData("SeriesTime", new StringMetaData(getDataElementValue(*studyIterator, "SeriesTime")));
                    origin.getMetaDataContainer().addMetaData("SeriesDate", new StringMetaData(getDataElementValue(*studyIterator,"SeriesDate")));

                    //add VolumeURL to result
                    result.push_back(origin);
                }
            }
        }
    }

    return result;
}*/

std::vector<CustomDicomDict> DicomNetworkConnector::getCustomDicts() const {
    std::vector<CustomDicomDict> specialDicts;

    std::string dirName = VoreenApplication::app()->getModulePath("gdcm") + "/dicts/CustomDicts";
    gdcm::Directory dir;
    dir.Load(dirName);
    std::vector<std::string> files = dir.GetFilenames();

    for (std::vector<std::string>::iterator it = files.begin(); it != files.end(); it++) {
        if (tgt::FileSystem::fileExtension(*it) == "xml") {
            CustomDicomDict s;
            try {
                s.loadFromFile(*it);
                specialDicts.push_back(s);
            }
            catch (tgt::Exception e) {
                LERROR("Couldn't load CustomDicomDict: " + (*it));
            }
        }
    }

    return specialDicts;
}


std::set<gdcm::Tag> DicomNetworkConnector::getRelevantTagsFromCustomDicts() const {
    std::map<std::string, gdcm::Tag> keywordMap;
    std::vector<CustomDicomDict> customDicts = getCustomDicts(); //get all custom dicts

    //iterate over custom dicts
    std::vector<CustomDicomDict>::iterator dictIterator;
    for (dictIterator = customDicts.begin(); dictIterator != customDicts.end(); ++dictIterator) {
        //get keywords for subdividing the images of one series
        const std::vector<std::string> * subdivisionKeywords = dictIterator->getSubdivisionKeywords();
        //iterate over these keywords and add the corresponding tag to the map
        std::vector<std::string>::const_iterator subIterator;
        for (subIterator = subdivisionKeywords->begin(); subIterator != subdivisionKeywords->end(); ++subIterator) {
            keywordMap.insert(std::make_pair(*subIterator,getTagFromDictEntry(dictIterator->getDict()->getDictEntryByKeyword(*subIterator))));
        }

        //get the conditions and add the keyword/tag combinations
        const std::vector<std::pair<std::string, std::vector<std::string> > >* conditions = dictIterator->getConditions();
        std::vector<std::pair<std::string, std::vector<std::string> > >::const_iterator condIter;
        for (condIter = conditions->begin(); condIter != conditions->end(); ++condIter) {
            keywordMap.insert(std::make_pair(condIter->first,getTagFromDictEntry(dictIterator->getDict()->getDictEntryByKeyword(condIter->first))));
        }
    }

    std::set<gdcm::Tag> result;

    //iterate over the map and insert the tags into the set
    std::map<std::string, gdcm::Tag>::iterator mapIter;
    for (mapIter = keywordMap.begin(); mapIter != keywordMap.end(); ++mapIter) {
        result.insert(mapIter->second);
    }

    return result;
}

std::vector<gdcm::DataSet> DicomNetworkConnector::subdivideSeriesImagesByCustomDict(std::vector<gdcm::DataSet> imageDataSets, CustomDicomDict customDict) const {
    std::map<std::string, gdcm::DataSet> groups;

    //iterate over all images and for each images over all subdivisionTags to get all value combinations in the series for the given images
    std::vector<std::string> subdivisionKeywords = *(customDict.getSubdivisionKeywords());

    std::vector<gdcm::DataSet>::const_iterator imageIterator;

    for (imageIterator = imageDataSets.begin(); imageIterator != imageDataSets.end(); ++imageIterator) {

        //iterate over keywords and concatenate one string out of the values of these keywords
        std::string s = "";

        std::vector<std::string>::iterator keywordIterator;
        for (keywordIterator = subdivisionKeywords.begin(); keywordIterator != subdivisionKeywords.end(); ++keywordIterator) {
            //get the value from the file
            std::string t = getDataElementValue(*imageIterator,*keywordIterator,customDict.getDict());
            t = trim(t, " ");

            if (!t.empty())
                s.append(t);
        }

        //add this value combination
        groups.insert(std::make_pair(s,*imageIterator));
    }

    //iterate over map and put into std::vector
    std::vector<gdcm::DataSet> groupDataSets;

    std::map<std::string, gdcm::DataSet>::const_iterator groupIter;
    for (groupIter = groups.begin(); groupIter != groups.end(); ++groupIter) {
        groupDataSets.push_back(groupIter->second);
    }

    return groupDataSets;
}


} //namespace

