/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#ifndef VRN_DICOMNETWORKCONNECTOR_H
#define VRN_DICOMNETWORKCONNECTOR_H

// include this before any GDCM header in order to make sure that C99 types are defined
#include "tgt/types.h"

#include "tgt/exception.h"

#include "voreen/core/utils/stringutils.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "../dicomdict.h"
#include "../customdicomdict.h"
#include "voreen/core/io/progressbar.h"
#include "voreen/core/voreenapplication.h"

//gdcm stuff
#include <gdcmFile.h>
#include <gdcmStringFilter.h>

//network functionality
#include <gdcmCompositeNetworkFunctions.h>
#include <gdcmPresentationContextGenerator.h>
#include <gdcmULConnectionManager.h>
#include <gdcmULWritingCallback.h>
//[...]

namespace voreen {

typedef struct dicomSeriesInfo_t {
    std::string seriesInstanceUID_;
    std::string seriesDescription_;
    std::string seriesDate_;
    std::string seriesTime_;
    int numberOfImages_;
} SeriesInfo;

typedef struct dicomStudyInfo_t {
    std::string studyInstanceUID_;
    std::string studyDescription_;
    std::string studyDate_;
    std::string studyTime_;
} StudyInfo;

typedef struct dicomPatientInfo_t {
    std::string patientID_;
    std::string patientName_;
} PatientInfo;

/**
 * DicomNetworkConnector implements the DICOM network functionality of GDCM that
 * is used in GdcmVolumeReader
 */
class VRN_CORE_API DicomNetworkConnector
{
public:
    DicomNetworkConnector(const DicomDict& dict);

    DicomNetworkConnector(const DicomDict& dict, const std::string& aeTitle, uint16_t incomingPort);

    void setIncomingPort(uint16_t port);

    void setAETitle(const std::string&);

    void setProgressBar(ProgressBar* progress);

    std::string getAETitle() const;


    /**
     * Try to contact a DICOM SCP and get every Data Set on this server.
     * Returns a list with a VolumeURL for each Data Set (mostly one for each SeriesInstanceUID).
     * If necessary, a series is subdivided using the available CustomDicomDicts.
     *
     * @param remote the remote server to contact
     * @param call the Application Entity Title of the called DICOM SCP
     * @param portno the port number to contact the server
     * @param patientID (optional) parameter to specify a PatientID
     * @param studyInstanceUID (optional) parameter to specify a StudyInstanceUID
     * @param seriesInstanceUID (optional) parameter to specify a Series InstanceUID
     *
     */
    //virtual std::vector<VolumeURL> listVolumes(const std::string& remote, const std::string& call, uint16_t portno, const std::string& patientID = "", const std::string& studyInstanceUID = "", const std::string& seriesInstanceUID = "") const;

    /**
     * Method to be called by GdcmVolumeReader to find patients
     *
     * @param remote the remote server to contact
     * @param call the Application Entity Title of the called DICOM SCP
     * @param portno the port number to contact the server
     */
    std::vector<PatientInfo> findNetworkPatients(const std::string& remote, const std::string& call, uint16_t portno) const;

    /**
     * Method to be called by GdcmVolumeReader to find studies for a patient
     *
     * @param remote the remote server to contact
     * @param call the Application Entity Title of the called DICOM SCP
     * @param portno the port number to contact the server
     * @param patientID parameter to specify the PatientID
     */
    std::vector<StudyInfo> findNetworkStudies(const std::string& remote, const std::string& call, uint16_t portno, const std::string& patientID) const;

    /**
     * Method to be called by GdcmVolumeReader to find series for a study of a patient
     *
     * @param remote the remote server to contact
     * @param call the Application Entity Title of the called DICOM SCP
     * @param portno the port number to contact the server
     * @param patientID parameter to specifies the PatientID
     * @param studyID parameter to specifies a StudyInstanceUID
     */
    std::vector<SeriesInfo> findNetworkSeries(const std::string& remote, const std::string& call, uint16_t portno, const std::string& patientID, const std::string& studyID) const;

    /**
     * Moves the images of a series from a DICOM SCP to a local path.
     *
     * @param origin VolumeURL has to contain the server's url, AE Title, port number as well as PatientID, StudyInstanceUID, and SeriesInstanceUID, as Search Parameters.
     * @param path the local path the images should be moved to
     */
    void moveSeries(VolumeURL origin, const std::string& path) const;

private:
    /**
     * C-FIND functionality... modified from the original method in the GDCM class CompositeNetworkFunctions
     */
    std::vector<gdcm::DataSet> myCFind(const std::string& remote, uint16_t portno, const gdcm::BaseRootQuery* query, const std::string& callee = "ANY-SCP") const;

    /**
     * C-MOVE functionality... modified from the original method in the GDCM class CompositeNetworkFunctions
     */
    bool myCMove( const std::string& remote, uint16_t portno, const gdcm::BaseRootQuery* query, const std::string& callee, const std::string& outputdir) const;

    /**
     * use C-MOVE to get a specific image file
     */
    void moveImage(const std::string& remote, uint16_t portno, const std::string& callee, const std::string& patientID,const std::string& studyInstanceUID,const std::string& seriesInstanceUID,const std::string& sopInstanceUID, const std::string& outputDirectory) const;

    /**
     * uses C-FIND functionality to find patients, optionally with a specific PatientID
     */
    std::vector<gdcm::DataSet> findPatients(const std::string& remote, uint16_t portno, const std::string& callee = "ANY-SCP", const std::string& patientID ="") const;

    /**
     * uses C-FIND functionality to find studies for a specific patient, optionally with a specific StudyInstanceUID
     */
    std::vector<gdcm::DataSet> findStudies(const std::string& remote, uint16_t portno, const std::string& patientID, const std::string& callee = "ANY-SCP", const std::string& studyInstanceUID ="") const;

    /**
     * uses C-FIND functionality to find series for a specific study of a patient, optionally with a specific SeriesInstanceUID
     */
    std::vector<gdcm::DataSet> findSeries(const std::string& remote, uint16_t portno, const std::string& patientID, const std::string& studyInstanceUID,const std::string& callee = "ANY-SCP", const std::string& seriesInstanceUID ="") const;

    /**
     * uses C-FIND functionality to find images for a specific study/series.
     * @param imageLevelTags tags that should be added (have to be accessible on image level)
     */
    std::vector<gdcm::DataSet> findImages(const std::string& remote, uint16_t portno, const std::string& studyInstanceUID, const std::string& seriesInstanceUID, std::set<gdcm::Tag> imageLevelTags, const std::string& callee = "ANY-SCP") const;

    gdcm::Tag getTagFromDictEntry(const DicomDictEntry &entry) const {
       return gdcm::Tag(entry.getGroupTagUint16(),entry.getElementTagUint16());
    }

    /**
     * Returns the value of a DataElement in a DataSet.
     *
     * @param set the DataSet
     * @param keyword specifies the DataElement by using the keyword to get the corresponding tag from the Dictionary
     * @param dict may be used for a specific Dictionary that should be used. Otherwise, the StandardDictionary is used.
     */
    std::string getDataElementValue(const gdcm::DataSet& set, const std::string& keyword, const DicomDict* dict = 0) const;

    /// get the custom dicts that are availabe
    std::vector<CustomDicomDict> getCustomDicts() const;

    /// get the tags that are necessary to test all CustomDicomDicts
    std::set<gdcm::Tag> getRelevantTagsFromCustomDicts() const;

    /// subdivides the images of one series and returns one image data set for each subdivision group (used for listing Volumes)
    std::vector<gdcm::DataSet> subdivideSeriesImagesByCustomDict(std::vector<gdcm::DataSet> imageDataSets, CustomDicomDict customDict) const;

    std::string aeTitle_;
    uint16_t incomingPort_;
    DicomDict dict_;
    ProgressBar* progressBar_;
    static const std::string loggerCat_;
};

} //namespace

#endif // DICOMNETWORKCONNECTOR_H
