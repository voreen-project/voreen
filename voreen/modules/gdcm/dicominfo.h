/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#ifndef VRN_DICOMINFO_H
#define VRN_DICOMINFO_H

#include "voreen/core/voreencoreapi.h"

// include this before any GDCM header in order to make sure that C99 types are defined
#include "tgt/types.h"
#include "tgt/vector.h"

#include <string>

/**
 * This class is used by GdcmVolumeReader to store some information about a DICOM volume that is loaded.
 * It contains basic Meta Information like StudyInstanceUID, SeriesInstanceUID, Patient's Name, Modality, etc.
 * as well as Volume-related information such as information about spacing.
 */
class VRN_CORE_API DicomInfo
{
public:
    DicomInfo();

    void setDx(int dx);
    void setDy(int dy);
    void setDz(int dz);
    void setNumberOfFrames(int frames);
    void setXSpacing(double x);
    void setYSpacing(double y);
    void setZSpacing(double z);
    void setXOrientationPatient(tgt::dvec3 orientation);
    void setYOrientationPatient(tgt::dvec3 orientation);
    void setSliceNormal(tgt::dvec3 normal);
    void setOffset(tgt::dvec3 offset);
    void setPixelRepresentation(unsigned short representation);
    void setBitsStored(int bits);
    void setSamplesPerPixel(int spp);
    void setBytesPerVoxel(int bytes);
    void setIntercept(float intercept);
    void setSlope(float slope);
    void setRescaleType(std::string type);
    void setModality(std::string modality);
    void setSeriesInstanceUID(std::string uid);
    void setStudyInstanceUID(std::string uid);
    void setSeriesDescription(std::string description);
    void setStudyDescription(std::string description);
    void setPatientName(std::string name);
    void setPatientId(std::string id);
    void setRwmDiffers(bool d);
    void setBaseType(std::string type);
    void setFormat(std::string format);


    int getDx() const;
    int getDy() const;
    int getDz() const;
    int getNumberOfFrames() const;
    int getBitsStored() const;
    int getSamplesPerPixel() const;
    int getBytesPerVoxel() const;
    unsigned short getPixelRepresentation() const;
    float getIntercept() const;
    float getSlope() const;
    double getXSpacing() const;
    double getYSpacing() const;
    double getZSpacing() const;
    tgt::dvec3 getXOrientationPatient() const;
    tgt::dvec3 getYOrientationPatient() const;
    tgt::dvec3 getSliceNormal() const;
    tgt::dvec3 getOffset() const;
    std::string getRescaleType() const;
    std::string getModality() const;
    std::string getSeriesInstanceUID() const;
    std::string getStudyInstanceUID() const;
    std::string getSeriesDescription() const;
    std::string getStudyDescription() const;
    std::string getPatientName() const;
    std::string getPatientId() const;
    bool rwmDiffers() const;
    std::string getBaseType() const;
    std::string getFormat() const;

private:

    int dx_, dy_, dz_;          ///< information about the size of the volume in x-, y- and z-direction
    int numberOfFrames_;        ///< for multiframe DICOM files

    //information about patient, study, modality etc.
    std::string patientId_;                 ///< Patient ID of the DICOM Meta Data
    std::string patientName_;               ///< Patient's Name
    std::string studyInstanceUID_;          ///< StudyInstanceUID of the DICOM Meta Data
    std::string studyDescription_;          ///< Study Description of the DICOM Meta Data
    std::string seriesInstanceUID_;         ///< SeriesInstanceUID of the DICOM Meta Data
    std::string seriesDescription_;         ///< Series Description of the DICOM Meta Data
    std::string modality_;                  ///< the Modality of the DICOM file(s)

    //spatial information about the volume
    double xSpacing_, ySpacing_, zSpacing_;     ///< spacing of the volume in x-, y- and z-direction
    tgt::dvec3 xOrientationPatient_;            ///< first vector of ImageOrientationPatient
    tgt::dvec3 yOrientationPatient_;            ///< second vector of ImageOrientationPatient
    tgt::dvec3 sliceNormal_;                    ///< direction in which the slices are arranged
    tgt::dvec3 offset_;                         ///< offset of the Volume's origin

    //pixel storage
    unsigned short pixelRepresentation_;    ///< 0 = unsigned, 1 = signed
    int bitsStored_;                        ///< bits stored per sample
    int bytesPerVoxel_;                     ///< bytes stored per Voxel (samplesPerPixel_ * bits per sample in bytes)

    //voxel format
    std::string baseType_;      ///< base type for format (e.g. "uint16")
    int samplesPerPixel_;       ///< samples stored per pixel (-> together with base type determines format)
    std::string format_;        ///< format string (e.g. "float", "Vector3(uint8)")

    //pixel rescale information
    bool rwmDiffers_;           ///< does the rescale slope and intercept of the image files differ?
    float intercept_;           ///< Rescale Intercept
    float slope_;               ///< Rescale Slope
    std::string rescaleType_;   ///< Rescale Type
};

#endif // DICOMINFO_H
