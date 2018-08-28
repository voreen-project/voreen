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

#include "gdcmvolumereader.h"

// include this before any GDCM header in order to make sure that C99 types are defined
#include "tgt/types.h"

#include <gdcmGlobal.h>
#include <gdcmReader.h>
#include <gdcmImageReader.h>
#include <gdcmXMLDictReader.h>
#include <gdcmStringFilter.h>
#include <gdcmFile.h>
#include <gdcmDirectory.h>
#include <gdcmSystem.h>
#include <gdcmRescaler.h>
//#include <gdcmIPPSorter.h> //currently not in use

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/io/progressbar.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/utils/stringutils.h"
#include "voreen/core/datastructures/meta/templatemetadata.h"
#include "voreen/core/datastructures/meta/filelistmetadata.h"
#include "voreen/core/datastructures/volume/volumefactory.h"

#include <algorithm>

#include "dicomdirparser.h"

namespace voreen {

using std::vector;
using std::pair;
using std::string;

const std::string GdcmVolumeReader::loggerCat_ = "voreen.gdcm.GdcmVolumeReader";

GdcmVolumeReader::GdcmVolumeReader(ProgressBar* progress)
    : VolumeReader(progress)
{
    extensions_.push_back("dcm");
    extensions_.push_back("ima");
    filenames_.push_back("DICOMDIR");
    protocols_.push_back("dicom");

    //DICOM network support
    protocols_.push_back("dicom-scp");

    dict_ = 0;
    lastBufferMod_ = DateTime::now();
}

GdcmVolumeReader::~GdcmVolumeReader() {
    if (dict_) {
        delete dict_;
        dict_ = 0;
    }
}

VolumeReader* GdcmVolumeReader::create(ProgressBar* progress) const {
    return new GdcmVolumeReader(progress);
}

VolumeURL GdcmVolumeReader::convertOriginToRelativePath(const VolumeURL& origin, const std::string& basePath) const {

    if (origin.getProtocol() == "dicom-scp")
        return origin;
    else
        return VolumeReader::convertOriginToRelativePath(origin,basePath);

}

VolumeURL GdcmVolumeReader::convertOriginToAbsolutePath(const VolumeURL& origin, const std::string& basePath) const {

    if (origin.getProtocol() == "dicom-scp")
        return origin;
    else
        return VolumeReader::convertOriginToAbsolutePath(origin,basePath);
}

gdcm::Tag GdcmVolumeReader::getTagFromDictEntry(const DicomDictEntry &entry) {
   return gdcm::Tag(entry.getGroupTagUint16(),entry.getElementTagUint16());
}

void GdcmVolumeReader::loadStandardDict() const {

    string fileName = VoreenApplication::app()->getModulePath("gdcm");
    fileName += "/dicts/StandardDictionary.xml";

    dict_ =  DicomDict::loadFromFile(fileName);
}

vector<string> GdcmVolumeReader::getFileNamesInDir(const string& dirName) const {
    gdcm::Directory dir;
    dir.Load(dirName);
    return dir.GetFilenames();
}

vector<CustomDicomDict> GdcmVolumeReader::getCustomDicts() const {
    vector<CustomDicomDict> specialDicts;

    vector<string> files = getFileNamesInDir(VoreenApplication::app()->getModulePath("gdcm") + "/dicts/CustomDicts");

    for (vector<string>::iterator it = files.begin(); it != files.end(); it++) {
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

/*
 * Sorts strings according to a double value interpreted as the distance of the Image Origin to the Volume Origin
 */
bool slices_cmp_dist(std::pair<string, double> a, std::pair<string, double> b) {
   return a.second < b.second;
}

/*
 * Removes slashes as well as backslaches within a std::string
 */
std::string removeSlashes(const std::string& s) {
    return strReplaceAll(strReplaceAll(s,"/",""),"\\","");
}

std::string GdcmVolumeReader::constructLocalPathFromNetworkOrigin(const VolumeURL& o) const {
    string path = VoreenApplication::app()->getTemporaryPath()+"/gdcm/";
    if (!trim(removeSlashes(o.getSearchParameter("PatientID"))).empty())
        path.append(trim(removeSlashes(o.getSearchParameter("PatientID")))+"/");
    if (!trim(removeSlashes(o.getSearchParameter("StudyInstanceUID"))).empty())
        path.append(trim(removeSlashes(o.getSearchParameter("StudyInstanceUID")))+"/");
    if (!trim(removeSlashes(o.getSearchParameter("SeriesInstanceUID"))).empty())
        path.append(trim(removeSlashes(o.getSearchParameter("SeriesInstanceUID")))+"/");
    return path;
}

bool GdcmVolumeReader::isDicomFile(const string &url) const {
    gdcm::Reader reader;
    reader.SetFileName(url.c_str());
    return reader.Read();
}

bool GdcmVolumeReader::isDicomDir(const string &url) const {
    if (!isDicomFile(url)) return false;

    gdcm::Reader reader;
    reader.SetFileName(url.c_str());
    reader.Read();

    gdcm::MediaStorage ms;
    ms.SetFromFile(reader.GetFile());

    return (ms == gdcm::MediaStorage::MediaStorageDirectoryStorage);
}

void GdcmVolumeReader::setMetaDataFromDict(MetaDataContainer* container, const DicomDict* dict, const std::string& file, bool setAll) {

    if ((!container) || (!dict))
        return;

    gdcm::Reader reader;
    reader.SetFileName(file.c_str());
    if ((!reader.Read()) || (!reader.GetFile().GetHeader().IsValid()))
        throw tgt::FileException("Cannot extract meta data from file " + file, file);

    gdcm::StringFilter sf;
    sf.SetFile(reader.GetFile());

    const vector<string> keys = dict->getKeywordVector();

    vector<string>::const_iterator keywordIterator;

    for (keywordIterator = keys.begin(); keywordIterator != keys.end(); ++keywordIterator) {
        //get entry and check if it should be set as MetaData
        DicomDictEntry entry = dict->getDictEntryByKeyword(*keywordIterator);
        if (entry.isMetaData() || setAll) {
            //get the value
            string s = trim(sf.ToString(getTagFromDictEntry(entry)), " ");
            if (!s.empty()) {
                //if value is not empty: construct MetaData and add to MetaDataContainer
                MetaDataBase* m = constructMetaData(entry, s);

                if (container->hasMetaData(*keywordIterator))
                    container->removeMetaData(*keywordIterator);
                container->addMetaData(*keywordIterator, m);
            }
        }
    }
}

std::string GdcmVolumeReader::getMetaDataFromFile(const std::string& filename, const DicomDict& dict, const std::string& keyword) const {
    std::map<string, MetaDataContainer>::iterator mapIterator;
    mapIterator = fileInfoBuffer_.find(filename);

    //if the file is not already present in the file buffer: create a new entry and add the meta information from the DicomDict
    if (mapIterator == fileInfoBuffer_.end()) {
        //put file and corresponding dict meta info into the buffer
        pair<string, MetaDataContainer> newPair;
        newPair.first = filename;
        //MetaDataContainer c;
        //newPair.second = c;
        setMetaDataFromDict(&newPair.second, &dict, filename, true);
        fileInfoBuffer_.insert(newPair);
        mapIterator = fileInfoBuffer_.find(filename);
        //set the timestamp to the current time because the buffer has been modified
        lastBufferMod_ = DateTime::now();
    }
    else if (!(mapIterator->second.hasMetaData(keyword))) {
        //if the filename is present, but does not have the right meta data: add the meta information from the DicomDict
        setMetaDataFromDict(&mapIterator->second, &dict, filename, true);
        //set the timestamp to the current time because the buffer has been modified
        lastBufferMod_ = DateTime::now();
    }

    //now the filename is present, the meta information is added and the iterator is set to the right position
    //extract the meta data and return it
    if (mapIterator->second.hasMetaData(keyword))
        return mapIterator->second.getMetaData(keyword)->toString();
    else
        return "";
}

MetaDataBase* GdcmVolumeReader::constructMetaData(const DicomDictEntry &entry, const std::string &valueString) {

    string vr = entry.getValueRepresentation();
    string vm = entry.getValueMultiplicity();
    string value = valueString;
    value = trim(value, " ");
    vr = trim(vr, " ");
    vm = trim(vm, " ");

    //select by ValueRepresentation and ValueMultiplicity
    if ((vr == "DA") && (vm == "1")) {
        //check for right length of string
        if (value.length() == 8) {
            //return DateTimeMetaData
            DateTime date = DateTime::createDate(stoi(value.substr(0,4)), stoi(value.substr(4,2)), stoi(value.substr(6,2)));
            DateTimeMetaData* d = new DateTimeMetaData(date);
            return d;
        }
        else if (value.length() == 10) {
            //return DateTimeMetaData
            DateTime date = DateTime::createDate(stoi(value.substr(0,4)), stoi(value.substr(5,2)), stoi(value.substr(8,2)));
            DateTimeMetaData* d = new DateTimeMetaData(date);
            return d;
        }
        else {
            //return StringMetaData
            LWARNING("Unknown Format in DICOM tag with ValueRepresentation DA: " + entry.getKeyword());
            StringMetaData* s = new StringMetaData(value);
            return s;
        }
    }
    else if ((vr == "DT") && (vm == "1")) {
        //check for right length of string
        if (value.length() == 14) {
            //return DateTimeMetaData
            DateTime dt(stoi(value.substr(0,4)), stoi(value.substr(4,2)), stoi(value.substr(6,2)), stoi(value.substr(8,2)), stoi(value.substr(10,2)), stoi(value.substr(12,2)), 0);
            DateTimeMetaData* d = new DateTimeMetaData(dt);
            return d;
        }
        else {
            //return StringMetaData
            LWARNING("Unknown Format in DICOM tag with ValueRepresentation DT: " + entry.getKeyword());
            StringMetaData* s = new StringMetaData(value);
            return s;
        }
    }
    else if ((vr == "TM") && (vm == "1")) {

        if (value.find_first_of(':') != string::npos) {  //older format
            //check for length of string
            if (value.length() >= 8) {
                //use DateTimeMetaData only with time
                DateTime time(DateTime::TIME);

                //check for fractional seconds
                int ms = 0;
                if (value.length() > 9) {
                    string fracString = value.substr(9);
                    if (fracString.length() > 3) {
                        ms = stoi(fracString.substr(0, 3));
                        if (stoi(fracString.substr(3,1)) >= 5)
                            ms++;
                    }
                    else if (fracString.length() == 3)
                        ms = stoi(fracString);
                    else if (fracString.length() == 2)
                        ms = stoi(fracString)*10;
                    else if (fracString.length() == 1)
                        ms = stoi(fracString)*100;
                }
                time.setTime(stoi(value.substr(0,2)), stoi(value.substr(3,2)), stoi(value.substr(6,2)), ms);
                DateTimeMetaData* d = new DateTimeMetaData(time);
                return d;
            }
            else if (value.length() == 5) {
                //use DateTimeMetaData only with time
                //only hours and minutes
                DateTime time(DateTime::TIME);
                time.setTime(stoi(value.substr(0,2)), stoi(value.substr(3,2)), 0, 0);
                DateTimeMetaData* d = new DateTimeMetaData(time);
                return d;
            }
            else {
                //return StringMetaData
                LWARNING("Unknown Format in DICOM tag with ValueRepresentation TM: " + entry.getKeyword());
                StringMetaData* s = new StringMetaData(value);
                return s;
            }
        }

        //check for length of string
        if (value.length() >= 6) {
            //use DateTimeMetaData only with time
            DateTime time(DateTime::TIME);

            //check for fractional seconds
            int ms = 0;
            if (value.length() > 7) {
                string fracString = value.substr(7);
                if (fracString.length() > 3) {
                    ms = stoi(fracString.substr(0, 3));
                    if (stoi(fracString.substr(3,1)) >= 5)
                        ms++;
                }
                else if (fracString.length() == 3)
                    ms = stoi(fracString);
                else if (fracString.length() == 2)
                    ms = stoi(fracString)*10;
                else if (fracString.length() == 1)
                    ms = stoi(fracString)*100;
            }

            time.setTime(stoi(value.substr(0,2)), stoi(value.substr(2,2)), stoi(value.substr(4,2)), ms);
            DateTimeMetaData* d = new DateTimeMetaData(time);
            return d;
        }
        else if (value.length() == 4) {
            //use DateTimeMetaData only with time
            //just hours and minutes
            DateTime time(DateTime::TIME);
            time.setTime(stoi(value.substr(0,2)), stoi(value.substr(2,2)), 0, 0);
            DateTimeMetaData* d = new DateTimeMetaData(time);
            return d;
        }
        else if (value.length() == 2) {
            //use DateTimeMetaData only with time
            //just hours
            DateTime time(DateTime::TIME);
            time.setTime(stoi(value.substr(0,2)), 0, 0, 0);
            DateTimeMetaData* d = new DateTimeMetaData(time);
            return d;
        }
        else {
            //return StringMetaData
            LWARNING("Unknown Format in DICOM tag with ValueRepresentation TM: " + entry.getKeyword());
            StringMetaData* s = new StringMetaData(value);
            return s;
        }
    }
    else if ((vr == "IS") && (vm == "1")) {
        //string to int
        IntMetaData* m = new IntMetaData();
        int i = stoi(value);
        m->setValue(i);
        return m;
    }
    else if ((vr == "IS") && (vm == "2")) {
        //string to tgt::ivec2
        IVec2MetaData* m = new IVec2MetaData();
        vector<string> values = strSplit(value, '\\');
        tgt::ivec2 v;
        v.x = stoi(values[0]);
        v.y = stoi(values[1]);
        m->setValue(v);
        return m;
    }
    else if ((vr == "IS") && (vm == "3")) {
        //string to tgt::ivec3
        IVec3MetaData* m = new IVec3MetaData();
        vector<string> values = strSplit(value, '\\');
        tgt::ivec3 v;
        v.x = stoi(values[0]);
        v.y = stoi(values[1]);
        v.z = stoi(values[2]);
        m->setValue(v);
        return m;
    }
    else if (((vr == "FL") || (vr == "DS") ||(vr == "OF")) && (vm == "1")) {
        //string to float
        FloatMetaData* m = new FloatMetaData();
        float f = stof(value);
        m->setValue(f);
        return m;
    }
    else if (((vr == "FL") || (vr == "DS") ||(vr == "OF")) && (vm == "2")) {
        //string to tgt::vec2
        Vec2MetaData* m = new Vec2MetaData();
        vector<string> values = strSplit(value, '\\');
        tgt::vec2 v;
        v.x = stof(values[0]);
        v.y = stof(values[1]);
        m->setValue(v);
        return m;
    }
    else if (((vr == "FL") || (vr == "DS") ||(vr == "OF")) && (vm == "3")) {
        //string to tgt::vec3
        Vec3MetaData* m = new Vec3MetaData();
        vector<string> values = strSplit(value, '\\');
        tgt::vec3 v;
        v.x = stof(values[0]);
        v.y = stof(values[1]);
        v.z = stof(values[2]);
        m->setValue(v);
        return m;
    }
    else if ((vr == "FD") && (vm == "1")) {
        //string to double
        DoubleMetaData* m = new DoubleMetaData();
        value = trim(value, " ");
        double d = stof(value);
        m->setValue(d);
        return m;
    }
    else if ((vr == "FD") && (vm == "2")) {
        //string to tgt::dvec2
        DVec2MetaData* m = new DVec2MetaData();
        vector<string> values = strSplit(value, '\\');
        tgt::dvec2 v;
        v.x = stod(values[0]);
        v.y = stod(values[1]);
        m->setValue(v);
        return m;
    }
    else if ((vr == "FD") && (vm == "3")) {
        //string to tgt::dvec3
        DVec3MetaData* m = new DVec3MetaData();
        vector<string> values = strSplit(value, '\\');
        tgt::dvec3 v;
        v.x = stod(values[0]);
        v.y = stod(values[1]);
        v.z = stod(values[2]);
        m->setValue(v);
        return m;
    }
    else {
        //other type: just use the string
        StringMetaData* m = new StringMetaData();
        m->setValue(value);
        return m;
    }

}



MetaDataBase* GdcmVolumeReader::constructVolumeDateTime(const DicomDict* dict, gdcm::StringFilter* sf) const {
    if ((!dict) || (!sf))
        return 0;

    string acquisitionDateTime  = trim(sf->ToString(getTagFromDictEntry(dict->getDictEntryByKeyword("AcquisitionDateTime"))));
    string acquisitionDate  = trim(sf->ToString(getTagFromDictEntry(dict->getDictEntryByKeyword("AcquisitionDate"))));
    string acquisitionTime  = trim(sf->ToString(getTagFromDictEntry(dict->getDictEntryByKeyword("AcquisitionTime"))));
    string seriesDate  = trim(sf->ToString(getTagFromDictEntry(dict->getDictEntryByKeyword("SeriesDate"))));
    string seriesTime  = trim(sf->ToString(getTagFromDictEntry(dict->getDictEntryByKeyword("SeriesTime"))));
    string studyDate  = trim(sf->ToString(getTagFromDictEntry(dict->getDictEntryByKeyword("StudyDate"))));
    string studyTime  = trim(sf->ToString(getTagFromDictEntry(dict->getDictEntryByKeyword("StudyTime"))));

    if (!acquisitionDateTime.empty()) {
        //use AcquisitionDateTime
        MetaDataBase* m = constructMetaData(dict_->getDictEntryByKeyword("AcquisitionDateTime"),acquisitionDateTime);
        return m;

    }
    else if (!acquisitionDate.empty() && !acquisitionTime.empty()) {
        //use AcquisitionDate and AcquisitionTime
        string year, month, day, hour, minute, second, millisecond;

        //check for right length of strings and get date
        if (acquisitionDate.length() == 8) {
            year = acquisitionDate.substr(0,4);
            month = acquisitionDate.substr(4,2);
            day = acquisitionDate.substr(6,2);
        }
        else if (acquisitionDate.length() == 10) {
            year = acquisitionDate.substr(0,4);
            month = acquisitionDate.substr(5,2);
            day = acquisitionDate.substr(8,2);
        }
        else {
            LWARNING("Could not create VolumeDateTime, because of unknown format in AcquisitionDate");
            return 0;
        }

        //get time
        if (acquisitionTime.find_first_of(':') != string::npos) {  //older format
            //check for length of string
            if (acquisitionTime.length() >= 8) {
                hour = acquisitionTime.substr(0,2);
                minute = acquisitionTime.substr(3,2);
                second = acquisitionTime.substr(6,2);

                millisecond = "0";

                //check for fractional seconds
                if (acquisitionTime.length() > 9) {
                    string fracString = acquisitionTime.substr(9);
                    if (fracString.length() > 3) {
                        millisecond = fracString.substr(0,3);
                    }
                    else if (fracString.length() == 3)
                        millisecond = fracString;
                    else if (fracString.length() == 2)
                        millisecond = fracString + "0";
                    else if (fracString.length() == 1)
                        millisecond = fracString + "00";
                }
            }
            else if (acquisitionTime.length() == 5) {
                //only hours and minutes
                hour = acquisitionTime.substr(0,2);
                minute = acquisitionTime.substr(3,2);
                second = "0";
                millisecond = "0";
            }
            else {
                LWARNING("Could not create VolumeDateTime, because of unknown format in AcquisitionTime");
                return 0;
            }
        }
        else {
            //check for length of string
            if (acquisitionTime.length() >= 6) {
                hour = acquisitionTime.substr(0,2);
                minute = acquisitionTime.substr(2,2);
                second = acquisitionTime.substr(4,2);

                millisecond = "0";

                //check for fractional seconds
                if (acquisitionTime.length() > 7) {
                    string fracString = acquisitionTime.substr(7);
                    if (fracString.length() > 3) {
                        millisecond = fracString.substr(0,3);
                    }
                    else if (fracString.length() == 3)
                        millisecond = fracString;
                    else if (fracString.length() == 2)
                        millisecond = fracString + "0";
                    else if (fracString.length() == 1)
                        millisecond = fracString + "00";
                }
            }
            else if (acquisitionTime.length() == 4) {
                //just hours and minutes
                hour = acquisitionTime.substr(0,2);
                minute = acquisitionTime.substr(2,2);
                second = "0";
                millisecond = "0";
            }
            else if (acquisitionTime.length() == 2) {
                //just hours
                hour = acquisitionTime.substr(0,2);
                minute = "0";
                second = "0";
                millisecond = "0";
            }
            else {
                LWARNING("Could not create VolumeDateTime, because of unknown format in AcquisitionTime");
                return 0;
            }
        }

        DateTime dateTime(stoi(year), stoi(month), stoi(day), stoi(hour), stoi(minute), stoi(second), stoi(millisecond));
        return new DateTimeMetaData(dateTime);

    }
    else if (!seriesDate.empty() && !seriesTime.empty()) {
        //use SeriesDate and SeriesTime
        string year, month, day, hour, minute, second, millisecond;

        //check for right length of strings and get date
        if (seriesDate.length() == 8) {
            year = seriesDate.substr(0,4);
            month = seriesDate.substr(4,2);
            day = seriesDate.substr(6,2);
        }
        else if (seriesDate.length() == 10) {
            year = seriesDate.substr(0,4);
            month = seriesDate.substr(5,2);
            day = seriesDate.substr(8,2);
        }
        else {
            LWARNING("Could not create VolumeDateTime, because of unknown format in SeriesDate");
            return 0;
        }

        //get time
        if (seriesTime.find_first_of(':') != string::npos) {  //older format
            //check for length of string
            if (seriesTime.length() >= 8) {
                hour = seriesTime.substr(0,2);
                minute = seriesTime.substr(3,2);
                second = seriesTime.substr(6,2);

                millisecond = "0";

                //check for fractional seconds
                if (seriesTime.length() > 9) {
                    string fracString = seriesTime.substr(9);
                    if (fracString.length() > 3) {
                        millisecond = fracString.substr(0,3);
                    }
                    else if (fracString.length() == 3)
                        millisecond = fracString;
                    else if (fracString.length() == 2)
                        millisecond = fracString + "0";
                    else if (fracString.length() == 1)
                        millisecond = fracString + "00";
                }
            }
            else if (seriesTime.length() == 5) {
                //only hours and minutes
                hour = seriesTime.substr(0,2);
                minute = seriesTime.substr(3,2);
                second = "0";
                millisecond = "0";
            }
            else {
                LWARNING("Could not create VolumeDateTime, because of unknown format in SeriesTime");
                return 0;
            }
        }
        else {
            //check for length of string
            if (seriesTime.length() >= 6) {
                hour = seriesTime.substr(0,2);
                minute = seriesTime.substr(2,2);
                second = seriesTime.substr(4,2);

                millisecond = "0";

                //check for fractional seconds
                if (seriesTime.length() > 7) {
                    string fracString = seriesTime.substr(7);
                    if (fracString.length() > 3) {
                        millisecond = fracString.substr(0,3);
                    }
                    else if (fracString.length() == 3)
                        millisecond = fracString;
                    else if (fracString.length() == 2)
                        millisecond = fracString + "0";
                    else if (fracString.length() == 1)
                        millisecond = fracString + "00";
                }
            }
            else if (seriesTime.length() == 4) {
                //just hours and minutes
                hour = seriesTime.substr(0,2);
                minute = seriesTime.substr(2,2);
                second = "0";
                millisecond = "0";
            }
            else if (seriesTime.length() == 2) {
                //just hours
                hour = seriesTime.substr(0,2);
                minute = "0";
                second = "0";
                millisecond = "0";
            }
            else {
                LWARNING("Could not create VolumeDateTime, because of unknown format in SeriesTime");
                return 0;
            }
        }

        DateTime dateTime(stoi(year), stoi(month), stoi(day), stoi(hour), stoi(minute), stoi(second), stoi(millisecond));
        return new DateTimeMetaData(dateTime);

    }
    else if (!studyDate.empty() && !studyTime.empty()) {
        //use StudyDate and StudyTime
        string year, month, day, hour, minute, second, millisecond;

        //check for right length of strings and get date
        if (studyDate.length() == 8) {
            year = studyDate.substr(0,4);
            month = studyDate.substr(4,2);
            day = studyDate.substr(6,2);
        }
        else if (studyDate.length() == 10) {
            year = studyDate.substr(0,4);
            month = studyDate.substr(5,2);
            day = studyDate.substr(8,2);
        }
        else {
            LWARNING("Could not create VolumeDateTime, because of unknown format in StudyDate");
            return 0;
        }

        //get time
        if (studyTime.find_first_of(':') != string::npos) {  //older format
            //check for length of string
            if (studyTime.length() >= 8) {
                hour = studyTime.substr(0,2);
                minute = studyTime.substr(3,2);
                second = studyTime.substr(6,2);

                millisecond = "0";

                //check for fractional seconds
                if (studyTime.length() > 9) {
                    string fracString = studyTime.substr(9);
                    if (fracString.length() > 3) {
                        millisecond = fracString.substr(0,3);
                    }
                    else if (fracString.length() == 3)
                        millisecond = fracString;
                    else if (fracString.length() == 2)
                        millisecond = fracString + "0";
                    else if (fracString.length() == 1)
                        millisecond = fracString + "00";
                }
            }
            else if (studyTime.length() == 5) {
                //only hours and minutes
                hour = studyTime.substr(0,2);
                minute = studyTime.substr(3,2);
                second = "0";
                millisecond = "0";
            }
            else {
                LWARNING("Could not create VolumeDateTime, because of unknown format in StudyTime");
                return 0;
            }
        }
        else {
            //check for length of string
            if (studyTime.length() >= 6) {
                hour = studyTime.substr(0,2);
                minute = studyTime.substr(2,2);
                second = studyTime.substr(4,2);

                millisecond = "0";

                //check for fractional seconds
                if (studyTime.length() > 7) {
                    string fracString = studyTime.substr(7);
                    if (fracString.length() > 3) {
                        millisecond = fracString.substr(0,3);
                    }
                    else if (fracString.length() == 3)
                        millisecond = fracString;
                    else if (fracString.length() == 2)
                        millisecond = fracString + "0";
                    else if (fracString.length() == 1)
                        millisecond = fracString + "00";
                }
            }
            else if (studyTime.length() == 4) {
                //just hours and minutes
                hour = studyTime.substr(0,2);
                minute = studyTime.substr(2,2);
                second = "0";
                millisecond = "0";
            }
            else if (studyTime.length() == 2) {
                //just hours
                hour = studyTime.substr(0,2);
                minute = "0";
                second = "0";
                millisecond = "0";
            }
            else {
                LWARNING("Could not create VolumeDateTime, because of unknown format in StudyTime");
                return 0;
            }
        }

        DateTime dateTime(stoi(year), stoi(month), stoi(day), stoi(hour), stoi(minute), stoi(second), stoi(millisecond));
        return new DateTimeMetaData(dateTime);

    }

    LWARNING("Could not construct VolumeDateTime, because no information was found.");
    return 0;
}

vector<string> GdcmVolumeReader::getFilesInSeries(vector<string> filenames, string seriesInstanceUID) const {
    vector<string> result;

    //show Progress Bar
    if (getProgressBar() && !filenames.empty())
        getProgressBar()->setTitle("Assigning files to SeriesInstanceUID: " + seriesInstanceUID);

    int progress = 1;

    //iterate over files
    vector<string>::iterator fileIterator;
    for (fileIterator = filenames.begin(); fileIterator != filenames.end(); fileIterator++) {
        if (getProgressBar()) {
            getProgressBar()->setProgressMessage("Reading files...");
            getProgressBar()->setProgress(static_cast<float>(progress) / static_cast<float>(filenames.size()));
            progress++;
        }

        //check SeriesInstanceUID of the file
        //if it fits: add this file to result
        try {
            string series = getMetaDataFromFile(*fileIterator, *dict_, "SeriesInstanceUID");
            if (series == seriesInstanceUID)
                result.push_back(*fileIterator);
        }
        catch(tgt::FileException) {
            //file is not a DICOM file... just ignore it
            LINFO("Ignoring file " + (*fileIterator) + " because it is not a DICOM file...");
        }
    }

    if (getProgressBar())
        getProgressBar()->hide();

    return result;
}

vector<vector<string> > GdcmVolumeReader::subdivideSeriesFilesByCustomDict(vector<string> fileNames, CustomDicomDict customDict) const {

    //iterate over all files and for each file over all subdivisionTags to get all value combinations in the series for the given files
    vector<string> subdivisionKeywords = *(customDict.getSubdivisionKeywords());

    //vector of value combinations and the corresponding files within this group
    vector<pair<string,vector<string> > > groups;

    std::vector<string>::const_iterator fileIterator;

    if (getProgressBar())
        getProgressBar()->setTitle("Subdividing files by CustomDict...");

    int itemused = 1;

    for (fileIterator = fileNames.begin(); fileIterator != fileNames.end(); ++fileIterator) {

        if (getProgressBar()) {
            getProgressBar()->setProgressMessage("Reading values of subdivision keywords from file \n" + (*fileIterator));
            getProgressBar()->setProgress(static_cast<float>(itemused) / static_cast<float>(fileNames.size()));
            itemused++;
        }

        //iterate over keywords and concatenate one string out of the values of these keywords
        string s = "";

        std::vector<string>::iterator keywordIterator;
        for (keywordIterator = subdivisionKeywords.begin(); keywordIterator != subdivisionKeywords.end(); ++keywordIterator) {
            //get the value from the file
            //string t = sf.ToString(getTagFromDictEntry(customDict.getDict()->getDictEntryByKeyword(*keywordIterator)));
            string t = getMetaDataFromFile(*fileIterator, *customDict.getDict(), *keywordIterator);
            t = trim(t, " ");

            if (!t.empty())
                s.append(t);
        }

        //check if the value combination is already present.
        //if so, insert the file into the corresponding vector.
        //else insert a new pair with the value combination and the file
        std::vector<pair<string, vector<string> > >::iterator combinationIterator;

        for (combinationIterator = groups.begin(); combinationIterator != groups.end(); ++combinationIterator) {
            if (combinationIterator->first == s) {
                combinationIterator->second.push_back(*fileIterator);
                break;
            }
        }

        if (combinationIterator == groups.end()) {
            pair<string, vector<string> > newPair;
            newPair.first = s;
            newPair.second.push_back(*fileIterator);
            groups.push_back(newPair);
        }

    }

    //now copy the groups into a result vector
    vector<vector<string> > resultVector;

    std::vector<pair<string, vector<string> > >::iterator combinationIterator;

    for (combinationIterator = groups.begin(); combinationIterator != groups.end(); ++combinationIterator) {
        //LERROR("Found " + itos(combinationIterator->second.size()) + "files for this group!");
        resultVector.push_back(combinationIterator->second);
    }

    if (getProgressBar())
        getProgressBar()->hide();

    return resultVector;
}

VolumeBase* GdcmVolumeReader::read(const VolumeURL& origin) {
    VolumeList* collection = read(origin.getURL());
    if (collection->size() > 1) {
        throw tgt::FileException("Could not load Volume, found more than one Volume: "+ origin.getPath(), origin.getPath());
    }
    else {
        if (getProgressBar())
            getProgressBar()->hide();
        return collection->first();
    }
}

VolumeList* GdcmVolumeReader::read(const string& url) {
    //check if fileInfoBuffer_ should be cleared (is the case when last buffer modification was at least 10 minutes ago)
    //get current time
    DateTime current = DateTime::now();
    double secondsSinceLastBufferMod = std::difftime(current.getTimestamp(),lastBufferMod_.getTimestamp());
    if ((!fileInfoBuffer_.empty()) && (secondsSinceLastBufferMod >= 600.0)) {
        LINFO("Clearing buffer (last modification was " + dtos(secondsSinceLastBufferMod / 60.0) + " minutes ago)");
        fileInfoBuffer_.clear();
        lastBufferMod_ = DateTime::now();
    }

    //if Standard Dictionary is not yet loaded: try to load
    if (!dict_) {
        LINFO("Trying to load Standard Dictionary...");
        loadStandardDict();
        LINFO("Standard Dictionary successfully loaded.");
    }

    VolumeURL origin;

    //check if this is a network path
    if (url.substr(0,12) == "dicom-scp://") {
        //get AE Title and Incoming Port Number
        string aet = dynamic_cast<const GdcmModule*>(VoreenApplication::app()->getModule("gdcm"))->getAeTitle();
        int incomingPort = dynamic_cast<const GdcmModule*>(VoreenApplication::app()->getModule("gdcm"))->getIncomingPortNumber();

        //create DicomNetworkConnector
        DicomNetworkConnector net(*dict_, aet, incomingPort);
        net.setProgressBar(getProgressBar());
        origin = VolumeURL(url);

        //construct path
        string path = constructLocalPathFromNetworkOrigin(origin);

        //create origin from local directory
        VolumeURL localOrigin("dicom",path,origin.getSearchString());

        VolumeBase* vHandle = 0;

        //check for files in tmp-directory and try to load from there
        if (!getFileNamesInDir(path).empty()) {
            LINFO("Found files in tmp-directory... trying to load from there.");
            try {
                vHandle = read(localOrigin);
            }
            catch (tgt::Exception t) {
                throw tgt::FileException("Couldn't load files in tmp-Directory: " + string(t.what()) + "... Please remove "+ path +" and try again.",path);
            }

            if (vHandle) {
                vHandle->setOrigin(origin);
                VolumeList* collection = new VolumeList();
                collection->add(vHandle);
                return collection;
            }
        }

        //else: create the directory and get the files by using C-MOVE
        LINFO("No local temporary version of Volume found... trying to use C-MOVE to retrieve images.");
        tgt::FileSystem::createDirectoryRecursive(path);

        //move the files
        try {
            net.moveSeries(origin,path);
        }
        catch (VoreenException v) {
            throw tgt::FileException(v.what(),"");
        }

        try {
            //try to load from there
            vHandle = read(localOrigin);
        }
        catch (tgt::Exception t) {
            throw tgt::FileException("Couldn't load files after C-MOVE: " + string(t.what()),path);
        }

        if (vHandle) {
            vHandle->setOrigin(origin);
            VolumeList* collection = new VolumeList();
            collection->add(vHandle);
            return collection;
        }
    }

    if (!(url.substr(0,8) == "dicom://")) {
        origin = VolumeURL("dicom", url);
    }
    else {
        origin = VolumeURL(url);
    }

    string fileName = origin.getPath();

    //Volume* volume;
    VolumeList* collection;

    //check if url is a file or a directory
    if (tgt::FileSystem::dirExists(fileName)) {
        //url is a directory
        collection = selectAndLoadDicomFiles(getFileNamesInDir(fileName), origin);
    }
    else if (tgt::FileSystem::fileExists(fileName)){
        //url is a file -> check if it is a Dicom file
        if (!isDicomFile(fileName)) {
            throw tgt::FileAccessException("Selected file is not a DICOM file or it is broken!", fileName);
        }
        else if (isDicomDir(fileName)) {
            //file is a DICOMDIR
            collection = readDicomDir(origin);
        }
        else {
            //file ist NOT a DICOMDIR -> read files in Directory
            LINFO("Loading files in " << tgt::FileSystem::dirName(fileName));

            //read SeriesInstanceUID from file if possible
            gdcm::Reader reader;
            reader.SetFileName(fileName.c_str());
            if (!reader.Read())
                throw tgt::FileAccessException("Could not read Dicom File", fileName);

            gdcm::StringFilter sf = gdcm::StringFilter();
            sf.SetFile(reader.GetFile());
            string seriesInstanceUID = sf.ToString(getTagFromDictEntry(dict_->getDictEntryByKeyword("SeriesInstanceUID")));
            seriesInstanceUID = trim(seriesInstanceUID, " ");

            //If there has not been a SeriesInstanceUID in the VolumeURL before: add it
            if (origin.getSearchParameter("SeriesInstanceUID").empty())
                origin.addSearchParameter("SeriesInstanceUID", seriesInstanceUID);

            //If the SeriesInstanceUID of the file and the one in the VolumeURL do not correspond: throw Exception
            if (!(seriesInstanceUID == origin.getSearchParameter("SeriesInstanceUID")))
                throw tgt::FileException("SeriesInstanceUID in VolumeURL does not correspond to file: " + origin.getPath(), origin.getPath());

            collection = selectAndLoadDicomFiles(getFileNamesInDir(tgt::FileSystem::dirName(fileName)), origin);
        }
    }
    else {
        //url is neither a file nor a directory -> error
        throw tgt::FileNotFoundException("GdcmVolumeReader: Unable to find ", fileName);
    }

    if (collection->empty()) {
       if (getProgressBar())
            getProgressBar()->hide();
       throw tgt::FileException("Found no volume to load! ("+origin.getPath()+")", origin.getPath());
    }
    if (getProgressBar())
        getProgressBar()->hide();
    return collection;
}

VolumeList* GdcmVolumeReader::selectAndLoadDicomFiles(const std::vector<std::string> &fileNames, const VolumeURL &vorigin) {

    VolumeURL origin(vorigin);

    if (fileNames.empty())
        throw tgt::FileNotFoundException("Directory is empty!", origin.getPath());

    //filter the given files, so that only readable DICOM files with one SeriesInstanceUID remain
    bool filter;
    if (origin.getSearchParameter("SeriesInstanceUID").empty())
        filter = false;
    else
        filter = true;

    vector<string> useFiles; //selected files to be used

    std::vector<std::string>::const_iterator i;
    int itemused = 1;

    if (getProgressBar())
        getProgressBar()->setTitle("Preselecting files");

    for (i = fileNames.begin(); i != fileNames.end(); i++) {

        if (getProgressBar()) {
            getProgressBar()->setProgressMessage("Checking SeriesInstanceUID of files...");
            getProgressBar()->setProgress(static_cast<float>(itemused) / static_cast<float>(fileNames.size()));
            itemused++;
        }

        if (!isDicomFile(*i)) {
            LINFO("Skipping file: " + (*i));
        }
        else {
            gdcm::Reader iReader;

            if (!filter) {
                //if there is no SeriesInstanceUID until now, get it from the first file
                //Read file
                iReader.SetFileName((*i).c_str());
                iReader.Read();

                //check if meta information is valid
                if (!iReader.GetFile().GetHeader().IsValid())
                    throw tgt::FileAccessException("GdcmVolumeReader: File Header not valid! ", (*i));

                //get SeriesInstanceUID
                string seriesInstanceUID = getMetaDataFromFile(*i, *dict_, "SeriesInstanceUID");

                if (seriesInstanceUID.empty())
                    LERROR("File has no SeriesInstanceUID: " << (*i));
                else {
                    seriesInstanceUID = trim(seriesInstanceUID, " ");
                    origin.addSearchParameter("SeriesInstanceUID", seriesInstanceUID);
                    filter = true;
                    useFiles.push_back(*i);
                    LINFO("Using first SeriesInstanceUID found: " + seriesInstanceUID);
                }
            }
            else {
                //check if this file belongs to the given series
                string seriesInstanceUID = getMetaDataFromFile(*i, *dict_, "SeriesInstanceUID");

                if (seriesInstanceUID == origin.getSearchParameter("SeriesInstanceUID"))
                    useFiles.push_back(*i);
            }
        }
    }

    if (getProgressBar())
        getProgressBar()->hide();

    //now all files in useFiles belong to the same SeriesInstanceUID
    if (useFiles.empty())
        throw tgt::FileNotFoundException("No file with corresponding SeriesInstanceUID could be found!", origin.getPath());

    return subdivideAndLoadDicomFiles(useFiles, origin);
}

VolumeList* GdcmVolumeReader::readDicomDir(const VolumeURL &origin) {

        string filepath = origin.getPath();
        string selectionSeriesInstanceUID = origin.getSearchParameter("SeriesInstanceUID", false);

        DicomDirParser parser(filepath);
        parser.open();

        LINFO("Loading DICOMDIR...");

        vector<string> filenames; //contains all image files
        vector<DicomInfo> seriesCollection; //contains all combinations of Patient -> Study -> Series in the DICOMDIR

        string patientId = ""; //current Patient ID
        string studyId = "";   //current Study ID
        string seriesId = "";  //current Series ID

        std::vector<VolumeList*> volumeLists; //contains the volumes to be returned

        if (getProgressBar())
            getProgressBar()->setTitle("Parsing DICOMDIR Structure");

        for (parser.setIteratorToFirstDirectoryRecordSequence(); !(parser.iteratorIsAtEnd()); parser.setIteratorToNextDirectoryRecordSequence()) {

            DirectoryRecordSequenceIterator sequenceIterator = parser.getPresentDirectoryRecordSequence();

            unsigned int itemused = 0;

            for (sequenceIterator.setIteratorToFirstElement(); !(sequenceIterator.iteratorIsAtEnd()); sequenceIterator.setIteratorToNextElement()) {

                if (getProgressBar()) {
                    getProgressBar()->setProgressMessage("Reading DirectoryRecordSequence ...");
                    getProgressBar()->setProgress(static_cast<float>(itemused) / static_cast<float>(sequenceIterator.getNumberOfItems()));
                    itemused++;
                }

                string itemName = sequenceIterator.getPresentItemName();
                itemName = trim(itemName, " ");

                if (itemName == "PATIENT") {
                    if (sequenceIterator.findTagInPresentItem(getTagFromDictEntry(dict_->getDictEntryByKeyword("PatientID"))))
                        patientId = sequenceIterator.getElementValue(getTagFromDictEntry(dict_->getDictEntryByKeyword("PatientID")));
                    patientId = trim(patientId, " ");
                }
                else if (itemName == "STUDY") {
                    if (sequenceIterator.findTagInPresentItem(getTagFromDictEntry(dict_->getDictEntryByKeyword("StudyInstanceUID"))))
                        studyId = sequenceIterator.getElementValue(getTagFromDictEntry(dict_->getDictEntryByKeyword("StudyInstanceUID")));
                    studyId = trim(studyId, " ");
                }
                else if (itemName == "SERIES") {
                    if (sequenceIterator.findTagInPresentItem(getTagFromDictEntry(dict_->getDictEntryByKeyword("SeriesInstanceUID"))))
                        seriesId = sequenceIterator.getElementValue(getTagFromDictEntry(dict_->getDictEntryByKeyword("SeriesInstanceUID")));

                    seriesId = trim(seriesId, " ");
                    //build new DicomInfo for this series
                    DicomInfo d;
                    d.setPatientId(patientId);
                    d.setStudyInstanceUID(studyId);
                    d.setSeriesInstanceUID(seriesId);

                    //get SERIES MODALITY
                    if (sequenceIterator.findTagInPresentItem(getTagFromDictEntry(dict_->getDictEntryByKeyword("Modality"))))
                        d.setModality(sequenceIterator.getElementValue(getTagFromDictEntry(dict_->getDictEntryByKeyword("Modality"))));

                    //only push back this Patient/Study/Series-Combination, if it fits to the given selectionSeriesInstanceUID:
                    //if a SelectionSeriesInstanceUID is given (ie. when loading a workspace), only this Series should be loaded
                    if (selectionSeriesInstanceUID.empty())
                        seriesCollection.push_back(d);
                    else if (selectionSeriesInstanceUID == d.getSeriesInstanceUID())
                        seriesCollection.push_back(d);
                }
                //PRIVATE and PRESENTATION tags should be handled in a different way
                //now they are simply ignored
                else if (itemName=="PRIVATE"){
                    //LINFO("Skipping Private Tag!");
                    //LINFO("--------------------------- ");
                }
                else if (itemName=="PRESENTATION"){
                    //LINFO("Skipping Presentation Tag!");
                    //LINFO("--------------------------- ");
                }
                else if (itemName=="IMAGE"){

                    //get PATH to the image
                    string filename;
                    if (sequenceIterator.findTagInPresentItem(getTagFromDictEntry(dict_->getDictEntryByKeyword("ReferencedFileID"))))
                        filename = sequenceIterator.getElementValue(getTagFromDictEntry(dict_->getDictEntryByKeyword("ReferencedFileID")));

                    filename = trim(filename, " ");

#ifndef WIN32
                size_t pos;
                while ((pos = filename.find_first_of('\\')) != string::npos) {
                    filename.replace(pos, 1, "/");
                }
#endif
                    //Get Path and try to add the file
                    string dir(string(filepath).substr(0, string(filepath).length() - string("DICOMDIR").length()));
                    dir.append(filename);

                    //remove space at the end of string
                    dir = trim(dir, " ");
                    //LINFO(dir);

                    if (!tgt::FileSystem::fileExists(dir)){
                        LERROR("File does not exist: " + dir);
                    }
                    else {
                        //files.push_back(dir);
                        filenames.push_back(dir);
                    }
                }
            } //end if(DirectoryRecordSequence)
        } //end iteration over DES

        //for each series: build volume
        std::vector<DicomInfo>::iterator seriesIterator;

        for(seriesIterator = seriesCollection.begin(); seriesIterator != seriesCollection.end(); ++seriesIterator){

            vector<string> seriesFilenames = getFilesInSeries(filenames, seriesIterator->getSeriesInstanceUID());

            if (getProgressBar())
                getProgressBar()->hide();

            //build the actual volume
            if (!seriesFilenames.empty()) {
                VolumeList* vc;

                try {
                    VolumeURL o(origin);
                    o.addSearchParameter("SeriesInstanceUID", seriesIterator->getSeriesInstanceUID());
                    vc = subdivideAndLoadDicomFiles(seriesFilenames, o);
                }
                catch(tgt::FileException f) {
                    string s(f.what());
                    LERROR(s);
                    vc = 0;
                    if (getProgressBar())
                        getProgressBar()->hide();
                }

                if (vc) {
                    volumeLists.push_back(vc);
                }
            }
        }

        //build VolumeList from all VolumeList-Objects
        VolumeList* volumeList = new VolumeList();

        vector<VolumeList*>::iterator iter = volumeLists.begin();
        while (iter != volumeLists.end()) {
              VolumeList* c = *iter;
              if (c) {
                  volumeList->add(c);
              }
              else {
                  throw tgt::FileException("GdcmVolumeReader: Unable to load " + origin.getPath(), origin.getPath());
              }
              iter++;
        }

        LINFO("Finished loading DICOMDIR.");

        if (getProgressBar())
            getProgressBar()->hide();

        return volumeList;
}



VolumeList* GdcmVolumeReader::subdivideAndLoadDicomFiles(const std::vector<std::string> &fileNames, const VolumeURL &origin) {

    VolumeList* vc = new VolumeList(); //the VolumeCollection to be returned

    //The given vector fileNames does now only contain the files of one SeriesInstanceUID, which is also a SearchParameter in the VolumeURL
    //Now look for CustomDicomDict-Files that may be relevant for the given files
    vector<CustomDicomDict> customDicts = getCustomDicts();

    //first file is used as reference
    string referenceFile = *fileNames.begin();

    if (getProgressBar())
        getProgressBar()->setTitle("Looking for CustomDicomDict-Files.");

    vector<CustomDicomDict>::iterator dictIt;
    int itemused = 1;

    for (dictIt = customDicts.begin(); dictIt != customDicts.end(); dictIt++) {

        if (getProgressBar()) {
            getProgressBar()->setProgressMessage("Checking CustomDicomDict-Files...");
            getProgressBar()->setProgress(static_cast<float>(itemused) / static_cast<float> (customDicts.size()));
            itemused++;
        }

        //check if the current CustomDicomDict is right for the given files
        const vector<pair<string, vector<string> > >* conditions = dictIt->getConditions();

        bool fits = true; //change if conditions do not apply

        //iterate over all conditions and check if they apply for the given files
        vector<pair<string, vector<string> > >::const_iterator condIterator;
        for (condIterator = conditions->begin(); condIterator != conditions->end(); condIterator++) {
            string value = getMetaDataFromFile(referenceFile, *dictIt->getDict(), condIterator->first);
            value = trim(value, " ");
            //iterate over possible values for the given tag
            fits = false;
            vector<string>::const_iterator valueIt;
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

        //if the CustomDicomDict fits for the given files: use it to subdivide the files into >= 1 volumes to be loaded separately
        if (fits) {
            if (getProgressBar())
                getProgressBar()->hide();

            //The Conditions of a specific CustomDicomDict match the given files -> use it for subdividing the series into several Volumes by using the subdivisionTags of the CustomDicomDict
            LINFO("Found CustomDicomDict for the given files...using it to subdivide the files of this SeriesInstanceUID...");

            vector<string> subdivisionKeywords = *(dictIt->getSubdivisionKeywords());
            const DicomDict customDict = *(dictIt->getDict());

            vector<vector<string> > divisionGroups = subdivideSeriesFilesByCustomDict(fileNames, *dictIt);

            vector<vector<string> >::iterator groupIterator;

            for (groupIterator = divisionGroups.begin(); groupIterator != divisionGroups.end(); ++groupIterator) {

                if (!groupIterator->empty()) {

                    bool ignore = false;

                    std::vector<string>::iterator keywordIterator;
                    for (keywordIterator = subdivisionKeywords.begin(); keywordIterator != subdivisionKeywords.end(); ++keywordIterator) {
                        //get the value from the first file of the group
                        string t = getMetaDataFromFile(*groupIterator->begin(), customDict, *keywordIterator);
                        t = trim(t, " ");

                        //if this value does not fit to the one specified in the VolumeURL: do not use this group
                        if (!(origin.getSearchParameter(*keywordIterator).empty() || (origin.getSearchParameter(*keywordIterator) == t))) {
                            ignore = true; //do not use this group
                            break;
                        }

                    }

                    if (ignore)
                        continue; //check out the next group

                    Volume* vh;

                    try {
                        vh = readDicomFiles(*groupIterator, origin);
                    }
                    catch(tgt::FileException f) {
                        string s(f.what());
                        LERROR(s);
                        vh = 0;
                    }

                    if (vh) {
                        //use the first file as a reference file to add additional MetaInformation and
                        //SearchParameters to the Volume and its VolumeURL as specified in the CustomDicomDict used
                        VolumeURL o(vh->getOrigin());

                        //add SearchParameter values to VolumeURL for subdivisionKeywords
                        vector<string>::iterator keyIterator;
                        for (keyIterator = subdivisionKeywords.begin(); keyIterator != subdivisionKeywords.end(); ++keyIterator) {
                            string value = getMetaDataFromFile(*groupIterator->begin(), customDict, *keyIterator);
                            value = trim(value, " ");
                            if (!value.empty())
                                o.addSearchParameter(*keyIterator, value);
                        }
                        vh->setOrigin(o);

                        LINFO("Additional Meta Information in CustomDicomDict:");
                        //add all Tags that are found in the Dictionary of this CustomDicomDict to the MetaInformation of the Volume, if the attribute metaData is set to true and if the value is not empty
                        setMetaDataFromDict(&(vh->getMetaDataContainer()), &customDict, *groupIterator->begin());

                        //print out Meta Information:
                        const vector<string> keys = customDict.getKeywordVector();
                        for (unsigned int i = 0; i < keys.size(); i++) {
                            MetaDataBase* m = vh->getMetaDataContainer().getMetaData(keys[i]);
                            if (m) {
                                LINFO(keys[i]+": "+ m->toString());
                            }
                        }

                        vc->add(vh);
                    }
                }
            }

            if (getProgressBar())
                getProgressBar()->hide();
            return vc;
        }
    }

    //if no CustomDicomDict is right: just load the files
    Volume* vh = readDicomFiles(fileNames, origin);
    if (vh) {
        vh->setOrigin(origin);
        vc->add(vh);
    }

    if (getProgressBar())
        getProgressBar()->hide();

    return vc;
}

Volume* GdcmVolumeReader::readDicomFiles(const vector<string> &fileNames, const VolumeURL &origin) {

    vector<pair<string, double> > slices; //filename, distance

    //take first file as reference, since all files belong to the same SeriesInstanceUID
    gdcm::ImageReader reader;
    reader.SetFileName(fileNames.begin()->c_str());
    reader.Read();
    gdcm::StringFilter sf;
    sf.SetFile(reader.GetFile());

    //Get Meta Information
    info_.setSeriesInstanceUID(sf.ToString(getTagFromDictEntry(dict_->getDictEntryByKeyword("SeriesInstanceUID"))));
    info_.setStudyInstanceUID(sf.ToString(getTagFromDictEntry(dict_->getDictEntryByKeyword("StudyInstanceUID"))));
    info_.setSeriesDescription(sf.ToString(getTagFromDictEntry(dict_->getDictEntryByKeyword("SeriesDescription"))));
    info_.setStudyDescription(sf.ToString(getTagFromDictEntry(dict_->getDictEntryByKeyword("StudyDescription"))));
    info_.setModality(sf.ToString(getTagFromDictEntry(dict_->getDictEntryByKeyword("Modality"))));
    info_.setPatientName(sf.ToString(getTagFromDictEntry(dict_->getDictEntryByKeyword("PatientName"))));
    info_.setPatientId(sf.ToString(getTagFromDictEntry(dict_->getDictEntryByKeyword("PatientID"))));

    info_.setIntercept(static_cast<float>(reader.GetImage().GetIntercept()));
    info_.setSlope(static_cast<float>(reader.GetImage().GetSlope()));
    info_.setRescaleType(sf.ToString(getTagFromDictEntry(dict_->getDictEntryByKeyword("RescaleType"))));

    scalarType_ = reader.GetImage().GetPixelFormat().GetScalarType();

    //for checks, if pixel data needs to be rescaled
    bool slopeDiffers = false;
    bool interceptDiffers = false;

    //get image related information
    //get image dimensions
    const unsigned int* dimensions = reader.GetImage().GetDimensions();

    if (dimensions) {
        info_.setDx(dimensions[0]);
        info_.setDy(dimensions[1]);
    }
    else
        throw tgt::FileException("No Dimensions found in Image File.", *fileNames.begin());

    //get bits stored
    info_.setBitsStored(reader.GetImage().GetPixelFormat().GetBitsStored());
    //get samples per pixel
    info_.setSamplesPerPixel(reader.GetImage().GetPixelFormat().GetSamplesPerPixel());

    LINFO("    Size: " << info_.getDx() << "x" << info_.getDy() << ", " << info_.getBitsStored()* info_.getSamplesPerPixel() << " bits");

    //get spacing for x and y
    const double* spacing = reader.GetImage().GetSpacing();
    info_.setXSpacing(spacing[0]);
    info_.setYSpacing(spacing[1]);

    //get ImageOrientation
    tgt::dvec3 xOrientationPatient;
    tgt::dvec3 yOrientationPatient;
    const double* directions = reader.GetImage().GetDirectionCosines();
    xOrientationPatient.x = directions[0];
    xOrientationPatient.y = directions[1];
    xOrientationPatient.z = directions[2];
    yOrientationPatient.x = directions[3];
    yOrientationPatient.y = directions[4];
    yOrientationPatient.z = directions[5];

    info_.setXOrientationPatient(xOrientationPatient);
    info_.setYOrientationPatient(yOrientationPatient);

    //calculate slice normal
    tgt::dvec3 sliceNormal = tgt::cross(xOrientationPatient,yOrientationPatient);

    info_.setSliceNormal(sliceNormal);

    //Read ImagePositionPatient from all files and calculate distance along slice normal
    vector<string>::const_iterator it_files = fileNames.begin();
    int i = 0;
    if (getProgressBar() && !fileNames.empty())
        getProgressBar()->setTitle("Loading DICOM Data Set");

    while (it_files != fileNames.end()) {
        if (getProgressBar()) {
            getProgressBar()->setProgressMessage("Reading slice '" + tgt::FileSystem::fileName(*it_files) + "' ...");
            getProgressBar()->setProgress(static_cast<float>(i) / static_cast<float>(fileNames.size()));
        }
        i++;

        //Read file
        gdcm::ImageReader reader;
        reader.SetFileName((*it_files).c_str());

        //check if file could be read -> should always be the case due to earlier checks
        if (!reader.Read()){
            if (getProgressBar())
                getProgressBar()->hide();
            throw tgt::FileAccessException("Could not read File which should be readable!", (*it_files));
        }

        //check if meta information is valid -> should always be the case due to earlier checks
        if (!reader.GetFile().GetHeader().IsValid()){
            if (getProgressBar())
                getProgressBar()->hide();
            throw tgt::FileAccessException("GdcmVolumeReader: File Header not valid! ", (*it_files));
        }

        gdcm::StringFilter sf = gdcm::StringFilter();
        sf.SetFile(reader.GetFile());

        //get position of image
        tgt::dvec3 imagePositionPatient;
        const double* origin = reader.GetImage().GetOrigin();
        imagePositionPatient.x = origin[0];
        imagePositionPatient.y = origin[1];
        imagePositionPatient.z = origin[2];

        //calculate distance
        double dist = tgt::dot(info_.getSliceNormal(),imagePositionPatient);

        //add image
        slices.push_back(make_pair(*it_files,dist));

        it_files++;
    }


    if (slices.size() == 0) {
        if (getProgressBar())
            getProgressBar()->hide();
        throw tgt::CorruptedFileException("Path does not contain any DICOM slices");
    }

    if (slices.size() == 1) {
        //check, if this file is a multiframe image
        if (reader.GetImage().GetNumberOfDimensions() == 3) {
            //get number of frames in this image
            string frames = sf.ToString(getTagFromDictEntry(dict_->getDictEntryByKeyword("NumberOfFrames")));
            if (!frames.empty())
                info_.setDz(stoi(frames));
            else {
                info_.setDz(1);
            }

            info_.setNumberOfFrames(info_.getDz());

            //get Z Spacing
            string sliceSpacing = sf.ToString(getTagFromDictEntry(dict_->getDictEntryByKeyword("SpacingBetweenSlices")));
            if (!sliceSpacing.empty())
                info_.setZSpacing(stod(sliceSpacing));
            else {
                info_.setZSpacing(1.0);
            }

        }
        else if (reader.GetImage().GetNumberOfDimensions() == 2){
            //not a multiframe image
            info_.setDz(1);
            info_.setZSpacing(1.0);
            info_.setNumberOfFrames(1);
        }
        else {
            if (getProgressBar())
                getProgressBar()->hide();
            throw tgt::FileException("Unexpected Number of Dimensions in Image File (Not supported): " + std::to_string(reader.GetImage().GetNumberOfDimensions()), slices.begin()->first);
        }

        //check samples per pixel and warn, if != 1
        /*if (info_.getSamplesPerPixel() != 1) {
            LWARNING("Unsupported Pixel Format: " + itos(info_.getSamplesPerPixel()) + " Samples per Pixel instead of 1! Might lead to unexpected results.");
        }*/


    }
    else
        info_.setNumberOfFrames(1);

    //sort slices by their distance from the origin, calculate Z spacing and do some additional checks
    if (slices.size() > 1) {
        //sort slices by distance from origin
        std::sort(slices.begin(), slices.end(), slices_cmp_dist);

        //check, if all images files are of dimension 2
        //also check if rescale intercept and slope are uniform and samples per pixel = 1 in all images
        int samplesPerPixel = 1;
        for (unsigned int i = 0; i < slices.size(); i++) {
            gdcm::ImageReader imR;
            imR.SetFileName(slices[i].first.c_str());
            imR.Read();
            if (imR.GetImage().GetNumberOfDimensions() != 2) {
                if (getProgressBar())
                    getProgressBar()->hide();
                throw tgt::FileException("Image file has unexpected Dimensions (Multiple slices are required to have Dimension 2): " + std::to_string(imR.GetImage().GetNumberOfDimensions()), slices[i].first);
            }

            //check if PixelRepresentation is uniform
            if (imR.GetImage().GetPixelFormat().GetScalarType() != scalarType_) {
                if (getProgressBar())
                    getProgressBar()->hide();
                throw tgt::FileException("Image files do not have uniform scalar type!");
            }

            //check if rescale slope and intercept are the same for alle images, otherwise: warning
            if (info_.getIntercept() != static_cast<float>(imR.GetImage().GetIntercept())) {
                interceptDiffers = true;
                //LERROR("Intercept: " + ftos(static_cast<float>(imR.GetImage().GetIntercept())));
            }
            if (info_.getSlope() != static_cast<float>(imR.GetImage().GetSlope())) {
                slopeDiffers = true;
                //LERROR("Slope: " + ftos(static_cast<float>(imR.GetImage().GetSlope())));
            }

            //check if samples per pixel are uniformly = 1
            if ((imR.GetImage().GetPixelFormat().GetSamplesPerPixel() != samplesPerPixel) && (imR.GetImage().GetPixelFormat().GetSamplesPerPixel() != 1)) {
                samplesPerPixel = imR.GetImage().GetPixelFormat().GetSamplesPerPixel();
                LWARNING("Found image files with unsupported Pixel Format: " + itos(samplesPerPixel) + " Samples per Pixel instead of 1! Might lead to unexpected results.");
                info_.setSamplesPerPixel(samplesPerPixel);
            }
        }

        if (interceptDiffers)
            LWARNING("Rescale Intercept differs within the image files!");

        if (slopeDiffers)
            LWARNING("Rescale Slope differs within the image files!");

        //calculate Z-Spacing
        info_.setZSpacing(slices[1].second - slices[0].second);
        if (info_.getZSpacing() == 0){
            if (getProgressBar())
                getProgressBar()->hide();
            throw tgt::FileException("Slice Spacing is 0: Found two or more Slices with the same Position! (Either not a Volume or Slices have to be subdivided by additional Attributes)");
        }

        //check, if slice spacing remains constant (with 10% tolerance)
        for (unsigned int i = 0; i < slices.size()-1; i++) {
            if (((slices[i+1].second - slices[i].second) < 0.9*info_.getZSpacing()) || ((slices[i+1].second - slices[i].second) > 1.1*info_.getZSpacing())) {
                //LERROR("Spacing between slices is not steady (10% tolerance)!");
                if ((slices[i+1].second - slices[i].second) == 0) {
                    if (getProgressBar())
                        getProgressBar()->hide();
                    throw tgt::FileException("Slice Spacing is 0: Found two or more Slices with the same Position! (Either not a Volume or Slices have to be subdivided by additional Attributes)");
                }
                if (getProgressBar())
                    getProgressBar()->hide();
                //either throw exception or display warning
                bool abort = dynamic_cast<const GdcmModule*>(VoreenApplication::app()->getModule("gdcm"))->ignoreSliceSpacing();
                if (abort)
                    throw tgt::FileException("Slice Spacing is not steady (differs > 10% Tolerance)! To ignore this set the option to ignore slice spacing differences in the module options and try to load the data set again.");
                else {
                    LWARNING("Slice Spacing is not steady (differs > 10% Tolerance)! The data set might be missing one or more slices!");
                    break;
                }
            }
        }

        info_.setDz(static_cast<int>(slices.size()));
    }

    //get position of first image to calculate offset
    gdcm::ImageReader yar;
    yar.SetFileName(slices[0].first.c_str());
    yar.Read();
    sf.SetFile(yar.GetFile()); //set file to gdcm::StringFilter for later use

    tgt::dvec3 offset;
    const double* volumeOrigin = yar.GetImage().GetOrigin();
    offset.x = volumeOrigin[0];
    offset.y = volumeOrigin[1];
    offset.z = volumeOrigin[2];

    info_.setOffset(offset);

    //get pixel representation
    info_.setPixelRepresentation(yar.GetImage().GetPixelFormat().GetPixelRepresentation());

    LINFO("We have " << info_.getDz() << " slices. [" << info_.getDx() << "x" << info_.getDy() << "]");
    LINFO("Spacing: (" << info_.getXSpacing() << "; " << info_.getYSpacing() << "; " << info_.getZSpacing() << ")");

    //Determine scalar base type of the data set if not already known (and if possible)
    if (scalarType_ == gdcm::PixelFormat::UNKNOWN) {

        LWARNING("Pixel Format: Scalar Type in DICOM files is UNKNOWN or format is not well defined... trying to compute the right format.");

        if (info_.getPixelRepresentation()) { //signed
            switch (info_.getBitsStored()) {
            case 8:
                scalarType_ = gdcm::PixelFormat::INT8;
                break;
            case 9:
            case 10:
            case 11:
            case 12:
            case 13:
            case 14:
            case 15:
            case 16:
                scalarType_ = gdcm::PixelFormat::INT16;
                break;
            case 24:
            case 32:
                scalarType_ = gdcm::PixelFormat::INT32;
                break;
            default:
                break;
            }
        }
        else {  //unsigned
            switch (info_.getBitsStored()) {
            case 8:
                scalarType_ = gdcm::PixelFormat::UINT8;
                break;
            case 9:
            case 10:
            case 11:
            case 12:
            case 13:
            case 14:
            case 15:
            case 16:
                scalarType_ = gdcm::PixelFormat::UINT16;
                break;
            case 24:
            case 32:
                scalarType_ = gdcm::PixelFormat::UINT32;
                break;
            default:
                break;
            }
        }
    }

    bool rwmDiffers = slopeDiffers || interceptDiffers;
    info_.setRwmDiffers(rwmDiffers);

    if (rwmDiffers) {
        if ((scalarType_ != gdcm::PixelFormat::UNKNOWN) && (info_.getSamplesPerPixel() == 1)) {
            LWARNING("Rescaling of pixel data necessary due to different scaling");
            //compute correct rescale slope and intercept values
            computeCorrectRescaleValues(slices);
        }
        else
            LWARNING("Rescaling of pixel data not possible due to unknown data type. Result may be broken.");
    }

    //write data base type as string to info object
    switch (scalarType_) {
        case gdcm::PixelFormat::UINT8:
            info_.setBaseType("uint8"); break;
        case gdcm::PixelFormat::INT8:
            info_.setBaseType("int8"); break;
        case gdcm::PixelFormat::UINT12:
            info_.setBaseType("uint16"); LWARNING("unsigned int (12 bit) is interpreted as 16 bit - might lead to unexpected results"); break;
        case gdcm::PixelFormat::INT12:
            info_.setBaseType("int16"); LWARNING("signed int (12 bit) is interpreted as 16 bit - might lead to unexpected results"); break;
        case gdcm::PixelFormat::UINT16:
            info_.setBaseType("uint16"); break;
        case gdcm::PixelFormat::INT16:
            info_.setBaseType("int16"); break;
        case gdcm::PixelFormat::UINT32:
            info_.setBaseType("uint32"); break;
        case gdcm::PixelFormat::INT32:
            info_.setBaseType("int32"); break;
        case gdcm::PixelFormat::FLOAT32:
            info_.setBaseType("float"); break;
        case gdcm::PixelFormat::FLOAT64:
            info_.setBaseType("double"); break;
        case gdcm::PixelFormat::FLOAT16:
            throw tgt::Exception("Data format: 16 bit float not supported");
        case gdcm::PixelFormat::UNKNOWN:
            throw tgt::Exception("Unknown data format!");
        default:
            throw tgt::Exception("Data format could not be computed!");
    }

    //determine format by using samples per pixel and base type
    if (info_.getSamplesPerPixel() == 1)
        info_.setFormat(info_.getBaseType());
    else if (info_.getSamplesPerPixel() < 1 || info_.getSamplesPerPixel() > 4)
        throw tgt::Exception("Unsupported format: " + itos(info_.getSamplesPerPixel()) + " samples per pixel");
    else {
        //build format as vector type
        std::stringstream s;
        s << "Vector" << info_.getSamplesPerPixel() << "(" << info_.getBaseType() << ")";
        info_.setFormat(s.str());
    }

    //determine bytes per voxel
    VolumeFactory volumeFac;
    info_.setBytesPerVoxel(volumeFac.getBytesPerVoxel(info_.getFormat()));

    //copy slices into new vector for disk representation
    std::vector<std::string> sliceFilenamesOnly(slices.size());
    for (size_t i = 0; i < slices.size(); ++i)
        sliceFilenamesOnly.at(i) = slices.at(i).first;

    Volume* vh;

    if (info_.getNumberOfFrames() > 1) {
        LWARNING("Multiframe DICOM file... loading VolumeRAM representation instead of VolumeDisk");
        VolumeRAM* ram = loadMultiframeDicomFile(info_, sliceFilenamesOnly);
        vh = new Volume(ram,
            tgt::vec3(static_cast<float>(info_.getXSpacing()),static_cast<float>(info_.getYSpacing()),static_cast<float>(info_.getZSpacing())), tgt::vec3(0.f));
    }
    else {
        //build volume disk representation
        VolumeDiskDicom* diskVolume = new VolumeDiskDicom(info_.getFormat(), tgt::svec3(info_.getDx(), info_.getDy(), info_.getDz()), info_, sliceFilenamesOnly);

        //build: VolumeDisk -> Volume
        vh = new Volume(diskVolume,
            tgt::vec3(static_cast<float>(info_.getXSpacing()),static_cast<float>(info_.getYSpacing()),static_cast<float>(info_.getZSpacing())), tgt::vec3(0.f));
    }

    //get volume origin in world coordinates:
    tgt::vec3 os(static_cast<float>(info_.getOffset().x), static_cast<float>(info_.getOffset().y), static_cast<float>(info_.getOffset().z));

    //construct PhysicalToWorld-Matrix for correct positioning of Volume
    tgt::dmat4 ptw(tgt::dvec4(info_.getXOrientationPatient(), os.x),
            tgt::dvec4(info_.getYOrientationPatient(), os.y), tgt::dvec4(info_.getSliceNormal(), os.z), tgt::dvec4(0,0,0,1));
    vh->setPhysicalToWorldMatrix(ptw);

    //set Modality and VolumeURL
    vh->setModality(Modality(info_.getModality()));
    vh->setOrigin(origin);

    //setting RealWorldMapping for denormalization and rescaling
    RealWorldMapping denormalize;
    try {
       denormalize = RealWorldMapping::createDenormalizingMapping(info_.getBaseType());
    } catch(tgt::UnsupportedFormatException) {
       LWARNING("Unknown format: " + info_.getBaseType());
    }

    RealWorldMapping rwm(info_.getSlope(), info_.getIntercept(), info_.getRescaleType());

    vh->setRealWorldMapping(RealWorldMapping::combine(denormalize, rwm));

    /*set Meta Information to Volume
      ------------------------------------*/

    LINFO("Setting Meta Information:");
    //add all Tags that are found in the Standard Dictionary to the MetaInformation of the Volume,
    //if the attribute metaData is set to true and if the value is not empty
    setMetaDataFromDict(&(vh->getMetaDataContainer()), dict_, slices[0].first);

    //Set VolumeDateTime, depending on the tags found in the file
    MetaDataBase* volumeDateTime = constructVolumeDateTime(dict_, &sf);
    if (volumeDateTime)
        vh->getMetaDataContainer().addMetaData("VolumeDateTime", volumeDateTime);


    //print out Meta Information:
    const vector<string> keys = dict_->getKeywordVector();
    for (unsigned int i = 0; i < keys.size(); i++) {
        MetaDataBase* m = vh->getMetaDataContainer().getMetaData(keys[i]);
        if (m) {
            LINFO(keys[i]+": "+ m->toString());
        }
    }

    //set file list meta data
    /*vector<string> fileList;
    std::vector<pair<string, double> >::iterator sliceIt;
    for (sliceIt = slices.begin(); sliceIt != slices.end(); ++sliceIt) {
        fileList.push_back(sliceIt->first);
    }
    FileListMetaData* fileListMetaData = new FileListMetaData(fileList);
    vh->getMetaDataContainer().addMetaData("FileList", fileListMetaData);*/

    if (getProgressBar())
        getProgressBar()->hide();

    return vh;
}

VolumeRAM* GdcmVolumeReader::loadMultiframeDicomFile(const DicomInfo& info, const std::vector<std::string>& sliceFiles) {
    if (sliceFiles.size() != 1)
        throw tgt::FileException("Multiple files containing multiframe DICOM data not supported");

    //now build the actual volume
    VolumeFactory volumeFac;
    VolumeRAM* dataset = 0;
    try {
        dataset = volumeFac.create(info.getFormat(), tgt::svec3(info.getDx(), info.getDy(), info.getDz()));
    }
    catch (std::exception& e) {
        LERROR(e.what());
        if (getProgressBar())
            getProgressBar()->hide();
        throw e;
    }

    // Now read the actual slices from the files
    LINFO("Building volume from multiframe DICOM file...");

    if (getProgressBar()) {
        getProgressBar()->setProgressMessage("Loading slice '" + tgt::FileSystem::fileName(sliceFiles.at(0)) + "' ...");
        getProgressBar()->setProgress(0.f);
    }

    size_t posScalar = 0;
    int slicesize = loadSlice(reinterpret_cast<char*>(dataset->getData()), sliceFiles.at(0), posScalar, info);

    if (getProgressBar())
        getProgressBar()->hide();

    if (slicesize == 0) {
        //obviously an error in loadSlice method
        delete dataset;
        throw tgt::FileException("Failed to read Pixel data.", sliceFiles.at(0));
    }

    LINFO("Building volume complete.");

    return dataset;
}


VolumeRAM* GdcmVolumeReader::loadDicomSlices(DicomInfo info, std::vector<std::string> sliceFiles) {

    if (sliceFiles.size() < 1)
        throw tgt::FileException("No slice files to build volume!");

    //LINFO("Building volume...");

    //build the actual volume
    VolumeFactory volumeFac;
    VolumeRAM* dataset = 0;
    try {
        //create data set for the slices
        dataset = volumeFac.create(info.getFormat(), tgt::svec3(info.getDx(), info.getDy(), sliceFiles.size()));
    }
    catch (std::exception& e) {
        LERROR(e.what());
        if (getProgressBar())
            getProgressBar()->hide();
        throw e;
    }

    //LINFO("Reading slice data from " << sliceFiles.size() << " files...");

    size_t posScalar = 0;

    std::vector<std::string>::iterator it_slices = sliceFiles.begin();
    size_t i = 0;
    while (it_slices != sliceFiles.end()) {
        if (getProgressBar()) {
            getProgressBar()->setProgressMessage("Loading slice '" + tgt::FileSystem::fileName(*it_slices) + "' ...");
            getProgressBar()->setProgress(static_cast<float>(i) / static_cast<float>(sliceFiles.size()));
            i++;
        }

        int slicesize = loadSlice(reinterpret_cast<char*>(dataset->getData()), *it_slices, posScalar, info);

        if (slicesize == 0) {
            //obviously an error in loadSlice method
            delete[] dataset;
            if (getProgressBar())
                getProgressBar()->hide();
            throw tgt::FileException("Failed to read Pixel data.", *it_slices);
        }

        posScalar += slicesize;
        it_slices++;
    }
    if (getProgressBar())
        getProgressBar()->hide();

    //LINFO("Building volume complete.");

    return dataset;
}

void GdcmVolumeReader::computeCorrectRescaleValues(std::vector<pair<string, double> > slices) {
    //for every slice, it is assumed that the whole domain of the scalar type values is used
    float rwmMin = std::numeric_limits<float>::max();
    float rwmMax = std::numeric_limits<float>::min();

    float dataTypeMin = 0.f;
    float dataTypeMax = 0.f;;

    //get max and min real world values
    std::vector<pair<string, double> >::iterator it_slices = slices.begin();

    int i = 0;
    while (it_slices != slices.end()) {
        if (getProgressBar()) {
            getProgressBar()->setProgressMessage("Calculating pixel rescaling...");
            getProgressBar()->setProgress(static_cast<float>(i) / static_cast<float>(slices.size()));
        }
        i++;

        gdcm::ImageReader reader;
        reader.SetFileName(it_slices->first.c_str());
        reader.Read();
        float slope = static_cast<float>(reader.GetImage().GetSlope());
        float intercept = static_cast<float>(reader.GetImage().GetIntercept());

        dataTypeMin = static_cast<float>(reader.GetImage().GetPixelFormat().GetMin());
        dataTypeMax = static_cast<float>(reader.GetImage().GetPixelFormat().GetMax());
        float sliceMin = dataTypeMin * slope + intercept;
        float sliceMax = dataTypeMax * slope + intercept;

        rwmMin = std::min(rwmMin, sliceMin);
        rwmMax = std::max(rwmMax, sliceMax);

        it_slices++;
    }

    //Calculate correct global slope and intercept values
    float globalSlope = (rwmMax - rwmMin) / (dataTypeMax - dataTypeMin + 1);
    float globalIntercept = rwmMin - dataTypeMin;

    info_.setSlope(globalSlope);
    info_.setIntercept(globalIntercept);

    if (getProgressBar())
        getProgressBar()->hide();
}


int GdcmVolumeReader::loadSlice(char* dataStorage, const std::string& fileName, size_t posScalar, DicomInfo info){

    gdcm::ImageReader reader;
    reader.SetFileName(fileName.c_str());

    if (!reader.Read()){
        LERROR("Error loading file " << fileName);
        return 0;
    }

    size_t dataLength = (static_cast<size_t>(info.getDx()) * static_cast<size_t>(info.getDy())
                        * static_cast<size_t>(info.getBytesPerVoxel()) * static_cast<size_t>(info.getNumberOfFrames()));

    if (reader.GetImage().GetBufferLength() != dataLength){
        LERROR("Failed to read Pixel data from file " << fileName << " because of unexpected Buffer Length!");
        return 0;
    }

    gdcm::PixelFormat scalarType = baseTypeStringToGdcm(info.getBaseType());

    //get pixel data
    if (info.rwmDiffers() && (scalarType != gdcm::PixelFormat::UNKNOWN) && (info.getSamplesPerPixel() == 1)) {
        //if rescale intercept and slope differ: recalculate the scalar values so that these fit the correct rescaling
        float slope = static_cast<float>(reader.GetImage().GetSlope());
        float intercept = static_cast<float>(reader.GetImage().GetIntercept());

        float nSlope = slope / info.getSlope();
        float nIntercept = intercept - info.getIntercept();

        //save the original buffer temporarily
        char* sliceScalars = new char[dataLength];
        reader.GetImage().GetBuffer(sliceScalars);


        //rescale and write buffer: use either gdcm rescaling or naive implementation (slower but might work for some data sets that the other method might not)
        bool useGdcmRescaling = dynamic_cast<const GdcmModule*>(VoreenApplication::app()->getModule("gdcm"))->useGdcmRescaling();
        if (useGdcmRescaling) {
            gdcm::Rescaler ir;
            ir.SetSlope(static_cast<double>(nSlope));
            ir.SetIntercept(static_cast<double>(nIntercept));
            ir.SetPixelFormat(scalarType);
            ir.SetUseTargetPixelType(true);
            ir.SetTargetPixelType(scalarType);
            ir.SetMinMaxForPixelType(static_cast<double>(gdcm::PixelFormat(scalarType).GetMin()), static_cast<double>(gdcm::PixelFormat(scalarType).GetMax()));
            ir.Rescale(&dataStorage[posScalar * info.getBytesPerVoxel()], sliceScalars, dataLength);
        }
        else {
            //do a dynamic cast depending on the data type and rescale every value
            switch (scalarType) {
                case gdcm::PixelFormat::UINT8:
                    {
                        uint8_t* typeScalars = reinterpret_cast<uint8_t*>(sliceScalars);
                        size_t length = dataLength;

                        for (size_t voxelIndex = 0; voxelIndex < length; ++voxelIndex) {
                            float valueToRescale = static_cast<float>(typeScalars[voxelIndex]);
                            float rescaledValue = ((valueToRescale * slope + intercept) - info.getIntercept()) / info.getSlope();
                            typeScalars[voxelIndex] = static_cast<uint8_t>(rescaledValue);
                        }
                    }
                    break;
                case gdcm::PixelFormat::INT8:
                    {
                        int8_t* typeScalars = reinterpret_cast<int8_t*>(sliceScalars);
                        size_t length = dataLength;

                        for (size_t voxelIndex = 0; voxelIndex < length; ++voxelIndex) {
                            float valueToRescale = static_cast<float>(typeScalars[voxelIndex]);
                            float rescaledValue = ((valueToRescale * slope + intercept) - info.getIntercept()) / info.getSlope();
                            typeScalars[voxelIndex] = static_cast<int8_t>(rescaledValue);
                        }
                    }
                    break;
                case gdcm::PixelFormat::UINT12:
                case gdcm::PixelFormat::UINT16:
                    {
                        uint16_t* typeScalars = reinterpret_cast<uint16_t*>(sliceScalars);
                        size_t length = dataLength / 2;

                        for (size_t voxelIndex = 0; voxelIndex < length; ++voxelIndex) {
                            float valueToRescale = static_cast<float>(typeScalars[voxelIndex]);
                            float rescaledValue = ((valueToRescale * slope + intercept) - info.getIntercept()) / info.getSlope();
                            typeScalars[voxelIndex] = static_cast<uint16_t>(rescaledValue);
                        }
                    }
                    break;
                case gdcm::PixelFormat::INT12:
                case gdcm::PixelFormat::INT16:
                    {
                        int16_t* typeScalars = reinterpret_cast<int16_t*>(sliceScalars);
                        size_t length = dataLength / 2;

                        for (size_t voxelIndex = 0; voxelIndex < length; ++voxelIndex) {
                            float valueToRescale = static_cast<float>(typeScalars[voxelIndex]);
                            float rescaledValue = ((valueToRescale * slope + intercept) - info.getIntercept()) / info.getSlope();
                            typeScalars[voxelIndex] = static_cast<int16_t>(rescaledValue);
                        }
                    }
                    break;
                case gdcm::PixelFormat::UINT32:
                    {
                        uint32_t* typeScalars = reinterpret_cast<uint32_t*>(sliceScalars);
                        size_t length = dataLength / 4;

                        for (size_t voxelIndex = 0; voxelIndex < length; ++voxelIndex) {
                            float valueToRescale = static_cast<float>(typeScalars[voxelIndex]);
                            float rescaledValue = ((valueToRescale * slope + intercept) - info.getIntercept()) / info.getSlope();
                            typeScalars[voxelIndex] = static_cast<uint32_t>(rescaledValue);
                        }
                    }
                    break;
                case gdcm::PixelFormat::INT32:
                    {
                        int32_t* typeScalars = reinterpret_cast<int32_t*>(sliceScalars);
                        size_t length = dataLength / 4;

                        for (size_t voxelIndex = 0; voxelIndex < length; ++voxelIndex) {
                            float valueToRescale = static_cast<float>(typeScalars[voxelIndex]);
                            float rescaledValue = ((valueToRescale * slope + intercept) - info.getIntercept()) / info.getSlope();
                            typeScalars[voxelIndex] = static_cast<int32_t>(rescaledValue);
                        }
                    }
                    break;
                case gdcm::PixelFormat::FLOAT16:
                case gdcm::PixelFormat::FLOAT32:
                    {
                        float* typeScalars = reinterpret_cast<float*>(sliceScalars);
                        size_t length = dataLength / 4;

                        for (size_t voxelIndex = 0; voxelIndex < length; ++voxelIndex) {
                            float valueToRescale = typeScalars[voxelIndex];
                            float rescaledValue = ((valueToRescale * slope + intercept) - info.getIntercept()) / info.getSlope();
                            typeScalars[voxelIndex] = rescaledValue;
                        }
                    }
                    break;
                case gdcm::PixelFormat::FLOAT64:
                    {
                        double* typeScalars = reinterpret_cast<double*>(sliceScalars);
                        size_t length = dataLength / 8;

                        for (size_t voxelIndex = 0; voxelIndex < length; ++voxelIndex) {
                            double valueToRescale = typeScalars[voxelIndex];
                            double rescaledValue = ((valueToRescale * slope + intercept) - info.getIntercept()) / info.getSlope();
                            typeScalars[voxelIndex] = rescaledValue;
                        }
                    }
                    break;
                default:
                    LERROR("Unexpected datatype while rescaling... no rescaling applied!");
            }

            //copy the rescaled values into the scalar buffer
            std::memcpy(&dataStorage[posScalar * info.getBytesPerVoxel()], sliceScalars, dataLength);
        }
        //delete temporary data
        delete[] sliceScalars;
    }
    else {
        reader.GetImage().GetBuffer(&dataStorage[posScalar * info.getBytesPerVoxel()]);
    }

    // Return number of voxels rendered
    return info.getDx() * info.getDy() * info.getNumberOfFrames();
}


std::vector<VolumeURL> GdcmVolumeReader::listVolumes(const std::string& url) const {

    VolumeURL urlOrigin(url);

    //check if this is a network path
    if (urlOrigin.getProtocol() == "dicom-scp") {

        vector<VolumeURL> networkOrigins;

        //first: check, if a local version is available
        //construct path
        string path = constructLocalPathFromNetworkOrigin(urlOrigin);

        if (!tgt::FileSystem::dirExists(path) || getFileNamesInDir(path).empty()) {
            //local version not available -> move Series
            //get AE Title and Incoming Port Number
            string aet = dynamic_cast<const GdcmModule*>(VoreenApplication::app()->getModule("gdcm"))->getAeTitle();
            int incomingPort = dynamic_cast<const GdcmModule*>(VoreenApplication::app()->getModule("gdcm"))->getIncomingPortNumber();

            //create DicomNetworkConnector
            DicomNetworkConnector net(*dict_, aet, incomingPort);
            net.setProgressBar(getProgressBar());

            tgt::FileSystem::createDirectoryRecursive(path);

            //move the files
            try {
                net.moveSeries(urlOrigin,path);
            }
            catch (VoreenException v) {
                if (getProgressBar())
                    getProgressBar()->hide();
                throw tgt::FileException(v.what(),"");
            }
        }

        //origin contains search parameters -> create origin from local directory
        //LINFO("Found local data...trying to list Volumes or load single volume directly...");
        VolumeURL localOrigin("dicom",path,urlOrigin.getSearchString());
        std::vector<VolumeURL> localOrigins = listVolumesDirectory(localOrigin);
        //iterate over local origins and change to network origins
        for (std::vector<VolumeURL>::iterator localIt = localOrigins.begin(); localIt != localOrigins.end(); localIt++) {
            VolumeURL nOrigin("dicom-scp",urlOrigin.getPath(),localIt->getSearchString());
            nOrigin.addSearchParameter("port",urlOrigin.getSearchParameter("port"));
            nOrigin.addSearchParameter("SCP-AET",urlOrigin.getSearchParameter("SCP-AET"));
            networkOrigins.insert(networkOrigins.end(),nOrigin);
        }

        if (networkOrigins.empty()) {
            if (getProgressBar())
                getProgressBar()->hide();
            throw tgt::FileException("Could not list Volumes in local directory. Remove local directory and try to load from server.", localOrigin.getPath());
        }
        else
            return networkOrigins;

        return networkOrigins;
    }

    // if URL already contains a seriesID, return it directly
    if (urlOrigin.getSearchParameter("SeriesInstanceUID", false) != "") {
        vector<VolumeURL> result;
        result.push_back(urlOrigin);
        if (getProgressBar())
            getProgressBar()->hide();
        return result;
    }

    //if Standard Dictionary is not yet loaded: try to load
    if (!dict_) {
        LINFO("Trying to load Standard Dictionary...");
        loadStandardDict();
        LINFO("Standard Dictionary successfully loaded.");
    }

    if (gdcm::System::FileIsDirectory(urlOrigin.getPath().c_str())) {
        //fetch all Volumes located in the directory
        return listVolumesDirectory(urlOrigin);
    }
    else if (gdcm::System::FileExists(urlOrigin.getPath().c_str())) {
        //select if it is a DICOMDIR or a single DICOM image
        if (isDicomDir(urlOrigin.getPath()))
            return listVolumesDicomDir(urlOrigin);
        else if (isDicomFile(urlOrigin.getPath())) {
            //return listVolumesSingleDicomImage(urlOrigin);
            //work around for some data sets: scan whole directory
            return listVolumesDirectory(tgt::FileSystem::dirName(urlOrigin.getPath()));
        }
        else
            throw tgt::FileAccessException("Not a DICOM file!", urlOrigin.getPath());
    }
    else {
        //file is neither a directory nor a file
        if (getProgressBar())
            getProgressBar()->hide();
        throw tgt::FileException("File(s) could not be read: " + urlOrigin.getPath(), urlOrigin.getPath());
    }
}

std::vector<VolumeURL> GdcmVolumeReader::listVolumesDirectory(const VolumeURL& origin) const {
    vector<VolumeURL> result;

    //get all SeriesInstanceUID values of the DICOM files in the Directory
    vector<string> filenames = getFileNamesInDir(origin.getPath());

    std::set<string> seriesInstanceUIDvalues;

    if (getProgressBar())
        getProgressBar()->setTitle("Reading SeriesInstanceUID from all files in Directory...");

    unsigned int itemused = 1;

    vector<string>::iterator fileIterator;
    for (fileIterator = filenames.begin(); fileIterator != filenames.end(); fileIterator++) {

        if (getProgressBar()) {
            getProgressBar()->setProgressMessage("Reading file " + (*fileIterator));
            getProgressBar()->setProgress(static_cast<float>(itemused)/static_cast<float>(filenames.size()));
            itemused++;
        }
        //try to read the file
        gdcm::Reader reader;
        reader.SetFileName(fileIterator->c_str());
        if (!reader.Read())
            continue;
        else {
            //check for valid File Header
            if (!reader.GetFile().GetHeader().IsValid()) {
                LERROR("File Header not valid: " << *fileIterator);
                continue;
            }

            //get the SeriesInstanceUID of the file
            gdcm::StringFilter sf;
            sf.SetFile(reader.GetFile());

            string seriesInstanceUID = sf.ToString(getTagFromDictEntry(dict_->getDictEntryByKeyword("SeriesInstanceUID")));
            seriesInstanceUID = trim(seriesInstanceUID, " ");

            if (seriesInstanceUID.empty()) {
                LERROR("File has empty SeriesInstanceUID: " << *fileIterator);
                continue;
            }
            else {
                //if the origin already contains a SeriesInstanceUID:
                //just use this one
                if (origin.getSearchParameter("SeriesInstanceUID").empty() || origin.getSearchParameter("SeriesInstanceUID") == seriesInstanceUID)
                    seriesInstanceUIDvalues.insert(seriesInstanceUID);
            }

        }
    }

    if (getProgressBar())
        getProgressBar()->hide();

    if (getProgressBar())
        getProgressBar()->setTitle("Checking possible Volumes...");

    unsigned int count = 1;

    //now for each SeriesInstanceUID: check if to subdivide by using CustomDicomDir
    std::set<string>::iterator seriesIterator;

    for (seriesIterator = seriesInstanceUIDvalues.begin(); seriesIterator != seriesInstanceUIDvalues.end(); ++seriesIterator) {

        if (getProgressBar()) {
            getProgressBar()->setProgressMessage("Checking Volume ...");
            getProgressBar()->setProgress(static_cast<float>(count-1) / static_cast<float>(seriesInstanceUIDvalues.size()));
            count++;
        }

        //check, if there is at least one Image file with this SeriesInstanceUID

        string imageFile = "";

        vector<string>::iterator fileIterator;
        for (fileIterator = filenames.begin(); fileIterator != filenames.end(); ++fileIterator) {
            //if this file is not a DICOM file: next file
            if (!isDicomFile(*fileIterator))
                continue;

            //get the SeriesInstanceUID of the DICOM file and compare it to the present one
            //if it does not match: next file
            string fileSeries = getMetaDataFromFile(*fileIterator, *dict_, "SeriesInstanceUID");
            fileSeries = trim(fileSeries, " ");

            if (!((*seriesIterator) == fileSeries))
                continue;
            else {
                //this image file can be used as a reference for the series now
                imageFile = (*fileIterator);
                break;
            }
        }

        if (!imageFile.empty()) {
            //there is at least one file in this Series that can be used as a reference
            //check if a CustomDicomDict needs to be used
            vector<CustomDicomDict> customDicts = getCustomDicts();

            vector<CustomDicomDict>::iterator dictIt;

            for (dictIt = customDicts.begin(); dictIt != customDicts.end(); dictIt++) {

                //check if the current CustomDicomDict is right for the given file
                const vector<pair<string, vector<string> > >* conditions = dictIt->getConditions();

                bool fits = true; //change if conditions do not apply

                //iterate over all conditions and check if they apply for the given files
                vector<pair<string, vector<string> > >::const_iterator condIterator;
                for (condIterator = conditions->begin(); condIterator != conditions->end(); condIterator++) {
                    string value = getMetaDataFromFile(imageFile, *dictIt->getDict(), condIterator->first);
                    value = trim(value, " ");
                    //iterate over possible values for the given tag
                    fits = false;
                    vector<string>::const_iterator valueIt;
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
                    //subdivide the files with the help of the CustomDicomDict
                    //and construct a VolumeOrgin for each group

                    //first: sort out files in another Series
                    vector<string> seriesFiles = getFilesInSeries(filenames, (*seriesIterator));
                    vector<vector<string> > groups = subdivideSeriesFilesByCustomDict(seriesFiles, *dictIt);

                    vector<vector<string> >::iterator groupIterator;

                    for (groupIterator = groups.begin(); groupIterator != groups.end(); ++groupIterator) {
                        if (!(groupIterator->empty())) {
                            VolumeURL o("dicom", origin.getPath());
                            o.addSearchParameter("SeriesInstanceUID", *seriesIterator);
                            o.addSearchParameter("PatientID", getMetaDataFromFile(*groupIterator->begin(), *dict_, "PatientID"));
                            o.addSearchParameter("StudyInstanceUID", getMetaDataFromFile(*groupIterator->begin(), *dict_, "StudyInstanceUID"));
                            //add SearchParameter values to VolumeURL for subdivisionKeywords
                            vector<string> subdivisionKeywords = *(dictIt->getSubdivisionKeywords());
                            const DicomDict customDict = *(dictIt->getDict());

                            vector<string>::iterator keyIterator;
                            for (keyIterator = subdivisionKeywords.begin(); keyIterator != subdivisionKeywords.end(); ++keyIterator) {
                                string value = getMetaDataFromFile(*groupIterator->begin(), customDict, *keyIterator);
                                value = trim(value, " ");

                                o.addSearchParameter(*keyIterator, value);
                            }

                            //now add Meta Information of StandardDictionary
                            setMetaDataFromDict(&(o.getMetaDataContainer()), dict_, *groupIterator->begin());
                            //now add Meta Information of CustomDicomDict
                            setMetaDataFromDict(&(o.getMetaDataContainer()), &customDict, *groupIterator->begin());

                            result.push_back(o);
                        }
                    }

                    break; //do not check the remaining CustomDicomDicts
                }

            }

            //if no CustomDicom could be found:
            if (dictIt == customDicts.end()) {
                //for now only the first file found is used
                VolumeURL v("dicom", origin.getPath());

                v.addSearchParameter("SeriesInstanceUID", *seriesIterator);
                v.addSearchParameter("PatientID", getMetaDataFromFile(imageFile, *dict_, "PatientID"));
                v.addSearchParameter("StudyInstanceUID", getMetaDataFromFile(imageFile, *dict_, "StudyInstanceUID"));

                //now add Meta Information of StandardDictionary
                setMetaDataFromDict(&(v.getMetaDataContainer()), dict_, imageFile);
                result.push_back(v);
            }
        }
    }

    if (getProgressBar())
        getProgressBar()->hide();

    return result;
}

std::vector<VolumeURL> GdcmVolumeReader::listVolumesDicomDir(const VolumeURL& origin) const {
    vector<VolumeURL> result;

    DicomDirParser parser(origin.getPath());
    parser.open();

    //parse the DICOMDIR to find every SeriesInstanceUID and all Image files
    vector<string> filenames; //contains all image files
    vector<string> seriesInstanceUIDs; //contains all SeriesInstanceUID values

    if (getProgressBar())
        getProgressBar()->setTitle("Parsing DICOMDIR Structure to find Volumes...");

    for (parser.setIteratorToFirstDirectoryRecordSequence(); !(parser.iteratorIsAtEnd()); parser.setIteratorToNextDirectoryRecordSequence()) {

        DirectoryRecordSequenceIterator sequence = parser.getPresentDirectoryRecordSequence();
        //iterate over the items
        unsigned int itemused = 0;

        for (sequence.setIteratorToFirstElement(); !(sequence.iteratorIsAtEnd()); sequence.setIteratorToNextElement()) {

            if (getProgressBar()) {
                getProgressBar()->setProgressMessage("Reading DirectoryRecordSequence ...");
                getProgressBar()->setProgress(static_cast<float>(itemused) / static_cast<float>(sequence.getNumberOfItems()));
                itemused++;
            }

            string itemName = sequence.getPresentItemName();
            itemName = trim(itemName, " ");

            if (itemName=="IMAGE"){
                //get path to the image
                string filename;

                if (sequence.findTagInPresentItem(getTagFromDictEntry(dict_->getDictEntryByKeyword("ReferencedFileID"))))
                   filename = sequence.getElementValue(getTagFromDictEntry(dict_->getDictEntryByKeyword("ReferencedFileID")));
                filename = trim(filename, " ");

#ifndef WIN32
                size_t pos;
                while ((pos = filename.find_first_of('\\')) != string::npos) {
                    filename.replace(pos, 1, "/");
                }
#endif
                //get path and try to add the file
                string dir(string(origin.getPath()).substr(0, string(origin.getPath()).length() - string("DICOMDIR").length()));
                dir.append(filename);
                //remove space at the end of string
                dir = trim(dir, " ");
                //add file if it is found
                if (tgt::FileSystem::fileExists(dir)){
                    filenames.push_back(dir);
                }
            }
            else if (itemName=="SERIES"){
                //get SeriesInstanceUID
                string seriesInstanceUID;
                if (sequence.findTagInPresentItem(getTagFromDictEntry(dict_->getDictEntryByKeyword("SeriesInstanceUID"))))
                   seriesInstanceUID = sequence.getElementValue(getTagFromDictEntry(dict_->getDictEntryByKeyword("SeriesInstanceUID")));
                seriesInstanceUID = trim(seriesInstanceUID, " ");

                seriesInstanceUIDs.push_back(seriesInstanceUID);
            }
        } //end (iteration over items)
    } //end of iteration over DES

    if (getProgressBar())
        getProgressBar()->hide();

    if (getProgressBar())
        getProgressBar()->setTitle("Checking possible Volumes...");

    unsigned int count = 1;

    //now for each SeriesInstanceUID: a) check if there are Image files, b) eventually subdivide by using CustomDicomDir
    vector<string>::iterator seriesIterator;

    for (seriesIterator = seriesInstanceUIDs.begin(); seriesIterator != seriesInstanceUIDs.end(); ++seriesIterator) {

        if (getProgressBar()) {
            getProgressBar()->setProgressMessage("Checking Volume ...");
            getProgressBar()->setProgress(static_cast<float>(count-1) / static_cast<float>(seriesInstanceUIDs.size()));
            count++;
        }

        //check, if there is at least one Image file with this SeriesInstanceUID

        string imageFile = "";

        vector<string>::iterator fileIterator;
        for (fileIterator = filenames.begin(); fileIterator != filenames.end(); ++fileIterator) {
            //if this file is not a DICOM file: next file
            if (!isDicomFile(*fileIterator))
                continue;

            //get the SeriesInstanceUID of the DICOM file and compare it to the present one
            //if it does not match: next file
            string fileSeries = getMetaDataFromFile(*fileIterator, *dict_, "SeriesInstanceUID");
            fileSeries = trim(fileSeries, " ");

            if (!((*seriesIterator) == fileSeries))
                continue;
            else {
                //this image file can be used as a reference for the series now
                imageFile = (*fileIterator);
                break;
            }
        }

        if (!imageFile.empty()) {
            //there is at least one file in this Series that can be used as a reference
            //check if a CustomDicomDict needs to be used
            vector<CustomDicomDict> customDicts = getCustomDicts();

            if (getProgressBar())
                getProgressBar()->setTitle("Looking for CustomDicomDict-Files.");

            vector<CustomDicomDict>::iterator dictIt;
            int itemused = 1;

            for (dictIt = customDicts.begin(); dictIt != customDicts.end(); dictIt++) {

                if (getProgressBar()) {
                    getProgressBar()->setProgressMessage("Checking CustomDicomDict-Files...");
                    getProgressBar()->setProgress(static_cast<float>(itemused) / static_cast<float> (customDicts.size()));
                    itemused++;
                }

                //check if the current CustomDicomDict is right for the given file
                const vector<pair<string, vector<string> > >* conditions = dictIt->getConditions();

                bool fits = true; //change if conditions do not apply

                //iterate over all conditions and check if they apply for the given files
                vector<pair<string, vector<string> > >::const_iterator condIterator;
                for (condIterator = conditions->begin(); condIterator != conditions->end(); condIterator++) {
                    string value = getMetaDataFromFile(imageFile, *dictIt->getDict(), condIterator->first);
                    value = trim(value, " ");
                    //iterate over possible values for the given tag
                    fits = false;
                    vector<string>::const_iterator valueIt;
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
                    //subdivide the files with the help of the CustomDicomDict
                    //and construct a VolumeOrgin for each group

                    //first: sort out files in another Series
                    vector<string> seriesFiles = getFilesInSeries(filenames, (*seriesIterator));

                    vector<vector<string> > groups = subdivideSeriesFilesByCustomDict(seriesFiles, *dictIt);

                    vector<vector<string> >::iterator groupIterator;

                    for (groupIterator = groups.begin(); groupIterator != groups.end(); ++groupIterator) {
                        if (!(groupIterator->empty())) {
                            VolumeURL o("dicom", origin.getPath());
                            o.addSearchParameter("SeriesInstanceUID", *seriesIterator);

                            //add SearchParameter values to VolumeURL for subdivisionKeywords
                            vector<string> subdivisionKeywords = *(dictIt->getSubdivisionKeywords());
                            const DicomDict customDict = *(dictIt->getDict());

                            vector<string>::iterator keyIterator;
                            for (keyIterator = subdivisionKeywords.begin(); keyIterator != subdivisionKeywords.end(); ++keyIterator) {
                                string value = getMetaDataFromFile(*groupIterator->begin(), customDict, *keyIterator);
                                value = trim(value, " ");

                                o.addSearchParameter(*keyIterator, value);
                            }

                            //now add Meta Information of StandardDictionary
                            setMetaDataFromDict(&(o.getMetaDataContainer()), dict_, *groupIterator->begin());
                            //now add Meta Information of CustomDicomDict
                            setMetaDataFromDict(&(o.getMetaDataContainer()), &customDict, *groupIterator->begin());

                            result.push_back(o);
                        }
                    }

                    break; //do not check the remaining CustomDicomDicts
                }
            }

            //if no CustomDicom could be found:
            if (dictIt == customDicts.end()) {
                //for now only the first file found is used
                VolumeURL v("dicom", origin.getPath());
                v.addSearchParameter("SeriesInstanceUID", *seriesIterator);
                //now add Meta Information of StandardDictionary
                setMetaDataFromDict(&(v.getMetaDataContainer()), dict_, imageFile);
                result.push_back(v);
            }
        }
    }

    if (getProgressBar())
        getProgressBar()->hide();

    return result;
}

std::vector<VolumeURL> GdcmVolumeReader::listVolumesSingleDicomImage(const VolumeURL& origin) const {
    //create new VolumeURL to be able to add protocol string and additional information
    VolumeURL newOrigin("dicom",origin.getPath());
    //get SeriesInstanceUID and add it to the VolumeURL
    string s = getMetaDataFromFile(newOrigin.getPath(), *dict_, "SeriesInstanceUID");
    s = trim(s, " ");
    newOrigin.addSearchParameter("SeriesInstanceUID", s);

    //now add Meta Information of StandardDictionary
    setMetaDataFromDict(&(newOrigin.getMetaDataContainer()), dict_, newOrigin.getPath());

    //Check, if there are CustomDicomDicts suitable for the given file
    vector<CustomDicomDict> customDicts = getCustomDicts();

    vector<CustomDicomDict>::iterator dictIt;

    for (dictIt = customDicts.begin(); dictIt != customDicts.end(); dictIt++) {

        //check if the current CustomDicomDict is right for the given file
        const vector<pair<string, vector<string> > >* conditions = dictIt->getConditions();

        bool fits = true; //change if conditions do not apply

        //iterate over all conditions and check if they apply for the given file
        vector<pair<string, vector<string> > >::const_iterator condIterator;
        for (condIterator = conditions->begin(); condIterator != conditions->end(); condIterator++) {
            string value = getMetaDataFromFile(newOrigin.getPath(),*dict_, condIterator->first);
            value = trim(value, " ");
            //iterate over possible values for the given tag
            fits = false;
            vector<string>::const_iterator valueIt;
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

            vector<string> subdivisionKeywords = *(dictIt->getSubdivisionKeywords());
            const DicomDict *customDict = dictIt->getDict();

            //add SearchParameter values to VolumeURL for subdivisionKeywords
            vector<string>::iterator keyIterator;
            for (keyIterator = subdivisionKeywords.begin(); keyIterator != subdivisionKeywords.end(); ++keyIterator) {
                string value = getMetaDataFromFile(newOrigin.getPath(), *dict_, *keyIterator);
                value = trim(value, " ");

                newOrigin.addSearchParameter(*keyIterator, value);
            }

            //now add Meta Information of CustomDicomDict
            setMetaDataFromDict(&(newOrigin.getMetaDataContainer()), customDict, newOrigin.getPath());

            break;
        }
    }

    //return a vector with this one VolumeURL
    vector<VolumeURL> v;
    v.push_back(newOrigin);
    if (getProgressBar())
        getProgressBar()->hide();
    return v;
}

gdcm::PixelFormat GdcmVolumeReader::baseTypeStringToGdcm(const std::string& type) const {
    if (type == "uint8")
        return gdcm::PixelFormat::UINT8;
    else if (type == "int8")
        return gdcm::PixelFormat::INT8;
    else if (type == "uint16")
        return gdcm::PixelFormat::UINT16;
    else if (type == "int16")
        return gdcm::PixelFormat::INT16;
    else if (type == "uint32")
        return gdcm::PixelFormat::UINT32;
    else if (type == "int32")
        return gdcm::PixelFormat::INT32;
    else if (type == "float")
        return gdcm::PixelFormat::FLOAT32;
    else if (type == "double")
        return gdcm::PixelFormat::FLOAT64;
    else
        return gdcm::PixelFormat::UNKNOWN;
}

std::vector<PatientInfo> GdcmVolumeReader::findNetworkPatients(const std::string& remote, const std::string& call, uint16_t portno) const {
    //if Standard Dictionary is not yet loaded: try to load
    if (!dict_) {
        LINFO("Trying to load Standard Dictionary...");
        loadStandardDict();
        LINFO("Standard Dictionary successfully loaded.");
    }

    //get AE Title and Incoming Port Number
    string aet = dynamic_cast<const GdcmModule*>(VoreenApplication::app()->getModule("gdcm"))->getAeTitle();
    int incomingPort = dynamic_cast<const GdcmModule*>(VoreenApplication::app()->getModule("gdcm"))->getIncomingPortNumber();

    //create DicomNetworkConnector
    DicomNetworkConnector net(*dict_, aet, incomingPort);
    net.setProgressBar(getProgressBar());

    try {
        return net.findNetworkPatients(remote,call,portno);
    }
    catch (VoreenException v) {
        throw tgt::FileException(v.what(),remote);
    }
}

std::vector<StudyInfo> GdcmVolumeReader::findNetworkStudies(const std::string& remote, const std::string& call, uint16_t portno, const std::string& patientID) const {
    //if Standard Dictionary is not yet loaded: try to load
    if (!dict_) {
        LINFO("Trying to load Standard Dictionary...");
        loadStandardDict();
        LINFO("Standard Dictionary successfully loaded.");
    }

    //get AE Title and Incoming Port Number
    string aet = dynamic_cast<const GdcmModule*>(VoreenApplication::app()->getModule("gdcm"))->getAeTitle();
    int incomingPort = dynamic_cast<const GdcmModule*>(VoreenApplication::app()->getModule("gdcm"))->getIncomingPortNumber();

    //create DicomNetworkConnector
    DicomNetworkConnector net(*dict_, aet, incomingPort);
    net.setProgressBar(getProgressBar());

    try {
        return net.findNetworkStudies(remote,call,portno,patientID);
    }
    catch (VoreenException v) {
        throw tgt::FileException(v.what(),remote);
    }
}

std::vector<SeriesInfo> GdcmVolumeReader::findNetworkSeries(const std::string& remote, const std::string& call, uint16_t portno, const std::string& patientID, const std::string& studyID) const {
    //if Standard Dictionary is not yet loaded: try to load
    if (!dict_) {
        LINFO("Trying to load Standard Dictionary...");
        loadStandardDict();
        LINFO("Standard Dictionary successfully loaded.");
    }

    //get AE Title and Incoming Port Number
    string aet = dynamic_cast<const GdcmModule*>(VoreenApplication::app()->getModule("gdcm"))->getAeTitle();
    int incomingPort = dynamic_cast<const GdcmModule*>(VoreenApplication::app()->getModule("gdcm"))->getIncomingPortNumber();

    //create DicomNetworkConnector
    DicomNetworkConnector net(*dict_, aet, incomingPort);
    net.setProgressBar(getProgressBar());

    try {
        return net.findNetworkSeries(remote,call,portno, patientID, studyID);
    }
    catch (VoreenException v) {
        throw tgt::FileException(v.what(),remote);
    }
}

} //namespace voreen
