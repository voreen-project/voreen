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

#include <limits>

#include "ext/tinyxml/tinyxml.h"

#include "tgt/init.h"
#include "tgt/vector.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/geometry/meshlistgeometry.h"
#include "voreen/core/io/datvolumewriter.h"
#include "voreen/core/io/datvolumereader.h"
#include "voreen/core/io/textfilereader.h"

//#define MAXDIM 512
//#define MAXDIM 256
//#define MAXDIM 128
#define MAXDIM 64
//#define MAXDIM 40
//#define MAXDIM 32

using namespace voreen;

struct ScalarData {
    float min_;
    float max_;
    std::vector<float> values_;
};

struct VecData {
    tgt::vec3 min_;
    tgt::vec3 max_;
    float minMag_;
    float maxMag_;
    std::vector<tgt::vec3> values_;
};

struct ScalarField {
    VecData coords_;
    ScalarData values_;
};

struct VecField {
    VecData coords_;
    VecData values_;
};

VecData loadVecValues(TiXmlNode* pParent, const std::string& tag) {

    VecData vd;

    if (!pParent)
        return vd;

    std::vector<tgt::vec3> values;
    tgt::vec3 max = -tgt::vec3(std::numeric_limits<float>::max());
    tgt::vec3 min = tgt::vec3(std::numeric_limits<float>::max());
    float maxMag = 0.f;
    float minMag = std::numeric_limits<float>::max();

    if(TiXmlNode* vtk = pParent->FirstChild("VTKFile")) {
        if(TiXmlNode* ugrid = vtk->FirstChild("UnstructuredGrid")) {
            if(TiXmlNode* piece = ugrid->FirstChild("Piece")) {
                if(TiXmlNode* ptData = piece->FirstChild(tag)) {
                    if(TiXmlNode* array = ptData->FirstChild("DataArray")) {
                        std::string strValues = array->ToElement()->GetText();
                        std::stringstream s(strValues);

                        bool done = false;
                        while(!done) {
                            tgt::vec3 out;
                            for(int i = 0; i < 3; i++) {
                                if(!(s >> out[i])) {
                                    done = true;
                                    break;
                                }
                            }
                            if(!done) {
                                if(out.x > max.x) max.x = out.x;
                                if(out.y > max.y) max.y = out.y;
                                if(out.z > max.z) max.z = out.z;
                                if(out.x < min.x) min.x = out.x;
                                if(out.y < min.y) min.y = out.y;
                                if(out.z < min.z) min.z = out.z;
                                if(length(out) > maxMag) maxMag = length(out);
                                if(length(out) < minMag) minMag = length(out);

                                values.push_back(out);
                            }
                        }
                    }
                }
            }
        }
    }


    vd.values_ = values;
    vd.min_ = min;
    vd.max_ = max;
    vd.minMag_ = minMag;
    vd.maxMag_ = maxMag;

    return vd;
}

ScalarData loadScalarValues(TiXmlNode* pParent, const std::string& tag) {

    ScalarData sd;

    if (!pParent)
        return sd;

    std::vector<float> values;
    float max = -std::numeric_limits<float>::max();
    float min = std::numeric_limits<float>::max();

    if(TiXmlNode* vtk = pParent->FirstChild("VTKFile")) {
        if(TiXmlNode* ugrid = vtk->FirstChild("UnstructuredGrid")) {
            if(TiXmlNode* piece = ugrid->FirstChild("Piece")) {
                if(TiXmlNode* ptData = piece->FirstChild(tag)) {
                    if(TiXmlNode* array = ptData->FirstChild("DataArray")) {
                        std::string strValues = array->ToElement()->GetText();
                        std::stringstream s(strValues);

                        bool done = false;
                        while(!done) {
                            float out;
                            if(!(s >> out))
                                done = true;

                            if(!done) {
                                if(out > max) max = out;
                                if(out < min) min = out;

                                values.push_back(out);
                            }
                        }
                    }
                }
            }
        }
    }

    sd.values_ = values;
    sd.min_ = min;
    sd.max_ = max;
    std::cout << "Loaded scalar field with " << values.size() << " values, minimum value: " << min << ", maximum value: " << max << std::endl;

    return sd;
}

void writeScalarFieldToVolume(ScalarField sf, const std::string& outName) {
    tgt::vec3 span = sf.coords_.max_ - sf.coords_.min_;
    tgt::ivec3 dimensions = tgt::ivec3(tgt::vec3(MAXDIM) * (span / tgt::max(span)));

    VolumeRAM_Float* target = new VolumeRAM_Float(dimensions);
    target->clear();

    for(size_t i = 0; i < sf.coords_.values_.size(); i++) {
        tgt::ivec3 loc = tgt::ivec3(((sf.coords_.values_.at(i) - sf.coords_.min_) / span) * (tgt::vec3(dimensions) - tgt::vec3(1)));

        float val = sf.values_.values_.at(i);
        target->voxel(loc) = val;
    }

    Volume* vh = new Volume(target, tgt::vec3(1.f), tgt::vec3(0.0f));
    oldVolumePosition(vh);
    DatVolumeWriter dvw;
    std::cout << "Saving dataset with dimensions " << dimensions << std::endl;
    dvw.write(outName, vh);

    delete vh;
}

void writeVecFieldToVolume(VecField vf, const std::string& outName) {
    tgt::vec3 span = vf.coords_.max_ - vf.coords_.min_;
    tgt::ivec3 dimensions = tgt::ivec3(tgt::vec3(MAXDIM) * (span / tgt::max(span)));

    VolumeRAM_3xFloat* target = new VolumeRAM_3xFloat(dimensions);
    target->clear();
    VolumeRAM_UInt8* counter = new VolumeRAM_UInt8(dimensions);
    counter->clear();

    for(size_t i = 0; i < vf.coords_.values_.size(); i++) {
        tgt::vec3 locN = (vf.coords_.values_.at(i) - vf.coords_.min_) / span;
        tgt::vec3 val = 0.1f * vf.values_.values_.at(i);

        //tgt::ivec3 loc = tgt::ivec3(locN * (tgt::vec3(dimensions) - tgt::vec3(1)));
        //target->voxel(loc) = val;

        tgt::vec3 step = tgt::vec3(1.f) / tgt::vec3(dimensions);
        int amountSteps = length(val) / tgt::min(step);
        tgt::vec3 sample = locN;
        tgt::vec3 dir = normalize(val);

        for(int j = 0; j < amountSteps; j++) {
            if(sample.x < 0.f || sample.x > 1.f || sample.y < 0.f || sample.y > 1.f || sample.z < 0.f || sample.z > 1.f)
                break;

            tgt::ivec3 loc = tgt::ivec3(sample * (tgt::vec3(dimensions) - tgt::vec3(1)));

            ////target->voxel(loc) = val;
            //target->voxel(loc) += val;
            ////target->voxel(loc) = val / vf.values_.maxMag_;
            //target->voxel(loc) += val / vf.values_.maxMag_;

            tgt::ivec3 tmpLoc;
            for(int k = -1; k <= 1; k++) {
                for(int l = -1; l <= 1; l++) {
                    for(int m = -1; m <= 1; m++) {
                        tgt::ivec3 offset = tgt::ivec3(k, l, m);
                        tmpLoc = loc + offset;
                        tmpLoc = clamp(tmpLoc, tgt::ivec3(0), dimensions - tgt::ivec3(1));
                        target->voxel(tmpLoc) += val;
                        counter->voxel(tmpLoc) += 1;
                    }
                }
            }

            sample += tgt::min(step) * dir;
        }
    }

    for(size_t i = 0; i < target->getNumVoxels(); i++) {
        if(counter->voxel(i) > 0) {
            target->voxel(i) /= float(counter->voxel(i));
            //target->voxel(i) = (target->voxel(i) - sf.values_.min_) / (sf.values_.max_ - sf.values_.min_);
            //target->voxel(i) = (target->voxel(i) - vf.values_.minMag_) / (vf.values_.maxMag_ - vf.values_.minMag_);
        }
    }

    Volume* vh = new Volume(target, tgt::vec3(1.f), tgt::vec3(0.0f));
    oldVolumePosition(vh);
    DatVolumeWriter dvw;
    std::cout << "Saving dataset with dimensions " << dimensions << std::endl;
    dvw.write(outName, vh);

    delete vh;
    delete counter;
}

void writeVecAndScalarFieldToVolume(VecField vf, ScalarField sf, const std::string& outName) {
    tgt::vec3 span = vf.coords_.max_ - vf.coords_.min_;
    tgt::ivec3 dimensions = tgt::ivec3(tgt::vec3(MAXDIM) * (span / tgt::max(span)));

    //VolumeRAM_4xFloat* target = new VolumeRAM_4xFloat(dimensions);
    VolumeRAM_Float* target = new VolumeRAM_Float(dimensions);
    VolumeRAM_UInt8* counter = new VolumeRAM_UInt8(dimensions);
    target->clear();
    counter->clear();

    for(size_t i = 0; i < vf.coords_.values_.size(); i++) {
        tgt::vec3 locN = (vf.coords_.values_.at(i) - vf.coords_.min_) / span;
        //tgt::vec3 val = 0.005f * vf.values_.values_.at(i);
        //tgt::vec3 val = 0.0025f * vf.values_.values_.at(i);
        tgt::vec3 val = 0.0015f * vf.values_.values_.at(i);
        //tgt::vec3 val = 0.0005f * vf.values_.values_.at(i);

        tgt::vec3 step = tgt::vec3(1.f) / tgt::vec3(dimensions);
        int amountSteps = length(val) / tgt::min(step);
        tgt::vec3 sample = locN;
        tgt::vec3 dir = normalize(val);

        for(int j = 0; j < amountSteps; j++) {
            if(sample.x < 0.f || sample.x > 1.f || sample.y < 0.f || sample.y > 1.f || sample.z < 0.f || sample.z > 1.f)
                break;
            tgt::ivec3 loc = tgt::ivec3(sample * (tgt::vec3(dimensions) - tgt::vec3(1)));
            //target->voxel(loc) = val;
            //target->voxel(loc).xyz() = val / vf.values_.maxMag_;
            //target->voxel(loc).w = (sf.values_.values_.at(i) - sf.values_.min_) / (sf.values_.max_ - sf.values_.min_);

            //target->voxel(loc) = (sf.values_.values_.at(i) - sf.values_.min_) / (sf.values_.max_ - sf.values_.min_);
            //target->voxel(loc) += (sf.values_.values_.at(i) - sf.values_.min_) / (sf.values_.max_ - sf.values_.min_);
            //target->voxel(loc) += 10.f * length(val);

            //target->voxel(loc) += sf.values_.values_.at(i);
            //counter->voxel(loc) += 1;
            //
            tgt::ivec3 tmpLoc;
            for(int k = -1; k <= 1; k++) {
                for(int l = -1; l <= 1; l++) {
                    for(int m = -1; m <= 1; m++) {
                        tgt::ivec3 offset = tgt::ivec3(k, l, m);
                        tmpLoc = loc + offset;
                        std::cout << "-------" << std::endl;
                        std::cout << tmpLoc << std::endl;
                        tmpLoc = clamp(tmpLoc, tgt::ivec3(0), dimensions - tgt::ivec3(1));
                        std::cout << tmpLoc << std::endl;
                        target->voxel(tmpLoc) += sf.values_.values_.at(i);
                        counter->voxel(tmpLoc) += 1;
                    }
                }
            }

            sample += tgt::min(step) * dir;
        }
    }

    for(size_t i = 0; i < target->getNumVoxels(); i++) {
        if(counter->voxel(i) > 0) {
            target->voxel(i) /= float(counter->voxel(i));
            target->voxel(i) = (target->voxel(i) - sf.values_.min_) / (sf.values_.max_ - sf.values_.min_);
            //target->voxel(i) = (target->voxel(i) - vf.values_.minMag_) / (vf.values_.maxMag_ - vf.values_.minMag_);
        }
    }

    Volume* vh = new Volume(target, tgt::vec3(1.f), tgt::vec3(0.0f));
    oldVolumePosition(vh);
    DatVolumeWriter dvw;
    std::cout << "Saving dataset with dimensions " << dimensions << std::endl;
    dvw.write(outName, vh);

    delete vh;
    delete counter;
}

std::vector<std::string> getSteps(const std::string& xmlName) {
    std::vector<std::string> fileNames;
    TiXmlDocument doc(xmlName);
    bool loadOkay = doc.LoadFile();
    if (!loadOkay)
        return fileNames;

    if(TiXmlNode* vtk = doc.FirstChild("VTKFile")) {
        if(TiXmlNode* coll = vtk->FirstChild("Collection")) {
            for(TiXmlNode* child = coll->FirstChild(); child; child = child->NextSibling()) {

                TiXmlElement* el = child->ToElement();
                const char* attr = el->Attribute("file");
                if(!attr)
                    return fileNames;

                std::string fn = tgt::FileSystem::dirName(xmlName);
                if(!fn.empty())
                    fn += "/";
                fileNames.push_back(fn + std::string(attr));
            }
        }
    }

    return fileNames;
}

ScalarField convertScalarData(const std::string& xmlName) {
    ScalarField sf;

    TiXmlDocument doc(xmlName);
    if(!doc.LoadFile()) {
        std::cerr << "Failed to load file " << xmlName << std::endl;
        return sf;
    }

    sf.coords_ = loadVecValues(&doc, "Points");
    sf.values_ = loadScalarValues(&doc, "PointData");

    return sf;
}

VecField convertVecData(const std::string& xmlName) {
    VecField vf;

    TiXmlDocument doc(xmlName);
    if(!doc.LoadFile()) {
        std::cerr << "Failed to load file " << xmlName << std::endl;
        return vf;
    }

    vf.coords_ = loadVecValues(&doc, "Points");
    vf.values_ = loadVecValues(&doc, "PointData");

    return vf;
}

void writeSeries(const std::vector<std::string>& fileNames) {

    if(fileNames.empty())
        return;

    DatVolumeReader dvr;
    VolumeBase* vhb = dvr.read(fileNames.at(0))->first();
    if(!vhb)
        return;

    std::stringstream datContent;
    tgt::ivec3 dims = vhb->getDimensions();
    tgt::vec3 spac = vhb->getSpacing();
    if(length(spac) == 0.f)
        spac = tgt::vec3(1.f);
    datContent << "Resolution:\t" << dims.x << " " << dims.y << " " << dims.z << std::endl;
    datContent << "SliceThickness:\t" << spac.x << " " << spac.y << " " << spac.z << std::endl;

    datContent << "Format:\t" << "FLOAT" << std::endl;

    if(vhb->getNumChannels() == 1)
        datContent << "ObjectModel:\t" << "I" << std::endl;
    else if(vhb->getNumChannels() == 3)
        datContent << "ObjectModel:\t" << "RGB" << std::endl;
    else if(vhb->getNumChannels() == 4)
        datContent << "ObjectModel:\t" << "RGBA" << std::endl;

    for(size_t i = 0; i < fileNames.size(); i++)
        datContent << "ObjectFileName:\t" << tgt::FileSystem::baseName(fileNames.at(i)).append(".raw") << std::endl;

    std::string sdatname = tgt::FileSystem::fullBaseName(fileNames.at(0)).append(".sdat");
    std::fstream sdatout(sdatname.c_str(), std::ios::out);

    if (!sdatout.is_open() || sdatout.bad())
        return;

    sdatout << datContent.str();
    if (sdatout.bad())
        std::cerr << "Error while writing volume series file." << std::endl;
    sdatout.close();
}

std::vector<tgt::vec3> convertGeometryDataFloat(const std::string& fileName) {
    TextFileReader tfr(fileName);

    std::string line;

    std::vector<tgt::vec3> result;

    tgt::vec3 max = -tgt::vec3(std::numeric_limits<float>::max());
    tgt::vec3 min = tgt::vec3(std::numeric_limits<float>::max());

    bool done = false;
    while(tfr.getNextLinePlain(line, false)) {
        std::stringstream s(line);
        tgt::vec3 out;
        for(int i = 0; i < 3; i++) {
            if(!(s >> out[i])) {
                done = true;
                break;
            }
        }
        if(done)
            break;

        result.push_back(out);
        if(out.x > max.x) max.x = out.x;
        if(out.y > max.y) max.y = out.y;
        if(out.z > max.z) max.z = out.z;
        if(out.x < min.x) min.x = out.x;
        if(out.y < min.y) min.y = out.y;
        if(out.z < min.z) min.z = out.z;
    }

    tgt::vec3 span = max - min;
    //tgt::vec3 cs = 2.f * span / tgt::max(span);
    tgt::vec3 cs = span / tgt::max(span);

    for(size_t i = 0; i < result.size(); i++) {
        //TODO doesn't yet work as intended
        //result.at(i) = cs * ((result.at(i) - min) / span) - 1.f;
        result.at(i) = 2.f * ((result.at(i) - min) / span) - 1.f;
        //result.at(i).y -= 0.25f;
        result.at(i) *= cs;
    }

    return result;
}

std::vector<tgt::ivec3> convertGeometryDataInt(const std::string& fileName) {
    TextFileReader tfr(fileName);

    std::string line;

    std::vector<tgt::ivec3> result;

    bool done = false;
    while(tfr.getNextLinePlain(line, false)) {
        tgt::ivec3 out;
        std::stringstream s(line);

        for(int i = 0; i < 3; i++) {
            if(!(s >> out[i])) {
                done = true;
                break;
            }
        }
        result.push_back(out);
        if(done)
            break;
    }

    return result;
}

void writeGeometry(const std::vector<tgt::vec3>& vf, const std::vector<tgt::ivec3>& indices, const std::string& outName) {

    MeshGeometry mesh;
    MeshListGeometry meshList;
    for(size_t i = 0; i < indices.size(); i++) {
        FaceGeometry face;
        for(int j = 0; j < 3; j++) {
            VertexGeometry vert;
            int index = indices.at(i)[j];
            vert.setCoords(vf[index]);
            vert.setColor(tgt::vec4(1.f, 1.f, 0.f, 1.0f));
            face.addVertex(vert);
        }

        tgt::vec3 v1 = face.getVertex(1).getCoords() - face.getVertex(0).getCoords();
        tgt::vec3 v2 = face.getVertex(2).getCoords() - face.getVertex(0).getCoords();

        tgt::vec3 normal = normalize(cross(v1, v2));
        face.setFaceNormal(normal);

        mesh.addFace(face);
    }

    meshList.addMesh(mesh);

    XmlSerializer s(outName);

    std::fstream f;
    f.open(outName.c_str(), std::ios::out);
    meshList.serialize(s);
    s.write(f);
    f.close();
}

int main(int argc, char* argv[]) {
    tgt::Singleton<tgt::LogManager>::init();

    if(argc < 2)
        return -1;

    if(argc == 3) {
        std::vector<std::string> fileNames = getSteps(argv[1]);
        std::vector<std::string> datNames;
        std::vector<VecField> vecFields;

        float globMaxMag = 0.f;
        float globMinMag = std::numeric_limits<float>::max();

        if(tgt::FileSystem::fileExtension(argv[1]) == "pvd" && tgt::FileSystem::fileExtension(argv[2]) == "pvd") {
            std::vector<std::string> fileNamesScalar = getSteps(argv[2]);
            std::vector<ScalarField> scalarFields;

            float globMaxScalar = -std::numeric_limits<float>::max();
            float globMinScalar = std::numeric_limits<float>::max();

            for(size_t i = 0; i < fileNames.size(); i++) {

                VecField vf = convertVecData(fileNames.at(i));
                ScalarField sf = convertScalarData(fileNamesScalar.at(i));
                if(vf.values_.maxMag_ > globMaxMag) globMaxMag = vf.values_.maxMag_;
                if(vf.values_.minMag_ < globMinMag) globMinMag = vf.values_.minMag_;
                if(sf.values_.max_ > globMaxScalar) globMaxScalar = sf.values_.max_;
                if(sf.values_.min_ < globMinScalar) globMinScalar = sf.values_.min_;

                vecFields.push_back(vf);
                scalarFields.push_back(sf);
                datNames.push_back(tgt::FileSystem::fullBaseName(fileNames.at(i)).append(".dat"));
            }

            std::cout << "Global minimum scalar value: " << globMinScalar << ", maximum value: " << globMaxScalar << std::endl;

            for(size_t i = 0; i < vecFields.size(); i++) {
                //vecFields.at(i).values_.minMag_ = globMinMag;
                //vecFields.at(i).values_.maxMag_ = globMaxMag;
                //scalarFields.at(i).values_.min_ = globMinScalar;
                //scalarFields.at(i).values_.max_ = globMaxScalar;
                writeVecAndScalarFieldToVolume(vecFields.at(i), scalarFields.at(i), datNames.at(i));
            }

            writeSeries(datNames);

        } else if(tgt::FileSystem::fileExtension(argv[1]) == "pvd") {

            for(size_t i = 0; i < fileNames.size(); i++) {
                VecField vf = convertVecData(fileNames.at(i));
                if(vf.values_.maxMag_ > globMaxMag) globMaxMag = vf.values_.maxMag_;
                if(vf.values_.minMag_ < globMinMag) globMinMag = vf.values_.minMag_;

                vecFields.push_back(vf);
                datNames.push_back(tgt::FileSystem::fullBaseName(fileNames.at(i)).append(".dat"));
            }

            for(size_t i = 0; i < vecFields.size(); i++) {
                //vecFields.at(i).values_.minMag_ = globMinMag;
                //vecFields.at(i).values_.maxMag_ = globMaxMag;
                writeVecFieldToVolume(vecFields.at(i), datNames.at(i));
                //writeVecAndScalarFieldToVolume(vecFields.at(i), ScalarField(), datNames.at(i));
            }

            writeSeries(datNames);
        } else if(tgt::FileSystem::fileExtension(argv[1]) == "vtu" && tgt::FileSystem::fileExtension(argv[2]) == "vtu") {
            VecField vf = convertVecData(argv[1]);
            ScalarField sf = convertScalarData(argv[2]);
            std::string outName = tgt::FileSystem::fullBaseName(argv[1]).append(".dat");
            writeVecAndScalarFieldToVolume(vf, sf , outName);

        } else if(tgt::FileSystem::fileExtension(argv[1]) == "txt") {
            std::vector<tgt::vec3> vf = convertGeometryDataFloat(argv[1]);
            std::vector<tgt::ivec3> indices  = convertGeometryDataInt(argv[2]);
            std::string outName = tgt::FileSystem::fullBaseName(argv[1]).append(".vge");

            writeGeometry(vf, indices, outName);
        }

    } else if(tgt::FileSystem::fileExtension(argv[1]) == "vtu") {
        VecField vf = convertVecData(argv[1]);
        std::string outName = tgt::FileSystem::fullBaseName(argv[1]).append(".dat");
        writeVecFieldToVolume(vf, outName);
    }
    else
        return -2;

    return 0;
}
