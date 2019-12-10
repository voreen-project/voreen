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

#include "quadhidacvolumereader.h"

#include <iostream>
#include <sstream>
#include <fstream>
#include <limits>

#include "tgt/exception.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatormirror.h"
#include "voreen/core/utils/stringutils.h"

/*
// excerpt of quadhidac file
Title=Dynamic PET Image Sequence.
;% Using List-mode Expectation Maximisation

Subsets=98
Rho=1.200000
ntf=1
dxyr=0.400000
dt=[900]

nbx=150.000000
nby=150.000000
nbz=276.000000

FrameTime=900

-------------
quadhidac files with title KS PET have a mixed x- and y-dimension
*/

namespace voreen {

const std::string QuadHidacVolumeReader::loggerCat_ = "voreen.io.VolumeReader.QuadHidac";

QuadHidacVolumeReader::QuadHidacVolumeReader(ProgressBar* progress) :
    VolumeReader(progress)
{
    extensions_.push_back("i4d");
    protocols_.push_back("i4d");
}

QuadHidacHeader QuadHidacVolumeReader::readHeader(const std::string& filename) const {
    QuadHidacHeader header;

    tgt::ivec3 croppedDims(0);

    std::fstream file(filename.c_str(), std::ios::in | std::ios::binary);

    if (!file)
        throw tgt::FileNotFoundException("reading; " + filename);

    LINFO(filename);
    char ch = 0;
    std::string typeStr;
    std::string argumentsStr;
    bool type = true;
    bool comment = false;
    int pos=-1;

    while ((ch != 0x0c) && (!file.eof())) {
        file.read(&ch, 1);

        if ((ch != 0x0c) && (ch != 0x0a)) {
            if (ch == '=') {
                type = false;
                LINFO("Type: " << typeStr);
            } else if (ch == ';') {
                comment = true;
                type = false;
            } else {
                if (type)
                    typeStr += ch;
                else
                    argumentsStr += ch;
            }
        } else {
            if (comment) {
                LINFO("Comment: " << argumentsStr);
            } else {
                std::stringstream argStream(argumentsStr);

                if (typeStr == "") {
                } else if (typeStr == "Title") {
                    header.title = argumentsStr;
                    pos = static_cast<int>(header.title.find("KS PET"));
                    LINFO("Value: " << header.title);
                } else if (typeStr == "Subsets") {
                    argStream >> header.subsets;
                    LINFO("Value: " << header.subsets);
                } else if (typeStr == "Rho") {
                    argStream >> header.rho;
                    LINFO("Value: " << header.rho);
                } else if (typeStr == "ntf") {
                    argStream >> header.ntf;
                    LINFO("Value: " << header.ntf);
                } else if (typeStr == "nbg") {
                    argStream >> header.nbg;
                    LINFO("Value: " << header.nbg);
                } else if (typeStr == "nbrg") {
                    argStream >> header.nbrg;
                    LINFO("Value: " << header.nbrg);
                } else if (typeStr == "dxyr") {
                    argStream >> header.dxyr;
                    LINFO("Value: " << header.dxyr);
                } else if (typeStr == "row1") {
                    argStream >> header.transformationMatrix.t00;
                    argStream >> header.transformationMatrix.t01;
                    argStream >> header.transformationMatrix.t02;
                    argStream >> header.transformationMatrix.t03;
                } else if (typeStr == "row2") {
                    argStream >> header.transformationMatrix.t10;
                    argStream >> header.transformationMatrix.t11;
                    argStream >> header.transformationMatrix.t12;
                    argStream >> header.transformationMatrix.t13;
                } else if (typeStr == "row3") {
                    argStream >> header.transformationMatrix.t20;
                    argStream >> header.transformationMatrix.t21;
                    argStream >> header.transformationMatrix.t22;
                    argStream >> header.transformationMatrix.t23;
                } else if (typeStr == "row4") {
                    argStream >> header.transformationMatrix.t30;
                    argStream >> header.transformationMatrix.t31;
                    argStream >> header.transformationMatrix.t32;
                    argStream >> header.transformationMatrix.t33;
                } else if (typeStr == "dt") {
                    // expected format: [float[,float]], e.g. '[30.0]' or '[2.0,2.0,30.0]'
                    std::string dtStr;
                    argStream >> dtStr;
                    LINFO("String: " << dtStr);

                    // split dt string into tokens
                    std::vector<std::string> dtTokens;
                    if (dtStr.size() < 2)
                        goto dt_malformed;
                    else if (dtStr.substr(0, 1) != "[" || dtStr.substr(dtStr.size()-1, 1) != "]")
                        goto dt_malformed;
                    dtTokens = strSplit(dtStr.substr(1, dtStr.size() - 2), ',');
                    if (dtTokens.empty())
                        goto dt_malformed;

                    // convert string tokens to float
                    for (size_t i=0; i<dtTokens.size(); i++)
                        header.dt.push_back(stof(dtTokens.at(i)));
                    goto dt_finished;

dt_malformed:       LWARNING("String malformed! Expected: '[float,...]'");
                    continue;

dt_finished:        LINFO("Parsed (" << header.dt.size() << " values" << "): " << "[" << strJoin(header.dt, ",") << "]");

                } else if (typeStr == "nbx") {
                    // if a KS PET file is read then pos = 0
                    if (pos == -1){
                        argStream >> header.dims.x;
                        LINFO("Value: " << header.dims.x);

                    } else { argStream >> header.dims.y;
                        LINFO("Value: " << header.dims.y);
                        LINFO("nbx and nby are twisted");
                    }

                } else if (typeStr == "nby") {
                    //argStream >> dims.y;
                    if (pos == -1){
                        argStream >> header.dims.y;
                        LINFO("Value: " << header.dims.y);
                    } else { argStream >> header.dims.x;
                        LINFO("Value: " << header.dims.x);
                    }
                } else if (typeStr == "nbz") {
                    argStream >> header.dims.z;
                    LINFO("Value: " << header.dims.z);

                } else if (typeStr == "croppednbx") {
                    // if a KS PET file is read then pos = 0
                    if (pos == -1) {
                        argStream >> croppedDims.x;
                    } else {
                        argStream >> croppedDims.y;
                    }

                } else if (typeStr == "croppednby") {
                    // if a KS PET file is read then pos = 0
                    if (pos == -1) {
                        argStream >> croppedDims.y;
                    } else {
                        argStream >> croppedDims.x;
                    }

                } else if (typeStr == "croppednbz") {
                    argStream >> croppedDims.z;

                } else if (typeStr == "FrameTime") {
                    argStream >> header.frameTime;
                    LINFO("Value: " << header.frameTime);
                } else {
                    LWARNING("Unknown parameter: " << typeStr.substr(0, 20) << " (value: " << argumentsStr << ")");
                }

                if (argStream.fail())
                    throw tgt::CorruptedFileException("Error during header parsing", filename);
            }

            typeStr = "";
            argumentsStr = "";
            type = true;
            comment = false;
        }
    }

    if (header.ntf != static_cast<int>(header.dt.size()))
        LWARNING("Number of time frames (" << header.ntf << ") does not match size of 'dt' vector (" << header.dt.size() << ")");

    // if croppedDims is available, use it, nicer version is in work
    if (croppedDims.z > 0)
        header.dims = croppedDims;

    header.dataOffset = static_cast<size_t>(file.tellg());

    file.close();

    return header;
}

Volume* QuadHidacVolumeReader::loadVolume(const std::string& filename, int i, const QuadHidacHeader& header) const {
    std::fstream file(filename.c_str(), std::ios::in | std::ios::binary);

    if (!file)
        throw tgt::FileNotFoundException("Reading vol " + itos(i) + " from " + filename);

    file.seekg(header.dataOffset + hmul(header.dims)*static_cast<size_t>(4*i), std::ios::beg);

    // read frame data (float)
    Volume* vh = 0;
    float* scalars = new float[hmul(header.dims)];
    try {
        file.read(reinterpret_cast<char*>(scalars), hmul(header.dims)*4);

        vh = new Volume(new VolumeRAM_Float(scalars, header.dims), tgt::vec3(header.dxyr), tgt::vec3(0.0f));
        centerVolume(vh);
        //Volume* t = VolumeOperatorMirrorX::APPLY_OP(vh);
        //delete vh;
        //vh = t;
    }
    catch (...) {
        delete[] scalars;
        file.close();
        LERROR("Failed to read frame data from file!");
        return 0;
    }

    vh->setPhysicalToWorldMatrix(header.transformationMatrix);
    /*
    tgt::mat4 m = tgt::mat4(0.0f, 1.0f, 0.0f, 0.0f,
                            1.0f, 0.0f, 0.0f, 0.0f,
                            0.0f, 0.0f, 1.0f, 0.0f,
                            0.0f, 0.0f, 0.0f, 1.0f);
    vh->setPhysicalToWorldMatrix(m);
    */

    if ((size_t)i < header.dt.size()) {
        float t = 0.0f;
        for(int j=0; j<i; j++) {
            t += header.dt[j];
        }
        vh->setTimestep(t);

        vh->setMetaDataValue<IntMetaData>("ActualFrameDuration", static_cast<int>(header.dt[i]));
        vh->setMetaDataValue<IntMetaData>("FrameTime", static_cast<int>(t));
    }

    VolumeURL frameOrigin("i4d", filename);
    frameOrigin.addSearchParameter("frame", itos(i));
    vh->setOrigin(frameOrigin);
    //vh->setRealWorldMapping(RealWorldMapping(1.0f, 0.0f, "BqCC")); // Not actually BqCC:
    vh->setRealWorldMapping(RealWorldMapping(1.0f, 0.0f, "CpsCC"));
    vh->setMetaDataValue<IntMetaData>("Frame", i);

    // clean up
    file.close();

    return vh;
}

VolumeList* QuadHidacVolumeReader::read(const std::string &url) {
    const VolumeURL origin(url);
    std::string filename = origin.getPath();

    QuadHidacHeader header = readHeader(filename);

    std::fstream file(filename.c_str(), std::ios::in | std::ios::binary);

    if (!file)
        throw tgt::FileNotFoundException("reading; " + url);

    // detect frame range to load
    int maxFrame = header.ntf-1;

    // result collection
    VolumeList* volumeList = new VolumeList();

    if (getProgressBar()) {
        getProgressBar()->setTitle("Quad-HIDAC Reader");
        getProgressBar()->show();
    }

    // load frames
    for (int frame=0; frame <= maxFrame; frame++) {

        LINFO("Reading frame " << frame+1 << "/" << header.ntf << " ...");
        if (getProgressBar()) {
            getProgressBar()->setProgressMessage("Loading frame " + itos(frame+1) + "/" + itos(header.ntf) + " ...");
            getProgressBar()->setProgress(static_cast<float>(frame)/(maxFrame + 1.f));
        }

        Volume* vh = loadVolume(filename, frame, header);
        volumeList->add(vh);
    }

    if (getProgressBar())
        getProgressBar()->hide();

    // clean up
    file.close();

    return volumeList;
}

VolumeBase* QuadHidacVolumeReader::read(const VolumeURL& origin) {
    std::string tmp = origin.getSearchParameter("frame");
    if (tmp == "") {
        LERROR("Origin does not specify frame id");
        return 0;
    }
    else {
        std::string filename = origin.getPath();
        QuadHidacHeader header = readHeader(filename);

        int frame = stoi(tmp);
        if(frame < header.ntf)
            return loadVolume(filename, frame, header);
        else
            return 0;
    }
}

std::vector<VolumeURL> QuadHidacVolumeReader::listVolumes(const std::string& url) const {
    VolumeURL urlOrigin(url);
    std::string fileName = urlOrigin.getPath();

    std::vector<VolumeURL> result;

    QuadHidacHeader header = readHeader(fileName);

    float t = 0.0f;
    for(int i=0; i<header.ntf; i++) {
        VolumeURL origin("i4d", fileName);
        origin.addSearchParameter("frame", itos(i));
        origin.getMetaDataContainer().addMetaData("frame", new IntMetaData(static_cast<int>(i)));
        //origin.getMetaDataContainer().addMetaData("timestep", new FloatMetaData(static_cast<int>(t)));
        result.push_back(origin);

        if((size_t)i < header.dt.size())
            t += header.dt[i];
    }

    return result;
}

VolumeReader* QuadHidacVolumeReader::create(ProgressBar* progress) const {
    return new QuadHidacVolumeReader(progress);
}

} // namespace voreen
