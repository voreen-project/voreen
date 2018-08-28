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

#include "commands_convert.h"
#include "voreen/core/io/volumeserializer.h"
#include "voreen/core/io/volumeserializerpopulator.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumecollection.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorconvert.h"
#include "tgt/vector.h"

#ifdef VRN_WITH_DEVIL
#include <IL/il.h>
#include <IL/ilu.h>
#endif

namespace voreen {

#ifdef VRN_WITH_DEVIL

CommandStackImg::CommandStackImg() :
    Command("--stackimg", "", "Stack multiple image files (PNG, JPEG, ...) to create a new volume data set.",
            "image [image ...] output", -1)
{
    loggerCat_ += "." + name_;
}

bool CommandStackImg::checkParameters(const std::vector<std::string>& parameters) {
    return (parameters.size() >= 2);
}

bool CommandStackImg::execute(const std::vector<std::string>& parameters) {
    VolumeRAM* targetDataset_;
    tgt::ivec3 dimensions;
    dimensions.z = parameters.size()-1;

    ilInit();
    iluInit();
    ILuint ImageName;

    //load first image to check dimensions...
    ilGenImages(1, &ImageName);
    ilBindImage(ImageName);
    ilLoadImage((ILstring)parameters[0].c_str());
    ILuint Width, Height, Format, Type;
    Width = ilGetInteger(IL_IMAGE_WIDTH);
    Height = ilGetInteger(IL_IMAGE_HEIGHT);
    Format = ilGetInteger(IL_IMAGE_FORMAT);
    Type = ilGetInteger(IL_IMAGE_TYPE);

    switch(Format) {
        case IL_RGB : LDEBUG("Format: IL_RGB");
            break;
        case IL_RGBA : LDEBUG("Format: IL_RGBA");
            break;
        case IL_BGR : LDEBUG("Format: IL_BGR");
            break;
        case IL_BGRA : LDEBUG("Format: IL_BGRA");
            break;
        case IL_LUMINANCE: LDEBUG("Format: IL_LUMINANCE");
            break;
    }

    switch(Type) {
        case IL_UNSIGNED_BYTE : LDEBUG("Type: IL_UNSIGNED_BYTE");
            break;
        case IL_BYTE : LDEBUG("Type: IL_BYTE");
            break;
        case IL_UNSIGNED_SHORT : LDEBUG("Type: IL_UNSIGNED_SHORT");
            break;
        case IL_SHORT : LDEBUG("Type: IL_SHORT");
            break;
        case IL_INT : LDEBUG("Type: IL_INT");
            break;
        case IL_UNSIGNED_INT : LDEBUG("Type: IL_UNSIGNED_INT");
            break;
    }

    ILenum Error;
    if ((Error = ilGetError()) != IL_NO_ERROR) {
        LERROR(Error << iluErrorString(Error));
        ilDeleteImages(1, &ImageName);
        return false;
    }
    dimensions.x = Width;
    dimensions.y = Height;
    LINFO("stacking " << dimensions.z << " images with dimensions: " << dimensions);
    ilDeleteImages(1, &ImageName);

    std::string fname;  //current filename
    targetDataset_ = new VolumeRAM_UInt8(dimensions);
    uint8_t* p = (uint8_t*)targetDataset_->getData();

    for (int i=0; i < dimensions.z; i++) {
        fname = parameters[i];
        ilGenImages(1, &ImageName);
        ilBindImage(ImageName);
        ilLoadImage((ILstring)fname.c_str());
        ilCopyPixels(0, 0, 0, dimensions.x, dimensions.y, 1, IL_LUMINANCE, IL_UNSIGNED_BYTE, p);
        p += dimensions.x * dimensions.y;
        ilDeleteImages(1, &ImageName);
    }

    if (targetDataset_) {
        VolumeSerializerPopulator volLoadPop;
        const VolumeSerializer* serializer = volLoadPop.getVolumeSerializer();
        serializer->write(parameters.back(), targetDataset_);
    }
    delete targetDataset_;
    return true;
}

#endif // VRN_WITH_DEVIL

//-----------------------------------------------------------------------------

CommandStackRaw::CommandStackRaw() /*:
    Command("--stackraw", "", "Stack raw-images given by sourceFilename(s) \n \n"
            "\t\tInput images with size DX*DY\n"
            "\t\tHEADERSIZE bytes at the beginning of each file are skipped.\n"
            "\t\tFormats:\n"
            "\t\tgs16: 16 bit grayscale\n"
            "\t\trgbni: 24 bit rgb, non interleaved",
            "DX DY HEADERSIZE <gs16|rgbni> IN [IN ...] OUT", -1) */
{
//    loggerCat_ += "." + name_;
}

bool CommandStackRaw::checkParameters(const std::vector<std::string>& parameters) {
    return (parameters.size() >= 6);
}

bool CommandStackRaw::execute(const std::vector<std::string>& parameters) {
    /*VolumeRAM* targetDataset_;

    tgt::ivec3 dimensions;
    dimensions.x = cast<int>(parameters[0]);
    dimensions.y = cast<int>(parameters[1]);

    dimensions.z = parameters.size()-5;

    int headersize;
    headersize = cast<int>(parameters[2]);

    std::string format = parameters[3];
    LINFO("stacking " << dimensions.z << " images with dimensions: " << dimensions);
    LINFO("skipping header with size " << headersize);
    LINFO("format: " << format);

    if (format == "gs16") {
        targetDataset_ = new VolumeRAM_UInt16(dimensions);
        uint8_t* p = (uint8_t*)targetDataset_->getData();
        std::ifstream fin;
        int readsize = dimensions.x * dimensions.y * 2;
        std::string fname;  //current filename

        for (int i=0; i < dimensions.z; i++) {
            fname = parameters[i+4];
            fin.open(fname.c_str(), std::ios::in | std::ios::binary);
            if (!fin.good()) {
                LERROR("failed to open " << fname);
                break;
            }
            else
                LINFO("opened: " << fname);
            //skip header:
            fin.seekg(headersize, std::ios::beg );
            fin.read(reinterpret_cast<char*>(p), readsize);
            LDEBUG("read " << fin.gcount() << " of " << readsize <<  " bytes");
            p += readsize;
            if ( fin.bad() ) {
                LERROR("Error reading data");
                exit( 0 );
            }
            fin.close();
        }
    }
    else if (format == "rgbni") {
        targetDataset_ = new VolumeRAM_4xUInt8(dimensions);
        uint8_t* p = (uint8_t*)targetDataset_->getData();
        std::ifstream fin;
        int readsize = dimensions.x * dimensions.y;
        uint8_t* buffer = new uint8_t[readsize];
        std::string fname;  //current filename

        for (int i=0; i < dimensions.z; i++) {
            fname = parameters[i+4];
            fin.open(fname.c_str(), std::ios::in | std::ios::binary);
            if (!fin.good()) {
                LERROR("failed to open " << fname);
                break;
            } else {
                LDEBUG("opened: " << fname);
            }
            //skip header:
            fin.seekg(headersize, std::ios::beg);
            //loop over rgb:
            for (int c=0; c<3; c++) {
                fin.read(reinterpret_cast<char*>(buffer), readsize);
                LDEBUG("read " << fin.gcount() << " of " << readsize <<  " bytes");
                if ( fin.bad() ) {
                    LERROR("Error reading data");
                    return false;
                }
                uint8_t* colorp = p;
                for (int j=0; j<readsize; j++) {
                    colorp[c] = buffer[j];
                    colorp += 4;
                }
            }
            //set alpha to opaque:
            uint8_t* colorp = p;
            for (int j=0; j<readsize; j++) {
                //calculate greyscale value for alpha channel:
                //colorp[3] = 0.3 * (float) colorp[0] + 0.59 * (float) colorp[1] + 0.11 * (float) colorp[2];
                colorp[3] = 255;
                colorp += 4;
            }
            fin.close();
            p += readsize*4;
        }
        delete[] buffer;
    }
    else
        throw tgt::Exception("Unknown format!");

    VolumeSerializerPopulator volLoadPop;
    const VolumeSerializer* serializer = volLoadPop.getVolumeSerializer();
    Volume targetHandle(targetDataset_, tgt::vec3(1.0f), tgt::vec3(0.0f));
    oldVolumePosition(&targetHandle);
    serializer->write(parameters.back(), &targetHandle);

    //delete targetDataset_; //gets deleted by handle */
    return true;
}

//-----------------------------------------------------------------------------

CommandConvert::CommandConvert() /*:
    Command("--convert", "",
            "Convert Volume Datasets\n"
            "\t\tTARGET:\n"
            "\t\t8: 8 bit DS\n"
            "\t\t12: 12 bit DS\n"
            "\t\t16: 16 bit DS",
            "<TARGET IN OUT>", 3) */
{
//    loggerCat_ += "." + name_;
}

bool CommandConvert::execute(const std::vector<std::string>& parameters) {
    Volume* targetDataset_ = 0;

    VolumeSerializerPopulator volLoadPop;
    const VolumeSerializer* serializer = volLoadPop.getVolumeSerializer();

    VolumeCollection* volumeCollection = serializer->read(parameters[1]);
    const VolumeBase* sourceDataset_ = volumeCollection->first();

    VolumeOperatorConvert voConvert;
    if (parameters[0] == "8") {
        //targetDataset_ = new VolumeRAM_UInt8(sourceDataset_->getDimensions(), sourceDataset_->getSpacing(), sourceDataset_->getTransformation());
        targetDataset_ = voConvert.apply<uint8_t>(sourceDataset_);
    }
    else if (parameters[0] == "16") {
        //targetDataset_ = new VolumeRAM_UInt16(sourceDataset_->getDimensions(), sourceDataset_->getSpacing(), sourceDataset_->getTransformation());
        targetDataset_ = voConvert.apply<uint16_t>(sourceDataset_);
    }
    else {
        delete sourceDataset_;
        throw tgt::Exception("Unknown target!");
    }

    serializer->write(parameters.back(), targetDataset_);

    delete volumeCollection;
    delete targetDataset_;
    return true;
}

//-----------------------------------------------------------------------------

CommandConvertFormat::CommandConvertFormat() /*:
    Command("--convertformat", "",
            "Convert Volume Datasets\n",
            "<IN OUT>", 2) */
{
//    loggerCat_ += "." + name_;
}

bool CommandConvertFormat::execute(const std::vector<std::string>& parameters) {
    /*VolumeSerializerPopulator volLoadPop;
    const VolumeSerializer* serializer = volLoadPop.getVolumeSerializer();

    VolumeCollection* volumeCollection = serializer->read(parameters[0]);
    const VolumeBase* sourceDataset = volumeCollection->first();

    serializer->write(parameters.back(), sourceDataset);

    delete volumeCollection;
    //delete serializer; */
    return true;
}

} // namespace voreen
