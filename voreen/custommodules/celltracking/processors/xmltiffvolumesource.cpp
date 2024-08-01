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

#include "xmltiffvolumesource.h"

#include "tgt/filesystem.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "tiffio.h"

namespace voreen {

XMLTiffVolumeSource::XMLTiffVolumeSource()
    : Processor()
    , outport_(Port::OUTPORT, "outport", "Outport")
    , selectedFile_("selectedFile", "File to load", "File to load", "")
    , timeStep_("timeStep", "Time step", 1, 1, 100)
    , channel_("channel", "Channel", 1, 1, 100)
    , autoload_("autoload", "Autoload volume", false)
    , loadButton_("loadButton", "Load file")
    , shouldLoadFileOnProcess_(false)
{
    addPort(outport_);
    addProperty(selectedFile_);
    addProperty(timeStep_);
    addProperty(channel_);
    addProperty(autoload_);
    addProperty(loadButton_);

    ON_CHANGE(selectedFile_, XMLTiffVolumeSource, hasSelectedFile);
    ON_CHANGE(loadButton_, XMLTiffVolumeSource, loadButtonPressed);

    
    hasSelectedFile();
}

Processor* XMLTiffVolumeSource::create() const {
    return new XMLTiffVolumeSource();
}

std::string XMLTiffVolumeSource::getClassName() const {
    return "XMLTiffVolumeSource";
}

std::string XMLTiffVolumeSource::getCategory() const {
    return "Input";
}

void XMLTiffVolumeSource::setDescriptions() {
    setDescription("Minimal sample processor that appends a user-defined prefix to a given text.");
}

void XMLTiffVolumeSource::process() {
    if (!shouldLoadFileOnProcess_ && !autoload_.get())
        return;
    Volume* vol = loadVolume();
    outport_.setData(vol);
    shouldLoadFileOnProcess_ = false;
}

void XMLTiffVolumeSource::hasSelectedFile(void)
{ 
    try{
        loadXMLFile(tgt::FileSystem::cleanupPath(selectedFile_.get()));
        loadButton_.setReadOnlyFlag(false);
    }catch(...){
        loadButton_.setReadOnlyFlag(true);
    }
}

void XMLTiffVolumeSource::loadButtonPressed(void)
{
    shouldLoadFileOnProcess_ = true;
}

void XMLTiffVolumeSource::loadXMLFile(const std::string path)
{   
    TiXmlDocument doc(path);   

    if(!doc.LoadFile())
        throw tgt::IOException(path + " could not be loaded!");
    TiXmlNode* curRoot = doc.FirstChildElement("ExportDocument");

    documentDesc_.clear();
    int maxchannel = -1;
    int maxz = -1;
    int maxtime = -1;

    for(TiXmlNode* node = curRoot->FirstChildElement("Image"); node; node = node->NextSibling()){
        std::string filename;
        int channel;
        int timeStep;
        int z;

        TiXmlElement* filenameNode = node->FirstChildElement("Filename");
        filename = filenameNode->GetText();

        TiXmlElement* boundsNode = node->FirstChildElement("Bounds");

        boundsNode->Attribute("StartZ", &z);
        boundsNode->Attribute("StartC", &channel);
        boundsNode->Attribute("StartT", &timeStep);
        boundsNode->Attribute("SizeX", &sizeX_);
        boundsNode->Attribute("SizeY", &sizeY_);
        
        maxchannel = std::max(maxchannel, channel);
        maxtime = std::max(maxtime, timeStep);
        maxz = std::max(maxz, z);

        documentDesc_[SliceDesc(z, channel, timeStep)] = filename;
    }

    if (maxz == -1 || maxchannel == -1 ||maxtime == -1){
        throw tgt::IOException("Bad file");
    }
    timeStep_.setMaxValue(maxtime+1);
    channel_.setMaxValue(maxchannel+1);
    sizeZ_ = maxz;
}

Volume* XMLTiffVolumeSource::loadVolume()
{
    VolumeRAM_UInt8* volumeram = new VolumeRAM_UInt8(tgt::svec3(sizeX_, sizeY_, sizeZ_));
    Volume *volume = new Volume(volumeram, tgt::vec3(1.0f), tgt::vec3(0.0f));

    int channel = channel_.get()-1;
    int timestep = timeStep_.get()-1;
    std::string path = tgt::FileSystem::dirName(selectedFile_.get())+"/";

    for(int slice = 0; slice != sizeZ_; slice++){
        std::string basename = documentDesc_[SliceDesc(slice, channel, timestep)];
        std::string filename = path+basename;
        TIFF* tif = TIFFOpen(filename.c_str(), "r");
        if (!tif) {
            LERROR("Failed to open TIFF stack");
            throw tgt::IOException("Failed to open TIFF stack", filename);
        }
        
        // load tiff file into buffer 'img'
        uint16_t depth, bps;
        TIFFGetField(tif, TIFFTAG_SAMPLESPERPIXEL, &depth);
        TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bps);
        size_t stripSize = TIFFStripSize (tif);
        int stripMax = TIFFNumberOfStrips (tif);
        std::vector<uint8_t> img(stripMax*stripSize);
        int offset = 0;
        for (int stripCount = 0; stripCount < stripMax; stripCount++) {
            int pixelsRead;
            if ((pixelsRead= static_cast<int>(TIFFReadEncodedStrip (tif, stripCount, img.data()+ offset, stripSize))) == -1){
                    LERROR("Read error on input strip number " << stripCount);
                    TIFFClose(tif);
                    throw tgt::CorruptedFileException("Read error on input strip number " + itos(stripCount));
            }
            offset+= pixelsRead;
        
        }

        // extract channel to volume
        uint8_t* out = static_cast<uint8_t*>(volumeram->getData())+sizeX_*sizeY_*slice;
        uint8_t* in = img.data();
        if (channel == 0)
            in+=1; // first channel is in green component
        if (channel == 1)
            in+=0; // second channel is in the red component
        for(int i = 0; i != sizeX_*sizeY_; i++){
            *out = *in;
            out+=1;
            in+=3;
        }

    }

    return volume;
}

XMLTiffVolumeSource::SliceDesc::SliceDesc(int z, int channel, int timeStep)
    : z_(z)
    , channel_(channel)
    , timeStep_(timeStep)
{
}

bool XMLTiffVolumeSource::SliceDesc::operator<(const SliceDesc & other) const
{
    if (z_ < other.z_)
        return true;
    if (z_ > other.z_)
        return false;
    if (channel_ < other.channel_)
        return true;
    if (channel_ > other.channel_)
        return false;
    if (timeStep_ < other.timeStep_)
        return true;
    else 
        return false;
}



 } // namespace
