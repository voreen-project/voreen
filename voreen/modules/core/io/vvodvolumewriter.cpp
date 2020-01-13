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

#include "vvodvolumewriter.h"
#include "vvodformat.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/voreenapplication.h"

#include "tgt/filesystem.h"
#include "tgt/matrix.h"

namespace voreen {

const std::string CACHE_SUBDIR =             "OctreeCreator";
const std::string OCTREE_FILENAME =          "octree.xml";
const std::string BRICK_BUFFER_SUBDIR =      "brickBuffer";
const std::string BRICK_BUFFER_FILE_PREFIX = "buffer_";


const std::string VvodVolumeWriter::loggerCat_("voreen.io.VvodVolumeWriter");

VvodVolumeWriter::VvodVolumeWriter(ProgressReporter* p)
    : progressBar_(p)
{
    extensions_.push_back("vvod");
}

std::string VvodVolumeWriter::getOctreeStoragePath(const VolumeOctree* octree) const {
    return VoreenApplication::app()->getCachePath(CACHE_SUBDIR + "/" + octree->getOctreeConfigurationHash());
}

void VvodVolumeWriter::write(const std::string& filename, const VolumeBase* volumeHandle) {
    tgtAssert(volumeHandle, "No volume");

    if (filename.empty()) {
        LWARNING("Empty filename");
        return;
    }

    // get octree representation
    const VolumeOctree* volume = volumeHandle->getRepresentation<VolumeOctree>();
    if (!volume) {
        LWARNING("No octree representation found, nothing to save");
        return;
    }

    // this should not happen as it is already checked by the octree save processor
    if (volumeHandle->getOrigin().getPath() == filename)
        throw tgt::IOException("Cannot overwrite volume with itself: " + filename);

    // determine the old octree file
    std::string octreeFile;
    if (tgt::FileSystem::fileExtension(volumeHandle->getOrigin().getPath()) == "vvod") {
        // volume is a stored octree
        octreeFile = tgt::FileSystem::cleanupPath(volumeHandle->getOrigin().getPath());
    }
    else {
        // determine cached octree.xml file
        octreeFile = tgt::FileSystem::cleanupPath(getOctreeStoragePath(volume) + "/" + OCTREE_FILENAME);
    }

    if (octreeFile.empty() || !tgt::FileSystem::fileExists(octreeFile)) {
        LWARNING("Found no cached or stored octree file: " + octreeFile);
        return;
    }
    std::string vvodName = tgt::FileSystem::cleanupPath(filename);

    //determine prefix
    std::string prefix = tgt::FileSystem::baseName(vvodName);
    std::string oldNodeBufferString, newNodeBufferString;
    std::string oldBrickBufferDir, newBrickBufferString;

    // create storage object for octree and copy the meta data into it
    VvodStorageObject storageObject;
    storageObject.volumeOctree_ = const_cast<VolumeOctree*>(volume);
    //if (const Volume* tmpVol = dynamic_cast<const Volume*>(volumeHandle))
    //    storageObject.metaData_ = MetaDataContainer(tmpVol->getMetaDataContainer());
    std::vector<std::string> keys = volumeHandle->getMetaDataKeys();
    for(size_t i=0; i<keys.size(); i++) {
        const MetaDataBase* md = volumeHandle->getMetaData(keys[i]);
        if(md) {
            storageObject.metaData_.addMetaData(keys[i], md->clone());
        }
    }


    // create file and serialize the storage object into it
    XmlSerializer serializer(vvodName);
    try {
        serializer.serialize("OctreeVolume", storageObject);
    }
    catch (SerializationException& e) {
        LWARNING("Failed to serialize octree volume: " << e.what());
        return;
    }

    // write serialization to stream
    std::ostringstream textStream;
    try {
        serializer.write(textStream);
        if (textStream.fail()) {
            LWARNING("Failed to write octree volume serialization to string stream");
            return;
        }
    }
    catch (std::exception& e) {
        LWARNING("Failed to write octree volume serialization to string stream: " << e.what());
        return;
    }

    // now we have a valid string stream containing the serialized octree
    // => open output file and write it to the file
    std::fstream fileStream(vvodName.c_str(), std::ios_base::out);
    if (fileStream.fail()) {
        LWARNING("Failed to open file '" << vvodName << "' for writing.");
        return;
    }

    try {
        fileStream << textStream.str();
    }
    catch (std::exception& e) {
        LWARNING("Failed to write serialization data stream to file '" << vvodName << "': " << std::string(e.what()));
    }
    fileStream.close();

    // now we should

    //tgt::FileSystem::copyFile(octreeFile, vvodName);

        //alter xml TODO: instead of altering the xml file, the VvodStorageObject should already contain the correct information !!!
        TiXmlDocument doc(vvodName);
        if(!doc.LoadFile())
            throw tgt::IOException(vvodName + " could not be loaded!");

        TiXmlHandle curRoot(&doc);
        TiXmlElement* pElem;
        pElem = curRoot.FirstChild("VoreenData").Element();
        if (!pElem) throw tgt::IOException("No VoreenData element in " + vvodName);
        else { curRoot = TiXmlHandle(pElem); }
        pElem = curRoot.FirstChild("OctreeVolume").Element();
        if (!pElem) throw tgt::IOException("No OctreeVolume element in " + vvodName);
        else { curRoot = TiXmlHandle(pElem); }
        pElem = curRoot.FirstChild("Octree").Element();
        if (!pElem) throw tgt::IOException("No Octree element in " + vvodName);
        else { curRoot = TiXmlHandle(pElem); }
        //rename nodebuffer
        pElem = curRoot.FirstChild("nodeBufferName").Element();
        if (!pElem) throw tgt::IOException("No nodeBufferName element in " + vvodName);
        else {
            oldNodeBufferString = std::string(pElem->Attribute("value"));
            newNodeBufferString = std::string(prefix + "_" + oldNodeBufferString);
            pElem->SetAttribute("value",newNodeBufferString);
        }
        //rename brickbuffer
        pElem = curRoot.FirstChild("brickPoolManager").Element();
        if (!pElem) throw tgt::IOException("No brickPoolManager element in " + vvodName);
        else { curRoot = TiXmlHandle(pElem); }
        pElem = curRoot.FirstChild("brickPoolPath").Element();
        if (!pElem) throw tgt::IOException("No brickPoolPath element in " + vvodName);
        else {
            if (tgt::FileSystem::isAbsolutePath(pElem->Attribute("value")))
                oldBrickBufferDir = tgt::FileSystem::cleanupPath(pElem->Attribute("value"), true);
            else
                oldBrickBufferDir = tgt::FileSystem::cleanupPath(tgt::FileSystem::parentDir(/*octreeFile*/vvodName) + "/" + pElem->Attribute("value"), true);
            newBrickBufferString = prefix + "_" + BRICK_BUFFER_SUBDIR;
            pElem->SetAttribute("value",/*"../"+ */newBrickBufferString);
        }
        //rename bricks
        pElem = curRoot.FirstChild("bufferFiles").Element();
        if (!pElem) throw tgt::IOException("No bufferFiles element in " + vvodName);
        else { curRoot = TiXmlHandle(pElem); }
        pElem = curRoot.FirstChild("item").Element();
        size_t pos = 0;
        if(pElem) {
            pos = std::string(pElem->Attribute("value")).find(BRICK_BUFFER_SUBDIR + "/" + BRICK_BUFFER_FILE_PREFIX);
            tgtAssert(pos != std::string::npos,"unknown buffer file name");
        }
        for( /*pElem*/; pElem != nullptr; pElem=pElem->NextSiblingElement()) {
            pElem->SetAttribute("value",/*"../"+*/ prefix + "_" + std::string(pElem->Attribute("value")).substr(pos));
        }
        //save file
        if(!doc.SaveFile())
            throw tgt::IOException(vvodName + " could not be saved!");

    // rename nodebuffer.raw
    std::string nodeBufferName = tgt::FileSystem::cleanupPath(tgt::FileSystem::parentDir(vvodName) + "/" + oldNodeBufferString);
    std::string newNodeBufferName = tgt::FileSystem::cleanupPath(tgt::FileSystem::parentDir(vvodName) + "/" + newNodeBufferString);
    // delete an already existing file if necessary
    if (tgt::FileSystem::fileExists(newNodeBufferName))
        tgt::FileSystem::deleteFile(newNodeBufferName);
    bool renameSuccessful = tgt::FileSystem::renameFile(nodeBufferName, newNodeBufferName);
    if (!renameSuccessful)
        throw tgt::IOException("Could not rename file " + nodeBufferName + " to " + newNodeBufferName);

    //copy brickBuffer
    std::string newBrickBufferDir = tgt::FileSystem::parentDir(vvodName) + "/" + newBrickBufferString;
    tgt::FileSystem::createDirectory(newBrickBufferDir);
    tgt::FileSystem::clearDirectory(newBrickBufferDir);
    if(tgt::FileSystem::dirExists(oldBrickBufferDir)) {
        std::vector<std::string> buffers = tgt::FileSystem::listFiles(oldBrickBufferDir);
        if (progressBar_)
            progressBar_->setProgress(0.f);
        for(size_t i = 0; i < buffers.size(); i++) {
            tgt::FileSystem::copyFile(oldBrickBufferDir + "/" + buffers[i], newBrickBufferDir + "/" + buffers[i]);
            if (progressBar_)
                progressBar_->setProgress(std::min(static_cast<float>(i+1)/static_cast<float>(buffers.size()), 0.99f));
        }
    } else
        throw tgt::IOException("Missing brick buffer directory: " + oldBrickBufferDir);

    if (progressBar_)
        progressBar_->setProgress(1.f);
    LINFO("Saving octree complete");
}

VolumeWriter* VvodVolumeWriter::create(ProgressBar* progress) const {
    return new VvodVolumeWriter(progress);
}

} // namespace voreen
