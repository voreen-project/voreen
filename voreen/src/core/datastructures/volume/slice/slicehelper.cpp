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

#include "voreen/core/datastructures/volume/slice/slicehelper.h"

#include "voreen/core/datastructures/geometry/trianglemeshgeometry.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/octree/volumeoctreebase.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorresample.h"
//#include "voreen/core/utils/hashing.h"
//#include <boost/thread/locks.hpp>
//#include <set>

namespace voreen {

const std::string SliceHelper::loggerCat_("voreen.SliceHelper");

SliceTexture* SliceHelper::getVolumeSlice(const VolumeBase* volume, SliceAlignment alignment, size_t sliceIndex, int* shiftArray /*= 0*/,
    size_t levelOfDetail /*= 0*/, clock_t timeLimit /*= 0*/, bool* complete /*= 0*/)
{
    tgtAssert(volume, "null pointer passed");
    tgtAssert(alignment != UNALIGNED_PLANE, "invalid alignment");
    tgtAssert(sliceIndex < volume->getDimensions()[alignment], "invalid slice index");

    const size_t sMinusOne = -1; // define to omit compiler warnings

    tgt::vec3 urb = volume->getURB();
    tgt::vec3 llf = volume->getLLF();
    tgt::vec3 sp = volume->getSpacing();
    tgt::Bounds b(llf, urb);
    tgt::svec3 dims = volume->getDimensions();

    tgt::vec3 bb_urb = b.getURB();
    tgt::vec3 bb_llf = b.getLLF();

    tgt::vec3 origin(0.0f);
    tgt::vec3 xVec(0.0f);
    tgt::vec3 yVec(0.0f);

    tgt::svec2 sliceDim;

    void* dataBuffer = 0;
    GLint textureFormat, internalFormat;
    GLenum textureDataType;

    size_t samplingLOD = 0;

    switch(alignment) {
        case YZ_PLANE: {
            // representation preference:
            // 1) use RAM volume, if present
            // 2) else: use octree, if present
            // 3) else: try to create RAM representation
            // note: disk representation is not usable for this alignment
            sliceDim = volume->getDimensions().yz();
            tgt::svec3 sliceDim3D(1, sliceDim);     //< expected dimension of the extracted slice
            const VolumeRAM* ramVolume3D = 0;       //< RAM representation
            const VolumeRAM* ramVolumeSlice = 0;    //< slice retrieved from octree (has to be deleted afterwards)
            if (volume->hasRepresentation<VolumeRAM>()) {
                ramVolume3D = volume->getRepresentation<VolumeRAM>();
            }
            else if (volume->hasRepresentation<VolumeOctreeBase>()) {
                try {
                    const VolumeOctreeBase* octree = volume->getRepresentation<VolumeOctreeBase>();
                    tgtAssert(octree, "no octree");
                    if(!shiftArray) {
                        ramVolumeSlice = octree->createSlice(YZ_PLANE, sliceIndex, levelOfDetail, timeLimit, complete);
                        samplingLOD = levelOfDetail;
                        tgtAssert(ramVolumeSlice, "null pointer returned (exception expected)");
                    } else {
                        //get all needed RAM slices
                        VolumeRAM** sliceArray = static_cast<VolumeRAM**>(malloc(sizeof(VolumeRAM*) * volume->getNumChannels()));
                        bool tmpComp = true, chComp = false; std::set<size_t> noDelete;
                        size_t firstSlice = -1;
                        for(size_t c = 0; c < volume->getNumChannels(); c++) {
                            sliceArray[c] = 0; //set null
                            //improvement
                                bool impro = false;
                                for(size_t d = 0; d < c; d++) {
                                    if(shiftArray[d] == shiftArray[c]) {
                                        sliceArray[c] = sliceArray[d];
                                        noDelete.insert(c);
                                        impro = true;
                                        break;
                                    }
                                }
                                if(impro)
                                    continue;
                            //get slice or set zero
                            int tmpSlice = (int)sliceIndex + shiftArray[c];
                            if( tmpSlice < 0 || tmpSlice >= static_cast<int>(dims.x))
                                sliceArray[c] = 0;
                            else {
                                sliceArray[c] = octree->createSlice(YZ_PLANE, tmpSlice, levelOfDetail, timeLimit/static_cast<clock_t>(volume->getNumChannels()), &chComp);
                                samplingLOD = levelOfDetail;
                                tmpComp &= chComp;
                                if(firstSlice == sMinusOne /*-1*/)
                                    firstSlice = c;
                            }
                        }
                        if(complete) *complete = tmpComp;
                        //compose slices
                        if(firstSlice == sMinusOne /*-1*/) {
                            //no slices have been created
                            switch (volume->getNumChannels()) {
                            case 1:
                                ramVolumeSlice = new VolumeRAM_UInt16(sliceDim3D, true);
                                break;
                            case 2:
                                ramVolumeSlice = new VolumeRAM_2xUInt16(sliceDim3D, true);
                                break;
                            case 3:
                                ramVolumeSlice = new VolumeRAM_3xUInt16(sliceDim3D, true);
                                break;
                            case 4:
                                ramVolumeSlice = new VolumeRAM_4xUInt16(sliceDim3D, true);
                                break;
                            default:
                                tgtAssert(false, "more than 4 channels");
                            }
                            std::fill_n(reinterpret_cast<uint16_t*>(const_cast<void*>(ramVolumeSlice->getData())),
                                        ramVolumeSlice->getNumVoxels()*ramVolumeSlice->getNumChannels(),(uint16_t)0);
                        } else {
                            //init with first slice
                            ramVolumeSlice = sliceArray[firstSlice];
                            uint16_t* modBuffer = reinterpret_cast<uint16_t*>(const_cast<void*>(ramVolumeSlice->getData()));
                            for(size_t c = 0; c < volume->getNumChannels(); c++) {
                                if(c == firstSlice)
                                    continue; // value is already been set
                                for(size_t x = 0; x < ramVolumeSlice->getDimensions().x; x++) {
                                    for(size_t y = 0; y < ramVolumeSlice->getDimensions().y; y++) {
                                        modBuffer[(y*ramVolumeSlice->getDimensions().x+x)*ramVolumeSlice->getNumChannels()+c] =
                                            (sliceArray[c] ? static_cast<uint16_t>(sliceArray[c]->getVoxelNormalized(x,y,0,c)*65535.f) : 0);
                                    }
                                }
                            }
                        }
                        //clean up
                        for(size_t c = 0; c < volume->getNumChannels(); c++) {
                            if(sliceArray[c] == sliceArray[firstSlice])
                                continue;
                            if(noDelete.find(c) == noDelete.end())
                                delete sliceArray[c];
                        }
                        free(sliceArray);
                    }
                }
                catch (tgt::Exception& e) {
                    LERROR("Failed to create YZ slice from octree: " << e.what());
                    return 0;
                }
            }
            else {
                ramVolume3D = volume->getRepresentation<VolumeRAM>();
                if (!ramVolume3D) {
                    LERROR("Failed to create YZ slice: neither RAM nor octree representation available");
                    return 0;
                }
            }

            float x = static_cast<float>(sliceIndex);
            float xcoord = llf.x + (x+0.5f) * sp.x; // We want our slice to be in the center of voxels

            origin = tgt::vec3(xcoord, bb_llf.y, bb_llf.z);
            xVec = tgt::vec3(xcoord, bb_urb.y, bb_llf.z);
            yVec = tgt::vec3(xcoord, bb_llf.y, bb_urb.z);

            tgtAssert(ramVolume3D || ramVolumeSlice, "no volume");
            try {
                if (ramVolume3D)
                    SliceHelper::extractAlignedSlicePixelDataHelper(ramVolume3D, alignment, sliceIndex, shiftArray, dataBuffer, textureFormat, internalFormat, textureDataType);
                else {
                    //shift Array has already been applied
                    SliceHelper::extractAlignedSlicePixelDataHelper(ramVolumeSlice, alignment, 0, 0, dataBuffer, textureFormat, internalFormat, textureDataType);
                }
            }
            catch (VoreenException& e) {
                LERROR(e.what());
                delete ramVolumeSlice;
                return 0;
            }

            delete ramVolumeSlice;
            ramVolumeSlice = 0;
        }
        break;

        case XZ_PLANE: {
            // representation preference:
            // 1) use RAM volume, if present
            // 2) else: use octree, if present
            // 3) else: try to create RAM representation
            // note: disk representation is not usable for this alignment
            sliceDim = tgt::svec2(volume->getDimensions().x, volume->getDimensions().z);
            tgt::svec3 sliceDim3D(sliceDim.x, 1, sliceDim.y);   //< expected slice dimension
            const VolumeRAM* ramVolume3D = 0;                   //< RAM representation
            const VolumeRAM* ramVolumeSlice = 0;                //< slice retrieved from octree (has to be deleted afterwards)
            if (volume->hasRepresentation<VolumeRAM>()) {
                ramVolume3D = volume->getRepresentation<VolumeRAM>();
            }
            else if (volume->hasRepresentation<VolumeOctreeBase>()) {
                try {
                   const VolumeOctreeBase* octree = volume->getRepresentation<VolumeOctreeBase>();
                    tgtAssert(octree, "no octree");
                    if(!shiftArray) {
                        ramVolumeSlice = octree->createSlice(XZ_PLANE, sliceIndex, levelOfDetail, timeLimit, complete);
                        samplingLOD = levelOfDetail;
                        tgtAssert(ramVolumeSlice, "null pointer returned (exception expected)");
                    } else {
                        //get all needed RAM slices
                        VolumeRAM** sliceArray = static_cast<VolumeRAM**>(malloc(sizeof(VolumeRAM*) * volume->getNumChannels()));
                        bool tmpComp = true, chComp = false; std::set<size_t> noDelete;
                        size_t firstSlice = -1;
                        for(size_t c = 0; c < volume->getNumChannels(); c++) {
                            sliceArray[c] = 0; //set null
                            //improvement
                                bool impro = false;
                                for(size_t d = 0; d < c; d++) {
                                    if(shiftArray[d] == shiftArray[c]) {
                                        sliceArray[c] = sliceArray[d];
                                        noDelete.insert(c);
                                        impro = true;
                                        break;
                                    }
                                }
                                if(impro)
                                    continue;
                            //get slice or set zero
                            int tmpSlice = (int)sliceIndex + shiftArray[c];
                            if( tmpSlice < 0 || tmpSlice >= static_cast<int>(dims.y))
                                sliceArray[c] = 0;
                            else {
                                sliceArray[c] = octree->createSlice(XZ_PLANE, tmpSlice, levelOfDetail, timeLimit/static_cast<clock_t>(volume->getNumChannels()), &chComp);
                                samplingLOD = levelOfDetail;
                                tmpComp &= chComp;
                                if(firstSlice == sMinusOne /*-1*/)
                                    firstSlice = c;
                            }
                        }
                        if(complete) *complete = tmpComp;
                        //compose slices
                        if(firstSlice == sMinusOne /*-1*/) {
                            //no slices have been created
                            switch (volume->getNumChannels()) {
                            case 1:
                                ramVolumeSlice = new VolumeRAM_UInt16(sliceDim3D, true);
                                break;
                            case 2:
                                ramVolumeSlice = new VolumeRAM_2xUInt16(sliceDim3D, true);
                                break;
                            case 3:
                                ramVolumeSlice = new VolumeRAM_3xUInt16(sliceDim3D, true);
                                break;
                            case 4:
                                ramVolumeSlice = new VolumeRAM_4xUInt16(sliceDim3D, true);
                                break;
                            default:
                                tgtAssert(false, "more than 4 channels");
                            }
                            std::fill_n(reinterpret_cast<uint16_t*>(const_cast<void*>(ramVolumeSlice->getData())),
                                        ramVolumeSlice->getNumVoxels()*ramVolumeSlice->getNumChannels(),(uint16_t)0);
                        } else {
                            //init with first slice
                            ramVolumeSlice = sliceArray[firstSlice];
                            uint16_t* modBuffer = reinterpret_cast<uint16_t*>(const_cast<void*>(ramVolumeSlice->getData()));
                            for(size_t c = 0; c < volume->getNumChannels(); c++) {
                                if(c == firstSlice)
                                    continue; // value is already been set
                                for(size_t x = 0; x < ramVolumeSlice->getDimensions().x; x++) {
                                    for(size_t y = 0; y < ramVolumeSlice->getDimensions().y; y++) {
                                        modBuffer[(y*ramVolumeSlice->getDimensions().x+x)*ramVolumeSlice->getNumChannels()+c] =
                                            (sliceArray[c] ? static_cast<uint16_t>(sliceArray[c]->getVoxelNormalized(x,y,0,c)*65535.f): 0);
                                    }
                                }
                            }
                        }
                        //clean up
                        for(size_t c = 0; c < volume->getNumChannels(); c++) {
                            if(sliceArray[c] == sliceArray[firstSlice])
                                continue;
                            if(noDelete.find(c) == noDelete.end())
                                delete sliceArray[c];
                        }
                        free(sliceArray);
                    }
                }
                catch (tgt::Exception& e) {
                    LERROR("Failed to create XZ slice from octree: " << e.what());
                    return 0;
                }
            }
            else {
                ramVolume3D = volume->getRepresentation<VolumeRAM>();
                if (!ramVolume3D) {
                    LERROR("Failed to create XZ slice: neither RAM nor octree representation available");
                    return 0;
                }
            }

            float y = static_cast<float>(sliceIndex);
            float ycoord = llf.y + (y+0.5f) * sp.y; // We want our slice to be in the center of voxels

            origin = tgt::vec3(bb_llf.x, ycoord, bb_llf.z);
            xVec = tgt::vec3(bb_urb.x, ycoord, bb_llf.z);
            yVec = tgt::vec3(bb_llf.x, ycoord, bb_urb.z);

            tgtAssert(ramVolume3D || ramVolumeSlice, "no volume");
            try {
                if (ramVolume3D)
                    SliceHelper::extractAlignedSlicePixelDataHelper(ramVolume3D, alignment, sliceIndex, shiftArray, dataBuffer, textureFormat, internalFormat, textureDataType);
                else
                    SliceHelper::extractAlignedSlicePixelDataHelper(ramVolumeSlice, alignment, 0, 0, dataBuffer, textureFormat, internalFormat, textureDataType);
            }
            catch (VoreenException& e) {
                LERROR(e.what());
                delete ramVolumeSlice;
                return 0;
            }

            delete ramVolumeSlice;
            ramVolumeSlice = 0;
        }
        break;

        case XY_PLANE: {
            // representation preference:
            // 1) use RAM volume, if present
            // 2) else: use octree, if present
            // 3) else: use disk volume, if present
            // 4) else: try to create RAM representation
            sliceDim = volume->getDimensions().xy();
            tgt::svec3 sliceDim3D(sliceDim, 1);    //< expected dimension of the extracted slice
            const VolumeRAM* ramVolume3D = 0;      //< RAM representation
            const VolumeRAM* ramVolumeSlice = 0;   //< slice retrieved from octree or disk volume (has to be deleted)
            if (volume->hasRepresentation<VolumeRAM>()) {
                ramVolume3D = volume->getRepresentation<VolumeRAM>();
            }
            else if (volume->hasRepresentation<VolumeOctreeBase>()) {
                try {
                   const VolumeOctreeBase* octree = volume->getRepresentation<VolumeOctreeBase>();
                    tgtAssert(octree, "no octree");
                    if(!shiftArray) {
                        ramVolumeSlice = octree->createSlice(XY_PLANE, sliceIndex, levelOfDetail, timeLimit, complete);
                        samplingLOD = levelOfDetail;
                        tgtAssert(ramVolumeSlice, "null pointer returned (exception expected)");
                    } else {
                        //get all needed RAM slices
                        VolumeRAM** sliceArray = static_cast<VolumeRAM**>(malloc(sizeof(VolumeRAM*) * volume->getNumChannels()));
                        bool tmpComp = true, chComp = false; std::set<size_t> noDelete;
                        size_t firstSlice = -1;
                        for(size_t c = 0; c < volume->getNumChannels(); c++) {
                            sliceArray[c] = 0; //set null
                            //improvement
                                bool impro = false;
                                for(size_t d = 0; d < c; d++) {
                                    if(shiftArray[d] == shiftArray[c]) {
                                        sliceArray[c] = sliceArray[d];
                                        noDelete.insert(c);
                                        impro = true;
                                        break;
                                    }
                                }
                                if(impro)
                                    continue;
                            //get slice or set zero
                            int tmpSlice = (int)sliceIndex + shiftArray[c];
                            if( tmpSlice < 0 || tmpSlice >= static_cast<int>(dims.z))
                                sliceArray[c] = 0;
                            else {
                                sliceArray[c] = octree->createSlice(XY_PLANE, tmpSlice, levelOfDetail, timeLimit/static_cast<clock_t>(volume->getNumChannels()), &chComp);
                                samplingLOD = levelOfDetail;
                                tmpComp &= chComp;
                                if(firstSlice == sMinusOne /*-1*/)
                                    firstSlice = c;
                            }
                        }
                        if(complete) *complete = tmpComp;
                        //compose slices
                        if(firstSlice == sMinusOne /*-1*/) {
                            //no slices have been created
                            switch (volume->getNumChannels()) {
                            case 1:
                                ramVolumeSlice = new VolumeRAM_UInt16(sliceDim3D, true);
                                break;
                            case 2:
                                ramVolumeSlice = new VolumeRAM_2xUInt16(sliceDim3D, true);
                                break;
                            case 3:
                                ramVolumeSlice = new VolumeRAM_3xUInt16(sliceDim3D, true);
                                break;
                            case 4:
                                ramVolumeSlice = new VolumeRAM_4xUInt16(sliceDim3D, true);
                                break;
                            default:
                                tgtAssert(false, "more than 4 channels");
                            }
                            std::fill_n(reinterpret_cast<uint16_t*>(const_cast<void*>(ramVolumeSlice->getData())),
                                        ramVolumeSlice->getNumVoxels()*ramVolumeSlice->getNumChannels(),(uint16_t)0);
                        } else {
                            //init with first slice
                            ramVolumeSlice = sliceArray[firstSlice];
                            uint16_t* modBuffer = reinterpret_cast<uint16_t*>(const_cast<void*>(ramVolumeSlice->getData()));
                            for(size_t c = 0; c < volume->getNumChannels(); c++) {
                                if(c == firstSlice)
                                    continue; // value is already been set
                                for(size_t x = 0; x < ramVolumeSlice->getDimensions().x; x++) {
                                    for(size_t y = 0; y < ramVolumeSlice->getDimensions().y; y++) {
                                        modBuffer[(y*ramVolumeSlice->getDimensions().x+x)*ramVolumeSlice->getNumChannels()+c] =
                                            (sliceArray[c] ? static_cast<uint16_t>(sliceArray[c]->getVoxelNormalized(x,y,0,c)*65535.f) : 0);
                                    }
                                }
                            }
                        }
                        //clean up
                        for(size_t c = 0; c < volume->getNumChannels(); c++) {
                            if(sliceArray[c] == sliceArray[firstSlice])
                                continue;
                            if(noDelete.find(c) == noDelete.end())
                                delete sliceArray[c];
                        }
                        free(sliceArray);
                    }
                }
                catch (tgt::Exception& e) {
                    LERROR("Failed to create XY slice from octeee: " << e.what());
                    return 0;
                }
            }
            else if (volume->hasRepresentation<VolumeDisk>() && !shiftArray) { //TODO: support channel shift
                try {
                    ramVolumeSlice = volume->getRepresentation<VolumeDisk>()->loadSlices(sliceIndex, sliceIndex);
                    tgtAssert(ramVolumeSlice, "null pointer returned (exception expected)");
                }
                catch (tgt::Exception& e) {
                    LERROR("Failed to load XY slice from disk volume: " << e.what());
                    return 0;
                }
            }
            else {
                ramVolume3D = volume->getRepresentation<VolumeRAM>();
                if (!ramVolume3D) {
                    LERROR("Failed to create XY slice: neither RAM nor octree, nor disk representation available");
                    return 0;
                }
            }

            float z = static_cast<float>(sliceIndex);
            float zcoord = llf.z + (z+0.5f) * sp.z; // We want our slice to be in the center of voxels

            origin = tgt::vec3(bb_llf.x, bb_llf.y, zcoord);
            xVec = tgt::vec3(bb_urb.x, bb_llf.y, zcoord);
            yVec = tgt::vec3(bb_llf.x, bb_urb.y, zcoord);

            tgtAssert(ramVolume3D || ramVolumeSlice, "no volume");
            try {
                if (ramVolume3D)
                    SliceHelper::extractAlignedSlicePixelDataHelper(ramVolume3D, alignment, sliceIndex, shiftArray, dataBuffer, textureFormat, internalFormat, textureDataType);
                else
                    SliceHelper::extractAlignedSlicePixelDataHelper(ramVolumeSlice, alignment, 0, 0, dataBuffer, textureFormat, internalFormat, textureDataType);
            }
            catch (VoreenException& e) {
                LERROR(e.what());
                delete ramVolumeSlice;
                return 0;
            }

            delete ramVolumeSlice;
            ramVolumeSlice = 0;
        }
        break;

        default:
            tgtAssert(false, "unknown/unsupported slice alignment");
            LERROR("Unknown/unsupported slice alignment: " << alignment);
            return 0;
    }

    tgtAssert(dataBuffer, "data buffer is empty");

    origin = volume->getPhysicalToWorldMatrix() * origin;
    xVec = volume->getPhysicalToWorldMatrix() * xVec;
    yVec = volume->getPhysicalToWorldMatrix() * yVec;
    xVec = xVec - origin;
    yVec = yVec - origin;

    return new SliceTexture(sliceDim, alignment, volume->getFormat(), volume->getBaseType(), origin, xVec, yVec,
        volume->getRealWorldMapping(), dataBuffer, textureFormat, internalFormat, textureDataType, samplingLOD);
}

SliceTexture* SliceHelper::getVolumeSlice(const VolumeBase* volume, tgt::plane pl, float samplingRate) {
    tgt::vec3 dims = volume->getDimensions();
    tgt::vec3 urb = volume->getURB();
    tgt::vec3 llf = volume->getLLF();
    tgt::vec3 center = (urb + llf) * 0.5f;

    tgt::vec3 xMax = center;
    xMax.x = urb.x;
    tgt::vec3 yMax = center;
    yMax.y = urb.y;
    tgt::vec3 zMax = center;
    zMax.z = urb.z;

    // check whether the plane normal matches one of the main directions of the volume:
    tgt::plane plVoxel = pl.transform(volume->getWorldToVoxelMatrix());
    if(fabs(fabs(dot(tgt::vec3(1.0f, 0.0f, 0.0f), plVoxel.n)) - 1.0f) < 0.001f) {
        float sliceNumber = std::round(plVoxel.d * plVoxel.n.x);
        if (sliceNumber >= 0.0f && sliceNumber < dims.x) {
            return SliceHelper::getVolumeSlice(volume, YZ_PLANE, static_cast<size_t>(sliceNumber));
        }
        return nullptr;
    }
    else if(fabs(fabs(dot(tgt::vec3(0.0f, 1.0f, 0.0f), plVoxel.n)) - 1.0f) < 0.001f) {
        float sliceNumber = std::round(plVoxel.d * plVoxel.n.y);
        if(sliceNumber >= 0.0f && sliceNumber < dims.y) {
            return SliceHelper::getVolumeSlice(volume, XZ_PLANE, static_cast<size_t>(sliceNumber));
        }
        return nullptr;
    }
    else if(fabs(fabs(dot(tgt::vec3(0.0f, 0.0f, 1.0f), plVoxel.n)) - 1.0f) < 0.001f) {
        float sliceNumber = std::round(plVoxel.d * plVoxel.n.z);
        if(sliceNumber >= 0.0f && sliceNumber < dims.z) {
            return SliceHelper::getVolumeSlice(volume, XY_PLANE, static_cast<size_t>(sliceNumber));
        }
        return nullptr;
    }
    const VolumeRAM* vol = volume->getRepresentation<VolumeRAM>();
    if(!vol)
        return nullptr;

    // transform to world coordinates:
    tgt::mat4 pToW = volume->getPhysicalToWorldMatrix();
    center = pToW * center;
    xMax = pToW * xMax;
    yMax = pToW * yMax;
    zMax = pToW * zMax;

    // project to plane:
    float d = pl.distance(center);
    center = center - (pl.n * d);
    d = pl.distance(xMax);
    xMax = xMax - (pl.n * d);
    d = pl.distance(yMax);
    yMax = yMax - (pl.n * d);
    d = pl.distance(zMax);
    zMax = zMax - (pl.n * d);

    // find max axis in plane:
    tgt::vec3 maxVec = xMax - center;
    if(distance(yMax, center) > length(maxVec))
        maxVec = yMax - center;
    if(distance(zMax, center) > length(maxVec))
        maxVec = zMax - center;

    maxVec = normalize(maxVec);
    tgt::vec3 temp = normalize(cross(maxVec, pl.n));

    // construct transformation to temporary system:
    tgt::mat4 m(maxVec.x, temp.x, pl.n.x, center.x,
                maxVec.y, temp.y, pl.n.y, center.y,
                maxVec.z, temp.z, pl.n.z, center.z,
                0.0f,     0.0f,   0.0f,   1.0f);
    tgt::mat4 mInv = tgt::mat4::identity;
    m.invert(mInv);

    // transform bounds to temp system in order to construct new coordinate frame
    tgt::Bounds b(volume->getLLF(), volume->getURB());
    b = b.transform(mInv*pToW);

    // construct new coordinate frame:
    tgt::vec3 origin = center;
    origin += b.getLLF().x * maxVec;
    origin += b.getLLF().y * temp;

    tgt::vec2 sp(tgt::min(volume->getSpacing()) / samplingRate);
    tgt::ivec2 res(tgt::iceil(b.diagonal().x / sp.x), tgt::iceil(b.diagonal().y / sp.y));
    int numChannels = static_cast<int>(volume->getNumChannels());

    tgt::vec3 xVec = maxVec * (sp.x * res.x);
    tgt::vec3 yVec = temp * (sp.y * res.y);

    LGL_ERROR;

    float* sliceData = new float[res.x*res.y*numChannels]; //SliceTexture gets ownership and deletes the array

    tgt::vec3 fetchX = normalize(xVec) * sp.x;
    tgt::vec3 fetchY = normalize(yVec) * sp.y;
    tgt::vec3 fetchOrigin = origin + (0.5f * fetchX) + (0.5f * fetchY);
    tgt::mat4 wToV = volume->getWorldToVoxelMatrix();
    for(int x=0; x<res.x; x++) {
        for(int y=0; y<res.y; y++) {
            tgt::vec3 pos = fetchOrigin + ((float) x * fetchX) + ((float) y * fetchY);
            pos = wToV * pos;
            for(int channel=0; channel < numChannels; channel++) {

                float valueFloat = 0.0f;
                if (hand(greaterThanEqual(pos, tgt::vec3(0.0f))) && hand(lessThanEqual(pos, dims)))
                    valueFloat = vol->getVoxelNormalizedLinear(pos, channel);

                sliceData[y * res.x * numChannels + x * numChannels + channel] = valueFloat;
            }
        }
    }

    GLint textureFormat, internalFormat;
    switch(numChannels) {
    case 1:
        textureFormat = GL_RED;
        internalFormat = GL_R32F;
        break;
    case 2:
        textureFormat = GL_RG;
        internalFormat = GL_RG32F;
        break;
    case 3:
        textureFormat = GL_RGB;
        internalFormat = GL_RGB32F;
        break;
    case 4:
        textureFormat = GL_RGBA;
        internalFormat = GL_RGBA32F;
        break;
    default:
        tgtAssert(false, "unsupported channel count");
        return nullptr;
    }

    //TODO: make dependent on input, add support for multiple channel
    SliceTexture* result = new SliceTexture(res, UNALIGNED_PLANE, volume->getFormat(), volume->getBaseType(),
                             origin, xVec, yVec, volume->getRealWorldMapping(),static_cast<void*>(sliceData), textureFormat, internalFormat, GL_FLOAT);

    return result;
}

template<typename T>
void copySliceData(const voreen::VolumeAtomic<T>* volume, T*& dataBuffer,
                   voreen::SliceAlignment sliceAlign, size_t sliceID, int* shiftArray, bool flipX, bool flipY) {

    tgtAssert(volume, "volume is null");

    tgt::ivec3 volDim = volume->getDimensions();
    size_t bytesPerChannel = volume->getBytesPerVoxel()/volume->getNumChannels();
    uint64_t zero = 0;
    if(!shiftArray) {
        switch (sliceAlign) {
        case voreen::YZ_PLANE:
            {
                tgt::ivec2 sliceDim = volDim.yz();
                tgtAssert(static_cast<int>(sliceID) < volDim.x, "invalid slice id");
                dataBuffer = new T[tgt::hmul(sliceDim)];
                for (size_t y=0; y<static_cast<size_t>(sliceDim.y); y++) {
                    for (size_t x=0; x<static_cast<size_t>(sliceDim.x); x++) {
                       dataBuffer[y*sliceDim.x + x] = volume->voxel(sliceID, !flipX ? x : (sliceDim.x-1) - x, !flipY ? y : (sliceDim.y-1) - y);
                    }
                }
            }
            break;
        case voreen::XZ_PLANE:
            {
                tgt::ivec2 sliceDim(volDim.x, volDim.z);
                tgtAssert(static_cast<int>(sliceID) < volDim.y, "invalid slice id");
                dataBuffer = new T[tgt::hmul(sliceDim)];
                for (size_t y=0; y<static_cast<size_t>(sliceDim.y); y++) {
                    for (size_t x=0; x<static_cast<size_t>(sliceDim.x); x++) {
                        dataBuffer[y*sliceDim.x + x] = volume->voxel(!flipX ? x : (sliceDim.x-1) - x, sliceID, !flipY ? y : (sliceDim.y-1) - y);
                    }
                }
            }
            break;
        case voreen::XY_PLANE:
            {
                tgt::ivec2 sliceDim(volDim.x, volDim.y);
                tgtAssert(static_cast<int>(sliceID) < volDim.z, "invalid slice id");
                dataBuffer = new T[tgt::hmul(sliceDim)];
                for (size_t y=0; y<static_cast<size_t>(sliceDim.y); y++) {
                    for (size_t x=0; x<static_cast<size_t>(sliceDim.x); x++) {
                        dataBuffer[y*sliceDim.x + x] = volume->voxel(!flipX ? x : (sliceDim.x-1) - x, !flipY ? y : (sliceDim.y-1) - y, sliceID);
                    }
                }
            }
            break;
        default:
            tgtAssert(false, "invalid slice alignment");
        }
    } else {
        for(size_t c = 0; c < volume->getNumChannels(); c++) {
            int tmpSlice = (int)sliceID + (shiftArray ? shiftArray[c] : 0);
            switch (sliceAlign) {
            case voreen::YZ_PLANE:
                {
                    tgt::ivec2 sliceDim = volDim.yz();
                    tgtAssert(static_cast<int>(sliceID) < volDim.x, "invalid slice id");
                    dataBuffer = new T[tgt::hmul(sliceDim)];
                    for (size_t y=0; y<static_cast<size_t>(sliceDim.y); y++) {
                        for (size_t x=0; x<static_cast<size_t>(sliceDim.x); x++) {
                            if(tmpSlice < 0 || tmpSlice >= volDim.x)
                                memcpy(&dataBuffer[y*sliceDim.x + x] + bytesPerChannel*c, &zero, bytesPerChannel);
                            else
                                memcpy(&dataBuffer[y*sliceDim.x + x] + bytesPerChannel*c,
                                &(volume->voxel((size_t)tmpSlice, !flipX ? x : (sliceDim.x-1) - x, !flipY ? y : (sliceDim.y-1) - y)) + bytesPerChannel*c,
                                bytesPerChannel);
                        }
                    }
                }
                break;
            case voreen::XZ_PLANE:
                {
                    tgt::ivec2 sliceDim(volDim.x, volDim.z);
                    tgtAssert(static_cast<int>(sliceID) < volDim.y, "invalid slice id");
                    dataBuffer = new T[tgt::hmul(sliceDim)];
                    for (size_t y=0; y<static_cast<size_t>(sliceDim.y); y++) {
                        for (size_t x=0; x<static_cast<size_t>(sliceDim.x); x++) {
                            if(tmpSlice < 0 || tmpSlice >= volDim.y)
                                memcpy(&dataBuffer[y*sliceDim.x + x] + bytesPerChannel*c, &zero, bytesPerChannel);
                            else
                                memcpy(&dataBuffer[y*sliceDim.x + x] + bytesPerChannel*c,
                                &(volume->voxel(!flipX ? x : (sliceDim.x-1) - x, (size_t)tmpSlice, !flipY ? y : (sliceDim.y-1) - y)) + bytesPerChannel*c,
                                bytesPerChannel);
                        }
                    }
                }
                break;
            case voreen::XY_PLANE:
                {
                    tgt::ivec2 sliceDim(volDim.x, volDim.y);
                    tgtAssert(static_cast<int>(sliceID) < volDim.z, "invalid slice id");
                    dataBuffer = new T[tgt::hmul(sliceDim)];
                    for (size_t y=0; y<static_cast<size_t>(sliceDim.y); y++) {
                        for (size_t x=0; x<static_cast<size_t>(sliceDim.x); x++) {
                            if(tmpSlice < 0 || tmpSlice >= volDim.z)
                                memcpy(&dataBuffer[y*sliceDim.x + x] + bytesPerChannel*c, &zero, bytesPerChannel);
                            else
                                memcpy(&dataBuffer[y*sliceDim.x + x] + bytesPerChannel*c,
                                &(volume->voxel(!flipX ? x : (sliceDim.x-1) - x, !flipY ? y : (sliceDim.y-1) - y, (size_t)tmpSlice)) + bytesPerChannel*c,
                                bytesPerChannel);
                        }
                    }
                }
                break;
            default:
                tgtAssert(false, "invalid slice alignment");
            }
        }
    }

}

void SliceHelper::extractAlignedSlicePixelDataHelper(const VolumeRAM* volumeRAM, SliceAlignment sliceAlign, size_t sliceID, int* shiftArray,
    void*& dataBuffer, GLint& textureFormat, GLint& internalFormat, GLenum& textureDataType) {
    tgtAssert(volumeRAM, "null pointer passed");
    tgtAssert(sliceID < volumeRAM->getDimensions()[sliceAlign], "invalid slice index");

    tgt::svec3 volDim = volumeRAM->getDimensions();

    bool flipX = false;
    bool flipY = false;
    tgt::svec3 sliceTexDim = volDim;
    if (sliceAlign == XY_PLANE) {
        sliceTexDim = tgt::svec3(volDim.xy(), 1);
    }
    else if (sliceAlign == XZ_PLANE) {
        sliceTexDim = tgt::svec3(volDim.x, volDim.z, 1);
    }
    else if (sliceAlign == YZ_PLANE) {
        sliceTexDim = tgt::svec3(volDim.y, volDim.z, 1);
    }
    else {
        tgtAssert(false, "unknown/unsupported slice alignment");
        LERROR("Unknown/unsupported slice alignment: " << sliceAlign);
        return;
    }

    //
    // allocate pixel data buffer and copy over pixel data from volume
    //

    try {
        // scalar
        if (dynamic_cast<const VolumeAtomic<uint8_t>*>(volumeRAM)) {
            textureFormat = GL_RED;
            internalFormat = GL_R8;
            textureDataType = GL_UNSIGNED_BYTE;
            copySliceData<uint8_t>(static_cast<const VolumeAtomic<uint8_t>*>(volumeRAM), reinterpret_cast<uint8_t*&>(dataBuffer),
                sliceAlign, sliceID, shiftArray, flipX, flipY);
        }
        else if (dynamic_cast<const VolumeAtomic<int8_t>*>(volumeRAM)) {
            textureFormat = GL_RED;
            internalFormat = GL_R8_SNORM;
            textureDataType = GL_BYTE;
            copySliceData<int8_t>(static_cast<const VolumeAtomic<int8_t>*>(volumeRAM), reinterpret_cast<int8_t*&>(dataBuffer),
                sliceAlign, sliceID, shiftArray, flipX, flipY);
        }
        else if (dynamic_cast<const VolumeAtomic<uint16_t>*>(volumeRAM)) {
            textureFormat = GL_RED;
            internalFormat = GL_R16;
            textureDataType = GL_UNSIGNED_SHORT;
            copySliceData<uint16_t>(static_cast<const VolumeAtomic<uint16_t>*>(volumeRAM), reinterpret_cast<uint16_t*&>(dataBuffer),
                sliceAlign, sliceID, shiftArray, flipX, flipY);
        }
        else if (dynamic_cast<const VolumeAtomic<int16_t>*>(volumeRAM)) {
            textureFormat = GL_RED;
            internalFormat = GL_R16_SNORM;
            textureDataType = GL_SHORT;
            copySliceData<int16_t>(static_cast<const VolumeAtomic<int16_t>*>(volumeRAM), reinterpret_cast<int16_t*&>(dataBuffer),
                sliceAlign, sliceID, shiftArray, flipX, flipY);
        }
        else if (dynamic_cast<const VolumeAtomic<uint32_t>*>(volumeRAM)) {
            textureFormat = GL_RED;
            internalFormat = GL_R32F;
            textureDataType = GL_UNSIGNED_INT;
            copySliceData<uint32_t>(static_cast<const VolumeAtomic<uint32_t>*>(volumeRAM), reinterpret_cast<uint32_t*&>(dataBuffer),
                sliceAlign, sliceID, shiftArray, flipX, flipY);
        }
        else if (dynamic_cast<const VolumeAtomic<int32_t>*>(volumeRAM)) {
            textureFormat = GL_RED;
            internalFormat = GL_R32F;
            textureDataType = GL_INT;
            copySliceData<int32_t>(static_cast<const VolumeAtomic<int32_t>*>(volumeRAM), reinterpret_cast<int32_t*&>(dataBuffer),
                sliceAlign, sliceID, shiftArray, flipX, flipY);
        }
        else if (dynamic_cast<const VolumeAtomic<uint64_t>*>(volumeRAM)) {
            throw VoreenException("Texture data type 'uint64' not supported by OpenGL");
        }
        else if (dynamic_cast<const VolumeAtomic<int64_t>*>(volumeRAM)) {
            throw VoreenException("Texture data type 'int64' not supported by OpenGL");
        }
        else if (dynamic_cast<const VolumeAtomic<float>*>(volumeRAM)) {
            textureFormat = GL_RED;
            internalFormat = GL_R32F;
            textureDataType = GL_FLOAT;
            copySliceData<float>(static_cast<const VolumeAtomic<float>*>(volumeRAM), reinterpret_cast<float*&>(dataBuffer),
                sliceAlign, sliceID, shiftArray, flipX, flipY);
        }
        else if (dynamic_cast<const VolumeAtomic<double>*>(volumeRAM)) {
            throw VoreenException("Texture data type 'double' not supported by OpenGL");
        }
        // vec2
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector2<uint8_t> >*>(volumeRAM)) {
            textureFormat = GL_RG;
            internalFormat = GL_RG8;
            textureDataType = GL_UNSIGNED_BYTE;
            copySliceData<tgt::Vector2<uint8_t> >(static_cast<const VolumeAtomic<tgt::Vector2<uint8_t> >*>(volumeRAM), reinterpret_cast<tgt::Vector2<uint8_t>*&>(dataBuffer),
                sliceAlign, sliceID, shiftArray, flipX, flipY);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector2<int8_t> >*>(volumeRAM)) {
            textureFormat = GL_RG;
            internalFormat = GL_RG8_SNORM;
            textureDataType = GL_BYTE;
            copySliceData<tgt::Vector2<int8_t> >(static_cast<const VolumeAtomic<tgt::Vector2<int8_t> >*>(volumeRAM), reinterpret_cast<tgt::Vector2<int8_t>*&>(dataBuffer),
                sliceAlign, sliceID, shiftArray, flipX, flipY);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector2<uint16_t> >*>(volumeRAM)) {
            textureFormat = GL_RG;
            internalFormat = GL_RG16;
            textureDataType = GL_UNSIGNED_SHORT;
            copySliceData<tgt::Vector2<uint16_t> >(static_cast<const VolumeAtomic<tgt::Vector2<uint16_t> >*>(volumeRAM), reinterpret_cast<tgt::Vector2<uint16_t>*&>(dataBuffer),
                sliceAlign, sliceID, shiftArray, flipX, flipY);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector2<int16_t> >*>(volumeRAM)) {
            textureFormat = GL_RG;
            internalFormat = GL_RG16_SNORM;
            textureDataType = GL_SHORT;
            copySliceData<tgt::Vector2<int16_t> >(static_cast<const VolumeAtomic<tgt::Vector2<int16_t> >*>(volumeRAM), reinterpret_cast<tgt::Vector2<int16_t>*&>(dataBuffer),
                sliceAlign, sliceID, shiftArray, flipX, flipY);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector2<uint32_t> >*>(volumeRAM)) {
            textureFormat = GL_RG;
            internalFormat = GL_RG;
            textureDataType = GL_UNSIGNED_INT;
            copySliceData<tgt::Vector2<uint32_t> >(static_cast<const VolumeAtomic<tgt::Vector2<uint32_t> >*>(volumeRAM), reinterpret_cast<tgt::Vector2<uint32_t>*&>(dataBuffer),
                sliceAlign, sliceID, shiftArray, flipX, flipY);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector2<int32_t> >*>(volumeRAM)) {
            textureFormat = GL_RG;
            internalFormat = GL_RG32F;
            textureDataType = GL_INT;
            copySliceData<tgt::Vector2<int32_t> >(static_cast<const VolumeAtomic<tgt::Vector2<int32_t> >*>(volumeRAM), reinterpret_cast<tgt::Vector2<int32_t>*&>(dataBuffer),
                sliceAlign, sliceID, shiftArray, flipX, flipY);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector2<uint64_t> >*>(volumeRAM)) {
            throw VoreenException("Texture data type 'uint64' not supported by OpenGL");
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector2<int64_t> >*>(volumeRAM)) {
            throw VoreenException("Texture data type 'int64' not supported by OpenGL");
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector2<float> >*>(volumeRAM)) {
            textureFormat = GL_RG;
            internalFormat = GL_RG32F;
            textureDataType = GL_FLOAT;
            copySliceData<tgt::Vector2<float> >(static_cast<const VolumeAtomic<tgt::Vector2<float> >*>(volumeRAM), reinterpret_cast<tgt::Vector2<float>*&>(dataBuffer),
                sliceAlign, sliceID, shiftArray, flipX, flipY);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector2<double> >*>(volumeRAM)) {
            throw VoreenException("Texture data type 'double' not supported by OpenGL");
        }
        // vec3
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector3<uint8_t> >*>(volumeRAM)) {
            textureFormat = GL_RGB;
            internalFormat = GL_RGB8;
            textureDataType = GL_UNSIGNED_BYTE;
            copySliceData<tgt::Vector3<uint8_t> >(static_cast<const VolumeAtomic<tgt::Vector3<uint8_t> >*>(volumeRAM), reinterpret_cast<tgt::Vector3<uint8_t>*&>(dataBuffer),
                sliceAlign, sliceID, shiftArray, flipX, flipY);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector3<int8_t> >*>(volumeRAM)) {
            textureFormat = GL_RGB;
            internalFormat = GL_RGB8_SNORM;
            textureDataType = GL_BYTE;
            copySliceData<tgt::Vector3<int8_t> >(static_cast<const VolumeAtomic<tgt::Vector3<int8_t> >*>(volumeRAM), reinterpret_cast<tgt::Vector3<int8_t>*&>(dataBuffer),
                sliceAlign, sliceID, shiftArray, flipX, flipY);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector3<uint16_t> >*>(volumeRAM)) {
            textureFormat = GL_RGB;
            internalFormat = GL_RGB16;
            textureDataType = GL_UNSIGNED_SHORT;
            copySliceData<tgt::Vector3<uint16_t> >(static_cast<const VolumeAtomic<tgt::Vector3<uint16_t> >*>(volumeRAM), reinterpret_cast<tgt::Vector3<uint16_t>*&>(dataBuffer),
                sliceAlign, sliceID, shiftArray, flipX, flipY);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector3<int16_t> >*>(volumeRAM)) {
            textureFormat = GL_RGB;
            internalFormat = GL_RGB16_SNORM;
            textureDataType = GL_SHORT;
            copySliceData<tgt::Vector3<int16_t> >(static_cast<const VolumeAtomic<tgt::Vector3<int16_t> >*>(volumeRAM), reinterpret_cast<tgt::Vector3<int16_t>*&>(dataBuffer),
                sliceAlign, sliceID, shiftArray, flipX, flipY);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector3<uint32_t> >*>(volumeRAM)) {
            textureFormat = GL_RGB;
            internalFormat = GL_RGB;
            textureDataType = GL_UNSIGNED_INT;
            copySliceData<tgt::Vector3<uint32_t> >(static_cast<const VolumeAtomic<tgt::Vector3<uint32_t> >*>(volumeRAM), reinterpret_cast<tgt::Vector3<uint32_t>*&>(dataBuffer),
                sliceAlign, sliceID, shiftArray, flipX, flipY);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector3<int32_t> >*>(volumeRAM)) {
            textureFormat = GL_RGB;
            internalFormat = GL_RGB32F;
            textureDataType = GL_INT;
            copySliceData<tgt::Vector3<int32_t> >(static_cast<const VolumeAtomic<tgt::Vector3<int32_t> >*>(volumeRAM), reinterpret_cast<tgt::Vector3<int32_t>*&>(dataBuffer),
                sliceAlign, sliceID, shiftArray, flipX, flipY);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector3<uint64_t> >*>(volumeRAM)) {
            throw VoreenException("Texture data type 'uint64' not supported by OpenGL");
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector3<int64_t> >*>(volumeRAM)) {
            throw VoreenException("Texture data type 'int64' not supported by OpenGL");
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector3<float> >*>(volumeRAM)) {
            textureFormat = GL_RGB;
            internalFormat = GL_RGB32F;
            textureDataType = GL_FLOAT;
            copySliceData<tgt::Vector3<float> >(static_cast<const VolumeAtomic<tgt::Vector3<float> >*>(volumeRAM), reinterpret_cast<tgt::Vector3<float>*&>(dataBuffer),
                sliceAlign, sliceID, shiftArray, flipX, flipY);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector3<double> >*>(volumeRAM)) {
            throw VoreenException("Texture data type 'double' not supported by OpenGL");
        }
        // vec4
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector4<uint8_t> >*>(volumeRAM)) {
            textureFormat = GL_RGBA;
            internalFormat = GL_RGBA8;
            textureDataType = GL_UNSIGNED_BYTE;
            copySliceData<tgt::Vector4<uint8_t> >(static_cast<const VolumeAtomic<tgt::Vector4<uint8_t> >*>(volumeRAM), reinterpret_cast<tgt::Vector4<uint8_t>*&>(dataBuffer),
                sliceAlign, sliceID, shiftArray, flipX, flipY);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector4<int8_t> >*>(volumeRAM)) {
            textureFormat = GL_RGBA;
            internalFormat = GL_RGBA8_SNORM;
            textureDataType = GL_BYTE;
            copySliceData<tgt::Vector4<int8_t> >(static_cast<const VolumeAtomic<tgt::Vector4<int8_t> >*>(volumeRAM), reinterpret_cast<tgt::Vector4<int8_t>*&>(dataBuffer),
                sliceAlign, sliceID, shiftArray, flipX, flipY);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector4<uint16_t> >*>(volumeRAM)) {
            textureFormat = GL_RGBA;
            internalFormat = GL_RGBA16;
            textureDataType = GL_UNSIGNED_SHORT;
            copySliceData<tgt::Vector4<uint16_t> >(static_cast<const VolumeAtomic<tgt::Vector4<uint16_t> >*>(volumeRAM), reinterpret_cast<tgt::Vector4<uint16_t>*&>(dataBuffer),
                sliceAlign, sliceID, shiftArray, flipX, flipY);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector4<int16_t> >*>(volumeRAM)) {
            textureFormat = GL_RGBA;
            internalFormat = GL_RGBA16_SNORM;
            textureDataType = GL_SHORT;
            copySliceData<tgt::Vector4<int16_t> >(static_cast<const VolumeAtomic<tgt::Vector4<int16_t> >*>(volumeRAM), reinterpret_cast<tgt::Vector4<int16_t>*&>(dataBuffer),
                sliceAlign, sliceID, shiftArray, flipX, flipY);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector4<uint32_t> >*>(volumeRAM)) {
            textureFormat = GL_RGBA;
            internalFormat = GL_RGBA;
            textureDataType = GL_UNSIGNED_INT;
            copySliceData<tgt::Vector4<uint32_t> >(static_cast<const VolumeAtomic<tgt::Vector4<uint32_t> >*>(volumeRAM), reinterpret_cast<tgt::Vector4<uint32_t>*&>(dataBuffer),
                sliceAlign, sliceID, shiftArray, flipX, flipY);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector4<int32_t> >*>(volumeRAM)) {
            textureFormat = GL_RGBA;
            internalFormat = GL_RGBA32F;
            textureDataType = GL_INT;
            copySliceData<tgt::Vector4<int32_t> >(static_cast<const VolumeAtomic<tgt::Vector4<int32_t> >*>(volumeRAM), reinterpret_cast<tgt::Vector4<int32_t>*&>(dataBuffer),
                sliceAlign, sliceID, shiftArray, flipX, flipY);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector4<uint64_t> >*>(volumeRAM)) {
            throw VoreenException("Texture data type 'uint64' not supported by OpenGL");
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector4<int64_t> >*>(volumeRAM)) {
            throw VoreenException("Texture data type 'int64' not supported by OpenGL");
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector4<float> >*>(volumeRAM)) {
            textureFormat = GL_RGBA;
            internalFormat = GL_RGBA32F;
            textureDataType = GL_FLOAT;
            copySliceData<tgt::Vector4<float> >(static_cast<const VolumeAtomic<tgt::Vector4<float> >*>(volumeRAM), reinterpret_cast<tgt::Vector4<float>*&>(dataBuffer),
                sliceAlign, sliceID, shiftArray, flipX, flipY);
        }
        else if (dynamic_cast<const VolumeAtomic<tgt::Vector4<double> >*>(volumeRAM)) {
            throw VoreenException("Texture data type 'double' not supported by OpenGL");
        }
        else {
            throw VoreenException("unknown or unsupported volume type");
        }
    }
    catch (std::exception& e) {
        throw VoreenException("Failed to extract slice pixel data: " + std::string(e.what()));
    }

}

TriangleMeshGeometryNormal* SliceHelper::getSliceGeometry(const VolumeBase* volume, SliceAlignment alignment, float sliceIndex, bool applyTransformation, const std::vector<const VolumeBase*> secondaryVolumes) {
    TriangleMeshGeometryNormal* slice = new TriangleMeshGeometryNormal();
    tgt::vec3 urb = volume->getURB();
    tgt::vec3 llf = volume->getLLF();
    tgt::vec3 sp = volume->getSpacing();
    tgt::Bounds b(llf, urb);

    tgt::mat4 wToP = volume->getWorldToPhysicalMatrix();
    for(size_t i=0; i<secondaryVolumes.size(); i++) {
        tgt::Bounds sb(secondaryVolumes[i]->getLLF(), secondaryVolumes[i]->getURB());
        tgt::Bounds sbTf = sb.transform(wToP * secondaryVolumes[i]->getPhysicalToWorldMatrix());

        b.addPoint(sbTf.getLLF());
        b.addPoint(sbTf.getURB());
    }

    tgt::vec3 bb_urb = b.getURB();
    tgt::vec3 bb_llf = b.getLLF();

    switch(alignment) {
        case YZ_PLANE: {
                           float x = sliceIndex;
                           float xcoord = llf.x + (x+0.5f) * sp.x; // We want our slice to be in the center of voxels

                           slice->addQuad(
                           VertexNormal(tgt::vec3(xcoord, bb_urb.y, bb_urb.z), tgt::vec3(xcoord, bb_urb.y, bb_urb.z)),
                           VertexNormal(tgt::vec3(xcoord, bb_urb.y, bb_llf.z), tgt::vec3(xcoord, bb_urb.y, bb_llf.z)),
                           VertexNormal(tgt::vec3(xcoord, bb_llf.y, bb_llf.z), tgt::vec3(xcoord, bb_llf.y, bb_llf.z)),
                           VertexNormal(tgt::vec3(xcoord, bb_llf.y, bb_urb.z), tgt::vec3(xcoord, bb_llf.y, bb_urb.z)));
                       }
                       break;
        case XZ_PLANE: {
                           float y = sliceIndex;
                           float ycoord = llf.y + (y+0.5f) * sp.y; // We want our slice to be in the center of voxels

                           slice->addQuad(
                           VertexNormal(tgt::vec3(bb_urb.x, ycoord, bb_urb.z), tgt::vec3(bb_urb.x, ycoord, bb_urb.z)),
                           VertexNormal(tgt::vec3(bb_urb.x, ycoord, bb_llf.z), tgt::vec3(bb_urb.x, ycoord, bb_llf.z)),
                           VertexNormal(tgt::vec3(bb_llf.x, ycoord, bb_llf.z), tgt::vec3(bb_llf.x, ycoord, bb_llf.z)),
                           VertexNormal(tgt::vec3(bb_llf.x, ycoord, bb_urb.z), tgt::vec3(bb_llf.x, ycoord, bb_urb.z)));
                       }
                       break;
        case XY_PLANE: {
                           float z = sliceIndex;
                           float zcoord = llf.z + (z+0.5f) * sp.z; // We want our slice to be in the center of voxels

                           slice->addQuad(
                           VertexNormal(tgt::vec3(bb_urb.x, bb_urb.y, zcoord), tgt::vec3(bb_urb.x, bb_urb.y, zcoord)),
                           VertexNormal(tgt::vec3(bb_urb.x, bb_llf.y, zcoord), tgt::vec3(bb_urb.x, bb_llf.y, zcoord)),
                           VertexNormal(tgt::vec3(bb_llf.x, bb_llf.y, zcoord), tgt::vec3(bb_llf.x, bb_llf.y, zcoord)),
                           VertexNormal(tgt::vec3(bb_llf.x, bb_urb.y, zcoord), tgt::vec3(bb_llf.x, bb_urb.y, zcoord)));
                       }
                       break;
        default: tgtAssert(false, "should not get here!");
    }
    tgt::mat4 m = volume->getPhysicalToTextureMatrix();

    // set coords to texture coordinates:
    for(size_t j=0; j<slice->getNumTriangles(); ++j) {
        Triangle<VertexNormal> t = slice->getTriangle(j);
        t.v_[0].pos_ = m * t.v_[0].pos_;
        t.v_[1].pos_ = m * t.v_[1].pos_;
        t.v_[2].pos_ = m * t.v_[2].pos_;
        t.v_[0].normal_ = t.v_[0].pos_;
        t.v_[1].normal_ = t.v_[1].pos_;
        t.v_[2].normal_ = t.v_[2].pos_;
        slice->setTriangle(t, j);
    }

    slice->setTransformationMatrix(volume->getTextureToWorldMatrix());
    return slice;
}

} // namespace voreen
