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

#include "commands_motion.h"
#include "voreen/core/io/volumeserializer.h"
#include "voreen/core/io/volumeserializerpopulator.h"
#ifndef VRN_SNAPSHOT
#include "modules/daoc/include/vq.h"
#endif

using namespace tgt;
using namespace voreen;

void cluster(voreen::VolumeRAM* vol) {
    const int numCodewords = 9;
    std::vector<vqVector<float, 3> > vectors;
    ivec3 stepSize(2,2,2);
    //ivec3 stepSize(3,3,3);
    ivec3 dim = vol->getDimensions();
    int num = 0;
    ivec3 pos;
    for (pos.z=0; pos.z<dim.z; pos.z+=stepSize.z) {
        for (pos.y=0; pos.y<dim.y; pos.y+=stepSize.y) {
            for (pos.x=0; pos.x<dim.x; pos.x+=stepSize.x) {
                vqVector<float, 3> elem;
                elem.data_[0] = vol->getVoxelNormalized(pos, 0)-0.5f;
                elem.data_[1] = vol->getVoxelNormalized(pos, 1)-0.5f;
                elem.data_[2] = vol->getVoxelNormalized(pos, 2)-0.5f;
                if (elem.data_[0] != 0 || elem.data_[1] != 0 || elem.data_[2] != 0) {
                    vectors.push_back(elem);
                    ++num;
                }
            }
        }
    }
    LINFOC("voreen.command.cluster", "Added vectors: " << num);
    vqCodebook<float, 3> codebook(numCodewords);
    codebook.learn(vectors);
    VolumeRAM_UInt8 indexDataset(dim);
    uint8_t* indexScalars = reinterpret_cast<uint8_t*>(indexDataset.getData());
    for (pos.z=0; pos.z<dim.z; ++pos.z) {
        for (pos.y=0; pos.y<dim.y; ++pos.y) {
            for (pos.x=0; pos.x<dim.x; ++pos.x) {
                vqVector<float, 3> elem/*, substitute*/;
                elem.data_[0] = vol->getVoxelNormalized(pos, 0)-0.5f;
                elem.data_[1] = vol->getVoxelNormalized(pos, 1)-0.5f;
                elem.data_[2] = vol->getVoxelNormalized(pos, 2)-0.5f;
                //std::cout << elem.data_[0] << ", " << elem.data_[1] << ", " << elem.data_[2] << std::endl;
                int bestCodeword = codebook.getClosestCodeword(elem);
                //substitute = codebook.getCodeword(bestCodeword);
                //voxel[0] = substitute.data_[0]+128;
                //voxel[1] = substitute.data_[1]+128;
                //voxel[2] = substitute.data_[2]+128;
                //voxel[3] = 250.0f/numCodewords * bestCodeword;
                //std::cout << (int)voxel[3] << " ";
                int index = (pos.z*dim.x*dim.y + pos.y*dim.x + pos.x);
                indexScalars[index] = bestCodeword*10;
                //std::cout << bestCodeword << " ";
            }
        }
        std::cout << ".";
    }
    //FIXME
    //VolumeSerializerPopulator volLoadPop;
    //const VolumeSerializer* serializer = volLoadPop.getVolumeSerializer();
    //serializer->save("indexData.dat", &indexDataset);//FIXME

    //indexDataset.save("indexData");
    // second step: cluster position


    //create dataset for one cluster
    VolumeRAM_4xUInt8 positionDataset(dim);
    uint8_t* posScalars = reinterpret_cast<uint8_t*>(positionDataset.getData());
    for (pos.z=0; pos.z<dim.z; ++pos.z) {
        for (pos.y=0; pos.y<dim.y; ++pos.y) {
            for (pos.x=0; pos.x<dim.x; ++pos.x) {
                //Voxel<uint16_t*> voxel = vol->getVoxel(pos);
                int index = (pos.z*dim.x*dim.y + pos.y*dim.x + pos.x);
                posScalars[index*4] = 127+pos.x;
                posScalars[index*4+1] = 127+pos.y;
                posScalars[index*4+2] = 127+pos.z;
                posScalars[index*4+3] = static_cast<uint8_t>(indexScalars[index]*255.f/numCodewords);
            }
        }
    }

    //serializer->save("position.dat", &positionDataset);//FIXME
    //positionDataset.save("position");
    num = 0;
    vectors.clear();
    vqVector<float, 3> nullElem;
    nullElem.data_[0] = 127;
    nullElem.data_[1] = 127;
    nullElem.data_[2] = 127;
    vectors.push_back(nullElem);
    for (int z=0; z<dim.z; z+=stepSize.z) {
        for (int y=0; y<dim.y; y+=stepSize.y) {
            for (int x=0; x<dim.x; x+=stepSize.x) {
                tgt::col4 voxel = positionDataset.voxel(x,y,z);
                tgt::vec3 elem(voxel[0], voxel[1], voxel[2]);
                if (voxel[0] != 127 || voxel[1] != 127 || voxel[2] != 127) {
                    vqVector<float, 3> normElem;
                    normElem.data_[0] = elem.x;
                    normElem.data_[1] = elem.y;
                    normElem.data_[2] = elem.z;
                    vectors.push_back(normElem);
                    ++num;
                }
            }
        }
    }
    LINFOC("voreen.command.cluster", "Added vectors: " << num);
    vqCodebook<float, 3> codebook2(numCodewords);
    codebook2.learn(vectors);
    VolumeRAM_UInt8 ergDataset(dim);
//     uint8_t* ergScalars = (uint8_t*)ergDataset.getData();
    for (int z=0; z<dim.z; ++z) {
        for (int y=0; y<dim.y; ++y) {
            for (int x=0; x<dim.x; ++x) {
                //tgt::col4 voxel = positionDataset.voxel(x,y,z);
                //vqVector<float, 3> elem;
                //elem.data_[0] = voxel[0];
                //elem.data_[1] = voxel[1];
                //elem.data_[2] = voxel[2];
//                 int bestCodeword = codebook2.getClosestCodeword(elem);
            }
        }
        std::cout << ".";
    }
    //serializer->save("test.dat", &ergDataset);//FIXME
    //ergDataset.save("test");
}

//-------------------------------------------------------------------------------------------------

namespace voreen {

using tgt::ivec3;
using tgt::vec3;

CommandCluster::CommandCluster() /*:
    Command("--cluster") */
{
//    parameterList_ = "IN OUT";
//    loggerCat_ += "." + name_;
}

bool CommandCluster::checkParameters(const std::vector<std::string>& /*parameters*/) {
    return false;
}

bool CommandCluster::execute(const std::vector<std::string>& /*parameters*/) {
    //TODO: repair or remove this (is this used?)...jennis?
    return false;
    /*
    DatVolumeReader reader;
    VolumeDataset32Bit* dataset = (VolumeDataset32Bit*)TexMgr.load(sourceFileName, tgt::LINEAR, false, true, false);
    cluster(dataset);
    */
#ifdef VRN_HAS_MATLAB
    MatVolumeReader reader;
#endif // VRN_HAS_MATLAB
//    VolumeContainer volConX, volConY, volConZ;
//     VolumeContainer* volCon = new VolumeContainer();
#ifdef VRN_HAS_MATLAB
    reader.readIntoContainer("../../../data/dawood/Datensatz1/AltekampXTo2Vx.mat", &volConX);
    reader.readIntoContainer("../../../data/dawood/Datensatz1/AltekampXTo2Vy.mat", &volConY);
    reader.readIntoContainer("../../../data/dawood/Datensatz1/AltekampXTo2Vz.mat", &volConZ);
#endif // VRN_HAS_MATLAB
    //no such function (WHAT TO DO HERE?)
    /*VolumeDatasetDirections* vol = new VolumeDatasetDirections(volConX.get(0), volConY.get(0), volConZ.get(0));
    cluster(vol);
    */
    /*
    MatVolumeReader reader;
    VolumeContainer volContainerX;
    VolumeContainer volContainerY;
    VolumeContainer volContainerZ;
    reader.readIntoContainer(argv[2], &volContainerX);
    reader.readIntoContainer(argv[3], &volContainerY);
    reader.readIntoContainer(argv[4], &volContainerZ);

    Dataset ds;
    ds.initialize(volContainerX.get(0), volContainerY.get(0), volContainerZ.get(0));
    cluster(&ds);
    */
    //dataset->saveToDisc(targetFileName);
    return true;
}

//-----------------------------------------------------------------------------

CommandCreateMotion::CommandCreateMotion() /*:
    Command("--createmotion", "", "", "<[1|2|3|4] SIZE OUT>", 3) */
{
//    loggerCat_ += "." + name_;
}

bool CommandCreateMotion::checkParameters(const std::vector<std::string>& parameters) {
    std::set<std::string> set;
    set.insert("1");
    set.insert("2");
    set.insert("3");
    set.insert("4");
//    return ((parameters.size() == 3) && is<int>(parameters[1]) && isValueInSet(parameters[0], set));
    return false;
}

bool CommandCreateMotion::execute(const std::vector<std::string>& parameters) {
/*    int size = cast<int>(parameters[1]);
    VolumeSerializerPopulator volLoadPop;
    const VolumeSerializer* serializer = volLoadPop.getVolumeSerializer();

    if (parameters[0] == "1") {
        VolumeRAM_4xUInt8* dataset = createMotionDS(size);
        Volume h(dataset, vec3(1.0f), vec3(0.0f));
        serializer->write(parameters.back(), &h);
    }
    else if (parameters[0] == "2") {
        VolumeRAM_4xUInt8* dataset = createMotionDS2(size);
        Volume h(dataset, vec3(1.0f), vec3(0.0f));
        serializer->write(parameters.back(), &h);
    }
    else if (parameters[0] == "3") {
        VolumeRAM_4xUInt8* dataset = createMotionDS3(size);
        Volume h(dataset, vec3(1.0f), vec3(0.0f));
        serializer->write(parameters.back(), &h);
    }
    else if (parameters[0] == "4") {
        VolumeRAM_4xUInt8* dataset = createMotionDS4(size);
        Volume h(dataset, vec3(1.0f), vec3(0.0f));
        serializer->write(parameters.back(), &h);
    } */

    return true;
}

VolumeRAM_4xUInt8* CommandCreateMotion::createMotionDS(int size) {
    VolumeRAM_4xUInt8* dataset = new VolumeRAM_4xUInt8(tgt::ivec3(size));
    createMotionBox(dataset, tgt::ivec3(size/2, size/2, size/2), size, 0, 127);
    return dataset;
}

VolumeRAM_4xUInt8* CommandCreateMotion::createMotionDS2(int size) {
    VolumeRAM_4xUInt8* dataset = new VolumeRAM_4xUInt8(tgt::ivec3(size));
    int localCenter = size/4;
    createMotionBox(dataset, tgt::ivec3(1*localCenter, 1*localCenter, 1*localCenter), size/2, 0, 15);
    createMotionBox(dataset, tgt::ivec3(1*localCenter, 1*localCenter, 3*localCenter), size/2, 16, 31);
    createMotionBox(dataset, tgt::ivec3(1*localCenter, 3*localCenter, 1*localCenter), size/2, 32, 47);
    createMotionBox(dataset, tgt::ivec3(1*localCenter, 3*localCenter, 3*localCenter), size/2, 48, 63);
    createMotionBox(dataset, tgt::ivec3(3*localCenter, 1*localCenter, 1*localCenter), size/2, 64, 79);
    createMotionBox(dataset, tgt::ivec3(3*localCenter, 1*localCenter, 3*localCenter), size/2, 80, 95);
    createMotionBox(dataset, tgt::ivec3(3*localCenter, 3*localCenter, 1*localCenter), size/2, 96, 111);
    createMotionBox(dataset, tgt::ivec3(3*localCenter, 3*localCenter, 3*localCenter), size/2, 112, 127);
    return dataset;
}

VolumeRAM_4xUInt8* CommandCreateMotion::createMotionDS3(int size) {
    VolumeRAM_4xUInt8* dataset = new VolumeRAM_4xUInt8(tgt::ivec3(size));
    ivec3 pos;
    for (pos.z=0; pos.z<size; ++pos.z) {
        for (pos.y=0; pos.y<size; ++pos.y) {
            for (pos.x=0; pos.x<size; ++pos.x) {
                dataset->voxel(pos)[0] = 127;
                dataset->voxel(pos)[1] = 127;
                dataset->voxel(pos)[2] = 127;
                dataset->voxel(pos)[3] = 0;
            }
        }
    }
    int localCenter = size/4;
    createMotionBox2(dataset, ivec3(1*localCenter+0, 1*localCenter, 1*localCenter-0), size/2, size/6, ivec3(128, 0, -64), 0, 15);
    createMotionBox2(dataset, ivec3(1*localCenter+0, 1*localCenter, 3*localCenter+0), size/2, size/6, ivec3(128, 0, 64), 16, 31);
    return dataset;
}

VolumeRAM_4xUInt8* CommandCreateMotion::createMotionDS4(int size) {
    VolumeRAM_4xUInt8* dataset = new VolumeRAM_4xUInt8(tgt::ivec3(size));
    createMotionBox3(dataset, ivec3(size/2, size/2, size/2), size, size/3, 0, 127);
    return dataset;
}

void CommandCreateMotion::createMotionBox(VolumeRAM_4xUInt8* dataset, ivec3 center, int size, int minValue/*=0*/, int maxValue/*=127*/) {
    ivec3 localCenter(size/2, size/2, size/2);
    float sizeFactor = (maxValue-minValue) / (sqrt(3.0f) * (size/2.0f));
    for (int z=0; z<size; ++z) {
        for (int y=0; y<size; ++y) {
            for (int x=0; x<size; ++x) {
                ivec3 localPos(x,y,z);
                ivec3 dir = localPos-localCenter;
                ivec3 pos = center+dir;
                dir.x = static_cast<int>(dir.x * sizeFactor);
                dir.y = static_cast<int>(dir.y * sizeFactor);
                dir.z = static_cast<int>(dir.z * sizeFactor);
                dataset->voxel(pos.x, pos.y, pos.z)[0] = dir.x + 127;
                dataset->voxel(pos.x, pos.y, pos.z)[1] = dir.y + 127;
                dataset->voxel(pos.x, pos.y, pos.z)[2] = dir.z + 127;
                dataset->voxel(pos.x, pos.y, pos.z)[3] = minValue + static_cast<uint8_t>(sqrt(static_cast<float>(dir.x*dir.x + dir.y*dir.y + dir.z*dir.z)));
                //std::cout << (int)voxel[3] << " ";
            }
        }
    }
}

void CommandCreateMotion::createMotionBox2(VolumeRAM_4xUInt8* dataset, ivec3 center, int size, int radius, ivec3 direction, int minValue/*=0*/, int maxValue/*=127*/) {
    ivec3 localCenter(size/2, size/2, size/2);
    float sizeFactor = (maxValue-minValue) / (sqrt(3.0f) * (size/2.0f));
    for (int z=0; z<size; ++z) {
        for (int y=0; y<size; ++y) {
            for (int x=0; x<size; ++x) {
                ivec3 localPos(x,y,z);
                ivec3 dir = localPos-localCenter;
                ivec3 pos = center+dir;
                if (pos.x < 0 || pos.y < 0 || pos.z < 0)
                    continue;
                if (pos.x >= 64 || pos.y >= 64 || pos.z >= 64)
                    continue;
                if (length(vec3(dir)) <= radius) {
                    dataset->voxel(pos.x, pos.y, pos.z)[0] = static_cast<uint8_t>(direction.x*sizeFactor + 127);
                    dataset->voxel(pos.x, pos.y, pos.z)[1] = static_cast<uint8_t>(direction.y*sizeFactor + 127);
                    dataset->voxel(pos.x, pos.y, pos.z)[2] = static_cast<uint8_t>(direction.z*sizeFactor + 127);
                    dataset->voxel(pos.x, pos.y, pos.z)[3] = minValue + static_cast<uint8_t>(sqrt(static_cast<float>(dir.x*dir.x + dir.y*dir.y + dir.z*dir.z)));
                }
                else {
                    dataset->voxel(pos.x, pos.y, pos.z)[0] = 127;
                    dataset->voxel(pos.x, pos.y, pos.z)[1] = 127;
                    dataset->voxel(pos.x, pos.y, pos.z)[2] = 127;
                    dataset->voxel(pos.x, pos.y, pos.z)[3] = minValue;
                }
//                std::cout << (int)dataset->voxel(pos.x, pos.y, pos.z)[0] << ", "
//                  << (int)dataset->voxel(pos.x, pos.y, pos.z)[1] << ", "
//                  << (int)dataset->voxel(pos.x, pos.y, pos.z)[2] << ", "
//                  << (int)dataset->voxel(pos.x, pos.y, pos.z)[3] << std::endl;
            }
        }
    }
}

void CommandCreateMotion::createMotionBox3(VolumeRAM_4xUInt8* dataset, ivec3 center, int size, int radius, int minValue/*=0*/, int maxValue/*=127*/) {
    ivec3 localCenter(size/2, size/2, size/2);
    float sizeFactor = (maxValue-minValue) / (sqrt(3.0f) * (size/2.0f));
    for (int z=0; z<size; ++z) {
        vec3 direction;
        float angle = 2*tgt::PIf * z/size;
        direction.x = 127 * sin(angle);
        direction.y = 127 * cos(angle);
        direction.z = 0;
            for (int y=0; y<size; ++y) {
            for (int x=0; x<size; ++x) {
                ivec3 localPos(x,y,z);
                ivec3 dir = localPos-localCenter;
                ivec3 pos = center+dir;
                if (pos.x < 0 || pos.y < 0 || pos.z < 0)
                    continue;
                if (pos.x >= 64 || pos.y >= 64 || pos.z >= 64)
                    continue;
                if (length(vec3(dir)) <= radius) {
                    dataset->voxel(pos.x, pos.y, pos.z)[0] = static_cast<uint8_t>(direction.x*sizeFactor + 127);
                    dataset->voxel(pos.x, pos.y, pos.z)[1] = static_cast<uint8_t>(direction.y*sizeFactor + 127);
                    dataset->voxel(pos.x, pos.y, pos.z)[2] = static_cast<uint8_t>(direction.z*sizeFactor + 127);
                    dataset->voxel(pos.x, pos.y, pos.z)[3] = minValue + static_cast<uint8_t>(sqrt(static_cast<float>(dir.x*dir.x + dir.y*dir.y + dir.z*dir.z)));
                }
                else {
                    dataset->voxel(pos.x, pos.y, pos.z)[0] = 127;
                    dataset->voxel(pos.x, pos.y, pos.z)[1] = 127;
                    dataset->voxel(pos.x, pos.y, pos.z)[2] = 127;
                    dataset->voxel(pos.x, pos.y, pos.z)[3] = minValue;
                }
//                std::cout << (int)dataset->voxel(pos.x, pos.y, pos.z)[0] << ", "
//                  << (int)dataset->voxel(pos.x, pos.y, pos.z)[1] << ", "
//                  << (int)dataset->voxel(pos.x, pos.y, pos.z)[2] << ", "
//                  << (int)dataset->voxel(pos.x, pos.y, pos.z)[3] << std::endl;
            }
        }
    }
}

}   //namespace voreen
