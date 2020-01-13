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

#include "commands_modify.h"
#include "voreen/core/io/volumeserializer.h"
#include "voreen/core/io/volumeserializerpopulator.h"
#include "voreen/core/datastructures/volume/volumecollection.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumeoperator.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatormirror.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorresample.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorsubset.h"

#include "tgt/vector.h"

namespace voreen {

CommandCutToPieces::CommandCutToPieces()/* :
    Command("--cuttopieces", "", "Devide dataset into pieces with dimensions (DX, DY, DZ)",
        "<DX DY DZ IN OUT>", 5) */
{
//    loggerCat_ += "." + name_;
}

bool CommandCutToPieces::checkParameters(const std::vector<std::string>& parameters) {
//    return (parameters.size() == 5) && is<int>(parameters[0]) && is<int>(parameters[1]) && is<int>(parameters[2]);
    return false;
}

bool CommandCutToPieces::execute(const std::vector<std::string>& parameters) {
    tgt::svec3 start, dimensions;
/*    dimensions.x = cast<int>(parameters[0]);
    dimensions.y = cast<int>(parameters[1]);
    dimensions.z = cast<int>(parameters[2]);

    VolumeSerializerPopulator volLoadPop;
    const VolumeSerializer* serializer = volLoadPop.getVolumeSerializer();

    VolumeCollection* volumeCollection = serializer->read(parameters[3]);
    const VolumeBase* sourceDataset_ = volumeCollection->first();

    size_t cx, cy, cz; // number of pieces in each dimension
    cx = sourceDataset_->getDimensions().x / dimensions.x;
    cy = sourceDataset_->getDimensions().y / dimensions.y;
    cz = sourceDataset_->getDimensions().z / dimensions.z;
    if ((sourceDataset_->getDimensions().x % dimensions.x) != 0)
        ++cx;
    if ((sourceDataset_->getDimensions().y % dimensions.y) != 0)
        ++cy;
    if ((sourceDataset_->getDimensions().z % dimensions.z) != 0)
        ++cz;
    LINFO("Cutting input dataset with dimensions: " << sourceDataset_->getDimensions());
    LINFO("Into " << cx << " * " << cy << " * " << cz << " = " << (cx * cy * cz) << " pieces with dimensions: " << dimensions);

    char string[500];
    for (size_t x = 0; x < cx; x++) {
        for (size_t y = 0; y < cy; y++) {
            for (size_t z = 0; z < cz; z++) {
                tgt::svec3 dim = dimensions;
                start.x = dimensions.x * x;
                start.y = dimensions.y * y;
                start.z = dimensions.z * z;

                if (start.x + dimensions.x > sourceDataset_->getDimensions().x)
                    dim.x = (sourceDataset_->getDimensions().x - start.x);
                if (start.y + dimensions.y > sourceDataset_->getDimensions().y)
                    dim.y = (sourceDataset_->getDimensions().y - start.y);
                if (start.z + dimensions.z > sourceDataset_->getDimensions().z)
                    dim.z = (sourceDataset_->getDimensions().z - start.z);

                Volume* targetDataset_ = VolumeOperatorSubset::APPLY_OP(sourceDataset_, start, dim);
                sprintf( string, "%s-%i-%i-%i", parameters.back().c_str(), (int)x, (int)y, (int)z);

                serializer->write(string, targetDataset_);
                //targetDataset_->save(string);
                //delete targetDataset_;
            }
        }
    } */
    return true;
}

//-----------------------------------------------------------------------------

CommandScale::CommandScale() /*:
    Command("--scale", "", "Rescale Datasets\n \
                           Downsample dataset to target dimensions dx dy dz\n \
                           FILTER:\n \
                           n|nearest: nearest filtering\n \
                           l|linear: linear filtering",
                           "<FILTER DX DY DZ IN OUT>", 6) */
{
//    loggerCat_ += "." + name_;
}

bool CommandScale::checkParameters(const std::vector<std::string>& parameters) {
    /*std::set<std::string> set;
    set.insert("n");
    set.insert("nearest");
    set.insert("l");
    set.insert("linear");
    return (parameters.size() == 6) && isValueInSet(parameters[0], set) && is<int>(parameters[1]) && is<int>(parameters[2]) && is<int>(parameters[3]); */
    return false;
}

bool CommandScale::execute(const std::vector<std::string>& parameters) {
/*    VolumeRAM::Filter filter = VolumeRAM::LINEAR;

    if ((parameters[0] == "n") || (parameters[0] == "nearest")) {
        filter = VolumeRAM::NEAREST;
        LINFO("using nearest filtering");
    }
    else if ((parameters[0] == "l") || (parameters[0] == "linear")) {
        filter = VolumeRAM::LINEAR;
        LINFO("using linear filtering");
    }

    tgt::ivec3 dimensions;
    dimensions.x = cast<int>(parameters[1]);
    dimensions.y = cast<int>(parameters[2]);
    dimensions.z = cast<int>(parameters[3]);

    if ( hor(lessThan(dimensions, tgt::ivec3(1))) )
        throw tgt::Exception("target dimensions must be greater equal one");

    VolumeSerializerPopulator volLoadPop;
    const VolumeSerializer* serializer = volLoadPop.getVolumeSerializer();

    VolumeCollection* volumeCollection = serializer->read(parameters[4]);
    const VolumeBase* sourceDataset_ = volumeCollection->first();

    Volume* targetDataset_ = VolumeOperatorResample::APPLY_OP(sourceDataset_, dimensions, filter);
    serializer->write(parameters.back(), targetDataset_);
    delete sourceDataset_;
    delete targetDataset_; */
    return true;
}

//-----------------------------------------------------------------------------

CommandMirrorZ::CommandMirrorZ() /*:
    Command("--mirrorz", "", "Mirror volume on z axis",
        "<IN OUT>", 2) */
{
//    loggerCat_ += "." + name_;
}

bool CommandMirrorZ::execute(const std::vector<std::string>& parameters) {
   /* VolumeSerializerPopulator volLoadPop;
    const VolumeSerializer* serializer = volLoadPop.getVolumeSerializer();

    VolumeCollection* volumeCollection = serializer->read(parameters[0]);
    const VolumeBase* sourceDataset_ = volumeCollection->first();

    Volume* targetDataset_ = VolumeOperatorMirrorZ::APPLY_OP(sourceDataset_);

    serializer->write(parameters.back(), targetDataset_);
    delete targetDataset_; */
    return true;
}

//-----------------------------------------------------------------------------

CommandSubSet::CommandSubSet() /*:
    Command("--subset", "", "Write subset of input volume starting at (X,Y,Z) with dimensions (DX,DY,DZ)",
        "<X Y Z DX DY DZ IN OUT>", 8) */
{
//    loggerCat_ += "." + name_;
}

bool CommandSubSet::checkParameters(const std::vector<std::string>& parameters) {
    //return (parameters.size() == 8) && is<int>(parameters[0]) && is<int>(parameters[1]) && is<int>(parameters[2]) &&
    //    is<int>(parameters[3]) && is<int>(parameters[4]) && is<int>(parameters[5]);
    return false;
}

bool CommandSubSet::execute(const std::vector<std::string>& parameters) {
    /*tgt::ivec3 start, dimensions;
    start.x = cast<int>(parameters[0]);
    start.y = cast<int>(parameters[1]);
    start.z = cast<int>(parameters[2]);
    dimensions.x = cast<int>(parameters[3]);
    dimensions.y = cast<int>(parameters[4]);
    dimensions.z = cast<int>(parameters[5]);

    VolumeSerializerPopulator volLoadPop;
    const VolumeSerializer* serializer = volLoadPop.getVolumeSerializer();

    VolumeCollection* volumeCollection = serializer->read(parameters[6]);
    const VolumeBase* sourceDataset_ = volumeCollection->first();

    Volume* targetDataset_ = VolumeOperatorSubset::APPLY_OP(sourceDataset_, start, dimensions);

    serializer->write(parameters.back(), targetDataset_);
    delete targetDataset_; */
    return true;
}

//-----------------------------------------------------------------------------

CommandScaleTexCoords::CommandScaleTexCoords() /*:
    Command("--scaleTC", "", "Sample dataset to target dimensions dx dy dz\n", "<DX DY DZ IN OUT>", 5) */
{
//    loggerCat_ += "." + name_;
}

bool CommandScaleTexCoords::execute(const std::vector<std::string>& parameters) {
   /* tgt::ivec3 newDims;
    newDims.x = cast<int>(parameters[0]);
    newDims.y = cast<int>(parameters[1]);
    newDims.z = cast<int>(parameters[2]);

    if( hor(lessThan(newDims, tgt::ivec3(1))) )
        throw tgt::Exception("target dimensions must be greater equal one");

    VolumeSerializerPopulator volLoadPop;
    const VolumeSerializer* serializer = volLoadPop.getVolumeSerializer();

    VolumeCollection* volumeCollection= serializer->read(parameters[3]);
    const VolumeRAM* sourceDataset_ = volumeCollection->first()->getRepresentation<VolumeRAM>();

    const VolumeRAM_4xUInt8* source = dynamic_cast<const VolumeRAM_4xUInt8*>(sourceDataset_);

    LINFO("resampling with sampledimensions from " << source->getDimensions() << " to " << newDims);

    // build target volume
    VolumeRAM_4xUInt8* v;
    try {
         v = new VolumeRAM_4xUInt8(newDims);
    }
    catch (std::bad_alloc) {
        throw; // throw it to the caller
    }

    tgt::vec3 ratio = tgt::vec3(source->getDimensions()) / tgt::vec3(newDims);
    //tgt::vec3 invDims = 1.f / tgt::vec3(source->getDimensions());

    tgt::ivec3 pos = tgt::ivec3::zero; // iteration variable
    tgt::vec3 nearest; // knows the new position of the target volume

        for (pos.z = 0; pos.z < newDims.z; ++pos.z) {
            nearest.z = static_cast<float>(pos.z) * ratio.z;

            for (pos.y = 0; pos.y < newDims.y; ++pos.y) {
                nearest.y = static_cast<float>(pos.y) * ratio.y;

                for (pos.x = 0; pos.x < newDims.x; ++pos.x) {
                    nearest.x = static_cast<float>(pos.x) * ratio.x;
                    tgt::vec3 p = nearest - floor(nearest); // get decimal part
                    tgt::ivec3 llb = tgt::ivec3(nearest);
                    tgt::ivec3 urf = tgt::ivec3(ceil(nearest));
                    urf = tgt::min(urf, tgt::ivec3(source->getDimensions()) - 1); // clamp so the lookups do not exceed the dimensions

                    //calculate target segment:
                    double segments[256];
                    for(int k=0; k<256; ++k)
                        segments[k] = 0.0f;

                    segments[source->voxel(llb.x, llb.y, llb.z).a] += ((1.f-p.x)*(1.f-p.y)*(1.f-p.z));  // llB
                    segments[source->voxel(urf.x, llb.y, llb.z).a] += ((    p.x)*(1.f-p.y)*(1.f-p.z));  // lrB
                    segments[source->voxel(urf.x, urf.y, llb.z).a] += ((    p.x)*(    p.y)*(1.f-p.z));  // urB
                    segments[source->voxel(llb.x, urf.y, llb.z).a] += ((1.f-p.x)*(    p.y)*(1.f-p.z));  // ulB
                    segments[source->voxel(llb.x, llb.y, urf.z).a] += ((1.f-p.x)*(1.f-p.y)*(    p.z));  // llF
                    segments[source->voxel(urf.x, llb.y, urf.z).a] += ((    p.x)*(1.f-p.y)*(    p.z));  // lrF
                    segments[source->voxel(urf.x, urf.y, urf.z).a] += ((    p.x)*(    p.y)*(    p.z));  // urF
                    segments[source->voxel(llb.x, urf.y, urf.z).a] += ((1.f-p.x)*(    p.y)*(    p.z));  // ulF

                    if(segments[0] < 1.0f) {
                        int seg = 1;
                        for(int k=0; k<256; ++k)
                            if(segments[k] > segments[seg])
                                seg = k;

                        tgt::dvec3 res = tgt::dvec3(.0f,.0f,.0f);
                        if(source->voxel(llb.x, llb.y, llb.z).a == seg) {
                            tgt::dvec3 value((double)source->voxel(llb.x, llb.y, llb.z).x, (double)source->voxel(llb.x, llb.y, llb.z).y, (double)source->voxel(llb.x, llb.y, llb.z).z);
                            res += value * (double)((1.f-p.x)*(1.f-p.y)*(1.f-p.z));  // llB
                        }
                        if(source->voxel(urf.x, llb.y, llb.z).a == seg) {
                            tgt::dvec3 value((double)source->voxel(urf.x, llb.y, llb.z).x, (double)source->voxel(urf.x, llb.y, llb.z).y, (double)source->voxel(urf.x, llb.y, llb.z).z);
                            res += value * (double)((    p.x)*(1.f-p.y)*(1.f-p.z));  // lrB
                        }
                        if(source->voxel(urf.x, urf.y, llb.z).a == seg) {
                            tgt::dvec3 value((double)source->voxel(urf.x, urf.y, llb.z).x, (double)source->voxel(urf.x, urf.y, llb.z).y, (double)source->voxel(urf.x, urf.y, llb.z).z);
                            res += value * double((    p.x)*(    p.y)*(1.f-p.z));  // urB
                        }
                        if(source->voxel(llb.x, urf.y, llb.z).a == seg) {
                            tgt::dvec3 value((double)source->voxel(llb.x, urf.y, llb.z).x, (double)source->voxel(llb.x, urf.y, llb.z).y, (double)source->voxel(llb.x, urf.y, llb.z).z);
                            res += value * double((1.f-p.x)*(    p.y)*(1.f-p.z));  // ulB
                        }
                        if(source->voxel(llb.x, llb.y, urf.z).a == seg) {
                            tgt::dvec3 value((double)source->voxel(llb.x, llb.y, urf.z).x, (double)source->voxel(llb.x, llb.y, urf.z).y, (double)source->voxel(llb.x, llb.y, urf.z).z);
                            res += value * double((1.f-p.x)*(1.f-p.y)*(    p.z));  // llF
                        }
                        if(source->voxel(urf.x, llb.y, urf.z).a == seg) {
                            tgt::dvec3 value((double)source->voxel(urf.x, llb.y, urf.z).x, (double)source->voxel(urf.x, llb.y, urf.z).y, (double)source->voxel(urf.x, llb.y, urf.z).z);
                            res += value * double((    p.x)*(1.f-p.y)*(    p.z));  // lrF
                        }
                        if(source->voxel(urf.x, urf.y, urf.z).a == seg) {
                            tgt::dvec3 value((double)source->voxel(urf.x, urf.y, urf.z).x, (double)source->voxel(urf.x, urf.y, urf.z).y, (double)source->voxel(urf.x, urf.y, urf.z).z);
                            res += value * double((    p.x)*(    p.y)*(    p.z));  // urF
                        }
                        if(source->voxel(llb.x, urf.y, urf.z).a == seg) {
                            tgt::dvec3 value((double)source->voxel(llb.x, urf.y, urf.z).x, (double)source->voxel(llb.x, urf.y, urf.z).y, (double)source->voxel(llb.x, urf.y, urf.z).z);
                            res += value * double((1.f-p.x)*(    p.y)*(    p.z));// ulF
                        }

                        res /= segments[seg];
                        v->voxel(pos).x = static_cast<uint8_t>(res.x);
                        v->voxel(pos).y = static_cast<uint8_t>(res.y);
                        v->voxel(pos).z = static_cast<uint8_t>(res.z);
                        v->voxel(pos).w = static_cast<uint8_t>(seg);
                    }
                    else {
                        v->voxel(pos).x = 0;
                        v->voxel(pos).y = 0;
                        v->voxel(pos).z = 0;
                        v->voxel(pos).w = 0;
                    }
                }
            }
        }
    Volume h(v, volumeCollection->first());
    serializer->write(parameters.back(), &h);
    delete sourceDataset_;
    //delete v; */
    return true;
}

}   //namespace voreen
