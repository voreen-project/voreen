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

#include "radarglyphvolumegenerator.h"
#include "voreen/core/datastructures/callback/memberfunctioncallback.h"
#include "voreen/core/datastructures/volume/volumelist.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

namespace voreen {

RadarGlyphVolumeGenerator::RadarGlyphVolumeGenerator()
    : Processor()
    , outport_(Port::OUTPORT,"volumecollection.outport")
    , optionProp_("optionProp","Data:")
    , generateProp_("generateProp","Generate")
    , generateNew_(true)
{
    addPort(outport_);
    addProperty(optionProp_);
        optionProp_.addOption("2_3_4_20","2_3_4_20",GO_2_3_4_20);
        optionProp_.addOption("50_50_50_50","50_50_50_50",GO_50_50_50_50);
        optionProp_.addOption("50_50_50_500","50_50_50_500",GO_50_50_50_500);
    addProperty(generateProp_);
    generateProp_.onChange(MemberFunctionCallback<RadarGlyphVolumeGenerator>(this, &RadarGlyphVolumeGenerator::generateOnChange));
}

RadarGlyphVolumeGenerator::~RadarGlyphVolumeGenerator() {
}

void RadarGlyphVolumeGenerator::generateOnChange(){
    generateNew_ = true;
    invalidate();
}

void RadarGlyphVolumeGenerator::process() {
    if(generateNew_){
        VolumeList* vc = new VolumeList();
        switch(optionProp_.getValue()) {
        case GO_2_3_4_20:
            //generate new volumes
            for(int i = 0; i < 20; ++i){
                setProgress(static_cast<float>(i)/19.f);
                VolumeRAM_3xFloat* v3f = new VolumeRAM_3xFloat(tgt::ivec3(2,3,4));
                for(size_t j = 1; j <= v3f->getNumVoxels(); ++j){
                    double value = (((double)(i)/20.0)*6.2831853);
                    /*if((i == 8) || (i == 12))
                        v3f->voxel(j-1) = tgt::vec3(0.f);
                    else*/
                        v3f->voxel(j-1) = tgt::vec3((float)std::cos(value),(float)std::sin(value),(float)i/20.f);
                }
                Volume* vh = new Volume(v3f,tgt::vec3(1.f),tgt::vec3(0.f));
                oldVolumePosition(vh);
                vc->add(vh);
            }
            break;
        case GO_50_50_50_50:
            //generate new volumes
            for(int i = 0; i < 50; ++i){
                setProgress((float)i/49.f);
                VolumeRAM_3xFloat* v3f = new VolumeRAM_3xFloat(tgt::ivec3(50,50,50));
                for(size_t j = 1; j <= v3f->getNumVoxels(); ++j){
                    double value = (((double)(i)/20.0)*6.2831853);
                    //if((i % 10) == 0)
                    //    v3f->voxel(j-1) = tgt::vec3(0.f);
                    //else
                        v3f->voxel(j-1) = tgt::vec3((float)std::cos(value),(float)std::sin(value),(float)i/50.f);
                }
                Volume* vh = new Volume(v3f,tgt::vec3(1.f),tgt::vec3(0.f));
                oldVolumePosition(vh);
                vc->add(vh);
            }
            break;
        case GO_50_50_50_500:
            //generate new volumes
            for(int i = 0; i < 500; ++i){
                setProgress((float)i/499.f);
                VolumeRAM_3xFloat* v3f = new VolumeRAM_3xFloat(tgt::ivec3(50,50,50));
                for(size_t j = 1; j <= v3f->getNumVoxels(); ++j){
                    double value = (((double)(i)/20.0)*6.2831853);
                    //if((i % 10) != 0)
                        v3f->voxel(j-1) = tgt::vec3((float)std::cos(value),(float)std::sin(value),(float)i/500.f);
                    //else
                    //    v3f->voxel(j-1) = tgt::vec3(0.f);
                }
                Volume* vh = new Volume(v3f,tgt::vec3(1.f),tgt::vec3(0.f));
                oldVolumePosition(vh);
                vc->add(vh);
            }
            break;
        default:
            tgtAssert(false,"Should not get here!");
            break;
        }
        outport_.setData(vc);
        //reset generateNew
        generateNew_ = false;
    }
}


}   // namespace
