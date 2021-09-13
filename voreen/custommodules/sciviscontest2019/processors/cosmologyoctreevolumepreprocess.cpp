/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2015 University of Muenster, Germany.                        *
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
#include "cosmologyoctreevolumepreprocess.h"
#include "../utils/cmmath.h"

#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumegl.h"
#include "voreen/core/datastructures/volume/volumeminmaxmagnitude.h"


namespace voreen{
CosmologyOctreeVolumePreprocess::CosmologyOctreeVolumePreprocess()
    : CachingVolumeProcessor()
    , outport_         (Port::OUTPORT, "volumehandle.output",   "Volume Output", false)
    , inport_          (Port::INPORT,  "volumehandle.input", "Particle Data Output")
    , volume_(nullptr)
{
    addPort(outport_);
    addPort(inport_);
    
}


void CosmologyOctreeVolumePreprocess::process(){
    if (!inport_.getData()) return;

    delete volume_;

    tgt::svec3 dim    = inport_.getData()->getDimensions();
    tgt::vec3 spacing = inport_.getData()->getSpacing();
    tgt::vec3 offset  = inport_.getData()->getOffset();
	
	VolumeRAM_UInt16* volumeData = new VolumeRAM_UInt16(dim, true);

    const VolumeRAM * volumeIn = inport_.getData()->getRepresentation<VolumeRAM>();
    
    const float* in = static_cast<const float*>(volumeIn->getData());
    uint16_t *out = volumeData->voxel();

	VolumeMinMaxMagnitude* vmm = inport_.getData()->getDerivedData<VolumeMinMaxMagnitude>();
	float min = vmm->getMinMagnitude();
	float max = vmm->getMaxMagnitude();


    size_t size = tgt::hmul(dim);
    for(size_t i = 0; i != size; i++){
        float f = *in++;
        f = (f-min)/(max-min);
        uint16_t v = static_cast<uint16_t>(65535*f);
        *out++=v;
    }

    volume_ = new Volume(volumeData, spacing, offset);
    outport_.setData(volume_);
}
}
