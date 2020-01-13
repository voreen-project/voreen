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

#include "volumefilter_itk.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

#include "modules/itk/utils/itkwrapper.h"

#include <itkMedianImageFilter.h>
#include <itkThresholdImageFilter.h>

#include <iostream>

namespace voreen {

const std::string VolumeFilterITK::loggerCat_("voreen.VolumeFilterITK");

VolumeFilterITK::VolumeFilterITK()
    : VolumeProcessor(),
    inport_(Port::INPORT, "input"),
    outport_(Port::OUTPORT, "output")
{
    addPort(inport_);
    addPort(outport_);
}

Processor* VolumeFilterITK::create() const {
    return new VolumeFilterITK();
}

std::string VolumeFilterITK::getProcessorInfo() const {
    return std::string("ITK test processor (median filter)");
}

void VolumeFilterITK::process() {
    const VolumeRAM* inputVolume = inport_.getData()->getRepresentation<VolumeRAM>();
    Volume* outputVolume = 0;

    if (dynamic_cast<const VolumeRAM_UInt8*>(inputVolume)) {
        typedef itk::Image<uint8_t, 3> ImageType;

        ImageType::Pointer p = voreenToITK<uint8_t>(inport_.getData());

        typedef itk::MedianImageFilter<ImageType, ImageType> FilterType;
        FilterType::Pointer filter = FilterType::New();

        ImageType::SizeType indexRadius;
        indexRadius.Fill(1);
        filter->SetRadius(indexRadius);

        filter->SetInput(p);
        filter->Update();
        outputVolume = ITKToVoreenCopy<uint8_t>(filter->GetOutput());
    }
    else {
        LERROR("Currently only VolumeRAM_UInt8 supported");
    }

    // assign computed volume to outport
    if (outputVolume)
        outport_.setData(outputVolume);
    else
        outport_.setData(0);

}

}   // namespace
