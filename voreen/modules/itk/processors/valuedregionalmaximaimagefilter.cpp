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

#include "valuedregionalmaximaimagefilter.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

#include "modules/itk/utils/itkwrapper.h"
#include <itkValuedRegionalMaximaImageFilter.h>

#include <iostream>

namespace voreen {

const std::string ValuedRegionalMaximaImageFilter::loggerCat_("voreen.ValuedRegionalMaximaImageFilter");

ValuedRegionalMaximaImageFilter::ValuedRegionalMaximaImageFilter()
    : ITKProcessor(),
    inport_(Port::INPORT, "input"),
    outport_(Port::OUTPORT, "output")
{
    addPort(inport_);
    addPort(outport_);
}

Processor* ValuedRegionalMaximaImageFilter::create() const {
    return new ValuedRegionalMaximaImageFilter();
}

std::string ValuedRegionalMaximaImageFilter::getProcessorInfo() const {
    return std::string("ITK ValuedRegionalMaximaImageFilter\nTransforms the image so that any pixel that is not a regional maxima is set to the minimum value for the pixel type.\nPixels that are regional maxima retain their value.\nRegional maxima are flat zones surrounded by pixels of lower value. A completely flat image will be marked as a regional maxima by this filter.");
}

template<class T>
Volume* ValuedRegionalMaximaImageFilter::applyFilter(const VolumeBase* inputVolume) {
    typedef itk::Image<T, 3> InputImageType;
    typedef itk::Image<T, 3> OutputImageType;

    typename InputImageType::Pointer p = voreenToITK<T>(inputVolume);

    typedef itk::ValuedRegionalMaximaImageFilter<InputImageType, OutputImageType>  FilterType;

    // Create the filter
    typename FilterType::Pointer filter = FilterType::New();

    filter->SetInput(p);

    observe(filter.GetPointer());

    // Run the filter
    filter->Update();

    return ITKToVoreenCopy<T>(filter->GetOutput());
}

void ValuedRegionalMaximaImageFilter::process() {
    const VolumeBase* inputVolume = inport_.getData();
    Volume* outputVolume = 0;

    if (dynamic_cast<const VolumeRAM_UInt8*>(inputVolume)) {
        outputVolume = applyFilter<uint8_t>(inputVolume);
    }
    else if (dynamic_cast<const VolumeRAM_Float*>(inputVolume)) {
        outputVolume = applyFilter<float>(inputVolume);
    }
    else {
        LERROR("Currently only VolumeRAM_UInt8 supported");
    }

    // assign computed volume to outport
    if (outputVolume) {
        transferRWM(inport_.getData(), outputVolume);
        transferTransformation(inport_.getData(), outputVolume);
        outport_.setData(outputVolume);
    }
    else
        outport_.setData(0);

}

}   // namespace
