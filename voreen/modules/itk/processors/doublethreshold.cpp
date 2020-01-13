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

#include "doublethreshold.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

#include "modules/itk/utils/itkwrapper.h"
#include <itkDoubleThresholdImageFilter.h>

#include <iostream>

namespace voreen {

const std::string DoubleThreshold::loggerCat_("voreen.DoubleThreshold");

DoubleThreshold::DoubleThreshold()
    : ITKProcessor(),
    inport_(Port::INPORT, "input"),
    outport_(Port::OUTPORT, "output"),
    threshold1_("threshold1", "Threshold 1", 0.0f, 0.0f, 65535.0f),
    threshold2_("threshold2", "Threshold 2", 20.0f, 0.0f, 65535.0f),
    threshold3_("threshold3", "Threshold 3", 200.0f, 0.0f, 65535.0f),
    threshold4_("threshold4", "Threshold 4", 255.0f, 0.0f, 65535.0f)
{
    addPort(inport_);
    addPort(outport_);

    addProperty(threshold1_);
    addProperty(threshold2_);
    addProperty(threshold3_);
    addProperty(threshold4_);
}

Processor* DoubleThreshold::create() const {
    return new DoubleThreshold();
}

std::string DoubleThreshold::getProcessorInfo() const {
    return std::string("ITK DoubleThreshold filter");
}

template<class T>
Volume* DoubleThreshold::doubleThreshold(const VolumeBase* inputVolume, T t1, T t2, T t3, T t4) {
    typedef itk::Image<T, 3> InputImageType;
    typedef itk::Image<uint8_t, 3> OutputImageType;

    typename InputImageType::Pointer p = voreenToITK<T>(inputVolume);

    // define the fillhole filter
    typedef itk::DoubleThresholdImageFilter<InputImageType, OutputImageType>  DoubleThresholdFilterType;

    // Create the filter
    typename DoubleThresholdFilterType::Pointer threshold = DoubleThresholdFilterType::New();

    // Setup the fillhole method
    threshold->SetInput(p);
    threshold->SetInsideValue(255);
    threshold->SetOutsideValue(0);
    threshold->SetThreshold1(t1);
    threshold->SetThreshold2(t2);
    threshold->SetThreshold3(t3);
    threshold->SetThreshold4(t4);

    observe(threshold.GetPointer());
    // Run the filter
    threshold->Update();

    return ITKToVoreenCopy<uint8_t>(threshold->GetOutput());
}

void DoubleThreshold::process() {
    const VolumeBase* inputVolume = inport_.getData();
    Volume* outputVolume = 0;

    if (dynamic_cast<const VolumeRAM_UInt8*>(inputVolume)) {
        outputVolume = doubleThreshold<uint8_t>(inputVolume, threshold1_.get(), threshold2_.get(), threshold3_.get(), threshold4_.get());
    }
    else if (dynamic_cast<const VolumeRAM_Float*>(inputVolume)) {
        outputVolume = doubleThreshold<float>(inputVolume, threshold1_.get()/65535.0f, threshold2_.get()/65535.0f, threshold3_.get()/65535.0f, threshold4_.get()/65535.0f);
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
