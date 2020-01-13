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

#include "vesselness.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

#include "voreen/core/datastructures/volume/operators/volumeoperatorconvert.h"
#include "modules/itk/utils/itkwrapper.h"

#include "itkHessian3DToVesselnessMeasureImageFilter.h"
#include "itkHessianRecursiveGaussianImageFilter.h"
#include "itkSymmetricSecondRankTensor.h"
#include "itkImage.h"

#include "itkHessian3DToVesselnessMeasureImageFilter.h"
#include "itkMultiScaleHessianBasedMeasureImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"

#include <iostream>

namespace voreen {

const std::string Vesselness::loggerCat_("voreen.Vesselness");

Vesselness::Vesselness()
    : ITKProcessor(),
    inport_(Port::INPORT, "input"),
    outport_(Port::OUTPORT, "output"),
    outportScale_(Port::OUTPORT, "outputScale"),
    volumeFloat_(0),
    sigmaMin_("sigmaMin", "Sigma min", 0.5f, 0.0f, 10.0f),
    sigmaMax_("sigmaMax", "Sigma max", 4.0f, 0.0f, 50.0f),
    sigmaSteps_("sigmaSteps", "Sigma steps", 5, 1, 100),
    alpha1_("alpha1", "Alpha 1", 0.1f,  0.0f, 10.0f),
    alpha2_("alpha2", "Alpha 2", 2.0f, 0.0f, 10.0f)
{
    addPort(inport_);
    addPort(outport_);
    addPort(outportScale_);

    addProperty(sigmaMin_);
    addProperty(sigmaMax_);
    addProperty(sigmaSteps_);
    addProperty(alpha1_);
    addProperty(alpha2_);
}

Processor* Vesselness::create() const {
    return new Vesselness();
}

std::string Vesselness::getProcessorInfo() const {
    return std::string("ITK test processor (median filter)");
}

bool Vesselness::isReady() const {
    if(inport_.isReady() && outport_.isReady())
        return true;
    else
        return false;
}

void Vesselness::process() {
    const VolumeRAM* inputVolume = inport_.getData()->getRepresentation<VolumeRAM>();
    Volume* outputVolume = 0;
    Volume* outputVolumeScale = 0;

    if(inport_.hasChanged()) {
        delete volumeFloat_;
        volumeFloat_ = 0;
    }

    if(!volumeFloat_) {
        if(dynamic_cast<const VolumeRAM_Float*>(inputVolume)) {
            volumeFloat_ = inport_.getData();
        }
        else {
            //volumeFloat_ = new VolumeRAM_Float(inputVolume->getDimensions(), inputVolume->getSpacing(), inputVolume->getTransformation());
            VolumeOperatorConvert voConvert;
            volumeFloat_ = voConvert.apply<float>(inport_.getData());
        }
    }

    // Define the dimension of the images
    typedef float PixelType;
    // Declare the types of the images
    typedef itk::Image<PixelType, 3> ImageType;
    typedef double InputPixelType;
    typedef itk::NumericTraits< InputPixelType >::RealType RealPixelType;
    typedef itk::SymmetricSecondRankTensor< RealPixelType, 3> HessianPixelType;
    typedef itk::Image< HessianPixelType, 3> HessianImageType;

    ImageType::Pointer p = voreenToITK<float>(volumeFloat_);

    typedef itk::RescaleIntensityImageFilter<ImageType> RescaleFilterType;
    // Declare the type of enhancement filter - use ITK’s 3D vesselness (Sato)
    typedef itk::Hessian3DToVesselnessMeasureImageFilter<float> VesselnessFilterType;
    // Declare the type of multiscale enhancement filter
    typedef itk::MultiScaleHessianBasedMeasureImageFilter<ImageType, HessianImageType> MultiScaleEnhancementFilterType;

    // Instantiate the multiscale filter and set the input image
    MultiScaleEnhancementFilterType::Pointer multiScaleEnhancementFilter = MultiScaleEnhancementFilterType::New();
    multiScaleEnhancementFilter->SetInput(p);
    multiScaleEnhancementFilter->SetSigmaMinimum(sigmaMin_.get());
    multiScaleEnhancementFilter->SetSigmaMaximum(sigmaMax_.get());
    multiScaleEnhancementFilter->SetNumberOfSigmaSteps(sigmaSteps_.get());

    // Get the vesselness filter and set the parameters
    VesselnessFilterType::Pointer vesselnessFilter = VesselnessFilterType::New();
    vesselnessFilter->SetAlpha1(alpha1_.get());
    vesselnessFilter->SetAlpha2(alpha2_.get());
    multiScaleEnhancementFilter->SetHessianToMeasureFilter(vesselnessFilter);

    multiScaleEnhancementFilter->SetGenerateScalesOutput(outportScale_.isReady());

    observe(multiScaleEnhancementFilter.GetPointer());
    // Now run the multiscale filter
    try
    {
        multiScaleEnhancementFilter->Update();
    }
    catch (itk::ExceptionObject &e)
    {
        LERROR(e);
    }
    outputVolume = ITKToVoreenCopy<float>(multiScaleEnhancementFilter->GetOutput());

    //if(outportScale_.isReady())
    outputVolumeScale = ITKToVoreenCopy<float>(multiScaleEnhancementFilter->GetScalesOutput());

    // assign computed volume to outport
    if (outputVolume) {
        transferTransformation(inport_.getData(), outputVolume);
        outport_.setData(outputVolume);
    }
    else
        outport_.setData(0);

    // assign computed volume to outport
    if (outputVolumeScale) {
        transferTransformation(inport_.getData(), outputVolumeScale);
        outportScale_.setData(outputVolumeScale);
    }
    else
        outportScale_.setData(0);

    if(dynamic_cast<const VolumeRAM_Float*>(inputVolume))
        volumeFloat_ = 0;
}

}   // namespace
