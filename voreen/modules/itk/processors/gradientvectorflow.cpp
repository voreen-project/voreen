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

#include "gradientvectorflow.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

#include "modules/itk/utils/itkwrapper.h"
#include <itkGradientRecursiveGaussianImageFilter.h>
#include <itkImageRegionIteratorWithIndex.h>
#include <itkCovariantVector.h>
#include <itkGradientImageFilter.h>
#include <itkGradientToMagnitudeImageFilter.h>
#include <itkDerivativeImageFilter.h>
#include <itkGradientVectorFlowImageFilter.h>
#include <itkLaplacianImageFilter.h>

#include "itkHessian3DToVesselnessMeasureImageFilter.h"

#include <iostream>

namespace voreen {

const std::string GradientVectorFlow::loggerCat_("voreen.GradientVectorFlow");

GradientVectorFlow::GradientVectorFlow()
    : ITKProcessor(),
    inport_(Port::INPORT, "input"),
    outport_(Port::OUTPORT, "output"),
    noiseLevel_("noiseLevel", "Noise Level", 500.0f, 0.0f, 10000.0f),
    timeStep_("timeStep", "Time Step", 0.001f, 0.0f, 1.0f),
    numSteps_("numSteps", "Num Steps", 20, 1, 1000),
    normalize_("normalize", "Normalize", false)
{
    addPort(inport_);
    addPort(outport_);

    addProperty(noiseLevel_);
    addProperty(timeStep_);
    addProperty(numSteps_);
    addProperty(normalize_);
}

Processor* GradientVectorFlow::create() const {
    return new GradientVectorFlow();
}

std::string GradientVectorFlow::getProcessorInfo() const {
    return std::string("ITK GradientVectorFlow filter");
}

template<class T>
Volume* GradientVectorFlow::gradientvectorflow(const VolumeBase* inputVolume) {
    typedef itk::Image<T, 3> InputImageType;
    typename InputImageType::Pointer p = voreenToITK<T>(inputVolume);

    // Declare gradient type
    typedef itk::CovariantVector<float, 3> myGradientType;
    typedef itk::Image<myGradientType, 3>   myGradientImageType;

    //-----------------------------------------------------
    //Calculate Gradient:
    typedef itk::GradientImageFilter<InputImageType, float, float> myGFilterType;
    typename myGFilterType::Pointer gfilter = myGFilterType::New();

    // Connect the input images
    gfilter->SetInput(p);

    //-----------------------------------------------------
    //Calculate Gradient Vector Flow:

    typedef itk::LaplacianImageFilter<InputImageType, InputImageType> myLaplacianFilterType;
    typename myLaplacianFilterType::Pointer m_LFilter = myLaplacianFilterType::New();

    typedef itk::GradientVectorFlowImageFilter<myGradientImageType, myGradientImageType, float> myGVFFilterType;
    myGVFFilterType::Pointer m_GVFFilter = myGVFFilterType::New();

    m_GVFFilter->SetInput(gfilter->GetOutput());
    m_GVFFilter->SetLaplacianFilter(m_LFilter);
    m_GVFFilter->SetNoiseLevel(noiseLevel_.get());
    m_GVFFilter->SetTimeStep(timeStep_.get());
    m_GVFFilter->SetIterationNum(numSteps_.get());

    // Run the filter
    observe(m_GVFFilter.GetPointer());
    m_GVFFilter->Update();

    //return ITKToVoreenCopy<float>(gtomfilter->GetOutput());

    //Convert to voreen vector volume:
    typename myGradientImageType::Pointer vol = m_GVFFilter->GetOutput();

    typename myGradientImageType::RegionType region;
    region = vol->GetBufferedRegion();

    typename myGradientImageType::SizeType size = region.GetSize();
//    typename myGradientImageType::IndexType start = region.GetIndex();
    tgt::ivec3 dim(size[0], size[1], size[2]);

    typename myGradientImageType::SpacingType sp = vol->GetSpacing();
    tgt::vec3 spacing(sp[0], sp[1], sp[2]);

    typename itk::ImageRegionConstIterator< myGradientImageType > it(vol, region);
    it.GoToBegin();

    VolumeAtomic<tgt::vec3>* out = new VolumeAtomic<tgt::vec3>(dim);
    tgt::vec3* data = (tgt::vec3*)out->getData();
    //VolumeAtomic<tgt::col3>* out = new VolumeAtomic<tgt::col3>(dim);
    //tgt::col3* data = (tgt::col3*)out->getData();

    while(!it.IsAtEnd())
    {
        tgt::vec3 v;
        v.x = it.Get()[0];
        v.y = it.Get()[1];
        v.z = it.Get()[2];
        if(normalize_.get()) {
            if(length(v) > 0.00001f)
                v = normalize(v);
        }

        //v *= 3.0f;
        //v.x = std::min(1.0f, v.x);
        //v.y = std::min(1.0f, v.y);
        //v.z = std::min(1.0f, v.z);
        //tgt::col3 c(v.x * 255.0f, v.y * 255.0f, v.z * 255.0f);
        //*data = c;

        *data = v;
        ++it;
        ++data;
    }
    return new Volume(out, inputVolume);

    // Declare the type for the vesselness filter
    //typedef itk::Hessian3DToVesselnessMeasureImageFilter<float>  myVesselnessFilterType;
    //typedef myVesselnessFilterType::OutputImageType myVesselnessImageType;
    //myVesselnessFilterType::Pointer filterVesselness = myVesselnessFilterType::New();

    //typename itk::Image<itk::SymmetricSecondRankTensor<double, 3>, 3>::Pointer hessianImage;

    //// Connect the input images
    //filterVesselness->SetInput( hessianImage );

    //// Execute the filter
    //filterVesselness->Update();

    //return ITKToVoreenCopy<float>(filterVesselness->GetOutput());
}

void GradientVectorFlow::process() {
    const VolumeBase* inputVolume = inport_.getData();
    Volume* outputVolume = 0;

    if (dynamic_cast<const VolumeRAM_Float*>(inputVolume)) {
        outputVolume = gradientvectorflow<float>(inputVolume);
    }
    else {
        LERROR("Currently only VolumeRAM_Float supported");
    }

    // assign computed volume to outport
    if (outputVolume)
        outport_.setData(outputVolume);
    else
        outport_.setData(0);

}

}   // namespace
