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

#include "narrowbandthresholdsegmentationlevelsetimagefilter.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include "modules/itk/utils/itkwrapper.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorconvert.h"
#include "itkImage.h"

#include "itkNarrowBandThresholdSegmentationLevelSetImageFilter.h"

#include <iostream>

namespace voreen {

const std::string NarrowBandThresholdSegmentationLevelSetImageFilterITK::loggerCat_("voreen.NarrowBandThresholdSegmentationLevelSetImageFilterITK");

NarrowBandThresholdSegmentationLevelSetImageFilterITK::NarrowBandThresholdSegmentationLevelSetImageFilterITK()
    : ITKProcessor(),
    inport1_(Port::INPORT, "InputImage"),
    inport2_(Port::INPORT, "FeatureImage"),
    outport1_(Port::OUTPORT, "OutputImage"),
    propagationScaling_("propagationScaling", "PropagationScaling", 2.0f, 0.0f, 100.0f),
    curvatureScaling_("curvatureScaling", "CurvatureScaling", 1.0f, 0.0f, 100.0f),
    advectionScaling_("advectionScaling", "AdvectionScaling", 1.0f, 0.0f, 100.0f),
    maximumRMSError_("maximumRMSError", "MaximumRMSError", 0.02f, 0.0f, 1.0f),
    numberOfIterations_("numberOfIterations", "NumberOfIterations", 800, 1, 10000)
{
    addPort(inport1_);
    PortConditionLogicalOr* orCondition1 = new PortConditionLogicalOr();
    orCondition1->addLinkedCondition(new PortConditionVolumeTypeFloat());
    inport1_.addCondition(orCondition1);
    addPort(inport2_);
    PortConditionLogicalOr* orCondition2 = new PortConditionLogicalOr();
    orCondition2->addLinkedCondition(new PortConditionVolumeTypeFloat());
    inport2_.addCondition(orCondition2);
    addPort(outport1_);

    addProperty(propagationScaling_);
    addProperty(curvatureScaling_);
    addProperty(advectionScaling_);
    addProperty(maximumRMSError_);
    addProperty(numberOfIterations_);

}

Processor* NarrowBandThresholdSegmentationLevelSetImageFilterITK::create() const {
    return new NarrowBandThresholdSegmentationLevelSetImageFilterITK();
}

template<class T, class S>
void NarrowBandThresholdSegmentationLevelSetImageFilterITK::narrowBandThresholdSegmentationLevelSetImageFilterITK() {


    typedef itk::Image<T, 3> InputImageType1;
    typedef itk::Image<S, 3> InputImageType2;
    typedef itk::Image<T, 3> OutputImageType1;

    typename InputImageType1::Pointer p1 = voreenToITK<T>(inport1_.getData());
    typename InputImageType2::Pointer p2 = voreenToITK<S>(inport2_.getData());


    //Filter define
    typedef itk::NarrowBandThresholdSegmentationLevelSetImageFilter<InputImageType1, InputImageType2> FilterType;
    typename FilterType::Pointer filter = FilterType::New();

    filter->SetInput(p1);
    filter->SetFeatureImage(p2);

    filter->SetPropagationScaling(propagationScaling_.get());
    filter->SetCurvatureScaling(curvatureScaling_.get());
    filter->SetAdvectionScaling(advectionScaling_.get());
    filter->SetMaximumRMSError(maximumRMSError_.get());
    filter->SetNumberOfIterations(numberOfIterations_.get());


    observe(filter.GetPointer());

    try
    {
        filter->Update();

    }
    catch (itk::ExceptionObject &e)
    {
        LERROR(e);
    }


    Volume* outputVolume1 = 0;
    outputVolume1 = ITKToVoreenCopy<T>(filter->GetOutput());

    if (outputVolume1) {
        transferRWM(inport1_.getData(), outputVolume1);
        transferTransformation(inport1_.getData(), outputVolume1);
        outport1_.setData(outputVolume1);
    } else
        outport1_.setData(0);



}




void NarrowBandThresholdSegmentationLevelSetImageFilterITK::process() {
    const VolumeBase* inputHandle1 = inport1_.getData();
    const VolumeRAM* inputVolume1 = inputHandle1->getRepresentation<VolumeRAM>();

    if (dynamic_cast<const VolumeRAM_Float*>(inputVolume1))  {
        volumeTypeSwitch1<float>();
    }
    else {
        LERROR("Inputformat of Volume 1 is not supported!");
    }

}

template <class T>
void NarrowBandThresholdSegmentationLevelSetImageFilterITK::volumeTypeSwitch1() {
    const VolumeBase* inputHandle2 = inport2_.getData();
    const VolumeRAM* inputVolume2 = inputHandle2->getRepresentation<VolumeRAM>();

    if (dynamic_cast<const VolumeRAM_Float*>(inputVolume2))  {
        narrowBandThresholdSegmentationLevelSetImageFilterITK<T, float>();
    }
    else {
        LERROR("Inputformat of Volume 2 is not supported!");
    }

}


}   // namespace
