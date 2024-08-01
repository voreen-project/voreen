/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "zerocrossingbasededgedetectionimagefilter.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include "modules/itk/utils/itkwrapper.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorconvert.h"
#include "itkImage.h"

#include "itkZeroCrossingBasedEdgeDetectionImageFilter.h"

#include <iostream>

namespace voreen {

const std::string ZeroCrossingBasedEdgeDetectionImageFilterITK::loggerCat_("voreen.ZeroCrossingBasedEdgeDetectionImageFilterITK");

ZeroCrossingBasedEdgeDetectionImageFilterITK::ZeroCrossingBasedEdgeDetectionImageFilterITK()
    : ITKProcessor(),
    inport1_(Port::INPORT, "InputImage"),
    outport1_(Port::OUTPORT, "OutputImage"),
    enableProcessing_("enabled", "Enable", false),
    maximumError_("maximumError", "MaximumError"),
    variance_("variance", "Variance"),
    backgroundValue_("backgroundValue", "BackgroundValue"),
    foregroundValue_("foregroundValue", "ForegroundValue")
{
    addPort(inport1_);
    PortConditionLogicalOr* orCondition1 = new PortConditionLogicalOr();
    orCondition1->addLinkedCondition(new PortConditionVolumeTypeFloat());
    inport1_.addCondition(orCondition1);
    addPort(outport1_);

    addProperty(enableProcessing_);
    addProperty(maximumError_);
    addProperty(variance_);
    addProperty(backgroundValue_);
    addProperty(foregroundValue_);

}

Processor* ZeroCrossingBasedEdgeDetectionImageFilterITK::create() const {
    return new ZeroCrossingBasedEdgeDetectionImageFilterITK();
}

template<class T>
void ZeroCrossingBasedEdgeDetectionImageFilterITK::zeroCrossingBasedEdgeDetectionImageFilterITK() {
    maximumError_.setVolume(inport1_.getData());
    variance_.setVolume(inport1_.getData());
    backgroundValue_.setVolume(inport1_.getData());
    foregroundValue_.setVolume(inport1_.getData());

    if (!enableProcessing_.get()) {
        outport1_.setData(inport1_.getData(), false);
        return;
    }

    typedef itk::Image<T, 3> InputImageType1;
    typedef itk::Image<T, 3> OutputImageType1;

    typename InputImageType1::Pointer p1 = voreenToITK<T>(inport1_.getData());


    //Filter define
    typedef itk::ZeroCrossingBasedEdgeDetectionImageFilter<InputImageType1, OutputImageType1> FilterType;
    typename FilterType::Pointer filter = FilterType::New();

    filter->SetInput(p1);

    filter->SetMaximumError(maximumError_.getValue<T>());
    filter->SetVariance(variance_.getValue<T>());
    filter->SetBackgroundValue(backgroundValue_.getValue<T>());
    filter->SetForegroundValue(foregroundValue_.getValue<T>());


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




void ZeroCrossingBasedEdgeDetectionImageFilterITK::process() {
    const VolumeBase* inputHandle1 = inport1_.getData();
    const VolumeRAM* inputVolume1 = inputHandle1->getRepresentation<VolumeRAM>();

    if (dynamic_cast<const VolumeRAM_Float*>(inputVolume1))  {
        zeroCrossingBasedEdgeDetectionImageFilterITK<float>();
    }
    else {
        LERROR("Inputformat of Volume 1 is not supported!");
    }

}


}   // namespace
