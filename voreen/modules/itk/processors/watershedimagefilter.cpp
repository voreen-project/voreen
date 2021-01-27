/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "watershedimagefilter.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include "modules/itk/utils/itkwrapper.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorconvert.h"
#include "itkImage.h"

#include "itkWatershedImageFilter.h"

#include <iostream>

namespace voreen {

const std::string WatershedImageFilterITK::loggerCat_("voreen.WatershedImageFilterITK");

WatershedImageFilterITK::WatershedImageFilterITK()
    : ITKProcessor(),
    inport1_(Port::INPORT, "input1"),
    outport1_(Port::OUTPORT, "output1"),
    enableProcessing_("enabled", "Enable", false),
    level_("level", "Level", 1.0f, 0.0f, 5000.0f),
    threshold_("threshold", "Threshold", 1.0f, 0.0f, 5000.0f)
{
    addPort(inport1_);
    PortConditionLogicalOr* orCondition1 = new PortConditionLogicalOr();
    orCondition1->addLinkedCondition(new PortConditionVolumeTypeFloat());
    inport1_.addCondition(orCondition1);
    addPort(outport1_);

    addProperty(enableProcessing_);
    addProperty(level_);
    addProperty(threshold_);

}

Processor* WatershedImageFilterITK::create() const {
    return new WatershedImageFilterITK();
}

template<class T>
void WatershedImageFilterITK::watershedImageFilterITK() {

    if (!enableProcessing_.get()) {
        outport1_.setData(inport1_.getData(), false);
        return;
    }

    typedef itk::Image<T, 3> InputImageType1;
    typedef itk::Image<VolumeRAM_UInt64, 3> OutputImageType1;


    typename InputImageType1::Pointer p1 = voreenToITK<T>(inport1_.getData());


    //Filter define
    typedef itk::WatershedImageFilter<InputImageType1> FilterType;
    typename FilterType::Pointer filter = FilterType::New();

    filter->SetInput(p1);

    filter->SetLevel(level_.get());
    filter->SetThreshold(threshold_.get());


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
    outputVolume1 = ITKToVoreenCopy<itk::IdentifierType>(filter->GetOutput());

    if (outputVolume1) {
        transferTransformation(inport1_.getData(), outputVolume1);
        outport1_.setData(outputVolume1);
    }
    else
        outport1_.setData(0);

}


void WatershedImageFilterITK::process() {
    const VolumeBase* inputHandle1 = inport1_.getData();
    const VolumeRAM* inputVolume1 = inputHandle1->getRepresentation<VolumeRAM>();

    if (dynamic_cast<const VolumeRAM_Float*>(inputVolume1))  {
        watershedImageFilterITK<float>();
    }
     else {
        LERROR("Inputformat of Volume 1 is not supported!");
    }

}


}   // namespace
