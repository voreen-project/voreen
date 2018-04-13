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

#include "hausdorffdistanceimagefilter.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include "modules/itk/utils/itkwrapper.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorconvert.h"
#include "itkImage.h"

#include "itkHausdorffDistanceImageFilter.h"

#include <iostream>

namespace voreen {

const std::string HausdorffDistanceImageFilterITK::loggerCat_("voreen.HausdorffDistanceImageFilterITK");

HausdorffDistanceImageFilterITK::HausdorffDistanceImageFilterITK()
    : ITKProcessor(),
    inport1_(Port::INPORT, "InputImage1"),
    inport2_(Port::INPORT, "InputImage2"),
    hausdorffDistance_("hausdorffDistance", "HausdorffDistance", 0.f, -FLT_MAX, FLT_MAX),
    averageHausdorffDistance_("averageHausdorffDistance", "AverageHausdorffDistance", 0.f, -FLT_MAX, FLT_MAX)
{
    addPort(inport1_);
    PortConditionLogicalOr* orCondition1 = new PortConditionLogicalOr();
    orCondition1->addLinkedCondition(new PortConditionVolumeTypeUInt8());
    orCondition1->addLinkedCondition(new PortConditionVolumeTypeInt8());
    orCondition1->addLinkedCondition(new PortConditionVolumeTypeUInt16());
    orCondition1->addLinkedCondition(new PortConditionVolumeTypeInt16());
    orCondition1->addLinkedCondition(new PortConditionVolumeTypeUInt32());
    orCondition1->addLinkedCondition(new PortConditionVolumeTypeInt32());
    orCondition1->addLinkedCondition(new PortConditionVolumeTypeFloat());
    orCondition1->addLinkedCondition(new PortConditionVolumeTypeDouble());
    inport1_.addCondition(orCondition1);
    addPort(inport2_);
    PortConditionLogicalOr* orCondition2 = new PortConditionLogicalOr();
    orCondition2->addLinkedCondition(new PortConditionVolumeTypeUInt8());
    orCondition2->addLinkedCondition(new PortConditionVolumeTypeInt8());
    orCondition2->addLinkedCondition(new PortConditionVolumeTypeUInt16());
    orCondition2->addLinkedCondition(new PortConditionVolumeTypeInt16());
    orCondition2->addLinkedCondition(new PortConditionVolumeTypeUInt32());
    orCondition2->addLinkedCondition(new PortConditionVolumeTypeInt32());
    orCondition2->addLinkedCondition(new PortConditionVolumeTypeFloat());
    orCondition2->addLinkedCondition(new PortConditionVolumeTypeDouble());
    inport2_.addCondition(orCondition2);

    addProperty(hausdorffDistance_);
    hausdorffDistance_.setReadOnlyFlag(true);
    addProperty(averageHausdorffDistance_);
    averageHausdorffDistance_.setReadOnlyFlag(true);

}

Processor* HausdorffDistanceImageFilterITK::create() const {
    return new HausdorffDistanceImageFilterITK();
}

template<class T, class S>
void HausdorffDistanceImageFilterITK::hausdorffDistanceImageFilterITK() {


    typedef itk::Image<T, 3> InputImageType1;
    typedef itk::Image<S, 3> InputImageType2;

    typename InputImageType1::Pointer p1 = voreenToITK<T>(inport1_.getData());
    typename InputImageType2::Pointer p2 = voreenToITK<S>(inport2_.getData());


    //Filter define
    typedef itk::HausdorffDistanceImageFilter<InputImageType1, InputImageType2> FilterType;
    typename FilterType::Pointer filter = FilterType::New();

    filter->SetInput1(p1);
    filter->SetInput2(p2);



    observe(filter.GetPointer());

    try
    {
        filter->Update();
        hausdorffDistance_.set(filter->GetHausdorffDistance());
        averageHausdorffDistance_.set(filter->GetAverageHausdorffDistance());

    }
    catch (itk::ExceptionObject &e)
    {
        LERROR(e);
    }




}




void HausdorffDistanceImageFilterITK::process() {
    const VolumeBase* inputHandle1 = inport1_.getData();
    const VolumeRAM* inputVolume1 = inputHandle1->getRepresentation<VolumeRAM>();

    if (dynamic_cast<const VolumeRAM_UInt8*>(inputVolume1))  {
        volumeTypeSwitch1<uint8_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Int8*>(inputVolume1))  {
        volumeTypeSwitch1<int8_t>();
    }
    else if (dynamic_cast<const VolumeRAM_UInt16*>(inputVolume1))  {
        volumeTypeSwitch1<uint16_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Int16*>(inputVolume1))  {
        volumeTypeSwitch1<int16_t>();
    }
    else if (dynamic_cast<const VolumeRAM_UInt32*>(inputVolume1))  {
        volumeTypeSwitch1<uint32_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Int32*>(inputVolume1))  {
        volumeTypeSwitch1<int32_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Float*>(inputVolume1))  {
        volumeTypeSwitch1<float>();
    }
    else if (dynamic_cast<const VolumeRAM_Double*>(inputVolume1))  {
        volumeTypeSwitch1<double>();
    }
    else {
        LERROR("Inputformat of Volume 1 is not supported!");
    }

}

template <class T>
void HausdorffDistanceImageFilterITK::volumeTypeSwitch1() {
    const VolumeBase* inputHandle2 = inport2_.getData();
    const VolumeRAM* inputVolume2 = inputHandle2->getRepresentation<VolumeRAM>();

    if (dynamic_cast<const VolumeRAM_UInt8*>(inputVolume2))  {
        hausdorffDistanceImageFilterITK<T, uint8_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Int8*>(inputVolume2))  {
        hausdorffDistanceImageFilterITK<T, int8_t>();
    }
    else if (dynamic_cast<const VolumeRAM_UInt16*>(inputVolume2))  {
        hausdorffDistanceImageFilterITK<T, uint16_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Int16*>(inputVolume2))  {
        hausdorffDistanceImageFilterITK<T, int16_t>();
    }
    else if (dynamic_cast<const VolumeRAM_UInt32*>(inputVolume2))  {
        hausdorffDistanceImageFilterITK<T, uint32_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Int32*>(inputVolume2))  {
        hausdorffDistanceImageFilterITK<T, int32_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Float*>(inputVolume2))  {
        hausdorffDistanceImageFilterITK<T, float>();
    }
    else if (dynamic_cast<const VolumeRAM_Double*>(inputVolume2))  {
        hausdorffDistanceImageFilterITK<T, double>();
    }
    else {
        LERROR("Inputformat of Volume 2 is not supported!");
    }

}


}   // namespace
