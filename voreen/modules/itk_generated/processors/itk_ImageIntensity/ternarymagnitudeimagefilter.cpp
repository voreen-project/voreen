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

#include "ternarymagnitudeimagefilter.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include "modules/itk/utils/itkwrapper.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorconvert.h"
#include "itkImage.h"

#include "itkTernaryMagnitudeImageFilter.h"

#include <iostream>

namespace voreen {

const std::string TernaryMagnitudeImageFilterITK::loggerCat_("voreen.TernaryMagnitudeImageFilterITK");

TernaryMagnitudeImageFilterITK::TernaryMagnitudeImageFilterITK()
    : ITKProcessor(),
    inport1_(Port::INPORT, "InputImage1"),
    inport2_(Port::INPORT, "InputImage2"),
    inport3_(Port::INPORT, "InputImage3"),
    outport1_(Port::OUTPORT, "OutputImage")
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
    addPort(inport3_);
    PortConditionLogicalOr* orCondition3 = new PortConditionLogicalOr();
    orCondition3->addLinkedCondition(new PortConditionVolumeTypeUInt8());
    orCondition3->addLinkedCondition(new PortConditionVolumeTypeInt8());
    orCondition3->addLinkedCondition(new PortConditionVolumeTypeUInt16());
    orCondition3->addLinkedCondition(new PortConditionVolumeTypeInt16());
    orCondition3->addLinkedCondition(new PortConditionVolumeTypeUInt32());
    orCondition3->addLinkedCondition(new PortConditionVolumeTypeInt32());
    orCondition3->addLinkedCondition(new PortConditionVolumeTypeFloat());
    orCondition3->addLinkedCondition(new PortConditionVolumeTypeDouble());
    inport3_.addCondition(orCondition3);
    addPort(outport1_);


}

Processor* TernaryMagnitudeImageFilterITK::create() const {
    return new TernaryMagnitudeImageFilterITK();
}

template<class T, class S, class R>
void TernaryMagnitudeImageFilterITK::ternaryMagnitudeImageFilterITK() {


    typedef itk::Image<T, 3> InputImageType1;
    typedef itk::Image<S, 3> InputImageType2;
    typedef itk::Image<R, 3> InputImageType3;
    typedef itk::Image<T, 3> OutputImageType1;

    typename InputImageType1::Pointer p1 = voreenToITK<T>(inport1_.getData());
    typename InputImageType2::Pointer p2 = voreenToITK<S>(inport2_.getData());
    typename InputImageType3::Pointer p3 = voreenToITK<R>(inport3_.getData());


    //Filter define
    typedef itk::TernaryMagnitudeImageFilter<InputImageType1, InputImageType2, InputImageType3, OutputImageType1> FilterType;
    typename FilterType::Pointer filter = FilterType::New();

    filter->SetInput1(p1);
    filter->SetInput2(p2);
    filter->SetInput3(p3);



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




void TernaryMagnitudeImageFilterITK::process() {
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
void TernaryMagnitudeImageFilterITK::volumeTypeSwitch1() {
    const VolumeBase* inputHandle2 = inport2_.getData();
    const VolumeRAM* inputVolume2 = inputHandle2->getRepresentation<VolumeRAM>();

    if (dynamic_cast<const VolumeRAM_UInt8*>(inputVolume2))  {
        volumeTypeSwitch2<T, uint8_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Int8*>(inputVolume2))  {
        volumeTypeSwitch2<T, int8_t>();
    }
    else if (dynamic_cast<const VolumeRAM_UInt16*>(inputVolume2))  {
        volumeTypeSwitch2<T, uint16_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Int16*>(inputVolume2))  {
        volumeTypeSwitch2<T, int16_t>();
    }
    else if (dynamic_cast<const VolumeRAM_UInt32*>(inputVolume2))  {
        volumeTypeSwitch2<T, uint32_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Int32*>(inputVolume2))  {
        volumeTypeSwitch2<T, int32_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Float*>(inputVolume2))  {
        volumeTypeSwitch2<T, float>();
    }
    else if (dynamic_cast<const VolumeRAM_Double*>(inputVolume2))  {
        volumeTypeSwitch2<T, double>();
    }
    else {
        LERROR("Inputformat of Volume 2 is not supported!");
    }

}

template <class T, class S>
void TernaryMagnitudeImageFilterITK::volumeTypeSwitch2() {
    const VolumeBase* inputHandle3 = inport3_.getData();
    const VolumeRAM* inputVolume3 = inputHandle3->getRepresentation<VolumeRAM>();

    if (dynamic_cast<const VolumeRAM_UInt8*>(inputVolume3))  {
        ternaryMagnitudeImageFilterITK<T, S, uint8_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Int8*>(inputVolume3))  {
        ternaryMagnitudeImageFilterITK<T, S, int8_t>();
    }
    else if (dynamic_cast<const VolumeRAM_UInt16*>(inputVolume3))  {
        ternaryMagnitudeImageFilterITK<T, S, uint16_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Int16*>(inputVolume3))  {
        ternaryMagnitudeImageFilterITK<T, S, int16_t>();
    }
    else if (dynamic_cast<const VolumeRAM_UInt32*>(inputVolume3))  {
        ternaryMagnitudeImageFilterITK<T, S, uint32_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Int32*>(inputVolume3))  {
        ternaryMagnitudeImageFilterITK<T, S, int32_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Float*>(inputVolume3))  {
        ternaryMagnitudeImageFilterITK<T, S, float>();
    }
    else if (dynamic_cast<const VolumeRAM_Double*>(inputVolume3))  {
        ternaryMagnitudeImageFilterITK<T, S, double>();
    }
    else {
        LERROR("Inputformat of Volume 3 is not supported!");
    }

}


}   // namespace
