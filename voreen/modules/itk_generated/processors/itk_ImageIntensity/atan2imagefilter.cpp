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

#include "atan2imagefilter.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include "modules/itk/utils/itkwrapper.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorconvert.h"
#include "itkImage.h"

#include "itkAtan2ImageFilter.h"

#include <iostream>

namespace voreen {

const std::string Atan2ImageFilterITK::loggerCat_("voreen.Atan2ImageFilterITK");

Atan2ImageFilterITK::Atan2ImageFilterITK()
    : ITKProcessor(),
    inport1_(Port::INPORT, "InputImage1"),
    inport2_(Port::INPORT, "InputImage2"),
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
    addPort(outport1_);


}

Processor* Atan2ImageFilterITK::create() const {
    return new Atan2ImageFilterITK();
}

template<class T, class S>
void Atan2ImageFilterITK::atan2ImageFilterITK() {


    typedef itk::Image<T, 3> InputImageType1;
    typedef itk::Image<S, 3> InputImageType2;
    typedef itk::Image<T, 3> OutputImageType1;

    typename InputImageType1::Pointer p1 = voreenToITK<T>(inport1_.getData());
    typename InputImageType2::Pointer p2 = voreenToITK<S>(inport2_.getData());


    //Filter define
    typedef itk::Atan2ImageFilter<InputImageType1, InputImageType2, OutputImageType1> FilterType;
    typename FilterType::Pointer filter = FilterType::New();

    filter->SetInput1(p1);
    filter->SetInput2(p2);



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




void Atan2ImageFilterITK::process() {
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
void Atan2ImageFilterITK::volumeTypeSwitch1() {
    const VolumeBase* inputHandle2 = inport2_.getData();
    const VolumeRAM* inputVolume2 = inputHandle2->getRepresentation<VolumeRAM>();

    if (dynamic_cast<const VolumeRAM_UInt8*>(inputVolume2))  {
        atan2ImageFilterITK<T, uint8_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Int8*>(inputVolume2))  {
        atan2ImageFilterITK<T, int8_t>();
    }
    else if (dynamic_cast<const VolumeRAM_UInt16*>(inputVolume2))  {
        atan2ImageFilterITK<T, uint16_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Int16*>(inputVolume2))  {
        atan2ImageFilterITK<T, int16_t>();
    }
    else if (dynamic_cast<const VolumeRAM_UInt32*>(inputVolume2))  {
        atan2ImageFilterITK<T, uint32_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Int32*>(inputVolume2))  {
        atan2ImageFilterITK<T, int32_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Float*>(inputVolume2))  {
        atan2ImageFilterITK<T, float>();
    }
    else if (dynamic_cast<const VolumeRAM_Double*>(inputVolume2))  {
        atan2ImageFilterITK<T, double>();
    }
    else {
        LERROR("Inputformat of Volume 2 is not supported!");
    }

}


}   // namespace
