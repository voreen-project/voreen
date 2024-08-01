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

#include "binaryerodeimagefilter.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include "modules/itk/utils/itkwrapper.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorconvert.h"
#include "itkImage.h"

#include "itkBinaryErodeImageFilter.h"
#include "itkBinaryBallStructuringElement.h"
#include "itkBinaryCrossStructuringElement.h"
#include "itkFlatStructuringElement.h"

#include <iostream>

namespace voreen {

const std::string BinaryErodeImageFilterITK::loggerCat_("voreen.BinaryErodeImageFilterITK");

BinaryErodeImageFilterITK::BinaryErodeImageFilterITK()
    : ITKProcessor(),
    inport1_(Port::INPORT, "InputImage"),
    outport1_(Port::OUTPORT, "OutputImage"),
    enableProcessing_("enabled", "Enable", false),
    erodeValue_("erodeValue", "ErodeValue"),
    structuringElement_("structuringElement", "Structuring Element"),
    shape_("shape", "Shape"),
    radius_("Radius", "Radius", tgt::ivec3(1),tgt::ivec3(1),tgt::ivec3(20))
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
    addPort(outport1_);


    structuringElement_.addOption("binaryBall", "BinaryBall");
    structuringElement_.addOption("binaryCross", "BinaryCross");
    structuringElement_.addOption("flat", "Flat");
    addProperty(structuringElement_);

    shape_.addOption("box","Box");
    shape_.addOption("ball","Ball");
    shape_.addOption("cross","Cross");
    shape_.addOption("annulus","Annulus");
    addProperty(shape_);
    shape_.setVisibleFlag(false);

    addProperty(radius_);
    addProperty(enableProcessing_);
    addProperty(erodeValue_);

}

Processor* BinaryErodeImageFilterITK::create() const {
    return new BinaryErodeImageFilterITK();
}

template<class T>
void BinaryErodeImageFilterITK::binaryErodeImageFilterITK() {
    erodeValue_.setVolume(inport1_.getData());

    if (!enableProcessing_.get()) {
        outport1_.setData(inport1_.getData(), false);
        return;
    }

    typedef itk::Image<T, 3> InputImageType1;
    typedef itk::Image<T, 3> OutputImageType1;


    typedef  T  PixelType;

    typename InputImageType1::Pointer p1 = voreenToITK<T>(inport1_.getData());


    if(structuringElement_.get() == "binaryBall"){
        shape_.setVisibleFlag(false);
        typedef itk::BinaryBallStructuringElement < PixelType, 3 > KernelType;
        typedef itk::BinaryErodeImageFilter<InputImageType1, OutputImageType1, KernelType> FilterType;
        typename FilterType::Pointer filter = FilterType::New();

    filter->SetInput(p1);

    filter->SetErodeValue(erodeValue_.getValue<T>());

        KernelType structuringElement;
        typename KernelType::SizeType radius;
        radius[0] = radius_.get().x;
        radius[1] = radius_.get().y;
        radius[2] = radius_.get().z;
        structuringElement.SetRadius(radius);
        structuringElement.CreateStructuringElement();
        filter->SetKernel(structuringElement);

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

    else if(structuringElement_.get() == "binaryCross"){
        shape_.setVisibleFlag(false);
        typedef itk::BinaryCrossStructuringElement < PixelType, 3 > KernelType;
        typedef itk::BinaryErodeImageFilter<InputImageType1, OutputImageType1, KernelType> FilterType;
        typename FilterType::Pointer filter = FilterType::New();

    filter->SetInput(p1);

    filter->SetErodeValue(erodeValue_.getValue<T>());

        KernelType structuringElement;
        typename KernelType::SizeType radius;
        radius[0] = radius_.get().x;
        radius[1] = radius_.get().y;
        radius[2] = radius_.get().z;
        structuringElement.SetRadius(radius);
        structuringElement.CreateStructuringElement();
        filter->SetKernel(structuringElement);

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

    else if(structuringElement_.get() == "flat"){
        shape_.setVisibleFlag(true);
        typedef itk::FlatStructuringElement < 3 > KernelType;
        typename KernelType::SizeType radius;
        radius[0] = radius_.get().x;
        radius[1] = radius_.get().y;
        radius[2] = radius_.get().z;

        if(shape_.get() == "box"){
            KernelType structuringElement = KernelType::Box(radius);
            typedef itk::BinaryErodeImageFilter<InputImageType1, OutputImageType1, KernelType> FilterType;
            typename FilterType::Pointer filter = FilterType::New();

    filter->SetInput(p1);

    filter->SetErodeValue(erodeValue_.getValue<T>());

            filter->SetKernel(structuringElement);

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

        else if(shape_.get() == "ball"){
            KernelType structuringElement = KernelType::Ball(radius);
            typedef itk::BinaryErodeImageFilter<InputImageType1, OutputImageType1, KernelType> FilterType;
            typename FilterType::Pointer filter = FilterType::New();

    filter->SetInput(p1);

    filter->SetErodeValue(erodeValue_.getValue<T>());

            filter->SetKernel(structuringElement);

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

        else if(shape_.get() == "cross"){
            KernelType structuringElement = KernelType::Cross(radius);
            typedef itk::BinaryErodeImageFilter<InputImageType1, OutputImageType1, KernelType> FilterType;
            typename FilterType::Pointer filter = FilterType::New();

    filter->SetInput(p1);

    filter->SetErodeValue(erodeValue_.getValue<T>());

            filter->SetKernel(structuringElement);


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
        else if(shape_.get() == "annulus"){
            KernelType structuringElement = KernelType::Annulus(radius);
            typedef itk::BinaryErodeImageFilter<InputImageType1, OutputImageType1, KernelType> FilterType;
            typename FilterType::Pointer filter = FilterType::New();

    filter->SetInput(p1);

    filter->SetErodeValue(erodeValue_.getValue<T>());

            filter->SetKernel(structuringElement);


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
    }
}




void BinaryErodeImageFilterITK::process() {
    const VolumeBase* inputHandle1 = inport1_.getData();
    const VolumeRAM* inputVolume1 = inputHandle1->getRepresentation<VolumeRAM>();

    if (dynamic_cast<const VolumeRAM_UInt8*>(inputVolume1))  {
        binaryErodeImageFilterITK<uint8_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Int8*>(inputVolume1))  {
        binaryErodeImageFilterITK<int8_t>();
    }
    else if (dynamic_cast<const VolumeRAM_UInt16*>(inputVolume1))  {
        binaryErodeImageFilterITK<uint16_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Int16*>(inputVolume1))  {
        binaryErodeImageFilterITK<int16_t>();
    }
    else if (dynamic_cast<const VolumeRAM_UInt32*>(inputVolume1))  {
        binaryErodeImageFilterITK<uint32_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Int32*>(inputVolume1))  {
        binaryErodeImageFilterITK<int32_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Float*>(inputVolume1))  {
        binaryErodeImageFilterITK<float>();
    }
    else if (dynamic_cast<const VolumeRAM_Double*>(inputVolume1))  {
        binaryErodeImageFilterITK<double>();
    }
    else {
        LERROR("Inputformat of Volume 1 is not supported!");
    }

}


}   // namespace
