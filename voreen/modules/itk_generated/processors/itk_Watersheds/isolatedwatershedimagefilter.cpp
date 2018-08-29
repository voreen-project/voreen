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

#include "isolatedwatershedimagefilter.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include "modules/itk/utils/itkwrapper.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorconvert.h"
#include "itkImage.h"

#include "itkIsolatedWatershedImageFilter.h"

#include <iostream>

namespace voreen {

const std::string IsolatedWatershedImageFilterITK::loggerCat_("voreen.IsolatedWatershedImageFilterITK");

IsolatedWatershedImageFilterITK::IsolatedWatershedImageFilterITK()
    : ITKProcessor(),
    inport1_(Port::INPORT, "InputImage"),
    outport1_(Port::OUTPORT, "OutputImage"),
    seedPointPort1_(Port::INPORT, "seedPointInput1"),
    seedPointPort2_(Port::INPORT, "seedPointInput2"),
    enableProcessing_("enabled", "Enable", false),
    numSeedPoint1_("numSeedPoint1", "Take Seed-Point", 0, 0, 0, Processor::VALID),
    seedPoint1_("SeedPoint1", "Seed-Point1", tgt::vec3(1),tgt::vec3(1),tgt::vec3(500)),
    numSeedPoint2_("numSeedPoint2", "Take Seed-Point", 0, 0, 0, Processor::VALID),
    seedPoint2_("SeedPoint2", "Seed-Point2", tgt::vec3(1),tgt::vec3(1),tgt::vec3(500)),
    threshold_("threshold", "Threshold", 0.0f, 0.0f, 1.0f),
    isolatedValueTolerance_("isolatedValueTolerance", "IsolatedValueTolerance", 0.001f, 0.0f, 1.0f),
    upperValueLimit_("upperValueLimit", "UpperValueLimit", 1.0f, 0.0f, 1.0f),
    replaceValue1_("replaceValue1", "ReplaceValue1"),
    replaceValue2_("replaceValue2", "ReplaceValue2")
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
    addPort(seedPointPort1_);
    addPort(seedPointPort2_);

    addProperty(enableProcessing_);
    addProperty(numSeedPoint1_);
    addProperty(seedPoint1_);
    numSeedPoint1_.onChange(MemberFunctionCallback<IsolatedWatershedImageFilterITK>(this, &IsolatedWatershedImageFilterITK::process));
    seedPoint1_.setReadOnlyFlag(true);
    addProperty(numSeedPoint2_);
    addProperty(seedPoint2_);
    numSeedPoint2_.onChange(MemberFunctionCallback<IsolatedWatershedImageFilterITK>(this, &IsolatedWatershedImageFilterITK::process));
    seedPoint2_.setReadOnlyFlag(true);
    addProperty(threshold_);
    addProperty(isolatedValueTolerance_);
    addProperty(upperValueLimit_);
    addProperty(replaceValue1_);
    addProperty(replaceValue2_);

}

Processor* IsolatedWatershedImageFilterITK::create() const {
    return new IsolatedWatershedImageFilterITK();
}

template<class T>
void IsolatedWatershedImageFilterITK::isolatedWatershedImageFilterITK() {
    replaceValue1_.setVolume(inport1_.getData());
    replaceValue2_.setVolume(inport1_.getData());

    if (!enableProcessing_.get()) {
        outport1_.setData(inport1_.getData(), false);
        return;
    }

    typedef itk::Image<T, 3> InputImageType1;
    typedef itk::Image<T, 3> OutputImageType1;

    typename InputImageType1::Pointer p1 = voreenToITK<T>(inport1_.getData());


    //Filter define
    typedef itk::IsolatedWatershedImageFilter<InputImageType1, OutputImageType1> FilterType;
    typename FilterType::Pointer filter = FilterType::New();

    filter->SetInput(p1);


    if (seedPointPort1_.hasChanged()) {
        const PointListGeometry<tgt::vec3>* pointList1 = dynamic_cast< const PointListGeometry<tgt::vec3>* >(seedPointPort1_.getData());
        if (pointList1) {
            seedPoints1 = pointList1->getData();
            if (!seedPoints1.empty()) {
                numSeedPoint1_.setMinValue(1);
                numSeedPoint1_.setMaxValue(seedPoints1.size());
                numSeedPoint1_.set(seedPoints1.size());
            }
            else {
                numSeedPoint1_.setMinValue(0);
                numSeedPoint1_.setMaxValue(0);
                numSeedPoint1_.set(0);
            }
        }
    }

    if(!seedPoints1.empty()) {
        seedPoint1_.set(seedPoints1[numSeedPoint1_.get()-1]);
    }
    else {
        seedPoint1_.set(tgt::vec3 (1));
    }
    typename InputImageType1::IndexType seed1;
    seed1[0] = seedPoint1_.get().x;
    seed1[1] = seedPoint1_.get().y;
    seed1[2] = seedPoint1_.get().z;

    filter->SetSeed1(seed1);

    if (seedPointPort2_.hasChanged()) {
        const PointListGeometry<tgt::vec3>* pointList2 = dynamic_cast< const PointListGeometry<tgt::vec3>* >(seedPointPort2_.getData());
        if (pointList2) {
            seedPoints2 = pointList2->getData();
            if (!seedPoints2.empty()) {
                numSeedPoint2_.setMinValue(1);
                numSeedPoint2_.setMaxValue(seedPoints2.size());
                numSeedPoint2_.set(seedPoints2.size());
            }
            else {
                numSeedPoint2_.setMinValue(0);
                numSeedPoint2_.setMaxValue(0);
                numSeedPoint2_.set(0);
            }
        }
    }

    if(!seedPoints2.empty()) {
        seedPoint2_.set(seedPoints2[numSeedPoint2_.get()-1]);
    }
    else {
        seedPoint2_.set(tgt::vec3 (1));
    }
    typename InputImageType1::IndexType seed2;
    seed2[0] = seedPoint2_.get().x;
    seed2[1] = seedPoint2_.get().y;
    seed2[2] = seedPoint2_.get().z;

    filter->SetSeed2(seed2);
    filter->SetThreshold(threshold_.get());
    filter->SetIsolatedValueTolerance(isolatedValueTolerance_.get());
    filter->SetUpperValueLimit(upperValueLimit_.get());
    filter->SetReplaceValue1(replaceValue1_.getValue<T>());
    filter->SetReplaceValue2(replaceValue2_.getValue<T>());


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




void IsolatedWatershedImageFilterITK::process() {
    const VolumeBase* inputHandle1 = inport1_.getData();
    const VolumeRAM* inputVolume1 = inputHandle1->getRepresentation<VolumeRAM>();

    if (dynamic_cast<const VolumeRAM_UInt8*>(inputVolume1))  {
        isolatedWatershedImageFilterITK<uint8_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Int8*>(inputVolume1))  {
        isolatedWatershedImageFilterITK<int8_t>();
    }
    else if (dynamic_cast<const VolumeRAM_UInt16*>(inputVolume1))  {
        isolatedWatershedImageFilterITK<uint16_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Int16*>(inputVolume1))  {
        isolatedWatershedImageFilterITK<int16_t>();
    }
    else if (dynamic_cast<const VolumeRAM_UInt32*>(inputVolume1))  {
        isolatedWatershedImageFilterITK<uint32_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Int32*>(inputVolume1))  {
        isolatedWatershedImageFilterITK<int32_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Float*>(inputVolume1))  {
        isolatedWatershedImageFilterITK<float>();
    }
    else if (dynamic_cast<const VolumeRAM_Double*>(inputVolume1))  {
        isolatedWatershedImageFilterITK<double>();
    }
    else {
        LERROR("Inputformat of Volume 1 is not supported!");
    }

}


}   // namespace
