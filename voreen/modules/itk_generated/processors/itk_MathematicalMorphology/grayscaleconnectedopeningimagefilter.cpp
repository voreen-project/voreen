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

#include "grayscaleconnectedopeningimagefilter.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include "modules/itk/utils/itkwrapper.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorconvert.h"
#include "itkImage.h"

#include "itkGrayscaleConnectedOpeningImageFilter.h"

#include <iostream>

namespace voreen {

const std::string GrayscaleConnectedOpeningImageFilterITK::loggerCat_("voreen.GrayscaleConnectedOpeningImageFilterITK");

GrayscaleConnectedOpeningImageFilterITK::GrayscaleConnectedOpeningImageFilterITK()
    : ITKProcessor(),
    inport1_(Port::INPORT, "InputImage"),
    outport1_(Port::OUTPORT, "OutputImage"),
    seedPointPort_(Port::INPORT, "seedPointInput"),
    enableProcessing_("enabled", "Enable", false),
    fullyConnected_("fullyConnected", "FullyConnected", false),
    numSeedPoint_("numSeedPoint", "Take Seed-Point", 0, 0, 0, Processor::VALID),
    seedPoint_("SeedPoint", "Seed-Point", tgt::vec3(1),tgt::vec3(1),tgt::vec3(500))
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
    addPort(seedPointPort_);

    addProperty(enableProcessing_);
    addProperty(fullyConnected_);
    addProperty(numSeedPoint_);
    addProperty(seedPoint_);
    numSeedPoint_.onChange(MemberFunctionCallback<GrayscaleConnectedOpeningImageFilterITK>(this, &GrayscaleConnectedOpeningImageFilterITK::process));
    seedPoint_.setReadOnlyFlag(true);

}

Processor* GrayscaleConnectedOpeningImageFilterITK::create() const {
    return new GrayscaleConnectedOpeningImageFilterITK();
}

template<class T>
void GrayscaleConnectedOpeningImageFilterITK::grayscaleConnectedOpeningImageFilterITK() {

    if (!enableProcessing_.get()) {
        outport1_.setData(inport1_.getData(), false);
        return;
    }

    typedef itk::Image<T, 3> InputImageType1;
    typedef itk::Image<T, 3> OutputImageType1;

    typename InputImageType1::Pointer p1 = voreenToITK<T>(inport1_.getData());


    //Filter define
    typedef itk::GrayscaleConnectedOpeningImageFilter<InputImageType1, OutputImageType1> FilterType;
    typename FilterType::Pointer filter = FilterType::New();

    filter->SetInput(p1);

    filter->SetFullyConnected(fullyConnected_.get());

    if (seedPointPort_.hasChanged()) {
        const PointListGeometry<tgt::vec3>* pointList = dynamic_cast< const PointListGeometry<tgt::vec3>* >(seedPointPort_.getData());
        if (pointList) {
            seedPoints = pointList->getData();
            if (!seedPoints.empty()) {
                numSeedPoint_.setMinValue(1);
                numSeedPoint_.setMaxValue(seedPoints.size());
                numSeedPoint_.set(seedPoints.size());
            }
            else {
                numSeedPoint_.setMinValue(0);
                numSeedPoint_.setMaxValue(0);
                numSeedPoint_.set(0);
            }
        }
    }

    if(!seedPoints.empty()) {
        seedPoint_.set(seedPoints[numSeedPoint_.get()-1]);
    }
    else {
        seedPoint_.set(tgt::vec3 (1));
    }
    typename InputImageType1::IndexType seed;
    seed[0] = seedPoint_.get().x;
    seed[1] = seedPoint_.get().y;
    seed[2] = seedPoint_.get().z;

    filter->SetSeed(seed);


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




void GrayscaleConnectedOpeningImageFilterITK::process() {
    const VolumeBase* inputHandle1 = inport1_.getData();
    const VolumeRAM* inputVolume1 = inputHandle1->getRepresentation<VolumeRAM>();

    if (dynamic_cast<const VolumeRAM_UInt8*>(inputVolume1))  {
        grayscaleConnectedOpeningImageFilterITK<uint8_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Int8*>(inputVolume1))  {
        grayscaleConnectedOpeningImageFilterITK<int8_t>();
    }
    else if (dynamic_cast<const VolumeRAM_UInt16*>(inputVolume1))  {
        grayscaleConnectedOpeningImageFilterITK<uint16_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Int16*>(inputVolume1))  {
        grayscaleConnectedOpeningImageFilterITK<int16_t>();
    }
    else if (dynamic_cast<const VolumeRAM_UInt32*>(inputVolume1))  {
        grayscaleConnectedOpeningImageFilterITK<uint32_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Int32*>(inputVolume1))  {
        grayscaleConnectedOpeningImageFilterITK<int32_t>();
    }
    else if (dynamic_cast<const VolumeRAM_Float*>(inputVolume1))  {
        grayscaleConnectedOpeningImageFilterITK<float>();
    }
    else if (dynamic_cast<const VolumeRAM_Double*>(inputVolume1))  {
        grayscaleConnectedOpeningImageFilterITK<double>();
    }
    else {
        LERROR("Inputformat of Volume 1 is not supported!");
    }

}


}   // namespace
