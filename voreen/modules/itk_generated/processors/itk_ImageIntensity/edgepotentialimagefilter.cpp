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

#include "edgepotentialimagefilter.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include "modules/itk/utils/itkwrapper.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorconvert.h"
#include "itkImage.h"

#include "itkEdgePotentialImageFilter.h"

#include <iostream>

namespace voreen {

const std::string EdgePotentialImageFilterITK::loggerCat_("voreen.EdgePotentialImageFilterITK");

EdgePotentialImageFilterITK::EdgePotentialImageFilterITK()
    : ITKProcessor(),
    inport1_(Port::INPORT, "InputImage"),
    outport1_(Port::OUTPORT, "OutputImage"),
    enableProcessing_("enabled", "Enable", false)
{
    addPort(inport1_);
    PortConditionLogicalOr* orCondition1 = new PortConditionLogicalOr();
    orCondition1->addLinkedCondition(new PortConditionVolumeType4xUInt8());
    orCondition1->addLinkedCondition(new PortConditionVolumeType4xInt8());
    orCondition1->addLinkedCondition(new PortConditionVolumeType4xUInt16());
    orCondition1->addLinkedCondition(new PortConditionVolumeType4xInt16());
    orCondition1->addLinkedCondition(new PortConditionVolumeType3xUInt8());
    orCondition1->addLinkedCondition(new PortConditionVolumeType3xInt8());
    orCondition1->addLinkedCondition(new PortConditionVolumeType3xUInt16());
    orCondition1->addLinkedCondition(new PortConditionVolumeType3xInt16());
    orCondition1->addLinkedCondition(new PortConditionVolumeType2xUInt8());
    orCondition1->addLinkedCondition(new PortConditionVolumeType2xInt8());
    orCondition1->addLinkedCondition(new PortConditionVolumeType2xUInt16());
    orCondition1->addLinkedCondition(new PortConditionVolumeType2xInt16());
    inport1_.addCondition(orCondition1);
    addPort(outport1_);

    addProperty(enableProcessing_);

}

Processor* EdgePotentialImageFilterITK::create() const {
    return new EdgePotentialImageFilterITK();
}

template<class T>
void EdgePotentialImageFilterITK::edgePotentialImageFilterVec2ITK() {

    if (!enableProcessing_.get()) {
        outport1_.setData(inport1_.getData(), false);
        return;
    }

    typedef itk::Image<itk::CovariantVector<T,2>, 3> InputImageType1;
    typedef itk::Image<itk::CovariantVector<T,2>, 3> OutputImageType1;

    typename InputImageType1::Pointer p1 = voreenVec2ToITKVec2<T>(inport1_.getData());


    //Filter define
    typedef itk::EdgePotentialImageFilter<InputImageType1, OutputImageType1> FilterType;
    typename FilterType::Pointer filter = FilterType::New();

    filter->SetInput(p1);



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
    outputVolume1 = ITKVec2ToVoreenVec2Copy<T>(filter->GetOutput());

    if (outputVolume1) {
        transferRWM(inport1_.getData(), outputVolume1);
        transferTransformation(inport1_.getData(), outputVolume1);
        outport1_.setData(outputVolume1);
    } else
        outport1_.setData(0);



}


template<class T>
void EdgePotentialImageFilterITK::edgePotentialImageFilterVec3ITK() {

    if (!enableProcessing_.get()) {
        outport1_.setData(inport1_.getData(), false);
        return;
    }

    typedef itk::Image<itk::CovariantVector<T,3>, 3> InputImageType1;
    typedef itk::Image<itk::CovariantVector<T,3>, 3> OutputImageType1;

    typename InputImageType1::Pointer p1 = voreenVec3ToITKVec3<T>(inport1_.getData());


    //Filter define
    typedef itk::EdgePotentialImageFilter<InputImageType1, OutputImageType1> FilterType;
    typename FilterType::Pointer filter = FilterType::New();

    filter->SetInput(p1);



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
    outputVolume1 = ITKVec3ToVoreenVec3Copy<T>(filter->GetOutput());

    if (outputVolume1) {
        transferRWM(inport1_.getData(), outputVolume1);
        transferTransformation(inport1_.getData(), outputVolume1);
        outport1_.setData(outputVolume1);
    } else
        outport1_.setData(0);



}


template<class T>
void EdgePotentialImageFilterITK::edgePotentialImageFilterVec4ITK() {

    if (!enableProcessing_.get()) {
        outport1_.setData(inport1_.getData(), false);
        return;
    }

    typedef itk::Image<itk::CovariantVector<T,4>, 3> InputImageType1;
    typedef itk::Image<itk::CovariantVector<T,4>, 3> OutputImageType1;

    typename InputImageType1::Pointer p1 = voreenVec4ToITKVec4<T>(inport1_.getData());


    //Filter define
    typedef itk::EdgePotentialImageFilter<InputImageType1, OutputImageType1> FilterType;
    typename FilterType::Pointer filter = FilterType::New();

    filter->SetInput(p1);



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
    outputVolume1 = ITKVec4ToVoreenVec4Copy<T>(filter->GetOutput());

    if (outputVolume1) {
        transferRWM(inport1_.getData(), outputVolume1);
        transferTransformation(inport1_.getData(), outputVolume1);
        outport1_.setData(outputVolume1);
    } else
        outport1_.setData(0);



}




void EdgePotentialImageFilterITK::process() {
    const VolumeBase* inputHandle1 = inport1_.getData();
    const VolumeRAM* inputVolume1 = inputHandle1->getRepresentation<VolumeRAM>();

    if (dynamic_cast<const VolumeRAM_4xUInt8*>(inputVolume1))  {
        edgePotentialImageFilterVec4ITK<uint8_t>();
    }
    else if (dynamic_cast<const VolumeRAM_4xInt8*>(inputVolume1))  {
        edgePotentialImageFilterVec4ITK<int8_t>();
    }
    else if (dynamic_cast<const VolumeRAM_4xUInt16*>(inputVolume1))  {
        edgePotentialImageFilterVec4ITK<uint16_t>();
    }
    else if (dynamic_cast<const VolumeRAM_4xInt16*>(inputVolume1))  {
        edgePotentialImageFilterVec4ITK<int16_t>();
    }
    else if (dynamic_cast<const VolumeRAM_3xUInt8*>(inputVolume1))  {
        edgePotentialImageFilterVec3ITK<uint8_t>();
    }
    else if (dynamic_cast<const VolumeRAM_3xInt8*>(inputVolume1))  {
        edgePotentialImageFilterVec3ITK<int8_t>();
    }
    else if (dynamic_cast<const VolumeRAM_3xUInt16*>(inputVolume1))  {
        edgePotentialImageFilterVec3ITK<uint16_t>();
    }
    else if (dynamic_cast<const VolumeRAM_3xInt16*>(inputVolume1))  {
        edgePotentialImageFilterVec3ITK<int16_t>();
    }
    else if (dynamic_cast<const VolumeRAM_2xUInt8*>(inputVolume1))  {
        edgePotentialImageFilterVec2ITK<uint8_t>();
    }
    else if (dynamic_cast<const VolumeRAM_2xInt8*>(inputVolume1))  {
        edgePotentialImageFilterVec2ITK<int8_t>();
    }
    else if (dynamic_cast<const VolumeRAM_2xUInt16*>(inputVolume1))  {
        edgePotentialImageFilterVec2ITK<uint16_t>();
    }
    else if (dynamic_cast<const VolumeRAM_2xInt16*>(inputVolume1))  {
        edgePotentialImageFilterVec2ITK<int16_t>();
    }
    else {
        LERROR("Inputformat of Volume 1 is not supported!");
    }

}


}   // namespace
