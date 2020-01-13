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

#include "vectorcurvatureanisotropicdiffusionimagefilter.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include "modules/itk/utils/itkwrapper.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorconvert.h"
#include "itkImage.h"

#include "itkVectorCurvatureAnisotropicDiffusionImageFilter.h"

#include <iostream>

namespace voreen {

const std::string VectorCurvatureAnisotropicDiffusionImageFilterITK::loggerCat_("voreen.VectorCurvatureAnisotropicDiffusionImageFilterITK");

VectorCurvatureAnisotropicDiffusionImageFilterITK::VectorCurvatureAnisotropicDiffusionImageFilterITK()
    : ITKProcessor(),
    inport1_(Port::INPORT, "InputImage"),
    outport1_(Port::OUTPORT, "OutputImage"),
    enableProcessing_("enabled", "Enable", false),
    numberOfIterations_("numberOfIterations", "NumberOfIterations", 1, 1, 1000),
    timeStep_("timeStep", "TimeStep", 0.0625f, 0.0f, 2.0f),
    conductanceParameter_("conductanceParameter", "ConductanceParameter", 0.1f, 0.0f, 10.0f)
{
    addPort(inport1_);
    PortConditionLogicalOr* orCondition1 = new PortConditionLogicalOr();
    orCondition1->addLinkedCondition(new PortConditionVolumeType4xFloat());
    orCondition1->addLinkedCondition(new PortConditionVolumeType4xDouble());
    orCondition1->addLinkedCondition(new PortConditionVolumeType3xFloat());
    orCondition1->addLinkedCondition(new PortConditionVolumeType3xDouble());
    orCondition1->addLinkedCondition(new PortConditionVolumeType2xFloat());
    orCondition1->addLinkedCondition(new PortConditionVolumeType2xDouble());
    inport1_.addCondition(orCondition1);
    addPort(outport1_);

    addProperty(enableProcessing_);
    addProperty(numberOfIterations_);
    addProperty(timeStep_);
    addProperty(conductanceParameter_);

}

Processor* VectorCurvatureAnisotropicDiffusionImageFilterITK::create() const {
    return new VectorCurvatureAnisotropicDiffusionImageFilterITK();
}

template<class T>
void VectorCurvatureAnisotropicDiffusionImageFilterITK::vectorCurvatureAnisotropicDiffusionImageFilterVec2ITK() {

    if (!enableProcessing_.get()) {
        outport1_.setData(inport1_.getData(), false);
        return;
    }

    typedef itk::Image<itk::CovariantVector<T,2>, 3> InputImageType1;
    typedef itk::Image<itk::CovariantVector<T,2>, 3> OutputImageType1;

    typename InputImageType1::Pointer p1 = voreenVec2ToITKVec2<T>(inport1_.getData());


    //Filter define
    typedef itk::VectorCurvatureAnisotropicDiffusionImageFilter<InputImageType1, OutputImageType1> FilterType;
    typename FilterType::Pointer filter = FilterType::New();

    filter->SetInput(p1);

    filter->SetNumberOfIterations(numberOfIterations_.get());
    filter->SetTimeStep(timeStep_.get());
    filter->SetConductanceParameter(conductanceParameter_.get());


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
void VectorCurvatureAnisotropicDiffusionImageFilterITK::vectorCurvatureAnisotropicDiffusionImageFilterVec3ITK() {

    if (!enableProcessing_.get()) {
        outport1_.setData(inport1_.getData(), false);
        return;
    }

    typedef itk::Image<itk::CovariantVector<T,3>, 3> InputImageType1;
    typedef itk::Image<itk::CovariantVector<T,3>, 3> OutputImageType1;

    typename InputImageType1::Pointer p1 = voreenVec3ToITKVec3<T>(inport1_.getData());


    //Filter define
    typedef itk::VectorCurvatureAnisotropicDiffusionImageFilter<InputImageType1, OutputImageType1> FilterType;
    typename FilterType::Pointer filter = FilterType::New();

    filter->SetInput(p1);

    filter->SetNumberOfIterations(numberOfIterations_.get());
    filter->SetTimeStep(timeStep_.get());
    filter->SetConductanceParameter(conductanceParameter_.get());


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
void VectorCurvatureAnisotropicDiffusionImageFilterITK::vectorCurvatureAnisotropicDiffusionImageFilterVec4ITK() {

    if (!enableProcessing_.get()) {
        outport1_.setData(inport1_.getData(), false);
        return;
    }

    typedef itk::Image<itk::CovariantVector<T,4>, 3> InputImageType1;
    typedef itk::Image<itk::CovariantVector<T,4>, 3> OutputImageType1;

    typename InputImageType1::Pointer p1 = voreenVec4ToITKVec4<T>(inport1_.getData());


    //Filter define
    typedef itk::VectorCurvatureAnisotropicDiffusionImageFilter<InputImageType1, OutputImageType1> FilterType;
    typename FilterType::Pointer filter = FilterType::New();

    filter->SetInput(p1);

    filter->SetNumberOfIterations(numberOfIterations_.get());
    filter->SetTimeStep(timeStep_.get());
    filter->SetConductanceParameter(conductanceParameter_.get());


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




void VectorCurvatureAnisotropicDiffusionImageFilterITK::process() {
    const VolumeBase* inputHandle1 = inport1_.getData();
    const VolumeRAM* inputVolume1 = inputHandle1->getRepresentation<VolumeRAM>();

    if (dynamic_cast<const VolumeRAM_4xFloat*>(inputVolume1))  {
        vectorCurvatureAnisotropicDiffusionImageFilterVec4ITK<float>();
    }
    else if (dynamic_cast<const VolumeRAM_4xDouble*>(inputVolume1))  {
        vectorCurvatureAnisotropicDiffusionImageFilterVec4ITK<double>();
    }
    else if (dynamic_cast<const VolumeRAM_3xFloat*>(inputVolume1))  {
        vectorCurvatureAnisotropicDiffusionImageFilterVec3ITK<float>();
    }
    else if (dynamic_cast<const VolumeRAM_3xDouble*>(inputVolume1))  {
        vectorCurvatureAnisotropicDiffusionImageFilterVec3ITK<double>();
    }
    else if (dynamic_cast<const VolumeRAM_2xFloat*>(inputVolume1))  {
        vectorCurvatureAnisotropicDiffusionImageFilterVec2ITK<float>();
    }
    else if (dynamic_cast<const VolumeRAM_2xDouble*>(inputVolume1))  {
        vectorCurvatureAnisotropicDiffusionImageFilterVec2ITK<double>();
    }
    else {
        LERROR("Inputformat of Volume 1 is not supported!");
    }

}


}   // namespace
