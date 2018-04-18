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

#include "anisotropicdiffusion.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"

#include "voreen/core/datastructures/volume/operators/volumeoperatorconvert.h"
#include "modules/itk/utils/itkwrapper.h"

#include "itkImage.h"
#include "itkGradientAnisotropicDiffusionImageFilter.h"

#include <iostream>

namespace voreen {

const std::string AnisotropicDiffusion::loggerCat_("voreen.AnisotropicDiffusion");

AnisotropicDiffusion::AnisotropicDiffusion()
    : ITKProcessor(),
    inport_(Port::INPORT, "input"),
    outport_(Port::OUTPORT, "output"),
    timeStep_("timeStep", "Time Step", 0.05f, 0.0f, 1.0f),
    conductance_("conductance", "Conductance", 2.0f, 0.0f, 10.0f),
    steps_("steps", "Steps", 5, 1, 100)
{
    addPort(inport_);

    //example port conditions:
    PortConditionLogicalOr* orCondition = new PortConditionLogicalOr();
    orCondition->addLinkedCondition(new PortConditionVolumeTypeFloat());
    //orCondition->addLinkedCondition(new PortConditionVolumeTypeDouble());
    inport_.addCondition(orCondition);

    addPort(outport_);

    addProperty(timeStep_);
    addProperty(conductance_);
    addProperty(steps_);
}

Processor* AnisotropicDiffusion::create() const {
    return new AnisotropicDiffusion();
}

std::string AnisotropicDiffusion::getProcessorInfo() const {
    return std::string("ITK test processor (median filter)");
}

bool AnisotropicDiffusion::isReady() const {
    if(inport_.isReady() && outport_.isReady())
        return true;
    else
        return false;
}

void AnisotropicDiffusion::process() {
    const VolumeBase* inputVolume = inport_.getData();
    Volume* outputVolume = 0;

    typedef    float    InputPixelType;
    typedef    float    OutputPixelType;

    typedef itk::Image< InputPixelType,  3 >   InputImageType;
    typedef itk::Image< OutputPixelType, 3 >   OutputImageType;

    typedef itk::GradientAnisotropicDiffusionImageFilter< InputImageType, OutputImageType >  FilterType;
    FilterType::Pointer filter = FilterType::New();

    InputImageType::Pointer p = voreenToITK<float>(inputVolume);

    filter->SetInput(p);

    filter->SetNumberOfIterations(steps_.get());
    filter->SetTimeStep( timeStep_.get() );
    filter->SetConductanceParameter( conductance_.get() );

    observe(filter.GetPointer());
    // Now run the filter
    try
    {
        filter->Update();
    }
    catch (itk::ExceptionObject &e)
    {
        LERROR(e);
    }

    outputVolume = ITKToVoreenCopy<float>(filter->GetOutput());
    outputVolume->setRealWorldMapping(inputVolume->getRealWorldMapping());
    transferTransformation(inputVolume, outputVolume);

    // assign computed volume to outport
    if (outputVolume)
        outport_.setData(outputVolume);
    else
        outport_.setData(0);
}

}   // namespace
