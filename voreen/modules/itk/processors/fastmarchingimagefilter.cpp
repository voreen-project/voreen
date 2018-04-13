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

#include "fastmarchingimagefilter.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include "modules/itk/utils/itkwrapper.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorconvert.h"
#include "itkImage.h"
#include "itkVectorImage.h"

#include "itkFastMarchingImageFilter.h"

#include <iostream>

namespace voreen {

const std::string FastMarchingImageFilterITK::loggerCat_("voreen.FastMarchingImageFilterITK");

FastMarchingImageFilterITK::FastMarchingImageFilterITK()
    : ITKProcessor(),
    inport1_(Port::INPORT, "InputImage"),
    outport1_(Port::OUTPORT, "OutputImage"),
    seedPointPort_(Port::INPORT, "seedPointInput"),
    enableProcessing_("enabled", "Enable", false),
    numSeedPoint_("numSeedPoint", "Take Seed-Point", 0, 0, 0, Processor::VALID),
    seedPoint_("SeedPoint", "Seed-Point", tgt::vec3(1),tgt::vec3(1),tgt::vec3(5000)),
    stoptime_("stoptime", "Stopping Value", 60.0f, 1.0f, 100000.0f),
    initialDistance_("initialDistance", "InitialDistance", 5.0f, 0.0f, 1000.0f)
{
    addPort(inport1_);
    PortConditionLogicalOr* orCondition1 = new PortConditionLogicalOr();
    orCondition1->addLinkedCondition(new PortConditionVolumeTypeFloat());
    orCondition1->addLinkedCondition(new PortConditionVolumeTypeDouble());
    inport1_.addCondition(orCondition1);
    addPort(outport1_);
    addPort(seedPointPort_);

    addProperty(enableProcessing_);

    //seed point
    addProperty(numSeedPoint_);
    addProperty(seedPoint_);
    numSeedPoint_.onChange(MemberFunctionCallback<FastMarchingImageFilterITK>(this, &FastMarchingImageFilterITK::process));
    seedPoint_.setReadOnlyFlag(!false);


    //properties for fast marching
    addProperty(stoptime_);
    addProperty (initialDistance_);
}

Processor* FastMarchingImageFilterITK::create() const {
    return new FastMarchingImageFilterITK();
}

template<class T>
void FastMarchingImageFilterITK::fastMarchingImageFilterITK() {

    if (!enableProcessing_.get()) {
        outport1_.setData(inport1_.getData(), false);
        return;
    }

    //typedefs images
    typedef itk::Image<T, 3> InputImageType;

    //input image p1
    typename InputImageType::Pointer p1 = voreenToITK<T>(inport1_.getData());

    typedef  itk::FastMarchingImageFilter< InputImageType, InputImageType >    FastMarchingFilterType;
    typename FastMarchingFilterType::Pointer  fastMarching = FastMarchingFilterType::New();

    typedef typename FastMarchingFilterType::NodeContainer   NodeContainer;
    typedef typename FastMarchingFilterType::NodeType        NodeType;

    typename NodeContainer::Pointer seeds = NodeContainer::New();

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


    typename InputImageType::IndexType  seedPosition;

    seedPosition[0] = seedPoint_.get().x;
    seedPosition[1] = seedPoint_.get().y;
    seedPosition[2] = seedPoint_.get().z;

    const double initialDistance = initialDistance_.get();

    NodeType node;

    const double seedValue = - initialDistance;

    node.SetValue( seedValue);
    node.SetIndex (seedPosition);

    seeds->Initialize();
    seeds->InsertElement( 0, node);

    fastMarching->SetTrialPoints( seeds );
    fastMarching->SetSpeedConstant( 1.0f );
    fastMarching->SetStoppingValue( stoptime_.get() );

    fastMarching->SetInput(p1);

    observe(fastMarching.GetPointer());

    try
    {
        fastMarching->Update();
    }
    catch (itk::ExceptionObject &e)
    {
        LERROR(e);
    }

    Volume* outputVolume1 = 0;
    outputVolume1 = ITKToVoreenCopy<T>(fastMarching->GetOutput());

    if (outputVolume1) {
        transferTransformation(inport1_.getData(), outputVolume1);
        outport1_.setData(outputVolume1);
    }
    else
        outport1_.setData(0);

}

void FastMarchingImageFilterITK::process() {
    const VolumeBase* inputHandle1 = inport1_.getData();
    const VolumeRAM* inputVolume1 = inputHandle1->getRepresentation<VolumeRAM>();

    if (dynamic_cast<const VolumeRAM_Float*>(inputVolume1))  {
        fastMarchingImageFilterITK<float>();
    }
    else {
        LERROR("Inputformat of Volume 1 is not supported!");
    }

}

}   // namespace
