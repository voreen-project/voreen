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

#include "narrowbandcurveslevelsetimagefilterworkflow.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/ports/conditions/portconditionvolumetype.h"
#include "modules/itk/utils/itkwrapper.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorconvert.h"
#include "itkImage.h"
#include "itkVectorImage.h"

#include "itkCurvatureAnisotropicDiffusionImageFilter.h"
#include "itkGradientMagnitudeRecursiveGaussianImageFilter.h"
#include "itkSigmoidImageFilter.h"
#include "itkDanielssonDistanceMapImageFilter.h"
#include "itkGradientImageFilter.h"
#include "itkEdgePotentialImageFilter.h"
#include "itkNarrowBandCurvesLevelSetImageFilter.h"
#include "itkFixedArray.h"
#include "itkLevelSetFunction.h"
#include "itkFastMarchingImageFilter.h"

#include <iostream>

namespace voreen {

const std::string NarrowBandCurvesLevelSetImageFilterWorkflowITK::loggerCat_("voreen.NarrowBandCurvesLevelSetImageFilterWorkflowITK");

NarrowBandCurvesLevelSetImageFilterWorkflowITK::NarrowBandCurvesLevelSetImageFilterWorkflowITK()
    : ITKProcessor(),
    inport1_(Port::INPORT, "InputImage"),
    outport1_(Port::OUTPORT, "OutputImage"),
    seedPointPort_(Port::INPORT, "seedPointInput"),
    enableProcessing_("enabled", "Enable", false),
    initLevelSet_("initialLevelSetMethod", "Initial LevelSet Method"),
    numSeedPoint_("numSeedPoint", "Take Seed-Point", 0, 0, 0, Processor::VALID),
    seedPoint_("SeedPoint", "Seed-Point", tgt::vec3(1),tgt::vec3(1),tgt::vec3(5000)),
    inputIsBinary_("inputIsBinary", "InputIsBinary", false),
    squaredDistance_("squaredDistance", "SquaredDistance", false),
    useImageSpacing_("useImageSpacing", "UseImageSpacing", false),
    timestep_("timestep", "Timestep", 0.0625f, 0.0f, 0.0625f),
    numberofiterations_ ("numberofiterations", "Number of Iterations", 5, 1, 65535),
    alpha_ ("alpha", "Alpha", -0.5f, -150.00f, 150.0f),
    beta_ ("beta", "Beta", 3.0f, -150.0f, 150.0f),
    sigma_("sigma", "Sigma", 0.5f, 0.0f, 50.0f),
    conductanceparameter_ ("conductanceparameter", "Conductance Parameter", 2.0f, 0.0f, 1000.0f),
    stoptime_("stoptime", "Stopping Value", 60.0f, 1.0f, 100000.0f),
    initialDistance_("initialDistance", "InitialDistance", 5.0f, 0.0f, 1000.0f),
    propagationScaling_("propagationScaling", "PropagationScaling", 2.0f, 0.0f, 100.0f)
{
    addPort(inport1_);
    PortConditionLogicalOr* orCondition1 = new PortConditionLogicalOr();
    orCondition1->addLinkedCondition(new PortConditionVolumeTypeFloat());
    orCondition1->addLinkedCondition(new PortConditionVolumeTypeDouble());
    inport1_.addCondition(orCondition1);
    addPort(outport1_);
    addPort(seedPointPort_);

    addProperty(enableProcessing_);
    initLevelSet_.addOption("danielsson", "SignedDanielssonDistanceMap");
    initLevelSet_.addOption("fastmarching", "FastMarching");
    addProperty(initLevelSet_);

    //seed point
    addProperty(numSeedPoint_);
    addProperty(seedPoint_);
    numSeedPoint_.onChange(MemberFunctionCallback<NarrowBandCurvesLevelSetImageFilterWorkflowITK>(this, &NarrowBandCurvesLevelSetImageFilterWorkflowITK::process));
    seedPoint_.setReadOnlyFlag(!false);

    //properties for DanielsonSignedDistance
    addProperty(inputIsBinary_);
    addProperty(squaredDistance_);
    addProperty(useImageSpacing_);

    inputIsBinary_.setGroupID("Danielsson");
    squaredDistance_.setGroupID("Danielsson");
    useImageSpacing_.setGroupID("Danielsson");
    setPropertyGroupGuiName("Danielsson", "Danielsson");

    //properties for smoothing, sigmoid, gradientMagnitude

    addProperty(timestep_);
    addProperty(numberofiterations_);
    addProperty(conductanceparameter_);
    timestep_.setGroupID("smooth");
    numberofiterations_.setGroupID("smooth");
    conductanceparameter_.setGroupID("smooth");
    setPropertyGroupGuiName("smooth", "Smoothing");

    addProperty(alpha_);
    addProperty(beta_);
    alpha_.setGroupID("sigmoid");
    beta_.setGroupID("sigmoid");
    setPropertyGroupGuiName("sigmoid", "Sigmoid");

    addProperty(sigma_);
    sigma_.setGroupID("magnitude");
    setPropertyGroupGuiName("magnitude", "Gradient Magnitude");

    //properties for fast marching
    addProperty(stoptime_);
    addProperty (initialDistance_);
    stoptime_.setGroupID("fastmarching");
    initialDistance_.setGroupID("fastmarching");
    setPropertyGroupGuiName("fastmarching", "Fast Marching");


    //properties for segmentation
    addProperty(propagationScaling_);
    propagationScaling_.setGroupID("segmentation");
    setPropertyGroupGuiName("segmentation", "GeodesicActiveContourLevelSet");
}

Processor* NarrowBandCurvesLevelSetImageFilterWorkflowITK::create() const {
    return new NarrowBandCurvesLevelSetImageFilterWorkflowITK();
}



template<class T>
void NarrowBandCurvesLevelSetImageFilterWorkflowITK::narrowBandCurvesLevelSetImageFilterWorkflowITK() {

    inputIsBinary_.setVisibleFlag(true);
    squaredDistance_.setVisibleFlag(true);
    useImageSpacing_.setVisibleFlag(true);

    if (!enableProcessing_.get()) {
        outport1_.setData(inport1_.getData(), false);
        return;
    }

    //typedefs images
    typedef itk::Image<T, 3> InputImageType;

    //input image p1
    typename InputImageType::Pointer p1 = voreenToITK<T>(inport1_.getData());

    //LevelsetFilter
    typedef itk::NarrowBandCurvesLevelSetImageFilter<InputImageType, InputImageType> LevelType;
    typename LevelType::Pointer levelfilter = LevelType::New();
    typename LevelType::Pointer levelfilter2 = LevelType::New();

    //smoothing, gradient, sigmoid
    typedef   itk::CurvatureAnisotropicDiffusionImageFilter<
                               InputImageType,
                               InputImageType >  SmoothingFilterType;

    typename SmoothingFilterType::Pointer smoothing = SmoothingFilterType::New();

    typedef   itk::GradientMagnitudeRecursiveGaussianImageFilter<
                               InputImageType,
                               InputImageType >  GradientFilterType;

    typedef   itk::SigmoidImageFilter<
                               InputImageType,
                               InputImageType >  SigmoidFilterType;

    typename GradientFilterType::Pointer  gradientMagnitude = GradientFilterType::New();

    typename SigmoidFilterType::Pointer sigmoid = SigmoidFilterType::New();
    sigmoid->SetOutputMinimum( 0.0 );
    sigmoid->SetOutputMaximum( 1.0 );

    smoothing->SetInput(p1);
    smoothing->SetTimeStep( timestep_.get() );
    smoothing->SetNumberOfIterations(  numberofiterations_.get() );
    smoothing->SetConductanceParameter( conductanceparameter_.get() );

    gradientMagnitude->SetInput( smoothing->GetOutput() );
    gradientMagnitude->SetSigma( sigma_.get());
    sigmoid->SetInput( gradientMagnitude->GetOutput() );
    sigmoid->SetAlpha( alpha_.get() );
    sigmoid->SetBeta( beta_.get() );
    sigmoid->SetOutputMinimum(  0.0  );
    sigmoid->SetOutputMaximum(  1.0f  );

    if (initLevelSet_.isSelected("danielsson")){

        //initial Level Set using DanielssonDistanceMapImageFiler
        typedef itk::DanielssonDistanceMapImageFilter<InputImageType, InputImageType> DanielssonType;
        typename DanielssonType::Pointer danielssonfilter = DanielssonType::New();


        danielssonfilter->SetInput(p1);
        danielssonfilter->SetInputIsBinary(inputIsBinary_.get());
        danielssonfilter->SetSquaredDistance(squaredDistance_.get());
        danielssonfilter->SetUseImageSpacing(useImageSpacing_.get());


        levelfilter->SetInput(danielssonfilter->GetOutput());
        levelfilter->SetFeatureImage (sigmoid->GetOutput());
        levelfilter->SetPropagationScaling(propagationScaling_.get());
        levelfilter->SetCurvatureScaling( 1.0 );
        levelfilter->SetAdvectionScaling( 1.0 );
        levelfilter->SetMaximumRMSError( 0.02 );
        levelfilter->SetNumberOfIterations( 800 );

        observe(levelfilter.GetPointer());

        try
        {
            levelfilter->Update();
        }
        catch (itk::ExceptionObject &e)
        {
            LERROR(e);
        }


        Volume* outputVolume1 = 0;
        outputVolume1 = ITKToVoreenCopy<T>(levelfilter->GetOutput());

        if (outputVolume1) {
            transferTransformation(inport1_.getData(), outputVolume1);
            outport1_.setData(outputVolume1);
        }
        else
            outport1_.setData(0);
    }

    else if (initLevelSet_.isSelected("fastmarching")){

        inputIsBinary_.setVisibleFlag(false);
        squaredDistance_.setVisibleFlag(false);
        useImageSpacing_.setVisibleFlag(false);

        //initial Level set using FastMarchingImageFilter
        typedef  itk::FastMarchingImageFilter<
                              InputImageType,
                              InputImageType >    FastMarchingFilterType;

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

        observe(levelfilter2.GetPointer());

        levelfilter2->SetInput(fastMarching->GetOutput());
        levelfilter2->SetFeatureImage (sigmoid->GetOutput());
        levelfilter2->SetPropagationScaling(propagationScaling_.get());
        levelfilter2->SetCurvatureScaling( 1.0 );
        levelfilter2->SetAdvectionScaling( 1.0 );
        levelfilter2->SetMaximumRMSError( 0.02 );
        levelfilter2->SetNumberOfIterations( 800 );


        try
        {
            levelfilter2->Update();
        }
        catch (itk::ExceptionObject &e)
        {
            LERROR(e);
        }

        Volume* outputVolume1 = 0;
        outputVolume1 = ITKToVoreenCopy<T>(levelfilter2->GetOutput());

        if (outputVolume1) {
            transferTransformation(inport1_.getData(), outputVolume1);
            outport1_.setData(outputVolume1);
        }
        else
            outport1_.setData(0);

    }
}



void NarrowBandCurvesLevelSetImageFilterWorkflowITK::process() {
    const VolumeBase* inputHandle1 = inport1_.getData();
    const VolumeRAM* inputVolume1 = inputHandle1->getRepresentation<VolumeRAM>();

    if (dynamic_cast<const VolumeRAM_Float*>(inputVolume1))  {
        narrowBandCurvesLevelSetImageFilterWorkflowITK<float>();
    }
    else {
        LERROR("Inputformat of Volume 1 is not supported!");
    }

}


}   // namespace
