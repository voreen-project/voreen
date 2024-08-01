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

#include "mutualinformationregistration.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

#include "voreen/core/datastructures/volume/operators/volumeoperatorconvert.h"
#include "voreen/core/datastructures/volume/volumedecorator.h"
#include "voreen/core/datastructures/meta/templatemetadata.h"
#include "modules/itk/utils/itkwrapper.h"

#include "itkMultiResolutionImageRegistrationMethod.h"

#include "itkMattesMutualInformationImageToImageMetric.h"
#include "itkMeanSquaresImageToImageMetric.h"
#include "itkNormalizedCorrelationImageToImageMetric.h"
#include "itkMeanReciprocalSquareDifferenceImageToImageMetric.h"
#include "itkMutualInformationImageToImageMetric.h"
#include "itkKullbackLeiblerCompareHistogramImageToImageMetric.h"
#include "itkMeanSquaresHistogramImageToImageMetric.h"
#include "itkCorrelationCoefficientHistogramImageToImageMetric.h"
#include "itkGradientDifferenceImageToImageMetric.h"

#include "itkLinearInterpolateImageFunction.h"
#include "itkMultiResolutionPyramidImageFilter.h"
#include "itkCenteredTransformInitializer.h"

#include <Eigen/SVD>
#include <Eigen/LU>

//#include <boost/date_time.hpp>

#include <iostream>

namespace voreen {

using tgt::mat4;
using tgt::vec3;

template <typename RegistrationType>
class RegistrationInterfaceCommand : public itk::Command
{
public:
    typedef  RegistrationInterfaceCommand   Self;
    typedef  itk::Command                   Superclass;
    typedef  itk::SmartPointer<Self>        Pointer;
    itkNewMacro( Self );
protected:
    float minStepLength_;
    float maxStepLength_;

    RegistrationInterfaceCommand() : minStepLength_(0.001f), maxStepLength_(0.2f) {};
public:
    typedef   itk::VersorRigid3DTransformOptimizer       OptimizerType;

    void setMinStepLength(float minStepLength) { minStepLength_ = minStepLength; }
    void setMaxStepLength(float maxStepLength) { maxStepLength_ = maxStepLength; }
    float getMinStepLength() { return minStepLength_; }
    float getMaxStepLength() { return maxStepLength_; }

    void Execute(itk::Object * object, const itk::EventObject & event)
    {
        if( !(itk::IterationEvent().CheckEvent( &event )) )
            return;

        RegistrationType* registration = dynamic_cast<RegistrationType*>(object);

        OptimizerType* optimizer = dynamic_cast<OptimizerType*>( registration->GetOptimizer() );

        std::cout << "-------------------------------------" << std::endl;
        std::cout << "MultiResolution Level : " << registration->GetCurrentLevel()  << std::endl;
        std::cout << std::endl;

        if ( registration->GetCurrentLevel() == 0 )
        {
            optimizer->SetMaximumStepLength(maxStepLength_);
            optimizer->SetMinimumStepLength(minStepLength_);
        }
        else
        {
            optimizer->SetMaximumStepLength( optimizer->GetMaximumStepLength() / 10.0 );
            optimizer->SetMinimumStepLength( optimizer->GetMinimumStepLength() / 2.0 );

            if(optimizer->GetStopCondition() == OptimizerType::Unknown)
                registration->StopRegistration();
        }
    }

    void Execute(const itk::Object * , const itk::EventObject & ) { return; }
};


//void MutualInformationRegistrationObserver::Execute(const itk::Object* object, const itk::EventObject & event)
void MutualInformationRegistrationObserver::Execute(itk::Object* object, const itk::EventObject& event)
{
    OptimizerType* optimizer = dynamic_cast<OptimizerType*>( object );
    if( !(itk::IterationEvent().CheckEvent( &event )) )
    {
        return;
    }
    std::cout << optimizer->GetCurrentIteration() << "   ";
    std::cout << optimizer->GetValue() << "   ";
    std::cout << optimizer->GetCurrentPosition() << std::endl;

    // Update transformation matrix:
    MutualInformationRegistration::ParametersType parameters = optimizer->GetCurrentPosition();
    MutualInformationRegistration::TransformType::Pointer transform = MutualInformationRegistration::TransformType::New();
    transform->SetParameters(parameters);

    try {
        boost::this_thread::sleep(boost::posix_time::seconds(0));
    }
    catch(boost::thread_interrupted&)
    {
        optimizer->StopOptimization();
    }

    processor_->lockMutex();
    processor_->tempResult_ = processor_->calculateVoreenTrafo(processor_->itkToVoreen(transform));
    processor_->tempResultUpdated_ = true;
    processor_->invalidate();
    processor_->unlockMutex();
}

const std::string MutualInformationRegistration::loggerCat_("voreen.MutualInformationRegistration");

MutualInformationRegistration::MutualInformationRegistration()
    : VolumeProcessor(),
    fixedVolumeInport_(Port::INPORT, "fixedVolumeInport", "Fixed Volume"),
    movingVolumeInport_(Port::INPORT, "movingVolumeInport", "Moving Volume"),
    transformationMatrix_("voreenTransformMatrix", "Voreen transformation matrix", tgt::mat4::identity, tgt::mat4(-2000.0), tgt::mat4(2000.0), VALID),
    numLevels_("numLevels", "Number of levels", 5, 1, 20, VALID),
    numIterations_("numIterations", "Number of iterations", 200, 1, 5000, VALID),
    minStepLength_("minStepLength", "Min Step Length", 0.001f, 0.0f, 10.0f, VALID),
    maxStepLength_("maxStepLength", "Max Step Length", 0.001f, 0.0f, 1000.0f, VALID),
    rotScale_("rotScale", "Rotation Scale", 0.001f, 0.0f, 1000.0f, VALID),
    metric_("metric", "Metric"),
    numHistogramBins_("numHistogramBins", "Histogram bins", 64, 1, 512, VALID),
    numSamples_("numSamples", "Number of samples [%]", 10, 1, 100, VALID),
    relaxationFactor_("relaxationFactor", "Relaxation factor", 0.9f, 0.0f, 2.0f, VALID),
    explicitPDF_("explicitPDF", "Use explicit PDF", false, VALID),
    lambda_("lambda", "Lambda", 1.0f, 0.0f, 100.0f),
    updateButton_("updateButton", "Update"),
    stopButton_("stopButton", "Stop"),
    initializeButton_("initializeButton", "Initialize"),
    tempResult_(mat4::identity),
    tempResultUpdated_(false)
{
    addPort(fixedVolumeInport_);
    addPort(movingVolumeInport_);

    addProperty(numLevels_);
    numLevels_.setGroupID("framework");
    setPropertyGroupGuiName("framework", "Framework");

    addProperty(rotScale_);
    rotScale_.setGroupID("optimizer");
    rotScale_.setStepping(0.001f);
    addProperty(minStepLength_);
    minStepLength_.setGroupID("optimizer");
    minStepLength_.setStepping(0.001f);
    addProperty(maxStepLength_);
    maxStepLength_.setGroupID("optimizer");
    maxStepLength_.setStepping(0.001f);
    addProperty(numIterations_);
    numIterations_.setGroupID("optimizer");
    addProperty(relaxationFactor_);
    relaxationFactor_.setGroupID("optimizer");
    setPropertyGroupGuiName("optimizer", "Optimizer");

    //metric_.addOption("violaWells", "Viola & Wells"); // needs more work...
    metric_.addOption("meanSquares", "Mean Squares");
    metric_.addOption("mattes", "Mattes");
    metric_.addOption("meanSquaresHistogram", "Mean Squares Histogram");
    metric_.addOption("normalizedCorrelation", "Normalized Correlation");
    //metric_.addOption("correlationCoefficient", "Correlation Coefficient Histogram"); // needs more work...
    metric_.addOption("gradientDifference", "Gradient Difference");
    //metric_.addOption("kullbackLeibler", "Kullback-Leibler"); // needs more work...
    metric_.onChange(MemberFunctionCallback<MutualInformationRegistration>(this, &MutualInformationRegistration::onMetricChange) );
    addProperty(metric_);
    metric_.setGroupID("metric");
    addProperty(numHistogramBins_);
    numHistogramBins_.setGroupID("metric");
    numHistogramBins_.setVisibleFlag(false);
    addProperty(numSamples_);
    numSamples_.setGroupID("metric");
    numSamples_.setVisibleFlag(false);
    addProperty(explicitPDF_);
    explicitPDF_.setGroupID("metric");
    explicitPDF_.setVisibleFlag(false);
    addProperty(lambda_);
    lambda_.setGroupID("metric");
    lambda_.setVisibleFlag(false);
    setPropertyGroupGuiName("metric", "Metric");

    addProperty(transformationMatrix_);
    addProperty(updateButton_);
    updateButton_.onChange(MemberFunctionCallback<MutualInformationRegistration>(this, &MutualInformationRegistration::updateRegistration));
    addProperty(stopButton_);
    stopButton_.onChange(MemberFunctionCallback<MutualInformationRegistration>(this, &MutualInformationRegistration::stopRegistration));
    addProperty(initializeButton_);
    initializeButton_.onChange(MemberFunctionCallback<MutualInformationRegistration>(this, &MutualInformationRegistration::initializeRegistration));
}


void MutualInformationRegistration::deinitialize() {
    VolumeProcessor::deinitialize();
}

Processor* MutualInformationRegistration::create() const {
    return new MutualInformationRegistration();
}

const VolumeBase* MutualInformationRegistration::getFixedVolume() const {
    return fixedVolumeInport_.getData();
}

const VolumeBase* MutualInformationRegistration::getMovingVolume() const {
    return fixedVolumeInport_.getData();
}

void MutualInformationRegistration::process() {
    // mutex is locked by NetworkEvaluator
    if(tempResultUpdated_) {
        transformationMatrix_.set(tempResult_);
        tempResultUpdated_ = false;
    }
}

void MutualInformationRegistration::onMetricChange() {
    numHistogramBins_.setVisibleFlag(false);
    numSamples_.setVisibleFlag(false);
    explicitPDF_.setVisibleFlag(false);
    lambda_.setVisibleFlag(false);

    if(metric_.get() == "mattes") {
        numHistogramBins_.setVisibleFlag(true);
        numSamples_.setVisibleFlag(true);
        explicitPDF_.setVisibleFlag(true);
    }
    else if(metric_.get() == "violaWells") {
    }
    else if(metric_.get() == "meanSquares") {
    }
    else if(metric_.get() == "meanSquaresHistogram") {
        numHistogramBins_.setVisibleFlag(true);
    }
    else if(metric_.get() == "normalizedCorrelation") {
        lambda_.setVisibleFlag(true);
    }
    else if(metric_.get() == "correlationCoefficient") {
    }
    else if(metric_.get() == "gradientDifference") {
    }
    else if(metric_.get() == "kullbackLeibler") {
    }
}

template<class fixedType, class movingType>
void MutualInformationRegistration::performRegistration(const VolumeBase* fixedVolume, const VolumeBase* movingVolume) {
    typedef itk::Image<fixedType, 3> FixedImageType;
    typedef itk::Image<movingType, 3> MovingImageType;

    typedef itk::VersorRigid3DTransformOptimizer OptimizerType;
    typedef OptimizerType::ScalesType       OptimizerScalesType;
    typedef itk::LinearInterpolateImageFunction<MovingImageType, double> InterpolatorType;
    typedef itk::MultiResolutionImageRegistrationMethod<FixedImageType, MovingImageType>   RegistrationType;

    typedef itk::MultiResolutionPyramidImageFilter<FixedImageType, FixedImageType>   FixedImagePyramidType; // second argument necessary?
    typedef itk::MultiResolutionPyramidImageFilter<MovingImageType, MovingImageType>   MovingImagePyramidType;

    typename OptimizerType::Pointer      optimizer     = OptimizerType::New();
    typename InterpolatorType::Pointer   interpolator  = InterpolatorType::New();
    typename RegistrationType::Pointer   registration  = RegistrationType::New();
    TransformType::Pointer      transform     = TransformType::New();

    typename FixedImagePyramidType::Pointer fixedImagePyramid = FixedImagePyramidType::New();
    typename MovingImagePyramidType::Pointer movingImagePyramid = MovingImagePyramidType::New();

    //---------------------------------------------------------------------------------------------
    // Metric setup:

    if(metric_.get() == "mattes") {
        typedef itk::MattesMutualInformationImageToImageMetric< FixedImageType, MovingImageType >   MetricType;
        typename MetricType::Pointer         metric        = MetricType::New();

        metric->SetNumberOfHistogramBins(numHistogramBins_.get());

        size_t numVoxels = hmul(fixedVolume->getDimensions());
        metric->SetNumberOfSpatialSamples(numVoxels * (static_cast<float>(numSamples_.get()) / 100.0f));

        metric->ReinitializeSeed( 76926294 );

        //// Define whether to calculate the metric derivative by explicitly
        //// computing the derivatives of the joint PDF with respect to the Transform
        //// parameters, or doing it by progressively accumulating contributions from
        //// each bin in the joint PDF.
        metric->SetUseExplicitPDFDerivatives(explicitPDF_.get());

        registration->SetMetric(metric);
    }
    else if(metric_.get() == "meanSquares") {
        typedef itk::MeanSquaresImageToImageMetric< FixedImageType, MovingImageType >   MetricType;
        typename MetricType::Pointer         metric        = MetricType::New();

        registration->SetMetric(metric);
    }
    else if(metric_.get() == "meanSquaresHistogram") {
        typedef itk::MeanSquaresHistogramImageToImageMetric< FixedImageType, MovingImageType >   MetricType;
        typename MetricType::Pointer         metric        = MetricType::New();
        itk::Array<unsigned long> ar(2);
        ar(0) = numHistogramBins_.get();
        ar(1) = numHistogramBins_.get();
        metric->SetHistogramSize(ar);

        registration->SetMetric(metric);
    }
    else if(metric_.get() == "normalizedCorrelation") {
        typedef itk::NormalizedCorrelationImageToImageMetric< FixedImageType, MovingImageType >   MetricType;
        typename MetricType::Pointer         metric        = MetricType::New();

        registration->SetMetric(metric);
    }
    else if(metric_.get() == "normalizedCorrelation") {
        typedef itk::MeanReciprocalSquareDifferenceImageToImageMetric< FixedImageType, MovingImageType >   MetricType;
        typename MetricType::Pointer         metric        = MetricType::New();
        metric->SetLambda(lambda_.get());

        registration->SetMetric(metric);
    }
    else if(metric_.get() == "violaWells") {
        typedef itk::MutualInformationImageToImageMetric< FixedImageType, MovingImageType >   MetricType;
        typename MetricType::Pointer         metric        = MetricType::New();

        registration->SetMetric(metric);
    }
    else if(metric_.get() == "kullbackLeibler") {
        typedef itk::KullbackLeiblerCompareHistogramImageToImageMetric< FixedImageType, MovingImageType >   MetricType;
        typename MetricType::Pointer         metric        = MetricType::New();

        registration->SetMetric(metric);
    }
    else if(metric_.get() == "correlationCoefficient") {
        typedef itk::CorrelationCoefficientHistogramImageToImageMetric< FixedImageType, MovingImageType >   MetricType;
        typename MetricType::Pointer         metric        = MetricType::New();

        registration->SetMetric(metric);
    }
    else if(metric_.get() == "gradientDifference") {
        typedef itk::GradientDifferenceImageToImageMetric< FixedImageType, MovingImageType >   MetricType;
        typename MetricType::Pointer         metric        = MetricType::New();

        registration->SetMetric(metric);
    }
    else {
        LERROR("Unknown metric!");
        return;
    }

    //---------------------------------------------------------------------------------------------

    registration->SetOptimizer(optimizer);
    registration->SetTransform(transform);
    registration->SetInterpolator(interpolator);
    registration->SetFixedImagePyramid(fixedImagePyramid);
    registration->SetMovingImagePyramid(movingImagePyramid);

    OptimizerScalesType optimizerScales( transform->GetNumberOfParameters() );

    optimizerScales[0] = rotScale_.get();
    optimizerScales[1] = rotScale_.get();
    optimizerScales[2] = rotScale_.get();
    optimizerScales[3] = 1.0f;
    optimizerScales[4] = 1.0f;
    optimizerScales[5] = 1.0f;
    optimizer->SetScales( optimizerScales );
    optimizer->SetMaximumStepLength(maxStepLength_.get());
    optimizer->SetMinimumStepLength(minStepLength_.get());
    //optimizer->MinimizeOn();

    typename FixedImageType::Pointer fixed = voreenToITK<fixedType>(fixedVolume);
    typename MovingImageType::Pointer moving = voreenToITK<movingType>(movingVolume);
    registration->SetFixedImage(fixed);
    registration->SetMovingImage(moving);
    registration->SetFixedImageRegion( fixed->GetBufferedRegion() );

    // we want to use the current matrix (in the property) as starting point:

    // get the matrix from prop and account for existing tranformation of fixed image:
    tgt::mat4 m = calculateITKTrafo();

    // ITK (i.e., VersorRigid3DTransformOptimizer) needs a orthogonal matrix and checks this with high precision
    // we therefore compute a polar decomposition of the existing matrix in double precision:
    Eigen::Matrix3d a; // copy matrix to Eigen
    a(0, 0) = m[0][0];
    a(0, 1) = m[0][1];
    a(0, 2) = m[0][2];

    a(1, 0) = m[1][0];
    a(1, 1) = m[1][1];
    a(1, 2) = m[1][2];

    a(2, 0) = m[2][0];
    a(2, 1) = m[2][1];
    a(2, 2) = m[2][2];

    Eigen::JacobiSVD<Eigen::Matrix3d> svd(a, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Vector3d rhs;
    svd.solve(rhs);
    Eigen::Matrix3d u = svd.matrixU();
    Eigen::Matrix3d v = svd.matrixV();

    Eigen::Matrix3d rot = u * v.transpose(); // rotation matrix with scaling/inaccuracies removed

    TransformType::MatrixType mITK; // copy this to ITK
    mITK[0][0] = rot(0, 0);
    mITK[0][1] = rot(0, 1);
    mITK[0][2] = rot(0, 2);

    mITK[1][0] = rot(1, 0);
    mITK[1][1] = rot(1, 1);
    mITK[1][2] = rot(1, 2);

    mITK[2][0] = rot(2, 0);
    mITK[2][1] = rot(2, 1);
    mITK[2][2] = rot(2, 2);

    transform->SetMatrix(mITK);

    TransformType::OutputVectorType trITK;
    trITK[0] = m[0][3];
    trITK[1] = m[1][3];
    trITK[2] = m[2][3];
    transform->SetTranslation(trITK);

    registration->SetInitialTransformParameters( transform->GetParameters() );

    optimizer->SetNumberOfIterations(numIterations_.get());
    optimizer->SetRelaxationFactor(relaxationFactor_.get());

    // Create the Command observer and register it with the optimizer.
    MutualInformationRegistrationObserver::Pointer observer = MutualInformationRegistrationObserver::New();
    observer->setProcessor(this);
    optimizer->AddObserver( itk::IterationEvent(), observer );

    typedef RegistrationInterfaceCommand<RegistrationType> CommandType;
    typename CommandType::Pointer command = CommandType::New();
    command->setMinStepLength(minStepLength_.get());
    command->setMaxStepLength(maxStepLength_.get());
    registration->AddObserver( itk::IterationEvent(), command );

    registration->SetNumberOfLevels(numLevels_.get());

    try
    {
        registration->StartRegistration();
        std::cout << "Optimizer stop condition: " << registration->GetOptimizer()->GetStopConditionDescription() << std::endl;
    }
    catch( itk::ExceptionObject & err )
    {
        LERROR("ITK ExceptionObject caught!");
        LERROR(err);
        return;
    }

    ParametersType finalParameters = registration->GetLastTransformParameters();
    transform->SetParameters(finalParameters);

    lockMutex();
    tempResult_ = calculateVoreenTrafo(itkToVoreen(transform));
    tempResultUpdated_ = true;
    invalidate();
    unlockMutex();
}

void MutualInformationRegistration::stopRegistration() {
    workerThread_.interrupt();
}

void MutualInformationRegistration::updateRegistration() {
    if(!isReady())
        return;

    if(!workerThread_.timed_join(boost::posix_time::seconds(0))) {
        LERROR("Thread is already running!");
        return;
    }

    const VolumeBase* fixedVolume = fixedVolumeInport_.getData();
    const VolumeBase* movingVolume = movingVolumeInport_.getData();

    if(dynamic_cast<const VolumeRAM_UInt16*>(fixedVolume->getRepresentation<VolumeRAM>()))
        typeSwitch<uint16_t>(fixedVolume, movingVolume);
    else if(dynamic_cast<const VolumeRAM_Int16*>(fixedVolume->getRepresentation<VolumeRAM>()))
        typeSwitch<int16_t>(fixedVolume, movingVolume);
    else if(dynamic_cast<const VolumeRAM_UInt8*>(fixedVolume->getRepresentation<VolumeRAM>()))
        typeSwitch<uint8_t>(fixedVolume, movingVolume);
    else if(dynamic_cast<const VolumeRAM_Int8*>(fixedVolume->getRepresentation<VolumeRAM>()))
        typeSwitch<int8_t>(fixedVolume, movingVolume);
    else if(dynamic_cast<const VolumeRAM_UInt32*>(fixedVolume->getRepresentation<VolumeRAM>()))
        typeSwitch<uint32_t>(fixedVolume, movingVolume);
    else if(dynamic_cast<const VolumeRAM_Int32*>(fixedVolume->getRepresentation<VolumeRAM>()))
        typeSwitch<int32_t>(fixedVolume, movingVolume);
    else if(dynamic_cast<const VolumeRAM_Float*>(fixedVolume->getRepresentation<VolumeRAM>()))
        typeSwitch<float>(fixedVolume, movingVolume);
    else if(dynamic_cast<const VolumeRAM_Double*>(fixedVolume->getRepresentation<VolumeRAM>()))
        typeSwitch<double>(fixedVolume, movingVolume);
    else {
        LERROR("Unsupported fixed volume type!");
    }

}

template<class fixedType>
void MutualInformationRegistration::typeSwitch(const VolumeBase* fixedVolume, const VolumeBase* movingVolume) {
    if(dynamic_cast<const VolumeRAM_UInt16*>(movingVolume->getRepresentation<VolumeRAM>()))
        workerThread_ = boost::thread(&MutualInformationRegistration::performRegistration<fixedType, uint16_t>, this, fixedVolume, movingVolume);
    else if(dynamic_cast<const VolumeRAM_Int16*>(movingVolume->getRepresentation<VolumeRAM>()))
        workerThread_ = boost::thread(&MutualInformationRegistration::performRegistration<fixedType, int16_t>, this, fixedVolume, movingVolume);
    else if(dynamic_cast<const VolumeRAM_UInt8*>(movingVolume->getRepresentation<VolumeRAM>()))
        workerThread_ = boost::thread(&MutualInformationRegistration::performRegistration<fixedType, uint8_t>, this, fixedVolume, movingVolume);
    else if(dynamic_cast<const VolumeRAM_Int8*>(movingVolume->getRepresentation<VolumeRAM>()))
        workerThread_ = boost::thread(&MutualInformationRegistration::performRegistration<fixedType, int8_t>, this, fixedVolume, movingVolume);
    else if(dynamic_cast<const VolumeRAM_UInt32*>(movingVolume->getRepresentation<VolumeRAM>()))
        workerThread_ = boost::thread(&MutualInformationRegistration::performRegistration<fixedType, uint32_t>, this, fixedVolume, movingVolume);
    else if(dynamic_cast<const VolumeRAM_Int32*>(movingVolume->getRepresentation<VolumeRAM>()))
        workerThread_ = boost::thread(&MutualInformationRegistration::performRegistration<fixedType, int32_t>, this, fixedVolume, movingVolume);
    else if(dynamic_cast<const VolumeRAM_Float*>(movingVolume->getRepresentation<VolumeRAM>()))
        workerThread_ = boost::thread(&MutualInformationRegistration::performRegistration<fixedType, float>, this, fixedVolume, movingVolume);
    else if(dynamic_cast<const VolumeRAM_Double*>(movingVolume->getRepresentation<VolumeRAM>()))
        workerThread_ = boost::thread(&MutualInformationRegistration::performRegistration<fixedType, double>, this, fixedVolume, movingVolume);
    else {
        LERROR("Unsupported moving volume type!");
    }
}

mat4 MutualInformationRegistration::itkToVoreen(TransformType::Pointer t) {
    mat4 m = mat4::identity;

    TransformType::MatrixType mITK;
    mITK = t->GetMatrix();
    mITK[0][0] = m[0][0];
    mITK[0][1] = m[0][1];
    mITK[0][2] = m[0][2];

    mITK[1][0] = m[1][0];
    mITK[1][1] = m[1][1];
    mITK[1][2] = m[1][2];

    mITK[2][0] = m[2][0];
    mITK[2][1] = m[2][1];
    mITK[2][2] = m[2][2];

    TransformType::OutputVectorType trITK;
    trITK = t->GetTranslation();
    m[0][3] = trITK[0];
    m[1][3] = trITK[1];
    m[2][3] = trITK[2];

    return m;
}

template<class fixedType>
void MutualInformationRegistration::typeSwitchInit(const VolumeBase* fixedVolume, const VolumeBase* movingVolume) {
    if(dynamic_cast<const VolumeRAM_UInt16*>(movingVolume->getRepresentation<VolumeRAM>()))
        performInitialization<fixedType, uint16_t>(fixedVolume, movingVolume);
    else if(dynamic_cast<const VolumeRAM_Int16*>(movingVolume->getRepresentation<VolumeRAM>()))
        performInitialization<fixedType, int16_t>(fixedVolume, movingVolume);
    else if(dynamic_cast<const VolumeRAM_UInt8*>(movingVolume->getRepresentation<VolumeRAM>()))
        performInitialization<fixedType, uint8_t>(fixedVolume, movingVolume);
    else if(dynamic_cast<const VolumeRAM_Int8*>(movingVolume->getRepresentation<VolumeRAM>()))
        performInitialization<fixedType, int8_t>(fixedVolume, movingVolume);
    else if(dynamic_cast<const VolumeRAM_UInt32*>(movingVolume->getRepresentation<VolumeRAM>()))
        performInitialization<fixedType, uint32_t>(fixedVolume, movingVolume);
    else if(dynamic_cast<const VolumeRAM_Int32*>(movingVolume->getRepresentation<VolumeRAM>()))
        performInitialization<fixedType, int32_t>(fixedVolume, movingVolume);
    else if(dynamic_cast<const VolumeRAM_Float*>(movingVolume->getRepresentation<VolumeRAM>()))
        performInitialization<fixedType, float>(fixedVolume, movingVolume);
    else if(dynamic_cast<const VolumeRAM_Double*>(movingVolume->getRepresentation<VolumeRAM>()))
        performInitialization<fixedType, double>(fixedVolume, movingVolume);
    else {
        LERROR("Unsupported moving volume type!");
    }
}

template<class fixedType, class movingType>
void MutualInformationRegistration::performInitialization(const VolumeBase* fixedVolume, const VolumeBase* movingVolume) {
    typedef itk::Image<fixedType, 3> FixedImageType;
    typedef itk::Image<movingType, 3> MovingImageType;

    typename FixedImageType::Pointer fixed = voreenToITK<fixedType>(fixedVolume);
    typename MovingImageType::Pointer moving = voreenToITK<movingType>(movingVolume);

    typedef itk::CenteredTransformInitializer< TransformType, FixedImageType, MovingImageType >  TransformInitializerType;

    typename TransformInitializerType::Pointer initializer = TransformInitializerType::New();
    TransformType::Pointer      transform     = TransformType::New();

    initializer->SetTransform(transform);
    initializer->SetFixedImage(fixed);
    initializer->SetMovingImage(moving);

    initializer->MomentsOn();
    initializer->InitializeTransform();

    typedef TransformType::VersorType  VersorType;
    typedef VersorType::VectorType     VectorType;
    VersorType     rotation;
    VectorType     axis;
    axis[0] = 0.0;
    axis[1] = 0.0;
    axis[2] = 1.0;
    const double angle = 0;
    rotation.Set(  axis, angle  );
    transform->SetRotation( rotation );

    transformationMatrix_.set(calculateVoreenTrafo(itkToVoreen(transform)));
}

void MutualInformationRegistration::initializeRegistration() {
    if(!isReady())
        return;

    const VolumeBase* fixedVolume = fixedVolumeInport_.getData();
    const VolumeBase* movingVolume = movingVolumeInport_.getData();

    if(dynamic_cast<const VolumeRAM_UInt16*>(fixedVolume->getRepresentation<VolumeRAM>()))
        typeSwitchInit<uint16_t>(fixedVolume, movingVolume);
    else if(dynamic_cast<const VolumeRAM_Int16*>(fixedVolume->getRepresentation<VolumeRAM>()))
        typeSwitchInit<int16_t>(fixedVolume, movingVolume);
    else if(dynamic_cast<const VolumeRAM_UInt8*>(fixedVolume->getRepresentation<VolumeRAM>()))
        typeSwitchInit<uint8_t>(fixedVolume, movingVolume);
    else if(dynamic_cast<const VolumeRAM_Int8*>(fixedVolume->getRepresentation<VolumeRAM>()))
        typeSwitchInit<int8_t>(fixedVolume, movingVolume);
    else if(dynamic_cast<const VolumeRAM_UInt32*>(fixedVolume->getRepresentation<VolumeRAM>()))
        typeSwitchInit<uint32_t>(fixedVolume, movingVolume);
    else if(dynamic_cast<const VolumeRAM_Int32*>(fixedVolume->getRepresentation<VolumeRAM>()))
        typeSwitchInit<int32_t>(fixedVolume, movingVolume);
    else if(dynamic_cast<const VolumeRAM_Float*>(fixedVolume->getRepresentation<VolumeRAM>()))
        typeSwitchInit<float>(fixedVolume, movingVolume);
    else if(dynamic_cast<const VolumeRAM_Double*>(fixedVolume->getRepresentation<VolumeRAM>()))
        typeSwitchInit<double>(fixedVolume, movingVolume);
    else {
        LERROR("Unsupported fixed volume type!");
    }
}

tgt::mat4 MutualInformationRegistration::calculateVoreenTrafo(const mat4& itkMatrix) {
    const VolumeBase* fixed = fixedVolumeInport_.getData();

    mat4 invertedITKTransform;
    bool success = itkMatrix.invert(invertedITKTransform);
    if(!success) {
        LERROR("Failed to invert ITK transformation matrix!");
        return tgt::mat4::identity;
    }
    mat4 f_physicalToWorld = fixed->getPhysicalToWorldMatrix();

    return f_physicalToWorld * invertedITKTransform;
}

tgt::mat4 MutualInformationRegistration::calculateITKTrafo() {
    const VolumeBase* fixed = fixedVolumeInport_.getData();
    const VolumeBase* moving = movingVolumeInport_.getData();

    mat4 t = transformationMatrix_.get();
    mat4 tInv;
    bool success = t.invert(tInv);
    if(!success) {
        LERROR("Failed to invert transformation matrix!");
        return tgt::mat4::identity;
    }

    return tInv * fixed->getPhysicalToWorldMatrix();
}

}   // namespace
