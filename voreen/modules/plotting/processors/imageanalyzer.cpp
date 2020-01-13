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

#include "imageanalyzer.h"


#include "voreen/core/utils/voreenqualitymode.h"
#include "../datastructures/aggregationfunction.h"
#include "../datastructures/plotdata.h"
#include "../datastructures/plotrow.h"

#include <limits>

namespace voreen {

const std::string ImageAnalyzer::loggerCat_("voreen.ImageAnalyzer");

ImageAnalyzer::ImageAnalyzer():
    ImageProcessor()
    , image_(Port::INPORT,"image.input")
    , inportMask_(Port::INPORT,"mask.input")
    , outport_(Port::OUTPORT,"PlotData.outport")
    , histogramOutPort_(Port::OUTPORT,"PlotDatahistogram.outport")
    , computeButton_("compute", "Compute information",Processor::VALID)
    , allOutPutZeroOne_("allOutPutZeroOne","All in [0,1] on PlotData.Outport",true,Processor::VALID)
    , interAction_("interActionModeAction","InterActionMode Action",Processor::VALID)
    , calculateStep_("calculatestep","Calculate Step",1,1,200,Processor::VALID)
    , histogramSteps_("histogramSteps","Histogram Steps",256,256,std::numeric_limits<int>::max(),Processor::VALID)
    , bytePerPixel_("Bytes Per Pixel","Bytes Per Pixel",0,0,64)
    , width_("Width","Width",0,0,std::numeric_limits<int>::max())
    , height_("Height","Height",0,0,std::numeric_limits<int>::max())
    , stepPerChannel_("stepperchannel","StepPerChannel",256,0,std::numeric_limits<int>::max())

    , minRed_("minRed","Min Red", 0.0f)
    , maxRed_("maxRed","Max Red", 0.0f)
    , medianRed_("medianRed","Median Red", 0.0f)
    , minGreen_("minGreen","Min Green", 0.0f)
    , maxGreen_("maxGreen","Max Green", 0.0f)
    , medianGreen_("medianGreen","Median Green", 0.0f)
    , minBlue_("minBlue","Min Blue", 0.0f)
    , maxBlue_("maxBlue","Max Blue", 0.0f)
    , medianBlue_("medianBlue","Median Blue", 0.0f)
    , minAlpha_("minAlpha","Min Alpha", 0.0f)
    , maxAlpha_("maxAlpha","Max Alpha", 0.0f)
    , medianAlpha_("medianAlpha","Median Alpha", 0.0f)

    , minGrey_("minGrey","Min Grey", 0.0f)
    , maxGrey_("maxGrey","Max Grey", 0.0f)
    , medianGrey_("medianGrey","Median Grey", 0.0f)
    , minLuminance_("minLuminance","Min Luminance", 0.0f)
    , maxLuminance_("maxLuminance","Max Luminance", 0.0f)
    , medianLuminance_("medianLuminance","Median Luminance", 0.0f)

    , minDepth_("minDepth", "Min Depth", 0.0f)
    , maxDepth_("maxDepth", "Max Depth", 0.0f)
    , medianDepth_("medianDepth_","Median Depth", 0.0f)

    , minC_("minC", "Min Cyan", 0.0f)
    , maxC_("maxC", "Max Cyan", 0.0f)
    , medianC_("medianC_","Median C", 0.0f)
    , minY_("minY", "Min Yellow", 0.0f)
    , maxY_("maxY", "Max Yellow", 0.0f)
    , medianY_("medianY_","Median Y", 0.0f)
    , minM_("minM", "Min Magenta", 0.0f)
    , maxM_("maxM", "Max Magenta", 0.0f)
    , medianM_("medianM_","Median M", 0.0f)
    , minK_("minK", "Min Key", 0.0f)
    , maxK_("maxK", "Max Key", 0.0f)
    , medianK_("medianK_","Median K", 0.0f)

    , minH_("minH", "Min H", 0.0f,0.0f,360.0f)
    , maxH_("maxH", "Max H", 0.0f,0.0f,360.0f)
    , medianH_("medianH_","Median H", 0.0f)
    , minS_("minS", "Min S", 0.0f,0.0f,1.0f)
    , maxS_("maxS", "Max S", 0.0f,0.0f,1.0f)
    , medianS_("medianS_","Median S", 0.0f)
    , minV_("minV", "Min V", 0.0f,0.0f,1.0f)
    , maxV_("maxV", "Max V", 0.0f,0.0f,1.0f)
    , medianV_("medianV_","Median V", 0.0f)

    , toHistoRed_("toHistRed","Add Red to Histogramm",true,Processor::VALID)
    , toHistoGreen_("toHistGreen","Add Green to Histogramm",true,Processor::VALID)
    , toHistoBlue_("toHistBlue","Add Blue to Histogramm",true,Processor::VALID)
    , toHistoAlpha_("toHistAlpha","Add Alpha to Histogramm",true,Processor::VALID)
    , toHistoGrey_("toHistGrey","Add Grey to Histogramm",true,Processor::VALID)
    , toHistoLuminance_("toHistLuminance","Add Luminance to Histogramm",true,Processor::VALID)
    , toHistoDepth_("toHistDepth","Add Depth to Histogramm",true,Processor::VALID)

    , toHistoC_("toHistC","Add C to Histogramm",true,Processor::VALID)
    , toHistoY_("toHistY","Add Y to Histogramm",true,Processor::VALID)
    , toHistoM_("toHistM","Add M to Histogramm",true,Processor::VALID)
    , toHistoK_("toHistK","Add K to Histogramm",true,Processor::VALID)

    , toHistoH_("toHistH","Add H to Histogramm",true,Processor::VALID)
    , toHistoS_("toHistS","Add S to Histogramm",true,Processor::VALID)
    , toHistoV_("toHistV","Add V to Histogramm",true,Processor::VALID)

    , pData_(0)
    , pDataHistogram_(0)

    , histoSteps_(0)

    , histogramRed_()
    , histogramGreen_()
    , histogramBlue_()
    , histogramAlpha_()
    , histogramGrey_()
    , histogramLuminance_()
    , histogramDepth_()

    , histogramC_()
    , histogramY_()
    , histogramM_()
    , histogramK_()

    , histogramH_()
    , histogramS_()
    , histogramV_()

{
    computeButton_.onChange(MemberFunctionCallback<ImageAnalyzer>(this, &ImageAnalyzer::calculate));
    calculateStep_.onChange(MemberFunctionCallback<ImageAnalyzer>(this, &ImageAnalyzer::activateComputeButton));
    allOutPutZeroOne_.onChange(MemberFunctionCallback<ImageAnalyzer>(this, &ImageAnalyzer::activateComputeButton));
    histogramSteps_.onChange(MemberFunctionCallback<ImageAnalyzer>(this, &ImageAnalyzer::activateComputeButton));
    toHistoRed_.onChange(MemberFunctionCallback<ImageAnalyzer>(this, &ImageAnalyzer::activateComputeButton));
    toHistoGreen_.onChange(MemberFunctionCallback<ImageAnalyzer>(this, &ImageAnalyzer::activateComputeButton));
    toHistoBlue_.onChange(MemberFunctionCallback<ImageAnalyzer>(this, &ImageAnalyzer::activateComputeButton));
    toHistoAlpha_.onChange(MemberFunctionCallback<ImageAnalyzer>(this, &ImageAnalyzer::activateComputeButton));
    toHistoGrey_.onChange(MemberFunctionCallback<ImageAnalyzer>(this, &ImageAnalyzer::activateComputeButton));
    toHistoLuminance_.onChange(MemberFunctionCallback<ImageAnalyzer>(this, &ImageAnalyzer::activateComputeButton));
    toHistoDepth_.onChange(MemberFunctionCallback<ImageAnalyzer>(this, &ImageAnalyzer::activateComputeButton));

    toHistoC_.onChange(MemberFunctionCallback<ImageAnalyzer>(this, &ImageAnalyzer::activateComputeButton));
    toHistoY_.onChange(MemberFunctionCallback<ImageAnalyzer>(this, &ImageAnalyzer::activateComputeButton));
    toHistoM_.onChange(MemberFunctionCallback<ImageAnalyzer>(this, &ImageAnalyzer::activateComputeButton));
    toHistoK_.onChange(MemberFunctionCallback<ImageAnalyzer>(this, &ImageAnalyzer::activateComputeButton));

    toHistoH_.onChange(MemberFunctionCallback<ImageAnalyzer>(this, &ImageAnalyzer::activateComputeButton));
    toHistoS_.onChange(MemberFunctionCallback<ImageAnalyzer>(this, &ImageAnalyzer::activateComputeButton));
    toHistoV_.onChange(MemberFunctionCallback<ImageAnalyzer>(this, &ImageAnalyzer::activateComputeButton));

    bytePerPixel_.setReadOnlyFlag(true);
    height_.setReadOnlyFlag(true);
    width_.setReadOnlyFlag(true);
    stepPerChannel_.setReadOnlyFlag(true);

    minLuminance_.setReadOnlyFlag(true);
    maxLuminance_.setReadOnlyFlag(true);
    medianLuminance_.setReadOnlyFlag(true);
    minDepth_.setReadOnlyFlag(true);
    maxDepth_.setReadOnlyFlag(true);
    medianDepth_.setReadOnlyFlag(true);
    minAlpha_.setReadOnlyFlag(true);
    maxAlpha_.setReadOnlyFlag(true);
    medianAlpha_.setReadOnlyFlag(true);
    minGrey_.setReadOnlyFlag(true);
    maxGrey_.setReadOnlyFlag(true);
    medianGrey_.setReadOnlyFlag(true);
    minRed_.setReadOnlyFlag(true);
    maxRed_.setReadOnlyFlag(true);
    medianRed_.setReadOnlyFlag(true);
    minGreen_.setReadOnlyFlag(true);
    maxGreen_.setReadOnlyFlag(true);
    medianGreen_.setReadOnlyFlag(true);
    minBlue_.setReadOnlyFlag(true);
    maxBlue_.setReadOnlyFlag(true);
    medianBlue_.setReadOnlyFlag(true);

    minC_.setReadOnlyFlag(true);
    maxC_.setReadOnlyFlag(true);
    medianC_.setReadOnlyFlag(true);
    minM_.setReadOnlyFlag(true);
    maxM_.setReadOnlyFlag(true);
    medianM_.setReadOnlyFlag(true);
    minY_.setReadOnlyFlag(true);
    maxY_.setReadOnlyFlag(true);
    medianY_.setReadOnlyFlag(true);
    minK_.setReadOnlyFlag(true);
    maxK_.setReadOnlyFlag(true);
    medianK_.setReadOnlyFlag(true);

    minH_.setReadOnlyFlag(true);
    maxH_.setReadOnlyFlag(true);
    medianH_.setReadOnlyFlag(true);
    minS_.setReadOnlyFlag(true);
    maxS_.setReadOnlyFlag(true);
    medianS_.setReadOnlyFlag(true);
    minV_.setReadOnlyFlag(true);
    maxV_.setReadOnlyFlag(true);
    medianV_.setReadOnlyFlag(true);

    interAction_.addOption("setzero","Set Null to Outport",0);
    interAction_.addOption("oldValue","Old Value",5);
    interAction_.addOption("Calculateon","Calculate On",10);
    interAction_.selectByValue(5);

    addProperty(&computeButton_);
    addProperty(&allOutPutZeroOne_);
    addProperty(&interAction_);

    addProperty(&calculateStep_);
    addProperty(&histogramSteps_);
    addProperty(&bytePerPixel_);
    addProperty(&width_);
    addProperty(&height_);
    addProperty(&stepPerChannel_);

    addProperty(&minRed_);
    addProperty(&maxRed_);
    addProperty(&medianRed_);
    addProperty(&minGreen_);
    addProperty(&maxGreen_);
    addProperty(&medianGreen_);
    addProperty(&minBlue_);
    addProperty(&maxBlue_);
    addProperty(&medianBlue_);
    addProperty(&minAlpha_);
    addProperty(&maxAlpha_);
    addProperty(&medianAlpha_);

    addProperty(&minGrey_);
    addProperty(&maxGrey_);
    addProperty(&medianGrey_);
    addProperty(&minLuminance_);
    addProperty(&maxLuminance_);
    addProperty(&medianLuminance_);

    addProperty(&minDepth_);
    addProperty(&maxDepth_);
    addProperty(&medianDepth_);

    addProperty(&minC_);
    addProperty(&maxC_);
    addProperty(&medianC_);
    addProperty(&minY_);
    addProperty(&maxY_);
    addProperty(&medianY_);
    addProperty(&minM_);
    addProperty(&maxM_);
    addProperty(&medianM_);
    addProperty(&minK_);
    addProperty(&maxK_);
    addProperty(&medianK_);

    addProperty(&minH_);
    addProperty(&maxH_);
    addProperty(&medianH_);
    addProperty(&minS_);
    addProperty(&maxS_);
    addProperty(&medianS_);
    addProperty(&minV_);
    addProperty(&maxV_);
    addProperty(&medianV_);

    addProperty(&toHistoRed_);
    addProperty(&toHistoGreen_);
    addProperty(&toHistoBlue_);
    addProperty(&toHistoAlpha_);
    addProperty(&toHistoGrey_);
    addProperty(&toHistoLuminance_);
    addProperty(&toHistoDepth_);

    addProperty(&toHistoC_);
    addProperty(&toHistoY_);
    addProperty(&toHistoM_);
    addProperty(&toHistoK_);

    addProperty(&toHistoH_);
    addProperty(&toHistoS_);
    addProperty(&toHistoV_);

    addPort(&image_);
    addPort(&outport_);
    addPort(&histogramOutPort_);

    computeButton_.setGroupID("Calculation Parameter");
    allOutPutZeroOne_.setGroupID("Calculation Parameter");
    interAction_.setGroupID("Calculation Parameter");
    calculateStep_.setGroupID("Calculation Parameter");
    histogramSteps_.setGroupID("Calculation Parameter");

    bytePerPixel_.setGroupID("Generic Information");
    height_.setGroupID("Generic Information");
    width_.setGroupID("Generic Information");
    stepPerChannel_.setGroupID("Generic Information");

    minRed_.setGroupID("RGBA-Values");
    maxRed_.setGroupID("RGBA-Values");
    medianRed_.setGroupID("RGBA-Values");
    minGreen_.setGroupID("RGBA-Values");
    maxGreen_.setGroupID("RGBA-Values");
    medianGreen_.setGroupID("RGBA-Values");
    minBlue_.setGroupID("RGBA-Values");
    maxBlue_.setGroupID("RGBA-Values");
    medianBlue_.setGroupID("RGBA-Values");
    minAlpha_.setGroupID("RGBA-Values");
    maxAlpha_.setGroupID("RGBA-Values");
    medianAlpha_.setGroupID("RGBA-Values");

    minGrey_.setGroupID("RGBA-Values");
    maxGrey_.setGroupID("RGBA-Values");
    medianGrey_.setGroupID("RGBA-Values");
    minLuminance_.setGroupID("RGBA-Values");
    maxLuminance_.setGroupID("RGBA-Values");
    medianLuminance_.setGroupID("RGBA-Values");

    minDepth_.setGroupID("RGBA-Values");
    maxDepth_.setGroupID("RGBA-Values");
    medianDepth_.setGroupID("RGBA-Values");

    minC_.setGroupID("CYMK-Values");
    maxC_.setGroupID("CYMK-Values");
    medianC_.setGroupID("CYMK-Values");
    minM_.setGroupID("CYMK-Values");
    maxM_.setGroupID("CYMK-Values");
    medianM_.setGroupID("CYMK-Values");
    minY_.setGroupID("CYMK-Values");
    maxY_.setGroupID("CYMK-Values");
    medianY_.setGroupID("CYMK-Values");
    minK_.setGroupID("CYMK-Values");
    maxK_.setGroupID("CYMK-Values");
    medianK_.setGroupID("CYMK-Values");

    minH_.setGroupID("HSV-Values");
    maxH_.setGroupID("HSV-Values");
    medianH_.setGroupID("HSV-Values");
    minS_.setGroupID("HSV-Values");
    maxS_.setGroupID("HSV-Values");
    medianS_.setGroupID("HSV-Values");
    minV_.setGroupID("HSV-Values");
    maxV_.setGroupID("HSV-Values");
    medianV_.setGroupID("HSV-Values");

    toHistoRed_.setGroupID("Histogramselection");
    toHistoGreen_.setGroupID("Histogramselection");
    toHistoBlue_.setGroupID("Histogramselection");
    toHistoAlpha_.setGroupID("Histogramselection");
    toHistoGrey_.setGroupID("Histogramselection");
    toHistoLuminance_.setGroupID("Histogramselection");
    toHistoDepth_.setGroupID("Histogramselection");

    toHistoC_.setGroupID("Histogramselection");
    toHistoY_.setGroupID("Histogramselection");
    toHistoM_.setGroupID("Histogramselection");
    toHistoK_.setGroupID("Histogramselection");

    toHistoH_.setGroupID("Histogramselection");
    toHistoS_.setGroupID("Histogramselection");
    toHistoV_.setGroupID("Histogramselection");
    Processor::setProgress(1.f);

    setPropertyGroupGuiName("Calculation Parameter","Calculation Parameter");
    setPropertyGroupGuiName("Generic Information","Generic Information");
    setPropertyGroupGuiName("RGBA-Values","RGBA-Values");
    setPropertyGroupGuiName("CYMK-Values","CYMK-Values");
    setPropertyGroupGuiName("HSV-Values","HSV-Values");
    setPropertyGroupGuiName("Histogramselection","Histogramselection");
}

Processor* ImageAnalyzer::create() const {
    return new ImageAnalyzer();
}

void ImageAnalyzer::initialize() {
    Processor::initialize();

    pData_ = new PlotData(1,1);
    pDataHistogram_ = new PlotData(1,1);
    histoSteps_ = histogramSteps_.get();
}

void ImageAnalyzer::deinitialize() {
    outport_.setData(0, false);
    histogramOutPort_.setData(0, false);
    delete pData_;
    delete pDataHistogram_;

    pData_ = 0;
    pDataHistogram_ = 0;

    Processor::deinitialize();
}

bool ImageAnalyzer::isReady() const {
    return image_.isConnected();
}

bool ImageAnalyzer::isEndProcessor() const {
    return ((!outport_.isConnected()) && (!histogramOutPort_.isConnected()));
}

void ImageAnalyzer::colorScale(tgt::Color* colorValue, unsigned int modValue, tgt::Vector4<unsigned int>* result){
    result->a = static_cast<unsigned int>(colorValue->a*modValue);
    result->r = static_cast<unsigned int>(colorValue->r*modValue);
    result->g = static_cast<unsigned int>(colorValue->g*modValue);
    result->b = static_cast<unsigned int>(colorValue->b*modValue);
}

void ImageAnalyzer::cymkScale(tgt::Color* colorValue, tgt::Vector4f* result){
    float valueC = 1.0f - colorValue->r;//( R / 255 )
    float valueM = 1.0f - colorValue->g;//( G / 255 )
    float valueY = 1.0f - colorValue->b;//( B / 255 )
    //CMYK and CMY values from 0 to 1

    float var_K = 1;

    if ( valueC < var_K )   var_K = valueC;
    if ( valueM < var_K )   var_K = valueM;
    if ( valueY < var_K )   var_K = valueY;
    if ( var_K == 1 ) { //Black
       valueC = 0;
       valueM = 0;
       valueY = 0;
    }
    else {
       valueC = ( valueC - var_K ) / ( 1.0f - var_K );
       valueM = ( valueM - var_K ) / ( 1.0f - var_K );
       valueY = ( valueY - var_K ) / ( 1.0f - var_K );
    }
    result->x = valueC;
    result->y = valueY;
    result->z = valueM;
    result->w = var_K;
}

void ImageAnalyzer::hsvScale(tgt::Color* colorValue, tgt::Vector3f* result) {
    float maxValue = std::max<float>(colorValue->r,colorValue->g);
    maxValue = std::max<float>(maxValue,colorValue->b);
    float minValue = std::min<float>(colorValue->r,colorValue->g);
    minValue = std::min<float>(minValue,colorValue->b);
    float valueH = 0.0f;
    float valueS = 0.0f;
    float valueV = 0.0f;
    if (maxValue == minValue)
        valueH = 0.0f;
    else if (maxValue == colorValue->r)
        valueH = 60.0f*(colorValue->g-colorValue->b)/(maxValue-minValue);
    else if (maxValue == colorValue->g)
        valueH = 120.0f + 60.0f*(colorValue->b-colorValue->r)/(maxValue-minValue);
    else if (maxValue == colorValue->b)
        valueH = 240.0f +  60.0f*(colorValue->r-colorValue->g)/(maxValue-minValue);
    if (valueH < 0)
        valueH = valueH + 360.0f;
    if (maxValue != 0)
        valueS = (maxValue-minValue)/maxValue;
    valueV = maxValue;
    result->x = valueH;
    result->y = valueS;
    result->z = valueV;
}

float ImageAnalyzer::greyScale(tgt::Color* colorVector){
    float resultred = (colorVector->r * 0.393f) + (colorVector->g * 0.769f) + (colorVector->b * 0.189f);
    float resultgreen = (colorVector->r * 0.349f) + (colorVector->g * 0.686f) + (colorVector->b * 0.168f);
    float resultblue = (colorVector->r * 0.272f) + (colorVector->g * 0.534f) + (colorVector->b * 0.131f);

    if (resultred > 1)
        resultred = 1;
    if (resultgreen > 1)
        resultgreen = 1;
    if (resultblue > 1)
        resultblue = 1;

    return (resultred + resultgreen + resultblue)/3;
}

float ImageAnalyzer::luminance(tgt::Color* colorVector){
    return (colorVector->r*0.2126f + colorVector->g*0.7152f + colorVector->b*0.0722f);
}

void ImageAnalyzer::activateComputeButton() {
    computeButton_.setReadOnlyFlag(false);
    histoSteps_ = histogramSteps_.get();
}

void ImageAnalyzer::process(){
    if (!image_.hasRenderTarget()) {
        return;
    }
    computeButton_.setReadOnlyFlag(false);

    width_.set(0);
    height_.set(0);
    bytePerPixel_.set(0);
    stepPerChannel_.set(0);
    minRedValue_ = 0;
    maxRedValue_ = 0;
    medianRedValue_ = 0;
    minGreenValue_ = 0;
    maxGreenValue_ = 0;
    medianGreenValue_ = 0;
    minBlueValue_ = 0;
    maxBlueValue_ = 0;
    medianBlueValue_ = 0;
    minAlphaValue_ = 0;
    maxAlphaValue_ = 0;
    medianAlphaValue_ = 0;

    minGreyValue_ = 0;
    maxGreyValue_ = 0;
    medianGreyValue_ = 0;

    minLuminanceValue_ = 0;
    maxLuminanceValue_ = 0;
    medianLuminanceValue_ = 0;

    minDepthValue_ = 0;
    maxDepthValue_ = 0;
    medianDepthValue_ = 0;

    //CYMK-Calculation
    minCValue_ = 0;
    maxCValue_ = 0;
    medianCValue_ = 0;
    minYValue_ = 0;
    maxYValue_ = 0;
    medianYValue_ = 0;
    minMValue_ = 0;
    maxMValue_ = 0;
    medianMValue_ = 0;
    minKValue_ = 0;
    maxKValue_ = 0;
    medianKValue_ = 0;

    //HSV-Calculation
    minHValue_ = 0;
    maxHValue_ = 0;
    medianHValue_ = 0;
    minSValue_ = 0;
    maxSValue_ = 0;
    medianSValue_ = 0;
    minVValue_ = 0;
    maxVValue_ = 0;
    medianVValue_ = 0;

    if (!QualityMode.isInteractionMode() || (interAction_.getValue() == 10))
    {
        calculate();
    }
    else if (interAction_.getValue() < 5){
        outport_.setData(0);
        histogramOutPort_.setData(0);
        setPropertyValues();
    }
    else
        setPropertyValues();
}

void ImageAnalyzer::setOutPortData(){

    if (!isInitialized())
        return;

    PlotData* pDataHistogramAll = new PlotData(1,14);
    std::vector<std::pair<int, AggregationFunction*> > implizitcells = std::vector<std::pair<int, AggregationFunction*> >(0);
    std::pair<int, AggregationFunction*> function;
    for (int i = 1; i <= 14; ++i){
        function = std::pair<int, AggregationFunction*>(i, new AggregationFunctionMedianHistogram());
        implizitcells.push_back(function);
    }
    pDataHistogramAll->setColumnLabel(0,"Index");
    pDataHistogramAll->setColumnLabel(1,"Red");
    pDataHistogramAll->setColumnLabel(2,"Green");
    pDataHistogramAll->setColumnLabel(3,"Blue");
    pDataHistogramAll->setColumnLabel(4,"Alpha");
    pDataHistogramAll->setColumnLabel(5,"Grey");
    pDataHistogramAll->setColumnLabel(6,"Luminance");
    pDataHistogramAll->setColumnLabel(7,"Depth");
    pDataHistogramAll->setColumnLabel(8,"C");
    pDataHistogramAll->setColumnLabel(9,"Y");
    pDataHistogramAll->setColumnLabel(10,"M");
    pDataHistogramAll->setColumnLabel(11,"K");
    pDataHistogramAll->setColumnLabel(12,"H");
    pDataHistogramAll->setColumnLabel(13,"S");
    pDataHistogramAll->setColumnLabel(14,"V");
    std::vector<PlotCellValue> cells = std::vector<PlotCellValue>(0);
    for (size_t i = 0; i < static_cast<size_t>(histoSteps_); ++i) {
        cells.push_back(PlotCellValue(static_cast<plot_t>(i)));
        cells.push_back(PlotCellValue(histogramRed_.at(i)));
        cells.push_back(PlotCellValue(histogramGreen_.at(i)));
        cells.push_back(PlotCellValue(histogramBlue_.at(i)));
        cells.push_back(PlotCellValue(histogramAlpha_.at(i)));
        cells.push_back(PlotCellValue(histogramGrey_.at(i)));
        cells.push_back(PlotCellValue(histogramLuminance_.at(i)));
        cells.push_back(PlotCellValue(histogramDepth_.at(i)));
        cells.push_back(PlotCellValue(histogramC_.at(i)));
        cells.push_back(PlotCellValue(histogramY_.at(i)));
        cells.push_back(PlotCellValue(histogramM_.at(i)));
        cells.push_back(PlotCellValue(histogramK_.at(i)));
        cells.push_back(PlotCellValue(histogramH_.at(i)));
        cells.push_back(PlotCellValue(histogramS_.at(i)));
        cells.push_back(PlotCellValue(histogramV_.at(i)));
        pDataHistogramAll->insert(cells);
        cells.clear();
    }
    Processor::setProgress(0.96);
    pDataHistogramAll->insertImplicit(implizitcells);
    for (int i = 1; i <= 14; ++i){
        function = std::pair<int, AggregationFunction*>(i, new AggregationFunctionMinHistogram());
        delete implizitcells[i-1].second;
        implizitcells[i-1] = function;
    }
    pDataHistogramAll->insertImplicit(implizitcells);
    for (int i = 1; i <= 14; ++i){
        function = std::pair<int, AggregationFunction*>(i, new AggregationFunctionMaxHistogram());
        delete implizitcells[i-1].second;
        implizitcells[i-1] = function;
    }
    pDataHistogramAll->insertImplicit(implizitcells);
    for (int i = 1; i <= 14; ++i){
        function = std::pair<int, AggregationFunction*>(i, new  AggregationFunctionCountHistogram());
        delete implizitcells[i-1].second;
        implizitcells[i-1] = function;
    }
    pDataHistogramAll->insertImplicit(implizitcells);
    for (int i = 1; i <= 14; ++i){
        function = std::pair<int, AggregationFunction*>(i, new  AggregationFunctionSumHistogram());
        delete implizitcells[i-1].second;
        implizitcells[i-1] = function;
    }
    pDataHistogramAll->insertImplicit(implizitcells);
    for (int i = 1; i <= 14; ++i){
        function = std::pair<int, AggregationFunction*>(i, new  AggregationFunctionAverage());
        delete implizitcells[i-1].second;
        implizitcells[i-1] = function;
    }
    pDataHistogramAll->insertImplicit(implizitcells);
    for (int i = 1; i <= 14; ++i){
        function = std::pair<int, AggregationFunction*>(i, new  AggregationFunctionMode());
        delete implizitcells[i-1].second;
        implizitcells[i-1] = function;
    }
    pDataHistogramAll->insertImplicit(implizitcells);
    for (size_t i = 0; i < implizitcells.size(); ++i) {
        delete implizitcells[i].second;
    }
    medianRedValue_ = static_cast<float>(pDataHistogramAll->getImplicitRow(0).getCellAt(1).getValue()/(histoSteps_ - 1.0f));
    medianGreenValue_ = static_cast<float>(pDataHistogramAll->getImplicitRow(0).getCellAt(2).getValue()/(histoSteps_ - 1.0f));
    medianBlueValue_ = static_cast<float>(pDataHistogramAll->getImplicitRow(0).getCellAt(3).getValue()/(histoSteps_ - 1.0f));
    medianAlphaValue_ = static_cast<float>(pDataHistogramAll->getImplicitRow(0).getCellAt(4).getValue()/(histoSteps_ - 1.0f));
    medianGreyValue_ = static_cast<float>(pDataHistogramAll->getImplicitRow(0).getCellAt(5).getValue()/(histoSteps_ - 1.0f));
    medianLuminanceValue_ = static_cast<float>(pDataHistogramAll->getImplicitRow(0).getCellAt(6).getValue()/(histoSteps_ - 1.0f));
    medianDepthValue_ = static_cast<float>(pDataHistogramAll->getImplicitRow(0).getCellAt(7).getValue()/(histoSteps_ - 1.0f));
    medianCValue_ = static_cast<float>(pDataHistogramAll->getImplicitRow(0).getCellAt(8).getValue()/(histoSteps_ - 1.0f));
    medianYValue_ = static_cast<float>(pDataHistogramAll->getImplicitRow(0).getCellAt(9).getValue()/(histoSteps_ - 1.0f));
    medianMValue_ = static_cast<float>(pDataHistogramAll->getImplicitRow(0).getCellAt(10).getValue()/(histoSteps_ - 1.0f));
    medianKValue_ = static_cast<float>(pDataHistogramAll->getImplicitRow(0).getCellAt(11).getValue()/(histoSteps_ - 1.0f));
    medianHValue_ = static_cast<float>(pDataHistogramAll->getImplicitRow(0).getCellAt(12).getValue()/(histoSteps_ - 1.0f));
    medianSValue_ = static_cast<float>(pDataHistogramAll->getImplicitRow(0).getCellAt(13).getValue()/(histoSteps_ - 1.0f));
    medianVValue_ = static_cast<float>(pDataHistogramAll->getImplicitRow(0).getCellAt(14).getValue()/(histoSteps_ - 1.0f));
    Processor::setProgress(0.97);
    if (outport_.isConnected()){
        pData_->reset(1,4);
        pData_->setColumnLabel(0,"Index");
        pData_->setColumnLabel(1,"Description");
        pData_->setColumnLabel(2,"Min");
        pData_->setColumnLabel(3,"Max");
        pData_->setColumnLabel(4,"Median");
        std::vector<PlotCellValue> cells = std::vector<PlotCellValue>(0);

        cells.push_back((PlotCellValue(1)));
        cells.push_back((PlotCellValue("Red")));
        cells.push_back((PlotCellValue(minRedValue_)));
        cells.push_back((PlotCellValue(maxRedValue_)));
        cells.push_back((PlotCellValue(medianRedValue_)));
        pData_->insert(cells);
        cells.clear();


        cells.push_back((PlotCellValue(2)));
        cells.push_back((PlotCellValue("Green")));
        cells.push_back((PlotCellValue(minGreenValue_)));
        cells.push_back((PlotCellValue(maxGreenValue_)));
        cells.push_back((PlotCellValue(medianGreenValue_)));
        pData_->insert(cells);
        cells.clear();


        cells.push_back((PlotCellValue(3)));
        cells.push_back((PlotCellValue("Blue")));
        cells.push_back((PlotCellValue(minBlueValue_)));
        cells.push_back((PlotCellValue(maxBlueValue_)));
        cells.push_back((PlotCellValue(medianBlueValue_)));
        pData_->insert(cells);
        cells.clear();


        cells.push_back((PlotCellValue(4)));
        cells.push_back((PlotCellValue("Luminance")));
        cells.push_back((PlotCellValue(minLuminanceValue_)));
        cells.push_back((PlotCellValue(maxLuminanceValue_)));
        cells.push_back((PlotCellValue(medianLuminanceValue_)));
        pData_->insert(cells);
        cells.clear();


        cells.push_back((PlotCellValue(5)));
        cells.push_back((PlotCellValue("Alpha")));
        cells.push_back((PlotCellValue(minAlphaValue_)));
        cells.push_back((PlotCellValue(maxAlphaValue_)));
        cells.push_back((PlotCellValue(medianAlphaValue_)));
        pData_->insert(cells);
        cells.clear();


        cells.push_back((PlotCellValue(6)));
        cells.push_back((PlotCellValue("Depth")));
        cells.push_back((PlotCellValue(minDepthValue_)));
        cells.push_back((PlotCellValue(maxDepthValue_)));
        cells.push_back((PlotCellValue(medianDepthValue_)));
        pData_->insert(cells);
        cells.clear();


        cells.push_back((PlotCellValue(7)));
        cells.push_back((PlotCellValue("Grey")));
        cells.push_back((PlotCellValue(minGreyValue_)));
        cells.push_back((PlotCellValue(maxGreyValue_)));
        cells.push_back((PlotCellValue(medianGreyValue_)));
        pData_->insert(cells);
        cells.clear();


        cells.push_back((PlotCellValue(8)));
        cells.push_back((PlotCellValue("C")));
        cells.push_back((PlotCellValue(minCValue_)));
        cells.push_back((PlotCellValue(maxCValue_)));
        cells.push_back((PlotCellValue(medianCValue_)));
        pData_->insert(cells);
        cells.clear();


        cells.push_back((PlotCellValue(9)));
        cells.push_back((PlotCellValue("Y")));
        cells.push_back((PlotCellValue(minYValue_)));
        cells.push_back((PlotCellValue(maxYValue_)));
        cells.push_back((PlotCellValue(medianYValue_)));
        pData_->insert(cells);
        cells.clear();


        cells.push_back((PlotCellValue(10)));
        cells.push_back((PlotCellValue("M")));
        cells.push_back((PlotCellValue(minMValue_)));
        cells.push_back((PlotCellValue(maxMValue_)));
        cells.push_back((PlotCellValue(medianMValue_)));
        pData_->insert(cells);
        cells.clear();


        cells.push_back((PlotCellValue(11)));
        cells.push_back((PlotCellValue("K")));
        cells.push_back((PlotCellValue(minKValue_)));
        cells.push_back((PlotCellValue(maxKValue_)));
        cells.push_back((PlotCellValue(medianKValue_)));
        pData_->insert(cells);
        cells.clear();


        cells.push_back((PlotCellValue(12)));
        cells.push_back((PlotCellValue("H")));
        if (allOutPutZeroOne_.get()) {
            cells.push_back((PlotCellValue(minHValue_/360.0f)));
            cells.push_back((PlotCellValue(maxHValue_/360.0f)));
        cells.push_back((PlotCellValue(medianHValue_)));
        }
        else {
            cells.push_back((PlotCellValue(minHValue_)));
            cells.push_back((PlotCellValue(maxHValue_)));
        cells.push_back((PlotCellValue(medianHValue_)));
        }
        pData_->insert(cells);
        cells.clear();


        cells.push_back((PlotCellValue(13)));
        cells.push_back((PlotCellValue("S")));
        cells.push_back((PlotCellValue(minSValue_)));
        cells.push_back((PlotCellValue(maxSValue_)));
        cells.push_back((PlotCellValue(medianSValue_)));
        pData_->insert(cells);
        cells.clear();


        cells.push_back((PlotCellValue(14)));
        cells.push_back((PlotCellValue("V")));
        cells.push_back((PlotCellValue(minVValue_)));
        cells.push_back((PlotCellValue(maxVValue_)));
        cells.push_back((PlotCellValue(medianVValue_)));
        pData_->insert(cells);
        cells.clear();
        Processor::setProgress(0.98);

        outport_.setData(pData_);
    }
    if (histogramOutPort_.isConnected()) {
        int j = toHistoRed_.get() + toHistoGreen_.get() + toHistoBlue_.get() + toHistoAlpha_.get() + toHistoGrey_.get()
            + toHistoLuminance_.get() + toHistoDepth_.get() + toHistoC_.get() + toHistoY_.get() + toHistoM_.get()
            + toHistoK_.get() + toHistoH_.get() + toHistoS_.get() + toHistoV_.get();
        std::vector<int> columns;
        columns.resize(j+1);
        columns[0] = 0;
        int k = 1;
        if (toHistoRed_.get()) {
            columns[k] = 1;
            ++k;
        }
        if (toHistoGreen_.get()) {
            columns[k] = 2;
            ++k;
        }
        if (toHistoBlue_.get()) {
            columns[k] = 3;
            ++k;
        }
        if (toHistoAlpha_.get()) {
            columns[k] = 4;
            ++k;
        }
        if (toHistoGrey_.get()) {
            columns[k] = 5;
            ++k;
        }
        if (toHistoLuminance_.get()) {
            columns[k] = 6;
            ++k;
        }
        if (toHistoDepth_.get()) {
            columns[k] = 7;
            ++k;
        }
        if (toHistoC_.get()) {
            columns[k] = 8;
            ++k;
        }
        if (toHistoY_.get()) {
            columns[k] = 9;
            ++k;
        }
        if (toHistoM_.get()) {
            columns[k] = 10;
            ++k;
        }
        if (toHistoK_.get()) {
            columns[k] = 11;
            ++k;
        }
        if (toHistoH_.get()) {
            columns[k] = 12;
            ++k;
        }
        if (toHistoS_.get()) {
            columns[k] = 13;
            ++k;
        }
        if (toHistoV_.get()) {
            columns[k] = 14;
            ++k;
        }
        PlotData* old;
        old = pDataHistogram_;
        if (k != pDataHistogramAll->getDataColumnCount()) {
            pDataHistogram_ = new PlotData(1,1);
            pDataHistogramAll->select(columns,1,j,*pDataHistogram_);
        }
        else {
            pDataHistogram_ = pDataHistogramAll;
        }
        Processor::setProgress(0.99);
        histogramOutPort_.setData(pDataHistogram_);
        delete old;
    }
    setPropertyValues();
    Processor::setProgress(1.f);
    if (pDataHistogram_ != pDataHistogramAll) {
        pDataHistogramAll->reset(0,0);
        delete pDataHistogramAll;
    }
}

void ImageAnalyzer::addToHistogram(std::vector<unsigned int>* values) {
    histogramRed_[values->at(0)]++;
    histogramGreen_[values->at(1)]++;
    histogramBlue_[values->at(2)]++;
    histogramAlpha_[values->at(3)]++;
    histogramGrey_[values->at(4)]++;
    histogramLuminance_[values->at(5)]++;
    histogramDepth_[values->at(6)]++;
    histogramC_[values->at(7)]++;
    histogramY_[values->at(8)]++;
    histogramM_[values->at(9)]++;
    histogramK_[values->at(10)]++;
    histogramH_[values->at(11)]++;
    histogramS_[values->at(12)]++;
    histogramV_[values->at(13)]++;
}

void ImageAnalyzer::setPropertyValues() {
    minRed_.set(minRedValue_);
    maxRed_.set(maxRedValue_);
    medianRed_.set(medianRedValue_);
    minGreen_.set(minGreenValue_);
    maxGreen_.set(maxGreenValue_);
    medianGreen_.set(medianGreenValue_);
    minBlue_.set(minBlueValue_);
    maxBlue_.set(maxBlueValue_);
    medianBlue_.set(medianBlueValue_);
    minAlpha_.set(minAlphaValue_);
    maxAlpha_.set(maxAlphaValue_);
    medianAlpha_.set(medianAlphaValue_);

    minGrey_.set(minGreyValue_);
    maxGrey_.set(maxGreyValue_);
    medianGrey_.set(medianGreyValue_);

    minLuminance_.set(minLuminanceValue_);
    maxLuminance_.set(maxLuminanceValue_);
    medianLuminance_.set(medianLuminanceValue_);

    minDepth_.set(minDepthValue_);
    maxDepth_.set(maxDepthValue_);
    medianDepth_.set(medianDepthValue_);

    //CYMK-Calculation
    minC_.set(minCValue_);
    maxC_.set(maxCValue_);
    medianC_.set(medianCValue_);
    minY_.set(minYValue_);
    maxY_.set(maxYValue_);
    medianY_.set(medianYValue_);
    minM_.set(minMValue_);
    maxM_.set(maxMValue_);
    medianM_.set(medianMValue_);
    minK_.set(minKValue_);
    maxK_.set(maxKValue_);
    medianK_.set(medianKValue_);

    //HSV-Calculation
    minH_.set(minHValue_);
    maxH_.set(maxHValue_);
    medianH_.set(medianHValue_);
    minS_.set(minSValue_);
    maxS_.set(maxSValue_);
    medianS_.set(medianSValue_);
    minV_.set(minVValue_);
    maxV_.set(maxVValue_);
    medianV_.set(medianVValue_);
}

tgt::Color ImageAnalyzer::getColorFromTexture(GLfloat* buffer, int position) {
    return tgt::Color(buffer[position],buffer[position+1],buffer[position+2],buffer[position+3]);
}

void ImageAnalyzer::calculate(){

    if (!isInitialized())
        return;

    if (!image_.hasRenderTarget())
        return;
    computeButton_.setReadOnlyFlag(true);
    image_.getColorTexture()->downloadTexture();

    int bytsPerPixel = image_.getColorTexture()->getBpt();
    bytePerPixel_.set(bytsPerPixel);

    GLfloat* textureBuffer = (GLfloat*)(image_.getColorTexture()->downloadTextureToBuffer(GL_RGBA,GL_FLOAT));//GL_RGBA16,GL_FLOAT


    unsigned int bitsright_ = (bytsPerPixel / 4)*8;
    unsigned int modValue_ = 1 << bitsright_;


    histogramRed_.clear();
    histogramGreen_.clear();
    histogramBlue_.clear();
    histogramAlpha_.clear();
    histogramGrey_.clear();
    histogramLuminance_.clear();
    histogramDepth_.clear();
    histogramC_.clear();
    histogramY_.clear();
    histogramM_.clear();
    histogramK_.clear();
    histogramH_.clear();
    histogramS_.clear();
    histogramV_.clear();

    histogramRed_.resize(histoSteps_);
    histogramGreen_.resize(histoSteps_);
    histogramBlue_.resize(histoSteps_);
    histogramAlpha_.resize(histoSteps_);
    histogramGrey_.resize(histoSteps_);
    histogramLuminance_.resize(histoSteps_);
    histogramDepth_.resize(histoSteps_);
    histogramC_.resize(histoSteps_);
    histogramY_.resize(histoSteps_);
    histogramM_.resize(histoSteps_);
    histogramK_.resize(histoSteps_);
    histogramH_.resize(histoSteps_);
    histogramS_.resize(histoSteps_);
    histogramV_.resize(histoSteps_);

    stepPerChannel_.set(modValue_);
    tgt::Color curColorVectorFloat = getColorFromTexture(textureBuffer,0);
    width_.set(image_.getColorTexture()->getWidth());
    height_.set(image_.getColorTexture()->getHeight());

    image_.getDepthTexture()->downloadTexture();
    float* fpixels = (float*)image_.getDepthTexture()->getCpuTextureData();

    float gray = greyScale(&curColorVectorFloat);
    float curLuminance = luminance(&curColorVectorFloat);
    std::vector<unsigned int> colorVectorInt;

    tgt::Vector4f curCYMKVector;
    cymkScale(&curColorVectorFloat,&curCYMKVector);

    tgt::Vector3f curHSVVector;
    hsvScale(&curColorVectorFloat,&curHSVVector);

    float curAlpha = curColorVectorFloat.a;
    float curRed = curColorVectorFloat.r;
    float curGreen = curColorVectorFloat.g;
    float curBlue = curColorVectorFloat.b;

    float curC = curCYMKVector.x;
    float curY = curCYMKVector.y;
    float curM = curCYMKVector.z;
    float curK = curCYMKVector.w;

    float curH = curHSVVector.x;
    float curS = curHSVVector.y;
    float curV = curHSVVector.z;


    minCValue_ = curC;
    maxCValue_ = curC;
    minYValue_ = curY;
    maxYValue_ = curY;
    minMValue_ = curM;
    maxMValue_ = curM;
    minKValue_ = curK;
    maxKValue_ = curK;

    minHValue_ = curH;
    maxHValue_ = curH;
    minSValue_ = curS;
    maxSValue_ = curS;
    minVValue_ = curV;
    maxVValue_ = curV;

    minLuminanceValue_ = curLuminance;
    maxLuminanceValue_ = curLuminance;

    minGreyValue_ = gray;
    maxGreyValue_ = gray;

    minAlphaValue_ = curAlpha;
    maxAlphaValue_ = curAlpha;

    float curDepth = *(fpixels);
    if (curDepth < 0)
        curDepth = 0;
    minDepthValue_ = curDepth;
    maxDepthValue_ = curDepth;


    minRedValue_ = curRed;
    maxRedValue_ = curRed;
    minGreenValue_ = curGreen;
    maxGreenValue_ = curGreen;
    minBlueValue_ = curBlue;
    maxBlueValue_ = curBlue;
//    buffer = textureBuffer;
    colorVectorInt.resize(14);
    for (int x = 0; x < image_.getSize().x; ++x) {
        Processor::setProgress(x*0.9f/(image_.getSize().x*1.f));
        for (int y = 0; y < image_.getSize().y; ++y) {

            if (((x*image_.getSize().y + y) % calculateStep_.get()) == 0){
                curDepth = *(++fpixels);
                if (curDepth < 0)
                    curDepth = 0;
                curColorVectorFloat = getColorFromTexture(textureBuffer, (x*image_.getSize().y + y)*4);//= image_.getColorTexture()->texelAsFloat(x,y);
                gray = greyScale(&curColorVectorFloat);
                curLuminance= luminance(&curColorVectorFloat);
                cymkScale(&curColorVectorFloat,&curCYMKVector);
                hsvScale(&curColorVectorFloat,&curHSVVector);

                curC = curCYMKVector.x;
                curY = curCYMKVector.y;
                curM = curCYMKVector.z;
                curK = curCYMKVector.w;
                curAlpha = curColorVectorFloat.a;
                curRed = curColorVectorFloat.r;
                curGreen = curColorVectorFloat.g;
                curBlue = curColorVectorFloat.b;
                curH = curHSVVector.x;
                curS = curHSVVector.y;
                curV = curHSVVector.z;

                colorVectorInt[0] = static_cast<unsigned int>(curColorVectorFloat.r*(histoSteps_ - 1));
                colorVectorInt[1] = static_cast<unsigned int>(curColorVectorFloat.g*(histoSteps_ - 1));
                colorVectorInt[2] = static_cast<unsigned int>(curColorVectorFloat.b*(histoSteps_ - 1));
                colorVectorInt[3] = static_cast<unsigned int>(curColorVectorFloat.a*(histoSteps_ - 1));
                colorVectorInt[4] = static_cast<unsigned int>(gray*(histoSteps_ - 1));
                colorVectorInt[5] = static_cast<unsigned int>(curLuminance*(histoSteps_ - 1));
                colorVectorInt[6] = static_cast<unsigned int>(curDepth*(histoSteps_ - 1));

                colorVectorInt[7] = static_cast<unsigned int>(curC*(histoSteps_ - 1));
                colorVectorInt[8] = static_cast<unsigned int>(curY*(histoSteps_ - 1));
                colorVectorInt[9] = static_cast<unsigned int>(curM*(histoSteps_ - 1));
                colorVectorInt[10] = static_cast<unsigned int>(curK*(histoSteps_ - 1));
                colorVectorInt[11] = static_cast<unsigned int>(curH/360.0f*(histoSteps_ - 1));
                colorVectorInt[12] = static_cast<unsigned int>(curS*(histoSteps_ - 1));
                colorVectorInt[13] = static_cast<unsigned int>(curV*(histoSteps_ - 1));
                addToHistogram(&colorVectorInt);


                if ((minDepthValue_ > curDepth))
                    minDepthValue_ = curDepth;
                else if ((maxDepthValue_ < curDepth))
                    maxDepthValue_ = curDepth;

/*
*       Colormodel    RGB     START
*/
                if ((minLuminanceValue_ > curLuminance))
                    minLuminanceValue_ = curLuminance;
                else if ((maxLuminanceValue_ < curLuminance))
                    maxLuminanceValue_ = curLuminance;

                if ((minGreyValue_ > gray))
                    minGreyValue_ = gray;
                else if ((maxGreyValue_ < gray))
                    maxGreyValue_ = gray;

                if ((minAlphaValue_ > curAlpha))
                    minAlphaValue_ = curAlpha;
                else if ((maxAlphaValue_ < curAlpha))
                    maxAlphaValue_ = curAlpha;

                if ((minRedValue_ > curRed))
                    minRedValue_ = curRed;
                else if ((maxRedValue_ < curRed))
                    maxRedValue_ = curRed;

                if ((minGreenValue_ > curGreen))
                    minGreenValue_ = curGreen;
                else if ((maxGreenValue_ < curGreen))
                    maxGreenValue_ = curGreen;

                if ((minBlueValue_ > curBlue))
                    minBlueValue_ = curBlue;
                else if ((maxBlueValue_ < curBlue))
                    maxBlueValue_ = curBlue;
/*
*       Colormodel      RGB     END
*
*/

/*
*       Colormodel      CMYK    START
*/
                if ((minCValue_ > curC))
                    minCValue_ = curC;
                else if ((maxCValue_ < curC))
                    maxCValue_ = curC;
                if ((minYValue_ > curY))
                    minYValue_ = curY;
                else if ((maxYValue_ < curY))
                    maxYValue_ = curY;
                if ((minMValue_ > curM))
                    minMValue_ = curM;
                else if ((maxMValue_ < curM))
                    maxMValue_ = curM;
                if ((minKValue_ > curK))
                    minKValue_ = curK;
                else if ((maxKValue_ < curK))
                    maxKValue_ = curK;
/*
*       Colormodel       CMYK     END
*
*/

/*
*       Colormodel      HSV    START
*/
                if ((minHValue_ > curH))
                    minHValue_ = curH;
                else if ((maxHValue_ < curH))
                    maxHValue_ = curH;
                if ((minSValue_ > curS))
                    minSValue_ = curS;
                else if ((maxSValue_ < curS))
                    maxSValue_ = curS;
                if ((minVValue_ > curV))
                    minVValue_ = curV;
                else if ((maxVValue_ < curV))
                    maxVValue_ = curV;
/*
*       Colormodel       HSV     END
*
*/
            }
        }
    }
    colorVectorInt.clear();
    image_.getColorTexture()->setCpuTextureData(0,false);
    image_.getDepthTexture()->setCpuTextureData(0,false);
    delete textureBuffer;
    setOutPortData();
}

} // namespace
