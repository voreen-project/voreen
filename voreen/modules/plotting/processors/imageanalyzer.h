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

#ifndef VRN_IMAGEANALYZER_H
#define VRN_IMAGEANALYZER_H

#include "voreen/core/processors/imageprocessor.h"

#include "../ports/plotport.h"
#include "voreen/core/ports/renderport.h"

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/floatproperty.h"

#include "../datastructures/plotcell.h"

namespace voreen {

class VRN_CORE_API ImageAnalyzer : public ImageProcessor {
public:
    ImageAnalyzer();
    virtual Processor* create() const;

    virtual std::string getCategory() const  { return "Utility"; }
    virtual std::string getClassName() const { return "ImageAnalyzer"; }
    virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL; }
    virtual bool isUtility() const           { return true; }
    //virtual bool usesExpensiveComputation() const { return true; }
    virtual bool isEndProcessor() const;
    virtual bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("Shows some information about the input image.");
    }

    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

    RenderPort image_;
    RenderPort inportMask_;
    PlotPort outport_;
    PlotPort histogramOutPort_;

    ButtonProperty computeButton_;
    BoolProperty allOutPutZeroOne_;

    IntOptionProperty interAction_;

    IntProperty calculateStep_;
    IntProperty histogramSteps_;
    IntProperty bytePerPixel_;
    IntProperty width_;
    IntProperty height_;
    IntProperty stepPerChannel_;

    //RGB - Calculations
    FloatProperty minRed_;
    FloatProperty maxRed_;
    FloatProperty medianRed_;
    FloatProperty minGreen_;
    FloatProperty maxGreen_;
    FloatProperty medianGreen_;
    FloatProperty minBlue_;
    FloatProperty maxBlue_;
    FloatProperty medianBlue_;
    FloatProperty minAlpha_;
    FloatProperty maxAlpha_;
    FloatProperty medianAlpha_;

    FloatProperty minGrey_;
    FloatProperty maxGrey_;
    FloatProperty medianGrey_;

    FloatProperty minLuminance_;
    FloatProperty maxLuminance_;
    FloatProperty medianLuminance_;

    FloatProperty minDepth_;
    FloatProperty maxDepth_;
    FloatProperty medianDepth_;

    //CYMK-Calculation
    FloatProperty minC_;
    FloatProperty maxC_;
    FloatProperty medianC_;
    FloatProperty minY_;
    FloatProperty maxY_;
    FloatProperty medianY_;
    FloatProperty minM_;
    FloatProperty maxM_;
    FloatProperty medianM_;
    FloatProperty minK_;
    FloatProperty maxK_;
    FloatProperty medianK_;

    //HSV-Calculation
    FloatProperty minH_;
    FloatProperty maxH_;
    FloatProperty medianH_;
    FloatProperty minS_;
    FloatProperty maxS_;
    FloatProperty medianS_;
    FloatProperty minV_;
    FloatProperty maxV_;
    FloatProperty medianV_;

    BoolProperty toHistoRed_;
    BoolProperty toHistoGreen_;
    BoolProperty toHistoBlue_;
    BoolProperty toHistoAlpha_;
    BoolProperty toHistoGrey_;
    BoolProperty toHistoLuminance_;
    BoolProperty toHistoDepth_;

    BoolProperty toHistoC_;
    BoolProperty toHistoY_;
    BoolProperty toHistoM_;
    BoolProperty toHistoK_;

    BoolProperty toHistoH_;
    BoolProperty toHistoS_;
    BoolProperty toHistoV_;

private:
    float greyScale(tgt::Color* colorVector);
    void colorScale(tgt::Color* colorValue, unsigned int modValue, tgt::Vector4<unsigned int>* result);
    float luminance(tgt::Color* colorVector);
    void cymkScale(tgt::Color* colorValue, tgt::Vector4f* result);
    void hsvScale(tgt::Color* colorValue,tgt::Vector3f* result);

    tgt::Color getColorFromTexture(GLfloat* buffer, int position);

    void setOutPortData();
    void calculate();
    void addToHistogram(std::vector<unsigned int>* values);
    void activateComputeButton();
    void setPropertyValues();

    PlotData* pData_;
    PlotData* pDataHistogram_;

    float minRedValue_;
    float maxRedValue_;
    float medianRedValue_;
    float minGreenValue_;
    float maxGreenValue_;
    float medianGreenValue_;
    float minBlueValue_;
    float maxBlueValue_;
    float medianBlueValue_;
    float minAlphaValue_;
    float maxAlphaValue_;
    float medianAlphaValue_;

    float minGreyValue_;
    float maxGreyValue_;
    float medianGreyValue_;

    float minLuminanceValue_;
    float maxLuminanceValue_;
    float medianLuminanceValue_;

    float minDepthValue_;
    float maxDepthValue_;
    float medianDepthValue_;

    //CYMK-Calculation
    float minCValue_;
    float maxCValue_;
    float medianCValue_;
    float minYValue_;
    float maxYValue_;
    float medianYValue_;
    float minMValue_;
    float maxMValue_;
    float medianMValue_;
    float minKValue_;
    float maxKValue_;
    float medianKValue_;

    //HSV-Calculation
    float minHValue_;
    float maxHValue_;
    float medianHValue_;
    float minSValue_;
    float maxSValue_;
    float medianSValue_;
    float minVValue_;
    float maxVValue_;
    float medianVValue_;

    int histoSteps_;

    std::vector<unsigned int> histogramRed_;
    std::vector<unsigned int> histogramGreen_;
    std::vector<unsigned int> histogramBlue_;
    std::vector<unsigned int> histogramAlpha_;
    std::vector<unsigned int> histogramGrey_;
    std::vector<unsigned int> histogramLuminance_;
    std::vector<unsigned int> histogramDepth_;

    std::vector<unsigned int> histogramC_;
    std::vector<unsigned int> histogramY_;
    std::vector<unsigned int> histogramM_;
    std::vector<unsigned int> histogramK_;

    std::vector<unsigned int> histogramH_;
    std::vector<unsigned int> histogramS_;
    std::vector<unsigned int> histogramV_;

    static const std::string loggerCat_;

};

}

#endif //VRN_IMAGEANALYZER_H
