/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2017 University of Muenster, Germany.                        *
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

#ifndef VRN_EnsembleSimilarityPlot_H
#define VRN_EnsembleSimilarityPlot_H

#include "voreen/core/processors/renderprocessor.h"

#include "voreen/core/io/volumeserializerpopulator.h"

#include "voreen/core/ports/volumeport.h"

#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/progressproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/datastructures/volume/histogram.h"

#include "../datastructures/fieldplotdata.h"
#include "../ports/ensembledatasetport.h"
#include "../ports/fieldplotdataport.h"

namespace voreen {

/**
 *
 */

class VRN_CORE_API EnsembleSimilarityPlot : public RenderProcessor {

public:
    EnsembleSimilarityPlot();
    virtual ~EnsembleSimilarityPlot();
    virtual Processor* create() const;

    virtual std::string getClassName() const      { return "EnsembleSimilarityPlot"; }
    virtual std::string getCategory() const       { return "Plotting";               }
    virtual CodeState getCodeState() const        { return CODE_STATE_EXPERIMENTAL;  }

protected:
    virtual void onEvent(tgt::Event* e);
    virtual void process();
    static const std::string loggerCat_;

protected:

    EnsembleDatasetPort inport_;
    RenderPort outport_;
    ButtonProperty calculateButton_;
    FloatProperty isoValueProperty_;
    OptionProperty<std::string> channelProperty_;
    ProgressProperty progressProperty_;
    IntProperty selectedTimeStep_;
    BoolProperty autoCalculateProperty_;

    /*
     * "run1" => [
     *    0 => [
     *       "channel1": 100%,
     *       "channel2": 0%,
     *       ...
     *    ],
     *    [ ... ]
     * ]
     * "run2" => [ ... ],
     *    ...
     * */
    std::map<std::string, std::vector<std::map<std::string, float>>> similarityPlot_;
    int sampleRate_;
    std::pair<float, float> similarityRange_;
    float timeStepPosition_;
    float viewPortWidth_;
    int maxNumTimeSteps_;

    void initProperties();
    void initIsoValueProperty();
    void createSimilarityPlot();
    float calculateSimilarity(const VolumeBase* timestepVolume, float isoValue, float normalizationFactor);
    void mouseEvent(tgt::MouseEvent* e);
    void drawTimeStepSelection();
    void drawSimilarityTimeSeries();
    void onIsoValueChanged();
    void onInportChanged();
};

} // namespace

#endif // VRN_EnsembleSimilarityPlot_H
