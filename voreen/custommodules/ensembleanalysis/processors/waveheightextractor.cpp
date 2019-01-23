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

#include "waveheightextractor.h"

#include "tgt/immediatemode/immediatemode.h"
#include "tgt/textureunit.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/utils/glsl.h"

#include "../datastructures/ensembledataset.h"

namespace voreen {

    const std::string WaveHeightExtractor::loggerCat_("voreen.viscontest2018.WaveHeightExtractor");

    WaveHeightExtractor::WaveHeightExtractor()
            : RenderProcessor(), inport_(Port::INPORT, "ensemble.inport", "Volume Input", false),
              volumeOutport_(Port::OUTPORT, "volumeout", "Volume Output", false),
              seaLevel_("seaLevel", "Sea Level", 5000.0f, 0.0f, 28000.0f),
              autoCompute_("autoCompute", "Auto-Compute", false),
              threshold_("threshold", "Iso-Value Threshold"), stepSize_("stepSize", "Step Size", 1.0f),
              recalculateButton_("calculate", "Calculate"), progressProperty_("progress", "Progress"),
              maxWaveHeight_("maxWaveHeight", "Max wave height") {
        addPort(inport_);
        addPort(volumeOutport_);
        addProperty(recalculateButton_);
        addProperty(progressProperty_);
        addProperty(seaLevel_);
        addProperty(threshold_);
        addProperty(autoCompute_);
        addProperty(stepSize_);
        addProperty(maxWaveHeight_);

        addProgressBar(&progressProperty_);

        ON_CHANGE(inport_, WaveHeightExtractor, onInportChanged);
        ON_CHANGE(recalculateButton_, WaveHeightExtractor, calculateWaveHeightImage);

        initProperties();

    }

    WaveHeightExtractor::~WaveHeightExtractor() {
    }

    void WaveHeightExtractor::initProperties() {
        stepSize_.setMinValue(0.0f);
        stepSize_.setMaxValue(10.0f);
        stepSize_.set(1.0f);

        seaLevel_.setMinValue(0.0f);
        seaLevel_.setMaxValue(10000.0f);
        seaLevel_.set(5000.0f);

        maxWaveHeight_.setMinValue(0);
        maxWaveHeight_.setMaxValue(10000);
        maxWaveHeight_.set(1000);

        threshold_.set(0.95f);
    }

    Processor *WaveHeightExtractor::create() const {
        return new WaveHeightExtractor();
    }

    void WaveHeightExtractor::process() {

        if (!inport_.hasData()) return;

    }

    void WaveHeightExtractor::calculateWaveHeightImage() {

        if (!inport_.hasData()) return;

        VolumeRAM_Float *output = new VolumeRAM_Float(tgt::vec3(inport_.getData()->getDimensions().x, inport_.getData()->getDimensions().y, 1));

        setProgress(0.0f);

        const VolumeBase* volumeBase = inport_.getData();
        const VolumeRAM_Float *volumeRepresentation = dynamic_cast<const VolumeRAM_Float *>(volumeBase->getRepresentation<VolumeRAM>());
        tgt::vec3 dimensions = volumeRepresentation->getDimensions();
        tgt::vec3 spacing = volumeBase->getSpacing();
        float totalSimulationHeight = static_cast<float>(spacing.y) * static_cast<float>(dimensions.y);
        float seaLevelMeters = seaLevel_.get();
        float meterStepDim = static_cast<float>(volumeRepresentation->getDimensions().y) / totalSimulationHeight;
        float seaLevelDim = seaLevelMeters * meterStepDim / 100.0f;
        float threshold = threshold_.get();


//#ifdef VRN_MODULE_OPENMP
//#pragma omp parallel for
//#endif
        for (size_t z = 0; z < static_cast<size_t>(volumeRepresentation->getDimensions().z); z++) {
            for (size_t x = 0; x < static_cast<size_t>(volumeRepresentation->getDimensions().x); x++) {

                int waveHeight = 0;
                for (int i = 0; i < maxWaveHeight_.get(); i++) {
                    float y = seaLevelDim + i * stepSize_.get() * meterStepDim;
                    if(y >= dimensions.y) break;
                    if (volumeRepresentation->getVoxelLinear(tgt::vec3(x,y,z)) >= threshold) waveHeight = i;
                    else break;
                }

                output->voxel(x, z, 0) = waveHeight * stepSize_.get() * meterStepDim;
            }
        }

        Volume *volume = new Volume(output, volumeBase->getSpacing(), tgt::vec3::zero);
        volumeOutport_.setData(volume, true);
    }

    void WaveHeightExtractor::onInportChanged() {
        if (!inport_.hasData() || !autoCompute_.get()) return;

        calculateWaveHeightImage();
    }

} // namespace voreen
