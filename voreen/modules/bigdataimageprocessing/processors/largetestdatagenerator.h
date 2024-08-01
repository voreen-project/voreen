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

#ifndef VRN_LARGETESTDATAGENERATOR_H
#define VRN_LARGETESTDATAGENERATOR_H

#include "voreen/core/processors/asynccomputeprocessor.h"

#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/progressproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/temppathproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"

#include "modules/hdf5/io/hdf5filevolume.h"

#include <random>


namespace voreen {

struct LargeTestDataGeneratorInput {
    enum Scenario {
        CELLS,
        VESSELS,
    };
    enum NoiseType {
        GAUSSIAN,
        POISSON,
    };
    enum LargeTestDataInvalidation {
        All,
        Seeds,
    };
    typedef std::mt19937 random_engine_type;

    Scenario scenario;
    LargeTestDataInvalidation invalidation;
    std::unique_ptr<HDF5FileVolume> outputVolumeNoisy;
    std::unique_ptr<HDF5FileVolume> outputVolumeGT;
    random_engine_type randomEngine;
    tgt::ivec3 dimensions;
    uint16_t foregroundMean;
    uint16_t backgroundMean;
    NoiseType noiseType;
    float gaussianNoiseSD;
    float density;
    tgt::ivec2 structureSizeRange;
    bool retainLabel;

    LargeTestDataGeneratorInput(
            Scenario scenario
            , LargeTestDataInvalidation invalidation
            , std::unique_ptr<HDF5FileVolume>&& outputVolumeNoisy
            , std::unique_ptr<HDF5FileVolume>&& outputVolumeGT
            , random_engine_type randomEngine
            , tgt::ivec3 dimensions
            , uint16_t foregroundMean
            , uint16_t backgroundMean
            , NoiseType noiseType
            , float gaussianNoiseSD
            , float density
            , tgt::ivec2 structureSizeRange
            , bool retainLabel
            )
        : scenario(scenario)
        , invalidation(invalidation)
        , outputVolumeNoisy(std::move(outputVolumeNoisy))
        , outputVolumeGT(std::move(outputVolumeGT))
        , randomEngine(randomEngine)
        , dimensions(dimensions)
        , foregroundMean(foregroundMean)
        , backgroundMean(backgroundMean)
        , noiseType(noiseType)
        , gaussianNoiseSD(gaussianNoiseSD)
        , density(density)
        , structureSizeRange(structureSizeRange)
        , retainLabel(retainLabel)
    {
    }

    LargeTestDataGeneratorInput(const LargeTestDataGeneratorInput&) = delete;
    LargeTestDataGeneratorInput(LargeTestDataGeneratorInput&& old)
        : scenario(old.scenario)
        , invalidation(old.invalidation)
        , outputVolumeNoisy(std::move(old.outputVolumeNoisy))
        , outputVolumeGT(std::move(old.outputVolumeGT))
        , randomEngine(old.randomEngine)
        , dimensions(old.dimensions)
        , foregroundMean(old.foregroundMean)
        , backgroundMean(old.backgroundMean)
        , noiseType(old.noiseType)
        , gaussianNoiseSD(old.gaussianNoiseSD)
        , density(old.density)
        , structureSizeRange(old.structureSizeRange)
        , retainLabel(old.retainLabel)
    {
    }
};
struct LargeTestDataGeneratorOutput {
    std::unique_ptr<PointSegmentListGeometryVec3> foregroundLabels;
    std::unique_ptr<PointSegmentListGeometryVec3> backgroundLabels;
    std::string outputVolumeNoisyFilePath;
    std::string outputVolumeGTFilePath;
};

class LargeTestDataGenerator : public AsyncComputeProcessor<LargeTestDataGeneratorInput, LargeTestDataGeneratorOutput> {
public:
    LargeTestDataGenerator();
    virtual ~LargeTestDataGenerator();

    virtual std::string getClassName() const         { return "LargeTestDataGenerator"; }
    virtual std::string getCategory() const       { return "Volume Processing"; }
    virtual VoreenSerializableObject* create() const;
    virtual void setDescriptions() {
        setDescription("Processor for generating test data potentially larger than RAM capacity.");
    }
    virtual CodeState getCodeState() const        { return CODE_STATE_EXPERIMENTAL;   }
    virtual bool isReady() const;

    virtual LargeTestDataGeneratorInput prepareComputeInput();
    virtual LargeTestDataGeneratorOutput compute(LargeTestDataGeneratorInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(LargeTestDataGeneratorOutput output);

    void invalidateAll();
    void invalidateSeeds();

    static const std::string loggerCat_;
private:
    // Ports
    VolumePort outportNoisy_;
    VolumePort outportGT_;
    GeometryPort foregroundLabelsPort_;
    GeometryPort backgroundLabelsPort_;

    // General properties
    TempPathProperty outputVolumeNoisyFilePath_;
    TempPathProperty outputVolumeGTFilePath_;
    IntProperty foregroundMean_;
    IntProperty backgroundMean_;
    OptionProperty<LargeTestDataGeneratorInput::NoiseType> noiseType_;
    FloatProperty gaussianNoiseSD_;
    FloatProperty density_;
    IntProperty seed_;
    IntVec3Property volumeDimensions_;
    IntIntervalProperty structureSizeRange_;
    OptionProperty<LargeTestDataGeneratorInput::Scenario> scenario_;
    BoolProperty retainLabel_;
    boost::optional<LargeTestDataGeneratorInput::LargeTestDataInvalidation> invalidation_;

    std::random_device randomDevice;
};
} // namespace voreen

#endif // VRN_LARGETESTDATAGENERATOR_H
