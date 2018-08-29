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

#ifndef VRN_CONNECTEDCOMPONENTANALYSIS_H
#define VRN_CONNECTEDCOMPONENTANALYSIS_H

#include "voreen/core/processors/asynccomputeprocessor.h"

#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/temppathproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/progressproperty.h"
#include "voreen/core/ports/volumeport.h"

#include "../algorithm/streamingcomponents.h"

#include "modules/hdf5/io/hdf5volumereader.h"
#include "modules/hdf5/io/hdf5filevolume.h"
#include "../util/csvwriter.h"

#include <functional>

namespace voreen {

class CCANodeMetaData {
public:
    CCANodeMetaData();
    CCANodeMetaData(tgt::svec2 yzPo, size_t lowerBound, size_t upperBound);
    ~CCANodeMetaData() {}
    CCANodeMetaData& operator+=(const CCANodeMetaData& rhs);

    size_t volume_;
    tgt::TemplateBounds<size_t> bounds_;
};

typedef CSVWriter<uint32_t, size_t, float, float, float, float, float, float> CCAWriterType;

enum CCANeighbourhoodMode {
    N_6 = 2,
    N_18 = 1,
    N_26 = 0,
};

struct CCAComputeInput {
    std::unique_ptr<CCAWriterType> statWriter;
    std::function<void(int uint32_t, const CCANodeMetaData&)> writeMetaData;
    const VolumeBase& inputVolume;
    std::unique_ptr<HDF5FileVolume> outputVolume;
    CCANeighbourhoodMode neighbourhoodMode;


    CCAComputeInput(std::unique_ptr<CCAWriterType>&& pStatWriter, std::function<void(int id, const CCANodeMetaData&)> pWriteMetaData, const VolumeBase& pInputVolume,
            std::unique_ptr<HDF5FileVolume>&& pOutputVolume, CCANeighbourhoodMode pNeighbourhoodMode)
        : statWriter(std::move(pStatWriter))
        , writeMetaData(pWriteMetaData)
        , inputVolume(pInputVolume)
        , outputVolume(std::move(pOutputVolume))
        , neighbourhoodMode(pNeighbourhoodMode)
    {}

    CCAComputeInput(const CCAComputeInput&) = delete;
    CCAComputeInput(CCAComputeInput&& old)
        : statWriter(old.statWriter.release())
        , writeMetaData(old.writeMetaData)
        , inputVolume(old.inputVolume)
        , outputVolume(old.outputVolume.release())
        , neighbourhoodMode(old.neighbourhoodMode)
    {}
};

struct CCAComputeOutput {
    StreamingComponentsStats stats;
};


class ConnectedComponentAnalysis : public AsyncComputeProcessor<CCAComputeInput, CCAComputeOutput> {
public:
    ConnectedComponentAnalysis();
    virtual ~ConnectedComponentAnalysis();

    virtual std::string getClassName() const         { return "ConnectedComponentAnalysis";      }
    virtual std::string getCategory() const       { return "Volume Processing"; }
    virtual bool isEndProcessor() const       { return true; }
    virtual VoreenSerializableObject* create() const;
    virtual void setDescriptions() {
        setDescription("Processor that performs connected component analysis on a binary input volume");
    }
    virtual CodeState getCodeState() const        { return CODE_STATE_EXPERIMENTAL;   }

    virtual bool isReady() const;

    virtual CCAComputeInput prepareComputeInput();
    virtual CCAComputeOutput compute(CCAComputeInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(CCAComputeOutput output);

protected:
    std::function<bool(const CCANodeMetaData&)> generateComponentConstraintTest(const VolumeBase& volume) const;
    virtual void adjustPropertiesToInput();

    template<int ADJACENCY>
    StreamingComponentsStats runCCA(const VolumeBase& input, HDF5FileVolume& output, std::function<void(uint32_t id, const CCANodeMetaData&)> writeMetaData, ProgressReporter& progressReporter) const;

private:
    // Ports
    VolumePort inport_;
    VolumePort outport_;

    // General properties
    TempPathProperty outputVolumeFilePath_;
    BoolProperty writeComponentStatFile_;
    FileDialogProperty componentStatFilePath_;
    IntProperty outputVolumeDeflateLevel_;
    OptionProperty<CCANeighbourhoodMode> neighbourhoodMode_;
    BoolProperty invertBinarization_;
    FloatProperty  binarizationThreshold_;
    FloatProperty minBoundsDiagonal_;
    FloatProperty minBoundsDiagonalRelative_;
    IntProperty minVoxelVolume_;
    BoolProperty applyLabeling_;

    static const std::string loggerCat_;
};

template<int ADJACENCY>
StreamingComponentsStats ConnectedComponentAnalysis::runCCA(const VolumeBase& input, HDF5FileVolume& output, std::function<void(uint32_t id, const CCANodeMetaData&)> writeMetaData, ProgressReporter& progressReporter) const {
    typedef StreamingComponents<ADJACENCY, CCANodeMetaData, CCAVoidLabel> SC;
    SC sc;
    typename SC::getClassFunc getClass;

    float binarizationThresholdNormalized;
    if(input.hasMetaData("RealWorldMapping")) {
        // If the input volume does not have a RealWorldMapping we need to convert the binarizationThreshold to a normalized value.
        binarizationThresholdNormalized = input.getRealWorldMapping().realWorldToNormalized(binarizationThreshold_.get());
    } else {
        // If the input volume does not have a RealWorldMapping we expect RW values to be normalized.
        binarizationThresholdNormalized = binarizationThreshold_.get();
    }

    if(invertBinarization_.get()) {
        getClass = [binarizationThresholdNormalized](const VolumeRAM& slice, tgt::svec3 pos) {
            return slice.getVoxelNormalized(pos) <= binarizationThresholdNormalized ? CCAVoidLabel::some() : boost::none;
            };
    } else {
        getClass = [binarizationThresholdNormalized](const VolumeRAM& slice, tgt::svec3 pos) {
            return slice.getVoxelNormalized(pos) > binarizationThresholdNormalized ? CCAVoidLabel::some() : boost::none;
            };
    }

    return sc.cca(input, output, writeMetaData, getClass, applyLabeling_.get(), generateComponentConstraintTest(input), progressReporter);
}

} // namespace voreen

#endif // VRN_CONNECTEDCOMPONENTANALYSIS_H
