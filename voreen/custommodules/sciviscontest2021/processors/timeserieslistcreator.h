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

#ifndef VRN_TIMESERIESLISTCREATOR_H
#define VRN_TIMESERIESLISTCREATOR_H

#include "voreen/core/processors/asynccomputeprocessor.h"
#include "modules/ensembleanalysis/ports/ensembledatasetport.h"
#include "custommodules/sciviscontest2021/ports/timeserieslistport.h"
#include "voreen/core/properties/string/stringlistproperty.h"

namespace voreen {

    struct TimeSeriesListCreatorInput {
        PortDataPointer<EnsembleDataset> ensemble;
        int numberOfSamples;
        float depth;
        std::vector<int> fieldIndices;
    };

    struct TimeSeriesListCreatorOutput {
        std::unique_ptr<TimeSeriesList> timeSeriesList;
    };

    /**
     * Creates an ensemble dataset with disk representations (VVD) from raw simulated data (local or cluster)
     * and adds measured data. The created ensemble can be used by the ensemble analysis module.
     */
    class VRN_CORE_API TimeSeriesListCreator : public AsyncComputeProcessor<TimeSeriesListCreatorInput, TimeSeriesListCreatorOutput> {
    public:
        TimeSeriesListCreator();
        virtual ~TimeSeriesListCreator();

        Processor* create() const;
        std::string getClassName() const { return "TimeSeriesListCreator"; }
        std::string getCategory() const { return "Time Series Processing"; }
        CodeState getCodeState() const { return CODE_STATE_EXPERIMENTAL; }

        virtual bool isReady() const;

        virtual ComputeInput prepareComputeInput();
        virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
        virtual void processComputeOutput(ComputeOutput output);

    protected:
        void setDescriptions() {
            setDescription(
                "Creates a time series list from an ensemble dataset."
            );
        }

    private:
        void adjustToEnsemble();

        EnsembleDatasetPort ensembleInport_;
        TimeSeriesListPort timeSeriesOutport_;
        VolumePort maskInport_;
        IntProperty numberOfSamples_;
        BoolProperty sampleOnSlice_;
        FloatProperty depth_;
        StringListProperty fields_;

        void onInportChange();

        static const std::string loggerCat_; ///< category used in logging
    };

}

#endif // VRN_TIMESERIESLISTCREATOR_H
