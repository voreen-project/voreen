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

#ifndef VRN_ENSEMBLEDATASET_H
#define VRN_ENSEMBLEDATASET_H

#include "voreen/core/voreencoreapi.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumelist.h"
#include "voreen/core/io/serialization/serializable.h"
#include "voreen/core/io/serialization/xmlserializer.h"
#include "voreen/core/io/serialization/xmldeserializer.h"
#include "voreen/core/ports/port.h"

#include "tgt/vector.h"

#include <map>

namespace voreen {

/**
 * Datastructure used to represent the structure of an ensemble dataset.
 */
class VRN_CORE_API EnsembleDataset : public DataInvalidationObservable {
public:

    struct TimeStep {
        std::string path_;
        float time_;
        float duration_;
        std::map<std::string, const VolumeBase*> channels_;
    };

    struct Run {
        std::string name_;
        tgt::vec3 color_;
        std::vector<TimeStep> timeSteps_;
    };

    /** Constructor */
    explicit EnsembleDataset();
    explicit EnsembleDataset(const EnsembleDataset& origin);
    explicit EnsembleDataset(const EnsembleDataset* const origin);
    /** Destructor */
    ~EnsembleDataset();

public:

    //----------------
    //  Access
    //----------------
    void addRun(const Run& run);
    const std::vector<Run>& getRuns() const;

    size_t getMinNumTimeSteps() const;
    size_t getMaxNumTimeSteps() const;
    size_t getTotalNumTimeSteps() const;

    float getMinTimeStepDuration() const;
    float getMaxTimeStepDuration() const;
    float getStartTime() const;
    float getEndTime() const;
    float getMaxTotalDuration() const;
    const tgt::vec2& getCommonTimeInterval() const;

    const tgt::vec2& getValueRange(const std::string& channel) const;

    const tgt::svec3& getDimensions() const;
    const tgt::vec3& getSpacing() const;

    const tgt::IntBounds& getRoi() const;
    void setRoi(const tgt::IntBounds& roi);

    const std::vector<std::string>& getCommonChannels() const;
    std::vector<const VolumeBase*> getVolumes() const;

    /**
     * This function takes a sample in voxel space
     *
     * @param volume
     * @param spacing
     * @param sample
     * @param filter
     * @return
     */
    float pickSample(const VolumeRAM_Float* volume, const tgt::vec3& spacing, tgt::vec3 sample, VolumeRAM::Filter filter = VolumeRAM::LINEAR) const;

    size_t pickTimeStep(size_t runIdx, float time) const;

private:

    //----------------
    //  Members
    //----------------
    std::vector<Run> runs_;
    std::vector<std::string> commonChannels_;

    std::map<std::string, tgt::vec2> valueRange_;

    size_t minNumTimeSteps_;
    size_t maxNumTimeSteps_;
    size_t totalNumTimeSteps_;

    float minTimeStepDuration_;
    float maxTimeStepDuration_;
    float startTime_;
    float endTime_;
    tgt::vec2 commonTimeInterval_;

    tgt::svec3 dimensions_;
    tgt::vec3 spacing_;

    tgt::IntBounds roi_;

};

}   // namespace

#endif //VRN_ENSEMBLEDATASET_H
