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

#ifndef VRN_ENSEMBLEDATASET_H
#define VRN_ENSEMBLEDATASET_H

#include "voreen/core/voreencoreapi.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/volume/volumelist.h"
#include "voreen/core/utils/statistics.h"

#include "tgt/vector.h"

#include <map>

namespace voreen {

/**
 * Datastructure used to represent the structure of an ensemble dataset.
 */
class VRN_CORE_API EnsembleDataset : public DataInvalidationObservable {
public:

    /**
     * This struct defines a single time step of a run.
     * It also holds the actual volume data.
     */
    struct TimeStep {
        std::string path_; ///< The time step's path containing the volume data.
        float time_;       ///< The point in time of this time step.
        float duration_;   ///< the duration of the time step.
        std::map<std::string, const VolumeBase*> fieldNames_; ///< Field names mapped to their volume data.
    };

    /**
     * This struct defines a run which is a unique member of the ensemble.
     */
    struct Run {
        std::string name_; ///< The run's name.
        tgt::vec3 color_;  ///< The run's distinct color.
        std::vector<TimeStep> timeSteps_; ///< List of time steps.
    };

    /** Constructor */
    EnsembleDataset();
    EnsembleDataset(const EnsembleDataset& origin);

public:

    /**
     * Add a new run to the ensemble.
     * This will update all ensemble meta data.
     */
    void addRun(const Run& run);

    /**
     * Returns all runs contained by the ensemble.
     */
    const std::vector<Run>& getRuns() const;

    /**
     * Returns the minimum number of time steps of all runs.
     */
    size_t getMinNumTimeSteps() const;

    /**
     * Returns the maximum number of time steps of all runs.
     */
    size_t getMaxNumTimeSteps() const;

    /**
     * Returns the total number of time steps of all runs combined.
     */
    size_t getTotalNumTimeSteps() const;

    /**
     * Returns the duration statistics for the specified run.
     */
    const Statistics& getTimeStepDurationStats(size_t runIdx) const;

    /**
     * Returns the color of the specified run.
     */
    const tgt::vec3& getColor(size_t runIdx) const;

    /**
     * Returns the minimum time step duration of all runs.
     */
    float getMinTimeStepDuration() const;

    /**
     * Returns the maximum time step duration of all runs.
     */
    float getMaxTimeStepDuration() const;
    
    /**
     * Returns the start time of the ensemble.
     * This is defined by the very first time step of all runs.
     */
    float getStartTime() const;

    /**
     * Returns the end time of the ensemble.
     * This is defined by the very last time step of all runs.
     */
    float getEndTime() const;

    /**
     * Returns the maximum total duration of the ensemble.
     * This is defined by getEndTime() - getStartTime().
     */
    float getMaxTotalDuration() const;

    /**
     * Returns the common time interval of all runs.
     * If non exist, tgt::vec2::zero is returned.
     */
    const tgt::vec2& getCommonTimeInterval() const;

    /**
     * Returns the scalar value range of the specified field.
     * @see VolumeMinMax
     */
    const tgt::vec2& getValueRange(const std::string& field) const;

    /**
     * Returns the magnitude value range of the specified field.
     * @see VolumeMinMaxMagnitude
     */
    const tgt::vec2& getMagnitudeRange(const std::string& field) const;

    /**
     * Returns the number of channels used by the specified field.
     */
    size_t getNumChannels(const std::string& field) const;

    /**
     * Returns the base type of the specified field.
     */
    const std::string& getBaseType(const std::string& field) const;

    /**
     * Returns the enclosing bounds of all time steps of all runs in physical coordinates.
     */
    const tgt::Bounds& getBounds() const;

    /**
     * Returns the bounds enclosing only the common space of all time steps of all runs in physical coordinates.
     * If there is no common space, the bounds are undefined.
     */
    const tgt::Bounds& getCommonBounds() const;

    /**
     * Returns the region of interest defined for the ensemble in physical coordinates.
     * Currently, only a single region defined by a bounding box is supported.
     */
    const tgt::Bounds& getRoi() const;

    /**
     * Sets the region of interest for the ensemble.
     * The region must be defined and intersect with the common bounds.
     * Otherwise, the region will not be set.
     */
    void setRoi(tgt::Bounds roi);

    /**
     * Returns a list of all unique field names of all runs.
     */
    const std::vector<std::string>& getUniqueFieldNames() const;

    /**
     * Returns a list of all common field names of all runs.
     */
    const std::vector<std::string>& getCommonFieldNames() const;

    /**
     * Returns all volumes contained in the ensemble. Order is not defined.
     * Hence, a filtered ensemble should be used before calling this function.
     * @see EnsembleFilter
     * @see EnsembleVolumeExtractor
     */
    std::vector<const VolumeBase*> getVolumes() const;

    /**
     * This utility function returns a time step index corresponding to the specified time of the specified run.
     * If time is behind the end of the specified run, the last time step index is returned.
     */
    size_t pickTimeStep(size_t runIdx, float time) const;

    /**
     * Creates a HTML string that provides a brief overview about the whole ensemble, including its parameters.
     * Parameters currently need to be stored as MetaData the key of which needs to contain the keyword 'Parameter'
     * in order to be listed inside the HTML file.
     * The file can be viewed using any browser.
     */
    std::string toHTML() const;

private:

    struct RunMetaData {
        //tgt::vec3 color_; // currently member of run but should be stored here
        Statistics timeStepDurationStats_{false};
        // could add simulation parameters here
    };

    struct FieldMetaData {
        tgt::vec2 valueRange_;
        tgt::vec2 magnitudeRange_;
        size_t numChannels_;
        std::string baseType_;
    };

    //----------------
    //  Members
    //----------------
    std::vector<Run> runs_;
    std::vector<std::string> uniqueFieldNames_;
    std::vector<std::string> commonFieldNames_;

    std::vector<RunMetaData> runMetaData_;
    std::map<std::string, FieldMetaData> fieldMetaData_;
    std::set<std::string> allParameters_;

    size_t minNumTimeSteps_;
    size_t maxNumTimeSteps_;
    size_t totalNumTimeSteps_;

    float minTimeStepDuration_;
    float maxTimeStepDuration_;
    float startTime_;
    float endTime_;
    tgt::vec2 commonTimeInterval_;

    tgt::Bounds bounds_;
    tgt::Bounds commonBounds_;
    tgt::Bounds roi_;

};

}   // namespace

#endif //VRN_ENSEMBLEDATASET_H
