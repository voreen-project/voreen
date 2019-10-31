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

#if 0

#ifndef VRN_STREAMLINECREATORBACKGROUNDTHREAD_H
#define VRN_STREAMLINECREATORBACKGROUNDTHREAD_H

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/utils/backgroundthread.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"

#include "modules/flowreen/processors/streamline/streamlinecreator.h"
#include "modules/flowreen/datastructures/streamlinelist.h"

#include <random>

namespace voreen {

/**
 * Background thread used to calculate the streamlines.
 */
class VRN_CORE_API StreamlineCreatorBackgroundThread : public ProcessorBackgroundThread<StreamlineCreator> {
public:
    /** Constructor */
    StreamlineCreatorBackgroundThread(StreamlineCreator* processor,
            int seedTime,
            const VolumeBase* flow,
            const VolumeBase* seedMask,
            StreamlineList* output,
            int maxNumStreamlines,
            tgt::ivec2 streamlineLengthThreshold,
            tgt::vec2 absoluteMagnitudeThreshold,
            int stopIntegrationAngleTreshold,
            StreamlineCreator::FilterMode filterMode);

protected:
    /** Main-Function used to calculate the streamlines */
    virtual void threadMain();

    /** Calcualtes a streamline */
    Streamline computeStreamlineRungeKutta(const tgt::vec3& start);
    /** Returns the vlocity in the selected filter mode. */
    tgt::vec3 getVelocityAt(const tgt::vec3& pos);
private:

    StreamlineList* output_;               ///< output, which will be used by the processor

    std::function<float()> rnd;

    const VolumeBase* flow_;               ///< input flow
    VolumeRAMRepresentationLock flowRepresentation_; ///< representation lock for input flow
    const VolumeBase* seedMask_;           ///< seed mask
    std::unique_ptr<VolumeRAMRepresentationLock> seedMaskRepresentation_; ///< representation lock for input flow

    size_t maxNumStreamlines_;             ///< maximal number of streamlines
    tgt::svec2 streamlineLengthThreshold_; ///< streamline length must be in this interval
    tgt::vec2 absoluteMagnitudeThreshold_; ///< only magnitudes in this interval are used
    float stopIntegrationAngleThreshold_;    ///< stop integration if angle is exceeded
    StreamlineCreator::FilterMode filterMode_; ///< filtering inside the volume
};

}   // namespace

#endif  // VRN_STREAMLINECREATORBACKGROUNDTHREAD_H

#endif