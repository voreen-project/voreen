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

#ifndef VRN_PATHLINECREATORBACKGROUNDTHREAD_H
#define VRN_PATHLINECREATORBACKGROUNDTHREAD_H

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/utils/backgroundthread.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"

#include "modules/flowreen/processors/streamline/pathlinecreator.h"
#include "modules/flowreen/datastructures/streamlinelist.h"

#include <random>

namespace voreen {

    /**
    * Background thread used to calculate the path lines.
    */
    class VRN_CORE_API PathlineCreatorBackgroundThread : public ProcessorBackgroundThread<PathlineCreator> {
    public:
        /** Constructor */
        PathlineCreatorBackgroundThread(PathlineCreator* processor, int seedTime, const VolumeList* flow, StreamlineList* output,
            int maxNumPathlines, tgt::ivec2 streamlineLengthThreshold, tgt::vec2 absoluteMagnitudeThreshold,
            PathlineCreator::FilterMode filterMode);
        /** Destructor */
        virtual ~PathlineCreatorBackgroundThread();
    protected:
        /** Main-Function used to calculate the streamlines */
        virtual void threadMain();

        //-----------------
        //  Helpers
        //-----------------
        /** Used to redefine a seed point. */
        void reseedPosition(tgt::vec3* seedingPositions, size_t currentPosition);
        /** Calculates one single Runge-Kutta step and returns true, if there is no need to continue. */
        bool calculateRungeKuttaStep(const VolumeRAM_3xFloat* volume, Streamline& pathline);
        /** Returns the velocity in the selected filter mode. */
        tgt::vec3 getVelocityAt(const VolumeRAM_3xFloat* volume, const tgt::vec3& pos) const;

    private:
        //-----------
        //  Members
        //-----------
        std::function<float()> rnd;
        //derived from properties
        const VolumeRAM_3xFloat* firstTimeStepFlow_;
        const VolumeList* flow_;               ///< input flow
        StreamlineList* output_;               ///< output, which will be used by the processor
        size_t maxNumPathlines_;                ///< maximal number of streamlines
        tgt::svec2 streamlineLengthThreshold_; ///< streamline length must be in this interval
        tgt::vec2 absoluteMagnitudeThreshold_; ///< only magnitudes in this intervall are used
        PathlineCreator::FilterMode filterMode_; ///< filtering inside the volume
    };

}   // namespace

#endif  // VRN_PATHLINECREATORBACKGROUNDTHREAD_H
