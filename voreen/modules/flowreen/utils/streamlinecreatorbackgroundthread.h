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
    * Background thread used to calculate the stream lines.
    */
    class VRN_CORE_API StreamlineCreatorBackgroundThread : public ProcessorBackgroundThread<StreamlineCreator> {
    public:
        /** Constructor */
        StreamlineCreatorBackgroundThread(StreamlineCreator* processor, int seedTime, const VolumeBase* flow, StreamlineList* output,
            int maxNumStreamlines, tgt::ivec2 streamlineLengthThreshold, tgt::vec2 absoluteMagnitudeThreshold,
            StreamlineCreator::FilterMode filterMode);
        /** Destructor */
        virtual ~StreamlineCreatorBackgroundThread();
    protected:
        /** Main-Function used to calculate the streamlines */
        virtual void threadMain();
        /** Used to clean up */
        virtual void handleInterruption();
        //-----------------
        //  Helpers
        //-----------------
        /** Used to redefine a seed point. */
        void reseedPosition(const size_t currentPosition);
        /** Calcualtes a streamline */
        Streamline computeStreamlineRungeKutta(const tgt::vec3& start);
        /** Returns the vlocity in the selected filter mode. */
        const tgt::vec3 getVelocityAt(const tgt::vec3& pos);
    private:
        //-----------
        //  Members
        //-----------
        std::function<float()> rnd;
        //derived from properties
        const VolumeBase* flow_;               ///< input flow
        VolumeRAMRepresentationLock representation_; ///< representation lock for input flow
        StreamlineList* output_;               ///< output, which will be used by the processor
        size_t maxNumStreamlines_;             ///< maximal number of streamlines
        tgt::svec2 streamlineLengthThreshold_; ///< streamline length must be in this interval
        tgt::vec2 absoluteMagnitudeThreshold_; ///< only magnitudes in this intervall are used
        StreamlineCreator::FilterMode filterMode_; ///< filtering inside the volume
        //created new
        tgt::vec3* seedingPositions_;          ///< used seeding points. calculated on the fly
    };

}   // namespace

#endif  // VRN_STREAMLINECREATORBACKGROUNDTHREAD_H
