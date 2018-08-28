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

#ifndef VRN_STREAMLINEBUNDLEDETECTORBACKGROUNDTHREAD_H
#define VRN_STREAMLINEBUNDLEDETECTORBACKGROUNDTHREAD_H

#include "voreen/core/utils/backgroundthread.h"

#include "modules/flowreen/processors/streamline/streamlinecreator.h"
#include "modules/flowreen/datastructures/streamlinelist.h"

namespace voreen {

    /**
    * Background thread used to detect streamline bundles.
    */
    class VRN_CORE_API StreamlineBundleDetectorBackgroundThread : public ProcessorBackgroundThread<StreamlineCreator> {
    public:
        /** Constructor */
        StreamlineBundleDetectorBackgroundThread(StreamlineCreator* processor,
                                                 StreamlineList* output,
                                                 float positionDeviation,
                                                 size_t noiseThreshold,
                                                 size_t resampleSize);
        /** Destructor */
        virtual ~StreamlineBundleDetectorBackgroundThread();

    protected:
        /** Main-Function used to detect the bundles. */
        virtual void threadMain();

        //-----------------
        //  Helpers
        //-----------------
        float calculateMDF(const Streamline& s, const Streamline& t) const;
        //float calculateMAM(const Streamline& s, const Streamline& t) const;

    private:
        //-----------
        //  Members
        //-----------
        StreamlineList* output_;    ///< output, which will be used by the processor

        // settings
        const float distanceThreshold_;
        const size_t noiseThreshold_;
        const size_t resampleSize_;
    };

}   // namespace

#endif  // VRN_STREAMLINEBUNDLEDETECTORBACKGROUNDTHREAD_H
