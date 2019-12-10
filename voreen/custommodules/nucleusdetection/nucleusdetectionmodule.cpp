/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "nucleusdetectionmodule.h"


// patch export, k-means import, caffe training data
#include "processors/patchextractor.h"
#include "processors/patchlistreader.h"

#include "processors/tileextractor.h"


#ifdef VRN_NUCLEUSDETECTION_CBLAS_FOUND
#include "processors/patchtrainingdataextractor.h"
#endif

// caffe classification
#ifdef VRN_NUCLEUSDETECTION_CAFFE_FOUND
#include "processors/patchcaffeclassifier.h"
#include "processors/tilecaffeclassifier.h"
#endif

#include "processors/roidetector.h"

namespace voreen {

NucleusDetectionModule::NucleusDetectionModule(const std::string& modulePath)
    : VoreenModule(modulePath)
{
    setID("nucleusdetection");
    setGuiName("Cell Nucleus Detection");

    registerProcessor(new PatchExtractor());
    registerProcessor(new PatchListReader());

    registerProcessor(new TileExtractor());

#ifdef VRN_NUCLEUSDETECTION_CBLAS_FOUND
    registerProcessor(new PatchTrainingDataExtractor());
#endif

#ifdef VRN_NUCLEUSDETECTION_CAFFE_FOUND
    registerProcessor(new PatchCaffeClassifier());
    registerProcessor(new TileCaffeClassifier());
#endif

    registerProcessor(new RoiDetector());
}

} // namespace
