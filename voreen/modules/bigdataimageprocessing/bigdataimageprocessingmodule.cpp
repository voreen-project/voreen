/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "bigdataimageprocessingmodule.h"

#include "processors/connectedcomponentanalysis.h"
#include "processors/largevolumedistancetransform.h"
#include "processors/largevolumeformatconversion.h"
#include "processors/volumeargmax.h"
#include "processors/volumefilterlist.h"
#include "processors/volumeresampletransformation.h"
#include "io/lz4slicevolumefilereader.h"

// nuclei cluster splitting
#include "operators/volumeoperatordistancetransform.h"
#include "operators/volumeoperatorwatershed.h"
#include "operators/volumeoperatorgradientdescent.h"
#include "operators/volumeoperatorfastvolumecombine.h"

#include "processors/nucleiclustersplitting.h"

namespace voreen {

BigDataImageProcessingModule::BigDataImageProcessingModule(const std::string& modulePath)
    : VoreenModule(modulePath)
{
    setID("Big Data Image Processing");
    setGuiName("Big Data Image Processing");

    registerProcessor(new ConnectedComponentAnalysis());
    registerProcessor(new LargeVolumeDistanceTransform());
    registerProcessor(new LargeVolumeFormatConversion());
    registerProcessor(new NucleiClusterSplitting());
    registerProcessor(new VolumeArgMax());
    registerProcessor(new VolumeFilterList());
    registerProcessor(new VolumeResampleTransformation());

    registerVolumeReader(new LZ4SliceVolumeFileReader());

    // instantiate volume operators (nuclei cluster splitting)
    INST_SCALAR_TYPES(VolumeOperatorSquaredEuclideanDistanceTransform, VolumeOperatorSquaredEuclideanDistanceTransformGeneric)
    INST_SCALAR_TYPES(VolumeOperatorEuclideanDistanceTransform, VolumeOperatorEuclideanDistanceTransformGeneric)
    //INST_SCALAR_TYPES(VolumeOperatorManhattanDistanceTransform, VolumeOperatorManhattanDistanceTransformGeneric)
    //INST_SCALAR_TYPES(VolumeOperatorChebychevDistanceTransform, VolumeOperatorChebychevDistanceTransformGeneric)
    INST_SCALAR_TYPES(VolumeOperatorWatershedTransform, VolumeOperatorWatershedTransformGeneric);
    INST_SCALAR_TYPES(VolumeOperatorGradientDescent, VolumeOperatorGradientDescentGeneric);
    INST_SCALAR_TYPES(VolumeOperatorFastVolumeCombine, VolumeOperatorFastVolumeCombineGeneric);
}

} // namespace
