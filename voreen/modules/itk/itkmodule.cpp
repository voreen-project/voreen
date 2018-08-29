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

#include "itkmodule.h"

//#include "io/itkvolumereader.h"
#include "processors/anisotropicdiffusion.h"
#include "processors/doublethreshold.h"
#include "processors/gradientvectorflow.h"
#include "processors/volumefilter_itk.h"
#include "processors/vesselness.h"
#include "processors/valuedregionalmaximaimagefilter.h"
#include "processors/mutualinformationregistration.h"

#include "processors/fastmarchingimagefilter.h"
#include "processors/curveslevelsetimagefilterworkflow.h"
#include "processors/geodesicactivecontourlevelsetimagefilterworkflow.h"
#include "processors/geodesicactivecontourshapepriorlevelsetimagefilterworkflow.h"
#include "processors/narrowbandcurveslevelsetimagefilterworkflow.h"
#include "processors/narrowbandthresholdsegmentationlevelsetimagefilterworkflow.h"
#include "processors/shapedetectionlevelsetimagefilterworkflow.h"
#include "processors/watershedimagefilter.h"
#include "processors/thresholdlabelerimagefilter.h"
#include "processors/labelmapoverlayimagefilter.h"
#include "processors/labelmapcontouroverlayimagefilter.h"
#include "processors/bayesianclassifierimagefilter.h"

namespace voreen {

ITKModule::ITKModule(const std::string& moduleName)
    : VoreenModule(moduleName)
{
    setID("ITK");
    setGuiName("ITK");

    registerProcessor(new AnisotropicDiffusion());
    registerProcessor(new DoubleThreshold());
    registerProcessor(new GradientVectorFlow());
    registerProcessor(new VolumeFilterITK());
    registerProcessor(new Vesselness());
    registerProcessor(new ValuedRegionalMaximaImageFilter());
    registerProcessor(new MutualInformationRegistration());

    registerProcessor(new BayesianClassifierImageFilterITK());
    registerProcessor(new CurvesLevelSetImageFilterWorkflowITK());
    registerProcessor(new FastMarchingImageFilterITK());
    registerProcessor(new GeodesicActiveContourLevelSetImageFilterWorkflowITK());
    registerProcessor(new GeodesicActiveContourShapePriorLevelSetImageFilterWorkflowITK());
    registerProcessor(new LabelMapContourOverlayImageFilterITK());
    registerProcessor(new ThresholdLabelerImageFilterITK());
    registerProcessor(new ShapeDetectionLevelSetImageFilterWorkflowITK());
    registerProcessor(new NarrowBandThresholdSegmentationLevelSetImageFilterWorkflowITK());
    registerProcessor(new NarrowBandCurvesLevelSetImageFilterWorkflowITK());
    registerProcessor(new LabelMapOverlayImageFilterITK());
    registerProcessor(new WatershedImageFilterITK());

    //registerVolumeReader(new ITKVolumeReader());
}

} // namespace
