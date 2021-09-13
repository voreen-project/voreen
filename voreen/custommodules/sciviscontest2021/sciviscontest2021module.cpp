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

//incude header file
#include "sciviscontest2021module.h"

// include classes to be registered
#include "processors/componentidselector.h"
#include "processors/connectedcomponenttracker.h"
#include "processors/diskseedpointcreator.h"
#include "processors/distancetooriginplot.h"
#include "processors/featureextractor.h"
#include "processors/sphericalcuttingplaneselector.h"
#include "processors/sphericalfakevolume.h"
#include "processors/sphericalgeometrytransformation.h"
#include "processors/sphericalstreamlinetransformation.h"
#include "processors/sphericalraycaster.h"
#include "processors/sphericalvolumelistproxy.h"
#include "processors/sphericalvolumeproxy.h"
#include "processors/starcoordinates.h"
#include "processors/timeserieslistcreator.h"
#include "processors/timeserieslistsave.h"
#include "processors/timeserieslistsource.h"
#include "processors/timeseriesplot.h"
#include "processors/volumeoverlaptracker.h"

//use voreen namespace
namespace voreen {

SciVisContest2021Module::SciVisContest2021Module(const std::string& modulePath)
    : VoreenModule(modulePath)
{
    setID("SciVisContest2021");
    setGuiName("SciVis Contest 2021");

    // adds the module glsl dir to the shader search path (if shaders are needed in the module)
    addShaderPath(getModulePath("glsl"));

    // processors
    registerProcessor(new ComponentIdSelector());
    registerProcessor(new ConnectedComponentTracker());
    registerProcessor(new DiskSeedPointCreator());
    registerProcessor(new DistanceToOriginPlot());
    registerProcessor(new FeatureExtractor());
    registerProcessor(new SphericalCuttingPlaneSelector());
    registerProcessor(new SphericalFakeVolume());
    registerProcessor(new SphericalGeometryTransformation());
    registerProcessor(new SphericalStreamlineTransformation());
    registerProcessor(new SphericalRaycaster());
    registerProcessor(new SphericalVolumeListProxy());
    registerProcessor(new SphericalVolumeProxy());
    registerProcessor(new StarCoordinates());
    registerProcessor(new TimeSeriesListCreator());
    registerProcessor(new TimeSeriesListSave());
    registerProcessor(new TimeSeriesListSource());
    registerProcessor(new TimeseriesPlot());
    registerProcessor(new VolumeOverlapTracker());
}


} // namespace
