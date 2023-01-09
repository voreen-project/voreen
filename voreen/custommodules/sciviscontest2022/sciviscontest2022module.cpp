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
#include "sciviscontest2022module.h"

#include "processors/binary_geometry/binarygeometrysave.h"
#include "processors/binary_geometry/binarygeometrysource.h"
#include "processors/binary_geometry/binarygeometrysequencesource.h"
#include "processors/geometry/geometryselector.h"
#include "processors/geometry/geometrysequencesource.h"
#include "processors/geometry/curvecreator.h"
#include "processors/vorticity/vorticityfieldcreator.h"
#include "processors/vorticity/pyrogenicvorticitymapper.h"
#include "processors/vorticity/divergencefieldcreator.h"
#include "processors/vorticity/materialderivative.h"
#include "processors/geometry/boundingboxsampler.h"
#include "processors/surface/streamsurfacecreator.h"
#include "processors/surface/pathsurfacecreator.h"
#include "processors/surface/streamsurfaceviewer.h"
#include "processors/corelinelengthvolumecreator.h"
#include "processors/similarityvolumecreator.h"
#include "processors/surface/pathsurfacerenderer.h"
#include "processors/volumeinterpolation.h"
#include "processors/seededstreamlinecreator.h"
#include "processors/seededpathlinecreator.h"
#include "processors/volumekernel.h"

#include "properties/link/linkevaluatorvectorelement.h"

namespace voreen {

SciVisContest2022Module::SciVisContest2022Module(const std::string& modulePath)
    : VoreenModule(modulePath)
{
    setID("SciVisContest2022");
    setGuiName("SciVis Contest 2022");
    
    // each module processor needs to be registered
	registerProcessor(new BinaryGeometrySave());
	registerProcessor(new BinaryGeometrySource());
	registerProcessor(new BinaryGeometrySequenceSource());
    registerProcessor(new GeometrySelector());
	registerProcessor(new GeometrySequenceSource());
	registerProcessor(new StreamSurfaceCreator());
    registerProcessor(new PathSurfaceCreator());
    registerProcessor(new PathSurfaceRenderer());
	registerProcessor(new StreamSurfaceViewer());
	registerProcessor(new CorelineLengthVolumeCreator());
	registerProcessor(new SimilarityVolumeCreator());
	registerProcessor(new VolumeInterpolation());
	registerProcessor(new BoundingBoxSampler());
    registerProcessor(new CurveCreator());
    registerProcessor(new SeededStreamlineCreator());
    registerProcessor(new SeededPathlineCreator());
    registerProcessor(new VorticityFieldCreator());
    registerProcessor(new PyrogenicVorticityMapper());
    registerProcessor(new VolumeKernel());
    registerProcessor(new DivergenceFieldCreator());
    registerProcessor(new MaterialDerivative());
    
    registerSerializableType(new LinkEvaluatorIntVectorElementFront());
    registerSerializableType(new LinkEvaluatorIntVectorElementBack());
}

}