/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "volumecurvature.h"
#include "voreen/core/datastructures/volume/volumeram.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/operators/volumeoperatorcurvature.h"

namespace voreen {

VolumeCurvature::VolumeCurvature()
    : VolumeProcessor()
    , curvatureType_("curvatureType", "Curvature Type")
    , inport_(Port::INPORT, "volumehandle.input")
    , outport_(Port::OUTPORT, "volumehandle.output", "volumehandle.output",false)
{
    curvatureType_.addOption("first", "First Principle");
    curvatureType_.addOption("second", "Second Principle");
    curvatureType_.addOption("mean", "Mean");
    curvatureType_.addOption("gauss", "Gaussian");
    addProperty(curvatureType_);

    addPort(inport_);
    addPort(outport_);
}

Processor* VolumeCurvature::create() const {
    return new VolumeCurvature();
}

void VolumeCurvature::process() {
    unsigned int curvatureType = -1;
    if (curvatureType_.get() == "first") curvatureType = 0;
    else if (curvatureType_.get() == "second") curvatureType = 1;
    else if (curvatureType_.get() == "mean") curvatureType = 2;
    else if (curvatureType_.get() == "gauss") curvatureType = 3;

    VolumeOperatorCurvature voOpCu;
    Volume* v = voOpCu.apply<float>(inport_.getData(), curvatureType);

    outport_.setData(v);
}

}   // namespace
