/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
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

#include "roiskeletonize.h"
#include "../pfskelmodule.h"
#include "../utils/pfskelwrapper.h"

#include "voreen/core/datastructures/volume/volumeatomic.h"
#include "voreen/core/datastructures/roi/roiraster.h"
#include "voreen/core/voreenapplication.h"

#include "modules/core/io/datvolumewriter.h"
#include "modules/core/io/datvolumereader.h"

#include "modules/roi/processors/roistorage.h"

namespace voreen {

using tgt::mat4;
using tgt::vec3;
using tgt::svec3;
using tgt::ivec3;

ROISkeletonize::ROISkeletonize()
    : ROISegmenter()
    , fieldSt_("fieldSt", "Field Strength", 4, 4, 10)
    , highDiv_("highDiv", "High Divergence %", 33, 0, 100)
    , branchTh_("branchTh", "Short Branch Threshold", 8, 2, 100)
{
    addProperty(fieldSt_);
    addProperty(highDiv_);
    addProperty(branchTh_);
}

bool ROISkeletonize::isCompatible(std::set<ROIBase*> rois) {
    return ( (rois.size() == 1) && dynamic_cast<ROIRaster*>(*rois.begin()) );
}

void ROISkeletonize::process() {
    if(!isActive())
        return;
}

void ROISkeletonize::segment() {
    ROIBase* roi = *rois_.begin();
    ROIRaster* rr = dynamic_cast<ROIRaster*>(roi);

    Graph g = PFSkelWrapper::calculateSkeleton(rr, fieldSt_.get(), highDiv_.get(), branchTh_.get());
    ROIGraph* rg = new ROIGraph(rr->getGrid(), g);
    rg->setGuiName("PFSkel(" + roi->getGuiName() + ", " + itos(fieldSt_.get()) + ", " + itos(highDiv_.get()) + ")");
    rg->setComment("Skeleton of " + roi->getGuiName() + " FieldSt: " + itos(fieldSt_.get()) + ", HighDiv:" + itos(highDiv_.get()) + ", BranchTh: " + itos(branchTh_.get()));
    rg->setColor(ROIStorage::getNextColor());
    storage_->addROI(rg);
}

Geometry* ROISkeletonize::getRasterMesh(tgt::plane /*pl*/) const {
    return 0;
}

void ROISkeletonize::activate() {
    //ROIBase* roi = *rois_.begin();
    ROISegmenter::activate();
}

void ROISkeletonize::deactivate() {
    ROISegmenter::deactivate();
}

} // namespace voreen
