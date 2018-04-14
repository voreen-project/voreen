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

#include "poipointsegmentgeometryexporter.h" 
#include <map>
#include <sstream>
#include <voreen/core/datastructures/geometry/pointsegmentlistgeometry.h>

namespace voreen {

POIPointSegmentExporter::POIPointSegmentExporter()
    : Processor()
    , cpPort_(Port::INPORT, "cpPort", "Coprocessors", false)
    , geometryport_(Port::OUTPORT, "geometryport", "Geometryport")
{
    addPort(cpPort_);
    addPort(geometryport_);
}

Processor* POIPointSegmentExporter::create() const {
    return new POIPointSegmentExporter();
}

std::string POIPointSegmentExporter::getClassName() const {
    return "POIPointSegmentExporter";
}

std::string POIPointSegmentExporter::getCategory() const {
    return "Points of Interest";
}

void POIPointSegmentExporter::setDescriptions() {
    setDescription("The POIPointListGeometryConverter converts the point to a PointSegmentGeometry, "
        "to allow their use with the usual voreen geometry pipeline. For each group, a segment "
        "with all its point is generated.");
}

void POIPointSegmentExporter::process() {
    POIStorage* storage = cpPort_.getConnectedProcessor();
    std::map<POIGroupID, std::vector<tgt::vec3> > groups;
    for(auto g: storage->getGroups()){
        POIGroupID gid = storage->getGroupID(g);
        groups[gid] = std::vector<tgt::vec3>();
    }
    for (auto p: storage->getPoints())
    {
        groups[p.group_].push_back(p.position_);
    }
    geometry_ = std::unique_ptr<voreen::PointSegmentListGeometryVec3>(new voreen::PointSegmentListGeometryVec3());
    for (auto g: groups){
        geometry_->addSegment(g.second);
    }
    geometryport_.setData(geometry_.get(), false);
}

} // namespace
