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

#include "poipointsegmentgeometryimporter.h" 
#include <map>
#include <sstream>
#include <voreen/core/datastructures/geometry/pointsegmentlistgeometry.h>

namespace voreen {

POIPointSegmentImporter::POIPointSegmentImporter()
    : Processor()
    , outport_(Port::OUTPORT, "outport", "POI Outport", true)
    , geometryport_(Port::INPORT, "geometryport", "Geometryport")
{
    addPort(outport_);
    addPort(geometryport_);
}

Processor* POIPointSegmentImporter::create() const {
    return new POIPointSegmentImporter();
}

std::string POIPointSegmentImporter::getClassName() const {
    return "POIPointSegmentImporter";
}

std::string POIPointSegmentImporter::getCategory() const {
    return "Points of Interest";
}

void POIPointSegmentImporter::setDescriptions() {
    setDescription("The POIPointListGeometryConverter converts the point to a PointSegmentGeometry, "
        "to allow their use with the usual voreen geometry pipeline. For each group, a segment "
        "with all its point is generated.");
}

void POIPointSegmentImporter::process() {
    const Geometry* geom = geometryport_.getData();
    const PointSegmentListGeometryVec3* psl = dynamic_cast<const PointSegmentListGeometryVec3*>(geom);
    if (!psl){
        LERROR("Geometry needs to be a PointSegmentList");
        return;
    }

    POIList *pois = new POIList;
    int segmentCount = psl->getNumSegments();
    for(int i = 0; i != segmentCount; i++){
        std::stringstream ss;
        ss << i;
        POIGroupID gid = pois->addGroup(ss.str());
        std::vector<tgt::vec3> segment = psl->getSegment(i);
        for(tgt::vec3 v : segment){
            pois->addPoint(v, gid);
        }
    }

    outport_.setData(pois, true);
}

} // namespace
