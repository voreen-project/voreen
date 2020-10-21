#include "seeds.h"

#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"
#include "voreen/core/datastructures/geometry/geometrysequence.h"

namespace voreen {

void getSeedListsFromPorts(std::vector<PortDataPointer<Geometry>>& geom, PointSegmentListGeometryVec3& seeds) {
    auto mat = tgt::mat4::createIdentity();
    for (size_t i=0; i<geom.size(); i++) {
        if(!seeds.collectSegmentsFromGeometry(*geom.at(i).get())) {
            LWARNINGC("voreen.RandomWalker.util", "Invalid geometry. PointSegmentListGeometry<vec3> or GeometrySequence expected.");
        }
    }
}

}
