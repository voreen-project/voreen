#include "seeds.h"

#include "voreen/core/datastructures/geometry/pointsegmentlistgeometry.h"
#include "voreen/core/datastructures/geometry/geometrysequence.h"

namespace voreen {

static void getSeedListsFromGeom(const Geometry* geom, const tgt::mat4& transform, PointSegmentListGeometry<tgt::vec3>& seeds);

static void getSeedListFromPointSegmentList(const PointSegmentListGeometry<tgt::vec3>& seedList, const tgt::mat4& transform, PointSegmentListGeometry<tgt::vec3>& seeds) {
    auto transformMat = transform * seedList.getTransformationMatrix();
    for (int j=0; j<seedList.getNumSegments(); j++) {
        std::vector<tgt::vec3> points;
        for(auto& vox : seedList.getSegment(j)) {
            points.push_back(transformMat.transform(vox));
        }
        seeds.addSegment(points);
    }
}
static void getSeedListFromSequence(const GeometrySequence& geoms, const tgt::mat4& transform, PointSegmentListGeometry<tgt::vec3>& seeds) {
    for (int j=0; j<geoms.getNumGeometries(); j++) {
        getSeedListsFromGeom(geoms.getGeometry(j), geoms.getTransformationMatrix(), seeds);
    }
}
static void getSeedListsFromGeom(const Geometry* geom, const tgt::mat4& transform, PointSegmentListGeometry<tgt::vec3>& seeds) {
    if(auto* seedList = dynamic_cast<const PointSegmentListGeometry<tgt::vec3>* >(geom)) {
        getSeedListFromPointSegmentList(*seedList, transform, seeds);
    } else if (auto* geoms = dynamic_cast<const GeometrySequence* >(geom)) {
        getSeedListFromSequence(*geoms, transform, seeds);
    } else {
        LWARNINGC("voreen.RandomWalker.util", "Invalid geometry. PointSegmentListGeometry<vec3> or GeometrySequence expected.");
    }
}
void getSeedListsFromPorts(std::vector<PortDataPointer<Geometry>>& geom, PointSegmentListGeometry<tgt::vec3>& seeds) {
    auto mat = tgt::mat4::createIdentity();
    for (size_t i=0; i<geom.size(); i++) {
        getSeedListsFromGeom(geom.at(i).get(), mat, seeds);
    }
}

}
