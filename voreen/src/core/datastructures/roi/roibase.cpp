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

#include "voreen/core/datastructures/roi/roibase.h"
#include "voreen/core/datastructures/geometry/trianglemeshgeometry.h"

#include "voreen/core/io/serialization/serialization.h"

namespace voreen {

using tgt::ivec3;
using tgt::vec3;
using tgt::vec4;
using tgt::mat4;
using tgt::plane;

//-----------------------------------------------------------------------------

Grid::Grid() : origin_(0.0f), xVec_(1.0f, 0.0f, 0.0f), yVec_(0.0f, 1.0f, 0.0f), zVec_(0.0f, 0.0f, 1.0f), spacing_(1.0f) {
}

Grid::Grid(tgt::vec3 origin, tgt::vec3 spacing, tgt::vec3 xVec, tgt::vec3 yVec, tgt::vec3 zVec) : origin_(origin), xVec_(xVec), yVec_(yVec), zVec_(zVec), spacing_(spacing) {
}

tgt::mat4 Grid::getPhysicalToWorldMatrix() const {
    mat4 m(xVec_.x, yVec_.x, zVec_.x, origin_.x,
           xVec_.y, yVec_.y, zVec_.y, origin_.y,
           xVec_.z, yVec_.z, zVec_.z, origin_.z,
           0.0f,    0.0f,    0.0f,   1.0f);
    return m;
}

tgt::mat4 Grid::getWorldToPhysicalMatrix() const {
    mat4 m = getPhysicalToWorldMatrix();
    mat4 inv;
    m.invert(inv);
    return inv;
}

void Grid::serialize(Serializer& s) const {
    s.serialize("origin", origin_);
    s.serialize("xVec", xVec_);
    s.serialize("yVec", yVec_);
    s.serialize("zVec", zVec_);
    s.serialize("spacing", spacing_);
}

void Grid::deserialize(Deserializer& s) {
    s.deserialize("origin", origin_);
    s.deserialize("xVec", xVec_);
    s.deserialize("yVec", yVec_);
    s.deserialize("zVec", zVec_);
    s.deserialize("spacing", spacing_);
}

//-----------------------------------------------------------------------------

Grid gridFromVolume(const VolumeBase* vol) {
    mat4 m = vol->getPhysicalToWorldMatrix();
    vec3 os = vol->getOffset();

    vec3 o = m * os;
    vec3 x = m * (os + vec3(1.0f, 0.0f, 0.0f));
    vec3 y = m * (os + vec3(0.0f, 1.0f, 0.0f));
    vec3 z = m * (os + vec3(0.0f, 0.0f, 1.0f));

    return Grid(o, vol->getSpacing(), x - o, y - o, z - o);
}

//-----------------------------------------------------------------------------

ROIMeshCache::ROIMeshCache() : mesh3D_(0), rasterMesh3D_(0) {
}

ROIMeshCache::~ROIMeshCache() {
    clear();
}

Geometry* ROIMeshCache::getMesh() const {
   return mesh3D_;
}

void ROIMeshCache::setMesh(Geometry* m) {
    delete mesh3D_;
    mesh3D_ = m;
}

Geometry* ROIMeshCache::getRasterMesh() const {
   return rasterMesh3D_;
}

void ROIMeshCache::setRasterMesh(Geometry* m) {
    delete rasterMesh3D_;
    rasterMesh3D_ = m;
}

Geometry* ROIMeshCache::getMesh(tgt::plane pl) const {
    for(std::list<CacheEntry>::iterator it=planeMeshes_.begin(); it != planeMeshes_.end(); it++) {
        if((*it).plane_.toVec4() == pl.toVec4()) {
            // move to front:
            CacheEntry ce = *it;
            planeMeshes_.erase(it);
            planeMeshes_.push_front(ce);
            return ce.mesh_;
        }
    }
    return 0;
}

void ROIMeshCache::setMesh(Geometry* m, tgt::plane pl) {
    if(getMesh(pl)) {
        //TODO: duplicate
        return;
    }
    else {
        planeMeshes_.push_front(CacheEntry(pl, m));

        // Limit size:
        while(planeMeshes_.size() > 10) {
            delete planeMeshes_.back().mesh_;
            planeMeshes_.pop_back();
        }
    }
}

Geometry* ROIMeshCache::getRasterMesh(tgt::plane pl) const {
    for(std::list<CacheEntry>::iterator it=planeRasterMeshes_.begin(); it != planeRasterMeshes_.end(); it++) {
        if((*it).plane_.toVec4() == pl.toVec4()) {
            // move to front:
            CacheEntry ce = *it;
            planeRasterMeshes_.erase(it);
            planeRasterMeshes_.push_front(ce);
            return ce.mesh_;
        }
    }
    return 0;
}

void ROIMeshCache::setRasterMesh(Geometry* m, tgt::plane pl) {
    if(getRasterMesh(pl)) {
        //TODO: duplicate
        return;
    }
    else {
        planeRasterMeshes_.push_front(CacheEntry(pl, m));

        // Limit size:
        while(planeRasterMeshes_.size() > 10) {
            delete planeRasterMeshes_.back().mesh_;
            planeRasterMeshes_.pop_back();
        }
    }
}

void ROIMeshCache::clear() {
    delete mesh3D_;
    mesh3D_ = 0;

    delete rasterMesh3D_;
    rasterMesh3D_ = 0;

    while(!planeMeshes_.empty()) {
        delete planeMeshes_.back().mesh_;
        planeMeshes_.pop_back();
    }

    while(!planeRasterMeshes_.empty()) {
        delete planeRasterMeshes_.back().mesh_;
        planeRasterMeshes_.pop_back();
    }
}

//-----------------------------------------------------------------------------

const std::string ROIBase::loggerCat_ = "voreen.ROIBase";

ROIBase::ROIBase() : PropertyOwner(),
    name_("name", "Name", "", METADATA_CHANGE),
    comment_("comment", "Comment", "", METADATA_CHANGE),
    isVisible_("isVisible", "Visible", true, METADATA_CHANGE),
    isSelected_("isSelected", "Selected", false, METADATA_CHANGE),
    snapToGrid_("snapToGrid", "Snap to Grid", true, METADATA_CHANGE),
    childrenVisible_("childrenVisible", "Children Visible", true, METADATA_CHANGE),
    //color_("color", "Color", vec4(1.0f), vec4(0.0f), vec4(1.0f), METADATA_CHANGE),
    color_("color", "Color", vec4(1.0f), vec4(0.0f), vec4(1.0f), ROI_CHANGE),
    transformation_("transformation", "Transformation", tgt::mat4::identity, tgt::mat4(-1000000.0f), tgt::mat4(100000.0f), ROI_CHANGE)
{
    addProperty(name_);
    addProperty(comment_);
    addProperty(isVisible_);
    addProperty(isSelected_);
    addProperty(snapToGrid_);
    addProperty(childrenVisible_);
    addProperty(color_);
    addProperty(transformation_);
}

ROIBase::ROIBase(Grid grid) : PropertyOwner(),
    name_("name", "Name", "", METADATA_CHANGE),
    comment_("comment", "Comment", "", METADATA_CHANGE),
    isVisible_("isVisible", "Visible", true, METADATA_CHANGE),
    isSelected_("isSelected", "Selected", false, METADATA_CHANGE),
    snapToGrid_("snapToGrid", "Snap to Grid", true, METADATA_CHANGE),
    childrenVisible_("childrenVisible", "Children Visible", true, METADATA_CHANGE),
    //color_("color", "Color", vec4(1.0f), vec4(0.0f), vec4(1.0f), METADATA_CHANGE),
    color_("color", "Color", vec4(1.0f), vec4(0.0f), vec4(1.0f), ROI_CHANGE),
    transformation_("transformation", "Transformation", tgt::mat4::identity, tgt::mat4(-1000000.0f), tgt::mat4(100000.0f), ROI_CHANGE),
    grid_(grid)
{
    addProperty(name_);
    addProperty(comment_);
    addProperty(isVisible_);
    addProperty(isSelected_);
    addProperty(snapToGrid_);
    addProperty(childrenVisible_);
    addProperty(color_);
    addProperty(transformation_);
}

ROIBase* ROIBase::clone() const {
    std::stringstream stream;

    // first serialize
    XmlSerializer s;
    s.serialize("this", this);
    s.write(stream);

    // then deserialize again
    XmlDeserializer d;
    d.read(stream);
    ROIBase* cloned = 0;
    d.deserialize("this", cloned);

    return cloned;
}

std::string ROIBase::getID() const {
    return name_.get();
}

std::string ROIBase::getGuiName() const {
    return name_.get();
}
void ROIBase::setGuiName(const std::string& name) {
    name_.set(name);
    id_ = name;
    guiName_ = name;
}

void ROIBase::setID(const std::string& id){
    name_.set(id);
    id_ = id;
    guiName_ = id;
}

VolumeRAM* ROIBase::rasterize(Grid /*g*/) const {
    return 0;
}

VolumeRAM* ROIBase::rasterize() const {
    return rasterize(getGrid());
}

Geometry* ROIBase::generateRasterMesh() const {
    return generateRasterMesh(getGrid());
}

Geometry* ROIBase::generateRasterMesh(const tgt::plane& pl) const
{
    return generateRasterMesh(pl, getGrid());
}

Geometry* ROIBase::generateRasterMesh(Grid /*g*/) const {
    //TODO
    return new TriangleMeshGeometrySimple();
}

Geometry* ROIBase::getMesh() const {
    Geometry* m = cache_.getMesh();
    if(!m) {
        m = generateMesh();
        cache_.setMesh(m);
    }
    return m;
}

Geometry* ROIBase::getRasterMesh() const {
    Geometry* m = cache_.getRasterMesh();
    if(!m) {
        m = generateRasterMesh();
        if(m)
            cache_.setMesh(m);
    }
    return m;
}

Geometry* ROIBase::getMesh(tgt::plane pl) const {
    Geometry* m = cache_.getMesh(pl);
    if(!m) {
        m = generateMesh(pl);
        if(m)
            cache_.setMesh(m, pl);
    }
    return m;
}

Geometry* ROIBase::getRasterMesh(tgt::plane pl) const {
    Geometry* m = cache_.getRasterMesh(pl);
    if(!m) {
        m = generateRasterMesh(pl);
        if(m)
            cache_.setRasterMesh(m, pl);
    }
    return m;
}

Geometry* ROIBase::generateRasterMesh(const tgt::plane& pl, Grid g) const {
    // Transform plane to grid coordinates:
    mat4 m = g.getWorldToPhysicalMatrix();
    const plane physicalPlane = pl.transform(m);

    tgt::Bounds bb = getBoundingBox(); // FIXME: BB is in own grid

    // Get BB in voxel coordinates:
    vec3 sp = g.getSpacing();
    ivec3 llfV(floor(bb.getLLF() / sp));
    ivec3 urbV(ceil(bb.getURB() / sp));

    // Optimize for orthogonal planes:
    int optDir = -1;
    float optDist = 0.0f;
    if(fabs(length(physicalPlane.n) - 1.0f) < 0.001f) {
        if(fabs(fabs(physicalPlane.n.x) - 1.0f) < 0.001f)
            optDir = 0;
        else if(fabs(fabs(physicalPlane.n.y) - 1.0f) < 0.001f)
            optDir = 1;
        else if(fabs(fabs(physicalPlane.n.z) - 1.0f) < 0.001f)
            optDir = 2;

        if(optDir >= 0) {
            optDist = physicalPlane.d * physicalPlane.n[optDir];
            llfV[optDir] = tgt::ifloor(optDist / sp[optDir]);
            urbV[optDir] = tgt::iceil(optDist / sp[optDir]);
        }
    }

    TriangleMeshGeometrySimple* planeFace = createQuad(physicalPlane);

    tgt::plane xm(-1.0f, 0.0f, 0.0f, 0.0f);
    tgt::plane xp(1.0f, 0.0f, 0.0f, 0.0f);
    tgt::plane ym(0.0f, -1.0f, 0.0f, 0.0f);
    tgt::plane yp(0.0f, 1.0f, 0.0f, 0.0f);
    tgt::plane zm(0.0f, 0.0f, -1.0f, 0.0f);
    tgt::plane zp(0.0f, 0.0f, 1.0f, 0.0f);

    TriangleMeshGeometrySimple* mg = new TriangleMeshGeometrySimple();
    int numFaces = 0;
    float d = max(sp);
    ivec3 i;
    for(i.x=llfV.x; i.x < urbV.x; i.x++) {
        xm.d = -(i.x * sp.x);
        xp.d = (i.x+1) * sp.x;
        for(i.y=llfV.y; i.y < urbV.y; i.y++) {
            ym.d = -(i.y * sp.y);
            yp.d = (i.y+1) * sp.y;
            for(i.z=llfV.z; i.z < urbV.z; i.z++) {
                zm.d = -(i.z * sp.z);
                zp.d = (i.z+1) * sp.z;

                // Calc bounds of this voxel:
                vec3 llf = vec3(i.x * sp.x, i.y * sp.y, i.z * sp.z);
                vec3 urb = llf + sp;

                vec3 c = (llf + urb) / 2.0f;

                if(fabs(dot(c, physicalPlane.n) - physicalPlane.d) > (2*d))
                    continue;

                if(!inROI(c))
                    continue;

                if(optDir == -1) {
                    TriangleMeshGeometrySimple* fg = createQuad(physicalPlane);
                    fg->clip(xm);
                    fg->clip(xp);
                    fg->clip(ym);
                    fg->clip(yp);
                    fg->clip(zm);
                    fg->clip(zp);

                    if(!fg->isEmpty()) {
                        mg->addMesh(fg);
                        numFaces++;
                    }
                }
                else if(optDir == 0) {
                    mg->addQuad(
                    VertexBase(vec3(optDist, llf.y, llf.z)),
                    VertexBase(vec3(optDist, llf.y, urb.z)),
                    VertexBase(vec3(optDist, urb.y, urb.z)),
                    VertexBase(vec3(optDist, urb.y, llf.z)));

                    numFaces++;
                }
                else if(optDir == 1) {
                    mg->addQuad(
                    VertexBase(vec3(llf.x, optDist, llf.z)),
                    VertexBase(vec3(llf.x, optDist, urb.z)),
                    VertexBase(vec3(urb.x, optDist, urb.z)),
                    VertexBase(vec3(urb.x, optDist, llf.z)));

                    numFaces++;
                }
                else if(optDir == 2) {
                    mg->addQuad(
                    VertexBase(vec3(llf.x, llf.y, optDist)),
                    VertexBase(vec3(llf.x, urb.y, optDist)),
                    VertexBase(vec3(urb.x, urb.y, optDist)),
                    VertexBase(vec3(urb.x, llf.y, optDist)));

                    numFaces++;
                }
            }
        }
    }

    mg->transform(g.getPhysicalToWorldMatrix());
    return mg;
}

Grid ROIBase::getGrid() const {
    return grid_;
}

void ROIBase::serialize(Serializer& s) const {
    PropertyOwner::serialize(s);
    s.serialize("Grid", grid_);
}

void ROIBase::deserialize(Deserializer& s) {
    PropertyOwner::deserialize(s);
    s.deserialize("Grid", grid_);
}

std::vector<std::string> ROIBase::getAvailableStatistics() const {
    std::vector<std::string> stats;
    for (std::map<std::string, MetaDataBase*>::const_iterator it = statistics_.begin(); it != statistics_.end(); ++it) {
        stats.push_back(it->first);
    }
    return stats;
}

const MetaDataBase* ROIBase::getStatistic(const std::string& id) const {
    if (statistics_.find(id) != statistics_.end())
        return statistics_.find(id)->second;
    else
        return 0;
}

void ROIBase::addStatistic(const std::string& id, MetaDataBase* stat) const {
    removeStatistic(id);
    statistics_[id] = stat;
}

void ROIBase::removeStatistic(const std::string& id) const {
    if (getStatistic(id) == 0)
        return;

    delete statistics_[id];
    statistics_.erase(statistics_.find(id));
}

void ROIBase::clearStatistics() const {
    for (std::map<std::string, MetaDataBase*>::iterator it = statistics_.begin(); it != statistics_.end(); ++it)
        delete it->second;

    statistics_.clear();
}

void ROIBase::invalidate(int inv) {
    if(inv >= ROI_CHANGE) {
        clearStatistics();
        cache_.clear();
    }

    PropertyOwner::invalidate(inv);
}

//-----------------------------------------------------------------------------

TriangleMeshGeometrySimple* createQuad(const tgt::plane& pl) {
    vec3 temp(1.0f, 1.0f, 0.0f);
    if(dot(temp, pl.n) > 0.9f)
        temp = vec3(0.0f, 1.0f, 1.0f);

    vec3 t1 = normalize(cross(pl.n, temp));
    vec3 t2 = normalize(cross(pl.n, t1));
    t1 *= 500.0f; //TODO: automatically use appropriate size
    t2 *= 500.0f;
    vec3 base = pl.d * pl.n;

    TriangleMeshGeometrySimple* face = new TriangleMeshGeometrySimple();
    VertexBase v1(base + t1 + t2);
    VertexBase v2(base + t1 - t2);
    VertexBase v3(base - t1 - t2);
    VertexBase v4(base - t1 + t2);
    face->addQuad(v1, v2, v3, v4);
    return face;
}

std::pair<tgt::ivec3, tgt::ivec3> getVoxelBoundingBox(const VolumeBase* vh, const ROIBase* roi) {
    tgt::mat4 m = tgt::mat4::createTranslation(-vh->getOffset()) * vh->getWorldToPhysicalMatrix() * roi->getGrid().getPhysicalToWorldMatrix();
    tgt::Bounds b = roi->getBoundingBox().transform(m);
    vec3 sp = vh->getSpacing();

    tgt::ivec3 llf = tgt::floor(b.getLLF() / sp);
    llf = tgt::max(llf, ivec3(0));
    tgt::ivec3 urb = tgt::ceil(b.getURB() / sp);
    urb = tgt::min(urb, ivec3(vh->getDimensions()));

    return std::pair<ivec3, ivec3>(llf, urb);
}

} // namespace
