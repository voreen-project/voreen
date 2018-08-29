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

#include "voreen/core/datastructures/roi/roiraster.h"
#include "voreen/core/datastructures/geometry/trianglemeshgeometry.h"

#include "voreen/core/io/serialization/serialization.h"

//TODO
#include <Eigen/Core>
#include <Eigen/Eigen>

namespace voreen {

using tgt::vec3;
using tgt::ivec3;
using tgt::svec3;
using tgt::vec4;
using tgt::mat4;

const std::string ROIRaster::loggerCat_ = "voreen.ROIRaster";

ROIRaster::ROIRaster() : ROISingle()
                         , llf_(0)
                         , dims_((size_t)0)
{
}

ROIRaster::ROIRaster(Grid grid, tgt::ivec3 llf, tgt::svec3 dims, bool value) : ROISingle(grid)
                                                                , llf_(llf)
                                                                , dims_(dims)
{
    size_t numVoxels = hmul(dims_);
    voxels_.assign(numVoxels, value);
}

ROIRaster::ROIRaster(Grid grid, const ROIBase* roi) : ROISingle(grid) {
    // Get BB in grid coordinates:
    tgt::mat4 m = grid.getWorldToPhysicalMatrix() * roi->getGrid().getPhysicalToWorldMatrix();
    tgt::Bounds b = roi->getBoundingBox().transform(m);

    vec3 sp = getGrid().getSpacing();
    llf_ = tgt::floor(b.getLLF() / sp);
    tgt::ivec3 urb = tgt::ceil(b.getURB() / sp);
    dims_ = urb - llf_;

    size_t numVoxels = hmul(dims_);
    //LINFO("nv: " << numVoxels << dims_);
    voxels_.assign(numVoxels, false);

    svec3 i;
    for(i.x=0; i.x<dims_.x; i.x++) {
        for(i.y=0; i.y<dims_.y; i.y++) {
            for(i.z=0; i.z<dims_.z; i.z++) {
                vec3 p = (vec3(ivec3(i)+llf_) + vec3(0.5f)) * sp;
                if(roi->inROI(p))
                    setVoxel(i, true);
            }
        }
    }
    optimize();
}

ROIRaster::ROIRaster(const ROIBase* roi) : ROISingle(roi->getGrid()) {
    tgt::Bounds b = roi->getBoundingBox();

    vec3 sp = getGrid().getSpacing();
    llf_ = tgt::floor(b.getLLF() / getGrid().getSpacing());
    tgt::ivec3 urb = tgt::ceil(b.getURB() / sp);
    dims_ = urb - llf_;

    size_t numVoxels = hmul(dims_);
    //LINFO("nv: " << numVoxels << dims_);
    voxels_.assign(numVoxels, false);

    svec3 i;
    for(i.x=0; i.x<dims_.x; i.x++) {
        for(i.y=0; i.y<dims_.y; i.y++) {
            for(i.z=0; i.z<dims_.z; i.z++) {
                vec3 p = (vec3(ivec3(i)+llf_) + vec3(0.5f)) * sp;
                if(roi->inROI(p))
                    setVoxel(i, true);
            }
        }
    }
    optimize();
}

void ROIRaster::serialize(Serializer& s) const {
    s.serialize("llf", ivec3(llf_));
    s.serialize("dims", ivec3(dims_));

    std::string tempString = "";

    for(size_t i=0; i<voxels_.size(); i++) {
        if(voxels_[i] == true)
            tempString += "1";
        else
            tempString += "0";
    }
    s.serialize("Voxels", tempString);

    ROISingle::serialize(s);
}

void ROIRaster::deserialize(Deserializer& s) {
    ivec3 llf;
    s.deserialize("llf", llf);
    llf_ = llf;
    ivec3 dims;
    s.deserialize("dims", dims);
    dims_ = dims;

    size_t numVoxels = hmul(dims_);
    voxels_.assign(numVoxels, false);

    std::string tempString = "";
    s.deserialize("Voxels", tempString);
    if(tempString.size() != numVoxels)
        LERROR("numVoxels mismatch! " << numVoxels << " / " << tempString.size());

    for(size_t i=0; i<voxels_.size(); i++) {
        if(tempString[i] == '1')
            voxels_[i] = true;
        else
            voxels_[i] = false;
    }

    ROISingle::deserialize(s);
}

bool ROIRaster::inROI(tgt::vec3 p) const {
    svec3 iVoxelCoords = physicalToVoxel(p);

    if(hand(tgt::greaterThanEqual(iVoxelCoords, svec3(size_t(0)))) && (hand(lessThan(iVoxelCoords, dims_))))
        return voxels_[calcPos(iVoxelCoords)];
    else
       return false;
}

void ROIRaster::setVoxel(const tgt::svec3& p, bool value) {
    if(hand(tgt::greaterThanEqual(p, svec3(size_t(0)))) && (hand(lessThan(p, dims_))))
        voxels_[calcPos(p)] = value;
}

bool ROIRaster::getVoxel(const tgt::svec3& p) const {
    if(hand(tgt::greaterThanEqual(p, svec3(size_t(0)))) && (hand(lessThan(p, dims_))))
        return voxels_[calcPos(p)];
    else
        return false;
}

tgt::ivec3 ROIRaster::worldToVoxel(const tgt::vec3& world) const {
    vec3 physical = getGrid().getWorldToPhysicalMatrix() * world;
    return physicalToVoxel(physical);
}

tgt::ivec3 ROIRaster::physicalToVoxel(const tgt::vec3& physical) const {
    vec3 voxelCoords = physical / getGrid().getSpacing();
    ivec3 iVoxelCoords = tgt::ifloor(voxelCoords);
    iVoxelCoords -= llf_;
    return iVoxelCoords;
}

Geometry* ROIRaster::generateMesh() const {
    return 0;
}

Geometry* ROIRaster::generateMesh(tgt::plane /*pl*/) const {
    return 0;
}

tgt::Bounds ROIRaster::getBoundingBox() const {
    vec3 llf = getGrid().getSpacing()*vec3(llf_);
    return tgt::Bounds(llf, llf+(vec3(dims_)*getGrid().getSpacing()));
}

std::vector<const ControlPoint*> ROIRaster::getControlPoints() const {
    std::vector<const ControlPoint*> cps = std::vector<const ControlPoint*>();
    //cps.push_back(new ControlPoint(this, ControlPoint::CENTER, getGrid().getPhysicalToWorldMatrix() * center_.get()));
    return cps;
}

std::vector<const ControlPoint*> ROIRaster::getControlPoints(tgt::plane /*pl*/) const {
    std::vector<const ControlPoint*> cps = std::vector<const ControlPoint*>();
    return cps;
}

bool ROIRaster::moveControlPoint(const ControlPoint* /*cp*/, tgt::vec3 /*offset*/) {
    return false;
}

void ROIRaster::setLLF(tgt::svec3 llf) {
    llf_ = llf;
    invalidate(ROIBase::ROI_CHANGE);
}

tgt::svec3 ROIRaster::getLLF() const {
    return llf_;
}

tgt::ivec3 ROIRaster::getDimensions() const {
    return dims_;
}

void ROIRaster::merge(const ROIRaster* m) {
    if(getGrid() != m->getGrid()) {
        float diff = getGrid().difference(m->getGrid());
        float th = 0.00000001f;
        if(diff > th) {
            LERROR("Failed to merge ROIRaster: Grids differ! " << diff);
            LERROR("Diff over threshold " << th);
            return;
        }
    }

    // Check if m is contained in this ROI:
    if(tgt::hand(tgt::greaterThanEqual(m->llf_, llf_)) && tgt::hand(lessThanEqual(m->llf_+ivec3(m->dims_), llf_+ivec3(dims_)))) {
        ivec3 offset = m->llf_ - llf_;
        svec3 i;
        for(i.x=0; i.x<m->dims_.x; i.x++) {
            for(i.y=0; i.y<m->dims_.y; i.y++) {
                for(i.z=0; i.z<m->dims_.z; i.z++) {
                    if(m->getVoxel(i))
                        setVoxel(offset+ivec3(i), true);
                }
            }
        }
    }
    else {
        // m is not in this => resize
        ivec3 llf = tgt::min(llf_, m->llf_);
        ivec3 urb = tgt::max(llf_+ivec3(dims_), m->llf_+ivec3(m->dims_));
        ivec3 dims = urb - llf;

        ROIRaster merged(getGrid(), llf, dims);

        ivec3 i;
        ivec3 offset = llf_ - llf;
        for(i.x=0; i.x<(int)dims_.x; i.x++) {
            for(i.y=0; i.y<(int)dims_.y; i.y++) {
                for(i.z=0; i.z<(int)dims_.z; i.z++) {
                    if(getVoxel(i))
                        merged.setVoxel(offset+i, true);
                }
            }
        }

        ivec3 offset2 = m->llf_ - llf;
        for(i.x=0; i.x<(int)m->dims_.x; i.x++) {
            for(i.y=0; i.y<(int)m->dims_.y; i.y++) {
                for(i.z=0; i.z<(int)m->dims_.z; i.z++) {
                    if(m->getVoxel(i))
                        merged.setVoxel(offset2+i, true);
                }
            }
        }

        llf_ = merged.llf_;
        dims_ = merged.dims_;
        voxels_ = merged.voxels_;
    }
}

void ROIRaster::erase(const ROIRaster* m) {
    if(getGrid() != m->getGrid()) {
        float diff = getGrid().difference(m->getGrid());
        float th = 0.00000001f;
        if(diff > th) {
            LERROR("Failed to merge ROIRaster: Grids differ! " << diff);
            LERROR("Diff over threshold " << th);
            return;
        }
    }

    //ivec3 llf = tgt::min(llf_, m->llf_);
    //ivec3 urb = tgt::max(llf_+ivec3(dims_), m->llf_+ivec3(m->dims_));
    //ivec3 dims = urb - llf;

    ivec3 i;
    //ivec3 offset = llf_ - llf;
    for(i.x=0; i.x<(int)dims_.x; i.x++) {
        for(i.y=0; i.y<(int)dims_.y; i.y++) {
            for(i.z=0; i.z<(int)dims_.z; i.z++) {
                if(getVoxel(i)) {
                    ivec3 voxelM = (i + llf_) - ivec3(m->getLLF());
                    if(m->getVoxel(voxelM))
                        setVoxel(i, false);
                }
            }
        }
    }
}

void ROIRaster::optimize() {
    ivec3 llf = ivec3(0);
    ivec3 urb = ivec3(0);
    bool voxelTrue = false;
    ivec3 i;
    for(i.x=0; i.x<(int)dims_.x; i.x++) {
        for(i.y=0; i.y<(int)dims_.y; i.y++) {
            for(i.z=0; i.z<(int)dims_.z; i.z++) {
                if(getVoxel(i)) {
                    if(voxelTrue) {
                        llf = min(llf, i);
                        urb = max(urb, i);
                        voxelTrue = true;
                    }
                    else {
                        //first voxel:
                        llf = i;
                        urb = i;
                        voxelTrue = true;
                    }
                }
            }
        }
    }

    llf += llf_;
    urb += llf_;

    ivec3 dims = (urb - llf) + ivec3(1);

    //LINFO("Optimizing from " << llf_ << " +> " << dims_ << " to " << llf << " +> " << dims << "(" << getName() << ")");

    if((svec3(dims) == dims_) && (llf == llf_))
        return; //nothing to optimize here...

    ROIRaster optimized(getGrid(), llf, dims, false);

    ivec3 offset = llf_ - llf;
    for(i.x=0; i.x<(int)dims_.x; i.x++) {
        for(i.y=0; i.y<(int)dims_.y; i.y++) {
            for(i.z=0; i.z<(int)dims_.z; i.z++) {
                if(getVoxel(i))
                    optimized.setVoxel(offset+i, true);
            }
        }
    }

    llf_ = optimized.llf_;
    dims_ = optimized.dims_;
    voxels_ = optimized.voxels_;
}

std::vector<vec3> computePCA(const std::vector<tgt::vec3>& points, vec3 center) {
    using namespace Eigen;

    unsigned int m = 3;              // dimension of each point
    unsigned int n = static_cast<unsigned int>(points.size());  // number of points

    MatrixXf DataPoints = MatrixXf::Zero(m, n);
    for (int i=0; i < DataPoints.cols(); i++) {
        for (int j=0; j < DataPoints.rows(); j++) {
            DataPoints(j, i) = points[i].elem[j] - center.elem[j];
        }
    }
    //std::cout << "Data:\n" << DataPoints << std::endl;

    // get the covariance matrix
    MatrixXf Covariance = MatrixXf::Zero(m, m);
    Covariance = (1.0 / (float) n) * DataPoints * DataPoints.transpose();
    //std::cout << "Covariance matrix:\n" << Covariance << "\n-----------------:\n";

    // compute the eigenvalue on the Cov Matrix
    EigenSolver<MatrixXf> m_solve(Covariance);
    VectorXf eigenvalues = VectorXf::Zero(m);
    eigenvalues = m_solve.eigenvalues().real();
    //std::cout << "Eigenvalues:\n" << eigenvalues << std::endl;

    MatrixXf eigenVectors = MatrixXf::Zero(n, m);  // matrix (n x m) (points, dims)
    eigenVectors = m_solve.eigenvectors().real();
    //std::cout << "Eigenvectors:\n"  << eigenVectors << std::endl;

    // sort and get the permutation indices
    typedef std::pair<double, int> myPair;
    typedef std::vector<myPair> PermutationIndices;
    PermutationIndices pi;
    for (int i=0 ; i<static_cast<int>(m); i++)
        pi.push_back(std::make_pair(eigenvalues(i), i));

    sort(pi.begin(), pi.end(), std::greater<myPair>());

    std::vector<vec3> evRet;
    for (unsigned int i=0; i<m ; i++) {
        int j = pi[i].second;
        vec3 ev = vec3(eigenVectors(0, j), eigenVectors(1, j), eigenVectors(2, j));
        ev = normalize(ev) * sqrtf(eigenvalues(j));
        std::cout << ev << std::endl;
        evRet.push_back(ev);
    }

    return evRet;
}

tgt::vec3 ROIRaster::getCenterOfMass() const {
    tgt::dvec3 com = tgt::dvec3(0.0);
    mat4 vToP = getVoxelToPhysicalMatrix();
    size_t numActive = 0;

    svec3 i;
    for(i.x=0; i.x<dims_.x; i.x++) {
        for(i.y=0; i.y<dims_.y; i.y++) {
            for(i.z=0; i.z<dims_.z; i.z++) {
                if(getVoxel(i)) {
                    tgt::dvec3 p = vToP * (vec3(i) + vec3(0.5f));
                    com += p;
                    numActive++;
                }
            }
        }
    }

    return (com / static_cast<double>(numActive));
}

std::vector<tgt::vec3> ROIRaster::getPrincipalComponents(tgt::vec3 center) const {
    mat4 vToP = getVoxelToPhysicalMatrix();

    svec3 i;
    std::vector<vec3> points;
    for(i.x=0; i.x<dims_.x; i.x++) {
        for(i.y=0; i.y<dims_.y; i.y++) {
            for(i.z=0; i.z<dims_.z; i.z++) {
                if(getVoxel(i)) {
                    tgt::vec3 p = vToP * (vec3(i) + vec3(0.5f));
                    points.push_back(p);
                }
            }
        }
    }

    return computePCA(points, center);
}

tgt::mat4 ROIRaster::getVoxelToPhysicalMatrix() const {
    return tgt::mat4::createScale(getGrid().getSpacing()) * tgt::mat4::createTranslation(vec3(llf_));
}

tgt::mat4 ROIRaster::getVoxelToWorldMatrix() const {
    return getGrid().getPhysicalToWorldMatrix() * getVoxelToPhysicalMatrix();
}

} // namespace
