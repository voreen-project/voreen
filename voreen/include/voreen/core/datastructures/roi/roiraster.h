/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#ifndef VRN_ROIRASTER_H
#define VRN_ROIRASTER_H

#include "voreen/core/datastructures/roi/roisingle.h"

namespace voreen {

class ROIRaster : public ROISingle {
public:
    ROIRaster(Grid grid, tgt::ivec3 llf, tgt::svec3 dims, bool value = false);
    /// Raster roi in given grid.
    ROIRaster(Grid grid, const ROIBase* roi);
    /// Raster roi in own grid.
    ROIRaster(const ROIBase* roi);
    ROIRaster();
    virtual ROIBase* create() const { return new ROIRaster(); }
    virtual std::string getClassName() const { return "ROIRaster"; }

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

    virtual bool inROI(tgt::vec3 p) const;

    virtual Geometry* generateMesh() const;
    virtual Geometry* generateMesh(tgt::plane pl) const;

    virtual tgt::Bounds getBoundingBox() const;

    virtual std::vector<const ControlPoint*> getControlPoints() const;
    virtual std::vector<const ControlPoint*> getControlPoints(tgt::plane pl) const;

    virtual bool moveControlPoint(const ControlPoint* cp, tgt::vec3 to);

    tgt::ivec3 worldToVoxel(const tgt::vec3& world) const;
    tgt::ivec3 physicalToVoxel(const tgt::vec3& physical) const;

    tgt::mat4 getVoxelToPhysicalMatrix() const;
    tgt::mat4 getVoxelToWorldMatrix() const;

    void setVoxel(const tgt::svec3& p, bool value);
    bool getVoxel(const tgt::svec3& p) const;

    void setLLF(tgt::svec3 llf);
    tgt::svec3 getLLF() const;
    tgt::ivec3 getDimensions() const;

    void merge(const ROIRaster* m);
    void erase(const ROIRaster* m);
    /// Shrinks ROI as far as possible
    void optimize();

    /// Computes the center of mass (in physical coordinates)
    tgt::vec3 getCenterOfMass() const;
    std::vector<tgt::vec3> getPrincipalComponents(tgt::vec3 center) const;
private:
    size_t calcPos(size_t x, size_t y, size_t z) const {
        return z*dims_.x*dims_.y + y*dims_.x + x;
    }

    size_t calcPos(const tgt::svec3& pos) const {
        return pos.z*dims_.x*dims_.y + pos.y*dims_.x + pos.x;
    }

    tgt::ivec3 llf_;
    tgt::svec3 dims_;
    std::vector<bool> voxels_;
    static const std::string loggerCat_;
};

} //namespace

#endif
