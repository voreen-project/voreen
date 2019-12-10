/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#ifndef VRN_GEOMETRYSEQUENCE_H
#define VRN_GEOMETRYSEQUENCE_H

#include "voreen/core/datastructures/geometry/geometry.h"

namespace voreen {

class VRN_CORE_API GeometrySequence : public Geometry {
public:
    GeometrySequence(bool takeOwnership_ = true);
    virtual ~GeometrySequence();

    virtual VoreenSerializableObject* create() const { return new GeometrySequence(); }

    virtual std::string getClassName() const { return "GeometrySequence"; }

    virtual void addGeometry(Geometry* g);
    //TODO: remove

    size_t getNumGeometries() const;
    const Geometry* getGeometry(size_t i) const;

    virtual bool equals(const Geometry* geometry, double epsilon = 1e-5) const;

    /**
     * Clips the GeometrySequence against the passed clipping plane.
     *
     * The clipping plane is represented by the equation
     * NX - c = 0 where N is the normal vector of the clipping plane
     * and corresponds to the x-, y-, and z-coordinate of the given vector,
     * and c is a scalar corresponding to the w-coordinate of the given
     * vector and X is any point on the plane.
     *
     * This function is supposed to be overridden by subclasses.
     * The default implementation is a no-op.
     *
     * @note It is not necessary that the clipping plane normal vector
     *       is normalized, but it may improve computational accuracy.
     *
     * @param clipPlane an arbitrary clipping plane
     * @param epsilon maximum distance at which two vertices are to be considered equal
     */
    virtual void clip(const tgt::plane& clipPlane, double epsilon = 1e-5);

    virtual void render() const;

    virtual tgt::Bounds getBoundingBox(bool transformed = true) const;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);
private:
    std::vector<Geometry*> geometries_;
    bool takeOwnership_;

    static const std::string loggerCat_;
};

} // namespace

#endif  //VRN_GEOMETRYSEQUENCE_H
