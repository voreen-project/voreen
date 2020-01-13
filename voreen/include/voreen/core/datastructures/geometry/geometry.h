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

#ifndef VRN_GEOMETRY_H
#define VRN_GEOMETRY_H

#include "voreen/core/voreenobject.h"
#include "voreen/core/datastructures/meta/metadatacontainer.h"
#include "voreen/core/datastructures/cloneable.h"

#include "tgt/vector.h"
#include "tgt/matrix.h"
#include "tgt/bounds.h"

#include <string>
#include <vector>

namespace voreen {


/**
 * Abstract base class for Geometry objects that
 * can be passed through GeometryPorts.
 *
 * @note In order to enable serialization (necessary for caching),
 *  a factory has to be provided for each concrete subtype.
 *
 * @see GeometryFactory
 * @see VoreenModule::registerSerializerFactory
 */
class VRN_CORE_API Geometry : public VoreenSerializableObject, public Cloneable<Geometry> {
public:
    Geometry();
    virtual ~Geometry();

    /**
     * Returns a copy of the geometry object.
     *
     * The default implementation uses XML serialization/deserialization
     * for cloning. For performance and accuracy reasons, however,
     * the cloning should be implemented directly.
     */
    virtual std::unique_ptr<Geometry> clone() const;

    /**
     * Returns whether the geometry is equal to the passed one.
     *
     * This function is supposed to be overridden by subclasses.
     * The default implementation compares the geometries' hashes.
     *
     * @param geometry Geometry object to compare to
     * @param epsilon maximum distance at which two vertices are to be considered equal
     */
    virtual bool equals(const Geometry* geometry, double epsilon = 1e-5) const;

    /// Transforms the geometry by multiplying the given transformation matrix to the existing one.
    virtual void transform(const tgt::mat4& m);

    tgt::mat4 getTransformationMatrix() const;
    tgt::mat4 getInvertedTransformationMatrix() const;

    void setTransformationMatrix(const tgt::mat4& m);

    /**
     * Clips the geometry against the passed clipping plane.
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

    /**
     * Renders the geometry.
     *
     * This function is supposed to be overridden by subclasses.
     * The default implementation is a no-op.
     */
    virtual void render() const;

    /**
     * Returns the axis-aligned bounding box of the geometry.
     *
     * This function is supposed to be overridden by subclasses.
     * The default implementation returns undefined bounds.
     *
     * @param transformed Apply transformation matrix?
     */
    virtual tgt::Bounds getBoundingBox(bool transformed = true) const;

    /**
     * Returns a hash of the geometry object,
     * which is used for the caching mechanism.
     *
     * The default implementation computes an MD5 hash
     * of the serialized XML string. Subclasses may
     * override this method, but are not required to do so.
     */
    virtual std::string getHash() const;

    /**
     * Indicates whether this geometry or parts of it are transparent.
     *
     * By default a geometry a geometry is considered opaque.
     * Subclasses may override this method,
     * but are not required to do so.
     */
    virtual bool isTransparent() const;

    /**
     * Returns the meta data container of this geometry.
     * Can be used to store texture names, for instance.
     *
     * @see MetaDataContainer
     */
    MetaDataContainer& getMetaDataContainer() const;

    /**
     * Supposed to be overridden by each concrete subclass.
     */
    virtual void serialize(Serializer& s) const;

    /**
     * Supposed to be overridden by each concrete subclass.
     */
    virtual void deserialize(Deserializer& s);

protected:
    static const std::string loggerCat_;
private:
    tgt::mat4 transformationMatrix_;    ///< Usually model to world

    /**
     * Contains the associated meta data.
     *
     */
    mutable MetaDataContainer metaDataContainer_;
};

} // namespace

#endif  //VRN_GEOMETRY_H
