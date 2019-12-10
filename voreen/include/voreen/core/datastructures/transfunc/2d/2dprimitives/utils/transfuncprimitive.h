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

#ifndef VRN_TRANSFUNCPRIMITIVE_H
#define VRN_TRANSFUNCPRIMITIVE_H

#include "voreen/core/io/serialization/serialization.h"

#include "tgt/vector.h"
#include "tgt/tgt_gl.h"

namespace voreen {

/**
 * Control point for a primitive. Is composed of a position and color (including opacity)
 */
class VRN_CORE_API TransFuncPrimitiveControlPoint : public VoreenSerializableObject {
public:
    tgt::vec2 position_;
    tgt::col4 color_;

    TransFuncPrimitiveControlPoint(tgt::vec2 pos = tgt::vec2(0.f), tgt::col4 color = tgt::col4(255))
        : position_(pos)
        , color_(color)
    {}

    virtual std::string getClassName() const   { return "TransFuncPrimitiveControlPoint";     }
    virtual TransFuncPrimitiveControlPoint* create() const { return new TransFuncPrimitiveControlPoint(); }

    tgt::vec4 getColorNormalized() {
        return tgt::vec4(color_) / tgt::vec4(255.f);
    }

    /**
     * @see serializable::serialize
     */
    virtual void serialize(Serializer& s) const {
        s.serialize("position", position_);
        s.serialize("color", color_);
    }

    /**
     * @see serializable::deserialize
     */
    virtual void deserialize(Deserializer& s) {
        s.deserialize("position", position_);
        s.deserialize("color", color_);
    }
};

/**
 * Abstract base class for all primitives that are used in 2D transfer functions.
 */
class VRN_CORE_API TransFuncPrimitive : public VoreenSerializableObject {
public:
    TransFuncPrimitive();

    /**
     * Paints the primitive.
     */
    virtual void paint() = 0;

    /**
     * Returns the number of control points this primitive consists of.
     */
    size_t getNumControlPoints() const;

    /**
     * Sets the color of the primitive (i.e. all control points!) to the given value.
     *
     * @param c desired color for the primitive
     * @note invalidateTexture() must be called from the transfer function to update the change
     */
    void setColor(const tgt::col4& c);

    /**
     * Returns a reference to the control point.
     *
     * @param index index of the control point, has to be within this range: 0 <= index < getNumControlPoints()
     */
    const TransFuncPrimitiveControlPoint& getControlPoint(size_t index) const;
    TransFuncPrimitiveControlPoint& getControlPoint(size_t index);

    /**
     * Sets the fuzziness of the primitive to the given value. With increasing fuzziness the
     * primitives get an alphagradient towards their border.
     *
     * @param f desired fuzziness
     * @note invalidateTexture() must be called from the transfer function to update the change
     */
    void setFuzziness(float f);

    /**
     * Returns the current fuzziness of the primitive.
     *
     * @return current fuzziness
     */
    float getFuzziness() const;

    /**
     * Returns the distance between pos and the closest control point.
     *
     * @param pos position the distance of control points is measured to
     * @param distance between pos and closest control point
     */
    virtual float getClosestControlPointDist(const tgt::vec2& pos) = 0;

    /**
     * Moves all of the primitive's control points by the given offset (without checking the coordinate range).
     *
     * @param offset offset the coordinates are moved by
     * @note invalidateTexture() must be called from the transfer function to update the change
     */
    virtual void move(const tgt::vec2& offset);

    /**
     * @see serializable::serialize
     */
    virtual void serialize(Serializer& s) const;

    /**
     * @see serializable::deserialize
     */
    virtual void deserialize(Deserializer& s);

    /**
     * Returns a copy of this object.
     */
    virtual TransFuncPrimitive* clone() const = 0;

protected:

    std::vector<TransFuncPrimitiveControlPoint> controlPoints_;     ///< the control points of the primitive
    float fuzziness_;   ///< fuzziness of the primitive
};

// ----------------------------------------------------------------------------

/**
 * A triangle primitive. It consists of 3 vertices that can be moved independently.
 */
class VRN_CORE_API TransFuncTriangle : public TransFuncPrimitive {
public:
    TransFuncTriangle();

    /**
     * Constructor
     *
     * @param a  coordinate of the first vertex
     * @param b  coordinate of the second vertex
     * @param c  coordinate of the third vertex
     * @param col color of the primitive
     */
    TransFuncTriangle(const tgt::vec2& a, const tgt::vec2& b, const tgt::vec2& c, const tgt::col4& col);

    virtual std::string getClassName() const   { return "TransFuncTriangle";    }
    virtual TransFuncTriangle* create() const { return new TransFuncTriangle();     }

    /**
     * Paints the triangle. The fuzziness factor is obeyed.
     */
    void paint();

    /**
     * Returns the distance between pos and the closest control point.
     *
     * @param pos position the distance of control points is measured to
     * @param distance between pos and closest control point
     */
    float getClosestControlPointDist(const tgt::vec2& pos);

    /**
     * Returns a copy of this object.
     */
    virtual TransFuncPrimitive* clone() const;

};


/**
 * A quad primitive. It consists of 4 vertices that can be moved independently.
 */
class VRN_CORE_API TransFuncQuad : public TransFuncPrimitive {
public:
    TransFuncQuad();

    /**
     * Constructor
     *
     * @param center center of the quad
     * @param size size of the quad
     * @param col color of the quad
     */
    TransFuncQuad(const tgt::vec2& center, float size, const tgt::col4& col);

    virtual std::string getClassName() const   { return "TransFuncQuad";    }
    virtual TransFuncQuad* create() const { return new TransFuncQuad();     }

    /**
     * Paints the quad. The fuzziness factor is obeyed.
     */
    void paint();

    /**
     * Returns the distance between pos and the closest control point.
     *
     * @param pos position the distance of control points is measured to
     * @param distance between pos and closest control point
     */
    float getClosestControlPointDist(const tgt::vec2& pos);

    /**
     * Returns a copy of this object.
     */
    virtual TransFuncPrimitive* clone() const;

};

// ----------------------------------------------------------------------------

/**
 * A banana primitive. It consists of 4 vertices that are connected by 2 splines.
 * The control points are arranged in the following way.
 *       1
 *
 *       2
 * 0          3
 * The first spline starts at point 0, goes through 1 and ends in 3.
 * The second spline starts at point 0, goes through 2 and ends in 3.
 */
class VRN_CORE_API TransFuncBanana : public TransFuncPrimitive {
public:
    TransFuncBanana();

    /**
     * Constructor
     *
     * @param a  coordinate of the left control point
     * @param b1 coordinate of the upper middle control point
     * @param b2 coordinate of the lower middle control point
     * @param c  coordinate of the right control point
     * @param col color of the primitive
     */
    TransFuncBanana(const tgt::vec2& a, const tgt::vec2& b1, const tgt::vec2& b2, const tgt::vec2& c, const tgt::col4& col);

    virtual std::string getClassName() const   { return "TransFuncBanana";  }
    virtual TransFuncBanana* create() const { return new TransFuncBanana(); }

    /**
     * Paints the banana. The fuzziness is obeyed.
     */
    void paint();

    /**
     * Returns the distance between pos and the closest control point.
     *
     * @param pos position the distance of control points is measured to
     * @param distance between pos and closest control point
     */
    float getClosestControlPointDist(const tgt::vec2& pos);

    /**
     * @see Serializable::serialize
     */
    virtual void serialize(Serializer& s) const;

    /**
     * @see Serializable::deserialize
     */
    virtual void deserialize(Deserializer& s);

    /**
     * Returns a copy of this object.
     */
    virtual TransFuncPrimitive* clone() const;

protected:
    /**
     * Paints the space between the both splines. steps_ triangles in a trianglestrip are used for that.
     */
    void paintInner();

    size_t steps_; ///< number of triangles used to fill the space between the both splines

};

} // namespace voreen

#endif // VRN_TRANSFUNCPRIMITIVE_H
