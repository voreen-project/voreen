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

#ifndef VRN_TRANSFUNC2D_H
#define VRN_TRANSFUNC2D_H

#include "voreen/core/datastructures/transfunc/transfuncbase.h"

namespace voreen {

/**
 * Two dimensional transfer function.
 *
 * Internally, it is represented by a two-dimensional RGBA texture of type GL_FLOAT.
 * It is used to handel the gamma value, domain and threshold.
 */
class VRN_CORE_API TransFunc2D : public TransFuncBase {
public:
    /**
     * Constructor
     *
     * @param width width of the texture of the transfer function
     * @param height height of the texture of the transfer function
     */
    TransFunc2D(int width = 256, int height = 256, DataType dataType = TF_FLOAT, tgt::Texture::Filter filter = tgt::Texture::LINEAR);

    /** Destructor */
    virtual ~TransFunc2D();

    /** @overridde */
    virtual void setMemberValuesFrom(const TransFuncBase* transfunc);
    /** @overridde */
    virtual bool compareTo(const TransFuncBase& tf) const;
    /** @overridde */
    virtual int getNumDimensions() const { return 2;}

    //--------------------------------------
    //  texture handling
    //--------------------------------------
protected:
    /** @overridde */
    virtual GLvoid* createTextureData();

    /** Has to be implemented in the sub-classes. */
    virtual tgt::Vector4<GLfloat> getMappingForValueFloat(float x, float y) = 0;
    /** Has to be implemented in the sub-classes. */
    virtual tgt::Vector4<GLubyte> getMappingForValueUByte(float x, float y) = 0;

    //--------------------------------------
    //  gamma handling
    //--------------------------------------
public:
    /** sets the gamma value*/
    void setGammaValue(tgt::vec2 gamma);
    /** retuns the alpha mode*/
    tgt::vec2 getGammaValue() const;
public:
    /** applies gamma to texture data */
    tgt::vec2 applyGammaToValues(float x, float y) const;

    tgt::vec2 gammaValues_;            ///< value for gamma correction

    //--------------------------------------
    //  domain handling
    //--------------------------------------
public:
    /**
     * Sets the transfer function's domain, i.e., the intensity range it covers,
     * for the specified dimension.
     */
    void setDomain(const tgt::vec2& domain, size_t dimension);
    /** @overload */
    void setDomain(float lower, float upper, size_t dimension);

    /**
     * Returns the transfer function's domain, i.e., the intensity range it covers,
     * for the specified dimension.
     */
    tgt::vec2 getDomain(size_t dimension) const;

    //--------------------------------------
    //  threshold handling
    //--------------------------------------
    /**
     * Sets the lower and upper intensity thresholds to given values. The thresholds have to be normalized
     * to the range [0,1]. The texture is not updated at this time.
     *
     * @param lower lower threshold
     * @param upper upper threshold
     *
     */
    void setThreshold(float lower, float upper, size_t dimension);
    /** @overload */
    void setThreshold(const tgt::vec2& thresholds, size_t dimension);

    /**
     * Returns the lower and upper intensity thresholds of the tranfer function.
     * The thresholds are normalized within the range [0,1].
     */
    tgt::vec2 getThreshold(size_t dimension) const;

    //--------------------------------------
    //  shader defines
    //--------------------------------------
public:
    /** @overridde */
    virtual std::string getSamplerType() const { return "sampler2D";}
    /** @overridde */
    virtual void setUniform(tgt::Shader* shader, const std::string& uniform, const std::string& uniformTex, const GLint texUnit);

    //--------------------------------------
    //  load and save
    //--------------------------------------
public:
    /** @see Serializable::serialize */
    virtual void serialize(Serializer& s) const;
    /** @see Serializable::deserialize  */
    virtual void deserialize(Deserializer& s);

    //--------------------------------------
    //  real world mapping
    //--------------------------------------
public:
    /**
     * Converts the passed real-world data value to a normalized value in the range [0.0,1.0],
     * with regard to the currently set domain.
     *
     * @param rw the real-world data value to normalize
     * @param dimension of the transfer function dimension to apply the mapping for
     */
    float realWorldToNormalized(float rw, int dimension = 0) const;
    /** @overload */
    tgt::vec2 realWorldToNormalized(tgt::vec2 rw) const;

    /**
     * Converts the passed normalized data value (range: [0.0,1.0]) to the corresponding real-world value,
     * with regard to the currently set transfer function domain.
     *
     * @param n the normalized data value to convert
     * @param dimension of the transfer function dimension to apply the mapping for
     */
    float normalizedToRealWorld(float n, int dimension) const;
    /** @overload */
    tgt::vec2 normalizedToRealWorld(tgt::vec2 n) const;

protected:
    tgt::vec2 domains_[2];                  ///< domain range
    tgt::vec2 thresholds_[2];               ///< thresholdrange
    static const std::string loggerCat_;    ///< the logger category
};

} // namespace voreen

#endif // VRN_TRANSFUNC2D_H
