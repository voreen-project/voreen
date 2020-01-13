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

#ifndef VRN_TRANSFUNC1D_H
#define VRN_TRANSFUNC1D_H

#include "voreen/core/datastructures/transfunc/transfuncbase.h"

#include "tgt/vector.h"

#include "voreen/core/datastructures/transfunc/1d/preintegrationtablemap.h"

namespace voreen {

class PreIntegrationTable;

/**
 * One dimensional transfer function. Supports pre-integration tables.
 *
 * Internally, it is represented by a one-dimensional RGBA texture of type GL_UNSIGNED_BYTE.
 * It is used to handel the gamma value, domain and threshold.
 */
class VRN_CORE_API TransFunc1D : public TransFuncBase {

    friend class PreIntegrationTable;

public:
    /**
     * Constructor
     *
     * @param width desired width of the transfer function
     */
    TransFunc1D(int width = 256, DataType dataType = TF_UBYTE, tgt::Texture::Filter filter = tgt::Texture::NEAREST);

    /**
     * Destructor - deletes the keys of the transfer function
     */
    virtual ~TransFunc1D();

    /** @overridde */
    virtual void setMemberValuesFrom(const TransFuncBase* transfunc);
    /** @overridde */
    virtual bool compareTo(const TransFuncBase& tf) const;
    /** @overridde */
    virtual int getNumDimensions() const { return 1;}

    //--------------------------------------
    //  texture handling
    //--------------------------------------
protected:
    /** @override */
    virtual void invalidateTexture();
    /** @overridde */
    virtual GLvoid* createTextureData();

    /** Has to be implemented in the sub-classes. */
    virtual tgt::Vector4<GLfloat> getMappingForValueFloat(float x) = 0;
    /** Has to be implemented in the sub-classes. */
    virtual tgt::Vector4<GLubyte> getMappingForValueUByte(float x) = 0;

    //--------------------------------------
    //  pre-integration handling
    //--------------------------------------
public:
    /**
     * Returns the pre-integration table of the transfer function with the specified attributes.
     * If necessary, the table is (re-)computed.
     *
     * @param dimension width of the pre-integration table, 0 chooses the width according to the bit depth of the volume (up to 1024)
     * @param samplingStepSize the segment length used for rendering, in texture (normalized) coordinates
     * @param useIntegral @see PreIntegrationTable
     * @param computeOnGPU @see PreIntegrationTable
     */
    const PreIntegrationTable* getPreIntegrationTable(float samplingStepSize = 1.f, size_t dimension = 0, bool useIntegral = true, bool computeOnGPU = false);
    //--------------------------------------
    //  gamma handling
    //--------------------------------------
public:
    /** sets the gamma value*/
    void setGammaValue(float gamma);
    /** retuns the alpha mode*/
    float getGammaValue() const;
public:
    /** applies gamma to texture data */
    float applyGammaToIntensity(float value) const;

    float gammaValue_;            ///< value for gamma correction

    //--------------------------------------
    //  domain handling
    //--------------------------------------
public:
    /**
     * Sets the transfer function's domain, i.e., the intensity range it covers,
     * for the specified dimension.
     */
    void setDomain(const tgt::vec2& domain);
    /** @overload */
    void setDomain(float lower, float upper);

    /**
     * Returns the transfer function's domain, i.e., the intensity range it covers,
     * for the specified dimension.
     */
    tgt::vec2 getDomain() const;

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
    void setThreshold(float lower, float upper);
    /** @overload */
    void setThreshold(const tgt::vec2& thresholds);

    /**
     * Returns the lower and upper intensity thresholds of the tranfer function.
     * The thresholds are normalized within the range [0,1].
     */
    tgt::vec2 getThreshold() const;

    //--------------------------------------
    //  shader defines
    //--------------------------------------
public:
    /** @overridde */
    virtual std::string getSamplerType() const { return "sampler1D";}
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
    float realWorldToNormalized(float rw) const;

    /**
     * Converts the passed normalized data value (range: [0.0,1.0]) to the corresponding real-world value,
     * with regard to the currently set transfer function domain.
     *
     * @param n the normalized data value to convert
     * @param dimension of the transfer function dimension to apply the mapping for
     */
    float normalizedToRealWorld(float n) const;

protected:
    tgt::vec2 domain_;                  ///< domain range
    tgt::vec2 threshold_;               ///< thresholdrange

    PreIntegrationTableMap preIntegrationTableMap_; ///< contains several pre-integration tables for the tf
    mutable tgt::Shader* preIntegrationProgram_; ///< shader program to compute pre-integration tables on the gpu

    static const std::string loggerCat_; ///< logger category
};

} // namespace voreen

#endif // VRN_TRANSFUNC1D_H
