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

#ifndef VRN_TRANSFUNCBASE_H
#define VRN_TRANSFUNCBASE_H

#include "voreen/core/io/serialization/serialization.h"
#include "voreen/core/datastructures/meta/metadatabase.h"
#include "voreen/core/voreenobject.h"

#include "voreen/core/voreencoreapi.h"
#include "tgt/texture.h"
#include "tgt/vector.h"

#include <vector>
#include <string>

namespace tgt {
    class Shader;
}

namespace voreen {

/**
 * Base class for transfer functions.
 *
 * The format is always GL_RGBA and the data type GL_UNSIGNED_BYTE or GL_FLOAT
 *
 * The base class is responsable for the texture generation ans the alpha mode.
 */
class VRN_CORE_API TransFuncBase : public VoreenSerializableObject {

    /// used to call resize()
    friend class TransFunc1DKeysProperty;
    friend class TransFunc1DGaussianProperty;
    friend class TransFunc2DPrimitivesProperty;

public:
    enum DataType {
        TF_UBYTE = GL_UNSIGNED_BYTE,
        TF_FLOAT = GL_FLOAT
    };

    /**
     * Constructor.
     *
     * @param width width of transfer function
     * @param height of transfer function. Pass one for a 1D transfer function.
     * @param depth of transfer function. Pass one for a 1D or 2D transfer function.
     * @param dataType internel data format
     * @param filter transfer function texture's filtering mode (see tgt::Texture)
     */
    TransFuncBase(int width = 1024, int height = 1, int depth = 1, DataType dataTyp = TF_UBYTE, tgt::Texture::Filter filter = tgt::Texture::NEAREST);

    /**
     * Destructor.
     */
    virtual ~TransFuncBase();

     /** Returns the dimensions of the transfer function's texture. */
    tgt::ivec3 getDimensions() const;
    /** Copies all member values into this. */
    virtual void setMemberValuesFrom(const TransFuncBase* transfunc);
    /** Compares all members */
    virtual bool compareTo(const TransFuncBase& tf) const;
    /** returns the dimensions of the TF. */
    virtual int getNumDimensions() const = 0;
    /** clones this */
    virtual TransFuncBase* clone() const = 0;
    /** Resets the transfer function to its default value. */
    virtual void reset() = 0;
protected:
    //--------------------------------------
    //  handle texture
    //--------------------------------------
public:
    /**
     * Returns the texture of the transfer function. If it has been marked invalid,
     * updateTexture() is called before.
     *
     * @note If the texture is not present when calling this function, it is created.
     *    Therefore, the caller has to make sure that a valid OpenGL context is active.
     */
    tgt::Texture* getTexture() const;

    /**
     * Returns the data type.
     */
    DataType getDataType() const;

    /** Returns if the current texture is valid. */
    virtual bool isTextureValid() const;
protected:
    /**
     * Resizes the transfer function.
     */
    void resize(int width, int height, int depth);

    /**
     * Sets textureInvalid to true.
     */
    virtual void invalidateTexture();

    /**
     * Implement in subclass
     *
     * @ return tex_ gets ownership of pixeldata
     */
    virtual GLvoid* createTextureData() = 0;

    /**
     * Generates the transfer function texture according to the specified parameters.
     */
    void createTex();

    /**
     * Updates the texture of the transfer function or creates it, if it is not present.
     * The base class implementation uploads the texture data to the GPU.
     *
     * @note calls apply (alpha/gamma/threshold)
     *
     * @note This function is called automatically on each texture access, when the texture is marked as invalid.
     *    You might call it directly in order to force an immediate update. But then you have to make sure that
     *    an valid OpenGL context is active.
     *
     */
    void updateTexture();

protected:
    tgt::Texture* tex_;           ///< the texture of the transfer function, is generated internally
    tgt::ivec3 dimensions_;       ///< dimensions of the transfer function texture
    DataType dataType_;           ///< data format of each texel
    tgt::Texture::Filter filter_; ///< filtering mode of the transfer function texture.
    bool textureInvalid_;         ///< indicates whether the transfer function texture has to be updated

    //--------------------------------------
    //  handle alpha
    //--------------------------------------
public:
    /**
     * Enum to determine, if the alpha value should be used or is constant zero or one.
     */
    enum AlphaMode {
        TF_ZERO_ALPHA,
        TF_USE_ALPHA,
        TF_ONE_ALPHA
    };
    /** sets the alpha mode*/
    void setAlphaMode(AlphaMode mode);
    /** returns the alpha mode*/
    AlphaMode getAlphaMode() const;
protected:
    /** applies alpha to texture data */
    void applyAlpha(GLvoid* data) const;

    AlphaMode alphaMode_;         ///< mode to determine the alpha value use.

    //--------------------------------------
    //  load and save
    //--------------------------------------
public:
    /**
     * Creates a transfer function out of the data contained in the file given by filename.
     *
     * @param filename The path to the file in which the data is stored
     * @return true when the load was successful, false otherwise
     */
    virtual bool load(const std::string& filename) = 0;

 /**
     * Saves the transfer function to a file. Any data in the file will be overwritten.
     * The supported extensions include:
     *
     * @param filename the name of the file the transfer function will be saved to
     * @return true, if the operation was successfull, false otherwise
     */
    virtual bool save(const std::string& filename) const = 0;

    /**
     * Returns a vector that contains the endings of suppported file formats for loading.
     *
     * @return vector with endings of supported file formats
     */
    virtual const std::vector<std::string> getLoadFileFormats() const = 0;

    /**
     * Returns a vector that contains the endings of supported file formats for saving.
     *
     * @return vector with endings of supported file formats
     */
    virtual const std::vector<std::string> getSaveFileFormats() const = 0;

    /**
     * Dummy. The lookup table is not serialized.
     *
     * @see Serializable::serialize
     */
    virtual void serialize(Serializer& s) const;

    /**
     * Dummy. The lookup table is not deserialized.
     *
     * @see Serializable::deserialize
     */
    virtual void deserialize(Deserializer& s);

    //--------------------------------------
    //  shader defines
    //--------------------------------------
public:
    /**
     * Returns a define for the usage of transfer functions in shaders.
     * For 1D transfer functions the define looks like: "#define TF_SAMPLER_TYPE sampler1D \n"
     * and with sampler2D for 2D transfer functions.
     *
     * @return define for usage of transfer functions in shaders
     */
    virtual std::string getShaderDefines(const std::string& defineName = "TF_SAMPLER_TYPE") const;

    /**
     * Sets the default uniforms for shders useing transfer functions.
     */
    virtual void setUniform(tgt::Shader* shader, const std::string& uniform, const std::string& uniformTex, const GLint texUnit) = 0;

    /**
     * Returns a string representation of the sampler type: "sampler1D" for 1D transfer
     * functions, "sampler2D" for 2D transfer functions, "sampler3D" for 3D transfer functions
     *
     * @return string representation of the sampler type used by the transfer function
     */
    virtual std::string getSamplerType() const = 0;

    //--------------------------------------
    //  real world mapping helper
    //--------------------------------------
protected:
    /**
     * Converts the passed normalized data value (range: [0.0,1.0]) to the corresponding real-world value,
     * with regard to the passed transfer function domain.
     *
     * @param n the normalized data value to convert
     * @param domain the transfer function domain to use for the transformation. domain.x must be smaller than domain.y.
     */
    static float normalizedToRealWorld(float n, const tgt::vec2& domain);
    /**
     * Converts the passed real-world data value to a normalized value in the range [0.0,1.0],
     * with regard to the passed transfer function domain.
     *
     * @param rw the real-world data value to normalize
     * @param domain transfer function domain to use. domain.x must be smaller than domain.y.
     */
    static float realWorldToNormalized(float rw, const tgt::vec2& domain);

    //--------------------------------------
    //  operators
    //--------------------------------------
public:
    bool operator==(const TransFuncBase& tf) { return  compareTo(tf); }
    bool operator!=(const TransFuncBase& tf) { return !compareTo(tf); }

protected:
    static const std::string loggerCat_; ///< logger category
};



///------------------------------------------------------------------------------------------------

/**
 * Metadata encapsulating a TransferFunction.
 * The meta data object takes ownership of the encapsulated transfer function
 * and deletes it on its own destruction.
 */
template<typename T>
class TransFuncMetaDataGeneric : public MetaDataBase {
public:
    TransFuncMetaDataGeneric();
    TransFuncMetaDataGeneric(T* transfunc);
    TransFuncMetaDataGeneric(const std::vector<T*>& transfunc);
    virtual ~TransFuncMetaDataGeneric();

    void setTransferFunction(T* transfunc, size_t channel = 0);
    void setTransferFunction(const std::vector<T*>& transfunc);

    T* getTransferFunction(size_t channel = 0) const;
    const std::vector<T*>& getValue() const;

    virtual std::string toString() const = 0;
    virtual std::string toString(const std::string& component) const;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

    virtual size_t getNumChannels() const;

protected:
    std::vector<T*> transFunc_;
};

template<typename T>
TransFuncMetaDataGeneric<T>::TransFuncMetaDataGeneric()
    : MetaDataBase()
{}

template<typename T>
TransFuncMetaDataGeneric<T>::TransFuncMetaDataGeneric(T* transfunc)
    : MetaDataBase()
{
    transFunc_.push_back(transfunc);
}

template<typename T>
TransFuncMetaDataGeneric<T>::TransFuncMetaDataGeneric(const std::vector<T*>& transfunc)
    : MetaDataBase()
{
    transFunc_.assign(transfunc.begin(),transfunc.end());
}

template<typename T>
TransFuncMetaDataGeneric<T>::~TransFuncMetaDataGeneric() {
    for(size_t i = 0; i < transFunc_.size(); i++)
        delete transFunc_[i];
}

template<typename T>
void TransFuncMetaDataGeneric<T>::setTransferFunction(T* transfunc, size_t channel) {
    if(channel < transFunc_.size()) {
        delete transFunc_[channel];
        transFunc_[channel] = transfunc;
    } else {
        while (channel >= transFunc_.size())
            transFunc_.push_back(0);
        transFunc_[channel] = transfunc;
    }
}

template<typename T>
void TransFuncMetaDataGeneric<T>::setTransferFunction(const std::vector<T*>& transfunc) {
    for(size_t i = 0; i < transFunc_.size(); i++)
        delete transFunc_[i];
    transFunc_.assign(transfunc.begin(),transfunc.end());
}

template<typename T>
T* TransFuncMetaDataGeneric<T>::getTransferFunction(size_t channel) const {
    if (channel < transFunc_.size())
        return transFunc_[channel];
    else
        return 0;
}

template<typename T>
const std::vector<T*>& TransFuncMetaDataGeneric<T>::getValue() const {
    return transFunc_;
}

template<typename T>
std::string TransFuncMetaDataGeneric<T>::toString(const std::string& /*component*/) const {
    return this->toString();
}

template<typename T>
void TransFuncMetaDataGeneric<T>::serialize(Serializer& s) const {
    if (!transFunc_.empty()) {
        try {
            s.serialize("transfunc", transFunc_);
        }
        catch (SerializationException& e) {
            LERRORC("voreen.TransFuncMetaData", std::string("Failed to serialize transfunc: ") + e.what());
        }
    }
}

template<typename T>
void TransFuncMetaDataGeneric<T>::deserialize(Deserializer& s) {
    for(size_t i = 0; i < transFunc_.size(); i++)
        delete transFunc_[i];
    transFunc_.clear();
    try {
        s.deserialize("transfunc", transFunc_);
    }
    catch (SerializationException& e) {
        LERRORC("voreen.TransFuncMetaData", std::string("Failed to deserialize transfunc: ") + e.what());
    }
}

template<typename T>
size_t TransFuncMetaDataGeneric<T>::getNumChannels() const {
    return transFunc_.size();
}

} // namespace voreen

#endif // VRN_TRANSFUNCBASE_H
