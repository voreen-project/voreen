/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2018 University of Muenster, Germany,           *
 * Department of Computer Science.                                    *
 *                                                                    *
 * This file is part of the tgt library. This library is free         *
 * software; you can redistribute it and/or modify it under the terms *
 * of the GNU Lesser General Public License version 2.1 as published  *
 * by the Free Software Foundation.                                   *
 *                                                                    *
 * This library is distributed in the hope that it will be useful,    *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the       *
 * GNU Lesser General Public License for more details.                *
 *                                                                    *
 * You should have received a copy of the GNU Lesser General Public   *
 * License in the file "LICENSE.txt" along with this library.         *
 * If not, see <http://www.gnu.org/licenses/>.                        *
 *                                                                    *
 **********************************************************************/

#ifndef TGT_TEXTURE_H
#define TGT_TEXTURE_H

#include <string>
#include "tgt/tgt_gl.h"
#include "tgt/types.h"
#include "tgt/vector.h"

namespace tgt {

/**
 * Wrapper of an OpenGL Texture.
 *
 * Note: There is a difference between the cpu and the gpu data of a texture.
 * A new created texture has NO cpu data and alloc(...) has to be called. If the
 * cpu data is no longer needed, simply call setCpuTexture(0,false).
 */
class TGT_API Texture {
public:
    friend class TextureManager;


    //--------------------------------------------------------------------------------------
    // Enum definitions
    //--------------------------------------------------------------------------------------

    /** Enum for texture filtering */
    enum Filter {
        NEAREST,
        LINEAR,
        MIPMAP,
        ANISOTROPIC
    };

    /** Enum for the wrapping at texture boarders */
    enum Wrapping {
        REPEAT = GL_REPEAT,
        CLAMP = GL_CLAMP_TO_EDGE,
        CLAMP_TO_EDGE = GL_CLAMP_TO_EDGE,
        CLAMP_TO_BORDER = GL_CLAMP_TO_BORDER,
        MIRRORED_REPEAT = GL_MIRRORED_REPEAT
    };

    //--------------------------------------------------------------------------------------
    // Constructors and Destructor
    //--------------------------------------------------------------------------------------
private:
    /** do not use */
    explicit Texture() {}
    /** do not use */
    explicit Texture(const Texture& t) {}

public:
    /** Constructor */
    Texture(const tgt::ivec3& dimensions, GLint format = GL_RGBA, GLint internalformat = GL_RGBA,
            GLenum dataType = GL_UNSIGNED_BYTE, Filter filter = LINEAR, Wrapping wrapping = REPEAT,
            GLubyte* data = 0, bool ownership = true);

    /**
    * The destructor deletes the Texture in OpenGL.
    * Handled by the texturemanager!
    */
    virtual ~Texture();

    //--------------------------------------------------------------------------------------
    // Getter and Setter
    //--------------------------------------------------------------------------------------
    const GLuint getId() const { return id_; }
    GLenum getType() const { return type_; }
    tgt::ivec3 getDimensions() const { return dimensions_;}
    int getWidth() const { return dimensions_.x; };
    int getHeight() const { return dimensions_.y; }
    int getDepth() const { return dimensions_.z; }
    GLint getGLFormat() const { return format_; }
    GLint getGLInternalFormat() const { return internalformat_; }
    Filter getFilter() const { return filter_; }
    Wrapping getWrapping() const { return wrapping_; }
    GLenum getGLDataType() const { return dataType_; }
    GLubyte getBpt() const { return bpt_; }
    size_t getNumChannels() const { return calcNumChannels(format_); }
    size_t getSizeOnCPU() const;
    int getSizeOnGPU() const;
    std::string getOptionalName() const { return name_; }
    GLubyte* getCpuTextureData() { return cpuTextureData_; }
    const GLubyte* getCpuTextureData() const { return cpuTextureData_; }
    bool getUseMipMap() const { return useMipMap_; }

    /** Sets the optinal name */
    void setOptionalName(const std::string& name) { name_ = name; }
    /** Updates the dimension. Can not change 1D, 2D or 3D. Old cpuData will get lost. No new cpuData is allocted. */
    void updateDimensions(const ivec3& dimensions, bool reuploadTexture);
    /** Sets the filtering for the texture. Binds the texture!!! When changing the filter to MIPMAP afterwards, you have to re-upload the texture!!! */
    void setFilter(Filter filter);
    /** Set texture wrapping mode. Binds the texture!!! */
    void setWrapping(Wrapping w);
    /** Sets how the values of RGBA appear to the shader.
     *  @param swizzleMask Defines the mapping. Example: GLenum swizzleMask[4] = {GL_ZERO, GL_ZERO, GL_ZERO, GL_RED}; */
    void setSwizzle(const GLint* swizzleMask);
    /** Updates the cpu texture data. If the texture has ownership, the previous data will be destroyed. */
    void setCpuTextureData(GLubyte* data, bool ownership);

    // removed because the flag is set implicitly by specifying the filtering
    /** Determines if the texture should use a mip map or not. You have to call upload texture to re-upload with the changes. */
    //void setUseMipMap(bool use) { useMipMap_ = use; }


    //--------------------------------------------------------------------------------------
    // Wrapper of OpenGL functions
    //--------------------------------------------------------------------------------------
protected:
    /** Generates a new OpenGL texture ID */
    void generateId();
public:
    /** Bind the texture to the active texture unit and target.
     *  Note: This does not enable texturing (use enable()). */
    void bind() const;
    /** Enable texturing on the active texture unit. */
    void enable() const;
    /** Disable texturing on the active texture unit. */
    void disable() const;
    /** Allocates an appropriate buffer for the texture on the CPU.
    *   Texture gets ownership of the buffer and deletes previous buffers. */
    GLubyte* alloc(bool ownership);

    //--------------------------------------------------------------------------------------
    // Upload and download texture
    //--------------------------------------------------------------------------------------

    /**
     * Upload Texture to graphics-card. Binds the texture.
     *
     * type_, format_, internalformat_, dimensions, dataType_ and the pixels_ pointer have to
     * be set before calling this method.
     */
    void uploadTexture();

    /**
     * Download Texture from graphics-card. Binds the texture.
     *
     * type_, format_, dimensions, dataType_ and the pixels_ pointer have to be set before
     * calling this method! Pixels will be allocated if pixels_ is a nullpointer.
     */
    void downloadTexture() const;

    /**
     * Download texture from the GPU to a newly allocated buffer, to which a
     * pointer is returned.  Binds the texture.
     *
     * type_, format_, dimensions, and dataType_ have to be set before
     * calling this method!
     */
    GLubyte* downloadTextureToBuffer() const;

    /**
     * Download texture from the GPU to a preallocated buffer. Binds the texture.
     *
     * type_, format_, dimensions, and dataType_ have to be set before
     * calling this method!
     */
     void downloadTextureToBuffer(GLubyte* pixels, size_t numBytesAllocated) const;

    /**
     * Download texture from the GPU to a newly allocated buffer with
     * the passed format/data type and the texture's dimensions.
     */
    GLubyte* downloadTextureToBuffer(GLint format, GLenum dataType) const;

    /**
     * Download texture from the GPU to a preallocated buffer with
     * the passed format/data type and the texture's dimensions.
     * Binds the texture.
     *
     * type_ and dimensions have to be set before calling this method!
     */
    void downloadTextureToBuffer(GLint format, GLenum dataType, GLubyte* pixels,
                                          size_t numBytesAllocated) const;

    //-----------------------------------------------------------------------------------------
    // Texel access (only if downloadTexture has been called or cpuTextureData_ is up to date)
    //-----------------------------------------------------------------------------------------
    template <class T>
    inline T& texel(size_t x, size_t y = 0, size_t z = 0) {
        return texel<T>(tgt::svec3(x, y, z));
    }
    template <class T>
    inline const T& texel(size_t x, size_t y = 0, size_t z = 0) const {
        return texel<T>(tgt::svec3(x, y, z));
    }
    template <class T>
    inline T& texel(const svec2& pos) {
        return texel<T>(tgt::svec3(pos, 0));
    }
    template <class T>
    inline const T& texel(const svec2& pos) const {
        return texel<T>(tgt::svec3(pos, 0));
    }
    template <class T>
    inline T& texel(const svec3& pos) {
        tgtAssert( sizeof(T) == bpt_, "sizeof(T) != bytes per texel here" );
        tgtAssert( type_ == GL_TEXTURE_3D || (type_ == GL_TEXTURE_2D && pos.z == 0) || (type_ == GL_TEXTURE_1D && pos.z == 0 && pos.y == 0), "Accessing more dimensions than available");
        tgtAssert( pos.z*dimensions_.x*dimensions_.y + pos.y*dimensions_.x + pos.x < size_t(hmul(dimensions_)), "index out of range" );
        tgtAssert( cpuTextureData_, "no cpu texture data. Call download texture or set texture data first.");
        return ((T*) cpuTextureData_)[pos.z*dimensions_.x*dimensions_.y + pos.y*dimensions_.x + pos.x];
    }
    template <class T>
    inline const T& texel(const svec3& pos) const {
        tgtAssert( sizeof(T) == bpt_, "sizeof(T) != bytes per texel here" );
        tgtAssert( type_ == GL_TEXTURE_3D || (type_ == GL_TEXTURE_2D && pos.z == 0) || (type_ == GL_TEXTURE_1D && pos.z == 0 && pos.y == 0), "Accessing more dimensions than available");
        tgtAssert( pos.z*dimensions_.x*dimensions_.y + pos.y*dimensions_.x + pos.x < size_t(hmul(dimensions_)), "index out of range" );
        tgtAssert( cpuTextureData_, "no cpu texture data. Call download texture or set texture data first.");
        return ((T*) cpuTextureData_)[pos.z*dimensions_.x*dimensions_.y + pos.y*dimensions_.x + pos.x];
    }

    ///Return texel as tgt::Color (slow!), downloadTexture() needs to be called first
    tgt::Color texelAsFloat(size_t x, size_t y=0) const;
    tgt::Color texelAsFloat(tgt::svec2 p) const { return texelAsFloat(p.x, p.y); }

    ///Sets texel as tgt::Color (slow!), downloadTexture() needs to be called first
    void texelFromFloat(const tgt::Color& color, size_t x, size_t y=0);
    void texelFromFloat(const tgt::Color& color, tgt::svec2 p) { texelFromFloat(color, p.x, p.y); }


    //--------------------------------------------------------------------------------------
    // Helper functions (mostly static)
    //--------------------------------------------------------------------------------------
    /** Calculates the bytes per texel from format and dataType */
    static int calcFormatBpt(GLint format, GLenum dataType);
    /** Calculates the bytes per texel from the internal format */
    static int calcInternalFormatBpt(GLint internalformat);
    /** Calculates the number of channels from the passed texture format */
    static int calcNumChannels(GLint format);
    /** Calculates the type_ (GL_TEXTURE_1D, GL_TEXTURE_2D, GL_TEXTURE_3D) from dimensions */
    static GLenum calcType(const ivec3& dimensions);

    //--------------------------------------------------------------------------------------
    // Members
    //--------------------------------------------------------------------------------------
protected:
    tgt::ivec3 dimensions_; ///< dimension of the texture
    GLint format_;          ///< GL_RGB...
    GLint internalformat_;  ///< GL_RGB...
    GLenum dataType_;       ///< GL_UNSIGNED_BYTE
    Filter filter_;         ///< current filter
    Wrapping wrapping_;     ///< current wrapping

    GLubyte* cpuTextureData_;   ///< pointer to the buffer containing the for loading image
    bool ownsCpuTextureData_;   ///< bool flag, if the texture owns the cpu data

    GLuint id_;             ///< OpenGL texture id (generated in the constructor)
    GLenum type_;           ///< 1D, 2D, 3D (calculated in the constructor)
    GLubyte bpt_;           ///< bytes per texel (calculated in the constructor)

    bool useMipMap_;        ///< when this is set to true, uploading the texture will call glGenerateMipMap

    std::string name_;      ///< optional, e.g. for storing texture file name
};

} // namespace tgt

#endif // TGT_TEXTURE_H
