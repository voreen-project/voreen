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

#ifndef VRN_SHADERPROPERTY_H
#define VRN_SHADERPROPERTY_H

#include "voreen/core/properties/templateproperty.h"
#include "tgt/shadermanager.h"

namespace voreen {

struct ShaderFileList {
    ShaderFileList(tgt::ShaderObject::ShaderType type, const std::string& name);
    ShaderFileList& operator() (tgt::ShaderObject::ShaderType type, const std::string& name);
    ShaderFileList& add(tgt::ShaderObject::ShaderType type, const std::string& name);

    std::vector<std::pair<tgt::ShaderObject::ShaderType, std::string>> files_;
};

class VRN_CORE_API ShaderSourceComponent : public Serializable {
    //We only want the ShaderProperty to update the obj_ file after compiling
    friend class ShaderProperty;

public:
    ShaderSourceComponent(); // For deserialization
    ShaderSourceComponent(tgt::ShaderObject::ShaderType type, const std::string& name);
    //ShaderSourceComponent(const ShaderSourceComponent& file) = default;
    //ShaderSourceComponent& operator=(const ShaderSourceComponent& file) = default;

    bool isExternal() const;
    bool isModified() const;
    tgt::ShaderObject::ShaderType getType() const;
    const std::string& getSource() const;
    void setSource(const std::string& source);

    void setExternalFilename(const std::string& externalFilename);

    /**
     * Get the _original_ file name
     */
    std::string getOriginalFileName() const;

    /**
     * Get the current file name (external or original)
     */
    std::string getCurrentFileName() const;

    /**
     * Get the complete path to the current file (external or original)
     * Note: May be empty if an interal source cannot be found
     */
    std::string getCurrentCompleteFilePath() const;

    /**
     * Get a short description of the current file, including:
     * type, internal/external, filename
     */
    std::string getDescription() const;

    /**
     * Reload the source from the current file (external or original)
     */
    void reloadSource();

    /**
     * Resets to original file
     */
    void reset();

    /*
     * Check if both file names and shader types in both entries match and if none are modified.
     * (According to previous implementation.)
     */
    bool operator==(const ShaderSourceComponent& other) const;
    bool operator!=(const ShaderSourceComponent& other) const;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);
    void deserializeFallback(Deserializer& s);

protected:

    const static std::string loggerCat_;

    tgt::ShaderObject::ShaderType type_;
    std::string originalFilename_;
    std::string externalFilename_;

    std::string source_;

    // Modified from the original file
    bool modified_;
};

class VRN_CORE_API ShaderSource : public Serializable {
public:
    typedef std::vector<ShaderSourceComponent>::iterator iterator;
    typedef std::vector<ShaderSourceComponent>::const_iterator const_iterator;
    ShaderSource();
    //ShaderSource(const ShaderSource& other);
    ShaderSource(const ShaderFileList& fileList);

    ShaderSourceComponent& operator[](size_t i);
    ShaderSourceComponent& at(size_t i);
    const ShaderSourceComponent& operator[](size_t i) const;
    const ShaderSourceComponent& at(size_t i) const;

    iterator begin();
    iterator end();
    const_iterator cbegin() const;
    const_iterator cend() const;

    size_t numComponents() const;

    /**
     * Operator to compare two ShaderSource objects.
     * Currently, the order of entries does matter. Should it?
     *
     * @param shaderSource ShaderSource instance that is compared to this instance
     * @return true when both ShaderSource instances are equal
     */
    bool operator==(const ShaderSource& other) const;
    /**
     * Operator to compare two ShaderSource objects.
     *
     * @param shaderSource ShaderSource instance that is compared to this instance
     * @return true when both ShaderSource instances are not equal
     */
    bool operator!=(const ShaderSource& shaderSource) const;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

protected:
    const static std::string loggerCat_;

    std::vector<ShaderSourceComponent> components_;
};


class VRN_CORE_API ShaderProperty : public TemplateProperty<ShaderSource> {
public:
    ShaderProperty(const std::string& id, const std::string& guiText,
                   const std::string& fragmentFileName,
                   const std::string& vertexFileName = "",
                   const std::string& geometryFileName = "",
                   int invalidationLevel=Processor::INVALID_PROGRAM,
                   Property::LevelOfDetail lod = Property::LOD_DEFAULT);

    ShaderProperty(const std::string& id, const std::string& guiText,
                   ShaderFileList shaderFiles,
                   int invalidationLevel=Processor::INVALID_PROGRAM,
                   Property::LevelOfDetail lod = Property::LOD_DEFAULT);

    ShaderProperty();
    virtual ~ShaderProperty();

    virtual Property* create() const;

    virtual std::string getClassName() const       { return "ShaderProperty"; }
    virtual std::string getTypeDescription() const { return "Shader"; }

    void initialize();
    void deinitialize();

    virtual void invalidate();

    /**
     * @see Property::serialize
     */
    virtual void serialize(Serializer& s) const;

    /**
     * @see Property::deserialize
     */
    virtual void deserialize(Deserializer& s);

    void setHeader(std::string header);
    std::string getHeader() const;

    tgt::Shader* getShader() const;
    bool rebuild();
    bool hasValidShader() const;

    /**
     * Tells whether the shader obtainable by getShader matches the source files and the header.
     */
    bool requiresRebuild() const;

private:

    std::string header_;

    std::unique_ptr<tgt::Shader> shader_;

    bool requiresRebuild_;

    static const std::string loggerCat_;
};
}   // namespace

#endif
