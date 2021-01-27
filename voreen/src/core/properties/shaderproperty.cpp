/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/utils/stringutils.h"

#include "voreen/core/utils/voreenfilepathhelper.h"

#include "tgt/filesystem.h"
#include "tgt/shadermanager.h"

namespace voreen {

ShaderFileList::ShaderFileList(tgt::ShaderObject::ShaderType type, const std::string& name)
    : files_()
{
    if(!name.empty()) {
        files_.emplace_back(type, name);
    }
}

ShaderFileList& ShaderFileList::operator()(tgt::ShaderObject::ShaderType type, const std::string& name) {
    return add(type, name);
}

ShaderFileList& ShaderFileList::add(tgt::ShaderObject::ShaderType type, const std::string& name) {
    if(!name.empty()) {
        files_.emplace_back(type, name);
    }
    return *this;
}

/// ------------------------------------------------------------------------
/// ShaderSourceComponent
/// ------------------------------------------------------------------------

static std::string shaderTypeToString(tgt::ShaderObject::ShaderType type) {
    switch(type) {
        case tgt::ShaderObject::ShaderType::VERTEX_SHADER:
            return "Vertex shader";
        case tgt::ShaderObject::ShaderType::FRAGMENT_SHADER:
            return "Fragment shader";
        case tgt::ShaderObject::ShaderType::GEOMETRY_SHADER:
            return "Geometry shader";
        case tgt::ShaderObject::ShaderType::COMPUTE_SHADER:
            return "Compute shader";
        default:
            tgtAssert(false, "Unimplemented shader type");
            return "";
    }
}

const std::string ShaderSourceComponent::loggerCat_("voreen.ShaderSourceComponent");

ShaderSourceComponent::ShaderSourceComponent()
    : type_(tgt::ShaderObject::ShaderType::UNDEFINED)
    , originalFilename_("")
    , externalFilename_("")
    , source_("")
    , modified_(false)
{
}
ShaderSourceComponent::ShaderSourceComponent(tgt::ShaderObject::ShaderType type, const std::string& name)
    : type_(type)
    , originalFilename_(name)
    , externalFilename_("")
    , source_("")
    , modified_(false)
{
}

std::string ShaderSourceComponent::getOriginalFileName() const {
    return originalFilename_;
}

std::string ShaderSourceComponent::getCurrentFileName() const {
    return isExternal() ? externalFilename_ : originalFilename_;
}
std::string ShaderSourceComponent::getCurrentCompleteFilePath() const {
    if(isExternal()) {
        return externalFilename_;
    } else {
        tgtAssert(tgt::Singleton<tgt::ShaderManager>::isInited(), "ShaderManager not initiated.");
        return ShdrMgr.completePath(originalFilename_);
    }
}
std::string ShaderSourceComponent::getDescription() const {
    std::string description = "";
    description += shaderTypeToString(type_);
    description += " (";
    description += getCurrentFileName();
    if(isExternal()) {
        description += " [external]";
    }
    description += ")";
    return description;
}

bool ShaderSourceComponent::operator==(const ShaderSourceComponent& other) const {
    if(    originalFilename_ != other.originalFilename_
        || externalFilename_ != other.externalFilename_
        || type_ != other.type_) {
        return false;
    }
    if(isModified() || other.isModified()) {
        return false;
    }
    return source_ == other.source_;
}

bool ShaderSourceComponent::operator!=(const ShaderSourceComponent& other) const {
    return !(*this == other);
}

void ShaderSourceComponent::reloadSource() {
    std::string filename = getCurrentFileName();
    tgtAssert(!filename.empty(), "filename empty");

    if(!tgt::Singleton<tgt::ShaderManager>::isInited() && !isExternal()) {
        LWARNING("ShaderManager not instantiated");
        source_ = "";
        return;
    }

    if(!tgt::Singleton<tgt::FileSystem>::isInited()) {
        LWARNING("FileSystem not instantiated");
        source_ = "";
        return;
    }

    std::string completeFilename = getCurrentCompleteFilePath();

    std::unique_ptr<tgt::File> file = std::unique_ptr<tgt::File>(FileSys.open(completeFilename));

    // check if file is open
    if(!file || !file->isOpen()) {
        LERROR("File not found: " << filename);
        source_ = "";
        return;
    }

#ifdef _WIN32
    source_ = convertNewlinesWindowsToUnix(file->getAsString());
#else
    source_ = file->getAsString();
#endif

    modified_ = false;

    file->close();
}
void ShaderSourceComponent::reset() {
    externalFilename_ = "";
    reloadSource();
}

void ShaderSourceComponent::setExternalFilename(const std::string& externalFilename) {
    externalFilename_ = externalFilename;
    reloadSource();
}
bool ShaderSourceComponent::isExternal() const {
    return !externalFilename_.empty();
}
bool ShaderSourceComponent::isModified() const {
    return modified_;
}
tgt::ShaderObject::ShaderType ShaderSourceComponent::getType() const {
    return type_;
}
const std::string& ShaderSourceComponent::getSource() const {
    return source_;
}
void ShaderSourceComponent::setSource(const std::string& source) {
    if(source != source_) {
        source_ = source;
        modified_ = true;
    }
}

// Make the compiler catch typos
static const char* SERIALIZATION_IDENTIFIER_SHADER_TYPE = "shaderType";
static const char* SERIALIZATION_IDENTIFIER_ORIGINAL_FILENAME = "originalFilename";
static const char* SERIALIZATION_IDENTIFIER_MODIFIED = "modified";
static const char* SERIALIZATION_IDENTIFIER_EXTERNAL_FILENAME_PATHS = "externalFilenamePaths";
static const char* SERIALIZATION_IDENTIFIER_SOURCE = "source";

void ShaderSourceComponent::serialize(Serializer& s) const {
    s.serialize(SERIALIZATION_IDENTIFIER_SHADER_TYPE, type_);

    s.serialize(SERIALIZATION_IDENTIFIER_ORIGINAL_FILENAME, originalFilename_);

    s.serialize(SERIALIZATION_IDENTIFIER_MODIFIED, modified_);

    // Only write external file name if there is one
    if(isExternal()) {
        s.serialize(SERIALIZATION_IDENTIFIER_EXTERNAL_FILENAME_PATHS, VoreenFilePathHelper(externalFilename_));
    }

    // Only write source if it has been modified
    if(modified_) {
        s.serialize(SERIALIZATION_IDENTIFIER_SOURCE, source_);
    }
}

void ShaderSourceComponent::deserialize(Deserializer& s) {
    try {
        unsigned int type;
        s.deserialize(SERIALIZATION_IDENTIFIER_SHADER_TYPE, type);
        type_ = static_cast<tgt::ShaderObject::ShaderType>(type);

        s.deserialize(SERIALIZATION_IDENTIFIER_ORIGINAL_FILENAME, originalFilename_);

        s.deserialize(SERIALIZATION_IDENTIFIER_MODIFIED, modified_);

        // Try to get external file name
        try {
            VoreenFilePathHelper pathHelper;
            s.deserialize(SERIALIZATION_IDENTIFIER_EXTERNAL_FILENAME_PATHS, pathHelper);
            externalFilename_ = pathHelper.getPath();
        } catch (SerializationNoSuchDataException&) {
            // No external file name => no external file
            s.removeLastError();
            externalFilename_ = "";
        }

        // If the source is not modified we can (and have to, see serialize) reload it from file instead
        if(modified_) {
            s.deserialize(SERIALIZATION_IDENTIFIER_SOURCE, source_);
        } else {
            reloadSource();
        }
    } catch (VoreenException& e) {
        LERROR(std::string("Error during deserialization: ") + e.what());
    }
}
// This is ugly, but because of VS2012 we cannot use initializer lists...
static std::string makeCamelCase(const std::string& first, const std::string& second, const std::string& third = "") {
    std::ostringstream output;
    output << first;
    if(!output.str().empty() && !second.empty()) {
        output << static_cast<unsigned char>(std::toupper(second[0]));
        output << second.substr(1);
    }
    if(!output.str().empty() && !third.empty()) {
        output << static_cast<unsigned char>(std::toupper(third[0]));
        output << third.substr(1);
    }
    return output.str();
}

void ShaderSourceComponent::deserializeFallback(Deserializer& s) {
    std::string typeName;
    switch(getType()) {
        case tgt::ShaderObject::VERTEX_SHADER:
            typeName = "vertex";
            break;
        case tgt::ShaderObject::GEOMETRY_SHADER:
            typeName = "geometry";
            break;
        case tgt::ShaderObject::FRAGMENT_SHADER:
            typeName = "fragment";
            break;
        case tgt::ShaderObject::COMPUTE_SHADER:
            typeName = ""; //? TODO check
            break;
        default:
            tgtAssert(false, "Unimplemented shader type");
    }
    try {
        s.deserialize(makeCamelCase(typeName, "modified"), modified_);
        bool isExternal;
        s.deserialize(makeCamelCase(typeName, "isExternal"), isExternal);

        if (isExternal) {
            VoreenFilePathHelper tmp;
            try {
                s.deserialize(makeCamelCase("external", typeName, "filenamePaths"), tmp);
                externalFilename_ = tmp.getPath();
            }
            catch (SerializationNoSuchDataException&) {
                //old deserialization
                s.removeLastError();
                s.deserialize(makeCamelCase("external", typeName, "filename"), externalFilename_);
                externalFilename_ = tgt::FileSystem::absolutePath(tgt::FileSystem::dirName(s.getDocumentPath()) + "/" + externalFilename_);
            }
            catch (VoreenException& ) {
                externalFilename_ = "";
            }
        } else {
            externalFilename_ = "";
        }

        // If the source is not modified we need to reload it from file instead
        if(modified_) {
            s.deserialize(makeCamelCase(typeName, "source"), source_);
        } else {
            reloadSource();
        }

    } catch (VoreenException& e) {
        LERROR(std::string("Error during (fallback) deserialization: ") + e.what());
    }
}

/// ------------------------------------------------------------------------
/// ShaderSource
/// ------------------------------------------------------------------------

const std::string ShaderSource::loggerCat_("voreen.ShaderSource");

ShaderSource::ShaderSource()
    : components_()
{
}
/*
ShaderSource::ShaderSource(const ShaderSource& other)
    : ShaderSource()
{
    for(auto& entry : other.components_) {
        components_.push_back(entry);
    }
}
*/
ShaderSource::ShaderSource(const ShaderFileList& fileList)
{
    for(auto& file : fileList.files_) {
        components_.emplace_back(file.first, file.second);
    }
}

ShaderSourceComponent& ShaderSource::operator[](size_t i) {
    return at(i);
}
ShaderSourceComponent& ShaderSource::at(size_t i) {
    tgtAssert(i < components_.size(), "Invalid index");
    return components_[i];
}

const ShaderSourceComponent& ShaderSource::operator[](size_t i) const {
    return at(i);
}

const ShaderSourceComponent& ShaderSource::at(size_t i) const {
    tgtAssert(i < components_.size(), "Invalid index");
    return components_[i];
}

ShaderSource::iterator ShaderSource::begin() {
    return components_.begin();
}

ShaderSource::iterator ShaderSource::end() {
    return components_.end();
}

ShaderSource::const_iterator ShaderSource::cbegin() const {
    return components_.cbegin();
}

ShaderSource::const_iterator ShaderSource::cend() const {
    return components_.cend();
}

size_t ShaderSource::numComponents() const {
    return components_.size();
}

bool ShaderSource::operator==(const ShaderSource& other) const {
    if(components_.size() != other.components_.size()) {
        return false;
    }
    for(size_t i=0; i<components_.size(); ++i) {
        if(components_[i] != other.components_[i]) {
            return false;
        }
    }
    return true;
}
bool ShaderSource::operator!=(const ShaderSource& other) const {
    return !(*this == other);
}
void ShaderSource::serialize(Serializer& s) const {
   s.serialize("components", components_);
}

void ShaderSource::deserialize(Deserializer& s) {
    bool newSerializationFailed = false;
    try {
        // First step: Deserialize all components into temporary vector
        std::vector<ShaderSourceComponent> deserializedComponents;
        s.deserialize("components", deserializedComponents);

        // Second step: Iterate over all deserialized components...
        for(auto& deserializedComponent: deserializedComponents) {
            // ... and find the corresponding component specified in the constructor.
            auto originalComponent = std::find_if(components_.begin(), components_.end(), [deserializedComponent] (ShaderSourceComponent& s) {
                        return s.getType() == deserializedComponent.getType() && s.getOriginalFileName() == deserializedComponent.getOriginalFileName();
                    });
            if(originalComponent != components_.end()) {
                // If there is such a component we can finish the deserialization for it
                *originalComponent = deserializedComponent;
            } else {
                // If not, we issue a warning and ignore that item.
                LWARNING("No " << shaderTypeToString(deserializedComponent.getType()) << " \"" << deserializedComponent.getOriginalFileName() << "\" found for deserialization.");
            }
        }
    } catch (VoreenException& ) {
        s.removeLastError();
        newSerializationFailed = true;
    }
    // If the new serialization was successfull, we are done and can return
    if(!newSerializationFailed) {
        return;
    }
    // Otherwise we try the fallback deserialization for a previous implementation, but issue some sanity checks first.
    if(std::count_if(components_.begin(), components_.end(),
                [] (const ShaderSourceComponent& c) { return c.getType() == tgt::ShaderObject::VERTEX_SHADER; }) > 1) {
        LWARNING("Fallback deserialization not applicable: Too many vertex shaders.");
        return;
    }
    if(std::count_if(components_.begin(), components_.end(),
                [] (const ShaderSourceComponent& c) { return c.getType() == tgt::ShaderObject::GEOMETRY_SHADER; }) > 1) {
        LWARNING("Fallback deserialization not applicable: Too many geometry shaders.");
        return;
    }
    if(std::count_if(components_.begin(), components_.end(),
                [] (const ShaderSourceComponent& c) { return c.getType() == tgt::ShaderObject::FRAGMENT_SHADER; }) > 1) {
        LWARNING("Fallback deserialization not applicable: Too many fragment shaders.");
        return;
    }
    if(std::count_if(components_.begin(), components_.end(),
                [] (const ShaderSourceComponent& c) { return c.getType() == tgt::ShaderObject::COMPUTE_SHADER; }) > 1) {
        LWARNING("Fallback deserialization not applicable: Too many compute shaders.");
        return;
    }

    for(ShaderSourceComponent& component : components_) {
        component.deserializeFallback(s);
    }
}

//---------------------------------------------------------------------------------------------------------------

const std::string ShaderProperty::loggerCat_("voreen.ShaderProperty");

ShaderProperty::ShaderProperty(const std::string& id, const std::string& guiText,
        const std::string& fragmentFilename,
        const std::string& vertexFilename,
        const std::string& geometryFilename,
        int invalidationLevel,
        Property::LevelOfDetail lod)
        : TemplateProperty<ShaderSource>(id, guiText,
            ShaderSource(ShaderFileList
                (tgt::ShaderObject::VERTEX_SHADER, vertexFilename)
                (tgt::ShaderObject::GEOMETRY_SHADER, geometryFilename)
                (tgt::ShaderObject::FRAGMENT_SHADER, fragmentFilename)),
            invalidationLevel, lod)
        , header_("")
        , shader_(nullptr)
        , requiresRebuild_(true)
{
}

ShaderProperty::ShaderProperty(const std::string& id, const std::string& guiText,
        ShaderFileList shaderFiles,
        int invalidationLevel,
        Property::LevelOfDetail lod)
    : TemplateProperty<ShaderSource>(id, guiText, ShaderSource(shaderFiles), invalidationLevel, lod)
    , header_("")
    , shader_(nullptr)
    , requiresRebuild_(true)
{}

ShaderProperty::ShaderProperty()
    : header_("")
    , shader_(nullptr)
{}

ShaderProperty::~ShaderProperty() {
    if (shader_) {
        LWARNING(getFullyQualifiedGuiName() << " has not been deinitialized before destruction.");
    }
}

Property* ShaderProperty::create() const {
    return new ShaderProperty();
}

void ShaderProperty::initialize() {
    TemplateProperty<ShaderSource>::initialize();

    shader_ = std::unique_ptr<tgt::Shader>(new tgt::Shader());
    LGL_ERROR;

    auto mut = getMutator();
    for(ShaderSourceComponent& component : *mut) {
        if(!component.isModified()) {
            component.reloadSource();
        }
    }
}

void ShaderProperty::deinitialize() {
    shader_ = nullptr;
    LGL_ERROR;

    TemplateProperty<ShaderSource>::deinitialize();
}

void ShaderProperty::invalidate() {
    requiresRebuild_ = true;
    TemplateProperty<ShaderSource>::invalidate();
}

void ShaderProperty::serialize(Serializer& s) const {
    Property::serialize(s);
    get().serialize(s);
}

void ShaderProperty::deserialize(Deserializer& s) {
    Property::deserialize(s);
    ShaderSource n = get();
    n.deserialize(s);
    set(n);
    invalidate();
    updateWidgets();
}

void ShaderProperty::setHeader(std::string header) {
    header_ = header;
    requiresRebuild_ = true;
}

std::string ShaderProperty::getHeader() const {
    return header_;
}


bool ShaderProperty::rebuild() {
    bool allSuccessful = true;

    // We can not simply delete the old shader since it might has been referenced using the getShader()
    // method which would lead to dangling pointers. Therefore, initialize once and detach each shader object.
    while (shader_->hasObjects()) {
        tgt::ShaderObject* object = shader_->getObjects().back();
        shader_->detachObject(object);
        delete object;
    }

    for(auto& entry : value_) {
        std::string filename = entry.getCurrentCompleteFilePath();
        bool validFilename = true;

        if(filename.empty() || !tgt::FileSystem::fileExists(filename)) {
            LWARNING(entry.getDescription() + " not found.");
            validFilename = false;
        }

        // should we allow the use of a modified source if the file does not exist in the file-system?
        if(validFilename || entry.isModified()) {
            tgt::ShaderObject* obj = new tgt::ShaderObject(filename, entry.getType());
            obj->setHeader(header_);
            obj->setSource(entry.getSource());

            obj->compileShader();
            shader_->attachObject(obj);

            if (!obj->isCompiled()) {
                allSuccessful = false;
                LWARNING("Failed to compile " + entry.getDescription() + ": " << obj->getCompilerLog());
            }
        }
    }

    if(allSuccessful) {
        if(shader_->hasObjects()) {
            shader_->linkProgram();
            if(!shader_->isLinked()) {
                LWARNING("Failed to link shader: " << shader_->getLinkerLog());
                allSuccessful = false;
            }
        } else {
            LWARNING("No objects have been attached to shader.");
            allSuccessful = false;
        }
    }

    requiresRebuild_ = false;

    updateWidgets();

    return allSuccessful;
}

tgt::Shader* ShaderProperty::getShader() const {
    return shader_.get();
}

bool ShaderProperty::hasValidShader() const {
    return (shader_ && shader_->isLinked());
}

bool ShaderProperty::requiresRebuild() const {
    return requiresRebuild_;
}

}   // namespace
