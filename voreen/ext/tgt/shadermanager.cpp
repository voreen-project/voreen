/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2024 University of Muenster, Germany,           *
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

#include "tgt/shadermanager.h"

#include "tgt/gpucapabilities.h"
#include "tgt/texturemanager.h"

#include <iostream>
#include <fstream>
#include <memory>

using std::string;

namespace tgt {

//------------------------------------------------------------------------------

namespace {

} // namespace

const std::string GLSLExtensionSet::loggerCat_("tgt.Shader.GLSLExtensionSet");

GLSLExtensionSet::GLSLExtensionSet()
    : extensions_()
{
}
void GLSLExtensionSet::set(const std::string& name, ExtensionStatus status) {
    auto extension = extensions_.find(name);
    if(extension == extensions_.end()) {
        extensions_[name] = status;
    } else {
        ExtensionStatus currentStatus = extensions_[name];
        if(currentStatus == status) {
            return;
        }

        int combinedStatus = currentStatus | status;
        if((combinedStatus & REQUIRED) && (combinedStatus & ENABLED)) {
            extensions_[name] = REQUIRED;
            return;
        }
        if((combinedStatus & REQUIRED) && (combinedStatus & WARN)) {
            LWARNING("Extension " << name << " is required, but also set as warn.");
            extensions_[name] = REQUIRED;
            return;
        }
        if((combinedStatus & REQUIRED) && (combinedStatus & DISABLED)) {
            LERROR("Extension " << name << " is both required and enabled!");
            extensions_[name] = REQUIRED;
            return;
        }
        if((combinedStatus & ENABLED) && (combinedStatus & WARN)) {
            extensions_[name] = WARN;
            return;
        }
        if((combinedStatus & ENABLED) && (combinedStatus & DISABLED)) {
            LWARNING("Extension " << name << " is both disabled and enabled. It will be disabled.");
            extensions_[name] = DISABLED;
            return;
        }
        if((combinedStatus & WARN) && (combinedStatus & DISABLED)) {
            LWARNING("Extension " << name << " is disabled, but also set als warn. It will be disabled.");
            extensions_[name] = DISABLED;
            return;
        }
    }
}
string GLSLExtensionSet::format(const std::string& name, ExtensionStatus estatus) const {
    std::string status;
    switch(estatus) {
        case ENABLED:
            status = "enable";
            break;
        case REQUIRED:
            status = "require";
            break;
        case WARN:
            status = "warn";
            break;
        case DISABLED:
            status = "disable";
            break;
    }
    return "#extension " + name + " : " + status;
}
string GLSLExtensionSet::generateGLSL() const {
    std::ostringstream output;
    for(auto pair : extensions_) {
        output << format(pair.first, pair.second) << std::endl;
    }
    return output.str();
}
void GLSLExtensionSet::generateGLSL(std::vector<std::string>& output) const {
    for(auto pair : extensions_) {
        output.push_back(format(pair.first, pair.second));
    }
}

//------------------------------------------------------------------------------

void LineTracker::addHeader(std::vector<std::string>& header) {
    header_ = header;
}

void LineTracker::addLineInfo(int lineNumber, const std::string& filename, int sourceLineNumber) {
    lineInfos_.emplace_back(lineNumber, filename, sourceLineNumber);
}

void LineTracker::reset() {
    header_.clear();
    lineInfos_.clear();
}

std::string LineTracker::resolveLineNumber(int number) const {
    // The error might be in the header:
    if(number <= static_cast<int>(header_.size())) {
        return "Header: " + header_[number-1];
    }

    // Ok, it isn't: Subtract the header size to get correct numbers from the LineInfo structure.
    number -= static_cast<int>(header_.size());
    int found = static_cast<int>(lineInfos_.size()) - 1;
    for (int i = static_cast<int>(lineInfos_.size()) - 1; i >= 0 ; i--) {
        if (lineInfos_[i].lineNumber_ <= number) {
            found = i;
            break;
        }
    }

    std::ostringstream result;
    //if (found >= 0) { //found is unsigned...
        result << FileSystem::fileName(lineInfos_[found].filename_) << ":"
               << (number - lineInfos_[found].lineNumber_ + lineInfos_[found].sourceLineNumber_);
    //}
    return result.str();
}


//------------------------------------------------------------------------------
//
const int GLSLVersionCollector::INVALID_VERSION = -1;

GLSLVersionCollector::GLSLVersionCollector()
    : version_(INVALID_VERSION)
    , additionalQualifier_(NONE)
{
}

void GLSLVersionCollector::addVersion(int version) {
    version_ = std::max(version, version_);
}

void GLSLVersionCollector::addAdditionalQualifier(AdditionalQualifier qualifier) {
    if(additionalQualifier_ == NONE || qualifier == additionalQualifier_) {
        additionalQualifier_ = qualifier;
    } else {
        LWARNINGC("GLSLVersionCollector", "Ignoring incompatible qualifier \"" << additionalQualifierAsString(qualifier) << "\" (incompatible to previously collected \"" << additionalQualifierAsString(additionalQualifier_) << "\")");
    }
}

bool GLSLVersionCollector::isInvalid() const {
    return version_ == INVALID_VERSION;
}

std::string GLSLVersionCollector::generateGLSL() const {
    std::string qualifier;
    return "#version " + std::to_string(version_) + " " + additionalQualifierAsString(additionalQualifier_) + "\n";
}

std::string GLSLVersionCollector::additionalQualifierAsString(AdditionalQualifier qualifier) {
    switch(qualifier) {
        case NONE:
            return "";
        case CORE:
            return "core";
        case COMPATIBILITY:
            return "compatibility";
        default:
            tgtAssert(false, "Unimplemented qualifier");
            return "";
    }
}

//------------------------------------------------------------------------------

const std::string ShaderPreprocessor::loggerCat_("tgt.Shader.ShaderPreprocessor");

ShaderPreprocessor::ShaderPreprocessor(ShaderObject* obj)
    : shd_(obj)
    , lineTracker_(obj->lineTracker_)
    , maxVersion_()
    , extensions_()
{
    parse();
}

static bool startsWith(const std::string& str, const std::string& start) {
    if(str.size() < start.size()) {
        return false;
    }
    for(string::size_type i=0; i < start.size(); ++i) {
        if(str[i] != start[i]) {
            return false;
        }
    }
    return true;
}

void ShaderPreprocessor::parsePart(const std::string& input, const std::string& name) {
    std::istringstream source(input);
    int locallinenumber = 0;

    lineTracker_.addLineInfo(activeLine_, name, 1);

    string line;
    bool inComment = false;
    while (std::getline(source, line)) {
        locallinenumber++;
        line = stripLine(line, inComment);

        if(startsWith(line, "#include ")) {
            string::size_type pos = line.find("\"", sizeof("#include"));
            string::size_type end = line.find("\"", pos + 1);
            string orig_filename(line, pos + 1, end - pos - 1);
            string filename = ShdrMgr.completePath(orig_filename);

            std::unique_ptr<File> file(FileSys.open(filename));
            string content;
            if ((!file) || (!file->isOpen())) {
                LERROR("Cannot open shader include '" << orig_filename << "' in " << lineTracker_.resolveLineNumber(activeLine_));
            }
            else {
                size_t len = file->size();
                // check if file is empty
                if (len == 0)
                    content = "";
                else
                    content = file->getAsString();
                file->close();

                outputComment("BEGIN INCLUDE " + filename, "BEGIN");

                if (!content.empty() && content[content.size() - 1] != '\n')
                    content += "\n";

                parsePart(content, filename);

                outputComment("END INCLUDE " + filename, "END");

                lineTracker_.addLineInfo(activeLine_, name, locallinenumber + 1);
            }
        } else if(startsWith(line, "#version ")) {
            std::istringstream versionStream(line.substr(sizeof("#version")));

            // Determine version
            int version;
            versionStream >> version;
            maxVersion_.addVersion(version);

            // Determine if core
            std::string rest;
            versionStream >> rest;

            if(rest.find("core") != string::npos) {
                maxVersion_.addAdditionalQualifier(GLSLVersionCollector::CORE);
            }

            if(rest.find("compatibility") != string::npos) {
                maxVersion_.addAdditionalQualifier(GLSLVersionCollector::COMPATIBILITY);
            }

            outputComment("PARSED VERSION: " + line);
        } else if(startsWith(line, "#extension ")) {
            std::istringstream extensionStream(line.substr(sizeof("#extension")));

            std::string extensionName;
            std::string delimiter;
            std::string extensionStatus;
            extensionStream >> extensionName;
            extensionStream >> delimiter;
            extensionStream >> extensionStatus;
            if(delimiter != ":") {
                LERROR("Could not find delimiter when parsing extension line: " << line);
                continue;
            }
            GLSLExtensionSet::ExtensionStatus status;
            if(extensionStatus == "enable") {
                status = GLSLExtensionSet::ExtensionStatus::ENABLED;
            } else if(extensionStatus == "require") {
                status = GLSLExtensionSet::ExtensionStatus::REQUIRED;
            } else if(extensionStatus == "warn") {
                status = GLSLExtensionSet::ExtensionStatus::WARN;
            } else if(extensionStatus == "disabled") {
                status = GLSLExtensionSet::ExtensionStatus::DISABLED;
            } else {
                LERROR("Could not determine extension status when parsing extension line: " << line);
                continue;
            }
            extensions_.set(extensionName, status);

            outputComment("PARSED EXTENSION: " + line);
        } else {
            result_ << line << "\n";
            activeLine_++;
        }
    }
}

void ShaderPreprocessor::parse() {
    activeLine_ = 1;
    lineTracker_.reset();
    result_.clear();

    if (!shd_->header_.empty()) {
        outputComment("BEGIN HEADER");
        parsePart(shd_->header_, "HEADER");
        outputComment("END HEADER");
    }

    parsePart(shd_->unparsedSource_, shd_->filename_);
    if(maxVersion_.isInvalid()) {
        LERROR("Could not find #version in " << shd_->filename_);
    }

    // Add version and extension header information to the line tracker
    std::vector<std::string> header;
    header.push_back(maxVersion_.generateGLSL());
    extensions_.generateGLSL(header);
    lineTracker_.addHeader(header);
}

string ShaderPreprocessor::stripLine(const string line, bool& inComment) {
    const char* pos = line.c_str();

    // True if we found text => Don't eat any more whitespaces!
    bool foundText = false;

    // The result accumulator
    string stripped = "";

    while(*pos != '\0') {
        if(inComment) {
            if(pos[0] == '*' && pos[1] == '/' /* the string is be null-terminated, so this is safe */) {
                // End of comment: Eat the '*/' and signal that the comment has ended
                pos += 2;
                inComment = false;
            } else {
                // Still in comment: Eat the char
                ++pos;
            }
        } else {
            if(pos[0] == '/' && pos[1] == '*' /* the string is be null-terminated, so this is safe */) {
                // Start of multiline comment: Eat '/*' and signal that a comment has started
                pos += 2;
                inComment = true;
            } else if(pos[0] == '/' && pos[1] == '/') {
                // Line comment. we are done!
                return stripped;
            } else {
                // Not in comment!
                if(foundText || (*pos != ' ' && *pos != '\t')) {
                    // Woohoo! We found a meaningful char! Write it to the result string
                    foundText = true;
                    stripped += *pos;
                    ++pos;
                } else {
                    // Remove leading whitespaces
                    ++pos;
                }
            }
        }
    }
    return stripped;
}

void ShaderPreprocessor::outputComment(const std::string& comment, const std::string& type) {
    result_ << "// " << comment << "\n";
    lineTracker_.addLineInfo(activeLine_, type, 0);
    activeLine_++;
}

std::string ShaderPreprocessor::getShaderDirective(const std::string& d) {
    string sourceStr(shd_->source_);
    string::size_type curPos = sourceStr.find(d + "(", 0);
    string::size_type length = d.length() + 1;
    if (curPos != string::npos) {
        string::size_type endPos = sourceStr.find(")", curPos);

        if (endPos != string::npos) {
            std::string ret = sourceStr.substr(curPos + length, endPos - curPos - length);
            // test for space, newline:
            if ((ret.find(" ", 0) == string::npos) && (ret.find("\n", 0) == string::npos) ) {
                LINFO("Directive " << d << ": " << ret);
                return ret;
            }
            else {
                LERROR("No spaces/newlines allowed inbetween directive brackets! Directive: " << d);
                return "";
            }
        }
        LERROR("Missing ending bracket for directive " << d);
        return "";
    }
    else {
        LWARNING("Could not locate directive " << d << "!");
        return "";
    }
}

std::string ShaderPreprocessor::getResult() const {
    return maxVersion_.generateGLSL() + extensions_.generateGLSL() + result_.str();
}

//------------------------------------------------------------------------------

const string ShaderObject::loggerCat_("tgt.Shader.ShaderObject");

static GLint getGLEnum(ShaderObject::ShaderType type) {
    switch(type) {
        case ShaderObject::VERTEX_SHADER:
            return GL_VERTEX_SHADER;
        case ShaderObject::FRAGMENT_SHADER:
            return GL_FRAGMENT_SHADER;
        case ShaderObject::GEOMETRY_SHADER:
            return GL_GEOMETRY_SHADER;
        case ShaderObject::COMPUTE_SHADER:
#ifdef GL_COMPUTE_SHADER
            return GL_COMPUTE_SHADER;
#else
            tgtAssert(false, "Compute shader not supported");
            return 0;
#endif
        case ShaderObject::UNDEFINED:
            tgtAssert(false, "Undefined shader type");
            return 0;
        default:
            tgtAssert(false, "Unimplemented shader type");
            return 0;
    }
}

ShaderObject::ShaderObject(const string& filename, ShaderType type)
    : filename_(filename)
    , shaderType_(type)
    , isCompiled_(false)
{
    id_ = glCreateShader(getGLEnum(shaderType_));
    if (id_ == 0)
        LERROR("ShaderObject(" + filename + ")::glCreateShader() returned 0");
}

ShaderObject::~ShaderObject() {
    glDeleteShader(id_);
}

void ShaderObject::loadSourceFromFile(const string& filename) {
    LDEBUG("Loading " << filename);
    File* file = FileSys.open(filename);

    // check if file is open
    if (!file || !file->isOpen()) {
        LDEBUG("File not found: " << filename);
        delete file;
        throw FileNotFoundException("", filename);
    }

    filename_ = filename;
    unparsedSource_ = file->getAsString();
    source_ = unparsedSource_;

    file->close();
    delete file;
}

void ShaderObject::uploadSource() {
    const GLchar* s = source_.c_str();
    glShaderSource(id_, 1,  &s, 0);
}

bool ShaderObject::compileShader() {
    isCompiled_ = false;

    ShaderPreprocessor p(this);
    source_ = p.getResult();

    uploadSource();

    glCompileShader(id_);
    GLint check = 0;
    glGetShaderiv(id_, GL_COMPILE_STATUS, &check);
    isCompiled_ = (check == GL_TRUE);
    return isCompiled_;
}

namespace {

int parseLogLineNumberNVIDIA(const std::string& message) {
    // Errors look like this:
    // 0(1397) : error C0000: syntax error, unexpected '=' at token "="
    std::istringstream ls(message);
    int id; // first number (probably used when multiple sources are bound), ignored
    if (ls >> id) {
        char c = 0; // should be an opening parenthesis
        if (ls >> c &&c == '(') {
            int num; // line number
            if (ls >> num)
                return num;
        }
    }
    return 0;
}

int parseLogLineNumberATI(const std::string& message) {
    // Errors look like this:
    // ERROR: 0:2785: 'frontPos' : undeclared identifier
    std::istringstream ls(message);
    std::string s;
    if (ls >> s && s == "ERROR:") {
        int id; // first number (probably used when multiple sources are bound), ignored
        if (ls >> id) {
            char c = 0; // should be a colon
            if (ls >> c && c == ':') {
                int num; // line number
                if (ls >> num)
                    return num;
            }
        }
    }

    return 0;
}

int parseLogLineNumberMesa(const std::string& message) {
    // Errors look like this:
    // 0:123(23): error: syntax error, unexpected NEW_IDENTIFIER
    std::istringstream ls(message);
    int id; // first number (probably used when multiple sources are bound), ignored
    if (ls >> id) {
        char c = 0; // should be a colon
        if (ls >> c &&c == ':') {
            int num; // line number
            if (ls >> num)
                return num;
        }
    }
    return 0;
}

/**
 * Tries to extract the line number for a given compile messages.
 * Returns 0 if no line number was found.
 */
int parseLogLineNumber(const std::string& message) {
    int n = 0;

    n = parseLogLineNumberNVIDIA(message);
    if (n > 0)
        return n;

    n = parseLogLineNumberATI(message);
    if (n > 0)
        return n;

    n = parseLogLineNumberMesa(message);
    if (n > 0)
        return n;

    return 0;
}

} // namespace

string ShaderObject::getCompilerLog() const {
    GLint len;
    glGetShaderiv(id_, GL_INFO_LOG_LENGTH , &len);

    if (len > 1) {
        GLchar* log = new GLchar[len];
        if (log == 0)
            return "Memory allocation for log failed!";
        GLsizei l;  // length returned
        glGetShaderInfoLog(id_, len, &l, log);
        std::istringstream str(log);
        delete[] log;

        std::ostringstream result;

        string line;
        while (getline(str, line)) {
            result << line;

            int num = parseLogLineNumber(line);
            if (num > 0)
                result << " [" << lineTracker_.resolveLineNumber(num) << "]";

            result << '\n';
        }

        return result.str();
    } else {
        return "";
    }
}

void ShaderObject::setHeader(const string& h) {
    header_ = h;
}

bool ShaderObject::rebuildFromFile() {
    try {
        loadSourceFromFile(filename_);
    }
    catch (const Exception& e) {
        LWARNING("Failed to load shader " << filename_ << ": " << e.what());
        return false;
    }

    uploadSource();

    if (!compileShader()) {
        LERROR("Failed to compile shader object " << filename_);
        LERROR("Compiler Log: \n" << getCompilerLog());
        return false;
    }

    return true;
}

//------------------------------------------------------------------------------

const string Shader::loggerCat_("tgt.Shader.Shader");

Shader::Shader()
    : isLinked_(false)
    , ignoreError_(false)
    , uniformsDirty_(true)
    , uniformsUnset_(0)
{
    id_ = glCreateProgram();
    if (id_ == 0)
        LERROR("Shader(): glCreateProgram() returned 0");
}

Shader::~Shader() {
    for (ShaderObjects::iterator iter = objects_.begin(); iter != objects_.end(); ++iter) {
        glDetachShader(id_, (*iter)->id_);
        delete (*iter);
    }
    glDeleteProgram(id_);
}

void Shader::attachObject(ShaderObject* obj) {
    glAttachShader(id_, obj->id_);
    objects_.push_back(obj);
    isLinked_ = false;
    uniformsDirty_ = true;
}

void Shader::detachObject(ShaderObject* obj) {
    glDetachShader(id_, obj->id_);
    objects_.remove(obj);
    isLinked_ = false;
    uniformsDirty_ = true;
}

void Shader::activate() {
    if (uniformsDirty_) {
        uniformNames_.clear();
        unsetUniforms_.clear();
    }

    if (isLinked_) {
        glUseProgram(id_);
        ShdrMgr.setActiveShader(this);

        if (uniformsDirty_) {
            GLsizei length = 0;
            GLint size = 0; // currently unused
            GLenum type = 0; // currently unused

            std::vector<GLchar> nameData(256); // Use constant size, since query of GL_ACTIVE_UNIFORM_MAX_LENGTH doesn't work with every driver.
            GLint numUniforms;
            glGetProgramiv(id_, GL_ACTIVE_UNIFORMS, &numUniforms);
            for (int i = 0; i < numUniforms; i++) {
                // Cache uniform location and name.
                glGetActiveUniform(id_, (GLuint)i, nameData.size(), &length, &size, &type, &nameData[0]);
                std::string name(static_cast<char*>(&nameData[0]), length);

                // Query location.
                GLint location = glGetUniformLocation(id_, name.c_str());
                if (location == -1)
                    continue;

                name = name.substr(0, name.find_first_of('[')); // Cut of array brackets.
                uniformNames_[location] = name;

                // Set uniform to unused, if not having gl_ prefix.
                if (name.find_first_of("gl_") != 0)
                    unsetUniforms_.insert(location);
            }

            uniformsDirty_ = false;
            uniformsUnset_ = unsetUniforms_.size();
        }
    } else {
        LERROR("::activate(): Shader is not linked");
    }
}

void Shader::deactivate() {

#ifdef TGT_DEBUG
    // Log unset uniforms.
    Shader* active = ShdrMgr.getActiveShader();
    if (active && active->uniformsUnset_ != active->unsetUniforms_.size()) {
        if (!active->unsetUniforms_.empty()) {

            std::string message = "\nUniforms not set:\n";
            for (GLint location : active->unsetUniforms_)
                message += " * " + active->uniformNames_.at(location) + "\n";

            message += "\nwithin shader originally built from the following files:\n";
            for (ShaderObject* object : active->getObjects()) {
                message += " * " + object->getOriginalFilename() + "\n";
            }

            // Consider to set a breakpoint here.
            LWARNING(message);
        }

        active->uniformsUnset_ = active->unsetUniforms_.size();
    }
#endif

    // Deactivate current program.
    glUseProgram(0);
    ShdrMgr.setActiveShader(0);
}

GLint Shader::getCurrentProgram() {
    GLint id;
    glGetIntegerv(GL_CURRENT_PROGRAM, &id);
    return id;
}

GLint Shader::getID() const {
    return id_;
}

bool Shader::hasObjects() const {
    return !objects_.empty();
}

const Shader::ShaderObjects& Shader::getObjects() const {
    return objects_;
}

bool Shader::isActivated() {
    GLint shader_nr;
    glGetIntegerv(GL_CURRENT_PROGRAM, &shader_nr);
    return (id_ == static_cast<GLuint>(shader_nr));
}

void Shader::detachObjectsByType(ShaderObject::ShaderType type) {
    for (ShaderObjects::iterator iter = objects_.begin(); iter != objects_.end(); ++iter) {
        if ((*iter)->getType() == type)
            detachObject(*iter);
        delete (*iter);
    }
    isLinked_ = false;
    uniformsDirty_ = true;
}

bool Shader::isLinked() const{
    return isLinked_;
}

bool Shader::linkProgram() {
    if (isLinked_) {
        // program is already linked: detach and re-attach everything
        for (ShaderObjects::iterator iter = objects_.begin(); iter != objects_.end(); ++iter) {
            glDetachShader(id_, (*iter)->id_);
            glAttachShader(id_, (*iter)->id_);
        }
    }

    isLinked_ = false;
    glLinkProgram(id_);
    GLint check = 0;
    glGetProgramiv(id_, GL_LINK_STATUS, &check);
    if (check)
        isLinked_ = true;

    uniformsDirty_ = true;

    return isLinked_;
}

string Shader::getLinkerLog() const {
    GLint len;
    glGetProgramiv(id_, GL_INFO_LOG_LENGTH , &len);

    if (len > 1) {
        GLchar* log = new GLchar[len];
        if (log == 0)
            return "Memory allocation for log failed!";
        GLsizei l;  // length returned
        glGetProgramInfoLog(id_, len, &l, log);
        string retStr(log);
        delete[] log;
        return retStr;
    }

    return "";
}

bool Shader::rebuild() {
    bool wasLinked = isLinked_;
    isLinked_ = false;
    uniformsDirty_ = true;

    if (wasLinked) {
        // program is already linked: detach and re-attach everything
        for (ShaderObjects::iterator iter = objects_.begin(); iter != objects_.end(); ++iter) {
            glDetachShader(id_, (*iter)->id_);
            (*iter)->uploadSource();
            if (!(*iter)->compileShader()) {
                LERROR("Failed to compile shader object.");
                LERROR("Compiler Log: \n" << (*iter)->getCompilerLog());
                return false;
            }

            glAttachShader(id_, (*iter)->id_);
        }
    }
    
    glLinkProgram(id_);
    GLint check = 0;
    glGetProgramiv(id_, GL_LINK_STATUS, &check);

    if (check) {
        isLinked_ = true;
        return true;
    } else {
        LERROR("Shader::rebuild(): Failed to link shader." );
        LERROR("Linker Log: \n" << getLinkerLog());
        return false;
    }
}

bool Shader::rebuildFromFile() {
    bool result = true;

    for (ShaderObjects::iterator iter = objects_.begin(); iter != objects_.end(); ++iter)
        result &= (*iter)->rebuildFromFile();

    result &= rebuild();

    return result;
}

void Shader::setHeaders(const string& customHeader) {
    for (ShaderObjects::iterator iter = objects_.begin(); iter != objects_.end(); ++iter) {
        (*iter)->setHeader(customHeader);
    }
}

void Shader::bindFragDataLocation(GLuint colorNumber, std::string name) {
    if (GpuCaps.getShaderVersion() >= GpuCapabilities::GlVersion::SHADER_VERSION_130) {
        glBindFragDataLocation(id_, colorNumber, name.c_str()); /* was ...EXT */
    }
}

void Shader::load(const string& filename, const string& customHeader) {
    return loadSeparate(filename + ".vert", "", filename + ".frag", customHeader);
}

void Shader::loadSeparate(const string& vert_filename, const string& geom_filename,
                          const string& frag_filename, const string& customHeader) {
    /*ShaderObject* frag = 0;
    ShaderObject* vert = 0;
    ShaderObject* geom = 0;

    if (!vert_filename.empty()) {
        vert = new ShaderObject(vert_filename, ShaderObject::VERTEX_SHADER);

        if (!customHeader.empty()) {
            vert->setHeader(customHeader);
        }

        try {
            vert->loadSourceFromFile(vert_filename);
        }
        catch (const Exception& e) {
            LDEBUG("Failed to load vertex shader " << vert_filename << ": " << e.what());
            delete vert;
            throw Exception("Failed to load vertex shader " + vert_filename + ": " + e.what());
        }

        vert->uploadSource();

        if (!vert->compileShader()) {
            LERROR("Failed to compile vertex shader " << vert_filename);
            LERROR("Compiler Log: \n" << vert->getCompilerLog());
            delete vert;
            throw Exception("Failed to compile vertex shader: " + vert_filename);
        }
    }

    if (!geom_filename.empty()) {
        geom = new ShaderObject(geom_filename, ShaderObject::GEOMETRY_SHADER);

        if (!customHeader.empty()) {
            geom->setHeader(customHeader);
        }

        try {
            geom->loadSourceFromFile(geom_filename);
        }
        catch (const Exception& e) {
            LDEBUG("Failed to load geometry shader " << geom_filename << ": " << e.what());
            delete vert;
            delete geom;
            throw Exception("Failed to load geometry shader " + geom_filename + ": " + e.what());
        }

        geom->uploadSource();
        if (!geom->compileShader()) {
            LERROR("Failed to compile geometry shader " << geom_filename);
            LERROR("Compiler Log: \n" << geom->getCompilerLog());
            delete vert;
            delete geom;
            throw Exception("Failed to compile geometry shader: " + geom_filename);
        }
    }

    if (!frag_filename.empty()) {
        frag = new ShaderObject(frag_filename, ShaderObject::FRAGMENT_SHADER);

        if (!customHeader.empty()) {
            frag->setHeader(customHeader);
        }

        try {
            frag->loadSourceFromFile(frag_filename);
        }
        catch (const Exception& e) {
            LDEBUG("Failed to load fragment shader " << frag_filename);
            delete frag;
            delete geom;
            delete vert;
            throw Exception("Failed to load fragment shader " + frag_filename + ": " + e.what());
        }

        if (GpuCaps.getShaderVersion() >= GpuCapabilities::GlVersion::SHADER_VERSION_130)
            bindFragDataLocation(0, "FragData0");

        frag->uploadSource();

        if (!frag->compileShader()) {
            LERROR("Failed to compile fragment shader " << frag_filename);
            LERROR("Compiler Log: \n" << frag->getCompilerLog());
            delete vert;
            delete geom;
            delete frag;
            throw Exception("Failed to compile fragment shader: " + frag_filename);
        }
    }

    // Attach ShaderObjects, dtor will take care of freeing them
    if (frag)
        attachObject(frag);
    if (vert)
        attachObject(vert);
    if (geom)
        attachObject(geom);

    if (!linkProgram()) {
        LERROR("Failed to link shader (" << vert_filename << ","  << frag_filename << "," << geom_filename << ")");
        if (vert) {
            LERROR(vert->filename_ << " Vertex shader compiler log: \n" << vert->getCompilerLog());
            detachObject(vert);
            delete vert;
        }
        if (geom) {
            LERROR(geom->filename_ << " Geometry shader compiler log: \n" << geom->getCompilerLog());
            detachObject(geom);
            delete geom;
        }
        if (frag) {
            LERROR(frag->filename_ << " Fragment shader compiler log: \n" << frag->getCompilerLog());
            detachObject(frag);
            delete frag;
        }

        LERROR("Linker Log: \n" << getLinkerLog());
        throw Exception("Failed to link shader (" + vert_filename + "," + frag_filename + "," + geom_filename + ")");
    }


    if (vert && vert->getCompilerLog().size() > 1) {
        LDEBUG("Vertex shader compiler log for file '" << vert_filename
               << "': \n" << vert->getCompilerLog());
    }
    if (geom && geom->getCompilerLog().size() > 1) {
        LDEBUG("Geometry shader compiler log for file '" << geom_filename
               << "': \n" << geom->getCompilerLog());
    }
    if (frag && frag->getCompilerLog().size() > 1) {
        LDEBUG("Fragment shader compiler log for file '" << frag_filename
               << "': \n" << frag->getCompilerLog());
    }

    if (getLinkerLog().size() > 1) {
        LDEBUG("Linker log for '" << vert_filename << "' and '"
               << frag_filename << "' and '"
               << geom_filename << "': \n" << getLinkerLog());
    }*/


    if (!geom_filename.empty()) {
        std::vector<std::string> filenames(3);
        filenames[0] = vert_filename;
        filenames[1] = geom_filename;
        filenames[2] = frag_filename;

        std::vector<ShaderObject::ShaderType> types(3);
        types[0] = ShaderObject::VERTEX_SHADER;
        types[1] = ShaderObject::GEOMETRY_SHADER;
        types[2] = ShaderObject::FRAGMENT_SHADER;

        return loadSeparate(filenames, types, customHeader);
    }
    else {
        std::vector<std::string> filenames(2);
        filenames[0] = vert_filename;
        filenames[1] = frag_filename;

        std::vector<ShaderObject::ShaderType> types(2);
        types[0] = ShaderObject::VERTEX_SHADER;
        types[1] = ShaderObject::FRAGMENT_SHADER;

        return loadSeparate(filenames, types, customHeader);
    }
}


void Shader::loadSeparate(const std::vector<string>& filenames, const std::vector<ShaderObject::ShaderType>& types,
                      const std::string& customHeader, bool setFragDataLocation) {

    if ((filenames.size() != types.size()) || filenames.empty()) {
        LDEBUG("Failed to load shaders, array size of filenames and types differs or is empty");
        throw Exception("Failes to load shaders, array size of filenames and types differs or is empty");
    }

    size_t numShaders = filenames.size();

    ShaderObject** shaderObjects = new ShaderObject*[numShaders];

    for (size_t i = 0; i < numShaders; ++i) {
        shaderObjects[i] = new ShaderObject(filenames[i], types[i]);

        if (!customHeader.empty()) {
            shaderObjects[i]->setHeader(customHeader);
        }

        try {
            shaderObjects[i]->loadSourceFromFile(filenames[i]);
        }
        catch (const Exception& e) {
            LDEBUG("Failed to load shader " << filenames[i] << ": " << e.what());
            for (size_t j = 0; j <= i; ++j)
                delete shaderObjects[j];
            delete[] shaderObjects;
            throw Exception("Failed to load shader " + filenames[i] + ": " + e.what());
        }

        if (types[i] == ShaderObject::FRAGMENT_SHADER && setFragDataLocation)
            if (GpuCaps.getShaderVersion() >= GpuCapabilities::GlVersion::SHADER_VERSION_130)
                bindFragDataLocation(0, "FragData0");

        shaderObjects[i]->uploadSource();

        if (!shaderObjects[i]->compileShader()) {
            LERROR("Failed to compile shader " << filenames[i]);
            LERROR("Compiler Log: \n" << shaderObjects[i]->getCompilerLog());
            for (size_t j = 0; j <= i; ++j)
                delete shaderObjects[j];
            delete[] shaderObjects;
            throw Exception("Failed to compile shader: " + filenames[i]);
        }
    }

    // Attach ShaderObjects, dtor will take care of freeing them
    for (size_t i = 0; i < numShaders; ++i) {
        if (shaderObjects[i])
            attachObject(shaderObjects[i]);
    }

    std::stringstream files;
    for (size_t i = 0; i < numShaders; ++i) {
        if (i != numShaders - 1)
            files << filenames[i] << ",";
        else
            files << filenames[i];
    }

    if (!linkProgram()) {
        LERROR("Failed to link shader (" << files.str() << ")");

        for (size_t i = 0; i < numShaders; ++i) {
            if (shaderObjects[i]) {
                LERROR(shaderObjects[i] << filenames[i] << " compiler log: \n" << shaderObjects[i]->getCompilerLog());
                detachObject(shaderObjects[i]);
                delete shaderObjects[i];
            }
        }

        delete[] shaderObjects;
        LERROR("Linker Log: \n" << getLinkerLog());
        throw Exception("Failed to link shader (" + files.str() + ")");
    }

    for (size_t i = 0; i < numShaders; ++i) {
        if (shaderObjects[i] && shaderObjects[i]->getCompilerLog().size() > 1) {
            LDEBUG("Shader compiler log for file '" << filenames[i]
               << "': \n" << shaderObjects[i]->getCompilerLog());
        }
    }

    if (getLinkerLog().size() > 1) {
        LDEBUG("Linker log for '" << files.str() << ": \n" << getLinkerLog());
    }

    delete[] shaderObjects;
}

std::vector<std::string> Shader::getUniformNames() const {
    std::vector<std::string> names;
    names.reserve(uniformNames_.size());
    for (const auto& iter : uniformNames_) {
        names.push_back(iter.second);
    }
    return names;
}

std::vector<std::string> Shader::getUnsetUniformNames() const {
    std::vector<std::string> names;
#ifdef TGT_DEBUG
    names.reserve(unsetUniforms_.size());
    for (GLint location : unsetUniforms_) {
        names.push_back(uniformNames_.at(location));
    }
#endif
    return names;
}

void Shader::setIgnoreUnsetUniform(const std::string& name) {
    GLint location = getUniformLocation(name);
    if (location != -1)
        unsetUniforms_.erase(location);
}

GLint Shader::useUniform(const std::string& name) {
#ifdef TGT_DEBUG
    if(!isLinked_) {
        return -1;
    }

    tgtAssert(isActivated(), "Shader not currently active");

    GLint l = getUniformLocation(name);
    if (l != -1) {
        unsetUniforms_.erase(l);
    }

    return l;
#else
    return getUniformLocation(name);
#endif
}

GLint Shader::getUniformLocation(const string& name) {
    GLint location = glGetUniformLocation(id_, name.c_str());
    if(location == -1 && !ignoreError_)
        LWARNING("Failed to locate uniform Location: " << name);
    return location;
}

void Shader::setIgnoreUniformLocationError(bool ignoreError) {
    ignoreError_ = ignoreError;
}

bool Shader::getIgnoreUniformLocationError() const {
    return ignoreError_;
}

// Floats
bool Shader::setUniform(const string& name, GLfloat value) {
    GLint l = useUniform(name);
    if (l == -1)
        return false;
    glUniform1f(l, value);
    return true;
}

bool Shader::setUniform(const string& name, GLfloat v1, GLfloat v2) {
    GLint l = useUniform(name);
    if (l == -1)
        return false;
    glUniform2f(l, v1, v2);
    return true;
}

bool Shader::setUniform(const string& name, GLfloat v1, GLfloat v2, GLfloat v3) {
    GLint l = useUniform(name);
    if (l == -1)
        return false;
    glUniform3f(l, v1, v2, v3);
    return true;
}

bool Shader::setUniform(const string& name, GLfloat v1, GLfloat v2, GLfloat v3, GLfloat v4) {
    GLint l = useUniform(name);
    if (l == -1)
        return false;
    glUniform4f(l, v1, v2, v3, v4);
    return true;
}

bool Shader::setUniform(const string& name, GLfloat* v, int count) {
    GLint l = useUniform(name);
    if (l == -1)
        return false;
    glUniform1fv(l, count, v);
    return true;
}

// Integers
bool Shader::setUniform(const string& name, GLint value) {
    GLint l = useUniform(name);
    if (l == -1)
        return false;
    glUniform1i(l, value);
    return true;
}

bool Shader::setUniform(const string& name, GLint v1, GLint v2) {
    GLint l = useUniform(name);
    if (l == -1)
        return false;
    glUniform2i(l, v1, v2);
    return true;
}

bool Shader::setUniform(const string& name, GLint v1, GLint v2, GLint v3) {
    GLint l = useUniform(name);
    if (l == -1)
        return false;
    glUniform3i(l, v1, v2, v3);
    return true;
}

bool Shader::setUniform(const string& name, GLint v1, GLint v2, GLint v3, GLint v4) {
    GLint l = useUniform(name);
    if (l == -1)
        return false;
    glUniform4i(l, v1, v2, v3, v4);
    return true;
}

bool Shader::setUniform(const string& name, GLint* v, int count) {
    GLint l = useUniform(name);
    if (l == -1)
        return false;
    glUniform1iv(l, count, v);
    return true;
}

bool Shader::setUniform(const string& name, bool value) {
    GLint l = useUniform(name);
    if (l == -1)
        return false;
    glUniform1i(l, static_cast<GLint>(value));
    return true;
}

bool Shader::setUniform(const string& name, bool v1, bool v2) {
    GLint l = useUniform(name);
    if (l == -1)
        return false;
    glUniform2i(l, static_cast<GLint>(v1), static_cast<GLint>(v2));
    return true;
}

bool Shader::setUniform(const string& name, bool v1, bool v2, bool v3) {
    GLint l = useUniform(name);
    if (l == -1)
        return false;
    glUniform3i(l, static_cast<GLint>(v1), static_cast<GLint>(v2), static_cast<GLint>(v3));
    return true;
}

bool Shader::setUniform(const string& name, bool v1, bool v2, bool v3, bool v4) {
    GLint l = getUniformLocation(name);
    if (l == -1)
        return false;
    glUniform4i(l, static_cast<GLint>(v1), static_cast<GLint>(v2), static_cast<GLint>(v3), static_cast<GLint>(v4));
    return true;
}

bool Shader::setUniform(const string& name, GLboolean* v, int count) {
    GLint l = useUniform(name);
    if (l == -1)
        return false;
    GLint* vector = new GLint[count];
    for (int i=0; i < count; i++)
        vector[i] = static_cast<GLint>( v[i] );
    glUniform1iv(l, count, vector);
    delete[] vector;
    return true;
}

// Vectors
bool Shader::setUniform(const string& name, const Vector2f& value) {
    GLint l = useUniform(name);
    if (l == -1)
        return false;
    glUniform2f(l, value.x, value.y);
    return true;
}

bool Shader::setUniform(const string& name, Vector2f* vectors, GLsizei count) {
    GLint l = useUniform(name);
    if (l == -1)
        return false;
    //TODO: use the adress directly, without copying, same below. joerg
    GLfloat* values = new GLfloat[2*count];
    for (int i=0; i < count; i++){
        values[2*i] = vectors[i].x;
        values[2*i+1] = vectors[i].y;
    }
    glUniform2fv(l, count, values);
    delete[] values;
    return true;
}

bool Shader::setUniform(const string& name, const Vector3f& value) {
    GLint l = useUniform(name);
    if (l == -1)
        return false;
    glUniform3f(l, value.x, value.y, value.z);
    return true;
}

bool Shader::setUniform(const string& name, Vector3f* vectors, GLsizei count) {
    GLint l = useUniform(name);
    if (l == -1)
        return false;
    GLfloat* values = new GLfloat[3*count];
    for (int i=0; i < count; i++) {
        values[3*i] = vectors[i].x;
        values[3*i+1] = vectors[i].y;
        values[3*i+2] = vectors[i].z;
    }
    glUniform3fv(l, count, values);
    delete[] values;
    return true;
}

bool Shader::setUniform(const string& name, const Vector4f& value) {
    GLint l = useUniform(name);
    if (l == -1)
        return false;
    glUniform4f(l, value.x, value.y, value.z, value.w);
    return true;
}

bool Shader::setUniform(const string& name, Vector4f* vectors, GLsizei count) {
    GLint l = useUniform(name);
    if (l == -1)
        return false;
    GLfloat* values = new GLfloat[4*count];
    for (int i=0; i < count; i++) {
        values[4*i] = vectors[i].x;
        values[4*i+1] = vectors[i].y;
        values[4*i+2] = vectors[i].z;
        values[4*i+3] = vectors[i].a;
    }
    glUniform4fv(l, count, values);
    delete[] values;
    return true;
}

bool Shader::setUniform(const string& name, const ivec2& value) {
    GLint l = useUniform(name);
    if (l == -1)
        return false;
    glUniform2i(l, value.x, value.y);
    return true;
}

bool Shader::setUniform(const string& name, ivec2* vectors, GLsizei count) {
    GLint l = useUniform(name);
    if (l == -1)
        return false;
    GLint* values = new GLint[2*count];
    for (int i=0; i < count; i++) {
        values[2*i] = vectors[i].x;
        values[2*i+1] = vectors[i].y;
    }
    glUniform2iv(l, count, values);
    delete[] values;
    return true;
}

bool Shader::setUniform(const string& name, const ivec3& value) {
    GLint l = useUniform(name);
    if (l == -1)
        return false;
    glUniform3i(l, value.x, value.y, value.z);
    return true;
}

bool Shader::setUniform(const string& name, ivec3* vectors, GLsizei count) {
    GLint l = useUniform(name);
    if (l == -1)
        return false;
    GLint* values = new GLint[3*count];
    for (int i=0; i < count; i++) {
        values[3*i] = vectors[i].x;
        values[3*i+1] = vectors[i].y;
        values[3*i+2] = vectors[i].z;
    }
    glUniform3iv(l, count, values);
    delete[] values;
    return true;
}

bool Shader::setUniform(const string& name, const ivec4& value) {
    GLint l = useUniform(name);
    if (l == -1)
        return false;
    glUniform4i(l, value.x, value.y, value.z, value.w);
    return true;
}

bool Shader::setUniform(const string& name, ivec4* vectors, GLsizei count) {
    GLint l = useUniform(name);
    if (l == -1)
        return false;
    GLint* values = new GLint[4*count];
    for (int i=0; i < count; i++) {
        values[4*i] = vectors[i].x;
        values[4*i+1] = vectors[i].y;
        values[4*i+2] = vectors[i].z;
        values[4*i+3] = vectors[i].a;
    }
    glUniform4iv(l, count, values);
    delete[] values;
    return true;
}

// Note: Matrix is transposed by OpenGL
bool Shader::setUniform(const string& name, const Matrix2f& value, bool transpose) {
    GLint l = useUniform(name);
    if (l == -1)
        return false;
    glUniformMatrix2fv(l, 1, !transpose, value.elem);
    return true;
}

bool Shader::setUniform(const string& name, const Matrix3f& value, bool transpose) {
    GLint l = useUniform(name);
    if (l == -1)
        return false;
    glUniformMatrix3fv(l, 1, !transpose, value.elem);
    return true;
}

bool Shader::setUniform(const string& name, const Matrix4f& value, bool transpose) {
    GLint l = useUniform(name);
    if (l == -1)
        return false;
    glUniformMatrix4fv(l, 1, !transpose, value.elem);
    return true;
}

// Attribute locations
void Shader::setAttributeLocation(GLuint index, const std::string& name) {
    glBindAttribLocation(id_, index, name.c_str());
}

GLint Shader::getAttributeLocation(const string& name) {
    GLint l;
    l = glGetAttribLocation(id_, name.c_str());
    if (l == -1)
        LERROR("Failed to locate attribute Location: " << name);
    return l;
}

//------------------------------------------------------------------------------

const string ShaderManager::loggerCat_("tgt.Shader.Manager");

ShaderManager::ShaderManager()
  : ResourceManager<Shader>(false)
  , activeShader_(0)
{}

ShaderManager::~ShaderManager() {
    for(Shader* shader : staticShaders_) {
        dispose(shader);
    }
}

Shader* ShaderManager::load(const string& filename, const string& customHeader,
                            bool activate) {
    return loadSeparate(filename + ".vert", filename + ".frag", customHeader, activate);
}

Shader* ShaderManager::loadSeparate(const string& vert_filename, const string& frag_filename,
                                    const string& customHeader, bool activate) {
    return loadSeparate(vert_filename, "", frag_filename, customHeader, activate);
}

Shader* ShaderManager::loadSeparate(const string& vert_filename, const string& geom_filename,
                                    const string& frag_filename,
                                    const string& customHeader, bool activate) {
    /*LDEBUG("Loading files " << vert_filename << " and " << frag_filename);
    if (!GpuCaps.areShadersSupported()) {
        LERROR("Shaders are not supported.");
        throw Exception("Shaders are not supported.");
    }

    // create a somewhat unique identifier for this shader triple
    string identifier = vert_filename + "#" +  frag_filename + "#" + geom_filename;

    if (isLoaded(identifier)) {
        LDEBUG("Shader already loaded. Increase usage count.");
        increaseUsage(identifier);
        return get(identifier);
    }

    Shader* shdr = new Shader();

    // searching in all paths for every shader
    string vert_completeFilename;
    if (!vert_filename.empty())
        vert_completeFilename = completePath(vert_filename);

    string geom_completeFilename;
    if (!geom_filename.empty())
        geom_completeFilename = completePath(geom_filename);

    string frag_completeFilename;
    if (!frag_filename.empty())
        frag_completeFilename = completePath(frag_filename);

    // loading and linking found shaders
    try {
        shdr->loadSeparate(vert_completeFilename, geom_completeFilename,
                           frag_completeFilename, customHeader);
        // register even when caching is disabled, needed for rebuildFromFile()
        reg(shdr, identifier);

        if (activate)
            shdr->activate();

        return shdr;
    }
    catch (const Exception& e) {
        delete shdr;
        throw;
    }*/

    if (!geom_filename.empty()) {
        std::vector<std::string> filenames(3);
        filenames[0] = vert_filename;
        filenames[1] = geom_filename;
        filenames[2] = frag_filename;

        std::vector<ShaderObject::ShaderType> types(3);
        types[0] = ShaderObject::VERTEX_SHADER;
        types[1] = ShaderObject::GEOMETRY_SHADER;
        types[2] = ShaderObject::FRAGMENT_SHADER;

        return loadSeparate(filenames, types, customHeader, activate);
    }
    else {
        std::vector<std::string> filenames(2);
        filenames[0] = vert_filename;
        filenames[1] = frag_filename;

        std::vector<ShaderObject::ShaderType> types(2);
        types[0] = ShaderObject::VERTEX_SHADER;
        types[1] = ShaderObject::FRAGMENT_SHADER;

        return loadSeparate(filenames, types, customHeader, activate);
    }
}

Shader* ShaderManager::loadSeparate(const std::vector<std::string>& filenames, const std::vector<ShaderObject::ShaderType>& types,
                                    const std::string& customHeader, bool activate, bool setFragDataLocation) {

    if ((filenames.size() != types.size()) || filenames.empty()) {
        LDEBUG("Failed to load shaders, array size of filenames and types differs or is empty");
        throw Exception("Failes to load shaders, array size of filenames and types differs or is empty");
    }

    size_t numShaders = filenames.size();

#ifdef TGT_DEBUG
    std::stringstream files;
    for (size_t i = 0; i < numShaders; ++i) {
        if (i != numShaders -1)
            files << filenames[i] << " and ";
        else
            files << filenames[i];
    }

    LDEBUG("Loading files " + files.str());
#endif

    if (!GpuCaps.areShadersSupported()) {
        LERROR("Shaders are not supported.");
        throw Exception("Shaders are not supported.");
    }

    // create a somewhat unique identifier for this shader combination
    std::stringstream identifierStream;
    for (size_t i = 0; i < numShaders; ++i) {
        if (i != numShaders -1)
            identifierStream << filenames[i] << "#";
        else
            identifierStream << filenames[i];
    }
    string identifier = identifierStream.str();

    if (isLoaded(identifier)) {
        LDEBUG("Shader already loaded. Increase usage count.");
        increaseUsage(identifier);
        return get(identifier);
    }

    Shader* shdr = new Shader();

    // searching in all paths for every shader
    std::vector<string> completeFilenames(numShaders);
    for (size_t i = 0; i < numShaders; ++i) {
        if (!filenames[i].empty())
            completeFilenames[i] = completePath(filenames[i]);
        else
            completeFilenames[i] = "";
    }

    try {
        shdr->loadSeparate(completeFilenames, types, customHeader, setFragDataLocation);
        // register even when caching is disabled, needed for rebuildFromFile()
        reg(shdr, identifier);

        if (activate)
            shdr->activate();

        return shdr;
    }
    catch (const Exception& /*e*/) {
        delete shdr;
        throw;
    }
}

void ShaderManager::registerStaticShader(Shader* shader) {
    staticShaders_.push_back(shader);
}

bool ShaderManager::rebuildAllShadersFromFile() {
    bool result = true;

    for (std::map<Shader*, ResourceManager<Shader>::Resource*>::iterator iter = resourcesByPtr_.begin();
         iter != resourcesByPtr_.end(); ++iter)
    {
        result &= iter->first->rebuildFromFile();
    }

    return result;
}

} // namespace
