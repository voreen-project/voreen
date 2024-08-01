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

#ifndef TGT_SHADERMANAGER_H
#define TGT_SHADERMANAGER_H

#include <list>
#include <string>
#include <unordered_set>
#include <unordered_map>

#include "tgt/exception.h"
#include "tgt/manager.h"
#include "tgt/matrix.h"
#include "tgt/tgt_gl.h"
#include "tgt/types.h"
#include "tgt/vector.h"

namespace tgt {

class LineTracker {
public:
    // Helper for resolving line number when includes are used in shader files
    struct LineInfo {
        LineInfo(int n, std::string s, int sn)
            : lineNumber_(n), filename_(s), sourceLineNumber_(sn) {}

        int lineNumber_;          //< line number in preprocessed file
        std::string filename_;    //< filename of included file
        int sourceLineNumber_;    //< line number in included file (needed when it itself
                                  //< includes another file)
    };
    /**
     * Add header information.
     */
    void addHeader(std::vector<std::string>& header);

    /**
     * Add header information.
     */
    void addLineInfo(int lineNumber, const std::string& filename, int sourceLineNumber);

    /**
     * Clear line info and header data.
     */
    void reset();

    /**
     * Resolve the line number to take into account the include directives.
     * Returns a string containing file name and line number in that file.
     */
    std::string resolveLineNumber(int number) const;

private:
    std::vector<LineInfo> lineInfos_; ///< keeps track of line numbers when includes are used
    std::vector<std::string> header_;
};

/**
 * Type of a shader object, can be vertex, fragment or geometry shader
 *
 * #include statements are allowed.
 *
 */
class TGT_API ShaderObject {
public:
    friend class Shader;
    friend class ShaderPreprocessor;

    enum ShaderType {
          VERTEX_SHADER = 0
        , FRAGMENT_SHADER = 1
        , GEOMETRY_SHADER = 2
        , COMPUTE_SHADER = 3
        , UNDEFINED = 255
    };


    /**
     * Creates a shader object of the specified type
     */
    ShaderObject(const std::string& filename, ShaderType type = VERTEX_SHADER);

    /**
     * Deletes the shader and source
     */
    ~ShaderObject();

    /**
     * Loads the shader source from the specified file.
     *
     * @throw Exception if loading failed.
     */
    void loadSourceFromFile(const std::string& filename);

    bool compileShader();

    bool isCompiled() const { return isCompiled_; }

    std::string getCompilerLog() const;

    bool rebuildFromFile();

    /**
     * Use h as header for shadersource (copies h)
     */
    void setHeader(const std::string& h);

    ShaderType getType() const { return shaderType_; }
    const std::string& getOriginalFilename() const { return filename_; }

    void setSource(std::string source) {
        source_ = source;
        unparsedSource_ = source;
    }
    const std::string getSource() { return unparsedSource_; }


protected:
    void uploadSource();

    std::string filename_;
    ShaderType shaderType_;

    GLuint id_;
    std::string source_;
    std::string unparsedSource_;
    std::string header_;
    bool isCompiled_;

    LineTracker lineTracker_;
    static const std::string loggerCat_;
};

//------------------------------------------------------------------------------

/**
 * Represents an OpenGL shader program, consisting of linked ShaderObjects.
 *
 * @note Convenient loading of shaders from file is provided by ShaderManager.
 */
class TGT_API Shader {
    friend class ShaderManager;
public:
    typedef std::list<ShaderObject*> ShaderObjects;

    Shader();

    /**
     * Detaches all shader objects, deletes them an disposes all textures
     */
    ~Shader();

    /**
     * Attach shader object to Shader
     */
    void attachObject(ShaderObject* obj);
    void detachObject(ShaderObject* obj);
    void detachObjectsByType(ShaderObject::ShaderType type);

    /**
     * Link all shader objects to one shader.
     * Will re-link already linked shaders.
     * @return true for success
     */
    bool linkProgram();

    bool rebuild();
    bool rebuildFromFile();

    void setHeaders(const std::string& customHeader);

    void bindFragDataLocation(GLuint colorNumber, std::string name);

    GLint getID() const;

    /**
     * Activates the shader
     */
    void activate();

    static void deactivate();

    static GLint getCurrentProgram();

    /**
     * Returns whether the Shader has at least one attached shader object.
     */
    bool hasObjects() const;

    /**
     * Returns a collection of objects currently attached to the shader
     */
    const ShaderObjects& getObjects() const;

    /**
     * Returns whether the Shader is currently activated
     */
    bool isActivated();

    bool isLinked() const;

    std::string getLinkerLog() const;

    //
    // Uniform stuff
    //

    /**
     * Returns all uniform names within the shader.
     * Returns an empty vector if program is not linked.
     */
    std::vector<std::string> getUniformNames() const;

    /**
     * Returns all names of the uniforms currently not being set.
     * Returns an empty vector if program is not linked or if not compiled in debug mode.
     */
    std::vector<std::string> getUnsetUniformNames() const;

    /**
     * Returns uniform location, or -1 on failure
     * Logs warning if location could not be found and 
     * setIgnoreUniformLocationError was not called before.
     */
    GLint getUniformLocation(const std::string& name);

    void setIgnoreUniformLocationError(bool ignoreError);
    bool getIgnoreUniformLocationError() const ;

    /**
     * In case the specified uniform has not been set via the Shader wrapper
     * it will be ignored in terms of unset uniforms.
     * This function can be used as workaround for explicit binding.
     * Does nothing, if the given uniform is not defined by any attached shader object.
     */
    void setIgnoreUnsetUniform(const std::string& name);
    
    // Floats
    bool setUniform(const std::string& name, GLfloat value);
    bool setUniform(const std::string& name, GLfloat v1, GLfloat v2);
    bool setUniform(const std::string& name, GLfloat v1, GLfloat v2, GLfloat v3);
    bool setUniform(const std::string& name, GLfloat v1, GLfloat v2, GLfloat v3, GLfloat v4);
    bool setUniform(const std::string& name, GLfloat* v, int count);

    // Integers
    bool setUniform(const std::string& name, GLint value);
    bool setUniform(const std::string& name, GLint v1, GLint v2);
    bool setUniform(const std::string& name, GLint v1, GLint v2, GLint v3);
    bool setUniform(const std::string& name, GLint v1, GLint v2, GLint v3, GLint v4);
    bool setUniform(const std::string& name, GLint* v, int count);

    // Booleans
    bool setUniform(const std::string& name, bool value);
    bool setUniform(const std::string& name, bool v1, bool v2);
    bool setUniform(const std::string& name, bool v1, bool v2, bool v3);
    bool setUniform(const std::string& name, bool v1, bool v2, bool v3, bool v4);
    bool setUniform(const std::string& name, GLboolean* v, int count);

    // Vectors
    bool setUniform(const std::string& name, const Vector2f& value);
    bool setUniform(const std::string& name, Vector2f* vectors, GLsizei count = 1);
    bool setUniform(const std::string& name, const Vector3f& value);
    bool setUniform(const std::string& name, Vector3f* vectors, GLsizei count = 1);
    bool setUniform(const std::string& name, const Vector4f& value);
    bool setUniform(const std::string& name, Vector4f* vectors, GLsizei count = 1);
    bool setUniform(const std::string& name, const ivec2& value);
    bool setUniform(const std::string& name, ivec2* vectors, GLsizei count = 1);
    bool setUniform(const std::string& name, const ivec3& value);
    bool setUniform(const std::string& name, ivec3* vectors, GLsizei count = 1);
    bool setUniform(const std::string& name, const ivec4& value);
    bool setUniform(const std::string& name, ivec4* vectors, GLsizei count = 1);

    // Note: Matrix is transposed by OpenGL
    bool setUniform(const std::string& name, const Matrix2f& value, bool transpose = false);
    bool setUniform(const std::string& name, const Matrix3f& value, bool transpose = false);
    bool setUniform(const std::string& name, const Matrix4f& value, bool transpose = false);
    
    // Attribute locations
    void setAttributeLocation(GLuint index, const std::string& name);
    GLint getAttributeLocation(const std::string& name);

protected:
    /**
     * Load filename.vert and filename.frag (vertex and fragment shader) and link shader.
     *
     * @param customHeader Header to be put in front of the shader source.
     *
     * @throw Exception if loading failed
     */
    void load(const std::string& filename, const std::string& customHeader = "");

    /**
     * Load vertex shader \p vertFilename, geometry shader \p geomFilename,
     * fragment shader \p fragFilename.
     *
     * @param customHeader header to be put in front of the shader source
     *
     * @throw Exception if loading failed
     */
    void loadSeparate(const std::string& vertFilename, const std::string& geomFilename,
        const std::string& fragFilename, const std::string& customHeader = "");

    /**
     * Load shaders.
     *
     * @param filenames shader code files
     * @param types the types of the shaders in the filenames
     * @param customHeader header to put in front of the shader sources
     * @param setFragDataLocation set frag data location? (do not set for GLSL version >= 4)
     * @param Exception if loading failes
     */
    void loadSeparate(const std::vector<std::string>& filenames, const std::vector<ShaderObject::ShaderType>& types,
                      const std::string& customHeader = "", bool setFragDataLocation = true);

    /**
     * Marks a uniform as properly initialized within the specified shader.
     * @param name uniform name to be used
     * @return the location of the uniform name or -1 if not defined.
     */
    GLint useUniform(const std::string& name);

    /// Members

    ShaderObjects objects_;

    GLuint id_;
    bool isLinked_;
    bool ignoreError_;

    std::unordered_map<GLint, std::string> uniformNames_;     ///< Maps locations to names
    std::unordered_set<GLint> unsetUniforms_;                 ///< Holds all unused/unset uniform locations
    bool uniformsDirty_;                                      ///< Determines if uniform locations need to be updated
    size_t uniformsUnset_;                                    ///< Number of uniforms being set right after the last call of 'deactivate'

    static const std::string loggerCat_;
};

//------------------------------------------------------------------------------

class ShaderManager;
#ifdef DLL_TEMPLATE_INST
template class TGT_API Singleton<ShaderManager>;
#endif
#ifdef DLL_TEMPLATE_INST
template class TGT_API ResourceManager<Shader>;
#endif

/**
 * Loads shaders from the file system, managing a shader search path.
 *
 * @see ResourceManager
 */
class TGT_API ShaderManager : public ResourceManager<Shader>, public Singleton<ShaderManager> {
    friend class Shader;
public:

    ShaderManager();

    virtual ~ShaderManager();

    /**
     * Load filename.vert and filename.frag (vertex and fragment shader), link shader and
     * activate it by default.
     *
     * @param customHeader Header to be put in front of the shader source
     * @param activate activate the shader after loading
     *
     * @return The loaded shader
     *
     * @throw Exception if loading failed
     */
    Shader* load(const std::string& filename, const std::string& customHeader = "",
                 bool activate = true);

    /**
     * Load vertex shader \p vertFilename and fragment shader \p fragFilename,
     * link shader and activate it by default.
     *
     * You have to pass the complete filenames, inclusive file extensions (".vert", ".frag").
     *
     * @param customHeader header to be put in front of the shader source
     * @param activate activate the shader after loading
     *
     * @return The loaded shader
     *
     * @throw Exception if loading failed
     */
    Shader* loadSeparate(const std::string& vertFilename, const std::string& fragFilename,
                         const std::string& customHeader = "", bool activate = true);

    /**
     * Load vertex shader \p vertFilename, geometry shader \p geomFilename,
     * fragment shader \p fragFilename, link shader and activate it by default.
     *
     * You have to pass the complete filenames, inclusive file extensions (".vert", ".geom", frag").
     *
     * @param customHeader header to be put in front of the shader source
     * @param activate activate the shader after loading
     *
     * @return The loaded shader
     *
     * @throw Exception if loading failed
     */
    Shader* loadSeparate(const std::string& vertFilename, const std::string& geomFilename,
                         const std::string& fragFilename,
                         const std::string& customHeader, bool activate = true);
     /**
     * Loads shaders from files, link shader program and activate it by default.
     *
     * You have to pass the complete filenames, inclusive file extensions (".vert", ".geom", frag").
     *
     * @param filenames the filenames of all the shaders
     * @param types the types for the corresponding file
     * @param customHeader header to be put in front of the shader source
     * @param activate activate the shader after loading
     * @param setFragDataLocation set frag data location? (do not set for GLSL version >= 4)
     *
     * @return The loaded shader
     *
     * @throw Exception if loading failed
     */
    Shader* loadSeparate(const std::vector<std::string>& filenames, const std::vector<ShaderObject::ShaderType>& types,
                         const std::string& customHeader = "", bool activate = true, bool setFragDataLocation = true);

    /**
     * Transfer ownership of a static shader, i.e. long lived shader to the shader manager.
     * The shader can still be used as long as the shader manager is alive,
     * but will be automatically disposed when the shader manager is deconstructed.
     *
     * This is useful for use cases where there is no explicit deinitialization
     * function that always will be called before the shader manager is deconstructed.
     */
    void registerStaticShader(Shader*);

    bool rebuildAllShadersFromFile();

    Shader* getActiveShader() const {
        return activeShader_;
    }

private:

    // Will be called by the Shader class itself, only!
    void setActiveShader(Shader* shader) {
        activeShader_ = shader;
    }

    Shader* activeShader_;

    /**
     * Long lived shaders that are automatically disposed
     * when the shader manager is deconstructed.
     */
    std::vector<Shader*> staticShaders_;

    static const std::string loggerCat_;
};

class GLSLExtensionSet {
public:
    enum ExtensionStatus {
        ENABLED  = (1 << 0),
        REQUIRED = (1 << 1),
        WARN     = (1 << 2),
        DISABLED = (1 << 3),
    };
    GLSLExtensionSet();
    ~GLSLExtensionSet() {}
    void set(const std::string& name, ExtensionStatus status);
    std::string generateGLSL() const;
    void generateGLSL(std::vector<std::string>& output) const;
private:
    std::string format(const std::string& name, ExtensionStatus estatus) const;
    std::map<std::string, ExtensionStatus> extensions_;

    static const std::string loggerCat_;
};
class GLSLVersionCollector {
public:
    enum AdditionalQualifier {
        NONE,
        CORE,
        COMPATIBILITY,
    };

    static const int INVALID_VERSION;

    GLSLVersionCollector();
    ~GLSLVersionCollector() {}

    void addVersion(int version);
    void addAdditionalQualifier(AdditionalQualifier qualifier);
    bool isInvalid() const;

    std::string generateGLSL() const;
private:
    static std::string additionalQualifierAsString(AdditionalQualifier);

    int version_;
    AdditionalQualifier additionalQualifier_;
};
/**
 * Parses #include statements and geometry shader settings
 */
class ShaderPreprocessor {
public:
    ShaderPreprocessor(ShaderObject* obj);

    // Returns the parsed result
    std::string getResult() const;

protected:

    void parse();
    void parsePart(const std::string& input, const std::string& name = "");

    /**
     * Removes comments and leading spaces from a line and returns the result.
     */
    std::string stripLine(const std::string line, bool& inComment);

    void outputComment(const std::string& comment, const std::string& type = "INFO");

    // not used atm, but may be useful at some point FL
    std::string getShaderDirective(const std::string& d);

    ShaderObject* shd_;
    LineTracker& lineTracker_; ///< keeps track of line numbers when includes are used
    int activeLine_;
    std::ostringstream result_;
    GLSLVersionCollector maxVersion_;
    GLSLExtensionSet extensions_;

    static const std::string loggerCat_;
};

} // namespace tgt

#define ShdrMgr tgt::Singleton<tgt::ShaderManager>::getRef()

#endif //TGT_SHADERMANAGER_H
