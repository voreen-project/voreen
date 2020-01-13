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

#ifndef VRN_VOREENMODULE_H
#define VRN_VOREENMODULE_H

#include "voreen/core/voreencoreapi.h"
#include "voreen/core/io/serialization/voreenserializableobjectfactory.h"
#include "voreen/core/properties/propertyowner.h"
#include "voreen/core/utils/exception.h"

#include <string>
#include <vector>

namespace voreen {

class VoreenApplication;

class Processor;
class ProcessorWidget;
class ProcessorWidgetFactory;

class Property;
class PropertyWidget;
class PropertyWidgetFactory;

class Port;
class SerializableFactory;
class VolumeReader;
class VolumeWriter;

/**
 * Abstract base class for module classes.
 *
 * A module class is responsible for registering the
 * module's resources at runtime:
 *  - Processors
 *  - Properties
 *  - Ports
 *  - VolumeReaders and VolumeWriters
 *  - SerializableFactories
 *  - Shader search path
 *
 * The registration is to be done in the constructor.
 *
 */
class VRN_CORE_API VoreenModule : public PropertyOwner, public VoreenSerializableObjectFactory {

    friend class VoreenApplication;

public:
    /**
     * Instantiate and register the module's classes
     * in the derived class's constructor.
     *
     * @param modulePath absolute path to the module directory
     */
    VoreenModule(const std::string& modulePath);

    /**
     * Deletes the registered objects.
     */
    virtual ~VoreenModule();

    virtual std::string getClassName() const { return "VoreenModule"; }
    /// Do not use!
    virtual VoreenModule* create() const { return 0; }

    /**
     * Returns the name of the module's base directory (without parent dir).
     */
    std::string getDirName() const;

    /**
     * Returns true, if the module has been successfully initialized,
     * and has not yet been deinitialized.
     */
    bool isInitialized() const;

    /**
     * Returns true, if the module has been successfully OpenGL initialized
     * and has not yet been deinitialized.
     */
    bool isInitializedGL() const;

    /**
     * Module documentation that is intended
     * to be presented to the user.
     */
    virtual std::string getDescription() const = 0;

    /// Sets the description
    void setDescription(const std::string& description) const;

    /**
     * Constructs an absolute path consisting of the module's base directory
     * and the given path suffix.
     */
    std::string getModulePath(const std::string& suffix = "") const;

    /**
     * Returns the module's processors.
     *
     * @see registerType
     */
    std::vector<const Processor*> getRegisteredProcessors() const;

    /**
     * Returns the module's properties.
     *
     * @see registerType
     */
    std::vector<const Property*> getRegisteredProperties() const;

    /**
     * Returns the module's VolumeReaders.
     *
     * @see registerVolumeReader
     */
    const std::vector<VolumeReader*>& getRegisteredVolumeReaders() const;

    /**
     * Returns the module's VolumeWriters.
     *
     * @see registerVolumeWriter
     */
    const std::vector<VolumeWriter*>& getRegisteredVolumeWriters() const;

    /**
     * Returns the module's ProcessorWidgetFactories.
     *
     * @see registerProcessorWidgetFactory
     */
    const std::vector<ProcessorWidgetFactory*>& getRegisteredProcessorWidgetFactories() const;

    /**
     * Returns the module's PropertyWidgetFactories.
     *
     * @see registerPropertyWidgetFactory
     */
    const std::vector<PropertyWidgetFactory*>& getRegisteredPropertyWidgetFactories() const;

    /**
     * Returns the GLSL shader search used by the module.
     *
     * @see registerShaderPath
     */
    const std::vector<std::string>& getRegisteredShaderPaths() const;

protected:
    /**
     * Time-consuming initializations that do not access
     * the OpenGL API should be performed here. Registration
     * of resources, however, should be done in the constructor.
     *
     * @throw tgt::Exception to indicate an initialization failure
     */
    virtual void initialize();

    /**
     * OpenGL-dependent initializations should be performed here.
     * However, it is usually not necessary to override this function.
     *
     * @throw tgt::Exception to indicate an initialization failure
     */
    virtual void initializeGL();

    /**
     * Time-consuming deinitializations that do not access
     * the OpenGL API should be performed here.
     * However, it is usually not necessary to override this function.
     *
     * @throw tgt::Exception to indicate a deinitialization failure
     */
    virtual void deinitialize();

    /**
     * OpenGL-dependent deinitializations should be performed here.
     * However, it is usually not necessary to override this function.
     *
     * @throw tgt::Exception to indicate an deinitialization failure
     */
    virtual void deinitializeGL();

    /**
     * Specifies the module's name. To be called in the
     * derived class's constructor.
     * GuiName and ID are the same!!!
     *
     * @note Setting a name is mandatory.
     */
    void setGuiName(const std::string& guiName);

    /**
     * Specifies the module's name. To be called in the
     * derived class's constructor.
     * GuiName and ID are the same!!!
     *
     * @note Setting a name is mandatory.
     */
    void setID(const std::string& id);

    /**
     * Registers the passed VoreenSerializableObject at the module, using its class name as type name.
     *
     * @note The modules takes ownership of the passed object.
     */
    virtual void registerSerializableType(const VoreenSerializableObject* type);

    /**
     * Registers the passed VoreenSerializableObject at the module.
     *
     * @note The module takes ownership of the passed object.
     */
    virtual void registerSerializableType(const std::string& typeName, const VoreenSerializableObject* type);

    /**
     * Registers the passed processor. To be called in the
     * derived class's constructor.
     *
     * @deprecated Use registerSerializableType() instead!
     */
    void registerProcessor(Processor* processor);

    /**
     * Registers the passed property. To be called in the
     * derived class's constructor.
     *
     * @deprecated Use registerSerializableType() instead!
*/
    void registerProperty(Property* property);

    /**
     * Registers the passed VolumeReader. To be called in the
     * derived class's constructor.
     */
    void registerVolumeReader(VolumeReader* reader);

    /**
     * Registers the passed VolumeWriter. To be called in the
     * derived class's constructor.
     */
    void registerVolumeWriter(VolumeWriter* writer);

    /**
     * Registers the passed ProcessorWidgetFactory. To be called in the
     * derived class's constructor.
     */
    void registerProcessorWidgetFactory(ProcessorWidgetFactory* factory);

    /**
     * Registers the passed PropertyWidgetFactory. To be called in the
     * derived class's constructor.
     */
    void registerPropertyWidgetFactory(PropertyWidgetFactory* factory);

    /**
     * Adds the passed directory to the shader search path.
     * To be called in the derived class's constructor.
     *
     * @note Use getModulePath() for determining the module's directory.
     */
    void addShaderPath(const std::string& path);

    bool initialized_;    ///< set by VoreenApplication after initialize() has been successfully called
    bool initializedGL_;  ///< set by VoreenApplication after initializeGL() has been successfully called

    static const std::string loggerCat_;

private:
    std::string dirName_;       //< name of the module's base directory (without parent dir)
    std::string modulePath_;    //< absolute path to module directory (set by registration header)

    std::vector<VolumeReader*> volumeReaders_;
    std::vector<VolumeWriter*> volumeWriters_;
    std::vector<std::string> shaderPaths_;

    std::vector<ProcessorWidgetFactory*> processorWidgetFactories_;
    std::vector<PropertyWidgetFactory*> propertyWidgetFactories_;
};

} // namespace

#endif // VRN_VOREENMODULE_H
