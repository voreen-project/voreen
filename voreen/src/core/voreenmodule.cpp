/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

#include "voreen/core/voreenmodule.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/ports/port.h"
#include "voreen/core/utils/stringutils.h"
#include "voreen/core/processors/processor.h"
#include "voreen/core/processors/processorwidgetfactory.h"
#include "voreen/core/properties/property.h"
#include "voreen/core/properties/propertywidget.h"
#include "voreen/core/properties/propertywidgetfactory.h"
#include "voreen/core/io/volumereader.h"
#include "voreen/core/io/volumewriter.h"
#include "voreen/core/io/serialization/serializablefactory.h"

#include "tgt/shadermanager.h"
#include "tgt/filesystem.h"

namespace voreen {

const std::string VoreenModule::loggerCat_("voreen.VoreenModule");

VoreenModule::VoreenModule(const std::string& modulePath)
    : PropertyOwner("<undefined>", "<undefined>")
    , initialized_(false)
    , initializedGL_(false)
    , dirName_("<undefined>")
    , modulePath_("<undefined>")
{
    if (modulePath == "")
        LWARNING("Module path is empty");
    else {
        if (!tgt::FileSystem::isAbsolutePath(modulePath))
            modulePath_ = VoreenApplication::app()->getBasePath(modulePath); //< module path is relative to base dir
        else
            modulePath_ = modulePath; //< module path lies outside the source tree

        dirName_ = tgt::FileSystem::fileName(modulePath);

        if (!tgt::FileSystem::dirExists(modulePath_))
            LDEBUG("Module path does not exist: " + modulePath);
    }
}

VoreenModule::~VoreenModule() {

    if (isInitializedGL()) {
        LWARNING("Module '" << getID() << "' has not been OpenGL deinitialized before destruction");
    }

    if (isInitialized()) {
        LWARNING("Module '" << getID() << "' has not been deinitialized before destruction");
    }

    // free resources
    for (size_t i=0; i<volumeReaders_.size(); i++)
        delete volumeReaders_[i];
    volumeReaders_.clear();

    for (size_t i=0; i<volumeWriters_.size(); i++)
        delete volumeWriters_[i];
    volumeWriters_.clear();

    for (size_t i=0; i<processorWidgetFactories_.size(); i++)
        delete processorWidgetFactories_[i];
    processorWidgetFactories_.clear();

    for (size_t i=0; i<propertyWidgetFactories_.size(); i++)
        delete propertyWidgetFactories_[i];
    propertyWidgetFactories_.clear();
}

std::string VoreenModule::getDirName() const {
    return dirName_;
}

std::vector<const Processor*> VoreenModule::getRegisteredProcessors() const {
    return getSerializableTypes<Processor>();
}

std::vector<const Property*> VoreenModule::getRegisteredProperties() const {
    return getSerializableTypes<Property>();
}

const std::vector<VolumeReader*>& VoreenModule::getRegisteredVolumeReaders() const {
    return volumeReaders_;
}

const std::vector<VolumeWriter*>& VoreenModule::getRegisteredVolumeWriters() const {
    return volumeWriters_;
}

const std::vector<ProcessorWidgetFactory*>& VoreenModule::getRegisteredProcessorWidgetFactories() const {
    return processorWidgetFactories_;
}

const std::vector<PropertyWidgetFactory*>& VoreenModule::getRegisteredPropertyWidgetFactories() const {
    return propertyWidgetFactories_;
}

const std::vector<std::string>& VoreenModule::getRegisteredShaderPaths() const {
    return shaderPaths_;
}


// protected

void VoreenModule::initialize() {
    if (getID() == "<undefined>" || getID().empty())
        LWARNING("module id is undefined");
    if (getModulePath() == "<undefined>" || getModulePath().empty())
        LWARNING("module path is undefined");
}

void voreen::VoreenModule::deinitialize() {
}

void VoreenModule::initializeGL() {
    for (size_t i=0; i<shaderPaths_.size(); i++) {
        ShdrMgr.addPath(shaderPaths_.at(i));
    }
}

void VoreenModule::deinitializeGL() {}

void VoreenModule::setGuiName(const std::string& guiName) {
    guiName_ = guiName;
}

void VoreenModule::setID(const std::string& id) {
    id_ = id;
    if (guiName_.empty())
        guiName_ = id;
}

void VoreenModule::registerSerializableType(const VoreenSerializableObject* type) {
    tgtAssert(type, "null pointer passed");
    // if type is a Processor, set module id
    if (const Processor* proc = dynamic_cast<const Processor*>(type))
        const_cast<Processor*>(proc)->setModuleName(getID());

    VoreenSerializableObjectFactory::registerSerializableType(type);
}

void VoreenModule::registerSerializableType(const std::string& typeName, const VoreenSerializableObject* type) {
    tgtAssert(type, "null pointer passed");
    // if type is a Processor, set module id
    if (const Processor* proc = dynamic_cast<const Processor*>(type))
        const_cast<Processor*>(proc)->setModuleName(getID());

    VoreenSerializableObjectFactory::registerSerializableType(typeName, type);
}

void VoreenModule::registerProcessor(Processor* processor) {
    registerSerializableType(processor);
}

void VoreenModule::registerProperty(Property* property) {
    registerSerializableType(property);
}

void VoreenModule::registerVolumeReader(VolumeReader* reader) {
    tgtAssert(reader, "null pointer passed");
    volumeReaders_.push_back(reader);
}

void VoreenModule::registerVolumeWriter(VolumeWriter* writer) {
    tgtAssert(writer, "null pointer passed");
    volumeWriters_.push_back(writer);
}

void VoreenModule::registerProcessorWidgetFactory(ProcessorWidgetFactory* factory) {
    tgtAssert(factory, "null pointer passed");
    processorWidgetFactories_.push_back(factory);
}

void VoreenModule::registerPropertyWidgetFactory(PropertyWidgetFactory* factory) {
    tgtAssert(factory, "null pointer passed");
    propertyWidgetFactories_.push_back(factory);
}

void VoreenModule::addShaderPath(const std::string& path) {
    shaderPaths_.push_back(path);
}

std::string VoreenModule::getModulePath(const std::string& suffix) const {
    return tgt::FileSystem::cleanupPath(modulePath_ + (suffix.empty() ? "" : "/" + suffix));
}

bool VoreenModule::isInitialized() const {
    return initialized_;
}

bool VoreenModule::isInitializedGL() const {
    return initializedGL_;
}

} // namespace
