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

#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include <exception>

#include "voreen/core/voreenapplication.h"
#include "tgt/filesystem.h"
#include "voreen/core/utils/stringutils.h"

#include "itk_wrapper.h"
#include "template.h"

using namespace voreen;

const std::string ITK_Wrapper::loggerCat_("voreen.ITK_Wrapper");

ITK_Wrapper::ITK_Wrapper() {}

ITK_Wrapper::ITK_Wrapper(std::vector<std::string> modules) {
    moduleXMLs_ = modules;
}

void ITK_Wrapper::initializeConstants() {
    codeStates_.push_back("EXPERIMENTAL");
    codeStates_.push_back("TESTING");
    codeStates_.push_back("STABLE");

    integerScalarTypes_.push_back("UInt8");
    integerScalarTypes_.push_back("Int8");
    integerScalarTypes_.push_back("UInt16");
    integerScalarTypes_.push_back("Int16");
    integerScalarTypes_.push_back("UInt32");
    integerScalarTypes_.push_back("Int32");
    //integerScalarTypes_.push_back("UInt64");
    //integerScalarTypes_.push_back("Int64");

    realScalarTypes_.push_back("Float");
    realScalarTypes_.push_back("Double");

    scalarTypes_ = integerScalarTypes_;
    scalarTypes_.push_back("Float");
    scalarTypes_.push_back("Double");

    integerVectorTypes_.push_back("4xUInt8");
    integerVectorTypes_.push_back("4xInt8");
    integerVectorTypes_.push_back("4xUInt16");
    integerVectorTypes_.push_back("4xInt16");
    integerVectorTypes_.push_back("3xUInt8");
    integerVectorTypes_.push_back("3xInt8");
    integerVectorTypes_.push_back("3xUInt16");
    integerVectorTypes_.push_back("3xInt16");
    integerVectorTypes_.push_back("2xUInt8");
    integerVectorTypes_.push_back("2xInt8");
    integerVectorTypes_.push_back("2xUInt16");
    integerVectorTypes_.push_back("2xInt16");

    realVectorTypes_.push_back("4xFloat");
    realVectorTypes_.push_back("4xDouble");
    realVectorTypes_.push_back("3xFloat");
    realVectorTypes_.push_back("3xDouble");
    realVectorTypes_.push_back("2xFloat");
    realVectorTypes_.push_back("2xDouble");

    vectorTypes_.push_back("4xUInt8");
    vectorTypes_.push_back("4xInt8");
    vectorTypes_.push_back("4xUInt16");
    vectorTypes_.push_back("4xInt16");
    vectorTypes_.push_back("4xFloat");
    vectorTypes_.push_back("4xDouble");
    vectorTypes_.push_back("3xUInt8");
    vectorTypes_.push_back("3xInt8");
    vectorTypes_.push_back("3xUInt16");
    vectorTypes_.push_back("3xInt16");
    vectorTypes_.push_back("3xFloat");
    vectorTypes_.push_back("3xDouble");
    vectorTypes_.push_back("2xUInt8");
    vectorTypes_.push_back("2xInt8");
    vectorTypes_.push_back("2xUInt16");
    vectorTypes_.push_back("2xInt16");
    vectorTypes_.push_back("2xFloat");
    vectorTypes_.push_back("2xDouble");

    allTypes_ = integerScalarTypes_;
    allTypes_.push_back("UInt64");
    allTypes_.push_back("Int64");
    allTypes_.push_back("Float");
    allTypes_.push_back("Double");
    for (size_t i=0; i<vectorTypes_.size(); ++i) {
        allTypes_.push_back(vectorTypes_[i]);
    }

    metaTypes_.push_back("Scalar");
    metaTypes_.push_back("IntegerScalar");
    metaTypes_.push_back("RealScalar");
    metaTypes_.push_back("Vector");
    metaTypes_.push_back("IntegerVector");
    metaTypes_.push_back("RealVector");

    argumentTypes_.push_back("Int");
    argumentTypes_.push_back("Float");
    argumentTypes_.push_back("Bool");
    argumentTypes_.push_back("KernelType");
    argumentTypes_.push_back("SizeType");
    argumentTypes_.push_back("IndexType");
    argumentTypes_.push_back("ArrayType");
    argumentTypes_.push_back("PixelType");
    argumentTypes_.push_back("AddSeedType");
    argumentTypes_.push_back("SetSeedType");
    argumentTypes_.push_back("VectorSeedType");
}

void ITK_Wrapper::wrapITKtoVoreen() {
    std::string filename;
    if (moduleXMLs_.empty()) {
        //take all modules
        moduleXMLs_ = tgt::FileSystem::listFiles(getXMLPath(), false);

        //delete all current files and directories
        tgt::FileSystem::deleteDirectoryRecursive(getProcessorOutputPath().c_str());
    }

    //deserializes all xml-files in moduleXMLs_ (and checks each module with checkFilters-Method)
    deserializeModuleXMLs();

    //delete old directories and create new
    for (size_t i=0; i<filtermodules_.size(); ++i) {
        tgt::FileSystem::deleteDirectoryRecursive(getProcessorOutputPath(filtermodules_[i].getName()).c_str());
        tgt::FileSystem::createDirectoryRecursive(getProcessorOutputPath(filtermodules_[i].getName()).c_str());
    }

    createITKGeneratedModule_CPP();
    createITKGenerated_CMake();
    createFilters();
}

void ITK_Wrapper::deserializeModuleXMLs() {
    ITK_Module itkModule;
    std::fstream stream;
    std::string filename;

    for (size_t i=0; i<moduleXMLs_.size(); ++i){
        XmlDeserializer de;
        de.setUseAttributes(true);
        stream.open(getXMLPath(moduleXMLs_[i]).c_str(), std::ios::in);
        de.read(stream);
        de.deserialize("ITK_Module", itkModule);
        checkFilters(itkModule);
        filtermodules_.push_back(itkModule);
        stream.clear();
        stream.close();
    }
}

void ITK_Wrapper::checkFilters(ITK_Module itkModule) {
    Filter filter;
    Ports port;
    Argument argument;
    std::string type;
    std::string errorMessage;

    // for each filter in module itkModule ...
    for (size_t i=0; i<itkModule.getFilters().size(); ++i) {
        filter = itkModule.getFilters()[i];
        errorMessage = "";
        errorMessage = errorMessage + "\nError in " + itkModule.getName() + ".xml!\n"
            + "filter: " + filter.getName() + "\n";

        if (!filter.isEnabled())
            continue;
        // ... ckeck codeState
        if (std::find(codeStates_.begin(), codeStates_.end(), filter.getCodeState()) == codeStates_.end()) {
            errorMessage = errorMessage + "codeState: " + filter.getCodeState() + "\n\n"
                + "CodeState is not supported! Supported CodeStates are:\n"
                + strJoin(codeStates_, ", ") + ".\n"
                + "Please change xml-File!";
            throw VoreenException(errorMessage);
        }

        // ... check for each inport ...
        for (size_t j=0; j<filter.getInports().size(); ++j) {
            port = filter.getInports()[j];
            // ... if it has a name ...
            if (port.getName().empty()) {
                errorMessage = errorMessage + "No inport name set!\n"
                    + "Please change xml-File!";
                throw VoreenException(errorMessage);
            }
            // ... and if its possibleTypes_ make sense
            for (size_t k=0; k<port.getPossibleTypes().size(); ++k) {
                type = port.getPossibleTypes()[k];
                if (std::find(allTypes_.begin(), allTypes_.end(), type) == allTypes_.end()
                    && std::find(metaTypes_.begin(), metaTypes_.end(), type) == metaTypes_.end()) {

                    errorMessage = errorMessage + "inport: " + port.getName() + "\n"
                        + "possibleType: " + type + "\n\n"
                        + "Volumetype is not supported! Supported types are:\n"
                        + strJoin(allTypes_, ", ") + ",\n"
                        + "or the metaTypes:\n"
                        + strJoin(metaTypes_, ", ") + ".\n"
                        + "Please change xml-File!";
                    throw VoreenException(errorMessage);
                }
            }
        }

        // ... check for each outport ...
        for (size_t j=0; j<filter.getOutports().size(); ++j) {
            port = filter.getOutports()[j];
            // ... if it has a name, ...
            if (port.getName().empty()) {
                errorMessage = errorMessage + "No ouport name set!\n"
                    + "Please change xml-File!";
                throw VoreenException(errorMessage);
            }
            if (!port.getPossibleTypes().empty()) {
                // ... if it has only one possible type ...
                if (port.getPossibleTypes().size() > 1) {
                    errorMessage = errorMessage + "outport: " + port.getName() + "\n"
                        + "Only one possibleType supported!\n"
                        + "Please change xml-File!";
                    throw VoreenException(errorMessage);
                }
                else {
                    type = port.getPossibleTypes()[0];
                    // ... and if this type makes sense
                    if (std::find(allTypes_.begin(), allTypes_.end(), type) == allTypes_.end()) {
                        errorMessage = errorMessage + "outport: " + port.getName()+ "\n"
                            + "possibleType: " + type + "\n\n"
                            + "Volumetype is not supported! Supported types are:\n"
                            + strJoin(allTypes_, ", ") + ". \n"
                            + "Please change xml-File!";
                        throw VoreenException(errorMessage);
                    }
                }
            }
        }

        // ... check for each argument ...
        for(size_t j=0; j<filter.getArguments().size(); ++j) {
            argument = filter.getArguments()[j];

            // ... if it has a name ...
            if (argument.getName().empty()) {
                errorMessage = errorMessage + "No argument name set!\n"
                    + "Please change xml-File!";
                throw VoreenException(errorMessage);
            }

            type = argument.getType();
            // ... and if the argument-type is supported
            if (std::find(argumentTypes_.begin(), argumentTypes_.end(), type) == argumentTypes_.end()) {
                errorMessage = errorMessage + "argument: " + argument.getName() + "\n"
                    + "type: " + type + "\n\n"
                    + "Type is not supported! Supported types are:\n"
                    + strJoin(argumentTypes_, ", ") + ".\n"
                    + "Please change xml-File!";
                throw VoreenException(errorMessage);
            }
            else if (type != "KernelType" && type != "PixelType" && type != "AddSeedType"
                && type != "SetSeedType" && type != "VectorSeedType" && argument.isInputArgument()) {
                // for some types a default value must be set
                    if (argument.getDefaultValue().empty()) {
                    errorMessage = errorMessage + "argument: " + argument.getName() + "\n\n"
                        + "No defaultValue set!\n"
                        + "Please change xml-File!";
                    throw VoreenException(errorMessage);
                }
                else if (type != "Bool") {
                    // for some types a min value must be set
                    if (argument.getMinValue().empty()) {
                        errorMessage = errorMessage + "argument: " + argument.getName() + "\n\n"
                            + "No minValue set!\n"
                            + "Please change xml-File!";
                        throw VoreenException(errorMessage);
                    }
                    // for some types a max value must be set
                    if (argument.getMaxValue().empty()) {
                        errorMessage = errorMessage + "argument: " + argument.getName() + "\n\n"
                            + "No maxValue set!\n"
                            + "Please change xml-File!";
                        throw VoreenException(errorMessage);
                    }
                }
            }
        }
    }
}

void ITK_Wrapper::createITKGeneratedModule_CPP() {
    ITK_Module itkModule;
    std::string moduleName;
    Filter filter;
    // replacement for §FILTERHEADER§ in template "itk_generatedmodule.txt"
    std::string rep_filterHeader;
    // replacement for §FILTERADD§ in template "itk_generatedmodule.txt"
    std::string rep_filterAdd;

    for (size_t f=0; f<filtermodules_.size(); ++f) {
        itkModule = filtermodules_[f];
        moduleName = itkModule.getName();

        for (size_t i=0; i<itkModule.getFilters().size(); ++i) {
            filter = itkModule.getFilters()[i];
            if (filter.isEnabled()) {
                rep_filterHeader = rep_filterHeader + "#include \"processors/" + moduleName
                    + "/" + filter.getFilename() + ".h\"\n";
                rep_filterAdd = rep_filterAdd + "    registerProcessor(new " + filter.getName() + "ITK());\n";
            }
        }
    }

    Template t(getTemplatePath("itk_generatedmodule.txt"));
    t.replace("FILTERHEADER", rep_filterHeader);
    t.replace("FILTERADD", rep_filterAdd);

    t.writeToFile(VoreenApplication::app()->getBasePath("modules/itk_generated/itk_generatedmodule.cpp"));
}

void voreen::ITK_Wrapper::createITKGenerated_CMake() {
    ITK_Module itkModule;
    std::string moduleName;
    Filter filter;
    // replacement for §FILTERCPP§ in template "itk_generated_core.txt"
    std::string rep_filterCPP;
    // replacement for §FILTERH§ in template "itk_generated_core.txt"
    std::string rep_filterH;

    for (size_t f=0; f<filtermodules_.size(); ++f) {
        itkModule = filtermodules_[f];
        moduleName = itkModule.getName();

        for (size_t i=0; i<itkModule.getFilters().size(); ++i) {
            filter = itkModule.getFilters()[i];
            if (filter.isEnabled()) {
                rep_filterCPP = rep_filterCPP
                    + "    ${MOD_DIR}/processors/" + moduleName
                    + "/" + filter.getFilename() + ".cpp \n";
                rep_filterH = rep_filterH
                    + "    ${MOD_DIR}/processors/" + moduleName
                    + "/" + filter.getFilename() + ".h \n";
            }
        }
    }

    Template t(getTemplatePath("itk_generated_cmake.txt"));
    t.replace("§FILTERCPP§", rep_filterCPP);
    t.replace("§FILTERH§", rep_filterH);

    t.writeToFile(VoreenApplication::app()->getBasePath("modules/itk_generated/itk_generated.cmake"));
}

void ITK_Wrapper::createFilters() {
    ITK_Module itkModule;
    Filter filter;
    size_t totalFilterCount = 0;
    size_t createdFilterCount = 0;
    size_t filterCount = 0;
    std::string infoMessage;
    std::string filename;

    for (size_t f=0; f<filtermodules_.size(); ++f) {
        itkModule = filtermodules_[f];
        filterCount = 0;
        infoMessage = "";
        infoMessage = infoMessage + "\nGenerate module " + itkModule.getName() + " ...\n";

        for (size_t i=0; i<itkModule.getFilters().size(); ++i) {
            filter = itkModule.getFilters()[i];
            if (filter.isEnabled()) {
                ++filterCount;
                // if a filter can be generated by the wrapper, initialize its ports and create
                // the .cpp-file and the .h-file
                if (filter.isAutoGenerated()) {
                    filter.initializeFilterPorts(integerScalarTypes_, realScalarTypes_
                        , scalarTypes_, integerVectorTypes_, realVectorTypes_, vectorTypes_);
                    filter.createCPP(itkModule.getName());
                    filter.createHeaderfile(itkModule.getName(), itkModule.getGroup());
                }
                // else "copy" its hand-coded .cpp-file and its .h-file from the specialFilters-
                // directory into its module directory
                else {
                    filename = "";
                    filename = filename + filter.getFilename() + ".cpp";
                    Template t1(getSpecialFilterPath(filename));
                    t1.writeToFile(filter.getProcessorOutputPath(itkModule.getName(), filename));
                    filename = "";
                    filename = filename + filter.getFilename() + ".h";
                    Template t2(getSpecialFilterPath(filename));
                    t2.writeToFile(filter.getProcessorOutputPath(itkModule.getName(), filename));
                }
            }
            else {
                infoMessage = infoMessage + "    - " + filter.getName() + " is disabled.\n"
                    + "      Description: " + filter.getDescription() + "\n";
            }
        }

        infoMessage = infoMessage + itos(filterCount) + " of " + itos(itkModule.getFilters().size())
            + " filters in this module generated. " + "\n\n";
        LINFO(infoMessage);
        totalFilterCount = totalFilterCount + itkModule.getFilters().size();
        createdFilterCount = createdFilterCount + filterCount;
    }
    infoMessage = "";
    infoMessage = infoMessage + "\nFinished to generate " + itos(createdFilterCount) + " of "
        + itos(totalFilterCount) + " filters. \n";
    LINFO(infoMessage);
}

std::string ITK_Wrapper::getXMLPath() {
    std::string path;
    path = path + VoreenApplication::app()->getBasePath() + "/apps/itk_wrapper/xml_Files";
    return path;
}

std::string ITK_Wrapper::getXMLPath(std::string xmlFilename) {
    std::string path;
    path = path + VoreenApplication::app()->getBasePath() + "/apps/itk_wrapper/xml_Files/" + xmlFilename;
    return path;
}

std::string ITK_Wrapper::getProcessorOutputPath() {
    std::string path;
    path = path + VoreenApplication::app()->getBasePath("modules/itk_generated/processors");
    return path;
}

std::string ITK_Wrapper::getProcessorOutputPath(std::string moduleName) {
    std::string path;
    path = path + VoreenApplication::app()->getBasePath("modules/itk_generated/processors/" + moduleName);
    return path;
}

std::string ITK_Wrapper::getTemplatePath(std::string templateFilename) {
    std::string path;
    path = path + VoreenApplication::app()->getBasePath() + "/apps/itk_wrapper/templates/" + templateFilename;
    return path;
}

std::string ITK_Wrapper::getSpecialFilterPath(std::string specialFilterFilename) {
    std::string path;
    path = path + VoreenApplication::app()->getBasePath() + "/apps/itk_wrapper/specialFilters/" + specialFilterFilename;
    return path;
}



/**
 * Generate Filters.
 */
int main(int argc, char** argv) {
    VoreenApplication app("itk_wrapper", "itk_wrapper", "Wraps the ITK library with Voreen processors, generating the module 'itk_generated'", argc, argv);
    app.initialize();
    std::cout << "Generate Filters ..." << std::endl << std::endl;

    // To generate only special modules:
    /*
    std::vector<std::string> modules;
    modules.push_back("moduleName.xml");
    ...
    ITK_Wrapper wrapper(modules);
    */

    ITK_Wrapper wrapper;
    wrapper.initializeConstants();
    try {
        wrapper.wrapITKtoVoreen();
    }
    catch (VoreenException e){
        LERRORC("voreen.ITK_Wrapper", e.what());
        std::cin.get();
        return 1;
    }
    std::cout << std::endl << "Press Enter to finish..." << std::endl;
    std::cin.get();
    app.deinitialize();
    return 0;
}
