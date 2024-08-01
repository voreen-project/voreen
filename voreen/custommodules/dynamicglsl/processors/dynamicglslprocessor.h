/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#ifndef VRN_DYNAMICGLSLPROCESSOR_H
#define VRN_DYNAMICGLSLPROCESSOR_H

#include "voreen/core/processors/volumeraycaster.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/matrixproperty.h"
#include "voreen/core/processors/processorwidgetfactory.h"

#include "voreen/core/ports/volumeport.h"

#include "voreen/core/utils/GLSLparser/glslprogram.h"
#include "voreen/core/utils/GLSLparser/glslannotation.h"
#include "voreen/core/utils/GLSLparser/glsl/glslsymbol.h"

#include "tgt/shadermanager.h"

namespace voreen {

class DynamicGLSLProcessor : public VolumeRaycaster {

public:
    DynamicGLSLProcessor();
    ~DynamicGLSLProcessor();

    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "DynamicGLSLProcessor";   }
    virtual std::string getCategory() const   { return "GLSL";      }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL; }

    ShaderProperty* getShader() { return shader_; }

protected:
    virtual void setDescriptions() {
        setDescription("This is a highly customizable processor. It can be customized by using the associated GLSL shader. Within the shader it is possible to define a varying number of inports for receiving image or volume data. Furthermore, up to eight render outports can be defined. To allow further customization, properties can be initialized and used as regular properties, e.g., linking is possible. To learn how to setup ports and properties have a look at the standard shader, which is associated with the processor.\
<p><b>Note:</b> To update ports and properties the refresh button has to be pressed.</p>");
    }

    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

    /// Calls rebuild() on the shader property.
    virtual void loadShader();

    /// Assigns the generated header to the property and rebuilds the shader.
    virtual void compile();

    virtual std::string generateHeaderFromVolume(Volume* volumeHandle = 0);

    virtual void setupPortGroup();
    virtual void initializePortsAndProperties();

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

private:

    struct InportStruct {
        InportStruct();
        InportStruct(const std::string& name, const std::string& mainSamplerIdentifier_,
            const std::string& depthSamplerIdentifier, const std::string& paramsIdentifier);
        InportStruct(const std::string& name, const std::string& mainSamplerIdentifier_,
            const std::string& depthSamplerIdentifier, const std::string& paramsIdentifier,
            GLint mainTexUnit, GLint depthTexUnit);

        std::string name_;
        std::string mainSamplerIdentifier_;
        std::string depthSamplerIdentifier_;
        std::string paramsIdentifier_;
        GLint mainTexUnit_;
        GLint depthTexUnit_;
    };

    std::vector<InportStruct> inportIDs_;
    std::vector<std::string> outportIDs_;
    std::vector<std::string> propertyIDs_;
    std::vector<std::string> baseClassProperties_;

    void addNewInport(glslparser::GLSLVariableSymbol* symbol);
    void removeOldPort(Port* port);
    void updateInports(glslparser::GLSLVariableSymbol* symbol);
    void updateOutports(glslparser::GLSLVariableSymbol* symbol);
    GLint getColorFormat(glslparser::GLSLVariableSymbol* symbol);

    void addNewProperty(glslparser::GLSLVariableSymbol* symbol);
    void removeOldProperty(Property* property);
    void updatePropertyValues(Property* property, glslparser::GLSLVariableSymbol* symbol);
    void updateProperties(glslparser::GLSLVariableSymbol* symbol);

    PortGroup portGroup_;
    ShaderProperty* shader_;

    static const std::string loggerCat_;
};

} // namespace

#endif // VRN_DYNAMICGLSLPROCESSOR_H
