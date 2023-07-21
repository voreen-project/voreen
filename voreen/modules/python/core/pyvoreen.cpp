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

// Ignore "dereferencing type-punned pointer will break strict-aliasing rules" warnings on gcc
// caused by Py_RETURN_TRUE and such. The problem lies in Python so we don't want the warning
// here.
#if (__GNUC__ > 4 || (__GNUC__ == 4 && __GNUC_MINOR__ >= 2))
#pragma GCC diagnostic ignored "-Wstrict-aliasing"
#endif

// Do this at very first
#include <Python.h>

#include "pyvoreenobjects.h"

#include "../pythonmodule.h"

#include "pyvoreen.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/version.h"

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/boundingboxproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/cameraproperty.h"
//#include "voreen/core/properties/colorproperty.h" // covered by FloatVec4Property
//#include "voreen/core/properties/filedialogproperty.h" // covered by StringProperty
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/fontproperty.h"
#include "voreen/core/properties/intproperty.h"
//#include "voreen/core/properties/lightsourceproperty.h" // covered by FloatVec4Property
#include "voreen/core/properties/matrixproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/planeproperty.h"
#include "voreen/core/properties/propertyvector.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/stringproperty.h"
//#include "voreen/core/properties/temppathproperty.h" // covered by FileDialogProperty
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/volumeurllistproperty.h"
//#include "voreen/core/properties/volumeurlproperty.h" // covered by StringProperty
//#include "voreen/core/properties/color/colorswitchproperty.h" // covered by ColorProperty
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"

#include "voreen/core/ports/volumeport.h"
#include "voreen/core/datastructures/volume/volumefactory.h"

#include "voreen/core/datastructures/transfunc/1d/1dkeys/transfunc1dkeys.h"
#include "voreen/core/network/processornetwork.h"

#include "voreen/core/processors/processor.h"
#include "voreen/core/network/networkevaluator.h"
#include "voreen/core/processors/processorwidget.h"

// core module is always available
#include "modules/core/processors/input/volumesource.h"
#include "modules/core/processors/input/volumelistsource.h"
#include "modules/core/processors/output/canvasrenderer.h"

#ifdef VRN_MODULE_BASE
#include "modules/base/processors/entryexitpoints/meshentryexitpoints.h"
#include "modules/base/processors/utility/clockprocessor.h"
#endif
#include "voreen/core/interaction/voreentrackball.h"

#include "tgt/init.h"
#include "tgt/glcontextmanager.h"


//-------------------------------------------------------------------------------------------------
// internal helper functions

namespace {



/**
 * Retrieves the current processor network.
 *
 * @param functionName name of the calling function, e.g. "setFloatProperty" (included in Python exception)
 */
voreen::ProcessorNetwork* getProcessorNetwork(const std::string& functionName);

/**
 * Retrieves the processor with the specified name from the network.
 */
voreen::Processor* getProcessor(const std::string& processorName, const std::string& functionName);

/**
 * Retrieves the processor with the specified name, if it matches the template parameter.
 */
template<typename T>
T* getTypedProcessor(const std::string& processorName, const std::string& processorTypeString,
                     const std::string& functionName);

/**
 * Retrieves a property with the specified ID of a certain processor.
 */
voreen::Property* getProperty(const std::string& processorName, const std::string& propertyID,
                              const std::string& functionName);

/**
 * Retrieves the property with the specified ID, if it matches the template parameter.
 */
template<typename T>
T* getTypedProperty(const std::string& processorName, const std::string& propertyID,
                    const std::string& propertyTypeString, const std::string& functionName);

/**
 * Assigns the passed value to the TemplateProperty with the specified ID
 * that is owned by the processor with the specified name.
 *
 * @tparam PropertyType the type of the property the value is to be assigned to
 * @tparam ValueType the type of the value to assign
 *
 * @param processorName the name of the processor that owns the property
 * @param propertyID the id of the property to modify
 * @param value the value to set
 * @param propertyTypeString the value's typename, e.g. "Float" (included in Python exception)
 * @param functionName name of the calling function, e.g. "setFloatProperty" (included in Python exception)
 *
 * Failure cases:
 *  - if processor or property do not exist, a PyExc_NameError is raised
 *  - if property type does not match, a PyExc_TypeError is raised
 *  - if property validation fails, a PyExc_ValueError with the corresponding validation message is raised
 *
 * @return true if property manipulation has been successful
 */
template<typename PropertyType, typename ValueType>
bool setPropertyValue(const std::string& processorName, const std::string& propertyID, const ValueType& value,
                      const std::string& propertyTypeString, const std::string& functionName);

/**
 * Assigns the passed value to the passed TemplateProperty.
 *
 * @tparam PropertyType the type of the property the value is to be assigned to
 * @tparam ValueType the type of the value to assign
 *
 * @param property the property to manipulate
 * @param value the value to set
 * @param functionName name of the calling function, e.g. "setFloatProperty" (included in Python exception)
 *
 * Failure cases:
 *  - if property validation fails, a PyExc_ValueError with the corresponding validation message is raised
 *
 * @return true if property manipulation has been successful
 */
template<typename PropertyType, typename ValueType>
bool setPropertyValue(PropertyType* property, const ValueType& value,
                      const std::string& functionName);

/**
 * Retrieves a port with the specified ID of a certain processor.
 */
voreen::Port* getPort(const std::string& processorName, const std::string& portID,
                      const std::string& functionName);

/**
 * Retrieves the port with the specified ID, if it matches the template parameter.
 */
template<typename T>
T* getTypedPort(const std::string& processorName, const std::string& portID,
                const std::string& portTypeString, const std::string& functionName);

/**
 * Uses the apihelper.py script to print documentation
 * about the module's functions.
 */
static PyObject* printModuleInfo(const std::string& moduleName, bool omitFunctionName = false,
                                 int spacing = 0, bool collapse = false, bool blanklines = false);

} // namespace anonymous


//-------------------------------------------------------------------------------------------------
// definitions of Python binding methods

using namespace voreen;

//
// Python module 'voreen.network'
//
static PyObject* voreen_setPropertyValue(PyObject* /*self*/, PyObject* args) {

    // check length of tuple
    if (PyTuple_Size(args) != 3) {
        std::ostringstream errStr;
        errStr << "setPropertyValue() takes exactly 3 arguments: processor name, property id, value";
        errStr << " (" << PyTuple_Size(args) << " given)";
        PyErr_SetString(PyExc_TypeError, errStr.str().c_str());
        return 0;
    }

    // check parameter 1 and 2, if they are strings
    if (!PyUnicode_Check(PyTuple_GetItem(args, 0)) || !PyUnicode_Check(PyTuple_GetItem(args, 1))) {
        PyErr_SetString(PyExc_TypeError, "setPropertyValue() arguments 1 and 2 must be strings");
        return 0;
    }

    // read processor name and property id
    std::string processorName_str = PyUnicodeAsString(PyTuple_GetItem(args, 0));
    std::string propertyID_str = PyUnicodeAsString(PyTuple_GetItem(args, 1));
    if (processorName_str.empty() || propertyID_str.empty()) {
        PyErr_SetString(PyExc_TypeError, "setPropertyValue() arguments 1 and 2 must be strings");
        return 0;
    }

    char* processorName;
    char* propertyID;

    // fetch property
    Property* property = getProperty(processorName_str, propertyID_str, "setPropertyValue");
    if (!property)
        return 0;

    // determine property type, convert and assign value
    if (BoolProperty* typedProp = dynamic_cast<BoolProperty*>(property)) {
        char value;
        if (!PyArg_ParseTuple(args, "ssb:setPropertyValue", &processorName, &propertyID, &value))
            return 0;
        if (setPropertyValue<BoolProperty, bool>(typedProp, (bool)value, "setPropertyValue"))
            Py_RETURN_NONE;
    }
    else if (ButtonProperty* typedProp = dynamic_cast<ButtonProperty*>(property)) {
        // directly trigger button property without reading passed value
        typedProp->clicked();
        Py_RETURN_NONE;
    }
    else if (CameraProperty* typedProp = dynamic_cast<CameraProperty*>(property)) {
        tgt::vec3 position, focus, up;
        if (!PyArg_ParseTuple(args, "ss((fff)(fff)(fff)):setPropertyValue",
                              &processorName, &propertyID,
                              &position.x, &position.y, &position.z,
                              &focus.x, &focus.y, &focus.z,
                              &up.x, &up.y, &up.z))
            return 0;

        typedProp->setPosition(position);
        typedProp->setFocus(focus);
        typedProp->setUpVector(up);
        Py_RETURN_NONE;
    }
    else if (FloatProperty* typedProp = dynamic_cast<FloatProperty*>(property)) {
        float value;
        if (!PyArg_ParseTuple(args, "ssf:setPropertyValue", &processorName, &propertyID, &value))
            return 0;
        if (setPropertyValue<FloatProperty, float>(typedProp, value, "setPropertyValue"))
            Py_RETURN_NONE;
    }
    else if (IntProperty* typedProp = dynamic_cast<IntProperty*>(property)) {
        int value;
        if (!PyArg_ParseTuple(args, "ssi:setPropertyValue", &processorName, &propertyID, &value))
            return 0;
        if (setPropertyValue<IntProperty, int>(typedProp, value, "setPropertyValue"))
            Py_RETURN_NONE;
    }
    else if (PlaneProperty* typedProp = dynamic_cast<PlaneProperty*>(property)) {
        tgt::plane value;
        if (!PyArg_ParseTuple(args, "ss(fff)f:setPropertyValue", &processorName, &propertyID,
                              &value.n.x, &value.n.y, &value.n.z, &value.d))
            return 0;
        if (setPropertyValue<PlaneProperty, tgt::plane>(typedProp, value, "setPropertyValue"))
            Py_RETURN_NONE;
    }
    else if (ShaderProperty* typedProp = dynamic_cast<ShaderProperty*>(property)) {
        char* vertexFilename;
        char* geometryFilename;
        char* fragmentFilename;

        if (!PyArg_ParseTuple(args, "sssss:setPropertyValue", &processorName, &propertyID,
                              &vertexFilename, &geometryFilename, &fragmentFilename))
            return 0;

        ShaderSource value(ShaderFileList
                           (tgt::ShaderObject::VERTEX_SHADER, vertexFilename)
                           (tgt::ShaderObject::GEOMETRY_SHADER, geometryFilename)
                           (tgt::ShaderObject::FRAGMENT_SHADER, fragmentFilename));

        if (setPropertyValue<ShaderProperty, ShaderSource>(typedProp, value, "setPropertyValue"))
            Py_RETURN_NONE;
    }
    else if (StringProperty* typedProp = dynamic_cast<StringProperty*>(property)) {
        char* value;
        if (!PyArg_ParseTuple(args, "sss:setPropertyValue", &processorName, &propertyID, &value))
            return 0;

        if (setPropertyValue<StringProperty, std::string>(typedProp, std::string(value), "setPropertyValue"))
            Py_RETURN_NONE;
    }
    else if (OptionPropertyBase* typedProp = dynamic_cast<OptionPropertyBase*>(property)) {
        char* value;
        if (!PyArg_ParseTuple(args, "sss:setPropertyValue", &processorName, &propertyID, &value))
            return 0;
        if (setPropertyValue<OptionPropertyBase, std::string>(typedProp, std::string(value), "setPropertyValue"))
            Py_RETURN_NONE;
    }
    else if (IntBoundingBoxProperty* typedProp = dynamic_cast<IntBoundingBoxProperty*>(property)) {
        tgt::ivec3 llf, urb;
        if (!PyArg_ParseTuple(args, "ss((iii)(iii)):setPropertyValue", &processorName, &propertyID,
                              &llf.x, &llf.y, &llf.z, &urb.x, &urb.y, &urb.z))
            return 0;

        tgt::IntBounds value(llf, urb);
        if (setPropertyValue<IntBoundingBoxProperty, tgt::IntBounds>(typedProp, value, "setPropertyValue"))
            Py_RETURN_NONE;
    }
    else if (FloatBoundingBoxProperty* typedProp = dynamic_cast<FloatBoundingBoxProperty*>(property)) {
        tgt::vec3 llf, urb;
        if (!PyArg_ParseTuple(args, "ss((fff)(fff)):setPropertyValue", &processorName, &propertyID,
                              &llf.x, &llf.y, &llf.z, &urb.x, &urb.y, &urb.z))
            return 0;

        tgt::Bounds value(llf, urb);
        if (setPropertyValue<FloatBoundingBoxProperty, tgt::Bounds>(typedProp, value, "setPropertyValue"))
            Py_RETURN_NONE;
    }
    else if (IntVec2Property* typedProp = dynamic_cast<IntVec2Property*>(property)) {
        tgt::ivec2 value;
        if (!PyArg_ParseTuple(args, "ss(ii):setPropertyValue", &processorName, &propertyID,
                &value.x, &value.y))
            return 0;
        if (setPropertyValue<IntVec2Property, tgt::ivec2>(typedProp, value, "setPropertyValue"))
            Py_RETURN_NONE;
    }
    else if (IntVec3Property* typedProp = dynamic_cast<IntVec3Property*>(property)) {
        tgt::ivec3 value;
        if (!PyArg_ParseTuple(args, "ss(iii):setPropertyValue", &processorName, &propertyID,
                &value.x, &value.y, &value.z))
            return 0;
        if (setPropertyValue<IntVec3Property, tgt::ivec3>(typedProp, value, "setPropertyValue"))
            Py_RETURN_NONE;
    }
    else if (IntVec4Property* typedProp = dynamic_cast<IntVec4Property*>(property)) {
        tgt::ivec4 value;
        if (!PyArg_ParseTuple(args, "ss(iiii):setPropertyValue", &processorName, &propertyID,
                &value.x, &value.y, &value.z, &value.w))
            return 0;
        if (setPropertyValue<IntVec4Property, tgt::ivec4>(typedProp, value, "setPropertyValue"))
            Py_RETURN_NONE;
    }
    else if (FloatVec2Property* typedProp = dynamic_cast<FloatVec2Property*>(property)) {
        tgt::vec2 value;
        if (!PyArg_ParseTuple(args, "ss(ff):setPropertyValue", &processorName, &propertyID,
                &value.x, &value.y))
            return 0;
        if (setPropertyValue<FloatVec2Property, tgt::vec2>(typedProp, value, "setPropertyValue"))
            Py_RETURN_NONE;
    }
    else if (FloatVec3Property* typedProp = dynamic_cast<FloatVec3Property*>(property)) {
        tgt::vec3 value;
        if (!PyArg_ParseTuple(args, "ss(fff):setPropertyValue", &processorName, &propertyID,
                &value.x, &value.y, &value.z))
            return 0;
        if (setPropertyValue<FloatVec3Property, tgt::vec3>(typedProp, value, "setPropertyValue"))
            Py_RETURN_NONE;
    }
    else if (FloatVec4Property* typedProp = dynamic_cast<FloatVec4Property*>(property)) {
        tgt::vec4 value;
        if (!PyArg_ParseTuple(args, "ss(ffff):setPropertyValue", &processorName, &propertyID,
                &value.x, &value.y, &value.z, &value.w))
            return 0;
        if (setPropertyValue<FloatVec4Property, tgt::vec4>(typedProp, value, "setPropertyValue"))
            Py_RETURN_NONE;
    }
    else if (FloatMat2Property* typedProp = dynamic_cast<FloatMat2Property*>(property)) {
        tgt::vec2 vec0, vec1;
        char *processorName, *propertyID;
        if (!PyArg_ParseTuple(args, "ss((ff)(ff)):setPropertyValue",
                &processorName, &propertyID,
                &vec0.x, &vec0.y, &vec1.x, &vec1.y))
            return 0;

        tgt::mat2 value(vec0, vec1);
        if (setPropertyValue<FloatMat2Property, tgt::mat2>(typedProp, value, "setPropertyValue"))
            Py_RETURN_NONE;
    }
    else if (FloatMat3Property* typedProp = dynamic_cast<FloatMat3Property*>(property)) {
        tgt::vec3 vec0, vec1, vec2;
        char *processorName, *propertyID;
        if (!PyArg_ParseTuple(args, "ss((fff)(fff)(fff)):setPropertyValue",
                &processorName, &propertyID,
                &vec0.x, &vec0.y, &vec0.z, &vec1.x, &vec1.y, &vec1.z,
                &vec2.x, &vec2.y, &vec2.z))
            return 0;

        tgt::mat3 value(vec0, vec1, vec2);
        if (setPropertyValue<FloatMat3Property, tgt::mat3>(typedProp, value, "setPropertyValue"))
            Py_RETURN_NONE;
    }
    else if (FloatMat4Property* typedProp = dynamic_cast<FloatMat4Property*>(property)) {
        tgt::vec4 vec0, vec1, vec2, vec3;
        char *processorName, *propertyID;
        if (!PyArg_ParseTuple(args, "ss((ffff)(ffff)(ffff)(ffff)):setPropertyValue",
                &processorName, &propertyID,
                &vec0.x, &vec0.y, &vec0.z, &vec0.w, &vec1.x, &vec1.y, &vec1.z, &vec1.w,
                &vec2.x, &vec2.y, &vec2.z, &vec2.w, &vec3.x, &vec3.y, &vec3.z, &vec3.w))
            return 0;

        tgt::mat4 value(vec0, vec1, vec2, vec3);
        if (setPropertyValue<FloatMat4Property, tgt::mat4>(typedProp, value, "setPropertyValue"))
            Py_RETURN_NONE;
    }


    // we only get here, if property value assignment has failed or
    // the property type is not supported at all

    if (!PyErr_Occurred()) {
        // no error so far => unknown property type
        std::ostringstream errStr;
        errStr << "setPropertyValue() Property '" << property->getFullyQualifiedID() << "'";
        errStr << " has unsupported type: '" << property->getClassName() << "'";
        PyErr_SetString(PyExc_ValueError, errStr.str().c_str());
    }

    return 0; //< indicates failure
}

static PyObject* voreen_getPropertyValue(PyObject* /*self*/, PyObject* args) {

    // Parse passed arguments: processor name, property ID
    char *processorName, *propertyID;
    PyArg_ParseTuple(args, "ss:getPropertyValue", &processorName, &propertyID);
    if (PyErr_Occurred())
        return 0;

    // fetch property
    Property* property = getProperty(std::string(processorName), std::string(propertyID), "getPropertyValue");
    if (!property)
        return 0;

    // determine property type and return value, if type compatible
    PyObject* result = (PyObject*)-1; //< to determine whether Py_BuildValue has been executed
    if (BoolProperty* typedProp = dynamic_cast<BoolProperty*>(property)) {
        result = Py_BuildValue("b", typedProp->get());
    }
    else if (CameraProperty* typedProp = dynamic_cast<CameraProperty*>(property)) {
        tgt::vec3 position = typedProp->get().getPosition();
        tgt::vec3 focus = typedProp->get().getFocus();
        tgt::vec3 upVector = typedProp->get().getUpVector();
        result = Py_BuildValue("([fff][fff][fff])",
                               position.x, position.y, position.z,
                               focus.x, focus.y, focus.z,
                               upVector.x, upVector.y, upVector.z
        );
    }
    else if (FloatProperty* typedProp = dynamic_cast<FloatProperty*>(property)) {
        result = Py_BuildValue("f", typedProp->get());
    }
    else if (IntProperty* typedProp = dynamic_cast<IntProperty*>(property)) {
        result = Py_BuildValue("i", typedProp->get());
    }
    else if (PlaneProperty* typedProp = dynamic_cast<PlaneProperty*>(property)) {
        tgt::vec3 normal = typedProp->get().n;
        float d = typedProp->get().d;
        result = Py_BuildValue("([fff]f)",
                               &normal.x, &normal.y, &normal.z,
                               d
        );
    }
    //else if (ShaderProperty* typedProp = dynamic_cast<ShaderProperty*>(property)) {
        // Not supported due to dynamic document count
    //}
    else if (StringProperty* typedProp = dynamic_cast<StringProperty*>(property)) {
        result = Py_BuildValue("s", typedProp->get().c_str());
    }
    else if (OptionPropertyBase* typedProp = dynamic_cast<OptionPropertyBase*>(property)) {
        result = Py_BuildValue("s", typedProp->get().c_str());
    }
    else if (IntBoundingBoxProperty* typedProp = dynamic_cast<IntBoundingBoxProperty*>(property)) {
        tgt::IntBounds value = typedProp->get();
        result = Py_BuildValue("[[iii][iii]]",
                               value.getLLF().x, value.getLLF().y, value.getLLF().z,
                               value.getURB().x, value.getURB().y, value.getURB().z
        );
    }
    else if (FloatBoundingBoxProperty* typedProp = dynamic_cast<FloatBoundingBoxProperty*>(property)) {
        tgt::Bounds value = typedProp->get();
        result = Py_BuildValue("[[fff][fff]]",
                               value.getLLF().x, value.getLLF().y, value.getLLF().z,
                               value.getURB().x, value.getURB().y, value.getURB().z
        );
    }
    else if (IntVec2Property* typedProp = dynamic_cast<IntVec2Property*>(property)) {
        tgt::ivec2 value = typedProp->get();
        result = Py_BuildValue("[ii]", value.x, value.y);
    }
    else if (IntVec3Property* typedProp = dynamic_cast<IntVec3Property*>(property)) {
        tgt::ivec3 value = typedProp->get();
        result = Py_BuildValue("[iii]", value.x, value.y, value.z);
    }
    else if (IntVec4Property* typedProp = dynamic_cast<IntVec4Property*>(property)) {
        tgt::ivec4 value = typedProp->get();
        result = Py_BuildValue("[iiii]", value.x, value.y, value.z, value.w);
    }
    else if (FloatVec2Property* typedProp = dynamic_cast<FloatVec2Property*>(property)) {
        tgt::vec2 value = typedProp->get();
        result = Py_BuildValue("[ff]", value.x, value.y);
    }
    else if (FloatVec3Property* typedProp = dynamic_cast<FloatVec3Property*>(property)) {
        tgt::vec3 value = typedProp->get();
        result = Py_BuildValue("[fff]", value.x, value.y, value.z);
    }
    else if (FloatVec4Property* typedProp = dynamic_cast<FloatVec4Property*>(property)) {
        tgt::vec4 value = typedProp->get();
        result = Py_BuildValue("[ffff]", value.x, value.y, value.z, value.w);
    }
    else if (FloatMat2Property* typedProp = dynamic_cast<FloatMat2Property*>(property)) {
        tgt::mat2 value = typedProp->get();
        result = Py_BuildValue("[[ff][ff]]",
                    value[0][0], value[0][1],
                    value[1][0], value[1][1]);
    }
    else if (FloatMat3Property* typedProp = dynamic_cast<FloatMat3Property*>(property)) {
        tgt::mat3 value = typedProp->get();
        result = Py_BuildValue("[[fff][fff][fff]]",
                    value[0][0], value[0][1], value[0][2],
                    value[1][0], value[1][1], value[1][2],
                    value[2][0], value[2][1], value[2][2]);
    }
    else if (FloatMat4Property* typedProp = dynamic_cast<FloatMat4Property*>(property)) {
        tgt::mat4 value = typedProp->get();
        result = Py_BuildValue("[[ffff][ffff][ffff][ffff]]",
                    value[0][0], value[0][1], value[0][2], value[0][3],
                    value[1][0], value[1][1], value[1][2], value[1][3],
                    value[2][0], value[2][1], value[2][2], value[2][3],
                    value[3][0], value[3][1], value[3][2], value[3][3]);
    }

    // if result is still -1, Py_BuildValue has not been executed
    if (result == (PyObject*)-1) {
        std::ostringstream errStr;
        errStr << "getPropertyValue() Property '" << property->getFullyQualifiedID() << "'";
        errStr << " has unsupported type: '" << property->getTypeDescription() << "'";
        PyErr_SetString(PyExc_ValueError, errStr.str().c_str());
        return 0;
    }

    return result;
}

static PyObject* voreen_setPropertyMinValue(PyObject* /*self*/, PyObject* args) {

    // check length of tuple
    if (PyTuple_Size(args) != 3) {
        std::ostringstream errStr;
        errStr << "setPropertyMinValue() takes exactly 3 arguments: processor name, property id, value";
        errStr << " (" << PyTuple_Size(args) << " given)";
        PyErr_SetString(PyExc_TypeError, errStr.str().c_str());
        return 0;
    }

    // check parameter 1 and 2, if they are strings
    if (!PyUnicode_Check(PyTuple_GetItem(args, 0)) || !PyUnicode_Check(PyTuple_GetItem(args, 1))) {
        PyErr_SetString(PyExc_TypeError, "setPropertyValue() arguments 1 and 2 must be strings");
        return 0;
    }

    // read processor name and property id
    std::string processorName_str = PyUnicodeAsString(PyTuple_GetItem(args, 0));
    std::string propertyID_str = PyUnicodeAsString(PyTuple_GetItem(args, 1));
    if (processorName_str.empty() || propertyID_str.empty()) {
        PyErr_SetString(PyExc_TypeError, "setPropertyMinValue() arguments 1 and 2 must be strings");
        return 0;
    }

    char* processorName;
    char* propertyID;

    // fetch property
    Property* property = getProperty(processorName_str, propertyID_str, "setPropertyMinValue");
    if (!property)
        return 0;

    // determine property type, convert and assign value
    if (IntProperty* typedProp = dynamic_cast<IntProperty*>(property)) {
        int value;
        if (!PyArg_ParseTuple(args, "ssi:setPropertyMinValue", &processorName, &propertyID, &value))
            return 0;
        typedProp->setMinValue(value);
        Py_RETURN_NONE;
    }
    else if (FloatProperty* typedProp = dynamic_cast<FloatProperty*>(property)) {
        float value;
        if (!PyArg_ParseTuple(args, "ssf:setPropertyMinValue", &processorName, &propertyID, &value))
            return 0;
        typedProp->setMinValue(value);
        Py_RETURN_NONE;
    }
    else if (IntBoundingBoxProperty* typedProp = dynamic_cast<IntBoundingBoxProperty*>(property)) {
        tgt::ivec3 value;
        if (!PyArg_ParseTuple(args, "ss(iii):setPropertyMinValue", &processorName, &propertyID,
                              &value.x, &value.y, &value.z))
            return 0;
        typedProp->setMinValue(value);
        Py_RETURN_NONE;
    }
    else if (FloatBoundingBoxProperty* typedProp = dynamic_cast<FloatBoundingBoxProperty*>(property)) {
        tgt::vec3 value;
        if (!PyArg_ParseTuple(args, "ss(fff):setPropertyMinValue", &processorName, &propertyID,
                              &value.x, &value.y, &value.z))
            return 0;
        typedProp->setMinValue(value);
        Py_RETURN_NONE;
    }
    else if (IntVec2Property* typedProp = dynamic_cast<IntVec2Property*>(property)) {
        tgt::ivec2 value;
        if (!PyArg_ParseTuple(args, "ss(ii):setPropertyMinValue", &processorName, &propertyID,
                &value.x, &value.y))
            return 0;
        typedProp->setMinValue(value);
        Py_RETURN_NONE;
    }
    else if (IntVec3Property* typedProp = dynamic_cast<IntVec3Property*>(property)) {
        tgt::ivec3 value;
        if (!PyArg_ParseTuple(args, "ss(iii):setPropertyMinValue", &processorName, &propertyID,
                &value.x, &value.y, &value.z))
            return 0;
        typedProp->setMinValue(value);
        Py_RETURN_NONE;
    }
    else if (IntVec4Property* typedProp = dynamic_cast<IntVec4Property*>(property)) {
        tgt::ivec4 value;
        if (!PyArg_ParseTuple(args, "ss(iiii):setPropertyMinValue", &processorName, &propertyID,
                &value.x, &value.y, &value.z, &value.w))
            return 0;
        typedProp->setMinValue(value);
        Py_RETURN_NONE;
    }
    else if (FloatVec2Property* typedProp = dynamic_cast<FloatVec2Property*>(property)) {
        tgt::vec2 value;
        if (!PyArg_ParseTuple(args, "ss(ff):setPropertyMinValue", &processorName, &propertyID,
                &value.x, &value.y))
            return 0;
        typedProp->setMinValue(value);
        Py_RETURN_NONE;
    }
    else if (FloatVec3Property* typedProp = dynamic_cast<FloatVec3Property*>(property)) {
        tgt::vec3 value;
        if (!PyArg_ParseTuple(args, "ss(fff):setPropertyMinValue", &processorName, &propertyID,
                &value.x, &value.y, &value.z))
            return 0;
        typedProp->setMinValue(value);
        Py_RETURN_NONE;
    }
    else if (FloatVec4Property* typedProp = dynamic_cast<FloatVec4Property*>(property)) {
        tgt::vec4 value;
        if (!PyArg_ParseTuple(args, "ss(ffff):setPropertyMinValue", &processorName, &propertyID,
                &value.x, &value.y, &value.z, &value.w))
            return 0;
        typedProp->setMinValue(value);
        Py_RETURN_NONE;
    }
    else if (FloatMat2Property* typedProp = dynamic_cast<FloatMat2Property*>(property)) {
        tgt::vec2 vec0, vec1;
        char *processorName, *propertyID;
        if (!PyArg_ParseTuple(args, "ss((ff)(ff)):setPropertyMinValue",
                &processorName, &propertyID,
                &vec0.x, &vec0.y, &vec1.x, &vec1.y))
            return 0;

        tgt::mat2 value(vec0, vec1);
        typedProp->setMinValue(value);
        Py_RETURN_NONE;
    }
    else if (FloatMat3Property* typedProp = dynamic_cast<FloatMat3Property*>(property)) {
        tgt::vec3 vec0, vec1, vec2;
        char *processorName, *propertyID;
        if (!PyArg_ParseTuple(args, "ss((fff)(fff)(fff)):setPropertyMinValue",
                &processorName, &propertyID,
                &vec0.x, &vec0.y, &vec0.z, &vec1.x, &vec1.y, &vec1.z,
                &vec2.x, &vec2.y, &vec2.z))
            return 0;

        tgt::mat3 value(vec0, vec1, vec2);
        typedProp->setMinValue(value);
        Py_RETURN_NONE;
    }
    else if (FloatMat4Property* typedProp = dynamic_cast<FloatMat4Property*>(property)) {
        tgt::vec4 vec0, vec1, vec2, vec3;
        char *processorName, *propertyID;
        if (!PyArg_ParseTuple(args, "ss((ffff)(ffff)(ffff)(ffff)):setPropertyMinValue",
                &processorName, &propertyID,
                &vec0.x, &vec0.y, &vec0.z, &vec0.w, &vec1.x, &vec1.y, &vec1.z, &vec1.w,
                &vec2.x, &vec2.y, &vec2.z, &vec2.w, &vec3.x, &vec3.y, &vec3.z, &vec3.w))
            return 0;

        tgt::mat4 value(vec0, vec1, vec2, vec3);
        typedProp->setMinValue(value);
        Py_RETURN_NONE;
    }


    // we only get here, if property value assignment has failed or
    // the property type is not supported at all

    if (!PyErr_Occurred()) {
        // no error so far => unknown property type
        std::ostringstream errStr;
        errStr << "setPropertyMinValue() Property '" << property->getFullyQualifiedID() << "'";
        errStr << " has unsupported type: '" << property->getTypeDescription() << "'";
        PyErr_SetString(PyExc_ValueError, errStr.str().c_str());
    }

    return 0; //< indicates failure
}

static PyObject* voreen_setPropertyMaxValue(PyObject* /*self*/, PyObject* args) {

    // check length of tuple
    if (PyTuple_Size(args) != 3) {
        std::ostringstream errStr;
        errStr << "setPropertyMaxValue() takes exactly 3 arguments: processor name, property id, value";
        errStr << " (" << PyTuple_Size(args) << " given)";
        PyErr_SetString(PyExc_TypeError, errStr.str().c_str());
        return 0;
    }

    // check parameter 1 and 2, if they are strings
    if (!PyUnicode_Check(PyTuple_GetItem(args, 0)) || !PyUnicode_Check(PyTuple_GetItem(args, 1))) {
        PyErr_SetString(PyExc_TypeError, "setPropertyValue() arguments 1 and 2 must be strings");
        return 0;
    }

    // read processor name and property id
    std::string processorName_str = PyUnicodeAsString(PyTuple_GetItem(args, 0));
    std::string propertyID_str = PyUnicodeAsString(PyTuple_GetItem(args, 1));
    if (processorName_str.empty() || propertyID_str.empty()) {
        PyErr_SetString(PyExc_TypeError, "setPropertyMaxValue() arguments 1 and 2 must be strings");
        return 0;
    }

    char* processorName;
    char* propertyID;

    // fetch property
    Property* property = getProperty(processorName_str, propertyID_str, "setPropertyMaxValue");
    if (!property)
        return 0;

    // determine property type, convert and assign value
    if (IntProperty* typedProp = dynamic_cast<IntProperty*>(property)) {
        int value;
        if (!PyArg_ParseTuple(args, "ssi:setPropertyMaxValue", &processorName, &propertyID, &value))
            return 0;
        typedProp->setMaxValue(value);
        Py_RETURN_NONE;
    }
    else if (FloatProperty* typedProp = dynamic_cast<FloatProperty*>(property)) {
        float value;
        if (!PyArg_ParseTuple(args, "ssf:setPropertyMaxValue", &processorName, &propertyID, &value))
            return 0;
        typedProp->setMaxValue(value);
        Py_RETURN_NONE;
    }
    else if (IntBoundingBoxProperty* typedProp = dynamic_cast<IntBoundingBoxProperty*>(property)) {
        tgt::ivec3 value;
        if (!PyArg_ParseTuple(args, "ss(iii):setPropertyMaxValue", &processorName, &propertyID,
                              &value.x, &value.y, &value.z))
            return 0;
        typedProp->setMaxValue(value);
        Py_RETURN_NONE;
    }
    else if (FloatBoundingBoxProperty* typedProp = dynamic_cast<FloatBoundingBoxProperty*>(property)) {
        tgt::vec3 value;
        if (!PyArg_ParseTuple(args, "ss(fff):setPropertyMaxValue", &processorName, &propertyID,
                              &value.x, &value.y, &value.z))
            return 0;
        typedProp->setMaxValue(value);
        Py_RETURN_NONE;
    }
    else if (IntVec2Property* typedProp = dynamic_cast<IntVec2Property*>(property)) {
        tgt::ivec2 value;
        if (!PyArg_ParseTuple(args, "ss(ii):setPropertyMaxValue", &processorName, &propertyID,
                &value.x, &value.y))
            return 0;
        typedProp->setMaxValue(value);
        Py_RETURN_NONE;
    }
    else if (IntVec3Property* typedProp = dynamic_cast<IntVec3Property*>(property)) {
        tgt::ivec3 value;
        if (!PyArg_ParseTuple(args, "ss(iii):setPropertyMaxValue", &processorName, &propertyID,
                &value.x, &value.y, &value.z))
            return 0;
        typedProp->setMaxValue(value);
        Py_RETURN_NONE;
    }
    else if (IntVec4Property* typedProp = dynamic_cast<IntVec4Property*>(property)) {
        tgt::ivec4 value;
        if (!PyArg_ParseTuple(args, "ss(iiii):setPropertyMaxValue", &processorName, &propertyID,
                &value.x, &value.y, &value.z, &value.w))
            return 0;
        typedProp->setMaxValue(value);
        Py_RETURN_NONE;
    }
    else if (FloatVec2Property* typedProp = dynamic_cast<FloatVec2Property*>(property)) {
        tgt::vec2 value;
        if (!PyArg_ParseTuple(args, "ss(ff):setPropertyMaxValue", &processorName, &propertyID,
                &value.x, &value.y))
            return 0;
        typedProp->setMaxValue(value);
        Py_RETURN_NONE;
    }
    else if (FloatVec3Property* typedProp = dynamic_cast<FloatVec3Property*>(property)) {
        tgt::vec3 value;
        if (!PyArg_ParseTuple(args, "ss(fff):setPropertyMaxValue", &processorName, &propertyID,
                &value.x, &value.y, &value.z))
            return 0;
        typedProp->setMaxValue(value);
        Py_RETURN_NONE;
    }
    else if (FloatVec4Property* typedProp = dynamic_cast<FloatVec4Property*>(property)) {
        tgt::vec4 value;
        if (!PyArg_ParseTuple(args, "ss(ffff):setPropertyMaxValue", &processorName, &propertyID,
                &value.x, &value.y, &value.z, &value.w))
            return 0;
        typedProp->setMaxValue(value);
        Py_RETURN_NONE;
    }
    else if (FloatMat2Property* typedProp = dynamic_cast<FloatMat2Property*>(property)) {
        tgt::vec2 vec0, vec1;
        char *processorName, *propertyID;
        if (!PyArg_ParseTuple(args, "ss((ff)(ff)):setPropertyMaxValue",
                &processorName, &propertyID,
                &vec0.x, &vec0.y, &vec1.x, &vec1.y))
            return 0;

        tgt::mat2 value(vec0, vec1);
        typedProp->setMaxValue(value);
        Py_RETURN_NONE;
    }
    else if (FloatMat3Property* typedProp = dynamic_cast<FloatMat3Property*>(property)) {
        tgt::vec3 vec0, vec1, vec2;
        char *processorName, *propertyID;
        if (!PyArg_ParseTuple(args, "ss((fff)(fff)(fff)):setPropertyMaxValue",
                &processorName, &propertyID,
                &vec0.x, &vec0.y, &vec0.z, &vec1.x, &vec1.y, &vec1.z,
                &vec2.x, &vec2.y, &vec2.z))
            return 0;

        tgt::mat3 value(vec0, vec1, vec2);
        typedProp->setMaxValue(value);
        Py_RETURN_NONE;
    }
    else if (FloatMat4Property* typedProp = dynamic_cast<FloatMat4Property*>(property)) {
        tgt::vec4 vec0, vec1, vec2, vec3;
        char *processorName, *propertyID;
        if (!PyArg_ParseTuple(args, "ss((ffff)(ffff)(ffff)(ffff)):setPropertyMaxValue",
                &processorName, &propertyID,
                &vec0.x, &vec0.y, &vec0.z, &vec0.w, &vec1.x, &vec1.y, &vec1.z, &vec1.w,
                &vec2.x, &vec2.y, &vec2.z, &vec2.w, &vec3.x, &vec3.y, &vec3.z, &vec3.w))
            return 0;

        tgt::mat4 value(vec0, vec1, vec2, vec3);
        typedProp->setMaxValue(value);
        Py_RETURN_NONE;
    }


    // we only get here, if property value assignment has failed or
    // the property type is not supported at all

    if (!PyErr_Occurred()) {
        // no error so far => unknown property type
        std::ostringstream errStr;
        errStr << "setPropertyMaxValue() Property '" << property->getFullyQualifiedID() << "'";
        errStr << " has unsupported type: '" << property->getTypeDescription() << "'";
        PyErr_SetString(PyExc_ValueError, errStr.str().c_str());
    }

    return 0; //< indicates failure
}

static PyObject* voreen_getPropertyMinValue(PyObject* /*self*/, PyObject* args) {

    // Parse passed arguments: processor name, property ID
    char *processorName, *propertyID;
    PyArg_ParseTuple(args, "ss:getPropertyMinValue", &processorName, &propertyID);
    if (PyErr_Occurred())
        return 0;

    // fetch property
    Property* property = getProperty(std::string(processorName), std::string(propertyID), "getPropertyMinValue");
    if (!property)
        return 0;

    // determine property type and return value, if type compatible
    PyObject* result = (PyObject*)-1; //< to determine whether Py_BuildValue has been executed
    if (IntProperty* typedProp = dynamic_cast<IntProperty*>(property))
        result = Py_BuildValue("i", typedProp->getMinValue());
    else if (FloatProperty* typedProp = dynamic_cast<FloatProperty*>(property))
        result = Py_BuildValue("f", typedProp->getMinValue());
    else if (IntVec2Property* typedProp = dynamic_cast<IntVec2Property*>(property)) {
        tgt::ivec2 value = typedProp->getMinValue();
        result = Py_BuildValue("[ii]", value.x, value.y);
    }
    else if (IntVec3Property* typedProp = dynamic_cast<IntVec3Property*>(property)) {
        tgt::ivec3 value = typedProp->getMinValue();
        result = Py_BuildValue("[iii]", value.x, value.y, value.z);
    }
    else if (IntVec4Property* typedProp = dynamic_cast<IntVec4Property*>(property)) {
        tgt::ivec4 value = typedProp->getMinValue();
        result = Py_BuildValue("[iiii]", value.x, value.y, value.z, value.w);
    }
    else if (FloatVec2Property* typedProp = dynamic_cast<FloatVec2Property*>(property)) {
        tgt::vec2 value = typedProp->getMinValue();
        result = Py_BuildValue("[ff]", value.x, value.y);
    }
    else if (FloatVec3Property* typedProp = dynamic_cast<FloatVec3Property*>(property)) {
        tgt::vec3 value = typedProp->getMinValue();
        result = Py_BuildValue("[fff]", value.x, value.y, value.z);
    }
    else if (FloatVec4Property* typedProp = dynamic_cast<FloatVec4Property*>(property)) {
        tgt::vec4 value = typedProp->getMinValue();
        result = Py_BuildValue("[ffff]", value.x, value.y, value.z, value.w);
    }
    else if (FloatMat2Property* typedProp = dynamic_cast<FloatMat2Property*>(property)) {
        tgt::mat2 value = typedProp->getMinValue();
        result = Py_BuildValue("[[ff][ff]]",
                    value[0][0], value[0][1],
                    value[1][0], value[1][1]);
    }
    else if (FloatMat3Property* typedProp = dynamic_cast<FloatMat3Property*>(property)) {
        tgt::mat3 value = typedProp->getMinValue();
        result = Py_BuildValue("[[fff][fff][fff]]",
                    value[0][0], value[0][1], value[0][2],
                    value[1][0], value[1][1], value[1][2],
                    value[2][0], value[2][1], value[2][2]);
    }
    else if (FloatMat4Property* typedProp = dynamic_cast<FloatMat4Property*>(property)) {
        tgt::mat4 value = typedProp->getMinValue();
        result = Py_BuildValue("[[ffff][ffff][ffff][ffff]]",
                    value[0][0], value[0][1], value[0][2], value[0][3],
                    value[1][0], value[1][1], value[1][2], value[1][3],
                    value[2][0], value[2][1], value[2][2], value[2][3],
                    value[3][0], value[3][1], value[3][2], value[3][3]);
    }

    // if result is still -1, Py_BuildValue has not been executed
    if (result == (PyObject*)-1) {
        std::ostringstream errStr;
        errStr << "getPropertyMinValue() Property '" << property->getFullyQualifiedID() << "'";
        errStr << " has unsupported type: '" << property->getTypeDescription() << "'";
        PyErr_SetString(PyExc_ValueError, errStr.str().c_str());
        return 0;
    }

    return result;
}

static PyObject* voreen_getPropertyMaxValue(PyObject* /*self*/, PyObject* args) {

    // Parse passed arguments: processor name, property ID
    char *processorName, *propertyID;
    PyArg_ParseTuple(args, "ss:getPropertyMaxValue", &processorName, &propertyID);
    if (PyErr_Occurred())
        return 0;

    // fetch property
    Property* property = getProperty(std::string(processorName), std::string(propertyID), "getPropertyMaxValue");
    if (!property)
        return 0;

    // determine property type and return value, if type compatible
    PyObject* result = (PyObject*)-1; //< to determine whether Py_BuildValue has been executed
    if (IntProperty* typedProp = dynamic_cast<IntProperty*>(property))
        result = Py_BuildValue("i", typedProp->getMaxValue());
    else if (FloatProperty* typedProp = dynamic_cast<FloatProperty*>(property))
        result = Py_BuildValue("f", typedProp->getMaxValue());
    else if (IntVec2Property* typedProp = dynamic_cast<IntVec2Property*>(property)) {
        tgt::ivec2 value = typedProp->getMaxValue();
        result = Py_BuildValue("[ii]", value.x, value.y);
    }
    else if (IntVec3Property* typedProp = dynamic_cast<IntVec3Property*>(property)) {
        tgt::ivec3 value = typedProp->getMaxValue();
        result = Py_BuildValue("[iii]", value.x, value.y, value.z);
    }
    else if (IntVec4Property* typedProp = dynamic_cast<IntVec4Property*>(property)) {
        tgt::ivec4 value = typedProp->getMaxValue();
        result = Py_BuildValue("[iiii]", value.x, value.y, value.z, value.w);
    }
    else if (FloatVec2Property* typedProp = dynamic_cast<FloatVec2Property*>(property)) {
        tgt::vec2 value = typedProp->getMaxValue();
        result = Py_BuildValue("[ff]", value.x, value.y);
    }
    else if (FloatVec3Property* typedProp = dynamic_cast<FloatVec3Property*>(property)) {
        tgt::vec3 value = typedProp->getMaxValue();
        result = Py_BuildValue("[fff]", value.x, value.y, value.z);
    }
    else if (FloatVec4Property* typedProp = dynamic_cast<FloatVec4Property*>(property)) {
        tgt::vec4 value = typedProp->getMaxValue();
        result = Py_BuildValue("[ffff]", value.x, value.y, value.z, value.w);
    }
    else if (FloatMat2Property* typedProp = dynamic_cast<FloatMat2Property*>(property)) {
        tgt::mat2 value = typedProp->getMaxValue();
        result = Py_BuildValue("[[ff][ff]]",
                    value[0][0], value[0][1],
                    value[1][0], value[1][1]);
    }
    else if (FloatMat3Property* typedProp = dynamic_cast<FloatMat3Property*>(property)) {
        tgt::mat3 value = typedProp->getMaxValue();
        result = Py_BuildValue("[[fff][fff][fff]]",
                    value[0][0], value[0][1], value[0][2],
                    value[1][0], value[1][1], value[1][2],
                    value[2][0], value[2][1], value[2][2]);
    }
    else if (FloatMat4Property* typedProp = dynamic_cast<FloatMat4Property*>(property)) {
        tgt::mat4 value = typedProp->getMaxValue();
        result = Py_BuildValue("[[ffff][ffff][ffff][ffff]]",
                    value[0][0], value[0][1], value[0][2], value[0][3],
                    value[1][0], value[1][1], value[1][2], value[1][3],
                    value[2][0], value[2][1], value[2][2], value[2][3],
                    value[3][0], value[3][1], value[3][2], value[3][3]);
    }

    // if result is still -1, Py_BuildValue has not been executed
    if (result == (PyObject*)-1) {
        std::ostringstream errStr;
        errStr << "getPropertyMaxValue() Property '" << property->getFullyQualifiedID() << "'";
        errStr << " has unsupported type: '" << property->getTypeDescription() << "'";
        PyErr_SetString(PyExc_ValueError, errStr.str().c_str());
        return 0;
    }

    return result;
}

static PyObject* voreen_setPortData(PyObject* /*self*/, PyObject* args) {

    // check length of tuple
    if (PyTuple_Size(args) != 3) {
        std::ostringstream errStr;
        errStr << "setPortData() takes exactly 3 arguments: processor name, port id, data";
        errStr << " (" << PyTuple_Size(args) << " given)";
        PyErr_SetString(PyExc_TypeError, errStr.str().c_str());
        return 0;
    }

    // check parameter 1 and 2, if they are strings
    if (!PyUnicode_Check(PyTuple_GetItem(args, 0)) || !PyUnicode_Check(PyTuple_GetItem(args, 1))) {
        PyErr_SetString(PyExc_TypeError, "setPortData() arguments 1 and 2 must be strings");
        return 0;
    }

    // read processor name and port id
    std::string processorName_str = PyUnicodeAsString(PyTuple_GetItem(args, 0));
    std::string portID_str = PyUnicodeAsString(PyTuple_GetItem(args, 1));
    if (processorName_str.empty() || portID_str.empty()) {
        PyErr_SetString(PyExc_TypeError, "setPortData() arguments 1 and 2 must be strings");
        return 0;
    }

    char* processorName;
    char* portID;

    // fetch port
    Port* port = getPort(processorName_str, portID_str, "setPortData");
    if (!port)
        return 0;

    if (port->isInport()) {
        PyErr_SetString(PyExc_TypeError, "setPortData() must be called for outgoing ports");
        return 0;
    }

    // Delete old data first.
    port->clear();

    // determine port type, convert and assign data
    if (VolumePort* typedPort = dynamic_cast<VolumePort*>(port)) {
        PyObject* object = NULL;
        if (!PyArg_ParseTuple(args, "ssO!:setPortData", &processorName, &portID, &VolumeObjectType, &object))
            return 0;
        tgtAssert(object, "parsing VolumeObject failed");

        const VolumeObject& volumeObject = *((VolumeObject*) object);

        unsigned int numVoxels = volumeObject.dimX * volumeObject.dimY * volumeObject.dimZ;
        Py_ssize_t size = PyList_Size(volumeObject.data);
        if(size != numVoxels) {
            std::string error = "Volume data has invalid size '" + std::to_string(size) + "', must be '" + std::to_string(numVoxels) + "') according to dimensions";
            PyErr_SetString(PyExc_ValueError, error.c_str());
            return 0;
        }

        std::string format = PyUnicodeAsString(volumeObject.format);
        std::unique_ptr<VolumeRAM> volume;
        try {
            volume.reset(VolumeFactory().create(format, tgt::svec3(volumeObject.dimX, volumeObject.dimY, volumeObject.dimZ)));
        } catch(const std::bad_alloc& error) {
            PyErr_SetString(PyExc_ValueError, "Could not allocate memory for volume");
            return 0;
        }
        if(!volume) {
            std::string error = "Volume of format '" + format + "' could not be created.";
            PyErr_SetString(PyExc_ValueError, error.c_str());
            return 0;
        }

        // Set voxel values.
        bool error = false;
        for(size_t i=0; i<volume->getNumVoxels() && !error; i++) {
            PyObject* p = PyList_GetItem(volumeObject.data, i);

            tgt::vec4 value = tgt::vec4::zero;
            switch(volume->getNumChannels()) {
            case 1:
                value.x = static_cast<float>(PyFloat_AsDouble(p));
                error |= PyErr_Occurred() != nullptr;
                break;
            case 2:
                error |= !PyArg_ParseTuple(p, "ff", &value.x, &value.y);
                break;
            case 3:
                error |= !PyArg_ParseTuple(p, "fff", &value.x, &value.y, &value.z);
                break;
            case 4:
                error |= !PyArg_ParseTuple(p, "ffff", &value.x, &value.y, &value.z, &value.w);
                break;
            default:
                tgtAssert(false, "unsupported channel count");
            }

            for(size_t channel = 0; channel < volume->getNumChannels(); channel++) {
                volume->setVoxelNormalized(value[channel], i, channel);
            }
        }

        if(error) {
            std::string message = "Volume contains invalid values.";
            PyErr_SetString(PyExc_ValueError, message.c_str());
            return 0;
        }

        // Set meta data.
        tgt::vec3 spacing(volumeObject.spacingX, volumeObject.spacingY, volumeObject.spacingZ);
        tgt::vec3 offset(volumeObject.offsetX, volumeObject.offsetY, volumeObject.offsetZ);
        RealWorldMapping rwm(volumeObject.rwmScale, volumeObject.rwmOffset, "");

        Volume* data = new Volume(volume.release(), spacing, offset);
        data->setRealWorldMapping(rwm);
        typedPort->setData(data, true);

        Py_RETURN_NONE;
    }
    else if (RenderPort* typedPort = dynamic_cast<RenderPort*>(port)) {
        PyObject* object = NULL;
        if (!PyArg_ParseTuple(args, "ssO!:setPortData", &processorName, &portID, &RenderTargetObjectType, &object))
            return 0;
        tgtAssert(object, "parsing RenderTargetObject failed");

        const RenderTargetObject& renderTargetObject = *((RenderTargetObject*) object);

        RenderTarget* renderTarget = typedPort->getRenderTarget();
        if(!renderTarget) {
            std::string error = "Port has no valid RenderTarget";
            PyErr_SetString(PyExc_ValueError, error.c_str());
            return 0;
        }

        unsigned int numPixels = renderTargetObject.width * renderTargetObject.height;
        Py_ssize_t colorTextureSize = PyList_Size(renderTargetObject.colorTexture);
        if(colorTextureSize != numPixels) {
            std::string error = "Color texture data has invalid size '" + std::to_string(colorTextureSize) + "', must be '" + std::to_string(numPixels) + "') according to dimensions";
            PyErr_SetString(PyExc_ValueError, error.c_str());
            return 0;
        }

        Py_ssize_t depthTextureSize = PyList_Size(renderTargetObject.depthTexture);
        if(depthTextureSize != numPixels) {
            std::string error = "Depth texture data has invalid size '" + std::to_string(depthTextureSize) + "', must be '" + std::to_string(numPixels) + "') according to dimensions";
            PyErr_SetString(PyExc_ValueError, error.c_str());
            return 0;
        }

        auto textureUpload = [] (tgt::Texture* texture, PyObject* target) {
            bool error = false;

            // Ensure cpu texture data is available.
            if(!texture->getCpuTextureData() && !texture->alloc(true))
                return false;

            for(int y = 0; y < texture->getDimensions().y; y++) {
                for(int x = 0; x < texture->getDimensions().x; x++) {

                    int index = y * texture->getDimensions().x + x;

                    PyObject *p = PyList_GetItem(target, index);
                    tgt::Color value;
                    switch (texture->getNumChannels()) {
                        case 1:
                            value.x = static_cast<float>(PyFloat_AsDouble(p));
                            error |= PyErr_Occurred() != nullptr;
                            break;
                        case 2:
                            error |= !PyArg_ParseTuple(p, "ff", &value.x, &value.y);
                            break;
                        case 3:
                            error |= !PyArg_ParseTuple(p, "fff", &value.x, &value.y, &value.z);
                            break;
                        case 4:
                            error |= !PyArg_ParseTuple(p, "ffff", &value.x, &value.y, &value.z, &value.w);
                            break;
                    }

                    if(error)
                        return false;

                    texture->texelFromFloat(value, x, y);
                }
            }

            // Upload data.
            texture->uploadTexture();

            return true;
        };

        tgt::GLConditionalContextStateGuard guard(tgt::isInitedGL());
        typedPort->activateTarget();
        typedPort->clearTarget();

        bool error = !textureUpload(renderTarget->getColorTexture(), renderTargetObject.colorTexture);
        if (!error) {
            error =  !textureUpload(renderTarget->getDepthTexture(), renderTargetObject.depthTexture);
        }

        renderTarget->deactivateTarget();

        if(error || PyErr_Occurred()) {
            std::string error = "Pixel data contains invalid values.";
            PyErr_SetString(PyExc_ValueError, error.c_str());
            return 0;
        }

        Py_RETURN_NONE;
    }

    // we only get here, if port data assignment has failed or
    // the port type is not supported at all

    if (!PyErr_Occurred()) {
        // no error so far => unknown port type
        std::ostringstream errStr;
        errStr << "setPortData() Port '" << port->getQualifiedName() << "'";
        errStr << " has unsupported type: '" << port->getClassName() << "'";
        PyErr_SetString(PyExc_ValueError, errStr.str().c_str());
    }

    return 0; //< indicates failure
}

static PyObject* voreen_getPortData(PyObject* /*self*/, PyObject* args) {

    // Parse passed arguments: processor name, property ID
    char *processorName, *portID;
    PyArg_ParseTuple(args, "ss:getPortData", &processorName, &portID);
    if (PyErr_Occurred())
        return 0;

    // fetch port
    Port* port = getPort(std::string(processorName), std::string(portID), "getPort");
    if (!port)
        return 0;

    if(!port->hasData()) {
        std::ostringstream errStr;
        errStr << "getPortData() Port '" << port->getQualifiedName() << "' has no data.";
        PyErr_SetString(PyExc_ValueError, errStr.str().c_str());
        return 0;
    }

    PyObject* result = (PyObject*)-1;
    if (VolumePort* typedPort = dynamic_cast<VolumePort*>(port)) {

        const VolumeBase* data = typedPort->getData();
        tgtAssert(data, "volume was null");
        const VolumeRAM* volume = data->getRepresentation<VolumeRAM>();
        if(!volume)
            return 0;

        // Create new volume object.
        VolumeObject* volumeObject = (VolumeObject*) VolumeObject_new(&VolumeObjectType, NULL, NULL);
        if(!volumeObject)
            return 0;

        Py_DECREF(volumeObject->format); // Format will be overwritten.
        volumeObject->format = PyUnicode_FromString(volume->getFormat().c_str());

        for(size_t i=0; i<volume->getNumVoxels(); i++) {
            tgt::vec4 value = tgt::vec4::zero;
            for(size_t channel = 0; channel < volume->getNumChannels(); channel++) {
                value[channel] = volume->getVoxelNormalized(i, channel);
            }

            PyObject* p = NULL;
            switch(volume->getNumChannels()) {
            case 1:
                p = PyFloat_FromDouble(value.x);
                break;
            case 2:
                p = Py_BuildValue("(ff)", value.x, value.y);
                break;
            case 3:
                p = Py_BuildValue("(fff)", value.x, value.y, value.z);
                break;
            case 4:
                p = Py_BuildValue("(ffff)", value.x, value.y, value.z, value.w);
                break;
            default:
                tgtAssert(false, "unsupported channel count");
            }

            tgtAssert(p, "value null");
            PyList_Append(volumeObject->data, p);

            // p is now owned by the list
            Py_XDECREF(p);
        }

        volumeObject->dimX = data->getDimensions().x;
        volumeObject->dimY = data->getDimensions().y;
        volumeObject->dimZ = data->getDimensions().z;

        volumeObject->spacingX = data->getSpacing().x;
        volumeObject->spacingY = data->getSpacing().y;
        volumeObject->spacingZ = data->getSpacing().z;

        volumeObject->offsetX = data->getOffset().x;
        volumeObject->offsetY = data->getOffset().y;
        volumeObject->offsetZ = data->getOffset().z;

        volumeObject->rwmScale  = data->getRealWorldMapping().getScale();
        volumeObject->rwmOffset = data->getRealWorldMapping().getOffset();

        result = (PyObject*) volumeObject;
    }
    else if (RenderPort* typedPort = dynamic_cast<RenderPort*>(port)) {

        const RenderTarget* renderTarget = typedPort->getRenderTarget();
        if(!renderTarget)
            return 0;

        // Create new render target object.
        RenderTargetObject* renderTargetObject = (RenderTargetObject*) RenderTargetObject_new(&RenderTargetObjectType, NULL, NULL);
        if(!renderTargetObject)
            return 0;

        renderTargetObject->internalColorFormat = renderTarget->getColorTexture()->getGLInternalFormat();
        renderTargetObject->internalDepthFormat = renderTarget->getDepthTexture()->getGLInternalFormat();
        renderTargetObject->width  = renderTarget->getSize().x;
        renderTargetObject->height = renderTarget->getSize().y;

        auto textureDownload = [] (const tgt::Texture* texture, PyObject* target) {

            // Ensure main context is active!
            tgt::GLConditionalContextStateGuard guard(tgt::isInitedGL());
            texture->downloadTexture();
            if(!texture->getCpuTextureData())
                return false;

            for(int y = 0; y < texture->getDimensions().y; y++) {
                for(int x = 0; x < texture->getDimensions().x; x++) {
                    tgt::vec4 value = texture->texelAsFloat(x, y);

                    PyObject* p = NULL;
                    switch(texture->getNumChannels()) {
                        case 1:
                            p = PyFloat_FromDouble(value.x);
                            break;
                        case 2:
                            p = Py_BuildValue("(ff)", value.x, value.y);
                            break;
                        case 3:
                            p = Py_BuildValue("(fff)", value.x, value.y, value.z);
                            break;
                        case 4:
                            p = Py_BuildValue("(ffff)", value.x, value.y, value.z, value.w);
                            break;
                        default:
                            tgtAssert(false, "unsupported channel count");
                    }

                    tgtAssert(p, "value null");
                    PyList_Append(target, p);
                    // p is now owned by the list
                    Py_XDECREF(p);
                }
            }

            return true;
        };

        bool error = !textureDownload(renderTarget->getColorTexture(), renderTargetObject->colorTexture);
        if (!error) {
            error =  !textureDownload(renderTarget->getDepthTexture(), renderTargetObject->depthTexture);
        }

        if(error || PyErr_Occurred()) {
            Py_CLEAR(renderTargetObject);
            PyErr_SetString(PyExc_ValueError, "Texture data could not be downloaded");
            return 0;
        }

        result = (PyObject*) renderTargetObject;
    }

    if (result == (PyObject*)-1) {
        std::ostringstream errStr;
        errStr << "getPortData() Port '" << port->getQualifiedName() << "'";
        errStr << " has unsupported type: '" << port->getClassName() << "'";
        PyErr_SetString(PyExc_ValueError, errStr.str().c_str());
        return 0;
    }

    return result;
}

static PyObject* voreen_setCameraPosition(PyObject* /*self*/, PyObject* args) {

    // parse arguments
    char* processorName = 0;
    char* propertyID = 0;
    tgt::vec3 position;
    PyArg_ParseTuple(args, "ss(fff):setCameraPosition", &processorName, &propertyID,
        &position.x, &position.y, &position.z);
    if (PyErr_Occurred())
        return 0;

    // find property and set position
    if (CameraProperty* camProp = getTypedProperty<CameraProperty>(
            std::string(processorName), std::string(propertyID), "Camera", "setCameraPosition")) {
        camProp->setPosition(position);
        camProp->invalidate();
        Py_RETURN_NONE;
    }

    // indicate failure
    return 0;
}

static PyObject* voreen_setCameraFocus(PyObject* /*self*/, PyObject* args) {

    // parse arguments
    char* processorName = 0;
    char* propertyID = 0;
    tgt::vec3 focus;
    PyArg_ParseTuple(args, "ss(fff):setCameraFocus", &processorName, &propertyID,
        &focus.x, &focus.y, &focus.z);
    if (PyErr_Occurred())
        return 0;

    // find property and set position
    if (CameraProperty* camProp = getTypedProperty<CameraProperty>(
            std::string(processorName), std::string(propertyID), "Camera", "setCameraFocus")) {
        camProp->setFocus(focus);
        camProp->invalidate();
        Py_RETURN_NONE;
    }

    // indicate failure
    return 0;
}

static PyObject* voreen_setCameraUpVector(PyObject* /*self*/, PyObject* args) {

    // parse arguments
    char* processorName = 0;
    char* propertyID = 0;
    tgt::vec3 up;
    PyArg_ParseTuple(args, "ss(fff):setCameraUpVector", &processorName, &propertyID,
        &up.x, &up.y, &up.z);
    if (PyErr_Occurred())
        return 0;

    // find property and set position
    if (CameraProperty* camProp = getTypedProperty<CameraProperty>(
            std::string(processorName), std::string(propertyID), "Camera", "setCameraUpVector")) {
        camProp->setUpVector(up);
        camProp->invalidate();
        Py_RETURN_NONE;
    }

    // indicate failure
    return 0;
}

static PyObject* voreen_loadVolume(PyObject* /*self*/, PyObject* args) {

    const char* filename = 0;
    const char* procStr = 0;
    if (!PyArg_ParseTuple(args, "s|s:loadVolume", &filename, &procStr))
        return 0;

    ProcessorNetwork* network = getProcessorNetwork("loadVolume");
    if (!network)
        return 0;

    VolumeSource* volumeSource = 0;
    if (!procStr) {
        // select first volumesource in network
        std::vector<VolumeSource*> sources = network->getProcessorsByType<VolumeSource>();
        if (sources.empty()) {
            PyErr_SetString(PyExc_RuntimeError, "loadVolume() Network does not contain a VolumeSource.");
            return 0;
        }
        volumeSource = sources.front();
    }
    else {
        // retrieve volumesource with given name from network
        volumeSource = getTypedProcessor<VolumeSource>(std::string(procStr), "VolumeSource", "loadVolume");
        if (!volumeSource)
            return 0;
    }
    tgtAssert(volumeSource, "no source proc");

    try {
        volumeSource->loadVolume(std::string(filename));
        Py_RETURN_NONE;
    }
    catch (std::exception &e) {
        PyErr_SetString(PyExc_RuntimeError, std::string("loadVolume() Failed to load data set '" +
            std::string(filename) + "': " + e.what()).c_str());
        return 0;
    }
}

static PyObject* voreen_loadVolumes(PyObject* /*self*/, PyObject* args) {

    const char* filename = 0;
    bool selected = true;
    bool clear = false;
    const char* procStr = 0;
    if (!PyArg_ParseTuple(args, "sbb|s:loadVolumes", &filename, &selected, &clear, &procStr))
        return 0;

    ProcessorNetwork* network = getProcessorNetwork("loadVolumes");
    if (!network)
        return 0;

    VolumeListSource* volumeListSource = 0;
    if (!procStr) {
        // select first volumesource in network
        std::vector<VolumeListSource*> sources = network->getProcessorsByType<VolumeListSource>();
        if (sources.empty()) {
            PyErr_SetString(PyExc_RuntimeError, "loadVolumes() Network does not contain a VolumeListSource.");
            return 0;
        }
        volumeListSource = sources.front();
    }
    else {
        // retrieve volumelistsource with given name from network
        volumeListSource = getTypedProcessor<VolumeListSource>(std::string(procStr), "VolumeListSource", "loadVolumes");
        if (!volumeListSource)
            return 0;
    }
    tgtAssert(volumeListSource, "no source proc");

    try {
        volumeListSource->loadVolumes(std::string(filename), selected, clear);
        Py_RETURN_NONE;
    }
    catch (std::exception &e) {
        PyErr_SetString(PyExc_RuntimeError, std::string("loadVolumes() Failed to load data set '" +
            std::string(filename) + "': " + e.what()).c_str());
        return 0;
    }
}

static PyObject* voreen_loadTransferFunction(PyObject* /*self*/, PyObject* args) {

    // parse arguments
    char* processorName = 0;
    char* propertyID = 0;
    char* filename = 0;
    PyArg_ParseTuple(args, "sss:loadTransferFunction", &processorName, &propertyID, &filename);
    if (PyErr_Occurred())
        return 0;

    // find property and set position
    if (TransFunc1DKeysProperty* property = getTypedProperty<TransFunc1DKeysProperty>(
            std::string(processorName), std::string(propertyID), "TransFunc", "loadTransferFunction")) {
        TransFunc1DKeys* transFunc = property->get();
        if (!transFunc) {
            PyErr_SetString(PyExc_SystemError, std::string("loadTransferFunction() Property '" +
                property->getFullyQualifiedID() + "' does not contain a transfer function").c_str());
            return 0;
        }
        else if (!transFunc->load(std::string(filename))) {
            PyErr_SetString(PyExc_ValueError, std::string("loadTransferFunction() Failed to load '" +
                std::string(filename) + "'").c_str());
            return 0;
        }

        property->invalidate();
        Py_RETURN_NONE;
    }

    // indicate failure
    return 0;
}

static PyObject* voreen_render(PyObject* /*self*/, PyObject* args) {
    using namespace voreen;

    int sync = 0;
    if (!PyArg_ParseTuple(args, "|i:render", &sync))
        return 0;

    if (VoreenApplication::app() && VoreenApplication::app()->getNetworkEvaluator()) {
        VoreenApplication::app()->getNetworkEvaluator()->process();
        if (sync)
            glFinish();
        Py_RETURN_NONE;
    }

    return 0;
}

static PyObject* voreen_repaint(PyObject* /*self*/, PyObject* args) {
    if (!PyArg_ParseTuple(args, ""))
        return NULL;

    using namespace voreen;

    ProcessorNetwork* network = getProcessorNetwork("repaint");
    if (!network)
        return 0;

    std::vector<CanvasRenderer*> canvases = network->getProcessorsByType<CanvasRenderer>();
    for (size_t i=0; i<canvases.size(); i++) {
        if (canvases[i]->getCanvas())
            canvases[i]->getCanvas()->repaint();
    }

    Py_RETURN_NONE;
}

static PyObject* voreen_setViewport(PyObject* /*self*/, PyObject* args) {
    int i1, i2;

    using namespace voreen;

    if (!PyArg_ParseTuple(args, "ii:setViewport", &i1, &i2))
        return 0;

    ProcessorNetwork* network = getProcessorNetwork("setViewport");
    if (!network)
        return 0;
    std::vector<CanvasRenderer*> canvases = network->getProcessorsByType<CanvasRenderer>();

    for (std::vector<CanvasRenderer*>::const_iterator iter = canvases.begin(); iter != canvases.end(); iter++) {
        //CanvasRenderer* canvasProc = dynamic_cast<CanvasRenderer*>(*iter);
        //if (canvasProc && canvasProc->getProcessorWidget())
            //canvasProc->getProcessorWidget()->setSize(i1, i2);
        IntVec2Property* sizeProp = dynamic_cast<IntVec2Property*>((*iter)->getProperty("canvasSize"));
        if (sizeProp)
            sizeProp->set(tgt::ivec2(i1, i2));
        else
            return 0;
    }

    Py_RETURN_NONE;
}

static PyObject* voreen_snapshot(PyObject* /*self*/, PyObject* args) {

    const char* filename = 0;
    const char* canvasStr = 0;
    if (!PyArg_ParseTuple(args, "s|s:snapshot", &filename, &canvasStr))
        return 0;

    ProcessorNetwork* network = getProcessorNetwork("snapshot");
    if (!network)
        return 0;

    CanvasRenderer* canvasProc = 0;
    if (!canvasStr) {
        // select first canvas in network
        std::vector<CanvasRenderer*> canvases = network->getProcessorsByType<CanvasRenderer>();
        if (canvases.empty()) {
            PyErr_SetString(PyExc_RuntimeError, "snapshot() Network does not contain a CanvasRenderer.");
            return 0;
        }
        canvasProc = canvases.front();
    }
    else {
        // retrieve canvas with given name from network
        canvasProc = getTypedProcessor<CanvasRenderer>(std::string(canvasStr), "CanvasRenderer", "snapshot");
        if (!canvasProc)
            return 0;
    }
    tgtAssert(canvasProc, "no canvas proc");

    // take snapshot
    bool success = canvasProc->renderToImage(filename);
    if (!success) {
        PyErr_SetString(PyExc_ValueError, (std::string("snapshot() ") + canvasProc->getRenderToImageError()).c_str());
        return 0;
    }

    Py_RETURN_NONE;
}

static PyObject* voreen_canvas_count(PyObject* /*self*/, PyObject* args) {
    if (!PyArg_ParseTuple(args, ":canvasCount"))
        return NULL;

    using namespace voreen;

    ProcessorNetwork* network = getProcessorNetwork("canvasCount");
    if (!network)
        return 0;
    std::vector<CanvasRenderer*> canvases = network->getProcessorsByType<CanvasRenderer>();

    PyObject* py_count = Py_BuildValue("i", canvases.size());
    return py_count;
}

static PyObject* voreen_canvas_snapshot(PyObject* /*self*/, PyObject* args) {
    int index;
    const char* filename = 0;
    if (!PyArg_ParseTuple(args, "is:canvasSnapshot", &index, &filename))
        return NULL;

    using namespace voreen;

    ProcessorNetwork* network = getProcessorNetwork("canvasSnapshot");
    if (!network)
        return 0;
    std::vector<CanvasRenderer*> canvases = network->getProcessorsByType<CanvasRenderer>();

    int count = 0;
    for (std::vector<CanvasRenderer*>::const_iterator iter = canvases.begin(); iter != canvases.end(); iter++) {
        if (count == index) {
            bool success = (*iter)->renderToImage(filename, (*iter)->getCanvas()->getSize());
            if (!success) {
                PyErr_SetString(PyExc_RuntimeError, (std::string("canvasSnapshot() renderToImage() failed: ")
                    + (*iter)->getRenderToImageError()).c_str());
                return 0;
            }
            break;
        }
        count++;
    }

    Py_RETURN_NONE;
}

static PyObject* voreen_rotateCamera(PyObject* /*self*/, PyObject* args) {
    using namespace voreen;

    // parse arguments
    char* processorName = 0;
    char* propertyID = 0;
    float f1, f2, f3, f4;
    PyArg_ParseTuple(args, "ssf(fff):rotateCamera", &processorName, &propertyID,
        &f1, &f2, &f3, &f4);
    if (PyErr_Occurred())
        return 0;

    // find property
    CameraProperty* camProp = getTypedProperty<CameraProperty>(
        std::string(processorName), std::string(propertyID), "Camera", "rotateCamera");
    if (!camProp)
        return 0;

    // rotate by trackball
    VoreenTrackball track(camProp);
    track.setCenter(camProp->get().getFocus());
    track.rotate(tgt::quat::createQuat(f1, tgt::vec3(f2, f3, f4)));
    camProp->invalidate();

    Py_RETURN_NONE;
}

static PyObject* voreen_invalidateProcessors(PyObject* /*self*/, PyObject* args) {
    if (!PyArg_ParseTuple(args, ""))
        return NULL;

    using namespace voreen;
    if (VoreenApplication::app() && VoreenApplication::app()->getNetworkEvaluator()) {
        VoreenApplication::app()->getNetworkEvaluator()->invalidateProcessors();
        Py_RETURN_NONE;
    }

    return 0;
}

static PyObject* voreen_tickClockProcessor(PyObject* /*self*/, PyObject* args) {

    // parse arguments
    char* processorName = 0;
    PyArg_ParseTuple(args, "s:tickClockProcessor", &processorName);
    if (PyErr_Occurred())
        return 0;

#ifdef VRN_MODULE_BASE
    // find clock processor
    ClockProcessor* clockProcessor = getTypedProcessor<ClockProcessor>(
        std::string(processorName), "ClockProcessor", "tickClockProcessor");
    if (!clockProcessor)
        return 0;

    // trigger timer event
    clockProcessor->timerEvent(0);
    Py_RETURN_NONE;
#endif

    return 0;
}

static PyObject* voreen_resetClockProcessor(PyObject* /*self*/, PyObject* args) {
    // parse arguments
    char* processorName = 0;
    PyArg_ParseTuple(args, "s:resetClockProcessor", &processorName);
    if (PyErr_Occurred())
        return 0;

#ifdef VRN_MODULE_BASE
    // find clock processor
    ClockProcessor* clockProcessor = getTypedProcessor<ClockProcessor>(
        std::string(processorName), "ClockProcessor", "resetClockProcessor");
    if (!clockProcessor)
        return 0;

    // reset
    clockProcessor->resetCounter();
    Py_RETURN_NONE;
#endif

    return 0;
}

static PyObject* voreen_getBasePath(PyObject* /*self*/, PyObject* /*args*/) {

    if (!VoreenApplication::app()) {
        PyErr_SetString(PyExc_SystemError, "getBasePath() VoreenApplication not instantiated.");
        return 0;
    }

    std::string basePath = VoreenApplication::app()->getBasePath();
    PyObject* arg = Py_BuildValue("s", basePath.c_str());
    if (!arg) {
        PyErr_SetString(PyExc_SystemError, "getBasePath() failed to create argument");
        return 0;
    }

    return arg;
}

static PyObject* voreen_getRevision(PyObject* /*self*/, PyObject* /*args*/) {
    return Py_BuildValue("s", VoreenVersion::getRevision().c_str());
}

static PyObject* voreen_info(PyObject* /*self*/, PyObject* /*args*/) {
    return printModuleInfo("voreen",  true, 0, false, true);
}

//------------------------------------------------------------------------------
// Python binding method tables

// module 'voreen'
static PyMethodDef voreen_methods[] = {
    {
        "setPropertyValue",
        voreen_setPropertyValue,
        METH_VARARGS,
        "setPropertyValue(processor name, property id, scalar or tuple)\n\n"
        "Assigns a value to a processor property. The value has to be passed\n"
        "as scalar or tuple, depending on the property's cardinality.\n"
        "Camera properties take a 3-tuple of 3-tuples, containing the position,\n"
        "focus and up vectors. Option properties expect an option key."
    },
    {
        "getPropertyValue",
        voreen_getPropertyValue,
        METH_VARARGS,
        "getPropertyValue(processor name, property id) -> scalar or tuple\n\n"
        "Returns the value of a processor property as scalar or tuple,\n"
        "depending on the property's cardinality. See: setPropertyValue"
    },
    {
        "setPropertyMinValue",
        voreen_setPropertyMinValue,
        METH_VARARGS,
        "setPropertyMinValue(processor name, property id, scalar or tuple)\n\n"
        "Defines the minimum value of a numeric property."
    },
    {
        "setPropertyMaxValue",
        voreen_setPropertyMaxValue,
        METH_VARARGS,
        "setPropertyMaxValue(processor name, property id, scalar or tuple)\n\n"
        "Defines the maximum value of a numeric property."
    },
    {
        "getPropertyMinValue",
        voreen_getPropertyMinValue,
        METH_VARARGS,
        "getPropertyMinValue(processor name, property id) -> scalar or tuple\n\n"
        "Returns the minimum value of a numeric property as scalar or tuple,\n"
        "depending on the property's cardinality."
    },
    {
        "getPropertyMaxValue",
        voreen_getPropertyMaxValue,
        METH_VARARGS,
        "getPropertyMaxValue(processor name, property id) -> scalar or tuple\n\n"
        "Returns the maximum value of a numeric property as scalar or tuple,\n"
        "depending on the property's cardinality."
    },
    {
        "setPortData",
        voreen_setPortData,
        METH_VARARGS,
        "setPortData(processor name, port id, data)\n\n"
        "Assigns data to a processor port.\n"
    },
    {
        "getPortData",
        voreen_getPortData,
        METH_VARARGS,
        "getPortData(processor name, port id) -> data\n\n"
        "Returns the data of a processor port,\n"
        "depending on the port's type. See: setPortData"
    },
    {
        "setCameraPosition",
        voreen_setCameraPosition,
        METH_VARARGS,
        "setCameraPosition(processor name, property id, (x, y, z))\n\n"
        "Convenience function for setting a camera's position.\n"
        "See also: setPropertyValue"
    },
    {
        "setCameraFocus",
        voreen_setCameraFocus,
        METH_VARARGS,
        "setCameraFocus(processor name, property id, (x, y, z))\n\n"
        "Convenience function for setting a camera's focus.\n"
        "See also: setPropertyValue"
    },
    {
        "setCameraUpVector",
        voreen_setCameraUpVector,
        METH_VARARGS,
        "setCameraUp(processor name, property id, (x, y, z))\n\n"
        "Convenience function for setting a camera's up vector.\n"
        "See also: setPropertyValue"
    },
    {
        "loadVolume",
        voreen_loadVolume,
        METH_VARARGS,
        "loadVolume(filename, [volume source])\n\n"
        "Loads a volume data set and assigns it to a VolumeSource processor.\n"
        "If no processor name is passed, the first volume source in the\n"
        "network is chosen."
    },
    {
        "loadVolumes",
        voreen_loadVolumes,
        METH_VARARGS,
        "loadVolumes(filename, selected, clear, [volume list source])\n\n"
        "Loads all volumes and assigns them to a VolumeListSource processor.\n"
        "If no processor name is passed, the first volume list source in the\n"
        "network is chosen."
    },
    {
        "loadTransferFunction",
        voreen_loadTransferFunction,
        METH_VARARGS,
        "loadTransferFunction(processor name, property id, filename)\n\n"
        "Loads a transfer function and assigns it to a transfer function property."
    },
    {
        "render",
        voreen_render,
        METH_VARARGS,
        "render([sync=0])\n\n"
        "Renders the current network by calling NetworkEvaluator::process().\n"
        "For sync=1, glFinish() is called afterwards."
    },
    {
        "repaint",
        voreen_repaint,
        METH_VARARGS,
        "repaint() Renders the network by forcing repaints of all canvases."
    },
    {
        "setViewport",
        voreen_setViewport,
        METH_VARARGS,
        "setViewport(width, height)\n\n"
        "Convenience function setting the canvas dimensions of\n"
        "of all CanvasRenderers in the network."
    },
    {
        "snapshot",
        voreen_snapshot,
        METH_VARARGS,
        "snapshot(filename, [canvas])\n\n"
        "Saves a snapshot of the specified canvas to the given file.\n"
        "If no canvas name is passed, the first canvas in the network is chosen."
    },
    {
        "snapshotCanvas",
        voreen_canvas_snapshot,
        METH_VARARGS,
        "snapshotCanvas(i, filename)\n\n"
        "Saves a snapshot of the ith canvas to the given file."
    },
    {
        "canvasCount",
        voreen_canvas_count,
        METH_VARARGS,
        "canvasCount() -> int\n\n"
        "Returns the number of canvases in the current network."
    },
    {
        "rotateCamera",
        voreen_rotateCamera,
        METH_VARARGS,
        "rotateCamera(processor name, property id, angle, (x,y,z))\n\n"
        "Rotates a camera by the specified angle around the specified axis."
    },
    {
        "invalidateProcessors",
        voreen_invalidateProcessors,
        METH_VARARGS,
        "invalidateProcessors() Invalidates all processors in the current network."
    },
    {
        "tickClockProcessor",
        voreen_tickClockProcessor,
        METH_VARARGS,
        "tickClockProcessor(processor name)\n\n"
        "Sends a timer event to a ClockProcessor."
    },
    {
        "resetClockProcessor",
        voreen_resetClockProcessor,
        METH_VARARGS,
        "resetClockProcessor(processor name)\n\n"
        "Resets the clock of a ClockProcessor."
    },
    {
        "getBasePath",
        voreen_getBasePath,
        METH_VARARGS,
        "getBasePath() -> path\n\n"
        "Returns the absolute Voreen base path."
    },
    {
        "getRevision",
        voreen_getRevision,
        METH_VARARGS,
        "getRevision() -> string\n\n"
        "Returns the revision of the Voreen binary."
    },
    {
        "info",
        voreen_info,
        METH_VARARGS,
        "info() Prints documentation of the module's functions."
    },
    { NULL, NULL, 0, NULL} // sentinal
};

static struct PyModuleDef voreenModuleDef =
{
    PyModuleDef_HEAD_INIT,
    "voreen",
    NULL,
    -1,
    voreen_methods
};

PyMODINIT_FUNC
PyInit_voreenModule(void)
{
    if (PyType_Ready(&VolumeObjectType) < 0)
        return NULL;

    if (PyType_Ready(&RenderTargetObjectType) < 0)
        return NULL;

    PyObject* m = PyModule_Create(&voreenModuleDef);
    if(m == NULL)
        return NULL;

    // Register custom objects.
    Py_INCREF(&VolumeObjectType);
    PyModule_AddObject(m, "Volume", (PyObject *) &VolumeObjectType);
    Py_INCREF(&RenderTargetObjectType);
    PyModule_AddObject(m, "RenderTarget", (PyObject *) &RenderTargetObjectType);

    return m;
}

namespace voreen {

const std::string PyVoreen::loggerCat_ = "voreen.Python.PyVoreen";

PyVoreen::PyVoreen() {
    if (!Py_IsInitialized()) {
        // initialize voreen module
        if(PyImport_AppendInittab("voreen", PyInit_voreenModule) == -1)
            LWARNING("Failed to init helper module 'voreen'");
    }
    else {
        LERROR("Python environment already initialized");
    }
}

} // namespace voreen


//-------------------------------------------------------------------------------------------------
// implementation of helper functions

namespace {

ProcessorNetwork* getProcessorNetwork(const std::string& functionName) {

    // retrieve evaluator from application
    if (!VoreenApplication::app()) {
        PyErr_SetString(PyExc_SystemError, std::string(functionName + "() VoreenApplication not instantiated").c_str());
        return 0;
    }
    else if (!VoreenApplication::app()->getNetworkEvaluator()) {
        PyErr_SetString(PyExc_SystemError, std::string(functionName + "() No network evaluator").c_str());
        return 0;
    }

    // get network from evaluator
    ProcessorNetwork* network = const_cast<ProcessorNetwork*>(VoreenApplication::app()->getNetworkEvaluator()->getProcessorNetwork());
    if (!network) {
        PyErr_SetString(PyExc_SystemError, std::string(functionName + "() No processor network").c_str());
        return 0;
    }

    return network;
}

Processor* getProcessor(const std::string& processorName, const std::string& functionName) {

    ProcessorNetwork* network = getProcessorNetwork(functionName);
    if (!network)
        return 0;

    // find processor
    Processor* processor = network->getProcessorByName(processorName);
    if (!processor) {
        PyErr_SetString(PyExc_NameError, std::string(functionName + "() Processor '" + processorName + "' not found").c_str());
        return 0;
    }

    return processor;
}

template<typename T>
T* getTypedProcessor(const std::string& processorName, const std::string& processorTypeString,
                     const std::string& functionName) {

    // fetch processor
    Processor* processor = getProcessor(processorName, functionName);
    if (!processor)
        return 0;

    // check type
    if (T* cProc = dynamic_cast<T*>(processor))
        return cProc;
    else {
        PyErr_SetString(PyExc_TypeError, std::string(functionName + "() Processor '" +
            processorName + "' is not of type " + processorTypeString).c_str());
        return 0;
    }
}

Property* getProperty(const std::string& processorName, const std::string& propertyID,
                      const std::string& functionName) {

    // fetch processor
    Processor* processor = getProcessor(processorName, functionName);
    if (!processor)
        return 0;

    // find property
    Property* property = processor->getProperty(propertyID);
    if (!property) {
        PyErr_SetString(PyExc_NameError, std::string(functionName + "() Processor '" +
            processorName + "' has no property '" + propertyID + "'").c_str());
        return 0;
    }

    return property;
}

template<typename T>
T* getTypedProperty(const std::string& processorName, const std::string& propertyID,
               const std::string& propertyTypeString, const std::string& functionName) {

    // fetch property
    Property* property = getProperty(processorName, propertyID, functionName);
    if (!property)
        return 0;

    // check type
    if (T* cProp = dynamic_cast<T*>(property))
        return cProp;
    else {
        PyErr_SetString(PyExc_TypeError, std::string(functionName + "() Property '" +
            property->getFullyQualifiedID() + "' is of type " + property->getTypeDescription() +
            ". Expected: " + propertyTypeString).c_str());
        return 0;
    }
}

template<typename PropertyType, typename ValueType>
bool setPropertyValue(const std::string& processorName, const std::string& propertyID, const ValueType& value,
                      const std::string& propertyTypeString, const std::string& functionName) {

    if (PropertyType* property = getTypedProperty<PropertyType>(
            processorName, propertyID, propertyTypeString, functionName)) {
        std::string errorMsg;
        if (property->isValidValue(value, errorMsg)) {
            property->set(value);
            return true;
        }
        else {
            // define Python exception
            PyErr_SetString(PyExc_ValueError, (functionName + std::string("() ") + errorMsg).c_str());
            return false;
        }
    }
    return false;
}

template<typename PropertyType, typename ValueType>
bool setPropertyValue(PropertyType* property, const ValueType& value,
                      const std::string& functionName) {
    tgtAssert(property, "Null pointer passed");

    std::string errorMsg;
    if (property->isValidValue(value, errorMsg)) {
        property->set(value);
        return true;
    }
    else {
        // define Python exception
        PyErr_SetString(PyExc_ValueError, (functionName + std::string("() ") + errorMsg).c_str());
        return false;
    }
}

Port* getPort(const std::string& processorName, const std::string& portID, const std::string& functionName) {

    // fetch processor
    Processor* processor = getProcessor(processorName, functionName);
    if (!processor)
        return 0;

    // find property
    Port* port = processor->getPort(portID);
    if (!port) {
        PyErr_SetString(PyExc_NameError, std::string(functionName + "() Processor '" +
                                                     processorName + "' has no port '" + portID + "'").c_str());
        return 0;
    }

    return port;
}

template<typename T>
T* getTypedPort(const std::string& processorName, const std::string& portID,
                const std::string& portTypeString, const std::string& functionName) {

    // fetch port
    Port* port = getPort(processorName, portID, functionName);
    if (!port)
        return 0;

    // check type
    if (T* cProp = dynamic_cast<T*>(port))
        return cProp;
    else {
        PyErr_SetString(PyExc_TypeError, std::string(functionName + "() Port '" +
                                                     port->getQualifiedName() + "' is not of type " +
                                                     portTypeString).c_str());
        return 0;
    }
}

static PyObject* printModuleInfo(const std::string& moduleName, bool omitFunctionName,
                                 int spacing, bool collapse, bool blanklines) {

    // import apihelper.py
    PyObject* apihelper = PyImport_ImportModule("apihelper");
    if (!apihelper) {
        PyErr_SetString(PyExc_SystemError,
            (moduleName + std::string(".info() apihelper module not found")).c_str());
        return 0;
    }

    // get reference to info function
    PyObject* func = PyDict_GetItemString(PyModule_GetDict(apihelper), "info");
    if (!func) {
        PyErr_SetString(PyExc_SystemError,
            (moduleName + std::string(".info() apihelper.info() not found")).c_str());
        Py_XDECREF(apihelper);
        return 0;
    }

    // get reference to module
    PyObject* mod = PyImport_AddModule(moduleName.c_str()); // const cast required for
                                                                               // Python 2.4
    if (!mod) {
        PyErr_SetString(PyExc_SystemError,
            (moduleName + std::string(".info() failed to access module ") + moduleName).c_str());
        Py_XDECREF(apihelper);
    }

    // build parameter tuple
    std::string docStr = "Module " + moduleName;
    PyObject* arg = Py_BuildValue("(O,s,i,i,i,i)", mod, docStr.c_str(),
        int(omitFunctionName), spacing, (int)collapse, (int)blanklines);
    if (!arg) {
        PyErr_SetString(PyExc_SystemError,
            (moduleName + std::string(".info() failed to create arguments")).c_str());
        Py_XDECREF(apihelper);
        return 0;
    }

    PyObject* callObj = PyObject_CallObject(func, arg);
    bool success = (callObj != 0);

    Py_XDECREF(callObj);
    Py_XDECREF(arg);
    Py_XDECREF(apihelper);

    if (success)
        Py_RETURN_NONE;
    else
        return 0;
}

} // namespace anonymous
