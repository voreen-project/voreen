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

#include "voreen/core/network/networkconverter.h"
#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/property.h"
#include "voreen/core/utils/stringutils.h"

#include "tinyxml/tinyxml.h"
#include "tgt/logmanager.h"

#include <map>

namespace voreen {

const std::string NetworkConverter::loggerCat_ = "voreen.NetworkConverter";

void NetworkConverter::changePropertyType(TiXmlElement* elem, const std::string& from, const std::string& to) {
    TiXmlElement* propertyLinkNode = elem->FirstChildElement("PropertyLinks");
    for (TiXmlElement* linkNode = propertyLinkNode->FirstChildElement("PropertyLink")
        ; linkNode != 0
        ; linkNode = linkNode->NextSiblingElement("PropertyLink"))
    {
        TiXmlElement* src = linkNode->FirstChildElement("SourceProperty");
        if(src && src->Attribute(std::string("type")) && !from.compare(*(src->Attribute(std::string("type")))))
            src->SetAttribute("type",to);
        TiXmlElement* dst = linkNode->FirstChildElement("DestinationProperty");
        if(dst && dst->Attribute(std::string("type")) && !from.compare(*(dst->Attribute(std::string("type")))))
            dst->SetAttribute("type",to);
    }
}

void NetworkConverter::changePropertyName(TiXmlElement* elem, const std::string& ofProcessor, const std::string& from, const std::string& to) {
    TiXmlElement* processorsNode = elem->FirstChildElement("Processors");
    for (TiXmlElement* procNode = processorsNode->FirstChildElement("Processor")
        ; procNode != 0
        ; procNode = procNode->NextSiblingElement("Processor"))
    {
        const std::string* type = procNode->Attribute(std::string("type"));
        if (ofProcessor == "" || (type && *type == ofProcessor)) {
            TiXmlElement* propertiesNode = procNode->FirstChildElement("Properties");
            for (TiXmlElement* propNode = propertiesNode->FirstChildElement("Property")
                ; propNode != 0
                ; propNode = propNode->NextSiblingElement("Property"))
            {
                const std::string* name = propNode->Attribute(std::string("name"));
                if (name && *name == from)
                    propNode->SetAttribute("name", to);
            }
        }
    }
}

void NetworkConverter::changePropertyValue(TiXmlElement* elem, const std::string& ofProcessor, const std::string& from, const std::string& to) {
    TiXmlElement* processorsNode = elem->FirstChildElement("Processors");
    for (TiXmlElement* procNode = processorsNode->FirstChildElement("Processor")
        ; procNode != 0
        ; procNode = procNode->NextSiblingElement("Processor"))
    {
        const std::string* type = procNode->Attribute(std::string("type"));
        if (ofProcessor == "" || (type && *type == ofProcessor)) {
            TiXmlElement* propertiesNode = procNode->FirstChildElement("Properties");
            for (TiXmlElement* propNode = propertiesNode->FirstChildElement("Property")
                ; propNode != 0
                ; propNode = propNode->NextSiblingElement("Property"))
            {
                const std::string* value = propNode->Attribute(std::string("value"));
                if (value && *value == from)
                    propNode->SetAttribute("value", to);
            }
        }
    }
}

void NetworkConverter::changeProcessorType(TiXmlElement* elem, const std::string& from, const std::string& to) {
    TiXmlElement* processorsNode = elem->FirstChildElement("Processors");
    for (TiXmlElement* node = processorsNode->FirstChildElement("Processor"); node != 0; node = node->NextSiblingElement("Processor")) {
        const std::string* type = node->Attribute(std::string("type"));
        if (type && *type == from) {
            node->SetAttribute("type", to);
            const std::string* name = node->Attribute(std::string("name"));
            if (name && *name == from)
                node->SetAttribute("name", to);
        }
    }
}

void NetworkConverter::changeLinkEvaluatorType(TiXmlElement* elem, const std::string& from, const std::string& to) {
    TiXmlElement* propertyLinkNode = elem->FirstChildElement("PropertyLinks");
    for (TiXmlElement* linkNode = propertyLinkNode->FirstChildElement("PropertyLink")
        ; linkNode != 0
        ; linkNode = linkNode->NextSiblingElement("PropertyLink"))
    {
        TiXmlElement* eval = linkNode->FirstChildElement("Evaluator");
        if(eval && eval->Attribute(std::string("type")) && !from.compare(*(eval->Attribute(std::string("type")))))
            eval->SetAttribute("type",to);
    }
}


//-------------------------------------------------------------------------------------------------

void NetworkConverter4to5::convert(TiXmlElement* elem) {
    //changePropertyName(elem, "", "OldPropertyDescriptor", "NewPropertyDescriptor");
    //changeProcessorType(elem, "OldProcessorClassname", "NewProcessorClassname");
    // AxialSliceProxyGeometry processor
    changePropertyName(elem, "", "set.Axis", "axis");
    changePropertyName(elem, "", "set.Begin", "begin");
    changePropertyName(elem, "", "set.Thickness", "thickness");
    // ClippingPlaneWidget processor
    changePropertyName(elem, "", "set.ClipSliceColor", "clipSliceColor");
    changePropertyName(elem, "", "set.polyOffsetFact", "polyOffsetFact");
    changePropertyName(elem, "", "set.polyOffsetUnit", "polyOffsetUnit");
    changePropertyName(elem, "", "switch.ClipXLock", "switchClipXLock");
    changePropertyName(elem, "", "switch.ClipYLock", "switchClipYLock");
    changePropertyName(elem, "", "switch.ClipZLock", "switchClipZLock");
    // ColorDepth processor
    changePropertyName(elem, "", "set.colordepth.factor", "factor");
    // Compositer processor
    changePropertyName(elem, "", "Compositer.set.blendMode", "blendMode");
    changePropertyValue(elem, "", "Alpha Compositing", "Alpha compositing");
    changePropertyValue(elem, "", "Maximum Alpha", "Maximum alpha");
    changePropertyName(elem, "", "Compositer.set.blendFactor", "blendFactor");
    // CubeProxyGeometry processor
    changePropertyName(elem, "", "bricking.region.priority", "regionPriority");
    changePropertyName(elem, "", "switch.virtualClipplane", "enableVirtualClipplane");
    changePropertyName(elem, "", "set.virtualClipplane", "virtualClipplane");
    // Gaussian processor
    changeProcessorType(elem, "Blur", "Gaussian");
    changePropertyName(elem, "", "set.blurDelta", "sigma");
    changePropertyName(elem, "", "set.blurRed", "blurRed");
    changePropertyName(elem, "", "set.blurGreen", "blurGreen");
    changePropertyName(elem, "", "set.blurBlue", "blurBlue");
    changePropertyName(elem, "", "set.blurAlpha", "blurAlpha");
    // GeomBoundingBox/SlicePositionRenderer processor
    changePropertyName(elem, "", "set.BoundingboxColor", "boundingboxColor");
    changePropertyName(elem, "", "set.BoundingBoxWidth", "boundingBoxWidth");
    changePropertyName(elem, "", "set.BoundingBoxStippleFactor", "boundingBoxStippleFactor");
    changePropertyName(elem, "", "set.BoundingBoxStipplePattern", "boundingBoxStipplePattern");
    // GeomLightWidget processor
    changePropertyName(elem, "", "set.showLightWidget", "showLightWidget");
    // GeomPlane processor
    changePropertyName(elem, "", "set.PlaneXSize", "planeXSize");
    changePropertyName(elem, "", "set.PlaneYSize", "planeYSize");
    // Grayscale processor
    changeProcessorType(elem, "Greyscale", "Grayscale");
    // EdgeDetect processor
    changePropertyName(elem, "", "set.edgeThreshold", "edgeThreshold");
    changePropertyName(elem, "", "set.showImage", "showImage");
    changePropertyName(elem, "", "set.edgeColor", "edgeColor");
    // IDRaycaster processor
    changePropertyName(elem, "", "set.penetrationDepth", "penetrationDepth");
    // ImageOverlay processor
    changePropertyName(elem, "", "set.imageoverlayLeft", "left");
    changePropertyName(elem, "", "set.imageoverlayWidth", "width");
    changePropertyName(elem, "", "set.imageoverlayTop", "top");
    changePropertyName(elem, "", "set.imageoverlayHeight", "height");
    changePropertyName(elem, "", "set.imageoverlayOpacity", "opacity");
    // Labeling processor
    changePropertyName(elem, "", "set.labelingWidget", "labelingWidget");
    changePropertyName(elem, "", "set.Layout", "layout");
    changePropertyName(elem, "", "set.LabelColor", "labelColor");
    changePropertyName(elem, "", "set.HaloColor", "haloColor");
    changePropertyName(elem, "", "set.Font", "font");
    changePropertyName(elem, "", "set.FontSize", "fontSize");
    changePropertyName(elem, "", "switch.drawHalo", "drawHalo");
    changePropertyName(elem, "", "set.SegmentDescriptionFile", "segmentDescriptionFile");
    changePropertyName(elem, "", "set.lockInternalFontSettings", "lockInternalFontSettings");
    changePropertyName(elem, "", "set.labelColorIntern", "labelColorIntern");
    changePropertyName(elem, "", "set.haloColorIntern", "haloColorIntern");
    changePropertyName(elem, "", "set.fontSizeIntern", "fontSizeIntern");
    changePropertyName(elem, "", "set.shape3D", "shape3D");
    changePropertyName(elem, "", "set.minSegmentSize", "minSegmentSize");
    changePropertyName(elem, "", "set.maxIterations", "maxIterations");
    changePropertyName(elem, "", "set.filterDelta", "filterDelta");
    changePropertyName(elem, "", "set.distanceMapStep", "distanceMapStep");
    changePropertyName(elem, "", "set.glyphAdvance", "glyphAdvance");
    changePropertyName(elem, "", "set.polynomialDegree", "polynomialDegree");
    changePropertyName(elem, "", "set.bezierHorzDegree", "bezierHorzDegree");
    changePropertyName(elem, "", "set.bezierVertDegree", "bezierVertDegree");
    changePropertyName(elem, "", "set.filterKernel", "filterKernel");
    changePropertyName(elem, "", "set.thresholdPPthreshold", "thresholdPPthreshold");
    changePropertyName(elem, "", "set.labelLayoutAsString", "labelLayoutAsString");
    changePropertyName(elem, "", "set.showLabels", "showLabels");
    // MultipleAxialSliceProxyGeometry processor
    changePropertyName(elem, "", "set.nSlicesPerRow", "nSlicesPerRow");
    changePropertyName(elem, "", "set.nSlicesPerCol", "nSlicesPerCol");
    // RenderProcessor processor
    changePropertyName(elem, "", "set.backgroundColor", "backgroundColor");
    // SliceProxyGeometry processor
    changePropertyName(elem, "", "set.slicePos", "slicePos");
    // Threshold processor
    changePropertyName(elem, "", "set.thresholdPPthreshold", "threshold");
    // Background processor
    changePropertyName(elem, "", "backgroundFirstColor", "color1");
    changePropertyName(elem, "", "backgroundSecondColor", "color2");
    changePropertyName(elem, "", "backgroundAngle", "angle");
    changePropertyName(elem, "", "backgroundfilenameAsString", "texture");
    changePropertyName(elem, "", "backgroundtile", "repeat");
}

void NetworkConverter5to6::convert(TiXmlElement* elem) {
    changeProcessorType(elem, "Compositer", "Compositor");
    changeProcessorType(elem, "LoopCompositer", "LoopCompositor");
    changeProcessorType(elem, "VVCompositer", "VVCompositor");
    changeProcessorType(elem, "SingleSliceRenderer", "SliceViewer");
    changeProcessorType(elem, "SliceSequenceRenderer", "SliceViewer");
    changeProcessorType(elem, "GeomLightWidget", "LightWidgetRenderer");
}


void NetworkConverter6to7::convert(TiXmlElement* elem) {
    changeProcessorType(elem, "VolumeConvert", "VolumeBitScale");
    changeProcessorType(elem, "VolumeMirrorZ", "VolumeMirror");
    changeProcessorType(elem, "TexturedGeometryRenderer", "GeometryRenderer");
    changeProcessorType(elem, "ClippingPlaneWidget", "PlaneWidgetProcessor");
    changeProcessorType(elem, "CylinderRenderer", "QuadricRenderer");
    changeProcessorType(elem, "SphereRenderer", "QuadricRenderer");
}

void NetworkConverter7to8::convert(TiXmlElement* elem) {
    changeProcessorType(elem, "VolumeSubSet", "VolumeCrop");
}

void NetworkConverter8to9::convert(TiXmlElement* elem) {
    changeProcessorType(elem, "FiberWriter", "FiberSave");
    changeProcessorType(elem, "RawTextureWriter", "RawTextureSave");
}

void NetworkConverter9to10::convert(TiXmlElement* elem) {
    changeProcessorType(elem, "TensorEllipsoidRenderer", "TensorGlyphRenderer");
    changeProcessorType(elem, "FiberTracker", "FiberTrackerFACT");
    changeProcessorType(elem, "TensorlineTracker", "FiberTrackerTensorline");
    changeProcessorType(elem, "VolumeBitScale", "VolumeFormatConversion");
    changeProcessorType(elem, "BoundingBox", "BoundingBoxRenderer");
}

void NetworkConverter10to11::convert(TiXmlElement* elem) {
    changeProcessorType(elem, "CubeMeshProxyGeometry", "CubeProxyGeometry");
}

// -------------------------------------------------------------------------------------

void NetworkConverter11to12::convert(TiXmlElement* elem) {
    changeProcessorType(elem, "CoordinateTransformation", "GeometryTransformationVolume");
    changeProcessorType(elem, "MeshClipping", "GeometryClipping");
    changeProcessorType(elem, "MeshSlabClipping", "GeometrySlabClipping");
    changeProcessorType(elem, "MeshClippingWidget", "GeometryClippingWidget");
    changeProcessorType(elem, "VolumeScaling", "VolumeSpacing");
    changeProcessorType(elem, "VolumeTranslation", "VolumeOffset");
    changeProcessorType(elem, "VolumeRegistration", "VolumeLandmarkRegistration");
}

void NetworkConverter11to12::convertVolumeContainer(TiXmlElement* workspaceNode) {
    if (!workspaceNode) {
        LERRORC("voreen.NetworkConverter11to12", "no workspace node");
        return;
    }

    // construct map from Volume-refID to URL
    std::map<std::string, std::string> refToUrl;
    TiXmlElement* volumeContainerNode = workspaceNode->FirstChildElement("VolumeContainer");
    if (!volumeContainerNode)
        return;
    LWARNINGC("voreen.NetworkConverter11to12", "converting volume container ...");
    TiXmlElement* volumeHandlesNode = volumeContainerNode->FirstChildElement("VolumeHandles");
    if (!volumeHandlesNode)
        return;
    for (TiXmlElement* handleNode = volumeHandlesNode->FirstChildElement("VolumeHandle"); handleNode != 0; handleNode = handleNode->NextSiblingElement("VolumeHandle")) {
        const std::string* refID = handleNode->Attribute(std::string("id"));
        if (!refID)
            continue;
        TiXmlElement* originNode = handleNode->FirstChildElement("Origin");
        if (!originNode)
            continue;
        const std::string* url = originNode->Attribute(std::string("url"));
        if (!url)
            url = originNode->Attribute(std::string("filename"));
        if (!url)
            continue;

        tgtAssert(refID, "no refID");
        tgtAssert(url, "no url");
        refToUrl.insert(std::pair<std::string, std::string>(*refID, *url));
    }

    // look for VolumeSource processors and replace their VolumeHandleProperty (with specific refID)
    // with a VolumeURLProperty with the corresponding volume URL
    TiXmlElement* networkNode = workspaceNode->FirstChildElement("ProcessorNetwork");
    if (!networkNode)
        return;
    TiXmlElement* processorsNode = networkNode->FirstChildElement("Processors");
    if (!processorsNode)
        return;
    for (TiXmlElement* procNode = processorsNode->FirstChildElement("Processor"); procNode != 0; procNode = procNode->NextSiblingElement("Processor")) {
        const std::string* type = procNode->Attribute(std::string("type"));
        if (!type || *type != "VolumeSource")
            continue;

        // it's a VolumeSource => extract VolumeHandleProperty node
        TiXmlElement* propertiesNode = procNode->FirstChildElement("Properties");
        if (!propertiesNode)
            continue;
        TiXmlElement* handlePropNode = propertiesNode->FirstChildElement("Property");
        if (!handlePropNode)
            continue;
        const std::string* propName = handlePropNode->Attribute(std::string("name"));
        if (!propName || *propName != "volumeHandle")
            continue;

        // convert VolumeHandleProperty to VolumeURLProperty
        handlePropNode->SetAttribute("name", "volumeURL");

        // retrieve Volume reference
        TiXmlElement* valueNode = handlePropNode->FirstChildElement("value");
        if (!valueNode)
            continue;
        const std::string* refID = valueNode->Attribute(std::string("ref"));
        if (!refID)
            continue;

        // insert referenced URL into converted VolumeHandleProperty
        if (refToUrl.find(*refID) == refToUrl.end())
            LWARNINGC("voreen.NetworkConverter11to12", "convertVolumeContainer(): unknown refID: " << *refID);
        else
            handlePropNode->SetAttribute("url", refToUrl[*refID]);
    }


}
//----------------------------------------------------------------------------

void NetworkConverter12to13::convert(TiXmlElement* elem) {
    changePropertyName(elem, "UnaryImageProcessor", "shader", "shader.program");
    changePropertyName(elem, "BinaryImageProcessor", "shader", "shader.program");
}

//----------------------------------------------------------------------------

void NetworkConverter13to14::convert(TiXmlElement* elem) {
    if(TiXmlElement* processorsNode = elem->FirstChildElement("Processors")) {
        for (TiXmlElement* procNode = processorsNode->FirstChildElement("Processor")
            ; procNode != 0
            ; procNode = procNode->NextSiblingElement("Processor"))
        {
            //add mapKey to processor->properties->property
            if(TiXmlElement* propertiesNode = procNode->FirstChildElement("Properties")) {
                for (TiXmlElement* propNode = propertiesNode->FirstChildElement("Property")
                     ; propNode != 0
                     ; propNode = propNode->NextSiblingElement("Property"))
                {
                    const std::string* name = propNode->Attribute(std::string("name"));
                    if (name)
                        propNode->SetAttribute("mapKey", name->c_str());
                }
            }
            //add mapKey to processor->inports->port->properties->property
            if(TiXmlElement* inportsNode = procNode->FirstChildElement("Inports")) {
                for (TiXmlElement* portNode = inportsNode->FirstChildElement("Port")
                    ; portNode != 0
                    ; portNode = portNode->NextSiblingElement("Port"))
                {
                    if(TiXmlElement* propertiesNode = portNode->FirstChildElement("Properties")) {
                        for (TiXmlElement* propNode = propertiesNode->FirstChildElement("Property")
                            ; propNode != 0
                            ; propNode = propNode->NextSiblingElement("Property"))
                        {
                            const std::string* name = propNode->Attribute(std::string("name"));
                            if (name)
                                propNode->SetAttribute("mapKey", name->c_str());
                        }
                    }
                }
            }
            //add mapKey to processor->outports->port->properties-<property
            if(TiXmlElement* outportsNode = procNode->FirstChildElement("Outports")) {
                for (TiXmlElement* portNode = outportsNode->FirstChildElement("Port")
                    ; portNode != 0
                    ; portNode = portNode->NextSiblingElement("Port"))
                {
                    if(TiXmlElement* propertiesNode = portNode->FirstChildElement("Properties")) {
                        for (TiXmlElement* propNode = propertiesNode->FirstChildElement("Property")
                            ; propNode != 0
                            ; propNode = propNode->NextSiblingElement("Property"))
                        {
                            const std::string* name = propNode->Attribute(std::string("name"));
                            if (name)
                                propNode->SetAttribute("mapKey", name->c_str());
                        }
                    }
                }
            }
            //add mapKey to processor->interactionhandlers->handler->properties->property
            if(TiXmlElement* interactionNode = procNode->FirstChildElement("InteractionHandlers")) {
                for (TiXmlElement* handlerNode = interactionNode->FirstChildElement("Handler")
                    ; handlerNode != 0
                    ; handlerNode = handlerNode->NextSiblingElement("Handler"))
                {
                    if(TiXmlElement* propertiesNode = handlerNode->FirstChildElement("Properties")) {
                        for (TiXmlElement* propNode = propertiesNode->FirstChildElement("Property")
                            ; propNode != 0
                            ; propNode = propNode->NextSiblingElement("Property"))
                        {
                            const std::string* name = propNode->Attribute(std::string("name"));
                            if (name)
                                propNode->SetAttribute("mapKey", name->c_str());
                        }
                    }
                }
            }
        }
    }
}

void NetworkConverter14to15::convert(TiXmlElement* elem) {
    changeProcessorType(elem, "VolumeCollectionSource",         "VolumeListSource");
    changeProcessorType(elem, "VolumeCollectionSave",           "VolumeListSave");
    changeProcessorType(elem, "VolumeCollectionModalityFilter", "VolumeListModalityFilter");
}


void NetworkConverter15to16::convert(TiXmlElement* networkElem) {
    // convert OptimizedProxyGeometry.modeString values
    TiXmlElement* processorsNode = networkElem->FirstChildElement("Processors");
    if (!processorsNode) {
        LWARNING("No <Processors> node found in <ProcessorNetwork>");
        return;
    }
    for (TiXmlElement* node = processorsNode->FirstChildElement("Processor"); node; node = node->NextSiblingElement("Processor")) {
        const std::string* type = node->Attribute(std::string("type"));
        if (type && *type == "OptimizedProxyGeometry") {
            TiXmlElement* propertiesNode = node->FirstChildElement("Properties");
            if (propertiesNode) {
                for (TiXmlElement* propertyNode = propertiesNode->FirstChildElement("Property"); propertyNode; propertyNode = propertyNode->NextSiblingElement("Property")) {
                    const std::string* propName = propertyNode->Attribute(std::string("name"));
                    if (propName && *propName == "modeString") {
                        const std::string* valueStr = propertyNode->Attribute(std::string("value"));
                        if (valueStr && *valueStr == "cube")
                            propertyNode->SetAttribute("value", "boundingbox");
                        else if (valueStr && *valueStr == "mincube")
                            propertyNode->SetAttribute("value", "minboundingbox");
                        else if (valueStr && *valueStr == "maximalcubes")
                            propertyNode->SetAttribute("value", "visiblebricks");
                        else if (valueStr && *valueStr == "octreemax")
                            propertyNode->SetAttribute("value", "octreebricks");
                    }
                }
            }
        }
    }
}

void NetworkConverter16to17::convert(TiXmlElement* networkElem) {
    //rename property type from TransFuncProperty to TransFunc1DKeysProperty
    changePropertyType(networkElem,"TransFuncProperty","TransFunc1DKeysProperty");
    //rename link evaluator type from LinkEvaluatorTransFuncId to LinkEvaluatorTransFunc1DKeysId
    changeLinkEvaluatorType(networkElem,"LinkEvaluatorTransFuncId","LinkEvaluatorTransFunc1DKeysId");
    //rename meta data type from type="SelectionMetaData::Processor" to type="SerializableVectorMetaData::Processor" isOwner="false"
    if(TiXmlElement* metaDataNode = networkElem->FirstChildElement("MetaData")) {
        for (TiXmlElement* itemNode = metaDataNode->FirstChildElement("MetaItem")
            ; itemNode != 0
            ; itemNode = itemNode->NextSiblingElement("MetaItem"))
        {
            if(itemNode->Attribute(std::string("type")) && !std::string("SelectionMetaData::Processor").compare(*(itemNode->Attribute(std::string("type"))))) {
                itemNode->SetAttribute("type","SerializableVectorMetaData::Processor");
                itemNode->SetAttribute("isOwner","false");
            }
        }
    }

    //rename meta data name from name="EditorWindow" to name="ToolWindow"
    TiXmlElement* processorsNode = networkElem->FirstChildElement("Processors");
    for (TiXmlElement* procNode = processorsNode->FirstChildElement("Processor")
        ; procNode != 0
        ; procNode = procNode->NextSiblingElement("Processor"))
    {
        TiXmlElement* propertiesNode = procNode->FirstChildElement("Properties");
        for (TiXmlElement* propNode = propertiesNode->FirstChildElement("Property")
            ; propNode != 0
            ; propNode = propNode->NextSiblingElement("Property"))
        {
            if(TiXmlElement* metaDataNode = propNode->FirstChildElement("MetaData")) {
                for (TiXmlElement* itemNode = metaDataNode->FirstChildElement("MetaItem")
                    ; itemNode != 0
                    ; itemNode = itemNode->NextSiblingElement("MetaItem"))
                {
                    const std::string* name = itemNode->Attribute(std::string("name"));
                    if (name && *name == "EditorWindow")
                        propNode->SetAttribute("name", "ToolWindow");
                }
            }
        }
    }
}

void NetworkConverter17to18::convert(TiXmlElement* networkElem) {
    TiXmlElement* processorsNode = networkElem->FirstChildElement("Processors");
    if (!processorsNode) {
        LWARNING("No <Processors> node found in <ProcessorNetwork>");
        return;
    }
    for (TiXmlElement* node = processorsNode->FirstChildElement("Processor"); node; node = node->NextSiblingElement("Processor")) {
        const std::string* type = node->Attribute(std::string("type"));
        if (type && *type == "StreamlineCreator") {
            TiXmlElement* propertiesNode = node->FirstChildElement("Properties");
            if (propertiesNode) {
                for (TiXmlElement* propertyNode = propertiesNode->FirstChildElement("Property"); propertyNode; propertyNode = propertyNode->NextSiblingElement("Property")) {
                    const std::string* propName = propertyNode->Attribute(std::string("name"));
                    if (!propName)
                        continue;
                    // Convert Vec2Property to IntervalProperty, therefore set appropriate attributes.
                    if (*propName == "absoluteMagnitudeThreshold") {
                        propertyNode->SetAttribute("minValue", "0");
                        propertyNode->SetAttribute("maxValue", "99999.99");
                        propertyNode->SetAttribute("minRange", "0");
                        propertyNode->SetAttribute("maxRange", "214748368");
                        propertyNode->SetAttribute("stepping", "0.01");
                    }
                    else if (*propName == "relativeMagnitudeThreshold") {
                        propertyNode->SetAttribute("minValue", "0");
                        propertyNode->SetAttribute("maxValue", "100");
                        propertyNode->SetAttribute("minRange", "0");
                        propertyNode->SetAttribute("maxRange", "214748368");
                        propertyNode->SetAttribute("stepping", "0.01");
                    }
                    else if (*propName == "streamlineLengthThresholdProp") {
                        propertyNode->SetAttribute("minValue", "2");
                        propertyNode->SetAttribute("maxValue", "1000");
                        propertyNode->SetAttribute("minRange", "0");
                        propertyNode->SetAttribute("maxRange", "214748368");
                        propertyNode->SetAttribute("stepping", "1");
                    }
                    // Set noise to 1%
                    else if (*propName == "minNumStreamlinesPerBundle") {
                        propertyNode->SetAttribute("value", "1");
                    }
                }
            }
        }
    }
}

void NetworkConverter18to19::convert(TiXmlElement* networkElem) {
    TiXmlElement* processorsNode = networkElem->FirstChildElement("Processors");
    if (!processorsNode) {
        LWARNING("No <Processors> node found in <ProcessorNetwork>");
        return;
    }

    std::string classificationName = "classification";

    std::set<std::string> raycasters;
    raycasters.insert("CPURaycaster");
    raycasters.insert("DynamicGLSLProcessor");
    raycasters.insert("IDRaycaster");
    raycasters.insert("IlluminationLineRaycaster");
    raycasters.insert("MultiVolumeRaycaster");
    raycasters.insert("RGBRaycaster");
    raycasters.insert("SegmentationRaycaster");
    raycasters.insert("SIMDRayCaster");
    raycasters.insert("SimpleRaycaster");
    raycasters.insert("SingleVolumeRaycaster");

    for (TiXmlElement* node = processorsNode->FirstChildElement("Processor"); node; node = node->NextSiblingElement("Processor")) {
        const std::string* type = node->Attribute(std::string("type"));
        if (type && raycasters.find(*type) != raycasters.end()) {
            TiXmlElement* propertiesNode = node->FirstChildElement("Properties");
            if (propertiesNode) {
                for (TiXmlElement* propertyNode = propertiesNode->FirstChildElement("Property"); propertyNode; propertyNode = propertyNode->NextSiblingElement("Property")) {
                    const std::string* propName = propertyNode->Attribute(std::string("name"));
                    if (!propName)
                        continue;

                    // Remove 'none' classification.
                    if (propName->find(classificationName) == 0) {
                        std::string suffix = propName->substr(classificationName.size());
                        // There might exist some old MultiVolumeRaycaster Workspaces
                        // which have serialized their Options (2, 3, 4) as well.
                        if (suffix == "" || suffix == "2" || suffix == "3" || suffix == "4") {
                            while (TiXmlElement* optionNode = propertyNode->FirstChildElement("Options")) {
                                propertyNode->RemoveChild(optionNode);
                            }
                        }
                    }
                }
            }
        }
    }
}

void NetworkConverter19to20::convert(TiXmlElement* networkElem) {
    TiXmlElement* processorsNode = networkElem->FirstChildElement("Processors");
    if (!processorsNode) {
        LWARNING("No <Processors> node found in <ProcessorNetwork>");
        return;
    }

    for (TiXmlElement* node = processorsNode->FirstChildElement("Processor"); node; node = node->NextSiblingElement("Processor")) {
        const std::string* type = node->Attribute(std::string("type"));
        if (type && *type == "OctreeCreator") {
            TiXmlElement* propertiesNode = node->FirstChildElement("Properties");
            if (propertiesNode) {
                std::vector<TiXmlElement*> deletedProperties;
                for (TiXmlElement* propertyNode = propertiesNode->FirstChildElement("Property"); propertyNode; propertyNode = propertyNode->NextSiblingElement("Property")) {
                    const std::string* propName = propertyNode->Attribute(std::string("name"));
                    if (!propName)
                        continue;

                    // AsyncComputProcessor brings it's own equivalents of the following properties.
                    if (*propName == "generateOctree" || *propName == "statusProperty" || *propName == "progressProperty" || *propName == "autogenerateOctree") {
                        deletedProperties.push_back(propertyNode);
                    }
                }
                // Delete selected nodes.
                for (TiXmlElement* node : deletedProperties) {
                    propertiesNode->RemoveChild(node);
                }
            }
        }
    }
}

} // namespace
