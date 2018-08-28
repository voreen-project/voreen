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

#ifndef VRN_ILLUMINATIONLINERAYCASTER_H
#define VRN_ILLUMINATIONLINERAYCASTER_H

#include "voreen/core/processors/volumeraycaster.h"
#include "voreen/core/datastructures/geometry/meshlistgeometry.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"

#define VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
//#define ERT_DOWNSAMPLE_SCHEME

namespace voreen {

/**
 * Performs a simple single pass raycasting using synchronized lines to generate entry/exit points, which enables the possibility of external lighting.
 */
class IlluminationLineRaycaster : public VolumeRaycaster {
public:
    IlluminationLineRaycaster();
    ~IlluminationLineRaycaster();
    virtual Processor* create() const;

    virtual std::string getCategory() const   { return "Raycasting";      }
    virtual std::string getClassName() const  { return "IlluminationLineRaycaster"; }

    //Specify member function pointer
#ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
    typedef bool (IlluminationLineRaycaster::*PostLineFunction)(int*);
#else
    typedef void (IlluminationLineRaycaster::*PostLineFunction)();
#endif

protected:
    virtual void setDescriptions() {
        setDescription("Performs a simple single pass raycasting using synchronized lines to generate entry/exit points, which enables the possibility of external lighting.");
    }

    virtual void beforeProcess();

#ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
    virtual void copyRenderPort(RenderPort *src, RenderPort *dst);
    virtual void visualizeProxyGeometrySize(RenderPort *dst, RenderPort *src, GLint src_unit, tgt::vec4 size);
    virtual void drawSquareWithLines(tgt::ivec2 v1, tgt::ivec2 v2, tgt::ivec2 v3, tgt::ivec2 v4);
    virtual bool drawSynchronizedLine(tgt::vec2 *startVertex, tgt::vec2 *endVertex, tgt::vec2 *startTexCoord, tgt::vec2 *endTexCoord, int *numLines);
    virtual bool postProcessLine(int *processed_lines);
    virtual bool postProcessLineGL4(int *processed_lines);
#else
    virtual void drawSynchronizedLine(tgt::vec2 *startVertex, tgt::vec2 *endVertex, tgt::vec2 *startTexCoord, tgt::vec2 *endTexCoord);
    virtual void postProcessLine();
    virtual void postProcessLineGL4();
#endif

    virtual void clearLightMaps();
    virtual void drawQuadIntoList();
    virtual bool drawLines(tgt::ivec4 coords, tgt::ivec2 light_coords, int lineWidth = 1);
    virtual bool drawLinesVBO(tgt::ivec4 coords, tgt::ivec2 light_coords, int lineWidth = 1);
    virtual bool drawZones(tgt::ivec4 coords, tgt::ivec2 light_coords, int lineWidth = 1);

    virtual tgt::vec4 calculateProxyGeometryScreenSize(MeshListGeometry *geom = 0, const tgt::mat4 *transform = 0, const tgt::vec2 pos = tgt::vec2(0.f,0.f), bool *isInside = 0);
    virtual void calculateERTmap(const GLint entryUnit, const GLint exitUnit);
    virtual void calculatePropagationZones(const tgt::vec2 light_pos, const GLint entryUnit, const GLint exitUnit);

    virtual void process();

    virtual void showAndHideProperties();
    virtual void screenSizeChanged();
    virtual void updatePropagation();

    virtual void initBufferObjects();

    virtual void initialize();
    virtual void deinitialize();

    virtual std::string generateHeader(const tgt::GpuCapabilities::GlVersion* version = 0);
    virtual void updateHeader();

private:
    //VolumePort volumePort_;
    GeometryPort geometryPort_;
    //RenderPort entryPort_;
    //RenderPort exitPort_;
#ifdef ERT_DOWNSAMPLE_SCHEME
    RenderPort ertPort_;
#endif
    //RenderPort outport_;
    RenderPort lightportOne_;
    RenderPort lightportTwo_;
    RenderPort lightportOnePos_;
    RenderPort lightportTwoPos_;
    RenderPort propagationZonesPort_;
    RenderPort propagationErtPort_;

    BoolProperty renderQuad_;
    BoolProperty autoLightCacheSize_;
    IntProperty     lightCacheWidth_;
    IntProperty     lightCacheHeight_;
    FloatProperty minLightDirLength_;
    FloatProperty diffCamLightDistance_;
    BoolProperty adaptiveNearPlane_;
    BoolProperty ertMod_;
#ifdef ERT_DOWNSAMPLE_SCHEME
    IntProperty     ertIncrement_;
    IntProperty     ertTextureSizeWidth_;
    IntProperty     ertTextureSizeHeight_;
#endif
    BoolProperty twoPassRenderingOn_;
    BoolProperty zoneRenderingOn_;
    OptionProperty<int> zoneOnOff_;
    IntProperty     zoneEdgeBlendingFactor_;
    BoolProperty accumulatedRayVis_;
    BoolProperty colorContribution_;
    FloatProperty shadowBrightness_;
    OptionProperty<int> interpolation_;
    BoolProperty piecewise_;
    BoolProperty adjustToProxyGeometry_;
#ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
    BoolProperty vboOn_;
    BoolProperty drawTrianglesAsLines_;
    BoolProperty linesOnlyAboveProxyGeometry_;
    BoolProperty visualizeProxySize_;
    BoolProperty pingPong_;
    IntProperty outputLineToDebugRT_;
#endif

    IntProperty lineWidth_;
    OptionProperty<tgt::GpuCapabilities::GlVersion> glslVersion_;

    tgt::Shader* raycastPrg_;         ///< The shader program used by this raycaster.
    tgt::Shader* copyImagePrg_;       ///< The shader program used to copy image.
    //tgt::Shader* blendPrg_;              ///< The shader program used to blend one rendering with a previous one. Currently not used
    tgt::Shader* propZonePrg_;        ///< The shader program used to define propagation zones.
    tgt::Shader* propErtPrg_;         ///< The shader program used when sweeping over the generated ERT texture.

    TransFunc1DKeysProperty transferFunc_;  ///< the property that controls the transfer-function

    // properties used for shading
    TransFunc1DKeysProperty albedo_;
    TransFunc1DKeysProperty surfaceness_;
    FloatProperty ambientIntensity_;

    PostLineFunction plfFunc_;          ///< function pointer to the function used after proccessing a synced line

    tgt::Camera lightCamera_;
    CameraProperty camera_;

    GLuint vboID_[2];                  ///< VBO for the vertices and texture coordinates, x+[0], y+[1]
    int currentVBO_;
    GLuint iboID_;                    ///< IBO for the indices, 0-MAX_DIM*2

    GLuint quadList_;                  ///< Geometry display list for a screen size quad.

    int lightMapUnits_[2];
    int lightMapPosUnits_[2];
    //float lightFrustumOffset_; //Currently not used

    bool glsl4On_;
    bool noLinesOverride_;
    bool verticalLines_;
    bool negativeDir_;

    tgt::ivec2 screenSize_;
    tgt::ivec2 propigationAxis_;
    tgt::ivec2 readWriteValues_;
    tgt::vec4  projLightPos_;
    tgt::vec4  projLightDirCam_;

    MeshListGeometry geometry_;

#ifndef VRN_ILLUMINATIONLINERAYCASTER_OPTIMIZED
    float triangleThickness_;
#endif
};


} // namespace

#endif // VRN_ILLUMINATIONLINERAYCASTER_H
