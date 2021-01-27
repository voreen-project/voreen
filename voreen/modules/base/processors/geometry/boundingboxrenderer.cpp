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

#include "boundingboxrenderer.h"

#include "voreen/core/voreenapplication.h"

#include "tgt/glmath.h"
#include "voreen/core/utils/stringutils.h"

#include "tgt/immediatemode/immediatemode.h"
#include "tgt/shadermanager.h"

namespace voreen {

using tgt::vec4;
using tgt::vec3;

const std::string BoundingBoxRenderer::loggerCat_("voreen.base.BoundingBoxRenderer");

BoundingBoxRenderer::BoundingBoxRenderer()
    : GeometryRendererBase()
    , volumeInport_(Port::INPORT, "volume", "Volume Input")
    , geometryInport_(Port::INPORT, "geometry", "Geometry Input")
    , enable_("enable", "Enable", true)
    , width_("boundingBoxWidth", "Line width", 1.0f, 1.0f, 10.0f)
    , bboxColor_("boundingboxColorSwitch", "Color", tgt::vec4(0.8f, 0.8f, 0.8f, 1.0f))
    , showGrid_("boundingboxGridShow", "Show Grid", false)
    , tilesSections_("boundingboxGridSize", "Min. Number of Grid Elements based on the longest side", 10, 2, 255)
    , tilesOpacity_("boundingboxGridOpacity", "Opacity Factor", 1.f, 0.f, 1.f)
    , showLegend_("showLegend", "Show Distance Legend", false)
    , legendScaleFactor_("legendScaleFactor", "Scale for legend", 1.0f, 1.0f, 4.0f)
    , fontProp_("fontProp", "Legend Font")
    , alignment_("layout", "Legend position:")
    , tileLength_(1.0f)
    , maxBoxLength_(1.0f)
    , legendDisplacement_("legendDisplacement", "Offset of legend: ", tgt::ivec2::zero, tgt::ivec2(-1000), tgt::ivec2(1000))
{
    //ports
    addPort(volumeInport_);
        volumeInport_.onChange(MemberFunctionCallback<BoundingBoxRenderer>(this, &BoundingBoxRenderer::onInputDataChange));
    addPort(geometryInport_);
        geometryInport_.onChange(MemberFunctionCallback<BoundingBoxRenderer>(this, &BoundingBoxRenderer::onInputDataChange));
    //enable
    addProperty(enable_);
    //line
    addProperty(width_);
        width_.setGroupID("line");
    addProperty(bboxColor_);
        bboxColor_.setGroupID("line");
    setPropertyGroupGuiName("line","Line Settings");
    //grid
    addProperty(showGrid_);
        showGrid_.setGroupID("grid");
    addProperty(tilesSections_);
        tilesSections_.onChange(MemberFunctionCallback<BoundingBoxRenderer>(this, &BoundingBoxRenderer::onTilesSectionChange));
        tilesSections_.setGroupID("grid");
    addProperty(tilesOpacity_);
        tilesOpacity_.setGroupID("grid");
    setPropertyGroupGuiName("grid","Grid Settings");
    //legend
    addProperty(showLegend_);
        showLegend_.setGroupID("legend");
    addProperty(legendScaleFactor_);
        legendScaleFactor_.setGroupID("legend");
   //addProperty(fontProp_);
        //fontProp_.setGroupID("legend");
        fontProp_.get()->setFontName(VoreenApplication::app()->getFontPath("Vera.ttf"));
        fontProp_.get()->setFontSize(12);
        fontProp_.get()->setTextAlignment(tgt::Font::TopCentered);

    addProperty(alignment_);
        alignment_.addOption("N",  "North", OPTION_N);
        alignment_.addOption("NE", "North-East", OPTION_NE);
        alignment_.addOption("E", "East", OPTION_E);
        alignment_.addOption("SE", "South-East", OPTION_SE);
        alignment_.addOption("S", "South", OPTION_S);
        alignment_.addOption("SW", "South-West", OPTION_SW);
        alignment_.addOption("W", "West", OPTION_W);
        alignment_.addOption("NW", "North-West", OPTION_NW);
        alignment_.addOption("CENTER", "Center", OPTION_CENTER);
        alignment_.selectByValue(OPTION_SE);
        alignment_.setGroupID("legend");
        
    addProperty(legendDisplacement_);
        legendDisplacement_.setGroupID("legend");
    setPropertyGroupGuiName("legend","Legend Settings");
}

Processor* BoundingBoxRenderer::create() const {
    return new BoundingBoxRenderer();
}

bool BoundingBoxRenderer::isReady() const {
    if (!isInitialized())
        return false;

    return (volumeInport_.isReady() || geometryInport_.isReady());
}

tgt::Bounds BoundingBoxRenderer::getBoundingBox() const {
    return boundingBox_.transform(inputToWorldTransformation_);
}

void BoundingBoxRenderer::initialize() {
    GeometryRendererBase::initialize();
    glGenQueries(6, occlusionQueries_);
}

void BoundingBoxRenderer::deinitialize() {
    glDeleteQueries(6, occlusionQueries_);
    GeometryRendererBase::deinitialize();
}

void BoundingBoxRenderer::render() {
    if (!enable_.get())
        return;

    // if both volume and geometry inports are ready, the volume has priority
    /*if (volumeInport_.isReady() && geometryInport_.isReady()) {
        LWARNING("Either volume inport or geometry inport may be connected");
        return;
    }*/

    // Extract corners.
    tgt::vec3 geomLLF = boundingBox_.getLLF();
    tgt::vec3 geomURB = boundingBox_.getURB();

    // Update transformation.
    MatStack.pushMatrix();
    MatStack.multMatrix(inputToWorldTransformation_);

    // render the outer bounding box
    IMode.color(bboxColor_.get());
    glLineWidth(width_.get());

    IMode.begin(tgt::ImmediateMode::LINES);

    // back face
    IMode.vertex(tgt::vec3(geomLLF[0], geomURB[1], geomURB[2]));
    IMode.vertex(tgt::vec3(geomLLF[0], geomLLF[1], geomURB[2]));

    IMode.vertex(tgt::vec3(geomLLF[0], geomLLF[1], geomURB[2]));
    IMode.vertex(tgt::vec3(geomURB[0], geomLLF[1], geomURB[2]));

    IMode.vertex(tgt::vec3(geomURB[0], geomLLF[1], geomURB[2]));
    IMode.vertex(tgt::vec3(geomURB[0], geomURB[1], geomURB[2]));

    IMode.vertex(tgt::vec3(geomURB[0], geomURB[1], geomURB[2]));
    IMode.vertex(tgt::vec3(geomLLF[0], geomURB[1], geomURB[2]));

    // front face
    IMode.vertex(tgt::vec3(geomLLF[0], geomURB[1], geomLLF[2]));
    IMode.vertex(tgt::vec3(geomLLF[0], geomLLF[1], geomLLF[2]));

    IMode.vertex(tgt::vec3(geomLLF[0], geomLLF[1], geomLLF[2]));
    IMode.vertex(tgt::vec3(geomURB[0], geomLLF[1], geomLLF[2]));

    IMode.vertex(tgt::vec3(geomURB[0], geomLLF[1], geomLLF[2]));
    IMode.vertex(tgt::vec3(geomURB[0], geomURB[1], geomLLF[2]));

    IMode.vertex(tgt::vec3(geomURB[0], geomURB[1], geomLLF[2]));
    IMode.vertex(tgt::vec3(geomLLF[0], geomURB[1], geomLLF[2]));

    // rest of the box
    IMode.vertex(tgt::vec3(geomLLF[0], geomURB[1], geomURB[2]));
    IMode.vertex(tgt::vec3(geomLLF[0], geomURB[1], geomLLF[2]));

    IMode.vertex(tgt::vec3(geomLLF[0], geomLLF[1], geomURB[2]));
    IMode.vertex(tgt::vec3(geomLLF[0], geomLLF[1], geomLLF[2]));

    IMode.vertex(tgt::vec3(geomURB[0], geomLLF[1], geomURB[2]));
    IMode.vertex(tgt::vec3(geomURB[0], geomLLF[1], geomLLF[2]));

    IMode.vertex(tgt::vec3(geomURB[0], geomURB[1], geomURB[2]));
    IMode.vertex(tgt::vec3(geomURB[0], geomURB[1], geomLLF[2]));

    IMode.end();

    if (showGrid_.get()) {

        IMode.color(tgt::vec4(bboxColor_.get().r, bboxColor_.get().g, bboxColor_.get().b, bboxColor_.get().a*tilesOpacity_.get()));

        // set opengl state for occlusion query
        glEnable(GL_CULL_FACE);
        glCullFace(GL_FRONT);
        glDepthMask(GL_FALSE);
        glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);

        // execute occlusion queries for each face
        glBeginQuery(GL_ANY_SAMPLES_PASSED, occlusionQueries_[0]);
        IMode.begin(tgt::ImmediateMode::QUADS);
            IMode.vertex(geomLLF);
            IMode.vertex(tgt::vec3(geomLLF[0], geomURB[1], geomLLF[2]));
            IMode.vertex(tgt::vec3(geomURB[0], geomURB[1], geomLLF[2]));
            IMode.vertex(tgt::vec3(geomURB[0], geomLLF[1], geomLLF[2]));
        IMode.end();
        glEndQuery(GL_ANY_SAMPLES_PASSED);

        glBeginQuery(GL_ANY_SAMPLES_PASSED, occlusionQueries_[1]);
        IMode.begin(tgt::ImmediateMode::QUADS);
            IMode.vertex(tgt::vec3(geomLLF[0], geomLLF[1], geomURB[2]));
            IMode.vertex(tgt::vec3(geomURB[0], geomLLF[1], geomURB[2]));
            IMode.vertex(tgt::vec3(geomURB[0], geomURB[1], geomURB[2]));
            IMode.vertex(tgt::vec3(geomLLF[0], geomURB[1], geomURB[2]));
        IMode.end();
        glEndQuery(GL_ANY_SAMPLES_PASSED);

        glBeginQuery(GL_ANY_SAMPLES_PASSED, occlusionQueries_[2]);
        IMode.begin(tgt::ImmediateMode::QUADS);
            IMode.vertex(tgt::vec3(geomLLF[0], geomLLF[1], geomLLF[2]));
            IMode.vertex(tgt::vec3(geomLLF[0], geomLLF[1], geomURB[2]));
            IMode.vertex(tgt::vec3(geomLLF[0], geomURB[1], geomURB[2]));
            IMode.vertex(tgt::vec3(geomLLF[0], geomURB[1], geomLLF[2]));
        IMode.end();
        glEndQuery(GL_ANY_SAMPLES_PASSED);

        glBeginQuery(GL_ANY_SAMPLES_PASSED, occlusionQueries_[3]);
        IMode.begin(tgt::ImmediateMode::QUADS);
            IMode.vertex(tgt::vec3(geomURB[0], geomLLF[1], geomURB[2]));
            IMode.vertex(tgt::vec3(geomURB[0], geomLLF[1], geomLLF[2]));
            IMode.vertex(tgt::vec3(geomURB[0], geomURB[1], geomLLF[2]));
            IMode.vertex(tgt::vec3(geomURB[0], geomURB[1], geomURB[2]));
        IMode.end();
        glEndQuery(GL_ANY_SAMPLES_PASSED);

        glBeginQuery(GL_ANY_SAMPLES_PASSED, occlusionQueries_[4]);
        IMode.begin(tgt::ImmediateMode::QUADS);
            IMode.vertex(tgt::vec3(geomLLF[0], geomLLF[1], geomLLF[2]));
            IMode.vertex(tgt::vec3(geomURB[0], geomLLF[1], geomLLF[2]));
            IMode.vertex(tgt::vec3(geomURB[0], geomLLF[1], geomURB[2]));
            IMode.vertex(tgt::vec3(geomLLF[0], geomLLF[1], geomURB[2]));
        IMode.end();
        glEndQuery(GL_ANY_SAMPLES_PASSED);

        glBeginQuery(GL_ANY_SAMPLES_PASSED, occlusionQueries_[5]);
        IMode.begin(tgt::ImmediateMode::QUADS);
            IMode.vertex(tgt::vec3(geomLLF[0], geomURB[1], geomLLF[2]));
            IMode.vertex(tgt::vec3(geomLLF[0], geomURB[1], geomURB[2]));
            IMode.vertex(tgt::vec3(geomURB[0], geomURB[1], geomURB[2]));
            IMode.vertex(tgt::vec3(geomURB[0], geomURB[1], geomLLF[2]));
        IMode.end();
        glEndQuery(GL_ANY_SAMPLES_PASSED);

        // set back opengl state
        glDepthMask(GL_TRUE);
        glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
        glDisable(GL_CULL_FACE);
        glCullFace(GL_BACK);

        // retrieve the result for each occlusion query
        int renderFaces[6];

        // get result of each query
        for (size_t i = 0; i < 6; ++i) {
            int queryReady = 0;
            while(!queryReady)
                glGetQueryObjectiv(occlusionQueries_[i], GL_QUERY_RESULT_AVAILABLE, &queryReady);
           glGetQueryObjectiv(occlusionQueries_[i], GL_QUERY_RESULT, &renderFaces[i]);
        }


        glEnable(GL_BLEND);
        // render the tiles for each face as lines, determine visibility using the occlusion query results
        IMode.begin(tgt::ImmediateMode::LINES);

        // render lines for creating front and back face tiles
        tgt::vec3 tileDim((geomURB[0] - geomLLF[0]), (geomURB[1] - geomLLF[1]), (geomURB[2] - geomLLF[2]));

        // start with vertical lines
        for (float x = tileLength_; x < tileDim.x; x += tileLength_) {
            if (renderFaces[0]) {
                IMode.vertex(tgt::vec3(geomLLF[0] + x, geomLLF[1], geomLLF[2]));
                IMode.vertex(tgt::vec3(geomLLF[0] + x, geomURB[1], geomLLF[2]));
            }
            if (renderFaces[1]) {
                IMode.vertex(tgt::vec3(geomLLF[0] + x, geomLLF[1], geomURB[2]));
                IMode.vertex(tgt::vec3(geomLLF[0] + x, geomURB[1], geomURB[2]));
            }
        }
        // horizontal lines
        for (float y = tileLength_; y < tileDim.y; y += tileLength_) {
            if (renderFaces[0]) {
                IMode.vertex(tgt::vec3(geomLLF[0], geomLLF[1] + y, geomLLF[2]));
                IMode.vertex(tgt::vec3(geomURB[0], geomLLF[1] + y, geomLLF[2]));
            }
            if (renderFaces[1]) {
                IMode.vertex(tgt::vec3(geomLLF[0], geomLLF[1] + y, geomURB[2]));
                IMode.vertex(tgt::vec3(geomURB[0], geomLLF[1] + y, geomURB[2]));
            }
        }

        // top and bottom face tiles
        for (float x = tileLength_; x < tileDim.x; x += tileLength_) {
            if (renderFaces[4]) {
                IMode.vertex(tgt::vec3(geomLLF[0] + x, geomLLF[1], geomLLF[2]));
                IMode.vertex(tgt::vec3(geomLLF[0] + x, geomLLF[1], geomURB[2]));
            }
            if (renderFaces[5]) {
                IMode.vertex(tgt::vec3(geomLLF[0] + x, geomURB[1], geomLLF[2]));
                IMode.vertex(tgt::vec3(geomLLF[0] + x, geomURB[1], geomURB[2]));
            }
        }
        for (float z = tileLength_; z < tileDim.z; z += tileLength_) {
            if (renderFaces[4]) {
                IMode.vertex(tgt::vec3(geomLLF[0], geomLLF[1], geomLLF[2] + z));
                IMode.vertex(tgt::vec3(geomURB[0], geomLLF[1], geomLLF[2] + z));
            }
            if (renderFaces[5]) {
                IMode.vertex(tgt::vec3(geomLLF[0], geomURB[1], geomLLF[2] + z));
                IMode.vertex(tgt::vec3(geomURB[0], geomURB[1], geomLLF[2] + z));
            }
        }

        // left and right face tiles
        for (float y = tileLength_; y < tileDim.y; y += tileLength_) {
            if (renderFaces[2]) {
                IMode.vertex(tgt::vec3(geomLLF[0], geomLLF[1] + y, geomLLF[2]));
                IMode.vertex(tgt::vec3(geomLLF[0], geomLLF[1] + y, geomURB[2]));
            }
            if (renderFaces[3]) {
                IMode.vertex(tgt::vec3(geomURB[0], geomLLF[1] + y, geomLLF[2]));
                IMode.vertex(tgt::vec3(geomURB[0], geomLLF[1] + y, geomURB[2]));
            }
        }
        for (float z = tileLength_; z < tileDim.z; z += tileLength_) {
            if (renderFaces[2]) {
                IMode.vertex(tgt::vec3(geomLLF[0], geomLLF[1], geomLLF[2] + z));
                IMode.vertex(tgt::vec3(geomLLF[0], geomURB[1], geomLLF[2] + z));
            }
            if (renderFaces[3]) {
                IMode.vertex(tgt::vec3(geomURB[0], geomLLF[1], geomLLF[2] + z));
                IMode.vertex(tgt::vec3(geomURB[0], geomURB[1], geomLLF[2] + z));
            }
        }

        IMode.end();
        glDisable(GL_BLEND);
    }

    if(showLegend_.get()) {
        renderLegend();
    }

    glLineWidth(1.f);

    IMode.color(tgt::vec4::one);
    MatStack.popMatrix();

    LGL_ERROR;
}

void BoundingBoxRenderer::renderLegend()
{
    float scale = legendScaleFactor_.get();

    //paint legend
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.pushMatrix();
    MatStack.loadIdentity();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();

    MatStack.loadMatrix(tgt::mat4::createOrtho(0, static_cast<float>(viewport_.x), 0, static_cast<float>(viewport_.y), -1, 1));

    float horizontalOffset = 12.f;
    float verticalOffset = 12.f;

    fontProp_.get()->setFontSize(static_cast<int>(12*scale));

    float lineWidth  = 0;
    if(showGrid_.get()) {
        lineWidth = std::min(scale*80.f,(float)viewport_.x);
    }
    else {
        lineWidth = std::min(scale*60.f,(float)viewport_.x);
    }
    float lineHeight = fontProp_.get()->getLineHeight();
    float height = scale*34+lineHeight;
    float width = lineWidth;


    float bottomOffset = static_cast<float>(legendDisplacement_.get().x);
    float leftOffset = static_cast<float>(legendDisplacement_.get().y);

    switch (alignment_.getValue()) {
    case OPTION_N:
        bottomOffset += viewport_.y - height - verticalOffset;
        leftOffset += viewport_.x / 2 - width / 2;
        break;
    case OPTION_NW:
        bottomOffset += viewport_.y - height - verticalOffset;
        leftOffset += horizontalOffset;
        break;
    case OPTION_NE:
        bottomOffset += viewport_.y - height - verticalOffset;
        leftOffset += viewport_.x - width - horizontalOffset;
        break;
    case OPTION_CENTER:
        bottomOffset += viewport_.y / 2 - height / 2;
        leftOffset += viewport_.x / 2 - width / 2;
        break;
    case OPTION_W:
        bottomOffset += viewport_.y / 2 - height / 2;
        leftOffset += horizontalOffset;
        break;
    case OPTION_E:
        bottomOffset += viewport_.y / 2 - height / 2;
        leftOffset += viewport_.x - width - horizontalOffset;
        break;
    case OPTION_S:
        bottomOffset += verticalOffset;
        leftOffset += viewport_.x / 2 - width / 2;
        break;
    case OPTION_SW:
        bottomOffset += verticalOffset;
        leftOffset += horizontalOffset;
        break;
    case OPTION_SE:
        bottomOffset += verticalOffset;
        leftOffset += viewport_.x - width - horizontalOffset;
        break;
    default:
        break;
    }

    IMode.color(bboxColor_.get());
    fontProp_.get()->setLineWidth(lineWidth);
    fontProp_.get()->setFontColor(bboxColor_.get());

    if(showGrid_.get()) {
        fontProp_.get()->render(tgt::vec3(leftOffset - 5.f, bottomOffset, 0), formatSpatialLength(tileLength_), viewport_);
        IMode.begin(tgt::ImmediateMode::LINE_LOOP);
        IMode.vertex(tgt::vec2(leftOffset+0.0f*lineWidth,bottomOffset+1.5f*lineHeight));
        IMode.vertex(tgt::vec2(leftOffset+0.8f*lineWidth,bottomOffset+1.5f*lineHeight));
        IMode.vertex(tgt::vec2(leftOffset+1.0f*lineWidth,bottomOffset+1.5f*lineHeight+scale*25.f));
        IMode.vertex(tgt::vec2(leftOffset+0.2f*lineWidth,bottomOffset+1.5f*lineHeight+scale*25.f));
        IMode.end();
    }
    else {
        fontProp_.get()->render(tgt::vec3(leftOffset, bottomOffset, 0), formatSpatialLength(maxBoxLength_), viewport_);
        IMode.begin(tgt::ImmediateMode::LINE_LOOP);
        IMode.vertex(tgt::vec2(leftOffset+0.2f*lineWidth,bottomOffset+1.5f*lineHeight));
        IMode.vertex(tgt::vec2(leftOffset+0.8f*lineWidth,bottomOffset+1.5f*lineHeight));
        IMode.vertex(tgt::vec2(leftOffset+0.8f*lineWidth,bottomOffset+1.5f*lineHeight+scale*25.f));
        IMode.vertex(tgt::vec2(leftOffset+0.2f*lineWidth,bottomOffset+1.5f*lineHeight+scale*25.f));
        IMode.end();
    }

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.popMatrix();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
}

//--------------------------------------------------------------------------------------
//      Callbacks
//--------------------------------------------------------------------------------------
void BoundingBoxRenderer::onInputDataChange() {

    if(!volumeInport_.isReady() && !geometryInport_.isReady())
        return;

    // Volume data has priority.
    if (volumeInport_.getData()) {
        boundingBox_ = volumeInport_.getData()->getBoundingBox(false).getBoundingBox(false);
        inputToWorldTransformation_ = volumeInport_.getData()->getPhysicalToWorldMatrix();
    }
    else if (geometryInport_.hasData()) {
        boundingBox_ = geometryInport_.getData()->getBoundingBox(false);
        inputToWorldTransformation_ = geometryInport_.getData()->getTransformationMatrix();
    }
    else {
        boundingBox_ = tgt::Bounds(); // undefined BoundingBox
        tgtAssert(false, "no input has data"); //< should never get here
    }

    // Once the input data changed, the tile-size needs to be recalculated as well.
    onTilesSectionChange();
}
void BoundingBoxRenderer::onTilesSectionChange() {
    // Extract corners.
    tgt::vec3 geomLLF = boundingBox_.getLLF();
    tgt::vec3 geomURB = boundingBox_.getURB();

    //compute spatial length
    tgt::vec3 lx = (inputToWorldTransformation_ * tgt::vec4(geomURB.x, geomLLF.y, geomLLF.z, 1.f)).xyz()
        - (inputToWorldTransformation_ * tgt::vec4(geomLLF, 1.f)).xyz();
    tgt::vec3 ly = (inputToWorldTransformation_ * tgt::vec4(geomLLF.x, geomURB.y, geomLLF.z, 1.f)).xyz()
        - (inputToWorldTransformation_ * tgt::vec4(geomLLF, 1.f)).xyz();
    tgt::vec3 lz = (inputToWorldTransformation_ * tgt::vec4(geomLLF.x, geomLLF.y, geomURB.z, 1.f)).xyz()
        - (inputToWorldTransformation_ * tgt::vec4(geomLLF, 1.f)).xyz();

    maxBoxLength_ = tgt::max(tgt::vec3(tgt::length(lx), tgt::length(ly), tgt::length(lz)));
    tileLength_ = maxBoxLength_ / tilesSections_.get();

    // Prevent the tileLength from being an inconvenient number,
    // therefore normalize according to how it's done in 'formatSpatialLength'.
    // Note that only one of the following loops is going to be executed.
    float power = 1.0f;
    float len = tileLength_;
    while (len < 1.0f) {
        power *= 1000.0f;
        len = tileLength_ * power;
    }
    while (len > 1000.0f) {
        power /= 1000.0f;
        len = tileLength_ * power;
    }

    // Cut of decimals. If this leads to undesired amount of grid cells,
    // cut of starting at the first decimal and so on, till the sweet spot has been found.
    float basePower = 1.0f / 10.0f; // Target multiples of 10.
    do {
        tileLength_ = std::floor(len * basePower) / (basePower * power);
        basePower *= 10.0f;
    } while (tileLength_ * (tilesSections_.get() + 1) <= maxBoxLength_);
}

}

