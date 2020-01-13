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

#include "volumelabeling.h"
#include "voreen/core/processors/volumeraycaster.h"

#include "tgt/types.h"
#include "tgt/tgt_gl.h"
#include "tgt/glmath.h"
#include "tgt/shadermanager.h"
#include "tgt/types.h"
#include "tgt/vector.h"
#include "tgt/gpucapabilities.h"

#include "voreen/core/utils/voreenpainter.h"
#include "voreen/core/utils/voreenqualitymode.h"
#include "voreen/core/voreenapplication.h"
#include "voreen/core/properties/cameraproperty.h"

#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <vector>
#include <algorithm>
#include <iostream>
#include <sstream>
#include <fstream>

namespace voreen {

using namespace tgt;
using namespace labeling;

using std::min;
using std::max;
using std::string;
using std::vector;
using std::swap;

const std::string VolumeLabeling::loggerCat_("voreen.volumelabeling.VolumeLabeling");

VolumeLabeling::VolumeLabeling()
    : ImageProcessor("volumelabeling")
    , labelingWidget_(0)
    , labelShader_(0)
    , showLabels_("showLabels", "Render Labels")
    , labelFile_("labelFile", "Segment Description File", "Path to Segment Description File",
        VoreenApplication::app()->getUserDataPath(""))
    , fontExtern_("fontExtern", "Font (extern)")
    , fontIntern_("fontIntern", "Font (intern)")
    , haloColorExtern_("haloColor", "Halo-Color", Color(0.8f, 0.8f, 0.8f, 1.0f))
    , lockInternalFontSettings_("lockInternalFontSettings", "Lock internal to external font settings", true)
    , haloColorIntern_("haloColorIntern", "Halo-Color (internal)", Color(0.8f, 0.8f, 0.8f, 1.0f))
    , shape3D_("shape3D", "3D shape fitting", true)
    , drawHalo_("drawHalo", "Draw Halo", true)
    , layout_("labelLayoutAsString", "Layout")
    , minSegmentSize_("minSegmentSize", "Minimal segment size", 10, 0, 500)
    , maxIterations_("maxIterations", "Iteration limit", 40, 1, 100)
    , filterKernel_("filterKernel", "Filter kernel")
    , filterDelta_("filterDelta", "Filter gap", 2, 1, 10)
    , distanceMapStep_("distanceMapStep", "Distance map step", 4, 1, 10)
    , polynomialDegree_("polynomialDegree", "Polynomial degree", 2, 1, 10)
    , bezierHorzDegree_("bezierHorzDegree", "Bezier horz. degree", 8, 1, 12)
    , bezierVertDegree_("bezierVertDegree", "Bezier vert. degree", 4, 1, 10)
    , valid_(false)
    , editMode_(false)
    , idMapPort_(Port::INPORT, "image.idmap", "ID Map Input")
    , volumePort_(Port::INPORT, "volumehandle.volumehandle", "Volume Input")
    , pickedLabel_(0)
    , drag_(false)
    , labelingPort_(Port::OUTPORT, "image.labeling", "Labeling Image Output", true)
    , camera_("camera", "Camera", tgt::Camera(vec3(0.f, 0.f, 3.5f), vec3(0.f, 0.f, 0.f), vec3(0.f, 1.f, 0.f)))
    , currentVolume_(0)

#ifdef labelDEBUG
    , showSegmentIDs_("showSegmentIDs", "Show segment IDs", false)
    , drawDebugMaps_("drawDebugMaps", "Draw debug maps", false)
    , drawConvexHull_("drawConvexHull", "Draw convex hull", false)
#endif

{


    // setup fonts
    fontExtern_.get()->setFontColor(tgt::vec4(0.0f, 0.0f, 1.0f, 1.0f));
    fontExtern_.get()->setTextAlignment(tgt::Font::TopLeft);
    fontExtern_.get()->setLineWidth(100000.0f);
    // "infinite" linewidth, because the Textrendering
    // does not work with linebreaks and so we can stop
    // the font rendering from breaking lines
    fontExtern_.invalidate();
    fontIntern_.get()->setFontColor(tgt::vec4(0.0f, 0.0f, 1.0f, 1.0f));
    fontIntern_.get()->setTextAlignment(tgt::Font::TopLeft);
    fontIntern_.get()->setLineWidth(100000.0f); // infinite linewidth. See above.
    fontIntern_.invalidate();


    lockInternalFontSettings_.onChange(MemberFunctionCallback<VolumeLabeling>(this, &VolumeLabeling::updateFontSettings));
    haloColorExtern_.onChange(MemberFunctionCallback<VolumeLabeling>(this, &VolumeLabeling::updateFontSettings));
    fontExtern_.onChange(MemberFunctionCallback<VolumeLabeling>(this, &VolumeLabeling::updateFontSettings));

    MemberFunctionCallback<VolumeLabeling> cmaUnlock(this, &VolumeLabeling::unlockInternalFontSettingsEvt);
    haloColorIntern_.onChange(cmaUnlock);
    fontIntern_.onChange(cmaUnlock);

    shape3D_.onChange(MemberFunctionCallback<VolumeLabeling>(this, &VolumeLabeling::setShape3DEvt));

    circle_ = generateCircle(10.f);
    image_.width = 0;
    image_.height = 0;

    showLabels_.addOption("all", "All");
    showLabels_.addOption("none", "None");
    showLabels_.addOption("external-only", "External only");
    showLabels_.addOption("internal-only", "Internal only");
    addProperty(showLabels_);

    labelFile_.onChange(MemberFunctionCallback<VolumeLabeling>(this, &VolumeLabeling::labelFileChanged));
    addProperty(labelFile_);

    layout_.addOption("silhouette", "Silhouette");
    layout_.addOption("left-right", "Left-Right");
    addProperty(layout_);

    addProperty(haloColorExtern_);
    addProperty(fontExtern_);
    addProperty( &haloColorIntern_ );
    addProperty(fontIntern_);
    addProperty(drawHalo_);

    createFilterKernels();
    for (unsigned int i=0; i<kernels_.size(); i++)
        filterKernel_.addOption(kernels_[i]->caption, kernels_[i]->caption, i);

    addProperty(filterKernel_);
    addProperty(filterDelta_);
    addProperty(distanceMapStep_);

    addProperty( &lockInternalFontSettings_ );

    addProperty( &shape3D_ );
    addProperty( &polynomialDegree_ );
    addProperty( &bezierHorzDegree_ );
    addProperty( &bezierVertDegree_ );

    addProperty(camera_);

#ifdef labelDEBUG

    segmentationTarget_ = -1;
    distanceMapTarget_ = -1;
    /* tc_->changeType(distanceMapTarget_+1, getTargetType(Identifier::Unknown));
    tc_->changeType(distanceMapTarget_+2, getTargetType(Identifier::Unknown));
    tc_->changeType(distanceMapTarget_+3, getTargetType(Identifier::Unknown));
    tc_->changeType(distanceMapTarget_+4, getTargetType(Identifier::Unknown));*/

    /*
    GroupProp* debugGroup = new GroupProp("debugGroup", "Debugging");
    debugGroup->addGroupedProp(&showSegmentIDs_);
    debugGroup->addGroupedProp(&drawConvexHull_);
    debugGroup->addGroupedProp(&drawDebugMaps_);
    addProperty(debugGroup);
    */
    addProperty(showSegmentIDs_);
    addProperty(drawConvexHull_);
    addProperty(drawDebugMaps_);

#endif

    addPort(volumePort_);
    addPort(idMapPort_);
    addPort(labelingPort_);
}

// destructor
VolumeLabeling::~VolumeLabeling() {
    while(!kernels_.empty()) {
        delete kernels_.back();
        kernels_.pop_back();
    }
}

void VolumeLabeling::initialize() {
    ImageProcessor::initialize();  //< loads the shader
    labelShader_ = program_;
    if(!labelShader_) {
        LERROR("Failed to load shaders!");
        processorState_ = PROCESSOR_STATE_NOT_INITIALIZED;
        throw VoreenException(getClassName() + ": Failed to load shader");
    }

    processorState_ = PROCESSOR_STATE_NOT_READY;

    labelFileChanged();
}

void VolumeLabeling::deinitialize() {
    for(size_t i = 0; i != kernels_.size(); i++){
        delete kernels_[i];
    }
    kernels_.resize(0);
    ImageProcessor::deinitialize();
}

// analyzes id-raycasting result and renders labels
void VolumeLabeling::process() {
    valid_ = false;

    if (!idMapPort_.isReady() || !labelingPort_.isReady() || !volumePort_.isReady())
        return;

    // if no user-interaction during the last 500ms and
    // not in dragging-mode, leave edit-mode
    if (editMode_ && !drag_ && time.getRuntime() > 500)
        editMode_ = false;

    // activate and clear render target
    labelingPort_.activateTarget();
    glPushAttrib(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    //glClearColor(backgroundColor_.get().r, backgroundColor_.get().g, backgroundColor_.get().b, 0.f);
    glClearDepth(0.f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glPopAttrib();

    if (volumePort_.hasChanged()) {
        currentVolume_ = volumePort_.getData();
    }

    // if not in coarseness-mode and label textures are available, perform labeling
    if (!QualityMode.isInteractionMode() && labelPersistentData_.size() > 0) {

        // if currently no user interaction, renew label layout

        if (!editMode_) {

            pickedLabel_ = NULL;
#ifdef labelDEBUG
            if (drawDebugMaps_.get())
                renderMaps();
#endif
        }

        deleteLabels();
        readImage();
        updateTextMetrics();
        findAnchors();
        labelLayout();

        renderLabels();

        valid_ = true;
    }

    labelingPort_.deactivateTarget();
    TextureUnit::setZeroUnit();
    glColor4f(1.f, 1.f, 1.f, 1.f);

    LGL_ERROR;

}

// ----------- labeling methods --------------

// performs the labeling
void VolumeLabeling::labelLayout() {
    if (layout_.get() == "silhouette")
        calcConvexHull();

    calcLabelPositions();
    appendUnsegmentedLabels();

    toWorld();
}



void VolumeLabeling::updateTextMetrics() {
    //FIXME(apv): this is useless
    for (size_t i = 0; i < labelPersistentData_.size(); i++) {
        std::string text = labelPersistentData_[i]->text.text;
        tgt::vec2 size = fontExtern_.get()->getSize(tgt::vec3::zero, text, idMapPort_.getSize());
        labelPersistentData_[i]->text.width  = static_cast<int>(size.x+0.5f);
        labelPersistentData_[i]->text.height = static_cast<int>(size.y+0.5f);
    }
}

namespace {
// used for creating distance maps: is called at every pixel position
void inline anchorAlgorithmCore(int* const &idBuffer, int* &distanceField,
                                int const &x, int const &y,
                                int const &imgWidth,
                                int &distance, int &newID, int &oldID, int &lastID,
                                int &lastDist, int &offset,
                                int const &step)
{
    newID = static_cast<int>(idBuffer[y * imgWidth + x]);
    if (newID == oldID && newID > 0) {
        distance += step;
        lastID = newID;
        lastDist = distance;
        offset = 0;
    }
    else if (newID == 0 && distance > 0) {
        offset += step;
        distance = 0; // modDist(lastDist, offset);
    }
    else {
        if (newID != lastID)
            distance = 0;
    }

    if (distanceField[y*imgWidth + x] == -1)
        distanceField[y*imgWidth + x] = distance;
    else
        distanceField[y*imgWidth + x] = min(distanceField[y*imgWidth + x], distance);

    oldID = newID;
}

} // anonymous

// - find anchor points in id-raycasting result
// - create a label for each anchor
// - calc bounding box of rendered object
void VolumeLabeling::findAnchors() {
    const int maxLabelID = 255;

    int* idBuffer = image_.idBuffer.getData();

    int* distanceField = new int[image_.width * image_.height];
    int* horzDistance = new int[image_.width*image_.height];
    int* vertDistance = new int[image_.width*image_.height];
    int* ascDiagDistance = new int[image_.width*image_.height];
    int* descDiagDistance = new int[image_.width*image_.height];

    // if pixel count exceeds limits, reduce number ob scanlines
    int step;
    if (image_.width*image_.height >= 2000000)
        step = 4;
    if (image_.width*image_.height >= 1000000)
        step = 3;
    else if (image_.width*image_.height >= 500000)
        step = 2;
    else
        step = 1;

    step = distanceMapStep_.get();

    // initialize distance maps
    for (int i=0; i<image_.width*image_.height; ++i) {
        horzDistance[i] = -1;
        vertDistance[i] = -1;
        ascDiagDistance[i] = -1;
        descDiagDistance[i] = -1;
    }

    //note(apv): i don't think this is a hash table
    // hash table for segments
    ivec2 anchors[maxLabelID+1];
    int segmentSize[maxLabelID+1];
    for (int i = 0; i <= maxLabelID; ++i) {
        anchors[i] = ivec2(-1,-1);
        segmentSize[i] = 0;
    }

    //
    // determine object's bounding box and silhouette points
    //
    // bottom-left and top-right of bounding box
    ivec2 bl = ivec2(image_.width,image_.height);
    ivec2 tr = ivec2(-1,-1);

    image_.silPoints.clear();
    // bottom: left to right
    for (int x = 0; x < image_.width-step; x += step) {
        int y = 0;
        while (y < (image_.height - step) && image_.idImage.getElem(x,y) == 0)
            y += step;

        if (image_.idImage.getElem(x,y) > 0) {
            image_.silPoints.push_back(ivec2(x,y));
            if (x < bl.x)
                bl.x = x;
            if (y < bl.y)
                bl.y = y;
            if (x > tr.x)
                tr.x = x;
            if (y > tr.y)
                tr.y = y;
        }
    }
    // right: bottom to top
    for (int y = bl.y; y < image_.height-step; y += step) {
        int x = tr.x;
        while (x > bl.x && image_.idImage.getElem(x,y) == 0)
            x -= step;

        if (image_.idImage.getElem(x,y) > 0) {
            image_.silPoints.push_back(ivec2(x,y));
            if (y > tr.y)
                tr.y = y;
        }
    }
    // top: right to left
    for (int x = tr.x; x > bl.x; x -= step) {
        int y = tr.y;
        while (y > bl.y && image_.idImage.getElem(x,y) == 0)
            y -= step;

        if (image_.idImage.getElem(x,y) > 0)
            image_.silPoints.push_back(ivec2(x,y));
    }
    // left: top to bottom
    for (int y = tr.y; y > bl.y; y -= step) {
        int x = bl.x;
        while (x < tr.x && image_.idImage.getElem(x,y) == 0)
            x += step;

        if (image_.idImage.getElem(x,y) > 0)
            image_.silPoints.push_back(ivec2(x,y));
    }

    //
    // creation of 4 direction-dependend distance maps:
    // horizontal, vertical, ascending diagonal, descending diagonal
    //

    int oldID;
    int newID;
    int lastID = -1;
    int distance;
    int lastDist;
    int offset;

    // find vertical distances
    // first pass (bottom to top)
    oldID = idBuffer[0];
    newID = idBuffer[0];
    for (int x = bl.x; x < tr.x; x += step) {
        distance = 0;
        for (int y = bl.y; y < tr.y; y += step) {
            anchorAlgorithmCore( idBuffer, vertDistance, x, y,
                image_.width, distance, newID, oldID, lastID,
                lastDist, offset, step);
        }
    }
    // second pass (top to bottom)
    oldID = 0;
    newID = 0;
    for (int x = bl.x; x < tr.x; x += step) {
        distance = 0;
        for (int y = tr.y; y > bl.y; y -= step) {
              anchorAlgorithmCore( idBuffer, vertDistance, x, y,
                image_.width, distance, newID, oldID, lastID,
                lastDist, offset, step);
        }
    }


    // find horizontal distances
    // first pass (left to right)
    oldID = 0;
    newID = 0;
    for (int y = bl.y; y < tr.y; y += step) {
        distance = 0;
        lastID = 0;
        lastDist = 0;
        offset = 0;
        oldID = 0;
        newID = 0;
        for (int x = bl.x; x < tr.x; x += step) {
            anchorAlgorithmCore( idBuffer, horzDistance, x, y,
                image_.width, distance, newID, oldID, lastID,
                lastDist, offset, step);
        }
    }
    // second pass (right to left)
    oldID = 0;
    newID = 0;
    for (int y = bl.y; y < tr.y; y += step) {
        distance = 0;
        lastID = 0;
        lastDist = 0;
        offset = 0;
        oldID = 0;
        newID = 0;
        for (int x = tr.x; x > bl.x; x -= step) {
            anchorAlgorithmCore( idBuffer, horzDistance, x, y,
                image_.width, distance, newID, oldID, lastID,
                lastDist, offset, step);
        }
    }

    //
    // diagonals
    //

    // descending diagonal distances

    // first pass (start from left)
    oldID = 0;
    newID = 0;
    for (int y_t = bl.y; y_t < tr.y; y_t += step) {
        int distance = 0;
        int y = y_t;
        for (int x = bl.x; x < tr.x && y > bl.y; ) {
            x += step;
            y -= step;
            anchorAlgorithmCore( idBuffer, descDiagDistance, x, y,
                image_.width, distance, newID, oldID, lastID,
                lastDist, offset, step);
        }
    }

    // second pass (start from top)
    oldID = 0;
    newID = 0;
    for (int x_t = bl.x; x_t <= tr.x; x_t += step) {
        int distance = 0;
        int x = x_t;
        for (int y = tr.y; y > bl.y && x < tr.x; ) {
            x += step;
            y -= step;
            anchorAlgorithmCore( idBuffer, descDiagDistance, x, y,
                image_.width, distance, newID, oldID, lastID,
                lastDist, offset, step);

        }
    }

    // third pass (start from right)
    oldID = 0;
    newID = 0;
    for (int y_t = bl.y; y_t <= tr.y; y_t += step) {
        int distance = 0;
        int y = y_t;
        for (int x = tr.x; x > bl.x && y < tr.y; ) {
            x -= step;
            y += step;
            anchorAlgorithmCore( idBuffer, descDiagDistance, x, y,
                image_.width, distance, newID, oldID, lastID,
                lastDist, offset, step);
        }
    }

    // fourth pass (start from bottom)
    oldID = 0;
    newID = 0;
    for (int x_t = bl.x; x_t < tr.x; x_t += step) {
        int distance = 0;
        int x = x_t;
        for (int y = bl.y; y < tr.y && x > bl.x; ) {
            x -= step;
            y += step;
            anchorAlgorithmCore( idBuffer, descDiagDistance, x, y,
                image_.width, distance, newID, oldID, lastID,
                lastDist, offset, step);
        }
    }


    //-----------------------------

    // ascending diagonal distances

    // first pass (start from left)
    oldID = 0;
    newID = 0;
    for (int y_t = bl.y; y_t < tr.y; y_t += step) {
        int distance = 0;
        int y = y_t;
        for (int x = bl.x; x < tr.x && y < tr.y; ) {
            x += step;
            y += step;
            anchorAlgorithmCore( idBuffer, ascDiagDistance, x, y,
                image_.width, distance, newID, oldID, lastID,
                lastDist, offset, step);
        }
    }

    // second pass (start from top)
    oldID = 0;
    newID = 0;
    for (int x_t = bl.x; x_t <= tr.x; x_t += step) {
        int distance = 0;
        int x = x_t;
        for (int y = tr.y; y > bl.y && x > bl.x; ) {
            x -= step;
            y -= step;
            anchorAlgorithmCore( idBuffer, ascDiagDistance, x, y,
                image_.width, distance, newID, oldID, lastID,
                lastDist, offset, step);

        }
    }

    // third pass (start from right)
    oldID = 0;
    newID = 0;
    for (int y_t = bl.y; y_t < tr.y; y_t += step) {
        int distance = 0;
        int y = y_t;
        for (int x = tr.x; x > bl.x && y > bl.y; ) {
            x -= step;
            y -= step;
            anchorAlgorithmCore( idBuffer, ascDiagDistance, x, y,
                image_.width, distance, newID, oldID, lastID,
                lastDist, offset, step);
        }
    }

    // fourth pass (start from bottom)
    oldID = 0;
    newID = 0;
    for (int x_t = bl.x + step; x_t < tr.x; x_t += step) {
        int distance = 0;
        int x = x_t;
        for (int y = bl.y; y < tr.y && x < tr.x; ) {
            x += step;
            y += step;
            anchorAlgorithmCore( idBuffer, ascDiagDistance, x, y,
                image_.width, distance, newID, oldID, lastID,
                lastDist, offset, step);
        }
    }


    //-----------------------------

    // find minimal distances by adding all 4 single-direction distance maps.
    // blur the resulting distance map.
    for (int y = bl.y; y < tr.y; y += step) {
        for (int x = bl.x; x < tr.x; x += step) {

            int id = idBuffer[y*image_.width + x];

            if (id > 0) {
                // blur the distancemap pixel
                int gap = (filterDelta_.get()/distanceMapStep_.get())*distanceMapStep_.get();
                int distance = 0;
                int* kernel = kernels_[filterKernel_.getValue()]->kernel;
                int kernelSize = kernels_[filterKernel_.getValue()]->size;
                for (int kX = 0; kX < kernelSize; kX++) {
                    for (int kY = 0; kY < kernelSize; kY++) {
                        int x_t = x + (kX - (kernelSize/2))*gap;
                        int y_t = y + (kY - (kernelSize/2))*gap;
                        if ( x_t >= bl.x && x_t <= tr.x && y_t >= bl.y && y_t <= tr.y &&
                             idBuffer[y_t*image_.width + x_t] == id) {
                            int coeff = kernel[kX*kernelSize + kY];
                            distance += horzDistance[y_t*image_.width + x_t] * coeff;
                            distance += vertDistance[y_t*image_.width + x_t] * coeff;
                            distance += ascDiagDistance[y_t*image_.width + x_t] * coeff;
                            distance += descDiagDistance[y_t*image_.width + x_t] * coeff;
                        }
                    }
                }
                distanceField[y*image_.width + x] = distance / kernels_[filterKernel_.getValue()]->coefficientsSum;

                // Set the anchor point.
                // if not set anchor or if set and the new distanceField point is larger than the old
                // it will set the new distanceField point. (anchors[id].x == -1) <=> (anchor not set)
                if ( ( (anchors[id].x > -1) &&
                    ( distanceField[anchors[id].y*image_.width + anchors[id].x] <
                    distanceField[y*image_.width + x]) )
                    || (anchors[id].x == -1) )
                        anchors[id] = ivec2(x,y);

            // set all distancemap pixel not belonging to any segment (id==0) to 0
            } else {
                distanceField[y*image_.width + x] = 0;
                horzDistance[y*image_.width + x] = 0;
                vertDistance[y*image_.width + x] = 0;
                ascDiagDistance[y*image_.width + x] = 0;
                descDiagDistance[y*image_.width + x] = 0;
            }
        }
    }


    // create label for each detected segment/anchor point
    for (int i = 0; i<maxLabelID; i++) {
        if ((anchors[i].x > -1) && (distanceField[anchors[i].y*image_.width + anchors[i].x] >
            (minSegmentSize_.get() / step) )) {
            LabelData* labelData = NULL;
            for (size_t j=0; (j<labelPersistentData_.size()) && !labelData; j++) {
                if (i == labelPersistentData_[j]->id)
                    labelData = labelPersistentData_[j];
            }
            if (labelData && labelData->text.text.length() > 0) {
                Label newLabel;
                // TODO: check if it's still the case (ab)
                // BezierPatch does not work with display lists on ATI :(
                if (GpuCaps.getVendor() == GpuCapabilities::GPU_VENDOR_ATI)
                    newLabel.bezierPatch = BezierPatch(false);
                newLabel.anchorPoint = anchors[i];
                newLabel.segmentSize = segmentSize[i];
                newLabel.labelData = labelData;
                newLabel.rotAngle = 0.f;
                newLabel.curve2D = NULL;
                newLabel.offLabel = false;
                labels_.push_back(newLabel);
            }
        }
    }

    // save bounding box to image_
    image_.bl = bl;
    image_.tr = tr;

    // save calculated distance fields into image_
    image_.distanceField.setData(distanceField, image_.width, image_.height, 1);
    image_.horzDistance.setData(horzDistance, image_.width, image_.height, 1);
    image_.vertDistance.setData(vertDistance, image_.width, image_.height, 1);
    image_.ascDiagDistance.setData(ascDiagDistance, image_.width, image_.height, 1);
    image_.descDiagDistance.setData(descDiagDistance, image_.width, image_.height, 1);
}

void VolumeLabeling::appendUnsegmentedLabels() {
    // create unsegmented labels
    LabelVec offLabels;
    for (size_t i=0; i<labelPersistentData_.size(); i++) {
        if (labelPersistentData_[i]->belongsToSegment)
            continue;
        Label newLabel;
        newLabel.labelData = labelPersistentData_[i];
        newLabel.intern = false;
        newLabel.offLabel = true;
        newLabel.rotAngle = 0.f;
        offLabels.push_back(newLabel);
    }

    // stack off-labels in the lower-left corner of the window
    //if (offLabels.size() > 1){
        placeOffLabels(offLabels);
        for (size_t i=0; i<offLabels.size(); i++) {
            Curve2DPolynomial* curve2D = new Curve2DPolynomial(1);
            vector<vec2> points;
            vec2 offsetVector = vec2(offLabels[i].labelData->text.width / 2.f, 0.f);
            points.push_back((vec2)offLabels[i].labelPos - offsetVector);
            points.push_back((vec2)offLabels[i].labelPos + offsetVector);
            curve2D->setCtrlPoints(points, static_cast<float>(offLabels[i].labelData->text.width));
            offLabels[i].curve2D = curve2D;
        }
    //}

    // append off-labels to global label vector
    for (size_t i = 0; i < offLabels.size(); ++i)
        labels_.push_back(offLabels[i]);
}

// calculate label positions and connection points for external labels
// depending on current layout-type
// and bezier patches for internal ones
void VolumeLabeling::calcLabelPositions() {

#ifdef labelDEBUG
    renderPoints_.clear();
#endif

    if (labels_.size() == 0)
        return;

    // separate intern from extern labels
    LabelVec internLabels, externLabels;
    for (size_t i = 0; i < labels_.size(); ++i) {
            Label label = labels_[i];
            if (label.labelData->internPreferred)
                internLabels.push_back(label);
            else {
                label.intern = false;
                externLabels.push_back(label);
            }
    }

    //
    // calculate control points for bezier patches of intern labels
    //
    for (size_t i = 0; i < internLabels.size(); i++) {

        // find 2D label path and and 3D bezier patch control points from it
        // order of the conditionals is important (shortcut-evaluation!)
        if ( findLabelPathBest(internLabels[i]) && findBezierPoints(internLabels[i]) )
            internLabels[i].intern = true;
        else {
            internLabels[i].intern = false;
            externLabels.push_back(internLabels[i]);
        }

    } // calculations of internal labels


    // determine layout of external labels
    if (layout_.get() == "left-right") {

    //---- left-right layout ---- //
        LabelVec leftLabels;
        LabelVec rightLabels;

        // sort labels by their horizontal position (ascending)
        for (size_t i = 0; i < externLabels.size(); ++i) {
            for (size_t j = i+1; j < externLabels.size(); ++j) {
                if (externLabels[i].anchorPoint.x > externLabels[j].anchorPoint.x ) {
                    swap(externLabels[i], externLabels[j]);
                }
            }
        }

        // left half of labels are put on the left side
        // right half on the right side
        // (but dont put a label to far away from its anchor)
        size_t mid = tgt::iround(externLabels.size() / 2.f);
        // calculate initial label layout
        int objectWidth = image_.tr.x-image_.bl.x;
        for (size_t i = 0; i < externLabels.size(); ++i) {
            int labelRelPos = externLabels[i].anchorPoint.x-image_.bl.x;
            if ( ( (i >= mid && labelRelPos > static_cast<int>(objectWidth*(1/3.f)) )  ||
                ( (labelRelPos > static_cast<int>(objectWidth*(2/3.f)) ))  ) ) {
                    externLabels[i].connectionPoint.x = image_.tr.x + gaps_.LR_BboxObject;
                    externLabels[i].connectionPoint.y = externLabels[i].anchorPoint.y;
                    rightLabels.push_back(externLabels[i]);
            }
            else {
                externLabels[i].connectionPoint.x = image_.bl.x - gaps_.LR_BboxObject;
                externLabels[i].connectionPoint.y = externLabels[i].anchorPoint.y;
                leftLabels.push_back(externLabels[i]);
            }
        }

        externLabels.clear();

        // sort left labels by their vertical label position (ascending)
        for (size_t i = 0; i < leftLabels.size(); ++i) {
            for (size_t j = i+1; j < leftLabels.size(); ++j) {
                if (leftLabels[i].connectionPoint.y > leftLabels[j].connectionPoint.y ) {
                    swap(leftLabels[i], leftLabels[j]);
                }
            }
        }

        // sort right labels by their vertical label position (ascending)
        for (size_t i = 0; i < rightLabels.size(); ++i) {
            for (size_t j = i+1; j < rightLabels.size(); ++j) {
                if (rightLabels[i].connectionPoint.y > rightLabels[j].connectionPoint.y ) {
                    swap(rightLabels[i], rightLabels[j]);
                }
            }
        }

        // resolve vertical label overlaps
        int bottom, top;
        if (leftLabels.size() > 1) {
            stackLabels(leftLabels, 0, static_cast<int>(leftLabels.size())-1, bottom, top);
        }
        if (rightLabels.size() > 1) {
            stackLabels(rightLabels, 0, static_cast<int>(rightLabels.size())-1, bottom, top);
        }

        // resolve line intersections
        resolveLineIntersections(leftLabels);
        resolveLineIntersections(rightLabels);

        // recombine left and right labels to one label vector
        for (size_t i = 0; i < leftLabels.size(); ++i) {
            externLabels.push_back(leftLabels[i]);
        }
        for (size_t i = 0; i < rightLabels.size(); ++i) {
            externLabels.push_back(rightLabels[i]);
        }

        // finally calculate label positions from connpoints
        for (size_t i = 0; i < externLabels.size(); ++i) {
            labelPosFromConnPoint(externLabels[i], externLabels[i].labelPos);

            // Fixme: labelPos is now the center of the label (before bottom-left corner)
            externLabels[i].labelPos += ivec2(externLabels[i].labelData->text.width/2,
                externLabels[i].labelData->text.height/2);
        }
    }
    else if (layout_.get() == "silhouette") {

    //---- silhouette layout ---- //

        // first find connection point closest to its anchor

        // number of projections per anchor point
        int stepCount = 60;

        float stepAngle = deg2rad(360.f / stepCount);

        // for each anchor point the distance to its nearest connpoint (so far) is saved
        int* nearest = new int[externLabels.size()];
        for (size_t i=0; i<externLabels.size(); ++i)
            nearest[i] = image_.width*image_.height;

        // dir is the direction of projection (DOP); it is rotated after each iteration
        // start with a normalized direction vector
        vec2 dir = vec2(0,1);
        mat2 rotationMatrix = mat2(cosf(stepAngle), sinf(stepAngle),
                                   -sinf(stepAngle), cosf(stepAngle));

        // for each anchor:
        // find projection of anchor point onto hull
        // regarding current DOP (vector dir).
        // if new connpoint is the closest for its anchor, it is saved
        for (int step=0; step < stepCount - 1; ++step) {
            for (size_t i=0; i < externLabels.size(); ++i) {

                ivec2 anchor = externLabels[i].anchorPoint;

                // project anchor point onto image boundary
                // using dir as DOP.
                // projected point is 'end'
                ivec2 end;
                float scale_x = (image_.width - 1 - anchor.x) / dir.x;
                if (scale_x < 0)
                    scale_x = (- anchor.x) / dir.x;

                float scale_y = (image_.height - 1 - anchor.y) / dir.y;
                if (scale_y < 0)
                    scale_y = (- anchor.y) / dir.y;

                if (scale_x < scale_y) {
                    end.x = anchor.x + tgt::iround(dir.x*scale_x);
                    end.y = anchor.y + tgt::iround(dir.y*scale_x);
                }
                else {
                    end.x = anchor.x + tgt::iround(dir.x*scale_y);
                    end.y = anchor.y + tgt::iround(dir.y*scale_y);
                }

                // find hull intersection of line anchor->end
                ivec2 intersect, normal;
                midPointLine(anchor, end, HullIntersect, intersect, normal);
                intersect.x += tgt::iround(dir.x*gaps_.SIL_LabelHull);
                intersect.y += tgt::iround(dir.y*gaps_.SIL_LabelHull);
                int len = dot(intersect-anchor, intersect-anchor);
                // if new connpoint is the closest, save it
                if ( len < nearest[i] ) {
                    nearest[i] = len;
                    externLabels[i].connectionPoint = intersect;
                    externLabels[i].normal = normal;
                }

            }
            // rotate dir-vector
            dir = rotationMatrix * dir;
        }

     // calculate initial label positions from connpoints
     for (size_t i = 0; i < externLabels.size(); ++i) {
         labelPosFromConnPoint(externLabels[i], externLabels[i].labelPos);
         correctHullGap(externLabels[i], externLabels[i].labelPos);
     }

     // resolve label / line overlaps
     bool finished = false;
     for (int i=0; (i < maxIterations_.get()) && !finished ; ++i) {
         finished = true;
         finished &= resolveLineIntersections(externLabels);
         finished &= resolveLabelOverlaps(externLabels);
     }

     // finally calculate connpoints from corrected labelpositions
     for (size_t i=0; i < externLabels.size(); ++i)
        connPointFromLabelPos(externLabels[i], externLabels[i].connectionPoint, externLabels[i].normal);

#ifdef labelDEBUG
     if (!finished)
         LINFO("Layout not finished! (iteration limit reached: " <<
                      maxIterations_.get() << ")" << std::endl );
#endif

     delete[] nearest;
    } // switch

    // recombine intern and extern labels
    labels_.clear();
    for (size_t i=0; i < externLabels.size(); i++)
        labels_.push_back(externLabels[i]);

    for (size_t i=0; i < internLabels.size(); i++)
        if (internLabels[i].intern)
            labels_.push_back(internLabels[i]);
}


// maps from texture (pixel) coords to world coords (NDC).
// if a label is supplied, only this label's coords are transformed,
// if not, all labels are transformed
void VolumeLabeling::toWorld(Label* pLabel) {
    const vec2 scale = vec2(2.f/(image_.width-0.5f), 2.f/(image_.height-0.5f));
    const vec2 shift = vec2(-1.f, -1.f);

    if (pLabel == NULL) {
        for (size_t i = 0; i < labels_.size(); i++)
            toWorld( &(labels_[i]) );
    }
    else {
        pLabel->anchorPointWorld = (vec2)pLabel->anchorPoint*scale + shift;
        pLabel->connectionPointWorld = (vec2)pLabel->connectionPoint*scale + shift;
        pLabel->labelPosWorld = (vec2)(pLabel->labelPos)*scale + shift;
        pLabel->labelData->text.widthWorld = static_cast<float>(pLabel->labelData->text.width) * scale.x;
        pLabel->labelData->text.heightWorld = static_cast<float>(pLabel->labelData->text.height) * scale.y;
        pLabel->labelData->text.widthInternWorld = static_cast<float>(pLabel->labelData->text.widthIntern) * scale.x;
        pLabel->labelData->text.heightInternWorld = static_cast<float>(pLabel->labelData->text.heightIntern) * scale.y;

        if ( !shape3D_.get() ) {
            BezierPatch patch = pLabel->bezierPatch;
            int dimS, dimT;
            vec3* ctrlPoints = patch.getCtrlPoints(dimS, dimT);
            for (int i=0; i<(dimS+1)*(dimT+1); ++i) {
                ctrlPoints[i].x = ctrlPoints[i].x*scale.x + shift.x;
                ctrlPoints[i].y = ctrlPoints[i].y*scale.y + shift.y;
            }
        }
    }
}


// renders labels, anchor points and connection lines
// AFTER calculation of anchorpoints, label positions, conn-points, bezier patches
//
// 2 rendering passes (the first only for halo):
//   1. pass: render connection lines and anchors
//(halo only) using halo color and thicker lines / points
//   2. pass: render connection lines and anchors with normal thickness
//            + render quads at label positions mapping a font texture onto them (external),
//              render bezier patches mapping a font texture onto them (internal)
void VolumeLabeling::renderLabels() {
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.pushMatrix();
    MatStack.loadIdentity();

    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();
    MatStack.loadIdentity();

    glDisable(GL_TEXTURE_2D);

    labelingPort_.activateTarget();

    glEnable(GL_ALPHA_TEST);
    glAlphaFunc(GL_GREATER, 0.f);
    //glDisable(GL_DEPTH_TEST);
    glDepthFunc(GL_ALWAYS);
    glDepthMask(false);

    for (int pass = 1; pass <= 2; ++pass) {
        // settings for different render passes
        if (pass == 1) {
            if ( drawHalo_.get() ) {
                glPointSize(4.0f);
                glLineWidth(3.0f);
                glColor4fv(haloColorExtern_.get().elem);
            }
            else
                continue;
        }
        else if (pass == 2) {
            glPointSize(2.5f);
            glLineWidth(1.0f);

            glColor4fv(fontExtern_.get()->getFontColor().elem);
        }

        LGL_ERROR;

        // anchor points and connection lines for external labels
        if ( showLabels_.get() == "all" || showLabels_.get() == "external-only" ) {
            // anchor points (only for extern labels)
            glBegin(GL_POINTS);
            for (size_t i = 0; i < labels_.size(); ++i) {
                if (    !labels_[i].intern && !(labels_[i].offLabel) &&
                    !(drag_ && pickedLabel_ && pickedLabel_->labelData->id == labels_[i].labelData->id)) {
                        glVertex3f(labels_[i].anchorPointWorld.x, labels_[i].anchorPointWorld.y, -1.0 );
                }
            }
            glEnd();

            // lines connecting anchors and connection points (only for extern labels)
            glBegin(GL_LINES);
            for (size_t i = 0; i < labels_.size(); ++i) {
                if (    !labels_[i].intern && !(labels_[i].offLabel) &&
                    !(drag_ && pickedLabel_ && pickedLabel_->labelData->id == labels_[i].labelData->id)) {
                        glVertex3f(labels_[i].anchorPointWorld.x, labels_[i].anchorPointWorld.y, -1 );
                        glVertex3f(labels_[i].connectionPointWorld.x, labels_[i].connectionPointWorld.y, -1);
                }
            }
            glEnd();
        }

        LGL_ERROR;

        // render anchor and connection-line of dragged label
        if (drag_ && pickedLabel_ && !pickedLabel_->intern && !pickedLabel_->offLabel) {

            glBegin(GL_POINTS);
            glVertex3f(pickedLabel_->anchorPointWorld.x, pickedLabel_->anchorPointWorld.y, -1.0 );
            glEnd();
            ivec2 bl = pickedLabel_->labelPos
                - ivec2(tgt::iround(pickedLabel_->labelData->text.width/2.f), tgt::iround(pickedLabel_->labelData->text.height/2.f) )
                - ivec2(gaps_.SIL_LabelLabel);
            ivec2 tr = pickedLabel_->labelPos
                + ivec2(tgt::iround(pickedLabel_->labelData->text.width/2.f), tgt::iround(pickedLabel_->labelData->text.height/2.f))
                + ivec2(gaps_.SIL_LabelLabel);
            if (!pointInBox(bl, tr, pickedLabel_->anchorPoint)) {
                glBegin(GL_LINES);
                glVertex3f(pickedLabel_->anchorPointWorld.x, pickedLabel_->anchorPointWorld.y, -1 );
                glVertex3f(pickedLabel_->connectionPointWorld.x, pickedLabel_->connectionPointWorld.y, -1);
                glEnd();
            }

        }

        LGL_ERROR;

    }

    for(size_t i = 0; i != labels_.size(); i++){
        tgt::Font *font = 0;

        if ((!labels_[i].intern ) &&
                    ( showLabels_.get() == "all" || showLabels_.get() == "external-only" ||
                      (drag_ && pickedLabel_->labelData->id == labels_[i].labelData->id )) ) {
            font = fontExtern_.get();
        }else if ( labels_[i].intern && (showLabels_.get() == "all" || showLabels_.get() == "internal-only")) {
            font = fontIntern_.get();
        }

        if (font){
            std::string lbl = labels_[i].labelData->text.text;
            tgt::vec3 pos = tgt::vec3(tgt::vec2(labels_[i].labelPos), 0.0f);
            float len = font->getSize(pos, lbl, idMapPort_.getSize()).x;
            font->render(pos-tgt::vec3(len*0.5f, 0, 0), lbl, idMapPort_.getSize());
        }
    }


    // reset settings
    TextureUnit::setZeroUnit();
    LGL_ERROR;
    glDepthFunc(GL_LESS);
    glDisable(GL_BLEND);
    glDisable(GL_ALPHA_TEST);
    glDepthMask(true);

#ifdef labelDEBUG
    // convex hull for debugging
    if ( (drawConvexHull_.get()) && (layout_->get() == SILHOUETTE) && (image_.hullPoints.size() > 0)) {
        glColor4fv(vec4::red.elem);
        glLineWidth(2.0);
        glBegin(GL_LINE_LOOP);
        for (size_t i = 0; i < image_.hullPoints.size()-1; ++i) {
            glVertex3f(image_.hullPoints[i].x*2.f/image_.width-1, image_.hullPoints[i].y*2.f/image_.height-1, -1 );
        }
        glEnd(); };

    glColor4f(1.f, 0.f, 0.f, 1.f);
    glPointSize(1.f);
    glBegin(GL_POINTS);
    for (size_t i=0; i<renderPoints_.size(); ++i) {
        glVertex3f( renderPoints_[i].x, renderPoints_[i].y, -1.f );
    }
    glEnd();
#endif

    glLineWidth(1.f);
    glPointSize(1.f);

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.popMatrix();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.popMatrix();

    labelingPort_.deactivateTarget();

    LGL_ERROR;

}

// clear current label vector
void VolumeLabeling::deleteLabels() {
    for (size_t i=0; i<labels_.size(); ++i) {
        if (labels_[i].curve2D != NULL) {
            delete labels_[i].curve2D;
        }
    }
    labels_.clear();
}

// calculates connection point from label's position (external labels only)
void VolumeLabeling::connPointFromLabelPos(const Label &pLabel, ivec2 &connPoint, ivec2 &normal, bool calcNormal) {
    Adjustments h_adjust;    // connPoint is placed left|right, up|down from label

    if (layout_.get() == "silhouette" || layout_.get() == "left-right"){

        // gradient thresholds for chosing one of the potential centers for the line connpoint.
        // after center is found, the connpoint is found by projecting the center onto the label's
        // textbox-hull in the direction of the anchor
        float l_thresh = 3.0;
        float u_thresh = 8.0;

        // label's adjustment regarding its connpoint,
        // depending on its relative position to anchor
        h_adjust = (pLabel.labelPos.x >= pLabel.anchorPoint.x) ? RIGHT : LEFT;

        // these are the potential centers
        ivec2 center_middle, center_out, center_semimiddle;

        // line would be connected to center of label
        center_middle.y = pLabel.labelPos.y;
        // line would be connected to the label's outside (left or right)
        center_out.y = pLabel.labelPos.y;
        // line would be connected to the center point between center_middle and center_out
        center_semimiddle.y = pLabel.labelPos.y;

        center_middle.x = pLabel.labelPos.x;
        if (h_adjust == LEFT) {
            center_out.x = pLabel.labelPos.x + tgt::iround(pLabel.labelData->text.width/2.f) -
                tgt::iround(pLabel.labelData->text.height/2.f);
            center_semimiddle.x = pLabel.labelPos.x + tgt::iround( (1/4.f)*pLabel.labelData->text.width );
        }
        else {
            center_out.x = pLabel.labelPos.x - tgt::iround(pLabel.labelData->text.width/2.f) +
                tgt::iround(pLabel.labelData->text.height/2.f);
            center_semimiddle.x = pLabel.labelPos.x - tgt::iround( (1/4.f)*pLabel.labelData->text.width );
        }

        // connection lines from anchor to potential connection points
        ivec2 t_middle = center_middle - pLabel.anchorPoint;
        ivec2 t_out = center_out - pLabel.anchorPoint;
        ivec2 t_semimiddle = center_semimiddle - pLabel.anchorPoint;

        // now choose one center by comparing the connection lines' gradients
        // to the thresholds
        ivec2 center;
        if ( (fabs(static_cast<float>(t_out.y) / t_out.x) < l_thresh) &&
                tgt::sign( (t_middle.x-pLabel.labelData->text.width/2)*(t_middle.x+pLabel.labelData->text.width/2)) == 1 )
            center = center_out;
        else if ( (fabs(static_cast<float>(t_middle.y) / t_middle.x) > u_thresh ) ||
            (tgt::sign(t_middle.x*t_semimiddle.x) == -1) )
            center = center_middle;
        else
            center = center_semimiddle;


        // finally we have to project the center onto the label's textbox-hull
        // (in the direction of the anchor point)
        ivec2 r = pLabel.anchorPoint - center;
        float s_y, s_x;
        if (r.y != 0)
            s_y = (pLabel.labelData->text.height/2.f + gaps_.LabelConnpoint) / abs(r.y);
        else
            s_y = static_cast<float>(image_.width*image_.height);

        if (r.x != 0) {
            if (h_adjust == LEFT)
                s_x = ( pLabel.labelData->text.width/2.f -
                    (center.x - pLabel.labelPos.x) + gaps_.LabelConnpoint) / static_cast<float>(r.x);
            else
                s_x = ( - pLabel.labelData->text.width/2.f +
                    ( pLabel.labelPos.x - center.x) - gaps_.LabelConnpoint) / static_cast<float>(r.x);
        }
        else
            s_x = static_cast<float>(image_.width*image_.height);

        // set connpoint
        float s = (s_x < s_y ? s_x : s_y);
        connPoint = center + (ivec2)(s*(vec2)r);

        // update hullintersection
        if (calcNormal) {
            ivec2 inter;
            midPointLine(connPoint, pLabel.anchorPoint, HullIntersect, inter, normal);
        }

    }
}


// calc position of one label using its connpoint (external labels only)
void VolumeLabeling::labelPosFromConnPoint(const Label &pLabel, ivec2 &labelPos) {
    int x_gap = gaps_.LabelConnpoint;

    Adjustments h_adjust, v_adjust;    // label is placed left|right, up|down from conn-point

    if (layout_.get() == "left-right") {
        h_adjust = (pLabel.anchorPoint.x < pLabel.connectionPoint.x) ? RIGHT : LEFT;

        if (h_adjust == RIGHT) {
            labelPos.x = pLabel.connectionPoint.x + x_gap;
        }
        else {
            labelPos.x = static_cast<int>(pLabel.connectionPoint.x - pLabel.labelData->text.width
                                          - x_gap);
        }
        labelPos.y = pLabel.connectionPoint.y - (pLabel.labelData->text.height-1) / 2;
    }
    else if (layout_.get() == "silhouette") {
        h_adjust = (pLabel.normal.x >= 0) ? RIGHT : LEFT;
        v_adjust = (pLabel.normal.y >= 0) ? UP : DOWN;

        float l_thresh = 3.0;
        float u_thresh = 8.0;

        ivec2 r = pLabel.connectionPoint - pLabel.anchorPoint;

        ivec2 shift;
        // connection-line has a small gradient
        if ((r.x != 0) && ((static_cast<float>(abs(r.y / r.x))) < l_thresh) ) {
            if (h_adjust == LEFT)
                shift.x = - tgt::iround(pLabel.labelData->text.width/2.f) - gaps_.LabelConnpoint;
            else
                shift.x = tgt::iround(pLabel.labelData->text.width/2.f) + gaps_.LabelConnpoint;
            shift.y = 0;
        // connection-line has a huge gradient
        }
        else if ((r.x == 0) || ((static_cast<float>(abs(r.y / r.x))) > u_thresh) ) {
            shift.x = 0;
            if (v_adjust == DOWN)
                shift.y = - (tgt::iround(pLabel.labelData->text.height/2.f) + gaps_.LabelConnpoint);
            else
                shift.y = (tgt::iround(pLabel.labelData->text.height/2.f) + gaps_.LabelConnpoint);
        // connection-line gradient between lower and upper threshold
        }
        else {
            shift.x = tgt::iround(pLabel.labelData->text.width/4.f * sign(r.x));
            if (v_adjust == DOWN)
                shift.y = - ( tgt::iround(pLabel.labelData->text.height/2.f) + gaps_.LabelConnpoint );
            else
                shift.y = ( tgt::iround(pLabel.labelData->text.height/2.f) + gaps_.LabelConnpoint );
        }

        labelPos = pLabel.connectionPoint + shift;
    }
}

// resolves intersections between two label-connecting lines depending on the current layout.
// both label positions and connpoints are corrected!
// @return false, if intersections have been found
bool VolumeLabeling::resolveLineIntersections(LabelVec &pLabels) {
    bool intersectionFound = false;
    bool finished = false;

    ivec2 intersection;
    bool onSegment;
    vec2 gapi, gapj;

    // while intersections are found, swap each two labels whose lines intersect
    for (int iter=0; (iter<maxIterations_.get()) && !finished; ++iter) {
        finished = true;

        // test each label against each other
        for (size_t i = 0; i < pLabels.size(); ++i) {
            for (size_t j = i+1; j < pLabels.size(); ++j) {

                // enlarge lines by line-gap
                gapi = normalize( static_cast<vec2>(pLabels[i].connectionPoint - pLabels[i].anchorPoint) );
                gapj = normalize( static_cast<vec2>(pLabels[j].connectionPoint - pLabels[j].anchorPoint) );
                gapi = gapi*static_cast<float>(gaps_.LineLine);
                gapj = gapj*static_cast<float>(gaps_.LineLine);

                // test, if lines intersect
                lineIntersection(
                    pLabels[i].connectionPoint+(ivec2)gapi,
                    pLabels[i].anchorPoint-(ivec2)gapi,
                    pLabels[j].connectionPoint+(ivec2)gapj,
                    pLabels[j].anchorPoint-(ivec2)gapj,
                    intersection, onSegment);
                if (onSegment) {
                    if (layout_.get() == "left-right") {
                        // swap connpoints and recalculate label positions afterwards
                        swap(pLabels[i].connectionPoint, pLabels[j].connectionPoint);
                        labelPosFromConnPoint(pLabels[i], pLabels[i].labelPos);
                        labelPosFromConnPoint(pLabels[j], pLabels[j].labelPos);

                    }
                    else if (layout_.get() == "silhouette") {
                        // swap connpoints
                        swap(pLabels[i].connectionPoint, pLabels[j].connectionPoint);

                        // recalculate label positions from new connpoints
                        labelPosFromConnPoint(pLabels[i], pLabels[i].labelPos);
                        labelPosFromConnPoint(pLabels[j], pLabels[j].labelPos);
                        correctHullGap(pLabels[i], pLabels[i].labelPos);
                        correctHullGap(pLabels[j], pLabels[j].labelPos);

                        // calculate final connpoints from label positions
                        connPointFromLabelPos(pLabels[i], pLabels[i].connectionPoint, pLabels[i].normal);
                        connPointFromLabelPos(pLabels[j], pLabels[j].connectionPoint, pLabels[i].normal);

                    }
                    finished = false;
                    intersectionFound = true;
                }
            }
        }
    }
    return !intersectionFound;
}


//  returns the intersection point of line through P1,Q1 and P2,Q2.
//  @intersection the calculated intersection point
//  @onSegments   true, if the intersection point lies on line segment between P1,Q1 and P2,Q2
void VolumeLabeling::lineIntersection(ivec2 P1, ivec2 Q1, ivec2 P2, ivec2 Q2,
                              ivec2 &intersection, bool &onSegments) {
    float s;
    float t;
    if (((Q2.y - P2.y) * (Q1.x - P1.x) == 0) && ((Q2.x - P2.x) * (Q1.y - P1.y) == 0)) {
        s = -1;
        t = -1;
    }
    else {
        s = static_cast<float>((Q2.x - P2.x) * (P1.y - P2.y) - (Q2.y - P2.y) * (P1.x - P2.x)) /
            ((Q2.y - P2.y) * (Q1.x - P1.x) - (Q2.x - P2.x) * (Q1.y - P1.y));
        t = static_cast<float>((Q1.x - P1.x) * (P1.y - P2.y) - (Q1.y - P1.y) * (P1.x - P2.x)) /
            ((Q2.y - P2.y) * (Q1.x - P1.x) - (Q2.x - P2.x) * (Q1.y - P1.y));
    }

    intersection = P1 + (ivec2)(s*(vec2)(Q1 - P1));
    onSegments = ((s >= 0) && (s <= 1) && (t >= 0) && (t <= 1));
}

// draws a one pixel thick halo with haloColor
// around all objects in bitmap with alpha-value greater than alphaTreshold.
// Egde halos get an alpha-value of 1.0, corner-halos 0.5 (=> rounded corners)
// @bitmap a four-channel RGBA bitmap
void VolumeLabeling::drawHalo(labeling::Bitmap<GLfloat> &bitmapParam, int width, int height,
                          tgt::Color haloColor, GLfloat alphaThreshold)
{
    GLfloat* result = new GLfloat[width*height*4];
    GLfloat* bitmap = bitmapParam.getData();

    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            GLfloat center = bitmap[(y*width+x)*4 + 3];

            GLfloat alpha=0.f;
            if (center <= alphaThreshold) {
                GLfloat east = bitmap[(y*width+min(x+1,width-1))*4 + 3];
                GLfloat west = bitmap[(y*width+max(x-1,0))*4 + 3];
                GLfloat north = bitmap[(min(y+1,height-1)*width+x)*4 + 3];
                GLfloat south = bitmap[(max(y-1,0)*width+x)*4 + 3];
                if (east > alphaThreshold || south > alphaThreshold ||
                    west > alphaThreshold || north > alphaThreshold) {
                        alpha = 1.f;
                }
                else {
                    GLfloat southeast = bitmap[(max(y-1,0)*width+min(x+1,width-1))*4 + 3];
                    GLfloat southwest = bitmap[(max(y-1,0)*width+max(x-1,0))*4 + 3];
                    GLfloat northwest = bitmap[(min(y+1,height-1)*width+max(x-1,0))*4 + 3];
                    GLfloat northeast = bitmap[(min(y+1,height-1)*width+min(x+1,width-1))*4 + 3];
                    if (southeast > alphaThreshold || southwest > alphaThreshold ||
                        northwest > alphaThreshold || northeast > alphaThreshold) {
                        alpha = 0.5f;
                    }
                }
            }

            if (alpha > 0.f) {
                result[(y*width+x)*4] = haloColor.r;
                result[(y*width+x)*4 + 1] = haloColor.g;
                result[(y*width+x)*4 + 2] = haloColor.b;
                result[(y*width+x)*4 + 3] = alpha;
            }
            else {
                result[(y*width+x)*4] = bitmap[(y*width+x)*4];
                result[(y*width+x)*4 + 1] = bitmap[(y*width+x)*4 + 1];
                result[(y*width+x)*4 + 2] = bitmap[(y*width+x)*4 + 2];
                result[(y*width+x)*4 + 3] = center;
            }
        }
    }

    bitmapParam.setData(result, width, height, 4);
}

// divide-and-conquer algorithm for resolving vertical label overlaps in left-right layout:
// 1. resolve overlaps in bottom half of labels
// 2. resolve overlaps in top half of labels
// 3. combine results: if bottom and top labeling regions overlap,
//                     shift labels downwards or upwards respectively
void VolumeLabeling::stackLabels(LabelVec& pLabels, int min, int max, int& bottom, int& top) {
    int y_gap = gaps_.LR_LabelLabel;

    if (min == max) {
        bottom = pLabels[min].connectionPoint.y - tgt::iround(pLabels[min].labelData->text.height/2.f);
        top = pLabels[min].connectionPoint.y + tgt::iround(pLabels[min].labelData->text.height/2.f);
        return;
    }
    else {
        // divide
        int mid = min + (max - min) / 2;
        int bottom1, top1, bottom2, top2;
        stackLabels(pLabels, min, mid, bottom1, top1);
        stackLabels(pLabels, mid+1, max, bottom2, top2);

        // combine results
        if (top1 > bottom2 - y_gap) {

            int overlap = top1 - bottom2 + y_gap;
            pLabels[mid].connectionPoint.y -= overlap / 2;
            pLabels[mid+1].connectionPoint.y += overlap / 2;

            for (int i = mid-1; i >= min; --i) {
                overlap = pLabels[i+1].connectionPoint.y - pLabels[i].connectionPoint.y -
                    tgt::iround(pLabels[i+1].labelData->text.height/2.f) - tgt::iround(pLabels[i].labelData->text.height/2.f) - y_gap;
                if (overlap < 0)
                    pLabels[i].connectionPoint.y += overlap;
            }

            for (int i = mid+2; i <= max; ++i) {
                overlap = pLabels[i].connectionPoint.y - pLabels[i-1].connectionPoint.y  -
                    tgt::iround(pLabels[i].labelData->text.height/2.f) - tgt::iround(pLabels[i-1].labelData->text.height/2.f) - y_gap;
                if (overlap < 0)
                    pLabels[i].connectionPoint.y += -overlap;
            }

            bottom = pLabels[min].connectionPoint.y - tgt::iround(pLabels[min].labelData->text.height/2.f);
            top = pLabels[max].connectionPoint.y + tgt::iround(pLabels[max].labelData->text.height/2.f);

        }
        else {
            bottom = bottom1;
            top = top2;
        }
    }
}


// "draws" a line from @start to @end using MidPoint-Line
// @action specifies the action that should be performed at a raster position:
//      * hullIntersect: find the intersection of line start->end with conex hull
//                       returns @intersection point and @normal at this point
//      * drawHull: draws a 3-pixel-thick line into image_.convexHull and calcs and sets normal
// @return true, if action==hullIntersect and intersection found
//         false, if action==hullIntersect and no intersection found
//         undefined otherwise
bool VolumeLabeling::midPointLine(ivec2 const& start, const ivec2& end, MidPointLineActions action,
                              ivec2& intersection, ivec2& normal)
{
     // clamp coords to image boundary
    int x0 = min(max(start.x,0),image_.width-1);
    int x1 = min(max(end.x,0),image_.width-1);
    int y0 = min(max(start.y,0),image_.height-1);
    int y1 = min(max(end.y,0),image_.height-1);

    int dist = -1;
    int grad_x = end.y - start.y;;
    int grad_y = -(end.x - start.x);;

    bool mirroringWH = false;
    bool mirroringXA = false;

    if (x1 < x0) {
        swap(x0,x1);
        swap(y0,y1);
    }

    if (y1 < y0) {
        y1 = y0 + (y0 - y1);
        mirroringXA = true;
    }

    int dx = x1 - x0;
    int dy = y1 - y0;
    if (dy > dx) {
        swap(x0,y0);
        swap(x1,y1);
        mirroringWH = true;
    }

    dx = x1 - x0;
    dy = y1 - y0;
    int d = 2*dy - dx;
    int deltaE = 2*dy;
    int deltaNE = 2*(dy - dx);
    int x = x0;
    int y = y0;

    while (x <= x1) {
        if (mirroringWH) {
            swap(x,y);
            swap(x0,y0);
        }

        //writePixel(x, y);
        int x_t;
        int y_t;
        if (mirroringXA) {
            x_t = x;
            y_t = 2*y0 - y;
        }
        else {
            x_t = x;
            y_t = y;
        }

        switch (action) {
            case HullIntersect:
                {
                    if ( (image_.convexHull.getElem(x_t, y_t, 0) == 1) ) {
                        int dist_t;
                        ivec2 res_t;
                        res_t.x = x_t;
                        res_t.y = y_t;
                        dist_t = dot(res_t - start, res_t - start);
                        if (dist_t > dist){
                            intersection = res_t;
                            dist = dist_t;
                            normal.x = image_.convexHull.getElem(x_t, y_t, 1);
                            normal.y = image_.convexHull.getElem(x_t, y_t, 2);
                        }
                    }
                }
                break;

            case DrawHull:
                {
                    int lineWidth = 3;

                    int min_x = max(x_t - lineWidth / 2, 0);
                    int max_x = min(x_t + lineWidth / 2, image_.width-1);
                    int min_y = max(y_t - lineWidth / 2, 0);
                    int max_y = min(y_t + lineWidth / 2, image_.height-1);

                    for (int x_new = min_x; x_new <= max_x; ++x_new) {
                        for (int y_new = min_y; y_new <= max_y; ++y_new) {
                            // hull
                            image_.convexHull.setElem(x_new, y_new, 1, 0);
                            // normal x-value
                            image_.convexHull.setElem(x_new, y_new, grad_x, 1) ;
                            // normal y-value
                            image_.convexHull.setElem(x_new, y_new, grad_y, 2);
                        }
                    }
                }
                break;

        }

        if (mirroringWH) {
            swap(x, y);
            swap(x0,y0);
        }

        ++x;
        if (d <= 0) {
            d += deltaE;
        }
        else {
            y++;
            d += deltaNE;
        }
    }

    if (action == HullIntersect)
        return (dist > -1);
    else
        return true;
}


// separates overlapping labels in SILHOUETTE! layout
// by shifting each two of them to opposite directions along hull
// both label positions and connpoints are corrected!
// @return false, if an overlap was found
bool VolumeLabeling::resolveLabelOverlaps(LabelVec& pLabels) {
    bool overlapFound = false;
    bool finished = false;

    int mask;
    int iter = 0;
    while ( (!finished) && (iter < maxIterations_.get()) ) {
        finished = true;

        // checks each label's bounding box against each other's
        for (size_t i=0; i< pLabels.size(); ++i) {
            for (size_t j=i+1; j<pLabels.size(); ++j) {
                ivec2 labelPos1 = pLabels[i].labelPos;
                ivec2 labelPos2 = pLabels[j].labelPos;
                LabelText text1 = pLabels[i].labelData->text;
                LabelText text2 = pLabels[j].labelData->text;

                // labels' bounding boxes (including gaps)
                ivec2 bl1 = labelPos1 - ivec2( tgt::iround(text1.width/2.f + gaps_.SIL_LabelLabel/2.f),
                    tgt::iround(text1.height/2.f + gaps_.SIL_LabelLabel/2.f) );
                ivec2 tr1 = labelPos1 + ivec2( tgt::iround(text1.width/2.f + gaps_.SIL_LabelLabel/2.f),
                    tgt::iround(text1.height/2.f + gaps_.SIL_LabelLabel/2.f) );
                ivec2 bl2 = labelPos2 - ivec2( tgt::iround(text2.width/2.f + gaps_.SIL_LabelLabel/2.f),
                    tgt::iround(text2.height/2.f + gaps_.SIL_LabelLabel/2.f) );
                ivec2 tr2 = labelPos2 + ivec2( tgt::iround(text2.width/2.f + gaps_.SIL_LabelLabel/2.f),
                    tgt::iround(text2.height/2.f + gaps_.SIL_LabelLabel/2.f) );

                // if bounding boxes intersect, separate labels by shifting them
                // to opposite directions parallel to hull
                if (intersectBoxBox(bl1, tr1, bl2, tr2, mask )) {
                    // interpolate hull normal from the two label normals
                    ivec3 n1 = ivec3(pLabels[i].normal, 0);
                    ivec3 n2 = ivec3(pLabels[j].normal, 0);
                    ivec3 normal = n1+n2;
                    if (normal == ivec3(0,0,0))
                        normal = cross(n1, ivec3(0,0,1));

                    // calculate tangent (= shift-vector) from normal
                    ivec3 t = cross(normal, ivec3(0,0,1));
                    ivec2 tangent = ivec2(t.x, t.y);

                    // center points of labels (= center of a label's textbox)
                    ivec2 center1 = labelPos1;
                    ivec2 center2 = labelPos2;
                    ivec2 c = center1 - center2;

                    // minimum distances between center points in both directions
                    int distx = ( tgt::iround(text1.width/2.f)
                        + tgt::iround(text2.width/2.f + gaps_.SIL_LabelLabel + 4) ) * sign(c.x) ;
                    int disty = ( tgt::iround(text1.height/2.f)
                        + tgt::iround(text2.height/2.f + gaps_.SIL_LabelLabel + 4) ) * sign(c.y);

                    // scale factor for shift vector (tangent vector)
                    float a,b,s;
                    a = ( tangent.x != 0 ? static_cast<float>(c.x-distx) / (2*tangent.x) : image_.width);
                    b = ( tangent.y != 0 ? static_cast<float>(c.y-disty) / (2*tangent.y) : image_.height);
                    s = (fabs(a) < fabs(b) ? a : b);

                    // create shift vector
                    tangent = (ivec2)((vec2)tangent*s);

                    // shift connection points to opposite directions
                    // by adding/subtracting shift vector
                    pLabels[i].labelPos -= tangent;
                    pLabels[j].labelPos += tangent;

                    // correct labels' distances to hull (shift in or out)
                    correctHullGap(pLabels[i], pLabels[i].labelPos);
                    correctHullGap(pLabels[j], pLabels[j].labelPos);

                    // recalculate labels' connpoints from new labelpositions
                    connPointFromLabelPos(pLabels[i], pLabels[i].connectionPoint, pLabels[i].normal );
                    connPointFromLabelPos(pLabels[j], pLabels[j].connectionPoint, pLabels[i].normal);

                    finished = false;
                    overlapFound = true;
                }
            }
        }
        iter++;
    }

    return !overlapFound;
}


// shifts label in or out to achieve correct distance from hull
void VolumeLabeling::correctHullGap(Label pLabel, ivec2& labelPos) {
    int gap = gaps_.SIL_LabelHull;
    float step = 2.0;       // must be greater zero
    int maxSteps = tgt::iround(max(image_.width, image_.height) / step);

    LabelText text = pLabel.labelData->text;

    ivec2 bl, tr;
    vec2 shift, labelPosNew;

    labelPosNew = (vec2)pLabel.labelPos;

    // create shift vector
    if (pLabel.anchorPoint != pLabel.connectionPoint)
        shift = (vec2)(pLabel.anchorPoint - pLabel.connectionPoint);
    else
        shift = (vec2)(pLabel.anchorPoint - pLabel.labelPos);

    shift = normalize(shift);
    shift = step*shift;

    bl = ivec2( tgt::iround(labelPosNew.x - text.width/2.f) - gap,
        tgt::iround(labelPosNew.y - text.height/2.f) - gap);
    tr = ivec2( tgt::iround(labelPosNew.x + text.width/2.f) + gap,
        tgt::iround(labelPosNew.y + text.height/2.f) + gap);

    // shift in until convex hull or anchor point are reached
    int iter = 0;
    while ( bl.x >= 0 && bl.y >= 0 && tr.x < image_.width && tr.y < image_.height &&
        !intersectBoxHull(bl,tr) && !pointInBox(bl, tr, pLabel.anchorPoint) && (iter < maxSteps) ) {
            labelPosNew += shift;
            bl = ivec2( tgt::iround(labelPosNew.x - text.width/2.f) - gap,
                tgt::iround(labelPosNew.y - text.height/2.f) - gap);
            tr = ivec2( tgt::iround(labelPosNew.x + text.width/2.f) + gap,
                tgt::iround(labelPosNew.y + text.height/2.f) + gap);
            iter++;
    }

    // shift out until convex hull has been left
    iter = 0;
    shift = -shift;
    while ( bl.x >= 0 && bl.y >= 0 && tr.x < image_.width && tr.y < image_.height &&
        ( intersectBoxHull(bl,tr) || pointInBox(bl, tr, pLabel.anchorPoint) ) && (iter < maxSteps) ) {
            labelPosNew += shift;
            bl = ivec2( tgt::iround(labelPosNew.x - text.width/2.f) - gap,
                tgt::iround(labelPosNew.y - text.height/2.f) - gap);
            tr = ivec2( tgt::iround(labelPosNew.x + text.width/2.f) + gap,
                tgt::iround(labelPosNew.y + text.height/2.f) + gap);
            iter++;
    }

    labelPos = (ivec2)labelPosNew;
}


// tests if the convex hull intersects an axis aligned 2D-box
// @bl  bottom left of the box
// @tr  top right of the box
bool VolumeLabeling::intersectBoxHull(ivec2 const &bl, ivec2 const &tr) {
    ivec2 br = ivec2(tr.x, bl.y);
    ivec2 tl = ivec2(bl.x, tr.y);
    ivec2 t1,t2;

    return ( (midPointLine(bl,br,HullIntersect,t1,t2)) ||
             (midPointLine(br,tr,HullIntersect,t1,t2)) ||
             (midPointLine(tr,tl,HullIntersect,t1,t2)) ||
             (midPointLine(tl,bl,HullIntersect,t1,t2)) );
}


// tests whether two axis-aligned 2D boxes intersect
// @param bl1 bottom-left of box 1
// @param tr1 top-right of box 1
// @param bl2 bottom-left of box 2
// @param tr2 top-right of box 2
// @param mask  bitmask for kind of intersection
//              (mask & 1000) == true  means top-left of box 2 lies in box 1
//              (mask & 0100) == true  means top-right of box 2 lies in box 1
//              (mask & 0010) == true  means bottom-right of box 2 lies in box 1
//              (mask & 0001) == true  means bottom-left of box 2 lies in box 1
// @return true if intersection found
bool VolumeLabeling::intersectBoxBox(ivec2 const& bl1, ivec2 const& tr1,
                                 ivec2 const& bl2, ivec2 const& tr2, int& mask) {
    mask = 0;

    ivec2 br2 = ivec2(tr2.x, bl2.y);
    ivec2 tl2 = ivec2(bl2.x, tr2.y);

    mask |= (((bl1.x <= tl2.x) && (tr1.x >= tl2.x) && (bl1.y <= tl2.y) && (tr1.y >= tl2.y)) ? 1000 : 0);
    mask |= (((bl1.x <= tr2.x) && (tr1.x >= tr2.x) && (bl1.y <= tr2.y) && (tr1.y >= tr2.y)) ? 0100 : 0);
    mask |= (((bl1.x <= br2.x) && (tr1.x >= br2.x) && (bl1.y <= br2.y) && (tr1.y >= br2.y)) ? 0010 : 0);
    mask |= (((bl1.x <= bl2.x) && (tr1.x >= bl2.x) && (bl1.y <= bl2.y) && (tr1.y >= bl2.y)) ? 0001 : 0);

    ivec2 br1 = ivec2(tr1.x, bl1.y);
    ivec2 tl1 = ivec2(bl1.x, tr1.y);

    mask |= (((bl2.x <= tl1.x) && (tr2.x >= tl1.x) && (bl2.y <= tl1.y) && (tr2.y >= tl1.y)) ? 1000 : 0);
    mask |= (((bl2.x <= tr1.x) && (tr2.x >= tr1.x) && (bl2.y <= tr1.y) && (tr2.y >= tr1.y)) ? 0100 : 0);
    mask |= (((bl2.x <= br1.x) && (tr2.x >= br1.x) && (bl2.y <= br1.y) && (tr2.y >= br1.y)) ? 0010 : 0);
    mask |= (((bl2.x <= bl1.x) && (tr2.x >= bl1.x) && (bl2.y <= bl1.y) && (tr2.y >= bl1.y)) ? 0001 : 0);

    return (mask != 0);
}

// returns true, if point lies inside the box specified by bottem-left bl and top-right tr
bool VolumeLabeling::pointInBox(ivec2 const &bl, ivec2 const &tr, ivec2 const &point) {
    return ( (point.x >= bl.x) && (point.x <= tr.x) && (point.y >= bl.y) && (point.y <= tr.y) );
}

bool VolumeLabeling::findLabelPathBest(Label &pLabel) {
    // step size = gap between control points in pixels
    float radius = 10.f;

    //
    // find several label paths for different start value configurations
    // and choose the best of them (the one with the greatest average distance map height)
    //

    int x = pLabel.anchorPoint.x;
    int y = pLabel.anchorPoint.y;
    float firstValue = 0.f, secondValue = 0.f;
    float firstAngle = 0.f, secondAngle = 0.f;

    float hValue = static_cast<float>(image_.horzDistance.getElem(x,y));
    float vValue = static_cast<float>(image_.vertDistance.getElem(x,y));
    float adValue = static_cast<float>(image_.ascDiagDistance.getElem(x,y));
    float ddValue = static_cast<float>(image_.descDiagDistance.getElem(x,y));

    // get initial direction
    if ( (hValue > vValue) && (hValue > adValue) && (hValue > ddValue) ) {
        firstValue = hValue;
        firstAngle = 0.f;
    }
    else if ( (vValue > adValue) && (vValue > ddValue) ) {
        firstValue = vValue;
        firstAngle = 90.f;
    }
    else if ( (adValue > ddValue)) {
        firstValue = adValue;
        firstAngle = 45.f;
    }
    else {
        firstValue = ddValue;
        firstAngle = 135.f;
    }

    if (firstAngle == 0.f) {
        secondAngle = (adValue > ddValue ? 45.f : -45.f);
        secondValue = max(adValue, ddValue); //(adValue > ddValue ? adValue : ddValue);
    }
    else if (firstAngle == 45.f) {
        secondAngle = (vValue > hValue ? 90.f : 0.f);
        secondValue = max(vValue, hValue); //(vValue > hValue ? vValue : hValue);
    }
    else if (firstAngle == 90.f) {
        secondAngle = (ddValue > adValue ? 135.f : 45.f);
        secondValue = max(ddValue, adValue); //(ddValue > adValue ? ddValue : adValue);
    }
    else if (firstAngle == 135.f) {
        secondAngle = (hValue > vValue ? 180.f : 90.f);
        secondValue = max(hValue, vValue); //(hValue > vValue ? hValue : vValue);
    }
    float angle = (firstAngle*firstValue + secondAngle*secondValue) / (firstValue + secondValue);
    if (angle < 0.f)
        angle += 360.f;

    float avgDist;
    int maxIndex=-1;
    float maxAvgDist=0.f;

    vector<ivec2> labelPath[6];

    if ( findLabelPathOne(pLabel, radius, angle,
        45.f, 15.f, 90.f, static_cast<float>(pLabel.labelData->text.width), labelPath[0]) ) {
            avgDist = 0;
            for (size_t point=0; point<labelPath[0].size(); point++)
                avgDist += image_.distanceField.getElem(labelPath[0][point].x, labelPath[0][point].y);

            avgDist /= labelPath[0].size();
            maxAvgDist = avgDist;
            maxIndex = 0;
        }

        if ( findLabelPathOne(pLabel, radius, firstAngle,
             45.f, 15.f, 90.f, static_cast<float>(pLabel.labelData->text.width), labelPath[1]) ) {
                avgDist = 0;
                for (size_t point=0; point<labelPath[1].size(); point++)
                    avgDist += image_.distanceField.getElem( labelPath[1][point].x, labelPath[1][point].y );

                avgDist /= labelPath[1].size();
                if (avgDist > maxAvgDist ){
                    maxIndex = 1;
                    maxAvgDist = avgDist;
                }
        }

        if ( findLabelPathOne(pLabel, radius, secondAngle,
             45.f, 15.f, 90.f, static_cast<float>(pLabel.labelData->text.width), labelPath[2]) ) {
                avgDist = 0;
                for (size_t point=0; point<labelPath[2].size(); point++)
                    avgDist += image_.distanceField.getElem( labelPath[2][point].x, labelPath[2][point].y );

                    avgDist /= labelPath[2].size();
                if (avgDist > maxAvgDist ) {
                    maxIndex = 2;
                    maxAvgDist = avgDist;
                }
        }

        if ( findLabelPathOne(pLabel, radius, angle,
             45.f, 15.f, 0.f, static_cast<float>(pLabel.labelData->text.width), labelPath[3]) ) {
                avgDist = 0;
                for (size_t point=0; point<labelPath[3].size(); point++)
                    avgDist += image_.distanceField.getElem( labelPath[3][point].x, labelPath[3][point].y );

                avgDist /= labelPath[3].size();
                if (avgDist > maxAvgDist ){
                    maxIndex = 3;
                    maxAvgDist = avgDist;
                }
        }

        if ( findLabelPathOne(pLabel, radius*1.5f, firstAngle,
            45.f, 15.f, 90.f, static_cast<float>(pLabel.labelData->text.width), labelPath[4]) ) {
                avgDist = 0;
                for (size_t point=0; point<labelPath[4].size(); point++)
                    avgDist += image_.distanceField.getElem( labelPath[4][point].x, labelPath[4][point].y );

                avgDist /= labelPath[4].size();
                if (avgDist > maxAvgDist ){
                    maxIndex = 4;
                    maxAvgDist = avgDist;
                }
        }

        if ( findLabelPathOne(pLabel, radius*2.f, secondAngle,
            45.f, 15.f, 90.f, static_cast<float>(pLabel.labelData->text.width), labelPath[5]) ) {
                avgDist = 0;
                for (size_t point=0; point<labelPath[5].size(); point++)
                    avgDist += image_.distanceField.getElem( labelPath[5][point].x, labelPath[5][point].y );

                avgDist /= labelPath[5].size();
                if (avgDist > maxAvgDist ){
                    maxIndex = 5;
                    maxAvgDist = avgDist;
                }
        }

        if (maxIndex > -1) {
            pLabel.controlPoints = labelPath[maxIndex];
            //float pathLength = 0;
            //for (size_t point=0; point<pLabel.controlPoints.size()-1; ++point) {
            //    pathLength += sqrtf(powf(pLabel.controlPoints[point+1].x -
            //                             static_cast<float>(pLabel.controlPoints[point].x), 2.f) +
            //                        powf(static_cast<float>(pLabel.controlPoints[point+1].y) -
            //                             pLabel.controlPoints[point].y, 2.f));
            //}
            pLabel.intern = pLabel.labelData->internPreferred;

            // invert control points if curve is running from right to left
            float xsum = 0.f;
            for (size_t j = 0; j < pLabel.controlPoints.size()-1; ++j)
                xsum += pLabel.controlPoints[j+1].x - pLabel.controlPoints[j].x;

            if (xsum < 0.f)
                for (size_t j = 0; j < (pLabel.controlPoints.size() / 2); ++j)
                    swap(pLabel.controlPoints[j], pLabel.controlPoints[pLabel.controlPoints.size()-1-j]);

            // calculate 2D fitting curve from control points in image space
            Curve2DPolynomial* curve2D = new Curve2DPolynomial(polynomialDegree_.get());
            vector<vec2> ctrlPoints2DCurve;
            for (size_t p=0; p < pLabel.controlPoints.size(); ++p) {
                vec2 ctrlPoint = vec2(static_cast<float>(pLabel.controlPoints[p].x), static_cast<float>(pLabel.controlPoints[p].y));
                ctrlPoints2DCurve.push_back(ctrlPoint);
            }
            if ( curve2D->setCtrlPoints(ctrlPoints2DCurve,
                static_cast<float>(pLabel.labelData->text.widthIntern)) ) {
                pLabel.curve2D = curve2D;
                return true;
            }
            else {
                delete curve2D;
                return false;
            }

        }
        else
            return false;
}


VolumeLabeling::Circle VolumeLabeling::generateCircle(float radius) {
    Circle circle;
    CirclePixel pixel;

    int d = static_cast<int>(1 - radius);
    int ddE = 3;
    int ddSE = static_cast<int>(5 - 2*radius);
    int x = 0;
    int y = static_cast<int>(radius);
    while (y >= x) {
        float angle = rad2deg(acosf(x / radius));

        // x,y
        pixel = CirclePixel();
        pixel.coords = ivec2(x,y);
        pixel.angle = angle;
        circle.pixels.push_back(pixel);

        // y,x
        pixel = CirclePixel();
        pixel.coords = ivec2(y,x);
        pixel.angle = 90.f-angle;
        circle.pixels.push_back(pixel);

        // x, -y
        pixel = CirclePixel();
        pixel.coords = ivec2(x,-y);
        pixel.angle = 360.f-angle;
        circle.pixels.push_back(pixel);

        // y, -x
        pixel = CirclePixel();
        pixel.coords = ivec2(y,-x);
        pixel.angle = 270.f+angle;
        circle.pixels.push_back(pixel);

        // -x, -y
        pixel = CirclePixel();
        pixel.coords = ivec2(-x,-y);
        pixel.angle = 180.f+angle;
        circle.pixels.push_back(pixel);

        // -y, -x
        pixel = CirclePixel();
        pixel.coords = ivec2(-y,-x);
        pixel.angle = 270.f-angle;
        circle.pixels.push_back(pixel);

        // -x, y
        pixel = CirclePixel();
        pixel.coords = ivec2(-x,y);
        pixel.angle = 180.f-angle;
        circle.pixels.push_back(pixel);

        // -y, x
        pixel = CirclePixel();
        pixel.coords = ivec2(-y,x);
        pixel.angle = 90.f+angle;
        circle.pixels.push_back(pixel);

        x++;
        // East
        if (d < 0) {
            d += ddE;
            ddE += 2;
            ddSE += 2;
        // SouthEast
        }
        else {
            y--;
            d += ddSE;
            ddE += 2;
            ddSE += 4;
        }
    }


    circle.radius = radius;
    circle.numPixels = static_cast<int>(circle.pixels.size());

    for (size_t i = 0 ; i < circle.pixels.size()-1 ; ++i) {
        for (size_t j = i+1 ; j < circle.pixels.size() ; ++j) {
            if (circle.pixels[i].angle > circle.pixels[j].angle) {
                CirclePixel tmp = circle.pixels[i];
                circle.pixels[i] = circle.pixels[j];
                circle.pixels[j] = tmp;
            }
        }
    }
    for (int i = 0 ; i < circle.numPixels ; ++i) {
        CirclePixel pixel;
        pixel.coords = circle.pixels[i].coords;
        pixel.angle = circle.pixels[i].angle + 360.f;
        circle.pixels.push_back(pixel);
    }

    return circle;
}

VolumeLabeling::CirclePixel VolumeLabeling::getNextCirclePixel(int id, Circle circle, ivec2 curPoint,
                                                   float minAngle, float maxAngle)
{
    int minIndex = max(0, iround( ((minAngle*0.9)/360.f)*circle.numPixels ) - 10);
        while (circle.pixels[minIndex].angle < minAngle)
            ++minIndex;

        int maxIndex = min(static_cast<int>(circle.numPixels*2-1), iround( ((maxAngle*1.1)/360.f)*circle.numPixels ) + 10);
        while (circle.pixels[maxIndex].angle > maxAngle && maxIndex > minIndex)
            --maxIndex;

        if (minIndex == maxIndex && minIndex > 0) {
            if ( (circle.pixels[minIndex].angle-circle.pixels[minIndex-1].angle / 2.f) > maxAngle ) {
                --minIndex;
                --maxIndex;
            }
        }
        // FIXME: Still neccessary? (ab)
        else if (minIndex > maxIndex) {
           // std::cout << "Fehler: minIndex > maxIndex" << std::endl;
        }
        else {
           // std::cout << "Pixels: " << maxIndex-minIndex << ", ";
        }

        CirclePixel candidate;
        candidate.coords.x = 0;
        candidate.coords.y = 0;
        int distValue = -1;
        int x,y;
        for (int i=minIndex; i<=maxIndex; ++i) {
            x = curPoint.x + circle.pixels[i].coords.x;
            y = curPoint.y + circle.pixels[i].coords.y;
            if (x >= 0 && x < image_.width && y >= 0 && y < image_.height &&
                image_.idBuffer.getElem(x,y) == id &&
                image_.distanceField.getElem(x, y) > distValue) {
                candidate = circle.pixels[i];
                distValue = image_.distanceField.getElem(x, y);
            }
        }
        return candidate;
}

bool VolumeLabeling::findLabelPathOne(const Label &pLabel, float step, float initAngle, float initMaxStepAngle,
                               float maxStepAngle, float maxAngle, float length, vector<ivec2> &result)
{
    float threshold = 0.3f;
    Circle circle = generateCircle(step);

    int id = pLabel.labelData->id;
    result.push_back(pLabel.anchorPoint);

    float curAngle = initAngle;
    float firstAngle=-1.f;
    ivec2 curPoint = pLabel.anchorPoint;
    int initValue = image_.distanceField.getElem(curPoint.x, curPoint.y);
    int distValue = initValue;
    float curMinAngle;
    float curMaxAngle;

    if (initAngle-maxAngle < 0.f) {
        initAngle += 360.f;
        curAngle += 360.f;
    }

    // find inital direction
    curMinAngle = max(curAngle-initMaxStepAngle, initAngle-maxAngle);
    curMaxAngle = min(curAngle+initMaxStepAngle, initAngle+maxAngle);

    CirclePixel nextPixel = getNextCirclePixel(id, circle, curPoint, curMinAngle, curMaxAngle);

    if (nextPixel.coords.x != 0 || nextPixel.coords.y != 0) {
        ivec2 nextPoint = curPoint + nextPixel.coords;
        result.push_back(nextPoint);
        curPoint = nextPoint;
        curAngle = nextPixel.angle;
        firstAngle = nextPixel.angle;
        distValue = image_.distanceField.getElem( nextPoint.x, nextPoint.y );
    }
    else
        distValue = -1;

    // find path to the right
    while (distValue > initValue*threshold) {
        curMinAngle = max(curAngle-maxStepAngle, initAngle-maxAngle);
        curMaxAngle = min(curAngle+maxStepAngle, initAngle+maxAngle);

        CirclePixel nextPixel = getNextCirclePixel(id, circle, curPoint, curMinAngle, curMaxAngle);

        if (nextPixel.coords.x != 0 || nextPixel.coords.y != 0) {
            ivec2 nextPoint = curPoint + nextPixel.coords;
            result.push_back(nextPoint);
            curPoint = nextPoint;
            curAngle = nextPixel.angle;
            distValue = image_.distanceField.getElem(nextPoint.x, nextPoint.y);
        }
        else
            distValue = -1;
    }

    // find path to the left
    distValue = initValue;
    curAngle = firstAngle + 180.f;
    if (curAngle > 360.f) curAngle -= 360.f;
    initAngle = curAngle;
    curPoint = pLabel.anchorPoint;
    int anchorIndex = 0;
    while (distValue > initValue*threshold) {
        curMinAngle = max(curAngle-maxStepAngle, initAngle-maxAngle);
        curMaxAngle = min(curAngle+maxStepAngle, initAngle+maxAngle);

        CirclePixel nextPixel = getNextCirclePixel(id, circle, curPoint, curMinAngle, curMaxAngle);

        if (nextPixel.coords.x != 0 || nextPixel.coords.y != 0) {
            ivec2 nextPoint = curPoint + nextPixel.coords;
            result.insert(result.begin(), nextPoint);
            ++anchorIndex;
            curPoint = nextPoint;
            curAngle = nextPixel.angle;
            distValue = image_.distanceField.getElem(nextPoint.x, nextPoint.y);
        }
        else
            distValue = -1;
    }

    // determine, how many control points are needed for the
    // specified path length, choose the ones with the highest
    // distance values and discard the rest (beginning at the anchor)
    if (length > 0.f) {
        int neededPoints = static_cast<int>( length / step + 1.f );
        if ( neededPoints > static_cast<int>(result.size()) )
            return false;
        else {
            neededPoints = min( neededPoints+2, static_cast<int>(result.size()) );
            int l = anchorIndex;
            int r = anchorIndex;
            int lvalue, rvalue;
            while (r-l+1 < neededPoints) {
                lvalue = (l > 0 ? image_.distanceField.getElem( result[l-1].x, result[l-1].y ) : -1);
                rvalue = (r < static_cast<int>(result.size())-1 ? image_.distanceField.getElem( result[r+1].x, result[r+1].y ) : -1);
                if (lvalue > rvalue)
                    --l;
                else
                    ++r;
            }
            if (r < static_cast<int>(result.size())-1)
                result.erase(result.begin()+r+1, result.end());
            if (l > 0)
                result.erase(result.begin(), result.begin()+l);
        }
    }

    return true;
}


// Calculates bezier control points for internal labels.
// In 3D mode the bezier patch is fit to its segment's 3D surface by using
// first-hit-positions.
// In 2D the internal labeling is performed in image space, without using
// any 3D shape information.
bool VolumeLabeling::findBezierPoints(Label &pLabel) {
    // -------------------------------------------------
    // 3D shape fitting (world space)
    // -------------------------------------------------
    if (shape3D_.get()) {

        int hDegree = bezierHorzDegree_.get();
        // round vertical degree to next even integral value.
        // this eases the following calculations.
        int vDegree = 2*((bezierVertDegree_.get()+1)/2);
        const int samplingRate = 3;
        const int samplings2D = 10;

        string text = pLabel.labelData->text.text;

        Curve2DPolynomial curve2D = *(pLabel.curve2D);

        // now sample the 2D curve, transform the sample points to object space using
        // first-hit-positions, and fit a 3D polynomial to the transformed control points
        Curve3DPolynomial curve(polynomialDegree_.get());
        vector<vec3> ctrlPoints;
        float paramOffset = 1.f / (samplings2D-1);
        for (int i=0; i < samplings2D; ++i) {
            vec2 ctrlpoint = curve2D.getCurvePoint( i*paramOffset);
            vec3 firstHitPoint = image_.firstHitPositions.getPixel3Ch(tgt::iround(ctrlpoint.x), tgt::iround(ctrlpoint.y));
            ctrlPoints.push_back(firstHitPoint);
        }

        // anchor-point in 2D screen coords
        ivec2 anchor2D = (!pLabel.offLabel ? pLabel.anchorPoint : ivec2(pLabel.curve2D->getCurvePoint(0.5f)));
        // anchor-point in object space
        vec3 anchorObjectSpace = image_.firstHitPositions.getPixel3Ch(anchor2D.x, anchor2D.y);
        // anchor-point in 3D screen coords
        vec3 anchorScreen = image_.firstHitPositions.projectToViewport(anchorObjectSpace);

        // calculate pixelScale: if the text would be put on a screen-parallel plane
        // "through" the anchor point in object-space, pixelScale is the distance in
        // object-coords that is projected to one pixel.
        // in other words: if you project two points lieing in this plane, whose distance is
        // pixelScale, the distance of their projections on the screen is one pixel
        vec3 tempVec1 = image_.firstHitPositions.inverseProject(anchorScreen);
        vec3 tempVec2 = image_.firstHitPositions.inverseProject(anchorScreen+vec3(0.f, 10.f, 0.f));
        float pixelScale = length( (tempVec2 - tempVec1) ) / 10.f;

        // h_offset is the gap between two bezier control point in horizontal direction
        // in object-space
        float h_offset = pLabel.labelData->text.widthIntern*pixelScale / hDegree;

        // calculate polynomial fitting function
        // the desired length in object space is provided -> curve internally calculates
        // start and end params from it
        if (!curve.setCtrlPoints(ctrlPoints, pLabel.labelData->text.widthIntern*pixelScale))
            return false;

        //
        // now calculate bezier points by sampling the fitting function
        //
        vec3* bezierPoints = new vec3[(hDegree+1)*(vDegree+1)];
        float param = 0.f;
        for (int k=0; k<=hDegree; ++k) {
            // sample curve
            vec3 curvePoint = curve.getCurvePoint(param);
            vec3 curvePointVP = image_.firstHitPositions.projectToViewport(curvePoint);
            if (curvePointVP.x < 0 || curvePointVP.x >= image_.width ||
                curvePointVP.y < 0 || curvePointVP.y >= image_.height ||
                image_.idImage.getElem(tgt::iround(curvePointVP.x), tgt::iround(curvePointVP.y)) == 0) {
                    delete[] bezierPoints;
                    return false;
            }

            // find firstHitPoint corresponding to curve sample point
            vec3 curSurfacePoint = image_.firstHitPositions.getPixel3Ch(tgt::iround(curvePointVP.x), tgt::iround(curvePointVP.y));
            vec3 curSurfacePointVP = image_.firstHitPositions.projectToViewport(curSurfacePoint);
            if (curSurfacePointVP.x < 0 || curSurfacePointVP.x >= image_.width ||
                curSurfacePointVP.y < 0 || curSurfacePointVP.y >= image_.height ||
                image_.idImage.getElem(tgt::iround(curSurfacePointVP.x), tgt::iround(curSurfacePointVP.y)) == 0  ) {
                    delete[] bezierPoints;
                    return false;
            }
            bezierPoints[(hDegree+1)*(vDegree/2)+k] = curSurfacePoint;

            // curve's tangent at parameter param
            vec3 tangent = normalize( curve.getTangent( param ) );

            float v_offset = pLabel.labelData->text.heightIntern*pixelScale / (vDegree*samplingRate);

            vec3 surfacePointMiddle = curSurfacePoint;
            vec3 surfacePointMiddleVP = curSurfacePointVP;

            // go from baseline up
            // vec3 curSurfacePoint = bezierPoint;
            // vec3 curSurfacePointVP = curvePointVP;
            for (int l=1; l<=vDegree*samplingRate/2; ++l) {
                // determine the normal at the surface point from the first hit positions
                // of the neighbourhood
                vec3 normal = image_.firstHitPositions.calcNormal( tgt::iround(curSurfacePointVP.x), tgt::iround(curSurfacePointVP.y) );
                vec3 up = normalize( cross(normal, tangent) )*v_offset;

                curSurfacePoint += up;
                curSurfacePointVP = image_.firstHitPositions.projectToViewport(curSurfacePoint);
                if (curSurfacePointVP.x < 0 || curSurfacePointVP.x >= image_.width ||
                    curSurfacePointVP.y < 0 || curSurfacePointVP.y >= image_.height ||
                    //image_.idBuffer.getElem(round(curSurfacePointVP.x), round(curSurfacePointVP.y)) != segmentID ){
                    image_.idImage.getElem(tgt::iround(curSurfacePointVP.x), tgt::iround(curSurfacePointVP.y)) == 0 ) {
                        delete[] bezierPoints;
                        return false;
                    }

                    if (l % samplingRate == 0)
                        bezierPoints[(hDegree+1)*(vDegree/2+l/samplingRate) + k] = curSurfacePoint;
            }

            // go from baseline down
            curSurfacePoint = surfacePointMiddle;
            curSurfacePointVP = surfacePointMiddleVP;
            for (int l=1; l<=vDegree*samplingRate/2; ++l) {
                // determine the normal at the surface point from the first hit positions
                // of the neighbourhood
                vec3 normal = image_.firstHitPositions.calcNormal( tgt::iround(curSurfacePointVP.x), tgt::iround(curSurfacePointVP.y) );
                vec3 up = normalize( cross(normal, tangent) )*v_offset;

                curSurfacePoint -= up;
                curSurfacePointVP = image_.firstHitPositions.projectToViewport(curSurfacePoint);
                if (curSurfacePointVP.x < 0 || curSurfacePointVP.x >= image_.width ||
                    curSurfacePointVP.y < 0 || curSurfacePointVP.y >= image_.height ||
                    image_.idImage.getElem(tgt::iround(curSurfacePointVP.x), tgt::iround(curSurfacePointVP.y)) == 0 ) {
                        delete[] bezierPoints;
                        return false;
                    }

                    if (l % samplingRate == 0)
                        bezierPoints[(hDegree+1)*(vDegree/2 - l/samplingRate) + k] = curSurfacePoint;
            }
            curve.getNextPoint(param, h_offset);
        }


        // set bezier control points
        pLabel.bezierPatch.setCtrlPoints(bezierPoints, hDegree, vDegree);

        // finally check calculated patch
        if ( !checkPatchQuality(pLabel) ) {
            //delete bezierPoints;
            return false;
        }


        // --------------------------------------------------
        // no 3D shape fitting (labeling only in image space)
        // --------------------------------------------------
    }
    else if (!shape3D_.get()) {
        // horizontal degree of bezier patch
        // vertical degree is always 1 here (patch is flat)
        int id = pLabel.labelData->id;
        int hDegree = bezierHorzDegree_.get();

        Curve2DPolynomial curve2D = *(pLabel.curve2D);
        vec3* bezierPoints = new vec3[(hDegree+1)*2];
        vec3 normal = vec3(0.f,0.f,1.f);

        float hOffset = static_cast<float>(pLabel.labelData->text.widthIntern) / (hDegree);
        float vOffset = static_cast<float>(pLabel.labelData->text.heightIntern) / 2.f;

        // now sample the 2D curve and calculate bezier control points in image space
        float param = 0.f;
        for (int i=0; i <= hDegree; ++i) {
            vec2 samplepoint = curve2D.getCurvePoint( param );
            vec3 samplepoint3D = vec3(samplepoint.x, samplepoint.y, 0.f);
            vec2 tangent = curve2D.getTangent( param );
            vec3 tangent3D = vec3(tangent.x, tangent.y, 0);
            vec3 updir = normalize(cross(normal, tangent3D));

            vec3 bezierpoint = samplepoint3D + updir*vOffset;
            if ( image_.idBuffer.getElem(tgt::iround(bezierpoint.x), tgt::iround(bezierpoint.y)) != id ) {
                delete[] bezierPoints;
                return false;
            }
            bezierPoints[(hDegree+1)*1 + i] = bezierpoint;

            bezierpoint = samplepoint3D - updir*vOffset;
            if ( image_.idBuffer.getElem(tgt::iround(bezierpoint.x), tgt::iround(bezierpoint.y)) != id ) {
                delete[] bezierPoints;
                return false;
            }
            bezierPoints[(hDegree+1)*0 + i] = bezierpoint;

            curve2D.getNextPoint(param, hOffset);
        }

        // set calculated control points
        pLabel.bezierPatch.setCtrlPoints(bezierPoints, hDegree, 1);
        delete[] bezierPoints;
    }
    return true;
}

// Returns false, if the bezier patch is too distorted
// (and thus not readable)
bool VolumeLabeling::checkPatchQuality(Label pLabel) {
    const float max_deviation = 16.f;

    // check if viewport-projection of bezier patch closely matches 2D-labelcurve
    vec2 leftPointCurve = pLabel.curve2D->getCurvePoint( 0.f );
    vec2 middlePointCurve = pLabel.curve2D->getCurvePoint( 0.5f  );
    vec2 rightPointCurve = pLabel.curve2D->getCurvePoint( 1.f );
    vec2 leftPointPatch = vec2( image_.firstHitPositions.projectToViewport(pLabel.bezierPatch.getPoint(0.f, 0.5f)).x ,
        image_.firstHitPositions.projectToViewport(pLabel.bezierPatch.getPoint(0.f, 0.5f)).y );
    vec2 middlePointPatch = vec2( image_.firstHitPositions.projectToViewport(pLabel.bezierPatch.getPoint(0.5f, 0.5f)).x ,
        image_.firstHitPositions.projectToViewport(pLabel.bezierPatch.getPoint(0.5f, 0.5f)).y );
    vec2 rightPointPatch = vec2( image_.firstHitPositions.projectToViewport(pLabel.bezierPatch.getPoint(1.f, 0.5f)).x ,
        image_.firstHitPositions.projectToViewport(pLabel.bezierPatch.getPoint(1.f, 0.5f)).y );

    float deviation = length(leftPointCurve - leftPointPatch) / 2.f +
                      length(middlePointCurve - middlePointPatch) +
                      length(rightPointCurve - rightPointPatch) / 2.f;

    //std::cout << deviation << std::endl;

    return (deviation <= max_deviation);
}


void VolumeLabeling::placeOffLabels(LabelVec &offLabels) {
    const int H_GAP = 10;
    const int V_GAP = 8;

    for (size_t i=0; i<offLabels.size(); ++i) {
        Label* newLabel = &(offLabels[i]);
        LabelData* labelData = newLabel->labelData;
        // FIXME(apv): WHAT THE FUCK?
        newLabel->labelPos = ivec2(labelData->text.width/2 + H_GAP,
            V_GAP + fontExtern_.get()->getFontSize()/2 + static_cast<int>(i)*(fontExtern_.get()->getFontSize()+V_GAP));
    }
}

void VolumeLabeling::addUnsegmentedLabelData(std::string text) {
    std::ostringstream os;

    int maxID = -1;
    bool found = false;
    for (size_t i=0; i<labelPersistentData_.size() && !found; ++i) {
        maxID = (labelPersistentData_[i]->id > maxID) ? labelPersistentData_[i]->id : maxID;
        found = (labelPersistentData_[i]->text.text == text);
    }
    if (found)
        return;

    maxID++;
    LabelData* labelData = new LabelData;
    labelData->text.text = text;
    labelData->xmlNode = NULL;
    labelData->id = maxID;
    labelData->internPreferred = true;
    labelData->belongsToSegment = false;
    os.clear();
    os << maxID;
    labelData->idstr = os.str();

    labelPersistentData_.push_back(labelData);

    // register segment as picking object
    //TODO: tc_
    //labelingPort_.getTextureContainer()->getIDManager().addNewPickObj(labelData->idstr);

    //std::cout << labelData->idstr<< std::endl;
}

void VolumeLabeling::removeUnsegmentedLabelData(std::string text) {
    int dataIndex = -1;
    for (size_t i = 0; (i < labelPersistentData_.size()) && (dataIndex == -1); ++i) {
        if (!labelPersistentData_[i]->belongsToSegment && labelPersistentData_[i]->text.text == text)
            dataIndex = static_cast<int>(i);
    }
    if (dataIndex == -1)
        return;

    int labelIndex = -1;
    for (size_t i = 0; (i < labels_.size()) && (labelIndex == -1); ++i) {
        if (labels_[i].labelData->text.text == text)
            labelIndex = static_cast<int>(i);
    }
    if (labelIndex > -1)
        labels_.erase(labels_.begin() + labelIndex);

    delete labelPersistentData_[dataIndex];
    labelPersistentData_.erase(labelPersistentData_.begin()+dataIndex);
}

// rotates a label's 2D curve by 'angle' degrees
void VolumeLabeling::rotateLabel(Label& pLabel, float angle) {
    const int NUM_SAMPLES = 10;

    float angleRad = deg2rad(angle);

    // find center of labels 2D curve
    vec2 curveCenter = vec2(0.f);
    float paramOffset = 1.f / (NUM_SAMPLES - 1);
    for (int i=0; i<NUM_SAMPLES; ++i) {
        curveCenter += pLabel.curve2D->getCurvePoint( i*paramOffset );
    }
    curveCenter = curveCenter * (1.f/NUM_SAMPLES);

    // rotate curve's control points
    mat2 rotationMatrix(cosf(angleRad), -sinf(angleRad), sinf(angleRad), cosf(angleRad));
    vector<vec2> ctrlPoints = pLabel.curve2D->getCtrlPoints();
    for (size_t i=0; i<ctrlPoints.size(); ++i) {
        ctrlPoints[i] = curveCenter + rotationMatrix*(ctrlPoints[i]-curveCenter);
    }
    pLabel.curve2D->setCtrlPoints(ctrlPoints, static_cast<float>(pLabel.labelData->text.width));
    pLabel.rotAngle += angle;
}

// reads result of id-raycasting into image
// - idBuffer: color coded segments (red channel)
// - idImage: greater zero, where intensity of volume between thresholds (=> volume visible)
//            (necessary for silhouette calculation, blue channel)
// - distanceField: distance field internally used by find anchors (can be used for debugging)
void VolumeLabeling::readImage() {
    image_.width = labelingPort_.getSize().x;
    image_.height = labelingPort_.getSize().y;

   // glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    //GLushort* buffer = new GLushort[image_.width*image_.height*4];
    //GLfloat* buffer = idMapPort_.getTargetAsFloats();
    //idMapPort_.getColorTexture()->downloadTexture();
    //GLfloat* buffer = (GLfloat*)idMapPort_.getColorTexture()->getPixelData();

    tgt::Texture* idTexture = idMapPort_.getColorTexture();
    idTexture->downloadTexture();

    int* idImage = new int[image_.width*image_.height];
    int* idBuffer = new int[image_.width*image_.height];
    float* positionBuffer = new float[image_.width*image_.height*3];
//    image_.depthField = new GLfloat[(image_.width+2) * (image_.height+2)];
    for (int y=0; y<image_.height; ++y) {
        for (int x=0; x<image_.width; ++x) {
            //int value = tgt::iround( buffer[4*i+3]*255);
            tgt::vec4 texel = idTexture->texelAsFloat(x,y);
            int value = tgt::iround(texel.a*255.f);
            int index = y*image_.width+x;
            if (value == 255)
                idBuffer[index] = 0;
            else
                idBuffer[index] = value;
            if ( value > 0 )
                idImage[index] = 255;
            else
                idImage[index] = 0;

            positionBuffer[3*index] = static_cast<float>(texel.x);
            positionBuffer[3*index+1] = static_cast<float>(texel.y);
            positionBuffer[3*index+2] = static_cast<float>(texel.z);
        }
    }

    image_.idBuffer.setData(idBuffer,image_.width, image_.height, 1);
    image_.idImage.setData(idImage, image_.width, image_.height, 1);

    image_.firstHitPositions.setData(positionBuffer, image_.width, image_.height, 3);
    image_.firstHitPositions.setVolumeSize(currentVolume_->getCubeSize() / 2.f);
    image_.firstHitPositions.calcTransformationMatrix(camera_.get().getViewMatrix(),
        camera_.get().getProjectionMatrix(labelingPort_.getSize()),
        ivec2(image_.width, image_.height));

    idMapPort_.getColorTexture()->setCpuTextureData(0,false);
}

// compare-function for graham-scan:
// P2 < P3 <=> polar angle of P2 is smaller than polar angle of P3
// (with regard to the origin)
bool angleComp(ivec2 P2, ivec2 P3) {
    return ( (P2.x*P3.y - P3.x*P2.y) > 0 );
}

// calculates convex hull according to prior found silhouette points
// Algorithm: Graham-Scan
void VolumeLabeling::calcConvexHull() {
    image_.hullPoints.clear();
    if (image_.silPoints.size() <= 3)
        return;

    // find most left point
    int index = 0;
    for (size_t i = 1; i < image_.silPoints.size(); ++i) {
        if ( (image_.silPoints[i].x < image_.silPoints[index].x) ||
             ((image_.silPoints[i].x == image_.silPoints[index].x) &&
             (image_.silPoints[i].y < image_.silPoints[index].y)) )
                index = static_cast<int>(i);
    }
    swap(image_.silPoints[0], image_.silPoints[index]);

    // sort by polar angle (ascending) with regard to the most left point
    ivec2 shift = image_.silPoints[0];
    for (size_t i=0; i<image_.silPoints.size(); ++i)
        image_.silPoints[i] += -shift;

    std::sort(image_.silPoints.begin(), image_.silPoints.end(), angleComp);
    for (size_t i=0; i<image_.silPoints.size(); ++i)
        image_.silPoints[i] += shift;

    // calculate convex hull points
    bool popped;
    image_.hullPoints.push_back(image_.silPoints[0]);
    image_.hullPoints.push_back(image_.silPoints[1]);
    ivec2 P2;
    ivec2 P3;
    for (size_t i = 2; i < image_.silPoints.size(); ++i) {
        popped = true;
        ivec2 P1 = image_.silPoints[i];
        while (popped && (image_.hullPoints.size() > 2)) {
            P2 = image_.hullPoints[image_.hullPoints.size()-1];
            image_.hullPoints.pop_back();
            P3 = image_.hullPoints[image_.hullPoints.size()-1];
            int select = (P2.x - P1.x)*(P3.y - P1.y) - (P3.x-P1.x)*(P2.y-P1.y);
            if (select < 0) {
                image_.hullPoints.push_back(P2);
                popped = false;
            }
        }
        image_.hullPoints.push_back(P1);
    }

    // draw convex hull
    image_.convexHull.setDim(image_.width, image_.height, 3);
    ivec2 inter, norm;
    for (size_t i = 0; i < image_.hullPoints.size()-1; ++i) {
        midPointLine(image_.hullPoints[i], image_.hullPoints[i+1], DrawHull, inter, norm);
    }
    midPointLine(image_.hullPoints[image_.hullPoints.size()-1], image_.hullPoints[0], DrawHull,  inter, norm);

}

// returns the label positioned at (x,y),
// or 0 if there is none
VolumeLabeling::Label* VolumeLabeling::getPickedLabel(int /*x*/, int /*y*/) {
    Label* pLabel = 0;

    for (size_t i = 0; (i < labels_.size()) && (pLabel == 0); ++i) {
        //TODO: tc_
        //if (labelingPort_.getTextureContainer()->getIDManager().isClicked(labels_[i].labelData->idstr, x, y)){
            //pLabel = &(labels_[i]);
        //}
    }
/*    for (size_t i=0; i<invisibleLabels_.size() && pLabel == 0; i++){
        if (idManager_.isClicked(invisibleLabels_[i].labelData->idstr, x, y)){
            pLabel = &(invisibleLabels_[i]);
        }
    } */
    return pLabel;
}

// Returns, if the user has moved a label into a segment
bool VolumeLabeling::catchedBySegment(const Label& pLabel, ivec2 mousePos) {
    //ivec2 bl = pLabel.labelPos - ivec2(gaps_.SIL_LabelLabel);
    //ivec2 tr = pLabel.labelPos + ivec2(pLabel.labelData->text.width, pLabel.labelData->text.height)
                    //+ ivec2(gaps_.SIL_LabelLabel);
    if ((mousePos.x >= 0) && (mousePos.x < image_.width) && (mousePos.y >= 0) && (mousePos.y < image_.height)
        && (image_.idBuffer.getElem(mousePos.x, mousePos.y) == pLabel.labelData->id) )
        return true;
    else
        return false;
}

// updates caption of the segment belonging to pLabel.
// update is written to xml-dom, xml-file, labels_ vector and segments_ vector
// @return true, if update was completely successful
bool VolumeLabeling::updateSegmentCaption(Label& pLabel, std::string const& newCaption) {
    if (!xmlDoc_.Error()) {
        // find segment-object corresponding to pLabel
        LabelData* segment = 0;
        for (size_t i=0; i<labelPersistentData_.size() && segment == 0; ++i) {
            if (labelPersistentData_[i]->idstr == pLabel.labelData->idstr)
                segment = labelPersistentData_[i];
        }
        if (segment == 0)
            return false;
        if (segment->text.text == newCaption)
            return true;

        // update segment caption
        TiXmlNode* pNode = segment->xmlNode;
        TiXmlNode* pChild;
        if ((pChild = pNode->FirstChild("caption")) != 0) {
            TiXmlNode* text = pChild->FirstChild();
            if ( (text) && (text->Type() == TiXmlNode::TINYXML_TEXT) ) {
                text->SetValue(newCaption);
                xmlDoc_.SaveFile();
                if (!xmlDoc_.Error()) {
                    pLabel.labelData->text.text = newCaption;
                    segment->text.text = newCaption;
                    return true;
                }
            }
        }

    }
    return false;
}



// reads labels' ids and captions from xml-file.
// the xml-dom is kept object-global for runtime-updates of the xml-file
//
// expected xml-file structure:
//  <segmentation>
//      <label>
//          <id>1</id>
//          <caption>caption1</caption>
//      </label>
//      <label>
//          <id>2</id>
//          <caption>caption2</caption>
//      </label>
//      ...
//  </segmentation>
//
bool VolumeLabeling::readSegmentData(std::string filename) {
    std::istringstream is;
    std::ostringstream os;

    labelPersistentData_.clear();

    if (xmlDoc_.LoadFile(filename.c_str())) {

        TiXmlNode* pParent = &xmlDoc_;

        // load segmented labels
        TiXmlElement* pLabel = pParent->FirstChildElement("segmentation")->FirstChildElement("label");
        TiXmlElement* pChild;
        while (pLabel) {
            LabelData* labelData = new LabelData;
          //  segment.text = new LabelText();

            pChild = pLabel->FirstChildElement("id");
            if (!pChild) {
                pLabel = pLabel->NextSiblingElement("label");
                continue;
            }
            is.clear();
            is.str(pChild->GetText());
            is >> labelData->id;
            labelData->idstr = pChild->GetText();

            pChild = pLabel->FirstChildElement("caption");
            if (!pChild) {
                pLabel = pLabel->NextSiblingElement("label");
                continue;
            }
            labelData->text.text = pChild->GetText();

            labelData->xmlNode = pLabel;
            labelData->internPreferred = false;
            labelData->belongsToSegment = true;

            labelPersistentData_.push_back(labelData);

            // register segment as picking object
            //TODO: tc_
            //labelingPort_.getTextureContainer()->getIDManager().addNewPickObj(labelData->idstr);

            pLabel = pLabel->NextSiblingElement("label");
        }

        // load unsegmented labels
        /*int cur_id = 1000;
        pLabel = pParent->FirstChildElement("unsegmentedLabels")->FirstChildElement("label");
        while (pLabel) {
            LabelData* labelData = new LabelData;

            pChild = pLabel->FirstChildElement("caption");
            if (!pChild){
                pLabel = pLabel->NextSiblingElement("label");
                continue;
            }
            labelData->text.text = pChild->GetText();
            labelData->text.textureExtern = -1;
            labelData->text.textureIntern = -1;

            labelData->xmlNode = pLabel;
            labelData->internPreferred = true;
            labelData->id = cur_id;
            os.clear();
            os << cur_id;
            labelData->idstr = os.str();
            cur_id++;
            labelData->belongsToSegment = false;

            labelPersistentData_.push_back(labelData);

            // register segment as picking object
            idManager_.addNewPickObj(labelData->idstr);

            pLabel = pLabel->NextSiblingElement("label");
        } */
    }
    else {
        LWARNING("Failed to load segment description file \"" << filename << "\"");
        return false;
    }

    LINFO("Read segment data from file \"" << filename << "\"");
    return true;
}

void VolumeLabeling::labelFileChanged() {
    if (isInitialized() && !labelFile_.get().empty()) {
        readSegmentData(labelFile_.get());
    }
}

void VolumeLabeling::createFilterKernels() {
    // no filter
    kernels_.push_back(new FilterKernel(1, "no filtering"));

    // average 3x3
    kernels_.push_back(new FilterKernel(3, "average 3x3"));

    // average 5x5
    kernels_.push_back(new FilterKernel(5, "average 5x5"));

    // average 7x7
    kernels_.push_back(new FilterKernel(7, "average 7x7"));

    // average 9x9
    kernels_.push_back(new FilterKernel(9, "average 9x9"));

    // gaussian blur 3x3
    kernels_.push_back(new FilterKernel(3, "gaussian 3x3"));
    kernels_.back()->coefficientsSum = 16;
    kernels_.back()->kernel[0] = 1;
    kernels_.back()->kernel[1] = 2;
    kernels_.back()->kernel[2] = 1;

    kernels_.back()->kernel[3] = 2;
    kernels_.back()->kernel[4] = 4;
    kernels_.back()->kernel[5] = 2;

    kernels_.back()->kernel[6] = 1;
    kernels_.back()->kernel[7] = 2;
    kernels_.back()->kernel[8] = 1;
}

#ifdef labelDEBUG
bool inline inList(std::vector<ivec2> list, ivec2 elem) {
    bool found = false;
    for (size_t i=0; i<list.size() && !found; i++) {
        if (list[i] == elem)
            found = true;
    }
    return found;
}

void VolumeLabeling::renderMaps() {
    if (image_.width == 0 || image_.height == 0)
        return;
/*
    if (segmentationTarget_ == -1) {
        segmentationTarget_ = tc_->allocTarget(getTargetType(Identifier::ttUnknown_), "Labeling::renderMaps (id map)");
        tc_->setPersistent(segmentationTarget_, true);
    }
    if (distanceMapTarget_ == -1) {
        distanceMapTarget_ = tc_->allocTarget(getTargetType(Identifier::ttUnknown_), "Labeling::renderMaps (distance map)");
        tc_->setPersistent(distanceMapTarget_, true);

     //   tc_->changeType(distanceMapTarget_+1, getTargetType(Identifier::Unknown));
     //   tc_->changeType(distanceMapTarget_+2, getTargetType(Identifier::Unknown));
     //   tc_->changeType(distanceMapTarget_+3, getTargetType(Identifier::Unknown));
     //   tc_->changeType(distanceMapTarget_+4, getTargetType(Identifier::Unknown));
    }
*/
    float* idField = new float[image_.width*image_.height];
    float* distanceField = new float[image_.width*image_.height];
    float* hDistanceField = new float[image_.width*image_.height];
    float* vDistanceField = new float[image_.width*image_.height];
    float* adDistanceField = new float[image_.width*image_.height];
    float* ddDistanceField = new float[image_.width*image_.height];

    int maxID = 0;
    int* maxDistance = new int[255];
    for (int i = 0; i < 255; i++)
        maxDistance[i] = 0;

    int maxhDistance = 0;
    int maxvDistance = 0;
    int maxadDistance = 0;
    int maxddDistance = 0;
    for (int i=0; i<image_.width*image_.height; i++) {
        int id = image_.idBuffer.getElem(i);
        if (image_.idBuffer.getElem(i) > maxID) maxID = image_.idBuffer.getElem(i);
        if (image_.distanceField.getElem(i) > maxDistance[id]) maxDistance[id] = image_.distanceField.getElem(i);
        if (image_.horzDistance.getElem(i) > maxhDistance) maxhDistance = image_.horzDistance.getElem(i);
        if (image_.vertDistance.getElem(i) > maxvDistance) maxvDistance = image_.vertDistance.getElem(i);
        if (image_.ascDiagDistance.getElem(i) > maxadDistance) maxadDistance = image_.ascDiagDistance.getElem(i);
        if (image_.descDiagDistance.getElem(i) > maxddDistance) maxddDistance = image_.descDiagDistance.getElem(i);
    }

    for (int x=0; x<image_.width; x++) {
        for (int y=0; y<image_.height; y++) {
            int id = image_.idBuffer.getElem(x,y);
            idField[x + y*image_.width] = image_.idBuffer.getElem(x,y) / (float)maxID;
            distanceField[x + y*image_.width] = image_.distanceField.getElem(x,y) / static_cast<float>(maxDistance[id]);
            hDistanceField[x + y*image_.width] = image_.horzDistance.getElem(x,y) / static_cast<float>(maxhDistance);
            vDistanceField[x + y*image_.width] = image_.vertDistance.getElem(x,y) / static_cast<float>(maxvDistance);
            adDistanceField[x + y*image_.width] = image_.ascDiagDistance.getElem(x,y) / static_cast<float>(maxadDistance);
            ddDistanceField[x + y*image_.width] = image_.descDiagDistance.getElem(x,y) / static_cast<float>(maxddDistance);
        }
    }

    glBindTexture(tc_->getGLTexTarget(segmentationTarget_), tc_->getGLTexID(segmentationTarget_));
    glTexImage2D(tc_->getGLTexTarget(segmentationTarget_), 0, GL_RGB, image_.width, image_.height, 0, GL_RED, GL_FLOAT, idField);

    glBindTexture(tc_->getGLTexTarget(distanceMapTarget_), tc_->getGLTexID(distanceMapTarget_));
    glTexImage2D(tc_->getGLTexTarget(distanceMapTarget_), 0, GL_RGB, image_.width, image_.height, 0, GL_RED, GL_FLOAT, distanceField);

  /*  glBindTexture(tc_->getGLTexTarget(distanceMapTarget_+1), tc_->getGLTexID(distanceMapTarget_+1));
    glTexImage2D(tc_->getGLTexTarget(distanceMapTarget_+1), 0, GL_RGB, image_.width, image_.height, 0, GL_RED, GL_FLOAT, hDistanceField);

    glBindTexture(tc_->getGLTexTarget(distanceMapTarget_+2), tc_->getGLTexID(distanceMapTarget_+2));
    glTexImage2D(tc_->getGLTexTarget(distanceMapTarget_+2), 0, GL_RGB, image_.width, image_.height, 0, GL_RED, GL_FLOAT, vDistanceField);

    glBindTexture(tc_->getGLTexTarget(distanceMapTarget_+3), tc_->getGLTexID(distanceMapTarget_+3));
    glTexImage2D(tc_->getGLTexTarget(distanceMapTarget_+3), 0, GL_RGB, image_.width, image_.height, 0, GL_RED, GL_FLOAT, adDistanceField);

    glBindTexture(tc_->getGLTexTarget(distanceMapTarget_+4), tc_->getGLTexID(distanceMapTarget_+4));
    glTexImage2D(tc_->getGLTexTarget(distanceMapTarget_+4), 0, GL_RGB, image_.width, image_.height, 0, GL_RED, GL_FLOAT, ddDistanceField); */

  /*  tc_->setActiveTarget(distanceMapTarget_);
    glDisable(GL_DEPTH_TEST);
    glPointSize(2.f);
    saveMatrixStacks();
    glBegin(GL_POINTS);
    for (size_t i=0; i<labels_.size(); ++i) {
        vector<ivec2> path = labels_[i].controlPoints;
        glColor3f(0.f, 0.f, 1.f);
        for (size_t j=0; j<path.size(); ++j) {
            glVertex2f(path[j].x / static_cast<float>(image_.width) * 2.f -1.f,
                path[j].y / static_cast<float>(image_.height) * 2.f -1.f);
        }
        glColor3f(0.f, 1.f, 0.f);
        glVertex2f(labels_[i].anchorPoint.x / static_cast<float>(image_.width) * 2.f -1.f,
            labels_[i].anchorPoint.y / static_cast<float>(image_.height) * 2.f -1.f);
    }
    glEnd();
    restoreMatrixStacks();
    glEnable(GL_DEPTH_TEST); */

    delete[] idField;
    delete[] distanceField;
    delete[] hDistanceField;
    delete[] vDistanceField;
    delete[] adDistanceField;
    delete[] ddDistanceField;
    delete[] maxDistance;
}

void VolumeLabeling::writeDistanceMapToFile() {
    std::ofstream myfile;
    std::ostringstream fileName;
    myfile.open ( "distancemap.dat" );
    myfile << "# labeling distance map " << std::endl;
    myfile << "# size: " << image_.tr.x - image_.bl.x << "x" << image_.tr.y - image_.bl.y << std::endl;
    myfile << "# step: " << distanceMapStep_.get() << std::endl;
    myfile << "# filter: " << kernels_[filterKernel_->get()].caption << std::endl;
    myfile << "# filter gap: " << filterDelta_.get() << std::endl;
    for (int x=image_.bl.x; x<image_.tr.x; ++x) {
        for (int y=image_.bl.y; y<image_.tr.y; ++y) {
            int value = image_.distanceField.getElem(x,y);
            myfile << x << " " << y << " " << value << std::endl;
        }
        myfile << std::endl;
    }
    myfile.close();
}
#endif

// propety methods invoked by onChange()
//


void VolumeLabeling::updateFontSettings() {
    if ( lockInternalFontSettings_.get() ){
        *fontIntern_.get() = *fontExtern_.get();
        fontIntern_.updateWidgets();
         haloColorIntern_.set( haloColorExtern_.get() );
    }
}

void VolumeLabeling::setShape3DEvt() {
    for (size_t i = 0; i < labels_.size(); ++i) {
        if (labels_[i].intern)
            labels_[i].intern = findBezierPoints(labels_[i]);
    }
    toWorld();
    pickedLabel_ = 0;
    time.reset();
    time.start();

    invalidate();
}

void VolumeLabeling::unlockInternalFontSettingsEvt() {
    lockInternalFontSettings_.set(false);
}

// ---- events ------ //

void VolumeLabeling::mousePressEvent(MouseEvent* e) {
    e->ignore();
    if (e->button() != MouseEvent::MOUSE_BUTTON_LEFT)
        return;

    // only accept mouse-events, if current labeling valid
    // e.g.: invalid in coarseness-mode
    if (valid_) {
        int x = e->x();
        int y = image_.height - e->y();
        Label* pLabel = getPickedLabel(x,y);

        // if a label is picked, accept mouse event and enter edit mode
        if (pLabel) {

            e->accept();
            editMode_ = true;
            if (pLabel->intern) {
                pLabel->labelPos.x = x;
                pLabel->labelPos.y = y;
            }

            pickedLabel_ = pLabel;
            drag_ = true;
            pickedPointOffset_ = ivec2(x,y) - pickedLabel_->labelPos;
            lastDragPoint_ = ivec2(x,y);

            time.reset();
            time.start();
        // if no label picked, stop edit mode
        }
        else {
            editMode_ = false;
            drag_ = false;
        }
    }
}

void VolumeLabeling::mouseDoubleClickEvent(MouseEvent* e) {
    e->ignore();

    if (labelingWidget_ == 0)
        return;

    if ( pickedLabel_ ) {
        int x = e->x();
        int y = image_.height - e->y();
        Label* pLabel = getPickedLabel(x,y);

        if ( pickedLabel_ && pLabel && pickedLabel_->labelData->id == pLabel->labelData->id ) {
            e->accept();

            // show input box
            string title = "Voreen";
            string label = "Please insert new segment caption \nfor segment";
            label.append(":    ");
            label.append(pLabel->labelData->text.text);
            string text = pLabel->labelData->text.text;
            string newcaption = labelingWidget_->showInputDialog(title, label, text);

            if (newcaption.length() > 0) {
                // if label was successfully updated, recalculate label layout
                // and request new rendering
                if (updateSegmentCaption(*pLabel, newcaption)) {
                    calcLabelPositions();
                    toWorld();
                    invalidate();
                }
                else
                    LERROR("Labeling: failed saving new segment caption to xml-file");
            }
            pickedLabel_ = 0;
        }
    }
}

void VolumeLabeling::mouseReleaseEvent(MouseEvent* e) {
    // only accept release-event if in edit mode
    // (= a label was picked before by a keypress-event)
    if (editMode_ && e->button() == MouseEvent::MOUSE_BUTTON_LEFT)
        e->accept();
    else
        e->ignore();

    // dragging mode ends with releasing mouse key
    if (drag_) {
        drag_ = false;

        // a label is considered "catched", if its connection line does not intersect the convex hull
        // and the current cursor pos is inside the volume
        if (pickedLabel_) {

            //if (pickedLabel_->offLabel){
            //    pickedLabel_->curve2D->fixShift();
            //}

            // check if label is "catched"
            if ( pickedLabel_->offLabel || catchedBySegment(*pickedLabel_, ivec2(e->x(), image_.height - e->y())) )
                pickedLabel_->labelData->internPreferred = true;
            else
                pickedLabel_->labelData->internPreferred = false;

            // synchronize picked label and its segment data.
            // if label is intern, calc its intern label pos
            // + redraw labels
            LabelData* labelData = pickedLabel_->labelData;
            if (labelData) {
                if (labelData->internPreferred && !pickedLabel_->intern) {
                    bool success = false;
                    // try to reuse 2D fitting curve by shifting it by the drag offset
                    if (pickedLabel_->curve2D != NULL)
                        success = true;
                    else
                        success = findLabelPathBest(*pickedLabel_);

                    // recalculate bezier control points
                    if ( success && findBezierPoints(*pickedLabel_) )
                        pickedLabel_->intern = true;
                    else
                        pickedLabel_->intern = false;

                    toWorld(pickedLabel_);
                }
            }
            time.reset();
            time.start();
            invalidate();
        }
        time.start();
    }

#ifdef labelDEBUG
//    std::cout << "Labeling: Mouse Released" << std::endl;
#endif
};

void VolumeLabeling::mouseMoveEvent(MouseEvent* e) {
    e->ignore();

    // if user is dragging a label, accept mouse event
    // and update the dragged label's position
    if (drag_) {

        e->accept();
        time.reset();
        time.start();
        editMode_ = true;

        // when dragging, update label position and redraw labels
        Label* label = pickedLabel_;
        if (label) {
            ivec2 mousePos = ivec2(e->x(), image_.height-e->y());
            if ( pickedLabel_->labelData->internPreferred && ( catchedBySegment(*pickedLabel_, mousePos) || pickedLabel_->offLabel)
                 && pickedLabel_->curve2D != NULL) {
                 // try to reuse 2D fitting curve by shifting it by the drag offset
                 vec2 shift = vec2(static_cast<float>(e->x()) - lastDragPoint_.x, static_cast<float>(image_.height-e->y()) - lastDragPoint_.y);
                 pickedLabel_->curve2D->shift( shift );
                 lastDragPoint_ = ivec2(e->x(), image_.height-e->y());
                 if ( findBezierPoints(*pickedLabel_) ) {
                    pickedLabel_->intern = true;
                 }
                 else {
                    pickedLabel_->intern = false;
                 }
            }
            else {
                pickedLabel_->intern = false;
            }
            if (!pickedLabel_->intern) {
                ivec2 normal;
                label->labelPos = ivec2(e->x(), image_.height - e->y()) - pickedPointOffset_;
                connPointFromLabelPos(*label, label->connectionPoint, normal, true);
            }
            toWorld(pickedLabel_);
            invalidate();
        }
    }

#ifdef labelDEBUG
 //   std::cout << "Labeling: Mouse Moved" << std::endl;
#endif
};

void VolumeLabeling::wheelEvent(tgt::MouseEvent* e) {
    e->ignore();

    if (pickedLabel_ && !pickedLabel_->labelData->belongsToSegment) {
        e->accept();
        float angle = 4.f * ( e->button() == e->MOUSE_WHEEL_DOWN ? 1 : -1);
        rotateLabel(*pickedLabel_, angle);
        if ( pickedLabel_->intern && findBezierPoints(*pickedLabel_) )
            pickedLabel_->intern = true;
        else
            pickedLabel_->intern = false;

        toWorld(pickedLabel_);
        invalidate();
    }
}
} // namespace voreen
