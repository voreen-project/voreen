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

#include "optimizedproxygeometry.h"

#include "voreen/core/datastructures/transfunc/transfuncbase.h"
#include "voreen/core/datastructures/transfunc/1d/transfunc1d.h"
#include "voreen/core/datastructures/transfunc/2d/transfunc2d.h"

#include "voreen/core/datastructures/geometry/trianglemeshgeometry.h"
#include <queue>

namespace {

/**
 * Helper class for iterating over a volume
 */
class VolumeIterator {
private:
    tgt::ivec3 llf_;
    tgt::ivec3 urb_;
    tgt::ivec3 pos_;
public:
    VolumeIterator(tgt::ivec3 llf, tgt::ivec3 urb) : llf_(llf), urb_(urb), pos_(llf) {}
    VolumeIterator(tgt::ivec3 size) : llf_(tgt::ivec3(0)), urb_(size-1), pos_(tgt::ivec3(0)) {}
    void next() {
        pos_.x++;
        if (pos_.x > urb_.x) {
            pos_.x = llf_.x;
            pos_.y++;
            if (pos_.y > urb_.y) {
                pos_.y = llf_.y;
                pos_.z++;
            }
        }
    }
    bool hasnext() {return pos_.x < urb_.x || pos_.y < urb_.y || pos_.z < urb_.z;}
    tgt::ivec3 value() {return pos_;}
    bool outofrange() {return pos_.z > urb_.z;}
    tgt::ivec3 getnext() {next(); return value();}
};

} // namespace anonymous

//-----------------------------------------------------------------------------

namespace voreen {

// OptimizedProxyGeometry
const std::string OptimizedProxyGeometry::loggerCat_("voreen.base.OptimizedProxyGeometry");

OptimizedProxyGeometry::OptimizedProxyGeometry()
    : Processor()
    , inport_(Port::INPORT, "volumehandle.volumehandle", "Volume Input")
    , outport_(Port::OUTPORT, "proxygeometry.geometry", "Proxy Geometry Output")
    , geometry_(0)
    , tmpGeometry_(0)
    , mode_("modeString", "Optimization Mode")
    , tfChannel0_("transferfunction", "Transfer Function 1D (Channel 1) (Link with Raycaster)")
    , resolutionMode_("resolutionMode", "Resolution Mode", Processor::VALID)
    , resolution_("resolution", "Resolution", 32, 1, 64)
    , resolutionVoxels_("resolutionvoxel", "Edge Length (Voxels)", 16, 1, 1024)
    , threshold_("threshold", "Visibility Threshold (*10e-4)", 1, 0, 100)
    , enableClipping_("useClipping", "Enable Clipping", false)
    , clipRegion_("clipRegion", "Clip Region", tgt::IntBounds(tgt::ivec3(0), tgt::ivec3(1)), tgt::ivec3(0), tgt::ivec3(1))
    , resetClipPlanes_("resetClipPlanes", "Reset Planes")
    , waitForOptimization_("waitForOptimization", "Wait for optimization", false, Processor::VALID)
    , geometryInvalid_(true)
    , structureInvalid_(true)
    , volStructureSize_(0,0,0)
    , backgroundThread_(0)
    , tfChannel1_("transferfunction2", "Transfer Function 1D (Channel 2) (Link with Raycaster)")
    , tfChannel2_("transferfunction3", "Transfer Function 1D (Channel 3) (Link with Raycaster)")
    , tfChannel3_("transferfunction4", "Transfer Function 1D (Channel 4) (Link with Raycaster)")
    , currentVolume_(0)
{
    //create mesh list geometry
    geometry_ = new TriangleMeshGeometryColorNormal();
    tmpGeometry_ = new TriangleMeshGeometryColorNormal();

    addPort(inport_);
    addPort(outport_);

    mode_.addOption("boundingbox",              "Bounding Box");
    mode_.addOption("minboundingbox",           "Minimal Visible Bounding Box");
    mode_.addOption("visiblebricks",            "Visible Bricks");
    mode_.addOption("outerfaces",               "Visible Bricks (Outer Faces)");
    mode_.addOption("volumeoctree",             "Volume Octree");
    mode_.addOption("volumeoctreeouterfaces",   "Volume Octree (Outer Faces)");
    mode_.set("boundingbox");
    addProperty(mode_);

    addProperty(tfChannel0_);
    tfChannel0_.onChange(MemberFunctionCallback<OptimizedProxyGeometry>(this, &OptimizedProxyGeometry::onTransFuncChange));

    addProperty(tfChannel1_);
    addProperty(tfChannel2_);
    addProperty(tfChannel3_);
    tfChannel1_.onChange(MemberFunctionCallback<OptimizedProxyGeometry>(this, &OptimizedProxyGeometry::onTransFuncChange));
    tfChannel2_.onChange(MemberFunctionCallback<OptimizedProxyGeometry>(this, &OptimizedProxyGeometry::onTransFuncChange));
    tfChannel3_.onChange(MemberFunctionCallback<OptimizedProxyGeometry>(this, &OptimizedProxyGeometry::onTransFuncChange));

    resolutionMode_.addOption("subdivide", "Subdivide Shortest Side");
    resolutionMode_.addOption("voxel", "Subdivide in Voxels");
    addProperty(resolutionMode_);
    resolutionMode_.onChange(MemberFunctionCallback<OptimizedProxyGeometry>(this, &OptimizedProxyGeometry::onResolutionModeChange));

    addProperty(resolution_);

    addProperty(resolutionVoxels_);

    addProperty(threshold_);
    threshold_.onChange(MemberFunctionCallback<OptimizedProxyGeometry>(this, &OptimizedProxyGeometry::onThresholdChange));

    clipRegion_.onChange(MemberFunctionCallback<OptimizedProxyGeometry>(this, &OptimizedProxyGeometry::onClipRegionChanged));
    enableClipping_.onChange(MemberFunctionCallback<OptimizedProxyGeometry>(this, &OptimizedProxyGeometry::adjustClipPropertiesVisibility));
    resetClipPlanes_.onChange(MemberFunctionCallback<OptimizedProxyGeometry>(this, &OptimizedProxyGeometry::resetClipPlanes));

    addProperty(enableClipping_);
    addProperty(clipRegion_);
    addProperty(resetClipPlanes_);
    addProperty(waitForOptimization_);

    clipRegion_.setGroupID("clipping");
    resetClipPlanes_.setGroupID("clipping");
    setPropertyGroupGuiName("clipping", "Clipping Planes");
    adjustClipPropertiesVisibility();

    oldVolumeDimensions_ = tgt::ivec3(0,0,0);

    updatePropertyVisibility();

    mode_.onChange(MemberFunctionCallback<OptimizedProxyGeometry>(this, &OptimizedProxyGeometry::onModeChange));
    tfChannel0_.onChange(MemberFunctionCallback<OptimizedProxyGeometry>(this, &OptimizedProxyGeometry::onTransFuncChange));
    resolution_.onChange(MemberFunctionCallback<OptimizedProxyGeometry>(this, &OptimizedProxyGeometry::onResolutionChange));
    resolutionVoxels_.onChange(MemberFunctionCallback<OptimizedProxyGeometry>(this, &OptimizedProxyGeometry::onResolutionVoxelChange));
}

OptimizedProxyGeometry::~OptimizedProxyGeometry() {

    if (backgroundThread_) {
        backgroundThread_->interrupt();
        delete backgroundThread_;
    }

    geometry_->clear();
    tmpGeometry_->clear();
    delete geometry_;
    delete tmpGeometry_;

    for (std::vector<TransFunc1DKeys*>::iterator i = tfCopies_.begin(); i != tfCopies_.end(); ++i) {
        delete *i;
    }
    tfCopies_.clear();
}

Processor* OptimizedProxyGeometry::create() const {
    return new OptimizedProxyGeometry();
}

void OptimizedProxyGeometry::adjustPropertiesToInput() {

    // interrupt background thread
    if (backgroundThread_) {
        backgroundThread_->interrupt();
        delete backgroundThread_;
        backgroundThread_ = 0;
    }

    // if we previously had another input volume which has not been deleted, we have to remove ourselves as observer
    if (currentVolume_)
        currentVolume_->Observable<VolumeObserver>::removeObserver(this);

    // if there is a new volume, we need to register as observer
    currentVolume_ = inport_.getData();
    if (currentVolume_)
        currentVolume_->Observable<VolumeObserver>::addObserver(this);

    // if we yhave new data, we need to adjust our properties
    if (inport_.getData())
        onVolumeChange();
}

void OptimizedProxyGeometry::volumeDelete(const VolumeBase* source) {
    tgtAssert(currentVolume_ == source, "Notified by other volume than current");

    // we do not want to observe the old volume anymore (although it is deleted anyway)
    source->Observable<VolumeObserver>::removeObserver(this);
    currentVolume_ = 0;

    if (backgroundThread_) {
        backgroundThread_->interrupt();
        delete backgroundThread_;
        backgroundThread_ = 0;
    }

    tfChannel0_.setVolume(0);
    tfChannel1_.setVolume(0);
    tfChannel2_.setVolume(0);
    tfChannel3_.setVolume(0);

    invalidate();
}

void OptimizedProxyGeometry::volumeRepresentationDelete(const VolumeBase* source, const VolumeRepresentation* rep) {
    tgtAssert(currentVolume_ == source, "Notified by other volume than current");
    if (dynamic_cast<const VolumeRAM*>(rep)) {

        if (backgroundThread_) {
            backgroundThread_->interrupt();
            delete backgroundThread_;
            backgroundThread_ = 0;
        }

        tfChannel0_.setVolume(0);
        tfChannel1_.setVolume(0);
        tfChannel2_.setVolume(0);
        tfChannel3_.setVolume(0);

        invalidate();
    }
}

void OptimizedProxyGeometry::volumeDataDelete(const VolumeBase* source) {
    tgtAssert(currentVolume_ == source, "Notified by other volume than current");

    if (backgroundThread_) {
        backgroundThread_->interrupt();
        delete backgroundThread_;
        backgroundThread_ = 0;
    }

    tfChannel0_.setVolume(0);
    tfChannel1_.setVolume(0);
    tfChannel2_.setVolume(0);
    tfChannel3_.setVolume(0);

    invalidate();
}


void OptimizedProxyGeometry::volumeChange(const VolumeBase* source) {
    tgtAssert(currentVolume_ == source, "Notified by other volume than current");

    if (backgroundThread_) {
        backgroundThread_->interrupt();
        delete backgroundThread_;
        backgroundThread_ = 0;
    }

    tfChannel0_.setVolume(0);
    tfChannel1_.setVolume(0);
    tfChannel2_.setVolume(0);
    tfChannel3_.setVolume(0);

    // adjust properties etc. if necessary
    onVolumeChange();

    invalidate();
}

void OptimizedProxyGeometry::process() {
    tgtAssert(inport_.getData(), "no input volume");

    //check if VolumeRAM is availabe (if necessary)
    if (!mode_.isSelected("boundingbox") &&
            (!mode_.hasKey("volumeoctree") || (!mode_.isSelected("volumeoctree") && !mode_.isSelected("volumeoctreeouterfaces")))) {

        inport_.getData()->getRepresentation<VolumeRAM>();
        bool hasVolumeRam = inport_.getData()->hasRepresentation<VolumeRAM>();

        if (!hasVolumeRam) {
            LWARNING("VolumeRAM not available. Falling back to bounding box.");
            mode_.set("boundingbox");
            invalidate();
            return;
        }
    }

    if (mode_.isSelected("boundingbox")) {
         processCube();
    }
    else {
        tgtAssert(tfChannel0_.get(), "no transfunc");

        //if background thread finished computation: do nothing (background thread invalidated processor, mesh geometry is valid)
        //else: compute new geometry in background thread and set temporary geometry to outport
        if (!backgroundThread_ || !backgroundThread_->isFinished()) {


        if (backgroundThread_) {
            backgroundThread_->interrupt();
            delete backgroundThread_;
            backgroundThread_ = 0;
        }

            //get number of channels in volume
            size_t numChannels = inport_.getData()->getNumChannels();

            //only 4 channels supported
            if (numChannels > 4) {
                LERROR("Currently only 4 channels supported. Falling back to Bounding Box mode.");
                mode_.set("boundingbox");
                invalidate();
                return;
            }

            //copy transfer functions
            if (tfCopies_.size() != numChannels) {
                for (std::vector<TransFunc1DKeys*>::iterator i = tfCopies_.begin(); i != tfCopies_.end(); ++i) {
                    delete *i;
                }
                tfCopies_.clear();

                tfCopies_.push_back(static_cast<TransFunc1DKeys*>(tfChannel0_.get()->clone()));
                if (numChannels > 1)
                    tfCopies_.push_back(static_cast<TransFunc1DKeys*>(tfChannel1_.get()->clone()));
                if (numChannels > 2)
                    tfCopies_.push_back(static_cast<TransFunc1DKeys*>(tfChannel2_.get()->clone()));
                if (numChannels > 3)
                    tfCopies_.push_back(static_cast<TransFunc1DKeys*>(tfChannel3_.get()->clone()));
            }

            // determine TF types and set up list of transfuncs
            std::vector<TransFunc1D*> tfVector;

            for (int i = 0; static_cast<size_t>(i) < tfCopies_.size(); ++i) {
                tfVector.push_back(tfCopies_[i]);
            }

            //get step size (for modes using bricks) and clipping parameters
            int stepSize = resolutionVoxels_.get();
            tgt::vec3 clipLlf = clipRegion_.get().getLLF();
            tgt::vec3 clipUrb = clipRegion_.get().getURB();


            //create background thread according to selected mode
            if (mode_.isSelected("minboundingbox"))
                backgroundThread_ = new MinCubeBackgroundThread(this, inport_.getData(), tfVector, static_cast<float>(threshold_.get()),
                    geometry_, &volumeStructure_, volStructureSize_, stepSize, false, enableClipping_.get(), clipLlf, clipUrb);
            else if (mode_.isSelected("visiblebricks"))
                backgroundThread_ = new VisibleBricksBackgroundThread(this, inport_.getData(), tfVector, static_cast<float>(threshold_.get()),
                    geometry_, &volumeStructure_, volStructureSize_, stepSize, false, enableClipping_.get(), clipLlf, clipUrb);
            else if (mode_.isSelected("outerfaces"))
                backgroundThread_ = new OuterFacesBackgroundThread(this, inport_.getData(), tfVector, static_cast<float>(threshold_.get()),
                    geometry_, &volumeStructure_, volStructureSize_, stepSize, false, enableClipping_.get(), clipLlf, clipUrb);
            else if (mode_.isSelected("volumeoctree")) {
                //select if the volume has an octree representation and start a VolumeOctreeBackgroundThread or use visible bricks as fallback mode
                if (inport_.getData()->hasRepresentation<VolumeOctreeBase>()) {
                    LDEBUG("VolumeOctree representation available");

                //create worker thread that computes a proxy geometry based on the volume octree
                backgroundThread_ = new VolumeOctreeBackgroundThread(this, inport_.getData(), tfVector,
                        static_cast<float>(threshold_.get()), geometry_, stepSize, enableClipping_.get(), clipLlf, clipUrb);
                }
                else {
                    LERROR("VolumeOctree representation not available (try to use OctreeCreator). Falling back to Visible Bricks mode.");
                    mode_.set("visiblebricks");
                    invalidate();
                    return;
                }
            }
            else if (mode_.isSelected("volumeoctreeouterfaces")) {
                //select if the volume has an octree representation and start a VolumeOctreeOuterFacesBackgroundThread or use outer faces as fallback mode
                if (inport_.getData()->hasRepresentation<VolumeOctreeBase>()) {
                    LDEBUG("VolumeOctree representation available");

                //create worker thread that computes a proxy geometry based on the volume octree
                backgroundThread_ = new VolumeOctreeOuterFacesBackgroundThread(this, inport_.getData(), tfVector,
                        static_cast<float>(threshold_.get()), geometry_, stepSize, enableClipping_.get(), clipLlf, clipUrb);
                }
                else {
                    LERROR("VolumeOctree representation not available (try to use OctreeCreator). Falling back to Outer Faces mode.");
                    mode_.set("outerfaces");
                    invalidate();
                    return;
                }
            }

            //start background computation
            backgroundThread_->run();

            if (waitForOptimization_.get()) {
                // wait for background thread to finish computation
                backgroundThread_->join();
            }
            else {
                //while background computation is not finished: use temporary bounding box geometry
                processTmpCube();
                outport_.setData(tmpGeometry_, false);

                return;
            }
        }
    }

    outport_.setData(geometry_, false);
}

void OptimizedProxyGeometry::processCube() {
    // clear the current geometry
    geometry_->clear();

    // get voxel space dimensions
    const VolumeBase* inputVolume = inport_.getData();
    tgt::vec3 voxelLLF = tgt::vec3(0.f);
    tgt::vec3 voxelURB = tgt::vec3(inputVolume->getDimensions() - tgt::svec3::one);

    // vertex coords of bounding box in voxel space without clipping (add a border of 0.5 voxels since a voxel is sampled at the center)
    tgt::vec3 coordLLF = voxelLLF - tgt::vec3(0.5f);
    tgt::vec3 coordURB = voxelURB + tgt::vec3(0.5f);

    // take into account the clipping
    if (enableClipping_.get()) {
        coordLLF = tgt::max(coordLLF, tgt::vec3(clipRegion_.get().getLLF()) - tgt::vec3(0.5f));
        coordURB = tgt::min(coordURB, tgt::vec3(clipRegion_.get().getURB()) + tgt::vec3(0.5f));
    }

    // compute the texture coordinates for the voxels
    tgt::vec3 texLLF = inputVolume->getVoxelToTextureMatrix() * coordLLF;
    tgt::vec3 texURB = inputVolume->getVoxelToTextureMatrix() * coordURB;

    // create output mesh
    geometry_->addCube(VertexNormal(texLLF, texLLF), VertexNormal(texURB, texURB));

    // coordinates are set in texture space -> use texture to world matrix as transformation
    geometry_->setTransformationMatrix(inputVolume->getTextureToWorldMatrix());
}

void OptimizedProxyGeometry::processTmpCube() {
    // clear the current temporary geometry
    tmpGeometry_->clear();

    // get voxel space dimensions
    const VolumeBase* inputVolume = inport_.getData();
    tgt::vec3 voxelLLF = tgt::vec3(0.f);
    tgt::vec3 voxelURB = tgt::vec3(inputVolume->getDimensions() - tgt::svec3::one);

    // vertex coords of bounding box in voxel space without clipping (add a border of 0.5 voxels since a voxel is sampled at the center)
    tgt::vec3 coordLLF = voxelLLF - tgt::vec3(0.5f);
    tgt::vec3 coordURB = voxelURB + tgt::vec3(0.5f);

    // take into account the clipping
    if (enableClipping_.get()) {
        coordLLF = tgt::max(coordLLF, tgt::vec3(clipRegion_.get().getLLF()) - tgt::vec3(0.5f));
        coordURB = tgt::min(coordURB, tgt::vec3(clipRegion_.get().getURB()) + tgt::vec3(0.5f));
    }

    // compute the texture coordinates for the voxels
    tgt::vec3 texLLF = inputVolume->getVoxelToTextureMatrix() * coordLLF;
    tgt::vec3 texURB = inputVolume->getVoxelToTextureMatrix() * coordURB;

    // create output mesh
    tmpGeometry_->addCube(VertexNormal(texLLF, texLLF), VertexNormal(texURB, texURB));

    // coordinates are set in texture space -> use texture to world matrix as transformation
    tmpGeometry_->setTransformationMatrix(inputVolume->getTextureToWorldMatrix());
}

void OptimizedProxyGeometry::addCubeMesh(TriangleMeshGeometryColorNormal* mesh, tgt::Bounds bounds, tgt::ivec3 dim) {
    // bounds for the block already contain the 0.5 voxel boundary so that we can directly use them
    tgt::vec3 coordllf = bounds.getLLF();
    tgt::vec3 coordurb = bounds.getURB();

    // these are the same computations as in voxelToTexture matrix
    tgt::vec3 texllf = (coordllf + tgt::vec3(0.5f)) / tgt::vec3(dim);
    tgt::vec3 texurb = (coordurb + tgt::vec3(0.5f)) / tgt::vec3(dim);

    // positions are set in texture coordinates
    mesh->addCube(VertexNormal(texllf, texllf), VertexNormal(texurb, texurb));
}

void OptimizedProxyGeometry::addCubeMeshClip(TriangleMeshGeometryColorNormal* mesh, tgt::Bounds bounds, tgt::ivec3 dim, tgt::Bounds clipBounds) {
    // bounds for the block already contain the 0.5 voxel boundary so that we can directly use them
    tgt::vec3 coordllf = bounds.getLLF();
    tgt::vec3 coordurb = bounds.getURB();

    // clip bounds are in voxels -> 0.5 offset necessary!
    tgt::vec3 clipLlf = clipBounds.getLLF() - tgt::vec3(0.5f);
    tgt::vec3 clipUrb = clipBounds.getURB() + tgt::vec3(0.5f);

    //clip the coordinates
    coordllf = tgt::max(coordllf, clipLlf);
    coordurb = tgt::min(coordurb, clipUrb);

    // these are the same computations as in voxelToTexture matrix
    tgt::vec3 texllf = (coordllf + tgt::vec3(0.5f)) / tgt::vec3(dim);
    tgt::vec3 texurb = (coordurb + tgt::vec3(0.5f)) / tgt::vec3(dim);

    // positions are set in texture coordinates
    mesh->addCube(VertexNormal(texllf, texllf), VertexNormal(texurb, texurb));
}


void OptimizedProxyGeometry::onClipRegionChanged(){
    geometryInvalid_ = true;

    if (backgroundThread_) {
        backgroundThread_->interrupt();
        delete backgroundThread_;
        backgroundThread_ = 0;
    }
}

void OptimizedProxyGeometry::resetClipPlanes() {
    tgt::ivec3 min = clipRegion_.getMinValue();
    tgt::ivec3 max = clipRegion_.getMaxValue();

    clipRegion_.set(tgt::IntBounds(min, max));
}

void OptimizedProxyGeometry::adjustResolutionPropertyRanges() {

    if (inport_.getData()) {
        //store old settings
        int oldVoxels = resolutionVoxels_.get();
        int oldDiv = resolution_.get();

        //get volume dimensions
        tgt::ivec3 dim = inport_.getData()->getDimensions();
        //determine shortest side
        int minDim = std::min(std::min(dim.x, dim.y), dim.z);

        //set max for voxel resolution
        resolutionVoxels_.setMaxValue(minDim);

        //set max for resolution
        resolution_.setMaxValue(minDim);

        /*int div = static_cast<int>(std::ceil(static_cast<float>(minDim) / static_cast<float>(resolutionVoxels_.get())));
        resolution_.set(div);*/
    }
}

void OptimizedProxyGeometry::adjustClipPropertiesRanges() {
    tgtAssert(inport_.getData() && inport_.getData(), "No input volume");

    tgt::ivec3 oldRange = clipRegion_.getMaxValue();
    tgt::ivec3 newRange = tgt::ivec3(inport_.getData()->getDimensions())-tgt::ivec3::one;

    if (oldRange != newRange){
        clipRegion_.setMaxValue(newRange);
        clipRegion_.set(tgt::IntBounds(tgt::ivec3::zero, newRange));
    }
}

void OptimizedProxyGeometry::adjustClipPropertiesVisibility() {
    bool clipEnabled = enableClipping_.get();
    setPropertyGroupVisible("clipping", clipEnabled);
    geometryInvalid_ = true;

    if (backgroundThread_) {
        backgroundThread_->interrupt();
        delete backgroundThread_;
        backgroundThread_ = 0;
    }
}

void OptimizedProxyGeometry::updatePropertyVisibility() {
    if (mode_.get() == "boundingbox") {
        tfChannel0_.setVisibleFlag(false);
        tfChannel1_.setVisibleFlag(false);
        tfChannel2_.setVisibleFlag(false);
        tfChannel3_.setVisibleFlag(false);
        resolutionMode_.setVisibleFlag(false);
        resolution_.setVisibleFlag(false);
        resolutionVoxels_.setVisibleFlag(false);
        threshold_.setVisibleFlag(false);
        enableClipping_.setVisibleFlag(true);
        adjustClipPropertiesVisibility();
    }
    else {

        //check number of channels and set transfer functions visible
        tfChannel0_.setVisibleFlag(true);

        if (inport_.getData()) {
            size_t numChannels = inport_.getData()->getNumChannels();

            if (numChannels > 1)
                tfChannel1_.setVisibleFlag(true);
            else
                tfChannel1_.setVisibleFlag(false);

            if (numChannels > 2)
                tfChannel2_.setVisibleFlag(true);
            else
                tfChannel2_.setVisibleFlag(false);

            if (numChannels > 3)
                tfChannel3_.setVisibleFlag(true);
            else
                tfChannel3_.setVisibleFlag(false);

            if (numChannels > 4)
                LWARNING("Number of channels in volume is greater than 4, only 4 channels currently supported!");
        }
        else {
            tfChannel1_.setVisibleFlag(false);
            tfChannel2_.setVisibleFlag(false);
            tfChannel3_.setVisibleFlag(false);
        }

        if (!mode_.hasKey("volumeoctree") || (!mode_.isSelected("volumeoctree") && !mode_.isSelected("volumeoctreeouterfaces"))) {
            resolutionMode_.setVisibleFlag(true);
            if (resolutionMode_.get() == "voxel") {
                resolution_.setVisibleFlag(false);
                resolutionVoxels_.setVisibleFlag(true);
            }
            else {
                resolution_.setVisibleFlag(true);
                resolutionVoxels_.setVisibleFlag(false);
            }
            resolution_.setVisibleFlag(true);
        }
        else {
            //volume octree mode does not support resolution
            resolutionMode_.setVisibleFlag(false);
            resolution_.setVisibleFlag(false);
            resolutionVoxels_.setVisibleFlag(false);
        }

        threshold_.setVisibleFlag(true);

        enableClipping_.setVisibleFlag(true);
        adjustClipPropertiesVisibility();
    }
}

void OptimizedProxyGeometry::adjustClippingToVolumeROI() {
    // adjust clipping sliders to volume ROI, if set
    if (inport_.getData()->hasMetaData("RoiPixelOffset") && inport_.getData()->hasMetaData("RoiPixelLength")) {
        const tgt::ivec3 volumeDim = static_cast<tgt::ivec3>(inport_.getData()->getDimensions());
        tgt::ivec3 roiLlf = inport_.getData()->getMetaDataValue<IVec3MetaData>("RoiPixelOffset", tgt::ivec3(0));
        tgt::ivec3 roiLength = inport_.getData()->getMetaDataValue<IVec3MetaData>("RoiPixelLength", volumeDim);
        if (tgt::hor(tgt::lessThan(roiLlf, tgt::ivec3::zero)) || tgt::hor(tgt::greaterThanEqual(roiLlf, volumeDim))) {
            LWARNING("Invalid ROI offset: " << roiLlf);
        }
        else if (tgt::hor(tgt::lessThanEqual(roiLength, tgt::ivec3::zero))) {
            LWARNING("Invalid ROI length: " << roiLength);
        }
        else {
            tgt::ivec3 roiUrb = tgt::min(roiLlf + roiLength - 1, volumeDim - tgt::ivec3::one);

            LINFO("Applying volume ROI: llf=" << roiLlf << ", urb=" << roiUrb);

            clipRegion_.set(tgt::IntBounds(roiLlf, roiUrb));
        }
    }
}

void OptimizedProxyGeometry::onThresholdChange() {
    geometryInvalid_ = true;

    if (backgroundThread_) {
        backgroundThread_->interrupt();
        delete backgroundThread_;
        backgroundThread_ = 0;
    }
}

void OptimizedProxyGeometry::onModeChange() {
    //structureInvalid_ = true;
    geometryInvalid_ = true;

    if (backgroundThread_) {
        backgroundThread_->interrupt();
        delete backgroundThread_;
        backgroundThread_ = 0;
    }

    updatePropertyVisibility();
}

void OptimizedProxyGeometry::onTransFuncChange() {

    geometryInvalid_ = true;

    if (backgroundThread_) {
        backgroundThread_->interrupt();
        delete backgroundThread_;
        backgroundThread_ = 0;
    }

    //transfunc copies are invalid
    for (std::vector<TransFunc1DKeys*>::iterator i = tfCopies_.begin(); i != tfCopies_.end(); ++i) {
        delete *i;
    }
    tfCopies_.clear();
}

void OptimizedProxyGeometry::onResolutionModeChange() {
    if (resolutionMode_.get() == "voxel") {
        resolution_.setVisibleFlag(false);
        resolutionVoxels_.setVisibleFlag(true);
    }
    else {
        //have to set the resolution according to the resolution voxel size
        if (inport_.getData()) {
            //get input data
            const VolumeBase* inputVolume = inport_.getData();
            tgt::ivec3 dim = inputVolume->getDimensions();
            //determine shortest side and for this side the number of subdivisions
            int minDim = std::min(std::min(dim.x, dim.y), dim.z);
            int div = static_cast<int>(std::ceil(static_cast<float>(minDim) / static_cast<float>(resolutionVoxels_.get())));
            resolution_.set(div);
        }

        resolution_.setVisibleFlag(true);
        resolutionVoxels_.setVisibleFlag(false);
    }
}

void OptimizedProxyGeometry::onVolumeChange() {

    //invalidate data structures
    structureInvalid_ = true;

    // adapt clipping plane properties on volume change
    adjustClipPropertiesRanges();

    // extract ROI from volume and adjust clipping sliders accordingly
    adjustClippingToVolumeROI();

    //clear region structure
    volumeStructure_.clear();

    //adjust the possible values for the resolution
    adjustResolutionPropertyRanges();

    //set the resolution in voxels, if the other mode is presently used
    if (resolutionMode_.get() == "subdivide") {
        //set voxel value according to previous resolution
        if (inport_.getData()) {
            //get input data
            const VolumeBase* inputVolume = inport_.getData();
            tgt::ivec3 dim = inputVolume->getDimensions();
            //determine shortest side and for this side the step size
            int minDim = std::min(std::min(dim.x, dim.y), dim.z);
            int stepSize = std::max(1, tgt::iround(static_cast<float>(minDim) / static_cast<float>(resolution_.get())));
            resolutionVoxels_.set(stepSize);
        }
    }

    size_t numChannels = inport_.getData()->getNumChannels();

    if (!mode_.isSelected("boundingbox")) {
        //check number of channels and set transfer functions visible

        if (numChannels > 1)
            tfChannel1_.setVisibleFlag(true);
        else
            tfChannel1_.setVisibleFlag(false);

        if (numChannels > 2)
            tfChannel2_.setVisibleFlag(true);
        else
            tfChannel2_.setVisibleFlag(false);

        if (numChannels > 3)
            tfChannel3_.setVisibleFlag(true);
        else
            tfChannel3_.setVisibleFlag(false);

        if (numChannels > 4)
            LWARNING("Number of channels in volume is greater than 4, only 4 channels currently supported!");
    }

    tfChannel0_.setVolume(inport_.getData());

    if (numChannels > 1)
        tfChannel1_.setVolume(inport_.getData(), 1);

    if (numChannels > 2)
        tfChannel2_.setVolume(inport_.getData(), 2);

    if (numChannels > 3)
        tfChannel3_.setVolume(inport_.getData(), 3);
}

void OptimizedProxyGeometry::onResolutionChange() {
    //do not invalidate the structure if the change happened because the resolution mode has been changed
    if (!resolution_.isVisibleFlagSet())
        return;

    //set voxel value according to previous resolution
    if (inport_.getData()) {
        //get input data
        const VolumeBase* inputVolume = inport_.getData();
        tgt::ivec3 dim = inputVolume->getDimensions();
        //determine shortest side and for this side the step size
        int minDim = std::min(std::min(dim.x, dim.y), dim.z);
        int stepSize = std::max(1, tgt::iround(static_cast<float>(minDim) / static_cast<float>(resolution_.get())));
        resolutionVoxels_.set(stepSize);
    }

    if (backgroundThread_) {
        backgroundThread_->interrupt();
        delete backgroundThread_;
        backgroundThread_ = 0;
    }
}

void OptimizedProxyGeometry::onResolutionVoxelChange() {
    //invalidate the structure
    structureInvalid_ = true;

    if (backgroundThread_) {
        backgroundThread_->interrupt();
        delete backgroundThread_;
        backgroundThread_ = 0;
    }
}

bool OptimizedProxyGeometry::structureInvalid() const {
    return structureInvalid_;
}

bool OptimizedProxyGeometry::geometryInvalid() const {
    return geometryInvalid_;
}

void OptimizedProxyGeometry::setStructureInvalid(bool value) {
    structureInvalid_ = value;
}

void OptimizedProxyGeometry::setGeometryInvalid(bool value) {
    geometryInvalid_ = value;
}

void OptimizedProxyGeometry::setVolStructureSize(tgt::ivec3 volStructureSize) {
    volStructureSize_ = volStructureSize;
}


//-------------------------------------------------------------------------------------------------
// background threads

OptimizedProxyGeometryBackgroundThread::OptimizedProxyGeometryBackgroundThread(OptimizedProxyGeometry* processor, const VolumeBase* volume,
          std::vector<TransFunc1D*> tfVector, float threshold, TriangleMeshGeometryColorNormal* geometry, int stepSize, bool debugOutput,
          bool clippingEnabled, tgt::vec3 clipLlf, tgt::vec3 clipUrb)
        : ProcessorBackgroundThread<OptimizedProxyGeometry>(processor)
        , volume_(volume)
        , tfCopyVector_(tfVector)
        , threshold_(threshold)
        , geometry_(geometry)
        , stepSize_(stepSize)
        , debugOutput_(debugOutput)
        , clippingEnabled_(clippingEnabled)
        , clipLlf_(clipLlf)
        , clipUrb_(clipUrb)
{}

void OptimizedProxyGeometryBackgroundThread::handleInterruption() {
    //nothing to handle
}

bool OptimizedProxyGeometryBackgroundThread::isRegionEmptyPi(float min, float max, const PreIntegrationTable* piTable) const {
    return (piTable->classify(min, max).a <= 0.0001 * static_cast<float>(threshold_));
}


//-------------------------------------------------------------------------------------------------
// StructureProxyGeometryBackgroundThread

StructureProxyGeometryBackgroundThread::StructureProxyGeometryBackgroundThread(OptimizedProxyGeometry* processor, const VolumeBase* volume,
          std::vector<TransFunc1D*> tfVector, float threshold, TriangleMeshGeometryColorNormal* geometry, std::vector<ProxyGeometryVolumeRegion>* volumeStructure,
          tgt::ivec3 volStructureSize, int stepSize, bool debugOutput, bool clippingEnabled, tgt::vec3 clipLlf, tgt::vec3 clipUrb)
        : OptimizedProxyGeometryBackgroundThread(processor, volume, tfVector, threshold, geometry, stepSize, debugOutput, clippingEnabled, clipLlf, clipUrb)
        , volumeStructure_(volumeStructure)
        , volStructureSize_(volStructureSize)
{}

void StructureProxyGeometryBackgroundThread::computeRegionStructure() {

    interruptionPoint();

    volumeStructure_->clear();

    if (debugOutput_) {
        stopWatch_.reset();
        stopWatch_.start();
    }

    tgt::ivec3 dim = tgt::ivec3(volume_->getDimensions());

    interruptionPoint();

    const VolumeRAM* vol = volume_->getRepresentation<VolumeRAM>();
    RealWorldMapping rwm = volume_->getRealWorldMapping();
    int numChannels = static_cast<int>(volume_->getNumChannels());

    interruptionPoint();

    const tgt::ivec3 step = tgt::ivec3(stepSize_);

    //determine size for this resolution
    const tgt::ivec3 size(
            static_cast<int>(std::ceil(static_cast<float>(dim.x) / static_cast<float>(stepSize_))),
            static_cast<int>(std::ceil(static_cast<float>(dim.y) / static_cast<float>(stepSize_))),
            static_cast<int>(std::ceil(static_cast<float>(dim.z) / static_cast<float>(stepSize_))));

    volStructureSize_ = size;
    processor_->setVolStructureSize(volStructureSize_);

    interruptionPoint();

    tgt::ivec3 pos, llf, urb;

    for (VolumeIterator it(size); !it.outofrange(); it.next()) {

        interruptionPoint();

        pos = it.value();
        llf = pos * step;
        urb = (pos+1) * step - tgt::ivec3::one;
        urb = tgt::min(urb, dim - tgt::ivec3::one);

        //since a voxel is the center, add/subtract 0.5 to/from the coordinates to get the bounding box
        tgt::Bounds regionBounds(tgt::vec3(llf) - tgt::vec3(0.5f), tgt::vec3(urb) + tgt::vec3(0.5f));

        interruptionPoint();

        // find min and max intensities for all channels (take into account one additional voxel as overlap)
        float minIntensity[4] = {std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max(), std::numeric_limits<float>::max()};
        float maxIntensity[4] = {std::numeric_limits<float>::min(), std::numeric_limits<float>::min(), std::numeric_limits<float>::min(), std::numeric_limits<float>::min()};
        llf = tgt::max(llf - tgt::ivec3(1),tgt::ivec3(0));
        llf = tgt::min(llf,dim - tgt::ivec3::one);
        urb = tgt::max(urb + tgt::ivec3(1),tgt::ivec3(0));
        urb = tgt::min(urb,dim - tgt::ivec3::one);

        //don't use macro because of interruption points within the loops...
        //VRN_FOR_EACH_VOXEL(pos, llf, urb) {
        for (pos = llf; pos.z <= urb.z; ++pos.z) {
            //interruption point after each slice
            interruptionPoint();

            for (pos.y = llf.y; pos.y <= urb.y; ++pos.y) {
                for (pos.x = llf.x; pos.x <= urb.x; ++pos.x) {

                    for (int i = 0; i < numChannels; ++i) {
                        float current = vol->getVoxelNormalized(pos, i);
                        //apply realworld mapping
                        current = rwm.normalizedToRealWorld(current);
                        minIntensity[i] = std::min(minIntensity[i],current);
                        maxIntensity[i] = std::max(maxIntensity[i],current);
                    }
                }
            }
        }

        //create list of min-max intensities
        std::vector<tgt::vec2> minMaxIntensities(numChannels);
        for (int i = 0; i < numChannels; ++i) {
            minMaxIntensities.at(i) = tgt::vec2(minIntensity[i], maxIntensity[i]);
        }

        //add region
        volumeStructure_->push_back(ProxyGeometryVolumeRegion(regionBounds,minMaxIntensities));

    }

    if (debugOutput_) {
        stopWatch_.stop();
        std::cout << "Computing region structure took " << stopWatch_.getRuntime() << " milliseconds" << std::endl;
    }
}

ProxyGeometryVolumeRegion& StructureProxyGeometryBackgroundThread::getVolumeRegion(tgt::ivec3 pos) {
    return volumeStructure_->at(pos.z * (volStructureSize_.x * volStructureSize_.y) + pos.y * volStructureSize_.x + pos.x);
}


//-------------------------------------------------------------------------------------------------
// MinCubeBackgroundThread

MinCubeBackgroundThread::MinCubeBackgroundThread(OptimizedProxyGeometry* processor, const VolumeBase* volume, std::vector<TransFunc1D*> tfVector,
                                                 float threshold, TriangleMeshGeometryColorNormal* geometry,
                                                 std::vector<ProxyGeometryVolumeRegion>* volumeStructure, tgt::ivec3 volStructureSize,
                                                 int stepSize, bool debugOutput, bool clippingEnabled, tgt::vec3 clipLlf, tgt::vec3 clipUrb)
        : StructureProxyGeometryBackgroundThread(processor, volume, tfVector, threshold, geometry, volumeStructure, volStructureSize, stepSize,
                                                 debugOutput, clippingEnabled, clipLlf, clipUrb)
{}

MinCubeBackgroundThread::~MinCubeBackgroundThread() {
    //wait for internal thread to finish
    join();
}

void MinCubeBackgroundThread::threadMain() {
    computeMinCube();
}

void MinCubeBackgroundThread::computeMinCube() {

    //get volume dimensions
    tgt::ivec3 dim = volume_->getDimensions();

    interruptionPoint();

    if (processor_->structureInvalid()) {
        //invalidate geometry
        processor_->setGeometryInvalid();

        //clear structure
        volumeStructure_->clear();

        interruptionPoint();

        //compute new structure
        computeRegionStructure();
        processor_->setStructureInvalid(false);
    }

    interruptionPoint();

    if (processor_->geometryInvalid()) {

        geometry_->clear();

        interruptionPoint();

        if (debugOutput_) {
            stopWatch_.reset();
            stopWatch_.start();
        }

        //compute pre-integration tables for every TF
        std::vector<const PreIntegrationTable*> piTables;
        for (int i = 0; static_cast<size_t>(i) < volume_->getNumChannels(); ++i) {
            piTables.push_back(tfCopyVector_.at(i)->getPreIntegrationTable(1.f, 256));
        }

        if (debugOutput_) {
            stopWatch_.stop();
            std::cout << "Fetching PreIntegration tables took " << stopWatch_.getRuntime() << " milliseconds" << std::endl;

            stopWatch_.reset();
            stopWatch_.start();
        }

        interruptionPoint();

        //scan through region structure, classify every region and add the bounds of non-transparent blocks to the minimal bounding box
        tgt::Bounds minBounds;

        std::vector<ProxyGeometryVolumeRegion>::const_iterator i;
        for (i = volumeStructure_->begin(); i != volumeStructure_->end(); ++i) {
            interruptionPoint();
            //check every channel if the region is opaque
            for (int channel = 0; static_cast<size_t>(channel) < volume_->getNumChannels(); ++channel) {
                //apply tf domain
                float minIntensity = tfCopyVector_.at(channel)->realWorldToNormalized(i->getMinIntensity(channel));
                float maxIntensity = tfCopyVector_.at(channel)->realWorldToNormalized(i->getMaxIntensity(channel));
                if (!isRegionEmptyPi(minIntensity,maxIntensity,piTables.at(channel))) {
                    minBounds.addVolume(i->getBounds());
                    break;
                }
            }
        }

        if (debugOutput_) {
            stopWatch_.stop();
            std::cout << "Determined (approximate) minimal bounding box in " << stopWatch_.getRuntime() << " milliseconds" << std::endl;
        }

        interruptionPoint();

        if (minBounds.isDefined()) {
            //now create the bounding box that has been found
            if (clippingEnabled_) {
                //get clipping planes
                tgt::Bounds clipBounds(clipLlf_, clipUrb_);
                // we have to use the actual clip bounds which extend 0.5 voxels more at the border for the checks, but use the original clip bounds for calling addCubeMeshClipped
                tgt::Bounds checkBounds(clipLlf_ - tgt::vec3(0.5f), clipUrb_ + tgt::vec3(0.5f));
                //only add and clip mesh if necessary
                if (checkBounds.containsVolume(minBounds))
                    OptimizedProxyGeometry::addCubeMesh(geometry_, minBounds, dim);
                else if (checkBounds.intersects(minBounds))
                    OptimizedProxyGeometry::addCubeMeshClip(geometry_, minBounds,dim,clipBounds);
            }
            else
                OptimizedProxyGeometry::addCubeMesh(geometry_, minBounds,dim);
        }

        // coordinates have been set in voxel space -> set voxelToWorld matrix as transformation
        geometry_->setTransformationMatrix(volume_->getTextureToWorldMatrix());
        processor_->setGeometryInvalid(false);

        interruptionPoint();

        if (debugOutput_)
            std::cout << "Created (approximate) minimal bounding box proxy geometry, volume "
                      << ftos(minBounds.volume()) << " (before clipping)" << std::endl;
    }
}

//-------------------------------------------------------------------------------------------------
// VisibleBricksBackgroundThread

VisibleBricksBackgroundThread::VisibleBricksBackgroundThread(OptimizedProxyGeometry* processor, const VolumeBase* volume, std::vector<TransFunc1D*> tfVector,
                                                            float threshold, TriangleMeshGeometryColorNormal* geometry,
                                                            std::vector<ProxyGeometryVolumeRegion>* volumeStructure,
                                                            tgt::ivec3 volStructureSize, int stepSize, bool debugOutput,
                                                            bool clippingEnabled, tgt::vec3 clipLlf, tgt::vec3 clipUrb)
        : StructureProxyGeometryBackgroundThread(processor, volume, tfVector, threshold, geometry, volumeStructure, volStructureSize, stepSize,
                                                debugOutput, clippingEnabled, clipLlf, clipUrb)
{}

VisibleBricksBackgroundThread::~VisibleBricksBackgroundThread() {
    //wait for internal thread to finish
    join();
}

void VisibleBricksBackgroundThread::threadMain() {
    computeMaximalBricks();
}

void VisibleBricksBackgroundThread::computeMaximalBricks() {

    //get volume dimensions
    tgt::ivec3 dim = volume_->getDimensions();

    interruptionPoint();

    if (processor_->structureInvalid()) {
        //invalidate geometry
        processor_->setGeometryInvalid();

        //clear structure
        volumeStructure_->clear();

        interruptionPoint();

        //compute new structure
        computeRegionStructure();
        processor_->setStructureInvalid(false);
    }

    interruptionPoint();

    if (processor_->geometryInvalid()) {

        geometry_->clear();

        interruptionPoint();

        if (debugOutput_) {
            stopWatch_.reset();
            stopWatch_.start();
        }

        //compute pre-integration tables for every TF
        std::vector<const PreIntegrationTable*> piTables;
        for (int i = 0; static_cast<size_t>(i) < volume_->getNumChannels(); ++i) {
            piTables.push_back(tfCopyVector_.at(i)->getPreIntegrationTable(1.f, 256));
        }

        if (debugOutput_) {
            stopWatch_.stop();
            std::cout << "Fetching PreIntegration tables took " << stopWatch_.getRuntime() << " milliseconds" << std::endl;
            stopWatch_.reset();
            stopWatch_.start();
        }

        setVolBound(tgt::ivec3(0),volStructureSize_-1,false);

        interruptionPoint();

        int numberOfCubes = 0; //number of created cubes
        float proxyVolume = 0.f;

        tgt::ivec3 pos, urbVol;
        tgt::vec3 llf, urb;
        for (VolumeIterator it(volStructureSize_); !it.outofrange(); it.next()) {
            interruptionPoint();
            pos = it.value();
            if (isVolNotEmptyPiNotBound(pos,tfCopyVector_, piTables)) {
                urbVol = getUrbPi(pos,tfCopyVector_, piTables);
                llf = getVolumeRegion(pos).getBounds().getLLF();
                urb = getVolumeRegion(urbVol).getBounds().getURB();

                interruptionPoint();

                tgt::Bounds cubeBounds(llf,urb);
                if (clippingEnabled_) {
                    tgt::Bounds clipBounds(clipLlf_, clipUrb_);
                    // we have to use the actual clip bounds which extend 0.5 voxels more at the border for the checks, but use the original clip bounds for calling addCubeMeshClipped
                    tgt::Bounds checkBounds(clipLlf_ - tgt::vec3(0.5f), clipUrb_ + tgt::vec3(0.5f));
                    //only add and clip mesh if necessary
                    if (checkBounds.containsVolume(cubeBounds))
                        OptimizedProxyGeometry::addCubeMesh(geometry_, cubeBounds, dim);
                    else if (checkBounds.intersects(cubeBounds))
                        OptimizedProxyGeometry::addCubeMeshClip(geometry_, cubeBounds,dim,clipBounds);
                }
                else
                    OptimizedProxyGeometry::addCubeMesh(geometry_, cubeBounds, dim);

                if (debugOutput_) {
                    proxyVolume += cubeBounds.volume();
                    numberOfCubes++;
                }

            }
        }

        if (debugOutput_) {
            stopWatch_.stop();
            std::cout << "Created maximal cubes proxy geometry in " << stopWatch_.getRuntime() << " milliseconds" << std::endl;
            std::cout << "Created Proxy Geometry consisting of " << numberOfCubes << " cubes using maximal cubes mode," << std::endl;
            std::cout << " volume < " << proxyVolume << " (before clipping)" << std::endl;
        }

        interruptionPoint();

        // geometry is specified in texture space -> set voxelToWorld matrix as transformation
        geometry_->setTransformationMatrix(volume_->getTextureToWorldMatrix());
        processor_->setGeometryInvalid(false);
    }
}


void VisibleBricksBackgroundThread::setVolBound(tgt::ivec3 llf, tgt::ivec3 urb, bool value) {
    for (int z=llf.z; z<=urb.z;z++) {
        interruptionPoint();
        for (int y=llf.y; y<=urb.y; y++) {
            for (int x=llf.x; x<=urb.x; x++) {
                ProxyGeometryVolumeRegion& v = getVolumeRegion(tgt::ivec3(x,y,z));
                v.setBound(value);
            }
        }
    }
}

bool VisibleBricksBackgroundThread::isVolNotEmptyPiNotBound(tgt::ivec3 pos, const std::vector<TransFunc1D*>& tfs, const std::vector<const PreIntegrationTable*>& piTables) {
    ProxyGeometryVolumeRegion& v = getVolumeRegion(pos);
    if (v.isBound())
        return false;

    //check every channel
    bool empty = true;

    for (int channel = 0; static_cast<size_t>(channel) < piTables.size(); ++channel) {
        //apply tf domain
        float minIntensity = tfs.at(channel)->realWorldToNormalized(v.getMinIntensity(channel));
        float maxIntensity = tfs.at(channel)->realWorldToNormalized(v.getMaxIntensity(channel));
        if (!isRegionEmptyPi(minIntensity,maxIntensity,piTables.at(channel))) {
            empty = false;
            break;
        }
     }

     return (!empty);
}

bool VisibleBricksBackgroundThread::isVolNotEmptyPiNotBound(tgt::ivec3 llf, tgt::ivec3 urb,
        const std::vector<TransFunc1D*>& tfs, const std::vector<const PreIntegrationTable*>& piTables)
{
    for (int z=llf.z; z<=urb.z;z++) {
        interruptionPoint();
        for (int y=llf.y; y<=urb.y; y++) {
            for (int x=llf.x; x<=urb.x; x++) {
                ProxyGeometryVolumeRegion& v = getVolumeRegion(tgt::ivec3(x,y,z));
                if (v.isBound())
                    return false;

                //check every channel
                bool empty = true;

                for (int channel = 0; static_cast<size_t>(channel) < piTables.size(); ++channel) {
                    //apply tf domain
                    float minIntensity = tfs.at(channel)->realWorldToNormalized(v.getMinIntensity(channel));
                    float maxIntensity = tfs.at(channel)->realWorldToNormalized(v.getMaxIntensity(channel));
                    if (!isRegionEmptyPi(minIntensity,maxIntensity,piTables.at(channel))) {
                        empty = false;
                        break;
                    }
                }

                if (empty)
                    return false;
            }
        }
    }

    return true;
}

tgt::ivec3 VisibleBricksBackgroundThread::getUrbPi(tgt::ivec3 llf, const std::vector<TransFunc1D*>& tfs,
        const std::vector<const PreIntegrationTable*>& piTables)
{

    interruptionPoint();

    const tgt::ivec3 size = volStructureSize_;
    tgt::bvec3 inc(true);
    tgt::ivec3 urb(llf);

    while(inc.x||inc.y||inc.z) {
        interruptionPoint();
        if (inc.x) {
            if (urb.x+1>size.x-1) {
                inc.x=false;
            } else {
                if (isVolNotEmptyPiNotBound(
                    tgt::ivec3(urb.x+1,llf.y,llf.z),
                    tgt::ivec3(urb.x+1,urb.y,urb.z),
                    tfs, piTables))
                {
                    urb.x += 1;
                } else {
                    inc.x=false;
                }
            }
        }
        if (inc.y) {
            if (urb.y+1>size.y-1) {
                inc.y=false;
            } else {
                if (isVolNotEmptyPiNotBound(
                    tgt::ivec3(llf.x,urb.y+1,llf.z),
                    tgt::ivec3(urb.x,urb.y+1,urb.z),
                    tfs, piTables)) {
                    urb.y += 1;
                } else {
                    inc.y=false;
                }
            }
        }
        if (inc.z) {
            if (urb.z+1>size.z-1) {
                inc.z=false;
            } else {
                if (isVolNotEmptyPiNotBound(
                    tgt::ivec3(llf.x,llf.y,urb.z+1),
                    tgt::ivec3(urb.x,urb.y,urb.z+1),
                    tfs, piTables)) {
                    urb.z += 1;
                } else {
                    inc.z=false;
                }
            }
        }
    }

    setVolBound(llf,urb);
    return urb;
}


//-------------------------------------------------------------------------------------------------
// OuterFacesBackgroundThread

OuterFacesBackgroundThread::OuterFacesBackgroundThread(OptimizedProxyGeometry* processor, const VolumeBase* volume, std::vector<TransFunc1D*> tfVector,
                                                            float threshold, TriangleMeshGeometryColorNormal* geometry,
                                                            std::vector<ProxyGeometryVolumeRegion>* volumeStructure,
                                                            tgt::ivec3 volStructureSize, int stepSize, bool debugOutput,
                                                            bool clippingEnabled, tgt::vec3 clipLlf, tgt::vec3 clipUrb)
        : StructureProxyGeometryBackgroundThread(processor, volume, tfVector, threshold, geometry, volumeStructure, volStructureSize, stepSize,
                                                debugOutput, clippingEnabled, clipLlf, clipUrb)
{}

OuterFacesBackgroundThread::~OuterFacesBackgroundThread() {
    //wait for internal thread to finish
    join();
}

void OuterFacesBackgroundThread::threadMain() {
    computeOuterFaces();
}

void OuterFacesBackgroundThread::computeOuterFaces() {

    interruptionPoint();

    if (processor_->structureInvalid()) {
        //invalidate geometry
        processor_->setGeometryInvalid();

        //clear structure
        volumeStructure_->clear();

        interruptionPoint();

        //compute new structure
        computeRegionStructure();
        processor_->setStructureInvalid(false);
    }

    interruptionPoint();

    if (processor_->geometryInvalid()) {

        geometry_->clear();

        interruptionPoint();

        if (debugOutput_) {
            stopWatch_.reset();
            stopWatch_.start();
        }

        //compute pre-integration tables for every TF
        std::vector<const PreIntegrationTable*> piTables;
        for (int i = 0; static_cast<size_t>(i) < volume_->getNumChannels(); ++i) {
            piTables.push_back(tfCopyVector_.at(i)->getPreIntegrationTable(1.f, 256));
        }

        if (debugOutput_) {
            stopWatch_.stop();
            std::cout << "Fetching PreIntegration table took " << stopWatch_.getRuntime() << " milliseconds" << std::endl;
            stopWatch_.reset();
            stopWatch_.start();
        }

        interruptionPoint();

        //first pass: scan through region structure and classify every region
        std::vector<ProxyGeometryVolumeRegion>::iterator i;
        for (i = volumeStructure_->begin(); i != volumeStructure_->end(); ++i) {
            interruptionPoint();

            //check every channel
            bool empty = true;

            for (int channel = 0; static_cast<size_t>(channel) < piTables.size(); ++channel) {
                //apply tf domain
                float minIntensity = tfCopyVector_.at(channel)->realWorldToNormalized(i->getMinIntensity(channel));
                float maxIntensity = tfCopyVector_.at(channel)->realWorldToNormalized(i->getMaxIntensity(channel));
                if (!isRegionEmptyPi(minIntensity,maxIntensity,piTables.at(channel))) {
                    empty = false;
                    break;
                }
            }

            i->setOpaque(!empty);
        }

        //second pass: scan through region structure
        //for every opaque region and for every side of the brick check the neighbors and add a face if the neighbor is transparent
        tgt::ivec3 pos;
        for (VolumeIterator it(volStructureSize_); !it.outofrange(); it.next()) {

            interruptionPoint();

            tgt::Bounds clipBounds(clipLlf_, clipUrb_);
            // we have to use the actual clip bounds which extend 0.5 voxels more at the border for the checks, but use the original clip bounds for calling addCubeMeshClipped
            tgt::Bounds checkBounds(clipLlf_ - tgt::vec3(0.5f), clipUrb_ + tgt::vec3(0.5f));

            //get position of current brick and a reference to the brick
            pos = it.value();
            ProxyGeometryVolumeRegion& current  = getVolumeRegion(pos);

            //if the current brick is transparent no faces are generated
            if (!current.isOpaque())
                continue;

            //if clipping is enabled and the current brick is completely outside the clipping area no faces are generated
            if (clippingEnabled_ && !checkBounds.containsVolume(current.getBounds()) && !checkBounds.intersects(current.getBounds()))
                continue;

            //outer loop: set sign of the normal from -1 to +1
            for (int sign = -1; sign <= 1; sign += 2) {
                //inner loop: iterate over the three dimensions
                for (int dim = 0; dim < 3; ++dim) {

                    interruptionPoint();

                    if (!clippingEnabled_ || checkBounds.containsVolume(current.getBounds())) {
                        //no clipping
                        //if direct neighbor in direction sign*dim is existing and opaque: continue
                        tgt::ivec3 neighborPos = pos;
                        neighborPos.elem[dim] += sign;

                        if ((neighborPos.elem[dim] >= 0) && (neighborPos.elem[dim] < volStructureSize_.elem[dim]) && getVolumeRegion(neighborPos).isOpaque())
                            continue;

                        //else: create a face in direction sign*dim

                        //get the coordinates of the current block and convert them to texture space
                        tgt::vec3 coordllf = volume_->getVoxelToTextureMatrix() * current.getBounds().getLLF();
                        tgt::vec3 coordurb = volume_->getVoxelToTextureMatrix() * current.getBounds().getURB();

                        // four vertices for the quad
                        tgt::vec3 vertices[4];

                        // set the coordinate of the face side according to dim and sign
                        for (int i = 0; i < 4; ++i) {
                            vertices[i].elem[dim] = (sign < 0) ? coordllf.elem[dim] : coordurb.elem[dim];
                        }

                        //set the other coordinates
                        for (int incr = 1; incr <= 2; ++incr) {
                            int curDim = (dim + incr) % 3;

                            vertices[0].elem[curDim] = coordllf.elem[curDim];
                            vertices[1].elem[curDim] = (incr == 1) ? coordurb.elem[curDim] : coordllf.elem[curDim];
                            vertices[2].elem[curDim] = coordurb.elem[curDim];
                            vertices[3].elem[curDim] = (incr == 1) ? coordllf.elem[curDim] : coordurb.elem[curDim];
                        }

                        // texture coordinates (second value) are computed like the transformation of the voxelToTexture-matrix (see VolumeBase)
                        VertexNormal ll(vertices[0], vertices[0]);
                        VertexNormal lr(vertices[1], vertices[1]);
                        VertexNormal ur(vertices[2], vertices[2]);
                        VertexNormal ul(vertices[3], vertices[3]);

                        //counter-clockwise order of vertices
                        if (sign > 0)
                            geometry_->addQuad(ll, lr, ur, ul);
                        else
                            geometry_->addQuad(ll, ul, ur, lr);

                    }
                    else {
                        //clipping needed

                        //if direct neighbor in direction sign*dim is existing and opaque AND is not completely outside the clipping area: continue
                        tgt::ivec3 neighborPos = pos;
                        neighborPos.elem[dim] += sign;

                        if ((neighborPos.elem[dim] >= 0) && (neighborPos.elem[dim] < volStructureSize_.elem[dim]) && getVolumeRegion(neighborPos).isOpaque()
                                && (checkBounds.containsVolume(getVolumeRegion(neighborPos).getBounds()) || checkBounds.intersects(getVolumeRegion(neighborPos).getBounds())))
                            continue;

                        //else: create a face in direction sign*dim

                        //get the coordinates of the current block
                        tgt::vec3 llf = current.getBounds().getLLF();
                        tgt::vec3 urb = current.getBounds().getURB();

                        //clip the coordinates and transform to texture coordinates
                        tgt::vec3 coordllf = volume_->getVoxelToTextureMatrix() * tgt::max(llf, clipLlf_ - tgt::vec3(0.5f));
                        tgt::vec3 coordurb = volume_->getVoxelToTextureMatrix() * tgt::min(urb, clipUrb_ + tgt::vec3(0.5f));

                        // four vertices for the quad
                        tgt::vec3 vertices[4];

                        // set the coordinate of the face side according to dim and sign
                        for (int i = 0; i < 4; ++i) {
                            vertices[i].elem[dim] = (sign < 0) ? coordllf.elem[dim] : coordurb.elem[dim];
                        }

                        //set the other coordinates
                        for (int incr = 1; incr <= 2; ++incr) {
                            int curDim = (dim + incr) % 3;

                            vertices[0].elem[curDim] = coordllf.elem[curDim];
                            vertices[1].elem[curDim] = (incr == 1) ? coordurb.elem[curDim] : coordllf.elem[curDim];
                            vertices[2].elem[curDim] = coordurb.elem[curDim];
                            vertices[3].elem[curDim] = (incr == 1) ? coordllf.elem[curDim] : coordurb.elem[curDim];
                        }

                        // texture coordinates (second value) are computed like the transformation of the voxelToTexture-matrix (see VolumeBase)
                        VertexNormal ll(vertices[0], vertices[0]);
                        VertexNormal lr(vertices[1], vertices[1]);
                        VertexNormal ur(vertices[2], vertices[2]);
                        VertexNormal ul(vertices[3], vertices[3]);

                        //counter-clockwise order of vertices
                        if (sign > 0)
                            geometry_->addQuad(ll, lr, ur, ul);
                        else
                            geometry_->addQuad(ll, ul, ur, lr);
                    }
                }
            }
        }

        if (debugOutput_) {
            stopWatch_.stop();
            std::cout << "Created outer faces proxy geometry in " << stopWatch_.getRuntime() << " milliseconds" << std::endl;
        }

        interruptionPoint();

        // geometry is specified in voxel space
        geometry_->setTransformationMatrix(volume_->getTextureToWorldMatrix());
        processor_->setGeometryInvalid(false);
    }
}

//-------------------------------------------------------------------------------------------------
// VolumeOctreeBackgroundThread

VolumeOctreeBackgroundThread::VolumeOctreeBackgroundThread(OptimizedProxyGeometry* processor, const VolumeBase* volume, std::vector<TransFunc1D*> tfVector,
            float threshold, TriangleMeshGeometryColorNormal* geometry, int stepSize, bool clippingEnabled, tgt::vec3 clipLlf, tgt::vec3 clipUrb)
        : OptimizedProxyGeometryBackgroundThread(processor, volume, tfVector, threshold, geometry, stepSize, false, clippingEnabled, clipLlf, clipUrb)
{}

VolumeOctreeBackgroundThread::~VolumeOctreeBackgroundThread() {
    //wait for internal thread to finish
    join();
}

void VolumeOctreeBackgroundThread::threadMain() {
    computeVolumeOctreeGeometry();
}

void VolumeOctreeBackgroundThread::computeVolumeOctreeGeometry() {
    interruptionPoint();

    //TODO: check if computation is necessary?!
    processor_->setGeometryInvalid(true);

    if (processor_->geometryInvalid()) {
        geometry_->clear();

        tgt::svec3 volumeDim = volume_->getDimensions();

        //compute pre-integration tables for every TF
        std::vector<const PreIntegrationTable*> piTables;
        for (int i = 0; static_cast<size_t>(i) < volume_->getNumChannels(); ++i) {
            piTables.push_back(tfCopyVector_.at(i)->getPreIntegrationTable(1.f, 256));
        }

        interruptionPoint();

        tgt::Bounds clipBounds(clipLlf_, clipUrb_);

        //get VolumeOctree representation and traverse it to create the geometry
        const VolumeOctreeBase* octree = volume_->getRepresentation<VolumeOctreeBase>();
        const VolumeOctreeNode* root = octree->getRootNode();
        int numChannels = static_cast<int>(volume_->getNumChannels());

        tgt::svec3 octreeDim = octree->getOctreeDim();

        traverseOctreeAndCreateGeometry(root, tgt::vec3(0.f), tgt::vec3(octreeDim - tgt::svec3::one), tfCopyVector_, numChannels, piTables, clipBounds, volumeDim);

        interruptionPoint();

        geometry_->setTransformationMatrix(volume_->getTextureToWorldMatrix());
        processor_->setGeometryInvalid(false);
    }

    interruptionPoint();
}

void VolumeOctreeBackgroundThread::traverseOctreeAndCreateGeometry(const VolumeOctreeNode* node, const tgt::vec3& nodeLlf, const tgt::vec3& nodeUrb,
            /*TransFunc1DKeys* tf, const PreIntegrationTable* piTable,*/
            std::vector<TransFunc1D*> tfVector, int numVolumeChannels, std::vector<const PreIntegrationTable*> piTables, tgt::Bounds clipBounds,
            tgt::ivec3 volumeDim)
{

    interruptionPoint();

    if (node->isLeaf()) {
        //compute intensity values
        RealWorldMapping rwm = volume_->getRealWorldMapping();

        //check if any of the channels is not transparent
        bool isOpaque = false;
        for (int i = 0; i < numVolumeChannels; ++i) {
            // transform normalized voxel intensities to real-world
            float minRealWorld = rwm.normalizedToRealWorld(node->getMinValue(i) / 65535.f);
            float maxRealWorld = rwm.normalizedToRealWorld(node->getMaxValue(i) / 65535.f);

            float minIntensity = tfVector.at(i)->realWorldToNormalized(minRealWorld);
            float maxIntensity = tfVector.at(i)->realWorldToNormalized(maxRealWorld);

            //check pre-integration heuristic
            if (!isRegionEmptyPi(minIntensity, maxIntensity, piTables.at(i))) {
                isOpaque = true;
                break;
            }
        }

        //do not render if all channels are transparent
        if (!isOpaque)
            return;

        //check if node is outside volume dimensions
        if (tgt::hor(tgt::greaterThanEqual(nodeLlf, tgt::vec3(volumeDim))))
            return;

        //clamp node to volume dimensions
        tgt::vec3 correctNodeUrb = tgt::min(nodeUrb, tgt::vec3(volumeDim - tgt::ivec3::one));

        //not transparent -> create cube (if necessery: clipping)
        //since the octree node does not take into account neighboring voxels we have to stretch the coordinates to one voxel further to account for trilinear interpolation at proxy geometry boundaries
        //we do not do this for the volume boundaries since it would throw off our bounding box if the volume is not empty at the border
        tgt::vec3 actualNodeLLF = tgt::max(nodeLlf - tgt::vec3(1.f), tgt::vec3(-0.5f));
        tgt::vec3 actualNodeURB = tgt::min(correctNodeUrb + tgt::vec3(1.f), tgt::vec3(volumeDim) - tgt::vec3(0.5f));
        tgt::Bounds nodeBounds(actualNodeLLF, actualNodeURB);    //get bounding box

        // account for offsets in clipping (see bounding box computation)
        tgt::Bounds checkBounds(clipBounds.getLLF() - tgt::vec3(0.5f), clipBounds.getURB() + tgt::vec3(0.5f));

        if (clippingEnabled_) {
            if (checkBounds.containsVolume(nodeBounds))
                OptimizedProxyGeometry::addCubeMesh(geometry_, nodeBounds, volumeDim);
            else if (checkBounds.intersects(nodeBounds))
                OptimizedProxyGeometry::addCubeMeshClip(geometry_, nodeBounds, volumeDim, clipBounds);
         }
         else
            OptimizedProxyGeometry::addCubeMesh(geometry_, nodeBounds, volumeDim);

    } else {

        //for every child: compute bounding box and traverse further
        tgt::vec3 nodeHalfDim = (nodeUrb-nodeLlf + tgt::vec3::one) / 2.f;
        for (size_t childID=0; childID<8; childID++) {
            const VolumeOctreeNode* childNode = node->children_[childID];
            if (childNode) {
                tgt::svec3 childCoord = linearCoordToCubic(childID, tgt::svec3::two);
                tgt::vec3 childLlf = nodeLlf + tgt::vec3(childCoord)*nodeHalfDim;
                tgt::vec3 childUrb = childLlf + nodeHalfDim - tgt::vec3::one;

                //traverseOctreeAndCreateGeometry(childNode, childLlf, childUrb, tf, piTable, clipBounds, volumeDim);
                traverseOctreeAndCreateGeometry(childNode, childLlf, childUrb, tfVector, numVolumeChannels, piTables, clipBounds, volumeDim);
            }
        }
    }
}

//-------------------------------------------------------------------------------------------------
// VolumeOctreeOuterFacesBackgroundThread

VolumeOctreeOuterFacesBackgroundThread::VolumeOctreeOuterFacesBackgroundThread(OptimizedProxyGeometry* processor, const VolumeBase* volume, std::vector<TransFunc1D*> tfVector,
            float threshold, TriangleMeshGeometryColorNormal* geometry, int stepSize, bool clippingEnabled, tgt::vec3 clipLlf, tgt::vec3 clipUrb)
        : OptimizedProxyGeometryBackgroundThread(processor, volume, tfVector, threshold, geometry, stepSize, false, clippingEnabled, clipLlf, clipUrb)
{}

VolumeOctreeOuterFacesBackgroundThread::~VolumeOctreeOuterFacesBackgroundThread() {
    //wait for internal thread to finish
    join();
}

void VolumeOctreeOuterFacesBackgroundThread::threadMain() {
    computeVolumeOctreeGeometry();
}

void VolumeOctreeOuterFacesBackgroundThread::computeVolumeOctreeGeometry() {
    interruptionPoint();

    //TODO: check if computation is necessary?!
    processor_->setGeometryInvalid(true);

    if (processor_->geometryInvalid()) {
        geometry_->clear();

        tgt::svec3 volumeDim = volume_->getDimensions();

        //compute pre-integration tables for every TF
        std::vector<const PreIntegrationTable*> piTables;
        for (int i = 0; static_cast<size_t>(i) < volume_->getNumChannels(); ++i) {
            piTables.push_back(tfCopyVector_.at(i)->getPreIntegrationTable(1.f, 256));
        }

        interruptionPoint();

        tgt::Bounds clipBounds(clipLlf_, clipUrb_);

        //get VolumeOctree representation and some meta information
        const VolumeOctreeBase* octree = volume_->getRepresentation<VolumeOctreeBase>();
        const VolumeOctreeNode* root = octree->getRootNode();
        int numChannels = static_cast<int>(volume_->getNumChannels());

        tgt::svec3 octreeDim = octree->getOctreeDim();
        octreeDepth_ = octree->getActualTreeDepth();

        //compute number of bricks in every direction
        brickDim_ = octree->getBrickDim();
        volDim_ = octree->getVolumeDim();
        brickStructureSize_ = volDim_ / brickDim_;
        if (volDim_.x % brickDim_.x != 0)
            brickStructureSize_.x += 1;
        if (volDim_.y % brickDim_.y != 0)
            brickStructureSize_.y += 1;
        if (volDim_.z % brickDim_.z != 0)
            brickStructureSize_.z += 1;

        //now allocate memory for the brick structure
        //TODO: allocate that in the processor, set flags -> do not re-compute the bricks if only the TF has changed?!
        brickStructure_ = std::vector<ProxyGeometryVolumeRegion>(brickStructureSize_.x * brickStructureSize_.y * brickStructureSize_.z,
                ProxyGeometryVolumeRegion());

        //set to zero to mark all bricks as being not opaque
        //std::memset(&brickStructure_, 0, brickStructure_.capacity() * sizeof(ProxyGeometryVolumeRegion));

        //do the actual traversal to create the brick structure
        traverseOctreeAndCreateBrickStructure(root, 0, tgt::svec3((size_t) 0), octreeDim, tfCopyVector_, numChannels, piTables, clipBounds, volumeDim);

        interruptionPoint();

        //create the outer faces geometry from the brick structure
        computeOuterFaces();

        geometry_->setTransformationMatrix(volume_->getTextureToWorldMatrix());
        processor_->setGeometryInvalid(false);
    }

    interruptionPoint();
}

void VolumeOctreeOuterFacesBackgroundThread::traverseOctreeAndCreateBrickStructure (const VolumeOctreeNode* node, size_t level,
        const tgt::svec3& nodeLlf, const tgt::svec3& nodeUrb, std::vector<TransFunc1D*> tfVector, int numVolumeChannels,
        std::vector<const PreIntegrationTable*> piTables, tgt::Bounds clipBounds, tgt::ivec3 volumeDim)
{

    interruptionPoint();

    if (node->isLeaf()) {
        //compute intensity values
        RealWorldMapping rwm = volume_->getRealWorldMapping();

        //check if any of the channels is not transparent
        bool isOpaque = false;
        for (int i = 0; i < numVolumeChannels; ++i) {
            // transform normalized voxel intensities to real-world
            float minRealWorld = rwm.normalizedToRealWorld(node->getMinValue(i) / 65535.f);
            float maxRealWorld = rwm.normalizedToRealWorld(node->getMaxValue(i) / 65535.f);

            float minIntensity = tfVector.at(i)->realWorldToNormalized(minRealWorld);
            float maxIntensity = tfVector.at(i)->realWorldToNormalized(maxRealWorld);

            //check pre-integration heuristic
            if (!isRegionEmptyPi(minIntensity, maxIntensity, piTables.at(i))) {
                isOpaque = true;
                break;
            }
        }

        //do not change the corresponding brick(s) if all channels are transparent
        if (!isOpaque)
            return;

        //check if node is outside volume dimensions and discard it
        if (tgt::hor(tgt::greaterThanEqual(nodeLlf, tgt::svec3(volumeDim))))
            return;

        //set bricks in vector opaque
        tgt::svec3 brickLlf = nodeLlf / brickDim_;
        tgt::svec3 brickUrb = nodeUrb / brickDim_ - tgt::svec3(1);
        setBricksVisible(brickLlf, brickUrb);

        //tgt::Bounds b(tgt::vec3(brickLlf), tgt::vec3(brickUrb));
        //std::cout << brickLlf << " " << brickUrb << std::endl;

    } else {

        //for every child: compute bounding box and traverse further
        tgt::svec3 tmpNodeDim = tgt::svec3(nodeUrb-nodeLlf);
        tgt::svec3 nodeHalfDim = tgt::svec3(tmpNodeDim.x / 2, tmpNodeDim.y / 2, tmpNodeDim.z / 2);
        for (size_t childID=0; childID<8; childID++) {
            const VolumeOctreeNode* childNode = node->children_[childID];
            if (childNode) {
                tgt::svec3 childCoord = linearCoordToCubic(childID, tgt::svec3::two);
                tgt::svec3 childLlf = nodeLlf + childCoord*nodeHalfDim;
                tgt::svec3 childUrb = childLlf + nodeHalfDim;

                //traverseOctreeAndCreateGeometry(childNode, childLlf, childUrb, tf, piTable, clipBounds, volumeDim);
                traverseOctreeAndCreateBrickStructure(childNode, level + 1, childLlf, childUrb, tfVector, numVolumeChannels, piTables, clipBounds, volumeDim);
            }
        }
    }
}

void VolumeOctreeOuterFacesBackgroundThread::setBricksVisible(tgt::svec3 llf, tgt::svec3 urb) {
    for (size_t z=llf.z; z<=urb.z && z < brickStructureSize_.z;z++) {
        interruptionPoint();
        for (size_t y=llf.y; y<=urb.y && y < brickStructureSize_.y; y++) {
            for (size_t x=llf.x; x<=urb.x && x < brickStructureSize_.x; x++) {
                // we need one voxel overlap because the octree does not check neighboring voxels
                tgt::vec3 rllf = tgt::vec3(brickDim_ * tgt::svec3(x,y,z)) - tgt::vec3(1.f);
                tgt::vec3 rurb = tgt::vec3((brickDim_ * tgt::svec3(x+1,y+1,z+1)) - tgt::svec3(1)) + tgt::vec3(1.f);
                tgt::Bounds regionBounds(rllf, rurb);
                //std::cout << rllf << " " << rurb << std::endl;
                std::vector<tgt::vec2> minMaxIntensities;       //<- not needed here, may be empty
                ProxyGeometryVolumeRegion& v = getVolumeRegion(tgt::ivec3((int)x,(int)y,(int)z));
                v = ProxyGeometryVolumeRegion(regionBounds,minMaxIntensities);
            }
        }
    }
}

ProxyGeometryVolumeRegion& VolumeOctreeOuterFacesBackgroundThread::getVolumeRegion(tgt::ivec3 pos) {
    return brickStructure_.at(pos.z * (brickStructureSize_.x * brickStructureSize_.y) + pos.y * brickStructureSize_.x + pos.x);
}

void VolumeOctreeOuterFacesBackgroundThread::computeOuterFaces() {

    //get volume dimensions
    tgt::ivec3 volDim = volume_->getDimensions();

    interruptionPoint();

    //second pass: scan through region structure
    //for every opaque region and for every side of the brick check the neighbors and add a face if the neighbor is transparent
    tgt::ivec3 pos;
    for (VolumeIterator it(brickStructureSize_); !it.outofrange(); it.next()) {

        interruptionPoint();

        tgt::Bounds clipBounds(clipLlf_, clipUrb_);
        // account for offsets in clipping (see bounding box computation)
        tgt::Bounds checkBounds(clipBounds.getLLF() - tgt::vec3(0.5f), clipBounds.getURB() + tgt::vec3(0.5f));

        //get position of current brick and a reference to the brick
        pos = it.value();
        ProxyGeometryVolumeRegion& current  = getVolumeRegion(pos);

        //if the current brick is transparent no faces are generated
        if (!current.isOpaque())
            continue;

        //if clipping is enabled and the current brick is completely outside the clipping area no faces are generated
        if (clippingEnabled_ && !checkBounds.containsVolume(current.getBounds()) && !checkBounds.intersects(current.getBounds()))
            continue;

        //outer loop: set sign of the normal from -1 to +1
        for (int sign = -1; sign <= 1; sign += 2) {
            //inner loop: iterate over the three dimensions
            for (int dim = 0; dim < 3; ++dim) {

                interruptionPoint();

                if (!clippingEnabled_ || checkBounds.containsVolume(current.getBounds())) {
                    //no clipping
                    //if direct neighbor in direction sign*dim is existing and opaque: continue
                    tgt::ivec3 neighborPos = pos;
                    neighborPos.elem[dim] += sign;

                    if ((neighborPos.elem[dim] >= 0) && (static_cast<size_t>(neighborPos.elem[dim]) < brickStructureSize_.elem[dim]) && getVolumeRegion(neighborPos).isOpaque())
                        continue;

                    //else: create a face in direction sign*dim

                    //get the coordinates of the current block
                    tgt::vec3 coordllf = current.getBounds().getLLF();
                    tgt::vec3 coordurb = current.getBounds().getURB();

                    // one additional voxel has already been added, but we do not do this for the volume boundaries since it would throw off our bounding box if the volume is not empty at the border
                    coordllf = tgt::max(coordllf, tgt::vec3(-0.5f));
                    coordurb = tgt::min(coordurb, tgt::vec3(volDim) - tgt::vec3(0.5f));

                    coordllf = volume_->getVoxelToTextureMatrix() * coordllf;
                    coordurb = volume_->getVoxelToTextureMatrix() * coordurb;

                    // four vertices for the quad
                    tgt::vec3 vertices[4];

                    // set the coordinate of the face side according to dim and sign
                    for (int i = 0; i < 4; ++i) {
                        vertices[i].elem[dim] = (sign < 0) ? coordllf.elem[dim] : coordurb.elem[dim];
                    }

                    //set the other coordinates
                    for (int incr = 1; incr <= 2; ++incr) {
                        int curDim = (dim + incr) % 3;

                        vertices[0].elem[curDim] = coordllf.elem[curDim];
                        vertices[1].elem[curDim] = (incr == 1) ? coordurb.elem[curDim] : coordllf.elem[curDim];
                        vertices[2].elem[curDim] = coordurb.elem[curDim];
                        vertices[3].elem[curDim] = (incr == 1) ? coordllf.elem[curDim] : coordurb.elem[curDim];
                    }

                    VertexNormal ll(vertices[0], vertices[0]);
                    VertexNormal lr(vertices[1], vertices[1]);
                    VertexNormal ur(vertices[2], vertices[2]);
                    VertexNormal ul(vertices[3], vertices[3]);

                    //counter-clockwise order of vertices
                    if (sign > 0)
                        geometry_->addQuad(ll, lr, ur, ul);
                    else
                        geometry_->addQuad(ll, ul, ur, lr);

                }
                else {
                    //clipping needed

                    //if direct neighbor in direction sign*dim is existing and opaque AND is not completely outside the clipping area: continue
                    tgt::ivec3 neighborPos = pos;
                    neighborPos.elem[dim] += sign;

                    if ((neighborPos.elem[dim] >= 0) && (static_cast<size_t>(neighborPos.elem[dim]) < brickStructureSize_.elem[dim]) && getVolumeRegion(neighborPos).isOpaque()
                            && (checkBounds.containsVolume(getVolumeRegion(neighborPos).getBounds()) || checkBounds.intersects(getVolumeRegion(neighborPos).getBounds())))
                       continue;

                    //else: create a face in direction sign*dim

                    //get the coordinates of the current block
                    tgt::vec3 coordllf = current.getBounds().getLLF();
                    tgt::vec3 coordurb = current.getBounds().getURB();

                    // one additional voxel has already been added, but we do not do this for the volume boundaries since it would throw off our bounding box if the volume is not empty at the border
                    coordllf = tgt::max(coordllf, tgt::vec3(-0.5f));
                    coordurb = tgt::min(coordurb, tgt::vec3(volDim) - tgt::vec3(0.5f));

                    //clip the coordinates
                    coordllf = tgt::max(coordllf, clipLlf_ - tgt::vec3(0.5f));
                    coordurb = tgt::min(coordurb, clipUrb_ + tgt::vec3(0.5f));

                    coordllf = volume_->getVoxelToTextureMatrix() * coordllf;
                    coordurb = volume_->getVoxelToTextureMatrix() * coordurb;

                    // four vertices for the quad
                    tgt::vec3 vertices[4];

                    // set the coordinate of the face side according to dim and sign
                    for (int i = 0; i < 4; ++i) {
                        vertices[i].elem[dim] = (sign < 0) ? coordllf.elem[dim] : coordurb.elem[dim];
                    }

                    //set the other coordinates
                    for (int incr = 1; incr <= 2; ++incr) {
                        int curDim = (dim + incr) % 3;

                        vertices[0].elem[curDim] = coordllf.elem[curDim];
                        vertices[1].elem[curDim] = (incr == 1) ? coordurb.elem[curDim] : coordllf.elem[curDim];
                        vertices[2].elem[curDim] = coordurb.elem[curDim];
                        vertices[3].elem[curDim] = (incr == 1) ? coordllf.elem[curDim] : coordurb.elem[curDim];
                    }

                    VertexNormal ll(vertices[0], vertices[0]);
                    VertexNormal lr(vertices[1], vertices[1]);
                    VertexNormal ur(vertices[2], vertices[2]);
                    VertexNormal ul(vertices[3], vertices[3]);

                    //counter-clockwise order of vertices
                    if (sign > 0)
                        geometry_->addQuad(ll, lr, ur, ul);
                    else
                        geometry_->addQuad(ll, ul, ur, lr);
                }
            }
        }
    }

    if (debugOutput_) {
        stopWatch_.stop();
        std::cout << "Created outer faces proxy geometry in " << stopWatch_.getRuntime() << " milliseconds" << std::endl;
    }

    interruptionPoint();

    geometry_->setTransformationMatrix(volume_->getTextureToWorldMatrix());
    processor_->setGeometryInvalid(false);
}

} // namespace
