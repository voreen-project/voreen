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

#ifndef VRN_OPTIMIZEDPROXYGEOMETRY_H
#define VRN_OPTIMIZEDPROXYGEOMETRY_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/geometryport.h"

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/boundingboxproperty.h"

#include "voreen/core/datastructures/transfunc/1d/preintegrationtable.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"

#include "voreen/core/datastructures/volume/volume.h"

#include "tgt/stopwatch.h"

#include "voreen/core/utils/backgroundthread.h"

#include "voreen/core/datastructures/octree/volumeoctree.h"
#include "voreen/core/datastructures/octree/octreeutils.h"


namespace voreen {

    class TriangleMeshGeometryColorNormal;

/**
 * A volume region (brick) used for spatial subdivision.
 */
class ProxyGeometryVolumeRegion {
public:
    ProxyGeometryVolumeRegion()
        : bounds_(tgt::vec3(0.f), tgt::vec3(0.f))
        , minMaxIntensities_()
        , isBound_(false)
        , isOpaque_(false)
    { }

    ProxyGeometryVolumeRegion(const tgt::Bounds& bounds, const std::vector<tgt::vec2>& minMaxIntensities)
        : bounds_(bounds)
        , minMaxIntensities_(minMaxIntensities)
        , isBound_(false)
        , isOpaque_(true)
    { }

    tgt::Bounds getBounds() const { return bounds_; }
    void setBounds(tgt::Bounds b) {bounds_ = b;}
    float getMinIntensity(int channel = 0) const { return minMaxIntensities_.at(channel).x; }
    float getMaxIntensity(int channel = 0) const { return minMaxIntensities_.at(channel).y; }
    bool isBound() const { return isBound_;}
    void setBound(bool value) {isBound_ = value;}
    bool isOpaque() const { return isOpaque_;}
    void setOpaque(bool value) {isOpaque_ = value;}

private:
    tgt::Bounds bounds_;        ///< bounding box in voxel coordinates
    std::vector<tgt::vec2> minMaxIntensities_; ///< contains min and max intensity values within the bounding box for every channel
    bool isBound_;              ///< is this region already bound by a maximal cube
    bool isOpaque_;             ///< is this region opaque? (depends on TF)
};


//-------------------------------------------------------------------------------------------------

/**
 * Generates an optimized proxy geometry, allowing for empty-space leaping.
 * Provides several modes that compute the necessary data structures in a background thread while the processor uses a temporary bounding box geometry.
 * The proxy geometry also provides axis-aligned clipping.
 *
 * @see EntryExitPoints
 * @see MeshClipping
 */
class OptimizedProxyGeometry : public Processor, public VolumeObserver {

    friend class OptimizedProxyGeometryBackgroundThread;
    friend class StructureProxyGeometryBackgroundThread;
    friend class MinCubeBackgroundThread;
    friend class VisibleBricksBackgroundThread;
    friend class OuterFacesBackgroundThread;
    friend class VolumeOctreeBackgroundThread;
    friend class VolumeOctreeOuterFacesBackgroundThread;

public:
    OptimizedProxyGeometry();
    virtual ~OptimizedProxyGeometry();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "OptimizedProxyGeometry"; }
    virtual std::string getCategory() const  { return "Volume Proxy Geometry";  }
    virtual CodeState getCodeState() const   { return CODE_STATE_STABLE;        }

    virtual void volumeDelete(const VolumeBase* source);
    virtual void volumeChange(const VolumeBase* source);
    virtual void volumeRepresentationDelete(const VolumeBase* source, const VolumeRepresentation* rep);
    virtual void volumeDataDelete(const VolumeBase* source);

protected:

    virtual void setDescriptions() {
        setDescription("Generates an optimized proxy geometry, allowing for empty-space leaping."
                       "<p>See also: <a href=\"http://voreen.uni-muenster.de/?q=optimized-proxy-geometry\">http://voreen.uni-muenster.de/?q=optimized-proxy-geometry</a></p>"
                       "<b>Note</b>: The optimization uses the transfer function to determine the currently visible parts of the volume. "
                       "Therefore, the transfer function has to be linked with the volume renderer that uses the generated proxy geometry!");

        //Modes, pre-integration heuristic
        mode_.setDescription("Selects how the proxy geometry is built. There are four modes available: \
                <ul> \
                <li> <b> Bounding Box </b><br> \
                This mode simply uses a cube containing the whole volume. No empty space leaping or optimization is applied here. </li> \
                <li> <b> Minimal Visible Bounding Box </b><br> \
                This mode subdivides the volume into several bricks, where the size of the bricks is adjustable. Afterwards an approximate minimal bounding box is computed by adding all the bricks that are not entirely transparent to an initially empty bounding box. Transparency is determined by a pre-integration heuristic that compares the opacity of the pre-integration table entry with the minimum and maximum intensity values within the brick to a threshold as described by [Knoll et al., 2011].</li> \
                <li> <b> Visible Bricks </b><br> \
                This mode subdivides the volume into several bricks, where the size of the bricks is adjustable. Afterwards the visibility of the data within the bricks is determined (ie. some bricks are marked as transparent regions by the pre-integration heuristic described above). Then a proxy geometry is built by combining neighboring bricks to larger boxes. </li>  \
                <li> <b> Visible Bricks (Outer Faces) </b><br> \
                This mode subdivides the volume into several bricks, where the size of the bricks is adjustable. Afterwards the visibility of the data within the bricks is determined (ie. some bricks are marked as transparent regions by the pre-integration heuristic described above). Then a proxy geometry is built by examining every non-transparent brick and adding each of its faces that has a transparent brick as its direct neighbor.  </li>  \
                <li> <b> Volume Octree </b><br> \
                This mode recursively traverses an octree representation of the volume (which is constructed by an OctreeCreator processor in the network) and builds a proxy geometry by determining the visibility of the leaf nodes by the pre-integration heuristic described above. If no octree representation of the volume is available, the visible bricks mode is used.</li> \
                <li> <b> Volume Octree (Outer Faces) </b><br> \
                Similar to the Volume Octree mode, but builds up a list of visible bricks (corresponding to octree nodes) while traversing the octree. Afterwards, only the outer faces of those bricks contribute to the geometry, similar to the Outer Faces mode. If no octree representation of the volume is available, the outer faces mode is used. </li></ul> \
               Note that the data structures that are built for the second, third, fourth and fifth mode only have to be newly computed once the volume or the resolution (ie. the size of the bricks) changes. The proxy geometry computed by these modes is built by a background thread. During background computation, a cube proxy geometry is temporarily used.");
        tfChannel0_.setDescription("The transfer function that is used to identify transparent regions within the volume (channel 0 in multi-channel volumes).");
        resolutionMode_.setDescription("Determines how the resolution for subdividing the volume into bricks is set, either by subdivisions of the shortest edge of the volume or by the edge length of the bricks in voxels.");
        resolution_.setDescription("Determines the size of the bricks by dividing the size of the shortest edge of the volume.");
        resolutionVoxels_.setDescription("Determines the edge length of the bricks in voxels.");
        threshold_.setDescription("Determines the threshold for classifying bricks as empty (multiplied by 10^(-4)). Increasing this value might lead to better performance, but choosing the visibility threshold too high might remove features of the data set that should be visible.");
        //checkHalfNodes_.setDescription("When selected, octree traversal checks if half of the children of a partially opaque octree node (ie. upper/lower, right/left, back/front children) are all opaque and may therefore be rendered directly instead of traversing further (might be faster).");
        enableClipping_.setDescription("Selects if axis-aligned clipping is enabled for the proxy geometry.");
        waitForOptimization_.setDescription("If selected, no temporary cube is used while computing the optimized proxy geometry, but the processor waits for the computation to finish. Does not have an effect in Cube mode.");
    }

    virtual void process();

    /// callback for changes on input volume port
    virtual void adjustPropertiesToInput();

    /// adjust properties to changed or new volume (only call if an input volume is present)
    void onVolumeChange();

    /// Is the volume region structure invalid (needed for visible bricks / minimal visible bounding box modes)?
    bool structureInvalid() const;

    /// Is the mesh list geometry invalid (indicates that the geometry has to be rebuilt)?
    bool geometryInvalid() const;

    void setStructureInvalid(bool value = true);

    void setGeometryInvalid(bool value = true);

    /// Sets the number of bricks in every direction for the current volume structure.
    void setVolStructureSize(tgt::ivec3 size);

    /// creates a cube mesh given the llf and urb in voxel coordinates, dim is getVolumeDimensions() without minus 1!!!
    static void addCubeMesh(TriangleMeshGeometryColorNormal* mesh, tgt::Bounds bounds, tgt::ivec3 dim);

    /// analogously to addCubeMesh, but takes into account the clipping parameters, dim is getVolumeDimensions() without minus 1!!!
    static void addCubeMeshClip(TriangleMeshGeometryColorNormal* mesh, tgt::Bounds bounds, tgt::ivec3 dim, tgt::Bounds clipBounds);

private:

    /// is called from process() to create a cube proxy geometry without optimization
    void processCube();

    /// while a background thread is computing the proxy geometry, this sets a temporary cube geometry
    void processTmpCube();

    /// Ensure Backgroundthread updating
    void onClipRegionChanged();

    /// Resets clipping plane parameters to extremal values.
    void resetClipPlanes();

    /// Adapt ranges of clip plane properties to the input volume's dimensions.
    void adjustClipPropertiesRanges();

    /// Adjust the resolution property ranges according to the volume dimensions
    void adjustResolutionPropertyRanges();

    /// Adjust visibility of the clipping plane properties according to enableClipping_ property.
    void adjustClipPropertiesVisibility();

    /// Adjust visibility of the properties according to the mode
    void updatePropertyVisibility();

    /// Extracts ROI from volume (if present) and adjusts clipping slider accordingly.
    void adjustClippingToVolumeROI();

    virtual void onThresholdChange();           ///< called when threshold_ has changed
    virtual void onModeChange();                ///< called when mode has been changed
    virtual void onTransFuncChange();           ///< called when transfer function has been changed
    virtual void onResolutionChange();          ///< called when resolution has been changed
    virtual void onResolutionVoxelChange();     ///< called when the resolution in voxels has been changed
    virtual void onResolutionModeChange();      ///< changes resolution_ property according to the mode

    VolumePort inport_;                     ///< Inport for the dataset the proxy is generated for
    GeometryPort outport_;                  ///< Outport for the generated proxy geometry

    StringOptionProperty mode_;             ///< Mode for building the proxy geometry
    TransFunc1DKeysProperty tfChannel0_;          ///< Transfer function for determining opacity of regions
    TransFunc1DKeysProperty tfChannel1_;          ///< Transfer function for determining opacity of regions
    TransFunc1DKeysProperty tfChannel2_;          ///< Transfer function for determining opacity of regions
    TransFunc1DKeysProperty tfChannel3_;          ///< Transfer function for determining opacity of regions

    StringOptionProperty resolutionMode_;   ///< Determines, if resolution_ sets a) the number of subdivisions of the shortest side, or b) the edge length of the bricks in voxels
    IntProperty resolution_;                ///< determines often the shortest edge of the volume should be subdivided to compute the edge length of the bricks
    IntProperty resolutionVoxels_;          ///< determines the edge length of the bricks in voxels

    IntProperty threshold_;                 ///< determine the threshold for classifying regions as empty (used with pre-integration heuristic)

    BoolProperty enableClipping_;           ///< Clipping enable / disable property.
    IntBoundingBoxProperty clipRegion_;          ///< Clipping plane
    ButtonProperty resetClipPlanes_;        ///< Reset clipping plane parameters to extremal values.
    BoolProperty waitForOptimization_;      ///< if set, no temporary cube is used while computing the optimized proxy geometry, but the processor waits for the computation to finish

    bool structureInvalid_;  ///< determines if the spatial subdivision of the volume has to be updated (e.g. volume or resolution changed)

    bool geometryInvalid_;   ///< determines if the mesh geometry has to be updated, e.g. because of a transfer function change

    /// Proxy geometry generated by this processor.
    TriangleMeshGeometryColorNormal* geometry_;

    TriangleMeshGeometryColorNormal* tmpGeometry_; ///< temporary geometry that is used while background thread computations are not finished yet

    /// data structure for spatial subdivision
    std::vector<ProxyGeometryVolumeRegion> volumeStructure_;

    /// size of VolumeRegion structure (ie. number of bricks in every direction)
    tgt::ivec3 volStructureSize_;

    const VolumeBase* currentVolume_;

    /// Category used for logging.
    static const std::string loggerCat_;

    ProcessorBackgroundThread<OptimizedProxyGeometry>* backgroundThread_; ///< background thread for proxy geometry computations

    std::vector<TransFunc1DKeys*> tfCopies_;      ///< contains copies of the TFs for background thread computations

}; // OptimizedProxyGeometry


//-------------------------------------------------------------------------------------------------
// Background threads performing the actual construction of the proxy geometry

/**
 * Abstract class derived from ProcessorBackgroundThread used as base class for proxy geometry background computations.
 */
class OptimizedProxyGeometryBackgroundThread : public ProcessorBackgroundThread<OptimizedProxyGeometry> {

public:
    OptimizedProxyGeometryBackgroundThread(OptimizedProxyGeometry* processor, const VolumeBase* volume,
        std::vector<TransFunc1D*> tfVector, float threshold, TriangleMeshGeometryColorNormal* geometry, int stepSize,
        bool debugOutput, bool clippingEnabled, tgt::vec3 clipLlf = tgt::vec3(0.f), tgt::vec3 clipUrb = tgt::vec3(1.f));

protected:
    virtual void handleInterruption();

    /**
     * This heuristically checks if a region of the volume can be treated as empty by comparing the opacity of the
     * pre-integration table (called with the min and max intensity of the region) to threshold_*1e-3
     */
    bool isRegionEmptyPi(float min, float max, const PreIntegrationTable* piTable) const;

    const VolumeBase* volume_;      ///< the volume for which the proxy geometry should be computed

    std::vector<TransFunc1D*> tfCopyVector_;        ///< TFs to identify empty (ie. completely transparent) regions

    float threshold_;               ///< threshold for the pre-integration heuristic (used by method isRegionEmptyPi)

    TriangleMeshGeometryColorNormal* geometry_;    ///< geometry of the processor

    int stepSize_;                  ///< step size (in voxels) for subdividing the volume into bricks

    bool clippingEnabled_;          ///< use axis-aligned clipping?

    tgt::vec3 clipLlf_;             ///< clipping region boundaries (in voxel space)
    tgt::vec3 clipUrb_;             ///< clipping region boundaries (in voxel space)

    bool debugOutput_;              ///< print out debug information?
    tgt::Stopwatch stopWatch_;      ///< for determining the building time of the data structures
};

//-------------------------------------------------------------------------------------------------

/**
 * Abstract class derived from ProxyGeometryBackgroundThread, used as base class for proxy geometry
 * background computations that use the volume region structure.
 */
class StructureProxyGeometryBackgroundThread : public OptimizedProxyGeometryBackgroundThread {

public:
    StructureProxyGeometryBackgroundThread(OptimizedProxyGeometry* processor, const VolumeBase* volume,
        std::vector<TransFunc1D*> tfVector, float threshold, TriangleMeshGeometryColorNormal* geometry, std::vector<ProxyGeometryVolumeRegion>* volumeStructure,
        tgt::ivec3 volStructureSize, int stepSize, bool debugOutput, bool clippingEnabled,
        tgt::vec3 clipLlf = tgt::vec3(0.f), tgt::vec3 clipUrb = tgt::vec3(1.f));

protected:
    void computeRegionStructure();               ///< determine the volume region data structurei

    ///get the volume region at a given position within the structure
    ProxyGeometryVolumeRegion& getVolumeRegion(tgt::ivec3 pos);

    std::vector<ProxyGeometryVolumeRegion>* volumeStructure_; ///< data structure for spatial subdivision

    tgt::ivec3 volStructureSize_;                ///< size of the spatial subdivision (ie. number of regions in every direction)
};

//-------------------------------------------------------------------------------------------------

/**
 * Background thread computing an approximate minimal bounding box as a proxy geometry
 */
class MinCubeBackgroundThread : public StructureProxyGeometryBackgroundThread {

public:
    MinCubeBackgroundThread(OptimizedProxyGeometry* processor, const VolumeBase* volume,
        std::vector<TransFunc1D*> tfVector, float threshold, TriangleMeshGeometryColorNormal* geometry, std::vector<ProxyGeometryVolumeRegion>* volumeStructure,
        tgt::ivec3 volStructureSize, int stepSize, bool debugOutput, bool clippingEnabled,
        tgt::vec3 clipLlf = tgt::vec3(0.f), tgt::vec3 clipUrb = tgt::vec3(1.f));

    ~MinCubeBackgroundThread();

protected:
    void threadMain();
    void computeMinCube();

};

//-------------------------------------------------------------------------------------------------

/**
 * Background thread computing a proxy geometry of (maximal) bricks using the volume region structure.
 */
class VisibleBricksBackgroundThread : public StructureProxyGeometryBackgroundThread {

public:

    VisibleBricksBackgroundThread(OptimizedProxyGeometry* processor, const VolumeBase* volume, std::vector<TransFunc1D*> tfVector,
        float threshold, TriangleMeshGeometryColorNormal* geometry, std::vector<ProxyGeometryVolumeRegion>* volumeStructure,
        tgt::ivec3 volStructureSize, int stepSize, bool debugOutput, bool clippingEnabled,
        tgt::vec3 clipLlf = tgt::vec3(0.f), tgt::vec3 clipUrb = tgt::vec3(1.f));

    ~VisibleBricksBackgroundThread();

protected:

    void threadMain();

    void computeMaximalBricks();

    ///set volume regions attribute bound_
    void setVolBound(tgt::ivec3 llf, tgt::ivec3 urb, bool value = true);

    ///check if one volume region is not empty (according to the pre-integration heuristic) and not bound
    bool isVolNotEmptyPiNotBound(tgt::ivec3 pos, const std::vector<TransFunc1D*>& tfs,  const std::vector<const PreIntegrationTable*>& piTables);

    ///check if the volume regions are not empty (according to the pre-integration heuristic) and not bound
    bool isVolNotEmptyPiNotBound(tgt::ivec3 llf, tgt::ivec3 urb, const std::vector<TransFunc1D*>& tfs,  const std::vector<const PreIntegrationTable*>& piTables);

    /// get the urb of a volume region by adding its non-empty (determinded by isRegionEmptyPi method)
    /// and not already bound neighbor regions, constructing a maximal cube
    tgt::ivec3 getUrbPi(tgt::ivec3 llf, const std::vector<TransFunc1D*>& tfs,  const std::vector<const PreIntegrationTable*>& piTables);

};

//-------------------------------------------------------------------------------------------------

/**
 * Background thread computing a proxy geometry of the outer faces of opaque neighbors in the volume region structure
 */
class OuterFacesBackgroundThread : public StructureProxyGeometryBackgroundThread {

public:

    OuterFacesBackgroundThread(OptimizedProxyGeometry* processor, const VolumeBase* volume, std::vector<TransFunc1D*> tfVector,
        float threshold, TriangleMeshGeometryColorNormal* geometry, std::vector<ProxyGeometryVolumeRegion>* volumeStructure,
        tgt::ivec3 volStructureSize, int stepSize, bool debugOutput, bool clippingEnabled,
        tgt::vec3 clipLlf = tgt::vec3(0.f), tgt::vec3 clipUrb = tgt::vec3(1.f));

    ~OuterFacesBackgroundThread();

protected:

    void threadMain();

    void computeOuterFaces();

};
//-------------------------------------------------------------------------------------------------

/**
 * Background thread computing a proxy geomtry using the VolumeOctree volume representation (if available).
 * If this representation does not exist, the visible bricks mode is used as a fallback mode.
 */
class VolumeOctreeBackgroundThread : public OptimizedProxyGeometryBackgroundThread {

public:

    VolumeOctreeBackgroundThread(OptimizedProxyGeometry* processor, const VolumeBase* volume, std::vector<TransFunc1D*> tfVector,
        float threshold, TriangleMeshGeometryColorNormal* geometry, int stepSize, bool clippingEnabled,
        tgt::vec3 clipLlf = tgt::vec3(0.f), tgt::vec3 clipUrb = tgt::vec3(1.f));

    ~VolumeOctreeBackgroundThread();

protected:

    /// compute an octree data structure and traverse it to create the proxy geometry
    void threadMain();

    /// uses an octree to create a proxygeometry of maximal cubes
    void computeVolumeOctreeGeometry();

    /// recursively traverses the octree and creates the geometry
    void traverseOctreeAndCreateGeometry(const VolumeOctreeNode* node, const tgt::vec3& nodeLlf, const tgt::vec3& nodeUrb,
            /*TransFunc1DKeys* tf, const PreIntegrationTable* piTable,*/
            std::vector<TransFunc1D*> tfVector, int numVolumeChannels, std::vector<const PreIntegrationTable*> piTables,
            tgt::Bounds clipBounds, tgt::ivec3 volumeDim);

};

/**
 * Background thread computing an outer faces proxy geomtry using the VolumeOctree volume representation (if available).
 * If this representation does not exist, the outer faces mode is used as a fallback mode.
 */
class VolumeOctreeOuterFacesBackgroundThread : public OptimizedProxyGeometryBackgroundThread {

public:

    VolumeOctreeOuterFacesBackgroundThread(OptimizedProxyGeometry* processor, const VolumeBase* volume, std::vector<TransFunc1D*> tfVector,
        float threshold, TriangleMeshGeometryColorNormal* geometry, int stepSize, bool clippingEnabled,
        tgt::vec3 clipLlf = tgt::vec3(0.f), tgt::vec3 clipUrb = tgt::vec3(1.f));

    ~VolumeOctreeOuterFacesBackgroundThread();

protected:

    /// compute an octree data structure and traverse it to create the proxy geometry
    void threadMain();

    /// uses an octree to create a proxygeometry of maximal cubes
    void computeVolumeOctreeGeometry();

    /// recursively traverses the octree and creates the geometry
    void traverseOctreeAndCreateBrickStructure(const VolumeOctreeNode* node, size_t level, const tgt::svec3& nodeLlf, const tgt::svec3& nodeUrb,
            std::vector<TransFunc1D*> tfVector, int numVolumeChannels, std::vector<const PreIntegrationTable*> piTables,
            tgt::Bounds clipBounds, tgt::ivec3 volumeDim);

    /// mark several bricks as visible
    void setBricksVisible(tgt::svec3 llf, tgt::svec3 urb);

    ProxyGeometryVolumeRegion& getVolumeRegion(tgt::ivec3 pos);

    void computeOuterFaces();

private:

    tgt::svec3 brickDim_;
    tgt::svec3 volDim_;
    tgt::svec3 brickStructureSize_;
    size_t octreeDepth_;      ///< actual depth of the optimized octree for determining if an octree leaf contains several bricks
    std::vector<ProxyGeometryVolumeRegion> brickStructure_; ///< data structure for spatial subdivision

};

} // namespace

#endif // VRN_OPTIMIZEDPROXYGEOMETRY_H
