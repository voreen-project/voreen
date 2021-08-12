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

#ifndef VRN_SLICEVIEWER_H
#define VRN_SLICEVIEWER_H

#include "voreen/core/processors/volumerenderer.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/matrixproperty.h"
#include "voreen/core/interaction/mwheelnumpropinteractionhandler.h"
#include "voreen/core/datastructures/volume/slice/slicecache.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/utils/voreenqualitymode.h"

namespace voreen {

template<typename P>
class OctreeSliceTexture {
public:
    typedef P Pixel;
    OctreeSliceTexture(GLint format, GLint datatype, tgt::Texture::Filter filter);
    void updateDimensions(tgt::svec2);
    void uploadTexture();
    void bindTexture();
    Pixel* buf();
    tgt::ivec2 dimensions();
    void clear();
private:
    std::vector<P> buf_;
    tgt::Texture texture_;
};

typedef OctreeSliceTexture<tgt::Vector4<uint16_t>> OctreeSliceTextureColor;
typedef OctreeSliceTexture<uint8_t> OctreeSliceTextureControl;

template<typename Pixel>
OctreeSliceTexture<Pixel>::OctreeSliceTexture(GLint format, GLint datatype, tgt::Texture::Filter filter)
    : buf_()
    , texture_(tgt::svec3(2,2,1), format, format, datatype, filter, tgt::Texture::CLAMP_TO_EDGE, nullptr, false)
{
}

template<typename Pixel>
void OctreeSliceTexture<Pixel>::updateDimensions(tgt::svec2 dim) {
    if(dim != tgt::svec2(texture_.getDimensions().xy())) {
        texture_.setCpuTextureData(nullptr, false);
        texture_.updateDimensions(tgt::ivec3(dim.x, dim.y, 1), false);

        buf_.resize(tgt::hmul(dim), Pixel(0));
        texture_.setCpuTextureData(reinterpret_cast<GLubyte*>(buf_.data()), false);
        LGL_ERROR;
    }
}

template<typename Pixel>
void OctreeSliceTexture<Pixel>::uploadTexture() {
    texture_.uploadTexture();
    LGL_ERROR;
}

template<typename Pixel>
void OctreeSliceTexture<Pixel>::bindTexture() {
    texture_.bind();
    LGL_ERROR;
}
template<typename Pixel>
Pixel* OctreeSliceTexture<Pixel>::buf() {
    return buf_.data();
}
template<typename Pixel>
tgt::ivec2 OctreeSliceTexture<Pixel>::dimensions() {
    return texture_.getDimensions().xy();
}
template<typename Pixel>
void OctreeSliceTexture<Pixel>::clear() {
    std::fill(buf_.begin(), buf_.end(), Pixel(0));
    uploadTexture();
}

struct OctreeSliceViewProgress {
    OctreeSliceViewProgress();
    size_t nextSlice_;
    size_t nextTileX_;
    size_t nextTileY_;
};

/**
 * Performs slice rendering of a single or multiple slices
 * along one of the three main axis of the volume.
 */
class VRN_CORE_API SliceViewer : public VolumeRenderer, public QualityModeObserver {

public:
    SliceViewer();
    virtual ~SliceViewer();
    virtual Processor* create() const;

    virtual std::string getClassName() const    { return "SliceViewer";      }
    virtual std::string getCategory() const     { return "Slice Rendering";  }
    virtual CodeState getCodeState() const      { return CODE_STATE_STABLE;  }

    struct VertexHelper {
        VertexHelper(tgt::vec2 pos, tgt::vec4 texCoord) : position_(pos), texCoord_(texCoord) {}
        tgt::vec2 position_;
        tgt::vec4 texCoord_;
    };

protected:

    virtual void initialize();
    virtual void deinitialize();

    virtual void beforeProcess();
    virtual void process();
    virtual void afterProcess();

    virtual void qualityModeChanged();

    virtual void adjustPropertiesToInput();

    void invalidateOctreeTexture();

    void renderFromVolumeTexture(tgt::mat4 toSliceCoordMatrix);
    void renderFromOctree();

    /// Generates the header for the shader depending on the choice of features to be used.
    std::string generateSliceShaderHeader(const tgt::GpuCapabilities::GlVersion* version = 0);

    /// Recompiles the shader.
    bool rebuildShader();

    /**
     * Adapts the min/max ranges of the respective properties to the
     * dimensions of the currently connected volume and adjusts the property visibility.
     */
    void updatePropertyConfiguration();

    void resetChannelShift();

    /**
     * Adapts the volume coord permutation and the states
     * of the slice properties when the slice alignment has changed.
     */
    void onSliceAlignmentChange();

    /**
     * Renders a GL_LINE_LOOP with white color and the current slice's
     * dimensions in order to visualize its boundary.
     */
    void renderSliceBoundaries() const;

    /**
     * Displays the cursor position in volume coordinates, the data value under the cursor,
     * and/or the current slice number. What data is actually shown, is determined by the
     * respective properties.
     *
     * @note requires freetype (module 'fontrendering')
     */
    void renderInfoTexts() const;

    /**
     * Draws the scale legend in the lower left corner.
     */
    void renderLegend();

    /**
     * Tiny helper function returning true when numSlicesPerRow_ and
     * numSlicesPerCol_ are set to one.
     */
    bool singleSliceMode() const;

    /**
     * Converts the given screen coordinates into the corresponding voxel
     * coordinates for the slice.
     */
    tgt::vec3 screenToVoxelPos(tgt::ivec2 screenPos) const;

    /**
     * Generates a matrix that maps from viewport coordinates
     * to volume coordinates.
     *
     * @note This matrix is only valid in single slice mode,
     *  since for multi slice mode, the mapping from viewport
     *  to volume coordinates cannot be expressed by a linear function.
     */
    tgt::mat4 generatePickingMatrix() const;

    /**
     * Is called by the respective event properties on mouse move or
     * mouse click events and saves the current cursor position.
     */
    void mouseLocalization(tgt::MouseEvent* e);

    /**
     * Is called by the mouseEventShift_ property when a slice shift
     * event has occurred.
     */
    void shiftEvent(tgt::MouseEvent* e);


    /**
     * Resets the zoomFactor and voxelOffset
     */
    void resetView();

protected:
    enum TextureMode {
        TEXTURE_2D,
        TEXTURE_3D,
        OCTREE,
    };

    enum InfoAlignment {
        ALIGNMENT_N,
        ALIGNMENT_NE,
        ALIGNMENT_E,
        ALIGNMENT_SE,
        ALIGNMENT_S,
        ALIGNMENT_SW,
        ALIGNMENT_W,
        ALIGNMENT_NW,
    };
    friend class OptionProperty<InfoAlignment>;

    virtual void setDescriptions() {
        setDescription("Displays 2D slices along one of the three main axis of the volume. Multiple slices can be viewed simultaneously.");
    }

    VolumePort inport_;
    RenderPort outport_;

    TransFunc1DKeysProperty transferFunc1_;
    TransFunc1DKeysProperty transferFunc2_;
    TransFunc1DKeysProperty transferFunc3_;
    TransFunc1DKeysProperty transferFunc4_;

    /// Property containing the available alignments: xy (axial), xz (coronal), yz (sagittal)
    OptionProperty<SliceAlignment> sliceAlignment_;
    IntProperty sliceIndex_;                ///< Property containing the currently selected slice
    IntProperty numGridRows_;               ///< Property containing the current row count of the displayed grid
    IntProperty numGridCols_;               ///< Property containing the current column count of the displayed grid
    BoolProperty selectCenterSliceOnInputChange_;   ///< if true, the center slice index is selected when the input volume changes

    IntVec3Property mouseCoord_;            ///< Property containing the current slice position of the mouse cursor

    BoolProperty renderSliceBoundaries_;    ///< Determines whether to render the slice boundaries
    ColorProperty boundaryColor_;           ///< Color to be used for rendering of the slice boundary
    IntProperty boundaryWidth_;             ///< line width of the boundary

    StringOptionProperty showCursorInfos_;  ///< Determines whether information about the cursor position is to be shown
    OptionProperty<InfoAlignment> infoAlignment_;   ///< alignment of the legend
    BoolProperty showSliceNumber_;          ///< Determines whether the slice number is to displayed on each slice
    IntProperty fontSize_;                  ///< Font size to be used for info texts
    BoolProperty showScaleLegend_;          ///< true, if the distance legend should be rendered
    FloatProperty legendLineLength_;        ///< length of the distance legend in pixels

    FloatVec2Property voxelOffset_;         ///< The 2D voxel position at which the view area is shifted
    FloatProperty zoomFactor_;              ///< Specifies the current slice zoom factor: the standard value of 1.0
                                            ///  causes the whole slice to be displayed
    ButtonProperty resetViewButton_;        ///< Button to reset zoomFactor_ and voxelOffset_ to it's defaults


    FloatMat4Property pickingMatrix_;       ///< Read-only (generated) matrix mapping from viewport to
                                            ///  volume coordinates.

    OptionProperty<TextureMode> texMode_;   ///< use 2D slice textures or 3D volume texture?
    IntProperty sliceLevelOfDetail_;        ///< level of detail used for 2D slice extraction
    IntProperty interactionLevelOfDetail_;  ///< level of detail during user interaction
    IntProperty sliceExtractionTimeLimit_;  ///< timelimit in milliseconds for 2D slice extraction
    IntProperty sliceCacheSize_;            ///< size of the 2D slice cache (0 means no cache).

    BoolProperty applyChannelShift_;
    FloatVec3Property channelShift1_;
    FloatVec3Property channelShift2_;
    FloatVec3Property channelShift3_;
    FloatVec3Property channelShift4_;

    ButtonProperty resetChannelShift_;

    EventProperty<SliceViewer>* mouseEventPress_;
    EventProperty<SliceViewer>* mouseEventMove_;
    EventProperty<SliceViewer>* mouseEventShift_;

    MWheelNumPropInteractionHandler<int> mwheelCycleHandler_;
    MWheelNumPropInteractionHandler<float> mwheelZoomHandler_;

    tgt::Shader* sliceShader_;
    tgt::Shader* octreeSliceTextureShader_;

    SliceCache sliceCache_;       ///< Cache for slices created in 2D texture mode.

    tgt::ivec3 voxelPosPermutation_;    ///< permutation of voxel coordinates according to current alignment
    tgt::vec2 sliceLowerLeft_;          ///< lower left of the whole slice plate in viewport coordinates
    tgt::vec2 sliceSize_;               ///< size of a single slice in viewport coordinates
    tgt::mat4 textureMatrix_;           ///< matrix that is currently applied to texture coordinates
    bool sliceComplete_;                ///< is set in process() and specifies whether the current slice has been created in full LOD,
                                        ///< if not, an invalidation is triggered in afterProcess().

    static const std::string fontName_; ///< path and name of the font used for text-rendering

    static const std::string loggerCat_;

private:
    void renderSliceGeometry(const tgt::vec4& t0, const tgt::vec4& t1, const tgt::vec4& t2, const tgt::vec4& t3) const;

    mutable tgt::ivec2 mousePosition_;          ///< Current mouse position
    mutable bool mouseIsPressed_;               ///< Is a mouse button currently pressed?
    mutable tgt::ivec3 lastPickingPosition_;    ///< Mouse position during previous interaction

    std::unique_ptr<OctreeSliceTextureColor> octreeTextureInteractive_;
    std::unique_ptr<OctreeSliceTextureControl> octreeTextureInteractiveControl_;
    std::unique_ptr<OctreeSliceTextureColor> octreeTexture_;
    std::unique_ptr<OctreeSliceTextureControl> octreeTextureControl_;
    OctreeSliceViewProgress octreeRenderProgress_;
};

} // namespace

#endif // VRN_SLICEVIEWER_H
