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

#ifndef VRN_FLOWSLICERENDERER_H
#define VRN_FLOWSLICERENDERER_H

#include "voreen/core/processors/renderprocessor.h"
#include "voreen/core/datastructures/volume/volume.h"
#include "modules/flowreen/utils/colorcodingability.h"
#include "flowreenprocessor.h"

#include "voreen/core/ports/renderport.h"
#include "voreen/core/ports/volumeport.h"

#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/optionproperty.h"

#include "voreen/core/datastructures/volume/volumegl.h"

namespace tgt { class Texture; }

namespace voreen {

class VolumeFlow3D;
class Flow2D;
class Flow3D;
template<typename T> class SimpleTexture;

/**
 * (Abstract) base class for all FlowSliceRenderer processors providing the
 * basic texture generation techniques.
 */
class FlowSliceRenderer : public RenderProcessor, public FlowreenProcessor {
public:
    FlowSliceRenderer();
    virtual ~FlowSliceRenderer();

    virtual void initialize();
    virtual void deinitialize();

protected:
    enum SliceAlignment {
        PLANE_XY = 2,   /** view from the volume's front to back (negative z-axis) */
        PLANE_XZ = 1,   /** view from the volume's bottom to top (negative y-axis) */
        PLANE_ZY = 0    /** view from the volume's right to left (negative x-axis) */
    };
    friend class OptionProperty<SliceAlignment>;

    enum RenderingTechnique {
        TECHNIQUE_COLOR_CODING,
        TECHNIQUE_INTEGRATE_DRAW,
        TECHNIQUE_FAST_LIC,
        TECHNIQUE_ARROW_PLOT_RAND,
        TECHNIQUE_ARROW_PLOT_GRID,
        TECHNIQUE_COLOR_CODING_PROJECTED,
        TECHNIQUE_INTEGRATE_DRAW_PROJECTED,
        TECHNIQUE_FAST_LIC_PROJECTED,
        TECHNIQUE_SPOTNOISE_PROJECTED,
        TECHNIQUE_SPOTNOISE
    };
    friend class OptionProperty<RenderingTechnique>;

    enum TextureZoom {
        ZOOM_FULL_SIZE = 0,
        ZOOM_1X = 1,
        ZOOM_2X = 2,
        ZOOM_4X = 4,
        ZOOM_8X = 8
    };
    friend class OptionProperty<TextureZoom>;

    typedef unsigned char BYTE;

    static tgt::ivec3 getCoordinatePermutation(const SliceAlignment& alignment);

protected:
    /**
     * Prepares the OpenGL machine for orthographic rendering by pushing
     * all attributes, setting up the ModelView and the Projection matrices,
     * setting the viewport and disabeling the depth test.
     * You MUST call <code>endOrthographicRendering()</code> after EVERY
     * call to <code>beginOrthographicRendering</code>.
     * This is basically used by any texture generation technique.
     *
     * @param   viewport    dimensions of the viewport (usually the size
     *                      of the texture to be generated)
     */
    void beginOrthographicRendering(const tgt::ivec2& viewport) const;

    /**
     * Restores all the settings which have been modified by calling
     * <code>beginOrthographicRendering()</code>. This concerns the
     * all attributes, the depth test, the ModelView and the Projecteion
     * matrix.
     */
    void endOrthographicRendering() const;

    /**
     * Prepares the texture containing 2D flow data for usage with OpenGL.
     * The texture is stored in flow2DTexture_ and a previously exisiting
     * one will be deleted and its memory will be freed.
     */
    void prepareFlowTexture2D(const Flow2D& flow2D);

    /**
     * Prepares the texture containing 3D flow data for usage with OpenGL.
     * The texture is stored in flow3DTexture_. If an error occurs (i.e.
     * if the supplied Volume is invalid), flow3DTexture_ will be
     * set to NULL.
     */
    void prepareFlowTexture3D();

    /**
     * Generates a tgt::Texture from the specified render target with the given texture size.
     * The texture contains 4 color channels (GL_RGBA) stored as floats (GL_FLOAT).
     * The caller has to free the texture using delete, afterwards.
     *
     * @param   rendeport   RenderPort object containing the image.
     * @param   textureSize Size of the texture. Note that this size must be less or equal to
     *                      the actual size of the render target!
     */
    tgt::Texture* copyTargetToTexture(RenderPort& renderport, const tgt::vec2& textureSize) const;

    /**
     * Generates a texture containing an arrow plot from the given flow volume data.
     *
     * @param   flow    The flow data from which the texture shall be generated.
     * @param   textureSize Target size of the output texture.
     * @param   tempPort    Renderport used for texture generation. This parameter
     *                      must be supplied by the caller and should be a private
     *                      port to avoid resource conflicts during rendering.
     * @param   thresholds  Lower and upper threshold for the flow magnitudes which
     *                      shall be rendered. The thresholds are absolute values
     *                      of the magnitudes from the input flow texture. Magnitudes
     *                      not within the range of the thresholds will be dropped
     *                      during texture generation.
     */
    tgt::Texture* renderArrowPlotTexture(const Flow2D& flow, const tgt::vec2& textureSize,
        RenderPort& tempPort, const tgt::vec2& thresholds);

    /**
     * Generates are color plot of a user-defined scalar quantity of the input flow.
     * The scalar quantity is determined by a property and the output texture can
     * be either generated by operating on a derived 2D texture by ignoring the
     * remaining component, or by computing values from the 3D texture and projecting
     * the result afterwards onto a 2D texture.
     *
     * @param   flow    The flow data from which the texture shall be generated.
     * @param   sliceNo The number of the slice for which the texture shall be
     *                  generated.
     * @param   textureSize Target size of the output texture.
     * @param   tempPort    Renderport used for texture generation. This parameter
     *                      must be supplied by the caller and should be a private
     *                      port to avoid resource conflicts during rendering.
     * @param   thresholds  Lower and upper threshold for the flow magnitudes which
     *                      shall be rendered. The thresholds are absolute values
     *                      of the magnitudes from the input flow texture. Magnitudes
     *                      not within the range of the thresholds will be dropped
     *                      during texture generation.
     * @param   projected   Determines whether the calculations are performed on
     *                      the 3D texture and are projected afterwards (in case
     *                      the value is true), or if a 2D texture is derived first
     *                      and calcuations are done on the 2D texture at once (case
     *                      of the value being false).
     */
    tgt::Texture* renderColorCodingTexture(const Flow3D& flow, const size_t sliceNo,
        const tgt::vec2& textureSize, RenderPort& tempPort, const tgt::vec2& threshold,
        const bool projected = false);

    /**
     * Renders a 2D texture containing flow images derived from the input flow volume
     * and by using the given technique. According to the technique, this methods
     * deligates to the other rendering methods.
     *
     * @param   flow    The flow data from which the texture shall be generated.
     * @param   sliceNo The number of the slice for which the texture shall be
     *                  generated.
     * @param   textureScaling  Determines the magnification of the output texture
     *                          according to the input flow texture respectively to
     *                          size of the slice taken from it. A value of 2 for
     *                          instance will generate a texture of double input size.
     *                          This parameter may not be negative.
     * @param   technique   Determines the technique to be used for texutre generation.
     * @param   thresholds  Lower and upper threshold for the flow magnitudes which
     *                      shall be rendered. The thresholds are absolute values
     *                      of the magnitudes from the input flow texture. Magnitudes
     *                      not within the range of the thresholds will be dropped
     *                      during texture generation.
     */
    tgt::Texture* renderFlowTexture(const Flow3D& flow, const size_t sliceNo,
        const int textureScaling, const RenderingTechnique& technique, const tgt::vec2& thresholds);

    /**
     * Generates a spot noise texture for flow visualization according the (enhaced)
     * spot noise technique by van Wijk and deLeeuw. This methods does not apply
     * flow magnitude thresholds.
     *
     * @param   flow    The flow data from which the texture shall be generated.
     * @param   sliceNo The number of the slice for which the texture shall be
     *                  generated.
     * @param   textureSize Target size of the output texture.
     * @param   viewportSize    The current size of the viewport. This is required
     *                      for internal calculations, e.g. the spot size.
     * @param   tempPorts   Renderports used for texture generation. This parameter
     *                      must be supplied by the caller and should be a vector of
     *                      private rendering ports to avoid resource conflicts during
     *                      rendering.
     * @param   projected   Determines whether the calculations are performed on
     *                      the 3D texture and are projected afterwards (in case
     *                      the value is true), or if a 2D texture is derived first
     *                      and calcuations are done on the 2D texture at once (case
     *                      of the value being false).
     */
    tgt::Texture* renderSpotNoiseTexture(const Flow3D& flow, const size_t sliceNo,
        const tgt::vec2& textureSize, const tgt::vec2& viewportSize,
        const std::vector<RenderPort*>& tempPorts, const bool projected = false);

    /**
     * Initializes the texture stored in noiseTexture_ with white noise and
     * frees a previously stored one. Used by FastLIC.
     */
    void initNoiseTexture(const tgt::ivec2& size);

    /**
     * Initializes the given number of random positions and random intensities.
     * Previously stored ones are freed and replaced. Used by Spot Noise.
     */
    void initRandomValues(const size_t count);

    /**
     * Utility method to facilitate re-seeding a random position which leads to
     * bad streamlines.
     *
     * @param   validPositions  the number of positions within the array randomPositions_
     *                          which are good (used for reseeding if > 0).
     *
     * @return  a new random position.
     */
    tgt::vec2 reseedRandomPosition(const size_t validPositions = 0);

    virtual std::string generateShaderHeader();

    bool loadShader(const std::string& vertexShaderName, const std::string& fragmentShaderName);

    bool rebuildShader();

private:
    void buildArrowDisplayList(const Flow2D& flow, const tgt::vec2& textureSize,
        const tgt::vec2& thresholds);

    void buildArrowGridDisplayList(const Flow2D& flow, const tgt::vec2& textureSize,
        const tgt::vec2& thresholds);

    /**
     * Performs the generation of an image of the specified slice through
     * the flow volume by using FastLIC (line integral convolution).
     * The returned pointer is stored in a tgt::Texture and deleted
     * by that dtor.
     */
    float* fastLIC(const Flow3D& flow, const size_t sliceNo, const int textureScaling);

    /**
     * Calculates and returns the initial intensity by evaluating a convolution
     * integral for projected FastLIC.
     *
     * @param   indexR0 index of the location on the streamline for which the
     *                  intensity has to be calculated
     * @param   kernelSize  size of the convolution kernel
     * @param   streamline  the streamline along which the convolution is
     *                      performed
     *
     * @return  the initial intensity for the FastLIC convolution
     */
    float initialFastLICIntensity(const int indexR0, const size_t kernelSize,
        const std::vector<tgt::ivec2>& streamline) const;

    /**
     * Performs the generation of an image of the specified slice
     * through the flow volume by using Integrate & Draw.
     * The returned pointer is stored in a tgt::Texture and deleted
     * by that dtor.
     */
    float* integrateDraw(const Flow3D& flow, const size_t sliceNo, const int textureScaling);

    /**
     * Prepares the texture to be re-created by indicating that it has become
     * invalid.
     */
    void invalidateTexture();

    void onColorCodingChange();

    /**
     * Action to be performed when the number of seeding positions changes
     * due to interaction with the corresponding Property (i.e.
     * re-initialization of random values)
     */
    void onNumberSeedingPositionsChange();

    /**
     * Action to be performed when the rendering techniques has been changed
     * due to Property changes.
     */
    void onTechniqueChange();

    /**
     * Prepares the shader in shader_ for rendering a color coded 2D flow slice.
     * The shader will be loaded, activated and the uniforms will be set.
     * The caller is responsible to deactivate that shader.
     */
    void prepareColorCoding2DShader(const Flow2D& flow2D, const tgt::vec2& thresholds);

    /**
     * Prepares the shader in shader_ for rendering a color coded slice using
     * a 3D input texture. This will result in a kind of "projection".
     * The shader will be loaded, activated and the uniforms will be set.
     * The caller is responsible to deactivate that shader.
     */
    void prepareColorCoding3DShader(const Flow3D& flow3D, const size_t sliceNo,
        const tgt::vec2& thresholds);

    /**
     * Renders a 2D arrow in XY-plane at [] with bounding dimensions consiting of
     * a quad and a triangle. This used by <code>buildArrowDisplayList()</code>.
     *
     * @param   transformation  Transformation matrix used to rotate, scale and
     *                          translate the arrow
     */
    void renderArrow(const tgt::mat3& transformation) const;

    /**
     * Renders a 2D quad at z = 0.0 with the given tessalation. The x- and y-values
     * of the vertices are within [0.0, 1.0].
     * The produced mesh is used by Spot Noise for application of the spot texture which
     * is generated by the corresponding shader. Used by <code>renderSpotNoiseTexture()</code>
     */
    void renderSpotMesh(const size_t tesselation = 1);

    /**
     * Toggles the visibility of the properties of this processor. The visibility
     * depends on the current rendering technique.
     */
    void toggleProperties();

protected:
    OptionProperty<RenderingTechnique>* techniqueProp_;
    IntProperty numSeedingPositionsProp_;
    IntProperty spotTesselationProp_;
    IntProperty spotScaleProp_;
    BoolProperty useSpotFilteringProp_; /** Determines whether to use filtring for spots when using Spot Noise */
    BoolProperty useSpotBendingProp_;
    BoolProperty useAdaptiveKernelSizeProp_;
    IntProperty kernelSizeProp_;
    IntProperty pixelSamplingProp_;
    IntProperty arrowSizeProp_;
    IntProperty arrowSpacingProp_;

    ColorCodingAbility colorCoding_;
    TextureZoom zoom_;

    tgt::ivec3 permutation_;        /** permutation of the components of the flow according to selected slice alignment */
    tgt::Shader* shader_;           /** all-purpose shader, used for several tasks */
    float* randomIntensities_;      /** array of random intensities used by Spot Noise */
    tgt::vec2* randomPositions_;    /** array of random positions for spots used by Spot Noise */
    SimpleTexture<unsigned char>* noiseTexture_;   /** "texture" containing white noise used by FastLIC */
    tgt::Texture* flow2DTexture_;
    const VolumeTexture* flow3DTexture_;
    bool rebuildTexture_;           /** indicates whether the texture containing the slice image needs to be rebuilt. */
    GLuint arrowList_;

    VolumePort volInport_;
    RenderPort imgOutport_;
    RenderPort privatePort1_;
    RenderPort privatePort2_;
    PortGroup portGroup_;
};

}   // namespace

#endif  // VRN_FLOWSLICERENDERER_H
