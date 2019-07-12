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

#ifndef VRN_ISOSURFACEEXTRACTOR_H
#define VRN_ISOSURFACEEXTRACTOR_H

#include "voreen/core/processors/volumeprocessor.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/progressproperty.h"


namespace voreen {

class IsosurfaceExtractor : public VolumeProcessor {
public:
    IsosurfaceExtractor();
    virtual ~IsosurfaceExtractor();

    virtual std::string getClassName() const         { return "IsosurfaceExtractor";      }
    virtual std::string getCategory() const       { return "Volume Processing"; }
    virtual VoreenSerializableObject* create() const;
    virtual void setDescriptions() {
        setDescription("Processor that extracts an isosurface mesh from a volume using the Marching Cubes algorithm "
                "(see Lorensen, W. E. and Cline, H. E. \"Marching cubes: A high resolution 3D surface construction algorithm (1987)\")");
        enabledProp_.setDescription("Enable or disable isosurface extraction.");
        surfaceTypeProp_.setDescription("Choose the type of surface to be extracted from the volume. "
                "Marching Cubes will extract a smooth approximation of an isosurface. "
                "Blocks will extract the Voronoi cell walls within the volume dimensions for all voxel neighbor "
                "pairs whose intensities enclose the isovalue. "
                "The outside of the volume is assumed to have an associated intensity smaller than the isovalue.");
        marchingCubesProcessingModeProp_.setDescription(
                "Choose between extracting the isosurface on the CPU or GPU. "
                "GPU computation should generally be faster."
                "Only available for the Marching Cubes surface type."
                );
        isoValueProp_.setDescription("The isovalue for the generated isosurface. (Note: in real world units)");
        meshBufferSizeProp_.setDescription(
                "Reserved size for the buffer that holds the isosurface mesh. "
                "The maximum value is large enough to always hold the entire isosurface. "
                "However, in most cases a much smaller size is sufficient, so that a user can "
                "by hand reduce the size if GPU memory is scarce."
                "If the buffer is to small, the mesh will be missing primitives, but a warning will be issued."
                );

        useSmoothingProp_.setDescription(
                "Enable or disable smoothing prior to extracting an isosurface. This is intended to be used on "
                "binary volumes. Note that the isovalue only accurately represents the volume for a normalized value of 0.5."
                "This is an implementation of a slightly modified version of the algorithm presented in "
                "\"Surface Extraction from Binary Volumes with Higher-Order Smoothness\" by Lempitsky, V. (2009)."
                );
        maxGradientProp_.setDescription(
                "Exponent (10<sup>x</sup>) for the maximum length of the gradient during smoothing of the input volume. "
                "Lower/Higher values will generate a smoother/rougher surface respectively. "
                "One of two methods to specify degree of smoothing."
                );
        maxSmoothingIterationsProp_.setDescription(
                "Maximum number of smoothing steps. "
                "Lower/Higher values will generate a rougher/smoother surface respectively. "
                "One of two methods to specify degree of smoothing."
                );
        smoothingStepSizeProp_.setDescription(
                "Initial step size for the smoothing operation."
                "Larger initial step sizes will lead to faster convergence. "
                "Beware: If the initial step size is too large the algorithm may step across a minimum "
                "and take time to readjust the step size or even diverge. "
                "Tip: Set voreen's log level to debug and observe changes in roughness and step size depending on the initial value. "
                "However, most of the time this value does not need to be changed at all."
                );

    }
    virtual CodeState getCodeState() const        { return CODE_STATE_TESTING;   }

    virtual void initialize();
    virtual void deinitialize();
    virtual void process();

protected:
    /**
     * Extract an isosurface from the given volume using the CPU.
     */
    Geometry* extractSurfaceMCCPU(const VolumeBase* vol, float normalizedIsoValue) const;

    /**
     * Extract an isosurface from the given volume using the GPU.
     */
    Geometry* extractSurfaceGPU(const VolumeBase* vol, tgt::Shader* shader, tgt::svec3 gridSize, float normalizedIsoValue);

    /**
     * Build the shader for extracting an isosurface using the GPU and the marching cubes algorithm.
     */
    void buildMarchingCubesShader();

    /**
     * Build the shader for extracting block surfaces using the GPU.
     */
    void buildBlocksShader();

    /**
     * Delete the smoothed volume and reset the smoothedVolume_ pointer.
     */
    void resetSmoothedVolume();

    /**
     * Adjust the visibility of properties associated with surface extration.
     */
    void adjustExtractionProperties();

    /**
     * Adjust the visibility of properties associated with smoothing a volume.
     */
    void adjustSmoothingProperties();

    /**
     * Return a volume smoother volume so that an isosurface for i=0.5 is still
     * a plausible isosurface for the original volume, but smoother than a direct
     * marching cubes approach.
     * Intended to be used on binary volumes.
     * @param input A binary volume to be smoothed.
     * @param normalizedBinarizationThreshold normalization threshold applied to the input volume
     * @return a smoothed version of input, or nullptr if smoothing failed.
     */
    VolumeRAM* smooth(const VolumeRAM* input, float normalizedBinarizationThreshold);


private:
    enum MarchingCubesProcessingMode {
        MODE_CPU,
        MODE_GPU,
    };

    enum SurfaceType {
        MARCHING_CUBES,
        BLOCKS
    };

    // Ports
    VolumePort inport_;
    GeometryPort outport_;

    // General properties
    BoolProperty enabledProp_;
    FloatProperty isoValueProp_;
    IntProperty meshBufferSizeProp_;
    OptionProperty<SurfaceType> surfaceTypeProp_;

    // Marching cubes properties
    OptionProperty<MarchingCubesProcessingMode> marchingCubesProcessingModeProp_;
    ShaderProperty marchingCubesShaderProp_;

    // Blocks properties
    ShaderProperty blocksShaderProp_;

    // Smoothing
    BoolProperty useSmoothingProp_;
    FloatProperty binarizationThresholdProp_;
    FloatProperty maxGradientProp_;
    IntProperty maxSmoothingIterationsProp_;
    FloatProperty smoothingStepSizeProp_;
    ButtonProperty recalculateSmoothingButtonProp_;
    ProgressProperty smoothingProgressProp_;


    // Buffers for extracting via OpenGL
    GLuint meshBuffer_;
    GLuint edgeIndicesBuffer_;
    GLuint triangleVertexIndicesBuffer_;
    GLuint vao_;

    // A smoothed version of the current volume, if one has been constructed
    std::unique_ptr<VolumeBase> smoothedVolume_;
    bool smoothingForced_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_ISOSURFACEEXTRACTOR_H
