/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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

/*
 * 4x4 float matrix
 */

typedef struct /*__attribute__((aligned(16)))*/ _mat4{
        float4 x;
        float4 y;
        float4 z;
        float4 w;
} mat4;

/*
 * This struct contains a volume and information about it, as
 * its dimensions and spacing. Additionally, the reciprocal
 * values of all parameters are available (suffix RCP) .
 */
typedef struct _VolumeParameters{

    mat4 physicalToWorldMatrix_;     // (see Volume)
    mat4 worldToPhysicalMatrix_;
    mat4 worldToTextureMatrix_;
    mat4 textureToWorldMatrix_;

    //Note: This is unfortunately not supported by all vendors:
    //sampler3D volume_;              // the actual dataset
    float4 cameraPositionPhysical_;        // camera position in volume object coordinates (see mod_shading.frag)
    float4 lightPositionPhysical_;         // light position in volume object coordinates (see mod_shading.frag)
    float4 datasetDimensions_;        // the dataset's resolution, e.g. [ 256.0, 128.0, 128.0]
    float4 datasetDimensionsRCP_;

    float4 datasetSpacing_;           // set dataset's voxel size, e.g. [ 0.02, 0.05, 0.05]
    float4 datasetSpacingRCP_;
    float4 volumeCubeSize_;           // the volume's size in physical coordinates, e.g. [ 1.0, 0.5, 0.5]
    float4 volumeCubeSizeRCP_;

    float4 volumeOffset_;             // see VolumeHandle::getOffset()
    int bitDepth_;                  // the volume's bit depth
    float rwmScale_;                // RealWorldMapping slope
    float rwmOffset_;               // RealWorldMapping intercept
    int numChannels_;
} VolumeParameters;

float4 texToPhysical(float4 samplePos, VolumeParameters volumeParams) {
    return ((samplePos*volumeParams.datasetDimensions_)*volumeParams.datasetSpacing_)+volumeParams.volumeOffset_;
}
