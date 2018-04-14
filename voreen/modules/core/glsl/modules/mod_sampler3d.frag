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

 #define MOD_SAMPLER3D

/*
 * This struct contains a volume and information about it, as
 * its dimensions and spacing. Additionally, the reciprocal
 * values of all parameters are available (suffix RCP) .
 * The values are set automatically by
 * VolumeRenderer::bindVolumes() if necessary.
 */
struct VolumeParameters {

    //Note: This is unfortunately not supported by all vendors:
    //sampler3D volume_;              // the actual dataset

    vec3 datasetDimensions_;        // the dataset's resolution, e.g. [ 256.0, 128.0, 128.0]
    vec3 datasetDimensionsRCP_;

    vec3 datasetSpacing_;           // set dataset's voxel size, e.g. [ 0.02, 0.05, 0.05]
    vec3 datasetSpacingRCP_;

    vec3 volumeCubeSize_;           // the volume's size in physical coordinates, e.g. [ 1.0, 0.5, 0.5]
    vec3 volumeCubeSizeRCP_;

    vec3 volumeOffset_;             // see VolumeHandle::getOffset()

    int bitDepth_;                  // the volume's bit depth

    float rwmScale_;                // RealWorldMapping slope
    float rwmOffset_;               // RealWorldMapping intercept

    int numChannels_;

    mat4 physicalToWorldMatrix_;     // (see Volume)
    mat4 worldToPhysicalMatrix_;

    mat4 worldToTextureMatrix_;
    mat4 textureToWorldMatrix_;

    vec3 cameraPositionPhysical_;        // camera position in volume object coordinates (see mod_shading.frag)
    vec3 lightPositionPhysical_;         // light position in volume object coordinates (see mod_shading.frag)
};

/*
 * Function for volume texture lookup. In addition to the texture coordinates
 * the corresponding VolumeParameters has to be passed.
 * Before returning the fetched value it is normalized to the interval [0,1], in order to deal
 * with 12 bit data sets.
 */
//vec4 textureLookup3D(VolumeParameters volumeStruct, vec3 texCoords) {
    //vec4 result;
//#if defined(GLSL_VERSION_130)
    //result = texture(volumeStruct.volume_, texCoords);
//#else
    //result = texture3D(volumeStruct.volume_, texCoords);
//#endif
    //result.a *= volumeStruct.rwmScale_;
    //result.a += volumeStruct.rwmOffset_;
    //return result;
//}

/*
 * DEPRECATED
 * Note: This function is re-activated until all vendors support samplers within GLSL structs.
 */
vec4 textureLookup3D(sampler3D volume, VolumeParameters volumeStruct, vec3 texCoords) {
    //return textureLookup3D(volumeStruct, texCoords);
    vec4 result;
#if defined(GLSL_VERSION_130)
    result = texture(volume, texCoords);
#else
    result = texture3D(volume, texCoords);
#endif
    result *= volumeStruct.rwmScale_;
    result += volumeStruct.rwmOffset_;
    return result;
}

/*
 * Function for volume texture lookup. In addition to the texture coordinates
 * the corresponding VolumeParameters has to be passed .
 * In contrast to textureLookup3D() this function does not normalize the intensity values,
 * in order to deal with 12 bit data sets.
 */
//vec4 textureLookup3DUnnormalized(VolumeParameters volumeStruct, vec3 texCoords) {
//#if defined(GLSL_VERSION_130)
    //return texture(volumeStruct.volume_, texCoords);
//#else
    //return texture3D(volumeStruct.volume_, texCoords);
//#endif
//}

/*
 * DEPRECATED
 * Note: This function is re-activated until all vendors support samplers within GLSL structs.
 */
vec4 textureLookup3DUnnormalized(sampler3D volume, VolumeParameters volumeStruct, vec3 texCoords) {
    //return textureLookup3DUnnormalized(volumeStruct, texCoords);
#if defined(GLSL_VERSION_130)
    return texture(volume, texCoords);
#else
    return texture3D(volume, texCoords);
#endif
}


/*
 * This function should be called by all raycasters in order to get the intensity from the volume.
 * In cases where volumeStruct indicates, that gradients are stored in the volume, these are
 * also fetched from the volume. Therefore, in addition to the volume and the texture coordinates
 * the corresponding VolumeParameters has to be passed.
 * Before returning the intensity value it is normalized to the interval [0,1], in order to deal
 * with 12 bit data sets.
 *
 * @return The result vec4 stores the intensity in the a channel and the gradient (if available)
 *         in the rgb channels.
 *
 */
//vec4 getVoxel(VolumeParameters volumeStruct, vec3 samplePos) {
    //return textureLookup3D(volumeStruct.volume_, volumeStruct, samplePos);
//}

/*
 * DEPRECATED
 * Note: This function is re-activated until all vendors support samplers within GLSL structs.
 */
vec4 getVoxel(sampler3D volume, VolumeParameters volumeStruct, vec3 samplePos) {
    return textureLookup3D(volume, volumeStruct, samplePos);
}


/*
 * This function should be called by all raycasters in order to get the intensity from the volume.
 * In cases where volumeStruct indicates, that gradients are stored in the volume, these are
 * also fetched from the volume. Therefore, in addition to the volume and the texture coordinates
 * the corresponding VolumeParameters has to be passed.
 * In contrast to getVoxel() this function does not normalize the intensity values, in order to deal
 * with 12 bit data sets.
 *
 * @return The result vec4 stores the intensity in the a channel and the gradient (if available)
 *         in the rgb channels.
 */
//vec4 getVoxelUnnormalized(VolumeParameters volumeStruct, vec3 samplePos) {
    //return textureLookup3DUnnormalized(volumeStruct.volume_, volumeStruct, samplePos);
//}

/*
 * DEPRECATED
 * Note: This function is re-activated until all vendors support samplers within GLSL structs.
 */
vec4 getVoxelUnnormalized(sampler3D volume, VolumeParameters volumeStruct, vec3 samplePos) {
    return textureLookup3DUnnormalized(volume, volumeStruct, samplePos);
}

vec3 texToPhysical(vec3 samplePos, VolumeParameters volumeParams) {
    return ((samplePos*volumeParams.datasetDimensions_ - vec3(0.5))*volumeParams.datasetSpacing_)+volumeParams.volumeOffset_;
}

vec3 physicalToTex(vec3 samplePos, VolumeParameters volumeStruct) {
    return (((samplePos - volumeStruct.volumeOffset_) * volumeStruct.datasetSpacingRCP_) + vec3(0.5)) * volumeStruct.datasetDimensionsRCP_;
}

bool inUnitCube(vec3 sample_val) {
    if(any(greaterThan(sample_val, vec3(1.0))) || any(lessThan(sample_val, vec3(0.0))))
        return false;
    return true;
}

vec3 worldToTex(vec3 worldSamplePos, VolumeParameters volumeStruct) {
    return (volumeStruct.worldToTextureMatrix_* vec4(worldSamplePos, 1.0)).xyz;
}

vec3 texToWorld(vec3 texSamplePos, VolumeParameters volumeStruct) {
    return (volumeStruct.textureToWorldMatrix_* vec4(texSamplePos, 1.0)).xyz;
}

