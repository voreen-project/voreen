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

#include "modules/mod_sampler3d.frag"           // contains struct VolumeParameters

uniform VolumeParameters volumeStruct_;
uniform sampler3D volume_;    // additional information about the volume the eep are generated for

out vec4 eep_;

/**
 * Similar to eep_simple.vert, but uses texture coordinates instead of vertex positions.
 * This allows for deformations when texture coordinate and vertex position differ. The
 * w-component of the texture coordinates is written to the alpha channel.
 */
void main() {
    eep_.xyz = physicalToTex(gl_MultiTexCoord0.xyz, volumeStruct_);
    eep_ = vec4(clamp(eep_.xyz, 0.0, 1.0), gl_MultiTexCoord0.w);

    // Brackets are important! Without them the compiler may do costly matrix multiplications
    // instead of a cheap matrix-vector multiplication.
    //
    // On NVIDIA 97.55 (Linux) it also triggers a compiler(?) bug, resulting in black rendering
    // of the entry params in TexCoordEntryExitPoints.
    gl_Position = gl_ProjectionMatrix * (gl_ModelViewMatrix * gl_Vertex);
}
