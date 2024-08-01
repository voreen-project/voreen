/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#include "modules/mod_sampler2d.frag"
#include "modules/mod_sampler3d.frag"
#include "modules/mod_transfunc.frag"

#include "modules/mod_raysetup.frag"

#include "modules/mod_depth.frag"
#include "modules/mod_compositing.frag"
#include "modules/mod_gradients.frag"
#include "modules/mod_shading.frag"

/*******************/
/* Declare Inports */
/*******************/
//$ @name = "inport", @depthtex = "inputTexDepth_", @params = "inputTexParams_"
uniform sampler2D inputTex_;
uniform sampler2D inputTexDepth_;
uniform TextureParameters inputTexParams_;

//$ @name = "volumeport", @params = "volParams_"
uniform VolumeParameters volumeStruct_;
uniform sampler3D volume_;

uniform float samplingStepSize_;

/********************/
/* Declare Outports */
/********************/

//$ @name = "outport0"
out vec4 FragData;

//$ @name = "outport1"
out vec4 FragData1;


/**********************/
/* Declare Properties */
/**********************/
//$ @name = "test", @label = "Test property"
uniform float test_;

/***
 * The main method.
 ***/
void main() {
    vec2 p = gl_FragCoord.xy * screenDimRCP_;
    FragData = vec4(0.0, test_, 0.0, 1.0);
    FragData1 = vec4(1.0, 1.0, 0.0, 1.0)*texture2D(inputTex_, p);
    gl_FragDepth = texture2D(inputTexDepth_, p).z;
}
