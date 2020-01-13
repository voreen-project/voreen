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

#version 430

layout(triangles) in;
layout(points, max_vertices = 8) out;

//position and normal in model coordinates
in VertexData {
    vec4 position;
#if defined(USE_COLOR)
    vec4 color;
#endif
#if defined(USE_NORMAL)
    vec3 normal;
#endif
#if defined(USE_TEXCOORD)
    vec2  texCoord;
    ivec2 texIndex;
#endif
} vin[];

// model transformation
layout(location = 0) uniform mat4 modelMatrix;

//clip plane (in world coordinates)
layout(location = 1) uniform vec4 plane;

//matrix for transformation to 2D
layout(location = 2) uniform mat4 xyTransform;

layout(stream = 0) out triangleData {
    vec4 position;
#if defined(USE_COLOR)
    vec4 color;
#endif
#if defined(USE_NORMAL)
    vec3 normal;
#endif
#if defined(USE_TEXCOORD)
    vec2  texCoord;
    ivec2 texIndex;
#endif
} tOut;

layout(stream = 1) out edgeData {
    vec4 position;
#if defined(USE_COLOR)
    vec4 color;
#endif
#if defined(USE_NORMAL)
    vec3 normal;
#endif
#if defined(USE_TEXCOORD)
    vec2  texCoord;
    ivec2 texIndex;
#endif
    vec2 pos2d;
} eOut;

//constant for numerical stability
//const float clipEpsilon = 0.00001;

void main() {
    // distances to the plane (compute using world coordinates)
    vec3 dist = vec3(dot(plane.xyz, gl_in[0].gl_Position.xyz) - plane.w, 
                     dot(plane.xyz, gl_in[1].gl_Position.xyz) - plane.w, 
                     dot(plane.xyz, gl_in[2].gl_Position.xyz) - plane.w);
  
    if (all(greaterThan/*Equal*/(dist,/*vec3(-clipEpsilon)*/vec3(0.0)))) {
        // Case 1 (none clipped) - output original triangle 
        // (in model coordinates)
        for(int i = 0; i < 3; i++) {
            tOut.position = vin[i].position;
#if defined(USE_COLOR)
            tOut.color = vin[i].color;
#endif
#if defined(USE_NORMAL)            
            tOut.normal = vin[i].normal;
#endif
#if defined(USE_TEXCOORD)
            tOut.texCoord = vin[i].texCoord;
            tOut.texIndex = vin[i].texIndex;
#endif
            EmitStreamVertex(0);
            EndStreamPrimitive(0);
        }
        return;
    }
    
    if (!any(greaterThan/*Equal*/(dist, /*vec3(clipEpsilon)*/vec3(0.0))))
        // Case 2 (all clipped) - no output
        return;

    ivec3 indices = ivec3(0,1,2);

    // There are either 1 or 2 vertices above the clipping plane.
    bvec3 above = greaterThanEqual(dist, vec3(0.0));
    bool nextIsAbove;
    
    // Find the CCW - most vertex above the plane.
    if (above[1] && !above[0]) {
        // Cycle once CCW .
        nextIsAbove = above[2];
        indices = ivec3(1,2,0);
        dist = dist.yzx ;
    } else if (above[2] && !above[1]) {
        // Cycle once CW.
        nextIsAbove = above[0];
        indices = ivec3(2,0,1);
        dist = dist.zxy ;
    } else 
        nextIsAbove = above[1];

    //values for clipped edge vertices
    vec4 ePos1, ePos2;
#if defined(USE_COLOR)    
    vec4 eColor1, eColor2;
#endif 
#if defined(USE_NORMAL)
    vec3 eNormal1, eNormal2;
#endif
#if defined(USE_TEXCOORD)
    vec2 eTexCoord1, eTexCoord2;
    ivec2 eTexIndex1, eTexIndex2;
#endif

    // We always need to clip v2 - v0.
    ePos2 = mix(vin[indices[0]].position, vin[indices[2]].position, 
	        dist[0] / (dist[0] - dist[2]));
#if defined(USE_COLOR)	        
    eColor2 = mix(vin[indices[0]].color, vin[indices[2]].color, 
		dist[0] / (dist[0] - dist[2]));
#endif
#if defined(USE_NORMAL)
    eNormal2 = mix(vin[indices[0]].normal, vin[indices[2]].normal, 
		dist[0] / (dist[0] - dist[2]));
#endif
#if defined(USE_TEXCOORD	)
    eTexCoord2 = mix(vin[indices[0]].texCoord, vin[indices[2]].texCoord, 
		dist[0] / (dist[0] - dist[2]));
    eTexIndex2 = ivec2(round(mix(vin[indices[0]].texIndex, vin[indices[2]].texIndex, 
		dist[0] / (dist[0] - dist[2]))));
#endif 

    if (nextIsAbove) {
        // geometric case 3      
        ePos1 = mix(vin[indices[1]].position, vin[indices[2]].position, 
		    dist[1] / (dist[1] - dist[2]));
#if defined(USE_COLOR)
	eColor1 = mix(vin[indices[1]].color, vin[indices[2]].color, 
		    dist[1] / (dist[1] - dist[2]));
#endif
#if defined(USE_NORMAL)
        eNormal1 = mix(vin[indices[1]].normal, vin[indices[2]].normal, 
		    dist[1] / (dist[1] - dist[2]));
#endif
#if defined(USE_TEXCOORD	)
	eTexCoord1 = mix(vin[indices[1]].texCoord, vin[indices[2]].texCoord, 
		dist[1] / (dist[1] - dist[2]));
	eTexIndex1 = ivec2(round(mix(vin[indices[1]].texIndex, vin[indices[2]].texIndex, 
		dist[1] / (dist[1] - dist[2]))));
#endif 


        //v0
        tOut.position = vin[indices[0]].position;
#if defined(USE_COLOR)
        tOut.color = vin[indices[0]].color;
#endif
#if defined(USE_NORMAL)
        tOut.normal = vin[indices[0]].normal;
#endif
#if defined(USE_TEXCOORD)
	tOut.texCoord = vin[indices[0]].texCoord;
	tOut.texIndex = vin[indices[0]].texIndex;
#endif
        EmitStreamVertex(0);
        EndStreamPrimitive(0);

        //v1
        tOut.position = vin[indices[1]].position;
#if defined(USE_COLOR)
        tOut.color = vin[indices[1]].color;
#endif
#if defined(USE_NORMAL)
        tOut.normal = vin[indices[1]].normal;
#endif
#if defined(USE_TEXCOORD)
	tOut.texCoord = vin[indices[1]].texCoord;
	tOut.texIndex = vin[indices[1]].texIndex;
#endif
        EmitStreamVertex(0);
        EndStreamPrimitive(0);

        //v2   
        tOut.position = ePos1; 
#if defined(USE_COLOR)
        tOut.color = eColor1;
#endif
#if defined(USE_NORMAL)
        tOut.normal = eNormal1;
#endif
#if defined(USE_TEXCOORD)
	tOut.texCoord = eTexCoord1;
	tOut.texIndex = eTexIndex1;
#endif
        EmitStreamVertex(0);
        EndStreamPrimitive(0);
 
        //v0
        tOut.position = vin[indices[0]].position;
#if defined(USE_COLOR)
        tOut.color = vin[indices[0]].color;
#endif
#if defined(USE_NORMAL)
        tOut.normal = vin[indices[0]].normal;
#endif
#if defined(USE_TEXCOORD)
	tOut.texCoord = vin[indices[0]].texCoord;
	tOut.texIndex = vin[indices[0]].texIndex;
#endif
        EmitStreamVertex(0);
        EndStreamPrimitive(0);

        //v2   
        tOut.position = ePos1; 
#if defined(USE_COLOR)
        tOut.color = eColor1;
#endif
#if defined(USE_NORMAL)
        tOut.normal = eNormal1;
#endif
#if defined(USE_TEXCOORD)
	tOut.texCoord = eTexCoord1;
	tOut.texIndex = eTexIndex1;
#endif
        EmitStreamVertex(0);
        EndStreamPrimitive(0);

        //v3
        tOut.position = ePos2; 
#if defined(USE_COLOR)
        tOut.color = eColor2;
#endif
#if defined(USE_NORMAL)
        tOut.normal = eNormal2;
#endif
#if defined(USE_TEXCOORD)
	tOut.texCoord = eTexCoord2;
	tOut.texIndex = eTexIndex2;
#endif
        EmitStreamVertex(0);
        EndStreamPrimitive(0);
    } else {
        // geometric case 4    
        ePos1 = mix(vin[indices[0]].position, vin[indices[1]].position, 
		    dist[0] / (dist[0] - dist[1]));
#if defined(USE_COLOR)
	eColor1 = mix(vin[indices[0]].color, vin[indices[1]].color, 
		    dist[0] / (dist[0] - dist[1]));
#endif
#if defined(USE_NORMAL)
        eNormal1 = mix(vin[indices[0]].normal, vin[indices[1]].normal, 
		    dist[0] / (dist[0] - dist[1]));
#endif
#if defined(USE_TEXCOORD	)
	eTexCoord1 = mix(vin[indices[0]].texCoord, vin[indices[1]].texCoord, 
		dist[0] / (dist[0] - dist[1]));
	eTexIndex1 = ivec2(round(mix(vin[indices[0]].texIndex, vin[indices[1]].texIndex, 
		dist[0] / (dist[0] - dist[1]))));
#endif 

        //v1   
        tOut.position = ePos1; 
#if defined(USE_COLOR)
        tOut.color = eColor1;
#endif
#if defined(USE_NORMAL)
        tOut.normal = eNormal1;
#endif
#if defined(USE_TEXCOORD)
        tOut.texCoord = eTexCoord1;
        tOut.texIndex = eTexIndex1;
#endif
        EmitStreamVertex(0);
        EndStreamPrimitive(0);

        //v2
        tOut.position = ePos2; 
#if defined(USE_COLOR)
        tOut.color = eColor2;
#endif
#if defined(USE_NORMAL)
        tOut.normal = eNormal2;
#endif
#if defined(USE_TEXCOORD)
        tOut.texCoord = eTexCoord2;
        tOut.texIndex = eTexIndex2;
#endif
        EmitStreamVertex(0);
        EndStreamPrimitive(0);

        //v0
        tOut.position = vin[indices[0]].position;
#if defined(USE_COLOR)
        tOut.color = vin[indices[0]].color;
#endif
#if defined(USE_NORMAL)
        tOut.normal = vin[indices[0]].normal;
#endif
#if defined(USE_TEXCOORD)
        tOut.texCoord = vin[indices[0]].texCoord;
        tOut.texIndex = vin[indices[0]].texIndex;
#endif
        EmitStreamVertex(0);
        EndStreamPrimitive(0);     
    }

    // write out the clipping edge to the other stream
    // correct orientation of edge is p2 -> p1
    eOut.position = ePos2;
#if defined(USE_COLOR)
    eOut.color = eColor2;
#endif
#if defined(USE_NORMAL)
    eOut.normal = eNormal2;
#endif
#if defined(USE_TEXCOORD)
    eOut.texCoord = eTexCoord2;
    eOut.texIndex = eTexIndex2;
#endif
    eOut.pos2d = (xyTransform * (modelMatrix * ePos2)).xy;
    EmitStreamVertex(1);
    EndStreamPrimitive(1);     

    eOut.position = ePos1; 
#if defined(USE_COLOR)
    eOut.color = eColor1;
#endif
#if defined(USE_NORMAL)
    eOut.normal = eNormal1;
#endif
#if defined(USE_TEXCOORD)
    eOut.texCoord = eTexCoord1;
    eOut.texIndex = eTexIndex1;
#endif
    eOut.pos2d = (xyTransform * (modelMatrix * ePos1)).xy;
    EmitStreamVertex(1);
    EndStreamPrimitive(1);
}
