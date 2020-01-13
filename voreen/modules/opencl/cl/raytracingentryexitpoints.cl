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

#pragma OPENCL EXTENSION cl_khr_byte_addressable_store: enable

const sampler_t smp = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_NEAREST;
const sampler_t smpLin = CLK_NORMALIZED_COORDS_FALSE | CLK_ADDRESS_CLAMP_TO_EDGE | CLK_FILTER_LINEAR;

void writeFloat(image2d_t img, int2 p, float4 value) {
    if(get_image_channel_data_type(img) == CLK_UNSIGNED_INT16) {
        uint4 write;
        write.x = value.x * 65535.0f;
        write.y = value.y * 65535.0f;
        write.z = value.z * 65535.0f;
        write.w = value.w * 65535.0f;
        write_imageui(img, p, write);
    }
    else {
        write_imagef(img, p, value);
    }
}

typedef struct _mat4 {
        float4 x,y,z,w;
} mat4;


typedef struct _Triangle {
    float4 a;
    float4 b;
    float4 c;
    float4 a_tc;
    float4 b_tc;
    float4 c_tc;
} Triangle;

typedef struct _Ray {
    float4 start;
    float4 dir;
} Ray;

typedef struct _IntersectionResult {
    bool intersect;
    float u;
    float v;
    float t;
} IntersectionResult;

__constant float epsilon = 0.00001f;

IntersectionResult rayIntersectsTriangle(Ray r, Triangle tri) {
    IntersectionResult res;
    res.intersect = false;

    float4 e1 = tri.b - tri.a;
    float4 e2 = tri.c - tri.a;
    float4 p = cross(r.dir, e2);
    float a = dot(e1,p);

    if ((a > -epsilon) && (a < epsilon))
        return res;

    float f = 1/a;
    float4 s = r.start - tri.a;
    res.u = f * dot(s,p);

    if ((res.u < 0.0) || (res.u > 1.0))
        return res;

    float4 q = cross(s, e1);
    res.v = f * dot(r.dir, q);
    if ((res.v < 0.0) || ((res.u + res.v) > 1.0))
        return res;

    // at this stage we can compute t to find out where
    // the intersection point is on the line
    res.t = f * dot(e2,q);
    //if (res.t > epsilon) { // ray intersection
        res.intersect = true;
        return res;
    //}
    //else // this means that there is a line intersection but not a ray intersection
        //return res;
}

//matrix vector n
float4 matrixVectorMult (__global mat4* matrix, float4 vector) {
    float4 result;
    result.x = dot(matrix->x, vector);
    result.y = dot(matrix->y, vector);
    result.z = dot(matrix->z, vector);
    result.w = dot(matrix->w, vector);
    return result;
}

float4 getTexCoord(Triangle t, float u, float v) {
    float4 tc;
    tc = ((1.0 - u - v) * t.a_tc) + (u * t.b_tc) + (v * t.c_tc);
    return tc;
}

__kernel void raytrace(write_only image2d_t entry, write_only image2d_t exit, __global Triangle* triangles, int numTriangles, float4 pos, float4 up, float4 look, float4 strafe, __global mat4* transformation_)
{
    int x = get_global_id(0);
    int y = get_global_id(1);

    Ray r;
    r.start = pos;
    r.dir = look;

    int width = get_image_width(entry);
    int height = get_image_height(exit);

    r.dir += ((((float)(x) / (float)(width)) * 2.0f) - 1.0f) * strafe;
    r.dir += ((((float)(y) / (float)(height)) * 2.0f) - 1.0f) * up;

    r.dir = normalize(r.dir);

    bool intersects = false;
    float4 entryColor = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
    float4 exitColor = (float4)(0.0f, 0.0f, 0.0f, 0.0f);
    float tEntry;
    float tExit;

    for(int i=0; i<numTriangles; i++) {
        Triangle tri = triangles[i];
        tri.a = matrixVectorMult(transformation_, tri.a);
        tri.b = matrixVectorMult(transformation_, tri.b);
        tri.c = matrixVectorMult(transformation_, tri.c);
        IntersectionResult intersectRes = rayIntersectsTriangle(r, tri);
        if(intersectRes.intersect) {
            if(intersects) {
                if(intersectRes.t < tEntry) {
                    entryColor = getTexCoord(tri, intersectRes.u, intersectRes.v);
                    tEntry = intersectRes.t;
                }

                if(intersectRes.t > tExit) {
                    exitColor = getTexCoord(tri, intersectRes.u, intersectRes.v);
                    tExit = intersectRes.t;
                }
            }
            else {
                intersects = true;
                entryColor = getTexCoord(tri, intersectRes.u, intersectRes.v);
                exitColor = getTexCoord(tri, intersectRes.u, intersectRes.v);
                tEntry = intersectRes.t;
                tExit = intersectRes.t;
            }
        }
    }

    if(intersects) {
        //camera inside volume:
        if(tEntry < 0.0f) {
            float tDiff = tExit - tEntry;
            entryColor = mix(entryColor, exitColor, -tEntry / tDiff);
        }

        entryColor.w = 1.0f;
        exitColor.w = 1.0f;
    }
    writeFloat(entry, (int2)(x,y), entryColor);
    writeFloat(exit, (int2)(x,y), exitColor);
}

