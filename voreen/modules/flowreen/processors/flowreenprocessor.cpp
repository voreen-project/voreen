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

#include "flowreenprocessor.h"

#include "voreen/core/voreenapplication.h"
#include "modules/flowreen/utils/flowmath.h"

#include "tgt/shadermanager.h"

namespace voreen {

const tgt::vec4 FlowreenProcessor::fancyColors_[NUM_COLORS] = {
    tgt::vec4(1.0f, 0.5f, 0.0f, 1.0f),  // orange
    tgt::vec4(1.0f, 0.0f, 0.8f, 1.0f),  // pink
    tgt::vec4(0.8f, 0.0f, 0.8f, 1.0f),  // violett
    tgt::vec4(0.8f, 0.8f, 0.0f, 1.0f),  // yellow
    tgt::vec4(0.0f, 0.8f, 0.8f, 1.0f),  // cyan
    tgt::vec4(0.8f, 0.0f, 0.5f, 1.0f),  // blackberry coloured? => fancy!
    tgt::vec4(0.0f, 0.0f, 0.45f, 1.0f), // dark blue
    tgt::vec4(0.0f, 0.45f, 0.0f, 1.0f), // dark green
    tgt::vec4(0.45f, 0.0f, 0.0f, 1.0f), // dark red
    tgt::vec4(0.3f, 0.3f, 0.3f, 1.0f)   // dark gray
};

// Vertices describe a unit circle in the XY-plane. It is used
// to generate streamtubes.
//
const tgt::vec3 FlowreenProcessor::planeCircleXY[NUM_CIRCLE_VERTICES] = {
    tgt::vec3(1.0f, 0.0f, 0.0f),
    tgt::vec3(0.984808f, 0.173648f, 0.0f),
    tgt::vec3(0.939693f, 0.34202f, 0.0f),
    tgt::vec3(0.866025f, 0.5f, 0.0f),
    tgt::vec3(0.766044f, 0.642788f, 0.0f),
    tgt::vec3(0.642788f, 0.766044f, 0.0f),
    tgt::vec3(0.5f, 0.866025f, 0.0f),
    tgt::vec3(0.34202f, 0.939693f, 0.0f),
    tgt::vec3(0.173648f, 0.984808f, 0.0f),
    tgt::vec3(0.0f, 1.0f, 0.0f),
    tgt::vec3(-0.173648f, 0.984808f, 0.0f),
    tgt::vec3(-0.34202f, 0.939693f, 0.0f),
    tgt::vec3(-0.5f, 0.866025f, 0.0f),
    tgt::vec3(-0.642788f, 0.766044f, 0.0f),
    tgt::vec3(-0.766044f, 0.642788f, 0.0f),
    tgt::vec3(-0.866025f, 0.5f, 0.0f),
    tgt::vec3(-0.939693f, 0.34202f, 0.0f),
    tgt::vec3(-0.984808f, 0.173648f, 0.0f),
    tgt::vec3(-1.0f, 0.0f, 0.0f),
    tgt::vec3(-0.984808f, -0.173648f, 0.0f),
    tgt::vec3(-0.939693f, -0.34202f, 0.0f),
    tgt::vec3(-0.866025f, -0.5f, 0.0f),
    tgt::vec3(-0.766044f, -0.642788f, 0.0f),
    tgt::vec3(-0.642788f, -0.766045f, 0.0f),
    tgt::vec3(-0.5f, -0.866025f, 0.0f),
    tgt::vec3(-0.34202f, -0.939693f, 0.0f),
    tgt::vec3(-0.173648f, -0.984808f, 0.0f),
    tgt::vec3(0.0f, -1.0f, 0.0f),
    tgt::vec3(0.173648f, -0.984808f, 0.0f),
    tgt::vec3(0.34202f, -0.939693f, 0.0f),
    tgt::vec3(0.5f, -0.866025f, 0.0f),
    tgt::vec3(0.642788f, -0.766044f, 0.0f),
    tgt::vec3(0.766044f, -0.642788f, 0.0),
    tgt::vec3(0.866025f, -0.5f, 0.0f),
    tgt::vec3(0.939693f, -0.34202f, 0.0f),
    tgt::vec3(0.984808f, -0.173648f, 0.0f)
};

bool FlowreenProcessor::flowreenInitialized_(false);

FlowreenProcessor::FlowreenProcessor()
    : maxStreamlineLengthProp_("maxStreamlineLength", "max. streamline length: ", 150.0f, 0.0f, 100000.0f),
    thresholdProp_("magnitudeThreshold", "threshold for magnitude (%): ", tgt::vec2(0.0f, 100.0f), tgt::vec2(0.0f), tgt::vec2(100.0f)),
    useBoundingBoxProp_("useBoundingBoxProp", "render bounding box:", true),
    boundingBoxColorProp_("boundingBoxColorProp", "bounding box color:", tgt::vec4(1.0f)),
    useCoordinateAxisProp_("useCoordinateAxisProp", "render coordinate axis:", true),
    currentVolume_(0),
    flowDimensions_(1),
    forceRemap_(true)
{
    maxStreamlineLengthProp_.setStepping(10.0f);
    thresholdProp_.setStepping(tgt::vec2(1.0f));
    thresholdProp_.onChange(MemberFunctionCallback<FlowreenProcessor>(this, &FlowreenProcessor::onThresholdChange));

}

FlowreenProcessor::~FlowreenProcessor() {}

// protected methods
//

std::vector<tgt::vec3> FlowreenProcessor::getTransformedCircle(const tgt::mat4& transformation) const
{
    // rotate all vertices of the plane circle in the origina so that it
    // becomes oriented along the streamline. then translate its center
    // to position given by the the position on the streamline.
    //
    tgt::vec4 r_trans;
    std::vector<tgt::vec3> circle;
    for (size_t i = 0; i < NUM_CIRCLE_VERTICES; ++i) {
        const tgt::vec3& r = planeCircleXY[i];
        r_trans = transformation * tgt::vec4(r, 1.0f);
        circle.push_back(r_trans._xyz);
    }
    return circle;
}

tgt::vec3 FlowreenProcessor::mapToFlowBoundingBox(const tgt::vec3& vertex) const {
    if(forceRemap_) {
        const tgt::vec3 dim(static_cast<tgt::vec3>(tgt::ivec3(flowDimensions_) - tgt::ivec3(1)));
        return (((2.0f * vertex) - dim) / tgt::max<float>(dim));
    }
    return vertex;
}

void FlowreenProcessor::renderArrow(const tgt::mat4& transformation, const float length) const {
    if (length == 0.0f)
        return;

    const float z1 = length / -2.0f;
    const float z2 = (0.75f * length) + z1;
    const float z3 = length / 2.0f;
    const float radius1 = 0.1f;
    const float radius2 = 0.2f;

    // cylinder along the z-axis
    //
    glBegin(GL_QUAD_STRIP);
    for (int k = NUM_CIRCLE_VERTICES; k >= 0; --k) {
        size_t a = k % NUM_CIRCLE_VERTICES;
        size_t b = (k + 1) % NUM_CIRCLE_VERTICES;

        tgt::vec4 r1(planeCircleXY[a].x * radius1, planeCircleXY[a].y * radius1, z1, 1.0f);
        tgt::vec4 r2(r1.x, r1.y, z2, 1.0f);
        tgt::vec4 r3(planeCircleXY[b].x * radius1, planeCircleXY[b].y * radius1, z1, 1.0f);
        r1 = (transformation * r1);
        r2 = (transformation * r2);
        r3 = (transformation * r3);

        tgt::vec3 e1 = FlowMath::normalize(r2.xyz() - r1.xyz());
        tgt::vec3 e2 = FlowMath::normalize(r3.xyz() - r1.xyz());
        tgt::vec3 normal = FlowMath::normalize(tgt::cross(e1, e2));

        glNormal3fv(normal.elem);
        glVertex3fv(mapToFlowBoundingBox(r2.xyz()).elem);
        glVertex3fv(mapToFlowBoundingBox(r1.xyz()).elem);
    }   // for (k
    glEnd();

    tgt::vec4 normal(0.0f, 0.0f, 1.0f, 0.0);
    const tgt::vec3 nz = FlowMath::normalize((transformation * normal).xyz());
    const tgt::vec3 _nz = -nz;

    // two disks covering the cylinder's bottom and
    // a bigger one covering the cone's bottom
    //
    for (size_t i = 0; i < 2; ++i) {
        float radius = radius1;
        float z = z1;
        if (i > 0) {
            radius = radius2;
            z = z2;
        }

        tgt::vec4 r0(0.0f, 0.0f, z, 1.0f);
        glBegin(GL_TRIANGLE_FAN);
            tgt::vec4 r_4(transformation * r0);

            glNormal3fv(_nz.elem);
            glVertex3fv(mapToFlowBoundingBox(r_4.xyz()).elem);

            for (size_t k = 0; k <= NUM_CIRCLE_VERTICES; ++k) {
                tgt::vec3 r(planeCircleXY[(k % NUM_CIRCLE_VERTICES)]);
                r.x *= radius;
                r.y *= radius;
                r.z = z;
                r_4 = transformation * tgt::vec4(r, 1.0f);

                glVertex3fv(mapToFlowBoundingBox(r_4.xyz()).elem);
            }
        glEnd();
    }

    // a cone with its peak at z3
    //
    tgt::vec4 r0 = transformation * tgt::vec4(0.0f, 0.0f, z3, 1.0f);
    glBegin(GL_TRIANGLE_FAN);
        glNormal3fv(nz.elem);
        glVertex3fv(mapToFlowBoundingBox(r0.xyz()).elem);

        for (int k = NUM_CIRCLE_VERTICES; k >= 0; --k) {
            size_t a = k % NUM_CIRCLE_VERTICES;
            size_t b = (k + 1) % NUM_CIRCLE_VERTICES;

            tgt::vec4 r1 = transformation * tgt::vec4(planeCircleXY[a].x * radius2, planeCircleXY[a].y * radius2, z2, 1.0f);
            tgt::vec4 r2 = transformation * tgt::vec4(planeCircleXY[b].x * radius2, planeCircleXY[b].y * radius2, z2, 1.0f);

            tgt::vec3 e1 = FlowMath::normalize(r0.xyz() - r1.xyz());
            tgt::vec3 e2 = FlowMath::normalize(r2.xyz() - r1.xyz());

            glNormal3fv(FlowMath::normalize(tgt::cross(e1, e2)).elem);
            glVertex3fv(mapToFlowBoundingBox(r1.xyz()).elem);
        }
    glEnd();
}

void FlowreenProcessor::renderBoundingBox(const tgt::vec4& color) const {
    tgt::vec3 urf, llb;
    if(forceRemap_) {
        urf = mapToFlowBoundingBox(static_cast<tgt::vec3>(tgt::ivec3(flowDimensions_) - tgt::ivec3(1)));
        llb = urf * -1.0f;
    } else {
        urf = static_cast<tgt::vec3>(tgt::ivec3(flowDimensions_) - tgt::ivec3(1));
        llb = tgt::vec3(0.f);
    }

    glPushAttrib(GL_POLYGON_BIT | GL_CURRENT_BIT);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glColor4fv(color.elem);

    glBegin(GL_QUADS);
        glNormal3f(0.0f, -1.0f, 0.0f);  // bottom
        glVertex3f(llb.x, llb.y, llb.z);
        glVertex3f(urf.x, llb.y, llb.z);
        glVertex3f(urf.x, llb.y, urf.z);
        glVertex3f(llb.x, llb.y, urf.z);

        glNormal3f(0.0f, 1.0f, 0.0f);  // top
        glVertex3f(llb.x, urf.y, llb.z);
        glVertex3f(llb.x, urf.y, urf.z);
        glVertex3f(urf.x, urf.y, urf.z);
        glVertex3f(urf.x, urf.y, llb.z);

        glNormal3f(-1.0f, 0.0f, 0.0f);  // left
        glVertex3f(llb.x, llb.y, llb.z);
        glVertex3f(llb.x, llb.y, urf.z);
        glVertex3f(llb.x, urf.y, urf.z);
        glVertex3f(llb.x, urf.y, llb.z);

        glNormal3f(1.0f, 0.0f, 0.0f);   // right
        glVertex3f(urf.x, llb.y, llb.z);
        glVertex3f(urf.x, urf.y, llb.z);
        glVertex3f(urf.x, urf.y, urf.z);
        glVertex3f(urf.x, llb.y, urf.z);

        glNormal3f(0.0f, 0.0f, 1.0f);   // front
        glVertex3f(llb.x, llb.y, urf.z);
        glVertex3f(urf.x, llb.y, urf.z);
        glVertex3f(urf.x, urf.y, urf.z);
        glVertex3f(llb.x, urf.y, urf.z);

        glNormal3f(0.0f, 0.0f, -1.0f);  // back
        glVertex3f(llb.x, llb.y, llb.z);
        glVertex3f(llb.x, urf.y, llb.z);
        glVertex3f(urf.x, urf.y, llb.z);
        glVertex3f(urf.x, llb.y, llb.z);
    glEnd();

    glPopAttrib();
}

void FlowreenProcessor::renderCoordinateAxis(const tgt::vec3& min, const tgt::vec3& max) const {
    glPushAttrib(GL_POLYGON_BIT | GL_CURRENT_BIT);
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    // render coordinate system
    //
    glBegin(GL_LINES);
        // x-axis
        glColor4f(1.0f, 0.0f, 0.0f, 1.0f);
        glVertex3f(min.x, 0.0f, 0.0f);
        glVertex3f(max.x, 0.0f, 0.0f);

        // y-axis
        glColor4f(0.0f, 1.0f, 0.0f, 1.0f);
        glVertex3f(0.0f, min.y, 0.0f);
        glVertex3f(0.0f, max.y, 0.0f);

        // z-axis
        glColor4f(0.0f, 0.0f, 1.0f, 1.0f);
        glVertex3f(0.0f, 0.0f, min.z);
        glVertex3f(0.0f, 0.0f, max.z);
    glEnd();

    glPopAttrib();
}

// private methods
//

void FlowreenProcessor::init() {
    if (flowreenInitialized_ == true)
        return;
    flowreenInitialized_ = true;
#ifdef VRN_NO_RANDOM
    srand(0);
#else
    srand(static_cast<unsigned int>(time(0)));
#endif

#if defined(_MSC_VER) && !defined(_WIN64)
    int i = _set_SSE2_enable(1);
    if (i != 0)
        std::cout << "SSE2 enabled.\n";
    else
        std::cout << "SSE2 not enabled; processor does not support SSE2.\n";
#endif
}

void FlowreenProcessor::onThresholdChange() {
    tgt::vec2 threshold(thresholdProp_.get());
    if (threshold.x > threshold.y) {
        threshold.y = threshold.x;
        thresholdProp_.set(threshold);
    }
}

}   // namespace
