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

#include "lp_graphic.h"
//include opengl
#include "tgt/tgt_gl.h"


namespace voreen {

void LP_Graphic::drawSphere(const tgt::vec3 position, float size, const tgt::vec4 color) {
    MatStack.pushMatrix();
    glPushAttrib(GL_COLOR_BUFFER_BIT);
        glColor4fv(color.elem);
        GLUquadric* quadric = gluNewQuadric();
        MatStack.translate(position.x, position.y, position.z);
        gluSphere(quadric,size,32,32);
        gluDeleteQuadric(quadric);
    glPopAttrib();
    MatStack.popMatrix();
}

void LP_Graphic::drawLineInterval(float start, float end, float length, float size, const tgt::vec4 color) {
    //test, if interval can be drawn
    if(start > length) return;
    if(end > length) end = length;
    if(start < 0) start = 0;
    //pre settings
    GLUquadric* quadric = gluNewQuadric();
    float tailSize = 3*size;
    float tailLength = 2*size;
    float headSize = 3*size;
    float headLength = 0.1f*length;

    MatStack.pushMatrix();
    glPushAttrib(GL_COLOR_BUFFER_BIT);
    glColor4fv(color.elem);
    //if interval is complete in arrowtail
    if(end < tailLength) {
        MatStack.translate(0,0,start);
        gluCylinder(quadric, tailSize, tailSize, end-start, 32, 32);
        gluDeleteQuadric(quadric);
        glPopAttrib();
        MatStack.popMatrix();
        return;
    }
    //if interval is complete in arrowhead
    if(start > (length-headLength)) {
        MatStack.translate(0,0,start);
        gluCylinder(quadric, headSize*(length-start)/headLength, headSize*(length-end)/headLength, end-start, 32, 32);
        gluDeleteQuadric(quadric);
        glPopAttrib();
        MatStack.popMatrix();
        return;
    }
    //if tail is in interval
    if(start <= tailLength) {
        MatStack.pushMatrix();
            MatStack.translate(0,0,start);
            gluCylinder(quadric, tailSize, tailSize, tailLength-start, 32, 32);
            MatStack.translate(0,0,tailLength-start);
            gluDisk(quadric, 0.0, tailSize, 32, 32);
            start = tailLength;
        MatStack.popMatrix();
    }
    //if head is in interval
    if(end >= (length-headLength)) {
        MatStack.pushMatrix();
            MatStack.translate(0,0,(length-headLength));
            gluCylinder(quadric, headSize,headSize*(length-end)/headLength, end-(length-headLength), 32, 32);
            gluDisk(quadric, 0.0, headSize, 32, 32);
            end = length-headLength;
        MatStack.popMatrix();
    }
    //darw middle part
    MatStack.translate(0,0,start);
    gluCylinder(quadric, size, size, end-start, 32, 32);

    gluDeleteQuadric(quadric);
    glPopAttrib();
    MatStack.popMatrix();
}

} // namespace voreen
