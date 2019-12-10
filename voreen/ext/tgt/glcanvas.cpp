/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2019 University of Muenster, Germany,           *
 * Department of Computer Science.                                    *
 *                                                                    *
 * This file is part of the tgt library. This library is free         *
 * software; you can redistribute it and/or modify it under the terms *
 * of the GNU Lesser General Public License version 2.1 as published  *
 * by the Free Software Foundation.                                   *
 *                                                                    *
 * This library is distributed in the hope that it will be useful,    *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of     *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the       *
 * GNU Lesser General Public License for more details.                *
 *                                                                    *
 * You should have received a copy of the GNU Lesser General Public   *
 * License in the file "LICENSE.txt" along with this library.         *
 * If not, see <http://www.gnu.org/licenses/>.                        *
 *                                                                    *
 **********************************************************************/

#include "tgt/glcanvas.h"
#include "tgt/camera.h"
#include "tgt/painter.h"
#include "tgt/glcontextmanager.h"
#include "tgt/tgt_gl.h"

#include <cstdlib>
#include <stdio.h>

namespace tgt {

GLCanvas::GLCanvas(const std::string& title,
                   const ivec2& size,
                   const GLCanvas::Buffers buffers)
    : GLContextBase(title)
    , title_(title)
    , size_(size)
    , buffers_(buffers)
    , rgbaSize_(0, 0, 0, 0)
    , depthSize_(0)
    , stencilSize_(0)
    , doubleBuffered_((buffers & DOUBLE_BUFFER) != 0)
    , quadBuffered_((buffers & QUAD_BUFFER) != 0)
    , fullscreen_(false)
    , autoFlush_(true)
    , initialized_(false)
    , painter_(0)
    , eventHandler_(new EventHandler())
{
}

GLCanvas::~GLCanvas() {
    eventHandler_->clear();
    delete eventHandler_;
    delete painter_;
}

void GLCanvas::sizeChanged(const ivec2& size) {
    size_ = size;
    if (painter_) {
        GLContextStateGuard guard(this);
        painter_->sizeChanged(size);
    }
}

void GLCanvas::paint() {
    if(initialized_ && painter_) {
        GLContextStateGuard guard(this);
        painter_->paint();
        if (autoFlush_) {
            finalizeRendering();
        }
    }
}

void GLCanvas::finalizeRendering() {
    if (doubleBuffered_)
        swap();
    else
        glFlush();
}

void GLCanvas::init() {
    tgtAssert(!initialized_,"canvas has been initialized twice!");
    if(painter_) {
        GLContextStateGuard guard(this);
        painter_->initialize();
    }
    initialized_ = true;
}

void GLCanvas::setPainter(Painter* p) {
    painter_ = p;
}

EventHandler* GLCanvas::getEventHandler() const {
    return eventHandler_;
}
void GLCanvas::setEventHandler(EventHandler* handler){
    delete eventHandler_;
    eventHandler_ = handler;
}

bool GLCanvas::getAutoFlush() {
    return autoFlush_;
}

bool GLCanvas::takeScreenshot(std::string fname) {
    GLContextStateGuard guard(this);

    //taken from gamedev.net forum:
    FILE *pFile;               // The file pointer.
    unsigned char uselessChar; // used for useless char.
    short int uselessInt;      // used for useless int.
    unsigned char imageType;   // Type of image we are saving.
    int index;                 // used with the for loop.
    unsigned char bits;    // Bit depth.
    long Size;                 // Size of the picture.
    int colorMode;
    unsigned char tempColors;
    short int width = getWidth();
    short int height = getHeight();

    // Open file for output.
    pFile = fopen(fname.c_str(), "wb");

   // Check if the file opened or not.
    if (!pFile) {
        return false;
    }

    unsigned char *image = new unsigned char[getWidth()*getHeight()*3];
    //read image from gl:
    glReadPixels( 0, 0, getWidth(), getHeight(), GL_RGB, GL_UNSIGNED_BYTE, image );
    // Set the image type, the color mode, and the bit depth.
    imageType = 2; colorMode = 3; bits = 24;

    // Set these two to 0.
    uselessChar = 0; uselessInt = 0;

    // Write useless data.
    size_t written = 0;
    written += fwrite(&uselessChar, sizeof(unsigned char), 1, pFile);
    written += fwrite(&uselessChar, sizeof(unsigned char), 1, pFile);

    // Now image type.
    written += fwrite(&imageType, sizeof(unsigned char), 1, pFile);

    // Write useless data.
    written += fwrite(&uselessInt, sizeof(short int), 1, pFile);
    written += fwrite(&uselessInt, sizeof(short int), 1, pFile);
    written += fwrite(&uselessChar, sizeof(unsigned char), 1, pFile);
    written += fwrite(&uselessInt, sizeof(short int), 1, pFile);
    written += fwrite(&uselessInt, sizeof(short int), 1, pFile);

    // Write the size that you want.
    written += fwrite(&width, sizeof(short int), 1, pFile);
    written += fwrite(&height, sizeof(short int), 1, pFile);
    written += fwrite(&bits, sizeof(unsigned char), 1, pFile);

    // Write useless data.
    written += fwrite(&uselessChar, sizeof(unsigned char), 1, pFile);

    // Get image size.
    Size = width * height * colorMode;

    // Now switch image from RGB to BGR.
    for (index = 0; index < Size; index += colorMode){
        tempColors = image[index];
        image[index] = image[index + 2];
        image[index + 2] = tempColors;
    }

    // Finally write the image.
    written += fwrite(image, sizeof(unsigned char), Size, pFile);

    // close the file.
    fclose(pFile);
    delete[] image;

    return (written > 0);
}

const ivec4& GLCanvas::getRgbaSize() const {
    return rgbaSize_;
}

int GLCanvas::getDepthSize() const {
    return depthSize_;
}

int GLCanvas::getStencilSize() const {
    return stencilSize_;
}

ivec2 GLCanvas::getSize() const {
    return size_;
}
ivec2 GLCanvas::getPhysicalSize() const {
    return getSize();
}
int GLCanvas::getPhysicalWidth() const {
    return getPhysicalSize().x;
}
int GLCanvas::getPhysicalHeight() const {
    return getPhysicalSize().y;
}

int GLCanvas::getWidth() const {
    return size_.x;
}

int GLCanvas::getHeight() const {
    return size_.y;
}

GLCanvas::Buffers GLCanvas::getBuffers() const {
    return buffers_;
}

} // namespace
