/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2013 Visualization and Computer Graphics Group, *
 * Department of Computer Science, University of Muenster, Germany.   *
 * <http://viscg.uni-muenster.de>                                     *
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

#include "maccorecontext.h"

namespace tgt {
    MacCoreContext::MacCoreContext(const QGLFormat& format, GLCanvas::Buffers buffers):
    QGLContext(format),
    buffers_(buffers)
    {
        
    }
    
    void* MacCoreContext::chooseMacVisual(GDHandle handle) {
        int max = 40;
        int i = 0;
        NSOpenGLPixelFormatAttribute attribs[max];
        
        // we always want core profile
        attribs[i++] = NSOpenGLPFAOpenGLProfile;
        attribs[i++] = NSOpenGLProfileVersion3_2Core;
        
        // wet set alpha buffer to 8 bit if one is required
        attribs[i++] = NSOpenGLPFAAlphaSize;
        if (buffers_ & GLCanvas::ALPHA_BUFFER) {
            attribs[i++] = 8;
        } else {
            attribs[i++] = 0;
        }
        
        // wet set depth buffer to 32 bit if one is required
        attribs[i++] = NSOpenGLPFADepthSize;
        if (buffers_ & GLCanvas::DEPTH_BUFFER) {
            attribs[i++] = 32;
        } else {
            attribs[i++] = 0;
        }
        
        if (buffers_ & GLCanvas::DOUBLE_BUFFER) {
            attribs[i++] = NSOpenGLPFADoubleBuffer;
        }
        
        // wet set stencil buffer to 32 bit if one is required
        attribs[i++] = NSOpenGLPFAStencilSize;
        if (buffers_ & GLCanvas::STENCIL_BUFFER) {
            attribs[i++] = 32;
        } else {
            attribs[i++] = 0;
        }
        
        // wet set accum buffer to 32 bit if one is required
        attribs[i++] = NSOpenGLPFAAccumSize;
        if (buffers_ & GLCanvas::ACCUM_BUFFER) {
            attribs[i++] = 32;
        } else {
            attribs[i++] = 0;
        }
        
        if (buffers_ & GLCanvas::MULTISAMPLING) {
            attribs[i++] = NSOpenGLPFAMultisample;
        }
        
        // list must be null terminated
        attribs[i] = 0;
        
        assert(i < max);
        
        return [[NSOpenGLPixelFormat alloc] initWithAttributes:attribs];
    }
}