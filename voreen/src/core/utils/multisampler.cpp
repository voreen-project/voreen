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

#include "voreen/core/utils/multisampler.h"

namespace voreen {

Multisampler::Multisampler(RenderPort& target, size_t numSamples)
    : tex_(0)
    , fbo_(0)
    , port_(target)
{
    glGenTextures(1, &tex_);
    glBindTexture(GL_TEXTURE_2D_MULTISAMPLE, tex_);
    size_t width = port_.getSize().x;
    size_t height = port_.getSize().y;
    glTexImage2DMultisample(GL_TEXTURE_2D_MULTISAMPLE, numSamples, port_.getColorTexture()->getGLFormat(), width, height, false);

    glGenFramebuffers(1, &fbo_);
    glBindFramebuffer(GL_FRAMEBUFFER, fbo_);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D_MULTISAMPLE, tex_, 0);

    glViewport(0, 0, width, height);

    glBindFramebuffer(GL_FRAMEBUFFER, fbo_);
}

Multisampler::~Multisampler() {
    size_t width = port_.getSize().x;
    size_t height = port_.getSize().y;

    port_.activateTarget();
    glBindFramebuffer(GL_READ_FRAMEBUFFER, fbo_); // Make sure your multisampled FBO is the read framebuffer
    //glDrawBuffer(GL_BACK);                        // Set the back buffer as the draw buffer
    glBlitFramebuffer(0, 0, width, height, 0, 0, width, height, GL_COLOR_BUFFER_BIT, GL_NEAREST);
    port_.deactivateTarget();

    glDeleteFramebuffers(1, &fbo_);
    glDeleteTextures(1, &tex_);
}

}
