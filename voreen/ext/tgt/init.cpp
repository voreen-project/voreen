/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2020 University of Muenster, Germany,           *
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

#include "tgt/init.h"

#include "tgt/tgt_gl.h"

#include "tgt/assert.h"
#include "tgt/cpucapabilities.h"
#include "tgt/singleton.h"
#include "tgt/glcontextmanager.h"
#include "tgt/gpucapabilities.h"
#include "tgt/matrixstack.h"
#include "tgt/tesselator.h"
#include "tgt/texturemanager.h"
#include "tgt/shadermanager.h"
#include "tgt/fontmanager.h"
#include "tgt/immediatemode/immediatemode.h"
#include "tgt/event/eventhandler.h"


#ifdef TGT_HAS_DEVIL
#include "tgt/texturereaderdevil.h"
#endif
#include "tgt/texturereadertga.h"

namespace tgt {

    static bool initedGL = false;

void init(InitFeature::Features featureset, LogLevel /*logLevel*/) {
    if (featureset & InitFeature::SHADER_MANAGER) {
        featureset = (InitFeature::Features) (featureset | InitFeature::GPU_PROPERTIES | InitFeature::FILE_SYSTEM);
    }

    if (featureset & InitFeature::TEXTURE_MANAGER) {
        featureset = (InitFeature::Features) (featureset | InitFeature::GPU_PROPERTIES | InitFeature::FILE_SYSTEM);
    }

    if (featureset & InitFeature::LOG_MANAGER) {
        LogManager::init();
    }

    if (featureset & InitFeature::FILE_SYSTEM)
        FileSystem::init();

    CPUCapabilities::init();

    // Despite the fact that the GLContextManager is for the use of OpenGL, only,
    // it must be initializied right before(!) calling initGL(), because
    // the initialization needs an active context which can't be registered
    // without having GLContextManager initialized.
    GLContextManager::init();
}

void deinit() {
    if (GLContextManager::isInited())
        GLContextManager::deinit();

    if (CPUCapabilities::isInited())
        CPUCapabilities::deinit();

    if (FileSystem::isInited())
        FileSystem::deinit();

    if (LogManager::isInited())
        LogManager::deinit();
}

void initGL(InitFeature::Features featureset) {
    tgtAssert(!initedGL, "tgt GL already initialized");

    if (featureset & InitFeature::SHADER_MANAGER) {
        featureset = (InitFeature::Features) (featureset | InitFeature::GPU_PROPERTIES | InitFeature::FILE_SYSTEM);
    }
    if (featureset & InitFeature::TEXTURE_MANAGER) {
        featureset = (InitFeature::Features) (featureset | InitFeature::GPU_PROPERTIES | InitFeature::FILE_SYSTEM);
    }

#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE
    GLenum err = glewInit();
    if (err != GLEW_OK) {
        // Problem: glewInit failed, something is seriously wrong.
        tgtAssert(false, "glewInit failed");
        std::cerr << "glewInit failed, error: " << glewGetErrorString(err) << std::endl;
        exit(EXIT_FAILURE);
    }
    LINFOC("tgt.init", "GLEW version:       " << glewGetString(GLEW_VERSION));
#else
    // init the glLoadGen
    ogl_LoadFunctions();
#endif

    if (featureset & InitFeature::GPU_PROPERTIES )
        GpuCapabilities::init();
#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE
    if (featureset & InitFeature::TESSELATOR)
        Tesselator::init();
#endif
    if (featureset & InitFeature::TEXTURE_MANAGER) {
        TextureManager::init();
#ifdef TGT_HAS_DEVIL
            TexMgr.registerReader(new TextureReaderDevil());
        //devil has tga support so we do not need the built-in reader:
#else
            TexMgr.registerReader(new TextureReaderTga());
#endif
    }

    ShaderManager::init();
    MatrixStack::init();
    ImmediateMode::init();
    FontManager::init();

    // Initialization was successful.
    initedGL = true;
}

void deinitGL() {
    tgtAssert(initedGL, "tgt GL not yet initialized");

    if (FontManager::isInited())
        FontManager::deinit();
    if (GpuCapabilities::isInited())
        GpuCapabilities::deinit();
    if (ImmediateMode::isInited())
        ImmediateMode::deinit();
    if (MatrixStack::isInited())
        MatrixStack::deinit();
    if (ShaderManager::isInited())
        ShaderManager::deinit();
    if (TextureManager::isInited())
        TextureManager::deinit();
#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE
    if (Tesselator::isInited())
        Tesselator::deinit();
#endif

    // Deinitialization was successful.
    initedGL = false;
}

bool isInitedGL() {
    return initedGL;
}

} // namespace
