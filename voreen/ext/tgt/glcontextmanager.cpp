/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2021 University of Muenster, Germany,           *
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

#include "glcontextmanager.h"
#include "glcontextbase.h"
#include "framebufferobject.h"
#include "logmanager.h"

namespace tgt {

const std::string GLContextManager::loggerCat_ = "tgt.GLContextManager";

GLContextManager::GLContextManager() 
    : mainContext_(nullptr)
{
}

GLContextManager::~GLContextManager() {
    tgtAssert(contexts_.empty(), "Context Stack not empty");
    tgtAssert(registered_.empty(), "Some contexts still registered");
}

GLContextBase* GLContextManager::getMainContext() const {
    return mainContext_;
}

void GLContextManager::registerContext(GLContextBase* context, const std::string& title) {
    tgtAssert(context, "Context was null");
    tgtAssert(!isRegisteredContext(context), "Context already registered");

    if (registered_.empty()) {
        tgtAssert(!mainContext_, "Main context needs to be registered first and only once.");
        mainContext_ = context;
        contexts_.push_back(mainContext_); // Push directly onto stack.
    }

    registered_.insert(std::make_pair(context, title));
#ifdef TGT_DEBUG
    LDEBUG("Registered Context: " << registered_[context]);
#endif
}

void GLContextManager::deregisterContext(GLContextBase* context) {
    tgtAssert(context, "Context was null");
    tgtAssert(registered_.find(context) != registered_.end(), "Context not registered");

#ifdef TGT_DEBUG
    LDEBUG("Deregistered Context: " << registered_[context]);
#endif

    registered_.erase(registered_.find(context));
    for (auto it = contexts_.begin(); it != contexts_.end();) {
        if (*it == context)
            it = contexts_.erase(it);
        else
            it++;
    }

    // Restore to valid state.
    if (!contexts_.empty() && !contexts_.back()->isActive()) {
        contexts_.back()->activate();
#ifdef TGT_DEBUG
        LDEBUG("Restored Context: " << registered_[contexts_.back()]);
#endif
    }

    if (registered_.empty())
        mainContext_ = nullptr;
}

bool GLContextManager::isRegisteredContext(GLContextBase* context) const {
    return registered_.find(context) != registered_.end();
}

void GLContextManager::pushContext(GLContextBase* context) {
    tgtAssert(context, "Context was null");
    tgtAssert(isRegisteredContext(context), "Context was not registered");

//#ifdef TGT_DEBUG
//    std::string indentation;
//    for(size_t i=0; i < contexts_.size(); i++) indentation += "* ";
//    LDEBUG(indentation << "Pushing Context: " << registered_[context]);
//#endif

    if (!context->isActive())
        context->activate();

    tgtAssert(context->isActive(), "Activating context failed!");

    contexts_.push_back(context);

#ifdef TGT_DEBUG
    // Enable Debug log output if available.
    if (GpuCaps.isExtensionSupported("KHR_debug")) {
        glEnable(GL_DEBUG_OUTPUT);
    }
#endif
}

void GLContextManager::popContext() {
    tgtAssert(!contexts_.empty(), "Context stack empty");
    contexts_.pop_back();

//#ifdef TGT_DEBUG
//    std::string indentation;
//    for(size_t i=0; i< contexts_.size(); i++) indentation += "* ";
//    LDEBUG(indentation << "Poping Context");
//#endif

    if (!contexts_.empty() && !contexts_.back()->isActive()) {
        contexts_.back()->activate();
        tgtAssert(contexts_.back()->isActive(), "Activating context failed!");
    }
}

GLContextBase* GLContextManager::getActiveContext() const {
    if (contexts_.empty())
        return nullptr;
    return contexts_.back();
}

bool GLContextManager::hasActiveContext() const { 
    return getActiveContext() != nullptr && getActiveContext()->isActive(); 
}

bool GLContextManager::isContextActive(GLContextBase* context) const { 
    return getActiveContext() == context && context->isActive();
}

bool GLContextManager::isMainContextActive() const {
    return isContextActive(mainContext_);
}

void GLContextManager::activateMainContext() {
    mainContext_->activate();
}

const std::string& GLContextManager::getDebugName(GLContextBase* context) const {
    tgtAssert(isRegisteredContext(context), "Context was not registered");
    return registered_.at(context);
}

///////////////////////////////////////////////////////////
//
// ContextStateGuard
//
///////////////////////////////////////////////////////////

const int GLContextStateGuard::SYNC_PASS_TIMEOUT = 1000; // in milliseconds

GLContextStateGuard::GLContextStateGuard(GLContextBase* context)
    : context_(context)
{
    // Flush rendering pipeline if context is about to be changed.
    if(context != GLContextMgr.getActiveContext()) {

        // Flush rendering pipeline without stalling, if possible.
#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE
        glFlush();
#else
        GLsync sync = glFenceSync(GL_SYNC_GPU_COMMANDS_COMPLETE, 0);
        glClientWaitSync(sync, GL_SYNC_FLUSH_COMMANDS_BIT, SYNC_PASS_TIMEOUT);
        glDeleteSync(sync);
#endif
    }

    // Switch contexts.
    GLContextMgr.pushContext(context);

#ifdef TGT_DEBUG
    // Check FBO completness.
    FramebufferObject::isComplete();
#endif
}

GLContextStateGuard::~GLContextStateGuard() {
    // Dont't pop stack if context was deleted within guarded block.
    if (GLContextMgr.isRegisteredContext(context_))
        GLContextMgr.popContext();
}

GLConditionalContextStateGuard::GLConditionalContextStateGuard(bool condition, GLContextBase* context) 
    : guard_(nullptr)
{
    if (condition)
        guard_.reset(new GLContextStateGuard(context ? context : GLContextMgr.getMainContext()));
}

} // namespace
