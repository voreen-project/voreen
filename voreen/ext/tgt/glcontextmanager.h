/**********************************************************************
 *                                                                    *
 * tgt - Tiny Graphics Toolbox                                        *
 *                                                                    *
 * Copyright (C) 2005-2018 Visualization and Computer Graphics Group, *
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

#ifndef TGT_GLCONTEXTMANAGER_H
#define TGT_GLCONTEXTMANAGER_H

#include "tgt/singleton.h"

#include <deque>
#include <unordered_map>
#include <memory>

namespace tgt {
    class GLContextBase;
    class GLContextManager;
    #ifdef DLL_TEMPLATE_INST
        template class TGT_API Singleton<GLContextManager>;
    #endif

    /**
     * Manager used to handle the different GLContexts in Voreen. It uses the GLContextBases, which are implemented as
     * GLUT or Qt version.
     *
     * As per contract, the first created GLContext object being created is the main context
     * which is used for the main rendering and resource management.
     *
     * @Note The manager takes NO OWNERSHIP of any context.
     */
class TGT_API GLContextManager : public Singleton<GLContextManager>{
    friend class Singleton<GLContextManager>;
public:
    // Convenient
    bool hasActiveContext() const;
    bool isContextActive(GLContextBase* context) const;
    bool isMainContextActive() const;
    void activateMainContext();
    GLContextBase* getMainContext() const;
    const std::string& getDebugName(GLContextBase* context) const;

    // Context registration
    void registerContext(GLContextBase* context, const std::string& title);
    void deregisterContext(GLContextBase* context);
    bool isRegisteredContext(GLContextBase* context) const;

    // Context activation.
    void pushContext(GLContextBase* context);
    void popContext();
    GLContextBase* getActiveContext() const;

private:
    /** Constructor */
    GLContextManager();
    //** Destructor */
    virtual ~GLContextManager();

private:
    GLContextBase* mainContext_;
    std::deque<GLContextBase*> contexts_;
    std::unordered_map<GLContextBase*, std::string> registered_;

    static const std::string loggerCat_;
};

/**
* This guard ensures that the given context is active when creating an instance.
* Nested instanciations are also possible.
* After the instance is destroyed, the former context will be activated, in case
* it hasn't been removed during lifetime.
*
* In case no context is passed, the Main context will be taken.
* @note: Also flushes the pipeline within the old context without stalling if possible.
*/
class TGT_API GLContextStateGuard {
    static const int SYNC_PASS_TIMEOUT;
public:
    GLContextStateGuard(GLContextBase* context = GLContextManager::getRef().getMainContext());
    ~GLContextStateGuard();
private:
    GLContextBase* context_;
};

/**
 * In case the guard is bound to a condition, you should use GLConditionalContextStateGuard.
 * This also can be useful when a guard needs to be placed in an environment,
 * where its possible that no GLContextManager is available (e.g. in tests).
 */
class TGT_API GLConditionalContextStateGuard {
public:
    GLConditionalContextStateGuard(bool condition, GLContextBase* context = nullptr);
private:
    std::unique_ptr<GLContextStateGuard> guard_;
};

}

#define GLContextMgr tgt::GLContextManager::getRef()

#endif //TGT_GLCONTEXTMANAGER_H
