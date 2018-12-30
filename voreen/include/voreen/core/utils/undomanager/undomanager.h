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

#ifndef VRN_UNDOMANAGER_H
#define VRN_UNDOMANAGER_H

#include "voreen/core/voreencoreapi.h"

#include "tgt/singleton.h"

#include <list>

#include <boost/shared_ptr.hpp>

namespace voreen {

class VRN_CORE_API UndoElement {
public:
    UndoElement() {};
    virtual ~UndoElement() {};

    virtual void undo() {};
};

template<class P, void* V>
class VRN_CORE_API UndoPropertyElement : public UndoElement {
public:
    UndoPropertyElement(P* prop, void* value) : property_(prop), value_(value) {};
    ~UndoPropertyElement() {};

    virtual void undo() {property_->set(V);}
private:
    P* property_;
    void* value_;
};


//-------------------------------------------------------------------------------------------------
class UndoManager;
#ifdef DLL_TEMPLATE_INST
template class TGT_API tgt::Singleton<UndoManager>;
#endif

class VRN_CORE_API UndoManager : public tgt::Singleton<UndoManager> {
public:
    UndoManager() {};
    virtual ~UndoManager() {};

    void addUndoElement(UndoElement elem);

    void undo();
    void redo();

private:
    std::list<UndoElement> undoList_;
    std::list<UndoElement> redoList_;

    static const int MAX_UNDOS = 10;
};

} //namespace

#endif // VRN_UNDOMANAGER_H
