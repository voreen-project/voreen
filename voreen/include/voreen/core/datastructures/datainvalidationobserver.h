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

#pragma once

#include "voreen/core/utils/observer.h"
namespace voreen {

class DataInvalidationObservable;
/**
 * An Observer for DataInvalidationObservable.
 *
 * @see DataInvalidationObservable for details.
 */
class VRN_CORE_API DataInvalidationObserver : public Observer {
public:
    virtual void dataAboutToInvalidate(const DataInvalidationObservable* data) = 0;
};


#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API Observable<DataInvalidationObserver>;
#endif

/*
 * A superclass for a type that notifies its observers if any of the ressources that
 * it shares are not safe to use anymore. This includes (but may not be limited to),
 * for example, deletion of VolumeRepresentations, deletion of the object itself, or
 * invalidation of exposed collections (e.g., due to reallocation).
 *
 * This means that a (const) pointer to  DataInvalidationObservable object is safe to
 * use in another thread, as long as the thread is joined/terminated before
 * notifyPendingDataInvalidation() finishes execution. The usual approach is to
 * implement DataInvalidationObserver and .join() the thread in dataAboutToInvalidate().
 *
 * @see DataInvalidationObserver
 * @see AsyncComputeProcessor for example usage
 * @see Volume for an example implementation
 */
class VRN_CORE_API DataInvalidationObservable : public Observable<DataInvalidationObserver> {
public:
    void notifyPendingDataInvalidation();
};

}
