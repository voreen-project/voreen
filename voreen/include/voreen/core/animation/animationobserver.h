/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#ifndef VRN_ANIMATIONOBSERVER_H
#define VRN_ANIMATIONOBSERVER_H

#include "voreen/core/utils/observer.h"

namespace voreen {

class AnimatedProcessor;

/**
 * Abstract class for an observer of the animationclass.
 * Changes given to the observer are added and removed animatedprocessors.
 */
class VRN_CORE_API AnimationObserver : public Observer {
public:

    /**
     * This function is called if there was an undo or a redo
     * -> all current animatedProcessors / propertytimelines are not valid anymore and have to be reloaded.
     * This function is not called after creation.
     */
    virtual void reloadAnimation() = 0;

    /**
     * This function is called if a new AnimatedProcessor is added.
     */
    virtual void animatedProcessorAdded(const AnimatedProcessor* processor) = 0;

    /**
     * This function is called if an AnimatedProcessor is deleted.
     */
    virtual void animatedProcessorRemoved(const AnimatedProcessor* processor) = 0;
};

} // namespace voreen

#endif
