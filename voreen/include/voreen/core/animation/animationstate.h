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

#ifndef VRN_ANIMATIONSTATE_H
#define VRN_ANIMATIONSTATE_H

#include <vector>
#include <string>
#include "voreen/core/animation/animatedprocessor.h"


namespace voreen {

class NetworkEvaluator;
class AnimationObserver;

/**
 * class to organize all animationsettings to a certain processornetwork
 */
class AnimationState {

public:

    /**
    * creator gets the pointer to the current NetworkEvaluator
    */
    AnimationState(NetworkEvaluator* network);
    /**
    * constructor
    */
    AnimationState();
    /**
    * constructor to load an animation from a file
    */
    AnimationState(std::string file);
    /**
    * destructor
    */
    ~AnimationState();
    /**
    * this function saves all animation settings to a file
    */
    void saveAnimation(std::string file);
    /**
    * returns a vector of all AnimatedProcessors
    */
    const std::vector<AnimatedProcessor*> getAnimatedProcessors() const;
    /**
    * renders the animation to a specific time
    */
    void renderAt(float time);
    /**
    * reruns a clone of this object
    */
    AnimationState* clone();
    /**
    * adds a new AnimatedProcessor to the given Processor
    */
    void addAnimatedProcessor(const Processor* processor);
    /**
    * remove the AnimatedProcessor to the corresponding Processor
    */
    void removeAnimatedProcessor(const Processor* processor);

    /**
    * adds an AnimationObserver
    */
    void addObserver(AnimationObserver* observer);

    /**
    * removes an AnimationObserver
    */
    void removeObserver(AnimationObserver* observer);

    /**
    * returns a vector of all Observers
    */
    std::vector<AnimationObserver*> getObservers();


private:
    std::vector<AnimationObserver*> observer_;

    NetworkEvaluator* network_;
    std::vector<AnimatedProcessor*> processors_;
};

} // namespace voreen

#endif
