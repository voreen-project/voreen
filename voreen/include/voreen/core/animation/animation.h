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

#ifndef VRN_ANIMATION_H
#define VRN_ANIMATION_H

#include <vector>
#include <string>
#include <queue>
#include "voreen/core/network/processornetwork.h"
#include "voreen/core/animation/animationobserver.h"
#include "modules/core/coremodule.h"

namespace voreen {

class AnimationObserver;
class ProcessorNetwork;
class AnimatedProcessor;
class PropertyTimeline;
class SerializableFactory;

#ifdef DLL_TEMPLATE_INST
template class VRN_CORE_API Observable<AnimationObserver>;
#endif

/**
 * Mainclass to create an animation.
 */
class VRN_CORE_API Animation : public Serializable, public ProcessorNetworkObserver, public Observable<AnimationObserver> {
public:
    /**
     * Constructor: creates a new animation for a given processornetwork.
     */
    Animation(ProcessorNetwork* network);

    /**
     * Destructor.
     */
    ~Animation();

    bool isEmpty() const;

    /**
     * This function returns the vector of all animatedprocessors.
     */
    const std::vector<AnimatedProcessor*>& getAnimatedProcessors() const;

    /**
    * This function renders with the settings at the given point of time by calling the 'renderAt'-function of all propertytimelines.
    */
    void renderAt(float time);

    /**
     * This function is called from the propertytimelines to inform of a change:
     * 1.the pointer to this propertytimeline is saved for undo / redo.
     * 2.change made -> delete all redostates
     * 3.save the pointer to the last changed propertytimeline
     * 4.if there are already too many undosteps -> delete oldest states
     */
    void animationChanged(PropertyTimeline* changedObject);

    /**
     * This function undos the last change by calling the undo-function of the last changed propertytimeline.
     */
    void undoLastChange();

    /**
     * This function redos the last undo by calling the redo-function of the corresponding propertytimeline.
     */
    void redoLastUndo();

    /**
     * Gets the frames per second number for video-output.
     */
    float getFPS() const;

    /**
     * Sets the frames per second number for video-output.
     */
    void setFPS(float fps);

    /**
     * Sets the number of steps that are stored for the undo-function.
     */
    void setUndoSteps(int steps);

    /**
     * Sets all settings of the current rendernetwork as keyvalues at the given time.
     *
     * If force is set to false, the keyvalues will only be set if the setting differs from the animated value at the given time.
     * If force is set to true, the keyvalues will always be set.
     */
    void setActualNetworkAsKeyvalues(float time, bool force = false);

    /**
     * Sets the interactionmode of all processors.
     */
    void setInteractionMode(bool interactionmode);

    /**
     * Sets the network (serialization).
     */
    void setNetwork(ProcessorNetwork* network);

    /**
     * This method is called on all modifications of the observed network
     * that do not match one of the more specialized modifications below.
     */
    virtual void networkChanged();

    /**
     * This method is called by the observed network after a processor has been added.
     */
    void processorAdded(const Processor* processor);

    /**
     * This method is called by the observed network immediately before
     * a processor is removed from the network and deleted.
     */
    void processorRemoved(const Processor* processor);

    /**
     * Returns the current time (the last rendered point of time).
     */
    float getCurrentTime() const;

    /**
     * Returns the duration of the animation.
     */
    float getDuration() const;

    /**
     * Sets the duration of the animation.
     * All ando and redosteps will be deleted,
     * all existing keyvalues in a timeline with a higher time will be deleted.
     */
    void setDuration(float duration);

    /**
     * @see Serializable::serialize
     */
    virtual void serialize(Serializer& s) const;

    /**
     * @see Serializable::deserialize
     */
    virtual void deserialize(Deserializer& s);

protected:
    static const std::string loggerCat_; ///< category used in logging

private:
    friend class CoreModule;
    friend class XmlDeserializer;

    /**
     * Default constructor.
     */
    Animation();

    /**
     * Construct all factories used for animation serialization.
     * Is called by the CoreModule.
     */
    static std::vector<SerializableFactory*> getSerializerFactories();

    static void deleteSerializerFactories();

    /**
     * Pointer to the corresponding rendernetwork.
     */
    ProcessorNetwork* network_;

    /**
     * Vector of all the animatedprocessors in the animation.
     */
    std::vector<AnimatedProcessor*> processors_;

    /**
     * Deque of pointers to the propertytimelines which were changed.
     */
    std::deque<PropertyTimeline*> lastChanges_;
    std::deque<PropertyTimeline*> lastUndos_;

    /**
     * Contains the fps-number for video-output.
     */
    float fps_;

    /**
     * Number of possible undo-steps.
     */
    int undoSteps_;

    /**
     * Contains the last rendered time.
     */
    float currentTime_;

    /**
     * Contains the duration of the animation in seconds.
     */
    float duration_;

    /**
     * Value if the function 'renderAt' is currently running
     * (only interesting if the program is multithreaded).
     */
    bool isRendering_;

    static std::vector<SerializableFactory*> factories_;

};

} // namespace voreen

#endif
