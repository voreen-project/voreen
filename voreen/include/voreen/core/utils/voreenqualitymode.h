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

#ifndef VRN_VOREENQUALITYMODE_H
#define VRN_VOREENQUALITYMODE_H

#include "voreen/core/voreencoreapi.h"
#include "voreen/core/utils/observer.h"

#include "tgt/singleton.h"

#include <set>

namespace voreen {

class QualityModeObserver;

/**
 * This Singleton class stores the quality settings of the rendering.
 * Every PropertyOwner can access this class to get or set the current rendering mode.
 */
class VRN_CORE_API VoreenQualityMode : public tgt::Singleton<VoreenQualityMode>, public Observable<QualityModeObserver> {
public:
    /** Determines the rendering quality. */
    enum RenderingQuality {
        RQ_INTERACTIVE = 0,     ///< Should be used to achive interactive frame rates by lower quality
        RQ_DEFAULT =1,          ///< The default rendering state
        RQ_HIGH = 2             ///< High quality with slow frame rates. Used for animation.
    };

    /** Constructor. */
    VoreenQualityMode();
    /** Destructor. */
    ~VoreenQualityMode();

    /** Returns the current quality. */
    RenderingQuality getQuality() const;
    /** true, if the current quality is RQ_INTERACTIVE. */
    bool isInteractionMode() const;
    /** true, if the current quality is RQ_DEFAULT. */
    bool isDefaultMode() const;
    /** true, if the current quality is RQ_HIGH. */
    bool isHighQualityMode() const;

    /**
     * Requests a new rendering mode. If the current mode is default, the mode will be switched.
     * If all sources have changed to RQ_DEFAULT, the quality switches back to default.
     */
    void requestQualityMode(RenderingQuality requestedQuality, void* source);

protected:

friend class AnimationExportWidget;
friend class CanvasRenderer;

    /**
     * Same as requestQualityMode(...) but always changes to the requested mode.
     *
     * @see requestQualityMode
     */
    void requestQualityModeForced(RenderingQuality requestedQuality, void* source);

    RenderingQuality quality_;      ///< current quality
    std::set<void*> activeSources_; ///< set of sources requesting the current quality mode


private:
    /** Sets the quality and calls notifyQualityChanged(). */
    void setQuality(RenderingQuality quality);
    /** Notifies all registered observers. */
    void notifyQualiyModeChanged();

};

/**
 * Observer interface for classes which need to react on quality changes.
 */
class VRN_CORE_API QualityModeObserver : public Observer {
    //to access qualityModeChanged() function
    friend class VoreenQualityMode;
public:
    QualityModeObserver() : processedInInteraction_(false) {}
protected:
    /** Called on mode change */
    virtual void qualityModeChanged() = 0;
    /** Used in renderer processors to trigger invalidation */
    bool processedInInteraction_;
};

/** Macro foe easy access to the singleton class. */
#define QualityMode tgt::Singleton<VoreenQualityMode>::getRef()

} // namespace voreen

#endif // VRN_VOREENQUALITYMODE_H
