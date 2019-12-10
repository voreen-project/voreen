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

#ifndef VRN_MEINPROCESSOR_H
#define VRN_MEINPROCESSOR_H

#include "voreen/core/processors/renderprocessor.h"

#include "voreen/core/ports/volumeport.h"

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/cameraproperty.h"

#include "tgt/shadermanager.h"

namespace voreen {

class CameraInteractionHandler;

/**
 * Calculates the entry and exit points for GPU raycasting and stores them in textures.
 */
class VRN_CORE_API MeinProcessor : public RenderProcessor {
public:
    MeinProcessor();
    virtual ~MeinProcessor();
    virtual Processor* create() const;

    virtual std::string getClassName() const    { return "MeinProcessor"; }
    virtual std::string getCategory() const     { return "Schülerakademie";   }
    virtual CodeState getCodeState() const      { return CODE_STATE_STABLE; }

    virtual bool isReady() const;

    virtual void adjustRenderOutportSizes();

protected:
    virtual void setDescriptions() {
        setDescription("Processor, der von den Schülern der Schülerakademie 2014 designed wurde.");
    }

    virtual void beforeProcess();
    virtual void process();

    void renderGeometry(const VolumeBase* volume, RenderPort& outport, GLenum depthFunc, float clearDepth, GLenum cullFace);
    void zeichneWuerfel(tgt::vec3 llf, tgt::vec3 urb);

    // ports
    VolumePort inport_;
    RenderPort entryPort_;
    RenderPort exitPort_;

    // properties
    CameraProperty camera_;  ///< camera used for rendering the proxy geometry

    // interaction handlers
    CameraInteractionHandler* cameraHandler_;

    /// Category used for logging.
    static const std::string loggerCat_;
};

} // namespace voreen

#endif //VRN_MESHENTRYEXITPOINTS_H
