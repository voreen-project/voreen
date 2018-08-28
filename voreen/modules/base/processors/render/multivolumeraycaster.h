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

#ifndef VRN_MULTIVOLUMERAYCASTER_H
#define VRN_MULTIVOLUMERAYCASTER_H

#include "voreen/core/processors/volumeraycaster.h"

#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/shaderproperty.h"

#include "voreen/core/ports/volumeport.h"

namespace voreen {

/**
 * Performs a single-pass multi-volume raycasting of up to four volumes,
 * providing several shading and compositing modes.
 *
 * @note The MeshEntryExitPoints must be set to WorldCoordinates
 *
 * @see MultiVolumeProxyGeometry, MeshEntryExitPoints
 */
class VRN_CORE_API MultiVolumeRaycaster : public VolumeRaycaster {
public:
    MultiVolumeRaycaster();
    Processor* create() const;

    std::string getClassName() const    { return "MultiVolumeRaycaster"; }
    std::string getCategory() const     { return "Raycasting"; }
    CodeState getCodeState() const      { return CODE_STATE_STABLE; }

    /**
     * One outport, one volume with associated entry and exit points and global
     * entry and exit points must be connected.
     */
    bool isReady() const;

protected:
    virtual void setDescriptions() {
        setDescription("Performs a single-pass multi-volume raycasting of up to four volumes, providing several shading and compositing modes. <br>" \
                       "It allows one to generate three output renderings, whereas only the first one provides depth values. <br>" \
                       "To use this processor, at least one volume, the associated entry and exit points, and global entry and exit points" \
                       "must be connected. <br>" \
                       "Use the template workspace <i>multivolumeraycasting.vws</i> to get an idea of how to use the raycaster. <br>" \
                       "<b>The MeshEntryExitPoints must be set to World-Coordiantes</b>" \
                       "<p>See MultiVolumeProxyGeometry, MeshEntryExitPoints.</p>");
    }

    void process();
    /// Compile shader etc.
    void beforeProcess();

    void initialize();
    void deinitialize();

    std::string generateHeader(const tgt::GpuCapabilities::GlVersion* version = 0);
    void compile();
private:
    void adjustPropertyVisibilities();

    VolumePort volumeInport2_;
    RenderPort entryPort2_;
    RenderPort exitPort2_;
    VolumePort volumeInport3_;
    RenderPort entryPort3_;
    RenderPort exitPort3_;
    VolumePort volumeInport4_;
    RenderPort entryPort4_;
    RenderPort exitPort4_;
    RenderPort entryPortGlobal_;
    RenderPort exitPortGlobal_;


    ShaderProperty shaderProp_;        ///< The shader property used by this raycaster.

    StringOptionProperty classificationMode2_;  ///< What type of transfer function should be used for classification
    StringOptionProperty classificationMode3_;  ///< What type of transfer function should be used for classification
    StringOptionProperty classificationMode4_;  ///< What type of transfer function should be used for classification

    StringOptionProperty shadeMode1_;       ///< shading mode to use for volume 1
    StringOptionProperty shadeMode2_;       ///< shading mode to use for volume 2
    StringOptionProperty shadeMode3_;       ///< shading mode to use for volume 3
    StringOptionProperty shadeMode4_;       ///< shading mode to use for volume 4

    TransFunc1DKeysProperty transferFunc1_;       ///< transfer function to apply to volume 1
    TransFunc1DKeysProperty transferFunc2_;       ///< transfer function to apply to volume 2
    TransFunc1DKeysProperty transferFunc3_;       ///< transfer function to apply to volume 3
    TransFunc1DKeysProperty transferFunc4_;       ///< transfer function to apply to volume 4

    CameraProperty camera_;                 ///< the camera used for lighting calculations

};


} // namespace voreen

#endif // VRN_MULTIVOLUMERAYCASTER_H
