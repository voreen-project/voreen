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

#ifndef VRN_SPHERICALVOLUMELISTPROXY_H
#define VRN_SPHERICALVOLUMELIStPROXY_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/numericproperty.h"

#include "voreen/core/ports/volumeport.h"

namespace voreen {

/**
 * This is the spherical VolumeProxy.
 *
 * @see CubeMeshProxyGeometry, MeshEntryExitPoints
 * @see MultiVolumeVolumeProxy
 */
class VRN_CORE_API SphericalVolumeListProxy : public Processor,
        public PortObserver,
        public DataInvalidationObserver {
public:
    SphericalVolumeListProxy();
    virtual ~SphericalVolumeListProxy();

    virtual Processor* create() const;

    virtual std::string getClassName() const { return "SphericalVolumeListProxy"; }

    virtual std::string getCategory() const { return "ProxyVolume"; }

    virtual CodeState getCodeState() const { return CODE_STATE_STABLE; }

    //! @see PortObserver
    virtual void afterConnectionAdded(const Port* source, const Port* connectedPort);
    virtual void beforeConnectionRemoved(const Port* source, const Port*);
    virtual void dataWillChange(const Port* source);
    virtual void dataHasChanged(const Port* source);

    //! @see DataInvalidationObserver
    virtual void dataAboutToInvalidate(const DataInvalidationObservable* data);

protected:
    virtual void setDescriptions() {
        setDescription("This is a proxy volume which gets a spherical volume as in-port and outputs a volume in cartesian coordinates. <p>See CubeProxyGeometry, MeshEntryExitPoints.</p>");
    }

    virtual void process();

private:

    VolumeListPort inport_;
    VolumeListPort outport_;

    IntProperty outputDimensions_;

    FloatProperty radiusMin_;
    FloatProperty radiusMax_;
    FloatProperty clampRadiusMin_;
    FloatProperty clampRadiusMax_;
    FloatProperty shiftFactor_;

    std::vector<std::unique_ptr<VolumeBase>> proxyVolumes_;

    static const std::string loggerCat_; ///< category used in logging
};


} // namespace voreen

#endif // VRN_SPHERICALVOLUMELISTPROXY_H
