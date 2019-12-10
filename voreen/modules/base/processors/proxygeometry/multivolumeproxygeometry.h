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

#ifndef VRN_MULTIVOLUMEPROXYGEOMETRY_H
#define VRN_MULTIVOLUMEPROXYGEOMETRY_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/geometryport.h"

namespace voreen {

/**
 * Provides a cube mesh proxy geometry for multi-volume raycasting.
 *
 * @see MeshEntryExitPoints
 * @see MeshClipping
 * @see MultiVolumeRaycaster
 */
class VRN_CORE_API MultiVolumeProxyGeometry : public Processor {

public:
    MultiVolumeProxyGeometry();
    virtual ~MultiVolumeProxyGeometry();
    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "MultiVolumeProxyGeometry"; }
    virtual std::string getCategory() const   { return "Volume Proxy Geometry"; }
    virtual CodeState getCodeState() const    { return CODE_STATE_STABLE; }

protected:
    virtual void setDescriptions() {
        setDescription("Provides a cube mesh proxy geometry for multi-volume raycasting.\
<p>See MeshEntryExitPoints, MultiVolumeRaycaster.</p>");
    }

    virtual void process();

    /**
     * Inport for the dataset.
     */
    VolumePort inport_;

    /**
     * Outport for the cube mesh proxy geometry.
     */
    GeometryPort outport_;
};

} // namespace

#endif // VRN_MULTIVOLUMEPROXYGEOMETRY_H
