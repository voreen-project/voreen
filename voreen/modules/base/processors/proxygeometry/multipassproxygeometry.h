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

#ifndef VRN_MULTIPASSPROXYGEOMETRY_H
#define VRN_MULTIPASSPROXYGEOMETRY_H

#include "voreen/core/processors/renderprocessor.h"
#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/renderport.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/ports/loopport.h"
#include "voreen/core/interaction/idmanager.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/cameraproperty.h"

#include <stack>

namespace voreen {

/**
 * Provides a ray casting proxy geometry that can be sectioned into bricks for exploded views.
 * The processor's loop port has to be connected to a ExplosionCompositor.
 *
 * Three bricking modes are offered:
 *
 * Regular Grid:
 * In this mode the proxy geometry of the volume can be split into a regular grid by choosing the number
 * of bricks on each axis. The bricks can "explode" by using the gap property of this mode.
 *
 * Manual Grid:
 * In this mode the proxy geometry of the volume can be split into a manual grid by manually adding
 * axis-aligned planes as boundaries for the bricks. As in regularGrid mode bricks can "explode" by using
 * the gap property of this mode.
 *
 * Custom:
 * In addition to manualGrid mode one or several bricks can be selected in this mode, so that added
 * axis-aligned planes only split the selected bricks and not the whole volume.
 * The "explosion" of the bricks no longer depends on a gap-value. Instead each brick can be translated everywhere
 * by using the translation property or the drag-and-drop translation mouse event of this mode.
 * Furthermore bricks can be hidden in this mode.
 *
 * The created proxy geometries of regular grid or manual grid mode can be used for custom mode by
 * using the "Use this for custom mode" property.
 *
 * @see ExplosionCompositor
 * @see MeshEntryExitPoints
 */
class VRN_CORE_API MultiPassProxyGeometry : public Processor {
public:
    MultiPassProxyGeometry();
    virtual ~MultiPassProxyGeometry();
    virtual Processor* create() const;

    virtual std::string getClassName() const    { return "MultiPassProxyGeometry"; }
    virtual std::string getCategory() const     { return "Volume Proxy Geometry";  }
    virtual CodeState getCodeState() const      { return CODE_STATE_TESTING;        }

protected:
    virtual void setDescriptions() {
        setDescription("<p>See ExplosionCompositor, MeshEntryExitPoints</p>");
    }

    virtual void beforeProcess();
    virtual void process();
    virtual void afterProcess();

    void performDepthSort();

    CameraProperty camera_;
    IntProperty test_;

    /// To be connected to the corresponding ExplosionCompositor's loop port.
    LoopPort loopInport_;

    GeometryPort inportProxyGeometry_;
    GeometryPort outportProxyGeometry_;

    std::vector<const Geometry*> sortedInput_;

private:
    /// Change listener of the canvasSize_ property. Calls RenderProcessor::portResized
    void canvasSizeChanged();
};

} // namespace

#endif // VRN_MULTIPASSPPROXYGEOMETRY_H
