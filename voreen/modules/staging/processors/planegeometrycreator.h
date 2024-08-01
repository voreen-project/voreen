/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#ifndef VRN_PLANEGEOMETRYCREATOR_H
#define VRN_PLANEGEOMETRYCREATOR_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/buttonproperty.h"

#include "voreen/core/ports/volumeport.h"
#include "voreen/core/ports/geometryport.h"

namespace voreen {

/**
 * Generates a plane geometry defined by the plane normal and the position in world space.
 * Best practice to adjust the plane settings is to link the plane normal and position to the <i>PlanManipulation</i> processor
 * for an easy handling. <br>
 * A common application case is to use this plane as input for the <i>MeshEntryExitPoints</i> processor and connect the entry/exit
 * texture with the <i>MultiSliceRenderer</i> to redner an arbitraty slice of the volume.
 * TODO: correkt handling of quad orientation
 */
class PlaneGeometryCreator : public Processor {
public:
    PlaneGeometryCreator();
    virtual ~PlaneGeometryCreator();

    virtual std::string getCategory() const override{ return "Volume"; }
    virtual std::string getClassName() const override{ return "PlaneGeometryCreator"; }
    virtual Processor::CodeState getCodeState() const override{ return CODE_STATE_TESTING; }
    virtual Processor* create() const override{ return new PlaneGeometryCreator(); }
    /** @see processor */
    virtual void process() override;
    /** VolumePort is optinal. */
    virtual bool isReady() const override;
protected:
    virtual void setDescriptions() {
        setDescription("Generates a plane geometry defined by the plane normal and the position in world space.. \
                        Best practice to adjust the plane settings is to link the plane normal and position to the <i>PlanManipulation</i> processor for an \
                        easy handling. Therefor link the plane normal with the plane normal of the <i>PlanManipulation</i> processor and the plane position with \
                        the center in world space of the <i>PlanManipulation</i> processor. <br> \
                        A common application case is to use this plane as input for the <i>MeshEntryExitPoints</i> processor and connect the entry texture with the \
                        <i>MultiSliceRenderer</i> to redner an arbitraty slice of the volume. Note: In the <i>MeshEntryExitPoints</i> processor \"Use culling\" and \
                        \"Camera inside volume technique\" must be disabled and the output coordinate system must be \"World Coordinates\". <br> \
                        <b>Important:</b>If no volume is connected, the maximal position and size is Float_Max which can not be adjusted by the GUI-Sliders. \
                        Use the textfiels instead.");
        volInport_.setDescription("This port is optional. If a volume is connected and the property option is toggled, the max range of the plane position \
                                   and size is adjusted to the input volumes dimensions.");
        geomOutport_.setDescription("Output of the created plane geometry.");
        adaptFromVolume_.setDescription("If true and a volume is connected, the max range of the plane position and size is adjusted to the input volumes dimensions.");
        normal_.setDescription("Plane normal. A value of (0,0,0) is not allowed.");
        position_.setDescription("Plane position in world space.");
        size_.setDescription("Size of each side of the created plane.");
        reset_.setDescription("Resets the plane settings to the default values. If the ranges are adapted to the input volumes dimensions, they are reset to the centered values.");
    }

    /**
     * Callback on volInport onChange.
     * Plane properties are adapted to the volume dimensions etc..
     */
    void onVolumePortChangeCallback();

    /**
     * Callback from reset_ property.
     * Resets all plane properties to default.
     */
    void onResetClickedCallback();

    //--------------//
    // Members      //
    //--------------//

    //properties
        //adapt
    BoolProperty adaptFromVolume_;  ///< if true, the plane parameters are adapted from the input volume
        //plane settings
    FloatVec3Property normal_;      ///< Clipping plane normal
    FloatVec3Property position_;    ///< Clipping plane position in world space
    FloatProperty size_;            ///< Size of the generated plane (plane)
    ButtonProperty reset_;          ///< reset values to default values

    //ports
    VolumePort volInport_;          ///< Optional volume inport, to adjust properties to useful default values
    GeometryPort geomOutport_;      ///< Output of the generated slice geometry

    static const std::string loggerCat_;    ///< meow
};

} // namespace voreen

#endif
