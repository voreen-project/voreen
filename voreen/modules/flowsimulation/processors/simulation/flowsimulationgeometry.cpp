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

#include "flowsimulationgeometry.h"

#include "include/voreen/core/datastructures/geometry/glmeshgeometry.h"

namespace voreen {

const std::string FlowSimulationGeometry::loggerCat_("voreen.flowsimulation.FlowSimulationGeometry");

FlowSimulationGeometry::FlowSimulationGeometry()
    : Processor()
    , flowParametrizationInport_(Port::INPORT, "flowparametrization.inport", "Flow Parametrization Input")
    , flowParametrizationOutport_(Port::OUTPORT, "flowparametrization.outport", "Flow Parametrization Output")
    , geometryPort_(Port::OUTPORT, "geometry", "Geometry Port")
    , geometryType_("geometryType", "Geometry Type")
    , ratio_("ratio", "Ratio", 1.0f, 0.1f, 10.0f)
    , radius_("radius", "Radius (mm)", 1.0f, 0.1f, 1000.0f)
    , length_("length", "Length (mm)", 1.0f, 0.1f, 1000.0f)
    , transformation_("transformation", "Transformation", tgt::mat4::identity, tgt::mat4(-1000.0f), tgt::mat4(+1000.0f))
    , flowProfile_("flowProfile", "Flow Profile")
    , inflowVelocity_("inflowVelocity", "Inflow Velocity (m/s)", 1.0f, 0.0f, 10.0f)
{
    addPort(flowParametrizationInport_);
    addPort(flowParametrizationOutport_);
    addPort(geometryPort_);

    addProperty(geometryType_);
    ON_CHANGE_LAMBDA(geometryType_, [this] {
        ratio_.setReadOnlyFlag(geometryType_.getValue() != FSGT_NARROWING);
    })
    geometryType_.addOption("channel", "Channel", FSGT_CHANNEL);
    geometryType_.addOption("cylinder", "Cylinder", FSGT_CYLINDER);
    geometryType_.addOption("narrowing", "Narrowing", FSGT_NARROWING);

    addProperty(ratio_);
    addProperty(radius_);
    addProperty(length_);
    addProperty(transformation_);
    addProperty(flowProfile_);
    ON_CHANGE_LAMBDA(flowProfile_, [this] {
        inflowVelocity_.setReadOnlyFlag(flowProfile_.getValue() == FlowProfile::FP_VOLUME);
    })
    flowProfile_.addOption("poiseuille", "Poiseuille", FlowProfile::FP_POISEUILLE);
    flowProfile_.addOption("powerlaw", "Powerlaw", FlowProfile::FP_POWERLAW);
    //flowProfile_.addOption("constant", "CONSTANT", FlowProfile::FP_CONSTANT);
    flowProfile_.addOption("volume", "Volume", FlowProfile::FP_VOLUME);
    addProperty(inflowVelocity_);
    inflowVelocity_.setNumDecimals(3);
}

void FlowSimulationGeometry::process() {

    // Add indicators.
    auto* config = new FlowSimulationConfig(*flowParametrizationInport_.getData());

    // Add a ramp up that lasts half of the simulated time.
    float rampUpTime = config->getSimulationTime() * 0.5f;

    FlowIndicator inlet;
    inlet.type_ = FIT_VELOCITY;
    inlet.flowProfile_ = flowProfile_.getValue();
    inlet.velocityCurve_ = VelocityCurve::createSinusoidalCurve(rampUpTime, inflowVelocity_.get());
    inlet.center_ = transformation_.get() * tgt::vec3(0.0f, 0.0f, 0.0f);
    inlet.normal_ = transformation_.get().getRotationalPart() * tgt::vec3(0.0f, 0.0f, 1.0f);
    inlet.radius_ = radius_.get();
    config->addFlowIndicator(inlet);

    FlowIndicator outlet;
    outlet.type_ = FIT_PRESSURE;
    outlet.center_ = transformation_.get() * tgt::vec3(0.0f, 0.0f, length_.get());
    outlet.normal_ = transformation_.get().getRotationalPart() * tgt::vec3(0.0f, 0.0f, 1.0f);
    outlet.radius_ = radius_.get();
    config->addFlowIndicator(outlet);

    flowParametrizationOutport_.setData(config);

    // Generate geometry.
    auto* geometry = new GlMeshGeometryUInt32Normal();
    geometry->setTransformationMatrix(transformation_.get());
    switch(geometryType_.getValue()) {
    case FSGT_CHANNEL:
        geometry->setCylinderGeometry(tgt::vec4::one, radius_.get(), radius_.get(), length_.get(), 4, 32, true, true);
        //geometry->setTransformationMatrix(transformation_.get() * tgt::mat4::createRotationZDegree(45.0f));
        break;
    case FSGT_CYLINDER:
        geometry->setCylinderGeometry(tgt::vec4::one, radius_.get(), radius_.get(), length_.get(), 32, 32, true, true);
        break;
    case FSGT_NARROWING: {
        const size_t slices = 32;
        const size_t tiles = 5;
        geometry->setCylinderGeometry(tgt::vec4::one, radius_.get(), radius_.get(), length_.get(), slices, tiles, true, true);
        for(size_t i=2; i<=3; i++) {
            size_t offset = (2 + i) * (slices + 1); // (cap + tiles)
            for (size_t j = 0; j <= slices; j++) {
                auto vertex = geometry->getVertex(offset + j);
                vertex.pos_.x *= ratio_.get() * radius_.get();
                vertex.pos_.y *= ratio_.get() * radius_.get();
                geometry->setVertex(offset + j, vertex);
            }
        }
        break;
    }
    default:
        tgtAssert(false, "Unhandled geometry type");
        break;
    }

    geometryPort_.setData(geometry);
}

}   // namespace
