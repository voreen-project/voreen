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

#include "flowindicatorfluid.h"


namespace voreen {

const std::string FlowIndicatorFluid::loggerCat_("voreen.flowsimulation.FlowIndicatorFluid");

FlowIndicatorFluid::FlowIndicatorFluid()
    : Processor()
    , configInport_(Port::INPORT, "flowConfig.inport", "Flow Config Input")
    , configOutport_(Port::OUTPORT, "flowConfig.outport", "Flow Config Output")
    , startDuration_("startDuration", "Start Duration", 0.0f, 0.0f, 10.0f)
{
    addPort(configInport_);
    addPort(configOutport_);

    addProperty(startDuration_);
}

void FlowIndicatorFluid::process() {

    // Here, we just set the output according to our currently set indicators.
    FlowSimulationConfig* config = new FlowSimulationConfig(*configInport_.getData());

    FlowIndicator indicator;
    indicator.type_ = FIT_VELOCITY;
    indicator.flowProfile_ = FP_VOLUME;
    indicator.id_ = MAT_FLUID;
    indicator.velocityCurve_ = VelocityCurve::createSinusoidalCurve(startDuration_.get(), 1.0f);

    config->addFlowIndicator(indicator, true);

    configOutport_.setData(config);
}

}   // namespace
