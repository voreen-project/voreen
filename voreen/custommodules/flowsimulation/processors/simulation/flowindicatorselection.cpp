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

#include "flowindicatorselection.h"

namespace voreen {

const std::string FlowIndicatorSelection::loggerCat_("voreen.flowreen.FlowIndicatorSelection");

FlowIndicatorSelection::FlowIndicatorSelection()
    : RenderProcessor()
    , renderInport_(Port::INPORT, "render.inport", "Render Inport")
    , renderOutport_(Port::OUTPORT, "render.outport", "Render Outport")
    , flowParametrizationPort_(Port::OUTPORT, "flowParametrization.outport", "Flow Parametrization")
    , ensembleName_("ensembleName", "Ensemble Name", "test_ensemble")
    , simulationTime_("simulationTime", "Simulation Time (s)", 2.0f, 0.1f, 20.0f)
    , temporalResolution_("temporalResolution", "Temporal Resolution (ms)", 3.1f, 1.0f, 200.0f)
    , numTimeSteps_("numTimeSteps", "Num. Output Time Steps", 50, 1, 1000)
    , pickingMatrix_("pickingMatrix", "Picking Matrix (link with SliceViewer)", tgt::mat4::createIdentity(), tgt::mat4(-1e6f), tgt::mat4(1e6f))
    , rebuildOutput_(true)
{
    addPort(renderInport_);
    addPort(renderOutport_);
    addPort(flowParametrizationPort_);

    addProperty(ensembleName_);
        ensembleName_.setGroupID("ensemble");
    addProperty(simulationTime_);
        simulationTime_.setGroupID("ensemble");
    addProperty(temporalResolution_);
        temporalResolution_.setGroupID("ensemble");
    addProperty(numTimeSteps_);
        numTimeSteps_.setGroupID("ensemble");
    setPropertyGroupGuiName("ensemble", "Ensemble");

    addProperty(pickingMatrix_);
}

bool FlowIndicatorSelection::isReady() const {
    // Ignore vessel port!
    return renderInport_.isReady();
}

void FlowIndicatorSelection::process() {

    if(rebuildOutput_) {
        FlowParametrizationList* parametrizationList = new FlowParametrizationList(ensembleName_.get());
        parametrizationList->setSimulationTime(simulationTime_.get());
        parametrizationList->setTemporalResolution(temporalResolution_.get() / 1000.0f); // Convert ms to s
        parametrizationList->setNumTimeSteps(numTimeSteps_.get());
        for(const FlowIndicator& flowIndicator : flowIndicators_) {
            parametrizationList->addFlowIndicator(flowIndicator);
        }

        flowParametrizationPort_.setData(parametrizationList);
        rebuildOutput_ = false;
    }

    renderOutport_.activateTarget();
    renderOutport_.clearTarget();

    for(const FlowIndicator& flowIndicator : flowIndicators_) {
    }

    tgt::TextureUnit colorTex, depthTex;
    renderInport_.bindTextures(colorTex.getEnum(), depthTex.getEnum());

    renderQuad();

    renderOutport_.deactivateTarget();

}

/*
FlowIndicatorSelection::Circle FlowIndicatorSelection::estimateCircle(const tgt::ivec2& texCoord) const {

    tgtAssert(slice_, "Slice null");
    tgtAssert(slice_->getCpuTextureData(), "Slice data null");

    // Use floodfill (4-neighborhood) to estimate radius.
    tgt::IntBounds bounds;
    std::stack<tgt::ivec2> next;
    next.push(texCoord);
    while(!next.empty()) {
        tgt::ivec2 pos = next.top();
        if(pos.x >= 0 && pos.y >= 0 &&
           pos.y < resolution_.y && pos.y < resolution_.y &&
           slice_->texel<float>(next.top()) > 0.0f)
        {
            bounds.addPoint(tgt::vec3(pos.x, pos.y, 0));

            next.push(pos + tgt::ivec2(-1,  0));
            next.push(pos + tgt::ivec2( 1,  0));
            next.push(pos + tgt::ivec2( 0, -1));
            next.push(pos + tgt::ivec2( 0,  1));
        }
        next.pop();
    }

    if(!bounds.isDefined()) {
        // Zero radius determines invalid pick.
        return Circle();
    }

    // DEBUG
    slice_->texel<float>(bounds.getLLF().xy()) = 10.0f;
    slice_->texel<float>(bounds.getURB().xy()) = 10.0f;
    // !DEBUG

    Circle circle;
    circle.radius_ = tgt::distance(screenToVoxelPos(bounds.getLLF().xy()), screenToVoxelPos(bounds.getURB().xy()))/2.0f;
    circle.center_ = screenToVoxelPos(bounds.center().xy());
    circle.normal_ = plane_.n;

    return circle;
}
*/
void FlowIndicatorSelection::selectRegion(tgt::MouseEvent* e) {
    if(e->modifiers() & tgt::MouseEvent::SHIFT) {
        tgt::ivec2 texCoords;
        texCoords.x = resolution_.x * e->x() / e->viewport().x;
        texCoords.y = resolution_.y * e->y() / e->viewport().y;

        //circles_.push_back(estimateCircle(texCoords));
    }
    else {
        //circles_.clear();
    }
}

}   // namespace
