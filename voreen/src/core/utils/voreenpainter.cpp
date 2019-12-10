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

#include "voreen/core/utils/voreenpainter.h"

#include "voreen/core/network/networkevaluator.h"
#include "voreen/core/network/processornetwork.h"

#include "modules/core/processors/output/canvasrenderer.h" //< core module is always available

#include "tgt/glcontextmanager.h"

namespace voreen {

const std::string VoreenPainter::loggerCat_ = "voreen.core.VoreenPainer";

VoreenPainter::VoreenPainter(tgt::GLCanvas* canvas, NetworkEvaluator* evaluator, CanvasRenderer* canvasRenderer)
    : tgt::Painter(canvas)
    , evaluator_(evaluator)
    , canvasRenderer_(canvasRenderer)
{
    tgtAssert(canvas, "No canvas");
    tgtAssert(evaluator_, "No network evaluator");
    tgtAssert(canvasRenderer_, "No canvas renderer");
}

VoreenPainter::~VoreenPainter() {
}

void VoreenPainter::initialize() {
}

void VoreenPainter::sizeChanged(const tgt::ivec2& size) {
    tgtAssert(canvasRenderer_, "No canvas renderer");
    canvasRenderer_->resizeCanvas(size);
}

void VoreenPainter::paint() {
    if (!canvas_) {
        LWARNING("No canvas assigned");
        return;
    }

    if (!evaluator_) {
        LWARNING("No network evaluator assigned");
        return;
    }

    evaluator_->process();
}

void VoreenPainter::repaint() {
    canvas_->repaint();
}

void VoreenPainter::renderToSnapshot(const std::string& fileName, const tgt::ivec2& size) {

    tgtAssert(canvasRenderer_, "No canvas renderer");
    canvasRenderer_->renderToImage(fileName, size);

}

NetworkEvaluator* VoreenPainter::getEvaluator() const {
    return evaluator_;
}

CanvasRenderer* VoreenPainter::getCanvasRenderer() const {
    return canvasRenderer_;
}

} // namespace
