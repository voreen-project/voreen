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

#include "plotpickingmanager.h"
#include "voreen/core/processors/processor.h"

#include "tgt/logmanager.h"
#include "tgt/immediatemode/immediatemode.h"

namespace voreen {

PlotPickingManager::PlotPickingManager() :
      rt_(0)
{}

void PlotPickingManager::initializeTarget() {
    if (!rt_)
        return;
    rt_->initialize(GL_RGB);
}

void PlotPickingManager::setGLColor(int row, int column) {
    tgt::col3 col = getColorFromCell(row, column);
    IMode.color(IMode.normalize(col));
}

tgt::ivec2 PlotPickingManager::getCellAtPos(tgt::ivec2 pos) {
    return getCellFromColor(getColorAtPos(pos));
}

tgt::col3 PlotPickingManager::getColorAtPos(tgt::ivec2 pos) {
    rt_->activateTarget();
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    tgt::col3 pixels;
    glReadPixels(pos.x, pos.y, 1, 1, GL_RGB, GL_UNSIGNED_BYTE, &pixels);
    rt_->deactivateTarget();
    return pixels;
}

tgt::col3 PlotPickingManager::getColorFromCell(int row, int column) {
    unsigned int color = row + 1; //we allow row = -1
    color = color << bitsReservedForColumns_;
    if (color > 0xFFFFFF) {
        tgtAssert(false, "PlotPickingManager: Row out of bounds!");
        return tgt::col3(0, 0, 0);
    }

    color += column + 1; //no cell is mapped to black (0,0,0)
    if (column >= numbersReservedForColumns_) {
        tgtAssert(false, "PlotPickingManager: Column out of bounds!");
        return tgt::col3(0, 0, 0);
    }
    //convert int to col3
    return tgt::col3(color >> 16, (color >> 8) % 256, color % 256);
}

tgt::ivec2 PlotPickingManager::getCellFromColor(tgt::col3 color) {
    //convert col3 to int
    unsigned int intColor = (color.r << 16) + (color.g << 8) + color.b;

    //return (row, column)
    return tgt::ivec2((intColor >> bitsReservedForColumns_) - 1, (intColor % numbersReservedForColumns_) - 1);
}

void PlotPickingManager::activateTarget(std::string debugLabel) {
    if (rt_) {
        rt_->activateTarget();
        rt_->setDebugLabel("ID target" + (debugLabel.empty() ? "" : " (" + debugLabel + ")"));
    }
    else
        LERRORC("voreen.plotpickingmanager", "No RenderTarget set!");
}

void PlotPickingManager::deactivateTarget() {
    if (rt_)
        rt_->deactivateTarget();
    else
        LERRORC("voreen.idmanager", "No RenderTarget set!");
}

void PlotPickingManager::clearTarget() {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}

void PlotPickingManager::setRenderTarget(RenderTarget* rt) {
    rt_ = rt;
}

RenderTarget* PlotPickingManager::getRenderTarget() {
    return rt_;
}

void PlotPickingManager::setColumnCount(int columnCount) {
    bitsReservedForColumns_ = static_cast<int>(floor(log(static_cast<float>(columnCount))/log(2.f))) + 1;
    numbersReservedForColumns_ = 1 << bitsReservedForColumns_;
}

tgt::Color PlotPickingManager::convertColor(tgt::col3 col) {
    return tgt::Color(static_cast<float>(col.r)/255.f,
        static_cast<float>(col.g)/255.f, static_cast<float>(col.b)/255.f,1.0);
}

} // namespace voreen
