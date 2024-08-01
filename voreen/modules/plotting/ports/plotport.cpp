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

#include "plotport.h"
#include "../datastructures/plotbase.h"

namespace voreen {

PlotPort::PlotPort(PortDirection direction, const std::string& name, const std::string& guiName,
        bool allowMultipleConnections, Processor::InvalidationLevel invalidationLevel)
 : GenericPort<PlotBase>(direction, name, guiName, allowMultipleConnections, invalidationLevel)
{
    if (isOutport())
        setData(0);
}

PlotPort::~PlotPort() {
    if ((portData_) && (isOutport()) /*&& (portData_->getColumnCount() > 0)*/) {
        delete portData_;
        portData_ = 0;
    }
}

bool PlotPort::isReady() const {
    if (isOutport())
        return isConnected();
    else
        return (hasData() && areConditionsMet());
}

void PlotPort::setData(const PlotBase* handle, bool deletePrevious) {
    tgtAssert(isOutport(), "called setData on inport!");
    const PlotBase* tempVol = portData_;
    if (handle != 0 || portData_ != 0) {
        portData_ = handle;
        invalidatePort();
    }
    else {
        portData_ = 0;
    }
    if (deletePrevious && tempVol && (tempVol != handle)) {
        delete tempVol;
    }
}

void PlotPort::deleteData() {
    tgtAssert(isOutport(), "deletePlotBase called on inport!");
    const PlotBase* tempVol = portData_;
    if (tempVol) {
        setData(0);
        delete tempVol;
    }
}

} // namespace
