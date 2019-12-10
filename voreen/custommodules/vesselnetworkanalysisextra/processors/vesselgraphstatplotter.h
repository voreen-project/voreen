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

#ifndef VRN_VESSELGRAPHSTATPLOTTER_H
#define VRN_VESSELGRAPHSTATPLOTTER_H

#include "voreen/core/processors/processor.h"

#include "modules/vesselnetworkanalysis/ports/vesselgraphport.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/buttonproperty.h"

#include "modules/plotting/ports/plotport.h"

#include <functional>

namespace voreen {

// Extracts featues from voxels in edges of VesselGraphs and makes them available via
// a PlotPort. Additionally per-voxel properties can be exported to a csv-file.
// Moreover, per-edge features are computed and availabe via read-only properties.
class VesselGraphStatPlotter : public Processor {
public:
    VesselGraphStatPlotter();
    virtual ~VesselGraphStatPlotter();
    virtual std::string getCategory() const { return "VesselGraph"; }
    virtual std::string getClassName() const { return "VesselGraphStatPlotter"; }
    virtual CodeState getCodeState() const { return Processor::CODE_STATE_EXPERIMENTAL; }
    virtual Processor* create() const { return new VesselGraphStatPlotter(); }

protected:
    virtual void setDescriptions() {
        setDescription("Plots stats for Vessel-Graphs.");
    }

    virtual void process();

    void adaptToNewInput();

    void exportToFile(const VesselGraphEdge& edge);

    // ports
    VesselGraphPort graphInport_;
    PlotPort plotOutport_;

    // properties
    FileDialogProperty exportFilePath_;
    ButtonProperty exportButton_;
    IntProperty activeEdgeID_;
    FloatProperty length_;
    FloatProperty distance_;
    FloatProperty curveness_;
    FloatProperty minRadiusAvg_;
    FloatProperty minRadiusStdDeviation_;
    FloatProperty avgRadiusAvg_;
    FloatProperty avgRadiusStdDeviation_;
    FloatProperty maxRadiusAvg_;
    FloatProperty maxRadiusStdDeviation_;
    FloatProperty roundnessAvg_;
    FloatProperty roundnessStdDeviation_;
    FloatProperty avgCrossSection_;
    FloatProperty volume_;
    IntProperty numEdgesLeftNode_;
    IntProperty numEdgesRightNode_;
    IntProperty numSkeletonVoxels_;

    bool exportForced_;

    static const std::string loggerCat_;

};

} // namespace voreen
#endif // VRN_VESSELGRAPHSTATPLOTTER_H
