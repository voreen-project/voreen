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

#ifndef VRN_VESSELGRAPHGLOBALSTATS_H
#define VRN_VESSELGRAPHGLOBALSTATS_H

#include "voreen/core/processors/processor.h"

#include "../ports/vesselgraphport.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/buttonproperty.h"

#include <functional>

namespace voreen {

class VesselGraphGlobalStats : public Processor {
public:
    VesselGraphGlobalStats();
    virtual ~VesselGraphGlobalStats();
    virtual std::string getCategory() const { return "Output"; }
    virtual std::string getClassName() const { return "VesselGraphGlobalStats"; }
    virtual CodeState getCodeState() const { return Processor::CODE_STATE_TESTING; }
    virtual Processor* create() const { return new VesselGraphGlobalStats(); }

protected:
    virtual void setDescriptions() {
        setDescription("This processor enables the export of per-node and per-edge properties as a csv-File for VesselGraphs.");
        segmentExportFilePath_.setDescription("Output path for the .csv file containing the edge data and properties.0");
        nodeExportFilePath_.setDescription("Output path for the .csv file containing the node data.");
        autoExport_.setDescription("Automatically overwrite specified files when the input VesselGraph changes. If not enabled, the export can be triggered manually using the 'Export'-Button.");
    }

    virtual void process();

    virtual void initialize();
    virtual void deinitialize();

    virtual bool isEndProcessor() const;

    void adaptToNewInput();

    // Export per-node and per-edge features to csv-files
    // (defined by segmentExportFilePath_, nodeExportFilePath_).
    void exportToFile(const VesselGraph& graph);

    // ports
    VesselGraphPort graphInport_;

    // properties
    FileDialogProperty segmentExportFilePath_;
    FileDialogProperty nodeExportFilePath_;
    BoolProperty autoExport_;
    ButtonProperty exportButton_;

    IntProperty numNodes_;
    IntProperty numEdges_;

    bool exportForced_;

    static const std::string loggerCat_;

};

} // namespace voreen
#endif // VRN_VESSELGRAPHGLOBALSTATS_H
