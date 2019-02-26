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

#ifndef VRN_FLOWSIMULATIONCLUSTER_H
#define VRN_FLOWSIMULATIONCLUSTER_H

#include "voreen/core/processors/asynccomputeprocessor.h"

#include "voreen/core/ports/geometryport.h"
#include "voreen/core/ports/volumeport.h"
#include "modules/flowreen/ports/flowparametrizationport.h"

#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/progressproperty.h"
#include "voreen/core/properties/temppathproperty.h"

namespace voreen {

/**
 * This processor performs simulations on the PALMAII cluster at WWU Muenster using a parameter set and as input.
 * TODO: Eventually, this processor should be able to react to emails send by the cluster and process the result automatically.
 */
class VRN_CORE_API FlowSimulationCluster : public Processor {
public:
    FlowSimulationCluster();
    virtual ~FlowSimulationCluster();
    virtual Processor* create() const         { return new FlowSimulationCluster();    }

    virtual std::string getClassName() const  { return "FlowSimulationCluster";        }
    virtual std::string getCategory() const   { return "Simulation";                   }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL;        }

    virtual bool isReady() const;
    virtual void process();

protected:
    virtual void setDescriptions() {
        setDescription("This processor performs simulations on the PALMAII cluster at WWU Muenster using a parameter set and as input.");
    }

private:

    void enqueueSimulations();
    void fetchResults();
    int executeCommand(const std::string& command) const;
    std::string generateEnqueueScript(const std::string& parametrizationPath) const;
    std::string generateSubmissionScript(const std::string& parametrizationName) const;

    GeometryPort geometryDataPort_;
    VolumeListPort measuredDataPort_;
    FlowParametrizationPort parameterPort_;

    StringProperty username_;
    StringProperty clusterAddress_;
    StringProperty simulationPath_;
    StringOptionProperty simulationType_;

    IntProperty configNodes_;
    IntProperty configTasks_;
    IntProperty configTasksPerNode_;
    IntProperty configCPUsPerTask_;
    IntProperty configMemory_;
    IntProperty configTime_;
    StringOptionProperty configPartition_;

    FileDialogProperty simulationResults_;
    FileDialogProperty uploadDataPath_;
    ButtonProperty triggerEnqueueSimulations_;
    ButtonProperty triggerFetchResults_;
    ProgressProperty progress_;

    static const std::string loggerCat_;
};

}   //namespace

#endif
