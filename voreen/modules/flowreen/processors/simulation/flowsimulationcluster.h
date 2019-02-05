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
#ifdef VRN_MODULE_ENSEMBLEANALYSIS
#include "custommodules/ensembleanalysis/ports/ensembledatasetport.h"
#endif

#include "voreen/core/properties/stringproperty.h"

#include "modules/hdf5/io/hdf5filevolume.h"

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

    // TODO: Create submit.cmd script in this processor and make setup configurable via UI.

    void fetchResults();
    int executeCommand(const std::string& command) const;

    GeometryPort geometryDataPort_;
    VolumeListPort measuredDataPort_; // TODO: Currently ignored.
    FlowParametrizationPort parameterPort_;
#ifdef VRN_MODULE_ENSEMBLEANALYSIS
    EnsembleDatasetPort ensemblePort_;
#endif

    StringProperty username_;
    StringProperty clusterAddress_;
    FileDialogProperty simulationResults_;
    ButtonProperty triggerFetchResults_;

    static const std::string loggerCat_;
};

}   //namespace

#endif
