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

#include "flowsimulationmodule.h"

// processors
#include "processors/geometry/geometryclose.h"
#include "processors/geometry/geometryoffsetremove.h"
#include "processors/render/unalignedsliceviewer.h"
#include "processors/simulation/flowcharacteristics.h"
#include "processors/simulation/flowensemblecreator.h"
//#include "processors/simulation/flowindicatorselection.h"
#include "processors/simulation/flowindicatorrenderer.h"
#include "processors/simulation/flowparametrization.h"
#include "processors/simulation/flowsimulationcluster.h"
#include "processors/simulation/flowsimulationgeometry.h"
#include "processors/volume/volumelistadapter.h"

#ifdef VRN_MODULE_VESSELNETWORKANALYSIS
#include "processors/simulation/flowindicatordetection.h"
#endif

#ifdef VRN_FLOWSIMULATION_USE_OPENLB
#include <olb3D.h>
#include "processors/simulation/flowsimulation.h"
#include "processors/geometry/implicitrepresentation.h"
#endif

namespace voreen {

FlowSimulationModule::FlowSimulationModule(const std::string& modulePath)
    : VoreenModule(modulePath)
{
    setID("FlowSimulation");
    setGuiName("FlowSimulation");

    // processors
    registerSerializableType(new GeometryClose());
    registerSerializableType(new GeometryOffsetRemove());
    registerSerializableType(new UnalignedSliceViewer());
    registerSerializableType(new FlowCharacteristics());
    registerSerializableType(new FlowEnsembleCreator());
    //registerSerializableType(new FlowIndicatorSelection());
    registerSerializableType(new FlowIndicatorRenderer());
    registerSerializableType(new FlowParametrization());
    registerSerializableType(new FlowSimulationCluster());
    registerSerializableType(new FlowSimulationGeometry());
    registerSerializableType(new VolumeListAdapter());
#ifdef VRN_MODULE_VESSELNETWORKANALYSIS
    registerSerializableType(new FlowIndicatorDetection());
#endif
#ifdef VRN_FLOWSIMULATION_USE_OPENLB
    registerSerializableType(new ImplicitRepresentation());
    registerSerializableType(new FlowSimulation());
#endif
}

void FlowSimulationModule::initialize() {
    VoreenModule::initialize();

#ifdef VRN_FLOWSIMULATION_USE_OPENLB
    olb::olbInit(nullptr, nullptr);
    olb::singleton::directories().setOutputDir(VoreenApplication::app()->getTemporaryPath("simulation")+"/");
#endif
}

void FlowSimulationModule::deinitialize() {
    VoreenModule::deinitialize();
}


} // namespace
