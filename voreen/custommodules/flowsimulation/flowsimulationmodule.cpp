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
#include "processors/features/wallshearstress.h"
#include "processors/geometry/geometryclose.h"
#include "processors/geometry/geometrysmoothnormals.h"
#include "processors/render/unalignedsliceviewer.h"
#include "processors/simulation/flowcharacteristics.h"
#include "processors/simulation/flowensemblecreator.h"
#include "processors/simulation/flowindicatorrenderer.h"
#include "processors/simulation/flowparametrizationensemble.h"
#include "processors/simulation/flowparametrizationrun.h"
#include "processors/simulation/flowsimulationcluster.h"
#include "processors/simulation/flowsimulationgeometry.h"
#include "processors/volume/connectedcomponentselector.h"
#include "processors/volume/flowtestdatagenerator.h"
#include "processors/volume/phaseunwrapping.h"
#include "processors/volume/vectordecompose.h"
#include "processors/volume/volumeapplyrealworldmapping.h"
#include "processors/volume/volumelistadapter.h"
#include "processors/volume/volumelistaggregate.h"
#include "processors/volume/volumelistmerger.h"
#include "processors/volume/volumelistoffset.h"
#include "processors/volume/volumelistrealworldmapping.h"
#include "processors/volume/volumelisttimestep.h"
#include "processors/volume/volumemerger.h"
#include "processors/volume/volumenoise.h"
#include "processors/volume/volumeselectormultichannel.h"

#ifdef VRN_MODULE_VESSELNETWORKANALYSIS
#include "processors/simulation/flowindicatordetection.h"
#endif

#ifdef VRN_MODULE_PLOTTING
#include "processors/plotting/roianalysis.h"
#include "processors/simulation/flowindicatoranalysis.h"
#endif

#ifdef VRN_FLOWSIMULATION_USE_OPENLB
#include <olb3D.h>
#include "processors/simulation/flowsimulation.h"
#include "processors/geometry/geometryinsidetest.h"
#endif

namespace voreen {

FlowSimulationModule::FlowSimulationModule(const std::string& modulePath)
    : VoreenModule(modulePath)
{
    setID("Flow Simulation");
    setGuiName("Flow Simulation");

    // processors
    registerProcessor(new ConnectedComponentSelector());
    registerProcessor(new GeometryClose());
    registerProcessor(new GeometrySmoothNormals());
    registerProcessor(new UnalignedSliceViewer());
    registerProcessor(new FlowCharacteristics());
    registerProcessor(new FlowEnsembleCreator());
    registerProcessor(new FlowIndicatorRenderer());
    registerProcessor(new FlowParametrizationEnsemble());
    registerProcessor(new FlowParametrizationRun());
    registerProcessor(new FlowSimulationCluster());
    registerProcessor(new FlowSimulationGeometry());
    registerProcessor(new FlowTestDataGenerator());
    registerProcessor(new PhaseUnwrapping());
    registerProcessor(new VectorDecompose());
    registerProcessor(new VolumeApplyRealWorldMapping());
    registerProcessor(new VolumeListAdapter());
    registerProcessor(new VolumeListAggregate());
    registerProcessor(new VolumeListMerger());
    registerProcessor(new VolumeListOffset());
    registerProcessor(new VolumeListRealWorldMapping());
    registerProcessor(new VolumeListTimeStep());
    registerProcessor(new VolumeMerger());
    registerProcessor(new VolumeNoise());
    registerProcessor(new VolumeSelectorMultiChannel());
    registerProcessor(new WallShearStress());
#ifdef VRN_MODULE_VESSELNETWORKANALYSIS
    registerProcessor(new FlowIndicatorDetection());
#endif
#ifdef VRN_MODULE_PLOTTING
    registerProcessor(new RoiAnalysis());
    registerProcessor(new FlowIndicatorAnalysis());
#endif
#ifdef VRN_FLOWSIMULATION_USE_OPENLB
    registerProcessor(new GeometryInsideTest());
    registerProcessor(new FlowSimulation());
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
