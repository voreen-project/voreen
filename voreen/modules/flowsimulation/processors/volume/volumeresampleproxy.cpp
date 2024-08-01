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

#include "volumeresampleproxy.h"

#include "voreen/core/datastructures/volume/volumehash.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/datastructures/volume/volumeminmaxmagnitude.h"
#include "voreen/core/utils/hashing.h"

#include "../../datastructures/volumeramremappingproxy.h"

#include "modules/ensembleanalysis/utils/utils.h"

namespace voreen {

const std::string VolumeResampleProxy::loggerCat_("voreen.VolumeResampleProxy");

VolumeResampleProxy::VolumeResampleProxy()
    : inport_(Port::INPORT, "volumehandle.volumehandle2", "Volume Inport")
    , outport_(Port::OUTPORT, "volumehandle.volumehandle3", "Volume Outport")
    , outputDimensions_("outputDimensions", "Output Resolution", tgt::ivec3(10, 10, 10), tgt::ivec3(2), tgt::ivec3(1000))
    , resetResolution_("resetResolution", "Reset Resolution")
    , autoResetResolution_("autoResetResolution", "Reset Resolution on Input Change", true)
{
    addPort(inport_);
    inport_.Observable<PortObserver>::addObserver(static_cast<PortObserver*>(this));
    addPort(outport_);

    // Output dimensions
    addProperty(outputDimensions_);
    outputDimensions_.setGroupID("output");
    addProperty(resetResolution_);
    ON_CHANGE(resetResolution_, VolumeResampleProxy, adjustPropertiesToInput);
    resetResolution_.setGroupID("output");
    addProperty(autoResetResolution_);
    autoResetResolution_.setGroupID("output");
    setPropertyGroupGuiName("output","Output");
}

VolumeResampleProxy::~VolumeResampleProxy() {
    inport_.Observable<PortObserver>::removeObserver(static_cast<PortObserver*>(this));
}

Processor* VolumeResampleProxy::create() const {
    return new VolumeResampleProxy();
}

void VolumeResampleProxy::adjustPropertiesToInput() {
    if(auto* volume = inport_.getData()) {
        tgt::ivec3 dim = volume->getDimensions();
        outputDimensions_.setMaxValue(dim*10);
        if(autoResetResolution_.get()) {
            outputDimensions_.set(dim);
        }
    }
}

void VolumeResampleProxy::process() {

    auto inputVolume = inport_.getData();

    tgt::vec3 originalDimensions(inputVolume->getDimensions());
    tgt::vec3 outputDimensions(outputDimensions_.get());

    tgt::vec3 scaling = originalDimensions / outputDimensions;

    auto remappingFunction = [scaling] (tgt::vec3& pos) -> bool {
        pos = tgt::vec3(pos) * scaling;
        return true;
    };

    VolumeRAMRemappingProxy* representation = new VolumeRAMRemappingProxy(outputDimensions_.get(),
                                                                          inputVolume,
                                                                          remappingFunction);

    tgt::vec3 spacing = inputVolume->getSpacing() * scaling;
    tgt::vec3 offset = inputVolume->getOffset();

    auto output = new Volume(representation, spacing, offset);

    // Copy derived data of original volume - the data didn't change and recalculation on the remapped data
    // both takes a long time and will call getData() internally, which might not fit into main memory..
    std::ostringstream configStr;
    configStr << outputDimensions_.get();
    std::string hash = VoreenHash::getHash(inputVolume->getHash() + "-" + configStr.str());
    output->addDerivedData(new VolumeHash(hash));
    output->addDerivedData(new VolumeMinMax(*inputVolume->getDerivedData<VolumeMinMax>()));
    //output->addDerivedData(new VolumeMinMaxMagnitude(*inputVolume->getDerivedData<VolumeMinMaxMagnitude>()));
    output->setRealWorldMapping(inputVolume->getRealWorldMapping());

    outport_.setData(output);
}

void VolumeResampleProxy::afterConnectionAdded(const Port* source, const Port* connectedPort) {
}
void VolumeResampleProxy::beforeConnectionRemoved(const Port* source, const Port*) {
}

void VolumeResampleProxy::dataWillChange(const Port* source) {
    outport_.clear();

    if(source->hasData()) {
        source->removeDataInvalidationObserver(static_cast<DataInvalidationObserver*>(this));
    }
}
void VolumeResampleProxy::dataHasChanged(const Port* source) {
    if(source->hasData()) {
        source->addDataInvalidationObserver(static_cast<DataInvalidationObserver*>(this));
    }
}

void VolumeResampleProxy::dataAboutToInvalidate(const DataInvalidationObservable* data) {
    outport_.clear();
}


} // namespace
