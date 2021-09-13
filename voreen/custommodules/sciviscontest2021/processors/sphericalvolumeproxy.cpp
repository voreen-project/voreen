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

#include "sphericalvolumeproxy.h"

#include "voreen/core/datastructures/volume/volumehash.h"
#include "voreen/core/datastructures/volume/volumeminmax.h"
#include "voreen/core/datastructures/volume/volumeminmaxmagnitude.h"
#include "voreen/core/utils/hashing.h"

#include "../datastructures/sphericalvolumeramproxy.h"

#include "modules/ensembleanalysis/utils/utils.h"

namespace voreen {

const std::string SphericalVolumeProxy::loggerCat_("voreen.SphericalVolumeProxy");

SphericalVolumeProxy::SphericalVolumeProxy()
    : inport_(Port::INPORT, "volumehandle.volumehandle2", "Volume Inport")
    , outport_(Port::OUTPORT, "volumehandle.volumehandle3", "Volume Outport")
    , radiusMin_("radiusMin", "Set Radius Min", 3485.0f, 0.0f, 10000.0f)
    , radiusMax_("radiusMax", "Set Radius Max", 6371.0f, 0.0f, 10000.0f)
    , clampRadiusMin_("clampingRadiusMin", "Set clamp Radius Min", 3485.0f, 0.0f, 10000.0f)
    , clampRadiusMax_("clampingRadiusMax", "Set clamp Radius Max", 6371.0f, 0.0f, 10000.0f)
    , shiftFactor_("shiftFactor", "Set shift factor", 100.0f, 0.0f, 1000.0f)
    , outputDimensions_("outputDimensions", "Set desired output dimensions", 100, 2, 2048)

{
    addPort(inport_);
    inport_.Observable<PortObserver>::addObserver(static_cast<PortObserver*>(this));
    addPort(outport_);

    //Radius
    addProperty(radiusMin_);
    addProperty(radiusMax_);
    addProperty(clampRadiusMin_);
    addProperty(clampRadiusMax_);
    addProperty(shiftFactor_);
    radiusMin_.setGroupID("radius");
    radiusMax_.setGroupID("radius");
    clampRadiusMin_.setGroupID("radius");
    clampRadiusMax_.setGroupID("radius");
    shiftFactor_.setGroupID("radius");
    setPropertyGroupGuiName("radius","Radius");

    // Output dimensions
    addProperty(outputDimensions_);
    outputDimensions_.setGroupID("output");
    setPropertyGroupGuiName("output","Output dimensions");
}

SphericalVolumeProxy::~SphericalVolumeProxy() {
    inport_.Observable<PortObserver>::removeObserver(static_cast<PortObserver*>(this));
}

Processor* SphericalVolumeProxy::create() const {
    return new SphericalVolumeProxy();
}

void SphericalVolumeProxy::process() {

    auto inputVolume = inport_.getData();

    SphericalVolumeRAMProxy* representation = new SphericalVolumeRAMProxy(outputDimensions_.get(),
                                                                          inputVolume,
                                                                          radiusMin_.get() / shiftFactor_.get(),
                                                                          radiusMax_.get() / shiftFactor_.get(),
                                                                          clampRadiusMin_.get() / shiftFactor_.get(),
                                                                          clampRadiusMax_.get() / shiftFactor_.get());

    tgt::vec3 spacing = tgt::vec3(radiusMax_.get() * 2.0f / shiftFactor_.get() / outputDimensions_.get());
    tgt::vec3 offset = -spacing * tgt::vec3(outputDimensions_.get() / 2);

    auto output = new Volume(representation, spacing, offset);

    // Copy derived data of original volume - the data didn't change and recalculation on the remapped data
    // both takes a long time and will call getData() internally, which might not fit into main memory..
    std::ostringstream configStr;
    configStr << radiusMin_.get() << "-" << radiusMax_.get() << "-"
              << clampRadiusMin_.get() << "-" << clampRadiusMax_.get() << "-"
              << shiftFactor_.get();
    std::string hash = VoreenHash::getHash(inputVolume->getHash() + "-" + configStr.str());
    output->addDerivedData(new VolumeHash(hash));
    output->addDerivedData(new VolumeMinMax(*inputVolume->getDerivedData<VolumeMinMax>()));
    output->addDerivedData(new VolumeMinMaxMagnitude(*inputVolume->getDerivedData<VolumeMinMaxMagnitude>()));
    output->setRealWorldMapping(inputVolume->getRealWorldMapping());

    outport_.setData(output);
}

void SphericalVolumeProxy::afterConnectionAdded(const Port* source, const Port* connectedPort) {
}
void SphericalVolumeProxy::beforeConnectionRemoved(const Port* source, const Port*) {
}

void SphericalVolumeProxy::dataWillChange(const Port* source) {
    outport_.clear();

    if(source->hasData()) {
        source->removeDataInvalidationObserver(static_cast<DataInvalidationObserver*>(this));
    }
}
void SphericalVolumeProxy::dataHasChanged(const Port* source) {
    if(source->hasData()) {
        source->addDataInvalidationObserver(static_cast<DataInvalidationObserver*>(this));
    }
}

void SphericalVolumeProxy::dataAboutToInvalidate(const DataInvalidationObservable* data) {
    outport_.clear();
}


} // namespace
