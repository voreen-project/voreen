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

#include "ensemblevolumeextractor.h"

#include "voreen/core/datastructures/volume/volumelist.h"

#include "../datastructures/ensembledataset.h"

namespace voreen {

const std::string EnsembleVolumeExtractor::loggerCat_("voreen.ensembleanalysis.EnsembleVolumeExtractor");

EnsembleVolumeExtractor::EnsembleVolumeExtractor()
    : Processor()
    , inport_(Port::INPORT, "input", "EnsembleDataset Input", false)
    , outport_(Port::OUTPORT, "output", "VolumeList Output", false)
{
    ON_CHANGE(inport_, EnsembleVolumeExtractor, updateOutput);

    addPort(inport_);
    addPort(outport_);
}

EnsembleVolumeExtractor::~EnsembleVolumeExtractor() {
}

Processor* EnsembleVolumeExtractor::create() const {
    return new EnsembleVolumeExtractor();
}

void EnsembleVolumeExtractor::process() {
    if(!inport_.hasData())
        outport_.setData(nullptr);
    else {
        VolumeList* volumeList = new VolumeList();

        for(const VolumeBase* volume : inport_.getData()->getVolumes())
            volumeList->add(const_cast<VolumeBase*>(volume));

        outport_.setData(volumeList, true);
    }
}

void EnsembleVolumeExtractor::updateOutput() {
    invalidate();
}

} // namespace
