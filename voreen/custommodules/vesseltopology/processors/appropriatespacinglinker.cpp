/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2015 University of Muenster, Germany.                        *
 * Visualization and Computer Graphics Group <http://viscg.uni-muenster.de>        *
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
#include "appropriatespacinglinker.h"

namespace voreen {

const std::string AppropriateSpacingLinker::loggerCat_("voreen.vesseltopology.appropriatespacinglinker");

AppropriateSpacingLinker::AppropriateSpacingLinker()
    : VolumeProcessor()
    , inport_(Port::INPORT, "volumethinning.inport", "Volume Input")
    , spacingDisplay_("spacingDisplay", "Resulting Spacing (mm)", tgt::vec3(1.0f), tgt::vec3(0.0f), tgt::vec3(1000.f))
    , multiplier_("multiplier", "Multiplier", 1.f, 0.000001f, 10.f)
{
    addPort(inport_);
    addProperty(multiplier_);
    addProperty(spacingDisplay_);
        spacingDisplay_.setReadOnlyFlag(true);
        spacingDisplay_.setNumDecimals(5);
}

AppropriateSpacingLinker::~AppropriateSpacingLinker() {
}

void AppropriateSpacingLinker::process() {
    const VolumeBase* invol = inport_.getData();
    if(!invol) {
        return;
    }
    tgt::vec3 current_spacing = invol->getSpacing();
    float wanted = tgt::min(current_spacing) * multiplier_.get();
    spacingDisplay_.set(tgt::vec3(wanted));
}
VoreenSerializableObject* AppropriateSpacingLinker::create() const {
    return new AppropriateSpacingLinker();
}

} // namespace voreen
