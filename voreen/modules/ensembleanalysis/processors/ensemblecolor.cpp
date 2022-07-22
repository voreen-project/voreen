/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2020 University of Muenster, Germany,                        *
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

#include "ensemblecolor.h"

#include "voreen/core/datastructures/callback/lambdacallback.h"

#include "../utils/ensemblehash.h"

#include "../modules/plotting/datastructures/colormap.h"

namespace voreen {

EnsembleColor::EnsembleColor()
    : Processor()
    , ensembleInport_(Port::INPORT, "ensembledatastructurein", "Ensemble Datastructure Input", false)
    , ensembleOutport_(Port::OUTPORT, "ensembledatastructureout", "Ensemble Datastructure Output", false)
    , members_("members", "Members", Processor::VALID)
    , color_("color", "Color")
    , applyGradient_("applyGradient", "Apply Gradient")
    , gradientColor_("gradientColor", "Gradient Color", tgt::vec4::zero)
    , apply_("apply", "Apply")
    , reset_("reset", "Reset")
    , needsProcess_(false)
{
    addPort(ensembleInport_);
    ON_CHANGE(ensembleInport_, EnsembleColor, adjustToEnsemble);
    addPort(ensembleOutport_);

    addProperty(members_);
    ON_CHANGE_LAMBDA(members_, [this] {
        std::vector<int> selection = members_.getSelectedRowIndices();
        if(selection.size() == 1) {
            color_.set(tgt::vec4(colors_[selection.front()], 1.0f));
        }
    });
    addProperty(color_);
    addProperty(applyGradient_);
    ON_CHANGE_LAMBDA(applyGradient_, [this] {
        gradientColor_.setReadOnlyFlag(!applyGradient_.get());
    });
    addProperty(gradientColor_);
    gradientColor_.setReadOnlyFlag(!applyGradient_.get());
    addProperty(apply_);
    ON_CHANGE(apply_, EnsembleColor, applyColors);
    addProperty(reset_);
    ON_CHANGE(reset_, EnsembleColor, adjustToEnsemble);
}

EnsembleColor::~EnsembleColor() {
}

Processor* EnsembleColor::create() const {
    return new EnsembleColor();
}

void EnsembleColor::process() {
    //if(needsProcess_)
    {
        ensembleOutport_.clear();

        members_.reset();
        std::unique_ptr<EnsembleDataset> ensemble(new EnsembleDataset());
        const auto& members = ensembleInport_.getData()->getMembers();
        for(size_t i=0; i<members.size(); i++) {
            members_.addRow(members[i].getName(), colors_[i]);
            ensemble->addMember(EnsembleMember(members[i].getName(), colors_[i], members[i].getTimeSteps()));
        }

        ensembleOutport_.setData(ensemble.release(), true);

        needsProcess_ = false;
    }
}

void EnsembleColor::adjustToEnsemble() {

    const auto* ensemble = ensembleInport_.getData();
    if(ensemble) {

        std::string hash = EnsembleHash(*ensemble).getHash();
        if(hash != hash_) {

            const auto& members = ensemble->getMembers();
            colors_.resize(members.size());

            for(size_t i=0; i<members.size(); i++) {
                colors_[i] = members[i].getColor();
            }

            hash_ = hash;
        }

        needsProcess_ = true;
    }
}

void EnsembleColor::applyColors() {
    std::vector<int> selection = members_.getSelectedRowIndices();
    if(selection.empty()) {
        return;
    }

    if(applyGradient_.get()) {
        std::vector<tgt::Color> colors;
        colors.push_back(color_.get());
        colors.push_back(gradientColor_.get());
        ColorMap colorMap = ColorMap::createFromVector(colors);

        ColorMap::InterpolationIterator colorIter = colorMap.getInterpolationIterator(selection.size());
        for(int idx : selection) {
            colors_[idx] = (*colorIter).xyz();
            ++colorIter;
        }
    }
    else {
        for(int idx : selection) {
            colors_[idx] = color_.get().xyz();
        }
    }

    needsProcess_ = true;
}

void EnsembleColor::serialize(Serializer& s) const {
    Processor::serialize(s);
    s.serialize("colors", colors_);
    s.serialize("hash", hash_);
}

void EnsembleColor::deserialize(Deserializer& s) {
    Processor::deserialize(s);
    s.optionalDeserialize("hash", hash_, std::string(""));
    s.deserialize("colors", colors_);
}

} // namespace
