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

#ifndef VRN_STREAMLINEPREDICATES_H
#define VRN_STREAMLINEPREDICATES_H

#include "voreen/core/processors/asynccomputeprocessor.h"

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"

#include "../../ports/streamlinelistport.h"
#include "voreen/core/ports/geometryport.h"

namespace voreen {

struct StreamlinePredicatesComputeInput {
    std::unique_ptr<StreamlineListBase> streamlines;
    PortDataPointer<VolumeBase> streamlinePredicateVolume;
    std::function<bool(const Streamline&, const std::function<bool(const tgt::vec3&)>&)> predicateFilter;
    tgt::vec2 physicalLengthRange;
    tgt::vec2 curvatureRange;
};

struct StreamlinePredicatesComputeOutput {
    std::unique_ptr<StreamlineListBase> streamlines;
};

/**
 * Used to filter streamlines according to (optionally) multiple predicates.
 */
class StreamlinePredicates : public AsyncComputeProcessor<StreamlinePredicatesComputeInput, StreamlinePredicatesComputeOutput> {
public:
    StreamlinePredicates();
    virtual ~StreamlinePredicates();

    virtual Processor* create() const { return new StreamlinePredicates(); }

    virtual std::string getCategory() const { return "Streamline Processing"; }
    virtual std::string getClassName() const { return "StreamlinePredicates"; }
    virtual CodeState getCodeState() const { return CODE_STATE_TESTING; }

    virtual ComputeInput prepareComputeInput();
    virtual ComputeOutput compute(ComputeInput input, ProgressReporter& progressReporter) const;
    virtual void processComputeOutput(ComputeOutput output);

protected:
    virtual void setDescriptions() {
        setDescription("Used to filter streamlines according a volume predicate.");
        //ports
        streamlineInport_.setDescription("Streamlines, which shall be filtered.");
        streamlinePredicateVolumeInport_.setDescription("(Optional) Predicate (aka mask) volume. Streamlines that do not intersect this volume will be discarded.");
        streamlineOutport_.setDescription("Filtered streamlines.");
        //properties
        predicateVolumeFiltering_.setDescription("A streamline will be kept if its discrete elements satisfy the predicate according to the selected option.");
    }

    virtual void adjustPropertiesToInput();
    virtual void dataWillChange(const Port* source);
    virtual bool isReady() const;

private:

    // ports
    StreamlineListPort streamlineInport_;
    VolumePort streamlinePredicateVolumeInport_;
    StreamlineListPort streamlineOutport_;

    // properties
    // enable
    BoolProperty enabled_;                       ///< toggles the processor on and off

    StringOptionProperty predicateVolumeFiltering_;
    FloatIntervalProperty physicalLengthRange_;
    IntIntervalProperty curvatureRange_;

    static const std::string loggerCat_;
};

}   // namespace

#endif
