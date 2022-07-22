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

#ifndef VRN_ENSEMBLECHANNELMERGER_H
#define VRN_ENSEMBLECHANNELMERGER_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "voreen/core/properties/string/stringlistproperty.h"

#include "modules/base/properties/interactivelistproperty.h"

#include "../ports/ensembledatasetport.h"

namespace voreen {

/**
 * This processor allows to merger channels of volumes, if these have homogeneous and identical resolutions.
 */
class VRN_CORE_API EnsembleChannelMerger : public Processor {
public:
    EnsembleChannelMerger();
    virtual ~EnsembleChannelMerger();

    virtual Processor* create() const;
    virtual std::string getClassName() const { return "EnsembleChannelMerger"; }
    virtual std::string getCategory() const  { return "Filter";                }
    virtual CodeState getCodeState() const   { return CODE_STATE_EXPERIMENTAL; }
    virtual bool isUtility() const           { return true;                    }

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

protected:

    void process();

    void merge();

    void adjustToEnsemble();

    /// Inport for the ensemble data structure.
    EnsembleDatasetPort ensembleInport_;

    /// The ensemble data
    EnsembleDatasetPort ensembleOutport_;

    // Properties
    InteractiveListProperty fields_;///< A list of all fields.
    ButtonProperty apply_;          ///< Apply settings.
    ButtonProperty reset_;          ///< Reset to input.

    /// Hash value of last valid data.
    std::string hash_;

    std::vector<std::unique_ptr<VolumeBase>> volumes_;
};

} // namespace

#endif
