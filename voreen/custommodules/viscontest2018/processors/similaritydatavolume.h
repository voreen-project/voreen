/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2017 University of Muenster, Germany.                        *
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

#ifndef VRN_SIMILARITYDATAVOLUME_H
#define VRN_SIMILARITYDATAVOLUME_H

#include "voreen/core/processors/renderprocessor.h"

#include "voreen/core/io/volumeserializerpopulator.h"

#include "voreen/core/ports/volumeport.h"

#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/progressproperty.h"
#include "voreen/core/properties/shaderproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"

#include "../datastructures/similaritydata.h"
#include "../ports/ensembledatasetport.h"
#include "../ports/similaritydataport.h"
#include "../properties/stringlistproperty.h"

namespace voreen {

/**
 *
 */
class VRN_CORE_API SimilartyDataVolume : public RenderProcessor {

public:
    SimilartyDataVolume();
    virtual ~SimilartyDataVolume();
    virtual Processor* create() const;

    virtual std::string getClassName() const      { return "SimilartyDataVolume"; }
    virtual std::string getCategory() const       { return "Plotting";                }
    virtual CodeState getCodeState() const        { return CODE_STATE_EXPERIMENTAL;   }

protected:

    virtual void onEvent(tgt::Event* e);

    virtual void initialize();
    virtual void deinitialize();

    virtual void process();

    virtual bool isReady() const;

    void onInportChange();
    void initSimilarityVolume();
    float calculateVariance(const std::vector<float> voxelData);
    float calculateMinMaxDiff(const std::vector<float> voxelData);
    const std::vector<float> applyGroupLogic(const std::vector<float> rawVoxelData);

protected:

    virtual void setDescriptions() {
        setDescription("");
    }

    SimilarityDataPort inport_;
    VolumePort outport_;
    StringOptionProperty similarityMethod_;

    VolumeBase* similarityVolume_;
    VolumeRAM_Float* similarityVolumeRepresentation_;

    StringListProperty group1_;
    StringListProperty group2_;

    static const std::string loggerCat_;
};

} // namespace

#endif
