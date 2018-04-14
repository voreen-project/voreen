/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
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

#ifndef VRN_RWMULTILABELLOOPINITIALIZER_H
#define VRN_RWMULTILABELLOOPINITIALIZER_H

#include "voreen/core/processors/volumeprocessor.h"
#include "voreen/core/datastructures/volume/volumeatomic.h"

#include "voreen/core/properties/optionproperty.h"

#include "voreen/core/ports/geometryport.h"
#include "voreen/core/ports/loopport.h"

#include <string>
#include "tgt/vector.h"


namespace voreen {

class Volume;

class RWMultiLabelLoopInitializer : public VolumeProcessor {

public:
    RWMultiLabelLoopInitializer();
    virtual ~RWMultiLabelLoopInitializer();
    virtual Processor* create() const;

    virtual std::string getCategory() const { return "Volume Processing"; }
    virtual std::string getClassName() const { return "RWMultiLabelLoopInitializer"; }
    virtual CodeState getCodeState() const { return CODE_STATE_EXPERIMENTAL; }

private:
    virtual void setDescriptions() {
        setDescription("Loop-wrapper for the RandomWalker processor for multi-label segmentations (experimental). <br/>"
                       "See: RWMultiLabelLoopFinalizer");
    }

    virtual void process();
    virtual void initialize();
    virtual void deinitialize();

    void startButtonClicked();

    GeometryPort seedsInport_;
    GeometryPort foregroundSeedsOutport_;
    GeometryPort backgroundSeedsOutport_;

    LoopPort loopInport_;

    ButtonProperty startButton_;
    ButtonProperty rwConnector_;

    bool running_;

    static const std::string loggerCat_; ///< category used in logging
};

} //namespace

#endif
