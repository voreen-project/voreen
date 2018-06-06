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

#ifndef VRN_FLOWREENADAPTER_H
#define VRN_FLOWREENADAPTER_H

#include "voreen/core/processors/processor.h"
#include "flowreenprocessor.h"
#include "voreen/core/ports/volumeport.h"
#include <string>

namespace voreen {

class Volume;

class FlowreenAdapter : public Processor, private FlowreenProcessor {
public:
    FlowreenAdapter();
    virtual ~FlowreenAdapter();
    virtual Processor* create() const { return new FlowreenAdapter(); }

    virtual std::string getClassName() const { return "FlowreenAdapter"; }
    virtual std::string getCategory() const { return "Flow Visualization"; }
    virtual CodeState getCodeState() const { return CODE_STATE_EXPERIMENTAL; }
    virtual bool isUtility() const { return true; }

protected:
    virtual void setDescriptions() {
        setDescription("Converts a gradient volume (3 channels, float) into a flow volume.");
    }

    virtual void process();
    void initialize();

private:
    void calculateMagnitudes();

private:
    Volume* processedVolume_;   /** The volume containing the magnitudes. */

    VolumePort volInport_;
    VolumePort volOutport_;
};

}   //namespace

#endif
