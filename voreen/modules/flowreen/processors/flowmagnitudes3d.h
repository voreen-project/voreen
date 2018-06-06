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

#ifndef VRN_FLOWMAGNITUDES3D_H
#define VRN_FLOWMAGNITUDES3D_H

#include "voreen/core/processors/processor.h"
#include "flowreenprocessor.h"
#include "voreen/core/ports/volumeport.h"

#include <string>

namespace voreen {

class Volume;

/**
 * A volume pre-processor which takes a volume containing flow data for input
 * and calculates an output volume containing only the flow magnitudes.
 * The resulting output then can be treated as any other volume and can be
 * rendered using Raycaster or SliceRenderer.
 */
class FlowMagnitudes3D : public Processor, private FlowreenProcessor {
public:
    FlowMagnitudes3D();
    virtual ~FlowMagnitudes3D();

    virtual Processor* create() const { return new FlowMagnitudes3D(); }
    virtual std::string getCategory() const { return "Flow Visualization"; }
    virtual std::string getClassName() const { return "FlowMagnitudes3D"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_STABLE; }
    virtual void process();

protected:
    virtual void setDescriptions() {
        setDescription("Calculates the magnitudes of the incomming Flow Volume and outputs them into a new Volume which can be used as a normal 3D texture.");
    }

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
