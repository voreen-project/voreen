/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany,                        *
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

#ifndef VRN_FLOWSTREAMLINES3D_H
#define VRN_FLOWSTREAMLINES3D_H

#include <string>
#include "voreen/core/processors/processor.h"
#include "flowreenprocessor.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/ports/volumeport.h"
namespace voreen {

class Volume;

class FlowStreamlinesTexture3D : public Processor, private FlowreenProcessor {
public:
    FlowStreamlinesTexture3D();
    virtual ~FlowStreamlinesTexture3D();

    virtual Processor* create() const { return new FlowStreamlinesTexture3D(); }
    virtual std::string getCategory() const { return "Flow Visualization"; }
    virtual std::string getClassName() const { return "StreamlineTexture3D"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_EXPERIMENTAL; }

    virtual void process();

private:
    virtual void setDescriptions() {
        setDescription("");
    }
    void calculateStreamlines();

private:
    Volume* processedVolume_;
    IntProperty voxelSamplingProp_;

    VolumePort volInport_;
    VolumePort volOutport_;
};

}   // namespace

#endif
