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

#ifndef VRN_SLICECOLORMAPPER_H
#define VRN_SLICECOLORMAPPER_H

#include "voreen/core/processors/renderprocessor.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"

namespace voreen {

class SliceColorMapper : public RenderProcessor {
public:
    SliceColorMapper();
    virtual ~SliceColorMapper();

    virtual std::string getCategory() const { return "Image Processing"; }
    virtual std::string getClassName() const { return "SliceColorMapper"; }
    virtual Processor::CodeState getCodeState() const { return CODE_STATE_EXPERIMENTAL; }
    virtual std::string getProcessorInfo() const;
    virtual Processor* create() const;

    virtual void initialize();
    virtual void deinitialize();
    virtual void process();
protected:
    virtual void setDescriptions() {
        setDescription("");
    }

    virtual void compile();

    TransFunc1DKeysProperty transferFunc_;

    RenderPort inport_;
    RenderPort outport_;

    tgt::Shader* shader_;
};


} // namespace

#endif // VRN_SLICECOLORMAPPER_H
