/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#ifndef VRN_VOLUMELISTMERGER_H
#define VRN_VOLUMELISTMERGER_H

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/genericport.h"

namespace voreen {

/**
 * This processor merges multiple lists into a single one.
 */
class VRN_CORE_API VolumeListMerger : public Processor {
public:
    VolumeListMerger();
    virtual ~VolumeListMerger();
    virtual Processor* create() const;

    virtual std::string getClassName() const  { return "VolumeListMerger";      }
    virtual std::string getCategory() const   { return "Volume Processing";     }
    virtual CodeState getCodeState() const    { return CODE_STATE_EXPERIMENTAL; }

protected:
    virtual void setDescriptions() {
        setDescription("This processor merges incoming Lists into one single List. <br>"
                       "Given two lists [a1, a2, a3] and [b1, b2], this results in: <br>"
                       "[a1, b1, a2, b2], since the length of the shortest list will be taken."
        );
    }

    virtual void process();

private:

    VolumeListPort inport_;
    VolumeListPort outport_;
};

}   //namespace

#endif
