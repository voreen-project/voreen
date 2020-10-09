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

#ifndef VRN_STREAMLINECOMBINE_H
#define VRN_STREAMLINECOMBINE_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/properties/optionproperty.h"

#include "../../ports/streamlinelistport.h"

namespace voreen {

    /**
     * Processor used to combine two StreamlineLists or pipe through on list.
     */
class VRN_CORE_API StreamlineCombine : public Processor {
public:
    StreamlineCombine();
    virtual ~StreamlineCombine();
    virtual Processor* create() const         {return new StreamlineCombine(); }
    virtual std::string getClassName() const  { return "StreamlineCombine";  }
    virtual std::string getCategory() const   { return "Process";          }
    virtual CodeState getCodeState() const    { return CODE_STATE_STABLE; }
    virtual bool isEndProcessor() const       { return true;              }

protected:
    virtual void setDescriptions() {
        setDescription("Processor used to combine two StreamlineLists or pipe through on list.");
    }

    virtual void process();

    //--------------
    //  Member
    //--------------
    /** Supported combine options */
    enum StreamlineCombineOptions {
        SCB_LEFT,
        SCB_RIGHT,
        SCB_COMBINE
    };

private:
    StreamlineListPort leftInport_;     ///< left inport used in process
    StreamlineListPort rightInport_;     ///< left inport used in process
    StreamlineListPort outport_;     ///< outport containing the combination

    OptionProperty<StreamlineCombineOptions> combineProp_;   ///< the choosen combination mode

    static const std::string loggerCat_;
};

}   //namespace

#endif
