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

#ifndef VRN_TEXT_H
#define VRN_TEXT_H

#include <string>
#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/textport.h"
#include "voreen/core/properties/stringproperty.h"

namespace voreen {

/**
 * Provides a textport output with content defined by the corresponding property.
 */
class VRN_CORE_API TextSource : public Processor {
public:
    TextSource();
    virtual ~TextSource();
    virtual Processor* create() const;

    virtual std::string getClassName() const { return "TextSource";      }
    virtual std::string getCategory() const  { return "Input";           }
    virtual CodeState getCodeState() const   { return CODE_STATE_STABLE; }

protected:
    virtual void setDescriptions() {
        setDescription("Provides a textport output with content defined by the corresponding property.");
    }

    virtual void process();

private:
    TextPort outport_;
    StringProperty text_;

    static const std::string loggerCat_; ///< category used in logging
};

}   //namespace

#endif
