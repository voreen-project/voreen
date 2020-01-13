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

#ifndef VRN_TEXTSERIESSOURCE_H
#define VRN_TEXTSERIESSOURCE_H

#include "voreen/core/processors/processor.h"

#include "voreen/core/ports/textport.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/filedialogproperty.h"

namespace voreen {

class TextSeriesSource : public Processor {
public:
    TextSeriesSource();
    virtual Processor* create() const;

    std::string getCategory() const     { return "Data Source";         }
    std::string getClassName() const    { return "TextSeriesSource";    }
    CodeState getCodeState() const      { return CODE_STATE_OBSOLETE;   }

protected:
    virtual void setDescriptions() {
        setDescription("");
    }

    void process();

    void openTextFile();
    void loadText();

    FileDialogProperty filename_;
    IntProperty step_;
    TextPort outport_;

    std::vector<std::string> texts_;
};

} // namespace

#endif // VRN_TEXTSERIESSOURCE_H
