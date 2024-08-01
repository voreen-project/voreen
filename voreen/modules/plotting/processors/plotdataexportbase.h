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

#ifndef VRN_PLOTDATAEXPORTBASE_H
#define VRN_PLOTDATAEXPORTBASE_H

#include "../ports/plotport.h"

#include "voreen/core/processors/processor.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/boolproperty.h"

namespace voreen {

class PlotData;
class PlotFunction;

class VRN_CORE_API PlotDataExportBase : public Processor {

public:
    PlotDataExportBase();
    ~PlotDataExportBase();
    virtual Processor* create() const;

    virtual std::string getCategory() const  { return "Export"; }
    virtual std::string getClassName() const { return "PlotDataExportBase"; }
    virtual CodeState getCodeState() const   { return CODE_STATE_TESTING; }
    virtual bool isEndProcessor() const;

protected:
    virtual void setDescriptions() {
        setDescription("");
    }
    virtual void process();
    virtual void deinitialize();

    void exportFile(std::string filename);

    PlotPort inPort_;

private:
    StringProperty separator_;
    BoolProperty valueBracket_;
    const PlotData* pData_;
    const PlotFunction* pPlotFunction_;

    static const std::string loggerCat_;
};

} //namespace
#endif

