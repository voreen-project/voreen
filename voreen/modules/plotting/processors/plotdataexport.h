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

#ifndef VRN_PLOTDATAEXPORT_H
#define VRN_PLOTDATAEXPORT_H

#include "../ports/plotport.h"

#include "plotdataexportbase.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/boolproperty.h"

namespace voreen {

class PlotData;
class PlotFunction;

class VRN_CORE_API PlotDataExport : public PlotDataExportBase {

public:
    PlotDataExport();
    virtual Processor* create() const;

    virtual std::string getCategory() const  { return "Export";          }
    virtual std::string getClassName() const { return "PlotDataExport";  }
    virtual CodeState getCodeState() const   { return CODE_STATE_STABLE; }

protected:
    virtual void setDescriptions() {
        setDescription("This Processor export a PlotDataObject to CSV-, TXT-, XML- or HTML-File.");
    }

    virtual void process();

private:
    void exportFile();

    FileDialogProperty outputFile_;   ///< filename of the file the value information should be saved in
    ButtonProperty rewrite_;
    BoolProperty continuousSave_;

    static const std::string loggerCat_;
};

} //namespace
#endif // VRN_PLOTDATAEXPORT_H

