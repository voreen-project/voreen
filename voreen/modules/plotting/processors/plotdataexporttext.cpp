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

#include "plotdataexporttext.h"

#include "../datastructures/plotdata.h"
#include "../datastructures/plotcell.h"
#include "../datastructures/plotrow.h"
#include "../datastructures/plotfunction.h"

#include "voreen/core/voreenapplication.h"

#include "tgt/tgt_math.h"
#include <iostream>
#include <fstream>

namespace voreen {

const std::string PlotDataExportText::loggerCat_("voreen.PlotDataExportText");

PlotDataExportText::PlotDataExportText():
    PlotDataExportBase()
    , filenamePort_(Port::INPORT, "filename")
    , rewrite_("rewrite","rewrite Export File",Processor::VALID)
{
    rewrite_.onChange(MemberFunctionCallback<PlotDataExportText>(this, &PlotDataExportText::exportFile));
    addProperty(rewrite_);
    addPort(filenamePort_);
}

Processor* PlotDataExportText::create() const {
    return new PlotDataExportText();
}

void PlotDataExportText::exportFile() {
    PlotDataExportBase::exportFile(filenamePort_.getData());
}


}
