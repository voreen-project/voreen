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

#include "plotdataexport.h"

#include "../datastructures/plotdata.h"
#include "../datastructures/plotcell.h"
#include "../datastructures/plotrow.h"
#include "../datastructures/plotfunction.h"

#include "voreen/core/voreenapplication.h"

#include "tgt/tgt_math.h"
#include <iostream>
#include <fstream>

namespace voreen {

const std::string PlotDataExport::loggerCat_("voreen.PlotDataExport");

PlotDataExport::PlotDataExport():
    PlotDataExportBase()
    , outputFile_("outputFile", "Save File", "Write CSV Output file", VoreenApplication::app()->getUserDataPath(),
        "Comma-Separated Files (*.csv);;Textfile (*.txt);;XML-File (*.xml);;HTML-File (*.htm *.html);;All TextFiles (*.csv *.xml *.txt *.htm *.html)",
        FileDialogProperty::SAVE_FILE,Processor::VALID)
    , rewrite_("rewrite","rewrite Export File",Processor::VALID)
    , continuousSave_("continuousSave", "Save continuously", false)
{
    outputFile_.onChange(MemberFunctionCallback<PlotDataExport>(this, &PlotDataExport::exportFile));
    rewrite_.onChange(MemberFunctionCallback<PlotDataExport>(this, &PlotDataExport::exportFile));
    addProperty(&outputFile_);
    addProperty(&rewrite_);
    addProperty(continuousSave_);

    outputFile_.setGroupID("FileSelection");
    rewrite_.setGroupID("FileSelection");
    continuousSave_.setGroupID("FileSelection");

    setPropertyGroupGuiName("FileSelection", "File-Selection");
}

Processor* PlotDataExport::create() const {
    return new PlotDataExport();
}

void PlotDataExport::process() {
    PlotDataExportBase::process();

    if (inPort_.hasChanged() && continuousSave_.get())
        exportFile();
}

void PlotDataExport::exportFile() {
    PlotDataExportBase::exportFile(outputFile_.get());
}

}
