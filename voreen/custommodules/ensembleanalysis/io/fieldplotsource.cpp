/***********************************************************************************
*                                                                                 *
* Voreen - The Volume Rendering Engine                                            *
*                                                                                 *
* Copyright (C) 2005-2017 University of Muenster, Germany.                        *
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

#include "fieldplotsource.h"

#include "voreen/core/voreenapplication.h"

#include "modules/hdf5/io/hdf5volumereader.h"

#include "tgt/filesystem.h"

namespace voreen {

const std::string FieldPlotSource::loggerCat_("voreen.viscontest2018.FieldPlotSource");

FieldPlotSource::FieldPlotSource()
    : Processor()
    // ports
    , fieldPlotOutport_(Port::OUTPORT, "fieldPlotOutport", "Field Plot Output", false)
    // properties
    , filenameProp_("filenameprop", "Load FPP File from", "Select file...", VoreenApplication::app()->getUserDataPath(), "field plot data (*.fpd)", FileDialogProperty::OPEN_FILE, Processor::INVALID_PATH)
    , loadButton_("loadButton", "Load")
    // members
    , loadFieldPlot_(false)
{
    addPort(fieldPlotOutport_);

    loadButton_.onChange(MemberFunctionCallback<FieldPlotSource>(this, &FieldPlotSource::loadFieldPlot));
    addProperty(filenameProp_);
    addProperty(loadButton_);
}

FieldPlotSource::~FieldPlotSource() {
}

bool FieldPlotSource::isReady() const {
    if (!isInitialized())
        return false;

    if (!fieldPlotOutport_.isReady())
        return false;

    return true;
}

void FieldPlotSource::invalidate(int inv) {
    Processor::invalidate(inv);

    if (inv == Processor::INVALID_PATH && isInitialized()) {
        loadFieldPlot_ = true;
    }
}

void FieldPlotSource::process() {
    if (loadFieldPlot_){
        loadFieldPlot();
        loadFieldPlot_ = false;
    }
}

//---------------------------------------------------------------------------------
//      Callbacks
//---------------------------------------------------------------------------------
void FieldPlotSource::loadFieldPlot() {
    if (!isInitialized())
        return;

    fieldPlotOutport_.setData(nullptr);

    if (filenameProp_.get().empty()) {
        LWARNING("no filename specified");
        return;
    }

    try {

        VolumeList* list = HDF5VolumeReader().read(filenameProp_.get());
        tgtAssert(list->size() == 1, "Only one volume expected");

        FieldPlotData* plotData = new FieldPlotData(static_cast<Volume*>(list->first()));
        fieldPlotOutport_.setData(plotData, true);

        // Delete volumes, since plot data creates a local copy.
        delete list->first();
        delete list;

    } catch(tgt::FileException e) {
        LERROR(e.what());
        filenameProp_.set("");
    }
}

}   // namespace
