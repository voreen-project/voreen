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

#include "fieldplotsave.h"

#include "voreen/core/voreenapplication.h"

#include "modules/hdf5/io/hdf5volumewriter.h"

#include "tgt/filesystem.h"

namespace voreen {

    const std::string FieldPlotSave::loggerCat_("voreen.viscontest2018.FieldPlotSave");

    FieldPlotSave::FieldPlotSave()
        : Processor()
        // ports
        , fieldPlotInport_(Port::INPORT, "fieldPlotInport", "Field Plot Input", false)
        // properties
        , filenameProp_("filenameprop", "Save File as", "Select file...", VoreenApplication::app()->getUserDataPath(), "field plot data (*.fpd)", FileDialogProperty::SAVE_FILE, Processor::INVALID_PATH)
        , saveButton_("saveButton", "Save")
        // members
        , saveFieldPlot_(false)
    {
        addPort(fieldPlotInport_);

        saveButton_.onChange(MemberFunctionCallback<FieldPlotSave>(this, &FieldPlotSave::saveFieldPlot));

        addProperty(filenameProp_);
        addProperty(saveButton_);
    }

    FieldPlotSave::~FieldPlotSave() {
    }

    bool FieldPlotSave::isReady() const {
        if (!isInitialized())
            return false;
        // only streamlineport must be ready
        if (!fieldPlotInport_.isReady())
            return false;

        return true;
    }

    void FieldPlotSave::invalidate(int inv) {
        Processor::invalidate(inv);

        if (inv == Processor::INVALID_PATH && isInitialized()) {
            saveFieldPlot_ = true;
        }
    }

    void FieldPlotSave::process() {
        if (saveFieldPlot_){
            saveFieldPlot();
            saveFieldPlot_ = false;
        }
    }

    //---------------------------------------------------------------------------------
    //      Callbacks
    //---------------------------------------------------------------------------------
    void FieldPlotSave::saveFieldPlot() {
        if (!isInitialized())
            return;
        if (!fieldPlotInport_.getData()) {
            LWARNING("no input field plot");
            return;
        }
        if (filenameProp_.get().empty()) {
            LWARNING("no filename specified");
            return;
        }
        std::string extension = FileSys.fileExtension(filenameProp_.get(), true);
        if (extension.compare("fpd")) {
            LWARNING("Selected file extension is not *.fpd");
            return;
        }

        try {
            HDF5VolumeWriter().write(filenameProp_.get(), fieldPlotInport_.getData()->getVolume());
        } catch(tgt::FileException e) {
            LERROR(e.what());
            filenameProp_.set("");
        }
    }

}   // namespace
