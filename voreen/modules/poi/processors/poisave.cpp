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

#include "poisave.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "voreen/core/io/serialization/serialization.h"
#include "voreen/core/io/serialization/serializable.h"

#include <fstream>
#include <sstream>
namespace voreen{

POISave::POISave()
    : inport_(Port::INPORT, "inport", "POI Inport", false)
    , fileProp_("fileProp", "Path to save", "Path for POI saving", ".", "Voreen Points of interesste (*.vpoi)", FileDialogProperty::SAVE_FILE)
    , autoSave_("autoSave", "Auto save on change")
    , saveSelection_("saveSelection", "Save selection", false)
    , saveButton_("saveButton", "Save")
    , shouldSaveOnProcess_(false)
{
    addPort(inport_);

    addProperty(fileProp_);
    addProperty(autoSave_);
    addProperty(saveSelection_);
    addProperty(saveButton_);

    ON_CHANGE_LAMBDA(saveButton_, [this]{
        shouldSaveOnProcess_ = true;
    });
}

Processor* POISave::create() const
{
    return new POISave();
}

std::string POISave::getClassName() const
{
    return "POISave";
}

std::string POISave::getCategory() const
{
    return "Points of Interest";
}

void POISave::setDescriptions()
{
    setDescription("The POISave processor is used to store the data from files in the vpoi format "
                   "If the Auto save on change property is checked, it will continiously save the data "
                   "to the selected file on every change. Otherwise it only stores the data when the Save "
                   "button is pressed.");
    fileProp_.setDefaultValue("File to store the data.");
    autoSave_.setDescription("Stores on every change of the data.");
    saveButton_.setDescription("Store data now.");
}

void POISave::process()
{
    const POIList* list = inport_.getData();
    if (!list)
        return;

    POIList newlist = *list;
    if (!saveSelection_.get()){
        newlist.clearSelection();
        tgtAssert(newlist.selectedPointCount() == 0, "Unexpected selected points found!");
    }

    if (shouldSaveOnProcess_ || autoSave_.get()){
        shouldSaveOnProcess_ = false;    
        XmlSerializer s;
        s.setUseAttributes(true);
        s.serialize("POIList", &newlist);
        std::stringstream ss;
        s.write(ss);
        std::ofstream file(fileProp_.get());
        file << ss.str();
    }
}
}
