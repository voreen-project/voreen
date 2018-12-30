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

#include "poisource.h"

#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "voreen/core/io/serialization/serializable.h"

#include <fstream>
namespace voreen{

POISource::POISource()
    : outport_(Port::OUTPORT, "outport", "POI Outport", true)
    , fileProp_("fileProp", "Path to load", "Path for POI loading", ".", "Voreen Points of interesste (*.vpoi)", FileDialogProperty::OPEN_FILE)
    , autoLoad_("autoLoad", "Auto load on change")
    , loadButton_("loadButton", "Load")
    , shouldLoadOnProcess_(false)
{
    addPort(outport_);

    addProperty(fileProp_);
    addProperty(autoLoad_);
    addProperty(loadButton_);

    ON_CHANGE_LAMBDA(loadButton_, [this]{
        shouldLoadOnProcess_ = true;
    });
}

Processor* POISource::create() const
{
    return new POISource();
}

std::string POISource::getClassName() const
{
    return "POISource";
}

std::string POISource::getCategory() const
{
    return "Points of Interest";
}

void POISource::setDescriptions()
{
    setDescription("The POISource is used to load data in the vpoi format.");
    fileProp_.setDescription("File to load data from.");
    autoLoad_.setDescription("Automatical data loading");
    loadButton_.setDescription("Load now.");
}

void POISource::process()
{
    try{
        if (autoLoad_.get() || shouldLoadOnProcess_){
            shouldLoadOnProcess_ = false;
            std::string fileName = fileProp_.get();
            POIList *pois;
            std::ifstream file(fileName);
            XmlDeserializer d(fileName);
            d.setUseAttributes(true);
            d.read(file);
            d.deserialize("POIList", pois);
            outport_.setData(pois, true);
        }
    }catch(...){
        LERROR("Could not load POIS");
    }
}
}