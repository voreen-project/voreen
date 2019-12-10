/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#include "poicsvexport.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "voreen/core/io/serialization/serialization.h"
#include "voreen/core/io/serialization/serializable.h"

#include <fstream>
#include <sstream>
namespace voreen{

/**
 * Exports points in CSV file format
 */
POICSVExport::POICSVExport()
    : inport_(Port::INPORT, "inport", "POI Inport", false)
    , volumeInport_(Port::INPORT, "volumeINport", "Volume port", false)
    , fileProp_("fileProp", "Path to save", "Path for POI saving", ".", "Comma seperated values (*.csv)", FileDialogProperty::SAVE_FILE)
    , autoSave_("autoSave", "Auto save on change")
    , saveButton_("saveButton", "Save")
    , shouldSaveOnProcess_(false)
{
    addPort(inport_);
    addPort(volumeInport_);

    addProperty(fileProp_);
    addProperty(autoSave_);
    addProperty(saveButton_);

    ON_CHANGE_LAMBDA(saveButton_, [this]{
        shouldSaveOnProcess_ = true;
    });
}

Processor* POICSVExport::create() const
{
    return new POICSVExport();
}

std::string POICSVExport::getClassName() const
{
    return "POICSVExport";
}

std::string POICSVExport::getCategory() const
{
    return "Points of Interest";
}

void POICSVExport::setDescriptions()
{
    setDescription("Stores Points in csv file.");
    fileProp_.setDefaultValue("File to export the data.");
    autoSave_.setDescription("Stores on every change of the data.");
    saveButton_.setDescription("Store data now.");
}


void POICSVExport::process()
{
    if (!volumeInport_.getData())
        return;
    const POIList* list = inport_.getData();
    if (!list)
        return;

    if (shouldSaveOnProcess_ || autoSave_.get()){
        shouldSaveOnProcess_ = false;    
        std::ofstream file(fileProp_.get());
        file << "id,x,y,z,voxel x,voxel y,voxel z,group"<<std::endl;
        tgt::mat4 voxelToWorld = volumeInport_.getData()->getVoxelToWorldMatrix();
        tgt::mat4 worldToVoxel;
        voxelToWorld.invert(worldToVoxel);
        for(auto p: list->getPoints()){
            POIGroup g = list->getGroup(p.group_);
            if (!g.enabled_)
                continue;

            tgt::vec3 pos_vx = (worldToVoxel*tgt::vec4(p.position_, 1.0)).xyz();
            file << p.id_ <<","<<p.position_.x<<","<<p.position_.y<<","<<p.position_.z<<","
                 << pos_vx.x << "," << pos_vx.y << "," << pos_vx.z << "," 
                 << g.name_<< std::endl;
        }

    }
}

}
