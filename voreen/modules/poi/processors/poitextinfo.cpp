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

#include "poitextinfo.h"
#include <map>
#include <set>
#include <sstream>

namespace voreen {

namespace{
std::string i2s(int i){
    std::stringstream ss;
    ss << i;
    return ss.str();
}
};

POITextInfo::POITextInfo()
    : Processor()
    , cpPort_(Port::INPORT, "cpPort", "Coprocessors", false)
    , textport_(Port::OUTPORT, "textport", "Textport")
    , showGroupInfo_("showGroupInfo", "Show group info", true)
    , textInfoFormat_("textInfoFormat", "Info Text Format", "")
{
    addPort(cpPort_);
    addPort(textport_);

    addProperty(showGroupInfo_);
    addProperty(textInfoFormat_);
        textInfoFormat_.addPlaceHolder("Id of active point", textInfoFormat_.makePlaceHolder("activeid"));
        textInfoFormat_.addPlaceHolder("Group of active point", textInfoFormat_.makePlaceHolder("activegroup"));
        textInfoFormat_.addPlaceHolder("Id of mouse point", textInfoFormat_.makePlaceHolder("mouseoverid"));
        textInfoFormat_.addPlaceHolder("Group of mouse point", textInfoFormat_.makePlaceHolder("mouseovergroup"));
        textInfoFormat_.addPlaceHolder("Distance", textInfoFormat_.makePlaceHolder("distance"));
        textInfoFormat_.setGroupID("info");
}

Processor* POITextInfo::create() const {
    return new POITextInfo();
}

std::string POITextInfo::getClassName() const {
    return "POITextInfo";
}

std::string POITextInfo::getCategory() const {
    return "Points of Interest";
}

void POITextInfo::setDescriptions() {
    setDescription("Minimal sample processor that appends a user-defined prefix to a given text.");
}

void POITextInfo::process() {
    POIStorage* storage = cpPort_.getConnectedProcessor();
    std::stringstream ss;
    if (showGroupInfo_.get()){
        std::map<POIGroupID, int> groups;
        for(auto g: storage->getGroups()){
            POIGroupID gid = storage->getGroupID(g);
            groups[gid] = 0;
        }
        for (auto p: storage->getPoints())
        {
            groups[p.group_]++;
        }
        std::vector<std::string> grouplist = storage->getGroups();
        std::sort(grouplist.begin(), grouplist.end());
    
        for (auto g: grouplist){
            if (storage->getGroup(g).enabled_){
                ss << g << " has " << groups[storage->getGroupID(g)] << " Elements" << std::endl;
            }
        }
    }

    POIPointID mouseoverid = storage->getMouseOverPoint();

    // Create info text
    std::map<std::string, std::string> infoReplacements;
    
    POIPoint mouseover_point;
    
    auto selInfo = selectionInfo();
    infoReplacements["activeid"] = std::get<0>(selInfo);
    infoReplacements["activegroup"] = std::get<1>(selInfo);


    if (mouseoverid != POI_NO_SUCH_POINT){
        mouseover_point = storage->getPointById(mouseoverid); 
        infoReplacements["mouseoverid"] = i2s(static_cast<int>(mouseoverid));
        infoReplacements["mouseovergroup"] = storage->getGroupName(mouseover_point.group_);
    }else{
        infoReplacements["mouseoverid"] = "-1";
        infoReplacements["mouseovergroup"] = "";
    }

    POIPointID activeid = POI_NO_SUCH_POINT;
    POIPoint active_point;
    if (storage->selectedPointCount() == 1){
        active_point = storage->getSelectedPoints().at(0);
        activeid = active_point.id_;
    }
    if (activeid != POI_NO_SUCH_POINT && mouseoverid != POI_NO_SUCH_POINT){
        float dist = tgt::distance(active_point.position_, mouseover_point.position_);
        infoReplacements["distance"] = formatSpatialLength(dist);
    }else{
        infoReplacements["distance"] = "";
    }
    std::string infotext = textInfoFormat_.replacePlaceHoldersInText(infoReplacements);

    ss << infotext;

    textport_.setData(ss.str());
}

std::tuple<std::string, std::string> POITextInfo::selectionInfo()
{
    std::ostringstream pointInfo;
    POIStorage* storage = cpPort_.getConnectedProcessor();
    auto selectedPoints = storage->getSelectedPoints();
    if (selectedPoints.size() == 0){
        pointInfo << "No point selected";
    }else if (selectedPoints.size() < 5){
        for(int i = 0; i != (int)selectedPoints.size(); i++){
            if (i != 0){
                pointInfo << ", ";
            }
            pointInfo << selectedPoints[i].id_;
        }
    }else{
        pointInfo << selectedPoints.size() << " points selected";
    }

    std::set<POIGroupID> gids;
    for(auto p: selectedPoints){
        gids.insert(p.group_);
    }
    std::ostringstream groupInfo;
    if (gids.size() == 0){
        groupInfo << "No groups selected";
    }else if (gids.size() < 4){
        std::vector<POIGroupID> selectedGroups(gids.begin(), gids.end());
        for(int i = 0; i != (int)selectedGroups.size(); i++){
            if (i != 0){
                groupInfo << ", ";
            }
            groupInfo << storage->getGroupName(selectedGroups[i]);
        }
    }else{
        groupInfo << gids.size() << " Groups selected";
    }
    return std::tuple<std::string, std::string>(pointInfo.str(), groupInfo.str());
}

} // namespace
