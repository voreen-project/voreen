/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2018 University of Muenster, Germany.                        *
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

#include "poilist.h"
#include <sstream>
#include "voreen/core/utils/exception.h"

namespace voreen {

POIList::POIList()
    : lastUsedGroupID_(POI_NO_SUCH_GROUP)
    , lastUsedPointID_(0)
{
}

void POIList::addPoint(tgt::vec3 position, POIGroupID group)
{
    POIPoint p;
    p.id_ = ++lastUsedPointID_;
    p.group_ = group;
    p.position_ = position;
    p.selected_ = false;
    points_.push_back(p);
}

void POIList::removePoint(POIPointID id)
{
    for(size_t i = 0; i != points_.size(); i++){
        if (points_[i].id_ == id){
            points_[i] = points_[points_.size()-1];
            points_.pop_back();
            return;
        }
    }
}

const std::vector<POIPoint>& POIList::getPoints() const
{
    return points_;
}

voreen::POIGroupID POIList::addGroup(std::string name, tgt::vec3 color, bool enabled)
{
    int idx = getGroupIndex(name);
    if (idx == -1){
        // group does not exist yet
        POIGroup group;
        group.color_ = color;
        group.enabled_ = enabled;
        group.name_ = name;
        group.id_ = ++lastUsedGroupID_;
        groups_.push_back(group);
        return group.id_;
    }else{
        return groups_[idx].id_;
    }
}

voreen::POIGroupID POIList::addGroup(std::string name, bool enabled)
{
    int idx = getGroupIndex(name);
    if (idx == -1){
        // group does not exist yet
        POIGroup group;
        
        group.enabled_ = enabled;
        group.name_ = name;
        group.id_ = ++lastUsedGroupID_;
        std::vector<tgt::vec3> colorMap;
        colorMap.push_back(tgt::vec3(255,0,0) / 255.f);
        colorMap.push_back(tgt::vec3(0,255,0) / 255.f);
        colorMap.push_back(tgt::vec3(0,0,255) / 255.f);
        colorMap.push_back(tgt::vec3(255,0,255) / 255.f);
        colorMap.push_back(tgt::vec3(0,255,255) / 255.f);
        colorMap.push_back(tgt::vec3(255,255,0) / 255.f);
        colorMap.push_back(tgt::vec3(255,100,20) / 255.f);
        colorMap.push_back(tgt::vec3(250,200,150) / 255.f);
        colorMap.push_back(tgt::vec3(150,200,250) / 255.f);
        colorMap.push_back(tgt::vec3(30,30,30) / 255.f);
        group.color_ = colorMap[lastUsedGroupID_%colorMap.size()];
        groups_.push_back(group);
        return group.id_;
    }else{
        return groups_[idx].id_;
    }
}

void POIList::removeGroup(std::string name)
{
    POIGroupID gid = getGroupID(name);
    size_t removed = 0;
    for(size_t i = 0; i < points_.size()-removed; i++){
        if (points_[i].group_ == gid){
            points_[i] = points_[points_.size()-1-removed];
            removed++;
            i--;
        }
    }
    points_.resize(points_.size()-removed);

    int idx = getGroupIndex(name);
    if (idx != -1){
        groups_[idx] = groups_[groups_.size()-1];
        groups_.pop_back();
    }
}

void POIList::removeGroup(POIGroupID gid)
{
    removeGroup(getGroupName(gid));
}

std::string POIList::getGroupName(POIGroupID gid) const
{
    int idx = getGroupIndex(gid);
    if (idx != -1){
        return groups_[idx].name_;
    }
    return "";
}

std::vector<std::string> POIList::getGroupNames() const
{
    std::vector<std::string> grouplist;
    for(auto group: groups_){
        grouplist.push_back(group.name_);
    }
    return grouplist;
}

POIGroupID POIList::getGroupID(std::string name) const
{
    int idx = getGroupIndex(name);
    if (idx != -1){
        return groups_[idx].id_;
    }else{
        return POI_NO_SUCH_GROUP;
    }
}

POIGroup POIList::getGroup(POIGroupID id) const
{
    int idx = getGroupIndexExcp(id);
    return groups_[idx];
}

POIGroup POIList::getGroup(std::string name) const
{
    int idx = getGroupIndexExcp(name);
    return groups_[idx];

}

int POIList::getGroupIndex(std::string name) const
{
    for(int i = 0; i != (int)groups_.size(); i++){
        if (groups_[i].name_ == name)
            return i;
    }
    return -1;
}

int POIList::getGroupIndex(POIGroupID gid) const
{
    for(int i = 0; i != (int)groups_.size(); i++){
        if (groups_[i].id_ == gid)
            return i;
    }
    return -1;
}

int POIList::getGroupIndexExcp(std::string name) const
{
    for(int i = 0; i != (int)groups_.size(); i++){
        if (groups_[i].name_ == name)
            return i;
    }
    throw VoreenException(std::string()+"POIList no group with name '"+ name+"'");
}

int POIList::getGroupIndexExcp(POIGroupID gid) const
{
    for(int i = 0; i != (int)groups_.size(); i++){
        if (groups_[i].id_ == gid)
            return i;
    }
    std::stringstream ss;
    ss << "POIList no group with id '"<< gid<<"'";
    throw VoreenException(ss.str());
}

void POIList::setGroupColor(POIGroupID id, tgt::vec3 color)
{
    int idx = getGroupIndexExcp(id);
    groups_[idx].color_ = color;
}

void POIList::setGroupEnabled(POIGroupID id, bool enabled)
{
    int idx = getGroupIndexExcp(id);
    groups_[idx].enabled_ = enabled;
}

void POIList::setGroupName(POIGroupID id, std::string name)
{
    int idx = getGroupIndexExcp(id);
    POIGroupID g = getGroupIndex(name);
    if (g != -1) throw VoreenException("POIList::setGroupName: Group name '"+name+"' is already used");
    groups_[idx].name_ = name;
}

void POIList::serialize(Serializer& s) const
{
    s.serialize("points", points_);
    s.serialize("groups", groups_);
    s.serialize("lastUsedGroupID", lastUsedGroupID_);
    s.serialize("lastUsedPointID", lastUsedPointID_);
}

void POIList::deserialize(Deserializer& s)
{
    s.deserialize("points", points_);
    s.deserialize("groups", groups_);
    s.deserialize("lastUsedGroupID", lastUsedGroupID_);
    s.deserialize("lastUsedPointID", lastUsedPointID_);
}

voreen::POIPoint POIList::getPointById(POIPointID id) const
{
    int idx = getPointIndexExcp(id);
    return points_[idx];
}

void POIList::setPointById(POIPoint p)
{
    int idx = getPointIndexExcp(p.id_);
    points_[idx] = p;
}

int POIList::getPointIndexExcp(POIPointID id) const
{
    for(int i = 0; i != (int)points_.size(); i++){
        if (points_[i].id_ == id)
            return i;
    }
    std::stringstream ss;
    ss << "POIList no point with id '"<< id<<"'";
    throw VoreenException(ss.str());
}

void POIList::removeSelectedPoints()
{
    size_t removed = 0;
    for(size_t i = 0; i < points_.size()-removed; i++){
        if (points_[i].selected_){
            points_[i] = points_[points_.size()-1-removed];
            removed++;
            i--;
        }
    }
    points_.resize(points_.size()-removed);
}

int POIList::selectedPointCount() const
{
    int n = 0;
    for(int i = 0; i != (int)points_.size(); i++){
        n += points_[i].selected_ == true;
    }
    return n;
}

std::vector<POIPoint> POIList::getSelectedPoints() const
{
    std::vector<POIPoint> selectedPoints;
    for(int i = 0; i != (int)points_.size(); i++){
        if (points_[i].selected_){
            selectedPoints.push_back(points_[i]);
        }
    }
    return selectedPoints;
}

void POIList::clearSelection()
{
    for(auto &p: points_)
        p.selected_ = false;
}

void POIPoint::serialize(Serializer& s) const
{
    s.serialize("position", position_);
    s.serialize("pointid", id_); // id is reserved in serialization
    s.serialize("group", group_);
    s.serialize("selected", selected_);
}

void POIPoint::deserialize(Deserializer& s)
{
    s.deserialize("position", position_);
    s.deserialize("pointid", id_); // id is reserved in serialization
    s.deserialize("group", group_);
    s.optionalDeserialize("selected", selected_, false);
}

void POIGroup::serialize(Serializer& s) const
{
    s.serialize("name", name_);
    s.serialize("groupid", id_);// id is reserved in serialization
    s.serialize("enabled", enabled_);
    s.serialize("color", color_);
}

void POIGroup::deserialize(Deserializer& s)
{
    s.deserialize("name", name_);
    s.deserialize("groupid", id_);// id is reserved in serialization
    s.deserialize("enabled", enabled_);
    s.deserialize("color", color_);
}
}