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

#include "poistorage.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"
namespace voreen {

POIStorage::POIStorage()
    : Processor()
    , cpPort_(Port::OUTPORT, "cpPort", "Coprocessors", true)
    , inport_(Port::INPORT, "inport", "POI inport", false)
    , outport_(Port::OUTPORT, "outport", "POI outport", true)
    , groupSelector_("groupSelector", "Selector of Groups", 1, 1)
    , removeGroup_("removeGroup", "Remove Group")
    , addGroup_("addGroup", "Add Group")
    , groupName_("groupName", "Name of Group")
    , groupEnabled_("groupEnabled", "Is Enabled")
    , groupColor_("groupColor", "Color of Group")
    , mouseOverPoint_("mouseOverPoint", "Mouseover Point", POI_NO_SUCH_POINT, -1, std::numeric_limits<int>::max())
    , clearButton_("clearButton", "Clear all Data")
    , enableInteraction_("isEnableInteraction", "Enable Interaction", true)
{
    addPort(cpPort_);
    addPort(inport_);
    addPort(outport_);

    groupName_.setInstantUpdate(false);
    addProperty(enableInteraction_);
    addProperty(mouseOverPoint_);
        mouseOverPoint_.setGroupID("active");
    setPropertyGroupGuiName("active", "Active Point");

    addProperty(groupSelector_);
        groupSelector_.setGroupID("group");
    addProperty(removeGroup_);
        removeGroup_.setGroupID("group");
    addProperty(addGroup_);
        addGroup_.setGroupID("group");
    addProperty(groupName_);
        groupName_.setGroupID("group");
    addProperty(groupEnabled_);
        groupEnabled_.setGroupID("group");
    addProperty(groupColor_);
        groupColor_.setGroupID("group");
    setPropertyGroupGuiName("group", "Group Configuration");

    addProperty(clearButton_);


    ON_CHANGE_LAMBDA(removeGroup_, [&]{
        std::string g = groupSelectorState_[groupSelector_.getSelectedRowIndex()];
        POIGroupID gid = list_.getGroupID(g);
        POIPointID mouseoverid = static_cast<POIPointID>(mouseOverPoint_.get());
        if (mouseoverid != POI_NO_SUCH_POINT && list_.getPointById(mouseoverid).group_ == gid)
            mouseOverPoint_.set(static_cast<int>(POI_NO_SUCH_POINT));
        list_.removeGroup(g);
        updateGroupSelector();
        updateGroupSettings();
        cpPort_.invalidate();
    });

    ON_CHANGE_LAMBDA(clearButton_, [this]{
        list_ = POIList();
        invalidatePOIs();
        updateGroupSelector();
        updateGroupSettings();
    });

    ON_CHANGE_LAMBDA(addGroup_, [this]{
        list_.addGroup("New Group");
        updateGroupSelector();
        for(int i = 0; i != (int)groupSelectorState_.size(); i++)
            if (groupSelectorState_[i] == "New Group")
                groupSelector_.setSelectedRowIndex(i);
        updateGroupSettings();
        cpPort_.invalidate();
    });

    std::function<void()> updateGroup=[&]{
        int row = groupSelector_.getSelectedRowIndex();
        if (row == -1 || row >= (int)groupSelectorState_.size())
            return;
        std::string g = groupSelectorState_[row];
        POIGroupID gid = list_.getGroupID(g);
        list_.setGroupColor(gid, groupColor_.get().xyz());
        list_.setGroupEnabled(gid, groupEnabled_.get());
        if (g != groupName_.get()){
            std::string gn = groupName_.get();
            std::vector<std::string> groups = list_.getGroupNames();
            bool found = false;
            for(auto g: groups){
                if (g == gn){
                    found = true;
                    break;
                }
            }
            if (found){
                LERRORC("POIStorage","Could not change group name to " << gn << ". Group with name already exists!");
            }else{
                list_.setGroupName(gid, groupName_.get());
                updateGroupSelector();
            }            
        }
            
        cpPort_.invalidate();
    };

    ON_CHANGE_LAMBDA(groupName_, updateGroup);
    ON_CHANGE_LAMBDA(groupColor_, updateGroup);
    ON_CHANGE_LAMBDA(groupEnabled_, updateGroup);
    ON_CHANGE_LAMBDA(groupSelector_, [&]{
        updateGroupSettings();
    });

    ON_CHANGE_LAMBDA(inport_, [&]{
        const POIList* list = inport_.getData();
        if (list){
            list_ = *list;
            updateGroupSelector();
            updateGroupSettings();
            //activePoint_.set(static_cast<int>(list_.getActivePoint()));
        }
    });
}


Processor* POIStorage::create() const {
    return new POIStorage();
}

std::string POIStorage::getClassName() const {
    return "POIStorage";
}

std::string POIStorage::getCategory() const {
    return "Points of Interest";
}

void POIStorage::setDescriptions() {
    setDescription("The POIStorage is the central processor for every network. It contains and manages "
                   "all data. The input port can be used to load data from different sources. If no data is "
                   "supplied in such a way, an empty list of points and groups is initialized. The Output "
                   "port is used to get the points as a whole object and for example to save them to the disk.");
    groupSelector_.setDescription("Used to select groups.");

    removeGroup_.setDescription("Remove selected group.");
    addGroup_.setDescription("Add a new group.");

    //activePoint_.setDescription("The active/selected point.");
    mouseOverPoint_.setDescription("The point under the mouse.");

    groupName_.setDescription("Name of the selected group.");
    groupEnabled_.setDescription("Enabled status of the active group");
    groupColor_.setDescription("Color of the selected group");
}

void POIStorage::process() 
{
    outport_.setData(&list_, false);
}

POIList& POIStorage::getPOIS()
{
    return list_;
}

void POIStorage::updateGroupSelector()
{
    int selectedRow = groupSelector_.getSelectedRowIndex();
    std::string selectedGroup = selectedRow != -1 && selectedRow < (int)groupSelectorState_.size()? groupSelectorState_[selectedRow]: "";

    int indexToSelect = -1;

    groupSelectorState_.clear();
    groupSelector_.reset();
    groupSelectorState_ = list_.getGroupNames();
    for(int i = 0; i != (int)groupSelectorState_.size(); i++){
        std::string group = groupSelectorState_[i];
        if (group == selectedGroup){
            indexToSelect = i;
        }
        std::vector<std::string> v;
        v.push_back(group);
        groupSelector_.addRow(v);
    }
    if (selectedRow != -1)
        groupSelector_.setSelectedRowIndex(indexToSelect);
    groupSelector_.invalidate();
}

void POIStorage::updateGroupSettings()
{
    int idx = groupSelector_.getSelectedRowIndex();
    if (idx == -1 || idx >= (int)groupSelectorState_.size()){
        removeGroup_.setReadOnlyFlag(true);
        groupName_.setVisibleFlag(false);
        groupEnabled_.setVisibleFlag(false);
        groupColor_.setVisibleFlag(false);
    }else{
        removeGroup_.setReadOnlyFlag(false);
        groupName_.setVisibleFlag(true);
        groupEnabled_.setVisibleFlag(true);
        groupColor_.setVisibleFlag(true);
        POIGroup group =list_.getGroup(groupSelectorState_[idx]);
        groupName_.set(group.name_);
        groupEnabled_.set(group.enabled_);
        groupColor_.set(tgt::vec4(group.color_, 1.0f));
    }
}

bool POIStorage::isReady() const
{
    return true;
}

void POIStorage::initialize() {
    Processor::initialize();
    groupSelector_.reset();
    updateGroupSelector();
    updateGroupSettings();
}

void POIStorage::invalidatePOIs()
{
    cpPort_.invalidate();
}

voreen::POIGroupID POIStorage::getActiveGroup()
{
    int idx = groupSelector_.getSelectedRowIndex();
    if (idx == -1 || idx >= (int)groupSelectorState_.size()){
        return POI_NO_SUCH_GROUP;
    }else{
        return list_.getGroupID(groupSelectorState_[idx]);
    }
}

void POIStorage::addPoint(tgt::vec3 position, POIGroupID group)
{
    list_.addPoint(position, group);
    invalidatePOIs();
}

void POIStorage::removePoint(POIPointID id)
{
    if (id == mouseOverPoint_.get()){
        mouseOverPoint_.set(static_cast<int>(POI_NO_SUCH_POINT));
    }
    list_.removePoint(id);
    invalidatePOIs();
}

const std::vector<POIPoint>& POIStorage::getPoints() const
{
    return list_.getPoints();
}

voreen::POIGroupID POIStorage::addGroup(std::string name, tgt::vec3 color/*=tgt::vec3::one*/, bool enabled/*=true*/)
{
    POIGroupID g = list_.addGroup(name, color, enabled);
    invalidatePOIs();
    return g;
}

void POIStorage::removeGroup(std::string name)
{
    list_.removeGroup(name);
    invalidatePOIs();
}

void POIStorage::removeGroup(POIGroupID gid)
{
    list_.removeGroup(gid);
    invalidatePOIs();
}

std::string POIStorage::getGroupName(POIGroupID gid) const
{
    return list_.getGroupName(gid);
}

std::vector<std::string> POIStorage::getGroups() const
{
    return list_.getGroupNames();
}

voreen::POIGroupID POIStorage::getGroupID(std::string name) const
{
    return list_.getGroupID(name);
}

voreen::POIGroup POIStorage::getGroup(POIGroupID id) const
{
    return list_.getGroup(id);
}

voreen::POIGroup POIStorage::getGroup(std::string name) const
{
    return list_.getGroup(name);
}

void POIStorage::setGroupColor(POIGroupID id, tgt::vec3 color)
{
    list_.setGroupColor(id, color);
    invalidatePOIs();
}

void POIStorage::setGroupEnabled(POIGroupID id, bool enabled)
{
    list_.setGroupEnabled(id, enabled);
    invalidatePOIs();
}

void POIStorage::setGroupName(POIGroupID id, std::string name)
{
    list_.setGroupName(id, name);
    invalidatePOIs();
}

voreen::POIPoint POIStorage::getPointById(POIPointID id) const
{
    return list_.getPointById(id);
}

void POIStorage::setPointById(POIPoint p)
{
    list_.setPointById(p);
    invalidatePOIs();
}

void POIStorage::deserialize(Deserializer& s)
{
    Processor::deserialize(s);
    mouseOverPoint_.set(static_cast<int>(POI_NO_SUCH_POINT));
}

voreen::POIPointID POIStorage::getMouseOverPoint()
{
    return static_cast<int>(mouseOverPoint_.get());
}

void POIStorage::setMouseOverPoint(POIPointID p)
{
    if (p != POI_NO_SUCH_POINT) list_.getPointById(p); // exception if invalid point!
    mouseOverPoint_.set(static_cast<int>(p));
    invalidatePOIs();
}

bool POIStorage::isInteractionEnabled() const
{
    return enableInteraction_.get();
}

void POIStorage::removeSelectedPoints() {
    list_.removeSelectedPoints();
    invalidatePOIs();
}

int POIStorage::selectedPointCount() const
{
    return list_.selectedPointCount();
}

std::vector<POIPoint> POIStorage::getSelectedPoints() const
{
    return list_.getSelectedPoints();
}

void POIStorage::clearSelection()
{
    list_.clearSelection();
    invalidatePOIs();
}

std::vector<std::string> POIStorage::getGroupNames() const
{
    return list_.getGroupNames();
}

void POIStorage::selectPoints( const std::vector<POIPoint> &points, SelectionMode sm )
{
    if (sm == ADD_TO_SELECTION){
        for(auto p: points){
            if (!list_.getGroup(p.group_).enabled_)
                continue;
            p.selected_ = true;
            list_.setPointById(p);
        }
    }else if (sm == TOGGLE_SELECTION){
        for(auto p: points){
            if (!list_.getGroup(p.group_).enabled_)
                continue;
            p.selected_ = !p.selected_;
            list_.setPointById(p);
        }
    }else if(sm == REPLACE_SELECTION){
        list_.clearSelection();
        for(auto p: points){
            if (!list_.getGroup(p.group_).enabled_)
                continue;
            p.selected_ = true;
            list_.setPointById(p);
        }
    }else if (sm == REMOVE_SELECTION){
        for(auto p: points){
            if (!list_.getGroup(p.group_).enabled_)
                continue;
            p.selected_ = false;
            list_.setPointById(p);
        }
    }
    invalidatePOIs();
}
} // namespace

