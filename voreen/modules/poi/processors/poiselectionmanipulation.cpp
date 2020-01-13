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

#include "poiselectionmanipulation.h"
#include "voreen/core/datastructures/callback/lambdacallback.h"
#include "voreen/core/io/serialization/serialization.h"
#include "voreen/core/io/serialization/serializable.h"

#include <fstream>
#include <sstream>
namespace voreen{
POISelectionManipulation::POISelectionManipulation()
    : cpPort_(Port::INPORT, "cpPort", "Coprocessors", false)
    , groupSelector_("groupSelector", "Selector of Groups", 1, 1)
    , moveToGroup_("moveToGroup", "Move selected points to group")
    , removeButton_("removeButton", "Delete points")
{
    addPort(cpPort_);
    addProperty(groupSelector_);
    addProperty(moveToGroup_);
    addProperty(removeButton_);

    ON_CHANGE(cpPort_, POISelectionManipulation, updateGroupSelector);
    ON_CHANGE_LAMBDA(removeButton_, [&]{
        POIStorage* storage = cpPort_.getConnectedProcessor();
        storage->removeSelectedPoints();
    });
    ON_CHANGE_LAMBDA(moveToGroup_, [&]{
        POIStorage* storage = cpPort_.getConnectedProcessor();
        std::string g = groupSelectorState_[groupSelector_.getSelectedRowIndex()];
        POIGroupID gid = storage->getGroupID(g);

        for (auto p: storage->getSelectedPoints()){
            p.group_ = gid;
            storage->setPointById(p);
        }
    });
}

Processor* POISelectionManipulation::create() const
{
    return new POISelectionManipulation();
}

std::string POISelectionManipulation::getClassName() const
{
    return "POISelectionManipulation";
}

std::string POISelectionManipulation::getCategory() const
{
    return "Points of Interest";
}

void POISelectionManipulation::setDescriptions()
{
    setDescription("Processors that allows working with the selected points.");
    groupSelector_.setDescription("Allows selection of a group to move the points into.");
    moveToGroup_.setDescription("Moves the selected points in the choosen group.");
    removeButton_.setDescription("Deletes all selected points permanantly.");
}

void POISelectionManipulation::process()
{

}

void POISelectionManipulation::updateGroupSelector()
{
    POIStorage* storage = cpPort_.getConnectedProcessor();
    if (!storage)
        return;
    int selectedRow = groupSelector_.getSelectedRowIndex();
    std::string selectedGroup = selectedRow != -1 && selectedRow < groupSelectorState_.size()? groupSelectorState_[selectedRow]: "";

    int indexToSelect = -1;

    groupSelectorState_.clear();
    groupSelector_.reset();
    groupSelectorState_ = storage->getGroupNames();
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

}

