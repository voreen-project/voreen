/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#ifndef VRN_POILIST_H
#define VRN_POILIST_H
#include "voreen/core/io/serialization/serializable.h"
#include "voreen/core/ports/genericport.h"
#include "tgt/vector.h"
#include <stdint.h>
#include <vector>
#include <string>
#include <map>

namespace voreen {
typedef int64_t POIPointID;
typedef int POIGroupID;

const static POIGroupID POI_NO_SUCH_GROUP = 0;
const static POIPointID POI_NO_SUCH_POINT =  -1;

/**
 * Point of interest
 */
struct POIPoint :public Serializable{
    tgt::vec3 position_;
    POIPointID id_;
    POIGroupID group_;
    bool selected_;

    virtual void serialize(Serializer& s) const override;
    virtual void deserialize(Deserializer& s) override;
};

/**
 * Group of points of interest
 */
struct POIGroup :public Serializable{
    std::string name_;
    POIGroupID id_;
    tgt::vec3 color_;
    bool enabled_;
    
    virtual void serialize(Serializer& s) const override;
    virtual void deserialize(Deserializer& s) override;
};

class POIList :public Serializable{
public:
    POIList();
    /**
     * Adds a Point to the list
     */
    void addPoint(tgt::vec3 position, POIGroupID group);
    /**
     * Removes the given point and throws a exception if it does not exist
     */
    void removePoint(POIPointID id);
    /**
     * Unsorted list of Points
     * This does not copy the points and vector should NOT be changed in any case
     */
    const std::vector<POIPoint>& getPoints() const;
    /**
     * Get a copy of a Point by it's id or throws an exception
     * 
     */
    POIPoint getPointById(POIPointID id) const;

    /**
     * Creates a new group if it does not exist or
     * returns an already existing group with the name.
     * In that second case it does NOT set color and enabled
     */
    POIGroupID addGroup(std::string name, tgt::vec3 color, bool enabled=true);
    /**
     * Creates a new group if it does not exist or
     * returns an already existing group with the name.
     * In that second case it does not set enabled
     */
    POIGroupID addGroup(std::string name, bool enabled=true);
    /**
     * Removes a group and all points belonging to it
     */
    void removeGroup(std::string name);
    /**
     * Removes a group and all points belonging to it
     */
    void removeGroup(POIGroupID);
    /**
     * Gets group name from id or throws exception if it does not exist
     */
    std::string getGroupName(POIGroupID gid) const;
    /**
     * Gets the names of all groups
     */
    std::vector<std::string> getGroupNames() const;
    /**
     * Gets the id of a group by name or throws exception if it does not exist
     */
    POIGroupID getGroupID(std::string name) const;
    /**
     * Gets a copy of a group or throws exception if it does not exist
     */
    POIGroup getGroup(POIGroupID id) const;
    /**
     * Gets a copy of a group or throws exception if it does not exist
     */
    POIGroup getGroup(std::string name) const;
    /**
     * Sets the color of a group or throws exception if it does not exist
     */
    void setGroupColor(POIGroupID id, tgt::vec3 color);
    /**
     * enables/disables a group or throws exception if it does not exist
     */
    void setGroupEnabled(POIGroupID id, bool enabled);
    /**
     * Sets the name of a group or throws exception if it does not exist or the name is already used
     */
    void setGroupName(POIGroupID id, std::string name);

    /**
     * Updates the point (identified by it's id) in the group list
     */
    void setPointById(POIPoint p);

    /**
     * Get number of selected point by the user.
     */ 
    int selectedPointCount() const;
    /**
     * Get all selected points.
     * Warning: This is extremely slow
     */
    std::vector<POIPoint> getSelectedPoints() const;

    /**
     * Makes all points not selected
     */
    void clearSelection();

    /**
     * Delete all selected points
     */
    void removeSelectedPoints();

    virtual void serialize(Serializer& s) const override;

    virtual void deserialize(Deserializer& s) override;

private:
    int getGroupIndex(std::string name) const;
    int getGroupIndex(POIGroupID gid) const;
    int getGroupIndexExcp(std::string name) const;
    int getGroupIndexExcp(POIGroupID gid) const;
    int getPointIndexExcp(POIPointID id) const;

    std::vector<POIPoint> points_; ///< unordered set of points
    std::vector<POIGroup> groups_; ///< unordered set of groups
    POIGroupID lastUsedGroupID_;
    POIPointID lastUsedPointID_;
};

class POIListPort : public GenericPort<POIList> {
public:
    POIListPort(PortDirection direction, const std::string& name, const std::string& guiName = "", bool allowMultipleConnections = false,
        Processor::InvalidationLevel invalidationLevel = Processor::INVALID_RESULT)
        : GenericPort<POIList>(direction, name, guiName, allowMultipleConnections, invalidationLevel) {}

    virtual Port* create(PortDirection direction, const std::string& id, const std::string& guiName = "") const {return new POIListPort(direction,id);}
    virtual std::string getClassName() const {return "POIListPort";}
};

} 

#endif
