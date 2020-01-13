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

#ifndef VRN_POISTORAGE_H
#define VRN_POISTORAGE_H

#include <stdint.h>

#include "voreen/core/processors/processor.h"
#include "voreen/core/ports/genericcoprocessorport.h"
#include "voreen/core/properties/string/stringtableproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/colorproperty.h"
#include "../datastructures/poilist.h"

namespace voreen {

/**
 * Processor that contains all the data for the POIPoints and the POIGroups
 */
class VRN_CORE_API POIStorage : public Processor {

public:
    POIStorage();

    virtual Processor* create() const;
    virtual std::string getClassName() const;
    virtual std::string getCategory() const;

    /**
     * Invalidate coprocessors
     */
    void invalidatePOIs();
    /**
     * Select the group that is selected in the gui
     */
    POIGroupID getActiveGroup();

    // see ../datastructures/poilist.h
    // these are just wrappers that handle invalidation
    void addPoint(tgt::vec3 position, POIGroupID group);
    void removePoint(POIPointID id);
    POIPoint getPointById(POIPointID id) const;
    const std::vector<POIPoint>& getPoints() const;
    void removeSelectedPoints();
    POIGroupID addGroup(std::string name, tgt::vec3 color=tgt::vec3::one, bool enabled=true);
    void removeGroup(std::string name);
    void removeGroup(POIGroupID);
    std::string getGroupName(POIGroupID gid) const;
    std::vector<std::string> getGroupNames() const;
    std::vector<std::string> getGroups() const;
    POIGroupID getGroupID(std::string name) const;
    POIGroup getGroup(POIGroupID id) const;
    POIGroup getGroup(std::string name) const;
    void setGroupColor(POIGroupID id, tgt::vec3 color);
    void setGroupEnabled(POIGroupID id, bool enabled);
    void setGroupName(POIGroupID id, std::string name);
    void setPointById(POIPoint p);
    int selectedPointCount() const;
    std::vector<POIPoint> getSelectedPoints() const;
    void clearSelection();

    POIPointID getMouseOverPoint();
    void setMouseOverPoint(POIPointID p);

    bool isInteractionEnabled() const;

    enum SelectionMode{
        ADD_TO_SELECTION  = 0,
        REPLACE_SELECTION = 1,
        TOGGLE_SELECTION  = 2,
        REMOVE_SELECTION  = 3
    };
    void selectPoints(const std::vector<POIPoint> &points, SelectionMode sm);


    /// @see Serializable::deserialize
    virtual void deserialize(Deserializer& s) override;
protected:
    virtual void setDescriptions();
    virtual void process();
    virtual bool isReady() const override;
    POIList& getPOIS();
    void updateGroupSelector();
    void updateGroupSettings();
    void initialize() override;
private:
    POIListPort inport_;
    POIListPort outport_;
    
    GenericCoProcessorPort<POIStorage> cpPort_;
    StringTableProperty groupSelector_;
    std::vector<std::string> groupSelectorState_;

    ButtonProperty removeGroup_;
    ButtonProperty addGroup_;

    IntProperty mouseOverPoint_;

    StringProperty groupName_;
    BoolProperty groupEnabled_;
    ColorProperty groupColor_;

    BoolProperty enableInteraction_;

    ButtonProperty clearButton_;

    POIList list_;
};

} // namespace

#endif // VRN_SAMPLEPROCESSOR_H
