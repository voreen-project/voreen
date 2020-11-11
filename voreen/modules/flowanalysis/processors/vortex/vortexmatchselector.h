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

#ifndef VRN_VORTEXMATCHSELECTOR_H
#define VRN_VORTEXMATCHSELECTOR_H

#include "voreen/core/processors/renderprocessor.h"
#include "voreen/core/ports/geometryport.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/properties/string/stringlistproperty.h"

#include "modules/flowanalysis/ports/vortexport.h"

namespace voreen {
class VortexMatchSelector : public RenderProcessor 	{
    enum class Visualization { eTree, eWorld };

public:
    VortexMatchSelector();

    Processor* create() const override
    {
        return new VortexMatchSelector();
    }
    std::string getClassName() const override
    {
        return "VortexMatchSelector";
    }
    std::string getCategory() const override
    {
        return "Vortex Processing";
    }

    bool isReady() const override
    {
        return _inportVortexCollection.isReady() && _outportRender.isReady();
    }

private:
    void process() override;
    void render( Visualization visualization );

    tgt::vec2 treeToClip( tgt::vec2 pos ) const;
    tgt::vec2 worldToClip( tgt::vec2 pos ) const;

    void onNewData();
    void onHoverEvent( tgt::MouseEvent* event );
    void onMouseEvent( tgt::MouseEvent* event );
    void onWheelEvent( tgt::MouseEvent* event );

    void updateProbabilities();
    void updateVortexInfos();
    void updateGeometry();

    void fillGroup( const int32_t group, const VortexCollection::VortexID& id );
    void fillSelection( const bool selected, const VortexCollection::VortexID& id );

    VortexCollectionPort _inportVortexCollection, _inportVortexCollectionMatch;
    RenderPort _outportRender;
    GeometryPort _outportGeometryClockwise, _outportGeometryCounterClockwise;

    std::unique_ptr<EventProperty<VortexMatchSelector>> _eventPropertyHover, _eventPropertyMouse, _eventPropertyWheel;

    FloatProperty _propertyMaximumMatchDistance;
    StringListProperty _propertySelectedMembers;
    IntProperty _propertyMinimumGroupSize;
    IntOptionProperty _propertyOrientation;
    IntIntervalProperty _propertyTimestepInterval;
    IntProperty _propertyCorelineLength;
    ButtonProperty _propertySelectAll, _propertyClearSelection;

    struct VortexInfo
    {
        tgt::vec2 positionTree = tgt::vec2::zero, positionWorld = tgt::vec2::zero;
        int32_t group = -1;
        float probability = 1.0f;
        bool visible = true, predecessor = false, highlighted = false, selected = false;
    };
    std::vector<std::vector<VortexInfo>> _vortexInfos;
    VortexCollection::VortexID _hoveredVortex = VortexCollection::VortexID::Invalid;

    struct GroupInfo
    {
        uint32_t members = 0, width = 0;
    };
    std::vector<GroupInfo> _groupInfos;

    float _zoomWorld = 0.95f;
    tgt::vec2 _offsetWorld = tgt::vec2::zero, _offsetTree = tgt::vec2::zero;
};

}

#endif // VRN_VORTEXMATCHSELECTOR_H
