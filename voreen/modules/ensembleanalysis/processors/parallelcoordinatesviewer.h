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

#ifndef VRN_PARALLELCOORDINATESVIEWER_H
#define VRN_PARALLELCOORDINATESVIEWER_H

#include "../ports/parallelcoordinatesaxesport.h"
#include "../properties/parallelcoordinatesselectionproperty.h"
#include "../properties/parallelcoordinatessectionsproperty.h"

#include "voreen/core/processors/renderprocessor.h"
#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/properties/numeric/intervalproperty.h"
#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"

namespace voreen {
class ParallelCoordinatesViewer : public RenderProcessor {
public:
    ParallelCoordinatesViewer();

    virtual Processor* create() const override { return new ParallelCoordinatesViewer(); }
    virtual std::string getClassName() const override { return "ParallelCoordinatesViewer"; }
    virtual std::string getCategory() const override { return "ParallelCoordinates"; }
    virtual CodeState getCodeState() const override { return CODE_STATE_EXPERIMENTAL; }

private:
    virtual void setDescriptions() {
        setDescription("Viewer for parallel coordinates created by ParallelCoordinatesCreator");
        _propertySelectedMember.setDescription("Currently displayed member");
        _propertyVisualizationMode.setDescription("Visualization mode for currently displayed member:<br>"
                                                  "<strong>Single Time Step</strong>: Each axis represents a field<br>"
                                                  "<strong>All Time Steps</strong>: Each axis represents a time step<br>");
        _propertySelectedTimestep.setDescription("If single Time Step visualization mode is selected, "
                                                 "this property determines the time step to be displayed");
        _propertySelectedTimestep.setDescription("If all Time Steps are visualized, this property determines<br>"
                                                 "the field to be displayed");
        _propertyDensityBlending.setDescription("Enable to use a density visualization for rendered lines");
        _propertyDensityVisibleSamples.setDescription("Adjust this property to enhance line density contrast of visible samples");
        _propertyDensitySelectedSamples.setDescription("Adjust this property to enhance line density contrast of selected samples");

        _propertyTransFuncField.front().setDescription("Assigns the respective transfer function to the selected field");
        _propertyTransFunc.front().setDescription("Transfer function to be used for the selected field");
    }

    virtual void initialize() override;
    virtual void deinitialize() override;
    virtual void process() override;

    virtual void serialize(Serializer& s) const;
    virtual void deserialize(Deserializer& s);

    void onNewInportData();
    void onHoverEvent( tgt::MouseEvent* event );
    void onMouseEvent( tgt::MouseEvent* event );

    void updateSampleStates();
    void updateTransferFunction( size_t index );

    float screenAxisToCoordinate( size_t axis ) const;

    // Ports
    ParallelCoordinatesAxesPort _axesport;
    RenderPort _renderport;

    // Event Properties
    std::unique_ptr<EventProperty<ParallelCoordinatesViewer>> _eventPropertyHover, _eventPropertyMouse;

    // Properties
    IntOptionProperty _propertySelectedMember;
    IntOptionProperty _propertyVisualizationMode;

    IntProperty _propertySelectedTimestep;
    IntOptionProperty _propertySelectedField;
    FloatIntervalProperty _propertyFieldInterval;

    BoolProperty _propertyDensityBlending;
    IntProperty _propertyDensityVisibleSamples;
    IntProperty _propertyDensitySelectedSamples;

    std::array<IntOptionProperty, 4> _propertyTransFuncField;
    std::array<TransFunc1DKeysProperty, 4> _propertyTransFunc;
    std::array<TransFunc1DKeysProperty, 4> _propertyTransFuncIntern;

    ParallelCoordinatesSelectionProperty _propertySelectedSamples;
    ParallelCoordinatesSectionsProperty _propertySections;

    // OpenGL
    GLuint _shaderProgram, _vertexArray, _indexBuffer;
    std::vector<GLuint> _indexBufferVec;
    std::vector<tgt::vec3> _uniformBufferVec;

    // Interaction Data
    enum class Interaction { eNone, eMoving, eSelection } _interaction = Interaction::eNone;
    size_t _hoveredScreenAxis = ~0u;
    tgt::vec2 _activeSection;
    std::vector<std::list<tgt::vec2>> _sections;
    std::vector<bool> _samplesVisiblity;
    std::vector<int> _samplesSelection;

    // Constants
    static const float X_LIMIT, Y_LIMIT;
};

}

#endif // VRN_PARALLELCOORDINATESVIEWER_H