/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2024 University of Muenster, Germany,                        *
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

#ifndef VRN_PLOTPROCESSOR_H
#define VRN_PLOTPROCESSOR_H

#include "voreen/core/processors/renderprocessor.h"


#include "../ports/plotport.h"

#include "../properties/plotdataproperty.h"
#include "../properties/plotentitiesproperty.h"
#include "../properties/plotselectionproperty.h"

#include "voreen/core/properties/eventproperty.h"
#include "voreen/core/properties/stringproperty.h"
#include "voreen/core/properties/intproperty.h"
#include "voreen/core/properties/floatproperty.h"
#include "voreen/core/properties/vectorproperty.h"
#include "voreen/core/properties/boolproperty.h"
#include "voreen/core/properties/cameraproperty.h"
#include "voreen/core/properties/filedialogproperty.h"
#include "voreen/core/properties/buttonproperty.h"
#include "voreen/core/properties/optionproperty.h"
#include "voreen/core/properties/colorproperty.h"

#include "../datastructures/plotfunction.h"
#include "../datastructures/plotdata.h"
#include "../utils/plotlibrary/plotlibrary.h"
#include "../utils/plotlibrary/plotlibraryopengl.h"
#include "../utils/plotlibrary/plotlibrarysvg.h"
#include "../utils/plotlibrary/plotlibrarylatex.h"
#include "../interaction/plotpickingmanager.h"

namespace voreen {

class PlotCameraInteractionHandler;
class CameraInteractionHandler;

class VRN_CORE_API PlotProcessor : public RenderProcessor {

public:
    PlotProcessor(PlotEntitySettings::Entity entity, bool threeDimensional);
    virtual ~PlotProcessor();

    virtual std::string getCategory() const { return "Plotting"; }

    void switchPlotLibrary();

    void setFile(PlotLibraryFileBase::PlotLibraryFileType fileType, const std::string& filename);

    std::string getXLabel() const;
    std::string getYLabel() const;
    std::string getZLabel() const;

protected:
    virtual void process();
    virtual void initialize();

    /// Adds Z-Label properties (not added by default)
    void addZLabelProperties();
    /// renders the plot
    virtual void render() = 0;
    /// renders the data
    virtual void renderData() = 0;
    /// renders axes including scales
    virtual void renderAxes() = 0;
    /// set the status of the plotlib
    virtual void setPlotStatus() = 0;
    /// reads plotdata/plotfunction from inport
    virtual void readFromInport() = 0;
    /// calculates the domains from the data
    virtual void calcDomains() = 0;
    /// toggles visibility of some properties
    virtual void toggleProperties();
    /// creates plotlabels using the PlotSelectionProperty
    virtual void createPlotLabels();
    /// selects a plotdata from plotfunction
    virtual void selectDataFromFunction();

    // functions called by mouse events
    virtual void select(tgt::MouseEvent* e, bool highlight, bool label, bool zoom, bool additive);
    virtual void highlight(tgt::MouseEvent* e);     ///< replacing highlighting
    virtual void label(tgt::MouseEvent* e);         ///< replacing zooming
    virtual void zoom(tgt::MouseEvent* e);          ///< replacing labelling
    //virtual void unselect(tgt::MouseEvent* e, bool highlight, bool label, bool zoom);
    virtual void highlightAdditive(tgt::MouseEvent* e);     ///< additive highlighting
    virtual void zoomAdditive(tgt::MouseEvent* e);          ///< additive zooming
    virtual void labelAdditive(tgt::MouseEvent* e);         ///< additive labelling
    virtual void updateMousePosition(tgt::MouseEvent* e);
    virtual void enterExitEvent(tgt::MouseEvent* e);

    /// renders the region selected by the mouse (2d)
    void renderSelectedRegion();
    /// renders the plot label
    void renderPlotLabel();
    /// render the selection planes (3d)
    void renderSelectionPlanes();
    /// starts 3d selection, returns true iff we are actually in selection mode (selectionEdgeIndex_ != -1)
    bool startSelection3D(tgt::dvec2 startPosition);
    /// calculated selectEnd_, returns true, if we are in selection mode
    bool calculateSelectEnd(tgt::dvec2 mousePressPosition);
    /// returns selected region specified by selectionEdgeIndex, selectStart_ and selectEnd_
    PlotSelection select3D();
    /// returns selected rectangle (2d)
    PlotSelection select2D(tgt::ivec2 first, tgt::ivec2 second);
    /// renders the position of the mouse
    void renderMousePosition();
    /// clips mouse position to margins
    std::pair<tgt::ivec2, tgt::ivec2> clipToPlotBounds(tgt::ivec2 first, tgt::ivec2 second);
    /// sets regenerate display list flags true
    void regenDisplayLists();
    /// renders legends
    void renderLegends();
    /// create a simple zoom state in selectionProp_
    void createSimpleZoomState();
    /// update File
    void updateFile();
    void xLabelOptionChanged();
    void yLabelOptionChanged();
    void zLabelOptionChanged();

    PlotPort inport_;           ///< inport serve a PlotData or PlotFunction
    RenderPort outport_;        ///< outport where rendered plot goes
    RenderPort pickingBuffer_;  ///< buffer used by plotPickingManager_
    PlotData data_;             ///< data table
    PlotFunction function_;     ///< function which can be evaluated to create a PlotData
    PlotLibrary* plotLib_;       ///< our library offering the render functions
    PlotLibrary* exchangePlotLib_;


    PlotPickingManager plotPickingManager_; ///< picking manager for current PlotData

    // interaction handlers
    EventProperty<PlotProcessor>* eventHighlight_;          ///< calls highlight
    EventProperty<PlotProcessor>* eventLabel_;              ///< calls label
    EventProperty<PlotProcessor>* eventZoom_;               ///< calls zoom
    EventProperty<PlotProcessor>* eventHighlightAdditive_;  ///< calls unhighlight
    EventProperty<PlotProcessor>* eventLabelAdditive_;      ///< calls unlabel
    EventProperty<PlotProcessor>* eventZoomAdditive_;       ///< calls unzoom
    EventProperty<PlotProcessor>* mousePositionUpdateEvent_;
    EventProperty<PlotProcessor>* mouseEventEnterExit_;

    // general properties
    FileDialogProperty FileProp_;               ///< dialog for file selection
    ButtonProperty updateFile_;                 ///< update file
    BoolProperty plotLibraryViewPortClipping_;  ///< choose if clipping is used
    PlotDataProperty dataProp_;             ///< Offering a Data Table
    PlotEntitiesProperty plotEntitiesProp_; ///< Selection of Plot Entities
    StringProperty plotLabel_;              ///< Label of the Plot
    IntProperty plotLabelSize_;             ///< Size of Plot Label
    ColorProperty plotLabelColor_;      ///< Color of Plot Label
    ColorProperty selectionPlaneColor_; ///< Color of the selection planes
    ColorProperty highlightColor_;      ///< highlight color in plot
    BoolProperty enablePicking_;            ///< enable object picking
    BoolProperty renderMousePosition_;      ///< enable render mouse position
    BoolProperty renderLegend_;             ///< enable render legends

    // Margin
    IntProperty marginLeft_;
    IntProperty marginRight_;
    IntProperty marginBottom_;
    IntProperty marginTop_;

    IntProperty discreteStep_;            /// number of steps used for discretization of functions

    PlotSelectionProperty selectionProp_; ///< handles selections

    // 3D camera properties and handlers
    CameraProperty camera_;
    PlotCameraInteractionHandler* plotCameraHandler_;
    CameraInteractionHandler* cameraHandler_;
    BoolProperty orthographicCamera_;

    // axis properties
    BoolProperty renderAxes_;
    BoolProperty renderScales_;
    FloatProperty axesWidth_;
    IntProperty xScaleStep_;
    IntProperty yScaleStep_;
    IntProperty zScaleStep_;
    BoolProperty renderXHelperLines_;
    BoolProperty renderYHelperLines_;
    BoolProperty renderZHelperLines_;

    tgt::ivec2 mousePressedPosition_; ///< viewport position of last MouseDown event
    tgt::ivec2 mousePosition_;        ///< current viewport mouse position if on canvas
    bool mousePressed_;               ///< flag if a MouseDown event was registered but no MouseReleased event yet
    bool mouseInsideCanvas_;

    PlotLibrary::SelectionEdge selectionEdge_;//< encapsulates information about the region which we select
    int selectionEdgeIndex_; ///< only for 3D plots: specifies the axis of select clipping planes, -1 if there are currently no select planes.
    plot_t selectStart_;     ///< only for 3D plots: plot coordinates for the first selection plane
    plot_t selectEnd_;       ///< only for 3D plots: plot coordinates for the second selection plane

    bool inportHasPlotFunction_; ///< true iff the inport serves a PlotFunction which has the right keycolumncount (dimension-1)

    bool threeDimensional_;      ///< flag whether the plot is three dimensional

    bool regenDataList_;      ///< true if the display list of the data must be regenerated
    bool regenPickingList_;   ///< true if the display list of the picking data must be regenerated
private:
    StringOptionProperty xLabelOption_;
    StringProperty xLabel_;                 ///< Label of x axis
    StringOptionProperty yLabelOption_;
    StringProperty yLabel_;                 ///< Label of y axis
    StringOptionProperty zLabelOption_;
    StringProperty zLabel_;                 ///< Label of z axis
};

}   //namespace

#endif // VRN_PLOTPROCESSOR_H
