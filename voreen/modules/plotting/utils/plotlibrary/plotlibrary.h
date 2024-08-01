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

#ifndef VRN_PLOTLIBRARY_H
#define VRN_PLOTLIBRARY_H

#include "../../datastructures/colormap.h"
#include "../../datastructures/plotdata.h"
#include "../smartlabel.h"
#include "../../datastructures/plotentitysettings.h"

#include "tgt/vector.h"
#include "tgt/font.h"
#include "tgt/camera.h"
#include "tgt/texture.h"

#include <list>
#include <valarray>

namespace voreen {

class PlotPickingManager;

/**
 *   This auxiliary library is a back end providing various rendering functions.
 *
 *   Like OpenGL this back end behaves as state machine, encapsulating drawing style
 *   and scales. The latter ones are calculated automatically by setting the dimensions
 *   and domain (i.e. min/max values) of the plot data. With the scales set, the rendering
 *   backend takes care of transforming the data values into plot coordinates. This library
 *   is normaly used as follows:
 *   1) call beforeRender
 *   2) set the status
 *   3) call setRenderStatus (sets in particular the transformation matrices)
 *   4a) call render methods
 *   4b) renew the status (e.g. drawing color)
 *   5) call resetRenderStatus (see https://dev.voreen.org/wiki/VoreenGl/Status)
 *   6) call afterRender
 *
 **/
class VRN_CORE_API PlotLibrary {
public:
    /// enum for bar grouping mode
    enum BarGroupingMode {
        GROUPED,
        STACKED,
        MERGED
    };

    /// enum for specifying axis
    enum Axis {
        X_AXIS = 0,
        Y_AXIS = 1,
        Z_AXIS = 2
    };

    /// enum for specifying dimension
    enum Dimension {
        TWO = 2,
        THREE = 3,
        FAKETHREE = 4
    };

    /**
     * \brief   Small struct encapsulating an outer projected axis edge and its orientation.
     *
     * \note    Outer projected means that this edge is part of the convex hull of all
     *          projected axis edges.
     **/
    struct SelectionEdge {
    public:
        SelectionEdge(Axis axis = X_AXIS, const tgt::dvec2& startVertex = tgt::dvec2(0,0), const tgt::dvec2& endVertex = tgt::dvec2(0,0),
                bool ascOrientation = true)
            : axis_(axis)
            , startVertex_(startVertex)
            , endVertex_(endVertex)
            , ascOrientation_(ascOrientation)
        {};

        Axis axis_;                 ///< axis of the edge

        tgt::dvec2 startVertex_;    ///< viewport coordinates of start vertex of the edge
        tgt::dvec2 endVertex_;      ///< viewport coordinates of end vertex of the edge

        bool ascOrientation_;       ///< orientation of this edge - true means ascending, false descending
    };

    PlotLibrary();

    virtual ~PlotLibrary();

    virtual void beforeRender() = 0;

    virtual void afterRender() = 0;


    ///sets the opengl status, i.e. color and projections matrix, should be called before rendering
    virtual bool setRenderStatus() = 0;

    /// resets the opengl status, should be called after rendering
    virtual void resetRenderStatus() = 0;

    /**
     * \brief   Renders a line into the plot using current line width and color.
     * \a indexX specifies the column containing x coordinates, \a indexY specifies the column
     * containing y coordinates.
     *
     * \note    This method is able to render pickable objects.
     *
     * \param   data        PlotData table containing the values
     * \param   indexX      column index of x coordinates
     * \param   indexY      column index of y coordinates
     **/
    virtual void renderLine(const PlotData& data, int indexX, int indexY) = 0;

    /**
     * \brief   Renders a spline into the plot using current line width and color.
     * \a indexX specifies the column containing x coordinates, \a indexY specifies the column
     * containing y coordinates.
     *
     * \note    This method is able to render pickable objects.
     *
     * \param   data        PlotData table containing the values
     * \param   indexX      column index of x coordinates
     * \param   indexY      column index of y coordinates
     **/
    virtual void renderSpline(const PlotData& data, int indexX, int indexY) = 0;

    /**
     * \brief   Renders a filled area behind the line visualizing the inaccurateness of the data using current fill color.
     * \a indexX specifies the column containing x coordinates, \a indexY specifies the column
     * containing y coordinates, \a indexError specifies the column containing the error data
     * which will be mapped to heigth of the error area.
     *
     * \note    This method is able to render pickable objects.
     *
     * \param   data        PlotData table containing the values
     * \param   indexX      column index of x coordinates
     * \param   indexY      column index of y coordinates
     * \param   indexError  column index of the error data
     **/
    virtual void renderErrorline(const PlotData& data, int indexX, int indexY, int indexError) = 0;

    /**
     * \brief   Renders a filled area behind the spline visualising the inaccurateness of the data using current fill color.
     * \a indexX specifies the column containing x coordinates, \a indexY specifies the column
     * containing y coordinates, \a indexError specifies the column containing the error data
     * which will be mapped to heigth of the error area.
     *
     * \note    This method is able to render pickable objects.
     *
     * \param   data        PlotData table containing the values
     * \param   indexX      column index of x coordinates
     * \param   indexY      column index of y coordinates
     * \param   indexError  column index of the error data
     **/
    virtual void renderErrorspline(const PlotData& data, int indexX, int indexY, int indexError) = 0;

    /**
     * \brief   Renders vertical error bars into the plot using the data in \a data.
     * \a indexX specifies the column containing x coordinates, \a indexY specifies the column
     * containing y coordinates, \a indexError specifies the column containing the error data
     * which will be mapped to the length of the error bars.
     *
     * \note    This method is able to render pickable objects.
     *
     * \param   data        PlotData table containing the values
     * \param   indexX      column index of x coordinates
     * \param   indexY      column index of y coordinates
     * \param   indexError  column index of the error data
     **/
    virtual void renderErrorbars(const PlotData& data, int indexX, int indexY, int indexError) = 0;

    /**
     * \brief   Renders a surface into the plot using the data in \a data.
     * \a indexX specifies the column containing x coordinates, \a indexY specifies the column
     * containing y coordinates, \a indexZ specifies the column containing z coordinates. This
     * method expects a valid delaunay triangulation given in \a triangleVertexIndices.
     *
     * \note    This method is able to render pickable objects.
     *          This method renders plot labels for highlighted cells.
     *
     * \param   data                    PlotData table containing the values
     * \param   triangleVertexIndices   Indices of the delaunay triangulation, each index corresponds to one row in \a data, every 3 indices build one triangle
     * \param   wire                    Flag whether to render only the wire frame
     * \param   indexX                  column index of x coordinates
     * \param   indexY                  column index of y coordinates
     * \param   indexZ                  column index of z coordinates
     * \param   indexCM                 column index of color data (default = -1 means current drawing color shall be used instead of current colormap)
     **/
    virtual void renderSurface(const PlotData& data, const std::vector<int>& triangleVertexIndices, bool wire, int indexX, int indexY,
            int indexZ, int indexCM = -1, bool wireonly = false) = 0;

    /**
     * \brief   Renders a heightmap into the plot using the data in \a data.
     * \a indexX specifies the column containing x coordinates, \a indexY specifies the column
     * containing y coordinates, \a indexZ specifies the column containing z coordinates. This
     * method expects a valid voronoi diagram given in \a voronoiRegions.
     *
     * \note    This method is able to render pickable objects.
     *          This method renders plot labels for highlighted cells.
     *
     * \param   data                PlotData table containing the values
     * \param   voronoiRegions      vector of voronoi regions each being a list of points (x,y-coordiantes) in counter-clockwise order
     * \param   wire                Flag whether to render only the wire frame
     * \param   indexX              column index of x coordinates
     * \param   indexY              column index of y coordinates
     * \param   indexZ              column index of z coordinates
     * \param   indexCM             column index of color data (-1 if current drawing color shall be used instead of current colormap )
     **/
    virtual void renderHeightmap(const voreen::PlotData& data, const std::vector< std::list< tgt::dvec2 > >& voronoiRegions, bool wire,
            int indexZ, int indexCM = -1, bool wireonly = false) = 0;

    /**
     * \brief   Renders a candle stick from data.
     *
     * \note    This method is able to render pickable objects.
     *
     * \param   data                PlotData table containing the values
     * \param   indexZ              column index of z coordinates
     * \param   stickTop            column index of stick top
     * \param   stickBottom         column index of stick bottom
     * \param   candleTop           column index of candle top
     * \param   candleBottom        column index of candle bottom
     **/
    virtual void renderCandlesticks(const PlotData& data, int indexX, int stickTop, int stickBottom, int candleTop, int candleBottom) = 0;

    /**
     * \brief   Renders a scatter plot from data.
     *
     * \note    This method is able to render pickable objects.
     *
     * \param   data                PlotData table containing the values
     * \param   indexX              column index of x coordinates
     * \param   indexY              column index of y coordinates
     * \param   indexZ              optional column index of z coordinates
     * \param   indexCM             column index of color data (-1 if the colormap isn't used)
     * \param   indexSize           column index of size data (-1 if only maxGlyphsize is used)
     **/
    virtual void renderScatter(const PlotData& data, int indexX, int indexY, int indexZ = -1, int indexCM = -1, int indexSize = -1) = 0;

    /**
     * \brief   Renders a bar plot of the data given in \a data.
     * The columns holding the values for the y axis values of the bars are
     * defined by \a indexesY. Bar colors are automatically chosen by current color map.
     *
     * \note    This method is able to render pickable objects.
     *          This method renders plot labels for highlighted cells.
     *
     * \param   data            PlotData table containing the values
     * \param   indexesY        column indexes of y coordinates of the bars (for each bar one index)
     **/
    virtual void renderBars(const PlotData& data, std::vector<int> indicesY) = 0;


    /**
     * \brief   Renders a node graph of the nodes and connections in \a graph.
     */
    virtual void renderNodeGraph(const PlotData& nodeData, const PlotData& connectionData, int indexX, int indexY, int indexDx, int indexDy) = 0;

    /**
     * \brief   Renders a legend for the current color map with given interval in the top right edge
     *
     * \param   data        PlotData table containing the values
     * \param   column      column which contains the values mapped to the current color map
     * \param   number      number of color map legend (optional, increment each time if you want to render more than one legend)
     **/
    virtual void renderColorMapLegend(const PlotData& data, int column, int number = 0)  = 0;

    /**
     * \brief   Renders the axes of the Plot using current axes width and color.
     **/
    virtual void renderAxes() = 0;

    /**
     * \brief   Renders the axis scales for axis \a axis and optionally helperlines.
     * The distance between each scale step will be automatically determined by according setting in
     * minimumScaleStep_ and be a power of 1, 2 or 5 (plot coordinates = values).
     * If \a helperLines is set, there will be rendered helper lines in the plot in current line size
     * and color for each scale.
     *
     * \param   axis            axis
     * \param   helperLines     flag whether to render helper lines
     **/
    virtual void renderAxisScales(Axis axis, bool helperLines, const std::string& label = "", plot_t offset = 0) = 0;

    /**
     * \brief   Renders the x axis scales using labels from \a data, and optionally helperlines.
     * The labels are read from the column \a indexLabel in \a data.
     * If \a helperLines is set, there will be rendered helper lines in the plot in current line size
     * and color for each scale.
     *
     * \param   data            data containing the labels
     * \param   indexLabel      column containing the labels
     * \param   helperLines     flag whether to render helper lines
     **/
    virtual void renderAxisLabelScales(const PlotData& data, int indexLabel, bool helperLines) = 0;

    /**
     * \brief   Labels the axis \a axis with the text in \a label using current font size and color.
     * X axis labels will be positioned below right the end of the axis.
     * Y axis labels will be positioned above left the end of the axis.
     * Please note that this method is just for 2d plots.
     *
     * \param   axis    the axis to label
     * \param   label   text for the label
     **/
    virtual void renderAxisLabel(Axis axis, const std::string& label);

    /**
     * \brief   Renders the label \a text aligned around the plot coordinates \a pos using current font size and color.
     *
     * \param   pos             position of label
     * \param   align           alignment of the label around given position
     * \param   text            the text to be rendered
     * \param   viewCoordinates indicates if pos is already in view coordinates
     * \param   padding         additional space in pixels (viewport coordinates) around label
     **/
    virtual void renderLabel(tgt::vec3 pos, const SmartLabel::Alignment align, const std::string& text, bool viewCoordinates = false,
            int padding = 10) = 0;

    /**
     * \brief   Renders the label \a text aligned around the viewport coordinates \a pos using current font size and color.
     *
     * \param   pos             position of label in viewport coordinates
     * \param   align           alignment of the label around given position
     * \param   text            the text to be rendered
     * \param   padding         additional space in pixels (viewport coordinates) around label
     **/
    virtual void renderLabel(tgt::dvec2 pos, const SmartLabel::Alignment align, const std::string& text, int padding = 10) = 0;


    /// returns the bounds of the plot area (surrounded by the axes)
    tgt::Bounds getBoundsPlot() const;
    /// returns the bounds of the area above the plot area
    tgt::Bounds getBoundsAbovePlot() const;
    /// returns the bounds of the area below the plot area
    tgt::Bounds getBoundsBelowPlot() const;
    /// returns the bounds of the area right of the plot area
    tgt::Bounds getBoundsRightOfPlot() const;
    /// returns the center above the plot area
    tgt::vec3 getCenterAbovePlot() const;

    /// resets the line labels, should be called before calling line render functions
    virtual void resetLineLabels() = 0;
    /// renders the line labels, should be called after calling line render functions
    virtual void renderLineLabels() = 0;

    /// resets the plot labels, should be called before calling plot render functions
    virtual void resetPlotLabels() = 0;
    /// renders the plot labels, should be called after calling plot render functions
    virtual void renderPlotLabels() = 0;

    /// renders a label showing the current mouse position
    virtual void renderMousePosition(tgt::ivec2 position);

    /// converts viewport coordinates to plot coordinates
    virtual tgt::dvec2 convertViewportToPlotCoordinates(tgt::ivec2 viewCoord) const = 0;

    void setDimension(Dimension dim);
    void setDomain(const Interval<plot_t>& interval, Axis axis);
//    void setPolarCoordinateFlag(bool polar);
    void setLogarithmicAxis(bool polar, Axis axis);
    void setCenterAxesFlag(bool value);

    void setDrawingColor(tgt::Color color);
    virtual void setWireColor(tgt::Color color) = 0;
    void setFillColor(tgt::Color color);
    void setFontColor(tgt::Color color);
    void setHighlightColor(tgt::Color color);
    void setColorMap(ColorMap cm);
    void setMinGlyphSize(float value);
    void setMaxGlyphSize(float value);
    void setFontSize(int size);
    void setLineWidth(float width);
    void setBarWidth(double width);
    void setAxesWidth(float width);
    void setLineStyle(PlotEntitySettings::LineStyle style);
    void setGlyphStyle(PlotEntitySettings::GlyphStyle style);
    void setMarginLeft(int margin);
    void setMarginRight(int margin);
    void setMarginTop(int margin);
    void setMarginBottom(int margin);
    void setWindowSize(const tgt::ivec2 windowSize);
    void setBarGroupingMode(BarGroupingMode mode);
    void setMinimumScaleStep(int value, Axis axis);
    void setCamera(const tgt::Camera& camera);
    void setShear(tgt::vec2 shear);
    void setSqueezeFactor(double sf);
    void setLightingFlag(bool value);
    void setOrthographicCameraFlag(bool value);
    void setTexture(tgt::Texture* texture, const std::string& texturepath);
    void setPlotPickingManager(PlotPickingManager* ppm);
    void setUsePlotPickingManager(bool value);
    void setViewPortClipping(bool value);

    /// Calculates the outer eges for each axis and saves them in according selectionEdges member.
    virtual void calculateSelectionEdges();

    const std::vector<SelectionEdge>& getSelectionEdgesX() const;
    const std::vector<SelectionEdge>& getSelectionEdgesY() const;
    const std::vector<SelectionEdge>& getSelectionEdgesZ() const;

    /// checks whether the polyline a-b-c is a left turn in the xy-plane
    static bool leftTurn(const tgt::dvec2& a, const tgt::dvec2& b, const tgt::dvec2& c);
    /// checks whether the polyline a-b-c is a right turn in the xy-plane
    static bool rightTurn(const tgt::dvec2& a, const tgt::dvec2& b, const tgt::dvec2& c);

    /// add a label to the plotLabelGroup_
    virtual void addPlotLabel(std::string text, tgt::vec3 position, tgt::Color color, int size, SmartLabel::Alignment align) = 0;

    /// add a label to the lineLabelGroup_
    virtual void addLineLabel(std::string text, tgt::vec3 position, tgt::Color color, int size, SmartLabel::Alignment align) = 0;

    /// converts the value from the PlotCells to plot coordinates
    virtual tgt::dvec2 convertPlotCoordinatesToViewport(const tgt::dvec3& plotCoordinates) const = 0;
    /// converts the value from the PlotCells to plot coordinates
    virtual tgt::dvec2 convertPlotCoordinatesToViewport(const plot_t x, const plot_t y, const plot_t z) const;
    /// converts the value from the PlotCells to plot coordinates
    virtual tgt::dvec3 convertPlotCoordinatesToViewport3(const tgt::dvec3& plotCoordinates) const = 0;
    /// converts the value from the PlotCells to plot coordinates
    virtual tgt::dvec3 convertPlotCoordinatesToViewport3(const plot_t x, const plot_t y, const plot_t z) const;

    inline virtual plot_t convertToLogCoordinates(plot_t value, Axis axis) const;

protected:


    /// compares pair<plot_t, color> but only the plot_t is used. We need it to order bars
    struct MergedBarSorter {
        bool operator() (std::pair<plot_t, tgt::Color> a, std::pair<plot_t, tgt::Color> b) const {
            if (a.first >= 0 || b.first >= 0)
                return (a.first < b.first);
            else
                return (b.first < a.first);
        };
    } mergedBarSorter;

    /// orientation predicate for the polyline a-b-c is a left turn in the xy-plane
    inline static double orientation(const tgt::dvec2& a, const tgt::dvec2& b, const tgt::dvec2& c);

    void renewPlotToViewportScale();

    inline plot_t convertFromLogCoordinates(plot_t value, Axis axis) const;
    inline plot_t round(plot_t value, plot_t roundingParameter) const;

    /// updates the scale steps for axis \a axis to powers of 1,2 or 5 and calculates
    /// by the way a rounding parameter
    tgt::dvec2 updateScaleSteps(Axis axis) const;

    /**
     * \brief helper function for \see PlotLibrary::calculateSelectionEdges()
     *
     * Finds all outer edges given by indices \a indices encoding vertices in \a vertices.
     * Outer edges will be written to \a out.
     **/
    void findEdges(const tgt::dvec2* vertices, const std::vector< std::valarray< int > >& indices, std::vector< SelectionEdge >& out, Axis axis) const;

    /// renders a single bar, the size can be squeezed
    virtual void renderSingleBar(plot_t left, plot_t right, plot_t bottom, plot_t top, tgt::Color c, plot_t squeeze = 1.0) = 0;

    /// renders a single glyph, the z value is only used for 3d plots
    virtual void renderGlyph(plot_t x, plot_t y, plot_t z = 0, plot_t size = 1) = 0;

    inline virtual tgt::Vector2<plot_t> logScale2dtoLogCoordinates(plot_t x, plot_t y) const;
    inline virtual tgt::Vector2<plot_t> logScale2dtoLogCoordinates(const tgt::Vector2<plot_t>& point) const;
    inline virtual tgt::Vector3<plot_t> logScale3dtoLogCoordinates(plot_t x, plot_t y, plot_t z) const;
    inline virtual tgt::Vector3<plot_t> logScale3dtoLogCoordinates(const tgt::Vector3<plot_t>& point) const;

    bool centerAxesFlag_;           ///< flag whether to center axes at origin or to render axes at bounds

    mutable tgt::Font labelFont_;
    tgt::ivec2 windowSize_;
    bool viewPortClipping_;

    Dimension dimension_;           ///< dimensionality of the plot (2, 3 or fake3d)
    Interval<plot_t> domain_[3];    ///< array of domain intervals for each axis
    bool logarithmicAxisFlags_[3];  ///< array of flags indicating if axis are logarithmic or not
    tgt::Color drawingColor_;       ///< current drawing color
    tgt::Color fillColor_;          ///< current filling color
    tgt::Color fontColor_;          ///< current font color
    tgt::Color highlightColor_;     ///< current drawing color for highlighted objects
    ColorMap colorMap_;             ///< current color map
    float minGlyphSize_;            ///< current min glyph size
    float maxGlyphSize_;            ///< current max glyph size
    float lineWidth_;               ///< current line width in pixel
    double  barWidth_;              ///< current bar width relative to maximal bar width
    float axesWidth_;               ///< current axes width in pixel
    int fontSize_;                  ///< current font size
    double squeezeFactor_;          ///< used for squeezing merged bars
    PlotEntitySettings::LineStyle lineStyle_;           ///< line style
    PlotEntitySettings::GlyphStyle glyphStyle_;         ///< glyph style
    BarGroupingMode barMode_;       ///< bar grouping mode
    int marginLeft_;                ///< left margin (in pixel)
    int marginRight_;               ///< right margin
    int marginBottom_;              ///< bottom margin
    int marginTop_;                 ///< top margin
    int minimumScaleStep_[3];       ///< minimum scale step for each axis in viewport coordinates
    tgt::Camera camera_;            ///< camera object for 3d plots
    tgt::vec2 shear_;               ///< shear used for fake 3d in bar plots
    bool lightingFlag_;             ///< flag whether lighting is enabled
    bool orthographicCameraFlag_;   ///< flag whether the camera uses orthographic projection
    tgt::Texture* texture_;         ///< texture used for glyphs
    std::string texturePath_;
    PlotPickingManager* ppm_;       ///< used to render pickable objects
    bool usePlotPickingManager_;    ///< if true, the color is set by PlotPickingManager if the render method is able to render pickable objects

    tgt::dvec2 plotToViewportScale_;

    // things for 3D zooming:
    std::vector<SelectionEdge> selectionEdgesX_;
    std::vector<SelectionEdge> selectionEdgesY_;
    std::vector<SelectionEdge> selectionEdgesZ_;

    ///// renders the SmartLabelGroup \a smg
    virtual void renderSmartLabelGroup(SmartLabelGroupBase* smg) = 0;

private:


    static const std::string loggerCat_;
};

//------------------------------ Implementation of the Inline-Functions --------------------------------------

plot_t PlotLibrary::convertFromLogCoordinates(plot_t value, Axis axis) const {
    //inverse of convertToLogCoordinates
    plot_t strechedLogValue = value - domain_[axis].getLeft();
    plot_t logValue = strechedLogValue * log(1 + domain_[axis].size()) / domain_[axis].size();
    return exp(logValue) + domain_[axis].getLeft() - 1;
}

plot_t PlotLibrary::convertToLogCoordinates(plot_t value, Axis axis) const {
    //normalize the interval to [1,1+intervalLength], apply log and fit it to the original domain
    plot_t logValue = log(value-domain_[axis].getLeft()+1);
    plot_t strechedLogValue = logValue * domain_[axis].size() / log(1+domain_[axis].size());
    return strechedLogValue + domain_[axis].getLeft();
}

double PlotLibrary::orientation(const tgt::dvec2& a, const tgt::dvec2& b, const tgt::dvec2& c) {
    return a.x*b.y + a.y*c.x +  b.x*c.y - c.x*b.y - c.y*a.x - b.x*a.y;
}

plot_t PlotLibrary::round(plot_t value, plot_t roundingParameter) const {
    return floor( value / roundingParameter + 0.5 ) * roundingParameter;
}

tgt::Vector2<plot_t> PlotLibrary::logScale2dtoLogCoordinates(plot_t x, plot_t y) const {
    tgt::Vector2<plot_t> copypoint = tgt::Vector2<plot_t>(x,y);
    if (logarithmicAxisFlags_[X_AXIS])
        copypoint.x = convertToLogCoordinates(x, X_AXIS);
    if (logarithmicAxisFlags_[Y_AXIS])
        copypoint.y = convertToLogCoordinates(y, Y_AXIS);
    return copypoint;
}

tgt::Vector2<plot_t> PlotLibrary::logScale2dtoLogCoordinates(const tgt::Vector2<plot_t>& point) const {
    return logScale2dtoLogCoordinates(point.x,point.y);
}

tgt::Vector3<plot_t> PlotLibrary::logScale3dtoLogCoordinates(plot_t x, plot_t y, plot_t z) const {
    tgt::Vector3<plot_t> copypoint = tgt::Vector3<plot_t>(x,y,z);
    if (logarithmicAxisFlags_[X_AXIS])
        copypoint.x = convertToLogCoordinates(x, X_AXIS);
    if (logarithmicAxisFlags_[Y_AXIS])
        copypoint.y = convertToLogCoordinates(y, Y_AXIS);
    if (logarithmicAxisFlags_[Z_AXIS])
        copypoint.z = convertToLogCoordinates(z, Z_AXIS);
    return copypoint;
}

tgt::Vector3<plot_t> PlotLibrary::logScale3dtoLogCoordinates(const tgt::Vector3<plot_t>& point) const {
    return logScale3dtoLogCoordinates(point.x,point.y,point.z);
}




} // namespace voreen

#endif // VRN_PLOTLIBRARY_H
