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

#ifndef VRN_PLOTLIBRARYSVG_H
#define VRN_PLOTLIBRARYSVG_H

#include "plotlibraryfilebase.h"

#include <sstream>
#include <string>

namespace voreen {

class PlotPickingManager;
class PlotLibrarySvgRender;



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
class VRN_CORE_API PlotLibrarySvg : public PlotLibraryFileBase {
public:

    PlotLibrarySvg();

    virtual ~PlotLibrarySvg();

    //set values before starting rendering
    void beforeRender();

    //set values after rendering i.e. write in File
    void afterRender();

    ///sets the SVG status, i.e. color and projections matrix, should be called before rendering
    bool setRenderStatus();

    /// resets the SVG status, should be called after rendering
    void resetRenderStatus();


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
    void renderLine(const PlotData& data, int indexX, int indexY);

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
    void renderSpline(const PlotData& data, int indexX, int indexY);

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
    void renderErrorline(const PlotData& data, int indexX, int indexY, int indexError);

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
    void renderErrorspline(const PlotData& data, int indexX, int indexY, int indexError);

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
    void renderErrorbars(const PlotData& data, int indexX, int indexY, int indexError);

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
    void renderSurface(const PlotData& data, const std::vector<int>& triangleVertexIndices, bool wire, int indexX, int indexY,
            int indexZ, int indexCM = -1, bool wireonly = false);

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
    void renderHeightmap(const voreen::PlotData& data, const std::vector< std::list< tgt::dvec2 > >& voronoiRegions, bool wire,
            int indexZ, int indexCM = -1, bool wireonly = false);

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
    void renderCandlesticks(const PlotData& data, int indexX, int stickTop, int stickBottom, int candleTop, int candleBottom);

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
    void renderScatter(const PlotData& data, int indexX, int indexY, int indexZ = -1, int indexCM = -1, int indexSize = -1);

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
    void renderBars(const PlotData& data, std::vector<int> indicesY);


    /**
     * \brief   Renders a node graph of the nodes and connections in \a graph.
     */
    void renderNodeGraph(const PlotData& nodeData, const PlotData& connectionData, int indexX, int indexY, int indexDx, int indexDy);

    /**
     * \brief   Renders a legend for the current color map with given interval in the top right edge
     *
     * \param   data        PlotData table containing the values
     * \param   column      column which contains the values mapped to the current color map
     * \param   number      number of color map legend (optional, increment each time if you want to render more than one legend)
     **/
    void renderColorMapLegend(const PlotData& data, int column, int number = 0);

    /**
     * \brief   Renders the axes of the Plot using current axes width and color.
     **/
    void renderAxes();

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
    void renderAxisScales(Axis axis, bool helperLines, const std::string& label = "", plot_t offset = 0);

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
    void renderAxisLabelScales(const PlotData& data, int indexLabel, bool helperLines);

    ///**
    // * \brief   Labels the axis \a axis with the text in \a label using current font size and color.
    // * X axis labels will be positioned below right the end of the axis.
    // * Y axis labels will be positioned above left the end of the axis.
    // * Please note that this method is just for 2d plots.
    // *
    // * \param   axis    the axis to label
    // * \param   label   text for the label
    // **/
    //void renderAxisLabel(Axis axis, const std::string& label);

    /**
     * \brief   Renders the label \a text aligned around the plot coordinates \a pos using current font size and color.
     *
     * \param   pos             position of label
     * \param   align           alignment of the label around given position
     * \param   text            the text to be rendered
     * \param   viewCoordinates indicates if pos is already in view coordinates
     * \param   padding         additional space in pixels (viewport coordinates) around label
     **/
    void renderLabel(tgt::vec3 pos, const SmartLabel::Alignment align, const std::string& text, bool viewCoordinates = false,
            int padding = 10);

    /**
     * \brief   Renders the label \a text aligned around the viewport coordinates \a pos using current font size and color.
     *
     * \param   pos             position of label in viewport coordinates
     * \param   align           alignment of the label around given position
     * \param   text            the text to be rendered
     * \param   padding         additional space in pixels (viewport coordinates) around label
     **/
    void renderLabel(tgt::dvec2 pos, const SmartLabel::Alignment align, const std::string& text, int padding = 10);


    /// resets the line labels, should be called before calling line render functions
    void resetLineLabels();
    /// renders the line labels, should be called after calling line render functions
    void renderLineLabels();

    /// resets the plot labels, should be called before calling plot render functions
    void resetPlotLabels();
    /// renders the plot labels, should be called after calling plot render functions
    void renderPlotLabels();

    /// add a label to the plotLabelGroup_
    void addPlotLabel(std::string text, tgt::vec3 position, tgt::Color color, int size, SmartLabel::Alignment align);

    /// add a label to the lineLabelGroup_
    void addLineLabel(std::string text, tgt::vec3 position, tgt::Color color, int size, SmartLabel::Alignment align);

protected:


    struct PlotLibrarySvgSorter {
    public:
        bool operator() (PlotLibrarySvgRender* a, PlotLibrarySvgRender* b) const;
    } plotLibrarySvgSorter;



    /// renders a single bar, the size can be squeezed
    void renderSingleBar(plot_t left, plot_t right, plot_t bottom, plot_t top, tgt::Color c, plot_t squeeze = 1.0);

    /// renders a single glyph, the z value is only used for 3d plots
    void renderGlyph(plot_t x, plot_t y, plot_t z = 0, plot_t size = 1);
    /// renders a Circle in SVG
    void svgCircle(const PlotLibraryFileBase::Projection_Coordinates& point,plot_t radius, tgt::Color fillColor, tgt::Color borderColor, int clipping_number, double lineWidth = 0, PlotEntitySettings::PolygonFillStyle fillstyle = PlotEntitySettings::NONZERO, PlotEntitySettings::LineStyle lineStyle = PlotEntitySettings::CONTINUOUS);
    /// renders a Line in SVG
    void svgLine(const PlotLibraryFileBase::Projection_Coordinates& startPoint, const PlotLibraryFileBase::Projection_Coordinates& endPoint, tgt::Color lineColor, int clipping_number, double lineWidth, PlotEntitySettings::LineStyle lineStyle = PlotEntitySettings::CONTINUOUS);
    /// renders a Polygon in SVG
    void svgPolygon(const std::vector< PlotLibraryFileBase::Projection_Coordinates >& points, tgt::Color fillColor, tgt::Color lineColor, int clipping_number, double lineWidth = 0, PlotEntitySettings::PolygonFillStyle fillstyle = PlotEntitySettings::EVENODD, PlotEntitySettings::LineStyle lineStyle = PlotEntitySettings::CONTINUOUS);
    /// renders a Polyline in SVG
    void svgPolyline(const std::vector< PlotLibraryFileBase::Projection_Coordinates >& points, tgt::Color fillColor, tgt::Color lineColor, int clipping_number, double lineWidth = 1, PlotEntitySettings::PolygonFillStyle fillstyle = PlotEntitySettings::FILLNONE, PlotEntitySettings::LineStyle lineStyle = PlotEntitySettings::CONTINUOUS);
    /// renders a Ellipse in SVG
    void svgEllipse(const PlotLibraryFileBase::Projection_Coordinates& mpoint, double xradius, double yradius, tgt::Color fillColor, tgt::Color lineColor, int clipping_number, double lineWidth = 0, PlotEntitySettings::PolygonFillStyle fillstyle = PlotEntitySettings::NONZERO, PlotEntitySettings::LineStyle lineStyle = PlotEntitySettings::CONTINUOUS);

    /// renders a Gradient in SVG
    void svgGradient(const tgt::Vector2<plot_t>& startPoint, const tgt::Vector2<plot_t>& endPoint, const std::string& id, tgt::Color startColor, tgt::Color endColor);
    /// renders a Rect in SVG
    void svgRect(const PlotLibraryFileBase::Projection_Coordinates& startPoint, double stepWidth, double height, const std::string& fillurl);

    /// renders a Text in SVG
    void svgText(const tgt::Vector3<plot_t>& point, const std::string& text, tgt::Font* font, double size, tgt::Color textColor, const SmartLabel::Alignment align);
    /// renders a Textarea in SVG 1.1
    void svgTextAreaSVG11(const tgt::Vector3<plot_t>& point, double width, double height, const std::string& text,  tgt::Font* font, double size, tgt::Color textColor);
    /// renders a Textarea in SVG 1.2
    void svgTextAreaSVG12(const tgt::Vector3<plot_t>& point, double width, double height, const std::string& text,  tgt::Font* font, double size, tgt::Color textColor);
    /// renders a Sphere in SVG
    void svgSphere(const tgt::Vector3<plot_t>& mpoint, double radius, tgt::Color fillColor, tgt::Color lineColor, int sections, int clipping_number, double lineWidth = 0,PlotEntitySettings::PolygonFillStyle fillstyle = PlotEntitySettings::NONZERO, PlotEntitySettings::LineStyle lineStyle = PlotEntitySettings::CONTINUOUS);
    /// renders a Quad in SVG
    void svgQuad(const tgt::Vector3<plot_t>& mpoint, double xsize, double ysize, double zsize, tgt::Color fillColor, tgt::Color lineColor, int clipping_number, double lineWidth = 0,PlotEntitySettings::PolygonFillStyle fillstyle = PlotEntitySettings::NONZERO, PlotEntitySettings::LineStyle lineStyle = PlotEntitySettings::CONTINUOUS);
    /// renders a Triangle in SVG
    void svgTriangle3D(const tgt::Vector3<plot_t>& mpoint, double height, double weight, tgt::Color fillColor, tgt::Color lineColor, int clipping_number, double lineWidth = 0,PlotEntitySettings::PolygonFillStyle fillstyle = PlotEntitySettings::NONZERO, PlotEntitySettings::LineStyle lineStyle = PlotEntitySettings::CONTINUOUS);

private:

    /// renders the SmartLabelGroup \a smg
    void renderSmartLabelGroup(SmartLabelGroupBase* smg);

    mutable SmartLabelGroupSvgNoLayoutWithBackground plotLabelGroup_;    ///< smartlabelgroup of plot labels
    mutable SmartLabelGroupSvgVerticalMoving         lineLabelGroup_;    ///< smartlabelgroup of linelabels
    mutable SmartLabelGroupSvgHorizontalMoving       xAxisLabelGroup_;   ///< smartlabelgroup of x axis labels
    mutable SmartLabelGroupSvgNoLayout               axisLabelGroup_;    ///< smartlabelgroup of 3d axis labels

    std::stringstream svgOutput_;
    tgt::Color svgColor_;

    int textureNumber_;
    int fountNumberCounter_;

    std::vector< PlotLibrarySvgRender* > renderVector_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_PLOTSVG_H
