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

#ifndef VRN_PLOTLIBRARYOPENGL_H
#define VRN_PLOTLIBRARYOPENGL_H

#include "plotlibrarynonefilebase.h"
#include "voreen/core/datastructures/geometry/glmeshgeometry.h"

#include "tgt/shadermanager.h"

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
class VRN_CORE_API PlotLibraryOpenGl : public PlotLibraryNoneFileBase {
public:

    PlotLibraryOpenGl();

    virtual ~PlotLibraryOpenGl();

    ///sets the opengl status, i.e. color and projections matrix, should be called before rendering
    bool setRenderStatus();

    /// resets the opengl status, should be called after rendering
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

    tgt::dvec2 convertPlotCoordinatesToViewport(const tgt::dvec3& plotCoordinates) const;
    tgt::dvec3 convertPlotCoordinatesToViewport3(const tgt::dvec3& plotCoordinates) const;

    /// converts viewport coordinates to plot coordinates using gluunproject
    tgt::dvec2 convertViewportToPlotCoordinates(tgt::ivec2 viewCoord) const;

    typedef GlMeshGeometryUInt32Color MeshType;

    /**
     * \brief   Creates a mesh useful to render a surface into the plot using the data in \a data.
     * \a indexX specifies the column containing x coordinates, \a indexY specifies the column
     * containing y coordinates, \a indexZ specifies the column containing z coordinates. This
     * method expects a valid delaunay triangulation given in \a triangleVertexIndices.
     *
     *
     * \param   data                    PlotData table containing the values
     * \param   triangleVertexIndices   Indices of the delaunay triangulation, each index corresponds to one row in \a data, every 3 indices build one triangle
     * \param   indexX                  column index of x coordinates
     * \param   indexY                  column index of y coordinates
     * \param   indexZ                  column index of z coordinates
     * \param   indexCM                 column index of color data (default = -1 means current drawing color shall be used instead of current colormap)
     **/
    MeshType* createSurfaceMesh(const PlotData& data, const std::vector<int>& triangleVertexIndices, int indexX, int indexY, int indexZ, int indexCM);

    /**
     * \brief   Creates a mesh useful to render the picking of a surface using the data in \a data.
     * \a indexX specifies the column containing x coordinates, \a indexY specifies the column
     * containing y coordinates, \a indexZ specifies the column containing z coordinates. This
     * method expects a valid delaunay triangulation given in \a triangleVertexIndices.
     *
     *
     * \param   data                    PlotData table containing the values
     * \param   triangleVertexIndices   Indices of the delaunay triangulation, each index corresponds to one row in \a data, every 3 indices build one triangle
     * \param   indexX                  column index of x coordinates
     * \param   indexY                  column index of y coordinates
     * \param   indexZ                  column index of z coordinates
     * \param   indexCM                 column index of color data (default = -1 means current drawing color shall be used instead of current colormap)
     **/
    MeshType* createSurfacePickingMesh(const PlotData& data, const std::vector<int>& triangleVertexIndices, int indexX, int indexY, int indexZ, int indexCM);

    /**
     * \brief   Creates a mesh useful to render a heightmap into the plot using the data in \a data.
     * \a indexZ specifies the column containing z coordinates. This
     * method expects a valid delaunay triangulation given in \a triangleVertexIndices.
     *
     *
     * \param   data                    PlotData table containing the values
     * \param   triangleVertexIndices   Indices of the delaunay triangulation, each index corresponds to one row in \a data, every 3 indices build one triangle
     * \param   indexZ                  column index of z coordinates
     * \param   indexCM                 column index of color data (default = -1 means current drawing color shall be used instead of current colormap)
     */
    MeshType* createHeightmapMesh(const voreen::PlotData& data, const std::vector< std::list< tgt::dvec2 > >& voronoiRegions, int indexZ, int indexCM = -1);

    /**
     * \brief   Creates a mesh useful to render the picking of a heightmap into the plot using the data in \a data.
     * \a indexZ specifies the column containing z coordinates. This
     * method expects a valid delaunay triangulation given in \a triangleVertexIndices.
     *
     *
     * \param   data                    PlotData table containing the values
     * \param   triangleVertexIndices   Indices of the delaunay triangulation, each index corresponds to one row in \a data, every 3 indices build one triangle
     * \param   indexZ                  column index of z coordinates
     * \param   indexCM                 column index of color data (default = -1 means current drawing color shall be used instead of current colormap)
     */
    MeshType* createHeightmapPickingMesh(const voreen::PlotData& data, const std::vector< std::list< tgt::dvec2 > >& voronoiRegions, int indexZ, int indexCM = -1);

protected:

    /// renders a single bar, the size can be squeezed
    void renderSingleBar(plot_t left, plot_t right, plot_t bottom, plot_t top, tgt::Color c, plot_t squeeze = 1.0);

    /// renders a single glyph, the z value is only used for 3d plots
    void renderGlyph(plot_t x, plot_t y, plot_t z = 0, plot_t size = 1);

    /// Used to render glyphs
    GlMeshGeometryUInt16TexCoord diskMesh_; /// No lighting => no normals
    GlMeshGeometryUInt16TexCoord triangleMesh_;
    GlMeshGeometryUInt16TexCoord rectMesh_;
    GlMeshGeometryUInt16NormalTexCoord sphereMesh_; /// Lighting => normals
    GlMeshGeometryUInt16NormalTexCoord tetrahedronMesh_;
    GlMeshGeometryUInt16NormalTexCoord quadMesh_;

    tgt::Shader& getSmoothPointShader();


private:

    inline void logGlVertex2d(plot_t x, plot_t y) const;
    inline void logGlVertex3d(plot_t x, plot_t y, plot_t z) const;


    /// renders the SmartLabelGroup \a smg
    void renderSmartLabelGroup(SmartLabelGroupBase* smg);


    mutable SmartLabelGroupOpenGlNoLayoutWithBackground plotLabelGroup_;    ///< smartlabelgroup of plot labels
    mutable SmartLabelGroupOpenGlVerticalMoving         lineLabelGroup_;    ///< smartlabelgroup of linelabels
    mutable SmartLabelGroupOpenGlHorizontalMoving       xAxisLabelGroup_;   ///< smartlabelgroup of x axis labels
    mutable SmartLabelGroupOpenGlNoLayout               axisLabelGroup_;    ///< smartlabelgroup of 3d axis labels

    tgt::Shader* smoothPointShader_;

    static const std::string loggerCat_;
};

} // namespace voreen

#endif // VRN_PLOTOPENGL_H
