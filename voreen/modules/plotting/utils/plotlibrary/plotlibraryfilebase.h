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

#ifndef VRN_PLOTLIBRARYFILEBASE_H
#define VRN_PLOTLIBRARYFILEBASE_H

#include "plotlibrary.h"


namespace voreen {
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
class VRN_CORE_API PlotLibraryFileBase : public PlotLibrary {
public:

    enum PlotLibraryFileType {
        SVG = 0,
        LATEX = 1
    };

    enum Coordinate_Type {
        ORIGIN = 1,
        MODELVIEW = 2,
        PROJECTION = 3
    };


    struct Projection_Coordinates {
    public:
        Projection_Coordinates(const tgt::Vector3<plot_t>& origin, const tgt::Vector3<plot_t>& afterModelview,
            const tgt::Vector3<plot_t>& afterProjection)
            : origin_(origin)
            , afterModelview_(afterModelview)
            , afterProjection_(afterProjection)
        {};

        Projection_Coordinates()
            : origin_(tgt::Vector3<plot_t>(0,0,0))
            , afterModelview_(tgt::Vector3<plot_t>(0,0,0))
            , afterProjection_(tgt::Vector3<plot_t>(0,0,0))
        {};

        tgt::Vector3<plot_t> get(Coordinate_Type type);

        tgt::Vector3<plot_t> origin_;
        tgt::Vector3<plot_t> afterModelview_;
        tgt::Vector3<plot_t> afterProjection_;
    };


    PlotLibraryFileBase();

    virtual ~PlotLibraryFileBase() {};
    /// to set the Wirecolor ist drawing ist devided in drawing area and lines
    virtual void setWireColor(tgt::Color color);
    /// the the file in which should be write
    void setOutputFile(const std::string& file);

    /// converts the value from the PlotCells to plot coordinates
    virtual tgt::dvec2 convertPlotCoordinatesToViewport(const tgt::dvec3& plotCoordinates) const;
    /// converts the value from the PlotCells to plot coordinates
    virtual tgt::dvec2 convertPlotCoordinatesToViewport(const plot_t x, const plot_t y, const plot_t z) const;
    /// converts the value from the PlotCells to plot coordinates
    virtual tgt::dvec3 convertPlotCoordinatesToViewport3(const tgt::dvec3& plotCoordinates) const;
    /// converts the value from the PlotCells to plot coordinates
    virtual tgt::dvec3 convertPlotCoordinatesToViewport3(const plot_t x, const plot_t y, const plot_t z) const;
    /// converts viewport coordinates to plot coordinates
    virtual tgt::dvec2 convertViewportToPlotCoordinates(tgt::ivec2 viewCoord) const;

    /// converts the value from the PlotCells to plot coordinates and all other coodinates during the calculation
    virtual Projection_Coordinates convertPlotCoordinatesToViewport3Projection(const tgt::dvec3& plotCoordinates) const;
    /// converts the value from the PlotCells to plot coordinates and all other coodinates during the calculation
    virtual Projection_Coordinates convertPlotCoordinatesToViewport3Projection(const plot_t x, const plot_t y, const plot_t z) const;


protected:

    struct Node {
    public:
        Node(double x,double y, double z, Node* next = 0,Node* prev = 0,Node* nextPoly = 0, Node* neighbor = 0,
            bool intersect = false, bool entry = false, bool visited = false, double alpha = 0)
            : x_(x)
            , y_(y)
            , z_(z)
            , next_(next)
            , prev_(prev)
            , nextPoly_(nextPoly)
            , neighbor_(neighbor)
            , intersect_(intersect)
            , entry_(entry)
            , visited_(visited)
            , alpha_(alpha)
        {};

        virtual ~Node() {};

        tgt::Vector3<plot_t> getPoint() {return tgt::Vector3<plot_t>(x_,y_,z_);};

        double x_,y_,z_;
        Node* next_;
        Node* prev_;
        Node* nextPoly_;
        Node* neighbor_;
        bool intersect_;
        bool entry_;
        bool visited_;
        double alpha_;
    };

    struct PlotLibraryPointSorter {
        bool operator() (PlotLibraryFileBase::Projection_Coordinates a, PlotLibraryFileBase::Projection_Coordinates b) const {
            return ((a.origin_.x > b.origin_.x) ||
                   (a.origin_.x == b.origin_.x && a.origin_.y > b.origin_.y) ||
                   (a.origin_.x == b.origin_.x && a.origin_.y == b.origin_.y && a.origin_.z > b.origin_.z));
        };
    } PlotLibraryPointSorter;


    int intersect_Triangle_Ray(const std::vector< tgt::Vector3<plot_t> >& ray, const std::vector< tgt::Vector3<plot_t> >& triangle, tgt::Vector3<plot_t>* result_point);

    static int intersect_Ray_Ray(const std::vector< tgt::Vector2<plot_t> >& ray1, const std::vector< tgt::Vector2<plot_t> >& ray2, tgt::Vector2<plot_t>* result_values);

    inline tgt::Matrix4<plot_t> plOrtho(double leftVal, double rightVal, double bottomVal, double topVal, double nearVal, double farVal);

    void clippolygons(const std::vector< tgt::Vector3<plot_t> >& subject, const std::vector< tgt::Vector3<plot_t> >& clip, std::vector< std::vector< tgt::Vector3<plot_t> > >* resultPolygons);

    bool inClippRegion(const tgt::Vector3<plot_t>& point);
    bool inClippRegion(const std::vector< PlotLibraryFileBase::Projection_Coordinates >& points);

    int tri_tri_intersect3D(const tgt::Vector3<plot_t>& c1, const tgt::Vector3<plot_t>& c2, const tgt::Vector3<plot_t>& c3,
        const tgt::Vector3<plot_t>& d1, const tgt::Vector3<plot_t>& d2, const tgt::Vector3<plot_t>& d3);

    double windingnumber(const std::vector< tgt::Vector3<plot_t> >& polygon, const tgt::Vector3<plot_t>& point);

    Projection_Coordinates projection_neu(const tgt::Vector3<plot_t>& point);
    Projection_Coordinates projection_neu(plot_t x, plot_t y, plot_t z);


    tgt::Color lineColor_;
    std::string outputFile_;
    tgt::Matrix4<plot_t> projectionsMatrix_;
    tgt::Matrix4<plot_t> modelviewMatrix_;
    tgt::Matrix4<plot_t> invertedmodelviewMatrix_;

    double left_;
    double right_;
    double bottom_;
    double top_;
    double near_;
    double far_;

    std::vector< tgt::Vector4<plot_t> > clippingPlanes_;


private:

    Node* createVertex(const tgt::Vector3d& startPoint, const tgt::Vector3d& endPoint, double alpha);
    bool intersect(const tgt::Vector2<plot_t>& startPoint1, const tgt::Vector2<plot_t>& endPoint1,
        const tgt::Vector2<plot_t>& startPoint2, const tgt::Vector2<plot_t>& endPoint2, double* alpha1, double* alpha2);

    bool clippOnPlane(const tgt::Vector4<plot_t>& clippplane, const tgt::Vector4<plot_t>& eyepoint);

    tgt::Vector2<plot_t> perpendicularVector(const tgt::Vector2<plot_t>& vector);

    Node* searchNode(Node* startNode, const tgt::Vector3<plot_t>& point1,const tgt::Vector3<plot_t>& point2, double alpha);

    static double det(const tgt::Vector2<plot_t>& column1, const tgt::Vector2<plot_t>& column2);

    static const std::string loggerCat_;
};



//------------------------------ Implementation of the Inline-Functions --------------------------------------

inline tgt::Matrix4<plot_t> PlotLibraryFileBase::plOrtho(double leftVal, double rightVal, double bottomVal, double topVal, double nearVal, double farVal) {
    return tgt::Matrix4<plot_t>(
        2/(rightVal-leftVal),        0,                 0,            -(rightVal+leftVal)/(rightVal-leftVal),
                0,         2/(topVal-bottomVal),        0,            -(topVal+bottomVal)/(topVal-bottomVal),
                0,               0,           -2/(farVal-nearVal),    -(farVal+nearVal)/(farVal-nearVal),
                0,               0,                 0,                        1
        );
}

} // namespace voreen

#endif // VRN_PLOTLIBRARYFILEBASE_H
