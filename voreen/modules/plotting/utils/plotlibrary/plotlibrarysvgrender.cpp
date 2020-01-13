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

#include "plotlibrarysvgrender.h"

#include "tgt/vector.h"

namespace voreen {

const std::string PlotLibrarySvgRender::loggerCat_("voreen.plotting.PlotLibrarySvgRender");

// PlotLibrarySvgRender methods -------------------------------------------------------

PlotLibrarySvgRender::PlotLibrarySvgRender(const std::vector< PlotLibraryFileBase::Projection_Coordinates >& points)
    : points_(points)
{
}

PlotLibrarySvgRender::~PlotLibrarySvgRender() {
}


const std::vector< PlotLibraryFileBase::Projection_Coordinates >& PlotLibrarySvgRender::getPoints() {
    return points_;
}

std::vector< tgt::Vector3<plot_t> > PlotLibrarySvgRender::getPoints(PlotLibraryFileBase::Coordinate_Type type) {
    std::vector< tgt::Vector3<plot_t> > points;
    points.resize(points_.size());
    for (size_t i = 0; i < points_.size(); ++i) {
        points[i] = points_[i].get(type);
    }
    return points;
}


plot_t PlotLibrarySvgRender::minValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis) {
    plot_t minV = points_[0].get(type).elem[axis];
    for (size_t i = 1; i < points_.size(); ++i) {
        minV = std::min(points_[i].get(type).elem[axis],minV);
    }
    return minV;
}

plot_t PlotLibrarySvgRender::maxValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis) {
    plot_t maxV = points_[0].get(type).elem[axis];
    for (size_t i = 1; i < points_.size(); ++i) {
        maxV = std::max(points_[i].get(type).elem[axis],maxV);
    }
    return maxV;
}

plot_t PlotLibrarySvgRender::avgValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis) {
    plot_t avg = 0;
    for (size_t i = 0; i < points_.size(); ++i) {
        avg += points_[i].get(type).elem[axis];
    }
    return avg/points_.size();
}


std::string PlotLibrarySvgRender::svgColor(tgt::Color color) {
    std::stringstream svg;
    svg << "\"rgb("
        << static_cast<int>(tgt::round(color.r*255)) << ","
        << static_cast<int>(tgt::round(color.g*255)) << ","
        << static_cast<int>(tgt::round(color.b*255))
        << ")\"";
    return svg.str();
}

std::string PlotLibrarySvgRender::svgTransparancy(tgt::Color color) {
    std::stringstream svg;
    svg << "\"" << color.a << "\"";
    return svg.str();
}

std::string PlotLibrarySvgRender::svgClipping(int clipping_number) {
    std::stringstream svg;
    svg << " clip-path=\"url(#cp" << clipping_number << ")\"";
    return svg.str();
}

std::string PlotLibrarySvgRender::svgLineStyle(PlotEntitySettings::LineStyle lineStyle) {
    std::stringstream svg;
    if (lineStyle!= PlotEntitySettings::CONTINUOUS) {
        if (lineStyle == PlotEntitySettings::DOTTED)
            svg << "stroke-dasharray=\"1%, 1%\"";
        else //DASHED
            svg << "stroke-dasharray=\"2%, 2%\"";
    }
    return svg.str();
}

// PlotLibrarySvgRenderLine methods -------------------------------------------------------

PlotLibrarySvgRenderLine::PlotLibrarySvgRenderLine(const PlotLibraryFileBase::Projection_Coordinates& startPoint, const PlotLibraryFileBase::Projection_Coordinates& endPoint, tgt::Color lineColor, int clipping_number, double lineWidth, PlotEntitySettings::LineStyle lineStyle)
: PlotLibrarySvgRender(std::vector< PlotLibraryFileBase::Projection_Coordinates >())
, lineColor_(lineColor)
, clipping_number_(clipping_number)
, lineWidth_(lineWidth)
, lineStyle_(lineStyle)
{
    points_.push_back(startPoint);
    points_.push_back(endPoint);
}

PlotLibrarySvgRenderLine::~PlotLibrarySvgRenderLine() {
}

std::string PlotLibrarySvgRenderLine::render() {
    std::stringstream svg;
    svg << "<line x1=\"" << points_[0].afterProjection_.x << "\" y1=\"" << points_[0].afterProjection_.y << "\" x2=\"" << points_[1].afterProjection_.x << "\" y2=\"" << points_[1].afterProjection_.y
        << "\" ";
    svg << PlotLibrarySvgRender::svgLineStyle(lineStyle_);
    if (clipping_number_ >= 0)
        svg << PlotLibrarySvgRender::svgClipping(clipping_number_);
    svg << " stroke=" << PlotLibrarySvgRender::svgColor(lineColor_);
    svg << " stroke-opacity=" << PlotLibrarySvgRender::svgTransparancy(lineColor_);
    svg << " stroke-width=\"" << lineWidth_ << "\" />";
    //svg << "//" << points_[0].afterProjection_.z << "," << points_[1].afterProjection_.z;
    //svg << "   " << avgValue(PlotLibraryFileBase::PROJECTION,PlotLibrary::Z_AXIS);
    svg << "\n";
    return svg.str();
}

std::vector< PlotLibrarySvgRender* > PlotLibrarySvgRenderLine::generateLines(const std::vector< PlotLibraryFileBase::Projection_Coordinates >& endPoints) {
    std::vector< PlotLibrarySvgRender* > result;
    for (size_t i = 0; i < endPoints.size() - 1; ++i) {
        result.push_back(new PlotLibrarySvgRenderLine(endPoints[i],endPoints[i+1],lineColor_,clipping_number_,lineWidth_,lineStyle_));
    }
    return result;
}

double PlotLibrarySvgRenderLine::getDeepValueAtAlpha(PlotLibraryFileBase::Coordinate_Type type, double alpha) {
    tgt::Vector3<plot_t> point = points_[0].get(type) + alpha* (points_[1].get(type)-points_[0].get(type));
    return point.z;
}

double PlotLibrarySvgRenderLine::getDeepValueAtValue(PlotLibraryFileBase::Coordinate_Type type, tgt::Vector2<plot_t> value) {
    tgt::Vector3<plot_t> point0 = points_[0].get(type);
    tgt::Vector3<plot_t> point1 = points_[1].get(type);
    tgt::Vector3<plot_t> rVector = point1 - point0;
    double alpha = 0;
    if (rVector.x != 0)
        alpha = (value.x - point0.x) / rVector.x;
    else if (rVector.y != 0)
        alpha = (value.y - point0.y) / rVector.y;
    else {
        return
            maxValue(type,PlotLibraryFileBase::Z_AXIS);
    }
    return
        getDeepValueAtAlpha(type,alpha);
}

// PlotLibrarySvgRenderCircle methods -------------------------------------------------------

PlotLibrarySvgRenderCircle::PlotLibrarySvgRenderCircle(const PlotLibraryFileBase::Projection_Coordinates& point, plot_t radius, tgt::Color fillColor, tgt::Color borderColor, int clipping_number, double lineWidth, int textureNumber, PlotEntitySettings::PolygonFillStyle fillstyle, PlotEntitySettings::LineStyle lineStyle)
: PlotLibrarySvgRender(std::vector< PlotLibraryFileBase::Projection_Coordinates >())
, radius_(radius)
, fillColor_(fillColor)
, lineColor_(borderColor)
, clipping_number_(clipping_number)
, lineWidth_(lineWidth)
, fillstyle_(fillstyle)
, lineStyle_(lineStyle)
, textureNumber_(textureNumber)
{
    points_.push_back(point);
}

PlotLibrarySvgRenderCircle::~PlotLibrarySvgRenderCircle() {
}

std::string PlotLibrarySvgRenderCircle::render() {
    std::stringstream svg;
    svg << "<circle cx=\""<< points_[0].afterProjection_.x << "\" cy=\""<<  points_[0].afterProjection_.y << "\" r=\"" << radius_ << "\" ";
    svg << PlotLibrarySvgRender::svgLineStyle(lineStyle_);
    if (clipping_number_ >= 0)
        svg << PlotLibrarySvgRender::svgClipping(clipping_number_);
    if (fillstyle_ == PlotEntitySettings::FILLNONE) {
        svg << " fill=\"none\"";
    }
    else if (fillstyle_ == PlotEntitySettings::TEXTURE) {
        svg << " fill=\"url(#texture"<< textureNumber_ << ")\"";
    }
    else {
        svg << " fill=" << PlotLibrarySvgRender::svgColor(fillColor_);
        svg << " fill-opacity=" << PlotLibrarySvgRender::svgTransparancy(fillColor_);
    }
    svg << " stroke=" << PlotLibrarySvgRender::svgColor(lineColor_);
    svg << " stroke-opacity=" << PlotLibrarySvgRender::svgTransparancy(lineColor_);
    svg << " stroke-width=\"" << lineWidth_ << "\" />" << "\n";
    return svg.str();
}


plot_t PlotLibrarySvgRenderCircle::minValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis) {
    if (axis == PlotLibrary::Z_AXIS) {
        return points_[0].get(type).z;
    }
    else {
        return points_[0].get(type).elem[axis] - radius_;
    }
}

plot_t PlotLibrarySvgRenderCircle::maxValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis) {
    if (axis == PlotLibrary::Z_AXIS) {
        return points_[0].get(type).z;
    }
    else {
        return points_[0].get(type).elem[axis] + radius_;
    }
}
plot_t PlotLibrarySvgRenderCircle::avgValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis) {
    return points_[0].get(type).elem[axis];
}




// PlotLibrarySvgRenderEllipse methods -------------------------------------------------------

PlotLibrarySvgRenderEllipse::PlotLibrarySvgRenderEllipse(const PlotLibraryFileBase::Projection_Coordinates& mpoint, double xradius, double yradius, tgt::Color fillColor, tgt::Color lineColor, int clipping_number, double lineWidth, int textureNumber, PlotEntitySettings::PolygonFillStyle fillstyle, PlotEntitySettings::LineStyle lineStyle)
: PlotLibrarySvgRender(std::vector< PlotLibraryFileBase::Projection_Coordinates >())
, xradius_(xradius)
, yradius_(yradius)
, fillColor_(fillColor)
, lineColor_(lineColor)
, clipping_number_(clipping_number)
, lineWidth_(lineWidth)
, fillstyle_(fillstyle)
, lineStyle_(lineStyle)
, textureNumber_(textureNumber)
{
    points_.push_back(mpoint);
}

PlotLibrarySvgRenderEllipse::~PlotLibrarySvgRenderEllipse() {
}

std::string PlotLibrarySvgRenderEllipse::render() {
    std::stringstream svg;
    svg << "<ellipse cx=\"" << points_[0].afterProjection_.x << "\" cy=\"" << points_[0].afterProjection_.y
        << "\" rx=\"" << xradius_ << "\" ry=\"" << yradius_ << "\" ";
    svg << PlotLibrarySvgRender::svgLineStyle(lineStyle_);
    if (clipping_number_ >= 0)
        svg << PlotLibrarySvgRender::svgClipping(clipping_number_);
    if (fillstyle_ == PlotEntitySettings::FILLNONE) {
        svg << " fill=\"none\"";
    }
    else if (fillstyle_ == PlotEntitySettings::TEXTURE) {
        svg << " fill=\"url(#texture"<< textureNumber_ << ")\"";
    }
    else {
        svg << " fill=" << PlotLibrarySvgRender::svgColor(fillColor_);
        svg << " fill-opacity=" << PlotLibrarySvgRender::svgTransparancy(fillColor_);
    }
    svg << " stroke=" << PlotLibrarySvgRender::svgColor(lineColor_);
    svg << " stroke-opacity=" << PlotLibrarySvgRender::svgTransparancy(lineColor_);
    svg << " stroke-width=\"" << lineWidth_ << "\" />" << "\n";
    return svg.str();
}

plot_t PlotLibrarySvgRenderEllipse::minValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis) {
    if (axis == PlotLibrary::Z_AXIS) {
        return points_[0].get(type).z;
    }
    else if (axis == PlotLibrary::X_AXIS){
        return points_[0].get(type).elem[axis] - xradius_;
    }
    else
        return points_[0].get(type).elem[axis] - yradius_;
}

plot_t PlotLibrarySvgRenderEllipse::maxValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis) {
    if (axis == PlotLibrary::Z_AXIS) {
        return points_[0].get(type).z;
    }
    else if (axis == PlotLibrary::X_AXIS) {
        return points_[0].get(type).elem[axis] + xradius_;
    }
    else
        return points_[0].get(type).elem[axis] + yradius_;
}
plot_t PlotLibrarySvgRenderEllipse::avgValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis) {
    return points_[0].get(type).elem[axis];
}


// PlotLibrarySvgRenderRect methods -------------------------------------------------------

PlotLibrarySvgRenderRect::PlotLibrarySvgRenderRect(const PlotLibraryFileBase::Projection_Coordinates& startPoint, double stepWidth, double height, const std::string& fillurl)
: PlotLibrarySvgRender(std::vector< PlotLibraryFileBase::Projection_Coordinates >())
, stepWidth_(stepWidth)
, height_(height)
, fillurl_(fillurl)
{
    points_.push_back(startPoint);
}


PlotLibrarySvgRenderRect::~PlotLibrarySvgRenderRect() {
}

std::string PlotLibrarySvgRenderRect::render() {
    std::stringstream svg;
    svg << "<rect fill=\"url(#" << fillurl_ << ")\" width=\"" << stepWidth_ << "\" height=\""
        << height_ << "\" x=\"" << points_[0].afterProjection_.x << "\" y=\"" << points_[0].afterProjection_.y << "\" />" << "\n";
    return svg.str();
}

plot_t PlotLibrarySvgRenderRect::minValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis) {
    return points_[0].get(type).elem[axis];;
}

plot_t PlotLibrarySvgRenderRect::maxValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis) {
    if (axis == PlotLibrary::Z_AXIS) {
        return points_[0].get(type).z;
    }
    else if (axis == PlotLibrary::X_AXIS) {
        return points_[0].get(type).elem[axis] + stepWidth_;
    }
    else
        return points_[0].get(type).elem[axis] + height_;
}

plot_t PlotLibrarySvgRenderRect::avgValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis) {
    return points_[0].get(type).elem[axis];
}



// PlotLibrarySvgRenderPolygon methods -------------------------------------------------------

PlotLibrarySvgRenderPolygon::PlotLibrarySvgRenderPolygon(const std::vector< PlotLibraryFileBase::Projection_Coordinates >& points, tgt::Color fillColor, tgt::Color lineColor, int clipping_number, double lineWidth, int textureNumber, PlotEntitySettings::PolygonFillStyle fillstyle, PlotEntitySettings::LineStyle lineStyle)
: PlotLibrarySvgRender(points)
, fillColor_(fillColor)
, lineColor_(lineColor)
, clipping_number_(clipping_number)
, lineWidth_(lineWidth)
, fillstyle_(fillstyle)
, lineStyle_(lineStyle)
, textureNumber_(textureNumber)
{
}

PlotLibrarySvgRenderPolygon::~PlotLibrarySvgRenderPolygon() {
}

std::string PlotLibrarySvgRenderPolygon::render() {
    std::stringstream svg;
    svg << "<polygon points=\"";
    for (size_t i = 0; i < points_.size(); ++i) {
        svg << points_.at(i).afterProjection_.x << "," << points_.at(i).afterProjection_.y << " ";
    }
    svg << "\" ";
    svg << PlotLibrarySvgRender::svgLineStyle(lineStyle_);
    if (clipping_number_ >= 0)
        svg << PlotLibrarySvgRender::svgClipping(clipping_number_);
    if (fillstyle_ == PlotEntitySettings::FILLNONE) {
        svg << " fill=\"none\"";
    }
    else if (fillstyle_ == PlotEntitySettings::TEXTURE) {
        svg << " fill=\"url(#texture"<< textureNumber_ << ")\"";
    }
    else {
        svg << " fill-rule=\"";
        if (fillstyle_ == PlotEntitySettings::EVENODD)
            svg << "evenodd";
        else
            svg << "nonzero";
        svg << "\" fill=" << PlotLibrarySvgRender::svgColor(fillColor_);
        svg << " fill-opacity=" << PlotLibrarySvgRender::svgTransparancy(fillColor_);
    }
    svg << " stroke=" << PlotLibrarySvgRender::svgColor(lineColor_);
    svg << " stroke-opacity=" << PlotLibrarySvgRender::svgTransparancy(lineColor_);
    svg << " stroke-width=\"" << lineWidth_ << "\" />";
    //svg << "//";
    //for (size_t i = 0; i < points_.size(); ++i) {
    //    svg << points_[i].afterProjection_.z << ",";
    //}
    //svg << "   " << avgValue(PlotLibraryFileBase::PROJECTION,PlotLibrary::Z_AXIS);
    svg << "\n";
    return svg.str();
}



// PlotLibrarySvgRenderPolyline methods -------------------------------------------------------

PlotLibrarySvgRenderPolyline::PlotLibrarySvgRenderPolyline(const std::vector< PlotLibraryFileBase::Projection_Coordinates >& points, tgt::Color fillColor, tgt::Color lineColor, int clipping_number, double lineWidth, int textureNumber, PlotEntitySettings::PolygonFillStyle fillstyle, PlotEntitySettings::LineStyle lineStyle)
: PlotLibrarySvgRender(points)
, fillColor_(fillColor)
, lineColor_(lineColor)
, clipping_number_(clipping_number)
, lineWidth_(lineWidth)
, fillstyle_(fillstyle)
, lineStyle_(lineStyle)
, textureNumber_(textureNumber)
{
}

PlotLibrarySvgRenderPolyline::~PlotLibrarySvgRenderPolyline() {
}

std::string PlotLibrarySvgRenderPolyline::render() {
    std::stringstream svg;
    svg << "<polyline points=\"";
    for (size_t i = 0; i < points_.size(); ++i) {
        svg << points_.at(i).afterProjection_.x << "," << points_.at(i).afterProjection_.y << " ";
    }
    svg << "\" ";
    svg << PlotLibrarySvgRenderPolyline::svgLineStyle(lineStyle_);
    if (clipping_number_ >= 0)
        svg << PlotLibrarySvgRenderPolyline::svgClipping(clipping_number_);
    if (fillstyle_ == PlotEntitySettings::FILLNONE) {
        svg << " fill=\"none\"";
    }
    else if (fillstyle_ == PlotEntitySettings::TEXTURE) {
        svg << " fill=\"url(#texture"<< textureNumber_ << ")\"";
    }
    else {
        svg << " fill-rule=\"";
        if (fillstyle_ == PlotEntitySettings::EVENODD)
            svg << "evenodd";
        else
            svg << "nonzero";
        svg << "\" fill=" << PlotLibrarySvgRenderPolyline::svgColor(fillColor_);
        svg << " fill-opacity=" << PlotLibrarySvgRenderPolyline::svgTransparancy(fillColor_);
    }
    svg << " stroke=" << PlotLibrarySvgRenderPolyline::svgColor(lineColor_);
    svg << " stroke-opacity=" << PlotLibrarySvgRenderPolyline::svgTransparancy(lineColor_);
    svg << " stroke-width=\"" << lineWidth_ << "\" />";
    //svg << "//";
    //for (size_t i = 0; i < points_.size(); ++i) {
    //    svg << points_[i].afterProjection_.z << ",";
    //}
    //svg << "   " << avgValue(PlotLibraryFileBase::PROJECTION,PlotLibrary::Z_AXIS);
    svg << "\n";
    return svg.str();
}


// PlotLibrarySvgRenderText methods -------------------------------------------------------

PlotLibrarySvgRenderText::PlotLibrarySvgRenderText(const tgt::Vector3<plot_t>& startPoint, const std::string& text)
: PlotLibrarySvgRender(std::vector< PlotLibraryFileBase::Projection_Coordinates >())
, startPoint_(startPoint)
, text_(text)
{
    points_.push_back(PlotLibraryFileBase::Projection_Coordinates(startPoint,startPoint,startPoint));
}


PlotLibrarySvgRenderText::~PlotLibrarySvgRenderText() {
}

std::string PlotLibrarySvgRenderText::render() {
    //std::stringstream ss;
    //ss << text_ <<  "\n";
    return text_;//ss.str();
}


plot_t PlotLibrarySvgRenderText::minValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis) {
    if (axis == PlotLibrary::Z_AXIS)
        return startPoint_.z-1.0;
    else
        return points_[0].get(type).elem[axis];
}

plot_t PlotLibrarySvgRenderText::maxValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis) {
    if (axis == PlotLibrary::Z_AXIS)
        return startPoint_.z-1.0;
    else
        return points_[0].get(type).elem[axis];
}

plot_t PlotLibrarySvgRenderText::avgValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis) {
    if (axis == PlotLibrary::Z_AXIS)
        return startPoint_.z-1.0;
    else
        return points_[0].get(type).elem[axis];
}

} // namespace
