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

#include "plotlibrarylatexrender.h"

#include "tgt/vector.h"

namespace voreen {

const std::string PlotLibraryLatexRender::loggerCat_("voreen.plotting.PlotLibraryLatexRender");

// PlotLibraryLatexRender methods -------------------------------------------------------

PlotLibraryLatexRender::PlotLibraryLatexRender(const std::vector< PlotLibraryFileBase::Projection_Coordinates >& points, int clippingNumber)
: points_(points)
, clippingNumber_(clippingNumber)
{
}

PlotLibraryLatexRender::~PlotLibraryLatexRender() {
}

const std::vector< PlotLibraryFileBase::Projection_Coordinates >& PlotLibraryLatexRender::getPoints() {
    return points_;
}

std::vector< tgt::Vector3<plot_t> > PlotLibraryLatexRender::getPoints(PlotLibraryFileBase::Coordinate_Type type) {
    std::vector< tgt::Vector3<plot_t> > points;
    points.resize(points_.size());
    for (size_t i = 0; i < points_.size(); ++i) {
        points[i] = points_[i].get(type);
    }
    return points;
}


plot_t PlotLibraryLatexRender::minValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis) {
    plot_t minV = points_[0].get(type).elem[axis];
    for (size_t i = 1; i < points_.size(); ++i) {
        minV = std::min(points_[i].get(type).elem[axis],minV);
    }
    return minV;
}

plot_t PlotLibraryLatexRender::maxValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis) {
    plot_t maxV = points_[0].get(type).elem[axis];
    for (size_t i = 1; i < points_.size(); ++i) {
        maxV = std::max(points_[i].get(type).elem[axis],maxV);
    }
    return maxV;
}

plot_t PlotLibraryLatexRender::avgValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis) {
    plot_t avg = 0;
    for (size_t i = 0; i < points_.size(); ++i) {
        avg += points_[i].get(type).elem[axis];
    }
    return avg/points_.size();
}

int PlotLibraryLatexRender::getClippingNumber() const {
    return clippingNumber_;
}

std::string PlotLibraryLatexRender::latexColor(tgt::Color color) {
    std::stringstream latex;
    latex << "" << color.r << " " << color.g << " " << color.b;
    return latex.str();
}

std::string PlotLibraryLatexRender::latexTransparancy(tgt::Color color) {
    std::stringstream latex;
    latex << "opacity=" << color.a;
    return latex.str();
}

std::string PlotLibraryLatexRender::latexLineStyle(PlotEntitySettings::LineStyle lineStyle, double lineWidth) {
    std::stringstream latex;
    if (lineStyle!= PlotEntitySettings::CONTINUOUS) {
        if (lineStyle == PlotEntitySettings::DOTTED)
            latex << "linestyle=dotted,dotsep=1";
        else //DASHED
            latex << "linestyle=dashed,dash=10 10";
    }
    else if (lineWidth != 0)
        latex << "linestyle=solid";
    else
        latex << "linestyle=none";
    return latex.str();
}

std::string PlotLibraryLatexRender::latexFillStyle(PlotEntitySettings::PolygonFillStyle fillStyle) {
    std::stringstream latex;
    if (fillStyle == PlotEntitySettings::FILLNONE)
        latex << "fillstyle=none";
    else if (fillStyle == PlotEntitySettings::EVENODD)
        latex << "fillstyle=solid";
    else if (fillStyle == PlotEntitySettings::TEXTURE)
        latex << "fillstyle=solid";
    else
        latex << "fillstyle=solid";
    return latex.str();
}

// PlotLibraryLatexRenderLine methods -------------------------------------------------------

PlotLibraryLatexRenderLine::PlotLibraryLatexRenderLine(const PlotLibraryFileBase::Projection_Coordinates& startPoint, const PlotLibraryFileBase::Projection_Coordinates& endPoint, tgt::Color lineColor, int clipping_number, double lineWidth, PlotEntitySettings::LineStyle lineStyle)
: PlotLibraryLatexRender(std::vector< PlotLibraryFileBase::Projection_Coordinates >(),clipping_number)
, lineColor_(lineColor)
, lineWidth_(lineWidth)
, lineStyle_(lineStyle)
{
    points_.push_back(startPoint);
    points_.push_back(endPoint);

}

PlotLibraryLatexRenderLine::~PlotLibraryLatexRenderLine() {
}

std::string PlotLibraryLatexRenderLine::render() {
    std::stringstream latex;
    latex << "\\newrgbcolor{lcolor}{"<< PlotLibraryLatexRenderLine::latexColor(lineColor_) <<"}\n";
    latex << "\\psline[linewidth="<< lineWidth_ << ",linecolor=lcolor"
        << ",stroke" << PlotLibraryLatexRenderLine::latexTransparancy(lineColor_)
        << "," << PlotLibraryLatexRenderLine::latexLineStyle(lineStyle_,lineWidth_)
        << "]"
        << "(" << points_.at(0).afterProjection_.x << "," << points_.at(0).afterProjection_.y << ")(" << points_.at(1).afterProjection_.x << "," << points_.at(1).afterProjection_.y << ")"
        << "\n";
    return latex.str();
}

std::vector< PlotLibraryLatexRender* > PlotLibraryLatexRenderLine::generateLines(const std::vector< PlotLibraryFileBase::Projection_Coordinates >& endPoints) {
    std::vector< PlotLibraryLatexRender* > result;
    for (size_t i = 0; i < endPoints.size() - 1; ++i) {
        result.push_back(new PlotLibraryLatexRenderLine(endPoints[i],endPoints[i+1],lineColor_,getClippingNumber(),lineWidth_,lineStyle_));
    }
    return result;
}



// PlotLibraryLatexRenderCircle methods -------------------------------------------------------

PlotLibraryLatexRenderCircle::PlotLibraryLatexRenderCircle(const PlotLibraryFileBase::Projection_Coordinates& point, plot_t radius, tgt::Color fillColor, tgt::Color borderColor, int clipping_number, double lineWidth, int textureNumber, PlotEntitySettings::PolygonFillStyle fillstyle, PlotEntitySettings::LineStyle lineStyle)
: PlotLibraryLatexRender(std::vector< PlotLibraryFileBase::Projection_Coordinates >(),clipping_number)
, radius_(radius)
, fillColor_(fillColor)
, lineColor_(borderColor)
, lineWidth_(lineWidth)
, fillstyle_(fillstyle)
, lineStyle_(lineStyle)
//, textureNumber_(textureNumber)
{
    points_.push_back(point);
}

PlotLibraryLatexRenderCircle::~PlotLibraryLatexRenderCircle() {
}

std::string PlotLibraryLatexRenderCircle::render() {
    std::stringstream latex;
    latex << "\\newrgbcolor{lcolor}{"<< PlotLibraryLatexRenderCircle::latexColor(lineColor_) <<"}\n";
    latex << "\\newrgbcolor{fcolor}{"<< PlotLibraryLatexRenderCircle::latexColor(fillColor_) <<"}\n";
    latex << "\\pscircle[linewidth="<< lineWidth_ << ",linecolor=lcolor"
        << ",stroke" << PlotLibraryLatexRenderCircle::latexTransparancy(lineColor_)
        << "," << PlotLibraryLatexRenderCircle::latexLineStyle(lineStyle_,lineWidth_)
        << "," << PlotLibraryLatexRenderCircle::latexFillStyle(fillstyle_);
        if (fillstyle_ != PlotEntitySettings::FILLNONE) {
            latex << ",fillcolor=fcolor"
            << "," << PlotLibraryLatexRenderPolyline::latexTransparancy(fillColor_);
        }
    latex << "]"
        << "(" << points_.at(0).afterProjection_.x << "," << points_.at(0).afterProjection_.y << "){" << radius_ << "}"
        << "\n";
    return latex.str();
}

plot_t PlotLibraryLatexRenderCircle::minValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis) {
    if (axis == PlotLibrary::Z_AXIS) {
        return points_[0].get(type).z;
    }
    else {
        return points_[0].get(type).elem[axis] - radius_;
    }
}

plot_t PlotLibraryLatexRenderCircle::maxValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis) {
    if (axis == PlotLibrary::Z_AXIS) {
        return points_[0].get(type).z;
    }
    else {
        return points_[0].get(type).elem[axis] + radius_;
    }
}
plot_t PlotLibraryLatexRenderCircle::avgValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis) {
    return points_[0].get(type).elem[axis];
}



// PlotLibraryLatexRenderEllipse methods -------------------------------------------------------

PlotLibraryLatexRenderEllipse::PlotLibraryLatexRenderEllipse(const PlotLibraryFileBase::Projection_Coordinates& mpoint, double xradius, double yradius, tgt::Color fillColor, tgt::Color lineColor, int clipping_number, double lineWidth, int textureNumber, PlotEntitySettings::PolygonFillStyle fillstyle, PlotEntitySettings::LineStyle lineStyle)
: PlotLibraryLatexRender(std::vector< PlotLibraryFileBase::Projection_Coordinates >(),clipping_number)
, xradius_(xradius)
, yradius_(yradius)
, fillColor_(fillColor)
, lineColor_(lineColor)
, lineWidth_(lineWidth)
, fillstyle_(fillstyle)
, lineStyle_(lineStyle)
//, textureNumber_(textureNumber)
{
    points_.push_back(mpoint);
}

PlotLibraryLatexRenderEllipse::~PlotLibraryLatexRenderEllipse() {
}

std::string PlotLibraryLatexRenderEllipse::render() {
    std::stringstream latex;
    latex << "\\newrgbcolor{lcolor}{"<< PlotLibraryLatexRenderCircle::latexColor(lineColor_) <<"}\n";
    latex << "\\newrgbcolor{fcolor}{"<< PlotLibraryLatexRenderCircle::latexColor(fillColor_) <<"}\n";
    latex << "\\psellipse[linewidth="<< lineWidth_ << ",linecolor=lcolor"
        << ",stroke" << PlotLibraryLatexRenderCircle::latexTransparancy(lineColor_)
        << "," << PlotLibraryLatexRenderCircle::latexLineStyle(lineStyle_,lineWidth_)
        << "," << PlotLibraryLatexRenderCircle::latexFillStyle(fillstyle_);
        if (fillstyle_ != PlotEntitySettings::FILLNONE) {
            latex << ",fillcolor=fcolor"
            << "," << PlotLibraryLatexRenderPolyline::latexTransparancy(fillColor_);
        }
    latex << "]"
        << "(" << points_.at(0).afterProjection_.x << "," << points_.at(0).afterProjection_.y << ")(" << xradius_ << "," << yradius_ << ")"
        << "\n";
    return latex.str();
}

plot_t PlotLibraryLatexRenderEllipse::minValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis) {
    if (axis == PlotLibrary::Z_AXIS) {
        return points_[0].get(type).z;
    }
    else if (axis == PlotLibrary::X_AXIS){
        return points_[0].get(type).elem[axis] - xradius_;
    }
    else
        return points_[0].get(type).elem[axis] - yradius_;
}

plot_t PlotLibraryLatexRenderEllipse::maxValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis) {
    if (axis == PlotLibrary::Z_AXIS) {
        return points_[0].get(type).z;
    }
    else if (axis == PlotLibrary::X_AXIS) {
        return points_[0].get(type).elem[axis] + xradius_;
    }
    else
        return points_[0].get(type).elem[axis] + yradius_;
}
plot_t PlotLibraryLatexRenderEllipse::avgValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis) {
    return points_[0].get(type).elem[axis];
}



// PlotLibraryLatexRenderRect methods -------------------------------------------------------

PlotLibraryLatexRenderRect::PlotLibraryLatexRenderRect(const PlotLibraryFileBase::Projection_Coordinates& startPoint, double stepWidth, double height, tgt::Color startColor, tgt::Color endColor, double gradangle, double gradmidpoint)
: PlotLibraryLatexRender(std::vector< PlotLibraryFileBase::Projection_Coordinates >(),-1)
, stepWidth_(stepWidth)
, height_(height)
, startColor_(startColor)
, endColor_(endColor)
, gradangle_(gradangle)
, gradmidpoint_(gradmidpoint)
{
    points_.push_back(startPoint);
}


PlotLibraryLatexRenderRect::~PlotLibraryLatexRenderRect() {
}

std::string PlotLibraryLatexRenderRect::render() {
    std::stringstream latex;
    latex << "\\newrgbcolor{scolor}{"<< PlotLibraryLatexRenderCircle::latexColor(startColor_) <<"}\n";
    latex << "\\newrgbcolor{ecolor}{"<< PlotLibraryLatexRenderCircle::latexColor(endColor_) <<"}\n";
    latex << "\\psframe[linestyle=none,fillstyle=gradient,gradbegin=scolor"
        << ",gradend=ecolor,gradangle=" << gradangle_<< ",gradmidpoint=" << gradmidpoint_
        << "]"
        << "(" << points_.at(0).afterProjection_.x-1 << "," << points_.at(0).afterProjection_.y-1 << ")(" << points_.at(0).afterProjection_.x+stepWidth_+0.7 << "," << points_.at(0).afterProjection_.y+height_+1 << ")"
        << "\n";
    return latex.str();
}

plot_t PlotLibraryLatexRenderRect::minValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis) {
    return points_[0].get(type).elem[axis];;
}

plot_t PlotLibraryLatexRenderRect::maxValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis) {
    if (axis == PlotLibrary::Z_AXIS) {
        return points_[0].get(type).z;
    }
    else if (axis == PlotLibrary::X_AXIS) {
        return points_[0].get(type).elem[axis] + stepWidth_;
    }
    else
        return points_[0].get(type).elem[axis] + height_;
}

plot_t PlotLibraryLatexRenderRect::avgValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis) {
    return points_[0].get(type).elem[axis];
}


// PlotLibraryLatexRenderPolygon methods -------------------------------------------------------

PlotLibraryLatexRenderPolygon::PlotLibraryLatexRenderPolygon(const std::vector< PlotLibraryFileBase::Projection_Coordinates >& points, tgt::Color fillColor, tgt::Color lineColor, int clipping_number, double lineWidth, int textureNumber, PlotEntitySettings::PolygonFillStyle fillstyle, PlotEntitySettings::LineStyle lineStyle)
: PlotLibraryLatexRender(points,clipping_number)
, fillColor_(fillColor)
, lineColor_(lineColor)
, lineWidth_(lineWidth)
, fillstyle_(fillstyle)
, lineStyle_(lineStyle)
//, textureNumber_(textureNumber)
{
}

PlotLibraryLatexRenderPolygon::~PlotLibraryLatexRenderPolygon() {
}

std::string PlotLibraryLatexRenderPolygon::render() {
    std::stringstream latex;
    latex << "\\newrgbcolor{lcolor}{"<< PlotLibraryLatexRenderPolygon::latexColor(lineColor_) <<"}\n";
    latex << "\\newrgbcolor{fcolor}{"<< PlotLibraryLatexRenderPolygon::latexColor(fillColor_) <<"}\n";
    latex << "\\pspolygon[linewidth="<< lineWidth_ << ",linecolor=lcolor"
        << ",stroke" << PlotLibraryLatexRenderPolygon::latexTransparancy(lineColor_)
        << "," << PlotLibraryLatexRenderPolygon::latexLineStyle(lineStyle_,lineWidth_)
        << "," << PlotLibraryLatexRenderPolygon::latexFillStyle(fillstyle_);
        if (fillstyle_ != PlotEntitySettings::FILLNONE) {
            latex << ",fillcolor=fcolor"
            << "," << PlotLibraryLatexRenderPolyline::latexTransparancy(fillColor_);
        }
    latex << "]";
    for (size_t i = 0; i < points_.size(); ++i) {
        latex << "(" << points_.at(i).afterProjection_.x << "," << points_.at(i).afterProjection_.y << ")";
    }
    latex << "\n";
    return latex.str();
}



// PlotLibraryLatexRenderPolyline methods -------------------------------------------------------

PlotLibraryLatexRenderPolyline::PlotLibraryLatexRenderPolyline(const std::vector< PlotLibraryFileBase::Projection_Coordinates >& points, tgt::Color fillColor, tgt::Color lineColor, int clipping_number, double lineWidth, int textureNumber, PlotEntitySettings::PolygonFillStyle fillstyle, PlotEntitySettings::LineStyle lineStyle)
: PlotLibraryLatexRender(points,clipping_number)
, fillColor_(fillColor)
, lineColor_(lineColor)
, lineWidth_(lineWidth)
, fillstyle_(fillstyle)
, lineStyle_(lineStyle)
//, textureNumber_(textureNumber)
{
}

PlotLibraryLatexRenderPolyline::~PlotLibraryLatexRenderPolyline() {
}

std::string PlotLibraryLatexRenderPolyline::render() {
    std::stringstream latex;
    latex << "\\newrgbcolor{lcolor}{"<< PlotLibraryLatexRenderPolyline::latexColor(lineColor_) <<"}\n";
    latex << "\\newrgbcolor{fcolor}{"<< PlotLibraryLatexRenderPolyline::latexColor(fillColor_) <<"}\n";
    latex << "\\psPline[linewidth="<< lineWidth_ << ",linecolor=lcolor"
        << ",stroke" << PlotLibraryLatexRenderPolyline::latexTransparancy(lineColor_)
        << "," << PlotLibraryLatexRenderPolyline::latexLineStyle(lineStyle_,lineWidth_)
        << "," << PlotLibraryLatexRenderPolyline::latexFillStyle(fillstyle_);
        if (fillstyle_ != PlotEntitySettings::FILLNONE) {
            latex << ",fillcolor=fcolor"
            << "," << PlotLibraryLatexRenderPolyline::latexTransparancy(fillColor_);
        }
    latex << "]";
    for (size_t i = 0; i < points_.size(); ++i) {
        latex << "(" << points_.at(i).afterProjection_.x << "," << points_.at(i).afterProjection_.y << ")";
    }
    latex << "\n";

    return latex.str();
}


// PlotLibraryLatexRenderRect methods -------------------------------------------------------

PlotLibraryLatexRenderText::PlotLibraryLatexRenderText(const tgt::Vector3<plot_t>& startPoint, const std::string& text)
: PlotLibraryLatexRender(std::vector< PlotLibraryFileBase::Projection_Coordinates >(),-1)
, text_(text)
{
    points_.push_back(PlotLibraryFileBase::Projection_Coordinates(startPoint,startPoint,startPoint));
}


PlotLibraryLatexRenderText::~PlotLibraryLatexRenderText() {
}

std::string PlotLibraryLatexRenderText::render() {
    //std::stringstream ss;
    //ss << text_ <<  "\n";
    return text_;//ss.str();
}


plot_t PlotLibraryLatexRenderText::minValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis) {
    if (axis == PlotLibrary::Z_AXIS)
        return startPoint_.z-1.0;
    else
        return points_[0].get(type).elem[axis];
}

plot_t PlotLibraryLatexRenderText::maxValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis) {
    if (axis == PlotLibrary::Z_AXIS)
        return startPoint_.z-1.0;
    else
        return points_[0].get(type).elem[axis];
}

plot_t PlotLibraryLatexRenderText::avgValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis) {
    if (axis == PlotLibrary::Z_AXIS)
        return startPoint_.z-1.0;
    else
        return points_[0].get(type).elem[axis];
}


} // namespace
