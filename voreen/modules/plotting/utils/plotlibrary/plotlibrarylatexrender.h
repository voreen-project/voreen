/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2021 University of Muenster, Germany,                        *
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

#ifndef VRN_PLOTLIBRARYLATEXRENDER_H
#define VRN_PLOTLIBRARYLATEXRENDER_H

#include "plotlibrarylatex.h"

#include <sstream>
#include <string>

namespace voreen {

class VRN_CORE_API PlotLibraryLatexRender  {
public:

    PlotLibraryLatexRender(const std::vector< PlotLibraryFileBase::Projection_Coordinates >& points, int clippingNumber);
    virtual ~PlotLibraryLatexRender();

    virtual std::string render() = 0;

    const std::vector< PlotLibraryFileBase::Projection_Coordinates >& getPoints();
    std::vector< tgt::Vector3<plot_t> > getPoints(PlotLibraryFileBase::Coordinate_Type type);
    virtual plot_t minValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis);
    virtual plot_t maxValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis);
    virtual plot_t avgValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis);
    int getClippingNumber() const;

    static std::string latexColor(tgt::Color color);
    static std::string latexTransparancy(tgt::Color color);
    static std::string latexLineStyle(PlotEntitySettings::LineStyle lineStyle, double lineWidth);
    static std::string latexFillStyle(PlotEntitySettings::PolygonFillStyle fillStyle);

protected:

    std::vector< PlotLibraryFileBase::Projection_Coordinates > points_;


private:
    int clippingNumber_;
    static const std::string loggerCat_;

};

class PlotLibraryLatexRenderLine : public PlotLibraryLatexRender {
public:

    PlotLibraryLatexRenderLine(const PlotLibraryFileBase::Projection_Coordinates& startPoint, const PlotLibraryFileBase::Projection_Coordinates& endPoint, tgt::Color lineColor, int clipping_number, double lineWidth, PlotEntitySettings::LineStyle lineStyle = PlotEntitySettings::CONTINUOUS);

    virtual ~PlotLibraryLatexRenderLine();

    virtual std::string render();

    std::vector< PlotLibraryLatexRender* > generateLines(const std::vector< PlotLibraryFileBase::Projection_Coordinates >& endPoints);

private:
    tgt::Color lineColor_;
    double lineWidth_;
    PlotEntitySettings::LineStyle lineStyle_;
};

class PlotLibraryLatexRenderCircle : public PlotLibraryLatexRender {
public:

    PlotLibraryLatexRenderCircle(const PlotLibraryFileBase::Projection_Coordinates& point, plot_t radius, tgt::Color fillColor, tgt::Color borderColor, int clipping_number, double linewidth, int textureNumber = 0, PlotEntitySettings::PolygonFillStyle fillstyle = PlotEntitySettings::NONZERO, PlotEntitySettings::LineStyle lineStyle = PlotEntitySettings::CONTINUOUS);

    virtual ~PlotLibraryLatexRenderCircle();

    virtual std::string render();
    virtual plot_t minValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis);
    virtual plot_t maxValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis);
    virtual plot_t avgValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis);

private:
    plot_t radius_;
    tgt::Color fillColor_;
    tgt::Color lineColor_;
    double lineWidth_;
    PlotEntitySettings::PolygonFillStyle fillstyle_;
    PlotEntitySettings::LineStyle lineStyle_;
    //int textureNumber_; //Currently not used
};

class PlotLibraryLatexRenderEllipse : public PlotLibraryLatexRender {
public:

    PlotLibraryLatexRenderEllipse(const PlotLibraryFileBase::Projection_Coordinates& mpoint, double xradius, double yradius, tgt::Color fillColor, tgt::Color lineColor, int clipping_number, double lineWidth, int textureNumber = 0, PlotEntitySettings::PolygonFillStyle fillstyle = PlotEntitySettings::NONZERO, PlotEntitySettings::LineStyle lineStyle = PlotEntitySettings::CONTINUOUS);

    virtual ~PlotLibraryLatexRenderEllipse();

    virtual std::string render();
    virtual plot_t minValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis);
    virtual plot_t maxValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis);
    virtual plot_t avgValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis);

private:

    double xradius_;
    double yradius_;
    tgt::Color fillColor_;
    tgt::Color lineColor_;
    double lineWidth_;
    PlotEntitySettings::PolygonFillStyle fillstyle_;
    PlotEntitySettings::LineStyle lineStyle_;
    //int textureNumber_; //Currently not used
};


class PlotLibraryLatexRenderRect : public PlotLibraryLatexRender {
public:

    PlotLibraryLatexRenderRect(const PlotLibraryFileBase::Projection_Coordinates& startPoint, double stepWidth, double height, tgt::Color startColor, tgt::Color endColor, double gradangle = 90, double gradmidpoint = 1);

    virtual ~PlotLibraryLatexRenderRect();

    virtual std::string render();
    virtual plot_t minValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis);
    virtual plot_t maxValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis);
    virtual plot_t avgValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis);

private:

    double stepWidth_;
    double height_;
    tgt::Color startColor_;
    tgt::Color endColor_;
    double gradangle_;
    double gradmidpoint_;
};



class PlotLibraryLatexRenderPolygon : public PlotLibraryLatexRender {
public:

    PlotLibraryLatexRenderPolygon(const std::vector< PlotLibraryFileBase::Projection_Coordinates >& points, tgt::Color fillColor, tgt::Color lineColor, int clipping_number, double lineWidth, int textureNumber = 0, PlotEntitySettings::PolygonFillStyle fillstyle = PlotEntitySettings::EVENODD, PlotEntitySettings::LineStyle lineStyle = PlotEntitySettings::CONTINUOUS);

    virtual ~PlotLibraryLatexRenderPolygon();

    virtual std::string render();

private:

    tgt::Color fillColor_;
    tgt::Color lineColor_;
    double lineWidth_;
    PlotEntitySettings::PolygonFillStyle fillstyle_;
    PlotEntitySettings::LineStyle lineStyle_;
    //int textureNumber_; //Currently not used
};

class PlotLibraryLatexRenderPolyline : public PlotLibraryLatexRender {
public:

    PlotLibraryLatexRenderPolyline(const std::vector< PlotLibraryFileBase::Projection_Coordinates >& points, tgt::Color fillColor, tgt::Color lineColor, int clipping_number, double lineWidth, int textureNumber = 0, PlotEntitySettings::PolygonFillStyle fillstyle = PlotEntitySettings::FILLNONE, PlotEntitySettings::LineStyle lineStyle = PlotEntitySettings::CONTINUOUS);

    virtual ~PlotLibraryLatexRenderPolyline();

    virtual std::string render();

private:

    tgt::Color fillColor_;
    tgt::Color lineColor_;
    double lineWidth_;
    PlotEntitySettings::PolygonFillStyle fillstyle_;
    PlotEntitySettings::LineStyle lineStyle_;
    //int textureNumber_; //Currently not used

};

class PlotLibraryLatexRenderText : public PlotLibraryLatexRender {
public:

    PlotLibraryLatexRenderText(const tgt::Vector3<plot_t>& startPoint, const std::string& text);

    virtual ~PlotLibraryLatexRenderText();

    virtual std::string render();
    virtual plot_t minValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis);
    virtual plot_t maxValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis);
    virtual plot_t avgValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis);

private:
    tgt::Vector3<plot_t> startPoint_;
    std::string text_;
};

} // namespace voreen

#endif // VRN_PLOTSVGRENDER_H
