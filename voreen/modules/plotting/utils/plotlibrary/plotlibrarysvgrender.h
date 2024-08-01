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

#ifndef VRN_PLOTLIBRARYSVGRENDER_H
#define VRN_PLOTLIBRARYSVGRENDER_H

#include "plotlibrarysvg.h"

#include <sstream>
#include <string>

namespace voreen {

class VRN_CORE_API PlotLibrarySvgRender  {
public:


    PlotLibrarySvgRender(const std::vector< PlotLibraryFileBase::Projection_Coordinates >& points);

    virtual ~PlotLibrarySvgRender();

    virtual std::string render() = 0;

    const std::vector< PlotLibraryFileBase::Projection_Coordinates >& getPoints();
    std::vector< tgt::Vector3<plot_t> > getPoints(PlotLibraryFileBase::Coordinate_Type type);
    virtual plot_t minValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis);
    virtual plot_t maxValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis);
    virtual plot_t avgValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis);

    static std::string svgColor(tgt::Color color);
    static std::string svgTransparancy(tgt::Color color);
    static std::string svgClipping(int clipping_number);
    static std::string svgLineStyle(PlotEntitySettings::LineStyle lineStyle);

protected:

    std::vector< PlotLibraryFileBase::Projection_Coordinates > points_;

private:

    static const std::string loggerCat_;
};

class VRN_CORE_API PlotLibrarySvgRenderLine : public PlotLibrarySvgRender {
public:

    PlotLibrarySvgRenderLine(const PlotLibraryFileBase::Projection_Coordinates& startPoint, const PlotLibraryFileBase::Projection_Coordinates& endPoint, tgt::Color lineColor, int clipping_number, double lineWidth, PlotEntitySettings::LineStyle lineStyle = PlotEntitySettings::CONTINUOUS);

    virtual ~PlotLibrarySvgRenderLine();

    virtual std::string render();

    std::vector< PlotLibrarySvgRender* > generateLines(const std::vector< PlotLibraryFileBase::Projection_Coordinates >& endPoints);

    double getDeepValueAtAlpha(PlotLibraryFileBase::Coordinate_Type type, double alpha);
    double getDeepValueAtValue(PlotLibraryFileBase::Coordinate_Type type, tgt::Vector2<plot_t> value);

private:
    tgt::Color lineColor_;
    int clipping_number_;
    double lineWidth_;
    PlotEntitySettings::LineStyle lineStyle_;
};

class VRN_CORE_API PlotLibrarySvgRenderCircle : public PlotLibrarySvgRender {
public:

    PlotLibrarySvgRenderCircle(const PlotLibraryFileBase::Projection_Coordinates& point, plot_t radius, tgt::Color fillColor, tgt::Color borderColor, int clipping_number, double linewidth, int textureNumber = 0, PlotEntitySettings::PolygonFillStyle fillstyle = PlotEntitySettings::NONZERO, PlotEntitySettings::LineStyle lineStyle = PlotEntitySettings::CONTINUOUS);

    virtual ~PlotLibrarySvgRenderCircle();

    virtual std::string render();
    virtual plot_t minValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis);
    virtual plot_t maxValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis);
    virtual plot_t avgValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis);


private:
    plot_t radius_;
    tgt::Color fillColor_;
    tgt::Color lineColor_;
    int clipping_number_;
    double lineWidth_;
    PlotEntitySettings::PolygonFillStyle fillstyle_;
    PlotEntitySettings::LineStyle lineStyle_;
    int textureNumber_;
};

class VRN_CORE_API PlotLibrarySvgRenderEllipse : public PlotLibrarySvgRender {
public:

    PlotLibrarySvgRenderEllipse(const PlotLibraryFileBase::Projection_Coordinates& mpoint, double xradius, double yradius, tgt::Color fillColor, tgt::Color lineColor, int clipping_number, double lineWidth, int textureNumber = 0, PlotEntitySettings::PolygonFillStyle fillstyle = PlotEntitySettings::NONZERO, PlotEntitySettings::LineStyle lineStyle = PlotEntitySettings::CONTINUOUS);

    virtual ~PlotLibrarySvgRenderEllipse();

    virtual std::string render();
    virtual plot_t minValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis);
    virtual plot_t maxValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis);
    virtual plot_t avgValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis);

private:
    double xradius_;
    double yradius_;
    tgt::Color fillColor_;
    tgt::Color lineColor_;
    int clipping_number_;
    double lineWidth_;
    PlotEntitySettings::PolygonFillStyle fillstyle_;
    PlotEntitySettings::LineStyle lineStyle_;
    int textureNumber_;
};


class VRN_CORE_API PlotLibrarySvgRenderRect : public PlotLibrarySvgRender {
public:

    PlotLibrarySvgRenderRect(const PlotLibraryFileBase::Projection_Coordinates& startPoint, double stepWidth, double height, const std::string& fillurl);

    virtual ~PlotLibrarySvgRenderRect();

    virtual std::string render();
    virtual plot_t minValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis);
    virtual plot_t maxValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis);
    virtual plot_t avgValue(PlotLibraryFileBase::Coordinate_Type type, PlotLibrary::Axis axis);

private:
    double stepWidth_;
    double height_;
    std::string fillurl_;
};



class VRN_CORE_API PlotLibrarySvgRenderPolygon : public PlotLibrarySvgRender {
public:

    PlotLibrarySvgRenderPolygon(const std::vector< PlotLibraryFileBase::Projection_Coordinates >& points, tgt::Color fillColor, tgt::Color lineColor, int clipping_number, double lineWidth, int textureNumber = 0, PlotEntitySettings::PolygonFillStyle fillstyle = PlotEntitySettings::EVENODD, PlotEntitySettings::LineStyle lineStyle = PlotEntitySettings::CONTINUOUS);

    virtual ~PlotLibrarySvgRenderPolygon();

    virtual std::string render();

private:
    tgt::Color fillColor_;
    tgt::Color lineColor_;
    int clipping_number_;
    double lineWidth_;
    PlotEntitySettings::PolygonFillStyle fillstyle_;
    PlotEntitySettings::LineStyle lineStyle_;
    int textureNumber_;
};

class VRN_CORE_API PlotLibrarySvgRenderPolyline : public PlotLibrarySvgRender {
public:

    PlotLibrarySvgRenderPolyline(const std::vector< PlotLibraryFileBase::Projection_Coordinates >& points, tgt::Color fillColor, tgt::Color lineColor, int clipping_number, double lineWidth, int textureNumber = 0, PlotEntitySettings::PolygonFillStyle fillstyle = PlotEntitySettings::FILLNONE, PlotEntitySettings::LineStyle lineStyle = PlotEntitySettings::CONTINUOUS);

    virtual ~PlotLibrarySvgRenderPolyline();

    virtual std::string render();

private:

    tgt::Color fillColor_;
    tgt::Color lineColor_;
    int clipping_number_;
    double lineWidth_;
    PlotEntitySettings::PolygonFillStyle fillstyle_;
    PlotEntitySettings::LineStyle lineStyle_;
    int textureNumber_;

};

class VRN_CORE_API PlotLibrarySvgRenderText : public PlotLibrarySvgRender {
public:

    PlotLibrarySvgRenderText(const tgt::Vector3<plot_t>& startPoint, const std::string& text);

    virtual ~PlotLibrarySvgRenderText();

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
