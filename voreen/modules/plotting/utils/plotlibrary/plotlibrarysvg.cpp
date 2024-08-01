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

#include "plotlibrarysvg.h"

//#include "functionlibrary.h"
#include "../../datastructures/plotrow.h"
//#include "../interaction/plotpickingmanager.h"
#include "plotlibrarysvgrender.h"

#include "voreen/core/voreenapplication.h"
#include "voreen/core/version.h"

//#include "tgt/glmath.h"
#include "tgt/tgt_math.h"
#include "tgt/spline.h"
//#include "tgt/quadric.h"
#include "tgt/matrix.h"
#include "tgt/vector.h"

#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

#include <iomanip>

namespace voreen {

//namespace {
///// compares pair<plot_t, color> but only the plot_t is used. We need it to order bars
//struct MergedBarSorter {
//    bool operator() (std::pair<plot_t, tgt::Color> a, std::pair<plot_t, tgt::Color> b) const {
//        if (a.first >= 0 || b.first >= 0)
//            return (a.first < b.first);
//        else
//            return (b.first < a.first);
//    };
//} mergedBarSorter;
//};

const std::string PlotLibrarySvg::loggerCat_("voreen.plotting.PlotLibrarySvg");

PlotLibrarySvg::PlotLibrarySvg()
    : PlotLibraryFileBase()
    , plotLabelGroup_(&labelFont_, 5, tgt::Bounds(), tgt::Color(1, 1, 1, 0.75))
    , lineLabelGroup_(&labelFont_, 6, tgt::Bounds())
    , xAxisLabelGroup_(&labelFont_, 10, tgt::Bounds())
    , axisLabelGroup_(&labelFont_, 6)
    , svgColor_(0.f,0.f,0.f,1.f)
    , textureNumber_(0)
    , fountNumberCounter_(0)
    , renderVector_()
    //, dimension_(TWO)
    //, drawingColor_(0.0f, 0.0f, 0.0f, 1.0f)
    //, fillColor_(0.0f, 0.0f, 0.0f, 0.5f)
    //, fontColor_(0.0f, 0.0f, 0.0f, 1.0f)
    //, highlightColor_(1.0f, 0.95f, 0.9f, 0.5f)
    //, colorMap_(ColorMap::createColdHot())
    //, minGlyphSize_(1.0f)
    //, maxGlyphSize_(1.0f)
    //, lineWidth_(1.0f)
    //, barWidth_(0.5)
    //, axesWidth_(1.0f)
    //, fontSize_(10)
    //, squeezeFactor_(1.0)
    //, lineStyle_(PlotEntitySettings::CONTINUOUS)
    //, glyphStyle_(PlotEntitySettings::POINT)
    //, barMode_(GROUPED)
    //, marginLeft_(60)
    //, marginRight_(100)
    //, marginBottom_(50)
    //, marginTop_(50)
    //, shear_(tgt::vec2(0.f,0.f))
    //, lightingFlag_(false)
    //, orthographicCameraFlag_(true)
    //, texture_(0)
    //, ppm_(NULL)
    //, usePlotPickingManager_(false)
{
    labelFont_.setFontSize(8);
    domain_[0] = Interval<plot_t>(0, 1, false, false);
    domain_[1] = Interval<plot_t>(0, 1, false, false);
    domain_[2] = Interval<plot_t>(0, 1, false, false);
    logarithmicAxisFlags_[0] = false;
    logarithmicAxisFlags_[1] = false;
    logarithmicAxisFlags_[2] = false;
    minimumScaleStep_[0] = 40;
    minimumScaleStep_[1] = 40;
    minimumScaleStep_[2] = 40;
}

PlotLibrarySvg::~PlotLibrarySvg() {
    for (size_t i = 0; i < renderVector_.size(); ++i)
        delete renderVector_[i];
    renderVector_.clear();
}

void PlotLibrarySvg::beforeRender() {
    svgOutput_.str("");
    svgOutput_.clear();
    fountNumberCounter_ = 0;
    svgColor_ = drawingColor_;
    svgOutput_ << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
        << "<!-- Created with Voreen " << VoreenVersion::getVersion() << " (http://www.voreen.org/) -->\n\n"
        << "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n"
        << "<svg xmlns=\"http://www.w3.org/2000/svg\"\n"
        << "  xmlns:xlink=\"http://www.w3.org/1999/xlink\" xmlns:ev=\"http://www.w3.org/2001/xml-events\"\n"
        << "  xmlns:voreen=\"http://www.voreen.org\"\n"
        << "  version=\"1.1\" baseProfile=\"full\"\n"
        << "  width=\"" << windowSize_.x << "px\" height=\"" << windowSize_.y << "px\" "
        << "  viewBox=\"0 0 " << windowSize_.x << " " << windowSize_.y << "\" \n"
        << "  voreen:version=\"" << VoreenVersion::getVersion() << "\"\n"
        << "  voreen:copyright=\"" << VoreenVersion::getCopyright() << "\""
        << " >\n"
        << "<title id=\"title\">" << VoreenApplication::app()->getGuiName() << " " << outputFile_ << " </title>\n"
        << "<desc>Voreen PlotSVG-Output-File</desc>\n"
        << "<defs>\n"
        << "  <clipPath id=\"cp0\">\n"
        << "    <rect x=\"" << -1 << "\" y=\"" << -1
        << "\" width=\"" << windowSize_.x + 2
        << "\" height=\"" << windowSize_.y + 2 << "\" />\n"
        << "  </clipPath>\n"
        << "  <clipPath id=\"cp1\">\n"
        << "    <rect x=\"" << marginLeft_ << "\" y=\"" << marginBottom_
        << "\" width=\"" << windowSize_.x - (marginLeft_+marginRight_)
        << "\" height=\"" << windowSize_.y - (marginTop_+marginBottom_) << "\" />\n"
        << "  </clipPath>\n"
        << "  <clipPath id=\"cp2\">\n"
        << "    <rect x=\"" << marginLeft_ - maxGlyphSize_/2.0 << "\" y=\"" << marginBottom_  - maxGlyphSize_/2.0
        << "\" width=\"" << windowSize_.x - (marginLeft_+marginRight_) + maxGlyphSize_
        << "\" height=\"" << windowSize_.y - (marginTop_+marginBottom_) + maxGlyphSize_  << "\" />\n"
        << "  </clipPath>\n"
        << "</defs>\n"
        << "<g transform=\" translate(0 " << windowSize_.y << ") scale(1 -1)\">\n";
     //   << "@font-face {\n"
     //   << "  font-family: 'VoreenFontFamily';\n"
     //   << "  src: url('font')  format('truetype');\n"
     //   << "}\n";
    if (viewPortClipping_) {
        std::vector< Projection_Coordinates > points;
        tgt::Vector3<plot_t> coor = tgt::Vector3<plot_t>(0,0,0);
        points.push_back(Projection_Coordinates(coor,coor,coor));
        coor = tgt::Vector3<plot_t>(windowSize_.x,0,0);
        points.push_back(Projection_Coordinates(coor,coor,coor));
        coor = tgt::Vector3<plot_t>(windowSize_.x,windowSize_.y,0);
        points.push_back(Projection_Coordinates(coor,coor,coor));
        coor = tgt::Vector3<plot_t>(0,windowSize_.y,0);
        points.push_back(Projection_Coordinates(coor,coor,coor));
        PlotLibrarySvgRenderPolygon polygon = PlotLibrarySvgRenderPolygon(points,tgt::Color(),tgt::Color(1,1,0.995,0.005),0,1,0,PlotEntitySettings::FILLNONE);
        svgOutput_ << polygon.render();
    }
}

bool PlotLibrarySvg::PlotLibrarySvgSorter::operator() (PlotLibrarySvgRender* a, PlotLibrarySvgRender* b) const {
    PlotLibrarySvgRenderLine* line1 = 0;
    PlotLibrarySvgRenderLine* line2 = 0;
    PlotLibrarySvgRenderPolygon* polygon1 = 0;
    PlotLibrarySvgRenderPolygon* polygon2 = 0;
    line1 = dynamic_cast<PlotLibrarySvgRenderLine*>(a);
    line2 = dynamic_cast<PlotLibrarySvgRenderLine*>(b);
    polygon1 = dynamic_cast<PlotLibrarySvgRenderPolygon*>(a);
    polygon2 = dynamic_cast<PlotLibrarySvgRenderPolygon*>(b);
    bool switched = false;
    if ((line2 != 0) && (polygon1 != 0)) {
        line1 = line2;
        line2 = 0;
        polygon2 = polygon1;
        polygon1 = 0;
        switched = true;
    }
    if ((line1 != 0) && (line2 != 0)) {
        return a->avgValue(PlotLibraryFileBase::PROJECTION,PlotLibrary::Z_AXIS) > b->avgValue(PlotLibraryFileBase::PROJECTION,PlotLibrary::Z_AXIS);
    }
    else if ((line1 != 0) && (polygon2 != 0)) {
        //return a->avgValue(PlotLibraryFileBase::PROJECTION,PlotLibrary::Z_AXIS) > b->avgValue(PlotLibraryFileBase::PROJECTION,PlotLibrary::Z_AXIS);
        std::vector< tgt::Vector3<plot_t> > linePoints = line1->getPoints(PlotLibraryFileBase::PROJECTION);
        std::vector< tgt::Vector2<plot_t> > linePoints2D, polygonPoints2D;
        tgt::Vector2<plot_t> intersect;
        tgt::Vector3<plot_t> intersectPointPolygon;
        linePoints2D.push_back(linePoints[0].xy());
        linePoints2D.push_back(linePoints[1].xy());
        std::vector< tgt::Vector3<plot_t> > polygonPoints = polygon2->getPoints(PlotLibraryFileBase::PROJECTION);
        for (size_t i = 0; i < polygonPoints.size(); ++i) {
            polygonPoints2D.clear();
            polygonPoints2D.push_back(polygonPoints[i].xy());
            polygonPoints2D.push_back(polygonPoints[(i+1)%(polygonPoints.size()-1)].xy());
            if (PlotLibraryFileBase::intersect_Ray_Ray(linePoints2D,polygonPoints2D,&intersect) == 1) {
                if (intersect.x > 0 && intersect.x < 1 && intersect.y > 0 && intersect.y < 1) {
                    intersectPointPolygon = polygonPoints[i] + intersect.x* (polygonPoints[(i+1)%(polygonPoints.size()-1)]-polygonPoints[i]);
                    if (switched)
                        return intersectPointPolygon.z > line1->getDeepValueAtAlpha(PlotLibraryFileBase::PROJECTION,intersect.y);
                    else
                        return intersectPointPolygon.z < line1->getDeepValueAtAlpha(PlotLibraryFileBase::PROJECTION,intersect.y);
                }
            }
        }
        //intersect = tgt::Vector2<plot_t>(polygon2->avgValue(PlotLibraryFileBase::PROJECTION,PlotLibrary::X_AXIS),
        //    polygon2->avgValue(PlotLibraryFileBase::PROJECTION,PlotLibrary::Y_AXIS));
        //if (intersect.x >= line1->minValue(PlotLibraryFileBase::PROJECTION,PlotLibrary::X_AXIS) &&
        //    intersect.x <= line1->maxValue(PlotLibraryFileBase::PROJECTION,PlotLibrary::X_AXIS) &&
        //    intersect.y >= line1->minValue(PlotLibraryFileBase::PROJECTION,PlotLibrary::Y_AXIS) &&
        //    intersect.y <= line1->maxValue(PlotLibraryFileBase::PROJECTION,PlotLibrary::Y_AXIS)) {
        //    return line1->getDeepValueAtValue(PlotLibraryFileBase::PROJECTION,intersect) > b->avgValue(PlotLibraryFileBase::PROJECTION,PlotLibrary::Z_AXIS);
        //}
    }
    return a->avgValue(PlotLibraryFileBase::PROJECTION,PlotLibrary::Z_AXIS) > b->avgValue(PlotLibraryFileBase::PROJECTION,PlotLibrary::Z_AXIS);
}




void PlotLibrarySvg::afterRender() {
    if (dimension_ == THREE) {
        PlotLibrarySvgRenderLine* line;
        PlotLibrarySvgRenderPolygon* polygon;
        std::vector< tgt::Vector3<plot_t> > linePoints, polygonpoints, trianglePoints;
        std::vector< PlotLibraryFileBase::Projection_Coordinates > linecutPoints;
        tgt::Vector3<plot_t> result = tgt::Vector3<plot_t>();
        std::vector< PlotLibrarySvgRender* >::iterator svgit;
        std::vector< PlotLibrarySvgRender* > newLines, cutLines;
        svgit = renderVector_.begin();
        while (svgit < renderVector_.end()) {
            line = dynamic_cast<PlotLibrarySvgRenderLine*>(*svgit);
            if (line != 0) {
                linePoints = line->getPoints(PlotLibraryFileBase::ORIGIN);
                for (size_t j = 0; j < renderVector_.size(); ++j) {
                    polygon = dynamic_cast<PlotLibrarySvgRenderPolygon*>(renderVector_[j]);
                    if (polygon != 0) {
                        polygonpoints = polygon->getPoints(PlotLibraryFileBase::ORIGIN);
                        trianglePoints.clear();
                        if (polygonpoints.size() == 3) {
                            trianglePoints = polygonpoints;
                            if (intersect_Triangle_Ray(linePoints,trianglePoints,&result) == 1)
                                linecutPoints.push_back(projection_neu(result));
                        }
                        else if (polygonpoints.size() == 4) {
                            trianglePoints.push_back(polygonpoints[0]);
                            trianglePoints.push_back(polygonpoints[1]);
                            trianglePoints.push_back(polygonpoints[2]);
                            if (intersect_Triangle_Ray(linePoints,trianglePoints,&result) == 1)
                                linecutPoints.push_back(projection_neu(result));
                            else {
                                trianglePoints.clear();
                                trianglePoints.push_back(polygonpoints[0]);
                                trianglePoints.push_back(polygonpoints[2]);
                                trianglePoints.push_back(polygonpoints[3]);
                                if (intersect_Triangle_Ray(linePoints,trianglePoints,&result) == 1)
                                    linecutPoints.push_back(projection_neu(result));
                            }
                        }
                    }

                }
                if (linecutPoints.size() > 0) {
                    linecutPoints.push_back(line->getPoints()[0]);
                    linecutPoints.push_back(line->getPoints()[1]);
                    std::stable_sort(linecutPoints.begin(),linecutPoints.end(),PlotLibraryPointSorter);
                    cutLines = line->generateLines(linecutPoints);
                    for (size_t i = 0; i < cutLines.size(); ++i) {
                        newLines.push_back(cutLines[i]);
                    }
                    cutLines.clear();
                    delete line;
                    line = 0;
                    svgit = renderVector_.erase(svgit);
                    linecutPoints.clear();
                }
                else
                    ++svgit;
            }
            else
                ++svgit;
        }
        for (size_t i = 0; i < newLines.size(); ++i) {
            renderVector_.push_back(newLines[i]);
        }
        newLines.clear();
        std::sort(renderVector_.begin(),renderVector_.end(),plotLibrarySvgSorter);
    }
    for (size_t i = 0; i < renderVector_.size(); ++i)
        svgOutput_ << renderVector_[i]->render();
    svgOutput_ << "</g></svg>";
    std::ofstream svgfile (outputFile_.c_str(),std::ios_base::out | std::ios_base::binary | std::ios_base::trunc);
    svgfile << svgOutput_.str();
    svgfile.close();
    for (size_t i = 0; i < renderVector_.size(); ++i)
        delete renderVector_[i];
    renderVector_.clear();
}

bool PlotLibrarySvg::setRenderStatus() {
    if (usePlotPickingManager_)
        return false;
    // helper domains
    invertedmodelviewMatrix_ = tgt::Matrix4<plot_t>::zero;
    plot_t xl = domain_[X_AXIS].getLeft(); plot_t xr = domain_[X_AXIS].getRight();
    plot_t yl = domain_[Y_AXIS].getLeft(); plot_t yr = domain_[Y_AXIS].getRight();
    plot_t zl = domain_[Z_AXIS].getLeft(); plot_t zr = domain_[Z_AXIS].getRight();
    clippingPlanes_.clear();
    renewPlotToViewportScale();
    // for 2D plots setup orthographic projection using margins
    if (dimension_ == TWO || dimension_ == FAKETHREE) {
        //check, if the canvas is big enough
        if (marginLeft_+marginRight_>=windowSize_.x || marginTop_+marginBottom_>=windowSize_.y)
            return false;
        left_ = xl-static_cast<double>(marginLeft_)/plotToViewportScale_.x;
        right_ = xr+static_cast<double>(marginRight_)/plotToViewportScale_.x;
        bottom_ = yl-static_cast<double>(marginBottom_)/plotToViewportScale_.y;
        top_ = yr+static_cast<double>(marginTop_)/plotToViewportScale_.y;
        near_ = zl-1;
        far_ = zr+1;
        projectionsMatrix_ = plOrtho(left_,right_,bottom_,top_,near_,far_);

        modelviewMatrix_ = tgt::Matrix4<plot_t>::createIdentity();
        if (dimension_ == TWO) {
            //set up clipping planes
            tgt::dvec2 leftBottom = convertViewportToPlotCoordinates(tgt::ivec2(marginLeft_,marginBottom_));
            tgt::dvec2 rightTop = convertViewportToPlotCoordinates(tgt::ivec2(windowSize_.x-marginRight_,windowSize_.y-marginTop_));
            clippingPlanes_.push_back(tgt::Vector4d(1.0, 0.0, 0.0, -leftBottom.x));
            clippingPlanes_.push_back(tgt::Vector4d(-1.0, 0.0, 0.0, rightTop.x));
            clippingPlanes_.push_back(tgt::Vector4d(0.0, 1.0, 0.0, -leftBottom.y));
            clippingPlanes_.push_back(tgt::Vector4d(0.0, -1.0, 0.0, rightTop.y));
        }
        else { //dimension_ == FAKETHREE
            tgt::Matrix4<plot_t> m = tgt::Matrix4<plot_t>(1,0,0,0,0,1,0,0,-shear_.x*static_cast<double>(xr-xl),
                -shear_.y*static_cast<double>(yr-yl),1,0,0,0,0,1);
            //translation to the center
            ////rescaling
            ////shearing for fake 3d effect
            modelviewMatrix_ = tgt::Matrix4<plot_t>::createTranslation(tgt::Vector3<plot_t>((xl+xr)/2,(yl+yr)/2,(zl+zr)/2))
                * tgt::transpose(m)
                * tgt::Matrix4<plot_t>::createScale(tgt::Vector3<plot_t>(1/(1+shear_.x*static_cast<double>(zr-zl)),1/(1+shear_.y*static_cast<double>(zr-zl)),1))
                * tgt::Matrix4<plot_t>::createTranslation(tgt::Vector3<plot_t>(-(xl+xr)/2,-(yl+yr)/2,-(zl+zr)/2));
        }
    }
    // for 3D plot setup projection using the camera
    else if (dimension_ == THREE) {
        //if (lightingFlag_ == true && usePlotPickingManager_ == false) {
        //}
        modelviewMatrix_ = camera_.getViewMatrix();
        modelviewMatrix_ = modelviewMatrix_
            * tgt::Matrix4<plot_t>::createTranslation(tgt::Vector3<plot_t>(-0.5,-0.5,-0.5))
            * tgt::Matrix4<plot_t>::createScale(tgt::Vector3<plot_t>(1/domain_[0].size(),1/domain_[1].size(),1/domain_[2].size()))
            * tgt::Matrix4<plot_t>::createTranslation(tgt::Vector3<plot_t>(-domain_[0].getLeft(),-domain_[1].getLeft(),-domain_[2].getLeft()));
        if (orthographicCameraFlag_) {
            left_ = -0.9;
            right_ = 0.9;
            bottom_ = -0.9;
            top_ = 0.9;
            near_ = -2;
            far_ = 10;
            projectionsMatrix_ = plOrtho(left_,right_,bottom_,top_,near_,far_);
        }
        else {
            projectionsMatrix_ = camera_.getProjectionMatrix(1.0f); //TODO?
            camera_.updateFrustum();
            left_ = camera_.getFrustLeft();
            right_ = camera_.getFrustRight();
            bottom_ = camera_.getFrustBottom();
            top_ = camera_.getFrustTop();
            near_ = camera_.getNearDist();
            far_ = camera_.getFarDist();
        }
        clippingPlanes_.push_back(tgt::Vector4d(1.0, 0.0, 0.0, -domain_[X_AXIS].getLeft()));
        clippingPlanes_.push_back(tgt::Vector4d(-1.0, 0.0, 0.0, domain_[X_AXIS].getRight()));
        clippingPlanes_.push_back(tgt::Vector4d(0.0, 1.0, 0.0, -domain_[Y_AXIS].getLeft()));
        clippingPlanes_.push_back(tgt::Vector4d(0.0, -1.0, 0.0, domain_[Y_AXIS].getRight()));
        clippingPlanes_.push_back(tgt::Vector4d(0.0, 0.0, 1.0, -domain_[Z_AXIS].getLeft()));
        clippingPlanes_.push_back(tgt::Vector4d(0.0, 0.0, -1.0, domain_[Z_AXIS].getRight()));

        // we need to calculate which are the outer edges (used for selection and
        // labeling the axes)
        calculateSelectionEdges();
    }
    svgColor_ = drawingColor_;
    textureNumber_ = 0;
    modelviewMatrix_.invert(invertedmodelviewMatrix_);

    return true;
}

void PlotLibrarySvg::renderLine(const PlotData& data, int indexX, int indexY) {
    if (usePlotPickingManager_)
        return;
    //row iterator
    std::vector<PlotRowValue>::const_iterator it = data.getRowsBegin();
    // check if only values or only tags in given cells
    bool tagsInX = (data.getColumnType(indexX) == PlotBase::STRING);
    bool lineIsHighlighted = data.isHighlighted(tgt::ivec2(-1,indexY));
    PlotLibraryFileBase::Projection_Coordinates point1, point2;
    int i = 0;
    double x = 0.0; double y = 0.0; //they are set in the loop
    //go to the first row with non null entries
    while (it->getCellAt(indexX).isNull() || it->getCellAt(indexY).isNull()) {
        ++it; ++i;
    }
    if (it == data.getRowsEnd() || it == --data.getRowsEnd())
        return;

    // draw the line
    double oldX = tagsInX ? i : it->getValueAt(indexX);
    double oldY = it->getValueAt(indexY);
    if (lineIsHighlighted)
        svgColor_ = highlightColor_;
    else
        svgColor_ = drawingColor_;
    for (++it, ++i; it != data.getRowsEnd(); ++it, ++i) {
        //we ignore rows with null entries
        if (it->getCellAt(indexX).isNull() || it->getCellAt(indexY).isNull()) {
            continue;
        }
        x = tagsInX ? i : it->getValueAt(indexX);
        y = it->getValueAt(indexY);
        point1 = convertPlotCoordinatesToViewport3Projection(oldX,oldY,0);
        point2 = convertPlotCoordinatesToViewport3Projection(x,y,0);
        svgLine(point1,point2,svgColor_,1,lineWidth_,lineStyle_);

        //if x is out of the interval, we leave the loop
        if (x>= domain_[X_AXIS].getRight())
            break;
        oldX = x;
        oldY = y;
    }

    // render the points
    for (it = data.getRowsBegin(), i = 0; it != data.getRowsEnd(); ++it, ++i) {
        //we ignore rows with null entries
        if (it->getCellAt(indexX).isNull() || it->getCellAt(indexY).isNull()) {
            continue;
        }
        x = tagsInX ? i : it->getValueAt(indexX);
        y = it->getValueAt(indexY);
        if (it->getCellAt(indexY).isHighlighted())
            svgColor_ = highlightColor_;
        else
            svgColor_ = drawingColor_;
        point2 = convertPlotCoordinatesToViewport3Projection(x,y,0);
        svgCircle(point2,maxGlyphSize_/2,svgColor_,svgColor_,2,1);
    }
}

void PlotLibrarySvg::renderSpline(const PlotData& data, int indexX, int indexY) {
    if (usePlotPickingManager_)
        return;
    std::vector<PlotRowValue>::const_iterator it = data.getRowsBegin(); // row iterator
    tgt::Spline spline;
    std::vector< PlotLibraryFileBase::Projection_Coordinates > points;
    // check if only values or only tags in given cells
    bool tagsInX = (data.getColumnType(indexX) == PlotBase::STRING);
    bool lineIsHighlighted = data.isHighlighted(tgt::ivec2(-1,indexY));
    int i = 0;
    //go to the first row with non null entrys
    while (it->getCellAt(indexX).isNull() || it->getCellAt(indexY).isNull()) {
        ++it; ++i;
    }
    if (it == data.getRowsEnd())
        return;

    // add the control points
    double x = tagsInX ? i : it->getValueAt(indexX);
    double y = it->getValueAt(indexY);
    spline.addControlPoint(tgt::dvec3(x,y,0));
    for (; it != data.getRowsEnd(); ++it, ++i) {
        //we ignore rows with null entries
        if (it->getCellAt(indexX).isNull() || it->getCellAt(indexY).isNull()) {
            continue;
        }
        x = tagsInX ? i : it->getValueAt(indexX);
        y = it->getValueAt(indexY);
        spline.addControlPoint(tgt::dvec3(x,y,0));
    }
    spline.addControlPoint(tgt::dvec3(x,y,0));

    //render spline, possibly with log coordinates
    GLfloat step = 1.f /spline.getStepCount();
    if (lineIsHighlighted)
        svgColor_ = highlightColor_;
    else
        svgColor_ = drawingColor_;
    for (GLfloat p = 0.f; p < 1.f; p+=step) {
        x = spline.getPoint(p).x;
        y = spline.getPoint(p).y;
        points.push_back(convertPlotCoordinatesToViewport3Projection(x,y,0));
        if (x >= domain_[X_AXIS].getRight())
            break;
    }
    svgPolyline(points,svgColor_,svgColor_,1,lineWidth_,PlotEntitySettings::FILLNONE,lineStyle_);

    // render the points
    for (it = data.getRowsBegin(), i = 0; it != data.getRowsEnd(); ++it, ++i) {
        //we ignore rows with null entries
        if (it->getCellAt(indexX).isNull() || it->getCellAt(indexY).isNull()) {
            continue;
        }
        x = tagsInX ? i : it->getValueAt(indexX);
        y = it->getValueAt(indexY);
        if (it->getCellAt(indexY).isHighlighted())
            svgColor_ = highlightColor_;
        else
            svgColor_ = drawingColor_;
        svgCircle(convertPlotCoordinatesToViewport3Projection(x,y,0),maxGlyphSize_/2,svgColor_,svgColor_,2);
    }
}

void PlotLibrarySvg::renderErrorline(const PlotData& data, int indexX, int indexY, int indexError) {
    if (usePlotPickingManager_)
        return;
    std::vector<PlotRowValue>::const_iterator it = data.getRowsBegin();
    std::vector< PlotLibraryFileBase::Projection_Coordinates > points;
    // check if only values or only tags in given cells
    bool tagsInX = (data.getColumnType(indexX) == PlotBase::STRING);

    int i = 0;
    //go to the first row with non null entrys
    while (it->getCellAt(indexX).isNull() || it->getCellAt(indexY).isNull()) {
        ++it; ++i;
    }
    if (it == data.getRowsEnd())
        return;

    double x             = tagsInX ? i : it->getValueAt(indexX);
    double errorTop      = it->getValueAt(indexY) + std::abs(it->getValueAt(indexError));
    double errorBottom   = it->getValueAt(indexY) - std::abs(it->getValueAt(indexError));
    double oldX = 0.0;
    double oldErrorTop = 0.0;
    double oldErrorBottom = 0.0;
    svgColor_ = fillColor_;
    // draw the errorline
    ++i;
    for (it = (++it); it != data.getRowsEnd(); ++it, ++i) {
        if (it->getCellAt(indexX).isNull() || it->getCellAt(indexY).isNull())
            continue;
        if (i > 1) {
            points.clear();
            points.push_back(convertPlotCoordinatesToViewport3Projection(oldX,oldErrorTop,0));
            points.push_back(convertPlotCoordinatesToViewport3Projection(x,errorTop,0));
            points.push_back(convertPlotCoordinatesToViewport3Projection(x,errorBottom,0));
            points.push_back(convertPlotCoordinatesToViewport3Projection(oldX,oldErrorBottom,0));
            svgPolygon(points,svgColor_,svgColor_,1,0);
        }
        oldX = x;
        oldErrorTop = errorTop;
        oldErrorBottom = errorBottom;
        x = tagsInX ? i : it->getValueAt(indexX);
        errorTop = it->getValueAt(indexY) + std::abs(it->getValueAt(indexError));
        errorBottom = it->getValueAt(indexY) - std::abs(it->getValueAt(indexError));
    }
    points.clear();
    points.push_back(convertPlotCoordinatesToViewport3Projection(oldX,oldErrorTop,0));
    points.push_back(convertPlotCoordinatesToViewport3Projection(x,errorTop,0));
    points.push_back(convertPlotCoordinatesToViewport3Projection(x,errorBottom,0));
    points.push_back(convertPlotCoordinatesToViewport3Projection(oldX,oldErrorBottom,0));
    svgPolygon(points,svgColor_,svgColor_,1,0);
}

void PlotLibrarySvg::renderErrorspline(const PlotData& data, int indexX, int indexY, int indexError) {
    if (usePlotPickingManager_)
        return;
    std::vector<PlotRowValue>::const_iterator it = data.getRowsBegin();
    tgt::Spline splineTop, splineBottom;
    std::vector< PlotLibraryFileBase::Projection_Coordinates > points;

    // check if only values or only tags in given cells
    bool tagsInX = (data.getColumnType(indexX) == PlotBase::STRING);
    svgColor_ = fillColor_;
    int i = 0;
    //go to the first row with non null entrys
    while (it->getCellAt(indexX).isNull() || it->getCellAt(indexY).isNull()) {
        ++it; ++i;
    }
    if (it == data.getRowsEnd())
        return;

    //fill the splines

    double x             = tagsInX ? i : it->getValueAt(indexX);
    double errorTop      = it->getValueAt(indexY) + std::abs(it->getValueAt(indexError));
    double errorBottom   = it->getValueAt(indexY) - std::abs(it->getValueAt(indexError));
    splineTop.addControlPoint(tgt::dvec3(x,errorTop,0));
    splineBottom.addControlPoint(tgt::dvec3(x,errorBottom,0));
    for (; it != data.getRowsEnd(); ++it, ++i) {
        if (it->getCellAt(indexX).isNull() || it->getCellAt(indexY).isNull())
            continue;

        x = tagsInX ? i : it->getValueAt(indexX);
        if (it->getCellAt(indexError).isNull()){
            errorTop = it->getValueAt(indexY);
            errorBottom = it->getValueAt(indexY);
        } else {
            errorTop = it->getValueAt(indexY) + std::abs(it->getValueAt(indexError));
            errorBottom = it->getValueAt(indexY) - std::abs(it->getValueAt(indexError));
        }
        splineTop.addControlPoint(tgt::dvec3(x,errorTop,0));
        splineBottom.addControlPoint(tgt::dvec3(x,errorBottom,0));
    }
    splineTop.addControlPoint(tgt::dvec3(x,errorTop,0));
    splineBottom.addControlPoint(tgt::dvec3(x,errorBottom,0));

    //draw the thick spline
    GLfloat step = 1.f /splineTop.getStepCount();
    for (GLfloat p = 0.f; p < 1.f-step; p+=step) {
        points.clear();
        points.push_back(convertPlotCoordinatesToViewport3Projection(splineTop.getPoint(p)));
        points.push_back(convertPlotCoordinatesToViewport3Projection(splineBottom.getPoint(p)));
        points.push_back(convertPlotCoordinatesToViewport3Projection(splineBottom.getPoint(p+step)));
        points.push_back(convertPlotCoordinatesToViewport3Projection(splineTop.getPoint(p+step)));
        svgPolygon(points,svgColor_,svgColor_,1,0);
    }
}

void PlotLibrarySvg::renderErrorbars(const PlotData& data, int indexX, int indexY, int indexError) {
    if (usePlotPickingManager_)
        return;
    PlotLibraryFileBase::Projection_Coordinates p11;
    PlotLibraryFileBase::Projection_Coordinates p12;
    PlotLibraryFileBase::Projection_Coordinates p21;
    PlotLibraryFileBase::Projection_Coordinates p22;
    PlotLibraryFileBase::Projection_Coordinates p31;
    PlotLibraryFileBase::Projection_Coordinates p32;
    std::vector<PlotRowValue>::const_iterator it;
    // check if only values or only tags in given cells
    bool tagsInX = (data.getColumnType(indexX) == PlotBase::STRING);

    float radius = static_cast<float>(domain_[X_AXIS].size()/(4.f*data.getRowsCount()));
    //float aspectRatio = static_cast<float>(windowSize_.x)/static_cast<float>(windowSize_.y)*
    //    static_cast<float>(domain_[Y_AXIS].size() / domain_[X_AXIS].size());

    int i = 0;
    // draw the errorbars
    for (it = data.getRowsBegin(); it != data.getRowsEnd(); ++it, ++i) {
        if (it->getCellAt(indexX).isNull() || it->getCellAt(indexY).isNull() || it->getCellAt(indexError).isNull()) {
            continue;
        }
        double x = tagsInX ? i : it->getValueAt(indexX);
        if (!domain_[X_AXIS].contains(x))
            continue;

        if (it->getCellAt(indexError).isHighlighted())
            svgColor_ = highlightColor_;
        else
            svgColor_ = drawingColor_;

        double y = it->getValueAt(indexY);
        double yTop = it->getValueAt(indexY) + it->getValueAt(indexError);
        double yBottom = it->getValueAt(indexY) - it->getValueAt(indexError);
        x = logarithmicAxisFlags_[X_AXIS] ? convertToLogCoordinates(x, X_AXIS) : x;
        y = logarithmicAxisFlags_[Y_AXIS] ? convertToLogCoordinates(y, Y_AXIS) : y;
        yTop = logarithmicAxisFlags_[Y_AXIS] ? convertToLogCoordinates(yTop, Y_AXIS) : yTop;
        yBottom = logarithmicAxisFlags_[Y_AXIS] ? convertToLogCoordinates(yBottom, Y_AXIS) : yBottom;

        p11 = projection_neu(x,yBottom,0);
        p12 = projection_neu(x,yTop,0);
        svgLine(p11,p12,svgColor_,1,lineWidth_);

        p21 = projection_neu(x-radius,yTop,0);
        p22 = projection_neu(x+radius,yTop,0);
        svgLine(p21,p22,svgColor_,0,lineWidth_);

        p31 = projection_neu(x-radius,yBottom,0);
        p32 = projection_neu(x+radius,yBottom,0);
        svgLine(p31,p32,svgColor_,0,lineWidth_);

        //tgt::Ellipse midpoint(tgt::dvec3(x, y, domain_[2].getLeft()), radius/2,
        //                      (radius*aspectRatio)/2, tgt::dvec3(0, 0, 1), tgt::dvec3(1, 0, 0 ), 32);
        svgEllipse(projection_neu(x,y,0),(p22.afterProjection_.x-p21.afterProjection_.x)*radius,(p22.afterProjection_.x-p21.afterProjection_.x)*radius,svgColor_,svgColor_,0,3,PlotEntitySettings::FILLNONE);
    }
}

void PlotLibrarySvg::renderCandlesticks(const PlotData& data, int indexX, int stickTop,
                                     int stickBottom, int candleTop, int candleBottom) {
    if (usePlotPickingManager_)
        return;
    std::vector<PlotRowValue>::const_iterator it;

    // check if only values or only tags in given cells
    bool tagsInX = (data.getColumnType(indexX) == PlotBase::STRING);
    std::vector< PlotLibraryFileBase::Projection_Coordinates > points;

    float width = static_cast<float>(domain_[X_AXIS].size()/(4.f*data.getRowsCount()));

    // draw the candlestick
    int i=0;
    double yStickTop     = 0.0;
    double yStickBottom  = 0.0;
    double yCandleTop    = 0.0;
    double yCandleBottom = 0.0;

    for (it = data.getRowsBegin(); it != data.getRowsEnd(); ++it, ++i) {
        if (it->getCellAt(indexX).isNull() || it->getCellAt(stickTop).isNull() || it->getCellAt(stickBottom).isNull()
            || it->getCellAt(candleTop).isNull() || it->getCellAt(candleBottom).isNull()) {
            continue;
        }
        double x = tagsInX ? i : it->getValueAt(indexX);
        yStickTop     = it->getValueAt(stickTop);
        yStickBottom  = it->getValueAt(stickBottom);
        yCandleTop    = it->getValueAt(candleTop);
        yCandleBottom = it->getValueAt(candleBottom);

        // we divide the stick and the candle in top and bottom half
        // draw stick
        if (it->getCellAt(stickTop).isHighlighted())
            svgColor_ = highlightColor_;
        else
            svgColor_ = drawingColor_;
        svgLine(convertPlotCoordinatesToViewport3Projection(x,yStickTop,0),convertPlotCoordinatesToViewport3Projection(x,(yStickTop+yStickBottom)/2.0,0),svgColor_,1,lineWidth_);
        if (it->getCellAt(stickBottom).isHighlighted())
            svgColor_ = highlightColor_;
        else
            svgColor_ = drawingColor_;
        svgLine(convertPlotCoordinatesToViewport3Projection(x,(yStickTop+yStickBottom)/2.0,0),convertPlotCoordinatesToViewport3Projection(x,yStickBottom,0),svgColor_,1,lineWidth_);
    //    //draw candle
        if (it->getCellAt(candleTop).isHighlighted())
            svgColor_ = highlightColor_;
        else
            svgColor_ = fillColor_;
        points.clear();
        points.push_back(convertPlotCoordinatesToViewport3Projection(x-width,yCandleTop,0));
        points.push_back(convertPlotCoordinatesToViewport3Projection(x-width,(yCandleBottom+yCandleTop)/2.0,0));
        points.push_back(convertPlotCoordinatesToViewport3Projection(x+width,(yCandleBottom+yCandleTop)/2.0,0));
        points.push_back(convertPlotCoordinatesToViewport3Projection(x+width,yCandleTop,0));
        svgPolygon(points,svgColor_,svgColor_,1,lineWidth_);
        if (it->getCellAt(candleBottom).isHighlighted())
            svgColor_ = highlightColor_;
        else
            svgColor_ = fillColor_;
        points.clear();
        points.push_back(convertPlotCoordinatesToViewport3Projection(x-width,(yCandleBottom+yCandleTop)/2.0,0));
        points.push_back(convertPlotCoordinatesToViewport3Projection(x-width,yCandleBottom,0));
        points.push_back(convertPlotCoordinatesToViewport3Projection(x+width,yCandleBottom,0));
        points.push_back(convertPlotCoordinatesToViewport3Projection(x+width,(yCandleBottom+yCandleTop)/2.0,0));
        svgPolygon(points,svgColor_,svgColor_,1,lineWidth_);
    }
}

void PlotLibrarySvg::renderSurface(const PlotData& data, const std::vector<int>& triangleVertexIndices, bool wire,
                                int indexX, int indexY, int indexZ, int indexCM, bool wireonly)  {
    if (usePlotPickingManager_ || !wire)
        return;
    // security check: if count of edge indices is not a multiple of 3 abort
    if (triangleVertexIndices.size() % 3 != 0)
        return;

    Interval<plot_t> colInterval(0, 0);

    if ( indexCM != -1 ) {
        colInterval = data.getInterval(indexCM);
        if (colInterval.size() == 0)
            indexCM = -1;
    }
    if ( indexCM == -1)
        svgColor_ = drawingColor_;

    // draw the triangles
    std::set<int> renderedHighlights; // we want to render each plot label // highlight only once
    PlotLibraryFileBase::Projection_Coordinates p1;
    PlotLibraryFileBase::Projection_Coordinates p2;
    PlotLibraryFileBase::Projection_Coordinates p3;
    PlotLibraryFileBase::Projection_Coordinates p4;
    std::vector< PlotLibraryFileBase::Projection_Coordinates > points;
    for (std::vector<int>::const_iterator it = triangleVertexIndices.begin(); it < triangleVertexIndices.end(); it += 3) {
        points.clear();
        svgColor_ = drawingColor_;
        for (int i=0; i<3; ++i) {
            const PlotRowValue& row = data.getRow(*(it+i));

            if (indexCM != -1 ) {
                float c = static_cast<float>((row.getValueAt(indexCM) - colInterval.getLeft()) / colInterval.size());
                svgColor_ = colorMap_.getColorAtPosition(c);
            }
            p1 = convertPlotCoordinatesToViewport3Projection(row.getValueAt(indexX), row.getValueAt(indexY), row.getValueAt(indexZ));
            if (!wire && row.getCellAt(indexZ).isHighlighted()) {
                points.push_back(p1);
            }
            else
                points.push_back(p1);
        }
        if (wireonly)
            svgPolygon(points,svgColor_,lineColor_,0,lineWidth_,PlotEntitySettings::FILLNONE);
        else
            svgPolygon(points,svgColor_,lineColor_,0,lineWidth_);
    }
}

void PlotLibrarySvg::renderHeightmap(const voreen::PlotData& data, const std::vector< std::list< tgt::dvec2 > >& voronoiRegions,
                                  bool wire, int indexZ, int indexCM, bool wireonly) {
    if (usePlotPickingManager_ || !wire)
        return;
    Interval<plot_t> colInterval(0, 0);
    plot_t yMin = data.getInterval(2).getLeft();
    PlotLibraryFileBase::Projection_Coordinates point1;
    PlotLibraryFileBase::Projection_Coordinates point2;
    PlotLibraryFileBase::Projection_Coordinates point3;
    PlotLibraryFileBase::Projection_Coordinates point4;
    std::vector< PlotLibraryFileBase::Projection_Coordinates > points;

    if ( indexCM != -1 ) {
        colInterval = data.getInterval(indexCM);
        if (colInterval.size() == 0)
            indexCM = -1;
    }

    int row = 0;
    std::vector< std::list< tgt::dvec2 > >::const_iterator rit = voronoiRegions.begin();
    std::vector<PlotRowValue>::const_iterator pit = data.getRowsBegin();
    for (; rit < voronoiRegions.end(); ++rit, ++pit, ++row){
        if (rit->empty())
            continue;

        if (!wire && pit->getCellAt(indexZ).isHighlighted())
            svgColor_ = highlightColor_;
        else if (indexCM != -1 ) {
            float c = static_cast<float>((pit->getValueAt(indexCM) - colInterval.getLeft()) / colInterval.size());
            svgColor_ = colorMap_.getColorAtPosition(c);
        }
        else {
            svgColor_ = drawingColor_;
        }

        plot_t height = pit->getValueAt(indexZ);
        //clip it:
        if (height > domain_[Z_AXIS].getRight())
            height = domain_[Z_AXIS].getRight();
        if (yMin < domain_[Z_AXIS].getLeft())
            yMin = domain_[Z_AXIS].getLeft();
        if (height > yMin) {
            // render the sides of the pillar
            int i = 0;
            for (std::list< tgt::dvec2 >::const_iterator eit = rit->begin(); eit != rit->end(); ++eit) {
                if (i == 0) {
                    point3 = convertPlotCoordinatesToViewport3Projection(eit->x, eit->y, height);
                    point4 = convertPlotCoordinatesToViewport3Projection(eit->x, eit->y, yMin);
                }
                else {
                    point2 = point3;
                    point1 = point4;
                    point3 = convertPlotCoordinatesToViewport3Projection(eit->x, eit->y, height);
                    point4 = convertPlotCoordinatesToViewport3Projection(eit->x, eit->y, yMin);
                    points.clear();
                    points.push_back(point1);
                    points.push_back(point2);
                    points.push_back(point3);
                    points.push_back(point4);
                    if (wireonly)
                        svgPolygon(points,svgColor_,lineColor_,0,lineWidth_,PlotEntitySettings::FILLNONE);
                    else
                        svgPolygon(points,svgColor_,lineColor_,0,lineWidth_);
                }
                ++i;
            }
            point1 = point4;
            point4 = convertPlotCoordinatesToViewport3Projection(rit->begin()->x, rit->begin()->y, yMin);
            point2 = point3;
            point3 = convertPlotCoordinatesToViewport3Projection(rit->begin()->x, rit->begin()->y, height);
            points.clear();
            points.push_back(point1);
            points.push_back(point2);
            points.push_back(point3);
            points.push_back(point4);
            if (wireonly)
                svgPolygon(points,svgColor_,lineColor_,0,lineWidth_,PlotEntitySettings::FILLNONE);
            else
                svgPolygon(points,svgColor_,lineColor_,0,lineWidth_);
            points.clear();
            for (std::list< tgt::dvec2 >::const_iterator eit = rit->begin(); eit != rit->end(); ++eit) {
                point1 = convertPlotCoordinatesToViewport3Projection(eit->x, eit->y, height);
                points.push_back(point1);
            }
            if (wireonly)
                svgPolygon(points,svgColor_,lineColor_,0,lineWidth_,PlotEntitySettings::FILLNONE);
            else
                svgPolygon(points,svgColor_,lineColor_,0,lineWidth_);
            points.clear();
            for (std::list< tgt::dvec2 >::const_iterator eit = rit->begin(); eit != rit->end(); ++eit) {
                point1 = convertPlotCoordinatesToViewport3Projection(eit->x, eit->y, yMin);
                points.push_back(point1);
            }
            if (wireonly)
                svgPolygon(points,svgColor_,lineColor_,0,lineWidth_,PlotEntitySettings::FILLNONE);
            else
                svgPolygon(points,svgColor_,lineColor_,0,lineWidth_);
        }
    }
}

void PlotLibrarySvg::renderBars(const PlotData& data, std::vector<int> indicesY) {
    if (usePlotPickingManager_)
        return;
    std::vector<PlotRowValue>::const_iterator rowIt = data.getRowsBegin();
    //stores y values and indices of a merged bar group, the indices are used to get the right color
    std::vector<std::pair<plot_t, tgt::Color> > mergedBars;
    //stores last y value, used for stacked bars
    plot_t lastY;
    //rowcounter
    plot_t row = 0;
    for (; rowIt < data.getRowsEnd(); ++rowIt) {
        lastY = 0;
        mergedBars.clear();
        // we do not use the iterator because we also iterate through the colormap
        for (size_t i = 0; i < indicesY.size(); ++i) {
            if (rowIt->getCellAt(indicesY.at(i)).isHighlighted())
                svgColor_ = highlightColor_;
            else
                svgColor_ = colorMap_.getColorAtIndex(static_cast<int>(i));
            if (barMode_ == STACKED) {
                if (rowIt->getCellAt(indicesY.at(i)).isNull() )
                    continue;
                plot_t newY = lastY+rowIt->getValueAt(indicesY.at(i));
                //negative stacked bars do not make any sense and are ignored
                if (newY < lastY)
                    continue;
                renderSingleBar(row-barWidth_/2.0,row+barWidth_/2.0, lastY, newY, svgColor_);
                lastY = newY;
            }
            else if (barMode_ == GROUPED) {
                if (rowIt->getCellAt(indicesY.at(i)).isNull())
                    continue;
                double singleBarWidth = barWidth_/(1.0*static_cast<double>(indicesY.size()));
                renderSingleBar(row-barWidth_/2.0+static_cast<double>(i)*singleBarWidth,
                    row-barWidth_/2.0+static_cast<double>(i+1)*singleBarWidth, 0, rowIt->getValueAt(indicesY.at(i)), svgColor_);
            }

            else { // MERGED
                // we can't skip null entries, so we set them 0
                if (rowIt->getCellAt(indicesY.at(i)).isNull())
                    mergedBars.push_back(std::pair<plot_t, tgt::Color>(0, svgColor_));
                // push the y value and the color index
                mergedBars.push_back(std::pair<plot_t, tgt::Color>(rowIt->getValueAt(indicesY.at(i)), svgColor_));
            }
        }
        if (barMode_ == MERGED) {
            // the values are stored in bars, but not yet drawn
            std::sort(mergedBars.begin(), mergedBars.end(), mergedBarSorter);
            std::vector<std::pair<plot_t, tgt::Color> >::const_iterator it;
            std::vector<std::pair<plot_t, tgt::Color> >::const_iterator lastit;
            double squeeze = 1.0;
            for (it = mergedBars.begin(); it < mergedBars.end(); ++it) {
                svgColor_ = it->second;
                if (it == mergedBars.begin()) {
                    renderSingleBar(row-barWidth_/2.0,row+barWidth_/2.0, 0, it->first, it->second, squeeze);
                }
                else
                    renderSingleBar(row-barWidth_/2.0,row+barWidth_/2.0, lastit->first, it->first, it->second, squeeze);
                squeeze = squeezeFactor_*squeeze;
                lastit = it;
            }
        }
        ++row;
    }
}

void PlotLibrarySvg::renderNodeGraph(const PlotData& nodeData, const PlotData& connectionData, int indexX, int indexY, int indexDx, int indexDy) {
    if (usePlotPickingManager_)
        return;
    std::vector<PlotRowValue>::const_iterator it;

    // render nodes
    plot_t glyphSize = (maxGlyphSize_ + minGlyphSize_)/2;
    int i = 0;
    for (it = nodeData.getRowsBegin(); it != nodeData.getRowsEnd(); ++it, ++i) {
        // render node
        svgColor_ = drawingColor_;
        renderGlyph(it->getValueAt(indexX), it->getValueAt(indexY), 0, glyphSize);

        // render force vector
        svgLine(convertPlotCoordinatesToViewport3Projection(it->getValueAt(indexX),it->getValueAt(indexY),0),
            convertPlotCoordinatesToViewport3Projection(it->getValueAt(indexX) + it->getValueAt(indexDx),
            it->getValueAt(indexY) + it->getValueAt(indexDy),0),
            fillColor_,1,lineWidth_);

        // render node label
        std::stringstream ss;
        ss << i;
        renderLabel(tgt::dvec3(it->getValueAt(indexX), it->getValueAt(indexY), 0), SmartLabel::CENTERED, ss.str(), false, 0);
    }

    std::vector<PlotCellValue> tester;
    tester.push_back(PlotCellValue(0));
    std::vector<PlotRowValue>::const_iterator firstIt, secondIt;

    // render connections
    for (it = connectionData.getRowsBegin(); it != connectionData.getRowsEnd(); ++it) {
        tester[0].setValue(static_cast<int>(it->getValueAt(0)));
        firstIt = nodeData.lower_bound(tester);

        tester[0].setValue(static_cast<int>(it->getValueAt(1)));
        secondIt = nodeData.lower_bound(tester);

        if (firstIt != nodeData.getRowsEnd() && secondIt != nodeData.getRowsEnd()) {
            svgLine(convertPlotCoordinatesToViewport3Projection(firstIt->getValueAt(indexX),firstIt->getValueAt(indexY),0),
                convertPlotCoordinatesToViewport3Projection(secondIt->getValueAt(indexX),secondIt->getValueAt(indexY),0),
                drawingColor_,1,lineWidth_);
        }
    }
}

void PlotLibrarySvg::renderColorMapLegend(const PlotData& data, int column, int number) {
    tgt::Color c1, c2;
    std::stringstream ss;
    // ColorMaps with less than 2 colors may not exist
    if (colorMap_.getColorCount() < 2)
        return;

    Interval<plot_t>   interval = data.getInterval(column);
    const std::string& label    = data.getColumnLabel(column);

    // switch to viewport coordinates

    // render legend
    int colorcount      = colorMap_.getColorCount();
    const double width  = 96;
    const double height = 16;
    const double xslide = 8;
    const double xStart = windowSize_.x - width - xslide;
    const double yStart = windowSize_.y - xslide - (number*32) - height;
    double stepWidth    = width / (colorcount - 1);

    // color map
    for (int i = 0; i < colorcount-1; ++i) {
        c1 = colorMap_.getColorAtIndex(i);
        c2 = colorMap_.getColorAtIndex(i+1);
        ss.clear();
        ss.str("");
        ss << "colormap" << i;
        svgGradient(tgt::Vector2<plot_t>(0,0),tgt::Vector2<plot_t>(1,0),ss.str(),c1,c2);
        tgt::Vector3<plot_t> startpoint = tgt::Vector3<plot_t>(xStart + (i*stepWidth),yStart,0);
        PlotLibraryFileBase::Projection_Coordinates startPoint = PlotLibraryFileBase::Projection_Coordinates(startpoint,startpoint,startpoint);
        svgRect(startPoint,stepWidth,height,ss.str());
    }

    svgColor_ = drawingColor_;
    std::vector< PlotLibraryFileBase::Projection_Coordinates > points;
    tgt::Vector3<plot_t> coor = tgt::Vector3<plot_t>(xStart,yStart,0);
    points.push_back(Projection_Coordinates(coor,coor,coor));
    coor = tgt::Vector3<plot_t>(xStart,yStart + height,0);
    points.push_back(Projection_Coordinates(coor,coor,coor));
    coor = tgt::Vector3<plot_t>(xStart + width,yStart + height,0);
    points.push_back(Projection_Coordinates(coor,coor,coor));
    coor = tgt::Vector3<plot_t>(xStart + width,yStart,0);
    points.push_back(Projection_Coordinates(coor,coor,coor));
    coor = tgt::Vector3<plot_t>(xStart,yStart,0);
    points.push_back(Projection_Coordinates(coor,coor,coor));
    svgPolyline(points,svgColor_,svgColor_,0,lineWidth_);

    // labels
    tgt::dvec3 position = tgt::dvec3(xStart, yStart + (height/2),0);
    SmartLabelGroupBaseSvg::renderSingleLabel(&labelFont_, label,fontSize_,SmartLabel::MIDDLERIGHT, 4,position);
    svgText(position,label,&labelFont_,fontSize_,fontColor_,SmartLabel::MIDDLERIGHT);

    ss.str("");
    ss.clear();
    ss << std::setprecision(4) << interval.getLeft();
    position = tgt::dvec3(xStart + xslide, yStart, 0);
    SmartLabelGroupBaseSvg::renderSingleLabel(&labelFont_, ss.str(),fontSize_, SmartLabel::BOTTOMRIGHT, 4,position);
    svgText(position,ss.str(),&labelFont_,fontSize_,fontColor_,SmartLabel::BOTTOMRIGHT);

    ss.str("");
    ss.clear();
    ss << std::setprecision(4) << interval.getRight();
    position = tgt::dvec3(xStart + xslide + width, yStart, 0);
    SmartLabelGroupBaseSvg::renderSingleLabel(&labelFont_, ss.str(),fontSize_, SmartLabel::BOTTOMRIGHT,4,position);
    svgText(position,ss.str(),&labelFont_,fontSize_,fontColor_,SmartLabel::BOTTOMRIGHT);

}

void PlotLibrarySvg::renderScatter(const PlotData& data, int indexX, int indexY, int indexZ, int indexCM, int indexSize) {
    if (usePlotPickingManager_)
        return;
    Interval<plot_t> colInterval(0, 0);
    if ( indexCM != -1 ) {
        colInterval = data.getInterval(indexCM);
        if (colInterval.size() == 0)
            indexCM = -1;
    }
    plot_t x, y, z;
    plot_t size = maxGlyphSize_;
    // if there is size data, we interpolate it to [minGlyphSize_,maxGlyphSize_]
    Interval<plot_t> sizeInterval(0, 0);
    if ( indexSize != -1 ) {
        sizeInterval = data.getInterval(indexSize);
        //we can only use an interval with positive size
        if (sizeInterval.size() == 0)
            indexSize = -1;
    }

    if (texture_ != 0) {
        ++textureNumber_;
        svgOutput_ << "<defs><pattern id=\"texture" << textureNumber_ << "\" x=\"0\" y=\"0\" width=\"" << texture_->getDimensions().x << "\" height=\"" << texture_->getDimensions().y << "\" patternUnits=\"userSpaceOnUse\">\n"
            << "  <image x=\"0\" y=\"0\" width=\"" << texture_->getDimensions().x << "\" height=\"" << texture_->getDimensions().y << "\" xlink:href=\"" << texturePath_ << "\" />\n"
            << "</pattern></defs>\n";
    }


    //row iterator
    int i = 0;
    std::vector<PlotRowValue>::const_iterator it = data.getRowsBegin();

    for (; it != data.getRowsEnd(); ++it, ++i) {
        //we ignore rows with null entries
        if (it->getCellAt(indexX).isNull() || it->getCellAt(indexY).isNull()) {
            continue;
        }
        x = it->getValueAt(indexX); y = it->getValueAt(indexY); z = (indexZ == -1 ? 0 : it->getValueAt(indexZ));
        //check if the point is inside the domains
        if (domain_[X_AXIS].contains(x) && domain_[Y_AXIS].contains(y) && (dimension_ == TWO || domain_[Z_AXIS].contains(z))) {
            // set color
            if ((indexZ != -1 && it->getCellAt(indexZ).isHighlighted())
                    || (indexZ == -1 && it->getCellAt(indexY).isHighlighted()))
                svgColor_ = highlightColor_;
            else if (indexCM != -1 ) {
                float c = static_cast<float>((it->getValueAt(indexCM) - colInterval.getLeft()) / colInterval.size());
                svgColor_ = colorMap_.getColorAtPosition(c);
            }
            else
                svgColor_ = drawingColor_;
            // set size
            if (indexSize != -1 ) {
                size = minGlyphSize_ + (maxGlyphSize_ - minGlyphSize_) *
                            (it->getValueAt(indexSize) - sizeInterval.getLeft()) / sizeInterval.size();
            }
            renderGlyph(x, y, z, size);
        }
    }
}

void PlotLibrarySvg::renderAxes() {
    // axes

     svgColor_ = drawingColor_;

    plot_t xl = domain_[0].getLeft();    plot_t xr = domain_[0].getRight();
    plot_t yl = domain_[1].getLeft();    plot_t yr = domain_[1].getRight();
    plot_t zl = domain_[2].getLeft();    plot_t zr = domain_[2].getRight();

    PlotLibraryFileBase::Projection_Coordinates p1;
    PlotLibraryFileBase::Projection_Coordinates p2;
    PlotLibraryFileBase::Projection_Coordinates p3;
    PlotLibraryFileBase::Projection_Coordinates p4;
    PlotLibraryFileBase::Projection_Coordinates p5;
    PlotLibraryFileBase::Projection_Coordinates p6;
    std::vector< PlotLibraryFileBase::Projection_Coordinates > points;

    if (! centerAxesFlag_) {
        //x and y axes
        p1 = projection_neu(xl, yr, zl);
        p2 = projection_neu(xl, yl, zl);
        p3 = projection_neu(xr, yl, zl);
        svgLine(p1,p2,svgColor_,0,axesWidth_);
        svgLine(p2,p3,svgColor_,0,axesWidth_);
    }

    if (dimension_ == TWO) {
    //draw arrows with viewport coordinates
        int arrowSize = 5; // in pixel
        tgt::Vector3<plot_t> p11 = tgt::Vector3<plot_t>(windowSize_.x-marginRight_,marginBottom_,0);
        tgt::Vector3<plot_t> p12 = tgt::Vector3<plot_t>(windowSize_.x-marginRight_+4*arrowSize,marginBottom_,0);
        p1 = PlotLibraryFileBase::Projection_Coordinates(p11,p11,p11);
        p2 = PlotLibraryFileBase::Projection_Coordinates(p12,p12,p12);
        svgLine(p1,p2,svgColor_,0,axesWidth_);

        p11 = tgt::Vector3<plot_t>(windowSize_.x-marginRight_+2*arrowSize,marginBottom_+arrowSize,0);
        p12 = tgt::Vector3<plot_t>(windowSize_.x-marginRight_+4*arrowSize,marginBottom_,0);
        p1 = PlotLibraryFileBase::Projection_Coordinates(p11,p11,p11);
        p2 = PlotLibraryFileBase::Projection_Coordinates(p12,p12,p12);
        svgLine(p1,p2,svgColor_,0,axesWidth_);

        p11 = tgt::Vector3<plot_t>(windowSize_.x-marginRight_+2*arrowSize,marginBottom_-arrowSize,0);
        p12 = tgt::Vector3<plot_t>(windowSize_.x-marginRight_+4*arrowSize,marginBottom_,0);
        p1 = PlotLibraryFileBase::Projection_Coordinates(p11,p11,p11);
        p2 = PlotLibraryFileBase::Projection_Coordinates(p12,p12,p12);
        svgLine(p1,p2,svgColor_,0,axesWidth_);

        p11 = tgt::Vector3<plot_t>(marginLeft_,windowSize_.y-marginTop_,0);
        p12 = tgt::Vector3<plot_t>(marginLeft_,windowSize_.y-marginTop_+4*arrowSize,0);
        p1 = PlotLibraryFileBase::Projection_Coordinates(p11,p11,p11);
        p2 = PlotLibraryFileBase::Projection_Coordinates(p12,p12,p12);
        svgLine(p1,p2,svgColor_,0,axesWidth_);

        p11 = tgt::Vector3<plot_t>(marginLeft_ + arrowSize,windowSize_.y-marginTop_+2*arrowSize,0);
        p12 = tgt::Vector3<plot_t>(marginLeft_,windowSize_.y-marginTop_+4*arrowSize,0);
        p1 = PlotLibraryFileBase::Projection_Coordinates(p11,p11,p11);
        p2 = PlotLibraryFileBase::Projection_Coordinates(p12,p12,p12);
        svgLine(p1,p2,svgColor_,0,axesWidth_);

        p11 = tgt::Vector3<plot_t>(marginLeft_ - arrowSize,windowSize_.y-marginTop_+2*arrowSize,0);
        p12 = tgt::Vector3<plot_t>(marginLeft_,windowSize_.y-marginTop_+4*arrowSize,0);
        p1 = PlotLibraryFileBase::Projection_Coordinates(p11,p11,p11);
        p2 = PlotLibraryFileBase::Projection_Coordinates(p12,p12,p12);
        svgLine(p1,p2,svgColor_,0,axesWidth_);
    }
    else if (dimension_ == THREE) {
    //draw cube mesh
        p1 = projection_neu(xr, yl, zl);
        p2 = projection_neu(xr, yr, zl);
        p3 = projection_neu(xl, yr, zl);
        p4 = projection_neu(xl, yr, zr);
        p5 = projection_neu(xr, yr, zr);
        p6 = projection_neu(xr, yr, zl);
        svgLine(p1,p2,svgColor_,0,axesWidth_);
        svgLine(p2,p3,svgColor_,0,axesWidth_);
        svgLine(p3,p4,svgColor_,0,axesWidth_);
        svgLine(p4,p5,svgColor_,0,axesWidth_);
        svgLine(p5,p6,svgColor_,0,axesWidth_);

        p1 = projection_neu(xr, yl, zl);
        p2 = projection_neu(xr, yl, zr);
        p3 = projection_neu(xl, yl, zr);
        p4 = projection_neu(xl, yl, zl);
        svgLine(p1,p2,svgColor_,0,axesWidth_);
        svgLine(p2,p3,svgColor_,0,axesWidth_);
        svgLine(p3,p4,svgColor_,0,axesWidth_);

        p1 = projection_neu(xr, yl, zr);
        p2 = projection_neu(xr, yr, zr);
        p3 = projection_neu(xl, yl, zr);
        p4 = projection_neu(xl, yr, zr);
        svgLine(p1,p2,svgColor_,0,axesWidth_);
        svgLine(p3,p4,svgColor_,0,axesWidth_);
    }
    else if (dimension_ == FAKETHREE) {
        //draw back
        p1 = projection_neu(xl, yr, zl);
        p2 = projection_neu(xl, yl, zl);
        p3 = projection_neu(xr, yl, zl);
        p4 = projection_neu(xr, yr, zl);
        points.clear();
        points.push_back(p1);
        points.push_back(p2);
        points.push_back(p3);
        points.push_back(p4);
        svgPolygon(points,svgColor_,svgColor_,0,axesWidth_,PlotEntitySettings::FILLNONE);

        //draw bottom
        p1 = projection_neu(xl, yl, zl);
        p2 = projection_neu(xl, yl, zr);
        svgLine(p1,p2,svgColor_,0,axesWidth_);

        p1 = projection_neu(xr, yl, zr);
        p2 = projection_neu(xr, yl, zl);
        svgLine(p1,p2,svgColor_,0,axesWidth_);

        //draw left
        p1 = projection_neu(xl, yr, zr);
        p2 = projection_neu(xl, yl, zr);
        p3 = projection_neu(xl, yl, zl);
        p4 = projection_neu(xl, yr, zl);
        points.clear();
        points.push_back(p1);
        points.push_back(p2);
        points.push_back(p3);
        points.push_back(p4);
        svgPolygon(points,svgColor_,svgColor_,0,axesWidth_,PlotEntitySettings::FILLNONE);

        //draw zero
        p1 = projection_neu(xl, 0, zl);
        p2 = projection_neu(xl, 0, zr);
        svgLine(p1,p2,svgColor_,0,axesWidth_);

        p1 = projection_neu(xr, 0, zr);
        p2 = projection_neu(xr, 0, zl);
        svgLine(p1,p2,svgColor_,0,axesWidth_);

        //the front is always above the plot
        p1 = projection_neu(xl, 0, zr);
        p2 = projection_neu(xr, 0, zr);
        svgLine(p1,p2,svgColor_,0,axesWidth_);

        p1 = projection_neu(xl, yl, zr);
        p2 = projection_neu(xr, yl, zr);
        svgLine(p1,p2,svgColor_,0,axesWidth_);
    }
    //glEnable(GL_DEPTH_TEST);

}

void PlotLibrarySvg::renderAxisScales(Axis axis, bool helperLines, const std::string& label, plot_t offset) {
    tgt::dvec2 step = updateScaleSteps(axis);
    svgColor_ = drawingColor_;
    plot_t xl = domain_[0].getLeft();    plot_t xr = domain_[0].getRight();
    plot_t yl = domain_[1].getLeft();    plot_t yr = domain_[1].getRight();
    plot_t zl = domain_[2].getLeft();    plot_t zr = domain_[2].getRight();

    std::stringstream stream;
    if (step.x < 1) {
        int precision = static_cast<int>(ceil(log10(1.0/step.x)));
        stream << std::fixed << std::setprecision(precision);
    }

    if (dimension_ == TWO || dimension_ == FAKETHREE) {
        xAxisLabelGroup_.reset();
        xAxisLabelGroup_.setBounds(getBoundsBelowPlot());
        // if respective scaleStep is e.g. 50 and interval begin e.g. 27.43, we want to have labels like
        // 50, 100, 150, ... instead of 27.43, 77.43, 127.43, ...
        // So we do some smart rounding:
        plot_t start = domain_[axis].getLeft();
        start = ceil(start / step.x) * step.x;

        for (plot_t i = start; i <= domain_[axis].getRight(); i += step.x) {
            plot_t pos = i;//(logarithmicAxisFlags_[axis] ? convertFromLogCoordinates(i, axis) : i);
            stream.str("");
            stream.clear();
            stream << round(pos + offset, step.y);
            if (axis == X_AXIS) {
                xAxisLabelGroup_.addLabel(stream.str(),
                                          convertPlotCoordinatesToViewport3(tgt::dvec3(i, yl, zr)),
                                          fontColor_, fontSize_, SmartLabel::TOPCENTERED);
                if (helperLines && dimension_ == TWO) {
                    PlotLibraryFileBase::Projection_Coordinates p1 = convertPlotCoordinatesToViewport3Projection(i, yl,0);
                    PlotLibraryFileBase::Projection_Coordinates p2 = convertPlotCoordinatesToViewport3Projection(i, yr,0);
                    svgLine(p1,p2,svgColor_,0,axesWidth_/2.f);
                }
            }
            else if (axis == Y_AXIS) {
                renderLabel(tgt::dvec3(xl, i, zr), SmartLabel::MIDDLERIGHT, stream.str());
                if (helperLines) {
                    PlotLibraryFileBase::Projection_Coordinates p2 = convertPlotCoordinatesToViewport3Projection(xl, i, zl);
                    PlotLibraryFileBase::Projection_Coordinates p3 = convertPlotCoordinatesToViewport3Projection(xr, i, zl);
                    if (dimension_ == FAKETHREE) {
                        PlotLibraryFileBase::Projection_Coordinates p1 = convertPlotCoordinatesToViewport3Projection(xl, i, zr);
                        svgLine(p1,p2,svgColor_,0,axesWidth_/2.f);
                    }
                    svgLine(p2,p3,svgColor_,0,axesWidth_/2.f);
                }
            }
        }
        renderSmartLabelGroup(&xAxisLabelGroup_);
    }
    else if (dimension_ == THREE) {
        // If we are inside the plot cube (or really close to it) we do not want to
        // render scales because it looks ugly it doesn't make any sense.
        if (!orthographicCameraFlag_ && tgt::distance(camera_.getPosition(), tgt::vec3(0,0,0)) < 1)
            return;

        axisLabelGroup_.reset();

        // avoid expensive copying by using iterators
        std::vector<SelectionEdge>::const_iterator minEdge = selectionEdgesX_.begin();
        std::vector<SelectionEdge>::const_iterator endEdge = selectionEdgesX_.end();
        if (axis == X_AXIS) {
            if (selectionEdgesX_.empty())
                return;
        }
        else if (axis == Y_AXIS) {
            if (selectionEdgesY_.empty())
                return;
            minEdge = selectionEdgesY_.begin();
            endEdge = selectionEdgesY_.end();
        }
        else if (axis == Z_AXIS) {
            if (selectionEdgesZ_.empty())
                return;
            minEdge = selectionEdgesZ_.begin();
            endEdge = selectionEdgesZ_.end();
        }

        // find edge with maximum length
        double length = tgt::length(minEdge->endVertex_ - minEdge->startVertex_);

        if (! orthographicCameraFlag_) {
            for (std::vector<SelectionEdge>::const_iterator it = ++minEdge; it < endEdge; ++it) {
                double val = tgt::length(it->endVertex_ - it->startVertex_);
                if (val > length) {
                    minEdge = it;
                    length = val;
                }
            }
        }

        // determine on which side the cube our axis is and with that the label alignment
        tgt::dvec2 edgeDirection = minEdge->endVertex_ - minEdge->startVertex_;
        tgt::dvec2 center(windowSize_.x/2, windowSize_.y/2);
        tgt::dvec2 ray((minEdge->startVertex_ + 0.5*edgeDirection) - center);
        ray = 1.0/tgt::length(ray) * ray;
        double angle = atan2(ray.y, ray.x);

        SmartLabel::Alignment align = SmartLabel::MIDDLELEFT;
        if (angle > tgt::PI/8 && angle <= 3.0*tgt::PI/8)
            align = SmartLabel::TOPLEFT;
        else if (angle > 3.0*tgt::PI/8 && angle <= 5.0*tgt::PI/8)
            align = SmartLabel::TOPCENTERED;
        else if (angle > 5.0*tgt::PI/8 && angle <= 7.0*tgt::PI/8)
            align = SmartLabel::TOPRIGHT;
        else if (angle < -7.0*tgt::PI/8 || angle > 7.0*tgt::PI/8)
            align = SmartLabel::MIDDLERIGHT;

        if (angle < -tgt::PI/8 && angle >= -3.0*tgt::PI/8)
            align = SmartLabel::BOTTOMLEFT;
        else if (angle < -3.0*tgt::PI/8 && angle >= -5.0*tgt::PI/8)
            align = SmartLabel::BOTTOMCENTERED;
        else if (angle < -5.0*tgt::PI/8 && angle >= -7.0*tgt::PI/8)
            align = SmartLabel::BOTTOMRIGHT;

        // render axis label first:
        if (label != "") {
            axisLabelGroup_.addLabel(label, minEdge->startVertex_ + (0.5 * edgeDirection) + (32.0 * ray), fontColor_, fontSize_, align);
        }

        // now render scales:
        // if respective scaleStep is e.g. 50 and interval begin e.g. 27.43, we want to have labels like
        // 50, 100, 150, ... instead of 27.43, 77.43, 127.43, ...
        // So we do some smart rounding:
        plot_t start = domain_[axis].getLeft();
        start = ceil(start / step.x) * step.x;
        plot_t domainSize = domain_[axis].size();

        for (plot_t i = start; i  <= domain_[axis].getRight(); i += step.x) {
            stream.str("");
            stream.clear();
            stream << round(i, step.y);
            if (minEdge->ascOrientation_)
                axisLabelGroup_.addLabel(stream.str(),
                        minEdge->startVertex_ + ((i - domain_[axis].getLeft())/domainSize)*edgeDirection,
                        fontColor_, fontSize_, align);
            else
                axisLabelGroup_.addLabel(stream.str(),
                        minEdge->endVertex_ - ((i - domain_[axis].getLeft())/domainSize)*edgeDirection,
                        fontColor_, fontSize_, align);

        }
        renderSmartLabelGroup(&axisLabelGroup_);

        // render helper lines
        if (helperLines) {
            PlotLibraryFileBase::Projection_Coordinates p1;
            PlotLibraryFileBase::Projection_Coordinates p2;
            PlotLibraryFileBase::Projection_Coordinates p3;
            PlotLibraryFileBase::Projection_Coordinates p4;
            std::vector< PlotLibraryFileBase::Projection_Coordinates > points;
            tgt::dvec2 xEdgeDirection;
            tgt::dvec2 yEdgeDirection;
            tgt::dvec2 zEdgeDirection;
            bool base, top, back, front, left, right;
            if (orthographicCameraFlag_) {
                xEdgeDirection = selectionEdgesX_.begin()->endVertex_ - selectionEdgesX_.begin()->startVertex_;
                yEdgeDirection = selectionEdgesY_.begin()->endVertex_ - selectionEdgesY_.begin()->startVertex_;
                //zEdgeDirection = selectionEdgesZ_.begin()->endVertex_ - selectionEdgesZ_.begin()->startVertex_;
                base = ((xEdgeDirection.x > 0.0  && xEdgeDirection.y > 0.0  && yEdgeDirection.x > 0.0  && yEdgeDirection.y < 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesY_.begin()->ascOrientation_ == false)) ||
                    (xEdgeDirection.x > 0.0  && xEdgeDirection.y < 0.0  && yEdgeDirection.x > 0.0  && yEdgeDirection.y > 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesY_.begin()->ascOrientation_ == true)) ||
                    (xEdgeDirection.x < 0.0  && xEdgeDirection.y < 0.0  && yEdgeDirection.x > 0.0  && yEdgeDirection.y < 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesY_.begin()->ascOrientation_ == true)) ||
                    (xEdgeDirection.x < 0.0 && xEdgeDirection.y > 0.0  && yEdgeDirection.x > 0.0  && yEdgeDirection.y > 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesY_.begin()->ascOrientation_ == false)));
                top = (!base);
                back = ((xEdgeDirection.x > 0.0  && xEdgeDirection.y < 0.0  && yEdgeDirection.x > 0.0  && yEdgeDirection.y > 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesY_.begin()->ascOrientation_ == true)) ||
                    (xEdgeDirection.x > 0.0  && xEdgeDirection.y > 0.0  && yEdgeDirection.x > 0.0  && yEdgeDirection.y < 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesY_.begin()->ascOrientation_ == false)) ||
                    (xEdgeDirection.x < 0.0  && xEdgeDirection.y > 0.0  && yEdgeDirection.x > 0.0  && yEdgeDirection.y > 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == false) && (selectionEdgesY_.begin()->ascOrientation_ == false)) ||
                    (xEdgeDirection.x < 0.0  && xEdgeDirection.y < 0.0  && yEdgeDirection.x > 0.0  && yEdgeDirection.y < 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == false) && (selectionEdgesY_.begin()->ascOrientation_ == true)));
                front = (!back);
                left = ((xEdgeDirection.x > 0.0  && xEdgeDirection.y < 0.0  && yEdgeDirection.x > 0.0  && yEdgeDirection.y > 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesY_.begin()->ascOrientation_ == true)) ||
                    (xEdgeDirection.x < 0.0  && xEdgeDirection.y < 0.0  && yEdgeDirection.x > 0.0  && yEdgeDirection.y < 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesY_.begin()->ascOrientation_ == true)) ||
                    (xEdgeDirection.x < 0.0  && xEdgeDirection.y < 0.0  && yEdgeDirection.x > 0.0  && yEdgeDirection.y < 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == false) && (selectionEdgesY_.begin()->ascOrientation_ == true)) ||
                    (xEdgeDirection.x > 0.0  && xEdgeDirection.y < 0.0  && yEdgeDirection.x > 0.0  && yEdgeDirection.y > 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == false) && (selectionEdgesY_.begin()->ascOrientation_ == true)));
                right = (!left);
            }
            else {
                if (selectionEdgesX_.size() > 0)
                    xEdgeDirection = selectionEdgesX_.begin()->endVertex_ - selectionEdgesX_.begin()->startVertex_;
                if (selectionEdgesY_.size() > 0)
                    yEdgeDirection = selectionEdgesY_.begin()->endVertex_ - selectionEdgesY_.begin()->startVertex_;
                if (selectionEdgesZ_.size() > 0)
                    zEdgeDirection = selectionEdgesZ_.begin()->endVertex_ - selectionEdgesZ_.begin()->startVertex_;
                base = ((xEdgeDirection.x > 0.0  && xEdgeDirection.y > 0.0  && yEdgeDirection.x > 0.0  && yEdgeDirection.y < 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesY_.begin()->ascOrientation_ == false)) ||
                    (xEdgeDirection.x > 0.0  && xEdgeDirection.y < 0.0  && yEdgeDirection.x > 0.0  && yEdgeDirection.y > 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesY_.begin()->ascOrientation_ == true)) ||
                    (xEdgeDirection.x < 0.0  && xEdgeDirection.y < 0.0  && yEdgeDirection.x > 0.0  && yEdgeDirection.y < 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesY_.begin()->ascOrientation_ == true)) ||
                    (xEdgeDirection.x < 0.0 && xEdgeDirection.y > 0.0  && yEdgeDirection.x > 0.0  && yEdgeDirection.y > 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesY_.begin()->ascOrientation_ == false)) ||
                    (xEdgeDirection.x > 0.0 && xEdgeDirection.y < 0.0  && yEdgeDirection.x < 0.0  && yEdgeDirection.y < 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesY_.begin()->ascOrientation_ == false)) ||
                    (xEdgeDirection.x > 0.0 && xEdgeDirection.y < 0.0  && zEdgeDirection.x < 0.0  && zEdgeDirection.y < 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesZ_.begin()->ascOrientation_ == false)) ||
                    (xEdgeDirection.x > 0.0 && xEdgeDirection.y > 0.0  && zEdgeDirection.x > 0.0  && zEdgeDirection.y < 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesZ_.begin()->ascOrientation_ == false)) ||
                    (xEdgeDirection.x < 0.0 && xEdgeDirection.y < 0.0  && zEdgeDirection.x > 0.0  && zEdgeDirection.y > 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesZ_.begin()->ascOrientation_ == true)) ||
                    (xEdgeDirection.x < 0.0 && xEdgeDirection.y < 0.0  && zEdgeDirection.x > 0.0  && zEdgeDirection.y < 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == false) && (selectionEdgesZ_.begin()->ascOrientation_ == true)) ||
                    (xEdgeDirection.x < 0.0 && xEdgeDirection.y > 0.0  && zEdgeDirection.x > 0.0  && zEdgeDirection.y > 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesZ_.begin()->ascOrientation_ == false)) ||
                    (xEdgeDirection.x > 0.0 && xEdgeDirection.y < 0.0  && zEdgeDirection.x > 0.0  && zEdgeDirection.y < 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == false) && (selectionEdgesZ_.begin()->ascOrientation_ == true)) ||
                    (xEdgeDirection.x < 0.0 && xEdgeDirection.y < 0.0  && zEdgeDirection.x < 0.0  && zEdgeDirection.y > 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesZ_.begin()->ascOrientation_ == false)) ||
                    (yEdgeDirection.x < 0.0 && yEdgeDirection.y < 0.0  && zEdgeDirection.x > 0.0  && zEdgeDirection.y < 0.0 && (selectionEdgesY_.begin()->ascOrientation_ == false) && (selectionEdgesZ_.begin()->ascOrientation_ == true)) ||
                    (yEdgeDirection.x < 0.0 && yEdgeDirection.y > 0.0  && zEdgeDirection.x < 0.0  && zEdgeDirection.y > 0.0 && (selectionEdgesY_.begin()->ascOrientation_ == true) && (selectionEdgesZ_.begin()->ascOrientation_ == false)) ||
                    (yEdgeDirection.x > 0.0 && yEdgeDirection.y > 0.0  && zEdgeDirection.x < 0.0  && zEdgeDirection.y > 0.0 && (selectionEdgesY_.begin()->ascOrientation_ == false) && (selectionEdgesZ_.begin()->ascOrientation_ == true)));
                top = ((xEdgeDirection.x > 0.0 && xEdgeDirection.y > 0.0  && zEdgeDirection.x > 0.0  && zEdgeDirection.y < 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesZ_.begin()->ascOrientation_ == false) && selectionEdgesY_.size() == 0) ||
                    (xEdgeDirection.x > 0.0 && xEdgeDirection.y < 0.0  && zEdgeDirection.x < 0.0  && zEdgeDirection.y < 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesZ_.begin()->ascOrientation_ == false)) ||
                    (xEdgeDirection.x < 0.0 && xEdgeDirection.y < 0.0  && zEdgeDirection.x < 0.0  && zEdgeDirection.y > 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesZ_.begin()->ascOrientation_ == false)) ||
                    (xEdgeDirection.x < 0.0 && xEdgeDirection.y > 0.0  && zEdgeDirection.x > 0.0  && zEdgeDirection.y > 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesZ_.begin()->ascOrientation_ == false)) ||
                    (xEdgeDirection.x > 0.0  && xEdgeDirection.y < 0.0  && yEdgeDirection.x > 0.0  && yEdgeDirection.y > 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesY_.begin()->ascOrientation_ == true)) ||
                    (xEdgeDirection.x < 0.0 && xEdgeDirection.y < 0.0  && zEdgeDirection.x > 0.0  && zEdgeDirection.y < 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == false) && (selectionEdgesZ_.begin()->ascOrientation_ == true)) ||
                    (yEdgeDirection.x < 0.0 && yEdgeDirection.y < 0.0  && zEdgeDirection.x > 0.0  && zEdgeDirection.y < 0.0 && (selectionEdgesY_.begin()->ascOrientation_ == false) && (selectionEdgesZ_.begin()->ascOrientation_ == true)));
                back = ((xEdgeDirection.x > 0.0  && xEdgeDirection.y < 0.0  && yEdgeDirection.x > 0.0  && yEdgeDirection.y > 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesY_.begin()->ascOrientation_ == true)) ||
                    (xEdgeDirection.x > 0.0  && xEdgeDirection.y > 0.0  && yEdgeDirection.x > 0.0  && yEdgeDirection.y < 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesY_.begin()->ascOrientation_ == false)) ||
                    (xEdgeDirection.x < 0.0  && xEdgeDirection.y > 0.0  && yEdgeDirection.x > 0.0  && yEdgeDirection.y > 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == false) && (selectionEdgesY_.begin()->ascOrientation_ == false)) ||
                    (xEdgeDirection.x < 0.0  && xEdgeDirection.y < 0.0  && yEdgeDirection.x > 0.0  && yEdgeDirection.y < 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == false) && (selectionEdgesY_.begin()->ascOrientation_ == true)) ||
                    (xEdgeDirection.x > 0.0 && xEdgeDirection.y > 0.0  && zEdgeDirection.x > 0.0  && zEdgeDirection.y < 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesZ_.begin()->ascOrientation_ == false)) ||
                    (xEdgeDirection.x > 0.0 && xEdgeDirection.y < 0.0  && zEdgeDirection.x < 0.0  && zEdgeDirection.y < 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesZ_.begin()->ascOrientation_ == false)) ||
                    (xEdgeDirection.x < 0.0 && xEdgeDirection.y < 0.0  && zEdgeDirection.x < 0.0  && zEdgeDirection.y > 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesZ_.begin()->ascOrientation_ == false)) ||
                    (xEdgeDirection.x < 0.0 && xEdgeDirection.y > 0.0  && zEdgeDirection.x > 0.0  && zEdgeDirection.y > 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesZ_.begin()->ascOrientation_ == false)) ||
                    (yEdgeDirection.x < 0.0 && yEdgeDirection.y < 0.0  && zEdgeDirection.x > 0.0  && zEdgeDirection.y < 0.0 && (selectionEdgesY_.begin()->ascOrientation_ == false) && (selectionEdgesZ_.begin()->ascOrientation_ == true)));
                front = ((xEdgeDirection.x < 0.0 && xEdgeDirection.y < 0.0  && zEdgeDirection.x > 0.0  && zEdgeDirection.y < 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == false) && (selectionEdgesZ_.begin()->ascOrientation_ == true)) ||
                    (xEdgeDirection.x > 0.0 && xEdgeDirection.y < 0.0  && zEdgeDirection.x > 0.0  && zEdgeDirection.y > 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == false) && (selectionEdgesZ_.begin()->ascOrientation_ == true)) ||
                    (yEdgeDirection.x < 0.0 && yEdgeDirection.y < 0.0  && zEdgeDirection.x > 0.0  && zEdgeDirection.y < 0.0 && (selectionEdgesY_.begin()->ascOrientation_ == false) && (selectionEdgesZ_.begin()->ascOrientation_ == true)));
                left = ((xEdgeDirection.x > 0.0  && xEdgeDirection.y < 0.0  && yEdgeDirection.x > 0.0  && yEdgeDirection.y > 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesY_.begin()->ascOrientation_ == true)) ||
                    (xEdgeDirection.x < 0.0  && xEdgeDirection.y < 0.0  && yEdgeDirection.x > 0.0  && yEdgeDirection.y < 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesY_.begin()->ascOrientation_ == true)) ||
                    (xEdgeDirection.x < 0.0  && xEdgeDirection.y < 0.0  && yEdgeDirection.x > 0.0  && yEdgeDirection.y < 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == false) && (selectionEdgesY_.begin()->ascOrientation_ == true)) ||
                    (xEdgeDirection.x > 0.0  && xEdgeDirection.y < 0.0  && yEdgeDirection.x > 0.0  && yEdgeDirection.y > 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == false) && (selectionEdgesY_.begin()->ascOrientation_ == true)) ||
                    (xEdgeDirection.x > 0.0 && xEdgeDirection.y < 0.0  && zEdgeDirection.x > 0.0  && zEdgeDirection.y > 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == false) && (selectionEdgesZ_.begin()->ascOrientation_ == true)) ||
                    (xEdgeDirection.x > 0.0 && xEdgeDirection.y > 0.0  && zEdgeDirection.x > 0.0  && zEdgeDirection.y < 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesZ_.begin()->ascOrientation_ == false)) ||
                    (xEdgeDirection.x > 0.0 && xEdgeDirection.y < 0.0  && zEdgeDirection.x < 0.0  && zEdgeDirection.y < 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesZ_.begin()->ascOrientation_ == false)) ||
                    (xEdgeDirection.x < 0.0 && xEdgeDirection.y < 0.0  && zEdgeDirection.x < 0.0  && zEdgeDirection.y > 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesZ_.begin()->ascOrientation_ == false)) ||
                    (xEdgeDirection.x < 0.0 && xEdgeDirection.y > 0.0  && zEdgeDirection.x > 0.0  && zEdgeDirection.y > 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesZ_.begin()->ascOrientation_ == false)) ||
                    (xEdgeDirection.x < 0.0 && xEdgeDirection.y < 0.0  && zEdgeDirection.x > 0.0  && zEdgeDirection.y < 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == false) && (selectionEdgesZ_.begin()->ascOrientation_ == true)) ||
                    (yEdgeDirection.x < 0.0 && yEdgeDirection.y < 0.0  && zEdgeDirection.x > 0.0  && zEdgeDirection.y < 0.0 && (selectionEdgesY_.begin()->ascOrientation_ == false) && (selectionEdgesZ_.begin()->ascOrientation_ == true)));
                right = ((xEdgeDirection.x > 0.0 && xEdgeDirection.y > 0.0  && zEdgeDirection.x > 0.0  && zEdgeDirection.y < 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesZ_.begin()->ascOrientation_ == false)) ||
                    (xEdgeDirection.x > 0.0 && xEdgeDirection.y < 0.0  && zEdgeDirection.x < 0.0  && zEdgeDirection.y < 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesZ_.begin()->ascOrientation_ == false) && selectionEdgesY_.size() == 0) ||
                    (xEdgeDirection.x < 0.0 && xEdgeDirection.y < 0.0  && zEdgeDirection.x < 0.0  && zEdgeDirection.y > 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesZ_.begin()->ascOrientation_ == false)) ||
                    (xEdgeDirection.x < 0.0 && xEdgeDirection.y > 0.0  && zEdgeDirection.x > 0.0  && zEdgeDirection.y > 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == true) && (selectionEdgesZ_.begin()->ascOrientation_ == false)) ||
                    (xEdgeDirection.x < 0.0 && xEdgeDirection.y < 0.0  && zEdgeDirection.x > 0.0  && zEdgeDirection.y < 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == false) && (selectionEdgesZ_.begin()->ascOrientation_ == true)) ||
                    (xEdgeDirection.x > 0.0 && xEdgeDirection.y < 0.0  && zEdgeDirection.x > 0.0  && zEdgeDirection.y > 0.0 && (selectionEdgesX_.begin()->ascOrientation_ == false) && (selectionEdgesZ_.begin()->ascOrientation_ == true)));
            }
            if (axis == X_AXIS) {
                for (plot_t i = start + step.x ; i <= domain_[axis].getRight(); i += step.x) {
                    if (base) {
                        p1 = projection_neu(i - step.x, yl, zl);
                        p2 = projection_neu(i, yl, zl);
                        p3 = projection_neu(i, yr, zl);
                        p4 = projection_neu(i - step.x, yr, zl);
                        svgLine(p1,p2,svgColor_,0,lineWidth_);
                        svgLine(p2,p3,svgColor_,0,lineWidth_);
                        svgLine(p3,p4,svgColor_,0,lineWidth_);
                        svgLine(p4,p1,svgColor_,0,lineWidth_);
                    }
                    if (top) {
                        p1 = projection_neu(i, yl, zr);
                        p2 = projection_neu(i - step.x, yl, zr);
                        p3 = projection_neu(i - step.x, yr, zr);
                        p4 = projection_neu(i, yr, zr);
                        svgLine(p1,p2,svgColor_,0,lineWidth_);
                        svgLine(p2,p3,svgColor_,0,lineWidth_);
                        svgLine(p3,p4,svgColor_,0,lineWidth_);
                        svgLine(p4,p1,svgColor_,0,lineWidth_);
                    }
                    if (front) {
                        p1 = projection_neu(i - step.x, yl, zl);
                        p2 = projection_neu(i - step.x, yl, zr);
                        p3 = projection_neu(i, yl, zr);
                        p4 = projection_neu(i, yl, zl);
                        svgLine(p1,p2,svgColor_,0,lineWidth_);
                        svgLine(p2,p3,svgColor_,0,lineWidth_);
                        svgLine(p3,p4,svgColor_,0,lineWidth_);
                        svgLine(p4,p1,svgColor_,0,lineWidth_);
                    }
                    if (back) {
                        p1 = projection_neu(i, yr, zl);
                        p2 = projection_neu(i, yr, zr);
                        p3 = projection_neu(i - step.x, yr, zr);
                        p4 = projection_neu(i - step.x, yr, zl);
                        svgLine(p1,p2,svgColor_,0,lineWidth_);
                        svgLine(p2,p3,svgColor_,0,lineWidth_);
                        svgLine(p3,p4,svgColor_,0,lineWidth_);
                        svgLine(p4,p1,svgColor_,0,lineWidth_);
                    }
                }
            }
            else if (axis == Y_AXIS) {
                for (plot_t i = start + step.x ; i <= domain_[axis].getRight(); i += step.x) {
                    if (base) {
                        p1 = projection_neu(xr, i - step.x, zl);
                        p2 = projection_neu(xr, i, zl);
                        p3 = projection_neu(xl, i, zl);
                        p4 = projection_neu(xl, i - step.x, zl);
                        svgLine(p1,p2,svgColor_,0,lineWidth_);
                        svgLine(p2,p3,svgColor_,0,lineWidth_);
                        svgLine(p3,p4,svgColor_,0,lineWidth_);
                        svgLine(p4,p1,svgColor_,0,lineWidth_);
                    }
                    if (top) {
                        p1 = projection_neu(xr, i, zr);
                        p2 = projection_neu(xr, i - step.x, zr);
                        p3 = projection_neu(xl, i - step.x, zr);
                        p4 = projection_neu(xl, i, zr);
                        svgLine(p1,p2,svgColor_,0,lineWidth_);
                        svgLine(p2,p3,svgColor_,0,lineWidth_);
                        svgLine(p3,p4,svgColor_,0,lineWidth_);
                        svgLine(p4,p1,svgColor_,0,lineWidth_);
                    }
                    if (left) {
                        p1 = projection_neu(xl, i - step.x, zr);
                        p2 = projection_neu(xl, i - step.x, zl);
                        p3 = projection_neu(xl, i, zl);
                        p4 = projection_neu(xl, i, zr);
                        svgLine(p1,p2,svgColor_,0,lineWidth_);
                        svgLine(p2,p3,svgColor_,0,lineWidth_);
                        svgLine(p3,p4,svgColor_,0,lineWidth_);
                        svgLine(p4,p1,svgColor_,0,lineWidth_);
                    }
                    if (right) {
                        p1 = projection_neu(xr, i, zr);
                        p2 = projection_neu(xr, i,  zl);
                        p3 = projection_neu(xr ,i - step.x, zl);
                        p4 = projection_neu(xr, i - step.x, zr);
                        svgLine(p1,p2,svgColor_,0,lineWidth_);
                        svgLine(p2,p3,svgColor_,0,lineWidth_);
                        svgLine(p3,p4,svgColor_,0,lineWidth_);
                        svgLine(p4,p1,svgColor_,0,lineWidth_);
                    }
                }
            }
            else if (axis == Z_AXIS) {
                for (plot_t i = start + step.x ; i <= domain_[axis].getRight(); i += step.x) {
                    if (back) {
                        p1 = projection_neu(xr, yr, i - step.x);
                        p2 = projection_neu(xr, yr, i);
                        p3 = projection_neu(xl, yr, i);
                        p4 = projection_neu(xl, yr, i - step.x);
                        svgLine(p1,p2,svgColor_,0,lineWidth_);
                        svgLine(p2,p3,svgColor_,0,lineWidth_);
                        svgLine(p3,p4,svgColor_,0,lineWidth_);
                        svgLine(p4,p1,svgColor_,0,lineWidth_);
                    }
                    if (front) {
                        p1 = projection_neu(xr, yl, i);
                        p2 = projection_neu(xr, yl, i - step.x);
                        p3 = projection_neu(xl, yl, i - step.x);
                        p4 = projection_neu(xl, yl, i);
                        svgLine(p1,p2,svgColor_,0,lineWidth_);
                        svgLine(p2,p3,svgColor_,0,lineWidth_);
                        svgLine(p3,p4,svgColor_,0,lineWidth_);
                        svgLine(p4,p1,svgColor_,0,lineWidth_);
                    }
                    if (left) {
                        p1 = projection_neu(xl, yl, i - step.x);
                        p2 = projection_neu(xl, yr, i - step.x);
                        p3 = projection_neu(xl, yr, i);
                        p4 = projection_neu(xl, yl, i);
                        svgLine(p1,p2,svgColor_,0,lineWidth_);
                        svgLine(p2,p3,svgColor_,0,lineWidth_);
                        svgLine(p3,p4,svgColor_,0,lineWidth_);
                        svgLine(p4,p1,svgColor_,0,lineWidth_);
                    }
                    if (right) {
                        p1 = projection_neu(xr, yl, i);
                        p2 = projection_neu(xr, yr, i);
                        p3 = projection_neu(xr, yr, i - step.x);
                        p4 = projection_neu(xr, yl,i - step.x);
                        svgLine(p1,p2,svgColor_,0,lineWidth_);
                        svgLine(p2,p3,svgColor_,0,lineWidth_);
                        svgLine(p3,p4,svgColor_,0,lineWidth_);
                        svgLine(p4,p1,svgColor_,0,lineWidth_);
                    }
                }
            }
        }
    }
}

void PlotLibrarySvg::renderAxisLabelScales(const PlotData& data, int indexLabel, bool helperLines) {
    std::string label;
    xAxisLabelGroup_.reset();
    xAxisLabelGroup_.setBounds(getBoundsBelowPlot());
    std::vector<PlotRowValue>::const_iterator it = data.getRowsBegin();
    plot_t x = 0;
    std::vector< PlotLibraryFileBase::Projection_Coordinates > points;
    for (;it!=data.getRowsEnd();++it) {
        if (!domain_[X_AXIS].contains(x) || it->getCellAt(indexLabel).isNull()){
            x += 1;
            continue;
        }
        if (data.getColumnType(indexLabel) == PlotBase::STRING)
            label = it->getTagAt(indexLabel);
        else {
            std::ostringstream stream;
            stream << it->getValueAt(indexLabel);
            label = stream.str();
        }
        xAxisLabelGroup_.addLabel(label,
                                  convertPlotCoordinatesToViewport3(tgt::dvec3(x, domain_[1].getLeft(), domain_[2].getRight())),
                                  fontColor_, fontSize_, SmartLabel::TOPCENTERED);
        if (helperLines && dimension_ == TWO) {
            points.clear();
            points.push_back(convertPlotCoordinatesToViewport3Projection(x,domain_[1].getLeft(),0));
            points.push_back(convertPlotCoordinatesToViewport3Projection(x,domain_[1].getRight(),0));
            svgLine(points.at(0),points.at(1),drawingColor_,0,axesWidth_/2.f);
        }
        x += 1;
    }
    renderSmartLabelGroup(&xAxisLabelGroup_);
}

void PlotLibrarySvg::renderLabel(tgt::vec3 pos, const SmartLabel::Alignment align, const std::string& text,
                              bool viewCoordinates, int padding) {
    if (!viewCoordinates)
        pos = convertPlotCoordinatesToViewport3(pos);

    tgt::Matrix4<plot_t> oldpjm = projectionsMatrix_;
    tgt::Matrix4<plot_t> oldmdlv = modelviewMatrix_;
    projectionsMatrix_ = plOrtho(0.0, windowSize_.x, 0.0, windowSize_.y, -1.0, 1.0);
    modelviewMatrix_ = tgt::Matrix4<plot_t>::createIdentity();

    tgt::dvec3 position = pos;
    SmartLabelGroupBaseSvg::renderSingleLabel(&labelFont_,  text,fontSize_, align, static_cast<double>(padding), position);
    svgText(position,text,&labelFont_,fontSize_,fontColor_,align);

    projectionsMatrix_ = oldpjm;
    modelviewMatrix_ = oldmdlv;
}

void PlotLibrarySvg::renderLabel(tgt::dvec2 pos, const SmartLabel::Alignment align, const std::string& text, int padding) {
    tgt::Matrix4<plot_t> oldpjm = projectionsMatrix_;
    tgt::Matrix4<plot_t> oldmdlv = modelviewMatrix_;
    projectionsMatrix_ = plOrtho(0.0, windowSize_.x, 0.0, windowSize_.y, -1.0, 1.0);
    modelviewMatrix_ = tgt::Matrix4<plot_t>::createIdentity();

    tgt::dvec3 position = tgt::dvec3(pos,0.0);
    SmartLabelGroupBaseSvg::renderSingleLabel(&labelFont_, text, fontSize_, align,static_cast<double>(padding), position);
    svgText(position,text,&labelFont_,fontSize_,fontColor_,align);

    projectionsMatrix_ = oldpjm;
    modelviewMatrix_ = oldmdlv;
}

void PlotLibrarySvg::addPlotLabel(std::string text, tgt::vec3 position, tgt::Color color,
                               int size, SmartLabel::Alignment align) {
    plotLabelGroup_.addLabel(text, position, color, size, align);
}

void PlotLibrarySvg::addLineLabel(std::string text, tgt::vec3 position, tgt::Color color,
                               int size, SmartLabel::Alignment align) {
    lineLabelGroup_.addLabel(text, position, color, size, align);
}

//
// helper functions
//

void PlotLibrarySvg::resetLineLabels() {
    lineLabelGroup_.reset();
    lineLabelGroup_.setBounds(getBoundsRightOfPlot());
}

void PlotLibrarySvg::renderLineLabels() {
    renderSmartLabelGroup(&lineLabelGroup_);
}

void PlotLibrarySvg::resetPlotLabels() {
    plotLabelGroup_.reset();
    plotLabelGroup_.setBounds(getBoundsPlot());
}

void PlotLibrarySvg::renderPlotLabels() {
    renderSmartLabelGroup(&plotLabelGroup_);
}

void PlotLibrarySvg::renderSmartLabelGroup(SmartLabelGroupBase* smg) {
    smg->performLayout();

    tgt::Matrix4<plot_t> oldpjm = projectionsMatrix_;
    tgt::Matrix4<plot_t> oldmdlv = modelviewMatrix_;
    projectionsMatrix_ = plOrtho(0.0, windowSize_.x, 0.0, windowSize_.y, -1.0, 1.0);
    modelviewMatrix_ = tgt::Matrix4<plot_t>::createIdentity();

    smg->render();
    if (dynamic_cast<SmartLabelGroupSvgNoLayoutWithBackground*>(smg) != 0) {
        SmartLabelGroupSvgNoLayoutWithBackground* slglabel = static_cast<SmartLabelGroupSvgNoLayoutWithBackground*>(smg);
        std::vector< std::vector< tgt::Vector3<double> > > polygonpoints = slglabel->getPolygonPoints();
        std::vector<SmartLabel> slabelvector = slglabel->getLabels();
        std::vector< tgt::Vector3<double> > points;
        std::vector< PlotLibraryFileBase::Projection_Coordinates > projPoints;
        points.resize(4);
        projPoints.resize(4);
        for (size_t i = 0; i < polygonpoints.size(); ++i) {
            for (int k = 0; k < 4; ++k) {
                points[k] = polygonpoints.at(i).at(k);
                projPoints[k] = PlotLibraryFileBase::Projection_Coordinates(points[k],points[k],points[k]);
            }
            svgPolygon(projPoints,slglabel->getBackgroundColor(),lineColor_,-1);
            svgTextAreaSVG11(tgt::Vector3<double>(points[0].x,points[0].y+(points[2].y-points[0].y),points[0].z),points[1].x-points[0].x,
                points[2].y-points[0].y,slabelvector[i].text_,
                slglabel->getFont(),slabelvector[i].size_,slabelvector[i].color_);
        }
    }
    else {
        SmartLabelGroupBaseSvg* slgb = dynamic_cast<SmartLabelGroupBaseSvg*>(smg);
        if (slgb != 0) {
            std::vector<SmartLabel> slabelvector = slgb->getLabels();
            tgt::Font* font = slgb->getFont();
            tgt::dvec3 pos;
            for (size_t i = 0; i < slabelvector.size(); ++i) {
                SmartLabel label = slabelvector.at(i);
                pos = label.position_;
                svgText(pos,label.text_,font,label.size_,label.color_,label.align_);
            }
        }
    }

    projectionsMatrix_ = oldpjm;
    modelviewMatrix_ = oldmdlv;
}

void PlotLibrarySvg::resetRenderStatus() {
    lineWidth_ = 1;
}

void PlotLibrarySvg::renderGlyph(plot_t x, plot_t y, plot_t z, plot_t size) {
    if (usePlotPickingManager_)
        return;
    PlotLibraryFileBase::Projection_Coordinates point = projection_neu(x,y,z);
    tgt::Vector3<plot_t> hpoint = point.afterProjection_;
    std::vector< PlotLibraryFileBase::Projection_Coordinates > points;
    if (glyphStyle_ == PlotEntitySettings::POINT) {
        svgCircle(point,size/2,svgColor_,svgColor_,0);
        return;
    }
    if (dimension_ == TWO){
        if (glyphStyle_ == PlotEntitySettings::CIRCLE){
            if (texture_ == 0) {
                svgCircle(point,size/2,svgColor_,svgColor_,0);
            }
            else
                svgCircle(point,size/2,svgColor_,svgColor_,0,0,PlotEntitySettings::TEXTURE);
        }
        else if (glyphStyle_ == PlotEntitySettings::TRIANGLE){
            points.clear();
            tgt::Vector3<plot_t> coor = tgt::Vector3<plot_t>(hpoint.x,hpoint.y + size/2,0);
            points.push_back(PlotLibraryFileBase::Projection_Coordinates(coor,coor,coor));
            coor = tgt::Vector3<plot_t>(hpoint.x - size/3,hpoint.y - size/3,0);
            points.push_back(PlotLibraryFileBase::Projection_Coordinates(coor,coor,coor));
            coor = tgt::Vector3<plot_t>(hpoint.x + size/3,hpoint.y - size/3,0);
            points.push_back(PlotLibraryFileBase::Projection_Coordinates(coor,coor,coor));
            if (texture_ == 0) {
                svgPolygon(points,svgColor_,svgColor_,0);
            }
            else
                svgPolygon(points,svgColor_,svgColor_,0,0,PlotEntitySettings::TEXTURE);
        }
        else if (glyphStyle_ == PlotEntitySettings::QUAD){
            points.clear();
            tgt::Vector3<plot_t> coor = tgt::Vector3<plot_t>(hpoint.x-size/2,hpoint.y - size/2,0);
            points.push_back(PlotLibraryFileBase::Projection_Coordinates(coor,coor,coor));
            coor = tgt::Vector3<plot_t>(hpoint.x+size/2,hpoint.y - size/2,0);
            points.push_back(PlotLibraryFileBase::Projection_Coordinates(coor,coor,coor));
            coor = tgt::Vector3<plot_t>(hpoint.x+size/2,hpoint.y + size/2,0);
            points.push_back(PlotLibraryFileBase::Projection_Coordinates(coor,coor,coor));
            coor = tgt::Vector3<plot_t>(hpoint.x+size/2,hpoint.y + size/2,0);
            points.push_back(PlotLibraryFileBase::Projection_Coordinates(coor,coor,coor));
            if (texture_ == 0) {
                svgPolygon(points,svgColor_,svgColor_,0);
            }
            else
                svgPolygon(points,svgColor_,svgColor_,0,0,PlotEntitySettings::TEXTURE);
        }
    }
    else if (dimension_ == THREE){
        if (glyphStyle_ == PlotEntitySettings::CIRCLE){
            if (texture_ == 0) {
                svgSphere(tgt::Vector3<plot_t>(x,y,z),size,svgColor_,tgt::Color(1,1,1,0.01f),20,0,1);
            }
            else
                svgSphere(tgt::Vector3<plot_t>(x,y,z),size,svgColor_,tgt::Color(1,1,1,0.01f),20,0,1,PlotEntitySettings::TEXTURE);
        }
        else if (glyphStyle_ == PlotEntitySettings::TRIANGLE){
            if (texture_ == 0) {
                svgTriangle3D(tgt::Vector3<plot_t>(x,y,z),size,size,svgColor_,tgt::Color(1,1,1,0.01f),0,1);
            }
            else
                svgTriangle3D(tgt::Vector3<plot_t>(x,y,z),size,size,svgColor_,tgt::Color(1,1,1,0.01f),0,1,PlotEntitySettings::TEXTURE);
        }
        else if (glyphStyle_ == PlotEntitySettings::QUAD){
            if (texture_ == 0) {
                svgQuad(tgt::Vector3<plot_t>(x,y,z),size,size,size,svgColor_,tgt::Color(1,1,1,0.01f),0,1);
            }
            else
                svgQuad(tgt::Vector3<plot_t>(x,y,z),size,size,size,svgColor_,tgt::Color(1,1,1,0.01f),0,1,PlotEntitySettings::TEXTURE);
        }
    }
}

void PlotLibrarySvg::renderSingleBar(plot_t left, plot_t right, plot_t bottom, plot_t top, tgt::Color c, plot_t squeeze) {
    if (usePlotPickingManager_)
        return;
    PlotLibraryFileBase::Projection_Coordinates p1;
    PlotLibraryFileBase::Projection_Coordinates p2;
    PlotLibraryFileBase::Projection_Coordinates p3;
    PlotLibraryFileBase::Projection_Coordinates p4;

    plot_t back = domain_[Z_AXIS].getLeft();
    plot_t front = domain_[Z_AXIS].getRight();
    //squeeze
    left = (left+squeeze*left+right-squeeze*right)/2.0;
    right = (left-squeeze*left+right+squeeze*right)/2.0;
    back = (back+squeeze*back+front-squeeze*front)/2.0;
    front = (back-squeeze*back+front+squeeze*front)/2.0;
    tgt::Color color;
    std::vector< PlotLibraryFileBase::Projection_Coordinates > points;
    points.resize(4);
    lineColor_ = tgt::Color(0.2f*c.r, 0.2f*c.g, 0.2f*c.b, 1.0f);

    //bottom
    color = tgt::Color(0.9f*c.r, 0.9f*c.g, 0.9f*c.b, 1.0f);
    p1 = convertPlotCoordinatesToViewport3Projection(left,bottom,front);
    p2 = convertPlotCoordinatesToViewport3Projection(left,bottom,back);
    p3 = convertPlotCoordinatesToViewport3Projection(right,bottom,back);
    p4 = convertPlotCoordinatesToViewport3Projection(right,bottom,front);
    points[0] = p1;
    points[1] = p2;
    points[2] = p3;
    points[3] = p4;
    svgPolygon(points,color,lineColor_,0,lineWidth_);

    //back
    color = tgt::Color(0.8f*c.r, 0.8f*c.g, 0.8f*c.b, 1.0f);
    p1 = convertPlotCoordinatesToViewport3Projection(right,bottom,back);
    p2 = convertPlotCoordinatesToViewport3Projection(right,top,back);
    p3 = convertPlotCoordinatesToViewport3Projection(left,top,back);
    p4 = convertPlotCoordinatesToViewport3Projection(left,bottom,back);
    points[0] = p1;
    points[1] = p2;
    points[2] = p3;
    points[3] = p4;
    svgPolygon(points,color,lineColor_,0,lineWidth_);

    //left
    color = tgt::Color(c.r, c.g, c.b, 1.0f);
    p1 = convertPlotCoordinatesToViewport3Projection(left,bottom,front);
    p2 = convertPlotCoordinatesToViewport3Projection(left,top,front);
    p3 = convertPlotCoordinatesToViewport3Projection(left,top,back);
    p4 = convertPlotCoordinatesToViewport3Projection(left,bottom,back);
    points[0] = p1;
    points[1] = p2;
    points[2] = p3;
    points[3] = p4;
    svgPolygon(points,color,lineColor_,0,lineWidth_);

    //right
    color = tgt::Color(0.8f*c.r, 0.8f*c.g, 0.8f*c.b, 1.0f);
    p1 = convertPlotCoordinatesToViewport3Projection(right,bottom,front);
    p2 = convertPlotCoordinatesToViewport3Projection(right,top,front);
    p3 = convertPlotCoordinatesToViewport3Projection(right,top,back);
    p4 = convertPlotCoordinatesToViewport3Projection(right,bottom,back);
    points[0] = p1;
    points[1] = p2;
    points[2] = p3;
    points[3] = p4;
    svgPolygon(points,color,lineColor_,0,lineWidth_);

    //top
    color = tgt::Color(0.9f*c.r, 0.9f*c.g, 0.9f*c.b, 1.0f);
    p1 = convertPlotCoordinatesToViewport3Projection(left,top,front);
    p2 = convertPlotCoordinatesToViewport3Projection(left,top,back);
    p3 = convertPlotCoordinatesToViewport3Projection(right,top,back);
    p4 = convertPlotCoordinatesToViewport3Projection(right,top,front);
    points[0] = p1;
    points[1] = p2;
    points[2] = p3;
    points[3] = p4;
    svgPolygon(points,color,lineColor_,0,lineWidth_);

    //front
    color = tgt::Color(c.r, c.g, c.b, 1.0f);
    p1 = convertPlotCoordinatesToViewport3Projection(left,bottom,front);
    p2 = convertPlotCoordinatesToViewport3Projection(left,top,front);
    p3 = convertPlotCoordinatesToViewport3Projection(right,top,front);
    p4 = convertPlotCoordinatesToViewport3Projection(right,bottom,front);
    points[0] = p1;
    points[1] = p2;
    points[2] = p3;
    points[3] = p4;
    svgPolygon(points,color,lineColor_,0,lineWidth_);
}

void PlotLibrarySvg::svgCircle(const PlotLibraryFileBase::Projection_Coordinates& point, plot_t radius, tgt::Color fillColor, tgt::Color borderColor, int clipping_number, double lineWidth, PlotEntitySettings::PolygonFillStyle fillStyle, PlotEntitySettings::LineStyle lineStyle) {
    if (clipping_number == 0 || inClippRegion(point.afterModelview_))
        renderVector_.push_back(new PlotLibrarySvgRenderCircle(point,radius,fillColor,borderColor,clipping_number,lineWidth,textureNumber_,fillStyle,lineStyle));
}

void PlotLibrarySvg::svgLine(const PlotLibraryFileBase::Projection_Coordinates& startPoint, const PlotLibraryFileBase::Projection_Coordinates& endPoint, tgt::Color lineColor, int clipping_number, double lineWidth, PlotEntitySettings::LineStyle lineStyle) {
    if (clipping_number == 0 || inClippRegion(startPoint.afterModelview_) || inClippRegion(endPoint.afterModelview_))
        renderVector_.push_back(new PlotLibrarySvgRenderLine(startPoint,endPoint,lineColor,clipping_number,lineWidth,lineStyle));
}

void PlotLibrarySvg::svgPolygon(const std::vector< PlotLibraryFileBase::Projection_Coordinates >& points, tgt::Color fillColor, tgt::Color lineColor, int clipping_number, double lineWidth, PlotEntitySettings::PolygonFillStyle fillStyle, PlotEntitySettings::LineStyle lineStyle) {
    if (clipping_number == 0 || inClippRegion(points))
        renderVector_.push_back(new PlotLibrarySvgRenderPolygon(points,fillColor,lineColor,clipping_number,lineWidth,textureNumber_,fillStyle,lineStyle));
}

void PlotLibrarySvg::svgPolyline(const std::vector< PlotLibraryFileBase::Projection_Coordinates >& points, tgt::Color fillColor, tgt::Color lineColor, int clipping_number, double lineWidth, PlotEntitySettings::PolygonFillStyle fillStyle, PlotEntitySettings::LineStyle lineStyle) {
    if (clipping_number == 0 || inClippRegion(points))
        renderVector_.push_back(new PlotLibrarySvgRenderPolyline(points,fillColor,lineColor,clipping_number,lineWidth,textureNumber_,fillStyle,lineStyle));
}

void PlotLibrarySvg::svgEllipse(const PlotLibraryFileBase::Projection_Coordinates& mpoint, double xradius, double yradius, tgt::Color fillColor, tgt::Color lineColor, int clipping_number, double lineWidth, PlotEntitySettings::PolygonFillStyle fillStyle, PlotEntitySettings::LineStyle lineStyle) {
    if (clipping_number == 0 || inClippRegion(mpoint.afterModelview_))
        renderVector_.push_back(new PlotLibrarySvgRenderEllipse(mpoint,xradius,yradius,fillColor,lineColor,clipping_number,lineWidth,textureNumber_,fillStyle,lineStyle));
}

void PlotLibrarySvg::svgGradient(const tgt::Vector2<plot_t>& startPoint, const tgt::Vector2<plot_t>& endPoint, const std::string& id, tgt::Color startColor, tgt::Color endColor) {
    svgOutput_ << "<defs><linearGradient id=\"" << id << "\" x1=\"" << startPoint.x << "\" y1=\"" << startPoint.y << "\" x2=\"" << endPoint.x << "\" y2=\"" << endPoint.y << "\">"
        << "<stop stop-color=";
    svgOutput_ << PlotLibrarySvgRender::svgColor(startColor);
    svgOutput_ << " offset=\"0\" />"
        << "<stop stop-color=";
    svgOutput_ << PlotLibrarySvgRender::svgColor(endColor);
    svgOutput_ << " offset=\"1\" />"
        << "</linearGradient></defs>\n";
}

void PlotLibrarySvg::svgRect(const PlotLibraryFileBase::Projection_Coordinates& startPoint, double stepWidth, double height, const std::string& fillurl) {
    renderVector_.push_back(new PlotLibrarySvgRenderRect(startPoint,stepWidth,height,fillurl));
}

void PlotLibrarySvg::svgText(const tgt::Vector3<plot_t>& point, const std::string& text,  tgt::Font* font, double size, tgt::Color textColor, const SmartLabel::Alignment align) {
    std::stringstream svgText;
    std::stringstream controlledText;
    ++fountNumberCounter_;
    svgText << "@font-face {\n"
       << "  font-family: 'VoreenFontFamily" << fountNumberCounter_ << "';\n"
       << "  src: url('" << font->getFontName() << "')  format('truetype');\n"
       << "}\n";
    svgText << "<g transform=\"translate(0 " << point.y << ") scale(1 -1) translate(0 " << -point.y << ")\"> ";
    svgText << "<text x=\"" << point.x << "\" y=\"" << point.y << "\" font-family=\"VoreenFontFamily" << fountNumberCounter_
        << "\" font-size=\"" << size << "\" font-style=\"normal\" fill=";
    svgText << PlotLibrarySvgRender::svgColor(textColor);
    svgText << " fill-opacity=";
    svgText << PlotLibrarySvgRender::svgTransparancy(textColor);
    switch (align) {
        case SmartLabel::TOPLEFT:
            svgText << " text-anchor=\"start\"";
            break;
        case SmartLabel::TOPCENTERED:
            svgText << " text-anchor=\"middle\"";
            break;
        case SmartLabel::TOPRIGHT:
            svgText << " text-anchor=\"end\"";
            break;
        case SmartLabel::MIDDLELEFT:
            svgText << " text-anchor=\"start\"";
            break;
        case SmartLabel::CENTERED:
            svgText << " text-anchor=\"middle\"";
            break;
        case SmartLabel::MIDDLERIGHT:
            svgText << " text-anchor=\"end\"";
            break;
        case SmartLabel::BOTTOMLEFT:
            svgText << " text-anchor=\"start\"";
            break;
        case SmartLabel::BOTTOMCENTERED:
            svgText << " text-anchor=\"middle\"";
            break;
        case SmartLabel::BOTTOMRIGHT:
            svgText << " text-anchor=\"end\"";
            break;
    }
    for (size_t i = 0; i < text.size(); ++i) {
        switch (text.at(i)) {
            case '&' : controlledText << "&amp;";
                break;
            case '<' : controlledText << "&lt;";
                break;
            case '>' : controlledText << "&gt;";
                break;
            case '\"' : controlledText << "&quot;";
                break;
            case '\'' : controlledText << "&apos;";
                break;
        default:
            controlledText << text.at(i);
        }
    }

    svgText << " >" << controlledText.str() << "</text> </g>\n";
    renderVector_.push_back(new PlotLibrarySvgRenderText(point,svgText.str()));
}

void PlotLibrarySvg::svgTextAreaSVG11(const tgt::Vector3<plot_t>& point, double /*width*/, double height, const std::string& text,  tgt::Font* font, double size, tgt::Color textColor) {
    std::vector<std::string > copytext;
    copytext.resize(1);
    int j =0;
    for (size_t i = 0; i < text.size(); ++i) {
        switch (text.at(i)) {
            case 'y':
            case 'z':
                ++j; copytext.push_back("");
            case '0':
            case '1':
            case '2':
            case '3':
            case '4':
            case '5':
            case '6':
            case '7':
            case '8':
            case '9':
            case '.':
            case ':':
            case ' ':
            case 'x':
                copytext[j] += text.at(i);
                break;
        }
    }
    for (size_t i = 0; i < copytext.size(); ++i) {
        svgText(tgt::Vector3<plot_t>(point.x+3,point.y-(i+0.76)*height/copytext.size(),point.z),copytext[i],font,size,textColor,SmartLabel::BOTTOMLEFT);
    }
}


void PlotLibrarySvg::svgTextAreaSVG12(const tgt::Vector3<plot_t>& point, double width, double height, const std::string& text,  tgt::Font* font, double size, tgt::Color textColor) {
    std::stringstream svgText;
    svgText << "<g transform=\"translate(0 " << point.y << ") scale(1 -1) translate(0 " << -point.y << ")\"> ";
    svgText << "<textArea x=\"" << point.x << "\" y=\"" << point.y << "\" width=\"" << width << "\" height=\"" << height << "\" "
        << "line-increment=\"auto\" font-size=\"" << size << "\" font-family=\"" << font->getFontName()
        << "\" fill=";
    svgText << PlotLibrarySvgRender::svgColor(textColor);
    svgText << " fill-opacity=" << PlotLibrarySvgRender::svgTransparancy(textColor);
    std::string copytext = "";
    for (size_t i = 0; i < text.size(); ++i) {
        switch (text.at(i)) {
            case 'y':
            case 'z':
                copytext += "<tbreak />";
            case '0':
            case '1':
            case '2':
            case '3':
            case '4':
            case '5':
            case '6':
            case '7':
            case '8':
            case '9':
            case '.':
            case ':':
            case ' ':
            case 'x':
                copytext += text.at(i);
                break;
        }
    }
    svgText << " >" << copytext << "</textArea> </g>\n";
    renderVector_.push_back(new PlotLibrarySvgRenderText(point,svgText.str()));
}

void PlotLibrarySvg::svgSphere(const tgt::Vector3<plot_t>& mpoint, double radius, tgt::Color fillColor, tgt::Color lineColor, int sections, int clipping_number, double lineWidth, PlotEntitySettings::PolygonFillStyle fillStyle, PlotEntitySettings::LineStyle lineStyle) {
    PlotLibraryFileBase::Projection_Coordinates point = convertPlotCoordinatesToViewport3Projection(mpoint);
    if (inClippRegion(point.afterModelview_)) {
        double xscale = domain_[X_AXIS].size()/static_cast<double>(2*windowSize_.x);
        double yscale = domain_[Y_AXIS].size()/static_cast<double>(2*windowSize_.y);
        double zscale = domain_[Z_AXIS].size()/1000.0;
        double a = 2*tgt::PI / static_cast<double>(sections);
        double theta, phi;
        std::vector< PlotLibraryFileBase::Projection_Coordinates > points;
        points.resize(4);
        float fakt = 0;
        for (int k = 0; k < sections; ++k) {
            for (int l = 0; l < sections/2.0; ++l) {
                theta = l*a;
                phi = k*a;
                PlotLibraryFileBase::Projection_Coordinates point1 = convertPlotCoordinatesToViewport3Projection(tgt::Vector3<plot_t>(std::sin(theta)*std::cos(phi)*radius*xscale+mpoint.x,std::sin(theta)*std::sin(phi)*radius*yscale+mpoint.y,std::cos(theta)*radius*zscale+mpoint.z));
                PlotLibraryFileBase::Projection_Coordinates point2 = convertPlotCoordinatesToViewport3Projection(tgt::Vector3<plot_t>(std::sin(theta)*std::cos(phi+a)*radius*xscale+mpoint.x,std::sin(theta)*std::sin(phi+a)*radius*yscale+mpoint.y,std::cos(theta)*radius*zscale+mpoint.z));
                PlotLibraryFileBase::Projection_Coordinates point3 = convertPlotCoordinatesToViewport3Projection(tgt::Vector3<plot_t>(std::sin(theta+a)*std::cos(phi+a)*radius*xscale+mpoint.x,std::sin(theta+a)*std::sin(phi+a)*radius*yscale+mpoint.y,std::cos(theta+a)*radius*zscale+mpoint.z));
                PlotLibraryFileBase::Projection_Coordinates point4 = convertPlotCoordinatesToViewport3Projection(tgt::Vector3<plot_t>(std::sin(theta+a)*std::cos(phi)*radius*xscale+mpoint.x,std::sin(theta+a)*std::sin(phi)*radius*yscale+mpoint.y,std::cos(theta+a)*radius*zscale+mpoint.z));
                points[0] = point1;
                points[1] = point2;
                points[2] = point3;
                points[3] = point4;
                fakt = static_cast<float>(0.4 + (0.6*std::cos(theta)+ 1)/2.0);
                svgPolygon(points,tgt::Color(fillColor.r*fakt,fillColor.g*fakt,fillColor.b*fakt,fillColor.a),lineColor,clipping_number,lineWidth,fillStyle,lineStyle);
            }
        }
    }
}

void PlotLibrarySvg::svgQuad(const tgt::Vector3<plot_t>& mpoint, double xsize, double ysize, double zsize, tgt::Color fillColor, tgt::Color lineColor, int clipping_number, double lineWidth,PlotEntitySettings::PolygonFillStyle fillStyle, PlotEntitySettings::LineStyle lineStyle) {
    PlotLibraryFileBase::Projection_Coordinates point = convertPlotCoordinatesToViewport3Projection(mpoint);
    if (inClippRegion(point.afterModelview_)) {
        double xscale = domain_[X_AXIS].size()/static_cast<double>(2*windowSize_.x);
        double yscale = domain_[Y_AXIS].size()/static_cast<double>(2*windowSize_.y);
        double zscale = domain_[Z_AXIS].size()/1000.0;
        PlotLibraryFileBase::Projection_Coordinates point1 = convertPlotCoordinatesToViewport3Projection(tgt::Vector3<plot_t>(mpoint.x - xsize*xscale,mpoint.y - ysize*yscale,mpoint.z - zsize*zscale));
        PlotLibraryFileBase::Projection_Coordinates point2 = convertPlotCoordinatesToViewport3Projection(tgt::Vector3<plot_t>(mpoint.x + xsize*xscale,mpoint.y - ysize*yscale,mpoint.z - zsize*zscale));
        PlotLibraryFileBase::Projection_Coordinates point3 = convertPlotCoordinatesToViewport3Projection(tgt::Vector3<plot_t>(mpoint.x + xsize*xscale,mpoint.y + ysize*yscale,mpoint.z - zsize*zscale));
        PlotLibraryFileBase::Projection_Coordinates point4 = convertPlotCoordinatesToViewport3Projection(tgt::Vector3<plot_t>(mpoint.x - xsize*xscale,mpoint.y + ysize*yscale,mpoint.z - zsize*zscale));
        PlotLibraryFileBase::Projection_Coordinates point5 = convertPlotCoordinatesToViewport3Projection(tgt::Vector3<plot_t>(mpoint.x - xsize*xscale,mpoint.y - ysize*yscale,mpoint.z + zsize*zscale));
        PlotLibraryFileBase::Projection_Coordinates point6 = convertPlotCoordinatesToViewport3Projection(tgt::Vector3<plot_t>(mpoint.x + xsize*xscale,mpoint.y - ysize*yscale,mpoint.z + zsize*zscale));
        PlotLibraryFileBase::Projection_Coordinates point7 = convertPlotCoordinatesToViewport3Projection(tgt::Vector3<plot_t>(mpoint.x + xsize*xscale,mpoint.y + ysize*yscale,mpoint.z + zsize*zscale));
        PlotLibraryFileBase::Projection_Coordinates point8 = convertPlotCoordinatesToViewport3Projection(tgt::Vector3<plot_t>(mpoint.x - xsize*xscale,mpoint.y + ysize*yscale,mpoint.z + zsize*zscale));
        std::vector< std::vector< PlotLibraryFileBase::Projection_Coordinates > > pointsvector;
        std::vector< PlotLibraryFileBase::Projection_Coordinates > points;
        std::vector<int> sortvector;
        std::vector<tgt::Color> fillcolorvector;
        sortvector.resize(6);
        points.resize(4);
        //Front
        points[0] = point5;
        points[1] = point6;
        points[2] = point7;
        points[3] = point8;
        pointsvector.push_back(points);
        fillcolorvector.push_back(fillColor);
        //Back
        points[0] = point1;
        points[1] = point2;
        points[2] = point3;
        points[3] = point4;
        pointsvector.push_back(points);
        fillcolorvector.push_back(tgt::Color(fillColor.r*0.8f,fillColor.g*0.8f,fillColor.b*0.8f,fillColor.a));
        //left
        points[0] = point1;
        points[1] = point4;
        points[2] = point8;
        points[3] = point5;
        pointsvector.push_back(points);
        fillcolorvector.push_back(tgt::Color(fillColor.r*0.7f,fillColor.g*0.7f,fillColor.b*0.7f,fillColor.a));
        //right
        points[0] = point2;
        points[1] = point3;
        points[2] = point7;
        points[3] = point6;
        pointsvector.push_back(points);
        fillcolorvector.push_back(tgt::Color(fillColor.r*0.8f,fillColor.g*0.8f,fillColor.b*0.8f,fillColor.a));
        //top
        points[0] = point3;
        points[1] = point4;
        points[2] = point8;
        points[3] = point7;
        pointsvector.push_back(points);
        fillcolorvector.push_back(tgt::Color(fillColor.r*0.9f,fillColor.g*0.9f,fillColor.b*0.9f,fillColor.a));
        //down
        points[0] = point1;
        points[1] = point2;
        points[2] = point6;
        points[3] = point5;
        pointsvector.push_back(points);
        fillcolorvector.push_back(tgt::Color(fillColor.r*0.9f,fillColor.g*0.9f,fillColor.b*0.9f,fillColor.a));

        for (size_t i = 0; i < fillcolorvector.size(); ++i) {
            svgPolygon(pointsvector[i],fillcolorvector[i],lineColor,clipping_number,lineWidth,fillStyle,lineStyle);
        }
    }
}

void PlotLibrarySvg::svgTriangle3D(const tgt::Vector3<plot_t>& mpoint, double height, double weight, tgt::Color fillColor, tgt::Color lineColor, int clipping_number, double lineWidth, PlotEntitySettings::PolygonFillStyle fillStyle, PlotEntitySettings::LineStyle lineStyle) {
    PlotLibraryFileBase::Projection_Coordinates point = convertPlotCoordinatesToViewport3Projection(mpoint);
    if (inClippRegion(point.afterModelview_)) {
        double xscale = domain_[X_AXIS].size()/static_cast<double>(2*windowSize_.x);
        double yscale = domain_[Y_AXIS].size()/static_cast<double>(2*windowSize_.y);
        double zscale = domain_[Z_AXIS].size()/1000.0;
        PlotLibraryFileBase::Projection_Coordinates point1 = convertPlotCoordinatesToViewport3Projection(tgt::Vector3<plot_t>(mpoint.x - weight*xscale/2.0,mpoint.y - weight*yscale/2.0,mpoint.z - height*zscale/2.0));
        PlotLibraryFileBase::Projection_Coordinates point2 = convertPlotCoordinatesToViewport3Projection(tgt::Vector3<plot_t>(mpoint.x + weight*xscale/2.0,mpoint.y - weight*yscale/2.0,mpoint.z - height*zscale/2.0));
        PlotLibraryFileBase::Projection_Coordinates point3 = convertPlotCoordinatesToViewport3Projection(tgt::Vector3<plot_t>(mpoint.x,mpoint.y + weight*yscale/2.0,mpoint.z - height*zscale/2.0));
        PlotLibraryFileBase::Projection_Coordinates point4 = convertPlotCoordinatesToViewport3Projection(tgt::Vector3<plot_t>(mpoint.x,mpoint.y,mpoint.z + height*zscale*3.0/2.0));
        std::vector< std::vector< PlotLibraryFileBase::Projection_Coordinates > > pointsvector;
        std::vector< PlotLibraryFileBase::Projection_Coordinates > points;
        std::vector< std::pair<double, int> > sortvector;
        std::vector<double> deepvector;
        points.resize(3);
        points[0] = point1;
        points[1] = point2;
        points[2] = point3;
        pointsvector.push_back(points);
        deepvector.push_back((point1.afterProjection_.z + point2.afterProjection_.z + point3.afterProjection_.z)/3);

        points[0] = point1;
        points[1] = point2;
        points[2] = point4;
        pointsvector.push_back(points);
        deepvector.push_back((point1.afterProjection_.z + point2.afterProjection_.z + point4.afterProjection_.z)/3);

        points[0] = point1;
        points[1] = point3;
        points[2] = point4;
        pointsvector.push_back(points);
        deepvector.push_back((point1.afterProjection_.z + point3.afterProjection_.z + point4.afterProjection_.z)/3);

        points[0] = point2;
        points[1] = point3;
        points[2] = point4;
        pointsvector.push_back(points);
        deepvector.push_back((point4.afterProjection_.z + point2.afterProjection_.z + point3.afterProjection_.z)/3);

        std::vector< std::pair<double, int> >::iterator it;
        bool stop;
        std::pair<double, int> item;
        for (size_t i = 0; i < deepvector.size(); ++i) {
            stop = false;
            if (sortvector.size() == 0) {
                sortvector.push_back(std::pair<double, int>(deepvector[i], static_cast<int>(i)));
            }
            else {
                for (it = sortvector.begin(); it < sortvector.end(); ++it) {
                    if (it->first <= deepvector[i]) {
                        item = std::pair<double, int>(deepvector[i], static_cast<int>(i));
                        sortvector.insert(it,item);
                        stop = true;
                        break;
                    }
                }
                if (!stop) {
                    sortvector.push_back(std::pair<double, int>(deepvector[i], static_cast<int>(i)));
                }
            }
        }
        float fakt;
        for(size_t i = 0; i < sortvector.size(); ++i) {
            fakt = 0.6f*(i+1)/sortvector.size() + 0.4f;
            svgPolygon(pointsvector[sortvector[i].second],tgt::Color(fillColor.r*fakt,fillColor.g*fakt,fillColor.b*fakt,fillColor.a),lineColor,clipping_number,lineWidth,fillStyle,lineStyle);
        }
    }
}

} // namespace
