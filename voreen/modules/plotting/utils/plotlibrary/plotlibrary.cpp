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

#include "plotlibrary.h"

#include "voreen/core/voreenapplication.h"

#include "tgt/texturemanager.h"

namespace voreen {

const std::string PlotLibrary::loggerCat_("voreen.plotting.PlotLibrary");

PlotLibrary::PlotLibrary()
    : centerAxesFlag_(false)
    , labelFont_(VoreenApplication::app()->getFontPath("Vera.ttf"))
    , windowSize_(256, 256)
    , viewPortClipping_(true)
    , dimension_(TWO)
    , drawingColor_(0.0f, 0.0f, 0.0f, 1.0f)
    , fillColor_(0.0f, 0.0f, 0.0f, 0.5f)
    , fontColor_(0.0f, 0.0f, 0.0f, 1.0f)
    , highlightColor_(1.0f, 0.95f, 0.9f, 0.5f)
    , colorMap_(ColorMap::createColdHot())
    , minGlyphSize_(1.0f)
    , maxGlyphSize_(1.0f)
    , lineWidth_(1.0f)
    , barWidth_(0.5)
    , axesWidth_(1.0f)
    , fontSize_(10)
    , squeezeFactor_(1.0)
    , lineStyle_(PlotEntitySettings::CONTINUOUS)
    , glyphStyle_(PlotEntitySettings::POINT)
    , barMode_(GROUPED)
    , marginLeft_(60)
    , marginRight_(100)
    , marginBottom_(50)
    , marginTop_(50)
    , shear_(tgt::vec2(0.f,0.f))
    , lightingFlag_(false)
    , orthographicCameraFlag_(true)
    , texture_(0)
    , texturePath_()
    , ppm_(NULL)
    , usePlotPickingManager_(false)
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

PlotLibrary::~PlotLibrary() {
}

void PlotLibrary::renderAxisLabel(Axis axis, const std::string& label) {
    switch (axis) {
        case X_AXIS:
            renderLabel(tgt::dvec3(domain_[0].getRight(), domain_[1].getLeft(),
                domain_[2].getRight()), SmartLabel::BOTTOMLEFT, label, false, 15);
            break;
        case Y_AXIS:
            renderLabel(tgt::dvec3(domain_[0].getLeft(), domain_[1].getRight(),
                domain_[2].getRight()), SmartLabel::TOPRIGHT, label, false, 15);
            break;
        case Z_AXIS:
            // should not be reached
            LWARNINGC("PlotLibrary", "this line PlotLibrary::renderAxisLabel should not be reached!");
            break;
    }
}

//// helper functions
////

void PlotLibrary::renderMousePosition(tgt::ivec2 position) {
    if (dimension_ == TWO) {
        tgt::dvec2 plotPosition = convertViewportToPlotCoordinates(position);
        std::stringstream ss;
        ss << "x: ";
        if (domain_[X_AXIS].contains(plotPosition.x))
            ss << plotPosition.x;
        else
            ss << "-";
        ss << " y: ";
        if (domain_[Y_AXIS].contains(plotPosition.y))
            ss << plotPosition.y;
        else
            ss << "-";
        renderLabel(tgt::ivec3(0, 0, 0), SmartLabel::TOPLEFT, ss.str(), true, 5);
    }
    else if (dimension_ == THREE) {
        // check if position was outside of the cube and which axis to zoom
        tgt::dvec2 center(windowSize_.x/2, windowSize_.y/2);
        tgt::dvec2 startPosition(position.x, position.y);

        int selectionEdgeIndex = -1; // -1 means no axis found
        SelectionEdge selectionEdge;
        std::stringstream ss;

        std::vector<SelectionEdge>::const_iterator it;
        // check against x-zoom-axes
        for (it = getSelectionEdgesX().begin(); it < getSelectionEdgesX().end(); ++it) {
            if (PlotLibrary::rightTurn(center, startPosition, it->startVertex_)
                && PlotLibrary::leftTurn(center, startPosition, it->endVertex_)
                && PlotLibrary::rightTurn(it->startVertex_, it->endVertex_, startPosition)) {
                selectionEdge = *it;
                selectionEdgeIndex = 0;
                ss << "x-axis: ";
            }
        }
        if (selectionEdgeIndex == -1) {
            // check against y-zoom-axes
            for (it = getSelectionEdgesY().begin(); it < getSelectionEdgesY().end(); ++it) {
                if (PlotLibrary::rightTurn(center, startPosition, it->startVertex_)
                    && PlotLibrary::leftTurn(center, startPosition, it->endVertex_)
                    && PlotLibrary::rightTurn(it->startVertex_, it->endVertex_, startPosition)) {
                    selectionEdge = *it;
                    selectionEdgeIndex = 1;
                    ss << "y-axis: ";
                }
            }
        }
        if (selectionEdgeIndex == -1) {
            // check against z-zoom-axes
            for (it = getSelectionEdgesZ().begin(); it < getSelectionEdgesZ().end(); ++it) {
                if (PlotLibrary::rightTurn(center, startPosition, it->startVertex_)
                    && PlotLibrary::leftTurn(center, startPosition, it->endVertex_)
                    && PlotLibrary::rightTurn(it->startVertex_, it->endVertex_, startPosition)) {
                    selectionEdge = *it;
                    selectionEdgeIndex = 2;
                    ss << "z-axis: ";
                }
            }
        }
        // if selectionEdgeIndex != -1, we are outside of the cube
        if (selectionEdgeIndex != -1) {
            tgtAssert(selectionEdgeIndex >= 0 && selectionEdgeIndex <= 2,
                    "PlotLibrary::renderMousePosition(): selectionEdgeIndex must be in [0, 2]");

            // determine position of the zoom clipping planes by orthogonal projection of
            // start-click onto the zoom edge (scaled dot product)
            tgt::dvec2 edge = selectionEdge.endVertex_ - selectionEdge.startVertex_;
            tgt::dvec2 click = startPosition - selectionEdge.startVertex_;
            double dot = (edge.x*click.x + edge.y*click.y) / (tgt::length(edge) * tgt::length(edge));

            // calculate plot coordinates of zoom clipping plane
            plot_t plotPosition;
            if (selectionEdge.ascOrientation_)
                plotPosition = domain_[selectionEdgeIndex].getLeft() + (dot*domain_[selectionEdgeIndex].size());
            else
                plotPosition = domain_[selectionEdgeIndex].getRight() - (dot*domain_[selectionEdgeIndex].size());

            if (domain_[selectionEdgeIndex].contains(plotPosition))
                ss << plotPosition;
            else
                ss << "-";
        }

        renderLabel(tgt::ivec3(0, 0, 0), SmartLabel::TOPLEFT, ss.str(), true, 5);
    }
}

namespace {
    /// tiny private function decoding array indices
    inline int decode(const std::valarray<int>& indices) {
        return 4*indices[0] + 2*indices[1] + indices[2];
    }
}

void PlotLibrary::findEdges(const tgt::dvec2* vertices, const std::vector< std::valarray< int > >& indices,
                            std::vector< SelectionEdge >& out, PlotLibrary::Axis axis) const {
    // for each edge
    for (int i=0; i<4; ++i) {
        // check if it is an outer edge (part of convex hull) by checking if it meets the leftTurn predicate
        // with the start points of all other edges
        if (   leftTurn(vertices[decode(indices[i])], vertices[decode(indices[i+4])], vertices[decode(indices[(i+1)%4])])
            && leftTurn(vertices[decode(indices[i])], vertices[decode(indices[i+4])], vertices[decode(indices[(i+2)%4])])
            && leftTurn(vertices[decode(indices[i])], vertices[decode(indices[i+4])], vertices[decode(indices[(i+3)%4])]))
            out.push_back(SelectionEdge(axis, vertices[decode(indices[i])], vertices[decode(indices[i+4])], true));

        // same with rightTurn predicate -> then the inverse edge is part of the convex hull
        if (   rightTurn(vertices[decode(indices[i])], vertices[decode(indices[i+4])], vertices[decode(indices[(i+1)%4])])
            && rightTurn(vertices[decode(indices[i])], vertices[decode(indices[i+4])], vertices[decode(indices[(i+2)%4])])
            && rightTurn(vertices[decode(indices[i])], vertices[decode(indices[i+4])], vertices[decode(indices[(i+3)%4])]))
            out.push_back(SelectionEdge(axis, vertices[decode(indices[i+4])], vertices[decode(indices[i])], false));
    }
}

void PlotLibrary::calculateSelectionEdges() {
    if (dimension_ != THREE)
        return;

    const plot_t& xmin = domain_[X_AXIS].getLeft();
    const plot_t& xmax = domain_[X_AXIS].getRight();
    const plot_t& ymin = domain_[Y_AXIS].getLeft();
    const plot_t& ymax = domain_[Y_AXIS].getRight();
    const plot_t& zmin = domain_[Z_AXIS].getLeft();
    const plot_t& zmax = domain_[Z_AXIS].getRight();

    // ----------------------------------------------------------------------
    // REMARK:
    //
    // We start a small battle of indices here to determine the outer edges
    // for each axis. At first sight this might look a bit odd, but doing
    // it this way save LOADS of lines of code! Otherwise we would
    // have to unroll each loop and hardcode the vertices. That would be
    // even more complex and harder to maintain.
    // ----------------------------------------------------------------------

    // index bits encoding which vertex to take:
    // most significant bit encodes x-axis, middle bit encodes y-axis, least significant bit z-axis
    // 0 is minimum, 1 is maximum,
    // decoding is done by local function decode(x, y, z) and gives according index of the vertices array
    int bits[24]  = {
     // x, y, z
        0, 0, 0,    // these index bits are aligned so that for each row i
        0, 0, 1,    // the vertices decoded by row i and row (i+4) form
        0, 1, 0,    // an x-axis edge.
        0, 1, 1,    // cyclically shifting each of these rows by one transforms
        1, 0, 0,    // the index bits so that they do the same for y-axis,
        1, 0, 1,    // and z-axis
        1, 1, 0,
        1, 1, 1
    };

    // get projected coordinates of all 8 vertices and put them into an array at the
    // right position so that they can be decoded later by the index bits
    tgt::dvec2 vertices[8] = {
        convertPlotCoordinatesToViewport(tgt::dvec3(xmin, ymin, zmin)),
        convertPlotCoordinatesToViewport(tgt::dvec3(xmin, ymin, zmax)),
        convertPlotCoordinatesToViewport(tgt::dvec3(xmin, ymax, zmin)),
        convertPlotCoordinatesToViewport(tgt::dvec3(xmin, ymax, zmax)),

        convertPlotCoordinatesToViewport(tgt::dvec3(xmax, ymin, zmin)),
        convertPlotCoordinatesToViewport(tgt::dvec3(xmax, ymin, zmax)),
        convertPlotCoordinatesToViewport(tgt::dvec3(xmax, ymax, zmin)),
        convertPlotCoordinatesToViewport(tgt::dvec3(xmax, ymax, zmax))
    };

    // fill valarrays with bits, each valarray gets one row
    std::vector<std::valarray<int> > indices;
    for (int i=0; i<8; ++i) {
        indices.push_back(std::valarray<int>(bits+(3*i), 3));
    }

    // find edges fox x-axis
    selectionEdgesX_.clear();
    findEdges(vertices, indices, selectionEdgesX_, X_AXIS);

    // cyclically shift each valarray 1 step right to get bits encoding y-axis comparisions
    for (std::vector<std::valarray<int> >::iterator it = indices.begin(); it < indices.end(); ++it)
        *it = it->cshift(-1);
    selectionEdgesY_.clear();
    findEdges(vertices, indices, selectionEdgesY_, Y_AXIS);

    // cyclically shift each valarray 1 step right to get bits encoding z-axis comparisions
    for (std::vector<std::valarray<int> >::iterator it = indices.begin(); it < indices.end(); ++it)
        *it = it->cshift(-1);
    selectionEdgesZ_.clear();
    findEdges(vertices, indices, selectionEdgesZ_, Z_AXIS);
}

const std::vector<PlotLibrary::SelectionEdge>& PlotLibrary::getSelectionEdgesX() const {
    return selectionEdgesX_;
}

const std::vector<PlotLibrary::SelectionEdge>& PlotLibrary::getSelectionEdgesY() const {
    return selectionEdgesY_;
}

const std::vector<PlotLibrary::SelectionEdge>& PlotLibrary::getSelectionEdgesZ() const {
    return selectionEdgesZ_;
}

bool PlotLibrary::leftTurn(const tgt::dvec2& a, const tgt::dvec2& b, const tgt::dvec2& c) {
    return (orientation(a, b, c) > 0);
}

bool PlotLibrary::rightTurn(const tgt::dvec2& a, const tgt::dvec2& b, const tgt::dvec2& c) {
    return (orientation(a, b, c) < 0);
}

tgt::Bounds PlotLibrary::getBoundsAbovePlot() const {
    plot_t yScale = domain_[1].size() / (windowSize_.y - marginTop_ - marginBottom_);
    return tgt::Bounds(tgt::dvec3(domain_[0].getLeft(), domain_[1].getRight(), domain_[2].getLeft()),
                       tgt::dvec3(domain_[0].getRight(), domain_[1].getRight() + yScale*marginTop_, domain_[2].getLeft()));
}

tgt::Bounds PlotLibrary::getBoundsPlot() const {
    if (dimension_ != THREE)
        return tgt::Bounds(tgt::ivec3(marginLeft_, marginBottom_, 1), tgt::ivec3(windowSize_.x-marginRight_, windowSize_.y-marginTop_, 0));
    else
        return tgt::Bounds(tgt::ivec3(0, 0, 1), tgt::ivec3(windowSize_.x, windowSize_.y, 0));
}

tgt::vec3 PlotLibrary::getCenterAbovePlot() const {
    return tgt::vec3(static_cast<float>(windowSize_.x/2), static_cast<float>(windowSize_.y-marginTop_/2), 0);
}

tgt::Bounds PlotLibrary::getBoundsBelowPlot() const {
    return tgt::Bounds(tgt::ivec3(0, 0, 0),
                       tgt::ivec3(windowSize_.x, marginBottom_, 0));
}


tgt::Bounds PlotLibrary::getBoundsRightOfPlot() const {
    return tgt::Bounds(tgt::ivec3(windowSize_.x - marginRight_, marginBottom_, 0),
                       tgt::ivec3(windowSize_.x, windowSize_.y, 0));
}

tgt::dvec2 PlotLibrary::convertPlotCoordinatesToViewport(const plot_t x, const plot_t y, const plot_t z) const {
    return convertPlotCoordinatesToViewport(tgt::dvec3(x,y,z));
}

tgt::dvec3 PlotLibrary::convertPlotCoordinatesToViewport3(const plot_t x, const plot_t y, const plot_t z) const {
    return convertPlotCoordinatesToViewport3(tgt::dvec3(x,y,z));
}

tgt::dvec2 PlotLibrary::updateScaleSteps(Axis axis) const {
    if (domain_[axis].size() <= 0)
        return tgt::dvec2(0,0);

    double viewportsize = 1;
    if (dimension_ == TWO || dimension_ == FAKETHREE) {
        viewportsize = ((axis == X_AXIS) ? windowSize_.x-marginLeft_-marginRight_ : windowSize_.y-marginTop_-marginBottom_);
    }
    else if (dimension_ == THREE) {
        tgt::dvec2 ray;
        if (axis == X_AXIS) {
            tgt::dvec2 start = convertPlotCoordinatesToViewport(tgt::dvec3(
                domain_[X_AXIS].getLeft(), domain_[Y_AXIS].getLeft(), domain_[Z_AXIS].getLeft()));
            tgt::dvec2 stop = convertPlotCoordinatesToViewport(tgt::dvec3(
                domain_[X_AXIS].getRight(), domain_[Y_AXIS].getLeft(), domain_[Z_AXIS].getLeft()));
            ray = stop - start;
        }
        else if (axis == Y_AXIS) {
            tgt::dvec2 start = convertPlotCoordinatesToViewport(tgt::dvec3(
                domain_[X_AXIS].getLeft(), domain_[Y_AXIS].getLeft(), domain_[Z_AXIS].getLeft()));
            tgt::dvec2 stop = convertPlotCoordinatesToViewport(tgt::dvec3(
                domain_[X_AXIS].getLeft(), domain_[Y_AXIS].getRight(), domain_[Z_AXIS].getLeft()));
            ray = stop - start;
        }
        else if (axis == Z_AXIS) {
            tgt::dvec2 start = convertPlotCoordinatesToViewport(tgt::dvec3(
                domain_[X_AXIS].getLeft(), domain_[Y_AXIS].getLeft(), domain_[Z_AXIS].getLeft()));
            tgt::dvec2 stop = convertPlotCoordinatesToViewport(tgt::dvec3(
                domain_[X_AXIS].getLeft(), domain_[Y_AXIS].getLeft(), domain_[Z_AXIS].getRight()));
            ray = stop - start;
        }

        viewportsize = tgt::length(ray);
    }
    // due to numerical erosion viewportsize can be 0 leading to infinite loops - we do not want that.
    if (viewportsize <= 0)
        viewportsize = .1;

    double maxStepCount = viewportsize / minimumScaleStep_[axis];
    double minStep = domain_[axis].size() / maxStepCount;
    double minBase = log10(minStep);
    if (minBase >= 0) {
        double exp = 0;
        double testStep = minStep;
        while (testStep > 7.5f) {
            testStep /= 10;
            ++exp;
        }
        if (testStep > 2.5f)
            return tgt::dvec2(5.0*pow(10, exp),pow(10, exp) / 100);
        else if (testStep > 1.5f)
            return tgt::dvec2(2.0*pow(10, exp),pow(10, exp) / 100);
        else
            return tgt::dvec2(1.0*pow(10, exp),pow(10, exp) / 100);
    }
    else {
        double exp = 0;
        double testStep = minStep;
        while (testStep < .75f) {
            testStep *= 10;
            --exp;
        }
        if (testStep < 1.5f)
            return tgt::dvec2(1.0*pow(10, exp),pow(10, exp) / 100);
        else if (testStep < 2.5f)
            return tgt::dvec2(2.0*pow(10, exp),pow(10, exp) / 100);
        else
            return tgt::dvec2(5.0*pow(10, exp),pow(10, exp) / 100);
    }
}


//
// getters and setters
//

void PlotLibrary::setDimension(Dimension dim){
    dimension_ = dim;
}

void PlotLibrary::setDomain(const Interval<plot_t>& interval, Axis axis) {
    if (interval.empty()) //this does not make much sense, so we take just [0,1]
        domain_[axis] = Interval<plot_t>(0.0, 1.0, false, false);
    else if (interval.size()==0) // not empty but interval [x, x], x \in \R
        domain_[axis] = Interval<plot_t>(interval.getLeft()-0.1, interval.getRight()+0.1, false, false);
    else
        domain_[axis] = interval;
    renewPlotToViewportScale();
}

void PlotLibrary::setLogarithmicAxis(bool polar, Axis axis){
    logarithmicAxisFlags_[axis] = polar;
}

void PlotLibrary::setCenterAxesFlag(bool value) {
    centerAxesFlag_ = value;
}

void PlotLibrary::setDrawingColor(tgt::Color color) {
    drawingColor_ = color;
}

void PlotLibrary::setFillColor(tgt::Color color) {
    fillColor_ = color;
}

void PlotLibrary::setFontColor(tgt::Color color) {
    fontColor_ = color;
}

void PlotLibrary::setHighlightColor(tgt::Color color) {
    highlightColor_ = color;
}

void PlotLibrary::setColorMap(ColorMap cm) {
    colorMap_ = cm;
}

void PlotLibrary::setLineWidth(float width) {
    if (width >= 0)
        lineWidth_ = width;
}

void PlotLibrary::setBarWidth(double width) {
    if (width >= 0 && width <= 100)
        barWidth_ = width;
}

void PlotLibrary::setAxesWidth(float width) {
    if (width >= 0)
        axesWidth_ = width;
}

void PlotLibrary::setFontSize(int size) {
    if (size > 0)
        fontSize_ = size;
}

void PlotLibrary::setLineStyle(PlotEntitySettings::LineStyle style) {
    lineStyle_ = style;
}

void PlotLibrary::setGlyphStyle(PlotEntitySettings::GlyphStyle style) {
    glyphStyle_ = style;
}


void PlotLibrary::setMarginLeft(int margin) {
    if (margin >= 0) {
        marginLeft_ = margin;
        renewPlotToViewportScale();
    }
}

void PlotLibrary::setMarginRight(int margin) {
    if (margin >= 0) {
        marginRight_ = margin;
        renewPlotToViewportScale();
    }
}

void PlotLibrary::setMarginTop(int margin) {
    if (margin >= 0) {
        marginTop_ = margin;
        renewPlotToViewportScale();
    }
}

void PlotLibrary::setMarginBottom(int margin) {
    if (margin >= 0) {
        marginBottom_ = margin;
        renewPlotToViewportScale();
    }
}

void PlotLibrary::setWindowSize(const tgt::ivec2 windowSize) {
    windowSize_ = windowSize;
    renewPlotToViewportScale();
}

void PlotLibrary::setBarGroupingMode(BarGroupingMode mode) {
    barMode_ = mode;
}

void PlotLibrary::renewPlotToViewportScale() {
    plotToViewportScale_.x = static_cast<double>(windowSize_.x-marginLeft_-marginRight_)/(domain_[0].size());
    plotToViewportScale_.y = static_cast<double>(windowSize_.y-marginBottom_-marginTop_)/(domain_[1].size());
}

void PlotLibrary::setMinimumScaleStep(int value, Axis axis) {
    if (value > 0)
        minimumScaleStep_[axis] = value;
}

void PlotLibrary::setCamera(const tgt::Camera& camera) {
    camera_ = camera;
}

void PlotLibrary::setShear(tgt::vec2 shear) {
    shear_ = shear;
}

void PlotLibrary::setSqueezeFactor(double sf) {
    squeezeFactor_ = sf;
}

void PlotLibrary::setLightingFlag(bool value) {
    lightingFlag_ = value;
}

void PlotLibrary::setOrthographicCameraFlag(bool value) {
    orthographicCameraFlag_ = value;
}

void PlotLibrary::setMinGlyphSize(float value) {
    minGlyphSize_ = value;
}

void PlotLibrary::setMaxGlyphSize(float value) {
    maxGlyphSize_ = value;
}

void PlotLibrary::setTexture(tgt::Texture* texture, const std::string& texturepath) {
    TexMgr.dispose(texture_);
    LGL_ERROR;
    texture_ = texture;
    texturePath_ = texturepath;
}


void PlotLibrary::setPlotPickingManager(PlotPickingManager* ppm) {
    ppm_ = ppm;
}

void PlotLibrary::setUsePlotPickingManager(bool value) {
    usePlotPickingManager_ = value;
}

void PlotLibrary::setViewPortClipping(bool value) {
    viewPortClipping_ = value;
}

} // namespace
