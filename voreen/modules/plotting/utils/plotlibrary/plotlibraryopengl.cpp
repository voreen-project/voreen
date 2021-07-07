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

#include "plotlibraryopengl.h"
#include "../../datastructures/plotrow.h"

#include "voreen/core/voreenapplication.h"
#include "../../interaction/plotpickingmanager.h"

#include "tgt/glmath.h"
#include "tgt/tgt_math.h"
#include "tgt/spline.h"
#include "tgt/quadric.h"
#include "tgt/immediatemode/immediatemode.h"

#include <iomanip>

namespace voreen {

const std::string PlotLibraryOpenGl::loggerCat_("voreen.plotting.PlotLibraryOpenGl");

PlotLibraryOpenGl::PlotLibraryOpenGl()
    : PlotLibraryNoneFileBase()
    , plotLabelGroup_(&labelFont_, 5, tgt::Bounds(), tgt::Color(1, 1, 1, 0.75), tgt::ivec2::zero)
    , lineLabelGroup_(&labelFont_, 6, tgt::Bounds(), tgt::ivec2::zero)
    , xAxisLabelGroup_(&labelFont_, 10, tgt::Bounds(), tgt::ivec2::zero)
    , axisLabelGroup_(&labelFont_, 6, tgt::ivec2::zero)
    , smoothPointShader_(nullptr)
    , diskMesh_()
    , triangleMesh_()
    , rectMesh_()
    , sphereMesh_()
    , tetrahedronMesh_()
    , quadMesh_()
{
    diskMesh_.setDiskGeometry(0.0f, 0.5f, 16);
    triangleMesh_.setTriangleGeometry(1);
    rectMesh_.setRectangleGeometry(1, 1);
    sphereMesh_.setSphereGeometry(0.5f, tgt::vec3::zero, tgt::vec4::one, 12);
    tetrahedronMesh_.setTetraHedronGeometry(1);
    quadMesh_.setCuboidGeometry(1, 1, 1);
}

PlotLibraryOpenGl::~PlotLibraryOpenGl() {
}

bool PlotLibraryOpenGl::setRenderStatus() {
    // helper domains
    plot_t xl = domain_[X_AXIS].getLeft(); plot_t xr = domain_[X_AXIS].getRight();
    plot_t yl = domain_[Y_AXIS].getLeft(); plot_t yr = domain_[Y_AXIS].getRight();
    plot_t zl = domain_[Z_AXIS].getLeft(); plot_t zr = domain_[Z_AXIS].getRight();
    // for 2D plots setup orthographic projection using margins
    if (dimension_ == TWO || dimension_ == FAKETHREE) {
        //check, if the canvas is big enough
        if (marginLeft_+marginRight_>=windowSize_.x || marginTop_+marginBottom_>=windowSize_.y)
            return false;
        MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
        MatStack.loadMatrix(tgt::mat4::createOrtho(xl-static_cast<double>(marginLeft_)/plotToViewportScale_.x
        , xr+static_cast<double>(marginRight_)/plotToViewportScale_.x
        , yl-static_cast<double>(marginBottom_)/plotToViewportScale_.y
        , yr+static_cast<double>(marginTop_)/plotToViewportScale_.y
        , zl-1, zr+1));
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        MatStack.loadIdentity();
        if (dimension_ == TWO) {
            //set up clipping planes
            tgt::dvec2 leftBottom = convertViewportToPlotCoordinates(tgt::ivec2(marginLeft_,marginBottom_));
            tgt::dvec2 rightTop = convertViewportToPlotCoordinates(tgt::ivec2(windowSize_.x-marginRight_,windowSize_.y-marginTop_));
            IMode.setModelSpaceClipPlaneEquation(0, tgt::plane(1.0, 0.0, 0.0, -leftBottom.x));
            IMode.setModelSpaceClipPlaneEquation(1, tgt::plane(-1.0, 0.0, 0.0, rightTop.x));
            IMode.setModelSpaceClipPlaneEquation(2, tgt::plane(0.0, 1.0, 0.0, -leftBottom.y));
            IMode.setModelSpaceClipPlaneEquation(3, tgt::plane(0.0, -1.0, 0.0, rightTop.y));
        }
        else { //dimension_ == FAKETHREE
            //translation to the center
            MatStack.translate((xl+xr)/2,(yl+yr)/2,(zl+zr)/2);
            //rescaling
            MatStack.scale(1/(1+shear_.x*(zr-zl)),1/(1+shear_.y*(zr-zl)),1);
            //shearing for fake 3d effect
            tgt::mat4 m = tgt::mat4(1,0,0,0,0,1,0,0,-shear_.x*static_cast<float>(xr-xl),
                -shear_.y*static_cast<float>(yr-yl),1,0,0,0,0,1);
            //tgt::multTransposeMatrix(m); //What was/is the purpose of this?
            MatStack.multMatrix(tgt::transpose(m));
            MatStack.translate(-(xl+xr)/2,-(yl+yr)/2,-(zl+zr)/2);
        }
    }
    // for 3D plot setup projection using the camera
    else if (dimension_ == THREE) {
        if (lightingFlag_ == true && usePlotPickingManager_ == false) {
            tgt::ImmediateMode::Material material;
            material.shininess = 0.0f;
            IMode.setMaterial(material);

            tgt::ImmediateMode::LightSource lightSource;
            lightSource.position = tgt::vec4(1.0f);
            lightSource.ambientColor = tgt::vec3(0.4f);
            lightSource.diffuseColor = tgt::vec3(0.7f);
            lightSource.specularColor = tgt::vec3(0.0f);
            IMode.setLightSource(lightSource);
        }
        MatStack.loadMatrix(camera_.getViewMatrix());
        //plot into [-0.5,0.5]^3
        MatStack.translate(-0.5,-0.5,-0.5);
        MatStack.scale(1/domain_[0].size(),1/domain_[1].size(),1/domain_[2].size());
        MatStack.translate(-domain_[0].getLeft(),-domain_[1].getLeft(),-domain_[2].getLeft());
        MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
        MatStack.loadIdentity();
        if (orthographicCameraFlag_)
           MatStack.loadMatrix(tgt::mat4::createOrtho(-.9,.9,-.9,.9,-2,10));
        else
            MatStack.loadMatrix(camera_.getProjectionMatrix(windowSize_));
        IMode.setModelSpaceClipPlaneEquation(0, tgt::plane(1.0, 0.0, 0.0, -domain_[X_AXIS].getLeft()));
        IMode.setModelSpaceClipPlaneEquation(1, tgt::plane(-1.0, 0.0, 0.0, domain_[X_AXIS].getRight()));
        IMode.setModelSpaceClipPlaneEquation(2, tgt::plane(0.0, 1.0, 0.0, -domain_[Y_AXIS].getLeft()));
        IMode.setModelSpaceClipPlaneEquation(3, tgt::plane(0.0, -1.0, 0.0, domain_[Y_AXIS].getRight()));
        IMode.setModelSpaceClipPlaneEquation(4, tgt::plane(0.0, 0.0, 1.0, -domain_[Z_AXIS].getLeft()));
        IMode.setModelSpaceClipPlaneEquation(5, tgt::plane(0.0, 0.0, -1.0, domain_[Z_AXIS].getRight()));
        // we need to calculate which are the outer edges (used for selection and
        // labeling the axes)
        calculateSelectionEdges();
    }
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    //some settings are not used for drawing pickable objects
    if (!usePlotPickingManager_) {
        //enable antialiasing
        glEnable (GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        //glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
        //glEnable(GL_POINT_SMOOTH);
        //glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
        //glEnable(GL_LINE_SMOOTH);
        //clear color buffer
        glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    }
    return true;
}

void PlotLibraryOpenGl::renderLine(const PlotData& data, int indexX, int indexY) {
    //row iterator
    std::vector<PlotRowValue>::const_iterator it = data.getRowsBegin();
    // check if only values or only tags in given cells
    bool tagsInX = (data.getColumnType(indexX) == PlotBase::STRING);
    //set up some opengl settings
    IMode.enableClipPlane(0); IMode.enableClipPlane(1);
    IMode.enableClipPlane(2); IMode.enableClipPlane(3);

    bool lineIsHighlighted = data.isHighlighted(tgt::ivec2(-1,indexY));
    glLineWidth(lineWidth_);
    glPointSize(maxGlyphSize_);
    glDisable(GL_DEPTH_TEST);
    //set line style
    if (lineStyle_!= PlotEntitySettings::CONTINUOUS) {
        LWARNING("Stippling is not supported by PlotLibraryOpenGl");
    }
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
    if (usePlotPickingManager_)
        ppm_->setGLColor(-1, indexY);
    else if (lineIsHighlighted)
        IMode.color(highlightColor_);
    else
        IMode.color(drawingColor_);
    for (++it, ++i; it != data.getRowsEnd(); ++it, ++i) {
        //we ignore rows with null entries
        if (it->getCellAt(indexX).isNull() || it->getCellAt(indexY).isNull()) {
            continue;
        }
        x = tagsInX ? i : it->getValueAt(indexX);
        y = it->getValueAt(indexY);
        IMode.begin(tgt::ImmediateMode::FAKE_LINES);
            logGlVertex2d(oldX, oldY);
            logGlVertex2d(x, y);
        IMode.end();
        //if x is out of the interval, we leave the loop
        if (x>= domain_[X_AXIS].getRight())
            break;
        oldX = x;
        oldY = y;
    }

    // render the points
    glPointSize(maxGlyphSize_);

    tgt::Shader& shader = getSmoothPointShader();
#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE
    glEnable(GL_POINT_SPRITE);
#endif
    shader.activate();
    IMode.begin(tgt::ImmediateMode::POINTS);
    for (it = data.getRowsBegin(), i = 0; it != data.getRowsEnd(); ++it, ++i) {
        //we ignore rows with null entries
        if (it->getCellAt(indexX).isNull() || it->getCellAt(indexY).isNull()) {
            continue;
        }
        x = tagsInX ? i : it->getValueAt(indexX);
        y = it->getValueAt(indexY);
        if (usePlotPickingManager_)
            ppm_->setGLColor(i, indexY);
        else if (it->getCellAt(indexY).isHighlighted())
            IMode.color(highlightColor_);
        else
            IMode.color(drawingColor_);
        logGlVertex2d(x, y);
    }
    IMode.end();
    shader.deactivate();
#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE
    glDisable(GL_POINT_SPRITE);
#endif

    IMode.disableClipPlane(0); IMode.disableClipPlane(1);
    IMode.disableClipPlane(2); IMode.disableClipPlane(3);
    glEnable(GL_DEPTH_TEST);
}

void PlotLibraryOpenGl::renderSpline(const PlotData& data, int indexX, int indexY) {
    std::vector<PlotRowValue>::const_iterator it = data.getRowsBegin(); // row iterator
    tgt::Spline spline;
    // check if only values or only tags in given cells
    bool tagsInX = (data.getColumnType(indexX) == PlotBase::STRING);
    bool lineIsHighlighted = data.isHighlighted(tgt::ivec2(-1,indexY));
    glLineWidth(lineWidth_);
    glDisable(GL_DEPTH_TEST);
    IMode.enableClipPlane(0); IMode.enableClipPlane(1);
    IMode.enableClipPlane(2); IMode.enableClipPlane(3);
    if (lineStyle_!= PlotEntitySettings::CONTINUOUS) { //set line style
        //glEnable(GL_LINE_STIPPLE);
        if (lineStyle_ == PlotEntitySettings::DOTTED)
            ;//glLineStipple(1, 0x0101);
        else //DASHED
            ;//glLineStipple(1, 0x00FF);
    }
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
    if (usePlotPickingManager_)
        ppm_->setGLColor(-1, indexY);
    else if (lineIsHighlighted)
        IMode.color(highlightColor_);
    else
        IMode.color(drawingColor_);
    IMode.begin(tgt::ImmediateMode::FAKE_LINE_STRIP);
    for (GLfloat p = 0.f; p < 1.f; p+=step) {
        x = spline.getPoint(p).x;
        y = spline.getPoint(p).y;
        logGlVertex2d(x, y);
        if (x>= domain_[X_AXIS].getRight())
            break;
    }
    IMode.end();

    // render the points
    glPointSize(maxGlyphSize_);
    IMode.begin(tgt::ImmediateMode::POINTS);
    for (it = data.getRowsBegin(), i = 0; it != data.getRowsEnd(); ++it, ++i) {
        //we ignore rows with null entries
        if (it->getCellAt(indexX).isNull() || it->getCellAt(indexY).isNull()) {
            continue;
        }
        x = tagsInX ? i : it->getValueAt(indexX);
        y = it->getValueAt(indexY);
        if (usePlotPickingManager_)
            ppm_->setGLColor(i, indexY);
        else if (it->getCellAt(indexY).isHighlighted())
            IMode.color(highlightColor_);
        else
            IMode.color(drawingColor_);
        logGlVertex2d(x, y);
    }
    IMode.end();

    IMode.disableClipPlane(0); IMode.disableClipPlane(1);
    IMode.disableClipPlane(2); IMode.disableClipPlane(3);
    glEnable(GL_DEPTH_TEST);
    //glDisable(GL_LINE_STIPPLE);
}

void PlotLibraryOpenGl::renderErrorline(const PlotData& data, int indexX, int indexY, int indexError) {
    std::vector<PlotRowValue>::const_iterator it = data.getRowsBegin();
    // check if only values or only tags in given cells
    bool tagsInX = (data.getColumnType(indexX) == PlotBase::STRING);
    glDisable(GL_DEPTH_TEST);
    IMode.enableClipPlane(0); IMode.enableClipPlane(1);
    IMode.enableClipPlane(2); IMode.enableClipPlane(3);

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
    if (usePlotPickingManager_)
        ppm_->setGLColor(-1, indexY);
    else
        IMode.color(fillColor_);
    // draw the errorline
    ++i;
    IMode.begin(tgt::ImmediateMode::TRIANGLE_STRIP);
        for (it = (++it); it != data.getRowsEnd(); ++it, ++i) {
            if (it->getCellAt(indexX).isNull() || it->getCellAt(indexY).isNull())
                continue;
            logGlVertex2d(x, errorBottom);
            logGlVertex2d(x, errorTop);

            x = tagsInX ? i : it->getValueAt(indexX);
            errorTop = it->getValueAt(indexY) + std::abs(it->getValueAt(indexError));
            errorBottom = it->getValueAt(indexY) - std::abs(it->getValueAt(indexError));
        }
        logGlVertex2d(x, errorBottom);
        logGlVertex2d(x, errorTop);
    IMode.end();
    glEnable(GL_DEPTH_TEST);
    IMode.disableClipPlane(0);  IMode.disableClipPlane(1);
    IMode.disableClipPlane(2);  IMode.disableClipPlane(3);
}

void PlotLibraryOpenGl::renderErrorspline(const PlotData& data, int indexX, int indexY, int indexError) {
    std::vector<PlotRowValue>::const_iterator it = data.getRowsBegin();
    tgt::Spline splineTop, splineBottom;

    // check if only values or only tags in given cells
    bool tagsInX = (data.getColumnType(indexX) == PlotBase::STRING);
    if (usePlotPickingManager_)
        ppm_->setGLColor(-1, indexY);
    else
        IMode.color(fillColor_);
    glDisable(GL_DEPTH_TEST);

    IMode.enableClipPlane(0); IMode.enableClipPlane(1);
    IMode.enableClipPlane(2); IMode.enableClipPlane(3);
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
    IMode.begin(tgt::ImmediateMode::TRIANGLE_STRIP);
    for (GLfloat p = 0.f; p < 1.f; p+=step) {
        logGlVertex2d(splineTop.getPoint(p).x, splineTop.getPoint(p).y);
        logGlVertex2d(splineBottom.getPoint(p).x, splineBottom.getPoint(p).y);
    }
    IMode.end();
    glEnable(GL_DEPTH_TEST);
    IMode.disableClipPlane(0);  IMode.disableClipPlane(1);
    IMode.disableClipPlane(2);  IMode.disableClipPlane(3);
}

void PlotLibraryOpenGl::renderErrorbars(const PlotData& data, int indexX, int indexY, int indexError) {
    std::vector<PlotRowValue>::const_iterator it;
    // check if only values or only tags in given cells
    bool tagsInX = (data.getColumnType(indexX) == PlotBase::STRING);
    glLineWidth(lineWidth_);
    glDisable(GL_DEPTH_TEST);
    IMode.enableClipPlane(2);
    IMode.enableClipPlane(3);

    float radius = static_cast<float>(domain_[X_AXIS].size()/(4.f*data.getRowsCount()));
    float aspectRatio = static_cast<float>(windowSize_.x)/static_cast<float>(windowSize_.y)*
        static_cast<float>(domain_[Y_AXIS].size() / domain_[X_AXIS].size());

    int i = 0;
    // draw the errorbars
    for (it = data.getRowsBegin(); it != data.getRowsEnd(); ++it, ++i) {
        if (it->getCellAt(indexX).isNull() || it->getCellAt(indexY).isNull() || it->getCellAt(indexError).isNull()) {
            continue;
        }
        double x = tagsInX ? i : it->getValueAt(indexX);
        if (!domain_[X_AXIS].contains(x))
            continue;

        if (usePlotPickingManager_)
            ppm_->setGLColor(i, indexError);
        else if (it->getCellAt(indexError).isHighlighted())
            IMode.color(highlightColor_);
        else
            IMode.color(drawingColor_);

        double y = it->getValueAt(indexY);
        double yTop = it->getValueAt(indexY) + it->getValueAt(indexError);
        double yBottom = it->getValueAt(indexY) - it->getValueAt(indexError);
        x = logarithmicAxisFlags_[X_AXIS] ? convertToLogCoordinates(x, X_AXIS) : x;
        y = logarithmicAxisFlags_[Y_AXIS] ? convertToLogCoordinates(y, Y_AXIS) : y;
        yTop = logarithmicAxisFlags_[Y_AXIS] ? convertToLogCoordinates(yTop, Y_AXIS) : yTop;
        yBottom = logarithmicAxisFlags_[Y_AXIS] ? convertToLogCoordinates(yBottom, Y_AXIS) : yBottom;
        IMode.begin(tgt::ImmediateMode::FAKE_LINES);
            IMode.vertex(x, yBottom);
            IMode.vertex(x, yTop);
        IMode.end();
        IMode.begin(tgt::ImmediateMode::FAKE_LINES);
            IMode.vertex(x-radius, yTop);
            IMode.vertex(x+radius, yTop);
        IMode.end();
        IMode.begin(tgt::ImmediateMode::FAKE_LINES);
            IMode.vertex(x-radius, yBottom);
            IMode.vertex(x+radius, yBottom);
        IMode.end();

        tgt::Ellipse midpoint(tgt::dvec3(x, y, domain_[2].getLeft()), radius/2,
                              (radius*aspectRatio)/2, tgt::dvec3(0, 0, 1), tgt::dvec3(1, 0, 0 ), 32);
        midpoint.render();
    }
    glEnable(GL_DEPTH_TEST);
    IMode.disableClipPlane(2);
    IMode.disableClipPlane(3);
}

void PlotLibraryOpenGl::renderCandlesticks(const PlotData& data, int indexX, int stickTop,
                                     int stickBottom, int candleTop, int candleBottom) {
    std::vector<PlotRowValue>::const_iterator it;

    // check if only values or only tags in given cells
    bool tagsInX = (data.getColumnType(indexX) == PlotBase::STRING);
    // openGL settings
    glDisable(GL_DEPTH_TEST);
    IMode.enableClipPlane(0); IMode.enableClipPlane(1);
    IMode.enableClipPlane(2); IMode.enableClipPlane(3);
    glLineWidth(lineWidth_);

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
        if (usePlotPickingManager_)
            ppm_->setGLColor(i, stickTop);
        else if (it->getCellAt(stickTop).isHighlighted())
            IMode.color(highlightColor_);
        else
            IMode.color(drawingColor_);
        IMode.begin(tgt::ImmediateMode::FAKE_LINES);
            logGlVertex2d(x, yStickTop);
            logGlVertex2d(x, (yStickTop+yStickBottom)/2.0);
        IMode.end();
        if (usePlotPickingManager_)
            ppm_->setGLColor(i, stickBottom);
        else if (it->getCellAt(stickBottom).isHighlighted())
            IMode.color(highlightColor_);
        else
            IMode.color(drawingColor_);
        IMode.begin(tgt::ImmediateMode::FAKE_LINES);
            logGlVertex2d(x, (yStickTop+yStickBottom)/2.0);
            logGlVertex2d(x, yStickBottom);
        IMode.end();
        //draw candle
        if (usePlotPickingManager_)
            ppm_->setGLColor(i, candleTop);
        else if (it->getCellAt(candleTop).isHighlighted())
            IMode.color(highlightColor_);
        else
            IMode.color(fillColor_);
        IMode.begin(tgt::ImmediateMode::TRIANGLE_FAN);
            logGlVertex2d(x-width, yCandleTop);
            logGlVertex2d(x-width, (yCandleBottom+yCandleTop)/2.0);
            logGlVertex2d(x+width, (yCandleBottom+yCandleTop)/2.0);
            logGlVertex2d(x+width, yCandleTop);
        IMode.end();
        if (usePlotPickingManager_)
            ppm_->setGLColor(i, candleBottom);
        else if (it->getCellAt(candleBottom).isHighlighted())
            IMode.color(highlightColor_);
        else
            IMode.color(fillColor_);
        IMode.begin(tgt::ImmediateMode::TRIANGLE_FAN);
            logGlVertex2d(x-width, (yCandleBottom+yCandleTop)/2.0);
            logGlVertex2d(x-width, yCandleBottom);
            logGlVertex2d(x+width, yCandleBottom);
            logGlVertex2d(x+width, (yCandleBottom+yCandleTop)/2.0);
        IMode.end();
    }
    glEnable(GL_DEPTH_TEST);
    IMode.disableClipPlane(0); IMode.disableClipPlane(1);
    IMode.disableClipPlane(2); IMode.disableClipPlane(3);
}

PlotLibraryOpenGl::MeshType* PlotLibraryOpenGl::createSurfaceMesh(const PlotData& data, const std::vector<int>& triangleVertexIndices, int indexX, int indexY, int indexZ, int indexCM) {
    // security check: if count of edge indices is not a multiple of 3 abort
    tgtAssert(triangleVertexIndices.size() % 3 == 0, "Count of edge indices is not a multiple of 3.");
    if (triangleVertexIndices.size() % 3 != 0) {
        LWARNING("Count of edge indices is not a multiple of 3. Aborting.");
        return nullptr;
    }

    MeshType* mesh = new MeshType();
    mesh->setPrimitiveType(GL_TRIANGLES);

    Interval<plot_t> colInterval(0, 0);
    if ( indexCM != -1 ) {
        colInterval = data.getInterval(indexCM);
        if (colInterval.size() == 0)
            indexCM = -1;
    }

    // No normal data available
    // TODO maybe get it from 3 vertices of the triangle?
    tgt::vec3 normal = tgt::vec3::zero;

    for (std::vector<int>::const_iterator it = triangleVertexIndices.begin(); it < triangleVertexIndices.end(); it += 3) {

        // Add a single triangle
        for (int i=0; i<3; ++i) {
            const PlotRowValue& row = data.getRow(*(it+i));

            // Determine color
            tgt::vec4 color;
            if (indexCM != -1 ) {
                float c = static_cast<float>((row.getValueAt(indexCM) - colInterval.getLeft()) / colInterval.size());
                color = colorMap_.getColorAtPosition(c);
            } else {
                color = drawingColor_;
            }

            if (row.getCellAt(indexZ).isHighlighted()) {
                color = highlightColor_;
            }

            // Add the vertex
            tgt::vec3 pos(row.getValueAt(indexX), row.getValueAt(indexY), row.getValueAt(indexZ));
            mesh->addVertex(pos, normal, color);
        }
    }
    return mesh;
}

PlotLibraryOpenGl::MeshType* PlotLibraryOpenGl::createSurfacePickingMesh(const PlotData& data, const std::vector<int>& triangleVertexIndices, int indexX, int indexY, int indexZ, int indexCM) {
    // security check: if count of edge indices is not a multiple of 3 abort
    tgtAssert(triangleVertexIndices.size() % 3 == 0, "Count of edge indices is not a multiple of 3.");
    if (triangleVertexIndices.size() % 3 != 0) {
        LWARNING("Count of edge indices is not a multiple of 3. Aborting.");
        return nullptr;
    }

    MeshType* mesh = new MeshType();
    mesh->setPrimitiveType(GL_TRIANGLES);

    // No normal for picking
    tgt::vec3 normal = tgt::vec3::zero;

    for (std::vector<int>::const_iterator it = triangleVertexIndices.begin(); it < triangleVertexIndices.end(); it += 3) {
        // to write out the PlotCell information to the ppm, we subdivide our triangle into three
        // quads, and render each with the according encoded color
        const PlotRowValue& rowA = data.getRow(*it);
        const PlotRowValue& rowB = data.getRow(*(it+1));
        const PlotRowValue& rowC = data.getRow(*(it+2));

        tgt::vec3 a(rowA.getValueAt(indexX), rowA.getValueAt(indexY), rowA.getValueAt(indexZ));
        tgt::vec3 b(rowB.getValueAt(indexX), rowB.getValueAt(indexY), rowB.getValueAt(indexZ));
        tgt::vec3 c(rowC.getValueAt(indexX), rowC.getValueAt(indexY), rowC.getValueAt(indexZ));

        // 3 linear interpolations
        tgt::vec3 ab = a + 0.5f*(b-a);
        tgt::vec3 bc = b + 0.5f*(c-b);
        tgt::vec3 ca = c + 0.5f*(a-c);


        tgt::dvec3 circumcenter = (a + b + c)*(1.0f/3.0f);

        tgt::vec4 color0(IMode.normalize(ppm_->getColorFromCell(*it, indexZ)), 1.0f);
        mesh->addVertex(a, normal, color0);
        mesh->addVertex(ab, normal, color0);
        mesh->addVertex(circumcenter, normal, color0);
        mesh->addVertex(ca, normal, color0);

        tgt::vec4 color1(IMode.normalize(ppm_->getColorFromCell(*it+1, indexZ)), 1.0f);
        mesh->addVertex(b, normal, color1);
        mesh->addVertex(bc, normal, color1);
        mesh->addVertex(circumcenter, normal, color1);
        mesh->addVertex(ab, normal, color1);

        tgt::vec4 color2(IMode.normalize(ppm_->getColorFromCell(*it+2, indexZ)), 1.0f);
        mesh->addVertex(c, normal, color2);
        mesh->addVertex(ca, normal, color2);
        mesh->addVertex(circumcenter, normal, color2);
        mesh->addVertex(bc, normal, color2);
    }

    return mesh;
}

PlotLibraryOpenGl::MeshType* PlotLibraryOpenGl::createHeightmapMesh(const voreen::PlotData& data, const std::vector< std::list< tgt::dvec2 > >& voronoiRegions, int indexZ, int indexCM) {
    Interval<plot_t> colInterval(0, 0);
    plot_t yMin = data.getInterval(2).getLeft();

    if ( indexCM != -1 ) {
        colInterval = data.getInterval(indexCM);
        if (colInterval.size() == 0)
            indexCM = -1;
    }
    MeshType* mesh = new MeshType();
    mesh->setPrimitiveType(GL_TRIANGLE_STRIP);
    mesh->enablePrimitiveRestart();

    // No normal data available
    // TODO maybe get it from 3 vertices of the triangle?
    tgt::vec3 normal = tgt::vec3::zero;

    int row = 0;
    std::vector< std::list< tgt::dvec2 > >::const_iterator rit = voronoiRegions.begin();
    std::vector<PlotRowValue>::const_iterator pit = data.getRowsBegin();
    for (; rit < voronoiRegions.end(); ++rit, ++pit, ++row){
        if (rit->empty())
            continue;

        tgt::vec4 color;
        if (pit->getCellAt(indexZ).isHighlighted())
            color = highlightColor_;
        else if (indexCM != -1 ) {
            float c = static_cast<float>((pit->getValueAt(indexCM) - colInterval.getLeft()) / colInterval.size());
            color = colorMap_.getColorAtPosition(c);
        }
        else
            color = drawingColor_;

        plot_t height = pit->getValueAt(indexZ);
        //clip it:
        if (height > domain_[Z_AXIS].getRight())
            height = domain_[Z_AXIS].getRight();
        if (yMin < domain_[Z_AXIS].getLeft())
            yMin = domain_[Z_AXIS].getLeft();
        if (height > yMin) {
            // render the sides of the pillar
            for (std::list< tgt::dvec2 >::const_iterator eit = rit->begin(); eit != rit->end(); ++eit) {
                mesh->addIndex(mesh->getNumVertices());
                mesh->addVertex(tgt::vec3(eit->x, eit->y, height), normal, color);
                mesh->addIndex(mesh->getNumVertices());
                mesh->addVertex(tgt::vec3(eit->x, eit->y, yMin), normal, color);
            }
            mesh->addIndex(mesh->getNumVertices());
            mesh->addVertex(tgt::vec3(rit->begin()->x, rit->begin()->y, height), normal, color);
            mesh->addIndex(mesh->getNumVertices());
            mesh->addVertex(tgt::vec3(rit->begin()->x, rit->begin()->y, yMin), normal, color);
            // Finish the pillar sides using the restart index
            mesh->addIndex(mesh->getPrimitiveRestartIndex());

            // render the top of the pillar
            for (std::list< tgt::dvec2 >::const_iterator eit = rit->begin(); eit != rit->end(); ++eit) {
                mesh->addIndex(mesh->getNumVertices());
                mesh->addVertex(tgt::vec3(eit->x, eit->y, height), normal, color);
            }
            // Finish the top of this pillar using the restart index
            mesh->addIndex(mesh->getPrimitiveRestartIndex());

            // render the bottom of the pillar
            for (std::list< tgt::dvec2 >::const_iterator eit = rit->begin(); eit != rit->end(); ++eit) {
                mesh->addIndex(mesh->getNumVertices());
                mesh->addVertex(tgt::vec3(eit->x, eit->y, yMin), normal, color);
            }
            // Finish the bottom of this pillar using the restart index
            mesh->addIndex(mesh->getPrimitiveRestartIndex());
        }
    }
    return mesh;
}

PlotLibraryOpenGl::MeshType* PlotLibraryOpenGl::createHeightmapPickingMesh(const voreen::PlotData& data, const std::vector< std::list< tgt::dvec2 > >& voronoiRegions, int indexZ, int indexCM) {
    Interval<plot_t> colInterval(0, 0);
    plot_t yMin = data.getInterval(2).getLeft();

    if ( indexCM != -1 ) {
        colInterval = data.getInterval(indexCM);
        if (colInterval.size() == 0)
            indexCM = -1;
    }
    MeshType* mesh = new MeshType();
    mesh->setPrimitiveType(GL_TRIANGLE_STRIP);
    mesh->enablePrimitiveRestart();

    // No normal data available
    // TODO maybe get it from 3 vertices of the triangle?
    tgt::vec3 normal = tgt::vec3::zero;

    int row = 0;
    std::vector< std::list< tgt::dvec2 > >::const_iterator rit = voronoiRegions.begin();
    std::vector<PlotRowValue>::const_iterator pit = data.getRowsBegin();
    for (; rit < voronoiRegions.end(); ++rit, ++pit, ++row){
        if (rit->empty())
            continue;

        tgt::vec4 color(IMode.normalize(ppm_->getColorFromCell(row, indexZ)), 1);

        plot_t height = pit->getValueAt(indexZ);
        //clip it:
        if (height > domain_[Z_AXIS].getRight())
            height = domain_[Z_AXIS].getRight();
        if (yMin < domain_[Z_AXIS].getLeft())
            yMin = domain_[Z_AXIS].getLeft();
        if (height > yMin) {
            // render the sides of the pillar
            for (std::list< tgt::dvec2 >::const_iterator eit = rit->begin(); eit != rit->end(); ++eit) {
                mesh->addIndex(mesh->getNumVertices());
                mesh->addVertex(tgt::vec3(eit->x, eit->y, height), normal, color);
                mesh->addIndex(mesh->getNumVertices());
                mesh->addVertex(tgt::vec3(eit->x, eit->y, yMin), normal, color);
            }
            mesh->addIndex(mesh->getNumVertices());
            mesh->addVertex(tgt::vec3(rit->begin()->x, rit->begin()->y, height), normal, color);
            mesh->addIndex(mesh->getNumVertices());
            mesh->addVertex(tgt::vec3(rit->begin()->x, rit->begin()->y, yMin), normal, color);
            // Finish the pillar sides using the restart index
            mesh->addIndex(mesh->getPrimitiveRestartIndex());

            // render the top of the pillar
            for (std::list< tgt::dvec2 >::const_iterator eit = rit->begin(); eit != rit->end(); ++eit) {
                mesh->addIndex(mesh->getNumVertices());
                mesh->addVertex(tgt::vec3(eit->x, eit->y, height), normal, color);
            }
            // Finish the top of this pillar using the restart index
            mesh->addIndex(mesh->getPrimitiveRestartIndex());

            // render the bottom of the pillar
            for (std::list< tgt::dvec2 >::const_iterator eit = rit->begin(); eit != rit->end(); ++eit) {
                mesh->addIndex(mesh->getNumVertices());
                mesh->addVertex(tgt::vec3(eit->x, eit->y, yMin), normal, color);
            }
            // Finish the bottom of this pillar using the restart index
            mesh->addIndex(mesh->getPrimitiveRestartIndex());
        }
    }
    return mesh;
}

void PlotLibraryOpenGl::renderSurface(const PlotData& data, const std::vector<int>& triangleVertexIndices, bool wire,
                                int indexX, int indexY, int indexZ, int indexCM, bool /*wireonly*/)  {
    // security check: if count of edge indices is not a multiple of 3 abort
    tgtAssert(triangleVertexIndices.size() % 3 == 0, "Count of edge indices is not a multiple of 3.");
    if (triangleVertexIndices.size() % 3 != 0) {
        LWARNING("Count of edge indices is not a multiple of 3. Aborting.");
        return;
    }

    IMode.enableClipPlane(0); IMode.enableClipPlane(1);
    IMode.enableClipPlane(2); IMode.enableClipPlane(3);
    IMode.enableClipPlane(4); IMode.enableClipPlane(5);

    Interval<plot_t> colInterval(0, 0);

    if (wire)
        glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
    else {
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(2.0, 1.0);
    }
    if ( indexCM != -1 ) {
        colInterval = data.getInterval(indexCM);
        if (colInterval.size() == 0)
            indexCM = -1;
    }
    if ( indexCM == -1)
        IMode.color(drawingColor_);

    glLineWidth(lineWidth_);
    // draw the triangles
    std::set<int> renderedHighlights; // we want to render each plot label // highlight only once
    for (std::vector<int>::const_iterator it = triangleVertexIndices.begin(); it < triangleVertexIndices.end(); it += 3) {
        // render plot picking data?
        if (usePlotPickingManager_) {
            // to write out the PlotCell information to the ppm, we subdivide our triangle into three
            // quads, and render each with the according encoded color
            const PlotRowValue& rowA = data.getRow(*it);
            const PlotRowValue& rowB = data.getRow(*(it+1));
            const PlotRowValue& rowC = data.getRow(*(it+2));

            tgt::dvec3 a(rowA.getValueAt(indexX), rowA.getValueAt(indexY), rowA.getValueAt(indexZ));
            tgt::dvec3 b(rowB.getValueAt(indexX), rowB.getValueAt(indexY), rowB.getValueAt(indexZ));
            tgt::dvec3 c(rowC.getValueAt(indexX), rowC.getValueAt(indexY), rowC.getValueAt(indexZ));

            // 3 linear interpolations
            tgt::vec3 ab = a + 0.5*(b-a);
            tgt::vec3 bc = b + 0.5*(c-b);
            tgt::vec3 ca = c + 0.5*(a-c);

            tgt::vec3 circumcenter = (a + b + c)*(1.0/3.0);

            ppm_->setGLColor(*it, indexZ);
            IMode.begin(tgt::ImmediateMode::TRIANGLE_STRIP);
                IMode.vertex(a);
                IMode.vertex(ab);
                IMode.vertex(circumcenter);
                IMode.vertex(ca);
            IMode.end();

            ppm_->setGLColor(*(it+1), indexZ);
            IMode.begin(tgt::ImmediateMode::TRIANGLE_STRIP);
                IMode.vertex(b);
                IMode.vertex(bc);
                IMode.vertex(circumcenter);
                IMode.vertex(ab);
            IMode.end();

            ppm_->setGLColor(*(it+2), indexZ);
            IMode.begin(tgt::ImmediateMode::TRIANGLE_STRIP);
                IMode.vertex(c);
                IMode.vertex(ca);
                IMode.vertex(circumcenter);
                IMode.vertex(bc);
            IMode.end();

        }
        // else render plot
        else {
            // render triangle first
            IMode.begin(tgt::ImmediateMode::TRIANGLE_STRIP);
            for (int i=0; i<3; ++i) {
                const PlotRowValue& row = data.getRow(*(it+i));

                if (indexCM != -1 ) {
                    float c = static_cast<float>((row.getValueAt(indexCM) - colInterval.getLeft()) / colInterval.size());
                    tgt::Color cc = colorMap_.getColorAtPosition(c);
                    IMode.color(cc);
                }
                if (!wire && row.getCellAt(indexZ).isHighlighted()) {
                    IMode.color(highlightColor_);
                    IMode.vertex(row.getValueAt(indexX), row.getValueAt(indexY), row.getValueAt(indexZ));
                    IMode.color(drawingColor_);
                }
                else {
                    IMode.vertex(row.getValueAt(indexX), row.getValueAt(indexY), row.getValueAt(indexZ));
                }
            }
            IMode.end();
        }
    }

    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
    glDisable(GL_POLYGON_OFFSET_FILL);
    IMode.disableClipPlane(0); IMode.disableClipPlane(1);
    IMode.disableClipPlane(2); IMode.disableClipPlane(3);
    IMode.disableClipPlane(4); IMode.disableClipPlane(5);
}

void PlotLibraryOpenGl::renderHeightmap(const voreen::PlotData& data, const std::vector< std::list< tgt::dvec2 > >& voronoiRegions, bool wire, int indexZ, int indexCM, bool /*wireonly*/) {
    Interval<plot_t> colInterval(0, 0);
    plot_t yMin = data.getInterval(2).getLeft();

    if (wire)
        glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
    else {
        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonOffset(2.0, 1.0);
    }
    if ( indexCM != -1 ) {
        colInterval = data.getInterval(indexCM);
        if (colInterval.size() == 0)
            indexCM = -1;
    }
    glLineWidth(lineWidth_);

    int row = 0;
    std::vector< std::list< tgt::dvec2 > >::const_iterator rit = voronoiRegions.begin();
    std::vector<PlotRowValue>::const_iterator pit = data.getRowsBegin();
    for (; rit < voronoiRegions.end(); ++rit, ++pit, ++row){
        if (rit->empty())
            continue;

        if (usePlotPickingManager_)
            ppm_->setGLColor(row, indexZ);
        else if (!wire && pit->getCellAt(indexZ).isHighlighted())
            IMode.color(highlightColor_);
        else if (indexCM != -1 ) {
            float c = static_cast<float>((pit->getValueAt(indexCM) - colInterval.getLeft()) / colInterval.size());
            IMode.color(colorMap_.getColorAtPosition(c));
        }
        else
            IMode.color(drawingColor_);

        plot_t height = pit->getValueAt(indexZ);
        //clip it:
        if (height > domain_[Z_AXIS].getRight())
            height = domain_[Z_AXIS].getRight();
        if (yMin < domain_[Z_AXIS].getLeft())
            yMin = domain_[Z_AXIS].getLeft();
        if (height > yMin) {
            // render the sides of the pillar
            IMode.begin(tgt::ImmediateMode::TRIANGLE_STRIP);
            for (std::list< tgt::dvec2 >::const_iterator eit = rit->begin(); eit != rit->end(); ++eit) {
                IMode.vertex(eit->x, eit->y, height);
                IMode.vertex(eit->x, eit->y, yMin);
            }
            IMode.vertex(rit->begin()->x, rit->begin()->y, height);
            IMode.vertex(rit->begin()->x, rit->begin()->y, yMin);
            IMode.end();
            // render the top of the pillar
            IMode.begin(tgt::ImmediateMode::TRIANGLE_STRIP);
            for (std::list< tgt::dvec2 >::const_iterator eit = rit->begin(); eit != rit->end(); ++eit) {
                IMode.vertex(eit->x, eit->y, height);
            }
            IMode.end();
            // render the bottom of the pillar
            IMode.begin(tgt::ImmediateMode::TRIANGLE_STRIP);
            for (std::list< tgt::dvec2 >::const_iterator eit = rit->begin(); eit != rit->end(); ++eit) {
                IMode.vertex(eit->x, eit->y, yMin);
            }
            IMode.end();
        }
    }
    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
    glDisable(GL_POLYGON_OFFSET_FILL);
}

void PlotLibraryOpenGl::renderBars(const PlotData& data, std::vector<int> indicesY) {
    std::vector<PlotRowValue>::const_iterator rowIt = data.getRowsBegin();
    //stores y values and indices of a merged bar group, the indices are used to get the right color
    std::vector<std::pair<plot_t, tgt::Color> > mergedBars;
    //stores last y value, used for stacked bars
    plot_t lastY;
    tgt::Color c;
    glLineWidth(lineWidth_);
    //rowcounter
    plot_t row = 0;
    for (; rowIt < data.getRowsEnd(); ++rowIt) {
        lastY = 0;
        mergedBars.clear();
        // we do not use the iterator because we also iterate through the colormap
        for (size_t i = 0; i < indicesY.size(); ++i) {
            if (rowIt->getCellAt(indicesY.at(i)).isHighlighted())
                c = highlightColor_;
            else
                c = colorMap_.getColorAtIndex(static_cast<int>(i));
            if (barMode_ == STACKED) {
                if (rowIt->getCellAt(indicesY.at(i)).isNull() )
                    continue;
                plot_t newY = lastY+rowIt->getValueAt(indicesY.at(i));
                //negative stacked bars do not make any sense and are ignored
                if (newY < lastY)
                    continue;
                if (usePlotPickingManager_)
                    ppm_->setGLColor(static_cast<int>(row), indicesY.at(i));
                renderSingleBar(row-barWidth_/2.0,row+barWidth_/2.0, lastY, newY, c);
                lastY = newY;
            }
            else if (barMode_ == GROUPED) {
                if (rowIt->getCellAt(indicesY.at(i)).isNull())
                    continue;
                double singleBarWidth = barWidth_/(1.0*static_cast<double>(indicesY.size()));
                if (usePlotPickingManager_)
                    ppm_->setGLColor(static_cast<int>(row), indicesY.at(i));
                renderSingleBar(row-barWidth_/2.0+static_cast<double>(i)*singleBarWidth,
                    row-barWidth_/2.0+static_cast<double>(i+1)*singleBarWidth, 0, rowIt->getValueAt(indicesY.at(i)), c);
            }

            else { // MERGED
                // we can't skip null entries, so we set them 0
                if (usePlotPickingManager_)
                    c = ppm_->convertColor(ppm_->getColorFromCell(static_cast<int>(row), indicesY.at(i)));
                if (rowIt->getCellAt(indicesY.at(i)).isNull())
                    mergedBars.push_back(std::pair<plot_t, tgt::Color>(0, c));
                // push the y value and the color index
                mergedBars.push_back(std::pair<plot_t, tgt::Color>(rowIt->getValueAt(indicesY.at(i)), c));
            }
        }
        if (barMode_ == MERGED) {
            // the values are stored in bars, but not yet drawn
            std::sort(mergedBars.begin(), mergedBars.end(), mergedBarSorter);
            std::vector<std::pair<plot_t, tgt::Color> >::const_iterator it;
            double squeeze = 1.0;
            for (it = mergedBars.begin(); it < mergedBars.end(); ++it) {
                IMode.color(it->second);
                renderSingleBar(row-barWidth_/2.0,row+barWidth_/2.0, 0, it->first, it->second, squeeze);
                squeeze = squeezeFactor_*squeeze;
            }
        }
        ++row;
    }
}

void PlotLibraryOpenGl::renderNodeGraph(const PlotData& nodeData, const PlotData& connectionData, int indexX, int indexY, int indexDx, int indexDy) {
    std::vector<PlotRowValue>::const_iterator it;

    // render nodes
    plot_t glyphSize = (maxGlyphSize_ + minGlyphSize_)/2;
    int i = 0;
    for (it = nodeData.getRowsBegin(); it != nodeData.getRowsEnd(); ++it, ++i) {
        // render node
        IMode.color(drawingColor_);
        renderGlyph(it->getValueAt(indexX), it->getValueAt(indexY), 0, glyphSize);

        // render force vector
        IMode.color(fillColor_);
        IMode.begin(tgt::ImmediateMode::FAKE_LINE_STRIP);
            logGlVertex2d(it->getValueAt(indexX), it->getValueAt(indexY));
            logGlVertex2d(it->getValueAt(indexX) + it->getValueAt(indexDx), it->getValueAt(indexY) + it->getValueAt(indexDy));
        IMode.end();

        // render node label
        std::stringstream ss;
        ss << i;
        renderLabel(tgt::dvec3(it->getValueAt(indexX), it->getValueAt(indexY), 0), SmartLabel::CENTERED, ss.str(), false, 0);
    }

    std::vector<PlotCellValue> tester;
    tester.push_back(PlotCellValue(0));
    std::vector<PlotRowValue>::const_iterator firstIt, secondIt;

    // render connections
    glLineWidth(lineWidth_);
    IMode.color(drawingColor_);
    for (it = connectionData.getRowsBegin(); it != connectionData.getRowsEnd(); ++it) {
        tester[0].setValue(static_cast<int>(it->getValueAt(0)));
        firstIt = nodeData.lower_bound(tester);

        tester[0].setValue(static_cast<int>(it->getValueAt(1)));
        secondIt = nodeData.lower_bound(tester);

        if (firstIt != nodeData.getRowsEnd() && secondIt != nodeData.getRowsEnd()) {
            IMode.begin(tgt::ImmediateMode::FAKE_LINE_STRIP);
                logGlVertex2d(firstIt->getValueAt(indexX), firstIt->getValueAt(indexY));
                logGlVertex2d(secondIt->getValueAt(indexX), secondIt->getValueAt(indexY));
            IMode.end();
        }
    }
}

void PlotLibraryOpenGl::renderColorMapLegend(const PlotData& data, int column, int number) {
    // ColorMaps with less than 2 colors may not exist
    if (colorMap_.getColorCount() < 2)
        return;

    Interval<plot_t>   interval = data.getInterval(column);
    const std::string& label    = data.getColumnLabel(column);

    // switch to viewport coordinates
    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.pushMatrix();
    MatStack.loadMatrix(tgt::mat4::createOrtho(0.0f, windowSize_.x, 0.f, windowSize_.y, -1.0f, 1.0f));

    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();
    MatStack.loadIdentity();

    glDisable(GL_DEPTH_TEST);

    // render legend
    int colorcount      = colorMap_.getColorCount();
    const double width  = 96;
    const double height = 16;
    const double xStart = windowSize_.x - width - 8;
    const double yStart = windowSize_.y - 8 - (number*32);
    double stepWidth    = width / (colorcount - 1);

    // color map
    IMode.color(colorMap_.getColorAtIndex(0));
    for (int i = 0; i < colorcount-1; ++i) {
        IMode.begin(tgt::ImmediateMode::TRIANGLE_FAN);
            IMode.vertex(xStart + (i*stepWidth), yStart);
            IMode.vertex(xStart + (i*stepWidth), yStart - height);
            IMode.color(colorMap_.getColorAtIndex(i + 1));
            IMode.vertex(xStart + ((i+1) * stepWidth), yStart - height);
            IMode.vertex(xStart + ((i+1) * stepWidth), yStart);
        IMode.end();
    }

    // bounding box
    IMode.color(drawingColor_);
    glLineWidth(lineWidth_);
    IMode.begin(tgt::ImmediateMode::FAKE_LINE_STRIP);
        IMode.vertex(xStart, yStart);
        IMode.vertex(xStart, yStart - height);
        IMode.vertex(xStart + width, yStart - height);
        IMode.vertex(xStart + width, yStart);
        IMode.vertex(xStart, yStart);
    IMode.end();

    // labels
    SmartLabelGroupBaseOpenGl::renderSingleLabel(&labelFont_, label, tgt::dvec3(xStart, yStart - (height/2), 0),
            fontColor_, fontSize_, SmartLabel::MIDDLERIGHT, 4, windowSize_);

    std::stringstream ss;
    ss << std::setprecision(4) << interval.getLeft();
    SmartLabelGroupBaseOpenGl::renderSingleLabel(&labelFont_, ss.str(), tgt::dvec3(xStart, yStart - height, 0),
            fontColor_, fontSize_, SmartLabel::BOTTOMCENTERED, 4, windowSize_);

    ss.str("");
    ss.clear();
    ss << std::setprecision(4) << interval.getRight();
    SmartLabelGroupBaseOpenGl::renderSingleLabel(&labelFont_, ss.str(), tgt::dvec3(xStart + width, yStart - height, 0),
            fontColor_, fontSize_, SmartLabel::BOTTOMCENTERED, 4, windowSize_);

    glEnable(GL_DEPTH_TEST);

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.popMatrix();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.popMatrix();
    LGL_ERROR;
}

void PlotLibraryOpenGl::renderScatter(const PlotData& data, int indexX, int indexY, int indexZ, int indexCM, int indexSize) {
    if (lightingFlag_ == true && usePlotPickingManager_ == false && glyphStyle_ != PlotEntitySettings::POINT) {
        IMode.enableLighting();
    }
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
            if (usePlotPickingManager_)
                ppm_->setGLColor(i,indexZ == -1 ? indexY : indexZ);
            else if ((indexZ != -1 && it->getCellAt(indexZ).isHighlighted())
                    || (indexZ == -1 && it->getCellAt(indexY).isHighlighted()))
                IMode.color(highlightColor_);
            else if (indexCM != -1 ) {
                float c = static_cast<float>((it->getValueAt(indexCM) - colInterval.getLeft()) / colInterval.size());
                tgt::Color cc = colorMap_.getColorAtPosition(c);
                IMode.color(cc);
            }
            else
                IMode.color(drawingColor_);
            // set size
            if (indexSize != -1 ) {
                size = minGlyphSize_ + (maxGlyphSize_ - minGlyphSize_) *
                            (it->getValueAt(indexSize) - sizeInterval.getLeft()) / sizeInterval.size();
            }
            renderGlyph(x, y, z, size);
        }
    }
    if (lightingFlag_ == true && usePlotPickingManager_ == false && glyphStyle_ != PlotEntitySettings::POINT) {
        IMode.disableLighting();
    }
}

void PlotLibraryOpenGl::renderAxes() {
    // axes
    if (dimension_ == TWO)
        glDisable(GL_DEPTH_TEST);

    glLineWidth(axesWidth_);
    IMode.color(drawingColor_);

    plot_t xl = domain_[0].getLeft();    plot_t xr = domain_[0].getRight();
    plot_t yl = domain_[1].getLeft();    plot_t yr = domain_[1].getRight();
    plot_t zl = domain_[2].getLeft();    plot_t zr = domain_[2].getRight();

    if (! centerAxesFlag_) {
        //x and y axes
        IMode.begin(tgt::ImmediateMode::FAKE_LINE_STRIP);
            IMode.vertex(xl, yr, zl);
            IMode.vertex(xl, yl, zl);
            IMode.vertex(xr, yl, zl);
        IMode.end();
    }

    if (dimension_ == TWO) {
    //draw arrows with viewport coordinates
        int arrowSize = 5; // in pixel
        MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
        MatStack.pushMatrix();
        MatStack.loadMatrix(tgt::mat4::createOrtho(0.0f, windowSize_.x, 0.f, windowSize_.y, -1.0f, 1.0f));
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
        MatStack.pushMatrix();
        MatStack.loadIdentity();
        IMode.begin(tgt::ImmediateMode::FAKE_LINES);
            IMode.vertex(windowSize_.x-marginRight_, marginBottom_);
            IMode.vertex(windowSize_.x-marginRight_+4*arrowSize, marginBottom_);
            IMode.vertex(windowSize_.x-marginRight_+2*arrowSize, marginBottom_+arrowSize);
            IMode.vertex(windowSize_.x-marginRight_+4*arrowSize, marginBottom_);
            IMode.vertex(windowSize_.x-marginRight_+2*arrowSize, marginBottom_-arrowSize);
            IMode.vertex(windowSize_.x-marginRight_+4*arrowSize, marginBottom_);

            IMode.vertex(marginLeft_, windowSize_.y-marginTop_);
            IMode.vertex(marginLeft_, windowSize_.y-marginTop_+4*arrowSize);
            IMode.vertex(marginLeft_ + arrowSize, windowSize_.y-marginTop_+2*arrowSize);
            IMode.vertex(marginLeft_, windowSize_.y-marginTop_+4*arrowSize);
            IMode.vertex(marginLeft_ - arrowSize, windowSize_.y-marginTop_+2*arrowSize);
            IMode.vertex(marginLeft_, windowSize_.y-marginTop_+4*arrowSize);
        IMode.end();
        MatStack.popMatrix();
        MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
        MatStack.popMatrix();
        MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    }
    else if (dimension_ == THREE) {
        if (centerAxesFlag_) {
            IMode.begin(tgt::ImmediateMode::FAKE_LINES);
                IMode.vertex(xl, 0, 0); IMode.vertex(xr, 0, 0);
                IMode.vertex(0, yl, 0); IMode.vertex(0, yr, 0);
                IMode.vertex(0, 0,  0); IMode.vertex(0, 0, zr);
            IMode.end();
        }
        else {
            //draw cube mesh
            IMode.begin(tgt::ImmediateMode::FAKE_LINE_STRIP);
                IMode.vertex(xr, yl, zl); IMode.vertex(xr, yr, zl);
                IMode.vertex(xl, yr, zl); IMode.vertex(xl, yr, zr);
                IMode.vertex(xr, yr, zr); IMode.vertex(xr, yr, zl);
            IMode.end();
            IMode.begin(tgt::ImmediateMode::FAKE_LINE_STRIP);
                IMode.vertex(xr, yl, zl); IMode.vertex(xr, yl, zr);
                IMode.vertex(xl, yl, zr); IMode.vertex(xl, yl, zl);
            IMode.end();
            IMode.begin(tgt::ImmediateMode::FAKE_LINES);
                IMode.vertex(xr, yl, zr); IMode.vertex(xr, yr, zr);
                IMode.vertex(xl, yl, zr); IMode.vertex(xl, yr, zr);
            IMode.end();
        }
    }
    else if (dimension_ == FAKETHREE) {
        //draw back
        IMode.begin(tgt::ImmediateMode::LINE_LOOP);
            IMode.vertex(xl, yr, zl); IMode.vertex(xl, yl, zl);
            IMode.vertex(xr, yl, zl); IMode.vertex(xr, yr, zl);
        IMode.end();

        //draw bottom
        IMode.begin(tgt::ImmediateMode::FAKE_LINES);
            IMode.vertex(xl, yl, zl); IMode.vertex(xl, yl, zr);
            IMode.vertex(xr, yl, zr); IMode.vertex(xr, yl, zl);
        IMode.end();

        //draw left
         IMode.begin(tgt::ImmediateMode::LINE_LOOP);
            IMode.vertex(xl, yr, zr); IMode.vertex(xl, yl, zr);
            IMode.vertex(xl, yl, zl); IMode.vertex(xl, yr, zl);
        IMode.end();

        //draw zero
        IMode.begin(tgt::ImmediateMode::FAKE_LINES);
            IMode.vertex(xl, 0, zl); IMode.vertex(xl, 0, zr);
            IMode.vertex(xr, 0, zr); IMode.vertex(xr, 0, zl);
        IMode.end();

        //the front is always above the plot
        glDisable(GL_DEPTH_TEST);
        IMode.begin(tgt::ImmediateMode::FAKE_LINES);
            IMode.vertex(xl, 0, zr);
            IMode.vertex(xr, 0, zr);
        IMode.end();
        IMode.begin(tgt::ImmediateMode::FAKE_LINES);
            IMode.vertex(xl, yl, zr);
            IMode.vertex(xr, yl, zr);
        IMode.end();
    }
    glEnable(GL_DEPTH_TEST);
}

void PlotLibraryOpenGl::renderAxisScales(Axis axis, bool helperLines, const std::string& label, plot_t offset) {
    tgt::dvec2 step = updateScaleSteps(axis);
    glLineWidth(axesWidth_/2.f);
    IMode.color(drawingColor_);
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
                    IMode.begin(tgt::ImmediateMode::FAKE_LINE_STRIP);
                        logGlVertex2d(i, yl);
                        logGlVertex2d(i, yr);
                    IMode.end();
                }
            }
            else if (axis == Y_AXIS) {
                renderLabel(tgt::dvec3(xl, i, zr), SmartLabel::MIDDLERIGHT, stream.str());
                if (helperLines) {
                    IMode.begin(tgt::ImmediateMode::FAKE_LINE_STRIP);
                        if (dimension_ == FAKETHREE)
                            logGlVertex3d(xl, i, zr);
                        logGlVertex3d(xl, i, zl);
                        logGlVertex3d(xr, i, zl);
                    IMode.end();
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
        std::vector<SelectionEdge>::const_iterator minEdge = getSelectionEdgesX().begin();
        std::vector<SelectionEdge>::const_iterator endEdge = getSelectionEdgesX().end();
        if (axis == X_AXIS) {
            if (getSelectionEdgesX().empty())
                return;
        }
        else if (axis == Y_AXIS) {
            if (getSelectionEdgesY().empty())
                return;
            minEdge = getSelectionEdgesY().begin();
            endEdge = getSelectionEdgesY().end();
        }
        else if (axis == Z_AXIS) {
            if (getSelectionEdgesZ().empty())
                return;
            minEdge = getSelectionEdgesZ().begin();
            endEdge = getSelectionEdgesZ().end();
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

            // set opengl state for occlusion query
            glEnable(GL_CULL_FACE);
            glCullFace(GL_FRONT);
            glDepthMask(GL_FALSE);
            glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);

            GLuint occlusionQueries[6];
            glGenQueries(6, occlusionQueries);

            tgt::vec3 geomLlf(xl, yl, zl);
            tgt::vec3 geomUrb(xr, yr, zr);


            // set opengl state for occlusion query
            glEnable(GL_CULL_FACE);
            glCullFace(GL_FRONT);
            glDepthMask(GL_FALSE);
            glColorMask(GL_FALSE, GL_FALSE, GL_FALSE, GL_FALSE);

            // execute occlusion queries for each face
            glBeginQuery(GL_ANY_SAMPLES_PASSED, occlusionQueries[0]);
            IMode.begin(tgt::ImmediateMode::QUADS);
                IMode.vertex(geomLlf);
                IMode.vertex(tgt::vec3(geomLlf[0], geomUrb[1], geomLlf[2]));
                IMode.vertex(tgt::vec3(geomUrb[0], geomUrb[1], geomLlf[2]));
                IMode.vertex(tgt::vec3(geomUrb[0], geomLlf[1], geomLlf[2]));
            IMode.end();
            glEndQuery(GL_ANY_SAMPLES_PASSED);

            glBeginQuery(GL_ANY_SAMPLES_PASSED, occlusionQueries[1]);
            IMode.begin(tgt::ImmediateMode::QUADS);
                IMode.vertex(tgt::vec3(geomLlf[0], geomLlf[1], geomUrb[2]));
                IMode.vertex(tgt::vec3(geomUrb[0], geomLlf[1], geomUrb[2]));
                IMode.vertex(tgt::vec3(geomUrb[0], geomUrb[1], geomUrb[2]));
                IMode.vertex(tgt::vec3(geomLlf[0], geomUrb[1], geomUrb[2]));
            IMode.end();
            glEndQuery(GL_ANY_SAMPLES_PASSED);

            glBeginQuery(GL_ANY_SAMPLES_PASSED, occlusionQueries[2]);
            IMode.begin(tgt::ImmediateMode::QUADS);
                IMode.vertex(tgt::vec3(geomLlf[0], geomLlf[1], geomLlf[2]));
                IMode.vertex(tgt::vec3(geomLlf[0], geomLlf[1], geomUrb[2]));
                IMode.vertex(tgt::vec3(geomLlf[0], geomUrb[1], geomUrb[2]));
                IMode.vertex(tgt::vec3(geomLlf[0], geomUrb[1], geomLlf[2]));
            IMode.end();
            glEndQuery(GL_ANY_SAMPLES_PASSED);

            glBeginQuery(GL_ANY_SAMPLES_PASSED, occlusionQueries[3]);
            IMode.begin(tgt::ImmediateMode::QUADS);
                IMode.vertex(tgt::vec3(geomUrb[0], geomLlf[1], geomUrb[2]));
                IMode.vertex(tgt::vec3(geomUrb[0], geomLlf[1], geomLlf[2]));
                IMode.vertex(tgt::vec3(geomUrb[0], geomUrb[1], geomLlf[2]));
                IMode.vertex(tgt::vec3(geomUrb[0], geomUrb[1], geomUrb[2]));
            IMode.end();
            glEndQuery(GL_ANY_SAMPLES_PASSED);

            glBeginQuery(GL_ANY_SAMPLES_PASSED, occlusionQueries[4]);
            IMode.begin(tgt::ImmediateMode::QUADS);
                IMode.vertex(tgt::vec3(geomLlf[0], geomLlf[1], geomLlf[2]));
                IMode.vertex(tgt::vec3(geomUrb[0], geomLlf[1], geomLlf[2]));
                IMode.vertex(tgt::vec3(geomUrb[0], geomLlf[1], geomUrb[2]));
                IMode.vertex(tgt::vec3(geomLlf[0], geomLlf[1], geomUrb[2]));
            IMode.end();
            glEndQuery(GL_ANY_SAMPLES_PASSED);

            glBeginQuery(GL_ANY_SAMPLES_PASSED, occlusionQueries[5]);
            IMode.begin(tgt::ImmediateMode::QUADS);
                IMode.vertex(tgt::vec3(geomLlf[0], geomUrb[1], geomLlf[2]));
                IMode.vertex(tgt::vec3(geomLlf[0], geomUrb[1], geomUrb[2]));
                IMode.vertex(tgt::vec3(geomUrb[0], geomUrb[1], geomUrb[2]));
                IMode.vertex(tgt::vec3(geomUrb[0], geomUrb[1], geomLlf[2]));
            IMode.end();
            glEndQuery(GL_ANY_SAMPLES_PASSED);

            // set back opengl state
            glDepthMask(GL_TRUE);
            glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
            glDisable(GL_CULL_FACE);
            glCullFace(GL_BACK);

            // retrieve the result for each occlusion query
            int renderFaces[6];

            // get result of each query
            for (size_t i = 0; i < 6; ++i) {
                int queryReady = false;
                while(!queryReady)
                    glGetQueryObjectiv(occlusionQueries[i], GL_QUERY_RESULT_AVAILABLE, &queryReady);
               glGetQueryObjectiv(occlusionQueries[i], GL_QUERY_RESULT, &renderFaces[i]);
            }

            IMode.begin(tgt::ImmediateMode::FAKE_LINES);
            if (axis == X_AXIS) {
                for (plot_t i = start; renderFaces[0] && i <= domain_[axis].getRight(); i += step.x) {
                    IMode.vertex(i, yl, zl);
                    IMode.vertex(i, yr, zl);
                }
                for (plot_t i = start; renderFaces[1] && i <= domain_[axis].getRight(); i += step.x) {
                    IMode.vertex(i, yl, zr);
                    IMode.vertex(i, yr, zr);
                }
                for (plot_t i = start; renderFaces[4] && i <= domain_[axis].getRight(); i += step.x) {
                    IMode.vertex(i, yl, zl);
                    IMode.vertex(i, yl, zr);
                }
                for (plot_t i = start; renderFaces[5] && i <= domain_[axis].getRight(); i += step.x) {
                    IMode.vertex(i, yr, zl);
                    IMode.vertex(i, yr, zr);
                }
            }
            else if (axis == Y_AXIS) {
                for (plot_t i = start; renderFaces[0] && i <= domain_[axis].getRight(); i += step.x) {
                    IMode.vertex(xl, i, zl);
                    IMode.vertex(xr, i, zl);
                }
                for (plot_t i = start; renderFaces[1] && i <= domain_[axis].getRight(); i += step.x) {
                    IMode.vertex(xl, i, zr);
                    IMode.vertex(xr, i, zr);
                }
                for (plot_t i = start; renderFaces[2] && i <= domain_[axis].getRight(); i += step.x) {
                    IMode.vertex(xl, i, zl);
                    IMode.vertex(xl, i, zr);
                }
                for (plot_t i = start; renderFaces[3] && i <= domain_[axis].getRight(); i += step.x) {
                    IMode.vertex(xr, i, zl);
                    IMode.vertex(xr, i, zr);
                }
            }
            else if (axis == Z_AXIS) {
                for (plot_t i = start; renderFaces[2] && i <= domain_[axis].getRight(); i += step.x) {
                    IMode.vertex(xl, yl, i);
                    IMode.vertex(xl, yr, i);
                }
                for (plot_t i = start; renderFaces[3] && i <= domain_[axis].getRight(); i += step.x) {
                    IMode.vertex(xr, yl, i);
                    IMode.vertex(xr, yr, i);
                }
                for (plot_t i = start; renderFaces[4] && i <= domain_[axis].getRight(); i += step.x) {
                    IMode.vertex(xl, yl, i);
                    IMode.vertex(xr, yl, i);
                }
                for (plot_t i = start; renderFaces[5] && i <= domain_[axis].getRight(); i += step.x) {
                    IMode.vertex(xl, yr, i);
                    IMode.vertex(xr, yr, i);
                }
            }
            IMode.end();

            glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
            glDisable(GL_CULL_FACE);
        }
    }
}

void PlotLibraryOpenGl::renderAxisLabelScales(const PlotData& data, int indexLabel, bool helperLines) {
    std::string label;
    xAxisLabelGroup_.reset();
    xAxisLabelGroup_.setBounds(getBoundsBelowPlot());
    std::vector<PlotRowValue>::const_iterator it = data.getRowsBegin();
    plot_t x = 0;
    glLineWidth(axesWidth_/2.f);
    IMode.color(drawingColor_);
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
            IMode.begin(tgt::ImmediateMode::FAKE_LINE_STRIP);
                logGlVertex2d(x, domain_[1].getLeft());
                logGlVertex2d(x, domain_[1].getRight());
            IMode.end();
        }
        x += 1;
    }
    renderSmartLabelGroup(&xAxisLabelGroup_);
}

void PlotLibraryOpenGl::renderLabel(tgt::vec3 pos, const SmartLabel::Alignment align, const std::string& text,
                              bool viewCoordinates, int padding) {
    if (!viewCoordinates)
        pos = convertPlotCoordinatesToViewport3(pos);

    glDisable(GL_DEPTH_TEST);
    SmartLabelGroupBaseOpenGl::renderSingleLabel(&labelFont_, text, pos, fontColor_, fontSize_, align, static_cast<float>(padding), windowSize_);
    glEnable(GL_DEPTH_TEST);
    LGL_ERROR;
}

void PlotLibraryOpenGl::renderLabel(tgt::dvec2 pos, const SmartLabel::Alignment align, const std::string& text, int padding) {

    glDisable(GL_DEPTH_TEST);
    SmartLabelGroupBaseOpenGl::renderSingleLabel(&labelFont_, text, tgt::vec3((float)pos.x, (float)pos.y, 0), fontColor_,
            fontSize_, align, static_cast<float>(padding), windowSize_);
    glEnable(GL_DEPTH_TEST);

    LGL_ERROR;
}

void PlotLibraryOpenGl::addPlotLabel(std::string text, tgt::vec3 position, tgt::Color color,
                               int size, SmartLabel::Alignment align) {
    plotLabelGroup_.addLabel(text, position, color, size, align);
}

void PlotLibraryOpenGl::addLineLabel(std::string text, tgt::vec3 position, tgt::Color color,
                               int size, SmartLabel::Alignment align) {
    lineLabelGroup_.addLabel(text, position, color, size, align);
}

//
// helper functions
//

void PlotLibraryOpenGl::resetLineLabels() {
    lineLabelGroup_.reset();
    lineLabelGroup_.setBounds(getBoundsRightOfPlot());
}

void PlotLibraryOpenGl::renderLineLabels() {
    renderSmartLabelGroup(&lineLabelGroup_);
}

void PlotLibraryOpenGl::resetPlotLabels() {
    plotLabelGroup_.reset();
    plotLabelGroup_.setBounds(getBoundsPlot());
}

void PlotLibraryOpenGl::renderPlotLabels() {
    renderSmartLabelGroup(&plotLabelGroup_);
}

void PlotLibraryOpenGl::renderSmartLabelGroup(SmartLabelGroupBase* smg) {
    smg->performLayout();

    // HACK(apv): opengl labels need the screensize for rendering
    SmartLabelGroupBaseOpenGl * smggl = dynamic_cast<SmartLabelGroupBaseOpenGl*>(smg);
    if (smggl){
        smggl->setScreenSize(windowSize_);
    }

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.pushMatrix();
    MatStack.loadMatrix(tgt::mat4::createOrtho(0.0f, windowSize_.x, 0.f, windowSize_.y, -1.0f, 1.0f));

    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.pushMatrix();
    MatStack.loadIdentity();

    glDisable(GL_DEPTH_TEST);
    smg->render();
    glEnable(GL_DEPTH_TEST);

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.popMatrix();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.popMatrix();
    LGL_ERROR;
}


void PlotLibraryOpenGl::resetRenderStatus() {
    glClearDepth(1.0);
    glClearColor(0.0, 0.0, 0.0, 0.0);
    glLineWidth(1.f);
    glEnable(GL_DEPTH_TEST);

    MatStack.matrixMode(tgt::MatrixStack::PROJECTION);
    MatStack.loadIdentity();
    MatStack.matrixMode(tgt::MatrixStack::MODELVIEW);
    MatStack.loadIdentity();

    glDisable(GL_BLEND);
    glBlendFunc(GL_ONE, GL_ZERO);
    glDisable(GL_LINE_SMOOTH);

    glPointSize(1.f);
    IMode.color(1.f,1.f,1.f,1.f);
}

void PlotLibraryOpenGl::logGlVertex2d(plot_t x, plot_t y) const {
    tgt::Vector2<plot_t> point = logScale2dtoLogCoordinates(x,y);
    IMode.vertex(point.x, point.y);
}

void PlotLibraryOpenGl::logGlVertex3d(plot_t x, plot_t y, plot_t z) const {
    tgt::Vector3<plot_t> point = logScale3dtoLogCoordinates(x,y,z);
    IMode.vertex(point.x, point.y, point.z);
}


namespace {
    /*
    // To avoid frequent con- and destruction of quadrics which might be very expensive
    // we instantiate here all quadrics which might be used to render glyphs locally.
    // This way each quadric exists only once during the whole lifetime of this class
    //
    // Hopefully it is ok to initialize them that early, but I haven't found any note
    // in the docs which said we need a valid OpenGL context.
    tgt::Disk        disk(0, 1, 16, 1, true, false);
    tgt::Triangle    triangle(1, true, false);
    tgt::Rect        rect(1, 1, true, false);
    tgt::Sphere      sphere(1, 12, 12, true, false);
    tgt::Tetrahedron tetrahedron(1, true, false);
    tgt::Quad        quad(1, 1, 1, true, false);
    */
}

void PlotLibraryOpenGl::renderGlyph(plot_t x, plot_t y, plot_t z, plot_t size) {
    if (glyphStyle_ == PlotEntitySettings::POINT) {
        glPointSize(static_cast<float>(size));

        tgt::Shader& shader = getSmoothPointShader();
#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE
        glEnable(GL_POINT_SPRITE);
#endif
        shader.activate();

        IMode.begin(tgt::ImmediateMode::POINTS);
            IMode.vertex(x,y,z);
        IMode.end();

        shader.deactivate();
#ifdef VRN_OPENGL_COMPATIBILITY_PROFILE
        glDisable(GL_POINT_SPRITE);
#endif
        return;
    }
    MatStack.pushMatrix();
    MatStack.translate(x,y,z);
    MatStack.scale(domain_[X_AXIS].size()/500, domain_[Y_AXIS].size()/500, domain_[Z_AXIS].size()/500);
    MatStack.scale(size, size, size);

    if (texture_ != 0 && !usePlotPickingManager_) {
        tgt::TextureUnit::setZeroUnit();
        texture_->bind();
        texture_->enable();
    }
    if (dimension_ == TWO){
        glDisable(GL_DEPTH_TEST);

        GlMeshGeometryUInt16TexCoord* mesh = nullptr;
        switch(glyphStyle_) {
            case PlotEntitySettings::CIRCLE:
                mesh = &diskMesh_;
                break;
            case PlotEntitySettings::TRIANGLE:
                mesh = &triangleMesh_;
                break;
            case PlotEntitySettings::QUAD:
                mesh = &rectMesh_;
                break;
            default:
                tgtAssert(false, "Invalid glyphStyle_");
                return;
        }
        if (texture_ != 0 && !usePlotPickingManager_) {
            mesh->render();
        } else {
            GlMeshGeometryUInt16Simple::renderDefault(*mesh);
        }

        glEnable(GL_DEPTH_TEST);
    }
    else if (dimension_ == THREE){
        GlMeshGeometryUInt16NormalTexCoord* mesh = nullptr;
        switch(glyphStyle_) {
            case PlotEntitySettings::CIRCLE:
                mesh = &sphereMesh_;
                break;
            case PlotEntitySettings::TRIANGLE:
                mesh = &tetrahedronMesh_;
                break;
            case PlotEntitySettings::QUAD:
                mesh = &quadMesh_;
                break;
            default:
                tgtAssert(false, "Invalid glyphStyle_");
                return;
        }
        if (texture_ != 0 && !usePlotPickingManager_) {
            mesh->render();
        } else {
            GlMeshGeometryUInt16Normal::renderDefault(*mesh);
        }
    }
    if (texture_ != 0 && !usePlotPickingManager_) {
        texture_->disable();
    }
    MatStack.popMatrix();
    LGL_ERROR;
}

void PlotLibraryOpenGl::renderSingleBar(plot_t left, plot_t right, plot_t bottom, plot_t top, tgt::Color c, plot_t squeeze) {
    plot_t back = domain_[Z_AXIS].getLeft();
    plot_t front = domain_[Z_AXIS].getRight();
    //squeeze
    left = (left+squeeze*left+right-squeeze*right)/2.0;
    right = (left-squeeze*left+right+squeeze*right)/2.0;
    back = (back+squeeze*back+front-squeeze*front)/2.0;
    front = (back-squeeze*back+front+squeeze*front)/2.0;

    ////
    // Draw filled quads
    ////
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonOffset(1.0, 1.0);

    if (!usePlotPickingManager_)
        IMode.color(0.8f*c.r, 0.8f*c.g, 0.8f*c.b, 1.0f);
    //right
    IMode.begin(tgt::ImmediateMode::QUADS);
        IMode.vertex(right, bottom, front); IMode.vertex(right, top, front);
        IMode.vertex(right, top, back);     IMode.vertex(right, bottom, back);
    IMode.end();

    //back
    IMode.begin(tgt::ImmediateMode::QUADS);
        IMode.vertex(right, bottom, back);  IMode.vertex(right, top, back);
        IMode.vertex(left, top, back);      IMode.vertex(left, bottom, back);
    IMode.end();

    //top
    if (!usePlotPickingManager_)
        IMode.color(0.9f*c.r, 0.9f*c.g, 0.9f*c.b, 1.0f);
    IMode.begin(tgt::ImmediateMode::QUADS);
        IMode.vertex(left, top, front);     IMode.vertex(left, top, back);
        IMode.vertex(right, top, back);     IMode.vertex(right, top, front);
    IMode.end();

    //bottom
    IMode.begin(tgt::ImmediateMode::QUADS);
        IMode.vertex(left, bottom, front);  IMode.vertex(left, bottom, back);
        IMode.vertex(right, bottom, back);  IMode.vertex(right, bottom, front);
    IMode.end();

    if (!usePlotPickingManager_)
        IMode.color(c.r, c.g, c.b, 1.0f);
    //front
    IMode.begin(tgt::ImmediateMode::QUADS);
        IMode.vertex(left, bottom, front);  IMode.vertex(left, top, front);
        IMode.vertex(right, top, front);    IMode.vertex(right, bottom, front);
    IMode.end();
    //left
    IMode.begin(tgt::ImmediateMode::QUADS);
        IMode.vertex(left, bottom, front);  IMode.vertex(left, top, front);
        IMode.vertex(left, top, back);      IMode.vertex(left, bottom, back);
    IMode.end();
    glDisable(GL_POLYGON_OFFSET_FILL);

    ////
    // Draw outlines
    ////
    if (!usePlotPickingManager_) {
        IMode.color(0.2f*c.r, 0.2f*c.g, 0.2f*c.b, 1.0f);

        //right
        IMode.begin(tgt::ImmediateMode::LINE_LOOP);
        IMode.vertex(right, bottom, front); IMode.vertex(right, top, front);
        IMode.vertex(right, top, back);     IMode.vertex(right, bottom, back);
        IMode.end();

        //back
        IMode.begin(tgt::ImmediateMode::LINE_LOOP);
        IMode.vertex(right, bottom, back);  IMode.vertex(right, top, back);
        IMode.vertex(left, top, back);      IMode.vertex(left, bottom, back);
        IMode.end();

        //top
        IMode.begin(tgt::ImmediateMode::LINE_LOOP);
        IMode.vertex(left, top, front);     IMode.vertex(left, top, back);
        IMode.vertex(right, top, back);     IMode.vertex(right, top, front);
        IMode.end();

        //bottom
        IMode.begin(tgt::ImmediateMode::LINE_LOOP);
        IMode.vertex(left, bottom, front);  IMode.vertex(left, bottom, back);
        IMode.vertex(right, bottom, back);  IMode.vertex(right, bottom, front);
        IMode.end();

        //front
        IMode.begin(tgt::ImmediateMode::LINE_LOOP);
        IMode.vertex(left, bottom, front);  IMode.vertex(left, top, front);
        IMode.vertex(right, top, front);    IMode.vertex(right, bottom, front);
        IMode.end();
        //left
        IMode.begin(tgt::ImmediateMode::LINE_LOOP);
        IMode.vertex(left, bottom, front);  IMode.vertex(left, top, front);
        IMode.vertex(left, top, back);      IMode.vertex(left, bottom, back);
        IMode.end();
        glDisable(GL_POLYGON_OFFSET_FILL);
    }



    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
}

tgt::Shader& PlotLibraryOpenGl::getSmoothPointShader() {
    if(!smoothPointShader_) {
        tgtAssert(ShdrMgr.isInited(), "ShaderManager not initialized");
        smoothPointShader_ = ShdrMgr.loadSeparate("pointsmoothing.vert", "", "pointsmoothing.frag", "", false);
        ShdrMgr.registerStaticShader(smoothPointShader_);
    }
    tgtAssert(smoothPointShader_, "Loading smooth point shader failed");
    return *smoothPointShader_;
}

tgt::dvec2 PlotLibraryOpenGl::convertPlotCoordinatesToViewport(const tgt::dvec3& plotCoordinates) const {
    tgt::dvec3 copy = logScale3dtoLogCoordinates(plotCoordinates);

    GLdouble x, y, z;
    tgt::Matrix4d mv(transpose(MatStack.getModelViewMatrix()));
    tgt::Matrix4d pr(transpose(MatStack.getProjectionMatrix()));
    GLint viewport[4]; glGetIntegerv(GL_VIEWPORT, viewport);
    tgt::gluProject(copy.x, copy.y, copy.z, mv.elem, pr.elem, viewport, &x, &y, &z);
    return tgt::dvec2(x, y);
}

tgt::dvec3 PlotLibraryOpenGl::convertPlotCoordinatesToViewport3(const tgt::dvec3& plotCoordinates) const {
    tgt::dvec3 copy = logScale3dtoLogCoordinates(plotCoordinates);

    GLdouble x, y, z;
    tgt::Matrix4d mv(transpose(MatStack.getModelViewMatrix()));
    tgt::Matrix4d pr(transpose(MatStack.getProjectionMatrix()));
    GLint viewport[4]; glGetIntegerv(GL_VIEWPORT, viewport);
    tgt::gluProject(copy.x, copy.y, copy.z, mv.elem, pr.elem, viewport, &x, &y, &z);
    return tgt::dvec3(x, y, z);
}

tgt::dvec2 PlotLibraryOpenGl::convertViewportToPlotCoordinates(tgt::ivec2 viewCoord) const {
    GLdouble x, y, z;
    tgt::Matrix4d mv(transpose(MatStack.getModelViewMatrix()));
    tgt::Matrix4d pr(transpose(MatStack.getProjectionMatrix()));
    GLint viewport[4]; glGetIntegerv(GL_VIEWPORT, viewport);
    tgt::gluUnProject(static_cast<double>(viewCoord.x), static_cast<double>(viewCoord.y), 0.0, mv.elem, pr.elem, viewport, &x, &y, &z);
    if (logarithmicAxisFlags_[X_AXIS])
        x = convertFromLogCoordinates(x, X_AXIS);
    if (logarithmicAxisFlags_[Y_AXIS])
        y = convertFromLogCoordinates(y, Y_AXIS);
    return tgt::dvec2(x, y);
}


} // namespace
