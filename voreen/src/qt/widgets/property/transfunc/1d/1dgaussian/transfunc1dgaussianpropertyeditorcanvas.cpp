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

#include "voreen/qt/widgets/property/transfunc/1d/1dgaussian/transfunc1dgaussianpropertyeditorcanvas.h"

#include "voreen/core/datastructures/volume/volume.h"
#include "voreen/core/datastructures/volume/volumedisk.h"
#include "voreen/core/datastructures/volume/volumedecorator.h"
#include "voreen/core/datastructures/volume/histogram.h"
#include "voreen/core/datastructures/transfunc/1d/1dgaussian/transfunc1dgaussian.h"
#include "voreen/core/datastructures/transfunc/1d/1dgaussian/utils/transfuncmappingcurve.h"

#include <QAction>
#include <QApplication>
#include <QColor>
#include <QColorDialog>
#include <QMouseEvent>
#include <QPainter>
#include <QString>
#include <QToolTip>
#include <QImage>

#include <iostream>

namespace  {
    const QString X_AXIS_TEXT = QString("intensity");
    const QString Y_AXIS_TEXT = QString("opacity");

    const int AXIS_OFFSET = 12;
    const int ARROW_LENGTH = 10;
    const int ARROW_WIDTH = 3;
    // sizes of the marker parts
    const int KEY_POINT_SIZE = 10;
    const int KEY_HIGHLIGHT_SIZE = 12;
    const float KEY_SPLIT_FACTOR = 1.5f;
}

namespace voreen {

    std::string TransFunc1DGaussianPropertyEditorCanvas::loggerCat_ = "TF1DGaussPropEditorCanvas";

    TransFunc1DGaussianPropertyEditorCanvas::TransFunc1DGaussianPropertyEditorCanvas(QWidget* parent, TransFunc1DGaussianProperty* prop)
        : QWidget(parent)
        //context
        , splitMergeAction_(0), deActivateAction_(0), unicolorAction_(0), multicolorAction_(0), deleteAction_(0)
        //interaction
        , selectedCurve_(0), selectedPart_(NO_KEY), isKeyBeingDraged_(false)
        //paint
        , gridSpacing_(0.1f, 0.1f), histogram_(0), histogramCache_(0), visibleHistogramRange_(0.f, 1.f)
        , showHistogram_(true), showTexture_(true)
        //general
        , tfProp_(prop)
    {
        setObjectName("TransFunc1DGaussianPropertyEditorCanvas");
        setMouseTracking(true);
        setFocusPolicy(Qt::StrongFocus);
        setFocus();

        // init the context menu
        initializeKeyContextMenu();
    }

    TransFunc1DGaussianPropertyEditorCanvas::~TransFunc1DGaussianPropertyEditorCanvas() {
        delete histogramCache_;
    }

    void TransFunc1DGaussianPropertyEditorCanvas::updateFromProperty() {

        tgtAssert(tfProp_ && tfProp_->get(), "no tf or tf property!");
        if (const VolumeBase* volume = tfProp_->getVolume()) {
            //get histogram
            if (volume->hasDerivedData<VolumeHistogramIntensity>()) {
                Histogram1D* histogram = &volume->getDerivedData<VolumeHistogramIntensity>()->getHistogram(tfProp_->getVolumeChannel());
                if (histogram_ != histogram) {
                    histogram_ = histogram;
                    delete histogramCache_;
                    histogramCache_ = 0;
                }
            }
            else {
                if (tfProp_->getComputeHistogram())
                    volume->getDerivedDataThreaded<VolumeHistogramIntensity>();
                histogram_ = 0;
                delete histogramCache_;
                histogramCache_ = 0;
            }
        }
        else {
            histogram_ = 0;
            delete histogramCache_;
            histogramCache_ = 0;
        }
        if (visibleHistogramRange_ != tfProp_->get()->getDomain()) {
            visibleHistogramRange_ = tfProp_->get()->getDomain();
            delete histogramCache_;
            histogramCache_ = 0;
        }
        update();
    }

    //-------------------------------------------------------------------------------------------------------------
    //      Context menu
    //-------------------------------------------------------------------------------------------------------------
    void TransFunc1DGaussianPropertyEditorCanvas::initializeKeyContextMenu() {
        QAction* colorChangeAction = new QAction(tr("Change color of key"), this);
        colorChangeAction->setIcon(QIcon(":/qt/icons/colorize.png"));
        keyContextMenu_.addAction(colorChangeAction);
        connect(colorChangeAction, SIGNAL(triggered()), this, SLOT(colorChangeActionSlot()));

        keyContextMenu_.addSeparator();

        QActionGroup* colorModeGroup = new QActionGroup(this);
        colorModeGroup->setExclusive(true);

        unicolorAction_ = new QAction(tr("unicolor"), this);
        unicolorAction_->setCheckable(true);
        connect(unicolorAction_, SIGNAL(triggered()), this, SLOT(unicolorSlot()));
        colorModeGroup->addAction(unicolorAction_);


        multicolorAction_ = new QAction(tr("multicolor"), this);
        multicolorAction_->setCheckable(true);
        connect(multicolorAction_, SIGNAL(triggered()), this, SLOT(multicolorSlot()));
        colorModeGroup->addAction(multicolorAction_);

        keyContextMenu_.addActions(colorModeGroup->actions());

        keyContextMenu_.addSeparator();

        splitMergeAction_ = new QAction(tr("split this curve"), this);
        splitMergeAction_->setCheckable(true);
        keyContextMenu_.addAction(splitMergeAction_);
        connect(splitMergeAction_, SIGNAL(triggered()), this, SLOT(splitMergeCurvesSlot()));

        keyContextMenu_.addSeparator();

        deActivateAction_ = new QAction(tr("deactivate this curve"), this);
        deActivateAction_->setCheckable(true);
        keyContextMenu_.addAction(deActivateAction_);
        connect(deActivateAction_, SIGNAL(triggered()), this, SLOT(deActivateCurveSlot()));

        keyContextMenu_.addSeparator();

        deleteAction_ = new QAction(tr("Delete this curve"), this);
        deleteAction_->setIcon(QIcon(":/qt/icons/eraser.png"));
        keyContextMenu_.addAction(deleteAction_);
        connect(deleteAction_, SIGNAL(triggered()), this, SLOT(deleteCurveSlot()));
    }

    void TransFunc1DGaussianPropertyEditorCanvas::showKeyContextMenu(QMouseEvent* event) {
        tgtAssert(tfProp_ && tfProp_->get(), "no tf or tf property!");

        // enable correct uni-/multicolor item if a curve is selected
        unicolorAction_->setChecked(selectedCurve_ && selectedCurve_->isUnicolor());
        multicolorAction_->setChecked(selectedCurve_ && !selectedCurve_->isUnicolor());

        // show a check mark if the curve is already split
        splitMergeAction_->setChecked(selectedCurve_ && selectedCurve_->isSplit());

        // show a check mark if the curve is already deactivated
        deActivateAction_->setChecked(selectedCurve_ && !selectedCurve_->isActive());

        // allow deletion of curves only if there are at least two curves
        deleteAction_->setEnabled(tfProp_->get()->getNumCurves() > 1);
        deleteAction_->setChecked(selectedCurve_ && !selectedCurve_->isActive());

        keyContextMenu_.popup(event->globalPos());
    }

    void TransFunc1DGaussianPropertyEditorCanvas::colorChangeActionSlot() {
        if (!selectedCurve_ || selectedPart_ == NO_KEY)
            return;

        QColor oldColor;
        // start with the color of the selected key
        switch (selectedPart_){
        case PEAK_LEFT:
            oldColor = Col2QColor(selectedCurve_->getColorL());
            break;
        case PEAK_RIGHT:
            oldColor = Col2QColor(selectedCurve_->getColorR());
            break;
        case BASE_LEFT:
            oldColor = Col2QColor(selectedCurve_->getBaseColorL());
            break;
        case BASE_RIGHT:
            oldColor = Col2QColor(selectedCurve_->getBaseColorR());
            break;
        default:
            tgtAssert(false, "no key of the curve is selected");
        }

        QColor newColor = QColorDialog::getColor(oldColor, 0);
        if (newColor.isValid())
            changeCurrentColorSlot(newColor);
    }

    void TransFunc1DGaussianPropertyEditorCanvas::changeCurrentColorSlot(const QColor& c) {
        if (!selectedCurve_ || !c.isValid())
            return;

        tgt::col4 tgtcolor = QColor2Col(c);
        bool changedColor = false;


        // left peak key selected
        if (selectedPart_ == PEAK_LEFT) {
            tgtcolor.a = selectedCurve_->getColorL().a;
            if (selectedCurve_->getColorL() != tgtcolor) {
                selectedCurve_->setColorL(tgtcolor);
                changedColor = true;
            }
        }
        // right peak key selected
        else if (selectedPart_ == PEAK_RIGHT) {
            tgtcolor.a = selectedCurve_->getColorR().a;
            if (selectedCurve_->getColorR() != tgtcolor) {
                selectedCurve_->setColorR(tgtcolor);
                changedColor = true;
            }
        }
        // left base key selected
        else if (selectedPart_ == BASE_LEFT) {
            tgtcolor.a = selectedCurve_->getBaseColorL().a;
            if (selectedCurve_->getBaseColorL() != tgtcolor) {
                selectedCurve_->setBaseColorL(tgtcolor);
                changedColor = true;
            }
        }
        // right key marker selected
        else if (selectedPart_ == BASE_RIGHT) {
            tgtcolor.a = selectedCurve_->getBaseColorR().a;
            if (selectedCurve_->getBaseColorR() != tgtcolor) {
                selectedCurve_->setBaseColorR(tgtcolor);
                changedColor = true;
            }
        }
        else {
            tgtAssert(false, "no key selected");
        }

        if (changedColor) {
            tfProp_->get()->updateCurve(selectedCurve_);
            tfProp_->invalidate();
            emit colorChangedSignal(c);
        }
    }

    void TransFunc1DGaussianPropertyEditorCanvas::unicolorSlot() {
        if (!selectedCurve_)
            return;

        selectedCurve_->setUnicolor(true, selectedPart_);

        // update the transfer function texture and invalidate the property
        tfProp_->get()->updateCurve(selectedCurve_);
        tfProp_->invalidate();
    }

    void TransFunc1DGaussianPropertyEditorCanvas::multicolorSlot() {
        if (!selectedCurve_)
            return;

        selectedCurve_->setUnicolor(false);

        // update the transfer function texture and invalidate the property
        tfProp_->get()->updateCurve(selectedCurve_);
        tfProp_->invalidate();
    }

    void TransFunc1DGaussianPropertyEditorCanvas::splitMergeCurvesSlot() {
        if (!selectedCurve_)
            return;

        // use the left side for merging if a marker on the left side is selected
        bool useLeft = (selectedPart_ == BASE_LEFT) || (selectedPart_ == PEAK_LEFT);
        selectedCurve_->setSplit(splitMergeAction_->isChecked(), useLeft);

        // update the transfer function texture and invalidate the property
        tfProp_->get()->updateCurve(selectedCurve_);
        tfProp_->invalidate();
    }

    void TransFunc1DGaussianPropertyEditorCanvas::deActivateCurveSlot() {
        if (!selectedCurve_)
            return;

        selectedCurve_->setActive(!deActivateAction_->isChecked());

        // update the transfer function texture and invalidate the property
        tfProp_->get()->updateCurve(selectedCurve_);
        tfProp_->invalidate();
    }

    void TransFunc1DGaussianPropertyEditorCanvas::deleteCurveSlot() {
        tgtAssert(tfProp_ && tfProp_->get(), "no tf or tf property!");

        if (!selectedCurve_ || tfProp_->get()->getNumCurves() < 2)
            return;

        tfProp_->get()->removeCurve(selectedCurve_);
        selectedCurve_ = 0;
        emit colorChangedSignal(QColor(0, 0, 0, 0));

        tfProp_->invalidate();
    }

    //-------------------------------------------------------------------------------------------------------------
    //      Mouse events and other
    //-------------------------------------------------------------------------------------------------------------
    void TransFunc1DGaussianPropertyEditorCanvas::mousePressEvent(QMouseEvent* event) {
        //activate interaction mode
        if (event->button() == Qt::LeftButton)
            emit toggleInteractionModeSignal(true);

        event->accept();

        tgt::vec2 sHit = tgt::vec2(event->x(), static_cast<float>(height()) - event->y());
        tgt::vec2 hit = stow(sHit);

        // see if a curve was selected at one of its three markers
        selectedCurve_ = 0;
        selectedPart_ = NO_KEY;
        for (int i = 0; i<tfProp_->get()->getNumCurves(); ++i) {
            TransFuncMappingCurve* curve = tfProp_->get()->getCurve(i);

            // a peak marker was clicked
            tgt::vec2 sp = wtos(tgt::vec2(curve->getIntensity(), curve->getAlphaL()));
            tgt::vec2 spr = wtos(tgt::vec2(curve->getIntensity(), curve->getAlphaR()));
            if (curve->isSplit()) {
                if (sHit.x > sp.x - KEY_SPLIT_FACTOR * KEY_POINT_SIZE && sHit.x <= sp.x &&
                    sHit.y > sp.y - KEY_POINT_SIZE && sHit.y < sp.y + KEY_POINT_SIZE)
                {
                    selectedCurve_ = curve;
                    selectedPart_ = PEAK_LEFT;
                }
                if (sHit.x >= spr.x && sHit.x < spr.x + KEY_SPLIT_FACTOR * KEY_POINT_SIZE &&
                    sHit.y > spr.y - KEY_POINT_SIZE && sHit.y < spr.y + KEY_POINT_SIZE)
                {
                    selectedCurve_ = curve;
                    selectedPart_ = PEAK_RIGHT;
                }
            }
            else {
                if (sHit.x > sp.x - KEY_POINT_SIZE && sHit.x < sp.x + KEY_POINT_SIZE &&
                    sHit.y > sp.y - KEY_POINT_SIZE && sHit.y < sp.y + KEY_POINT_SIZE)
                {
                    selectedCurve_ = curve;
                    // in a non-splitted curve PEAK_LEFT and PEAK_RIGHT have equal meaning
                    selectedPart_ = PEAK_LEFT;
                }
            }

            // left base marker
            float baseX = 3 * sqrt(curve->getWidth());
            sp = wtos(tgt::vec2(curve->getIntensity() - baseX, curve->getOpacityAt(curve->getIntensity() - baseX)));
            if (sHit.x > sp.x - KEY_POINT_SIZE && sHit.x < sp.x + KEY_POINT_SIZE &&
                sHit.y > sp.y - KEY_POINT_SIZE && sHit.y < sp.y + KEY_POINT_SIZE)
            {
                selectedCurve_ = curve;
                selectedPart_ = BASE_LEFT;
            }

            // right base marker
            sp = wtos(tgt::vec2(curve->getIntensity() + baseX, curve->getOpacityAt(curve->getIntensity() + baseX)));
            if (sHit.x > sp.x - KEY_POINT_SIZE && sHit.x < sp.x + KEY_POINT_SIZE &&
                sHit.y > sp.y - KEY_POINT_SIZE && sHit.y < sp.y + KEY_POINT_SIZE)
            {
                selectedCurve_ = curve;
                selectedPart_ = BASE_RIGHT;
            }
        }

        // if a curve was selected, re-sort all curves so that this curve
        // will be drawn on top of the other curves
        if (selectedCurve_)
            tfProp_->get()->sortCurves(selectedCurve_);

        // if the right mouse button is pressed, open the context menu for the selected key
        if (event->button() == Qt::RightButton) {
            if (selectedCurve_) {
                showKeyContextMenu(event);

                // update the color shown in the color picker
                if (selectedPart_ == PEAK_LEFT) {
                    emit colorChangedSignal(Col2QColor(selectedCurve_->getColorL()));
                }
                else if (selectedPart_ == PEAK_RIGHT) {
                    emit colorChangedSignal(Col2QColor(selectedCurve_->getColorR()));
                }
                else if (selectedPart_ == BASE_LEFT) {
                    emit colorChangedSignal(Col2QColor(selectedCurve_->getBaseColorL()));
                }
                else if (selectedPart_ == BASE_RIGHT) {
                    emit colorChangedSignal(Col2QColor(selectedCurve_->getBaseColorR()));
                }
            }
            else {
                emit colorChangedSignal(QColor(0, 0, 0, 0));
            }
            update();
            return;
        }
        // if the left mouse button is pressed, the key can be dragged
        if (selectedCurve_ != 0 && event->button() == Qt::LeftButton) {
            isKeyBeingDraged_ = true;
            //keep values within valid range
            hit = tgt::clamp(hit, 0.f, 1.f);
            updateToolTipCoordinates(event->pos(), hit);

            // update the color shown in the color picker
            if (selectedPart_ == PEAK_LEFT) {
                emit colorChangedSignal(Col2QColor(selectedCurve_->getColorL()));
            }
            else if (selectedPart_ == PEAK_RIGHT) {
                emit colorChangedSignal(Col2QColor(selectedCurve_->getColorR()));
            }
            else if (selectedPart_ == BASE_LEFT) {
                emit colorChangedSignal(Col2QColor(selectedCurve_->getBaseColorL()));
            }
            else if (selectedPart_ == BASE_RIGHT) {
                emit colorChangedSignal(Col2QColor(selectedCurve_->getBaseColorR()));
            }
            update();
            return;
        }

        // if no curve was selected, insert a new curve
        if (hit.x >= 0.f && hit.x <= 1.f &&
            hit.y >= 0.f && hit.y <= 1.f &&
            event->button() == Qt::LeftButton)
        {
            insertNewCurve(hit); //calls invalidate
            isKeyBeingDraged_ = true;
            updateToolTipCoordinates(event->pos(), hit);
            update();
            emit colorChangedSignal(Col2QColor(selectedCurve_->getColorL()));
        }
    }

    void TransFunc1DGaussianPropertyEditorCanvas::mouseMoveEvent(QMouseEvent* event) {
        event->accept();
        mousePos_ = event->pos();

        tgt::vec2 sHit = tgt::vec2(event->x(), static_cast<float>(height()) - event->y());
        tgt::vec2 hit = stow(sHit);

        // return when no key is dragged
        if (!isKeyBeingDraged_)
            return;

        // keep location within valid texture coord range
        hit = tgt::clamp(hit, 0.f, 1.f);
        // if a base marker is moved, keep it on its correct side of the curve
        if ((selectedPart_ == BASE_LEFT && hit.x > selectedCurve_->getIntensity())
            || (selectedPart_ == BASE_RIGHT && hit.x < selectedCurve_->getIntensity()))
            hit.x = selectedCurve_->getIntensity();

        // marker gets dragged
        if (selectedCurve_ != 0) {
            updateToolTipCoordinates(event->pos(), hit);

            // a peak marker is dragged
            if (selectedPart_ == PEAK_LEFT || selectedPart_ == PEAK_RIGHT) {
                // x-position determines the curve's intentisy
                if (event->modifiers() != Qt::ShiftModifier) {
                    selectedCurve_->setIntensity(hit.x);
                }
                // y-position determines the intensity of the selected key
                if (event->modifiers() != Qt::ControlModifier) {
                    if (selectedPart_ == PEAK_LEFT)
                        selectedCurve_->setAlphaL(hit.y);
                    else
                        selectedCurve_->setAlphaR(hit.y);
                }
            }
            // a base marker, left or right, is dragged
            else if (selectedPart_ == BASE_LEFT || selectedPart_ == BASE_RIGHT) {

                // x-position determines the variance of the Gauss curve
                if (event->modifiers() != Qt::ShiftModifier) {
                    float diffX = (selectedCurve_->getIntensity() - hit.x) / 3.f;
                    if (diffX != 0)
                        selectedCurve_->setWidth(diffX * diffX);
                }

                // y-position determines the base value of the curve
                if (event->modifiers() != Qt::ControlModifier) {
                    selectedCurve_->setBaseValue(hit.y);
                }
            }

            tfProp_->get()->updateCurve(selectedCurve_);
            update();
            tfProp_->invalidate();
        }
    }

    void TransFunc1DGaussianPropertyEditorCanvas::mouseReleaseEvent(QMouseEvent* event) {
        event->accept();
        if (event->button() == Qt::LeftButton) {
            isKeyBeingDraged_ = false;
            QToolTip::hideText();
            update();
            emit toggleInteractionModeSignal(false);
        }
    }

    void TransFunc1DGaussianPropertyEditorCanvas::mouseDoubleClickEvent(QMouseEvent *event) {
        event->accept();
        if (event->button() == Qt::LeftButton)
            colorChangeActionSlot();
    }

    void TransFunc1DGaussianPropertyEditorCanvas::keyReleaseEvent(QKeyEvent* event) {
        if (event->key() == Qt::Key_Delete && selectedCurve_ != 0) {
            event->accept();
            deleteCurveSlot();
        }
    }

    void TransFunc1DGaussianPropertyEditorCanvas::updateToolTipCoordinates(QPoint pos, tgt::vec2 values) {
        std::ostringstream os;
        os.precision(2);
        os.setf(std::ios::fixed, std::ios::floatfield);

        float intensity = values.x;
        if (tfProp_->get()) {
            intensity = tfProp_->get()->getDomain().x + (tfProp_->get()->getDomain().y - tfProp_->get()->getDomain().x) * intensity;
        }
        os << intensity << " / " << values.y; // intensity / alpha
        QToolTip::showText(mapToGlobal(pos), QString(os.str().c_str()));
    }

    //-------------------------------------------------------------------------------------------------------------
    //      Paint Functions
    //-------------------------------------------------------------------------------------------------------------
    void TransFunc1DGaussianPropertyEditorCanvas::paintEvent(QPaintEvent* event) {
        tgtAssert(tfProp_ && tfProp_->get(), "no tf or property!");

        //the histogram is automatically painted onto this widget
        //we do not need to call the paintevent for the Histogrampainter directly
        event->accept();

        QPainter paint(this);

        // put origin in lower lefthand corner
        QMatrix m;
        m.translate(0.0, static_cast<float>(height()) - 1);
        m.scale(1.f, -1.f);
        paint.setMatrix(m);
        //draw white field
        paint.setMatrixEnabled(true);
        paint.setRenderHint(QPainter::Antialiasing, false);
        paint.setPen(Qt::NoPen);
        paint.setBrush(Qt::white);
        paint.drawRect(0, 0, width() - 1, height() - 1);

        //draw the grid
        drawCheckBoard(&paint);
        //draw transfunc texture
        if (showTexture_)
            drawTransFunc(&paint);
        //draw histogram
        if (showHistogram_)
            drawHistogram(&paint);
        //draw axis
        drawAxes(&paint);
        //draw threshold
        drawThreshold(&paint);
        //draw mapping curves
        drawMappingCurves(&paint);

        paint.setRenderHint(QPainter::Antialiasing, false);
        paint.setPen(Qt::lightGray);
        paint.setBrush(Qt::NoBrush);
        paint.drawRect(0, 0, width() - 1, height() - 1);
    }

    void TransFunc1DGaussianPropertyEditorCanvas::drawCheckBoard(QPainter* painter) {
        // draw grid
        painter->setPen(QColor(220, 220, 220));
        painter->setBrush(QColor(220, 220, 220));
        painter->setRenderHint(QPainter::Antialiasing, false);

        tgt::vec2 minCoord = wtos(tgt::vec2(0.f, 0.f));
        tgt::vec2 maxCoord = wtos(tgt::vec2(1.f, 1.f));

        if (showTexture_) {
            int cX = 0, cY = 0;
            tgt::vec2 size = wtos(tgt::vec2(gridSpacing_.x, gridSpacing_.y)) - minCoord;
            for (float x = 0.f; x < 1.f - gridSpacing_.x*0.5; x += gridSpacing_.x) {
                for (float y = 0.f; y < 1.f - gridSpacing_.y*0.5; y += gridSpacing_.y) {
                    if ((cX + cY) % 2 == 0) {
                        tgt::vec2 start = wtos(tgt::vec2(x, y));
                        tgt::vec2 tmp = wtos(tgt::vec2(x + gridSpacing_.x, y + gridSpacing_.y));
                        painter->drawRect(start.x, start.y, size.x, size.y);
                    }
                    cY++;
                }
                cY = 0;
                cX++;
            }
            painter->drawLine(maxCoord.x, minCoord.y, maxCoord.x, maxCoord.y);
            painter->drawLine(minCoord.x, maxCoord.y, maxCoord.x, maxCoord.y);
        }
        else {
            for (float x = gridSpacing_.x; x < 1.f + gridSpacing_.x*0.5; x += gridSpacing_.x) {
                tgt::vec2 tmp = wtos(tgt::vec2(x, 0.f));
                painter->drawLine(tmp.x, minCoord.y, tmp.x, maxCoord.y);
            }
            for (float y = gridSpacing_.y; y < 1.f + gridSpacing_.y*0.5; y += gridSpacing_.y) {
                tgt::vec2 tmp = wtos(tgt::vec2(0.f, y));
                painter->drawLine(minCoord.x, tmp.y, maxCoord.x, tmp.y);
            }
        }
    }

    void TransFunc1DGaussianPropertyEditorCanvas::drawTransFunc(QPainter* painter) {
        QPixmap icon(256, 2);
        icon.fill(Qt::transparent);
        //icon.set
        QPainter pmPainter(&icon);
        for (int i = 0; i < 256; i++) {
            tgt::Color color = tfProp_->get()->getTexture()->texelAsFloat(static_cast<size_t>((static_cast<float>(i) / 255.f)*(tfProp_->get()->getDimensions().x - 1)));
            pmPainter.setPen(QColor(color.r * 255, color.g * 255, color.b * 255, color.a * 255));
            pmPainter.drawLine(i, 0, i, 1);
        }
        tgt::vec2 start = wtos(tgt::vec2::zero);
        tgt::vec2 end = wtos(tgt::vec2::one);
        painter->drawPixmap(start.x, start.y, end.x - start.x, end.y - start.y, icon);
    }

    void TransFunc1DGaussianPropertyEditorCanvas::drawAxes(QPainter* painter) {
        // draw x and y axes
        painter->setRenderHint(QPainter::Antialiasing, true);
        painter->setPen(Qt::gray);
        painter->setBrush(Qt::gray);

        // draw axes independently from visible range
        tgt::vec2 origin = wtos(tgt::vec2(0.f, 0.f));
        origin.x = floor(origin.x) + 0.5f;
        origin.y = floor(origin.y) - 0.5f;

        painter->setRenderHint(QPainter::Antialiasing, true);

        painter->drawLine(QPointF(AXIS_OFFSET, origin.y),
            QPointF(width() - AXIS_OFFSET, origin.y));

        painter->drawLine(QPointF(origin.x, AXIS_OFFSET),
            QPointF(origin.x, height() - AXIS_OFFSET));

        QPointF arrow[3];
        arrow[0] = QPointF(origin.x, height() - AXIS_OFFSET);
        arrow[1] = QPointF(origin.x + ARROW_WIDTH, height() - AXIS_OFFSET - ARROW_LENGTH);
        arrow[2] = QPointF(origin.x - ARROW_WIDTH, height() - AXIS_OFFSET - ARROW_LENGTH);

        painter->drawConvexPolygon(arrow, 3);

        arrow[0] = QPointF(width() - AXIS_OFFSET, origin.y);
        arrow[1] = QPointF(width() - AXIS_OFFSET - ARROW_LENGTH, origin.y - ARROW_WIDTH);
        arrow[2] = QPointF(width() - AXIS_OFFSET - ARROW_LENGTH, origin.y + ARROW_WIDTH);

        painter->drawConvexPolygon(arrow, 3);

        painter->scale(-1.f, 1.f);
        painter->rotate(180.f);
        painter->drawText(static_cast<int>(width() - painter->fontMetrics().width(X_AXIS_TEXT) - 2.78f * AXIS_OFFSET), static_cast<int>(-1 * (origin.y + 1.0f - 0.8f * AXIS_OFFSET)), X_AXIS_TEXT);
        painter->drawText(static_cast<int>(1.6f * AXIS_OFFSET), static_cast<int>(-1 * (height() - 1.85f * AXIS_OFFSET)), Y_AXIS_TEXT);
        painter->rotate(180.f);
        painter->scale(-1.f, 1.f);
    }

    void TransFunc1DGaussianPropertyEditorCanvas::drawHistogram(QPainter* painter) {
        if (!histogram_) {
            painter->setMatrixEnabled(false);
            painter->setPen(Qt::red);
            painter->drawText(QRectF(0, 7, width() - 1, height() - 8), tr("No volume or calculating histogram"), QTextOption(Qt::AlignHCenter));
            painter->setMatrixEnabled(true);
        }
        else {
            if (histogramCache_ == 0 || histogramCache_->rect() != rect()) {
                delete histogramCache_;
                histogramCache_ = new QPixmap(rect().size());
                histogramCache_->fill(Qt::transparent);

                QPainter paint(histogramCache_);

                if (histogram_) {
                    // draw histogram
                    paint.setPen(Qt::NoPen);
                    paint.setBrush(QColor(255, 135, 135, 255)); //200 0 0 120
                    paint.setRenderHint(QPainter::Antialiasing, true);

                    int histogramWidth = static_cast<int>(histogram_->getNumBuckets());
                    tgt::vec2 p;

                    // Qt can't handle polygons that have more than 65536 points
                    // so we have to split the polygon
                    int maxSize = 65536; //max size of polygon
                    std::vector<QPointF*> points;
                    int vi = 0; //iterator in the points vector
                    points.push_back(new QPointF[maxSize]);
                    int count = 0;

                    for (int x = 0; x < histogramWidth; ++x) {
                        float xpos = static_cast<float>(x) / histogramWidth;
                        xpos = histogram_->getMinValue() + (histogram_->getMaxValue() - histogram_->getMinValue()) * xpos;
                        // Do some simple clipping here, as the automatic clipping of drawPolygon()
                        // gets very slow if lots of polygons have to be clipped away, e.g. when
                        // zooming to small part of the histogram.
                        if (xpos >= visibleHistogramRange_.x && xpos <= visibleHistogramRange_.y) {
                            //open new list, if old one is full
                            if (count == maxSize - 2){
                                count = 0;
                                points.push_back(new QPointF[maxSize]);
                                vi++;
                                //copy last point to connect two polygons
                                points[vi][count].rx() = p.x;
                                points[vi][count].ry() = p.y;
                                count++;
                            }
                            float value = (true ? histogram_->getBucketLogNormalized(x) : histogram_->getBucketNormalized(x));
                            p = histoWtos(tgt::vec2(xpos, value));

                            // optimization: if the y-coord has not changed from the two last points
                            // then just update the last point's x-coord to the current one
                            if ((count >= 2) && (points[vi][count - 2].ry() == p.y) && (points[vi][count - 1].ry() == p.y) && (count >= 2)){
                                points[vi][count - 1].rx() = p.x;
                            }
                            else {
                                points[vi][count].rx() = p.x;
                                points[vi][count].ry() = p.y;
                                count++;
                            }
                        }
                    }

                    for (size_t i = 0; i < points.size(); ++i){
                        if (count > 0) {
                            if (i == vi){
                                // needed for a closed polygon
                                p = histoWtos(tgt::vec2(0.f, 0.f));
                                points[i][count].rx() = points[i][count - 1].rx();
                                points[i][count].ry() = p.y;
                                count++;
                                points[i][count].rx() = points[i][0].rx();
                                points[i][count].ry() = p.y;
                                count++;

                                paint.drawPolygon(points[i], count);
                            }
                            else {
                                // needed for a closed polygon
                                p = histoWtos(tgt::vec2(0.f, 0.f));
                                points[i][maxSize - 2].rx() = points[i][maxSize - 3].rx();
                                points[i][maxSize - 2].ry() = p.y;
                                points[i][maxSize - 1].rx() = points[i][0].rx();
                                points[i][maxSize - 1].ry() = p.y;

                                paint.drawPolygon(points[i], maxSize);
                            }
                        }
                    }

                    for (size_t i = 0; i < points.size(); ++i)
                        delete[] points[i];
                }
            }
            painter->drawPixmap(0, 0, *histogramCache_);
        }
    }

    void TransFunc1DGaussianPropertyEditorCanvas::drawMappingCurves(QPainter* painter) {

        // draw the outline of all Gauss functions represented by the mapping curves
        QPen pen = QPen(Qt::darkRed);

        for (int i = 0; i < tfProp_->get()->getNumCurves(); ++i) {

            TransFuncMappingCurve *curve = tfProp_->get()->getCurve(i);

            // deactivated curves are greyed out
            if (!curve->isActive()) {
                pen.setColor(Qt::darkGray);
                pen.setWidthF(1.5f);
            }
            // a selected curve is bigger than the rest
            else if (curve == selectedCurve_) {
                pen.setColor(Qt::red);
                pen.setWidthF(1.8f);
            }
            else {
                pen.setColor(Qt::darkRed);
                pen.setWidthF(1.5f);
            }
            painter->setPen(pen);

            const float RESOLUTION = 0.002f;

            // left from the curve
            float x = curve->getIntensity() - (curve->isSplit() ? 0.01f : 0.f);
            float y = curve->getOpacityAt(x);
            tgt::vec2 old = wtos(tgt::vec2(x, y));
            tgt::vec2 cur;
            while (x > 0) {
                x -= RESOLUTION;
                if (x < 0.f)
                    x = 0.f;
                y = curve->getOpacityAt(x);
                cur = wtos(tgt::vec2(x, y));

                if (y > 0.001f)
                    painter->drawLine(QPointF(old.x, old.y), QPointF(cur.x, cur.y));
                old = cur;
            }


            // right from the curve
            x = curve->getIntensity() + (curve->isSplit() ? 0.01f : 0.f);
            y = curve->getOpacityAt(x);
            old = wtos(tgt::vec2(x, y));
            while (x < 1.f) {
                x += RESOLUTION;
                if (x > 1.f)
                    x = 1.f;
                y = curve->getOpacityAt(x);
                cur = wtos(tgt::vec2(x, y));

                if (y > 0.001f)
                    painter->drawLine(QPointF(old.x, old.y), QPointF(cur.x, cur.y));
                old = cur;
            }

            // paint the markers for the curve's keys
            paintCurveMarkers(curve, *painter);
        }
    }

    void TransFunc1DGaussianPropertyEditorCanvas::drawThreshold(QPainter* painter) {
        // grey out threshold area
        tgt::vec2 origin = wtos(tgt::vec2(0.f));
        painter->setBrush(QBrush(QColor(192, 192, 192, 230), Qt::SolidPattern));
        painter->setPen(Qt::NoPen);
        tgt::vec2 upperRight = wtos(tgt::vec2(1.f));
        tgt::vec2 lowerLeft = wtos(tgt::vec2(0.f));
        int w = static_cast<int>(upperRight.x - lowerLeft.x);
        int h = static_cast<int>(upperRight.y - lowerLeft.y);

        if (tfProp_->get()->getThreshold().x > 0.f) {
            painter->drawRect(static_cast<int>(origin.x), static_cast<int>(origin.y),
                static_cast<int>(tfProp_->get()->getThreshold().x * w + 1), h);
        }
        if (tfProp_->get()->getThreshold().y < 1.f) {
            painter->drawRect(static_cast<int>(origin.x + floor(tfProp_->get()->getThreshold().y * w)),
                static_cast<int>(origin.y), static_cast<int>((1 - tfProp_->get()->getThreshold().y) * w + 1), h);
        }
    }

    void TransFunc1DGaussianPropertyEditorCanvas::toggleHistogram(bool state) {
        showHistogram_ = state;
        update();
    }

    void TransFunc1DGaussianPropertyEditorCanvas::toggleTexture(bool state) {
        showTexture_ = state;
        update();
    }

    //-------------------------------------------------------------------------------------------------------------
    //      Mapping Curve Functions
    //-------------------------------------------------------------------------------------------------------------
    void TransFunc1DGaussianPropertyEditorCanvas::insertNewCurve(tgt::vec2& hit) {

        if (!tfProp_->get())
            return;

        hit = tgt::clamp(hit, 0.f, 1.f);

        TransFuncMappingCurve* curve = new TransFuncMappingCurve(hit.x, 0.01f, 0, QColor2Col(Qt::lightGray));

        tfProp_->get()->addCurve(curve);
        tgt::col4 curveColor = tgt::col4(255, 255, 255, 255);

        // set the peak as well as the base colors of the new curve
        curve->setColorL(curveColor);
        curve->setBaseColorL(curveColor);
        curve->setBaseColorR(curveColor);

        // overwrite alpha value with clicked position
        curve->setAlphaL(hit.y);

        // select the peak key of the new curve
        selectedCurve_ = curve;
        selectedPart_ = PEAK_LEFT;
        tfProp_->get()->sortCurves(selectedCurve_);

        // ivalidate the property and update the underlying transfer function texture
        tfProp_->get()->updateCurve(selectedCurve_);
        tfProp_->invalidate();
    }


    void TransFunc1DGaussianPropertyEditorCanvas::paintCurveMarkers(TransFuncMappingCurve* curve, QPainter& paint) {

        if (!curve)
            return;

        // properties of the currently drawn marker
        int props;

        // main marker at the peak of the curve -------------- //
        tgt::vec2 pos = wtos(tgt::vec2(curve->getIntensity(), curve->getColorL().a / 255.0));

        // the peak marker can be split: draw two markers
        if (curve->isSplit()) {
            // left part
            props = getMarkerProps(curve, PEAK_LEFT);
            drawMarker(paint, curve->getColorL(), pos, props);

            // right part
            pos = wtos(tgt::vec2(curve->getIntensity(), curve->getColorR().a / 255.0));
            props = getMarkerProps(curve, PEAK_RIGHT);
            drawMarker(paint, curve->getColorR(), pos, props);
        }
        // ..if it's not, draw a full marker
        else {
            // since the curve is not split, it is irrelevant which part of the peak key is passed
            props = getMarkerProps(curve, PEAK_LEFT);
            drawMarker(paint, curve->getColorL(), pos, props);
        }

        // base markers ------------------------------------ //
        // The base markers can leave the histogram area.
        // Therefore, its x coordinates have to be checked before drawing.

        // the base marker's positions are 3*sigma away from the peak
        float baseX = 3 * sqrt(curve->getWidth());

        // left base marker
        float x = curve->getIntensity() - baseX;
        if (x >= 0.f) {
            pos = wtos(tgt::vec2(x, curve->getOpacityAt(x)));
            props = getMarkerProps(curve, BASE_LEFT);
            drawMarker(paint, curve->getBaseColorL(), pos, props);
        }

        // right base marker
        x = curve->getIntensity() + baseX;
        if (x <= 1.f) {
            pos = wtos(tgt::vec2(x, curve->getOpacityAt(x)));
            props = getMarkerProps(curve, BASE_RIGHT);
            drawMarker(paint, curve->getBaseColorR(), pos, props);
        }
    }

    int TransFunc1DGaussianPropertyEditorCanvas::getMarkerProps(TransFuncMappingCurve* curve, CurvePart marker) {

        int props = 0;

        // selected curves and markers will be highlighted
        if (curve == selectedCurve_) {
            props |= MARKER_HIGHLIGHT;
            if (selectedPart_ == marker)
                props |= MARKER_SELECTED;
        }
        // markers of deactivated curves will be grayed out
        if (!curve->isActive())
            props |= MARKER_GRAYOUT;


        // the peak marker can be split..
        if ((marker == PEAK_LEFT || marker == PEAK_RIGHT) && curve->isSplit()) {
            if (marker == PEAK_LEFT)
                props |= MARKER_LEFT;
            else
                props |= MARKER_RIGHT;

            return props;
        }
        // ..if it's not a peak marker or not splitted, draw a full marker
        else {
            props |= MARKER_NORMAL;
            return props;
        }
    }

    void TransFunc1DGaussianPropertyEditorCanvas::drawMarker(QPainter& paint, const tgt::col4& tgtcolor, const tgt::vec2& p, int props) {

        // markers can be greyed out
        if (props & MARKER_GRAYOUT) {
            uint8_t sum = static_cast<uint8_t>(0.213f*tgtcolor.r + 0.715f*tgtcolor.g + 0.072f*tgtcolor.b);
            paint.setBrush(Col2QColor(tgt::col4(sum, sum, sum, 255)));
        }
        else {
            paint.setBrush(Col2QColor(tgtcolor));
        }

        QPen pen(QBrush(Qt::darkGray), Qt::SolidLine);

        // highlighted markers are bigger
        float size;
        if (props & MARKER_HIGHLIGHT) {
            size = KEY_HIGHLIGHT_SIZE;
            pen.setWidth(1);
        }
        else
            size = KEY_POINT_SIZE;

        // selected markers have a bigger outline
        if (props & MARKER_SELECTED)
            pen.setWidth(3);

        paint.setPen(pen);

        // draw the left part, right part or the whole marker
        if (props & MARKER_LEFT) {
            paint.drawPie(QRectF(p.x - KEY_SPLIT_FACTOR * size / 2, p.y - size / 2,
                KEY_SPLIT_FACTOR * size, size),
                90 * 16, 180 * 16);
        }
        else if (props & MARKER_RIGHT) {
            paint.drawPie(QRectF(p.x - KEY_SPLIT_FACTOR * size / 2, p.y - size / 2,
                KEY_SPLIT_FACTOR * size, size),
                270 * 16, 180 * 16);
        }
        else {
            paint.drawEllipse(QRectF(p.x - size / 2, p.y - size / 2,
                size, size));
        }
    }

    //-------------------------------------------------------------------------------------------------------------
    //      Helper Functions
    //-------------------------------------------------------------------------------------------------------------
    tgt::vec2 TransFunc1DGaussianPropertyEditorCanvas::wtos(tgt::vec2 p) {
        float sx = p.x * (static_cast<float>(width()) - 2.f * AXIS_OFFSET - 1.5f * ARROW_LENGTH) + AXIS_OFFSET;
        float sy = p.y * (static_cast<float>(height()) - 2.f * AXIS_OFFSET - 1.5f * ARROW_LENGTH) + AXIS_OFFSET;
        return tgt::vec2(sx, sy);
    }

    tgt::vec2 TransFunc1DGaussianPropertyEditorCanvas::histoWtos(tgt::vec2 p) {
        float sx = (p.x - visibleHistogramRange_.x) / (visibleHistogramRange_.y - visibleHistogramRange_.x) * (static_cast<float>(width()) - 2.f * AXIS_OFFSET - 1.5f * ARROW_LENGTH) + AXIS_OFFSET;
        float sy = p.y * (static_cast<float>(height()) - 2.f * AXIS_OFFSET - 1.5f * ARROW_LENGTH) + AXIS_OFFSET;
        return tgt::vec2(sx, sy);
    }

    tgt::vec2 TransFunc1DGaussianPropertyEditorCanvas::stow(tgt::vec2 p) {
        float wx = (p.x - AXIS_OFFSET) / (static_cast<float>(width()) - 2.f * AXIS_OFFSET - 1.5f * ARROW_LENGTH);
        float wy = (p.y - AXIS_OFFSET) / (static_cast<float>(height()) - 2.f * AXIS_OFFSET - 1.5f * ARROW_LENGTH);
        return tgt::vec2(wx, wy);
    }

    QColor TransFunc1DGaussianPropertyEditorCanvas::Col2QColor(const tgt::col4& color) {
        return QColor(color.r, color.g, color.b); // ignore alpha
    }

    tgt::col4 TransFunc1DGaussianPropertyEditorCanvas::QColor2Col(const QColor& color) {
        return tgt::col4(color.red(), color.green(), color.blue(), 255); // ignore alpha
    }

} // namespace voreen
