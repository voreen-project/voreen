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

#ifndef VRN_TRANSFUNC1DPROPERTYWIDGETPAINTER_H
#define VRN_TRANSFUNC1DPROPERTYWIDGETPAINTER_H

#include "voreen/qt/widgets/property/transfunc/transfuncpropertywidgetpainterbase.h"

#include "voreen/qt/voreenqtapi.h"

namespace voreen {

class TransFunc1D;
class Histogram1D;

/**
 * The painter used in the TransFuncPropertyWidget to render the display canvas.
 */
class VRN_QT_API TransFunc1DPropertyWidgetPainter : public TransFuncPropertyWidgetPainterBase {
    Q_OBJECT
public:
    /**
     * @param canvas canvas that belongs to this painter
     */
    TransFunc1DPropertyWidgetPainter(tgt::GLCanvas* canvas);
    ~TransFunc1DPropertyWidgetPainter();

    /**
     * Paints the texture of the transfer function. A checkerboard is displayed
     * in the background.
     */
    void paint();

    /**
     * Initializes the painter, e.g. the projection and modelview matrix are set.
     */
    void initialize();

    /**
     * This method is called whenever the size of the widget changes.
     * It sets the viewport to the new size and updates the projection matrix.
     *
     * @param size new size of the widget
     */
    void sizeChanged(const tgt::ivec2& size);

    //--------------------
    //  getter and setter
    //--------------------
    /**
     * Sets the transfer function, which will be displayed.
     * Also controls the zoom level. If min = 1.f and max = 0.f the zoom will be reset.
     */
    void setTransFunc(TransFuncBase* tf);
    /** Sets the histogram, which will be displayed. */
    void setHistogram(const Histogram* histogram);
private:
    /** Help function used in paint. */
    void drawCheckBoard();
    /** Help function used in paint. */
    void drawTransferFunction();
    /** Help function used in paint. */
    void drawHistogram();
    /** Help function used in paint. */
    void drawThreshold(); //not used

    //--------------------
    //  mouse events
    //--------------------
    /** @see tgt::EventListener::mousePressEvent */
    void mousePressEvent(tgt::MouseEvent* e);
    /** @see tgt::EventListener::mouseMoveEvent */
    void mouseMoveEvent(tgt::MouseEvent* e);
    /** @see tgt::EventListener::mouseReleaseEvent */
    void mouseReleaseEvent(tgt::MouseEvent* e);
    /** @see tgt::EventListener::mouseDoubleClickEvent */
    void mouseDoubleClickEvent(tgt::MouseEvent* e);
    /** emits showInfoToolTip. Uses pressedSlider_ to create right tooltip. */
    void createInfoToolTip(QPoint mousePos);
protected:
    /** Enum used to communicate between mouse events. */
    enum MoveSlider {
        NO_SLIDER,
        GAMMA_SLIDER,
        DOMAIN_LEFT_SLIDER,
        DOMAIN_RIGHT_SLIDER,
        DOMAIN_BOTH_SLIDER
    } pressedSlider_;               ///< stores, which slider was clicked during mouse press.
    int mousePressedPosition_;      ///< position of the mouse press event. Used to move both sliders.
    int domainLeftPressedPosition_;        ///< position of the left domain slider during press
    int domainRightPressedPosition_;       ///< position of the right domain slider during press

protected slots:
    /** Slot connected to the gamma slider signal. Emits changed(). */
    void gammaSlot(float gamma);
    /** Slot connected to the domain slider signal. Emits changed(). */
    void domainSlot(tgt::vec2 domain);

    //--------------------
    //  zoom functions
    //--------------------
public:
    /** Zooms in by changing the minimal and maximal domain values. */
    void zoomIn();
    /** Zooms in by changing the minimal and maximal domain values. */
    void zoomOut();
    /** Fits min/maxDomainValue to tf_ domain. */
    void resetZoom();

    /** Returns the minimal domain value. It is the most left possible position of the domain slider. */
    float getMinVisibleDomainValue() const;
    /** Returns the minimal domain value. It is the most left possible position of the domain slider. */
    float getMaxVisibleDomainValue() const;

    /** Sets the visible array of the domain sliders. */
    void setVisibleDomainValues(tgt::vec2 visibleDomain);

private:
    //----------------
    //  Member
    //----------------
        //basic member
    TransFunc1D* tf_;                 ///< the displayed transfer function
    const Histogram1D* histogram_;  ///< the displayed histogram

        //slider
    TransFuncPropertyWidgetPainterGammaSlider* gammaSlider_; ///< slider to modify the gamma value
    TransFuncPropertyWidgetPainterDomainSlider* domainSlider_; ///< slider to modify the domain values

        //helper
    bool isInitialized_;            ///< painter must be initialized by calling initialize()

    bool renderHistogram_;          ///< if true, the histogram will be displayed
    bool renderGammaSlider_;        ///< if true, the gamma slider will be displayed
    bool renderDomainSlider_;       ///< if true, the domain slider will be displayed

    bool logarithmic_;              ///< if true, the histogram scale will be logarithmic

    tgt::vec2  visibleDomainValues_;  ///< most left value of the texture canvas
};

} // namespace voreen

#endif // VRN_TRANSFUNC1DPROPERTYWIDGETPAINTER_H
