/***********************************************************************************
 *                                                                                 *
 * Voreen - The Volume Rendering Engine                                            *
 *                                                                                 *
 * Copyright (C) 2005-2019 University of Muenster, Germany,                        *
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

#ifndef VRN_DOUBLESLIDER_H
#define VRN_DOUBLESLIDER_H

#include <QWidget>
#include "voreen/qt/voreenqtapi.h"

class QMouseEvent;

namespace voreen {

/**
 * DoubleSlider is a QWidget providing two sliders on one ridge to define a range.
 * The sliders can be moved by mouse interaction or by \a setMaxValue() and \a setMinValue()
 * method. The range of the values is from 0.0 to 1.0 and stored in float variables.
 * To avoid an overlay of the two sliders the active one shunts the inactive slider.
 */
class VRN_QT_API DoubleSlider : public QWidget {
    Q_OBJECT
public:
    /**
     * Constructor - inits member variables to default
     *
     * @param parent the parent widget
     * @param vertical defines if the slider should be vertical
     */
    DoubleSlider(QWidget* parent = 0, bool vertical = false);

    /**
     * The border to left,right or bottom,top are set to given values.
     *
     * @param leftTop left or bottom border
     * @param rightBottom right or top border
     */
    void setOffsets(int leftTop, int rightBottom);

    /**
     * Sets the minimal allowed distance between both sliders to given value.
     *
     * @param dist minimal allowed distance between both sliders
     */
    void setMinimalAllowedSliderDistance(float dist);

    /**
     * Returns the minimum of the range marked out by the two sliders, e.g.
     * the x-coordinate of the left slider in a normalized coordinate system.
     */
    float getMinValue();

    /**
     * Returns the maximum of the range marked out by the two sliders, e.g.
     * the x-coordinate of the right slider in a normalized coordinate system.
     */
    float getMaxValue();

    /// Sets mapping to real values (only used for tooltip)
    void setMapping(float min, float max);
    void setUnit(std::string unit);
    void showToolTip(bool stt);

public slots:
     /**
      * Sets the maximal value to max and emits valuesChanged(float, float).
      *
      * @param max new maximal value
      */
     void setMaxValue(float max);

     /**
      * Sets the minimal value to min and emits valuesChanged(float, float).
      *
      * @param min new minimal value
      */
     void setMinValue(float min);

     /**
      * Sets the minimal and maximal value. Emits valuesChanged(float, float).
      *
      * @param min new minimal value
      * @param max new maximal value
      */
     void setValues(float min, float max);

     /**
      * Sets the width of a slider. Default value is 5 pixels.
      *
      * @param sliderW width of a slider
      */
     void setSliderWidthOrHeight(int sliderWH);

signals:
     /**
      * This signal is emitted when the user drags a slider.
      *
      * @param min the current max value
      * @param max the current min value
      */
     void valuesChanged(float min, float max);

     /**
      * This signal is emitted when the user drags a slider. It will switch
      * the coarseness mode of volume rendering on and off.
      *
      * @param on should coarseness mode switched on or off?
      */
     void toggleInteractionMode(bool on);

protected:
    /**
     * This method is called by Qt whenever a repaint of the widget is necessary.
     * It will draw a horizontal line and two sliders.
     *
     * @param event the paint event emitted by Qt
     */
    void paintEvent(QPaintEvent* event);

    /**
     * Calculates the slider that is next to mouse position, sets it active
     * and calls \a moveSlider(float). Emits valuesChanged(float, float).
     *
     * @param event the mouse event
     */
    void mousePressEvent(QMouseEvent* event);

    /**
     * Moves the active slider by calling moveSlider() and emits
     * valuesChanged(float, float).
     *
     * @param event the mouse event
     */
    void mouseMoveEvent(QMouseEvent* event);

    /**
     * Only emits the toggleInteractionMode signal, so that coarseness mode is siwtched off.
     *
     * @param event the mouse event
     */
    void mouseReleaseEvent(QMouseEvent* event);

    /// Show a tooltip at the current mouse position
    void showToolTip(std::string text);

    std::string generateToolTipText(bool minBold, bool maxBold);

    float getMappedValue(float norm);

private:
    /**
     * Internal Function to move the active slider(s). Called by \a mouseMoveEvent() and \a mousePressEvent().
     *
     * @param pos new slider position
     */
    void moveSlider(float pos);

    float minValue_;           ///< relative x or y coordinate of position of left slider
    float maxValue_;           ///< relative x or y coordinate of position of right slider

    float normalizedMousePos_; ///< position where the drag of a slider started

    float mV1_;                ///< x or y coordinate of left slider when dragging started
    float mV2_;                ///< x or y coordinate of right slider when dragging started

    float minimalAllowedSliderDistance_;    ///< mininmal allowed distance between both sliders
    int sliderWidthHeight_;                 ///< pixel width or height of sliders

    int leftTopOffset_;          ///< offset to left or top before slider and line will start
    int rightBottomOffset_;            ///< offset to right or bottom before slider and line will end

    bool leftTopSliderActive_;   ///< true when left or top slider was clicked, false otherwise
    bool rightBottomSliderActive_;     ///< true when right or bottom slider was clicked, false otherwise

    bool vertical_;                 ///< defines the orientation
    //Tooltip stuff:
    float minRWValue_;
    float maxRWValue_;
    std::string unit_;
    bool showToolTip_;
    QPoint globalMousePos_;
};

} // namespace voreen

#endif // VRN_DOUBLESLIDER_H
