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

#ifndef VRN_QUADRUPLESLIDER_H
#define VRN_QUADRUPLESLIDER_H

#include <QWidget>
#include "voreen/qt/voreenqtapi.h"

class QMouseEvent;

namespace voreen {

/**
 * QuadrupleSlider is a QWidget providing four sliders on one ridge to define two ranges at a time.
 * One range (main) is dominant and must enclose the other range (extra).
 * The sliders can be moved by mouse interaction or by \a setMainMaxValue() and \a setMainMinValue()
 * method. The range of the values is from 0.0 to 1.0 and stored in float variables.
 * To avoid an overlay of the two main sliders the active one shunts the inactive slider.
 * Same for the extra sliders. Also the main sliders shunt the extra sliders to avoid an overlay.
 */
class VRN_QT_API QuadrupleSlider : public QWidget {
    Q_OBJECT
public:
    /**
     * Constructor - inits member variables to default
     *
     * @param parent the parent widget
     */
    QuadrupleSlider(QWidget* parent = 0);

    /**
     * The border to left and right are set to given values.
     *
     * @param left left border
     * @param right right border
     */
    void setOffsets(int left, int right);

    /**
     * Sets the minimal allowed distance between both sliders to given value.
     *
     * @param dist minimal allowed distance between both sliders
     */
    void setMinimalAllowedSliderDistance(float dist);

    /**
     * Returns the minimum of the range marked out by the two main sliders, e.g.
     * the x-coordinate of the left main slider in a normalized coordinate system.
     */
    float getMainMinValue();


    /**
     * Returns the maximum of the range marked out by the two main sliders, e.g.
     * the x-coordinate of the right main slider in a normalized coordinate system.
     */
    float getMainMaxValue();

    /**
     * Returns the minimum of the range marked out by the two extra sliders, e.g.
     * the x-coordinate of the left extra slider in a normalized coordinate system.
     */
    float getExtraMinValue();

    /**
     * Returns the maximum of the range marked out by the two extra sliders, e.g.
     * the x-coordinate of the right extra slider in a normalized coordinate system.
     */
    float getExtraMaxValue();

    /**
     * If extra sliders are disabled, they are not rendered or pickable.
     * The values can still be setted or getted and will emit siganls.
     */
    void setUseExtraSlider(bool b);

    /// Sets mapping to real values (only used for tooltip)
    void setMapping(float min, float max);
    void setUnit(std::string unit);
    void showToolTip(bool stt);

public slots:
     /**
      * Sets the main maximal value to max and emits maxValuesChanged(float, float).
      *
      * @param max new main maximal value
      */
     void setMainMaxValue(float max);

     /**
      * Sets the extra maximal value to max and emits extraValuesChanged(float, float).
      *
      * @param max new extra maximal value
      */
     void setExtraMaxValue(float max);

     /**
      * Sets the main minimal value to min and emits mainValuesChanged(float, float).
      *
      * @param min new main minimal value
      */
     void setMainMinValue(float min);

     /**
      * Sets the extra minimal value to min and emits extraValuesChanged(float, float).
      *
      * @param min new extra minimal value
      */
     void setExtraMinValue(float max);

     /**
      * Sets the main minimal and main maximal value. Emits mainValuesChanged(float, float).
      *
      * @param min new main minimal value
      * @param max new main maximal value
      */
     void setMainValues(float min, float max);

     /**
      * Sets the extra minimal and extra maximal value. Emits extraValuesChanged(float, float).
      *
      * @param min new extra minimal value
      * @param max new extra maximal value
      */
     void setExtraValues(float min, float max);

     /**
      * Sets the width of the main sliders. Default value is 9 pixels.
      *
      * @param sliderW width of main sliders
      */
     void setSliderMainWidth(int sliderW);

     /**
      * Sets the width of the extra sliders. Default value is 6 pixels.
      *
      * @param sliderW width of extra sliders
      */
     void setSliderExtraWidth(int sliderW);
signals:
     /**
      * This signal is emitted when the user drags a main slider.
      *
      * @param min the current main max value
      * @param max the current main min value
      */
     void mainValuesChanged(float min, float max);

     /**
      * This signal is emitted when the user drags an extra slider.
      *
      * @param min the current extra max value
      * @param max the current extra min value
      */
     void extraValuesChanged(float min, float max);

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
     * It will draw a horizontal line and four sliders.
     *
     * @param event the paint event emitted by Qt
     */
    void paintEvent(QPaintEvent* event);

    /**
     * Calculates the slider that is next to mouse position, sets it active
     * and calls \a moveSlider(float). Emits mainValuesChanged(float, float) and extraValuesChanged(float, float).
     *
     * @param event the mouse event
     */
    void mousePressEvent(QMouseEvent* event);

    /**
     * Moves the active slider by calling moveSlider() and emits
     * mainValuesChanged(float, float) and extraValuesChanged(float, float).
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

    std::string generateToolTipText(bool minMainBold, bool maxMainBold, bool minExtraBold, bool maxExtraBold);

    float getMappedValue(float norm);

private:
    /**
     * Internal Function to move the active slider(s). Called by \a mouseMoveEvent() and \a mousePressEvent().
     *
     * @param pos new slider position
     */
    void moveSlider(float pos);

    float minValueMain_;           ///< relative x coordinate of position of left main slider (e.g. TF domain)
    float maxValueMain_;           ///< relative x coordinate of position of right main slider (e.g. TF domain)
    float minValueExtra_;           ///< relative x coordinate of position of left extra slider (e.g. TF threshold)
    float maxValueExtra_;           ///< relative x coordinate of position of right extra slider (e.g. TF threshold)

    float normalizedMousePos_; ///< position where the drag of a slider started

    float mV1_;                ///< x coordinate of left main slider when dragging started
    float mV2_;                ///< x coordinate of right main slider when dragging started

    float minimalAllowedSliderDistance_; ///< mininmal allowed distance between both sliders
    int sliderMainWidth_;          ///< pixel width of main sliders
    int sliderExtraWidth_;         ///< pixel width of extra sliders

    int leftOffset_;           ///< offset to left before slider and line will start
    int rightOffset_;          ///< offset to right before slider and line will end

    bool leftMainSliderActive_;    ///< true when left main slider was clicked, false otherwise
    bool rightMainSliderActive_;   ///< true when right main slider was clicked, false otherwise
    bool leftExtraSliderActive_;    ///< true when left extra slider was clicked, false otherwise
    bool rightExtraSliderActive_;   ///< true when right extra slider was clicked, false otherwise

    bool useExtraSlider_;       ///< if false, the sliders are not rendered or selectable

    //Tooltip stuff:
    float minRWValue_;
    float maxRWValue_;
    std::string unit_;
    bool showToolTip_;
    QPoint globalMousePos_;
};

} // namespace voreen

#endif // VRN_QUADRUPLESLIDER_H
