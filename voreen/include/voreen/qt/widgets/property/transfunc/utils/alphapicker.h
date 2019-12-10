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

#ifndef VRN_ALPHAPICKER_H
#define VRN_ALPHAPICKER_H

#include <QWidget>

class QColor;
class QMouseEvent;
class QPaintEvent;

namespace voreen {

/**
 * A widget that allows the user to selected the tranzparency for a given color.
 * A colorbar is diplayed that shows a colorgradient with all possible alpha values.
 * Furthermore there is an arrow that shows the current alpha. The user can adjust the
 * alpha by dragging the arrow or clicking on the colorbar.
 * @Note: The ColorPicker uses the HSV color space and ONLY takes the H+S component into account.
 *        Use the ColorLuminancePicker to adjust the V and the AlphaPicker to adjust an optional A.
 */
class AlphaPicker : public QWidget {
    Q_OBJECT

public:
    /** Constructor */
    AlphaPicker(QWidget* parent=0);
    /** Destructor */
    ~AlphaPicker();

    //---------------------------
    //  Interface Signals/Slots
    //---------------------------
public slots:
    /**
     * Slot that is called when the selected color in the colorpicker changes.
     * @param h new hue value
     * @param s new saturation value
     * @param v new value value
     */
    void updateHSVSlot(int h, int s, int v);

    /**
     * Slot that is called when the selected color in the colorpicker changes.
     * @param c selected color
     */
    void setHSVAByColorSlot(const QColor& c);

signals:
    /**
     * This signal is emitted when the alpha changed by the user.
     */
    void newHSVASignal(int h, int s, int v, int a);

    /**
     * This signal is emitted when the user changes the luminance. It will switch
     * the coarseness mode of volume rendering on and off.
     *
     * @param on should coarseness mode switched on or off?
     */
    void toggleInteractionModeSignal(bool b);

    //---------------------------
    //  Mouse and Paint Events
    //---------------------------
protected:
    /**
     * Paints the luminancepicker.
     *
     * @param event the paint event
     */
    void paintEvent(QPaintEvent* event);

    /**
     * Sets the arrow that indicates the luminance to the y-position of the mouse.
     *
     * @param event the mouse event
     */
    void mouseMoveEvent(QMouseEvent* event);

    /**
     * Sets the arrow that indicates the luminance to the y-position of the mouse.
     *
     * @param event the mouse event
     */
    void mousePressEvent(QMouseEvent* event);

    /**
     * Only emits the toggleInteractionMode signal, so that coarseness mode is siwtched off.
     *
     * @param event the mouse event
     */
    void mouseReleaseEvent(QMouseEvent* event);

    //---------------------------
    //  Convert Functions
    //---------------------------
private:
    /**
     * Converts the given height (y-value of mouse position) to alpha.
     * @param y y-position of mouse
     * @return alpha value
     */
    int yPositionToAlpha(int y);

    /**
     * Converts the alpha value to y-position of the mouse.
     * @param val value of hsv-colorspace
     * @return y-position of the mouse
     */
    int alphaToYPosition(int val);

    /**
     * This method sets the displayed color to the given values. The arrow that
     * indicates the alpha is set to a position that represents the value a.
     *
     * @param h hue value
     * @param s saturation value
     * @param v value in hsv-colorspace
     * @param a alpha value
     */
    void updateAllMembers(int h, int s, int v, int a);

    /**
     * Sets the alpha value to the given value.
     * @param a new alpha value
     */
    void updateAlpha(int a);

    //---------------------------
    //  Members
    //---------------------------
    static const int F_OFF_ = 0; ///< offet of the frame
    static const int C_OFF_ = 1; ///< offset of the content

    int hue_; ///< hue that is represented by the alphapicker.
    int sat_; ///< saturation that is represented by the alphapicker.
    int val_; ///< value in hsv-colorspace that is represented by the alphapicker.
    int alpha_; // alpha that is represented by the alphapicker.

    QPixmap* pix_; ///< pixmap that displays a bar with current choosen color and all possible luminance values.

};

 } // namespace voreen

#endif // VRN_COLORLUMINANCEPICKER_H
