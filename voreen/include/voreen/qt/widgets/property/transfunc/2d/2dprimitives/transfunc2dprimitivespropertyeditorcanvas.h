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

#ifndef VRN_TRANSFUNC2DPRIMITIVESPROPERTYEDITORCANVAS_H
#define VRN_TRANSFUNC2DPRIMITIVESPROPERTYEDITORCANVAS_H

#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"
#include "voreen/core/properties/transfunc/2d/2dprimitives/transfunc2dprimitivesproperty.h"

#include "voreen/qt/voreenqtapi.h"
#include "tgt/vector.h"

#include <QWidget>
#include <QMenu>
#include <QPixmap>

class QAction;
class QColor;
class QMouseEvent;

namespace voreen {

// Forward Declarations
class Histogram2D;
class TransFuncMappingKey;
class TransFunc1DKeys;

// ------------------------------------------------------------------------- //

/**
 * A widget that provides a canvas to edit the keys of an 2d primitive transfer function.
 */
class VRN_QT_API TransFunc2DPrimitivesPropertyEditorCanvas : public QWidget {
    Q_OBJECT
public:
    /**
     * Constructor
     * @param parent the parent widget
     * @param tf the transfer function that is displayed in this widget
     */
    TransFunc2DPrimitivesPropertyEditorCanvas(QWidget* parent, TransFunc2DPrimitivesProperty* prop);
    /** Destructor */
    virtual ~TransFunc2DPrimitivesPropertyEditorCanvas();

    /** Sets the histogram and all needed information. */
    void updateFromProperty();
    /** Sets the histogram. */
    void setHistogram(const Histogram2D* histogram);
    /** Sets the lower and upper threshold to the given values. */
    void setThreshold(tgt::vec2 thresholdX, tgt::vec2 thresholdY);

    //------------------
    //  Context Menu
    //------------------
protected:
    /** Initializes the key context menu. */
    void initializeKeyContextMenu();
    /** Initializes the key context menu. */
    void initializePrimitiveContextMenu();

    QMenu keyContextMenu_;   ///< context menu for right mouse click when a key is selected
    QMenu primitiveContextMenu_;   ///< context menu for right mouse click when a primitive is selected

public slots:
    /**
     * Changes the color of the selected key to the given value.
     * @note called from the editor/colorpicker
     * @param c new color of the selected key
     */
    void changeCurrentColorSlot(const QColor& c);
    /** Called from editor on histogramBrightnessSlider changes. */
    void setHistogramBrightness(int value);
    /**
     * Slot called to remove the selected primitive. It is called from the editor, context menu
     * and delete key. If no primitive or key is selected, the slot does nothing.
     */
    void deletePrimitiveSlot();
    /**
     * Changes the fuzziness of the selected primitive.
     * If no key or primitive is selected, this function does nothing
     */
    void updateFuzzinessSlot(int value);
    /**
     * Slot triggered, if the editor adds a new primitive to the tf.
     * @param prim the new primitive.
     */
    void newPrimitiveAddedSlot(TransFuncPrimitive* prim);

protected slots:
    /** Opens a colordialog for choosing the color of the current selected key. */
    void keyColorChangeActionSlot();
    /** Opens a colordialog for choosing the color of the current selected primitive. */
    void primitiveColorChangeActionSlot();

    //------------------
    //  Mouse events
    //------------------
public:
    /**
     * Coarseness mode is turned on, a new key is inserted when no key is at the mouse position
     * or the context menu is opened when right mousebutton was pressed.
     *
     * @param event the mouse event
     */
    virtual void mousePressEvent(QMouseEvent* event);

    /**
     * Switches coarseness mode off and hides the tooltip.
     *
     * @param event the mouse event
     */
    virtual void mouseReleaseEvent(QMouseEvent* event);

    /**
     * Moves the selected key to new mouse position or moves both keys at the
     * ends of the dragged line. Nothing is done when no key is selected nor a line is dragged.
     *
     * @param event the mouse event
     */
    virtual void mouseMoveEvent(QMouseEvent* event);

    /**
     * Opens a colordialog to change the color of the selected key. Nothing happens when no key
     * is selected.
     *
     * @param event the mouse event
     */
    virtual void mouseDoubleClickEvent(QMouseEvent* event);

    //virtual void keyPressEvent(QKeyEvent* event);

    /**
     * Unsets the cursor and deletes the selected key when del on keyboard was pressed.
     *
     * @param event the key event
     */
    virtual void keyReleaseEvent(QKeyEvent* event);

protected:
    //used during interaction
    QPoint lastMousePos_;                  ///< current position of the mouse

    TransFuncPrimitiveControlPoint* selectedKey_; ///< key that was selected by the user
    TransFuncPrimitive* selectedPrimitive_;       ///< primitive that was selected by the user

    bool isBeingDraged_;           ///< is the user dragging a key or primitive?

private:

    /**
     * Displays a tooltip at position pos with given values.
     *
     * @param pos position of the tooltip
     * @param values values that are displayed in the tooltip
     */
    void updateToolTipCoordinates(QPoint pos, tgt::vec2 values);

    //------------------
    //  Paint functions
    //------------------
protected:
    /**
     * Paints the current transfer function in a coordinate system
     * with a grid and a caption.
     *
     * @param event the paint event
     */
    virtual void paintEvent(QPaintEvent* event);
public:
    /** Called from editor. */
    void toggleHistogram(bool state);
    /** Called from editor. */
    void toggleTexture(bool state);
private:
    //helper
    void drawCheckBoard(QPainter* painter);
    void drawAxes(QPainter* painter);
    void drawTransFunc(QPainter* painter);
    void drawHistogram(QPainter* painter);
    void drawMappingKeys(QPainter* painter);
        void drawControlPoint(QPainter& paint, const tgt::col4& color, const tgt::vec2& p, bool selected = false);
    void drawThreshold(QPainter* painter);

    tgt::vec2 gridSpacing_;     ///< width and height of the underlying grid
    Histogram2D* histogram_;    ///< pointer to the histogram
    QPixmap* histogramCache_;   ///< image cache of the histogram
    tgt::vec2 visibleHistogramRangeX_; ///< visible histogram range
    tgt::vec2 visibleHistogramRangeY_; ///< visible histogram range
    bool showHistogram_;        ///< if true, the histogram is shown
    float histogramBrightness_; ///< factor for the histogram brightness
    bool showTexture_;          ///< if true, the texture is shown

    //------------------
    //  Helper Functions
    //------------------
private:
    /**
     * Helper function for calculation of pixel coordinates from relative coordinates.
     *
     * @param p relative coordinates in the interval [0,1]
     * @return pixel coordinates
     */
    tgt::vec2 normalizedToPixel(tgt::vec2 p);
    /**
     * Helper function for calculation of relative coordinates from pixel coordinates.
     *
     * @param p pixel coordinates in the interval [0,width]x[0,height]
     * @return relative coordinates
     */
    tgt::vec2 pixelToNormalized(tgt::vec2 p);

    //------------------
    //  Member
    //------------------
protected:
    TransFunc2DPrimitivesProperty* tfProp_;           ///< pointer to the transfer function property

signals:
    /**
     * Signal that is emitted when the user drags a key or a line. It turns coarseness mode on or off.
     * @param on should coarseness mode switched on or off?
     */
    void toggleInteractionModeSignal(bool on);
    /**
     * Signal that is emitted when the color of the selected key changed (via ColorDialog).
     * @param c new color of the selected key
     */
    void colorChangedSignal(const QColor& c);
    /**
     * Signal used to update the slider in the editor.
     * It is called, if a primitive is selected.
     * @param fuzziness it is between 0 and 100. -1 meens, no seleteced primitive.
     */
    void fuzzinessChangedSignal(int fuzziness);

};

} // namespace voreen

#endif // VRN_TRANSFUNC2DPRIMITIVESPROPERTYEDITORCANVAS_H
