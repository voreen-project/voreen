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

#ifndef VRN_TRANSFUNC1DKEYSPROPERTYEDITORCANVAS_H
#define VRN_TRANSFUNC1DKEYSPROPERTYEDITORCANVAS_H

#include "voreen/core/properties/transfunc/1d/1dkeys/transfunc1dkeysproperty.h"

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
class Histogram1D;
class TransFuncMappingKey;
class TransFunc1DKeys;

// ------------------------------------------------------------------------- //

/**
 * A widget that provides a canvas to edit the keys of an intensity transfer function.
 * The user can insert new keys by clicking at the desired location and reposition keys with holding
 * down left mouse button. Furthermore keys can be splitted, merged and deleted. The color of a key
 * can also be changed.
 */
class VRN_QT_API TransFunc1DKeysPropertyEditorCanvas : public QWidget {
    Q_OBJECT
public:
    /**
     * Constructor
     * @param parent the parent widget
     * @param tf the transfer function that is displayed in this widget
     */
    TransFunc1DKeysPropertyEditorCanvas(QWidget* parent, TransFunc1DKeysProperty* prop);
    /** Destructor */
    virtual ~TransFunc1DKeysPropertyEditorCanvas();

    /** Sets the histogram and all needed information. */
    void updateFromProperty();
    /** Sets the histogram. */
    void setHistogram(const Histogram1D* histogram);
    /** Sets the lower and upper threshold to the given values. */
    void setThreshold(tgt::vec2 threshold);

    //------------------
    //  Context Menu
    //------------------
protected:
    /** Initializes the context menu. */
    void initializeKeyContextMenu();
    /**
     * Diplays the context menu at the given mouseposition
     * for the case of a keyselection.
     *
     * @param event the mouse event
     */
    void showKeyContextMenu(QMouseEvent* event);

    QMenu keyContextMenu_;   ///< context menu for right mouse click when a key is selected
    QAction* splitMergeAction_; ///< action for split/merge context menu entry
    QAction* zeroAction_;       ///< action for zero to right context menu entry
    QAction* deleteAction_;     ///< action for delete key context menu entry

public slots:
    /**
     * Changes the color of the selected key to the given value.
     * @note called from the editor/colorpicker
     * @param c new color of the selected key
     */
    void changeCurrentColorSlot(const QColor& c);

protected slots:
    /** Opens a colordialog for choosing the color of the current selected key. */
    void colorChangeActionSlot();
    /** Splits or merges the current selected key. */
    void splitMergeKeysSlot();
    /** Sets the left or right part of the current selected key to zero. */
    void zeroKeySlot();
    /** Deletes the current selected key. */
    void deleteKeySlot();

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

    /**
     * Sets the cursor to vertical size cursor when the mouse is on a line of the transfer function
     * and shift was pressed.
     *
     * @param event the key event
     */
    virtual void keyPressEvent(QKeyEvent* event);

    /**
     * Unsets the cursor and deletes the selected key when del on keyboard was pressed.
     *
     * @param event the key event
     */
    virtual void keyReleaseEvent(QKeyEvent* event);

protected:
    //used during interaction
    QPoint mousePos_;                  ///< current position of the mouse

    TransFuncMappingKey* selectedKey_; ///< key that was selected by the user
    bool isLeftPartSelected_;          ///< when selected key is split, was the left part selected?
    bool isKeyBeingDraged_;           ///< is the user dragging a key?

    int dragLine_;                     ///< number that indicates the line that was dragged using the shift modifier
    int dragLineStartY_;               ///< y position where the drag of the line started
    float dragLineAlphaLeft_;          ///< floating alpha value of the left key of the dragged line
    float dragLineAlphaRight_;         ///< floating alpha value of the right key of the dragged line
private:
    /**
     * Returns the number of the key that is left to the mouse position when
     * the position lies on the line between 2 keys. Otherwise -1 is returned.
     *
     * @param p mouse position in pixel coordinates
     * @return number of the nearest left key when on the line between two keys or -1 otherwise
     */
    int isLineHit(const tgt::vec2& pixelCoordinates);

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
    void drawThreshold(QPainter* painter);

    tgt::vec2 gridSpacing_;     ///< width and height of the underlying grid
    Histogram1D* histogram_;    ///< pointer to the histogram
    QPixmap* histogramCache_;   ///< image cache of the histogram
    tgt::vec2 visibleHistogramRange_; ///< visible histogram range
    bool showHistogram_;        ///< if true, the histogram is shown
    bool showTexture_;          ///< if true, the texture is shown
    //--------------------------
    //  Mapping Key Functions
    //--------------------------
protected:
    // enum for the status of a key
    enum MarkerProps {
        MARKER_NORMAL   =  0, ///< key is not selected and not split
        MARKER_LEFT     =  1, ///< left part of a key
        MARKER_RIGHT    =  2, ///< right part of a key
        MARKER_SELECTED =  4  ///< key is selected
    };

    /**
     * Creates a new key at the given position.
     *
     * @param hit position of the key in relative coordinates
     */
    void insertNewKey(tgt::vec2& hit);

    /**
     * Returns the nearest left or the nearest right key of the given key.
     * If no key exists at left or right 0 is returned.
     *
     * @param selectedKey key which nearest neighbour is searched
     * @param selectedLeftPart should the left neighbour be returned?
     * @return key that is the left or right neighbour of the given key. If no neighbour exists 0 is returned.
     */
    TransFuncMappingKey* getOtherKey(TransFuncMappingKey* selectedKey, bool selectedLeftPart);

    /**
     * Paints all keys of the transfer function.
     *
     * @param paint the painter
     */
    void paintKeys(QPainter& paint);

    /**
     * Draws the marker at the keys of the transfer function.
     *
     * @param paint the painter
     * @param color the color of the key
     * @param p the position of the key in pixel coordinates
     * @param props the status of the key
     */
    void drawMarker(QPainter& paint, const tgt::col4& color, const tgt::vec2& p,
                    int props = 0);

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
    tgt::vec2 wtos(tgt::vec2 p);
    /** wtos used for histogram calculation. */
    tgt::vec2 histoWtos(tgt::vec2 p);
    /**
     * Helper function for calculation of relative coordinates from pixel coordinates.
     *
     * @param p pixel coordinates in the interval [0,width]x[0,height]
     * @return relative coordinates
     */
    tgt::vec2 stow(tgt::vec2 p);
    /** Converts a tgt::col4 into a QColor. */
    QColor Col2QColor(const tgt::col4& color);
    /** Converts a QColor into a tgt::col4. */
    tgt::col4 QColor2Col(const QColor& color);

    //------------------
    //  Member
    //------------------
protected:
    TransFunc1DKeysProperty* tfProp_;           ///< pointer to the transfer function property

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

};

} // namespace voreen

#endif // VRN_TRANSFUNC1DKEYSPROPERTYEDITORCANVAS_H
