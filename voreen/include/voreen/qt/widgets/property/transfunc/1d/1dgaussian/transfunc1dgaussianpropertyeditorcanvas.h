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

#ifndef VRN_TRANSFUNC1DGAUSSIANPROPERTYEDITORCANVAS_H
#define VRN_TRANSFUNC1DGAUSSIANPROPERTYEDITORCANVAS_H

#include "voreen/core/properties/transfunc/1d/1dgaussian/transfunc1dgaussianproperty.h"

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
class TransFuncMappingCurve;
class TransFunc1DGaussian;

// ------------------------------------------------------------------------- //

/**
 * A widget that provides a canvas to edit the curves of an intensity transfer function.
 * The user can insert new curves by clicking at the desired location and reposition curves with holding
 * down left mouse button. Furthermore curves can be splitted, merged and deleted. The color of a curve
 * can also be changed.
 */
class VRN_QT_API TransFunc1DGaussianPropertyEditorCanvas : public QWidget {
    Q_OBJECT

public:

    /**
     * Constructor
     * @param parent the parent widget
     * @param tf the transfer function that is displayed in this widget
     */
    TransFunc1DGaussianPropertyEditorCanvas(QWidget* parent, TransFunc1DGaussianProperty* prop);
    /** Destructor */
    virtual ~TransFunc1DGaussianPropertyEditorCanvas();

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

    QMenu keyContextMenu_;   ///< context menu for right mouse click when a curve is selected
    QAction* unicolorAction_; ///< action for curve unicolor context menu entry
    QAction* multicolorAction_; ///< action for curve multicolor context menu entry
    QAction* splitMergeAction_; ///< action for split/merge context menu entry
    QAction* deActivateAction_; ///< action for curve (de-)activation context menu entry
    QAction* deleteAction_;     ///< action for delete curve context menu entry

public slots:
    /**
     * Changes the color of the selected key to the given value.
     * @note called from the editor/colorpicker
     * @param c new color of the selected key
     */
    void changeCurrentColorSlot(const QColor& c);

protected slots:
    /** Opens a colordialog for choosing the color of the current selected curve. */
    void colorChangeActionSlot();
    /** Sets the curves behavior to unicolor */
    void unicolorSlot();
    /** Sets the curves behavior to multicolor */
    void multicolorSlot();
    /** Splits or merges the current selected curve. */
    void splitMergeCurvesSlot();
    /** Deactivates/Reactivates the current selected curve */
    void deActivateCurveSlot();
    /** Deletes the current selected curve. */
    void deleteCurveSlot();

    //------------------
    //  Mouse events
    //------------------
public:
    /**
     * Coarseness mode is turned on, a new curve is inserted when no marker is at the mouse position
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
     * Moves the selected marker to new mouse position.
     * Nothing is done when no marker is selected.
     *
     * @param event the mouse event
     */
    virtual void mouseMoveEvent(QMouseEvent* event);

    /**
     * Opens a colordialog to change the color of the selected marker. Nothing happens when no marker
     * is selected.
     *
     * @param event the mouse event
     */
    virtual void mouseDoubleClickEvent(QMouseEvent* event);


    /**
     * Deletes the selected curve when del on keyboard was pressed.
     *
     * @param event the key event
     */
    virtual void keyReleaseEvent(QKeyEvent* event);

protected:
    //used during interaction
    QPoint mousePos_;                  ///< current position of the mouse

    TransFuncMappingCurve* selectedCurve_; ///< curve that was selected by the user
    CurvePart selectedPart_;               ///< selected key of the curve (BASE_LEFT/-RIGHT or PEAK_LEFT/-RIGHT)
    bool isKeyBeingDraged_;                ///< is the user dragging a key?

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
    // draw helper
    void drawCheckBoard(QPainter* painter);
    void drawAxes(QPainter* painter);
    void drawTransFunc(QPainter* painter);
    void drawHistogram(QPainter* painter);
    void drawMappingCurves(QPainter* painter);
    void drawThreshold(QPainter* painter);

    tgt::vec2 gridSpacing_;     ///< width and height of the underlying grid
    Histogram1D* histogram_;    ///< pointer to the histogram
    QPixmap* histogramCache_;   ///< image cache of the histogram
    tgt::vec2 visibleHistogramRange_; ///< visible histogram range
    bool showHistogram_;        ///< if true, the histogram is shown
    bool showTexture_;          ///< if true, the texture is shown

    //--------------------------
    //  Mapping Curve Functions
    //--------------------------
protected:
    // enum for the parts and attributes of a marker icon
    enum MarkerProps {
        MARKER_NORMAL   =  0, ///< marker is not selected and not split
        MARKER_LEFT     =  1, ///< left part of a key
        MARKER_RIGHT    =  2, ///< right part of a key
        MARKER_SELECTED =  4, ///< this explicit marker is selected
        MARKER_HIGHLIGHT = 8, ///< to highlight markers of selected curves
        MARKER_GRAYOUT  = 16  ///< to gray out markers of deactivated curves
    };

    /**
     * Creates a new curve at the given position.
     *
     * @param hit position of the key in relative coordinates
     */
    void insertNewCurve(tgt::vec2& hit);

    /**
     * Paints the three markers for the corresponding keys of the given curve.
     *
     * @param curve the curve whose markers will be drawn
     * @param paint the painter
     */
    void paintCurveMarkers(TransFuncMappingCurve* curve, QPainter& paint);

    /**
     * Draws a marker for a key of a curve of the transfer function.
     *
     * @param paint the painter
     * @param color the color of the key
     * @param p the position of the key in pixel coordinates
     * @param props the status of the key
     */
    void drawMarker(QPainter& paint, const tgt::col4& color, const tgt::vec2& p,
                    int props = 0);

    /**
    * Returns the property value for drawing the specified marker
    * (Peak, L/R-Bottom) of the curve.
    *
    * @param curve the curve which is drawn
    * @param marker which curve's key is drawn
    * @return the property value for drawing this marker
    */
    int getMarkerProps(TransFuncMappingCurve* curve, CurvePart marker);

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
    TransFunc1DGaussianProperty* tfProp_;           ///< pointer to the transfer function property
    static std::string loggerCat_;

signals:
    /**
     * Signal that is emitted when the user drags a key. It turns coarseness mode on or off.
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

#endif // VRN_TRANSFUNC1DGAUSSIANPROPERTYEDITORCANVAS_H
