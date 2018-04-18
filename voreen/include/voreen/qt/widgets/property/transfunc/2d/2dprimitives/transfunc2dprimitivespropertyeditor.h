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

#ifndef VRN_TRANSFUNC2DPRIMITIVESEDITOR_H
#define VRN_TRANSFUNC2DPRIMITIVESEDITOR_H

#include "voreen/qt/widgets/property/transfunc/2d/transfunc2dpropertyeditor.h"
#include "voreen/core/datastructures/transfunc/2d/2dprimitives/utils/transfuncprimitive.h"

#include "voreen/qt/widgets/property/transfunc/utils/doubleslider.h"
#include "voreen/qt/widgets/property/transfunc/utils/alphapicker.h"

class QCheckBox;

namespace voreen {
class TransFunc2DPrimitivesPropertyEditorCanvas;
class TransFunc2DPrimitivesProperty;

/**
 * TransFuncEditorIntensityGradient is an editor widget for 2D transfer functions.
 * It provides a canvas where the transfer function can be edited and the histogram
 * of the rendered dataset is shown. Furthermore it provides several control elements
 * to adjust the appearance of the histogram and the primitives.
 */
class TransFunc2DPrimitivesPropertyEditor : public TransFunc2DPropertyEditor {
    Q_OBJECT
public:
    /** Constructor
     * @param prop the transfer function property
     * @param parent the parent widget
     */
    TransFunc2DPrimitivesPropertyEditor(TransFunc2DPrimitivesProperty* prop, QWidget* parent = 0);
    /** Destructor */
    ~TransFunc2DPrimitivesPropertyEditor();

    /** @override TransFuncPropertyEditorBase */
    virtual void updateFromProperty();
    /** @override TransFunc2DPropertyEditor */
    virtual std::string getXGuiName() {return "Intensity";}
    /** @override TransFunc2DPropertyEditor */
    virtual std::string getYGuiName() {return "Gradient";}
protected:
    /** @override TransFunc2DPropertyEditor */
    virtual void updateVolumeDataBoundsFromProperty();
    //----------------------
    //  Layout functions
    //----------------------
protected:
    /**
     * Adds the alpha picker to the color selection
     * @override TransFuncpropertyEditorBase
     */
    virtual QGroupBox* createColorPickerBox();
    /** @override TransFuncpropertyEditorBase */
    virtual QWidget* layoutLeftComponents();

private:
    AlphaPicker* alphaPicker_;                                    ///< picker to select the alpha value of a primitive
    TransFunc2DPrimitivesPropertyEditorCanvas* mappingCanvas_;    ///< mapping canvas
    DoubleSlider* thresholdSliderX_;                              ///< 2 slider for adjusting the thresholds
    DoubleSlider* thresholdSliderY_;                              ///< 2 slider for adjusting the thresholds

    QCheckBox* showHistogramCB_;    ///< toggles the histogram
    QSlider*  histogramBrightnessSlider_; ///< multiplies a factor to the alpha value
    QCheckBox* showTextureCB_;      ///< toggles the background texture

    QCheckBox* computeHistogram_;   ///< auto computation of histogram?

    QToolButton* addTriangleButton_;///< adds a triangle
    QToolButton* addQuadButton_;    ///< adds a quad
    QToolButton* addBananaButton_;  ///< adds a banana

    QSlider* primitiveFuzzinessSlider_; ///< manipulates the fuzziness of the selected primitive
    QToolButton* removePrimitiveButton_; ///< removes the selected primitive

    //----------------------
    //  Layout slots
    //----------------------
protected slots:
    /** Applies threshold changes to the tf */
    void sliderThresholdChanged();
    /** Toggles the texture in the mapping canvas. */
    void showHistogramToggled(int state);
    /** Toggles the histogram in the mapping canvas. */
    void showTextureToggled(int state);
    /** Adds a triangle to the TF. */
    void addTriangleButtonClicked();
    /** Adds a quad to the TF. */
    void addQuadButtonClicked();
    /** Adds a banana to the TF. */
    void addBananaButtonClicked();
    /** Changes the fuzziness of the selected priitive. */
    void fuzzinessChangedSlot(int value);

    void computeHistogramToggled(int state);

signals:
    void newPrimitiveAddedSignal(TransFuncPrimitive* prim);

    //----------------------
    //  Helper Functions
    //----------------------
protected slots:
    /**
     * This slot is called whenever the color of a key was changed.
     * The new color is propagated to the mapping canvas.
     *
     * @param h new hue of the key
     * @param s new saturation of the key
     * @param v new intensity of the key
     * @param a new alpha of the key
     */
    void markerColorChanged(int h, int s, int v, int a);

    //----------------------
    // Member
    //----------------------
private:
    TransFunc2DPrimitivesProperty* transFunc2DPrimitivesProp_;
    static const std::string loggerCat_;  ///< the logger category






    /**
     * Creates the whole layout for this widget.
     */
//    void createWidgets();

    /**
     * Creates the connections for all control elements.
     */
//    void createConnections();

    /**
     * Tests whether the volume associated with the transFuncProperty has changed and updates
     * the editor to the new volume when necessary.
     */
//    void updateFromProperty();

    /**
     * Propagates the new volume to the painter for the mapping canvas and updates the control
     * elements for the histogram.
     */
//    void volumeChanged();

    /**
     * Tests whether the value inside the transFuncProp is the same as transFuncGradient_.
     * A new transfer function object is created if they do not match. The transfer function
     * is reset if they match. Furthermore a repaint for the volume rendering is triggered.
     */
//    void resetEditor();

public slots:
    /**
     * Opens a filedialog and loads the selected transfer function.
     */
    //void loadTransferFunction();

    /**
     * Opens a filedialog and saves the current transfer function to the desired file.
     */
    //void saveTransferFunction();

    /**
     * Enables or disables the grid in the background.
     */
    //void toggleShowGrid();

    /**
     * Enables or disables the display of the histogram texture.
     */
    //void toggleShowHistogram();

    /**
     * Slot that is called when a primitive was selected. It enables the transparency
     * and fuzzines slider.
     */
    //void primitiveSelected();

    /**
     * Slot that is called when a primitive was deselected. The transparency and fuzzines
     * slider are disabled.
     */
    //void primitiveDeselected();

    /**
     * Sets the transparency slider to given value.
     *
     * @param trans new value for transparency slider
     */
    //void setTransparency(int trans);

    /**
     * Activates the coarseness mode.
     */
    //void startTracking();

    /**
     * Deactivates the coarseness mode.
     */
    //void stopTracking();

    /**
     * Informs the transfer function property about a change in transfer function
     * and emits transferFunctionChanged() signal.
     */
    //void repaintSignal();

    //void domainChanged();

    //void fitToDomain();

private:

    /**
     * Creates the layout with mapping canvas.
     *
     * @return the layout with mapping canvas
     */
    //QLayout* createMappingLayout();

    /**
     * Creates the layout with all standard buttons.
     *
     * @return the layout with all standard buttons
     */
    //QLayout* createButtonLayout();

    /**
     * Creates the layout with all buttons for manipulation of primitives.
     */
    //QLayout* createPrimitivesButtonLayout();

    /**
     * Creates the layout with transparency and fuzzines slider and control elements for the histogram.
     *
     * @return layout with 3 sliders and a checkbox
     */
    //QLayout* createSliderLayout();

  /*  TransFunc2DPrimitives* transFuncGradient_; ///< the 2D transfer function
    tgt::QtCanvas* transCanvas_;                    ///< canvas in which the tf can be edited
    TransFunc2DPrimitivesPropertyWidgetPainter* painter_;    ///< painter for the canvas

    QToolButton* loadButton_;             ///< button for loading of a 2D tf
    QToolButton* saveButton_;             ///< button for saving of a 2D tf
    QToolButton* clearButton_;            ///< button for resetting of the tf
    QToolButton* gridEnabledButton_;      ///< button for enabling and disabling the grid
    QToolButton* histogramEnabledButton_; ///< button for enabling and disabling the histogram
    QToolButton* fitToDomainButton_;      ///< button for enabling and disabling the histogram

    QToolButton* quadButton_;             ///< button for adding a quad primitive
    QToolButton* bananaButton_;           ///< button for adding a banana primitive
    QToolButton* deleteButton_;           ///< button for deletion of current selected primitive
    QToolButton* colorButton_;            ///< button for opening of a colordialog for the current selected primitive
    QCheckBox* histogramLog_;             ///< checkbox for enabling and disabling the logarithmic display of the histogram
    QSlider* histogramBrightness_;        ///< slider for adjusting the brightness of the histogram
    QSlider* fuzziness_;                  ///< slider for adjusting the fuzzines of a primitive
    QSlider* transparency_;               ///< slider for adjusting the transparency of a primitive
    QDoubleSpinBox* lowerIntensityDomainSpin_;         ///< spinbox for lower mapping value
    QDoubleSpinBox* upperIntensityDomainSpin_;         ///< spinbox for upper mapping value
    QDoubleSpinBox* lowerGradientDomainSpin_;          ///< spinbox for lower mapping value
    QDoubleSpinBox* upperGradientDomainSpin_;          ///< spinbox for upper mapping value
    Qt::Orientation orientation_;         ///< orientation of the widget*/


};

} // namespace voreen

#endif // VRN_TRANSFUNC2DPRIMITIVESEDITOR_H
